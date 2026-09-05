#include <Rcpp.h>
#include <zlib.h>

#include "vendor/rad_core/include/kseq/kseq.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <unordered_map>
#include <vector>

KSEQ_INIT(gzFile, gzread)

namespace {

constexpr std::uint8_t kSequenceBits = 2;
constexpr std::uint8_t kQualityBits = 7;
constexpr std::uint64_t kFrontCodingRestart = 64;
constexpr std::int64_t kInteger64NA = std::numeric_limits<std::int64_t>::min();

SEXP packed_store_tag() {
  return Rf_install("rad::PackedReadStore");
}

SEXP fastq_stream_tag() {
  return Rf_install("rad::FastqStream");
}

std::uint32_t next_store_id() {
  static std::atomic<std::uint32_t> counter{1U};
  const std::uint32_t value = counter.fetch_add(1U, std::memory_order_relaxed);
  if (value == 0U || value > 0x7fffffffU) {
    throw std::runtime_error("RAD packed-store identifier space is exhausted");
  }
  return value;
}

class BitBuffer {
 public:
  void append(std::uint8_t value, std::uint8_t width) {
    if (width == 0 || width > 8) {
      throw std::invalid_argument("bit width must be between 1 and 8");
    }
    const std::uint64_t mask = (std::uint64_t{1} << width) - 1;
    value = static_cast<std::uint8_t>(value & mask);

    const std::uint64_t word_index = bit_size_ >> 6;
    const unsigned int word_offset = static_cast<unsigned int>(bit_size_ & 63U);
    if (word_index == words_.size()) {
      words_.push_back(0);
    }
    words_[word_index] |= static_cast<std::uint64_t>(value) << word_offset;
    if (word_offset + width > 64U) {
      words_.push_back(static_cast<std::uint64_t>(value) >> (64U - word_offset));
    }
    bit_size_ += width;
  }

  std::uint8_t get(std::uint64_t symbol_index, std::uint8_t width) const {
    const std::uint64_t bit_position = symbol_index * width;
    if (bit_position + width > bit_size_) {
      throw std::out_of_range("packed value index is out of range");
    }
    const std::uint64_t word_index = bit_position >> 6;
    const unsigned int word_offset = static_cast<unsigned int>(bit_position & 63U);
    std::uint64_t value = words_[word_index] >> word_offset;
    if (word_offset + width > 64U) {
      value |= words_[word_index + 1] << (64U - word_offset);
    }
    return static_cast<std::uint8_t>(value & ((std::uint64_t{1} << width) - 1));
  }

  std::uint64_t bit_size() const { return bit_size_; }
  std::uint64_t payload_bytes() const { return (bit_size_ + 7U) / 8U; }
  std::uint64_t allocated_bytes() const {
    return static_cast<std::uint64_t>(words_.capacity()) * sizeof(std::uint64_t);
  }

 private:
  std::vector<std::uint64_t> words_;
  std::uint64_t bit_size_ = 0;
};

void append_varint(std::vector<std::uint8_t>& output, std::uint64_t value) {
  while (value >= 0x80U) {
    output.push_back(static_cast<std::uint8_t>((value & 0x7fU) | 0x80U));
    value >>= 7U;
  }
  output.push_back(static_cast<std::uint8_t>(value));
}

std::uint64_t read_varint(const std::vector<std::uint8_t>& input, std::uint64_t& offset) {
  std::uint64_t value = 0;
  unsigned int shift = 0;
  while (offset < input.size() && shift <= 63U) {
    const std::uint8_t byte = input[offset++];
    value |= static_cast<std::uint64_t>(byte & 0x7fU) << shift;
    if ((byte & 0x80U) == 0) {
      return value;
    }
    shift += 7U;
  }
  throw std::runtime_error("corrupt variable-length integer in packed store");
}

bool is_lower_hex(char value) {
  return (value >= '0' && value <= '9') || (value >= 'a' && value <= 'f');
}

bool is_canonical_lower_uuid(const std::string& value) {
  if (value.size() != 36) {
    return false;
  }
  for (std::size_t i = 0; i < value.size(); ++i) {
    if (i == 8 || i == 13 || i == 18 || i == 23) {
      if (value[i] != '-') return false;
    } else if (!is_lower_hex(value[i])) {
      return false;
    }
  }
  return true;
}

std::uint8_t hex_value(char value) {
  if (value >= '0' && value <= '9') {
    return static_cast<std::uint8_t>(value - '0');
  }
  return static_cast<std::uint8_t>(value - 'a' + 10);
}

char hex_digit(std::uint8_t value) {
  static constexpr char digits[] = "0123456789abcdef";
  return digits[value & 0x0fU];
}

std::uint64_t stable_string_handle(const std::string& value) {
  // FNV-1a followed by a Murmur-style avalanche. The sign bit is cleared so
  // the payload remains a positive bit64 integer. Collisions are checked
  // against the original bytes while each store is built.
  std::uint64_t hash = 14695981039346656037ULL;
  for (unsigned char byte : value) {
    hash ^= static_cast<std::uint64_t>(byte);
    hash *= 1099511628211ULL;
  }
  hash ^= hash >> 33U;
  hash *= 0xff51afd7ed558ccdULL;
  hash ^= hash >> 33U;
  hash *= 0xc4ceb9fe1a85ec53ULL;
  hash ^= hash >> 33U;
  hash &= 0x7fffffffffffffffULL;
  return hash == 0U ? 1U : hash;
}

class PackedStringStore {
 public:
  enum class Codec : std::uint8_t { kNone = 0, kFront = 1, kUuidLower = 2 };

  void append(const std::string& value, bool present = true) {
    if (finalized_) {
      throw std::logic_error("cannot append to a finalized string store");
    }
    values_.push_back(present ? value : std::string());
    present_.push_back(present ? 1U : 0U);
    if (present) raw_bytes_ += value.size();
  }

  void finalize() {
    if (finalized_) return;
    const bool any_present = std::any_of(
      present_.begin(), present_.end(), [](std::uint8_t value) { return value != 0; });
    if (!any_present) {
      codec_ = Codec::kNone;
      values_.clear();
      values_.shrink_to_fit();
      finalized_ = true;
      return;
    }
    build_handles();
    bool uuid_compatible = !values_.empty();
    for (std::size_t i = 0; i < values_.size(); ++i) {
      if (!present_[i] || !is_canonical_lower_uuid(values_[i])) {
        uuid_compatible = false;
        break;
      }
    }

    if (uuid_compatible) {
      codec_ = Codec::kUuidLower;
      encode_uuids();
    } else {
      codec_ = Codec::kFront;
      encode_front();
    }
    values_.clear();
    values_.shrink_to_fit();
    finalized_ = true;
  }

  std::size_t size() const { return present_.size(); }
  bool present(std::size_t index) const {
    check_index(index);
    return present_[index] != 0;
  }

  std::string decode(std::size_t index) const {
    check_index(index);
    if (!present_[index]) return std::string();
    if (!finalized_) {
      throw std::logic_error("string store has not been finalized");
    }
    if (codec_ == Codec::kNone) {
      throw std::logic_error("cannot decode a missing packed string");
    }
    return codec_ == Codec::kUuidLower ? decode_uuid(index) : decode_front(index);
  }

  std::int64_t handle(std::size_t index) const {
    check_index(index);
    if (!present_[index]) return kInteger64NA;
    if (index >= handles_.size()) {
      throw std::logic_error("packed string handle index is unavailable");
    }
    return static_cast<std::int64_t>(handles_[index]);
  }

  std::string decode_handle(std::int64_t handle_value) const {
    if (handle_value == kInteger64NA) return std::string();
    if (handle_value <= 0) throw std::out_of_range("invalid packed string handle");
    const std::uint64_t target = static_cast<std::uint64_t>(handle_value);
    const auto found = std::lower_bound(sorted_handles_.begin(), sorted_handles_.end(), target);
    if (found == sorted_handles_.end() || *found != target) {
      throw std::out_of_range("packed string handle is absent from its backing RAD store");
    }
    const std::size_t position = static_cast<std::size_t>(found - sorted_handles_.begin());
    return decode(handle_indices_[position]);
  }

  const char* codec_name() const {
    if (codec_ == Codec::kNone) return "none";
    return codec_ == Codec::kUuidLower ? "uuid16" : "front";
  }
  std::uint64_t raw_bytes() const { return raw_bytes_; }
  std::uint64_t payload_bytes() const { return payload_.size(); }
  std::uint64_t index_bytes() const {
    return restart_offsets_.size() * sizeof(std::uint64_t) + present_.size() +
      sorted_handles_.size() * sizeof(std::uint64_t) +
      handle_indices_.size() * sizeof(std::uint32_t) +
      handles_.size() * sizeof(std::uint64_t);
  }

 private:
  void build_handles() {
    handles_.resize(values_.size(), 0U);
    std::unordered_map<std::uint64_t, std::size_t> unique_handles;
    for (std::size_t i = 0; i < values_.size(); ++i) {
      if (!present_[i]) continue;
      const std::uint64_t handle_value = stable_string_handle(values_[i]);
      const auto found = unique_handles.find(handle_value);
      if (found != unique_handles.end() && values_[found->second] != values_[i]) {
        throw std::runtime_error(
          "64-bit read-ID fingerprint collision; retain IDs as character data for this input");
      }
      if (found == unique_handles.end()) unique_handles[handle_value] = i;
      handles_[i] = handle_value;
    }
    std::vector<std::pair<std::uint64_t, std::uint32_t>> ordered;
    ordered.reserve(unique_handles.size());
    for (const auto& entry : unique_handles) {
      if (entry.second > std::numeric_limits<std::uint32_t>::max()) {
        throw std::overflow_error("packed string index exceeds 32-bit descriptor capacity");
      }
      ordered.emplace_back(entry.first, static_cast<std::uint32_t>(entry.second));
    }
    std::sort(ordered.begin(), ordered.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });
    sorted_handles_.reserve(ordered.size());
    handle_indices_.reserve(ordered.size());
    for (const auto& entry : ordered) {
      sorted_handles_.push_back(entry.first);
      handle_indices_.push_back(entry.second);
    }
  }

  void check_index(std::size_t index) const {
    if (index >= present_.size()) {
      throw std::out_of_range("string handle is out of range");
    }
  }

  void encode_uuids() {
    payload_.reserve(values_.size() * 16U);
    for (const std::string& value : values_) {
      std::uint8_t high = 0;
      bool have_high = false;
      for (char character : value) {
        if (character == '-') continue;
        if (!have_high) {
          high = static_cast<std::uint8_t>(hex_value(character) << 4U);
          have_high = true;
        } else {
          payload_.push_back(static_cast<std::uint8_t>(high | hex_value(character)));
          have_high = false;
        }
      }
    }
  }

  void encode_front() {
    std::string previous;
    for (std::size_t i = 0; i < values_.size(); ++i) {
      if (i % kFrontCodingRestart == 0) {
        restart_offsets_.push_back(payload_.size());
        previous.clear();
      }
      const std::string& current = values_[i];
      std::size_t common = 0;
      const std::size_t limit = std::min(previous.size(), current.size());
      while (common < limit && previous[common] == current[common]) ++common;
      append_varint(payload_, common);
      append_varint(payload_, current.size() - common);
      payload_.insert(payload_.end(), current.begin() + common, current.end());
      previous = current;
    }
  }

  std::string decode_uuid(std::size_t index) const {
    static constexpr std::array<int, 4> hyphen_positions{{8, 13, 18, 23}};
    const std::size_t offset = index * 16U;
    if (offset + 16U > payload_.size()) {
      throw std::runtime_error("corrupt UUID block in packed store");
    }
    std::string result;
    result.reserve(36);
    std::size_t nibble_index = 0;
    for (int output_index = 0; output_index < 36; ++output_index) {
      if (std::find(hyphen_positions.begin(), hyphen_positions.end(), output_index) !=
          hyphen_positions.end()) {
        result.push_back('-');
        continue;
      }
      const std::uint8_t byte = payload_[offset + nibble_index / 2U];
      const std::uint8_t nibble = nibble_index % 2U == 0 ? byte >> 4U : byte & 0x0fU;
      result.push_back(hex_digit(nibble));
      ++nibble_index;
    }
    return result;
  }

  std::string decode_front(std::size_t index) const {
    const std::size_t restart = index / kFrontCodingRestart;
    if (restart >= restart_offsets_.size()) {
      throw std::runtime_error("corrupt front-coded string index");
    }
    std::uint64_t offset = restart_offsets_[restart];
    const std::size_t start = restart * kFrontCodingRestart;
    std::string previous;
    for (std::size_t i = start; i <= index; ++i) {
      const std::uint64_t common = read_varint(payload_, offset);
      const std::uint64_t suffix_size = read_varint(payload_, offset);
      if (common > previous.size() || offset + suffix_size > payload_.size()) {
        throw std::runtime_error("corrupt front-coded string block");
      }
      std::string current = previous.substr(0, common);
      current.append(reinterpret_cast<const char*>(payload_.data() + offset), suffix_size);
      offset += suffix_size;
      previous.swap(current);
    }
    return previous;
  }

  Codec codec_ = Codec::kFront;
  bool finalized_ = false;
  std::uint64_t raw_bytes_ = 0;
  std::vector<std::string> values_;
  std::vector<std::uint8_t> present_;
  std::vector<std::uint8_t> payload_;
  std::vector<std::uint64_t> restart_offsets_;
  std::vector<std::uint64_t> handles_;
  std::vector<std::uint64_t> sorted_handles_;
  std::vector<std::uint32_t> handle_indices_;
};

std::uint8_t encode_base(char base, bool& exceptional) {
  exceptional = false;
  switch (base) {
    case 'A': return 0;
    case 'C': return 1;
    case 'T': return 2;
    case 'G': return 3;
    default:
      exceptional = true;
      return 0;
  }
}

char decode_base(std::uint8_t code) {
  switch (code & 3U) {
    case 0: return 'A';
    case 1: return 'C';
    case 2: return 'T';
    default: return 'G';
  }
}

class PackedReadStore {
 public:
  explicit PackedReadStore(int phred_offset)
      : phred_offset_(phred_offset), store_id_(next_store_id()) {
    if (phred_offset < 0 || phred_offset > 126) {
      throw std::invalid_argument("phred_offset must be between 0 and 126");
    }
    sequence_offsets_.push_back(0);
    sequence_exception_offsets_.push_back(0);
    quality_offsets_.push_back(0);
  }

  void append(const std::string& id,
              bool id_present,
              const std::string& comment,
              bool comment_present,
              const std::string& sequence,
              const std::string* quality) {
    if (finalized_) {
      throw std::logic_error("cannot append to a finalized read store");
    }
    ids_.append(id, id_present);
    comments_.append(comment, comment_present);

    std::uint64_t previous_exception_position = 0;
    bool first_exception = true;
    for (std::size_t position = 0; position < sequence.size(); ++position) {
      bool exceptional = false;
      sequence_bits_.append(encode_base(sequence[position], exceptional), kSequenceBits);
      if (exceptional) {
        const std::uint64_t delta = first_exception
          ? static_cast<std::uint64_t>(position)
          : static_cast<std::uint64_t>(position) - previous_exception_position;
        append_varint(sequence_exceptions_, delta);
        sequence_exceptions_.push_back(static_cast<std::uint8_t>(sequence[position]));
        previous_exception_position = position;
        first_exception = false;
      }
    }
    sequence_offsets_.push_back(sequence_offsets_.back() + sequence.size());
    sequence_exception_offsets_.push_back(sequence_exceptions_.size());
    raw_sequence_bytes_ += sequence.size();

    double sum = 0.0;
    int minimum = std::numeric_limits<int>::max();
    std::uint64_t at_least_q20 = 0;
    if (quality != nullptr) {
      if (quality->size() != sequence.size()) {
        throw std::invalid_argument("sequence and quality lengths differ");
      }
      for (unsigned char character : *quality) {
        if (character > 126U) {
          std::ostringstream message;
          message << "quality byte " << static_cast<int>(character)
                  << " is outside the supported FASTQ ASCII range";
          throw std::invalid_argument(message.str());
        }
        const int score = static_cast<int>(character) - phred_offset_;
        if (score < 0 || score > 93) {
          std::ostringstream message;
          message << "quality byte " << static_cast<int>(character)
                  << " is outside Q0-Q93 for Phred+" << phred_offset_;
          throw std::invalid_argument(message.str());
        }
        quality_bits_.append(static_cast<std::uint8_t>(score), kQualityBits);
        sum += score;
        minimum = std::min(minimum, score);
        if (score >= 20) ++at_least_q20;
      }
      quality_present_.push_back(1U);
      quality_offsets_.push_back(quality_offsets_.back() + quality->size());
      raw_quality_bytes_ += quality->size();
      if (quality->empty()) {
        mean_quality_.push_back(R_NaReal);
        minimum_quality_.push_back(NA_INTEGER);
        fraction_q20_.push_back(R_NaReal);
      } else {
        mean_quality_.push_back(sum / static_cast<double>(quality->size()));
        minimum_quality_.push_back(minimum);
        fraction_q20_.push_back(
          static_cast<double>(at_least_q20) / static_cast<double>(quality->size()));
      }
    } else {
      quality_present_.push_back(0U);
      quality_offsets_.push_back(quality_offsets_.back());
      mean_quality_.push_back(R_NaReal);
      minimum_quality_.push_back(NA_INTEGER);
      fraction_q20_.push_back(R_NaReal);
    }
  }

  void finalize() {
    if (finalized_) return;
    ids_.finalize();
    comments_.finalize();
    finalized_ = true;
  }

  std::size_t size() const { return quality_present_.size(); }
  std::uint32_t store_id() const { return store_id_; }
  std::int64_t handle(std::size_t index) const {
    check_index(index);
    if (index >= std::numeric_limits<std::uint32_t>::max()) {
      throw std::overflow_error("a packed RAD store cannot exceed 2^32 - 1 reads");
    }
    const std::uint64_t value = (static_cast<std::uint64_t>(store_id_) << 32U) |
      static_cast<std::uint64_t>(index + 1U);
    return static_cast<std::int64_t>(value);
  }
  bool has_id(std::size_t index) const { return ids_.present(index); }
  bool has_comment(std::size_t index) const { return comments_.present(index); }
  bool has_quality(std::size_t index) const {
    check_index(index);
    return quality_present_[index] != 0;
  }
  std::uint64_t sequence_length(std::size_t index) const {
    check_index(index);
    return sequence_offsets_[index + 1] - sequence_offsets_[index];
  }
  std::uint64_t quality_length(std::size_t index) const {
    check_index(index);
    return quality_offsets_[index + 1] - quality_offsets_[index];
  }
  std::string id(std::size_t index) const { return ids_.decode(index); }
  std::string comment(std::size_t index) const { return comments_.decode(index); }
  std::int64_t id_handle(std::size_t index, bool comment = false) const {
    return comment ? comments_.handle(index) : ids_.handle(index);
  }
  std::string string_from_handle(std::int64_t value, bool comment = false) const {
    return comment ? comments_.decode_handle(value) : ids_.decode_handle(value);
  }

  std::string sequence(std::size_t index) const {
    check_index(index);
    const std::uint64_t start = sequence_offsets_[index];
    const std::uint64_t length = sequence_length(index);
    std::string result(length, 'A');
    for (std::uint64_t i = 0; i < length; ++i) {
      result[i] = decode_base(sequence_bits_.get(start + i, kSequenceBits));
    }

    std::uint64_t exception_offset = sequence_exception_offsets_[index];
    const std::uint64_t exception_end = sequence_exception_offsets_[index + 1];
    std::uint64_t position = 0;
    bool first = true;
    while (exception_offset < exception_end) {
      const std::uint64_t delta = read_varint(sequence_exceptions_, exception_offset);
      position = first ? delta : position + delta;
      first = false;
      if (exception_offset >= exception_end || position >= length) {
        throw std::runtime_error("corrupt sequence exception lane");
      }
      result[position] = static_cast<char>(sequence_exceptions_[exception_offset++]);
    }
    return result;
  }

  std::vector<int> quality_values(std::size_t index) const {
    check_index(index);
    if (!has_quality(index)) return {};
    const std::uint64_t start = quality_offsets_[index];
    const std::uint64_t length = quality_length(index);
    std::vector<int> result(length);
    for (std::uint64_t i = 0; i < length; ++i) {
      result[i] = quality_bits_.get(start + i, kQualityBits);
    }
    return result;
  }

  std::string quality_ascii(std::size_t index, int output_offset) const {
    check_index(index);
    if (!has_quality(index)) return std::string();
    if (output_offset < 0 || output_offset > 126) {
      throw std::invalid_argument("quality output offset must be between 0 and 126");
    }
    const std::uint64_t start = quality_offsets_[index];
    const std::uint64_t length = quality_length(index);
    std::string result(length, '!');
    for (std::uint64_t i = 0; i < length; ++i) {
      const int value = quality_bits_.get(start + i, kQualityBits) + output_offset;
      if (value > 126) {
        throw std::invalid_argument("quality score cannot be represented with requested offset");
      }
      result[i] = static_cast<char>(value);
    }
    return result;
  }

  double mean_quality(std::size_t index) const {
    check_index(index);
    return mean_quality_[index];
  }
  int minimum_quality(std::size_t index) const {
    check_index(index);
    return minimum_quality_[index];
  }
  double fraction_at_least(std::size_t index, int threshold) const {
    check_index(index);
    if (!has_quality(index) || quality_length(index) == 0) return R_NaReal;
    const std::uint64_t start = quality_offsets_[index];
    const std::uint64_t length = quality_length(index);
    std::uint64_t passing = 0;
    for (std::uint64_t i = 0; i < length; ++i) {
      if (quality_bits_.get(start + i, kQualityBits) >= threshold) ++passing;
    }
    return static_cast<double>(passing) / static_cast<double>(length);
  }

  Rcpp::List statistics() const {
    const std::uint64_t index_bytes =
      sequence_offsets_.capacity() * sizeof(std::uint64_t) +
      sequence_exception_offsets_.capacity() * sizeof(std::uint64_t) +
      quality_offsets_.capacity() * sizeof(std::uint64_t) +
      quality_present_.capacity() +
      mean_quality_.capacity() * sizeof(double) +
      minimum_quality_.capacity() * sizeof(int) +
      fraction_q20_.capacity() * sizeof(double) +
      ids_.index_bytes() + comments_.index_bytes();

    return Rcpp::List::create(
      Rcpp::Named("reads") = static_cast<double>(size()),
      Rcpp::Named("store_id") = static_cast<double>(store_id_),
      Rcpp::Named("id_codec") = ids_.codec_name(),
      Rcpp::Named("id_raw_bytes") = static_cast<double>(ids_.raw_bytes()),
      Rcpp::Named("id_packed_bytes") = static_cast<double>(ids_.payload_bytes()),
      Rcpp::Named("comment_codec") = comments_.codec_name(),
      Rcpp::Named("comment_raw_bytes") = static_cast<double>(comments_.raw_bytes()),
      Rcpp::Named("comment_packed_bytes") = static_cast<double>(comments_.payload_bytes()),
      Rcpp::Named("sequence_bases") = static_cast<double>(sequence_offsets_.back()),
      Rcpp::Named("sequence_raw_bytes") = static_cast<double>(raw_sequence_bytes_),
      Rcpp::Named("sequence_packed_bytes") = static_cast<double>(sequence_bits_.payload_bytes()),
      Rcpp::Named("sequence_exception_bytes") = static_cast<double>(sequence_exceptions_.size()),
      Rcpp::Named("quality_scores") = static_cast<double>(quality_offsets_.back()),
      Rcpp::Named("quality_raw_bytes") = static_cast<double>(raw_quality_bytes_),
      Rcpp::Named("quality_packed_bytes") = static_cast<double>(quality_bits_.payload_bytes()),
      Rcpp::Named("index_allocated_bytes") = static_cast<double>(index_bytes),
      Rcpp::Named("sequence_buffer_allocated_bytes") =
        static_cast<double>(sequence_bits_.allocated_bytes()),
      Rcpp::Named("quality_buffer_allocated_bytes") =
        static_cast<double>(quality_bits_.allocated_bytes()),
      Rcpp::Named("phred_offset") = phred_offset_);
  }

  const std::vector<double>& means() const { return mean_quality_; }
  const std::vector<int>& minima() const { return minimum_quality_; }
  const std::vector<double>& fractions_q20() const { return fraction_q20_; }

 private:
  void check_index(std::size_t index) const {
    if (index >= size()) {
      throw std::out_of_range("packed read handle is out of range");
    }
  }

  int phred_offset_;
  std::uint32_t store_id_;
  bool finalized_ = false;
  PackedStringStore ids_;
  PackedStringStore comments_;
  BitBuffer sequence_bits_;
  BitBuffer quality_bits_;
  std::vector<std::uint64_t> sequence_offsets_;
  std::vector<std::uint64_t> quality_offsets_;
  std::vector<std::uint64_t> sequence_exception_offsets_;
  std::vector<std::uint8_t> sequence_exceptions_;
  std::vector<std::uint8_t> quality_present_;
  std::vector<double> mean_quality_;
  std::vector<int> minimum_quality_;
  std::vector<double> fraction_q20_;
  std::uint64_t raw_sequence_bytes_ = 0;
  std::uint64_t raw_quality_bytes_ = 0;
};

void set_integer64_value(Rcpp::NumericVector& output, R_xlen_t index, std::int64_t value) {
  static_assert(sizeof(double) == sizeof(std::int64_t), "integer64 transport requires 64-bit double");
  std::memcpy(output.begin() + index, &value, sizeof(value));
}

std::int64_t get_integer64_value(const Rcpp::NumericVector& input, R_xlen_t index) {
  std::int64_t value;
  std::memcpy(&value, input.begin() + index, sizeof(value));
  return value;
}

Rcpp::NumericVector make_handles(const Rcpp::XPtr<PackedReadStore>& store,
                                 const std::string& kind,
                                 const std::string& field = "") {
  const std::size_t count = store->size();
  Rcpp::NumericVector result(count);
  for (std::size_t i = 0; i < count; ++i) {
    bool present = true;
    if (kind == "rad_id") {
      present = field == "comment" ? store->has_comment(i) : store->has_id(i);
    } else if (kind == "rad_phred") {
      present = store->has_quality(i);
    }
    const std::int64_t value = kind == "rad_id"
      ? store->id_handle(i, field == "comment")
      : (present ? store->handle(i) : kInteger64NA);
    set_integer64_value(result, i, value);
  }
  result.attr("class") = Rcpp::CharacterVector::create(kind, "integer64");
  result.attr(".rad_store") = store;
  if (!field.empty()) result.attr(".rad_field") = field;
  return result;
}

Rcpp::NumericVector make_read_keys(std::size_t count) {
  Rcpp::NumericVector result(count);
  for (std::size_t i = 0; i < count; ++i) {
    set_integer64_value(result, i, static_cast<std::int64_t>(i + 1U));
  }
  result.attr("class") = "integer64";
  return result;
}

Rcpp::XPtr<PackedReadStore> store_from_handles(SEXP handles) {
  Rcpp::NumericVector vector(handles);
  if (!vector.hasAttribute(".rad_store")) {
    Rcpp::stop("packed handle has lost its backing RAD store");
  }
  SEXP pointer = vector.attr(".rad_store");
  if (TYPEOF(pointer) != EXTPTRSXP || R_ExternalPtrAddr(pointer) == nullptr) {
    Rcpp::stop("packed RAD store is no longer available");
  }
  if (R_ExternalPtrTag(pointer) != packed_store_tag()) {
    Rcpp::stop("external pointer is not a RAD packed-read store");
  }
  return Rcpp::XPtr<PackedReadStore>(pointer);
}

std::vector<std::size_t> indices_from_handles(SEXP handles,
                                              const PackedReadStore& store,
                                              bool allow_missing = true) {
  Rcpp::NumericVector vector(handles);
  std::vector<std::size_t> result(vector.size(), std::numeric_limits<std::size_t>::max());
  for (R_xlen_t i = 0; i < vector.size(); ++i) {
    const std::int64_t value = get_integer64_value(vector, i);
    if (value == kInteger64NA) {
      if (!allow_missing) Rcpp::stop("missing packed handle is not permitted here");
      continue;
    }
    if (value <= 0) {
      Rcpp::stop("packed handle points outside its backing RAD store");
    }
    const std::uint64_t bits = static_cast<std::uint64_t>(value);
    const std::uint32_t handle_store = static_cast<std::uint32_t>(bits >> 32U);
    const std::uint32_t ordinal = static_cast<std::uint32_t>(bits & 0xffffffffU);
    if (handle_store != store.store_id()) {
      Rcpp::stop("packed handles from different RAD stores were combined; use rad_rbind_reads()");
    }
    if (ordinal == 0U || static_cast<std::uint64_t>(ordinal) > store.size()) {
      Rcpp::stop("packed handle points outside its backing RAD store");
    }
    result[i] = static_cast<std::size_t>(ordinal - 1U);
  }
  return result;
}

Rcpp::DataFrame table_from_store(std::unique_ptr<PackedReadStore> store_owner,
                                 bool include_comment) {
  store_owner->finalize();
  if (store_owner->size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    Rcpp::stop("one rad_reads table cannot exceed R's compact row-name limit");
  }
  Rcpp::XPtr<PackedReadStore> store(store_owner.get(), true);
  R_SetExternalPtrTag(store, packed_store_tag());
  store_owner.release();
  const std::size_t count = store->size();
  Rcpp::NumericVector lengths(count);
  for (std::size_t i = 0; i < count; ++i) {
    lengths[i] = static_cast<double>(store->sequence_length(i));
  }
  Rcpp::NumericVector means = Rcpp::wrap(store->means());
  Rcpp::IntegerVector minima = Rcpp::wrap(store->minima());
  Rcpp::NumericVector fractions = Rcpp::wrap(store->fractions_q20());

  Rcpp::List columns;
  columns["read_key"] = make_read_keys(count);
  columns["seq_id"] = make_handles(store, "rad_id", "id");
  if (include_comment) {
    columns["seq_comment"] = make_handles(store, "rad_id", "comment");
  }
  columns["seq"] = make_handles(store, "rad_seq");
  columns["qual"] = make_handles(store, "rad_phred");
  columns["length"] = lengths;
  columns["mean_q"] = means;
  columns["min_q"] = minima;
  columns["frac_q20"] = fractions;
  columns.attr("class") = Rcpp::CharacterVector::create("data.frame");
  columns.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -static_cast<int>(count));
  columns.attr(".rad_store") = store;
  return Rcpp::DataFrame(columns);
}

class FastqStream {
 public:
  FastqStream(const std::string& path, int phred_offset, bool keep_comment)
      : path_(path), phred_offset_(phred_offset), keep_comment_(keep_comment) {
    file_ = gzopen(path.c_str(), "rb");
    if (file_ == nullptr) {
      throw std::runtime_error("could not open FASTQ file: " + path);
    }
    record_ = kseq_init(file_);
    if (record_ == nullptr) {
      gzclose(file_);
      file_ = nullptr;
      throw std::runtime_error("could not initialize FASTQ parser: " + path);
    }
  }

  ~FastqStream() { close(); }
  FastqStream(const FastqStream&) = delete;
  FastqStream& operator=(const FastqStream&) = delete;

  void close() {
    if (record_ != nullptr) {
      kseq_destroy(record_);
      record_ = nullptr;
    }
    if (file_ != nullptr) {
      gzclose(file_);
      file_ = nullptr;
    }
    closed_ = true;
  }

  bool closed() const { return closed_; }
  bool done() const { return eof_; }
  bool keep_comment() const { return keep_comment_; }

  std::unique_ptr<PackedReadStore> next(std::size_t count) {
    if (failed_) {
      throw std::runtime_error(
        "FASTQ stream cannot be resumed after a previous read failure");
    }
    if (closed_) throw std::runtime_error("FASTQ stream is closed");
    if (count == 0) throw std::invalid_argument("chunk size must be positive");
    try {
      std::unique_ptr<PackedReadStore> store(new PackedReadStore(phred_offset_));
      std::size_t loaded = 0;
      while (loaded < count) {
        const std::int64_t status = kseq_read(record_);
        if (status == -1) {
          int zlib_code = Z_OK;
          const char* detail = gzerror(file_, &zlib_code);
          if (zlib_code != Z_OK && zlib_code != Z_STREAM_END) {
            const std::string message = detail == nullptr ? std::string() : detail;
            throw std::runtime_error(
              "compressed FASTQ ended before its gzip stream was complete: " + path_ +
              (message.empty() ? std::string() : " (" + message + ")"));
          }
          eof_ = true;
          break;
        }
        if (status == -2) {
          throw std::runtime_error("truncated or length-mismatched FASTQ quality string in " + path_);
        }
        if (status == -3) {
          throw std::runtime_error("error while reading FASTQ stream " + path_);
        }
        if (status < 0) {
          throw std::runtime_error("unknown FASTQ parser error in " + path_);
        }
        if (!record_->is_fastq) {
          throw std::runtime_error("encountered a FASTA record in FASTQ stream " + path_);
        }
        const std::string id(record_->name.s, record_->name.l);
        const bool has_comment = keep_comment_ && record_->comment.l > 0;
        const std::string comment = has_comment
          ? std::string(record_->comment.s, record_->comment.l)
          : std::string();
        const std::string sequence(record_->seq.s, record_->seq.l);
        const std::string quality(record_->qual.s, record_->qual.l);
        store->append(id, true, comment, has_comment, sequence, &quality);
        ++loaded;
      }
      if (loaded == 0) return nullptr;
      return store;
    } catch (...) {
      failed_ = true;
      close();
      throw;
    }
  }

 private:
  std::string path_;
  int phred_offset_;
  bool keep_comment_;
  bool eof_ = false;
  bool closed_ = false;
  bool failed_ = false;
  gzFile file_ = nullptr;
  kseq_t* record_ = nullptr;
};

Rcpp::XPtr<FastqStream> stream_pointer(SEXP pointer) {
  if (TYPEOF(pointer) != EXTPTRSXP || R_ExternalPtrAddr(pointer) == nullptr) {
    Rcpp::stop("RAD stream pointer is no longer available");
  }
  if (R_ExternalPtrTag(pointer) != fastq_stream_tag()) {
    Rcpp::stop("external pointer is not a RAD FASTQ stream");
  }
  return Rcpp::XPtr<FastqStream>(pointer);
}

void require_handle_class(SEXP value, const char* class_name) {
  if (!Rf_inherits(value, class_name)) {
    Rcpp::stop(std::string("expected a ") + class_name + " handle column");
  }
}

void validate_handle_alignment(SEXP id_handles,
                               SEXP sequence_handles,
                               SEXP quality_handles,
                               SEXP comment_handles) {
  require_handle_class(id_handles, "rad_id");
  require_handle_class(sequence_handles, "rad_seq");
  require_handle_class(quality_handles, "rad_phred");
  const bool have_comments = comment_handles != R_NilValue;
  if (have_comments) require_handle_class(comment_handles, "rad_id");

  Rcpp::NumericVector ids(id_handles);
  Rcpp::NumericVector sequences(sequence_handles);
  Rcpp::NumericVector qualities(quality_handles);
  Rcpp::NumericVector comments;
  if (ids.size() != sequences.size() || ids.size() != qualities.size()) {
    Rcpp::stop("id, sequence, and quality columns must have equal length");
  }
  if (have_comments) {
    comments = Rcpp::NumericVector(comment_handles);
    if (comments.size() != ids.size()) Rcpp::stop("comment column has the wrong length");
  }

  Rcpp::XPtr<PackedReadStore> id_store = store_from_handles(id_handles);
  Rcpp::XPtr<PackedReadStore> sequence_store = store_from_handles(sequence_handles);
  Rcpp::XPtr<PackedReadStore> quality_store = store_from_handles(quality_handles);
  Rcpp::XPtr<PackedReadStore> comment_store = have_comments
    ? store_from_handles(comment_handles) : id_store;
  if (id_store.get() != sequence_store.get() || quality_store.get() != sequence_store.get() ||
      (have_comments && comment_store.get() != sequence_store.get())) {
    Rcpp::stop("rad_reads columns do not share one backing store");
  }

  const std::vector<std::size_t> sequence_indices =
    indices_from_handles(sequence_handles, *sequence_store, false);
  const std::vector<std::size_t> quality_indices =
    indices_from_handles(quality_handles, *quality_store, true);
  for (R_xlen_t i = 0; i < sequences.size(); ++i) {
    const std::size_t record = sequence_indices[i];
    const std::int64_t id_value = get_integer64_value(ids, i);
    if (id_value != sequence_store->id_handle(record, false)) {
      Rcpp::stop("seq_id and sequence handles are not aligned at row " +
                 std::to_string(i + 1));
    }
    if (have_comments) {
      const std::int64_t comment_value = get_integer64_value(comments, i);
      if (comment_value != sequence_store->id_handle(record, true)) {
        Rcpp::stop("seq_comment and sequence handles are not aligned at row " +
                   std::to_string(i + 1));
      }
    }
    if (sequence_store->has_quality(record)) {
      if (quality_indices[i] == std::numeric_limits<std::size_t>::max() ||
          quality_indices[i] != record) {
        Rcpp::stop("quality and sequence handles are not aligned at row " +
                   std::to_string(i + 1));
      }
    } else if (quality_indices[i] != std::numeric_limits<std::size_t>::max()) {
      Rcpp::stop("quality presence does not match the sequence record at row " +
                 std::to_string(i + 1));
    }
  }
}

bool ends_with(const std::string& value, const std::string& suffix) {
  return value.size() >= suffix.size() &&
    value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool contains_line_break(const std::string& value) {
  return value.find('\n') != std::string::npos ||
    value.find('\r') != std::string::npos;
}

void gz_write_all(gzFile file, const std::string& data) {
  std::size_t offset = 0;
  while (offset < data.size()) {
    const std::size_t remaining = data.size() - offset;
    const unsigned int amount = static_cast<unsigned int>(std::min<std::size_t>(
      remaining, static_cast<std::size_t>(std::numeric_limits<unsigned int>::max())));
    const int written = gzwrite(file, data.data() + offset, amount);
    if (written <= 0) throw std::runtime_error("failed while writing gzip FASTQ output");
    offset += static_cast<std::size_t>(written);
  }
}

}  // namespace

// [[Rcpp::export]]
SEXP packed_reads_build_cpp(Rcpp::CharacterVector ids,
                            Rcpp::CharacterVector sequences,
                            SEXP qualities_sexp = R_NilValue,
                            SEXP comments_sexp = R_NilValue,
                            int phred_offset = 33) {
  const R_xlen_t count = sequences.size();
  if (ids.size() != count) Rcpp::stop("ids and sequences must have equal length");

  const bool have_qualities = qualities_sexp != R_NilValue;
  const bool have_comments = comments_sexp != R_NilValue;
  Rcpp::CharacterVector qualities;
  Rcpp::CharacterVector comments;
  if (have_qualities) {
    qualities = Rcpp::CharacterVector(qualities_sexp);
    if (qualities.size() != count) Rcpp::stop("qualities and sequences must have equal length");
  }
  if (have_comments) {
    comments = Rcpp::CharacterVector(comments_sexp);
    if (comments.size() != count) Rcpp::stop("comments and sequences must have equal length");
  }

  std::unique_ptr<PackedReadStore> store(new PackedReadStore(phred_offset));
  for (R_xlen_t i = 0; i < count; ++i) {
    if (sequences[i] == NA_STRING) Rcpp::stop("sequences cannot contain NA");
    const bool id_present = ids[i] != NA_STRING;
    const std::string id = id_present
      ? std::string(Rf_translateCharUTF8(STRING_ELT(ids, i))) : std::string();
    const bool comment_present = have_comments && comments[i] != NA_STRING;
    const std::string comment = comment_present
      ? std::string(Rf_translateCharUTF8(STRING_ELT(comments, i))) : std::string();
    const std::string sequence = Rcpp::as<std::string>(sequences[i]);
    std::string quality;
    const std::string* quality_pointer = nullptr;
    if (have_qualities && qualities[i] != NA_STRING) {
      quality = Rcpp::as<std::string>(qualities[i]);
      quality_pointer = &quality;
    }
    try {
      store->append(id, id_present, comment, comment_present, sequence, quality_pointer);
    } catch (const std::exception& error) {
      std::ostringstream message;
      message << "cannot pack read " << (i + 1) << ": " << error.what();
      Rcpp::stop(message.str());
    }
  }
  return table_from_store(std::move(store), have_comments);
}

// [[Rcpp::export]]
Rcpp::CharacterVector packed_id_decode_cpp(SEXP handles) {
  Rcpp::NumericVector vector(handles);
  Rcpp::XPtr<PackedReadStore> store = store_from_handles(handles);
  const bool comment = vector.hasAttribute(".rad_field") &&
    Rcpp::as<std::string>(vector.attr(".rad_field")) == "comment";
  Rcpp::CharacterVector result(vector.size());
  for (R_xlen_t i = 0; i < vector.size(); ++i) {
    const std::int64_t handle_value = get_integer64_value(vector, i);
    if (handle_value == kInteger64NA) {
      result[i] = NA_STRING;
    } else {
      try {
        const std::string decoded = store->string_from_handle(handle_value, comment);
        if (decoded.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
          Rcpp::stop("packed identifier exceeds R's string-length limit");
        }
        result[i] = Rf_mkCharLenCE(
          decoded.data(), static_cast<int>(decoded.size()), CE_UTF8);
      } catch (const std::exception& error) {
        Rcpp::stop(error.what());
      }
    }
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::CharacterVector packed_sequence_decode_cpp(SEXP handles) {
  Rcpp::NumericVector vector(handles);
  Rcpp::XPtr<PackedReadStore> store = store_from_handles(handles);
  const std::vector<std::size_t> indices = indices_from_handles(handles, *store, false);
  Rcpp::CharacterVector result(vector.size());
  for (R_xlen_t i = 0; i < vector.size(); ++i) result[i] = store->sequence(indices[i]);
  return result;
}

// [[Rcpp::export]]
Rcpp::CharacterVector packed_quality_decode_cpp(SEXP handles, int phred_offset = 33) {
  Rcpp::NumericVector vector(handles);
  Rcpp::XPtr<PackedReadStore> store = store_from_handles(handles);
  const std::vector<std::size_t> indices = indices_from_handles(handles, *store);
  Rcpp::CharacterVector result(vector.size());
  for (R_xlen_t i = 0; i < vector.size(); ++i) {
    if (indices[i] == std::numeric_limits<std::size_t>::max()) {
      result[i] = NA_STRING;
    } else {
      result[i] = store->quality_ascii(indices[i], phred_offset);
    }
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::List packed_quality_values_cpp(SEXP handles) {
  Rcpp::NumericVector vector(handles);
  Rcpp::XPtr<PackedReadStore> store = store_from_handles(handles);
  const std::vector<std::size_t> indices = indices_from_handles(handles, *store);
  Rcpp::List result(vector.size());
  for (R_xlen_t i = 0; i < vector.size(); ++i) {
    if (indices[i] == std::numeric_limits<std::size_t>::max()) {
      result[i] = R_NilValue;
    } else {
      result[i] = Rcpp::wrap(store->quality_values(indices[i]));
    }
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector packed_quality_mean_cpp(SEXP handles) {
  Rcpp::NumericVector vector(handles);
  Rcpp::XPtr<PackedReadStore> store = store_from_handles(handles);
  const std::vector<std::size_t> indices = indices_from_handles(handles, *store);
  Rcpp::NumericVector result(vector.size(), R_NaReal);
  for (R_xlen_t i = 0; i < vector.size(); ++i) {
    if (indices[i] != std::numeric_limits<std::size_t>::max()) {
      result[i] = store->mean_quality(indices[i]);
    }
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::IntegerVector packed_quality_min_cpp(SEXP handles) {
  Rcpp::NumericVector vector(handles);
  Rcpp::XPtr<PackedReadStore> store = store_from_handles(handles);
  const std::vector<std::size_t> indices = indices_from_handles(handles, *store);
  Rcpp::IntegerVector result(vector.size(), NA_INTEGER);
  for (R_xlen_t i = 0; i < vector.size(); ++i) {
    if (indices[i] != std::numeric_limits<std::size_t>::max()) {
      result[i] = store->minimum_quality(indices[i]);
    }
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector packed_quality_fraction_cpp(SEXP handles, int threshold) {
  if (threshold < 0 || threshold > 93) Rcpp::stop("threshold must be between Q0 and Q93");
  Rcpp::NumericVector vector(handles);
  Rcpp::XPtr<PackedReadStore> store = store_from_handles(handles);
  const std::vector<std::size_t> indices = indices_from_handles(handles, *store);
  Rcpp::NumericVector result(vector.size(), R_NaReal);
  for (R_xlen_t i = 0; i < vector.size(); ++i) {
    if (indices[i] != std::numeric_limits<std::size_t>::max()) {
      result[i] = store->fraction_at_least(indices[i], threshold);
    }
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::CharacterVector packed_sequence_format_cpp(SEXP handles, int width = 24) {
  if (width < 4) width = 4;
  Rcpp::CharacterVector decoded = packed_sequence_decode_cpp(handles);
  Rcpp::CharacterVector result(decoded.size());
  for (R_xlen_t i = 0; i < decoded.size(); ++i) {
    const std::string value = Rcpp::as<std::string>(decoded[i]);
    if (static_cast<int>(value.size()) <= width) {
      result[i] = value;
    } else {
      result[i] = value.substr(0, static_cast<std::size_t>(width - 3)) + "...";
    }
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::CharacterVector packed_quality_format_cpp(SEXP handles) {
  Rcpp::NumericVector vector(handles);
  Rcpp::XPtr<PackedReadStore> store = store_from_handles(handles);
  const std::vector<std::size_t> indices = indices_from_handles(handles, *store);
  Rcpp::CharacterVector result(vector.size());
  for (R_xlen_t i = 0; i < vector.size(); ++i) {
    if (indices[i] == std::numeric_limits<std::size_t>::max()) {
      result[i] = "<no quality>";
      continue;
    }
    std::ostringstream label;
    label.setf(std::ios::fixed);
    label.precision(1);
    label << "<Q:" << store->quality_length(indices[i])
          << "; mean=" << store->mean_quality(indices[i]) << ">";
    result[i] = label.str();
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::List packed_store_stats_cpp(SEXP handles) {
  return store_from_handles(handles)->statistics();
}

// [[Rcpp::export]]
bool packed_validate_handles_cpp(SEXP id_handles,
                                 SEXP sequence_handles,
                                 SEXP quality_handles,
                                 SEXP comment_handles = R_NilValue) {
  validate_handle_alignment(id_handles, sequence_handles, quality_handles, comment_handles);
  return true;
}

// [[Rcpp::export]]
Rcpp::LogicalVector packed_handle_is_na_cpp(SEXP handles) {
  Rcpp::NumericVector vector(handles);
  Rcpp::LogicalVector result(vector.size());
  for (R_xlen_t i = 0; i < vector.size(); ++i) {
    result[i] = get_integer64_value(vector, i) == kInteger64NA;
  }
  return result;
}

// [[Rcpp::export]]
SEXP packed_fastq_stream_open_cpp(std::string path,
                                  int phred_offset = 33,
                                  bool keep_comment = true) {
  try {
    Rcpp::XPtr<FastqStream> pointer(new FastqStream(path, phred_offset, keep_comment), true);
    R_SetExternalPtrTag(pointer, fastq_stream_tag());
    pointer.attr("class") = "rad_stream_ptr";
    return pointer;
  } catch (const std::exception& error) {
    Rcpp::stop(error.what());
  }
  return R_NilValue;
}

// [[Rcpp::export]]
SEXP packed_fastq_stream_next_cpp(SEXP pointer, int chunk_size = 50000) {
  if (chunk_size <= 0) Rcpp::stop("chunk_size must be positive");
  Rcpp::XPtr<FastqStream> stream = stream_pointer(pointer);
  try {
    std::unique_ptr<PackedReadStore> store = stream->next(static_cast<std::size_t>(chunk_size));
    if (!store) return R_NilValue;
    return table_from_store(std::move(store), stream->keep_comment());
  } catch (const std::exception& error) {
    Rcpp::stop(error.what());
  }
  return R_NilValue;
}

// [[Rcpp::export]]
Rcpp::LogicalVector packed_fastq_stream_state_cpp(SEXP pointer) {
  Rcpp::XPtr<FastqStream> stream = stream_pointer(pointer);
  return Rcpp::LogicalVector::create(
    Rcpp::Named("closed") = stream->closed(),
    Rcpp::Named("done") = stream->done());
}

// [[Rcpp::export]]
void packed_fastq_stream_close_cpp(SEXP pointer) {
  Rcpp::XPtr<FastqStream> stream = stream_pointer(pointer);
  stream->close();
}

// [[Rcpp::export]]
void packed_write_fastq_cpp(SEXP id_handles,
                           SEXP sequence_handles,
                           SEXP quality_handles,
                           SEXP comment_handles,
                           std::string path,
                           int phred_offset = 33,
                           bool append = false) {
  validate_handle_alignment(id_handles, sequence_handles, quality_handles, comment_handles);
  Rcpp::NumericVector ids(id_handles);
  Rcpp::NumericVector sequences(sequence_handles);
  Rcpp::NumericVector qualities(quality_handles);
  if (ids.size() != sequences.size() || ids.size() != qualities.size()) {
    Rcpp::stop("id, sequence, and quality columns must have equal length");
  }
  const bool have_comments = comment_handles != R_NilValue;
  Rcpp::NumericVector comments;
  if (have_comments) {
    comments = Rcpp::NumericVector(comment_handles);
    if (comments.size() != ids.size()) Rcpp::stop("comment column has the wrong length");
  }

  Rcpp::XPtr<PackedReadStore> id_store = store_from_handles(id_handles);
  Rcpp::XPtr<PackedReadStore> sequence_store = store_from_handles(sequence_handles);
  Rcpp::XPtr<PackedReadStore> quality_store = store_from_handles(quality_handles);
  Rcpp::XPtr<PackedReadStore> comment_store = have_comments
    ? store_from_handles(comment_handles) : id_store;
  const std::vector<std::size_t> sequence_indices =
    indices_from_handles(sequence_handles, *sequence_store, false);
  const std::vector<std::size_t> quality_indices =
    indices_from_handles(quality_handles, *quality_store, false);

  // Validate the complete input before opening the destination. Opening in
  // overwrite mode truncates immediately, so a bad ID, comment, sequence, or
  // requested output encoding must never be allowed to destroy an existing
  // file before the error is reported.
  for (R_xlen_t i = 0; i < ids.size(); ++i) {
    const std::int64_t id_handle = get_integer64_value(ids, i);
    if (id_handle == kInteger64NA) {
      Rcpp::stop("FASTQ output cannot contain missing IDs");
    }
    const std::string id = id_store->string_from_handle(id_handle, false);
    const std::string sequence = sequence_store->sequence(sequence_indices[i]);
    const std::string quality =
      quality_store->quality_ascii(quality_indices[i], phred_offset);
    if (sequence.size() != quality.size()) {
      Rcpp::stop("sequence and quality lengths differ while writing FASTQ");
    }
    if (contains_line_break(id) || contains_line_break(sequence) ||
        contains_line_break(quality)) {
      Rcpp::stop("FASTQ fields cannot contain carriage returns or newlines");
    }
    if (have_comments) {
      const std::int64_t comment_handle = get_integer64_value(comments, i);
      if (comment_handle != kInteger64NA) {
        const std::string comment =
          comment_store->string_from_handle(comment_handle, true);
        if (contains_line_break(comment)) {
          Rcpp::stop("FASTQ fields cannot contain carriage returns or newlines");
        }
      }
    }
  }

  const bool gzip = ends_with(path, ".gz");
  gzFile gzip_file = nullptr;
  std::ofstream plain_file;
  if (gzip) {
    gzip_file = gzopen(path.c_str(), append ? "ab" : "wb");
    if (gzip_file == nullptr) Rcpp::stop("could not open output file: " + path);
  } else {
    plain_file.open(path, append ? (std::ios::out | std::ios::app) : std::ios::out);
    if (!plain_file.is_open()) Rcpp::stop("could not open output file: " + path);
  }

  try {
    for (R_xlen_t i = 0; i < ids.size(); ++i) {
      const std::int64_t id_handle = get_integer64_value(ids, i);
      if (id_handle == kInteger64NA) throw std::runtime_error("FASTQ output cannot contain missing IDs");
      const std::string id = id_store->string_from_handle(id_handle, false);
      const std::string sequence = sequence_store->sequence(sequence_indices[i]);
      const std::string quality = quality_store->quality_ascii(quality_indices[i], phred_offset);
      if (sequence.size() != quality.size()) {
        throw std::runtime_error("sequence and quality lengths differ while writing FASTQ");
      }
      if (contains_line_break(id) || contains_line_break(sequence) ||
          contains_line_break(quality)) {
        throw std::runtime_error("FASTQ fields cannot contain carriage returns or newlines");
      }
      std::string record;
      record.reserve(id.size() + sequence.size() + quality.size() + 8U);
      if (id.empty() || id[0] != '@') record.push_back('@');
      record += id;
      const std::int64_t comment_handle = have_comments
        ? get_integer64_value(comments, i) : kInteger64NA;
      if (have_comments && comment_handle != kInteger64NA) {
        const std::string comment = comment_store->string_from_handle(comment_handle, true);
        if (contains_line_break(comment)) {
          throw std::runtime_error("FASTQ fields cannot contain carriage returns or newlines");
        }
        record.push_back(' ');
        record += comment;
      }
      record += '\n';
      record += sequence;
      record += "\n+\n";
      record += quality;
      record.push_back('\n');
      if (gzip) {
        gz_write_all(gzip_file, record);
      } else {
        plain_file.write(record.data(), static_cast<std::streamsize>(record.size()));
        if (!plain_file.good()) throw std::runtime_error("failed while writing FASTQ output");
      }
    }
  } catch (...) {
    if (gzip_file != nullptr) gzclose(gzip_file);
    if (plain_file.is_open()) plain_file.close();
    throw;
  }
  if (gzip_file != nullptr && gzclose(gzip_file) != Z_OK) {
    Rcpp::stop("failed while closing gzip FASTQ output");
  }
  if (plain_file.is_open()) {
    plain_file.close();
    if (plain_file.fail()) Rcpp::stop("failed while closing FASTQ output");
  }
}
