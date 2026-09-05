#pragma once
#include "rad_headers.h"
#include <boost/math/distributions/poisson.hpp>
#include <cerrno>

inline constexpr double kDefaultScanWlMaxErrorRatio = 0.3;

inline bool scan_wl_has_suffix_ci(const std::string& value,
                                  const std::string& suffix) {
    if (value.size() < suffix.size()) return false;
    const size_t offset = value.size() - suffix.size();
    for (size_t index = 0; index < suffix.size(); ++index) {
        if (std::tolower(static_cast<unsigned char>(value[offset + index])) !=
            std::tolower(static_cast<unsigned char>(suffix[index]))) {
            return false;
        }
    }
    return true;
}

struct extracted_bc {
    std::string read_id;
    std::string sequence;
    int position;
    bool is_rc;
    std::string barcode_1;
    std::string barcode_2;
};

using barcode_count_map = std::unordered_map<std::string, uint64_t>;

// scan-wl only needs the frequency of each extracted barcode.  Keeping one
// extracted_bc (four std::strings) per matching read made memory scale with the
// total read count, even though the per-read metadata was discarded as soon as
// the whitelist statistics were calculated.
struct barcode_count_result {
    barcode_count_map counts;
    uint64_t reads_processed = 0;
    uint64_t total_extractions = 0;
    uint64_t reference_whitelist_size = 0;
    bool succeeded = true;
    bool filtered_to_whitelist = false;

    uint64_t size() const {
        return total_extractions;
    }
};

// A scan-wl barcode of 1..31 canonical bases fits in one uint64_t while still
// retaining its length.  The leading 1 bit is a sentinel; two bits are then
// appended for each base using the same encoding as int64_seq.
using packed_scan_barcode = uint64_t;
using packed_reference_set =
    phmap::flat_hash_set<packed_scan_barcode>;
using packed_observed_count_map =
    phmap::flat_hash_map<packed_scan_barcode, uint64_t>;

static inline bool pack_scan_barcode(const char* data, size_t length,
                                     packed_scan_barcode& packed) {
    if (length == 0 || length > 31) {
        return false;
    }

    packed_scan_barcode value = 1;
    for (size_t i = 0; i < length; ++i) {
        uint64_t base = 0;
        switch (data[i]) {
            case 'A': base = 0; break;
            case 'C': base = 1; break;
            case 'T': base = 2; break;
            case 'G': base = 3; break;
            default: return false;
        }
        value = (value << 2) | base;
    }
    packed = value;
    return true;
}

static inline bool pack_scan_barcode(const std::string& sequence,
                                     packed_scan_barcode& packed) {
    return pack_scan_barcode(sequence.data(), sequence.size(), packed);
}

static inline std::string unpack_scan_barcode(packed_scan_barcode packed) {
    if (packed == 0) {
        return "";
    }

    packed_scan_barcode sentinel = packed;
    size_t length = 0;
    while (sentinel > 1) {
        sentinel >>= 2;
        ++length;
    }
    if (sentinel != 1 || length == 0 || length > 31) {
        return "";
    }

    std::string sequence(length, 'A');
    for (size_t i = length; i > 0; --i) {
        switch (packed & 3U) {
            case 0: sequence[i - 1] = 'A'; break;
            case 1: sequence[i - 1] = 'C'; break;
            case 2: sequence[i - 1] = 'T'; break;
            case 3: sequence[i - 1] = 'G'; break;
        }
        packed >>= 2;
    }
    return sequence;
}

// Reference-guided scan-wl keeps its large whitelist as a compact membership
// set and allocates counters only for barcodes actually observed in the input.
// Canonical barcodes longer than 31 bases use the existing int64_seq lookup
// path, which supports multiple 64-bit chunks.
struct scan_whitelist_filter {
    packed_reference_set packed_barcodes;
    packed_observed_count_map observed_counts;
    std::unordered_set<int64_seq> long_barcodes;

    bool empty() const {
        return packed_barcodes.empty() && long_barcodes.empty();
    }

    size_t size() const {
        return packed_barcodes.size() + long_barcodes.size();
    }
};

static inline void merge_packed_whitelist_counts(
    std::vector<packed_scan_barcode>& barcodes,
    scan_whitelist_filter& whitelist_filter) {
    if (barcodes.empty()) {
        return;
    }

    std::sort(barcodes.begin(), barcodes.end());
    for (size_t first = 0; first < barcodes.size();) {
        size_t last = first + 1;
        while (last < barcodes.size() &&
               barcodes[last] == barcodes[first]) {
            ++last;
        }
        const auto key = barcodes[first];
        if (whitelist_filter.packed_barcodes.find(key) !=
            whitelist_filter.packed_barcodes.end()) {
            auto observed =
                whitelist_filter.observed_counts.try_emplace(key, 0);
            observed.first->second += static_cast<uint64_t>(last - first);
        }
        first = last;
    }
    barcodes.clear();
}

static scan_whitelist_filter load_scan_whitelist(
    const std::string& spec, uint16_t default_length, bool verbose) {
    scan_whitelist_filter out;
    const std::string path = whitelist_utils::kit_to_path(spec);
    const bool is_bitlist =
        whitelist_utils::check_if_bitlist(path, verbose, /*N=*/10);

    std::unique_ptr<std::istream> in;
    igzstream* gzip_in = nullptr;
    if (scan_wl_has_suffix_ci(path, ".gz")) {
        auto stream = std::make_unique<igzstream>(path.c_str());
        gzip_in = stream.get();
        in = std::move(stream);
    } else {
        in = std::make_unique<std::ifstream>(path);
    }
    if (!in || !*in) {
        throw std::runtime_error("Failed to open whitelist: " + path);
    }

    uint64_t malformed = 0;
    std::string line;
    while (std::getline(*in, line)) {
        if (line.empty()) {
            continue;
        }
        const auto delim = line.find_first_of(",\t");
        std::string token = seq_utils::trim(
            delim == std::string::npos ? line : line.substr(0, delim));
        if (token.empty()) {
            continue;
        }

        if (is_bitlist) {
            size_t consumed = 0;
            uint64_t raw_bits = 0;
            try {
                raw_bits = std::stoull(token, &consumed);
            } catch (...) {
                ++malformed;
                continue;
            }
            if (consumed != token.size() || default_length == 0) {
                ++malformed;
                continue;
            }

            if (default_length <= 31) {
                const uint64_t max_bits =
                    (uint64_t{1} << (2 * default_length)) - 1;
                if (raw_bits > max_bits) {
                    ++malformed;
                    continue;
                }
                const packed_scan_barcode key =
                    (uint64_t{1} << (2 * default_length)) | raw_bits;
                out.packed_barcodes.insert(key);
            } else if (default_length == 32) {
                out.long_barcodes.emplace(
                    static_cast<int64_t>(raw_bits), default_length);
            } else {
                // A one-column int64 bitlist cannot encode more than 32 bases.
                ++malformed;
            }
            continue;
        }

        packed_scan_barcode key = 0;
        if (pack_scan_barcode(token, key)) {
            out.packed_barcodes.insert(key);
            continue;
        }

        if (token.size() > 31 &&
            token.size() <= std::numeric_limits<uint16_t>::max()) {
            int64_seq encoded(token);
            if (encoded.is_valid() &&
                encoded.length == static_cast<uint16_t>(token.size())) {
                out.long_barcodes.insert(std::move(encoded));
                continue;
            }
        }
        ++malformed;
    }

    if (gzip_in) {
        const int zlib_code = gzip_in->zlib_error_code();
        const std::string zlib_message = gzip_in->zlib_error_message();
        gzip_in->close();
        if (zlib_code != Z_OK && zlib_code != Z_STREAM_END) {
            throw std::runtime_error(
                "scan whitelist gzip read failed: " + path +
                " (zlib code " + std::to_string(zlib_code) +
                (zlib_message.empty() ? ")"
                                      : ", " + zlib_message + ")"));
        }
        if (gzip_in->bad()) {
            throw std::runtime_error(
                "scan whitelist gzip close failed: " + path);
        }
    } else if (in->bad()) {
        throw std::runtime_error("scan whitelist read failed: " + path);
    }

    if (verbose && malformed > 0) {
        RAD_COUT << "[load_scan_whitelist] Skipped " << malformed
                  << " malformed whitelist row"
                  << (malformed == 1 ? "" : "s") << "\n";
    }
    return out;
}

static const EdlibEqualityPair kWildcardEqualities[8] = {
                {'N', 'A'}, 
                {'N', 'C'}, 
                {'N', 'G'}, 
                {'N', 'T'},
                {'A', 'N'}, 
                {'C', 'N'}, 
                {'G', 'N'}, 
                {'T', 'N'}
};

// Compute allowed edit distance for adapter search given a ratio.
// - Ratio <= 0 -> exact match only (0 edits).
// - Positive ratios are rounded up so small adapters still allow >=1 edit.
// - Clamp to int range for safety.
static int compute_max_edit_distance(size_t adapter_len, double max_edit_distance_ratio) {
    if (adapter_len == 0 || !std::isfinite(max_edit_distance_ratio) || max_edit_distance_ratio <= 0.0) {
        return 0;
    }
    double raw = std::floor(static_cast<double>(adapter_len) * max_edit_distance_ratio);
    raw = std::min(raw, static_cast<double>(std::numeric_limits<int>::max()));
    int max_edits = static_cast<int>(raw);
    return std::max(1, max_edits);
}

static inline bool barcode_bounds(
    std::string_view sequence, int adapter_start, int adapter_end,
    int bc_length, int m_left, int m_right, bool is_rc,
    int& start_pos, int& end_pos) {
    if (adapter_start < 0 || adapter_end < adapter_start || bc_length <= 0 ||
        m_left < 0 || m_right < 0 ||
        sequence.length() >
            static_cast<size_t>(std::numeric_limits<int>::max())) {
        return false;
    }
    const int64_t sequence_length =
        static_cast<int64_t>(sequence.length());
    if (adapter_end > sequence_length) {
        return false;
    }

    // `left_margin` is the gap between the adapter and the extracted window;
    // `right_margin` extends that window away from the adapter.  Do not clamp
    // either read boundary: a scan result must contain the complete requested
    // barcode (and the complete optional right-margin extension).
    const int64_t extraction_length =
        static_cast<int64_t>(bc_length) + m_right;
    int64_t start = 0;
    int64_t end = 0;
    if (is_rc) {
        // Mirrored forward window in the original read: skip left-margin
        // bases immediately upstream of the reverse-complement adapter, then
        // extend farther upstream by barcode_length + right_margin.
        end = static_cast<int64_t>(adapter_start) - m_left;
        start = end - extraction_length;
    } else {
        start = static_cast<int64_t>(adapter_end) + m_left;
        end = start + extraction_length;
    }
    if (start < 0 || end < 0 || start >= end || end > sequence_length ||
        start > std::numeric_limits<int>::max() ||
        end > std::numeric_limits<int>::max()) {
        return false;
    }
    start_pos = static_cast<int>(start);
    end_pos = static_cast<int>(end);
    return true;
}

// Extract a barcode relative to the matched adapter interval [start, end).
std::string extract_barcode(
    const std::string& sequence, int adapter_start, int adapter_end,
    int bc_length, int m_left, int m_right, bool is_rc) {
    int start_pos = 0;
    int end_pos = 0;
    if (!barcode_bounds(
            sequence, adapter_start, adapter_end, bc_length,
            m_left, m_right, is_rc, start_pos, end_pos)) {
        return "";
    }
    if (is_rc) {
        std::string bc_seq = sequence.substr(start_pos, end_pos - start_pos);
        return seq_utils::revcomp(bc_seq);
    }
    return sequence.substr(start_pos, end_pos - start_pos);
}

// Encode directly from the read so the reference-guided hot path does not
// build a temporary barcode string and then scan that string a second time.
// `extracted` remains true for non-canonical or >31-base sequences so they
// still contribute to total_extractions even though they cannot match the
// compact reference.
static inline bool extract_packed_barcode(
    std::string_view sequence, int adapter_start, int adapter_end,
    int bc_length, int m_left, int m_right, bool is_rc,
    packed_scan_barcode& packed, bool& extracted) {
    int start_pos = 0;
    int end_pos = 0;
    extracted = barcode_bounds(
        sequence, adapter_start, adapter_end, bc_length,
        m_left, m_right, is_rc, start_pos, end_pos);
    if (!extracted || end_pos - start_pos > 31) {
        return false;
    }

    packed_scan_barcode value = 1;
    auto append_base = [&](char base) {
        uint64_t code = 0;
        switch (base) {
            case 'A': code = is_rc ? 2 : 0; break;
            case 'C': code = is_rc ? 3 : 1; break;
            case 'T': code = is_rc ? 0 : 2; break;
            case 'G': code = is_rc ? 1 : 3; break;
            default: return false;
        }
        value = (value << 2) | code;
        return true;
    };

    if (is_rc) {
        for (int i = end_pos; i > start_pos; --i) {
            if (!append_base(sequence[static_cast<size_t>(i - 1)])) {
                return false;
            }
        }
    } else {
        for (int i = start_pos; i < end_pos; ++i) {
            if (!append_base(sequence[static_cast<size_t>(i)])) {
                return false;
            }
        }
    }
    packed = value;
    return true;
}

// Structure to hold a read for processing
struct read_chunk {
    std::string read_id;
    std::string sequence;
};

struct adapter_match {
    int start;
    int end;
    bool is_rc;
};

static inline std::optional<size_t> find_zero_edit_adapter_match(
    std::string_view sequence, std::string_view adapter) {
    const auto literal_position = sequence.find(adapter);
    if (adapter.find('N') == std::string_view::npos) {
        if (literal_position == 0) {
            return 0;
        }
        const size_t relevant_prefix =
            literal_position == std::string_view::npos
                ? sequence.size()
                : std::min(
                      sequence.size(),
                      literal_position + adapter.size() - 1);
        if (sequence.substr(0, relevant_prefix).find('N') ==
            std::string_view::npos) {
            return literal_position == std::string_view::npos
                       ? std::nullopt
                       : std::optional<size_t>(literal_position);
        }
    }
    if (adapter.size() > sequence.size()) {
        return std::nullopt;
    }

    auto equal_with_n = [](char lhs, char rhs) {
        if (lhs == rhs) {
            return true;
        }
        auto canonical = [](char base) {
            return base == 'A' || base == 'C' ||
                   base == 'G' || base == 'T';
        };
        return (lhs == 'N' && canonical(rhs)) ||
               (rhs == 'N' && canonical(lhs));
    };

    const size_t final_start = sequence.size() - adapter.size();
    for (size_t start = 0; start <= final_start; ++start) {
        size_t offset = 0;
        while (offset < adapter.size() &&
               equal_with_n(sequence[start + offset], adapter[offset])) {
            ++offset;
        }
        if (offset == adapter.size()) {
            return start;
        }
    }
    return std::nullopt;
}

static inline std::optional<adapter_match> find_adapter_match(
    std::string_view sequence,
    const std::string& adapter_seq,
    const std::string& adapter_seq_rc,
    bool exact_match_only,
    const EdlibAlignConfig& forward_config,
    const EdlibAlignConfig& rc_config) {
    // Edlib's HW/LOC empty-target result reports one location but does not
    // allocate startLocations. Treat an empty read as having no adapter.
    if (sequence.empty()) {
        return std::nullopt;
    }

    if (exact_match_only) {
        const auto forward_pos =
            find_zero_edit_adapter_match(sequence, adapter_seq);
        if (forward_pos &&
            *forward_pos <=
                static_cast<size_t>(std::numeric_limits<int>::max()) &&
            adapter_seq.size() <=
                static_cast<size_t>(std::numeric_limits<int>::max()) -
                    *forward_pos) {
            return adapter_match{
                static_cast<int>(*forward_pos),
                static_cast<int>(*forward_pos + adapter_seq.size()), false};
        }
        const auto rc_pos =
            find_zero_edit_adapter_match(sequence, adapter_seq_rc);
        if (rc_pos &&
            *rc_pos <=
                static_cast<size_t>(std::numeric_limits<int>::max()) &&
            adapter_seq_rc.size() <=
                static_cast<size_t>(std::numeric_limits<int>::max()) -
                    *rc_pos) {
            return adapter_match{
                static_cast<int>(*rc_pos),
                static_cast<int>(*rc_pos + adapter_seq_rc.size()), true};
        }
        return std::nullopt;
    }

    // An exact hit is also the globally optimal semi-global alignment.  Check
    // it before invoking Edlib so the common error-free case avoids allocating
    // and running the fuzzy matcher.  Forward remains authoritative over
    // reverse-complement matches, matching the existing search order.
    const auto forward_pos =
        find_zero_edit_adapter_match(sequence, adapter_seq);
    if (forward_pos &&
        *forward_pos <=
            static_cast<size_t>(std::numeric_limits<int>::max()) &&
        adapter_seq.size() <=
            static_cast<size_t>(std::numeric_limits<int>::max()) -
                *forward_pos) {
        return adapter_match{
            static_cast<int>(*forward_pos),
            static_cast<int>(*forward_pos + adapter_seq.size()), false};
    }

    EdlibAlignResult result = edlibAlign(
        adapter_seq.c_str(), adapter_seq.length(),
        sequence.data(), sequence.length(), forward_config);
    if (result.status == EDLIB_STATUS_OK && result.numLocations > 0 &&
        result.startLocations != nullptr && result.endLocations != nullptr &&
        result.startLocations[0] >= 0 && result.endLocations[0] >=
            result.startLocations[0] &&
        result.endLocations[0] < std::numeric_limits<int>::max()) {
        const int start = result.startLocations[0];
        const int end = result.endLocations[0] + 1;
        edlibFreeAlignResult(result);
        return adapter_match{start, end, false};
    }
    edlibFreeAlignResult(result);

    const auto rc_pos =
        find_zero_edit_adapter_match(sequence, adapter_seq_rc);
    if (rc_pos &&
        *rc_pos <= static_cast<size_t>(std::numeric_limits<int>::max()) &&
        adapter_seq_rc.size() <=
            static_cast<size_t>(std::numeric_limits<int>::max()) - *rc_pos) {
        return adapter_match{
            static_cast<int>(*rc_pos),
            static_cast<int>(*rc_pos + adapter_seq_rc.size()), true};
    }

    result = edlibAlign(
        adapter_seq_rc.c_str(), adapter_seq_rc.length(),
        sequence.data(), sequence.length(), rc_config);
    if (result.status == EDLIB_STATUS_OK && result.numLocations > 0 &&
        result.startLocations != nullptr && result.endLocations != nullptr &&
        result.startLocations[0] >= 0 && result.endLocations[0] >=
            result.startLocations[0] &&
        result.endLocations[0] < std::numeric_limits<int>::max()) {
        const int start = result.startLocations[0];
        const int end = result.endLocations[0] + 1;
        edlibFreeAlignResult(result);
        return adapter_match{start, end, true};
    }
    edlibFreeAlignResult(result);
    return std::nullopt;
}

// Process a chunk of reads while retaining only extracted barcode strings.
std::vector<std::string> process_read_chunk(const std::vector<read_chunk>& chunk,
                                            const std::string& adapter_seq,
                                            const std::string& adapter_seq_rc,
                                            int bc_length, int m_left, int m_right,
                                            double max_edit_distance_ratio) {
    std::vector<std::string> chunk_sequences;
    chunk_sequences.reserve(chunk.size() / 5);

    const int max_edits_forward = compute_max_edit_distance(adapter_seq.length(), max_edit_distance_ratio);
    const int max_edits_rc = compute_max_edit_distance(adapter_seq_rc.length(), max_edit_distance_ratio);
    const bool exact_match_only = (max_edits_forward == 0 && max_edits_rc == 0);

    EdlibAlignConfig forward_config{};
    EdlibAlignConfig rc_config{};
    if (!exact_match_only) {
        forward_config = edlibNewAlignConfig(
            max_edits_forward,
            EDLIB_MODE_HW,
            EDLIB_TASK_LOC,
            kWildcardEqualities,
            8);
        rc_config = edlibNewAlignConfig(
            max_edits_rc,
            EDLIB_MODE_HW,
            EDLIB_TASK_LOC,
            kWildcardEqualities,
            8);
    }

    const int worker_count = std::max(1, omp_get_max_threads());
    std::vector<std::vector<std::string>> per_thread(worker_count);
    for (auto& values : per_thread) {
        values.reserve(chunk.size() /
                           (10 * static_cast<size_t>(worker_count)) +
                       1);
    }
    std::atomic<bool> processing_failed{false};
    std::exception_ptr processing_error;
    std::mutex processing_error_mutex;

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        
        #pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < chunk.size(); ++i) {
            if (processing_failed.load(std::memory_order_acquire)) {
                continue;
            }
            try {
                const auto& read = chunk[i];
                const auto match = find_adapter_match(
                    read.sequence, adapter_seq, adapter_seq_rc,
                    exact_match_only, forward_config, rc_config);
                if (!match) continue;

                std::string barcode = extract_barcode(
                    read.sequence, match->start, match->end,
                    bc_length, m_left, m_right, match->is_rc);
                if (!barcode.empty()) {
                    per_thread[tid].push_back(std::move(barcode));
                }
            } catch (...) {
                {
                    std::lock_guard<std::mutex> lock(processing_error_mutex);
                    if (!processing_error) {
                        processing_error = std::current_exception();
                    }
                }
                processing_failed.store(true, std::memory_order_release);
            }
        }
    }
    if (processing_error) std::rethrow_exception(processing_error);
    for (auto& values : per_thread) {
        chunk_sequences.insert(
            chunk_sequences.end(),
            std::make_move_iterator(values.begin()),
            std::make_move_iterator(values.end()));
    }
    return chunk_sequences;
}

struct packed_chunk_result {
    std::vector<packed_scan_barcode> barcodes;
    uint64_t total_extractions = 0;
};

static packed_chunk_result process_read_chunk_packed(
    const std::vector<read_chunk>& chunk,
    const std::string& adapter_seq,
    const std::string& adapter_seq_rc,
    int bc_length, int m_left, int m_right,
    double max_edit_distance_ratio) {
    packed_chunk_result out;
    out.barcodes.reserve(chunk.size() / 5);

    const int max_edits_forward =
        compute_max_edit_distance(adapter_seq.length(), max_edit_distance_ratio);
    const int max_edits_rc =
        compute_max_edit_distance(adapter_seq_rc.length(), max_edit_distance_ratio);
    const bool exact_match_only =
        max_edits_forward == 0 && max_edits_rc == 0;

    EdlibAlignConfig forward_config{};
    EdlibAlignConfig rc_config{};
    if (!exact_match_only) {
        forward_config = edlibNewAlignConfig(
            max_edits_forward, EDLIB_MODE_HW, EDLIB_TASK_LOC,
            kWildcardEqualities, 8);
        rc_config = edlibNewAlignConfig(
            max_edits_rc, EDLIB_MODE_HW, EDLIB_TASK_LOC,
            kWildcardEqualities, 8);
    }

    const int worker_count = std::max(1, omp_get_max_threads());
    std::vector<std::vector<packed_scan_barcode>> per_thread(worker_count);
    std::vector<uint64_t> per_thread_extractions(worker_count, 0);
    for (auto& values : per_thread) {
        values.reserve(chunk.size() /
                           (10 * static_cast<size_t>(worker_count)) +
                       1);
    }
    std::atomic<bool> processing_failed{false};
    std::exception_ptr processing_error;
    std::mutex processing_error_mutex;

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();

        #pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < chunk.size(); ++i) {
            if (processing_failed.load(std::memory_order_acquire)) {
                continue;
            }
            try {
                const auto& sequence = chunk[i].sequence;
                const auto match = find_adapter_match(
                    sequence, adapter_seq, adapter_seq_rc,
                    exact_match_only, forward_config, rc_config);
                if (!match) continue;

                packed_scan_barcode packed = 0;
                bool extracted = false;
                const bool packable = extract_packed_barcode(
                    sequence, match->start, match->end,
                    bc_length, m_left, m_right, match->is_rc,
                    packed, extracted);
                if (extracted) ++per_thread_extractions[tid];
                if (packable) per_thread[tid].push_back(packed);
            } catch (...) {
                {
                    std::lock_guard<std::mutex> lock(processing_error_mutex);
                    if (!processing_error) {
                        processing_error = std::current_exception();
                    }
                }
                processing_failed.store(true, std::memory_order_release);
            }
        }
    }
    if (processing_error) std::rethrow_exception(processing_error);
    for (int tid = 0; tid < worker_count; ++tid) {
        out.total_extractions += per_thread_extractions[tid];
        out.barcodes.insert(
            out.barcodes.end(),
            std::make_move_iterator(per_thread[tid].begin()),
            std::make_move_iterator(per_thread[tid].end()));
    }
    return out;
}

// Parse FASTQ file (regular or gzipped) and extract barcodes
barcode_count_result process_fastq(const std::string& input_path,
                                   const std::string& adapter_seq,
                                   int bc_length,
                                   int m_left,
                                   int m_right,
                                   int max_reads,
                                   double max_edit_distance_ratio =
                                       kDefaultScanWlMaxErrorRatio,
                                   int chunk_size = 10000,
                                   int num_threads = 0,
                                   scan_whitelist_filter* whitelist_filter = nullptr,
                                   std::function<void()> boundary_poll = {}) {
    barcode_count_result results;
    if (chunk_size <= 0) {
        RAD_CERR << "Error: Chunk size must be greater than zero\n";
        results.succeeded = false;
        return results;
    }
    results.filtered_to_whitelist =
        whitelist_filter != nullptr && !whitelist_filter->empty();
    if (results.filtered_to_whitelist) {
        results.reference_whitelist_size = whitelist_filter->size();
        whitelist_filter->observed_counts.clear();
    }

    if (!std::isfinite(max_edit_distance_ratio) || max_edit_distance_ratio < 0.0) {
        RAD_CERR << "Warning: invalid max_edit_distance_ratio (" << max_edit_distance_ratio
                  << "); clamping to 0 (exact matches only)\n";
        max_edit_distance_ratio = 0.0;
    } else if (max_edit_distance_ratio == 0.0) {
        RAD_COUT << "Max edit distance ratio set to 0 - requiring exact adapter matches\n";
    }
    
    // Set number of OpenMP threads
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
    const int effective_threads = std::max(1, omp_get_max_threads());

    RAD_COUT << "Using " << effective_threads
              << " threads for parallel processing" << std::endl;

    std::unique_ptr<file_streaming> files_ptr;
    std::unique_ptr<read_streaming> reader_ptr;

    try {
        const int decompression_threads = effective_threads;
        files_ptr = std::make_unique<file_streaming>(
            input_path, decompression_threads,
            /*use_pigz=*/decompression_threads > 1);
        reader_ptr = std::make_unique<read_streaming>(*files_ptr);
    } catch (const std::exception& e) {
        RAD_CERR << "Error: Could not open input " << input_path << ": " << e.what() << "\n";
        results.succeeded = false;
        return results;
    }
    
    std::string adapter_seq_rc = seq_utils::revcomp(adapter_seq);
    uint64_t reads_processed = 0;
    uint64_t chunks_processed = 0;
    phmap::flat_hash_map<packed_scan_barcode, uint32_t> packed_chunk_counts;
    barcode_count_map long_chunk_counts;
    const bool use_direct_packed_extraction =
        results.filtered_to_whitelist &&
        !whitelist_filter->packed_barcodes.empty() &&
        whitelist_filter->long_barcodes.empty();
    const bool use_serial_packed_extraction =
        use_direct_packed_extraction && effective_threads == 1;
    if (results.filtered_to_whitelist) {
        if (!whitelist_filter->packed_barcodes.empty() &&
            !use_direct_packed_extraction) {
            packed_chunk_counts.reserve(static_cast<size_t>(chunk_size));
        }
        if (!whitelist_filter->long_barcodes.empty()) {
            long_chunk_counts.reserve(static_cast<size_t>(chunk_size));
        }
    }
    
    if (use_serial_packed_extraction) {
        const int max_edits_forward = compute_max_edit_distance(
            adapter_seq.length(), max_edit_distance_ratio);
        const int max_edits_rc = compute_max_edit_distance(
            adapter_seq_rc.length(), max_edit_distance_ratio);
        const bool exact_match_only =
            max_edits_forward == 0 && max_edits_rc == 0;

        EdlibAlignConfig forward_config{};
        EdlibAlignConfig rc_config{};
        if (!exact_match_only) {
            forward_config = edlibNewAlignConfig(
                max_edits_forward, EDLIB_MODE_HW, EDLIB_TASK_LOC,
                kWildcardEqualities, 8);
            rc_config = edlibNewAlignConfig(
                max_edits_rc, EDLIB_MODE_HW, EDLIB_TASK_LOC,
                kWildcardEqualities, 8);
        }

        std::vector<packed_scan_barcode> packed_chunk;
        packed_chunk.reserve(static_cast<size_t>(chunk_size));
        auto flush_packed_chunk = [&]() {
            merge_packed_whitelist_counts(
                packed_chunk, *whitelist_filter);
        };

        while (max_reads <= 0 ||
               reads_processed < static_cast<uint64_t>(max_reads)) {
            if (reads_processed % static_cast<uint64_t>(chunk_size) == 0 &&
                boundary_poll) {
                boundary_poll();
            }
            auto sequence = reader_ptr->next_sequence_view();
            if (!sequence) {
                break;
            }
            ++reads_processed;

            const auto match = find_adapter_match(
                *sequence, adapter_seq, adapter_seq_rc,
                exact_match_only, forward_config, rc_config);
            if (match) {
                packed_scan_barcode packed = 0;
                bool extracted = false;
                const bool packable = extract_packed_barcode(
                    *sequence, match->start, match->end,
                    bc_length, m_left, m_right, match->is_rc,
                    packed, extracted);
                if (extracted) {
                    ++results.total_extractions;
                }
                if (packable) {
                    packed_chunk.push_back(packed);
                }
            }

            if (reads_processed % static_cast<uint64_t>(chunk_size) == 0) {
                flush_packed_chunk();
            }
            if (reads_processed %
                    (static_cast<uint64_t>(chunk_size) * 100U) == 0) {
                chunks_processed =
                    reads_processed / static_cast<uint64_t>(chunk_size);
                RAD_COUT << "Processed " << reads_processed << " reads in "
                          << chunks_processed << " chunks, found "
                          << results.total_extractions
                          << " barcodes so far" << std::endl;
            }
        }
        flush_packed_chunk();
    } else {
    while (max_reads <= 0 ||
           reads_processed < static_cast<uint64_t>(max_reads)) {
        if (boundary_poll) boundary_poll();
        // Read a chunk of data
        std::vector<read_chunk> chunk;
        chunk.reserve(chunk_size);
        
        for (int i = 0; i < chunk_size; ++i) {
            if (max_reads > 0 &&
                reads_processed >= static_cast<uint64_t>(max_reads)) {
                break;
            }
            auto sequence = reader_ptr->next_sequence_view();
            if (!sequence) {
                break;
            }
            reads_processed++;
            // Own only the sequence across the parallel chunk; read ID,
            // comment, and quality are not consumed by whitelist scanning.
            chunk.push_back({"", std::string(*sequence)});

        }
        
        if (chunk.empty()) break;
        
        chunks_processed++;
        if (chunks_processed % 100 == 0) {
            RAD_COUT << "Processed " << reads_processed << " reads in " << chunks_processed << " chunks, found " << results.total_extractions
                      << " barcodes so far" << std::endl;
        }
        
        if (use_direct_packed_extraction) {
            auto packed_results = process_read_chunk_packed(
                chunk, adapter_seq, adapter_seq_rc, bc_length,
                m_left, m_right, max_edit_distance_ratio);
            results.total_extractions += packed_results.total_extractions;
            merge_packed_whitelist_counts(
                packed_results.barcodes, *whitelist_filter);
        } else {
            auto chunk_results = process_read_chunk(
                chunk, adapter_seq, adapter_seq_rc, bc_length,
                m_left, m_right, max_edit_distance_ratio);
            results.total_extractions +=
                static_cast<uint64_t>(chunk_results.size());

            // Mixed or long reference whitelists retain the generic string
            // fallback; de-novo mode also needs the original sequence keys.
            if (results.filtered_to_whitelist) {
                packed_chunk_counts.clear();
                long_chunk_counts.clear();
                for (auto& barcode : chunk_results) {
                    packed_scan_barcode key = 0;
                    if (!whitelist_filter->packed_barcodes.empty() &&
                        pack_scan_barcode(barcode, key)) {
                        ++packed_chunk_counts[key];
                    } else if (!whitelist_filter->long_barcodes.empty() &&
                               barcode.size() > 31) {
                        auto inserted =
                            long_chunk_counts.try_emplace(std::move(barcode), 0);
                        ++inserted.first->second;
                    }
                }
                for (const auto& kv : packed_chunk_counts) {
                    if (whitelist_filter->packed_barcodes.find(kv.first) !=
                        whitelist_filter->packed_barcodes.end()) {
                        auto observed =
                            whitelist_filter->observed_counts.try_emplace(
                                kv.first, 0);
                        observed.first->second += kv.second;
                    }
                }
                for (const auto& kv : long_chunk_counts) {
                    int64_seq encoded(kv.first);
                    if (whitelist_filter->long_barcodes.find(encoded) !=
                        whitelist_filter->long_barcodes.end()) {
                        auto inserted =
                            results.counts.try_emplace(kv.first, 0);
                        inserted.first->second += kv.second;
                    }
                }
            } else {
                for (auto& barcode : chunk_results) {
                    auto inserted =
                        results.counts.try_emplace(std::move(barcode), 0);
                    ++inserted.first->second;
                }
            }
        }
    };
    }

    results.reads_processed = reads_processed;
    if (reader_ptr->read_failed()) {
        RAD_CERR << "Error: sequence read failed (kseq code "
                  << reader_ptr->read_error_code() << ") in "
                  << reader_ptr->read_error_path() << "\n";
        results.succeeded = false;
        return results;
    }
    if (boundary_poll) boundary_poll();

    if (results.filtered_to_whitelist) {
        std::vector<packed_scan_barcode> observed_keys;
        observed_keys.reserve(whitelist_filter->observed_counts.size());
        for (const auto& kv : whitelist_filter->observed_counts) {
            observed_keys.push_back(kv.first);
        }
        std::sort(observed_keys.begin(), observed_keys.end());
        for (const auto key : observed_keys) {
            const auto found = whitelist_filter->observed_counts.find(key);
            std::string sequence = unpack_scan_barcode(key);
            if (!sequence.empty()) {
                results.counts.emplace(std::move(sequence), found->second);
            }
        }
    }
    
    RAD_COUT << "Processing complete: " << reads_processed 
              << " reads processed, " << results.total_extractions << " barcodes extracted"
              << " (" << results.counts.size() << " unique"
              << (results.filtered_to_whitelist ? " whitelist matches retained" : "")
              << ")" << std::endl;
    
    return results;
}

// -----------------------------------------------------------------------
// Two-part barcode (BC1 + BC2) support
// -----------------------------------------------------------------------

// Load a two-part barcode whitelist CSV. Text inputs contain bare barcode
// sequences (with an optional leading header). Numeric bitlists require an
// explicit sequence length so their integer payload cannot be mistaken for a
// literal barcode. Returns the sequence set plus a sorted vector of unique
// lengths found in the whitelist.
struct split_whitelist {
    std::unordered_set<std::string> seqs;
    std::vector<int> lengths;
    std::vector<std::string> ordered;              // insertion-order list for index mapping
    std::unordered_map<std::string, size_t> idx;   // sequence -> index
};

split_whitelist load_split_whitelist_csv(
    const std::string& path,
    std::optional<uint16_t> packed_length = std::nullopt) {
    split_whitelist wl;

    const bool is_bitlist =
        whitelist_utils::check_if_bitlist(path, /*verbose=*/false, /*N=*/10);
    if (is_bitlist && !packed_length) {
        throw std::invalid_argument(
            "numeric split-barcode whitelist requires an explicit barcode "
            "length: " + path);
    }
    if (packed_length && (*packed_length == 0 || *packed_length > 32)) {
        throw std::invalid_argument(
            "split-barcode bitlist length must be between 1 and 32");
    }

    std::unique_ptr<std::istream> input;
    igzstream* gzip_input = nullptr;
    if (scan_wl_has_suffix_ci(path, ".gz")) {
        auto stream = std::make_unique<igzstream>(path.c_str());
        gzip_input = stream.get();
        input = std::move(stream);
    } else {
        input = std::make_unique<std::ifstream>(path, std::ios::binary);
    }
    if (!input || !*input) {
        throw std::runtime_error(
            "could not open split-barcode whitelist: " + path);
    }

    auto is_unsigned_integer = [](const std::string& value) {
        return !value.empty() &&
               std::all_of(value.begin(), value.end(),
                           [](unsigned char ch) { return std::isdigit(ch); });
    };
    auto is_canonical_barcode = [](const std::string& value) {
        return !value.empty() &&
               std::all_of(value.begin(), value.end(), [](char base) {
                   return base == 'A' || base == 'C' ||
                          base == 'G' || base == 'T';
               });
    };
    auto lower = [](std::string value) {
        std::transform(value.begin(), value.end(), value.begin(),
                       [](unsigned char ch) {
                           return static_cast<char>(std::tolower(ch));
                       });
        return value;
    };

    bool first_content_line = true;
    size_t line_number = 0;
    std::string line;
    while (std::getline(*input, line)) {
        ++line_number;
        line = seq_utils::trim(line);
        if (line.empty()) continue;

        const auto delimiter = line.find_first_of(",\t");
        std::string token = seq_utils::trim(
            delimiter == std::string::npos ? line
                                           : line.substr(0, delimiter));
        if (first_content_line && token.size() >= 3 &&
            static_cast<unsigned char>(token[0]) == 0xef &&
            static_cast<unsigned char>(token[1]) == 0xbb &&
            static_cast<unsigned char>(token[2]) == 0xbf) {
            token.erase(0, 3);
            token = seq_utils::trim(token);
        }

        if (first_content_line) {
            first_content_line = false;
            const std::string header = lower(token);
            if (header == "whitelist_bcs" || header == "barcode" ||
                (is_bitlist && !is_unsigned_integer(token))) {
                continue;
            }
        }
        if (token.empty()) continue;

        std::string barcode;
        if (is_bitlist) {
            if (!is_unsigned_integer(token)) {
                throw std::runtime_error(
                    "invalid packed split-barcode value at line " +
                    std::to_string(line_number) + " in " + path + ": " +
                    token);
            }
            size_t consumed = 0;
            uint64_t raw_bits = 0;
            try {
                raw_bits = std::stoull(token, &consumed);
            } catch (const std::exception&) {
                throw std::runtime_error(
                    "invalid packed split-barcode value at line " +
                    std::to_string(line_number) + " in " + path + ": " +
                    token);
            }
            if (consumed != token.size()) {
                throw std::runtime_error(
                    "invalid packed split-barcode value at line " +
                    std::to_string(line_number) + " in " + path + ": " +
                    token);
            }
            if (*packed_length < 32) {
                const uint64_t max_bits =
                    (uint64_t{1} << (2 * *packed_length)) - 1;
                if (raw_bits > max_bits) {
                    throw std::runtime_error(
                        "packed split-barcode value exceeds the declared "
                        "length at line " + std::to_string(line_number) +
                        " in " + path);
                }
            }
            barcode = int64_seq(static_cast<int64_t>(raw_bits),
                                *packed_length).bits_to_sequence();
        } else {
            if (!is_canonical_barcode(token)) {
                throw std::runtime_error(
                    "invalid DNA split barcode at line " +
                    std::to_string(line_number) + " in " + path + ": " +
                    token);
            }
            barcode = std::move(token);
        }

        if (wl.seqs.insert(barcode).second) {
            wl.idx.emplace(barcode, wl.ordered.size());
            wl.lengths.push_back(static_cast<int>(barcode.size()));
            wl.ordered.push_back(std::move(barcode));
        }
    }

    if (gzip_input) {
        const int zlib_code = gzip_input->zlib_error_code();
        const std::string zlib_message = gzip_input->zlib_error_message();
        gzip_input->close();
        if (zlib_code != Z_OK && zlib_code != Z_STREAM_END) {
            throw std::runtime_error(
                "split-barcode whitelist gzip read failed: " + path +
                " (zlib code " + std::to_string(zlib_code) +
                (zlib_message.empty() ? ")"
                                      : ", " + zlib_message + ")"));
        }
        if (gzip_input->bad()) {
            throw std::runtime_error(
                "split-barcode whitelist gzip close failed: " + path);
        }
    } else if (input->bad()) {
        throw std::runtime_error(
            "split-barcode whitelist read failed: " + path);
    }

    if (wl.seqs.empty()) {
        throw std::runtime_error(
            "split-barcode whitelist contains no valid barcodes: " + path);
    }

    std::sort(wl.lengths.begin(), wl.lengths.end());
    wl.lengths.erase(std::unique(wl.lengths.begin(), wl.lengths.end()), wl.lengths.end());
    return wl;
}

// Search for a valid BC1+BC2 pair in a read window that begins at
// adapter_end + umi_length.  The extra search offset range
// [search_extra_min, search_extra_max] handles the small positional
// jitter expected in two-part barcode protocols.
// Returns split barcode fields when a valid pair is found.
struct split_barcode_hit {
    std::string bc1;
    std::string bc2;
    std::string combined() const { return bc1 + bc2; }
};

// Split scans need aggregate pair frequencies, not one heap-heavy record for
// every matching read. Keeping this map bounded by the number of distinct
// observed BC1/BC2 pairs makes whole-file scans practical for large datasets.
struct split_barcode_count_result {
    std::unordered_map<std::string, uint64_t> pair_counts;
    uint64_t reads_processed = 0;
    uint64_t total_extractions = 0;
};

inline std::string split_pair_key(const std::string& bc1,
                                  const std::string& bc2) {
    std::string key;
    key.reserve(bc1.size() + bc2.size() + 1);
    key.append(bc1);
    key.push_back('\t');
    key.append(bc2);
    return key;
}

std::optional<split_barcode_hit> extract_split_barcode(
    const std::string& seq,
    int adapter_end,
    int umi_length,
    const std::unordered_set<std::string>& bc1_set,
    const std::vector<int>& bc1_lengths,
    const std::unordered_set<std::string>& bc2_set,
    const std::vector<int>& bc2_lengths,
    int search_extra_min,
    int search_extra_max)
{
    if (seq.size() >
        static_cast<size_t>(std::numeric_limits<int>::max())) {
        return std::nullopt;
    }
    const int64_t seq_len = static_cast<int64_t>(seq.size());
    const int64_t window_base =
        static_cast<int64_t>(adapter_end) + umi_length;
    const int64_t first_extra = std::max<int64_t>(
        search_extra_min, -window_base);
    const int64_t last_extra = std::min<int64_t>(
        search_extra_max, seq_len - 1 - window_base);

    for (int64_t extra = first_extra; extra <= last_extra; ++extra) {
        const int64_t bc1_start = window_base + extra;

        for (int len1 : bc1_lengths) {
            const int64_t bc1_end = bc1_start + len1;
            if (bc1_end > seq_len) continue;

            if (bc1_set.count(seq.substr(
                    static_cast<size_t>(bc1_start),
                    static_cast<size_t>(len1)))) {
                for (int len2 : bc2_lengths) {
                    if (bc1_end + len2 > seq_len) continue;
                    std::string bc2 = seq.substr(
                        static_cast<size_t>(bc1_end),
                        static_cast<size_t>(len2));
                    if (bc2_set.count(bc2)) {
                        split_barcode_hit hit;
                        hit.bc1 = seq.substr(
                            static_cast<size_t>(bc1_start),
                            static_cast<size_t>(len1));
                        hit.bc2 = std::move(bc2);
                        return hit;
                    }
                }
            }
        }
    }
    return std::nullopt;
}

std::vector<extracted_bc> process_read_chunk_split_barcode(
    const std::vector<read_chunk>& chunk,
    const std::string& adapter_seq,
    const std::string& adapter_seq_rc,
    int umi_length,
    const std::unordered_set<std::string>& bc1_set,
    const std::vector<int>& bc1_lengths,
    const std::unordered_set<std::string>& bc2_set,
    const std::vector<int>& bc2_lengths,
    int search_extra_min,
    int search_extra_max,
    double max_edit_distance_ratio)
{
    std::vector<extracted_bc> chunk_results;
    chunk_results.reserve(chunk.size() / 5);

    const int max_edits_fwd = compute_max_edit_distance(adapter_seq.size(), max_edit_distance_ratio);
    const int max_edits_rc  = compute_max_edit_distance(adapter_seq_rc.size(), max_edit_distance_ratio);
    const bool exact_only   = (max_edits_fwd == 0 && max_edits_rc == 0);

    EdlibAlignConfig fwd_cfg{}, rc_cfg{};
    if (!exact_only) {
        fwd_cfg = edlibNewAlignConfig(max_edits_fwd, EDLIB_MODE_HW, EDLIB_TASK_LOC, kWildcardEqualities, 8);
        rc_cfg  = edlibNewAlignConfig(max_edits_rc,  EDLIB_MODE_HW, EDLIB_TASK_LOC, kWildcardEqualities, 8);
    }

    const int worker_count = std::max(1, omp_get_max_threads());
    std::vector<std::vector<extracted_bc>> per_thread(worker_count);
    for (auto& values : per_thread) {
        values.reserve(chunk.size() /
                           (5 * static_cast<size_t>(worker_count)) +
                       1);
    }
    std::atomic<bool> processing_failed{false};
    std::exception_ptr processing_error;
    std::mutex processing_error_mutex;

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();

        #pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < chunk.size(); ++i) {
            if (processing_failed.load(std::memory_order_acquire)) {
                continue;
            }
            try {
                const auto& read = chunk[i];
                auto& thread_results = per_thread[tid];

                auto try_pair = [&](const std::string& seq, int adapter_end,
                                    bool is_rc) -> bool {
                    auto hit = extract_split_barcode(
                        seq, adapter_end, umi_length,
                        bc1_set, bc1_lengths, bc2_set, bc2_lengths,
                        search_extra_min, search_extra_max);
                    if (!hit) return false;
                    thread_results.push_back({
                        read.read_id,
                        hit->combined(),
                        adapter_end,
                        is_rc,
                        hit->bc1,
                        hit->bc2
                    });
                    return true;
                };

                if (exact_only) {
                    auto pos = find_zero_edit_adapter_match(
                        read.sequence, adapter_seq);
                    if (pos &&
                        *pos + adapter_seq.size() <=
                            static_cast<size_t>(
                                std::numeric_limits<int>::max())) {
                        try_pair(read.sequence,
                                 static_cast<int>(*pos + adapter_seq.size()),
                                 false);
                    } else {
                        pos = find_zero_edit_adapter_match(
                            read.sequence, adapter_seq_rc);
                        if (pos && read.sequence.size() <=
                                       static_cast<size_t>(
                                           std::numeric_limits<int>::max())) {
                            std::string rc_seq =
                                seq_utils::revcomp(read.sequence);
                            try_pair(rc_seq,
                                     static_cast<int>(rc_seq.size() - *pos),
                                     true);
                        }
                    }
                    continue;
                }

                // Fuzzy forward
                bool found = false;
                EdlibAlignResult res = edlibAlign(
                    adapter_seq.c_str(), adapter_seq.size(),
                    read.sequence.c_str(), read.sequence.size(), fwd_cfg);
                if (res.status == EDLIB_STATUS_OK && res.numLocations > 0) {
                    found = try_pair(
                        read.sequence, res.endLocations[0] + 1, false);
                }
                edlibFreeAlignResult(res);

                if (!found) {
                    res = edlibAlign(
                        adapter_seq_rc.c_str(), adapter_seq_rc.size(),
                        read.sequence.c_str(), read.sequence.size(), rc_cfg);
                    if (res.status == EDLIB_STATUS_OK && res.numLocations > 0) {
                        std::string rc_seq =
                            seq_utils::revcomp(read.sequence);
                        // In the RC'd read the adapter ends at
                        // (read_len - start_loc).
                        try_pair(
                            rc_seq,
                            static_cast<int>(rc_seq.size() -
                                             res.startLocations[0]),
                            true);
                    }
                    edlibFreeAlignResult(res);
                }
            } catch (...) {
                {
                    std::lock_guard<std::mutex> lock(processing_error_mutex);
                    if (!processing_error) {
                        processing_error = std::current_exception();
                    }
                }
                processing_failed.store(true, std::memory_order_release);
            }
        }
    }
    if (processing_error) std::rethrow_exception(processing_error);
    for (auto& values : per_thread) {
        chunk_results.insert(
            chunk_results.end(),
            std::make_move_iterator(values.begin()),
            std::make_move_iterator(values.end()));
    }
    return chunk_results;
}

split_barcode_count_result process_fastq_split_barcode(
    const std::string& input_path,
    const std::string& adapter_seq,
    int umi_length,
    const std::unordered_set<std::string>& bc1_set,
    const std::vector<int>& bc1_lengths,
    const std::unordered_set<std::string>& bc2_set,
    const std::vector<int>& bc2_lengths,
    int search_extra_min,
    int search_extra_max,
    int max_reads,
    double max_edit_distance_ratio = kDefaultScanWlMaxErrorRatio,
    int chunk_size = 10000,
    int num_threads = 0,
    std::function<void()> boundary_poll = {})
{
    split_barcode_count_result results;

    if (chunk_size <= 0) {
        throw std::invalid_argument(
            "split-barcode scan chunk size must be greater than zero");
    }
    if (!std::isfinite(max_edit_distance_ratio) ||
        max_edit_distance_ratio < 0.0) {
        throw std::invalid_argument(
            "split-barcode scan max-error ratio must be non-negative and finite");
    }
    if (search_extra_min > search_extra_max) {
        throw std::invalid_argument(
            "split-barcode scan offset minimum exceeds its maximum");
    }
    if (bc1_set.empty() || bc2_set.empty() || bc1_lengths.empty() ||
        bc2_lengths.empty()) {
        throw std::invalid_argument(
            "split-barcode scan requires non-empty BC1 and BC2 whitelists");
    }

    if (num_threads > 0) omp_set_num_threads(num_threads);
    const int effective_threads = std::max(1, omp_get_max_threads());
    RAD_COUT << "Using " << effective_threads << " threads\n";

    std::unique_ptr<file_streaming> files_ptr;
    std::unique_ptr<read_streaming> reader_ptr;
    try {
        files_ptr = std::make_unique<file_streaming>(
            input_path, effective_threads,
            // file_streaming cannot currently surface a non-zero exit status
            // from its external pigz subprocess.  Use zlib directly here so
            // truncated/corrupt gzip input is distinguishable from clean EOF.
            /*use_pigz=*/false);
        reader_ptr = std::make_unique<read_streaming>(*files_ptr);
    } catch (const std::exception& e) {
        throw std::runtime_error(
            "could not open split-barcode input " + input_path + ": " +
            e.what());
    }

    std::string adapter_seq_rc = seq_utils::revcomp(adapter_seq);
    uint64_t chunks_processed = 0;

    while (max_reads <= 0 ||
           results.reads_processed < static_cast<uint64_t>(max_reads)) {
        if (boundary_poll) boundary_poll();
        std::vector<read_chunk> chunk;
        chunk.reserve(static_cast<size_t>(chunk_size));
        for (int i = 0; i < chunk_size; ++i) {
            if (max_reads > 0 &&
                results.reads_processed >=
                    static_cast<uint64_t>(max_reads)) {
                break;
            }
            auto rec = reader_ptr->next_sequence();
            if (!rec) break;
            ++results.reads_processed;
            chunk.push_back({rec->id, rec->seq});
        }
        if (chunk.empty()) break;

        ++chunks_processed;
        if (chunks_processed % 100 == 0) {
            RAD_COUT << "Processed " << results.reads_processed
                      << " reads, " << results.total_extractions
                      << " BC pairs found\n";
        }

        auto chunk_res = process_read_chunk_split_barcode(
            chunk, adapter_seq, adapter_seq_rc,
            umi_length, bc1_set, bc1_lengths, bc2_set, bc2_lengths,
            search_extra_min, search_extra_max, max_edit_distance_ratio);
        results.total_extractions +=
            static_cast<uint64_t>(chunk_res.size());
        for (const auto& hit : chunk_res) {
            if (hit.barcode_1.empty() || hit.barcode_2.empty()) continue;
            ++results.pair_counts[
                split_pair_key(hit.barcode_1, hit.barcode_2)];
        }
    }

    if (reader_ptr->read_failed()) {
        std::string detail = reader_ptr->read_error_detail();
        throw std::runtime_error(
            "split-barcode sequence read failed (kseq code " +
            std::to_string(reader_ptr->read_error_code()) + ") in " +
            reader_ptr->read_error_path() +
            (detail.empty() ? "" : ": " + detail));
    }
    if (boundary_poll) boundary_poll();

    RAD_COUT << "Done: " << results.reads_processed << " reads, "
              << results.total_extractions << " BC pairs extracted ("
              << results.pair_counts.size() << " unique)\n";
    return results;
}

// -----------------------------------------------------------------------

struct batch_entry {
    std::string input_fastq;
    std::string output_prefix;
};

static bool ensure_output_directory(const std::string& prefix) {
    std::filesystem::path p(prefix);
    std::filesystem::path dir = p.parent_path();
    if (dir.empty()) {
        return true; // current directory
    }

    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    if (ec) {
        RAD_CERR << "Error: Could not create output directory " << dir << ": " << ec.message() << "\n";
        return false;
    }
    return true;
}

bool write_split_pair_counts(const split_barcode_count_result& extracted_barcodes,
                             const std::string& output_csv) {
    std::vector<std::tuple<std::string, std::string, uint64_t>> rows;
    rows.reserve(extracted_barcodes.pair_counts.size());
    for (const auto& kv : extracted_barcodes.pair_counts) {
        auto tab = kv.first.find('\t');
        if (tab == std::string::npos) continue;
        rows.emplace_back(
            kv.first.substr(0, tab),
            kv.first.substr(tab + 1),
            kv.second
        );
    }

    std::sort(rows.begin(), rows.end(),
              [](const auto& a, const auto& b) {
                  if (std::get<2>(a) != std::get<2>(b)) return std::get<2>(a) > std::get<2>(b);
                  if (std::get<0>(a) != std::get<0>(b)) return std::get<0>(a) < std::get<0>(b);
                  return std::get<1>(a) < std::get<1>(b);
              });

    std::ofstream out(output_csv);
    if (!out.is_open()) {
        RAD_CERR << "Error: could not write split pair counts to "
                  << output_csv << "\n";
        return false;
    }
    out << "bc1,bc2,pair,count\n";
    for (const auto& row : rows) {
        const auto& bc1 = std::get<0>(row);
        const auto& bc2 = std::get<1>(row);
        const auto cnt = std::get<2>(row);
        out << bc1 << "," << bc2 << "," << bc1 << bc2 << "," << cnt << "\n";
    }
    out.flush();
    if (!out) {
        RAD_CERR << "Error: failed while writing split pair counts to "
                  << output_csv << "\n";
        return false;
    }
    out.close();
    if (out.fail()) {
        RAD_CERR << "Error: failed while closing split pair counts at "
                  << output_csv << "\n";
        return false;
    }
    return true;
}

/// Write a dense NxM binary matrix CSV (0s and 1s) from the observed pairs.
/// Row = BC1 index, Col = BC2 index.  The file has no header.
bool write_spat_mask_csv(
    const split_barcode_count_result& extracted_barcodes,
    const split_whitelist& bc1_wl,
    const split_whitelist& bc2_wl,
    const std::string& output_path)
{
    const size_t nrows = bc1_wl.ordered.size();
    const size_t ncols = bc2_wl.ordered.size();
    if (nrows == 0 || ncols == 0) {
        RAD_CERR << "Error: cannot write a spatial mask for empty split "
                     "whitelists\n";
        return false;
    }
    if (nrows > std::numeric_limits<size_t>::max() / ncols) {
        RAD_CERR << "Error: split whitelist dimensions overflow the spatial "
                     "mask size\n";
        return false;
    }

    // Build the bit grid
    std::vector<bool> grid(nrows * ncols, false);
    size_t pairs_set = 0;
    for (const auto& pair : extracted_barcodes.pair_counts) {
        const auto tab = pair.first.find('\t');
        if (tab == std::string::npos) continue;
        const std::string barcode_1 = pair.first.substr(0, tab);
        const std::string barcode_2 = pair.first.substr(tab + 1);
        auto it1 = bc1_wl.idx.find(barcode_1);
        auto it2 = bc2_wl.idx.find(barcode_2);
        if (it1 == bc1_wl.idx.end() || it2 == bc2_wl.idx.end()) continue;
        size_t pos = static_cast<size_t>(it1->second) * ncols + it2->second;
        if (!grid[pos]) { grid[pos] = true; ++pairs_set; }
    }

    // Write row by row
    std::ofstream out(output_path);
    if (!out.is_open()) {
        RAD_CERR << "Error: could not write spat_mask to " << output_path
                  << "\n";
        return false;
    }
    for (size_t r = 0; r < nrows; ++r) {
        for (size_t c = 0; c < ncols; ++c) {
            if (c > 0) out << ',';
            out << (grid[r * ncols + c] ? '1' : '0');
        }
        out << '\n';
    }
    out.flush();
    if (!out) {
        RAD_CERR << "Error: failed while writing spat_mask to " << output_path
                  << "\n";
        return false;
    }
    out.close();
    if (out.fail()) {
        RAD_CERR << "Error: failed while closing spat_mask at " << output_path
                  << "\n";
        return false;
    }
    RAD_COUT << "[spat_mask] Written " << nrows << "x" << ncols
              << " matrix (" << pairs_set << " valid pairs) to " << output_path << "\n";
    return true;
}

std::vector<batch_entry> read_batch_csv(const std::string& filename) {
    std::vector<batch_entry> entries;
    std::ifstream in(filename);
    if (!in) {
        RAD_CERR << "Error: Could not open batch CSV file " << filename << "\n";
        return entries;
    }

    std::string line;
    size_t line_no = 0;
    while (std::getline(in, line)) {
        line_no++;
        if (line.empty()) continue;

        std::stringstream ss(line);
        std::string fastq_path, prefix;
        std::getline(ss, fastq_path, ',');
        std::getline(ss, prefix);

        fastq_path = seq_utils::trim(fastq_path);
        prefix = seq_utils::trim(prefix);

        if (fastq_path.empty() || prefix.empty()) {
            RAD_CERR << "Warning: Skipping line " << line_no << " in " << filename
                      << " (need two non-empty columns)\n";
            continue;
        }

        entries.push_back({fastq_path, prefix});
    }

    if (entries.empty()) {
        RAD_CERR << "Warning: No valid entries found in batch CSV " << filename << "\n";
    }
    return entries;
}

// Barcode statistics structure
struct bc_wl_stats {
    std::string sequence;
    uint64_t count;
    double ncpm;
    double log1p_ncpm;
    double log1p_ncpm_ztpois;
    
    // Calculate NCPM per-barcode
    void calculate_bc_ncpm(uint64_t barcode_count, double total_reads) {
        if (total_reads > 0.0) {
            ncpm = (static_cast<double>(barcode_count) / total_reads) * 1e6;
        } else {
            ncpm = 0.0;
        }
    }
    
    // Calculate log1p_ncpm per-barcode
    void calculate_bc_log1p_ncpm() {
        log1p_ncpm = std::log1p(ncpm);
    }
    
    void calculate_bc_ztpois_pct(double k, double lambda) {
        if (k <= 0.0) {
            log1p_ncpm_ztpois = 0.0;
            return; // no zeroes allowed
        }
        if (lambda <= 0.0) {
            log1p_ncpm_ztpois = 0.0;
            return;
        }
        
        // Convert to integer for discrete distribution
        int k_int = static_cast<int>(std::floor(k));
        if (k_int <= 0) {
            log1p_ncpm_ztpois = 0.0;
            return;
        }
        
        // Pre-calculate zero probability for truncation
        double prob_zero = std::exp(-lambda);
        double truncation_denom = 1.0 - prob_zero;
        
        // Calculate regular Poisson CDF: P(X <= k)
        double regular_cdf = 0.0;
        for (int j = 0; j <= k_int; j++) {
            // Poisson PMF: P(X = j) = (lambda^j * exp(-lambda)) / j!
            double log_pmf = j * std::log(lambda) - lambda - std::lgamma(j + 1.0);
            regular_cdf += std::exp(log_pmf);
        }
        
        // Zero-truncated CDF: P(X <= k | X > 0) = [P(X <= k) - P(X = 0)] / [1 - P(X = 0)]
        double zt_cdf = (regular_cdf - prob_zero) / truncation_denom;
        
        // Clamp to [0, 1] and convert to percentage
        zt_cdf = std::max(0.0, std::min(1.0, zt_cdf));
        // Scale to percentage
        log1p_ncpm_ztpois = zt_cdf * 100.0;
    }

    // Raw barcode counts can be much larger than the transformed observations
    // used by the legacy path.  Evaluate their Poisson CDF with Boost's
    // incomplete-gamma implementation instead of looping from zero through k.
    void calculate_bc_ztpois_pct(uint64_t k, double lambda) {
        if (k == 0 || !(lambda > 0.0) || !std::isfinite(lambda)) {
            log1p_ncpm_ztpois = 0.0;
            return;
        }

        const boost::math::poisson_distribution<double> distribution(lambda);
        const double k_as_double = static_cast<double>(k);
        const double prob_zero = std::exp(-lambda);
        const double truncation_denom = -std::expm1(-lambda);
        const double upper_tail = boost::math::cdf(
            boost::math::complement(distribution, k_as_double));

        double zt_cdf = 0.0;
        if (upper_tail <= truncation_denom / 2.0) {
            zt_cdf = 1.0 - upper_tail / truncation_denom;
        } else {
            const double regular_cdf =
                boost::math::cdf(distribution, k_as_double);
            zt_cdf = (regular_cdf - prob_zero) / truncation_denom;
        }

        zt_cdf = std::max(0.0, std::min(1.0, zt_cdf));
        log1p_ncpm_ztpois = zt_cdf * 100.0;
    }

};

enum class scan_whitelist_selection {
    high_specificity,
    above_floor
};

inline const char* scan_whitelist_selection_name(
    scan_whitelist_selection selection) {
    return selection == scan_whitelist_selection::above_floor
               ? "above_floor"
               : "high_specificity";
}

struct whitelist_scan_stats {
    uint64_t reads_processed = 0;
    uint64_t total_extractions = 0;
    size_t unique_sequences = 0;
    uint64_t total_perfect_matches = 0;
    size_t unique_perfect_matches = 0;
    double match_rate_percent = 0.0;
    double floor = 0.0;
    double threshold = 0.0;
    std::string threshold_rule;
    size_t above_floor_barcodes = 0;
    size_t final_barcodes = 0;
    size_t selected_barcodes = 0;
    scan_whitelist_selection selection =
        scan_whitelist_selection::high_specificity;
};

inline bool write_whitelist_scan_log(const std::string& output_path,
                                     const std::string& input_path,
                                     const whitelist_scan_stats& stats,
                                     double wall_time_seconds) {
    std::ofstream out(output_path);
    if (!out) {
        RAD_CERR << "Error: could not open whitelist scan log for writing: "
                  << output_path << "\n";
        return false;
    }

    out << "RAD whitelist scan summary\n"
        << "input=" << input_path << "\n"
        << "selection=" << scan_whitelist_selection_name(stats.selection)
        << "\n"
        << "reads_processed=" << stats.reads_processed << "\n"
        << "total_extractions=" << stats.total_extractions << "\n"
        << "unique_sequences=" << stats.unique_sequences << "\n"
        << "total_perfect_matches=" << stats.total_perfect_matches << "\n"
        << "unique_perfect_matches=" << stats.unique_perfect_matches << "\n"
        << std::fixed << std::setprecision(6)
        << "match_rate_percent=" << stats.match_rate_percent << "\n"
        << "floor=" << stats.floor << "\n"
        << "threshold=" << stats.threshold << "\n"
        << "threshold_rule=" << stats.threshold_rule << "\n"
        << "above_floor_barcodes=" << stats.above_floor_barcodes << "\n"
        << "high_specificity_barcodes=" << stats.final_barcodes << "\n"
        << "selected_barcodes=" << stats.selected_barcodes << "\n"
        << "wall_time_seconds=" << wall_time_seconds << "\n";

    out.flush();
    if (!out) {
        RAD_CERR << "Error: failed while writing whitelist scan log: "
                  << output_path << "\n";
        return false;
    }
    out.close();
    if (out.fail()) {
        RAD_CERR << "Error: failed while closing whitelist scan log: "
                  << output_path << "\n";
        return false;
    }
    return true;
}

// Silverman's bandwidth calculation
double calculate_silverman_bandwidth(const std::vector<double>& data) {
    if (data.size() < 2) return 0.1;
    
    // Calculate standard deviation
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    double sq_sum = 0.0;
    for (double val : data) {
        sq_sum += (val - mean) * (val - mean);
    }
    double std_dev = std::sqrt(sq_sum / (data.size() - 1));
    
    // Calculate IQR
    std::vector<double> sorted_data = data;
    std::sort(sorted_data.begin(), sorted_data.end());
    double q1 = sorted_data[sorted_data.size() / 4];
    double q3 = sorted_data[3 * sorted_data.size() / 4];
    double iqr = q3 - q1;
    
    // Silverman's rule: 0.9 * min(sd, IQR/1.34) * n^(-1/5)
    double n = static_cast<double>(data.size());
    return 0.9 * std::min(std_dev, iqr / 1.34) * std::pow(n, -0.2);
}

// Simple kernel density estimation
std::vector<std::pair<double, double>> kde(const std::vector<double>& data, double bandwidth, int n_points = 512) {
    if (data.empty()) return {};
    
    double min_val = *std::min_element(data.begin(), data.end()) - 2 * bandwidth;
    double max_val = *std::max_element(data.begin(), data.end()) + 2 * bandwidth;
    double step = (max_val - min_val) / (n_points - 1);
    
    std::vector<std::pair<double, double>> density;
    
    for (int i = 0; i < n_points; ++i) {
        double x = min_val + i * step;
        double y = 0.0;
        
        // Gaussian kernel
        for (double data_point : data) {
            double z = (x - data_point) / bandwidth;
            y += std::exp(-0.5 * z * z);
        }
        
        y /= (data.size() * bandwidth * std::sqrt(2 * M_PI));
        density.push_back({x, y});
    }
    
    return density;
}

// Find peaks in density curve
std::vector<double> find_peaks(const std::vector<std::pair<double, double>>& density,
                               double min_height = 0.2, // default is 20%
                               bool debug = false) {
    std::vector<std::pair<double, double>> potential_peaks; // (x, height)

    // 1) Max and threshold
    double max_density = 0.0;
    for (const auto& point : density) max_density = std::max(max_density, point.second);
    double min_height_threshold = max_density * min_height; 
    if (debug) {
        RAD_CERR << "[find_peaks] max_density=" << std::setprecision(6) << max_density
                  << "  min_height_threshold=" << min_height_threshold << " (20%)\n";
    }

    // 2) Collect local maxima; print all local max w/ PASS/FAIL
    for (size_t i = 1; i + 1 < density.size(); ++i) {
        const double y0 = density[i-1].second;
        const double y1 = density[i].second;
        const double y2 = density[i+1].second;

        const bool is_local_max = (y1 > y0 && y1 > y2);
        if (debug && is_local_max) {
            RAD_CERR << "  local_max at x=" << std::fixed << std::setprecision(3) << density[i].first
                      << " y=" << std::setprecision(6) << y1
                      << " -> " << (y1 >= min_height_threshold ? "PASS" : "FAIL") << "\n";
        }
        if (is_local_max && y1 >= min_height_threshold) {
            potential_peaks.push_back({density[i].first, y1});
        }
    }

    // 3) Non-maximum suppression by separation
    std::vector<double> final_peaks;
    const double log_dist = 0.5; // doubled because floor doubled as well
    const double eps = 1e-9; // for x-comparison

    for (const auto& peak : potential_peaks) {
        bool keep_peak = true;
        for (const auto& other_peak : potential_peaks) {
            if (std::abs(peak.first - other_peak.first) < eps) continue; // same peak
            double distance = std::abs(peak.first - other_peak.first);
            if (distance <= log_dist && other_peak.second > peak.second) {
                keep_peak = false;
                break;
            }
        }
        if (keep_peak) final_peaks.push_back(peak.first);
    }

    if (debug) {
        RAD_CERR << "  kept peaks: ";
        for (double px : final_peaks) RAD_CERR << std::setprecision(3) << px << " ";
        RAD_CERR << "\n";
    }

    return final_peaks;
}

static inline std::vector<std::pair<double,double>>
kde_on_grid(
    const std::vector<double>& data, 
    double bandwidth,
    double min_val,
    double max_val, 
    int n_points = 512)
{
    std::vector<std::pair<double,double>> density;
    if (data.empty() || bandwidth <= 0.0 || n_points < 2) {
        return density;
    }
    density.reserve(n_points);
    const double step = (max_val - min_val) / (n_points - 1);
    const double inv_norm = 1.0 / (data.size() * bandwidth * std::sqrt(2.0 * M_PI));

    for (int i = 0; i < n_points; ++i) {
        const double x = min_val + i * step;
        double y = 0.0;
        for (double z0 : data) {
            const double z = (x - z0) / bandwidth;
            y += std::exp(-0.5 * z * z);
        }
        density.push_back({x, y * inv_norm});
    }
    return density;
}

static inline double
tau_intersection_on_grid(const std::vector<double>& xs, const std::vector<double>& LRw,
                         const std::vector<double>& BGw, double tau)
{
    const size_t n = xs.size();
    if (n == 0){
        return 0.0;
    } 
    const double c = (tau == 0.5) ? 1.0 : (1.0 - tau) / tau;

    auto diff_at = [&](size_t i){ 
        return LRw[i] - c * BGw[i]; 
    };

    // find first sign change
    for (size_t i = 0; i + 1 < n; ++i) {
        const double d0 = diff_at(i), d1 = diff_at(i+1);
        if ((d0 <= 0.0 && d1 >= 0.0) || (d0 >= 0.0 && d1 <= 0.0)) {
            // linear interp
            const double t = xs[i] + (xs[i+1] - xs[i]) * (0.0 - d0) / (d1 - d0 + 1e-300);
            return t;
        }
    }
    // no crossing: pick closest approach to zero
    size_t k = 0;
    double best = std::fabs(diff_at(0));
    for (size_t i = 1; i < n; ++i) {
        const double v = std::fabs(diff_at(i));
        if (v < best) { 
            best = v; k = i; 
        }
    }
    return xs[k];
}

static inline std::vector<double>
right_tail_area(const std::vector<double>& xs, const std::vector<double>& ys)
{
    const size_t n = xs.size();
    std::vector<double> tail(n, 0.0);
    if (n < 2) return tail;
    // segment areas
    std::vector<double> seg(n - 1, 0.0);
    for (size_t i = 0; i + 1 < n; ++i) {
        seg[i] = 0.5 * (ys[i] + ys[i+1]) * (xs[i+1] - xs[i]);
    }
    // cumulative from right
    double acc = 0.0;
    for (size_t i = n - 1; i-- > 0; ) {
        acc += seg[i];
        tail[i] = acc;
    }
    tail[n-1] = 0.0;
    return tail;
}

// ---------- Expanded-context upper-tail hazard diagnostic ----------
//
// The steepest log-count/log-rank slope has the following continuous
// equivalent in RAD's x = log1p(NCPM) coordinate:
//
//   d log(count) / d log(rank)
//       = -S(x) / ((1 - exp(-x)) f(x))
//
// where f is the density of unique-barcode scores and S is its upper-tail
// survival.  Its strongest cliff therefore corresponds to a local minimum of
//
//   H(x) = (1 - exp(-x)) f(x) / S(x).
//
// The hazard curve is also used below to decide whether RAD's preliminary gate
// has hidden a stronger whole-whitelist boundary.  A prominent minimum at or
// below that gate may own the final rank-like cutoff; minima above the gate
// remain diagnostic because they would only tighten the existing callset.
struct weighted_score {
    double x = 0.0;
    size_t multiplicity = 0;  // number of unique barcodes, never read depth
};

struct hazard_candidate {
    double x = 0.0;
    double hazard = 0.0;
    size_t estimated_rank = 0;  // empirical number of barcodes at or above x
    double prominence = 0.0;    // relative depth versus nearby hazard flanks
};

struct upper_tail_hazard_result {
    bool ok = false;
    double bw = 0.0;
    double grid_step = 0.0;
    size_t n_barcodes = 0;
    size_t n_distinct_scores = 0;
    size_t n_prominent_minima = 0;
    std::vector<hazard_candidate> prominent_minima;
    std::optional<hazard_candidate> below_floor;
    std::optional<hazard_candidate> associated_with_gate;
    std::optional<hazard_candidate> above_gate;
};

// A whole-whitelist hazard may propose moving left of an already credible
// saddle.  Measure that expansion against the finite barcode design rather
// than an expected cell count:
//
//   q0  = (M - n_positive) / M
//   rho = (K_hazard - K_preliminary) / (M - K_preliminary)
//
// The proposal is inside the finite-population trust region when rho <= q0.
// Cross multiplication avoids division; long-double products are exact for
// RAD's supported whitelist sizes while avoiding uint64_t multiplication
// overflow.  This helper only evaluates population sizes; the caller is
// responsible for invoking it solely for a credible preliminary saddle and
// count-tied callsets.
enum class finite_population_hazard_decision {
    invalid_universe,
    invalid_observed_count,
    invalid_preliminary_count,
    invalid_hazard_count,
    not_an_expansion,
    accepted,
    rejected
};

static inline const char* finite_population_hazard_decision_name(
    finite_population_hazard_decision decision)
{
    switch (decision) {
        case finite_population_hazard_decision::invalid_universe:
            return "invalid_universe";
        case finite_population_hazard_decision::invalid_observed_count:
            return "invalid_observed_count";
        case finite_population_hazard_decision::invalid_preliminary_count:
            return "invalid_preliminary_count";
        case finite_population_hazard_decision::invalid_hazard_count:
            return "invalid_hazard_count";
        case finite_population_hazard_decision::not_an_expansion:
            return "not_an_expansion";
        case finite_population_hazard_decision::accepted:
            return "accepted";
        case finite_population_hazard_decision::rejected:
            return "rejected";
    }
    return "unknown";
}

struct finite_population_hazard_result {
    bool valid = false;
    bool accept = false;
    uint64_t universe_size = 0;
    uint64_t observed_positive = 0;
    uint64_t zero_count = 0;
    uint64_t preliminary_count = 0;
    uint64_t hazard_count = 0;
    uint64_t expansion_count = 0;
    double zero_fraction = std::numeric_limits<double>::quiet_NaN();
    double residual_expansion_fraction =
        std::numeric_limits<double>::quiet_NaN();
    double finite_population_ratio =
        std::numeric_limits<double>::quiet_NaN();
    finite_population_hazard_decision decision =
        finite_population_hazard_decision::invalid_universe;
};

static inline finite_population_hazard_result
evaluate_finite_population_hazard_expansion(
    uint64_t universe_size,
    uint64_t observed_positive,
    uint64_t preliminary_count,
    uint64_t hazard_count)
{
    finite_population_hazard_result result;
    result.universe_size = universe_size;
    result.observed_positive = observed_positive;
    result.preliminary_count = preliminary_count;
    result.hazard_count = hazard_count;

    if (universe_size == 0) {
        result.decision =
            finite_population_hazard_decision::invalid_universe;
        return result;
    }
    if (observed_positive > universe_size) {
        result.decision =
            finite_population_hazard_decision::invalid_observed_count;
        return result;
    }
    if (preliminary_count > observed_positive) {
        result.decision =
            finite_population_hazard_decision::invalid_preliminary_count;
        return result;
    }
    if (hazard_count > observed_positive ||
        hazard_count < preliminary_count) {
        result.decision =
            finite_population_hazard_decision::invalid_hazard_count;
        return result;
    }

    result.valid = true;
    result.zero_count = universe_size - observed_positive;
    result.expansion_count = hazard_count - preliminary_count;
    result.zero_fraction =
        static_cast<double>(result.zero_count) /
        static_cast<double>(universe_size);

    if (result.expansion_count == 0) {
        result.accept = true;
        result.residual_expansion_fraction = 0.0;
        result.finite_population_ratio = 0.0;
        result.decision =
            finite_population_hazard_decision::not_an_expansion;
        return result;
    }

    const uint64_t residual_universe =
        universe_size - preliminary_count;
    // A positive expansion implies preliminary_count < observed_positive <= M,
    // so residual_universe is non-zero after the validation above.
    result.residual_expansion_fraction =
        static_cast<double>(result.expansion_count) /
        static_cast<double>(residual_universe);
    result.finite_population_ratio =
        result.zero_fraction > 0.0
            ? result.residual_expansion_fraction / result.zero_fraction
            : std::numeric_limits<double>::infinity();

    // delta / (M - K0) <= (M - n+) / M
    // Products of two uint64_t values may overflow, so compare in long double.
    const long double lhs =
        static_cast<long double>(result.expansion_count) *
        static_cast<long double>(universe_size);
    const long double rhs =
        static_cast<long double>(result.zero_count) *
        static_cast<long double>(residual_universe);
    result.accept = lhs <= rhs;
    result.decision = result.accept
        ? finite_population_hazard_decision::accepted
        : finite_population_hazard_decision::rejected;
    return result;
}

static inline upper_tail_hazard_result
compute_upper_tail_adjusted_hazard(
    const std::vector<weighted_score>& input_scores,
    double bandwidth,
    double floor,
    double rad_threshold,
    double upper_mode_bound = std::numeric_limits<double>::infinity(),
    double search_min = 1.0,
    size_t min_survivors = 50,
    double min_prominence = 0.05)
{
    upper_tail_hazard_result result;
    result.bw = bandwidth;
    if (!(bandwidth > 0.0) || !std::isfinite(bandwidth) ||
        !std::isfinite(floor) || !std::isfinite(rad_threshold) ||
        !std::isfinite(search_min) || input_scores.empty()) {
        return result;
    }

    // Keep every finite score.  A single source more than four bandwidths
    // below search_min has negligible influence, but a large collection of
    // such sources need not.  Score-level compression keeps this affordable
    // while avoiding a population-size-dependent boundary artefact.
    std::vector<weighted_score> scores;
    scores.reserve(input_scores.size());
    for (const auto& score : input_scores) {
        if (score.multiplicity == 0 || !std::isfinite(score.x)) {
            continue;
        }
        scores.push_back(score);
    }
    if (scores.empty()) {
        return result;
    }

    std::sort(scores.begin(), scores.end(),
              [](const weighted_score& a, const weighted_score& b) {
                  return a.x < b.x;
              });

    // Callers normally provide one row per raw-count level.  Merge again here
    // so tests and future callers cannot accidentally turn runtime back into
    // O(number of barcodes * number of grid points).
    std::vector<weighted_score> merged;
    merged.reserve(scores.size());
    for (const auto& score : scores) {
        if (!merged.empty() && score.x == merged.back().x) {
            merged.back().multiplicity += score.multiplicity;
        } else {
            merged.push_back(score);
        }
    }
    scores.swap(merged);

    size_t n_barcodes = 0;
    for (const auto& score : scores) {
        n_barcodes += score.multiplicity;
    }
    result.n_barcodes = n_barcodes;
    result.n_distinct_scores = scores.size();
    if (n_barcodes < min_survivors || scores.size() < 2) {
        return result;
    }

    const double grid_min = search_min - 4.0 * bandwidth;
    const double grid_max = scores.back().x + 4.0 * bandwidth;
    if (!(grid_max > grid_min)) {
        return result;
    }

    // Keep at least twelve grid intervals per bandwidth while bounding the
    // diagnostic cost.  RAD's x range is small enough that the upper clamp is
    // reached only for unusually tiny bandwidths.
    const double target_step = bandwidth / 12.0;
    const double desired_points =
        std::ceil((grid_max - grid_min) / target_step) + 1.0;
    size_t n_points = 4096;
    if (std::isfinite(desired_points) && desired_points < 4096.0) {
        n_points = std::max<size_t>(
            512, static_cast<size_t>(desired_points));
    }
    const double step =
        (grid_max - grid_min) / static_cast<double>(n_points - 1);
    result.grid_step = step;

    std::vector<double> xs(n_points, 0.0);
    std::vector<double> density(n_points, 0.0);
    const double inv_norm =
        1.0 / (static_cast<double>(n_barcodes) * bandwidth *
               std::sqrt(2.0 * M_PI));
    for (size_t i = 0; i < n_points; ++i) {
        const double x = grid_min + static_cast<double>(i) * step;
        xs[i] = x;
        double kernel_sum = 0.0;
        for (const auto& score : scores) {
            const double z = (x - score.x) / bandwidth;
            kernel_sum += static_cast<double>(score.multiplicity) *
                          std::exp(-0.5 * z * z);
        }
        density[i] = kernel_sum * inv_norm;
    }

    const std::vector<double> survival = right_tail_area(xs, density);
    if (survival.empty() || !(survival.front() > 0.0)) {
        return result;
    }

    std::vector<size_t> suffix(scores.size(), 0);
    size_t suffix_total = 0;
    for (size_t i = scores.size(); i-- > 0;) {
        suffix_total += scores[i].multiplicity;
        suffix[i] = suffix_total;
    }

    std::vector<double> hazard(n_points,
                               std::numeric_limits<double>::quiet_NaN());
    std::vector<size_t> empirical_rank(n_points, 0);
    for (size_t i = 0; i < n_points; ++i) {
        const double x = xs[i];
        // Evaluate enough of the left context to measure a candidate's full
        // prominence window.  Candidate reporting still begins at search_min;
        // this merely prevents that reporting bound from creating a false
        // one-sided minimum.
        if (!(x > 0.0) || !(survival[i] > 1e-15)) {
            continue;
        }
        const auto level = std::lower_bound(
            scores.begin(), scores.end(), x,
            [](const weighted_score& score, double value) {
                return score.x < value;
            });
        const size_t level_index =
            static_cast<size_t>(level - scores.begin());
        const size_t rank =
            (level_index < suffix.size()) ? suffix[level_index] : 0;
        empirical_rank[i] = rank;
        if (rank < min_survivors) {
            continue;
        }

        const double adjustment = -std::expm1(-x);
        const double value = adjustment * density[i] / survival[i];
        if (std::isfinite(value) && value >= 0.0) {
            hazard[i] = value;
        }
    }

    const size_t prominence_radius = std::max<size_t>(
        3, static_cast<size_t>(std::ceil(2.0 * bandwidth / step)));
    const double gate_slack = 0.5 * bandwidth;

    auto consider = [](std::optional<hazard_candidate>& slot,
                       const hazard_candidate& candidate) {
        if (!slot ||
            candidate.hazard < slot->hazard ||
            (candidate.hazard == slot->hazard &&
             candidate.prominence > slot->prominence)) {
            slot = candidate;
        }
    };

    for (size_t i = 1; i + 1 < n_points; ++i) {
        if (xs[i] < search_min) {
            continue;
        }
        if (!std::isfinite(hazard[i - 1]) ||
            !std::isfinite(hazard[i]) ||
            !std::isfinite(hazard[i + 1])) {
            continue;
        }
        if (!(hazard[i] <= hazard[i - 1] &&
              hazard[i] < hazard[i + 1])) {
            continue;
        }

        double left_max = hazard[i];
        double right_max = hazard[i];
        const size_t left_begin =
            (i > prominence_radius) ? i - prominence_radius : 0;
        const size_t right_end =
            std::min(n_points - 1, i + prominence_radius);
        for (size_t j = left_begin; j <= i; ++j) {
            if (std::isfinite(hazard[j])) {
                left_max = std::max(left_max, hazard[j]);
            }
        }
        for (size_t j = i; j <= right_end; ++j) {
            if (std::isfinite(hazard[j])) {
                right_max = std::max(right_max, hazard[j]);
            }
        }
        const double flank = std::min(left_max, right_max);
        const double prominence =
            (flank > 0.0)
                ? std::max(0.0, 1.0 - hazard[i] / flank)
                : 0.0;
        if (prominence < min_prominence) {
            continue;
        }

        const hazard_candidate candidate{
            xs[i], hazard[i], empirical_rank[i], prominence};
        ++result.n_prominent_minima;
        result.prominent_minima.push_back(candidate);
        if (candidate.x < floor) {
            consider(result.below_floor, candidate);
        } else if (candidate.x <= rad_threshold + gate_slack) {
            consider(result.associated_with_gate, candidate);
        } else if (!std::isfinite(upper_mode_bound) ||
                   candidate.x <= upper_mode_bound) {
            consider(result.above_gate, candidate);
        }
    }

    result.ok = true;
    return result;
}

enum class adaptive_floor_decision {
    not_evaluated,
    no_interior_minimum,
    no_stable_minimum,
    boundary_minimum,
    no_observed_level,
    insufficient_side_support,
    insufficient_gain,
    candidate_ready
};

static inline const char*
adaptive_floor_decision_name(adaptive_floor_decision decision)
{
    switch (decision) {
        case adaptive_floor_decision::not_evaluated:
            return "not_evaluated";
        case adaptive_floor_decision::no_interior_minimum:
            return "no_interior_minimum";
        case adaptive_floor_decision::no_stable_minimum:
            return "no_stable_minimum";
        case adaptive_floor_decision::boundary_minimum:
            return "boundary_minimum";
        case adaptive_floor_decision::no_observed_level:
            return "no_observed_level";
        case adaptive_floor_decision::insufficient_side_support:
            return "insufficient_side_support";
        case adaptive_floor_decision::insufficient_gain:
            return "insufficient_gain";
        case adaptive_floor_decision::candidate_ready:
            return "candidate_ready";
    }
    return "unknown";
}

struct adaptive_floor_result {
    bool evaluated = false;
    bool has_three_bandwidth_candidates = false;
    bool stable = false;
    bool sufficient_side_support = false;
    bool material_gain = false;
    bool candidate_ready = false;

    double legacy_floor = 2.0;
    double bandwidth = 0.0;
    double candidate_x_08 = std::numeric_limits<double>::quiet_NaN();
    double candidate_x_10 = std::numeric_limits<double>::quiet_NaN();
    double candidate_x_125 = std::numeric_limits<double>::quiet_NaN();
    double candidate_spread = std::numeric_limits<double>::quiet_NaN();
    double candidate_hazard = std::numeric_limits<double>::quiet_NaN();
    double candidate_prominence = std::numeric_limits<double>::quiet_NaN();
    double snapped_floor = std::numeric_limits<double>::quiet_NaN();

    size_t legacy_above_floor = 0;
    size_t candidate_above_floor = 0;
    size_t left_interval_barcodes = 0;
    size_t added_barcodes = 0;
    size_t required_gain = 0;
    adaptive_floor_decision decision =
        adaptive_floor_decision::not_evaluated;
};

// Find a reproducible adjusted-hazard candidate strictly between search_min
// and the legacy floor.  This function establishes only the candidate's
// stability, boundary clearance, observed score level, and support.  The
// downstream hazard/saddle selector decides whether it may replace RAD's
// legacy saddle.
static inline adaptive_floor_result
compute_adaptive_hazard_floor(
    const std::vector<weighted_score>& input_scores,
    double bandwidth,
    double legacy_floor = 2.0,
    double search_min = 1.0,
    size_t min_side_barcodes = 50,
    size_t min_absolute_gain = 50,
    double min_fractional_gain = 0.01,
    double max_stability_spread_h = 0.5,
    double min_prominence = 0.05)
{
    adaptive_floor_result result;
    result.legacy_floor = legacy_floor;
    result.bandwidth = bandwidth;
    result.evaluated =
        bandwidth > 0.0 && std::isfinite(bandwidth) &&
        std::isfinite(legacy_floor) && std::isfinite(search_min) &&
        std::isfinite(min_fractional_gain) &&
        min_fractional_gain >= 0.0 &&
        legacy_floor > search_min && !input_scores.empty();
    if (!result.evaluated) {
        return result;
    }

    // Evaluate exactly the same floor interval at three bandwidths.  We retain
    // every prominent minimum so that the same feature can be matched across
    // smoothing scales; independently choosing each curve's deepest minimum
    // can compare three unrelated troughs.
    const upper_tail_hazard_result hazard_08 =
        compute_upper_tail_adjusted_hazard(
            input_scores, 0.8 * bandwidth, legacy_floor, legacy_floor,
            legacy_floor, search_min, min_side_barcodes, min_prominence);
    const upper_tail_hazard_result hazard_10 =
        compute_upper_tail_adjusted_hazard(
            input_scores, bandwidth, legacy_floor, legacy_floor,
            legacy_floor, search_min, min_side_barcodes, min_prominence);
    const upper_tail_hazard_result hazard_125 =
        compute_upper_tail_adjusted_hazard(
            input_scores, 1.25 * bandwidth, legacy_floor, legacy_floor,
            legacy_floor, search_min, min_side_barcodes, min_prominence);

    auto is_interior = [&](const upper_tail_hazard_result& hazard,
                           const hazard_candidate& candidate) {
        if (!hazard.ok) return false;
        const double edge_margin = 2.0 * hazard.grid_step;
        return candidate.x > search_min + edge_margin &&
               candidate.x < legacy_floor - edge_margin;
    };

    std::vector<const hazard_candidate*> candidates_08;
    std::vector<const hazard_candidate*> candidates_10;
    std::vector<const hazard_candidate*> candidates_125;
    for (const auto& candidate : hazard_08.prominent_minima) {
        if (candidate.x < legacy_floor &&
            is_interior(hazard_08, candidate)) {
            candidates_08.push_back(&candidate);
        }
    }
    for (const auto& candidate : hazard_10.prominent_minima) {
        if (candidate.x < legacy_floor &&
            is_interior(hazard_10, candidate)) {
            candidates_10.push_back(&candidate);
        }
    }
    for (const auto& candidate : hazard_125.prominent_minima) {
        if (candidate.x < legacy_floor &&
            is_interior(hazard_125, candidate)) {
            candidates_125.push_back(&candidate);
        }
    }
    if (candidates_08.empty() || candidates_10.empty() ||
        candidates_125.empty()) {
        result.decision = adaptive_floor_decision::no_interior_minimum;
        return result;
    }

    const hazard_candidate* selected_08 = nullptr;
    const hazard_candidate* selected_10 = nullptr;
    const hazard_candidate* selected_125 = nullptr;
    double selected_spread = std::numeric_limits<double>::infinity();
    const double maximum_spread =
        max_stability_spread_h * bandwidth;
    for (const hazard_candidate* candidate_10 : candidates_10) {
        const hazard_candidate* matched_08 = nullptr;
        const hazard_candidate* matched_125 = nullptr;
        double matched_spread = std::numeric_limits<double>::infinity();
        for (const hazard_candidate* candidate_08 : candidates_08) {
            for (const hazard_candidate* candidate_125 :
                 candidates_125) {
                const double min_x = std::min(
                    {candidate_08->x, candidate_10->x,
                     candidate_125->x});
                const double max_x = std::max(
                    {candidate_08->x, candidate_10->x,
                     candidate_125->x});
                const double spread = max_x - min_x;
                if (spread <= maximum_spread &&
                    spread < matched_spread) {
                    matched_08 = candidate_08;
                    matched_125 = candidate_125;
                    matched_spread = spread;
                }
            }
        }
        if (!matched_08 || !matched_125) continue;

        result.has_three_bandwidth_candidates = true;
        const bool better =
            !selected_10 ||
            candidate_10->hazard < selected_10->hazard ||
            (candidate_10->hazard == selected_10->hazard &&
             candidate_10->prominence > selected_10->prominence) ||
            (candidate_10->hazard == selected_10->hazard &&
             candidate_10->prominence == selected_10->prominence &&
             candidate_10->x > selected_10->x);
        if (better) {
            selected_08 = matched_08;
            selected_10 = candidate_10;
            selected_125 = matched_125;
            selected_spread = matched_spread;
        }
    }
    if (!selected_10) {
        result.decision = adaptive_floor_decision::no_stable_minimum;
        return result;
    }
    result.stable = true;
    result.candidate_x_08 = selected_08->x;
    result.candidate_x_10 = selected_10->x;
    result.candidate_x_125 = selected_125->x;
    result.candidate_spread = selected_spread;
    result.candidate_hazard = selected_10->hazard;
    result.candidate_prominence = selected_10->prominence;

    // The globally deepest stable minimum wins before any acceptance guard is
    // applied.  In particular, do not discard a boundary/count-stratum
    // minimum and fall through to a shallower minimum farther right.
    const bool clear_of_search_boundary =
        selected_08->x > search_min + 0.5 * hazard_08.bw &&
        selected_10->x > search_min + 0.5 * hazard_10.bw &&
        selected_125->x > search_min + 0.5 * hazard_125.bw;
    if (!clear_of_search_boundary) {
        result.decision = adaptive_floor_decision::boundary_minimum;
        return result;
    }

    // Merge repeated score levels, then snap the continuous KDE minimum to the
    // first observed score on its right.  This makes the selected barcode rank
    // deterministic and preserves every barcode at that raw-count level.
    std::vector<weighted_score> scores;
    scores.reserve(input_scores.size());
    for (const auto& score : input_scores) {
        if (score.multiplicity > 0 && std::isfinite(score.x)) {
            scores.push_back(score);
        }
    }
    std::sort(scores.begin(), scores.end(),
              [](const weighted_score& a, const weighted_score& b) {
                  return a.x < b.x;
              });
    std::vector<weighted_score> merged;
    merged.reserve(scores.size());
    for (const auto& score : scores) {
        if (!merged.empty() && score.x == merged.back().x) {
            merged.back().multiplicity += score.multiplicity;
        } else {
            merged.push_back(score);
        }
    }
    scores.swap(merged);
    if (scores.empty()) {
        return result;
    }

    std::vector<size_t> suffix(scores.size(), 0);
    size_t total_barcodes = 0;
    for (size_t i = scores.size(); i-- > 0;) {
        total_barcodes += scores[i].multiplicity;
        suffix[i] = total_barcodes;
    }

    const auto legacy_level = std::lower_bound(
        scores.begin(), scores.end(), legacy_floor,
        [](const weighted_score& score, double value) {
            return score.x < value;
        });
    const size_t legacy_index =
        static_cast<size_t>(legacy_level - scores.begin());
    result.legacy_above_floor =
        (legacy_index < suffix.size()) ? suffix[legacy_index] : 0;

    const auto snapped_level = std::lower_bound(
        scores.begin(), scores.end(), result.candidate_x_10,
        [](const weighted_score& score, double value) {
            return score.x < value;
        });
    if (snapped_level == scores.end()) {
        result.decision = adaptive_floor_decision::no_observed_level;
        return result;
    }
    const size_t snapped_index =
        static_cast<size_t>(snapped_level - scores.begin());
    result.snapped_floor = snapped_level->x;
    if (!(result.snapped_floor > search_min) ||
        !(result.snapped_floor < legacy_floor)) {
        result.decision = adaptive_floor_decision::no_observed_level;
        return result;
    }

    result.candidate_above_floor = suffix[snapped_index];
    result.added_barcodes =
        (result.candidate_above_floor > result.legacy_above_floor)
            ? result.candidate_above_floor - result.legacy_above_floor
            : 0;
    const auto search_level = std::lower_bound(
        scores.begin(), scores.end(), search_min,
        [](const weighted_score& score, double value) {
            return score.x < value;
        });
    const size_t search_index =
        static_cast<size_t>(search_level - scores.begin());
    for (size_t i = search_index; i < snapped_index; ++i) {
        result.left_interval_barcodes += scores[i].multiplicity;
    }
    const size_t fractional_gain = static_cast<size_t>(
        std::ceil(min_fractional_gain *
                  static_cast<double>(result.legacy_above_floor)));
    result.required_gain =
        std::max(min_absolute_gain, fractional_gain);
    result.sufficient_side_support =
        result.left_interval_barcodes >= min_side_barcodes &&
        result.added_barcodes >= min_side_barcodes;
    result.material_gain =
        result.added_barcodes >= result.required_gain;
    if (!result.sufficient_side_support) {
        result.decision =
            adaptive_floor_decision::insufficient_side_support;
        return result;
    }
    if (!result.material_gain) {
        result.decision = adaptive_floor_decision::insufficient_gain;
        return result;
    }
    result.candidate_ready = true;
    result.decision = adaptive_floor_decision::candidate_ready;
    return result;
}

// ---------- AF saddle cut ----------
struct saddle_cut_result {
    bool   ok              = false;  // found >= 2 peaks and a valley
    double bw              = 0.0;
    int    n_peaks         = 0;
    double cut             = 0.0;    // raw valley x
    double final_cut       = 0.0;    // widened plateau cut (clamped by left_bound)
    bool   used_flat_widen = false;

    double valley_height   = 0.0;
    double sole_peak       = std::numeric_limits<double>::quiet_NaN();
    double left_peak       = 0.0;
    double right_peak      = 0.0;
    double plateau_left    = 0.0;    // widened plateau interval (for logging)
    double plateau_right   = 0.0;
};

static inline bool should_force_all_af_for_single_peak_ztpois(
    const saddle_cut_result& saddle,
    double ztpois_cut,
    double floor)
{
    // A sole peak immediately above the truncation floor can be a boundary
    // pile-up rather than a cell mode (as in the SN controls).  Likewise, a
    // ZTP cutoff many bandwidths into the tail is not splitting the modal
    // core.  Only override ZTP when both the peak and cutoff geometry show
    // that the fallback line actually bisects an interior unimodal lobe.
    const double modal_radius = 2.0 * saddle.bw;
    return !saddle.ok &&
           saddle.n_peaks == 1 &&
           saddle.bw > 0.0 &&
           std::isfinite(saddle.bw) &&
           std::isfinite(saddle.sole_peak) &&
           std::isfinite(ztpois_cut) &&
           std::isfinite(floor) &&
           saddle.sole_peak - floor >= modal_radius &&
           ztpois_cut >= saddle.sole_peak &&
           ztpois_cut - saddle.sole_peak <= modal_radius;
}

static inline saddle_cut_result compute_saddle_cut_af(const std::vector<double>& af_x, double bw = -1.0, int n_points = 512,
                      double min_height_ratio = 0.20, bool debug = false,
                      double left_bound = -std::numeric_limits<double>::infinity()
                      ) {
    saddle_cut_result res;
    if (af_x.size() < 10) {
        return res;
    }
    // 1) Bandwidth (Silverman fallback)
    double use_bw = (bw > 0.0) ? bw : calculate_silverman_bandwidth(af_x);
    if (!(use_bw > 0.0)) use_bw = 0.1;
    res.bw = use_bw;

    // 2) KDE on AF range
    const double xmin = *std::min_element(af_x.begin(), af_x.end());
    const double xmax = *std::max_element(af_x.begin(), af_x.end());
    auto d = kde_on_grid(af_x, use_bw, xmin, xmax, n_points);
    if (d.size() < 3) return res;

    // 3) Peaks on AF KDE
    auto pkx = find_peaks(d, /*min_height=*/min_height_ratio, /*debug=*/debug);
    res.n_peaks = static_cast<int>(pkx.size());
    if (pkx.size() == 1) {
        res.sole_peak = pkx.front();
    }
    if (pkx.size() < 2) {
        if (debug) RAD_CERR << "[saddle] < 2 peaks → no saddle\n";
        return res;
    }

    auto idx_of_x = [&](double x){
        auto it = std::lower_bound(
            d.begin(), d.end(), std::make_pair(x, -INFINITY),
            [](const auto& a, const auto& b){ return a.first < b.first; });
        size_t i = (it == d.end()) ? (d.size()-1) : static_cast<size_t>(it - d.begin());
        if (i > 0 && std::fabs(d[i].first - x) > std::fabs(d[i-1].first - x)) --i;
        return i;
    };

    struct peak_h { 
        double x; 
        double h; 
        size_t i; 
    };
    std::vector<peak_h> ph;
    ph.reserve(pkx.size());
    for (double x : pkx) {
        size_t i = idx_of_x(x);
        ph.push_back({x, d[i].second, i});
    }
    std::sort(ph.begin(), ph.end(), [](const peak_h& a, const peak_h& b){ return a.h > b.h; });

    const peak_h p1 = ph[0];
    const peak_h p2 = ph[1];
    size_t i_left  = std::min(p1.i, p2.i);
    size_t i_right = std::max(p1.i, p2.i);

    // 4) Raw valley minimum between the two tallest peaks
    size_t i_valley = i_left;
    double y_min = d[i_left].second;
    for (size_t i = i_left + 1; i < i_right; ++i) {
        if (d[i].second < y_min) { 
            y_min = d[i].second; i_valley = i; 
        }
    }

    res.ok            = (i_right > i_left + 1);
    res.cut           = d[i_valley].first;
    res.valley_height = d[i_valley].second;
    res.left_peak     = d[i_left].first;
    res.right_peak    = d[i_right].first;

    if (!res.ok) return res;

    // 5) “Keep walking” across flat valley -> widen both directions,
    //    but: do not move left of left_bound; do not cross the right peak.
    const double valley_h = res.valley_height;

    // addressing peak switchoff--we expect p1 to be greater than p2 due to background
    // but in cases with spatial data, or whether a peak needs to be pruned more than it needs to be preserved,
    // the tolerance needs to be adjusted to shift left or right to whatever peak size is dominant.
    // Get heights of spatially left and right peaks
    double left_peak_height = d[i_left].second;
    double right_peak_height = d[i_right].second;

    double climb_left = left_peak_height - valley_h;
    double climb_right = right_peak_height - valley_h;

    // Stop threshold: the density level at which we stop walking
    // Whichever is lower: half peak height, or halfway up from valley
    double stop_at_left = std::min(0.5 * left_peak_height, valley_h + 0.5 * climb_left);
    double stop_at_right = std::min(0.5 * right_peak_height, valley_h + 0.5 * climb_right);

    // Left expansion - stop when density exceeds threshold
    const double left_peak_bound = std::max(left_bound, d[i_left].first);
    size_t L = i_valley;
    while (L > i_left + 1 && d[L-1].second <= stop_at_left && d[L-1].first > left_peak_bound) {
        --L;
    }

    // Right expansion - stop when density exceeds threshold
    size_t R = i_valley;
    while (R + 1 < d.size() && (R + 1) < i_right && d[R+1].second <= stop_at_right) {
        ++R;
    }


    res.plateau_left  = std::max(left_peak_bound, d[L].first);
    res.plateau_right = d[R].first;

    // Choose the leftmost point in the flat valley (“keep walking” left)
    res.final_cut = std::max(left_peak_bound, res.plateau_left);
    res.used_flat_widen = (res.final_cut != res.cut);

    if (debug) {
        RAD_CERR << "[saddle] peaks @ " << p1.x << ", " << p2.x
                  << "  raw valley = " << res.cut << " (h=" << valley_h << ")"
                  << "  plateau=[" << res.plateau_left << ", " << res.plateau_right << "]"
                  << "  final_cut = " << res.final_cut << "  (left_bound=" << left_peak_bound
                  << ", bw = " << res.bw << ")\n";
    }
    return res;
}

enum class adaptive_saddle_decision {
    not_evaluated,
    candidate_unavailable,
    no_expanded_saddle,
    invalid_hazard_saddle_order,
    rescued_background,
    accepted
};

static inline const char*
adaptive_saddle_decision_name(adaptive_saddle_decision decision)
{
    switch (decision) {
        case adaptive_saddle_decision::not_evaluated:
            return "not_evaluated";
        case adaptive_saddle_decision::candidate_unavailable:
            return "candidate_unavailable";
        case adaptive_saddle_decision::no_expanded_saddle:
            return "no_expanded_saddle";
        case adaptive_saddle_decision::invalid_hazard_saddle_order:
            return "invalid_hazard_saddle_order";
        case adaptive_saddle_decision::rescued_background:
            return "rescued_background";
        case adaptive_saddle_decision::accepted:
            return "accepted";
    }
    return "unknown";
}

struct adaptive_saddle_result {
    bool evaluated = false;
    bool applied = false;
    double hazard_floor = std::numeric_limits<double>::quiet_NaN();
    size_t rescued_total = 0;
    size_t rescued_lr = 0;
    size_t rescued_bg = 0;
    saddle_cut_result expanded_saddle;
    adaptive_saddle_decision decision =
        adaptive_saddle_decision::not_evaluated;
};

static inline bool
valid_adaptive_hazard_saddle_order(double hazard, double left_peak,
                                   double saddle, double right_peak)
{
    return hazard < left_peak &&
           left_peak < saddle &&
           saddle < right_peak;
}

static inline bool
in_adaptive_rescued_band(double score, double hazard, double saddle)
{
    return score >= hazard && score < saddle;
}

// Decide whether a stable pre-floor hazard is the foreground boundary rather
// than merely context for RAD's saddle.  The candidate floor is evaluated
// provisionally: rejected candidates never alter the legacy x>=2 population.
// A hazard may replace the saddle only when the expanded KDE has the topology
//
//     hazard < left peak < saddle < right peak
//
// and every barcode rescued between the hazard and saddle remains on the
// legacy ZT-Poisson LR side.  The rescued interval is [hazard, saddle), so all
// raw-count ties at the selected hazard are kept deterministically.
static inline adaptive_saddle_result
evaluate_adaptive_hazard_saddle(
    const adaptive_floor_result& candidate,
    const std::vector<bc_wl_stats>& barcode_stats,
    double ztpois_rule = 95.0)
{
    adaptive_saddle_result result;
    if (!candidate.candidate_ready ||
        !std::isfinite(candidate.snapped_floor)) {
        result.decision =
            adaptive_saddle_decision::candidate_unavailable;
        return result;
    }

    result.evaluated = true;
    result.hazard_floor = candidate.snapped_floor;

    std::vector<double> expanded_af;
    expanded_af.reserve(barcode_stats.size());
    for (const auto& stats : barcode_stats) {
        if (stats.log1p_ncpm >= result.hazard_floor) {
            expanded_af.push_back(stats.log1p_ncpm);
        }
    }
    if (expanded_af.size() < 10) {
        result.decision =
            adaptive_saddle_decision::no_expanded_saddle;
        return result;
    }

    result.expanded_saddle = compute_saddle_cut_af(
        expanded_af, /*bw=*/-1.0, /*n_points=*/512,
        /*min_height_ratio=*/0.20, /*debug=*/false);
    if (!result.expanded_saddle.ok) {
        result.decision =
            adaptive_saddle_decision::no_expanded_saddle;
        return result;
    }

    const double left_peak = result.expanded_saddle.left_peak;
    const double saddle = result.expanded_saddle.final_cut;
    const double right_peak = result.expanded_saddle.right_peak;
    if (!valid_adaptive_hazard_saddle_order(
            result.hazard_floor, left_peak, saddle, right_peak)) {
        result.decision =
            adaptive_saddle_decision::invalid_hazard_saddle_order;
        return result;
    }

    for (const auto& stats : barcode_stats) {
        if (!in_adaptive_rescued_band(
                stats.log1p_ncpm, result.hazard_floor, saddle)) {
            continue;
        }
        ++result.rescued_total;
        if (stats.log1p_ncpm_ztpois >= ztpois_rule) {
            ++result.rescued_lr;
        } else {
            ++result.rescued_bg;
        }
    }
    if (result.rescued_bg != 0) {
        result.decision =
            adaptive_saddle_decision::rescued_background;
        return result;
    }

    result.applied = true;
    result.decision = adaptive_saddle_decision::accepted;
    return result;
}

bool count_perfect_matches_with_stats(const barcode_count_result& extracted_barcodes,
                                      const std::unordered_set<int64_seq>& whitelist_set,
                                      const std::string& output_csv, const std::string text_out,
                                      bool verbose = false,
                                      scan_whitelist_selection selection =
                                          scan_whitelist_selection::high_specificity,
                                      whitelist_scan_stats* summary = nullptr,
                                      const std::function<void()>& boundary_poll = {})
{
    constexpr size_t INTERRUPT_POLL_INTERVAL = 4096;
    size_t work_since_poll = 0;
    auto poll_now = [&]() {
        if (boundary_poll) boundary_poll();
        work_since_poll = 0;
    };
    auto poll_iteration = [&]() {
        if (boundary_poll &&
            ++work_since_poll >= INTERRUPT_POLL_INTERVAL) {
            poll_now();
        }
    };
    poll_now();

    // Observations were collapsed online while the FASTQ was scanned, so this
    // stage no longer allocates a second read-sized container.
    const auto& extracted_counts = extracted_barcodes.counts;
    RAD_COUT << "Using streaming barcode counts: " << extracted_counts.size()
              << " unique sequences from " << extracted_barcodes.total_extractions
              << " total extractions";
    if (extracted_barcodes.filtered_to_whitelist) {
        RAD_COUT << " (non-whitelist sequences discarded per chunk)";
    }
    RAD_COUT << "\n";

    // 2) Whitelist lookup (optional - if empty, analyze all sequences)
    bool use_whitelist =
        extracted_barcodes.filtered_to_whitelist || !whitelist_set.empty();
    if (use_whitelist) {
        const uint64_t whitelist_size =
            extracted_barcodes.reference_whitelist_size > 0
                ? extracted_barcodes.reference_whitelist_size
                : static_cast<uint64_t>(whitelist_set.size());
        RAD_COUT << "Using whitelist with " << whitelist_size << " barcodes\n";
    } else {
        RAD_COUT << "No whitelist provided - analyzing all sequences\n";
    }

    // 3) Build stats (WL barcodes if provided, otherwise all)
    std::vector<bc_wl_stats> barcode_stats;
    barcode_stats.reserve(extracted_counts.size());

    double total_reads = static_cast<double>(extracted_barcodes.total_extractions);
    uint64_t total_perfect_matches = 0;
    size_t unique_matches = 0;

    for (const auto& kv : extracted_counts) {
        poll_iteration();
        const std::string& seq = kv.first;
        uint64_t cnt = kv.second;
        
        // Skip if whitelist provided and sequence not in whitelist
        if (use_whitelist && !extracted_barcodes.filtered_to_whitelist) {
            int64_seq encoded(seq);
            if (whitelist_set.find(encoded) == whitelist_set.end()) {
                continue;
            }
        }
        bc_wl_stats s;
        s.sequence = seq;
        s.count    = cnt;
        s.calculate_bc_ncpm(cnt, total_reads);
        s.calculate_bc_log1p_ncpm();
        barcode_stats.push_back(s);

        total_perfect_matches += cnt;
        unique_matches++;
    }

    // Hash-map iteration order is intentionally unspecified and can differ
    // with chunk scheduling. Stabilize all subsequent floating-point
    // accumulation, statistics, and artifact rows by barcode sequence.
    poll_now();
    std::sort(barcode_stats.begin(), barcode_stats.end(),
              [](const bc_wl_stats& first, const bc_wl_stats& second) {
                  return first.sequence < second.sequence;
              });
    poll_now();

    // 4) zt_poisson on raw barcode counts
    double lambda = 0.0;
    if (!barcode_stats.empty()) {
        double sum = 0.0;
        for (const auto& s : barcode_stats) {
            poll_iteration();
            sum += static_cast<double>(s.count);
        }
        lambda = sum / barcode_stats.size();
    }
    RAD_COUT << "Calculated lambda for zero-truncated Poisson (raw count): "
              << lambda << "\n";

    // ZT-Poisson percentile for each barcode
    for (auto& s : barcode_stats) {
        poll_iteration();
        s.calculate_bc_ztpois_pct(s.count, lambda);
    }

    // 5) Above-floor population (log1p scale).  The legacy floor remains the
    // default.  It moves left only when the upper-tail hazard has the same
    // prominent interior minimum across three bandwidths and that move adds a
    // material, well-supported population.
    const double LEGACY_FLOOR = 2.0;
    const double ZTPOIS_RULE = 95.0; // fallback percentile if saddle/50% fail

    std::vector<double> legacy_af_x;
    legacy_af_x.reserve(barcode_stats.size());
    std::unordered_map<uint64_t, size_t> multiplicity_by_count;
    for (const auto& s : barcode_stats) {
        poll_iteration();
        if (s.log1p_ncpm >= LEGACY_FLOOR) {
            legacy_af_x.push_back(s.log1p_ncpm);
        }
        ++multiplicity_by_count[s.count];
    }

    std::vector<weighted_score> hazard_scores;
    hazard_scores.reserve(multiplicity_by_count.size());
    if (total_reads > 0.0) {
        for (const auto& count_level : multiplicity_by_count) {
            poll_iteration();
            const double ncpm =
                (static_cast<double>(count_level.first) / total_reads) *
                1e6;
            hazard_scores.push_back(
                {std::log1p(ncpm), count_level.second});
        }
    }

    double adaptive_bandwidth = 0.0;
    adaptive_floor_result adaptive_floor;
    if (legacy_af_x.size() >= 10 && !hazard_scores.empty()) {
        adaptive_bandwidth =
            calculate_silverman_bandwidth(legacy_af_x);
        if (!(adaptive_bandwidth > 0.0) ||
            !std::isfinite(adaptive_bandwidth)) {
            adaptive_bandwidth = 0.1;
        }
        adaptive_floor = compute_adaptive_hazard_floor(
            hazard_scores, adaptive_bandwidth, LEGACY_FLOOR,
            /*search_min=*/1.0, /*min_side_barcodes=*/50,
            /*min_absolute_gain=*/50,
            /*min_fractional_gain=*/0.01,
            /*max_stability_spread_h=*/0.5,
            /*min_prominence=*/0.05);
    }
    poll_now();
    const adaptive_saddle_result adaptive_saddle =
        evaluate_adaptive_hazard_saddle(
            adaptive_floor, barcode_stats, ZTPOIS_RULE);
    poll_now();
    double FLOOR = adaptive_saddle.applied
        ? adaptive_saddle.hazard_floor
        : LEGACY_FLOOR;

    RAD_COUT << "[hazard floor candidate] legacy=" << LEGACY_FLOOR
              << "  bw=" << adaptive_bandwidth;
    if (adaptive_floor.stable) {
        RAD_COUT << "  candidates="
                  << adaptive_floor.candidate_x_08 << "/"
                  << adaptive_floor.candidate_x_10 << "/"
                  << adaptive_floor.candidate_x_125
                  << "  spread=" << adaptive_floor.candidate_spread;
        if (std::isfinite(adaptive_floor.snapped_floor)) {
            RAD_COUT << "  snapped=" << adaptive_floor.snapped_floor;
        } else {
            RAD_COUT << "  snapped=unavailable";
        }
    } else {
        RAD_COUT << "  candidate=none";
    }
    RAD_COUT << "  left_support="
              << adaptive_floor.left_interval_barcodes
              << "  gain=" << adaptive_floor.added_barcodes
              << "  required=" << adaptive_floor.required_gain
              << "  decision="
              << adaptive_floor_decision_name(adaptive_floor.decision)
              << "\n";

    RAD_COUT << "[hazard/saddle selector] ";
    if (adaptive_saddle.evaluated) {
        RAD_COUT << "hazard=" << adaptive_saddle.hazard_floor;
        if (adaptive_saddle.expanded_saddle.ok) {
            RAD_COUT << "  peaks="
                      << adaptive_saddle.expanded_saddle.left_peak
                      << "/"
                      << adaptive_saddle.expanded_saddle.right_peak
                      << "  saddle="
                      << adaptive_saddle.expanded_saddle.final_cut;
        } else {
            RAD_COUT << "  saddle=unavailable";
        }
        RAD_COUT << "  rescued=" << adaptive_saddle.rescued_total
                  << "  LR_rescued=" << adaptive_saddle.rescued_lr
                  << "  BG_rescued=" << adaptive_saddle.rescued_bg;
    } else {
        RAD_COUT << "hazard=unavailable";
    }
    RAD_COUT << "  decision="
              << adaptive_saddle_decision_name(adaptive_saddle.decision)
              << "  effective_floor=" << FLOOR << "\n";

    std::vector<double> af_x; // log1p_ncpm for AF
    af_x.reserve(barcode_stats.size());
    for (const auto& s : barcode_stats){
        poll_iteration();
        if (s.log1p_ncpm >= FLOOR) {
            af_x.push_back(s.log1p_ncpm);
        }
    }

    // 6) PRIMARY THRESHOLD: AF-KDE saddle point on AF distribution
    double threshold = FLOOR;                // numeric gate (log1p_ncpm)
    std::string rule = "af_kde_saddle";
    bool have_threshold = false;
    bool force_all_af_single_peak = false;
    bool force_all_af_adaptive_hazard = false;
    bool force_final_hazard = false;

    double t_saddle = std::numeric_limits<double>::quiet_NaN();
    bool   have_saddle = false;
    saddle_cut_result sd;

    if (adaptive_saddle.applied) {
        sd = adaptive_saddle.expanded_saddle;
        t_saddle = sd.final_cut;
        threshold = FLOOR;
        have_threshold = true;
        have_saddle = true;
        force_all_af_adaptive_hazard = true;
        rule = "adaptive_hazard_bg0";
        RAD_COUT << "\n[AF-KDE saddle] peaks @ " << sd.left_peak
                  << " & " << sd.right_peak
                  << "  → saddle cut=" << sd.final_cut
                  << " (bw = " << sd.bw
                  << "); stable pre-saddle hazard accepted at "
                  << threshold;
    } else if (af_x.size() >= 10) {
        sd = compute_saddle_cut_af(af_x, /*bw=*/-1.0, /*n_points=*/512, /*min_height_ratio=*/0.20, /*debug=*/verbose);
        if (sd.ok) {
            t_saddle = sd.final_cut;
            threshold = std::max(sd.final_cut, FLOOR);
            have_threshold = true;
            have_saddle = true;
            RAD_COUT << "\n[AF-KDE saddle] peaks @ " << sd.left_peak << " & " << sd.right_peak << "  → saddle cut=" << sd.final_cut << " (bw = " << sd.bw << ")";
        } else {
            RAD_COUT << "\n[AF-KDE saddle] No stable valley; will try fallback.\n";
        }
    } else {
        RAD_COUT << "\n[AF-KDE saddle] Not enough AF points (<10); will try fallback.\n";
    }

    // 7) FALLBACK THRESHOLD: ZT-Poisson 95th percentile → smallest x meeting it
    if (!have_threshold) {
        rule = "ztpois_95pct";
        double best = std::numeric_limits<double>::infinity();
        for (const auto& s : barcode_stats) {
            poll_iteration();
            if (s.log1p_ncpm >= FLOOR && s.log1p_ncpm_ztpois >= ZTPOIS_RULE) {
                if (s.log1p_ncpm < best) best = s.log1p_ncpm;
            }
        }
        if (std::isfinite(best)) {
            have_threshold = true;
            RAD_COUT << "[ZT-Poisson fallback] threshold = " << best << " (first x with ≥ "
                      << ZTPOIS_RULE << "%)\n";
            // A cutoff on or beyond the only accepted AF mode is not a
            // between-population boundary.  Keep RAD's hard floor in that
            // case instead of discarding the modal above-floor population.
            if (should_force_all_af_for_single_peak_ztpois(
                    sd, best, FLOOR)) {
                threshold = FLOOR;
                force_all_af_single_peak = true;
                rule = "single_peak_ztpois_guard";
                RAD_COUT << "[ZT-Poisson single-peak guard] candidate " << best
                          << " is at/right of sole AF-KDE peak " << sd.sole_peak
                          << " — using FLOOR = " << FLOOR << " (return all AF)\n";
            } else {
                threshold = best;
            }
        } else {
            threshold = FLOOR;
            RAD_COUT << "[ZT-Poisson fallback] No items ≥ " << ZTPOIS_RULE
                      << "% — using FLOOR = " << FLOOR << "\n";
        }
    }

    // 7b) WHOLE-WHITELIST HAZARD SELECTOR
    //
    // The preliminary saddle/fallback gate supplies the bandwidth and an
    // upper search context, but it does not own the result when the expanded
    // whole-whitelist hazard contains a stronger pre-gate boundary.  Prefer a
    // minimum associated with the preliminary gate.  A below-floor minimum
    // must either match RAD's stable three-bandwidth feature or lie within
    // half a hazard bandwidth of the current floor, where endpoint smoothing
    // can erase a real boundary.  Above-gate minima remain diagnostic:
    // selecting them could only remove barcodes already accepted by the
    // preliminary selector.
    const double preliminary_threshold = threshold;
    const double preliminary_floor = FLOOR;
    upper_tail_hazard_result hazard_result;
    if (af_x.size() >= 10 && sd.bw > 0.0 &&
        std::isfinite(sd.bw)) {
        double upper_mode_bound =
            std::numeric_limits<double>::infinity();
        if (sd.ok && std::isfinite(sd.right_peak)) {
            upper_mode_bound = sd.right_peak;
        } else if (std::isfinite(sd.sole_peak)) {
            upper_mode_bound = sd.sole_peak;
        } else if (!af_x.empty()) {
            upper_mode_bound =
                *std::max_element(af_x.begin(), af_x.end());
        }

        poll_now();
        hazard_result = compute_upper_tail_adjusted_hazard(
            hazard_scores, sd.bw, preliminary_floor,
            preliminary_threshold, upper_mode_bound,
            /*search_min=*/1.0, /*min_survivors=*/50,
            /*min_prominence=*/0.05);
        poll_now();
    }

    auto print_hazard_candidate =
        [](const char* label,
           const std::optional<hazard_candidate>& candidate) {
            RAD_COUT << "[hazard " << label << "] ";
            if (!candidate) {
                RAD_COUT << "none\n";
                return;
            }
            RAD_COUT << "x=" << candidate->x
                      << "  H=" << candidate->hazard
                      << "  empirical_rank="
                      << candidate->estimated_rank
                      << "  prominence="
                      << candidate->prominence << "\n";
        };

    if (verbose) {
        RAD_COUT
            << "\n[upper-tail hazard diagnostic] "
            << (hazard_result.ok ? "computed" : "unavailable")
            << " (preliminary RAD gate " << preliminary_threshold << ")"
            << "\n         source_barcodes="
            << hazard_result.n_barcodes
            << "  distinct_scores="
            << hazard_result.n_distinct_scores
            << "  bw=" << hazard_result.bw
            << "  grid_step=" << hazard_result.grid_step
            << "  prominent_minima="
            << hazard_result.n_prominent_minima << "\n";
        print_hazard_candidate(
            "below-floor", hazard_result.below_floor);
        print_hazard_candidate(
            "gate-associated", hazard_result.associated_with_gate);
        print_hazard_candidate(
            "above-gate", hazard_result.above_gate);
    }

    std::optional<hazard_candidate> selected_hazard;
    const char* selected_hazard_source = nullptr;
    const char* selected_hazard_acceptance = nullptr;
    bool below_floor_considered = false;
    bool below_floor_stable_match = false;
    bool below_floor_adjacent = false;
    double below_floor_gap = std::numeric_limits<double>::quiet_NaN();
    double below_floor_adjacency_limit =
        std::numeric_limits<double>::quiet_NaN();
    if (hazard_result.associated_with_gate &&
        hazard_result.associated_with_gate->x <
            preliminary_threshold) {
        selected_hazard = hazard_result.associated_with_gate;
        selected_hazard_source = "gate-associated";
        selected_hazard_acceptance = "gate-associated";
    } else if (hazard_result.below_floor &&
               hazard_result.below_floor->x <
                   preliminary_threshold) {
        below_floor_considered = true;
        below_floor_gap =
            preliminary_floor - hazard_result.below_floor->x;
        below_floor_adjacency_limit = 0.5 * hazard_result.bw;
        below_floor_adjacent =
            below_floor_gap >= 0.0 &&
            below_floor_adjacency_limit > 0.0 &&
            below_floor_gap <= below_floor_adjacency_limit;

        // candidate_x_10 is the same-bandwidth member of the stable
        // 0.8x/1.0x/1.25x triplet.  Require the selected whole-whitelist
        // minimum to match that feature rather than allowing an unrelated
        // stable trough elsewhere in the sub-floor interval to authorize it.
        const double stable_match_limit =
            0.5 * adaptive_bandwidth;
        below_floor_stable_match =
            adaptive_floor.candidate_ready &&
            stable_match_limit > 0.0 &&
            std::isfinite(adaptive_floor.candidate_x_10) &&
            std::fabs(hazard_result.below_floor->x -
                      adaptive_floor.candidate_x_10) <=
                stable_match_limit;

        if (below_floor_stable_match || below_floor_adjacent) {
            selected_hazard = hazard_result.below_floor;
            selected_hazard_source = "below-floor";
            selected_hazard_acceptance =
                below_floor_stable_match
                    ? "three-bandwidth-stable"
                    : "floor-adjacent";
        }
    }

    // A credible KDE saddle already defines a count-tied preliminary
    // population.  Before a whole-whitelist hazard is allowed to expand that
    // population, account for saturation of the finite barcode design.  This
    // is deliberately not applied to a ZT-Poisson/no-saddle fallback: there is
    // no competing saddle population in that case.
    finite_population_hazard_result finite_population_guard;
    if (selected_hazard && have_saddle && sd.ok &&
        std::isfinite(t_saddle) &&
        selected_hazard->x < preliminary_threshold) {
        uint64_t preliminary_count = 0;
        uint64_t hazard_count = 0;
        for (const auto& stats : barcode_stats) {
            poll_iteration();
            if (stats.log1p_ncpm >= preliminary_threshold) {
                ++preliminary_count;
            }
            if (stats.log1p_ncpm >= selected_hazard->x) {
                ++hazard_count;
            }
        }
        finite_population_guard =
            evaluate_finite_population_hazard_expansion(
                extracted_barcodes.reference_whitelist_size,
                static_cast<uint64_t>(barcode_stats.size()),
                preliminary_count, hazard_count);
        RAD_COUT
            << "\n[finite-whitelist hazard correction] M="
            << finite_population_guard.universe_size
            << "  n_positive="
            << finite_population_guard.observed_positive
            << "  n_zero=" << finite_population_guard.zero_count
            << "  K_preliminary="
            << finite_population_guard.preliminary_count
            << "  K_hazard=" << finite_population_guard.hazard_count
            << "  delta=" << finite_population_guard.expansion_count;
        if (finite_population_guard.valid) {
            RAD_COUT
                << "  q0=" << finite_population_guard.zero_fraction
                << "  rho="
                << finite_population_guard.residual_expansion_fraction
                << "  ratio="
                << finite_population_guard.finite_population_ratio;
        }
        RAD_COUT
            << "  decision="
            << finite_population_hazard_decision_name(
                   finite_population_guard.decision)
            << "\n";

        // Invalid/missing M fails open for backwards compatibility.  Normal
        // scan-wl runs with a supplied whitelist populate M; replay/tests must
        // pass it explicitly to exercise this correction.
        if (finite_population_guard.valid &&
            !finite_population_guard.accept) {
            selected_hazard.reset();
        }
    }

    if (selected_hazard) {
        FLOOR = selected_hazard->x;
        threshold = selected_hazard->x;
        have_threshold = true;
        force_final_hazard = true;
        rule = "guarded_upper_tail_hazard_minimum";

        // The accepted hazard can sit on either side of the old hard floor.
        // Rebuild the population used by every downstream diagnostic so the
        // public above-floor set and final whitelist describe the same cutoff.
        af_x.clear();
        for (const auto& s : barcode_stats) {
            poll_iteration();
            if (s.log1p_ncpm >= FLOOR) {
                af_x.push_back(s.log1p_ncpm);
            }
        }

        RAD_COUT
            << "[upper-tail hazard selector] accepted "
            << selected_hazard_source
            << " x=" << selected_hazard->x
            << "  empirical_rank=" << selected_hazard->estimated_rank
            << "  previous_floor=" << preliminary_floor
            << "  previous_gate=" << preliminary_threshold
            << "  acceptance=" << selected_hazard_acceptance
            << "\n";
    } else {
        if (below_floor_considered) {
            RAD_COUT
                << "[upper-tail hazard guard] rejected below-floor x="
                << hazard_result.below_floor->x
                << "  three_bandwidth_match="
                << (below_floor_stable_match ? "true" : "false")
                << "  floor_gap=" << below_floor_gap
                << "  adjacency_limit="
                << below_floor_adjacency_limit << "\n";
        }
        RAD_COUT
            << "[upper-tail hazard selector] no eligible pre-gate minimum; "
            << "floor/gate remain " << FLOOR << "/" << threshold << "\n";
    }

    // 8) DENSITY CALCS : AF mixture using a provisional LR/BG split (ZT-Poisson 95)
    //    Also compute the 50% posterior intersection (t_purity50).
    std::vector<double> lr_x, bg_x;
    lr_x.reserve(af_x.size());
    bg_x.reserve(af_x.size());
    for (const auto& s : barcode_stats) {
        poll_iteration();
        if (s.log1p_ncpm < FLOOR) {
            continue;
        }
        if (s.log1p_ncpm_ztpois >= ZTPOIS_RULE) {
            lr_x.push_back(s.log1p_ncpm);
        } else {
            bg_x.push_back(s.log1p_ncpm);
        }
    }

    double t_purity50 = std::numeric_limits<double>::quiet_NaN();
    bool   have_tpurity50 = false;

    if(af_x.size() >= 10 && lr_x.size() >= 10 && bg_x.size() >= 10) {
        double bw = calculate_silverman_bandwidth(lr_x);

        if (!(bw > 0.0)){
            bw = 0.1;
        }

        const int KDE_N = 512;
        const double xmin = *std::min_element(af_x.begin(), af_x.end());
        const double xmax = *std::max_element(af_x.begin(), af_x.end());

        auto dR = kde_on_grid(lr_x, bw, xmin, xmax, KDE_N);  // fR|AF
        auto dB = kde_on_grid(bg_x, bw, xmin, xmax, KDE_N);  // fB|AF

        std::vector<double> xs, fR, fB;
        xs.reserve(dR.size()); fR.reserve(dR.size()); fB.reserve(dB.size());
        for (size_t i = 0; i < dR.size(); ++i) {
            xs.push_back(dR[i].first);
            fR.push_back(dR[i].second);
            fB.push_back(dB[i].second);
        }

        const double piR = static_cast<double>(lr_x.size()) /
                           static_cast<double>(lr_x.size() + bg_x.size());

        // mixture-weighted components and AF
        std::vector<double> LRw(fR.size()), BGw(fB.size()), AFw(fR.size());
        for (size_t i = 0; i < fR.size(); ++i) {
            LRw[i] = piR * fR[i];
            BGw[i] = (1.0 - piR) * fB[i];
            AFw[i] = LRw[i] + BGw[i];
        }

        RAD_COUT << "\n[AF mixture diagnostics] bw = " << bw
                  << "  |AF|=" << af_x.size()
                  << "  |LR_in_AF|=" << lr_x.size()
                  << "  |BG_in_AF|=" << bg_x.size()
                  << "  piR=" << std::setprecision(6) << piR << "\n";

        // 50% posterior intersection
        t_purity50 = tau_intersection_on_grid(xs, LRw, BGw, /*tau=*/0.5);
        have_tpurity50 = std::isfinite(t_purity50);
        if (have_tpurity50) {
            RAD_COUT << "[t_purity50] 50% boundary at x = " << t_purity50 << "\n";
        } else {
            RAD_COUT << "[t_purity50] could not compute.\n";
        }
        if ((force_all_af_single_peak ||
             force_all_af_adaptive_hazard ||
             force_final_hazard) &&
            have_tpurity50) {
            RAD_COUT << "[t_purity50] diagnostic ignored because "
                      << (force_final_hazard
                              ? "the whole-whitelist hazard owns the cutoff"
                              : (force_all_af_adaptive_hazard
                                    ? "the accepted pre-saddle hazard"
                                    : "the single-peak ZT-Poisson guard"))
                      << " forced an above-floor result\n";
            have_tpurity50 = false;
        } else if (t_purity50 > sd.right_peak) {
            RAD_COUT << "[t_purity50] " << t_purity50 << " exceeds right peak " 
                  << sd.right_peak << " — invalid, using saddle only\n";
            have_tpurity50 = false;
        }

        // A purity boundary derived from the provisional ZT-Poisson split is
        // not independent evidence and must not override a stable empirical
        // saddle.  Keep it as a diagnostic when the saddle exists; retain it
        // only as a selector for the no-saddle case.
        if (have_tpurity50 && have_saddle) {
            RAD_COUT << "[t_purity50] diagnostic only; stable saddle "
                      << t_saddle << " remains the final threshold\n";
        } else if (have_tpurity50 && !have_saddle) {
            double prev = threshold;
            threshold = std::max({t_purity50, FLOOR, prev});
            RAD_COUT << "[gate] prev=" << prev << "  -> chosen=" << threshold << "\n";
        }
    } else {
        RAD_COUT << "\n[AF mixture diagnostics] " << "  |AF|=" << af_x.size()
                  << "  |LR_in_AF|=" << lr_x.size() << "  |BG_in_AF|=" << bg_x.size() <<   "\n";
        RAD_COUT << "\n[AF mixture diagnostics] Skipped (need ≥10 in AF/LR/BG).\n";
    }

    // Prepare the "between" band for annotation (only when both lines exist)
    double band_lo = std::numeric_limits<double>::infinity();
    double band_hi = -std::numeric_limits<double>::infinity();
    bool   have_band = false;
    if (have_saddle && have_tpurity50) {
        band_lo = std::min(t_saddle, t_purity50);
        band_hi = std::max(t_saddle, t_purity50);
        have_band = true;
    }

    // 9) COLLAPSE DECISION: combine left-tail fraction + spread check
    //    (a) left-tail %: if AF below threshold is small (≤10%), return ALL AF
    //    (b) spread: if both sides are narrow (≤1.0) OR full range ≤ 2.0, return ALL AF
    // 9b) LEFT_BUDGET fallback: if collapse did not fire but left_frac is only
    //     marginally over the collapse threshold (10% < left_frac ≤ 20%),
    //     rescue the top LEFT_TAIL_MAX_FRAC * af_total barcodes from AF_left
    //     (sorted by count desc). This recovers real cells that fall just
    //     under the ZT-Poisson threshold without admitting the deep noise tail.
    //     Above LEFT_TAIL_BUDGET_MAX_FRAC (20%) the threshold is treated as
    //     unreliable / the data is too noisy and we fall back to the strict
    //     "drop AF_left" behavior.
    const double LEFT_TAIL_MAX_FRAC        = 0.10;  // 10%
    const double LEFT_TAIL_BUDGET_MAX_FRAC = 0.20;  // 20% — upper bound for budget fallback
    const double SPREAD_EPS = 1.0;           // in log1p_ncpm units
    const int    MIN_SIDE_N = 5;             // avoid tiny/empty groups

    size_t af_total = 0, af_left = 0, af_right = 0;
    double left_min  =  std::numeric_limits<double>::infinity();
    double right_min =  std::numeric_limits<double>::infinity();
    double left_max  = -std::numeric_limits<double>::infinity();
    double right_max = -std::numeric_limits<double>::infinity();

    for (const auto& s : barcode_stats) {
        poll_iteration();
        if (s.log1p_ncpm < FLOOR) continue;
        af_total++;
        if (s.log1p_ncpm < threshold) {
            af_left++;
            left_min  = std::min(left_min,  s.log1p_ncpm);
            left_max  = std::max(left_max,  s.log1p_ncpm);
        } else {
            af_right++;
            right_min = std::min(right_min, s.log1p_ncpm);
            right_max = std::max(right_max, s.log1p_ncpm);
        }
    }

    const double left_frac   = (af_total > 0) ? (static_cast<double>(af_left) / af_total) : 1.0;
    const bool   have_left   = (af_left  >= MIN_SIDE_N);
    const bool   have_right  = (af_right >= MIN_SIDE_N);
    const double left_range  = (have_left  ? (left_max  - left_min)  : 0.0);
    const double right_range = (have_right ? (right_max - right_min) : 0.0);
    const double full_range  = (std::isfinite(left_min) && std::isfinite(right_max))
                                 ? (right_max - left_min)
                                 : 0.0;

    const bool narrow_both   = have_left && have_right &&
                               (left_range <= SPREAD_EPS) && (right_range <= SPREAD_EPS);
    const bool narrow_full   = (full_range > 0.0) && (full_range <= 2.0 * SPREAD_EPS);
    const bool collapse_all_af =
        !force_final_hazard &&
        (force_all_af_single_peak ||
         force_all_af_adaptive_hazard ||
         ((af_total >= 10) &&
          ((left_frac <= LEFT_TAIL_MAX_FRAC) ||
           narrow_both || narrow_full)));

    // Budget rescue: only fires when collapse didn't, and left_frac is in the
    // borderline (LEFT_TAIL_MAX_FRAC, LEFT_TAIL_BUDGET_MAX_FRAC] window.
    const bool use_left_budget = !force_final_hazard
                              && (af_total >= 10) && !collapse_all_af
                              && (left_frac >  LEFT_TAIL_MAX_FRAC)
                              && (left_frac <= LEFT_TAIL_BUDGET_MAX_FRAC);
    std::unordered_set<std::string> budget_rescue;
    size_t budget_n = 0;
    if (use_left_budget) {
        std::vector<size_t> af_left_idx;
        af_left_idx.reserve(af_left);
        for (size_t i = 0; i < barcode_stats.size(); ++i) {
            poll_iteration();
            const auto& s = barcode_stats[i];
            if (s.log1p_ncpm >= FLOOR && s.log1p_ncpm < threshold) {
                af_left_idx.push_back(i);
            }
        }
        budget_n = static_cast<size_t>(LEFT_TAIL_MAX_FRAC * static_cast<double>(af_total));
        if (budget_n > af_left_idx.size()) budget_n = af_left_idx.size();
        if (budget_n > 0) {
            poll_now();
            std::partial_sort(
                af_left_idx.begin(),
                af_left_idx.begin() + budget_n,
                af_left_idx.end(),
                [&](size_t a, size_t b) {
                    if (barcode_stats[a].count != barcode_stats[b].count) {
                        return barcode_stats[a].count > barcode_stats[b].count;
                    }
                    return barcode_stats[a].sequence <
                           barcode_stats[b].sequence;
                });
            poll_now();
            for (size_t i = 0; i < budget_n; ++i) {
                poll_iteration();
                budget_rescue.insert(barcode_stats[af_left_idx[i]].sequence);
            }
        }
    }

    RAD_COUT << "\n[AF split] floor=" << FLOOR
              << "  threshold=" << threshold
              << "  rule=" << rule
              << "\n         AF_total=" << af_total
              << "  AF_left=" << af_left << " (" << left_frac * 100.0 << "%)"
              << "  AF_right=" << af_right
              << "\n[Spread]  left_range="  << left_range
              << "  right_range=" << right_range
              << "  full_range="  << full_range
              << "\n[Collapse] "
              << (collapse_all_af
                    ? "YES → return ALL AF"
                    : (use_left_budget
                          ? ("NO; LEFT_BUDGET fallback → +"
                             + std::to_string(budget_n)
                             + " AF_left rescues (top "
                             + std::to_string(static_cast<int>(LEFT_TAIL_MAX_FRAC * 100.0))
                             + "% of AF_total by count)")
                          : "NO"))
              << "\n";

    // 10) Write CSV — explicit outputs + annotation band
    //     final_barcode = TRUE iff:
    //       - collapse_all_af: (log1p_ncpm >= FLOOR)
    //       - else:            (log1p_ncpm >= threshold)
    std::ofstream csv_out(output_csv);
    if (!csv_out) {
        RAD_CERR << "Error: could not open statistics CSV for writing: "
                  << output_csv << "\n";
        return false;
    }
    std::vector<std::string> selected_barcodes;

    csv_out << "barcode,count,ncpm,log1p_ncpm,ztpois_percentile,above_floor,over_threshold,"
               "final_barcode,final_bc_annotation\n";

    size_t n_final = 0, n_above_floor = 0;
    for (const auto& s : barcode_stats) {
        poll_iteration();
        const bool above_floor   = (s.log1p_ncpm >= FLOOR);
        if (above_floor) n_above_floor++;

        const bool over_threshold = (s.log1p_ncpm >= threshold);
        const bool budget_pass    = use_left_budget &&
                                    above_floor && !over_threshold &&
                                    (budget_rescue.find(s.sequence) != budget_rescue.end());
        const bool final_barcode  = collapse_all_af
                                      ? above_floor
                                      : (over_threshold || budget_pass);
        if (final_barcode) n_final++;

        // Annotation:
        // - collapse mode: AF that pass are "high_confidence"
        // - otherwise:
        //     * ≥ chosen threshold: "high_confidence"
        //     * budget rescue (10-20% borderline window): "high_sensitivity"
        //     * between t_saddle and t_purity50: "high_sensitivity"
        std::string ann;
        if (collapse_all_af) {
            if (final_barcode){
                ann = "high_confidence";
            }
        } else {
            if (over_threshold && above_floor) {
                ann = "high_confidence";
            } else if (budget_pass) {
                ann = "high_sensitivity";
            } else if (have_band && above_floor && s.log1p_ncpm >= band_lo && s.log1p_ncpm < band_hi) {
                ann = "high_sensitivity";
            } else {
                ann = "low_confidence";
            }
        }

        const bool selected =
            selection == scan_whitelist_selection::above_floor
                ? above_floor
                : final_barcode;
        if (selected) {
            selected_barcodes.push_back(s.sequence);
        }

        csv_out << s.sequence << ","
                << s.count    << ","
                << std::fixed << std::setprecision(6) << s.ncpm << ","
                << s.log1p_ncpm << ","
                << s.log1p_ncpm_ztpois << ","
                << (above_floor   ? "TRUE" : "FALSE") << ","
                << (over_threshold? "TRUE" : "FALSE") << ","
                << (final_barcode ? "TRUE" : "FALSE") << ","
                << ann
                << "\n";
    }
    csv_out.close();
    if (!csv_out) {
        RAD_CERR << "Error: failed while writing statistics CSV: "
                  << output_csv << "\n";
        return false;
    }

    // 11) Console summary
    RAD_COUT << "\nStatistical analysis complete!\n"
              << "Total extracted sequences: " << extracted_barcodes.total_extractions << "\n";
    if (extracted_barcodes.filtered_to_whitelist) {
        RAD_COUT << "Unique whitelist-matching sequences retained: "
                  << extracted_counts.size() << "\n";
    } else {
        RAD_COUT << "Unique extracted sequences: " << extracted_counts.size() << "\n";
    }
    RAD_COUT << "Total perfect matches: " << total_perfect_matches << "\n"
              << "Unique barcodes with perfect matches: " << unique_matches << "\n"
              << "Match rate: " << std::fixed << std::setprecision(2)
              << (100.0 * static_cast<double>(total_perfect_matches) /
                  std::max<double>(1.0, static_cast<double>(extracted_barcodes.total_extractions))) << "%\n"
              << "Above-floor barcodes: " << n_above_floor << "\n"
              << "Final barcodes (TRUE): " << n_final << "\n"
              << "Selected barcodes (" << scan_whitelist_selection_name(selection)
              << "): " << selected_barcodes.size() << "\n"
              << "Results written to: " << output_csv << "\n";

    if (summary) {
        summary->reads_processed = extracted_barcodes.reads_processed;
        summary->total_extractions = extracted_barcodes.total_extractions;
        summary->unique_sequences = extracted_counts.size();
        summary->total_perfect_matches = total_perfect_matches;
        summary->unique_perfect_matches = unique_matches;
        summary->match_rate_percent =
            100.0 * static_cast<double>(total_perfect_matches) /
            std::max<double>(
                1.0,
                static_cast<double>(extracted_barcodes.total_extractions));
        summary->floor = FLOOR;
        summary->threshold = threshold;
        summary->threshold_rule = rule;
        summary->above_floor_barcodes = n_above_floor;
        summary->final_barcodes = n_final;
        summary->selected_barcodes = selected_barcodes.size();
        summary->selection = selection;
    }

    if (!text_out.empty()) {
        std::ofstream hc_out(text_out);
        if (!hc_out) {
            RAD_CERR << "Error: could not open whitelist output for writing: "
                      << text_out << "\n";
            return false;
        }
        for (const auto& bc : selected_barcodes) {
            poll_iteration();
            hc_out << bc << "\n";
        }
        hc_out.close();
        if (!hc_out) {
            RAD_CERR << "Error: failed while writing whitelist output: "
                      << text_out << "\n";
            return false;
        }
        RAD_COUT << "Selected whitelist written to: " << text_out << "\n";
    }
    poll_now();
    return true;
}

// Two-part barcode mode still needs the observed BC1/BC2 pairs for its pair
// matrix outputs.  Collapse its combined sequences before entering the shared
// statistics path and, critically, avoid reserving hash buckets for every read.
bool count_perfect_matches_with_stats(const split_barcode_count_result& extracted_barcodes,
                                      const std::unordered_set<int64_seq>& whitelist_set,
                                      const std::string& output_csv, const std::string text_out,
                                      bool verbose = false,
                                      scan_whitelist_selection selection =
                                          scan_whitelist_selection::high_specificity,
                                      whitelist_scan_stats* summary = nullptr,
                                      const std::function<void()>& boundary_poll = {})
{
    barcode_count_result counts;
    counts.reads_processed = extracted_barcodes.reads_processed;
    counts.total_extractions = extracted_barcodes.total_extractions;
    counts.counts.reserve(
        std::min<size_t>(extracted_barcodes.pair_counts.size(), 1u << 20));
    size_t pair_index = 0;
    if (boundary_poll) boundary_poll();
    for (const auto& pair : extracted_barcodes.pair_counts) {
        if (boundary_poll && (++pair_index % 4096) == 0) boundary_poll();
        const auto tab = pair.first.find('\t');
        if (tab == std::string::npos) continue;
        const std::string combined =
            pair.first.substr(0, tab) + pair.first.substr(tab + 1);
        counts.counts[combined] += pair.second;
    }
    return count_perfect_matches_with_stats(
        counts, whitelist_set, output_csv, text_out, verbose, selection,
        summary, boundary_poll);
}

struct scan_wl_rescan_input {
    barcode_count_result counts;
    uint64_t input_rows = 0;
};

inline std::string scan_wl_rescan_trim(const std::string& value) {
    const auto first = std::find_if_not(
        value.begin(), value.end(),
        [](unsigned char ch) { return std::isspace(ch); });
    if (first == value.end()) {
        return "";
    }
    const auto last = std::find_if_not(
                          value.rbegin(), value.rend(),
                          [](unsigned char ch) { return std::isspace(ch); })
                          .base();
    return std::string(first, last);
}

inline std::string scan_wl_rescan_lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char ch) {
                       return static_cast<char>(std::tolower(ch));
                   });
    return value;
}

inline bool scan_wl_rescan_ends_with(const std::string& value,
                                     const std::string& suffix) {
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) ==
               0;
}

inline std::vector<std::string> scan_wl_rescan_split_delimited(
    const std::string& line, char delimiter) {
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (size_t i = 0; i < line.size(); ++i) {
        const char ch = line[i];
        if (quoted) {
            if (ch == '"') {
                if (i + 1 < line.size() && line[i + 1] == '"') {
                    field.push_back('"');
                    ++i;
                } else {
                    quoted = false;
                }
            } else {
                field.push_back(ch);
            }
        } else if (ch == '"') {
            quoted = true;
        } else if (ch == delimiter) {
            fields.push_back(scan_wl_rescan_trim(field));
            field.clear();
        } else if (ch != '\r') {
            field.push_back(ch);
        }
    }
    if (quoted) {
        throw std::runtime_error("unterminated quoted field");
    }
    fields.push_back(scan_wl_rescan_trim(field));
    return fields;
}

inline std::vector<std::string> scan_wl_rescan_split(
    const std::string& line, char& delimiter) {
    if (delimiter == '\0') {
        if (line.find(',') != std::string::npos) {
            delimiter = ',';
        } else if (line.find('\t') != std::string::npos) {
            delimiter = '\t';
        } else {
            delimiter = ' ';
        }
    }
    if (delimiter != ' ') {
        return scan_wl_rescan_split_delimited(line, delimiter);
    }

    std::vector<std::string> fields;
    std::istringstream input(line);
    std::string field;
    while (input >> field) {
        fields.push_back(std::move(field));
    }
    return fields;
}

inline uint64_t scan_wl_rescan_parse_count(const std::string& text,
                                           size_t line_number) {
    if (text.empty() || text.front() == '-') {
        throw std::runtime_error(
            "line " + std::to_string(line_number) +
            " has an invalid positive count: " + text);
    }
    errno = 0;
    char* end = nullptr;
    const unsigned long long value = std::strtoull(text.c_str(), &end, 10);
    if (errno == ERANGE || end == text.c_str() || *end != '\0' || value == 0) {
        throw std::runtime_error(
            "line " + std::to_string(line_number) +
            " has an invalid positive count: " + text);
    }
    return static_cast<uint64_t>(value);
}

inline scan_wl_rescan_input read_scan_wl_rescan_counts(
    const std::string& input_path,
    const std::function<void()>& boundary_poll = {}) {
    std::unique_ptr<std::istream> input;
    igzstream* gzip_input = nullptr;
    if (scan_wl_rescan_ends_with(
            scan_wl_rescan_lower(input_path), ".gz")) {
        auto stream = std::make_unique<igzstream>(input_path.c_str());
        gzip_input = stream.get();
        input = std::move(stream);
    } else {
        input = std::make_unique<std::ifstream>(input_path, std::ios::binary);
    }
    if (!input || !*input) {
        throw std::runtime_error("could not open rescan input: " + input_path);
    }

    scan_wl_rescan_input result;
    result.counts.filtered_to_whitelist = true;
    char delimiter = '\0';
    bool columns_known = false;
    size_t barcode_column = 0;
    size_t count_column = 1;
    size_t line_number = 0;
    std::string line;

    if (boundary_poll) boundary_poll();
    while (std::getline(*input, line)) {
        ++line_number;
        if (boundary_poll && (line_number % 4096) == 0) boundary_poll();
        const std::string trimmed = scan_wl_rescan_trim(line);
        if (trimmed.empty() || trimmed.front() == '#') {
            continue;
        }

        std::vector<std::string> fields;
        try {
            fields = scan_wl_rescan_split(trimmed, delimiter);
        } catch (const std::exception& error) {
            throw std::runtime_error(
                "line " + std::to_string(line_number) + ": " + error.what());
        }
        if (fields.size() < 2) {
            throw std::runtime_error(
                "line " + std::to_string(line_number) +
                " has fewer than two columns");
        }

        if (!columns_known) {
            if (fields.front().size() >= 3 &&
                static_cast<unsigned char>(fields.front()[0]) == 0xef &&
                static_cast<unsigned char>(fields.front()[1]) == 0xbb &&
                static_cast<unsigned char>(fields.front()[2]) == 0xbf) {
                fields.front().erase(0, 3);
            }
            size_t named_barcode = fields.size();
            size_t named_count = fields.size();
            for (size_t i = 0; i < fields.size(); ++i) {
                const std::string name = scan_wl_rescan_lower(fields[i]);
                if (name == "barcode") {
                    named_barcode = i;
                } else if (name == "count") {
                    named_count = i;
                }
            }
            if (named_barcode != fields.size() || named_count != fields.size()) {
                if (named_barcode == fields.size() || named_count == fields.size()) {
                    throw std::runtime_error(
                        "header must contain both barcode and count columns");
                }
                barcode_column = named_barcode;
                count_column = named_count;
                columns_known = true;
                continue;
            }
            columns_known = true;
        }

        const size_t required_column = std::max(barcode_column, count_column);
        if (fields.size() <= required_column) {
            throw std::runtime_error(
                "line " + std::to_string(line_number) +
                " is missing a barcode or count value");
        }
        const std::string barcode = scan_wl_rescan_trim(fields[barcode_column]);
        if (barcode.empty()) {
            throw std::runtime_error(
                "line " + std::to_string(line_number) +
                " has an empty barcode");
        }
        const uint64_t count =
            scan_wl_rescan_parse_count(fields[count_column], line_number);
        uint64_t& aggregate = result.counts.counts[barcode];
        if (count > std::numeric_limits<uint64_t>::max() - aggregate ||
            count > std::numeric_limits<uint64_t>::max() -
                        result.counts.total_extractions) {
            throw std::runtime_error(
                "count overflow at line " + std::to_string(line_number));
        }
        aggregate += count;
        result.counts.total_extractions += count;
        ++result.input_rows;
    }

    if (gzip_input) {
        const int zlib_code = gzip_input->zlib_error_code();
        const std::string zlib_message = gzip_input->zlib_error_message();
        gzip_input->close();
        if (zlib_code != Z_OK && zlib_code != Z_STREAM_END) {
            throw std::runtime_error(
                "rescan gzip read failed: " + input_path +
                " (zlib code " + std::to_string(zlib_code) +
                (zlib_message.empty() ? ")"
                                      : ", " + zlib_message + ")"));
        }
        if (gzip_input->bad()) {
            throw std::runtime_error(
                "rescan gzip close failed: " + input_path);
        }
    } else if (!input->eof()) {
        throw std::runtime_error(
            "failed while reading rescan input: " + input_path);
    }
    if (result.counts.counts.empty()) {
        throw std::runtime_error(
            "rescan input contains no barcode/count rows: " + input_path);
    }
    result.counts.reference_whitelist_size =
        static_cast<uint64_t>(result.counts.counts.size());
    if (boundary_poll) boundary_poll();
    return result;
}

inline std::string default_scan_wl_rescan_prefix(std::string input_path) {
    if (scan_wl_rescan_ends_with(
            scan_wl_rescan_lower(input_path), ".gz")) {
        input_path.resize(input_path.size() - 3);
    }
    const std::string::size_type slash = input_path.find_last_of("/\\");
    const std::string::size_type dot = input_path.find_last_of('.');
    if (dot != std::string::npos &&
        (slash == std::string::npos || dot > slash)) {
        input_path.resize(dot);
    }
    return input_path + "_rescan";
}


void usage_scan_wl(const char* program_name) {
    RAD_COUT << "Usage: " << program_name << " [options]\n"
              << "Options:\n"
              << "  -i, --input FILE            Input FASTQ file (use this for single runs)\n"
              << "  -b, --batch-csv FILE        CSV with FASTQ path and output prefix (two columns, optional)\n"
              << "  -o, --output-prefix PREFIX  Output prefix for generated files (.txt and .csv)\n"
              << "  -p, --adapter_seq SEQ       Adapter/primer sequence to search for\n"
              << "  -n, --barcode-length INT    Number of bases to extract (barcode length)\n"
              << "  -l, --left-margin INT       Bases to include on left side [default: 0]\n"
              << "  -r, --right-margin INT      Bases to include on right side [default: 0]\n"
              << "  -m, --max-reads INT         Maximum number of reads to process [default: all]\n"
              << "  -e, --max-error FLOAT       Maximum edit distance ratio [default: "
              << kDefaultScanWlMaxErrorRatio << "]\n"
              << "  -w, --whitelist FILE        Barcode whitelist kit key or path (optional)\n"
              << "  -t, --threads INT           Number of threads for parallel processing [default: auto]\n"
              << "  -k, --chunk-size INT        Chunk size for parallel processing [default: 10000]\n"
              << "      --af-bcs                Write above-floor barcodes to .txt\n"
              << "      --hs-bcs                Write high-specificity final calls to .txt [default]\n"
              << "      --rescan                Treat --input as barcode/count data and rerun only\n"
              << "                              the statistical selector (no FASTQ scan)\n"
              << "  -v, --verbose               Enable verbose/debug output\n"
              << "  -h, --help                  Show this help message\n"
              << "\nTwo-part barcode options (BC1+BC2 mode):\n"
              << "  -1, --bc1-whitelist FILE    BC1 whitelist kit key or path (activates two-part mode)\n"
              << "  -2, --bc2-whitelist FILE    BC2 whitelist kit key or path\n"
              << "  -u, --umi-length INT        UMI length between primer and BC1 [default: 9]\n"
              << "      --offset-min INT        Min extra bases after UMI to start BC1 search [default: 0]\n"
              << "      --offset-max INT        Max extra bases after UMI to start BC1 search [default: 3]\n";
}

inline int cmd_scan_wl(int argc, char* argv[]) {
    optind = 1; // reset getopt state for subcommand dispatch
    std::string input_file, output_prefix, adapter_seq, whitelist_file, batch_csv_file;
    std::string bc1_whitelist_file, bc2_whitelist_file;
    int bc_length = 16, m_left = 0, m_right = 0, max_reads = 0, num_threads = 0, chunk_size = 10000;
    int umi_length = 9, offset_min = 0, offset_max = 3;
    double max_error = kDefaultScanWlMaxErrorRatio;
    bool verbose = false, af_bcs = false, hs_bcs = false, rescan = false;

    static struct option long_options[] = {
        {"input",          required_argument, 0, 'i'},
        {"batch-csv",      required_argument, 0, 'b'},
        {"output-prefix",  required_argument, 0, 'o'},
        {"adapter_seq",    required_argument, 0, 'p'},
        {"barcode-length", required_argument, 0, 'n'},
        {"left-margin",    required_argument, 0, 'l'},
        {"right-margin",   required_argument, 0, 'r'},
        {"max-reads",      required_argument, 0, 'm'},
        {"max-error",      required_argument, 0, 'e'},
        {"whitelist",      required_argument, 0, 'w'},
        {"threads",        required_argument, 0, 't'},
        {"chunk-size",     required_argument, 0, 'k'},
        {"bc1-whitelist",  required_argument, 0, '1'},
        {"bc2-whitelist",  required_argument, 0, '2'},
        {"umi-length",     required_argument, 0, 'u'},
        {"offset-min",     required_argument, 0,  3 },
        {"offset-max",     required_argument, 0,  4 },
        {"af-bcs",          no_argument,       0,  5 },
        {"hs-bcs",          no_argument,       0,  6 },
        {"rescan",          no_argument,       0,  7 },
        // Backward-compatible aliases
        {"hd-offset-min",  required_argument, 0,  3 },
        {"hd-offset-max",  required_argument, 0,  4 },
        {"verbose",        no_argument,       0, 'v'},
        {"help",           no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "i:b:o:p:n:l:r:m:e:w:t:k:1:2:u:hv", long_options, NULL)) != -1) {
        switch (opt) {
            case 'i': input_file = optarg; break;
            case 'b': batch_csv_file = optarg; break;
            case 'o': output_prefix = optarg; break;
            case 'p': adapter_seq = optarg; break;
            case 'n': bc_length = std::atoi(optarg); break;
            case 'l': m_left = std::atoi(optarg); break;
            case 'r': m_right = std::atoi(optarg); break;
            case 'm': max_reads = std::atoi(optarg); break;
            case 'e': max_error = std::atof(optarg); break;
            case 'w': whitelist_file = optarg; break;
            case 't': num_threads = std::atoi(optarg); break;
            case 'k': chunk_size = std::atoi(optarg); break;
            case '1': bc1_whitelist_file = optarg; break;
            case '2': bc2_whitelist_file = optarg; break;
            case 'u': umi_length = std::atoi(optarg); break;
            case  3 : offset_min = std::atoi(optarg); break;
            case  4 : offset_max = std::atoi(optarg); break;
            case  5 : af_bcs = true; break;
            case  6 : hs_bcs = true; break;
            case  7 : rescan = true; break;
            case 'v': verbose = true; break;
            case 'h': usage_scan_wl(argv[0]); return 0;
            default: usage_scan_wl(argv[0]); return 1;
        }
    }
    
    if (!rescan && adapter_seq.empty()) {
        RAD_CERR << "Error: Missing required adapter sequence\n";
        usage_scan_wl(argv[0]);
        return 1;
    }
    if (chunk_size <= 0) {
        RAD_CERR << "Error: Chunk size must be greater than zero\n";
        return 1;
    }
    if (af_bcs && hs_bcs) {
        RAD_CERR << "Error: --af-bcs and --hs-bcs are mutually exclusive\n";
        return 1;
    }

    const bool using_single = !input_file.empty();
    const bool using_batch  = !batch_csv_file.empty();

    if (!using_single && !using_batch) {
        RAD_CERR << "Error: Provide either --input (single run) or --batch-csv (batch mode)\n";
        usage_scan_wl(argv[0]);
        return 1;
    }

    if (using_batch && using_single) {
        RAD_COUT << "Batch CSV provided - ignoring single --input value\n";
    }

    if (rescan) {
        if (using_batch) {
            RAD_CERR << "Error: --rescan does not support --batch-csv\n";
            return 1;
        }
        if (!bc1_whitelist_file.empty() || !bc2_whitelist_file.empty()) {
            RAD_CERR << "Error: --rescan is not compatible with split-barcode options\n";
            return 1;
        }

        const auto rescan_started = std::chrono::steady_clock::now();
        const std::string resolved_prefix =
            output_prefix.empty()
                ? default_scan_wl_rescan_prefix(input_file)
                : output_prefix;
        if (!ensure_output_directory(resolved_prefix)) {
            return 1;
        }
        const std::string csv_output = resolved_prefix + ".csv";
        const std::string text_output = resolved_prefix + ".txt";
        const std::string scan_log_output =
            resolved_prefix + "_scan_wl.log";
        if (csv_output == input_file || text_output == input_file) {
            RAD_CERR << "Error: --rescan output would overwrite its input\n";
            return 1;
        }

        scan_wl_rescan_input retained;
        try {
            retained = read_scan_wl_rescan_counts(input_file);
        } catch (const std::exception& error) {
            RAD_CERR << "Error: " << error.what() << "\n";
            return 1;
        }
        RAD_COUT << "Rescanning " << retained.input_rows
                  << " barcode/count rows ("
                  << retained.counts.counts.size() << " unique barcodes)\n"
                  << "Total extractions from count sum: "
                  << retained.counts.total_extractions << "\n"
                  << "Whitelist universe from unique input barcodes: "
                  << retained.counts.reference_whitelist_size << "\n";

        whitelist_scan_stats scan_stats;
        const scan_whitelist_selection selection =
            af_bcs ? scan_whitelist_selection::above_floor
                   : scan_whitelist_selection::high_specificity;
        const std::unordered_set<int64_seq> no_secondary_filter;
        if (!count_perfect_matches_with_stats(
                retained.counts, no_secondary_filter, csv_output, text_output,
                verbose, selection, &scan_stats)) {
            return 1;
        }
        const double rescan_wall_seconds =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - rescan_started)
                .count();
        if (!write_whitelist_scan_log(
                scan_log_output, input_file, scan_stats,
                rescan_wall_seconds)) {
            return 1;
        }
        RAD_COUT << "Whitelist rescan summary written to: "
                  << scan_log_output << "\n";
        return 0;
    }

    const bool split_barcode_mode = !bc1_whitelist_file.empty() && !bc2_whitelist_file.empty();
    bool run_failed = false;

    auto run_single = [&](const std::string& fastq_path, const std::string& output_prefix) {
        const auto scan_started = std::chrono::steady_clock::now();
        std::string resolved_prefix = output_prefix.empty() ? "barcodes" : output_prefix;
        if (!ensure_output_directory(resolved_prefix)) {
            run_failed = true;
            return;
        }

        std::string csv_output  = resolved_prefix + ".csv";
        std::string text_output = resolved_prefix + ".txt";
        std::string scan_log_output = resolved_prefix + "_scan_wl.log";
        whitelist_scan_stats scan_stats;
        const scan_whitelist_selection selection =
            af_bcs
                ? scan_whitelist_selection::above_floor
                : scan_whitelist_selection::high_specificity;

        if (split_barcode_mode) {
            std::string bc1_path, bc2_path;
            try {
                bc1_path = whitelist_utils::kit_to_path(bc1_whitelist_file);
                bc2_path = whitelist_utils::kit_to_path(bc2_whitelist_file);
            } catch (const std::exception& e) {
                RAD_CERR << "Error: could not resolve split whitelist path(s): " << e.what() << "\n";
                run_failed = true;
                return;
            }

            RAD_COUT << "Two-part barcode mode: " << fastq_path << "\n"
                      << "  Adapter: " << adapter_seq << "  UMI length: " << umi_length
                      << "  BC1 search offset: [" << offset_min << ", " << offset_max << "]\n";
            if (verbose) {
                RAD_COUT << "  BC1 whitelist: " << bc1_path << "\n"
                          << "  BC2 whitelist: " << bc2_path << "\n";
            }

            split_whitelist bc1_wl;
            split_whitelist bc2_wl;
            try {
                bc1_wl = load_split_whitelist_csv(bc1_path);
                bc2_wl = load_split_whitelist_csv(bc2_path);
            } catch (const std::exception& e) {
                RAD_CERR << "Error: failed to load two-part barcode "
                             "whitelists: " << e.what() << "\n";
                run_failed = true;
                return;
            }
            RAD_COUT << "Loaded " << bc1_wl.seqs.size() << " BC1 / " << bc2_wl.seqs.size() << " BC2 sequences\n";

            split_barcode_count_result barcodes;
            try {
                barcodes = process_fastq_split_barcode(
                    fastq_path, adapter_seq,
                    umi_length, bc1_wl.seqs, bc1_wl.lengths,
                    bc2_wl.seqs, bc2_wl.lengths,
                    offset_min, offset_max,
                    max_reads, max_error, chunk_size, num_threads);
            } catch (const std::exception& e) {
                RAD_CERR << "Error: split-barcode FASTQ scan failed: "
                          << e.what() << "\n";
                run_failed = true;
                return;
            }

            std::string pair_counts_output = resolved_prefix + "_valid_pairs.csv";
            if (!write_split_pair_counts(barcodes, pair_counts_output)) {
                run_failed = true;
                return;
            }
            RAD_COUT << "Raw valid BC1+BC2 pair counts written to: " << pair_counts_output << "\n";

            std::string spat_mask_output = resolved_prefix + "_spat_mask.csv";
            if (!write_spat_mask_csv(
                    barcodes, bc1_wl, bc2_wl, spat_mask_output)) {
                run_failed = true;
                return;
            }

            std::unordered_set<int64_seq> empty_wl;
            if (!count_perfect_matches_with_stats(
                    barcodes, empty_wl, csv_output, text_output, verbose,
                    selection, &scan_stats)) {
                run_failed = true;
                return;
            }
        } else {
            RAD_COUT << "Processing " << fastq_path << "...\n";
            RAD_COUT << "Adapter: " << adapter_seq << "\n";
            RAD_COUT << "Extracting " << bc_length << " bases with margins: left=" << m_left << ", right=" << m_right << "\n";
            RAD_COUT << "Chunk size: " << chunk_size << " reads per chunk\n";

            scan_whitelist_filter wl_filter;
            if (!whitelist_file.empty()) {
                std::string resolved_wl;
                try {
                    resolved_wl = whitelist_utils::kit_to_path(whitelist_file);
                    RAD_COUT << "\nLoading whitelist for streaming validation...\n";
                    if (verbose) {
                        RAD_COUT << "Whitelist source: " << whitelist_file << "\n"
                                  << "Resolved path: " << resolved_wl << "\n";
                    }
                    wl_filter = load_scan_whitelist(
                        resolved_wl, static_cast<uint16_t>(bc_length), verbose);
                } catch (const std::exception& e) {
                    RAD_CERR << "Error: could not load whitelist: "
                              << e.what() << "\n";
                    run_failed = true;
                    return;
                }
                RAD_COUT << "Loaded " << wl_filter.size()
                          << " whitelist barcodes\n";
                if (wl_filter.empty()) {
                    RAD_CERR << "Error: whitelist contains no valid barcodes: "
                              << resolved_wl << "\n";
                    run_failed = true;
                    return;
                }
            }

            auto barcodes = process_fastq(
                fastq_path, adapter_seq, bc_length, m_left, m_right,
                max_reads, max_error, chunk_size, num_threads,
                wl_filter.empty() ? nullptr : &wl_filter);
            if (!barcodes.succeeded) {
                run_failed = true;
                return;
            }

            RAD_COUT << "Found " << barcodes.size() << " barcodes\n";

            if (!whitelist_file.empty()) {
                if (!count_perfect_matches_with_stats(
                        barcodes, wl_filter.long_barcodes,
                        csv_output, text_output, verbose, selection,
                        &scan_stats)) {
                    run_failed = true;
                    return;
                }
            } else {
                std::unordered_set<int64_seq> empty_whitelist;
                RAD_COUT << " No whitelist provided; generating whitelist de novo.\n";
                if (!count_perfect_matches_with_stats(
                        barcodes, empty_whitelist,
                        csv_output, text_output, verbose, selection,
                        &scan_stats)) {
                    run_failed = true;
                    return;
                }
            }
        }

        const double scan_wall_seconds =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - scan_started)
                .count();
        if (!write_whitelist_scan_log(
                scan_log_output, fastq_path, scan_stats,
                scan_wall_seconds)) {
            run_failed = true;
            return;
        }
        RAD_COUT << "Whitelist scan summary written to: "
                  << scan_log_output << "\n";
    };

    if (using_batch) {
        auto batch_entries = read_batch_csv(batch_csv_file);
        if (batch_entries.empty()) {
            return 1;
        }
        RAD_COUT << "Running batch of " << batch_entries.size() << " FASTQ files from " << batch_csv_file << "\n";
        for (size_t idx = 0; idx < batch_entries.size(); ++idx) {
            const auto& entry = batch_entries[idx];
            RAD_COUT << "\n[Batch " << (idx + 1) << "/" << batch_entries.size() << "] Prefix: "
                      << entry.output_prefix << "\n";
            run_single(entry.input_fastq, entry.output_prefix);
        }
    } else {
        run_single(input_file, output_prefix);
    }
    return run_failed ? 1 : 0;
}
