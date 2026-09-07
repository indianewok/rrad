#pragma once

#include "rad_headers.h"

// Reusable, in-process implementation of RAD's `reformat` command.  The
// public entry point intentionally accepts typed options rather than CLI
// arguments so hosts such as R can validate their own interfaces and call the
// engine without constructing argv or starting a subprocess.

struct reformat_options {
    std::string input_path;
    std::string split_output_dir;
    bool split_by_barcode = false;
    bool reformat_header = false;
    bool presto_header = false;
    bool to_fasta = false;
    bool parse_collapsed_id = false;
    char delimiter = '_';
    std::optional<std::string> coordinate_mode;
    int bin_size_um = 2;
    std::string output_path;
    int threads = 1;
    size_t chunk_size = 5000;
    bool verbose = false;
    bool overwrite = false;
    size_t max_open_split_writers = 64;
};

struct reformat_stats {
    uint64_t records_read = 0;
    uint64_t records_written = 0;
    uint64_t records_skipped_missing_cb = 0;
    uint64_t records_reformatted = 0;
    uint64_t records_converted_to_fasta = 0;
    uint64_t coordinate_mapped = 0;
    uint64_t coordinate_unmapped = 0;
    uint64_t empty_sequences_skipped = 0;
    uint64_t split_barcodes = 0;
    uint64_t gzip_members_written = 0;
};

struct reformat_result {
    std::string input_path;
    std::string output_path;
    std::string split_output_dir;
    std::vector<std::pair<std::string, std::string>> split_paths;
    reformat_stats stats;
};

static_assert(std::is_nothrow_move_constructible<reformat_result>::value,
              "reformat_result must remain nothrow-movable after publication");

namespace reformat_detail {

namespace fs = std::filesystem;

inline bool ends_with(const std::string& value, const std::string& suffix) {
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

inline std::string lower_copy(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) {
                       return static_cast<char>(std::tolower(c));
                   });
    return value;
}

inline bool has_fasta_suffix(const std::string& path) {
    return ends_with(path, ".fa") || ends_with(path, ".fasta") ||
           ends_with(path, ".fa.gz") || ends_with(path, ".fasta.gz");
}

inline std::optional<std::string> extract_tag_value(
    const std::string& source, const std::string& tag) {
    size_t start = 0;
    while (start < source.size()) {
        start = source.find_first_not_of(" \t", start);
        if (start == std::string::npos) break;
        size_t end = source.find_first_of(" \t", start);
        if (end == std::string::npos) end = source.size();
        const size_t token_size = end - start;
        if (token_size >= tag.size() &&
            source.compare(start, tag.size(), tag) == 0) {
            return source.substr(start + tag.size(), token_size - tag.size());
        }
        start = end;
    }
    return std::nullopt;
}

inline std::optional<std::string> extract_cb(const std::string& id,
                                             const std::string& comment) {
    static const std::string tag = "CB:Z:";
    if (auto value = extract_tag_value(comment, tag)) return value;
    return extract_tag_value(id, tag);
}

inline std::optional<std::string> extract_ub(const std::string& comment) {
    static const std::string tag = "UB:Z:";
    return extract_tag_value(comment, tag);
}

inline std::optional<std::string> extract_bc(const std::string& comment) {
    static const std::string tag = "BC:Z:";
    return extract_tag_value(comment, tag);
}

inline std::string collapse_id(const std::string& qname,
                               const std::optional<std::string>& cb,
                               const std::optional<std::string>& ub,
                               char delimiter) {
    std::string output = qname;
    if (cb && !cb->empty()) {
        output.push_back(delimiter);
        output.append(*cb);
    }
    if (ub && !ub->empty()) {
        output.push_back(delimiter);
        output.append(*ub);
    }
    return output;
}

inline void validate_presto_component(const std::string& value,
                                      const char* label) {
    const bool unsafe = std::any_of(
        value.begin(), value.end(), [](unsigned char character) {
            return std::isspace(character) || character == '|' ||
                   character == '=' || character == ',';
        });
    if (unsafe) {
        throw std::invalid_argument(
            std::string("pRESTO ") + label +
            " contains whitespace or a reserved delimiter (|, =, or ,)");
    }
}

inline std::string presto_id(const std::string& qname,
                             const std::optional<std::string>& barcode,
                             const std::optional<std::string>& ub) {
    validate_presto_component(qname, "QNAME");
    std::string output = qname;
    if (barcode && !barcode->empty()) {
        validate_presto_component(*barcode, "BARCODE value");
        output += "|BARCODE=";
        output += *barcode;
    }
    if (ub && !ub->empty()) {
        validate_presto_component(*ub, "UB value");
        output += "|UMI=";
        output += *ub;
    }
    return output;
}

inline std::optional<std::tuple<std::string, std::string, std::string>>
parse_collapsed_id(const std::string& id, char delimiter) {
    // Parse from the right so delimiters in the original QNAME are retained.
    const auto second = id.rfind(delimiter);
    if (second == std::string::npos) return std::nullopt;
    const auto first =
        id.rfind(delimiter,
                 second == 0 ? std::string::npos : second - 1);
    if (first == std::string::npos) return std::nullopt;

    std::string qname = id.substr(0, first);
    std::string cb = id.substr(first + 1, second - first - 1);
    std::string ub = id.substr(second + 1);
    if (qname.empty() || cb.empty() || ub.empty()) return std::nullopt;
    return std::make_tuple(std::move(qname), std::move(cb), std::move(ub));
}

inline std::string visium_hd_comment(const std::string& spatial_barcode,
                                     const std::string& cb,
                                     const std::string& ub,
                                     const std::string& read_group,
                                     const std::string& original_comment) {
    std::string comment;
    comment.reserve(spatial_barcode.size() + cb.size() + ub.size() +
                    read_group.size() + original_comment.size() + 40);
    comment += "BC:Z:";
    comment += spatial_barcode;
    comment += "\tCB:Z:";
    comment += cb;
    comment += "\tUB:Z:";
    comment += ub;
    comment += "\tRG:Z:";
    comment += read_group;

    // Coordinate conversion owns these four SAM-style fields, but unrelated
    // metadata remains part of the read. Re-emit non-coordinate fields after
    // the canonical replacements instead of discarding the whole comment.
    std::istringstream fields(original_comment);
    std::string field;
    while (fields >> field) {
        const bool replaced =
            field.rfind("BC:", 0) == 0 || field.rfind("CB:", 0) == 0 ||
            field.rfind("UB:", 0) == 0 || field.rfind("RG:", 0) == 0;
        if (!replaced) {
            comment.push_back('\t');
            comment += field;
        }
    }
    return comment;
}

inline std::string zero_pad(int value, int width) {
    if (value < 0) value = 0;
    std::ostringstream output;
    output << std::setw(width) << std::setfill('0') << value;
    return output.str();
}

inline std::vector<std::string> split_string(const std::string& value,
                                             char delimiter) {
    std::vector<std::string> output;
    size_t start = 0;
    while (start <= value.size()) {
        size_t position = value.find(delimiter, start);
        if (position == std::string::npos) position = value.size();
        output.push_back(value.substr(start, position - start));
        start = position + 1;
        if (position == value.size()) break;
    }
    return output;
}

struct vizhd_axis_reference {
    std::unordered_map<std::string, int> x_index_by_bc;
    std::unordered_map<std::string, int> y_index_by_bc;
    std::vector<int> x_lengths;
    std::vector<int> y_lengths;
};

inline std::vector<std::string> load_barcode_csv_first_column(
    const std::string& path) {
    auto lines = streaming_utils::import_text(path);
    std::vector<std::string> output;
    output.reserve(lines.size());
    bool first = true;
    for (auto& line : lines) {
        auto trimmed = seq_utils::trim(line);
        if (trimmed.empty()) continue;
        if (first) {
            first = false;
            if (lower_copy(trimmed) == "whitelist_bcs") continue;
        }
        const auto comma = trimmed.find(',');
        if (comma != std::string::npos) trimmed = trimmed.substr(0, comma);
        trimmed = seq_utils::trim(trimmed);
        if (!trimmed.empty()) output.push_back(trimmed);
    }
    return output;
}

inline vizhd_axis_reference load_vizhd_v1_axis_reference() {
    vizhd_axis_reference reference;
    const auto x_barcodes = load_barcode_csv_first_column(
        whitelist_utils::kit_to_path("visium_hd_bc1"));
    const auto y_barcodes = load_barcode_csv_first_column(
        whitelist_utils::kit_to_path("visium_hd_bc2"));
    if (x_barcodes.empty() || y_barcodes.empty()) {
        throw std::runtime_error(
            "Visium HD coordinate barcode resources are empty");
    }

    reference.x_index_by_bc.reserve(x_barcodes.size());
    reference.y_index_by_bc.reserve(y_barcodes.size());
    for (size_t index = 0; index < x_barcodes.size(); ++index) {
        reference.x_index_by_bc.emplace(
            x_barcodes[index], static_cast<int>(index));
        reference.x_lengths.push_back(
            static_cast<int>(x_barcodes[index].size()));
    }
    for (size_t index = 0; index < y_barcodes.size(); ++index) {
        reference.y_index_by_bc.emplace(
            y_barcodes[index], static_cast<int>(index));
        reference.y_lengths.push_back(
            static_cast<int>(y_barcodes[index].size()));
    }

    std::sort(reference.x_lengths.begin(), reference.x_lengths.end());
    reference.x_lengths.erase(
        std::unique(reference.x_lengths.begin(), reference.x_lengths.end()),
        reference.x_lengths.end());
    std::sort(reference.y_lengths.begin(), reference.y_lengths.end());
    reference.y_lengths.erase(
        std::unique(reference.y_lengths.begin(), reference.y_lengths.end()),
        reference.y_lengths.end());
    return reference;
}

inline std::optional<std::pair<std::string, std::string>>
resolve_xy_barcodes(const std::string& cb_tag,
                    const vizhd_axis_reference& reference) {
    if (cb_tag.empty()) return std::nullopt;

    const auto direct = split_string(cb_tag, '-');
    if (direct.size() >= 2) {
        const std::string& first = direct[0];
        const std::string& second = direct[1];
        if (reference.x_index_by_bc.find(first) !=
                reference.x_index_by_bc.end() &&
            reference.y_index_by_bc.find(second) !=
                reference.y_index_by_bc.end()) {
            return std::make_pair(first, second);
        }
        if (reference.x_index_by_bc.find(second) !=
                reference.x_index_by_bc.end() &&
            reference.y_index_by_bc.find(first) !=
                reference.y_index_by_bc.end()) {
            return std::make_pair(second, first);
        }
    }

    std::vector<std::pair<std::string, std::string>> candidates;
    for (int x_length : reference.x_lengths) {
        for (int y_length : reference.y_lengths) {
            if (static_cast<int>(cb_tag.size()) != x_length + y_length) {
                continue;
            }

            std::string x = cb_tag.substr(0, static_cast<size_t>(x_length));
            std::string y = cb_tag.substr(static_cast<size_t>(x_length),
                                          static_cast<size_t>(y_length));
            if (reference.x_index_by_bc.find(x) !=
                    reference.x_index_by_bc.end() &&
                reference.y_index_by_bc.find(y) !=
                    reference.y_index_by_bc.end()) {
                candidates.emplace_back(x, y);
            }

            y = cb_tag.substr(0, static_cast<size_t>(y_length));
            x = cb_tag.substr(static_cast<size_t>(y_length),
                              static_cast<size_t>(x_length));
            if (reference.x_index_by_bc.find(x) !=
                    reference.x_index_by_bc.end() &&
                reference.y_index_by_bc.find(y) !=
                    reference.y_index_by_bc.end()) {
                candidates.emplace_back(x, y);
            }
        }
    }
    if (candidates.size() == 1) return candidates.front();
    return std::nullopt;
}

inline int project_axis_for_bin_size(int coordinate_2um, int bin_size_um) {
    const long long scaled = static_cast<long long>(coordinate_2um) * 2LL;
    return static_cast<int>(scaled / static_cast<long long>(bin_size_um));
}

inline std::string bin_read_group(int bin_size_um) {
    std::ostringstream output;
    output << 's' << std::setw(3) << std::setfill('0') << bin_size_um << "um";
    return output.str();
}

inline std::string spatial_mapping_id(int bin_size_um, int x, int y) {
    std::ostringstream output;
    output << "s_" << std::setw(3) << std::setfill('0') << bin_size_um
           << "um_" << zero_pad(x, 5) << '_' << zero_pad(y, 5) << "-1";
    return output.str();
}

inline std::string coordinate_read_id(const std::string& base_id, int x,
                                      int y) {
    return base_id + '_' + zero_pad(x, 5) + '_' + zero_pad(y, 5);
}

inline bool is_safe_barcode_filename(const std::string& barcode) {
    if (barcode.empty() || barcode == "." || barcode == "..") return false;
    return std::all_of(barcode.begin(), barcode.end(), [](unsigned char c) {
        const bool ascii_alphanumeric =
            (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
            (c >= '0' && c <= '9');
        return ascii_alphanumeric || c == '_' || c == '-' || c == '.' ||
               c == '+';
    });
}

inline fs::path absolute_normal(const fs::path& path) {
    return fs::absolute(path).lexically_normal();
}

inline bool existing_paths_equivalent(const fs::path& first,
                                      const fs::path& second) {
    std::error_code first_error;
    std::error_code second_error;
    const bool first_exists = fs::exists(first, first_error);
    const bool second_exists = fs::exists(second, second_error);
    if (first_error || second_error || !first_exists || !second_exists) {
        return false;
    }
    std::error_code equivalent_error;
    const bool equivalent = fs::equivalent(first, second, equivalent_error);
    if (equivalent_error) {
        throw std::runtime_error("could not verify reformat path identity");
    }
    return equivalent;
}

inline uint64_t next_stage_nonce() {
    static std::atomic<uint64_t> counter{1};
    const auto ticks = static_cast<uint64_t>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
    return ticks ^ counter.fetch_add(1, std::memory_order_relaxed);
}

class stage_directory {
public:
    explicit stage_directory(const fs::path& parent) {
        std::error_code error;
        if (!fs::is_directory(parent, error) || error) {
            throw std::runtime_error("reformat staging parent is not a directory: " +
                                     parent.string());
        }
        for (int attempt = 0; attempt < 1000; ++attempt) {
            path_ = parent /
                    (".rrad-reformat-stage-" +
                     std::to_string(next_stage_nonce()) + '-' +
                     std::to_string(attempt));
            error.clear();
            if (fs::create_directory(path_, error)) return;
            if (error && error != std::errc::file_exists) {
                throw std::runtime_error(
                    "could not create reformat staging directory: " +
                    error.message());
            }
        }
        throw std::runtime_error(
            "could not allocate a unique reformat staging directory");
    }

    stage_directory(const stage_directory&) = delete;
    stage_directory& operator=(const stage_directory&) = delete;

    ~stage_directory() noexcept {
        try {
            std::error_code error;
            fs::remove_all(path_, error);
        } catch (...) {
            // Cleanup is best-effort and must never turn a completed atomic
            // publication into a reported failure or terminate the process.
        }
    }

    const fs::path& path() const { return path_; }

private:
    fs::path path_;
};

class output_writer {
public:
    output_writer(const fs::path& path, bool append)
        : path_(path), gzip_(ends_with(path.string(), ".gz")) {
        if (gzip_) {
            gzip_file_ = gzopen(path.string().c_str(), append ? "ab" : "wb");
            if (!gzip_file_) {
                throw std::runtime_error("could not open gzip output: " +
                                         path.string());
            }
            if (gzbuffer(gzip_file_, 1 << 20) != 0) {
                (void)gzclose(gzip_file_);
                gzip_file_ = nullptr;
                throw std::runtime_error(
                    "could not allocate gzip output buffer: " + path.string());
            }
        } else {
            const auto mode = std::ios::out | std::ios::binary |
                              (append ? std::ios::app : std::ios::trunc);
            plain_file_.open(path, mode);
            if (!plain_file_) {
                throw std::runtime_error("could not open output: " +
                                         path.string());
            }
        }
    }

    output_writer(const output_writer&) = delete;
    output_writer& operator=(const output_writer&) = delete;

    ~output_writer() { close_noexcept(); }

    void write(const std::string& data) {
        if (closed_) {
            throw std::runtime_error("write attempted after output close: " +
                                     path_.string());
        }
        if (gzip_) {
            size_t offset = 0;
            while (offset < data.size()) {
                const size_t remaining = data.size() - offset;
                const unsigned int amount = static_cast<unsigned int>(
                    std::min<size_t>(remaining,
                                     static_cast<size_t>(INT_MAX)));
                const int written =
                    gzwrite(gzip_file_, data.data() + offset, amount);
                if (written <= 0) {
                    throw std::runtime_error("failed while writing gzip output: " +
                                             path_.string());
                }
                offset += static_cast<size_t>(written);
            }
        } else {
            plain_file_.write(data.data(),
                              static_cast<std::streamsize>(data.size()));
            if (!plain_file_) {
                throw std::runtime_error("failed while writing output: " +
                                         path_.string());
            }
        }
    }

    void close() {
        if (closed_) return;
        closed_ = true;
        if (gzip_file_) {
            gzFile file = gzip_file_;
            gzip_file_ = nullptr;
            const int status = gzclose(file);
            if (status != Z_OK) {
                throw std::runtime_error("failed while closing gzip output: " +
                                         path_.string() + " (zlib code " +
                                         std::to_string(status) + ')');
            }
        } else if (plain_file_.is_open()) {
            plain_file_.flush();
            const bool flush_failed = plain_file_.fail();
            plain_file_.close();
            if (flush_failed || plain_file_.fail()) {
                throw std::runtime_error("failed while closing output: " +
                                         path_.string());
            }
        }
    }

private:
    void close_noexcept() noexcept {
        if (closed_) return;
        closed_ = true;
        if (gzip_file_) {
            (void)gzclose(gzip_file_);
            gzip_file_ = nullptr;
        }
        if (plain_file_.is_open()) plain_file_.close();
    }

    fs::path path_;
    bool gzip_ = false;
    bool closed_ = false;
    gzFile gzip_file_ = nullptr;
    std::ofstream plain_file_;
};

class gzip_input_guard {
public:
    explicit gzip_input_guard(gzFile file) : file_(file) {}
    gzip_input_guard(const gzip_input_guard&) = delete;
    gzip_input_guard& operator=(const gzip_input_guard&) = delete;

    ~gzip_input_guard() {
        if (file_) (void)gzclose(file_);
    }

    gzFile get() const { return file_; }

    int close() {
        if (!file_) return Z_OK;
        gzFile file = file_;
        file_ = nullptr;
        return gzclose(file);
    }

private:
    gzFile file_ = nullptr;
};

inline std::string serialize_record(const read_streaming::sequence& record,
                                    bool to_fasta) {
    std::string output;
    if (record.is_fastq && !to_fasta) {
        output.reserve(record.id.size() + record.comment.size() +
                       record.seq.size() + record.plus.size() +
                       record.qual.size() + 16);
        output.push_back('@');
        output += record.id;
        if (!record.comment.empty()) {
            output.push_back('\t');
            output += record.comment;
        }
        output += '\n';
        output += record.seq;
        output += "\n+";
        output += record.plus;
        output.push_back('\n');
        output += record.qual;
        output.push_back('\n');
    } else {
        output.reserve(record.id.size() + record.comment.size() +
                       record.seq.size() + 8);
        output.push_back('>');
        output += record.id;
        if (!record.comment.empty()) {
            output.push_back('\t');
            output += record.comment;
        }
        output.push_back('\n');
        output += record.seq;
        output.push_back('\n');
    }
    return output;
}

inline void verify_gzip(const fs::path& path,
                        const std::function<void()>& boundary_poll) {
    gzip_input_guard input(gzopen(path.string().c_str(), "rb"));
    if (!input.get()) {
        throw std::runtime_error("could not reopen gzip output: " +
                                 path.string());
    }
    std::array<unsigned char, 64 * 1024> buffer{};
    size_t bytes_since_poll = 0;
    int count = 0;
    while ((count = gzread(input.get(), buffer.data(),
                           static_cast<unsigned int>(buffer.size()))) > 0) {
        bytes_since_poll += static_cast<size_t>(count);
        if (boundary_poll && bytes_since_poll >= (32U << 20)) {
            boundary_poll();
            bytes_since_poll = 0;
        }
    }
    if (count < 0) {
        int zlib_code = Z_OK;
        const char* detail = gzerror(input.get(), &zlib_code);
        const std::string message = detail ? detail : "";
        (void)input.close();
        throw std::runtime_error("gzip output failed validation: " +
                                 path.string() +
                                 (message.empty() ? "" : " (" + message + ')'));
    }
    const int status = input.close();
    if (status != Z_OK) {
        throw std::runtime_error("failed while closing verified gzip output: " +
                                 path.string() + " (zlib code " +
                                 std::to_string(status) + ')');
    }
}

inline void verify_output(const fs::path& path,
                          const std::function<void()>& boundary_poll) {
    std::error_code error;
    if (!fs::is_regular_file(path, error) || error) {
        throw std::runtime_error("reformat output is missing: " +
                                 path.string());
    }
    if (ends_with(path.string(), ".gz")) verify_gzip(path, boundary_poll);
}

inline void atomic_replace(const fs::path& staged, const fs::path& final_path) {
    // rrad currently targets Unix. C rename(2) atomically replaces an existing
    // non-directory destination without deleting it first; on failure the old
    // destination remains intact.
    if (std::rename(staged.string().c_str(), final_path.string().c_str()) != 0) {
        const int error = errno;
        throw std::runtime_error("could not atomically install reformat output " +
                                 final_path.string() + ": " +
                                 std::strerror(error));
    }
}

inline void atomic_install_no_clobber(const fs::path& staged,
                                      const fs::path& final_path) {
    // The stage is created beside the destination, so a hard-link operation
    // is both atomic and same-filesystem. Unlike rename(), create_hard_link()
    // cannot replace a destination created by a concurrent process.
    std::error_code error;
    fs::create_hard_link(staged, final_path, error);
    if (error) {
        if (error == std::errc::file_exists) {
            throw std::invalid_argument(
                "reformat output already exists: " + final_path.string());
        }
        throw std::runtime_error(
            "could not atomically install reformat output " +
            final_path.string() + ": " + error.message());
    }
    // The hard link above is the publication commit point. Cleanup after it
    // must be non-throwing: both names identify the same completed inode, and
    // the staging-directory guard can retry removal.
    try {
        error.clear();
        fs::remove(staged, error);
    } catch (...) {
        return;
    }
}

inline std::string input_read_error(const read_streaming& reader) {
    std::string message = "sequence read failed (kseq code " +
                          std::to_string(reader.read_error_code()) + ") in " +
                          reader.read_error_path();
    if (!reader.read_error_detail().empty()) {
        message += " (" + reader.read_error_detail() + ')';
    }
    return message;
}

struct split_writer_state {
    std::string barcode;
    fs::path staged_path;
    fs::path final_path;
    std::unique_ptr<output_writer> writer;
    bool has_data = false;
    uint64_t last_used = 0;
    uint64_t gzip_members = 0;
};

inline size_t split_commit_fail_after_install() noexcept {
#ifdef RRAD_EMBEDDED
    // Test-only allocation failpoint. Parsing is deliberately allocation-free
    // because it is consulted while preparing transactional bookkeeping.
    const char* raw = std::getenv(
        "RRAD_TEST_REFORMAT_FAIL_AFTER_SPLIT_INSTALL");
    if (!raw || !*raw) return 0;
    char* end = nullptr;
    errno = 0;
    const unsigned long long value = std::strtoull(raw, &end, 10);
    if (errno != 0 || end == raw || *end != '\0' || value == 0 ||
        value > static_cast<unsigned long long>(
                    std::numeric_limits<size_t>::max())) {
        return 0;
    }
    return static_cast<size_t>(value);
#else
    return 0;
#endif
}

inline void close_split_writers(
    std::unordered_map<std::string, std::unique_ptr<split_writer_state>>& states) {
    std::exception_ptr first_error;
    for (auto& item : states) {
        auto& state = *item.second;
        if (!state.writer) continue;
        try {
            state.writer->close();
        } catch (...) {
            if (!first_error) first_error = std::current_exception();
        }
        state.writer.reset();
    }
    if (first_error) std::rethrow_exception(first_error);
}

inline void commit_split_outputs(
    const std::vector<split_writer_state*>& ordered,
    const fs::path& stage_path, const fs::path& input_path, bool overwrite) {
    const fs::path backup_dir = stage_path / "previous";
    struct split_commit_state {
        split_writer_state* output;
        fs::path final_path;
        fs::path backup_path;
        bool had_existing_output;
        bool backup_installed = false;
        bool output_installed = false;
    };
    std::vector<split_commit_state> commit_states;
    commit_states.reserve(ordered.size());
    const size_t fail_after_install = split_commit_fail_after_install();

    // Preflight the full set and allocate every path and rollback slot before
    // moving any existing output. Once publication starts, bookkeeping is
    // reduced to setting booleans and cannot itself throw.
    for (size_t index = 0; index < ordered.size(); ++index) {
        auto* state = ordered[index];
        const fs::path target = absolute_normal(state->final_path);
        if (target == input_path || existing_paths_equivalent(target, input_path)) {
            throw std::invalid_argument(
                "split output would overwrite the reformat input: " +
                target.string());
        }
        std::error_code error;
        const bool exists = fs::exists(target, error);
        if (error) {
            throw std::runtime_error("could not inspect split output: " +
                                     target.string());
        }
        if (exists && fs::is_directory(target, error)) {
            throw std::invalid_argument(
                "split output identifies a directory: " + target.string());
        }
        if (error) {
            throw std::runtime_error("could not inspect split output: " +
                                     target.string());
        }
        if (exists && !overwrite) {
            throw std::invalid_argument("split output already exists: " +
                                        target.string());
        }
        commit_states.push_back(split_commit_state{
            state,
            target,
            backup_dir / (std::to_string(index) + ".fq.gz"),
            exists
        });
    }

    try {
        if (overwrite) {
            std::error_code error;
            if (!fs::create_directory(backup_dir, error) && error) {
                throw std::runtime_error(
                    "could not create split-output backup directory: " +
                    error.message());
            }
            for (auto& state : commit_states) {
                if (!state.had_existing_output) continue;
                fs::rename(state.final_path, state.backup_path);
                state.backup_installed = true;
            }
        }

        size_t installed_count = 0;
        for (auto& state : commit_states) {
            if (overwrite) {
                fs::rename(state.output->staged_path, state.final_path);
            } else {
                atomic_install_no_clobber(
                    state.output->staged_path, state.final_path);
            }
            state.output_installed = true;
            ++installed_count;
            if (fail_after_install == installed_count) {
                throw std::bad_alloc();
            }
        }
    } catch (...) {
        for (auto it = commit_states.rbegin();
             it != commit_states.rend(); ++it) {
            if (!it->output_installed) continue;
            std::error_code error;
            fs::remove(it->final_path, error);
        }
        for (auto it = commit_states.rbegin();
             it != commit_states.rend(); ++it) {
            if (!it->backup_installed) continue;
            std::error_code error;
            fs::rename(it->backup_path, it->final_path, error);
        }
        throw;
    }
}

}  // namespace reformat_detail

inline reformat_result reformat_fastx(
    const reformat_options& options,
    std::function<void()> boundary_poll = {},
    std::function<void(const reformat_result&)> before_publish = {}) {
    namespace detail = reformat_detail;
    namespace fs = std::filesystem;

    if (options.input_path.empty()) {
        throw std::invalid_argument("reformat input_path is required");
    }
    if (!options.split_by_barcode && !options.reformat_header &&
        !options.presto_header && !options.coordinate_mode.has_value() &&
        !options.to_fasta) {
        throw std::invalid_argument(
            "reformat requires split_by_barcode, header formatting, coordinate_mode, or to_fasta");
    }
    if (options.reformat_header && options.presto_header) {
        throw std::invalid_argument(
            "collapsed and pRESTO header formatting cannot be combined");
    }
    if (options.reformat_header && options.coordinate_mode.has_value()) {
        throw std::invalid_argument(
            "reformat_header cannot be combined with coordinate_mode");
    }
    if (options.presto_header && options.coordinate_mode.has_value()) {
        throw std::invalid_argument(
            "pRESTO header formatting cannot be combined with coordinate_mode");
    }
    if (options.presto_header && options.parse_collapsed_id) {
        throw std::invalid_argument(
            "pRESTO header formatting requires SAM tags and cannot parse collapsed IDs");
    }
    if (options.split_by_barcode && options.split_output_dir.empty()) {
        throw std::invalid_argument(
            "split_output_dir is required when split_by_barcode is enabled");
    }
    if (options.split_by_barcode && !options.output_path.empty()) {
        throw std::invalid_argument(
            "output_path cannot be combined with split_by_barcode");
    }
    if (options.to_fasta && !options.split_by_barcode &&
        options.output_path.empty()) {
        throw std::invalid_argument(
            "output_path is required for aggregate FASTA conversion");
    }
    if (options.to_fasta && !options.split_by_barcode &&
        !detail::has_fasta_suffix(options.output_path)) {
        throw std::invalid_argument(
            "FASTA output_path must end in .fa, .fasta, .fa.gz, or .fasta.gz");
    }
    if (options.delimiter == '\0' ||
        std::isspace(static_cast<unsigned char>(options.delimiter))) {
        throw std::invalid_argument(
            "reformat delimiter cannot be NUL or whitespace");
    }
    if (options.chunk_size == 0) {
        throw std::invalid_argument("reformat chunk_size must be positive");
    }
    if (options.max_open_split_writers == 0) {
        throw std::invalid_argument(
            "max_open_split_writers must be positive");
    }
    if (options.threads < 1) {
        throw std::invalid_argument("reformat threads must be positive");
    }

    std::error_code error;
    fs::path input_path = detail::absolute_normal(options.input_path);
    if (!fs::is_regular_file(input_path, error) || error) {
        throw std::invalid_argument("reformat input is not a regular file: " +
                                    input_path.string());
    }
    input_path = fs::canonical(input_path);
    const std::string split_file_suffix = options.to_fasta
        ? ".fa.gz"
        : path_utils::get_fastqa_type(input_path.string()) + ".gz";

    std::optional<detail::vizhd_axis_reference> spatial_reference;
    if (options.coordinate_mode) {
        if (options.coordinate_mode->empty()) {
            throw std::invalid_argument("coordinate_mode cannot be empty");
        }
        const std::string mode = detail::lower_copy(*options.coordinate_mode);
        if (!(mode == "vizhd-v1" || mode == "visium-hd-v1" ||
              mode == "visium_hd_v1")) {
            throw std::invalid_argument(
                "unsupported coordinate_mode: " + *options.coordinate_mode +
                " (supported: vizHD-v1)");
        }
        if (options.bin_size_um <= 0) {
            throw std::invalid_argument("bin_size_um must be positive");
        }
        spatial_reference = detail::load_vizhd_v1_axis_reference();
        if (options.verbose) {
            RAD_COUT << "[reformat] coordinate mode: "
                      << *options.coordinate_mode
                      << ", bin-size: " << options.bin_size_um << "um"
                      << ", x barcodes: "
                      << spatial_reference->x_index_by_bc.size()
                      << ", y barcodes: "
                      << spatial_reference->y_index_by_bc.size() << '\n';
        }
    }

    reformat_result result;
    result.input_path = input_path.string();

    fs::path final_output;
    std::unique_ptr<detail::stage_directory> stage;
    std::unique_ptr<detail::output_writer> single_writer;
    std::unordered_map<std::string,
                       std::unique_ptr<detail::split_writer_state>>
        split_states;
    size_t open_split_writers = 0;
    uint64_t use_counter = 0;

    if (options.split_by_barcode) {
        fs::path split_dir = detail::absolute_normal(options.split_output_dir);
        if (fs::exists(split_dir, error)) {
            if (error || !fs::is_directory(split_dir, error) || error) {
                throw std::invalid_argument(
                    "split_output_dir is not a directory: " +
                    split_dir.string());
            }
        } else if (!fs::create_directories(split_dir, error) || error) {
            throw std::runtime_error(
                "could not create split_output_dir: " + split_dir.string());
        }
        split_dir = fs::canonical(split_dir);
        result.split_output_dir = split_dir.string();
        stage = std::make_unique<detail::stage_directory>(split_dir);
    } else {
        final_output = options.output_path.empty()
                           ? input_path
                           : detail::absolute_normal(options.output_path);
        fs::path parent = final_output.parent_path();
        if (parent.empty()) parent = fs::current_path();
        if (!fs::exists(parent, error)) {
            if (!fs::create_directories(parent, error) || error) {
                throw std::runtime_error(
                    "could not create reformat output directory: " +
                    parent.string());
            }
        }
        parent = fs::canonical(parent);
        final_output = parent / final_output.filename();

        const bool same_path = final_output == input_path;
        if (!same_path && detail::existing_paths_equivalent(final_output,
                                                             input_path)) {
            throw std::invalid_argument(
                "output_path is a second name for the reformat input");
        }
        if (fs::exists(final_output, error)) {
            if (error || fs::is_directory(final_output, error) || error) {
                throw std::invalid_argument(
                    "reformat output identifies a directory: " +
                    final_output.string());
            }
            if (!same_path && !options.overwrite) {
                throw std::invalid_argument("reformat output already exists: " +
                                            final_output.string());
            }
        } else if (error) {
            throw std::runtime_error("could not inspect reformat output: " +
                                     final_output.string());
        }

        stage = std::make_unique<detail::stage_directory>(parent);
        const fs::path staged_output =
            stage->path() / final_output.filename();
        single_writer =
            std::make_unique<detail::output_writer>(staged_output, false);
        result.output_path = final_output.string();
    }

    auto acquire_split_writer = [&](const std::string& barcode)
        -> detail::split_writer_state& {
        if (!detail::is_safe_barcode_filename(barcode)) {
            throw std::invalid_argument(
                "CB tag is unsafe for a split-output filename: " + barcode);
        }

        auto found = split_states.find(barcode);
        if (found == split_states.end()) {
            auto state = std::make_unique<detail::split_writer_state>();
            state->barcode = barcode;
            state->staged_path =
                stage->path() / (barcode + split_file_suffix);
            state->final_path =
                fs::path(result.split_output_dir) /
                (barcode + split_file_suffix);
            std::error_code staged_error;
            const bool staged_exists =
                fs::exists(state->staged_path, staged_error);
            if (staged_error) {
                throw std::runtime_error(
                    "could not inspect split-output staging path: " +
                    state->staged_path.string());
            }
            if (staged_exists) {
                throw std::invalid_argument(
                    "distinct CB tags resolve to the same split-output filename "
                    "under this filesystem's case rules: " + barcode);
            }
            if (detail::absolute_normal(state->final_path) == input_path ||
                detail::existing_paths_equivalent(state->final_path,
                                                  input_path)) {
                throw std::invalid_argument(
                    "split output would overwrite the reformat input: " +
                    state->final_path.string());
            }
            found = split_states.emplace(barcode, std::move(state)).first;
        }

        auto& state = *found->second;
        if (!state.writer) {
            if (open_split_writers >= options.max_open_split_writers) {
                detail::split_writer_state* oldest = nullptr;
                for (auto& item : split_states) {
                    auto* candidate = item.second.get();
                    if (candidate->writer &&
                        (!oldest || candidate->last_used < oldest->last_used)) {
                        oldest = candidate;
                    }
                }
                if (!oldest) {
                    throw std::logic_error(
                        "split writer accounting is inconsistent");
                }
                oldest->writer->close();
                oldest->writer.reset();
                --open_split_writers;
            }
            state.writer = std::make_unique<detail::output_writer>(
                state.staged_path, state.has_data);
            ++open_split_writers;
            ++state.gzip_members;
        }
        state.last_used = ++use_counter;
        return state;
    };

    // Disable pigz explicitly: the embedded API performs no subprocesses and
    // its gzip validation/error behavior must not depend on the host PATH.
    file_streaming files(input_path.string(), 1, false);
    read_streaming reader(files, true);
    uint64_t records_since_poll = 0;

    struct prepared_record {
        std::optional<std::string> barcode;
        std::string serialized;
        std::string unmapped_id;
        bool skip_missing_barcode = false;
        bool reformatted = false;
        bool converted_to_fasta = false;
        bool coordinate_mapped = false;
        bool coordinate_unmapped = false;
    };

    auto prepare_public_result = [&]() {
        // Logging and host result construction may allocate. They must happen
        // before either publication path makes its staged output visible.
        if (options.verbose) {
            RAD_COUT << "[reformat] wrote " << result.stats.records_written
                     << " reads\n";
            if (options.split_by_barcode) {
                RAD_COUT << "[reformat] barcodes: "
                         << result.stats.split_barcodes << '\n';
            }
        }
        if (before_publish) before_publish(result);
    };

    try {
        bool reached_end = false;
        while (!reached_end) {
            if (boundary_poll) boundary_poll();

            std::vector<read_streaming::sequence> chunk;
            chunk.reserve(options.chunk_size);
            while (chunk.size() < options.chunk_size) {
                auto next = reader.next_sequence();
                if (!next) {
                    reached_end = true;
                    break;
                }
                chunk.push_back(std::move(*next));
            }
            if (chunk.empty()) break;

            result.stats.records_read +=
                static_cast<uint64_t>(chunk.size());
            records_since_poll += static_cast<uint64_t>(chunk.size());

            std::vector<prepared_record> prepared(chunk.size());
            std::vector<std::exception_ptr> errors(chunk.size());
#pragma omp parallel for num_threads(options.threads) schedule(static)
            for (std::ptrdiff_t raw_index = 0;
                 raw_index < static_cast<std::ptrdiff_t>(chunk.size());
                 ++raw_index) {
                const std::size_t index =
                    static_cast<std::size_t>(raw_index);
                try {
                    auto record = std::move(chunk[index]);
                    auto cb = options.presto_header
                                  ? detail::extract_tag_value(
                                        record.comment, "CB:Z:")
                                  : detail::extract_cb(
                                        record.id, record.comment);
                    auto ub = detail::extract_ub(record.comment);
                    std::optional<std::string> presto_barcode;
                    if (options.presto_header) {
                        presto_barcode = detail::extract_bc(record.comment);
                    }
                    if (options.presto_header &&
                        (!presto_barcode || presto_barcode->empty())) {
                        presto_barcode = cb;
                    }
                    std::string qname = record.id;
                    // Explicit SAM-style tags are authoritative. Only infer
                    // CB/UB from a previously collapsed ID when neither tag is
                    // present, and do not mutate the record merely to inspect
                    // it. Split-only and failed coordinate conversions must
                    // preserve the original read name byte-for-byte.
                    if (options.parse_collapsed_id && !cb && !ub) {
                        if (auto parsed = detail::parse_collapsed_id(
                                record.id, options.delimiter)) {
                            auto [parsed_qname, parsed_cb, parsed_ub] = *parsed;
                            qname = std::move(parsed_qname);
                            cb = std::move(parsed_cb);
                            ub = std::move(parsed_ub);
                        }
                    }

                    if (spatial_reference) {
                        bool mapped = false;
                        if (cb && ub && !cb->empty() && !ub->empty()) {
                            const auto xy = detail::resolve_xy_barcodes(
                                *cb, *spatial_reference);
                            if (xy) {
                                const auto x =
                                    spatial_reference->x_index_by_bc.find(
                                        xy->first);
                                const auto y =
                                    spatial_reference->y_index_by_bc.find(
                                        xy->second);
                                if (x !=
                                        spatial_reference->x_index_by_bc.end() &&
                                    y != spatial_reference->y_index_by_bc.end()) {
                                    const int x_bin =
                                        detail::project_axis_for_bin_size(
                                            x->second, options.bin_size_um);
                                    const int y_bin =
                                        detail::project_axis_for_bin_size(
                                            y->second, options.bin_size_um);
                                    record.id = detail::coordinate_read_id(
                                        qname, x_bin, y_bin);
                                    record.comment = detail::visium_hd_comment(
                                        detail::spatial_mapping_id(
                                            options.bin_size_um, x_bin, y_bin),
                                        *cb, *ub,
                                        detail::bin_read_group(
                                            options.bin_size_um),
                                        record.comment);
                                    mapped = true;
                                    prepared[index].coordinate_mapped = true;
                                    prepared[index].reformatted = true;
                                }
                            }
                        }
                        if (!mapped) {
                            prepared[index].coordinate_unmapped = true;
                            prepared[index].unmapped_id = record.id;
                        }
                    } else if (options.presto_header) {
                        record.id = detail::presto_id(
                            qname, presto_barcode, ub);
                        record.comment.clear();
                        prepared[index].reformatted = true;
                    } else if (options.reformat_header) {
                        record.id = detail::collapse_id(
                            qname, cb, ub, options.delimiter);
                        record.comment.clear();
                        prepared[index].reformatted = true;
                    }

                    prepared[index].barcode = std::move(cb);
                    prepared[index].converted_to_fasta =
                        options.to_fasta && record.is_fastq;
                    if (options.split_by_barcode &&
                        (!prepared[index].barcode ||
                         prepared[index].barcode->empty())) {
                        prepared[index].skip_missing_barcode = true;
                    } else {
                        prepared[index].serialized =
                            detail::serialize_record(record, options.to_fasta);
                    }
                } catch (...) {
                    errors[index] = std::current_exception();
                }
            }

            for (const auto& error_ptr : errors) {
                if (error_ptr) std::rethrow_exception(error_ptr);
            }

            for (auto& record : prepared) {
                if (record.reformatted) {
                    ++result.stats.records_reformatted;
                }
                if (record.coordinate_mapped) {
                    ++result.stats.coordinate_mapped;
                }
                if (record.coordinate_unmapped) {
                    ++result.stats.coordinate_unmapped;
                    if (options.verbose &&
                        result.stats.coordinate_unmapped <= 10) {
                        RAD_COUT
                            << "[reformat] skipped coordinate mapping for read "
                            << record.unmapped_id << '\n';
                    } else if (options.verbose &&
                               result.stats.coordinate_unmapped == 11) {
                        RAD_COUT
                            << "[reformat] additional unmapped-read diagnostics "
                               "suppressed\n";
                    }
                }
                if (record.skip_missing_barcode) {
                    ++result.stats.records_skipped_missing_cb;
                    continue;
                }
                if (record.converted_to_fasta) {
                    ++result.stats.records_converted_to_fasta;
                }

                if (options.split_by_barcode) {
                    auto& state = acquire_split_writer(*record.barcode);
                    state.writer->write(record.serialized);
                    state.has_data = true;
                } else {
                    single_writer->write(record.serialized);
                }
                ++result.stats.records_written;

                if (options.verbose &&
                    result.stats.records_written % 500000 == 0) {
                    RAD_COUT << "[reformat] processed "
                              << result.stats.records_written << " reads\n";
                }
            }

            if (boundary_poll &&
                records_since_poll >= options.chunk_size) {
                boundary_poll();
                records_since_poll = 0;
            }
        }

        if (reader.read_failed()) {
            throw std::runtime_error(detail::input_read_error(reader));
        }
        result.stats.empty_sequences_skipped = reader.empty_sequences_skipped();
        if (boundary_poll) boundary_poll();

        if (options.split_by_barcode) {
            detail::close_split_writers(split_states);
            open_split_writers = 0;

            std::vector<detail::split_writer_state*> ordered;
            ordered.reserve(split_states.size());
            for (auto& item : split_states) ordered.push_back(item.second.get());
            std::sort(ordered.begin(), ordered.end(), [](const auto* first,
                                                         const auto* second) {
                return first->barcode < second->barcode;
            });
            for (const auto* state : ordered) {
                if (boundary_poll) boundary_poll();
                detail::verify_output(state->staged_path, boundary_poll);
                result.stats.gzip_members_written += state->gzip_members;
            }

            // Result construction can allocate. Complete it before the first
            // split file is made visible so a late allocation failure cannot
            // turn a successful publication into a reported failed call.
            result.split_paths.reserve(ordered.size());
            for (const auto* state : ordered) {
                result.split_paths.emplace_back(
                    state->barcode,
                    detail::absolute_normal(state->final_path).string());
            }
            result.stats.split_barcodes =
                static_cast<uint64_t>(ordered.size());
            prepare_public_result();
            // No potentially throwing work follows a successful commit: the
            // result and host return object are already complete.
            detail::commit_split_outputs(
                ordered, stage->path(), input_path, options.overwrite);
        } else {
            single_writer->close();
            single_writer.reset();
            const fs::path staged_output =
                stage->path() / final_output.filename();
            detail::verify_output(staged_output, boundary_poll);
            if (detail::ends_with(staged_output.string(), ".gz")) {
                result.stats.gzip_members_written = 1;
            }

            if (final_output == input_path) {
                const auto permissions = fs::status(input_path).permissions();
                fs::permissions(staged_output, permissions,
                                fs::perm_options::replace);
            }
            prepare_public_result();
            // Atomic publication is the last potentially throwing operation.
            if (final_output == input_path || options.overwrite) {
                detail::atomic_replace(staged_output, final_output);
            } else {
                detail::atomic_install_no_clobber(
                    staged_output, final_output);
            }
        }
    } catch (...) {
        if (options.split_by_barcode) {
            try {
                detail::close_split_writers(split_states);
            } catch (...) {
                // Preserve the primary read/write/interrupt failure.
            }
        } else if (single_writer) {
            try {
                single_writer->close();
            } catch (...) {
                // Preserve the primary read/write/interrupt failure.
            }
        }
        throw;
    }
    return result;
}
