#pragma once
#include "rad_headers.h"

namespace basic_stats {
    double mean(const std::vector<double>& x) {
        if (x.empty()) return 0.0;
        return std::accumulate(x.begin(), x.end(), 0.0) / x.size();
    }

    double sd(const std::vector<double>& x) {
        if (x.size() <= 1) return 0.0;
        double m = mean(x);
        double sum = 0.0;
        for (const auto& val : x) {
            double diff = val - m;
            sum += diff * diff;
        }
        return std::sqrt(sum / (x.size() - 1));
    }

    double kmer_match_prob(int str1_len, int str2_len, int k) {
        if (k > str1_len || k > str2_len || k <= 0) {
            return 0.0;
        }
        // Number of k-mers in each string
        int n1 = str1_len - k + 1;  // k-mers in string1
        int n2 = str2_len - k + 1;  // k-mers in string2
        // Total possible k-mers (4^k for DNA)
        double total_possible_kmers = std::pow(4.0, k);
        // Expected number of matches
        double expected_matches = (static_cast<double>(n1) * n2) / total_possible_kmers;
        // Probability of at least one match (1 - P(no matches))
        return 1.0 - std::exp(-expected_matches);
    }

    int min_match_len(int str1_len, int str2_len, double max_prob) {
        if (str1_len <= 0 || str2_len <= 0) {
            return -1; // Invalid input
        }
        // Start with k=1 and find first k where probability drops below threshold
        for (int k = 1; k <= std::min(str1_len, str2_len); ++k) {
            double prob = kmer_match_prob(str1_len, str2_len, k);
            if (prob < max_prob) {
                return k;
            }
        }
        // If we never reach the threshold, return the maximum possible k
        return std::min(str1_len, str2_len);
    }

};

namespace memory_utils {

    inline double to_mib(std::size_t bytes) {
        return double(bytes) / 1024.0 / 1024.0;
    }

    inline std::optional<double> current_rss_gib() {
        #if defined(__APPLE__)
            mach_task_basic_info info;
            mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
            if (task_info(mach_task_self(),
                          MACH_TASK_BASIC_INFO,
                          reinterpret_cast<task_info_t>(&info),
                          &count) == KERN_SUCCESS) {
                return double(info.resident_size) /
                       (1024.0 * 1024.0 * 1024.0);
            }
        #elif defined(__linux__)
            const long page_size = sysconf(_SC_PAGESIZE);
            if (page_size <= 0) {
                return std::nullopt;
            }
            std::ifstream statm("/proc/self/statm");
            long total_pages = 0;
            long resident_pages = 0;
            if (statm && (statm >> total_pages >> resident_pages)) {
                return (double(resident_pages) * double(page_size)) /
                       (1024.0 * 1024.0 * 1024.0);
            }
        #endif
        return std::nullopt;
    }

    inline std::optional<double> peak_rss_gib() {
        struct rusage usage {};
        if (getrusage(RUSAGE_SELF, &usage) != 0) {
            return std::nullopt;
        }
        #if defined(__APPLE__)
            // macOS reports ru_maxrss in bytes.
            return double(usage.ru_maxrss) /
                   (1024.0 * 1024.0 * 1024.0);
        #elif defined(__linux__)
            // Linux reports ru_maxrss in KiB.
            return double(usage.ru_maxrss) / (1024.0 * 1024.0);
        #else
            return std::nullopt;
        #endif
    }

    // pointer mem size
    inline constexpr std::size_t get_pointer_mem() {
        return sizeof(void*);
    }

    inline void get_rss() {
        const auto rss = current_rss_gib();
        if (rss) {
            RAD_COUT << "[memory_utils] RSS: " << *rss << " GiB\n";
        } else {
            RAD_CERR << "[memory_utils] RSS: unavailable\n";
        }
    }
};

namespace path_utils {
    namespace bfs = std::filesystem;
  
    inline bfs::path find_resource_root() {
      // 0) Explicit override.
      if (const char* env_root = std::getenv("RAD_RESOURCES"); env_root && *env_root) {
        std::error_code ec;
        bfs::path p(env_root);
        if (p.is_relative()) p = bfs::absolute(p, ec);
        if (!ec && bfs::exists(p) && bfs::is_directory(p)) {
          // Accept either <root>/resources or direct resources path.
          bfs::path cand = (p.filename() == "resources") ? p : (p / "resources");
          if (bfs::exists(cand) && bfs::is_directory(cand)) {
            return bfs::canonical(cand);
          }
        }
      }

      // 1) Compile-time source directory (developer builds only).  Installed
      // R packages pass their private resource directory through
      // RAD_RESOURCES, so executable-location probing is deliberately omitted.
#ifdef RAD_SOURCE_DIR
      {
        bfs::path src_res = bfs::path(RAD_SOURCE_DIR) / "resources";
        if (bfs::exists(src_res) && bfs::is_directory(src_res)) {
          return bfs::canonical(src_res);
        }
      }
#endif

      // 2) Fallback: ./resources in your CWD
      bfs::path alt = bfs::current_path() / "resources";
      if (bfs::exists(alt) && bfs::is_directory(alt)) {
        return bfs::canonical(alt);
      }

      throw std::runtime_error("Cannot locate 'resources' folder up the tree or in CWD");
    }

    // returns either ".fa" or ".fq" based on the file extension it's passed originally
    inline std::string get_fastqa_type(const std::string& file_path) {
        bfs::path p(file_path);
        // Get the extension
        std::string ext = p.extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(),
                       [](unsigned char value) {
                           return static_cast<char>(std::tolower(value));
                       });
        // Handle .gz files - check the extension of the stem
        if (ext == ".gz") {
            ext = p.stem().extension().string();
            std::transform(ext.begin(), ext.end(), ext.begin(),
                           [](unsigned char value) {
                               return static_cast<char>(std::tolower(value));
                           });
        }
        if (ext == ".fa" || ext == ".fasta") {
            return ".fa";
        }
        if (ext == ".fq" || ext == ".fastq") {
            return ".fq";
        }
        return ".fq"; //default to fq
    }
};

namespace seq_utils {

    std::string revcomp(const std::string& seq) {
        if (seq.empty()) return seq;
        std::string rc = seq;
        std::reverse(rc.begin(), rc.end());
        std::transform(rc.begin(), rc.end(), rc.begin(), [](char c) {
            switch (c) {
                case 'A': return 'T';
                case 'T': return 'A';
                case 'G': return 'C';
                case 'C': return 'G';
                case 'N': return 'N';
                case 'n': return 'n';
                default: return c;
            }
        });
        return rc;
    }
    
    std::string trim(const std::string& s) {
        size_t start = s.find_first_not_of(" \t\n\r\f\v");
        if (start == std::string::npos)
            return "";
        size_t end = s.find_last_not_of(" \t\n\r\f\v");
        return s.substr(start, end - start + 1);
    }

    bool is_rc(const std::string& str) {
        return str.size() >= 3 && str.substr(0, 3) == "rc_";
    }

    std::string remove_rc(const std::string& str) {
        if (str.size() >= 3 && str.substr(0, 3) == "rc_") {
            return str.substr(3); // Return string without the "rc_" prefix
        }
        return str; // Return original string if no prefix found
    }

    /*
    std::string substr_w_padding(const std::string& str, int start, int stop, int pad) {
        int start_pos = std::max(1, start - pad) -1 ;
        int stop_pos = std::min(static_cast<int>(str.length()), stop + pad);
        int length = stop_pos - start_pos;
        if(length <= 0) {
            return ""; // Return empty string if length is non-positive
        }
        return str.substr(start_pos, length);
    }*/
   std::string substr_w_padding(const std::string& str, int start, int stop, int pad) {
    if (str.empty()) return "";

    // Safety: clamp start/stop into [1, str.length()]
    int n = static_cast<int>(str.length());
    start = std::max(1, std::min(start, n));
    stop  = std::max(1, std::min(stop,  n));

    // Ensure stop >= start
    if (stop < start) std::swap(start, stop);

    int start_pos = std::max(1, start - pad) - 1; // convert to 0-based
    int stop_pos  = std::min(n, stop + pad);
    int length    = stop_pos - start_pos;

    if (length <= 0) return "";
    return str.substr(start_pos, length);
}

    static std::vector<std::pair<size_t,std::string>>extract(const std::string& seq, const std::string& pattern, size_t length, int offset){
        std::vector<std::pair<size_t,std::string>> results;
        if (pattern.empty() || length == 0 || seq.empty())
            return results;
        size_t patLen = pattern.size();
        // search for every match of pattern
        for (size_t pos = seq.find(pattern);
            pos != std::string::npos;
            pos = seq.find(pattern, pos + 1))
        {
            // compute start index (can be before or after the match)
            ptrdiff_t start = offset >= 0 ? static_cast<ptrdiff_t>(pos + patLen) + offset  : static_cast<ptrdiff_t>(pos) + offset;

            // ensure we stay in bounds
            if (start >= 0
            && static_cast<size_t>(start) + length <= seq.size())
            {
                results.emplace_back(
                    pos,
                    seq.substr(static_cast<size_t>(start), length)
                );
            }
        }
        return results;
    }

    int int_kmerize(const std::string& sequence, int k) {
        if (k <= 0 || k > static_cast<int>(sequence.length())) {
            return 0;
        }
        std::unordered_set<std::string> unique_kmers;
        for (int i = 0; i <= static_cast<int>(sequence.length()) - k; ++i) {
            std::string kmer = sequence.substr(i, k);
            unique_kmers.insert(kmer);
        }
        return static_cast<int>(unique_kmers.size());
    }

    std::vector<std::string> kmerize(const std::string& sequence, int k) {
        std::vector<std::string> kmers;
        if (k <= 0 || k > static_cast<int>(sequence.length())) {
            return kmers;  // Return empty vector if invalid k
        }
        kmers.reserve(sequence.length() - k + 1);
        // Generate all k-mers
        for (size_t i = 0; i <= sequence.length() - k; ++i) {
            kmers.push_back(sequence.substr(i, k));
        }
        return kmers;
    }
    
    std::vector<std::string> circ_kmerize(const std::string& seq, size_t k = 16) {
        if (k == 0) 
            throw std::invalid_argument("k must be > 0");
        size_t n = seq.size();
        if (n == 0) 
            return {};  // nothing to do on an empty sequence
        // Case 1: normal k‑mers
        if (n >= k) {
            std::vector<std::string> kmers;
            kmers.reserve(n - k + 1);
            for (size_t i = 0; i + k <= n; ++i) {
                kmers.emplace_back(seq.substr(i, k));
            }
            return kmers;
        }
    
        // Case 2: k > n -> circular wrap
        std::vector<std::string> kmers;
        kmers.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            std::string kmer;
            kmer.reserve(k);
            for (size_t j = 0; j < k; ++j) {
                // wrap index via modulo
                kmer.push_back(seq[(i + j) % n]);
            }
            kmers.emplace_back(std::move(kmer));
        }
        return kmers;
    }

};

namespace term_print_utils {
    inline const std::vector<std::string>& palette() {
        static const std::vector<std::string> colors = {
            "\033[38;5;31m",  // teal
            "\033[38;5;27m",  // blue
            "\033[38;5;35m",  // purple
            "\033[38;5;166m", // orange
            "\033[38;5;100m", // slate
            "\033[38;5;64m",  // green
            "\033[38;5;160m", // red
            "\033[38;5;214m", // amber
            "\033[38;5;45m",  // cyan
            "\033[38;5;94m"   // brownish
        };
        return colors;
    }

    inline std::string colorize(const std::string& text, const std::string& color_code) {
        static const char* reset = "\033[0m";
        return color_code + text + reset;
    }

    inline const std::string& direction_label_color() {
        static const std::string color_dir = "\033[38;5;247m"; // muted gray
        return color_dir;
    }

    inline const std::string& class_color(const std::string& class_id, const std::string& seq = "") {
        static std::unordered_map<std::string, std::string> cache;
        std::string key = seq.empty() ? class_id : seq;
        if (!seq.empty()) {
            std::string rc = seq_utils::revcomp(seq);
            if (rc < key) key = rc; // canonicalize to group fwd/rc together
        }
        auto it = cache.find(key);
        if (it != cache.end()) return it->second;
        size_t idx = std::hash<std::string>{}(key) % palette().size();
        cache[key] = palette()[idx];
        return cache[key];
    }
}

namespace streaming_utils { 
    std::vector<std::string> import_text(const std::string &path, size_t max_lines = SIZE_MAX) {
        // 1) open the right kind of stream
        std::unique_ptr<std::istream> in;
        igzstream* gzip_in = nullptr;
        std::string lower_path = path;
        std::transform(lower_path.begin(), lower_path.end(),
                       lower_path.begin(), [](unsigned char value) {
                           return static_cast<char>(std::tolower(value));
                       });
        if (lower_path.size() > 3 &&
            lower_path.substr(lower_path.size() - 3) == ".gz") {
            auto stream = std::make_unique<igzstream>(path.c_str());
            gzip_in = stream.get();
            in = std::move(stream);
        } else {
            in = std::make_unique<std::ifstream>(path);
        }
        if (!in || !*in) {
            throw std::runtime_error("Failed to open file: " + path);
        }
        // 2) read lines
        std::vector<std::string> out;
        out.reserve(std::min(max_lines, size_t(10)));
        std::string line;
        size_t count = 0;
        while (count < max_lines && std::getline(*in, line)) {
            out.push_back(std::move(line));
            ++count;
        }
        if (gzip_in) {
            const int zlib_code = gzip_in->zlib_error_code();
            const std::string zlib_message = gzip_in->zlib_error_message();
            gzip_in->close();
            if (zlib_code != Z_OK && zlib_code != Z_STREAM_END) {
                throw std::runtime_error(
                    "gzip text read failed: " + path + " (zlib code " +
                    std::to_string(zlib_code) +
                    (zlib_message.empty() ? ")"
                                          : ", " + zlib_message + ")"));
            }
            if (gzip_in->bad()) {
                throw std::runtime_error(
                    "gzip text close failed: " + path);
            }
        } else if (in->bad()) {
            throw std::runtime_error("text read failed: " + path);
        }
        return out;
    }
};

namespace whitelist_utils {
    inline const std::unordered_map<std::string, std::string>& default_kit_wl_paths() {
        static const std::unordered_map<std::string, std::string> defaults = {
        //splitseq barcodes
        {"splitseq_bc1", "resources/wl/splitseq_bc1_bitlist.csv.gz"},
        {"splitseq_bc2", "resources/wl/splitseq_bc2_bitlist.csv.gz"},
        //cellranger kits and whitelists--this is just the total list that I could pull together

        //10x 3' kits
        { "10x_3v1", "resources/wl/737K-august-2016_bitlist.csv.gz" },
        { "10x_3v2", "resources/wl/737K-august-2016_bitlist.csv.gz" },
        { "10x_3v3", "resources/wl/3M-february-2018-3v3.txt_bitlist.csv.gz" },
        { "10x_3v3.1", "resources/wl/3M-february-2018-3v3.txt_bitlist.csv.gz" },
        { "10x_3HTv3.1","resources/wl/3M-february-2018-3v3.txt_bitlist.csv.gz" },
        { "10x_3v4", "resources/wl/3M-3pgex-may-2023.txt_bitlist.csv.gz" },

        //10x 5' kits
        { "10x_5v1", "resources/wl/737K-august-2016_bitlist.csv.gz" },
        { "10x_5v2", "resources/wl/737K-august-2016_bitlist.csv.gz" },
        { "10x_5HTv2", "resources/wl/737K-august-2016_bitlist.csv.gz" },
        { "10x_5v3", "resources/wl/3M-5pgex-jan-2023.txt_bitlist.csv.gz" },

        //10x Visium kits (and Xenium? need to check)
        { "10x_Vis_V1", "resources/wl/visium-v1_v2_bitlist.csv.gz" },
        { "10x_Vis_V2", "resources/wl/visium-v1_v2_bitlist.csv.gz" },
        { "10x_Vis_V3", "resources/wl/visium-v3_v4_bitlist.csv.gz" },
        { "10x_Vis_V4", "resources/wl/visium-v3_v4_bitlist.csv.gz" },
        { "10x_Vis_V5", "resources/wl/visium-v5_bitlist.csv.gz" },

        //10x Visium HD
        { "visium_hd_bc1", "resources/wl/visium_hd_bc1.csv.gz" },
        { "visium_hd_bc2", "resources/wl/visium_hd_bc2.csv.gz" }
        };
        return defaults;
    }

    // Effective map (defaults + runtime overrides - tombstoned removals)
    static std::unordered_map<std::string, std::string> kit_wl_paths = default_kit_wl_paths();

    inline std::filesystem::path whitelist_override_path() {
        if (const char* home = std::getenv("HOME"); home && *home) {
            return std::filesystem::path(home) / ".rad" / "whitelist_overrides.tsv";
        }
        return std::filesystem::current_path() / ".rad_whitelist_overrides.tsv";
    }

    inline std::unordered_set<std::string>& removed_whitelist_keys() {
        static std::unordered_set<std::string> removed;
        return removed;
    }

    inline bool& whitelist_overrides_loaded() {
        static bool loaded = false;
        return loaded;
    }

    inline void load_whitelist_overrides_once() {
        // rrad's embedded backend must be reproducible from its bundled
        // resources.  Do not let per-user RAD CLI state change that backend.
        if (const char* embedded_only = std::getenv("RAD_EMBEDDED_ONLY");
            embedded_only && *embedded_only) {
            return;
        }
        if (whitelist_overrides_loaded()) return;
        whitelist_overrides_loaded() = true;

        std::ifstream in(whitelist_override_path().string());
        if (!in.is_open()) return;

        std::string line;
        while (std::getline(in, line)) {
            std::string trimmed_line = seq_utils::trim(line);
            if (trimmed_line.empty() || trimmed_line[0] == '#') continue;

            auto tab = line.find('\t');
            if (tab == std::string::npos) continue;

            std::string key = seq_utils::trim(line.substr(0, tab));
            std::string value = seq_utils::trim(line.substr(tab + 1));
            if (key.empty()) continue;

            if (value.empty()) {
                kit_wl_paths.erase(key);
                removed_whitelist_keys().insert(key);
            } else {
                kit_wl_paths[key] = value;
                removed_whitelist_keys().erase(key);
            }
        }
    }

    inline bool persist_whitelist_overrides() {
        auto out_path = whitelist_override_path();
        std::error_code ec;
        std::filesystem::create_directories(out_path.parent_path(), ec);
        if (ec) return false;

        std::ofstream out(out_path.string());
        if (!out.is_open()) return false;

        const auto& defaults = default_kit_wl_paths();
        std::vector<std::pair<std::string, std::string>> rows;
        rows.reserve(kit_wl_paths.size());
        for (const auto& kv : kit_wl_paths) {
            auto dit = defaults.find(kv.first);
            if (dit != defaults.end() && dit->second == kv.second) {
                continue; // not an override; omit
            }
            rows.emplace_back(kv.first, kv.second);
        }
        std::sort(rows.begin(), rows.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });

        std::vector<std::string> tombstones;
        tombstones.reserve(removed_whitelist_keys().size());
        for (const auto& key : removed_whitelist_keys()) {
            if (defaults.find(key) != defaults.end()) {
                tombstones.push_back(key);
            }
        }
        std::sort(tombstones.begin(), tombstones.end());

        if (rows.empty() && tombstones.empty()) {
            out.close();
            std::error_code rm_ec;
            std::filesystem::remove(out_path, rm_ec);
            return true;
        }

        for (const auto& kv : rows) {
            out << kv.first << '\t' << kv.second << '\n';
        }
        for (const auto& key : tombstones) {
            out << key << '\t' << '\n'; // tombstone removal marker
        }
        return true;
    }

    inline const std::unordered_map<std::string, std::string>& list_whitelist_paths() {
        load_whitelist_overrides_once();
        return kit_wl_paths;
    }

    // Split "kit:path", "kit:kit", or "path" into one or two specs.
    // Empty components and additional separators are malformed input: silently
    // dropping them can leave callers indexing an empty result or resolve a
    // different whitelist than the one the user supplied.
    static std::vector<std::string> parse_whitelist_specs(std::string const &field) {
        const std::string value = seq_utils::trim(field);
        if (value.empty()) {
            throw std::invalid_argument(
                "Malformed whitelist specification: value is empty");
        }

        const auto separator = value.find(':');
        if (separator == std::string::npos) {
            return {value};
        }
        if (value.find(':', separator + 1) != std::string::npos) {
            throw std::invalid_argument(
                "Malformed whitelist specification: expected at most two "
                "non-empty components separated by one ':'");
        }

        std::string first = seq_utils::trim(value.substr(0, separator));
        std::string second = seq_utils::trim(value.substr(separator + 1));
        if (first.empty() || second.empty()) {
            throw std::invalid_argument(
                "Malformed whitelist specification: both ':'-separated "
                "components must be non-empty");
        }

        std::vector<std::string> specs;
        specs.reserve(2);
        specs.push_back(std::move(first));
        specs.push_back(std::move(second));
        return specs;
    }

    inline std::string kit_to_path(const std::string &field) {
        namespace bfs = std::filesystem;
        load_whitelist_overrides_once();
    
        // 1) Split into up to two specs
        auto specs = whitelist_utils::parse_whitelist_specs(field);
        std::vector<std::string> resolved;
        resolved.reserve(specs.size());
    
        for (auto &spec : specs) {
            // if it’s a known kit key, map into resources/wl/…
            if (auto it = kit_wl_paths.find(spec); it != kit_wl_paths.end()) {
                bfs::path mapped(it->second); // e.g. "resources/wl/foo.csv.gz" or absolute path

                const auto default_it = default_kit_wl_paths().find(spec);
                const bool is_bundled_resource =
                    default_it != default_kit_wl_paths().end() &&
                    default_it->second == it->second;

                // explicit/custom filesystem path override
                if (!is_bundled_resource &&
                    (mapped.is_absolute() ||
                     (it->second.find(bfs::path::preferred_separator) != std::string::npos && bfs::exists(mapped)))) {
                    if (!bfs::exists(mapped))
                        throw std::runtime_error("Whitelist resource not found: " + mapped.string());
                    resolved.push_back(bfs::canonical(mapped).string());
                    continue;
                }
    
                // strip leading "resources/" if present
                auto iter = mapped.begin();
                if (iter!=mapped.end() && *iter=="resources") ++iter;
    
                // recombine the remainder
                bfs::path sub;
                for (; iter!=mapped.end(); ++iter) sub /= *iter;
    
                // now prepend the true resources folder
                bfs::path root = path_utils::find_resource_root();
                bfs::path full = root / sub;

                // Installed resources are authoritative.  RRAD_RESOURCE_CACHE
                // is only an overlay for default kit assets that were omitted
                // from the package (for example, large CRAN-external files).
                // Its root mirrors the resource tree, so `sub` is resolved as
                // <cache-root>/wl/<filename>.
                if (bfs::exists(full)) {
                    resolved.push_back(bfs::canonical(full).string());
                    continue;
                }

                if (is_bundled_resource) {
                    if (const char* cache_root = std::getenv("RRAD_RESOURCE_CACHE");
                        cache_root && *cache_root) {
                        std::error_code ec;
                        bfs::path cache(cache_root);
                        if (cache.is_relative()) cache = bfs::absolute(cache, ec);
                        if (!ec) {
                            const bfs::path cached = cache / sub;
                            if (bfs::exists(cached) && bfs::is_regular_file(cached)) {
                                resolved.push_back(bfs::canonical(cached).string());
                                continue;
                            }
                        }
                    }

                    // This prefix and resource-relative POSIX path are a
                    // machine-readable contract with rrad's R download layer.
                    throw std::runtime_error(
                        "RRAD_RESOURCE_MISSING:" + sub.generic_string());
                }

                throw std::runtime_error("Whitelist resource not found: " + full.string());
            }
            else {
                // 2b) otherwise treat spec as a path
                bfs::path p(spec);
                if (p.is_relative()) p = bfs::absolute(p);
                if (!bfs::exists(p))
                    throw std::runtime_error("Whitelist file not found: " + p.string());
                resolved.push_back(bfs::canonical(p).string());
            }
        }
    
        // 3) re-join with ":" if there were two pieces
        if (resolved.size() == 2) {
            return resolved[0] + ":" + resolved[1];
        } else if (resolved.size() == 1) {
            return resolved.front();
        }
        throw std::runtime_error(
            "Malformed whitelist specification: no whitelist was resolved");
    }

    // Check if it's a kit
    bool is_kit(const std::string &name){
        load_whitelist_overrides_once();
        if(kit_wl_paths.find(name) != kit_wl_paths.end()){
            return true;
        }
        return false;
    }

    //check if the whitelist we're importing is a bitlist or a whitelist
    bool check_if_bitlist(const std::string &field, bool verbose, size_t N = 10) {
        const std::string path = kit_to_path(field);

        // Inspect the first rows without assuming that a header is present.
        // Most bundled bitlists are headerless, while a few use a text header;
        // dropping row zero unconditionally also misclassified one-row lists.
        auto lines = streaming_utils::import_text(path, N + 1);
        if (lines.empty()) {
            if(verbose) RAD_CERR << "[check_if_bitlist] ERROR: no lines read from " << path << "\n";
            return false;
        }
        if(verbose) RAD_COUT << "[check_if_bitlist] debug: first " << lines.size() << " raw lines from " << path << ":\n";

        for (size_t i = 0; i < lines.size(); ++i) {
            if(verbose) RAD_COUT << "  [" << i << "]: " << lines[i] << "\n";
        }

        // Extract up to N tokens (before first comma/tab). A single leading
        // non-numeric token is accepted as a header only when the remaining
        // rows are numeric.
        std::vector<std::string> tokens;
        tokens.reserve(N);
        for (auto &ln : lines) {
          if (ln.empty()) continue;
          auto p = ln.find_first_of(",\t");
          std::string token = seq_utils::trim(
              p == std::string::npos ? ln : ln.substr(0, p));
          if (token.empty()) continue;
          tokens.push_back(std::move(token));
        }
    
        if (tokens.empty()) {
            return false;
        }

        auto is_unsigned_integer = [](const std::string& token) {
            return !token.empty() &&
                   std::all_of(token.begin(), token.end(),
                               [](unsigned char c) { return std::isdigit(c); });
        };

        size_t first_content = 0;
        if (!is_unsigned_integer(tokens.front())) {
            if (tokens.size() == 1) {
                return false;
            }
            first_content = 1;
        }
        const size_t content_end = std::min(tokens.size(), first_content + N);
        for (size_t i = first_content; i < content_end; ++i) {
            if (!is_unsigned_integer(tokens[i])) {
                return false;
            }
        }
        return first_content < content_end;
      }

    inline const std::string& get_whitelist_path(const std::string &kit) {
        load_whitelist_overrides_once();
        auto it = kit_wl_paths.find(kit);
        if (it == kit_wl_paths.end())
            throw std::invalid_argument("Unknown kit name: " + kit);
        return it->second;
    }

    inline void set_whitelist_path(const std::string &kit, const std::string &path) {
        load_whitelist_overrides_once();
        removed_whitelist_keys().erase(kit);
        kit_wl_paths[kit] = path;
        if (!persist_whitelist_overrides()) {
            RAD_CERR << "[WARNING] Failed to persist whitelist overrides to "
                      << whitelist_override_path().string() << "\n";
        }
    }

    inline bool remove_whitelist_path(const std::string &kit) {
        load_whitelist_overrides_once();
        bool removed = false;
        if (kit_wl_paths.erase(kit) > 0) {
            removed = true;
        }
        if (default_kit_wl_paths().find(kit) != default_kit_wl_paths().end()) {
            removed_whitelist_keys().insert(kit);
            removed = true;
        }
        if (removed && !persist_whitelist_overrides()) {
            RAD_CERR << "[WARNING] Failed to persist whitelist overrides to "
                      << whitelist_override_path().string() << "\n";
        }
        return removed;
    }

};

namespace config_utils {

    //config for read layout
    inline const std::unordered_map<std::string, std::string>& default_layout_files() {
        static const std::unordered_map<std::string, std::string> defaults = {
        {"five_prime", "resources/read_layout/five_prime_read_layout.csv" },
        {"sctagger", "resources/read_layout/sctagger_sim_read_layout.csv" },
        {"three_prime", "resources/read_layout/three_prime_read_layout.csv" },
        {"splitseq", "resources/read_layout/splitseq_read_layout.csv" },
        {"visium", "resources/read_layout/visium_three_prime_read_layout.csv" },
        {"visium_hd", "resources/read_layout/visium_hd_read_layout.csv" },
        {"nanopore_rapid_bc", "resources/read_layout/nanopore_bulk_rapid_bc_read_layout.csv" }
        };
        return defaults;
    }

    static std::unordered_map<std::string, std::string> layout_files = default_layout_files();

    inline std::filesystem::path layout_override_path() {
        if (const char* home = std::getenv("HOME"); home && *home) {
            return std::filesystem::path(home) / ".rad" / "layout_overrides.tsv";
        }
        return std::filesystem::current_path() / ".rad_layout_overrides.tsv";
    }

    inline std::unordered_set<std::string>& removed_layout_keys() {
        static std::unordered_set<std::string> removed;
        return removed;
    }

    inline bool& layout_overrides_loaded() {
        static bool loaded = false;
        return loaded;
    }

    inline void load_layout_overrides_once() {
        // rrad's embedded backend must be reproducible from its bundled
        // resources.  Do not let per-user RAD CLI state change that backend.
        if (const char* embedded_only = std::getenv("RAD_EMBEDDED_ONLY");
            embedded_only && *embedded_only) {
            return;
        }
        if (layout_overrides_loaded()) return;
        layout_overrides_loaded() = true;

        std::ifstream in(layout_override_path().string());
        if (!in.is_open()) return;

        std::string line;
        while (std::getline(in, line)) {
            std::string trimmed_line = seq_utils::trim(line);
            if (trimmed_line.empty() || trimmed_line[0] == '#') continue;

            auto tab = line.find('\t');
            if (tab == std::string::npos) continue;

            std::string key = seq_utils::trim(line.substr(0, tab));
            std::string value = seq_utils::trim(line.substr(tab + 1));
            if (key.empty()) continue;

            if (value.empty()) {
                layout_files.erase(key);
                removed_layout_keys().insert(key);
            } else {
                layout_files[key] = value;
                removed_layout_keys().erase(key);
            }
        }
    }

    inline bool persist_layout_overrides() {
        auto out_path = layout_override_path();
        std::error_code ec;
        std::filesystem::create_directories(out_path.parent_path(), ec);
        if (ec) return false;

        std::ofstream out(out_path.string());
        if (!out.is_open()) return false;

        const auto& defaults = default_layout_files();
        std::vector<std::pair<std::string, std::string>> rows;
        rows.reserve(layout_files.size());
        for (const auto& kv : layout_files) {
            auto dit = defaults.find(kv.first);
            if (dit != defaults.end() && dit->second == kv.second) {
                continue; // not an override; omit
            }
            rows.emplace_back(kv.first, kv.second);
        }
        std::sort(rows.begin(), rows.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });

        std::vector<std::string> tombstones;
        tombstones.reserve(removed_layout_keys().size());
        for (const auto& key : removed_layout_keys()) {
            if (defaults.find(key) != defaults.end()) {
                tombstones.push_back(key);
            }
        }
        std::sort(tombstones.begin(), tombstones.end());

        if (rows.empty() && tombstones.empty()) {
            out.close();
            std::error_code rm_ec;
            std::filesystem::remove(out_path, rm_ec);
            return true;
        }

        for (const auto& kv : rows) {
            out << kv.first << '\t' << kv.second << '\n';
        }
        for (const auto& key : tombstones) {
            out << key << '\t' << '\n'; // tombstone removal marker
        }
        return true;
    }

    inline const std::unordered_map<std::string, std::string>& list_read_layouts() {
        load_layout_overrides_once();
        return layout_files;
    }

    // check if the layout is a custom one
    inline bool check_if_custom_rl(const std::string &type) {
        load_layout_overrides_once();
        return layout_files.find(type) == layout_files.end();
    }
    
   inline std::string get_read_layout(const std::string &type) {
    namespace bfs = std::filesystem;
    load_layout_overrides_once();
    // 1) If `type` looks like an absolute or contains a path‐sep, treat it as a real file:
    bfs::path p(type);
    if (p.is_absolute() ||
        type.find(bfs::path::preferred_separator) != std::string::npos)
    {
      if (!bfs::exists(p))
        throw std::runtime_error("Layout file not found: " + type);
      return bfs::canonical(p).string();
    }

    // 2) Otherwise look up our map
    auto it = layout_files.find(type);
    if (it == layout_files.end())
      throw std::invalid_argument("Unknown layout key: " + type);

    // 3) If mapped value is a direct path and exists, use it directly.
    bfs::path mapped(it->second);
    if (mapped.is_absolute() ||
        (it->second.find(bfs::path::preferred_separator) != std::string::npos && bfs::exists(mapped))) {
        if (!bfs::exists(mapped))
            throw std::runtime_error("Missing layout resource: " + mapped.string());
        return bfs::canonical(mapped).string();
    }

    // 4) Take the stored relative‐to‐resources path
    //    e.g. it->second == "splitseq_read_layout.csv"   OR
    //         it->second == "read_layout/splitseq_read_layout.csv" OR
    //         it->second == "resources/read_layout/splitseq_read_layout.csv"
    
    bfs::path orig_rel(it->second);

    // strip any leading "resources/" component:
    auto iter = orig_rel.begin();
    if (iter!=orig_rel.end() && *iter == "resources") {
      ++iter;
    }
    // now build a new path from whatever remains
    bfs::path rel;
    for (; iter!=orig_rel.end(); ++iter) {
      rel /= *iter;
    }

    // 5) Join against the real resources folder
    bfs::path root = path_utils::find_resource_root();
    bfs::path full = root / rel;

    if (!bfs::exists(full))
      throw std::runtime_error("Missing layout resource: " + full.string());
    return bfs::canonical(full).string();
  }

    inline void save_read_layout(const std::string &type, const std::string &path) {
        load_layout_overrides_once();
        removed_layout_keys().erase(type);
        layout_files[type] = path;
        if (!persist_layout_overrides()) {
            RAD_CERR << "[WARNING] Failed to persist layout overrides to "
                      << layout_override_path().string() << "\n";
        }
    }
    
    inline bool remove_read_layout(const std::string &type) {
        load_layout_overrides_once();
        bool removed = false;
        if (layout_files.erase(type) > 0) {
            removed = true;
        }
        if (default_layout_files().find(type) != default_layout_files().end()) {
            removed_layout_keys().insert(type);
            removed = true;
        }
        if (removed && !persist_layout_overrides()) {
            RAD_CERR << "[WARNING] Failed to persist layout overrides to "
                      << layout_override_path().string() << "\n";
        }
        return removed;
    }

};

namespace plotting_utils {
    struct peak {
        size_t idx;   // index in y
        double x;     // x[idx]
        double y;     // y[idx]
    };

    inline std::vector<double> moving_average(
        const std::vector<double>& y, int w = 5) {
        if (w <= 1 || (int)y.size() <= w) return y;
        std::vector<double> out(y.size());
        double acc = 0.0;
        for (int i = 0; i < (int)y.size(); ++i) {
            acc += y[i];
            if (i >= w) acc -= y[i - w];
            if (i >= w-1) out[i - (w-1)/2] = acc / w;
        }
        int half = w/2;
        for (int i = 0; i < half; ++i) out[i] = y[i];
        for (int i = (int)y.size()-half; i < (int)y.size(); ++i) out[i] = y[i];
        return out;
    }

    // strict local maxima + simple spacing and height threshold
    inline std::vector<peak> find_peaks(
        const std::vector<double>& x, const std::vector<double>& y_in,
        double min_height = -INFINITY, int min_distance = 1, bool smooth = true, int ma_window = 7) {
        const std::vector<double> y = smooth ? moving_average(y_in, ma_window) : y_in;
        std::vector<size_t> cand;
        if (y.size() >= 3) {
            for (size_t i = 1; i + 1 < y.size(); ++i) {
                if (y[i] > y[i-1] && y[i] > y[i+1] && y[i] >= min_height)
                    cand.push_back(i);
            }
        }
        if (cand.empty()) return {};

        // enforce min_distance by greedy “keep highest in window”
        std::vector<size_t> kept;
        size_t j = 0;
        while (j < cand.size()) {
            size_t best = cand[j];
            size_t k = j + 1;
            while (k < cand.size() && (int)(cand[k] - cand[j]) < min_distance) {
                if (y[cand[k]] > y[best]) best = cand[k];
                ++k;
            }
            kept.push_back(best);
            j = k;
        }

        std::vector<peak> out;
        out.reserve(kept.size());
        for (auto i : kept) out.push_back(peak{i, x[i], y[i]});
        return out;
    }

    // first “flat” local minimum after a given peak (index in y)
    inline std::optional<peak> first_flat_min_after(
        const std::vector<double>& x, const std::vector<double>& y,
        size_t start_idx, int slope_window = 7, double slope_tol = 1e-3) {
        auto avg_abs_slope = [&](size_t i)->double {
            int half = std::max(1, slope_window/2);
            int L = std::max<int>(1, (int)i - half);
            int R = std::min<int>((int)y.size()-2, (int)i + half);
            double s = 0.0; int c = 0;
            for (int k = L; k <= R; ++k) { s += std::fabs(y[k+1]-y[k]); ++c; }
            return (c>0) ? s / c : 0.0;
        };

        for (size_t i = start_idx + 1; i + 1 < y.size(); ++i) {
            if (y[i] <= y[i-1] && y[i] <= y[i+1]) {
                if (avg_abs_slope(i) <= slope_tol)
                    return peak{i, x[i], y[i]};
            }
        }
        return std::nullopt;
    }
};

namespace gmm_utils {
    struct gmm_comp {
        double weight; // π_k
        double mean;   // μ_k
        double var;    // σ_k^2
    };

    struct gmm_fit {
        std::vector<gmm_comp> comps;
        int K = 0;
        double bic = std::numeric_limits<double>::infinity();
        // evaluated on x grid for easy overlay
        std::vector<std::vector<double>> comp_y; // K x x.size()
        std::vector<double> mix_y;               // sum_k comp_y[k]
    };

    struct gmm_sweep {
        std::vector<gmm_fit> fits;
        int best_idx = -1;
    };

    // normalize y to sum=1; return weights
    inline std::vector<double> normalize_weights(const std::vector<double>& y) {
        double s = 0.0;
        for (double v : y) s += v;
        std::vector<double> w(y.size());
        if (s <= 0) return w; // all zeros -> all zeros
        for (size_t i = 0; i < y.size(); ++i) w[i] = y[i] / s;
        return w;
    }

    // effective N for BIC from weights
    inline double effective_n(const std::vector<double>& w) {
        double s1 = 0.0, s2 = 0.0;
        for (double wi : w) { s1 += wi; s2 += wi*wi; }
        if (s2 == 0.0) return (double)w.size();
        return (s1*s1) / s2;
    }

    inline double normal_pdf_1d(double x, double m, double v) {
        if (v <= 0) v = 1e-12;
        constexpr double two_pi = boost::math::constants::two_pi<double>();
        double inv = 1.0 / std::sqrt(two_pi * v);
        double z = (x - m);
        return inv * std::exp(-0.5 * z * z / v);
    }

    // simple weighted quantiles in [0,1] of x under weights w (assumes x sorted ascending)
    inline double weighted_quantile(const std::vector<double>& x, const std::vector<double>& w, double q
    ) {
        double acc = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            acc += w[i];
            if (acc >= q) return x[i];
        }
        return x.back();
    }

    inline gmm_fit fit_gmm_weighted(const std::vector<double>& x, const std::vector<double>& y,
                                    int Kmin = 1, int Kmax = 5, int max_iters = 200, double tol = 1e-6,
                                    double min_var = 1e-6) {
        // Preconditions
        gmm_fit best;
        if (x.size() != y.size() || x.size() < 3 || Kmin < 1 || Kmax < Kmin) return best;

        // Ensure x is ascending (your sampling already is)
        // Build normalized weights (sum to 1)
        std::vector<double> w = normalize_weights(y);
        const size_t N = x.size();
        const double Neff = effective_n(w);

        // global weighted mean/var (for init)
        double mu0 = 0.0, s0 = 0.0;
        for (size_t i = 0; i < N; ++i) mu0 += w[i] * x[i];
        for (size_t i = 0; i < N; ++i) { double d = x[i] - mu0; s0 += w[i] * d * d; }
        s0 = std::max(s0, min_var);

        auto run_em = [&](int K)->gmm_fit {
            // ---- init by weighted quantiles for means; equal weights; same var ----
            std::vector<gmm_comp> comps(K);
            for (int k = 0; k < K; ++k) {
                double q = (k + 0.5) / K; // mid-quantiles
                comps[k].mean = weighted_quantile(x, w, q);
                comps[k].weight = 1.0 / K;
                comps[k].var = s0;
            }

            std::vector<std::vector<double>> resp(K, std::vector<double>(N)); // responsibilities r_ik
            double prev_ll = -std::numeric_limits<double>::infinity();

            for (int iter = 0; iter < max_iters; ++iter) {
                // E-step: r_ik = π_k N(x_i|μ_k,σ^2_k) / sum_j
                for (size_t i = 0; i < N; ++i) {
                    double denom = 0.0;
                    for (int k = 0; k < K; ++k) {
                        double v = comps[k].weight * normal_pdf_1d(x[i], comps[k].mean, comps[k].var);
                        resp[k][i] = v;
                        denom += v;
                    }
                    if (denom <= 0.0) {
                        // if all underflowed, assign uniform
                        for (int k = 0; k < K; ++k) resp[k][i] = 1.0 / K;
                    } else {
                        double inv = 1.0 / denom;
                        for (int k = 0; k < K; ++k) resp[k][i] *= inv;
                    }
                }

                // M-step using weighted x-weights: w_i are sample weights for KDE grid
                for (int k = 0; k < K; ++k) {
                    // effective component mass
                    double Nk = 0.0;
                    for (size_t i = 0; i < N; ++i) Nk += w[i] * resp[k][i];
                    Nk = std::max(Nk, 1e-12);

                    // mean
                    double mk = 0.0;
                    for (size_t i = 0; i < N; ++i) mk += w[i] * resp[k][i] * x[i];
                    mk /= Nk;

                    // variance
                    double vk = 0.0;
                    for (size_t i = 0; i < N; ++i) {
                        double d = x[i] - mk;
                        vk += w[i] * resp[k][i] * d * d;
                    }
                    vk = std::max(vk / Nk, min_var);

                    comps[k].mean = mk;
                    comps[k].var  = vk;
                    comps[k].weight = Nk; // temporary; normalize below
                }
                // normalize weights
                double sum_pi = 0.0;
                for (int k = 0; k < K; ++k) sum_pi += comps[k].weight;
                if (sum_pi <= 0) {
                    for (int k = 0; k < K; ++k) comps[k].weight = 1.0 / K;
                } else {
                    for (int k = 0; k < K; ++k) comps[k].weight /= sum_pi;
                }

                // log-likelihood of weighted grid: sum_i w_i log sum_k π_k N(x_i)
                double ll = 0.0;
                for (size_t i = 0; i < N; ++i) {
                    double s = 0.0;
                    for (int k = 0; k < K; ++k)
                        s += comps[k].weight * normal_pdf_1d(x[i], comps[k].mean, comps[k].var);
                    if (s <= 0) s = 1e-300;
                    ll += w[i] * std::log(s);
                }
                if (std::abs(ll - prev_ll) < tol) break;
                prev_ll = ll;
            }

            // BIC: -2*LL + p*log(Neff), with p = (K-1) + K*2 (weights-1 + mean + var per comp)
            // (mixture weights have K-1 free params)
            double ll = 0.0;
            for (size_t i = 0; i < N; ++i) {
                double s = 0.0;
                for (int k = 0; k < K; ++k)
                    s += comps[k].weight * normal_pdf_1d(x[i], comps[k].mean, comps[k].var);
                if (s <= 0) s = 1e-300;
                ll += w[i] * std::log(s);
            }
            int p = (K - 1) + 2 * K;
            double bic = -2.0 * ll + p * std::log(std::max(Neff, 1.0));

            // build comp_y and mix_y on the same grid (scale to density units)
            std::vector<std::vector<double>> comp_y(K, std::vector<double>(N));
            std::vector<double> mix_y(N, 0.0);
            for (int k = 0; k < K; ++k) {
                for (size_t i = 0; i < N; ++i) {
                    double v = comps[k].weight * normal_pdf_1d(x[i], comps[k].mean, comps[k].var);
                    comp_y[k][i] = v;
                    mix_y[i] += v;
                }
            }

            return gmm_fit{ std::move(comps), K, bic, std::move(comp_y), std::move(mix_y) };
        };

        // sweep K and keep best BIC
        for (int K = Kmin; K <= Kmax; ++K) {
            gmm_fit cur = run_em(K);
            if (cur.bic < best.bic) best = std::move(cur);
        }
        return best;
    }

    // overload that accepts initial means (<= K seeds)
    inline gmm_fit fit_gmm_weighted(const std::vector<double>& x, const std::vector<double>& y,
                                int Kmin, int Kmax, int max_iters, double tol, double min_var,
                                const std::vector<double>* init_means) {
        // nullptr => quantile init
        gmm_fit best;
        // basic guards
        if (x.size() != y.size() || x.size() < 3 || Kmin < 1 || Kmax < Kmin) return best;

        // weights from density y (sum to 1)
        std::vector<double> w = normalize_weights(y);
        const size_t N = x.size();
        const double Neff = effective_n(w);

        // global weighted stats (for variance init)
        double mu0 = 0.0, s0 = 0.0;
        for (size_t i = 0; i < N; ++i) mu0 += w[i] * x[i];
        for (size_t i = 0; i < N; ++i) { double d = x[i] - mu0; s0 += w[i] * d * d; }
        s0 = std::max(s0, min_var);

        auto run_em = [&](int K)->gmm_fit {
            std::vector<gmm_comp> comps(K);

            // ---- initialization ----
            if (init_means && !init_means->empty()) {
                const int S = std::min<int>((int)init_means->size(), K);
                for (int k = 0; k < S; ++k) comps[k].mean = (*init_means)[k];
                for (int k = S; k < K; ++k) {
                    double q = (k + 0.5) / K;                  // mid-quantiles
                    comps[k].mean = weighted_quantile(x, w, q);
                }
            } else {
                for (int k = 0; k < K; ++k) {
                    double q = (k + 0.5) / K;
                    comps[k].mean = weighted_quantile(x, w, q);
                }
            }
            for (int k = 0; k < K; ++k) { comps[k].weight = 1.0 / K; comps[k].var = s0; }

            // responsibilities r_ik
            std::vector<std::vector<double>> resp(K, std::vector<double>(N));
            double prev_ll = -std::numeric_limits<double>::infinity();

            // ---- EM iterations ----
            for (int iter = 0; iter < max_iters; ++iter) {
                // E-step
                for (size_t i = 0; i < N; ++i) {
                    double denom = 0.0;
                    for (int k = 0; k < K; ++k) {
                        double v = comps[k].weight * normal_pdf_1d(x[i], comps[k].mean, comps[k].var);
                        resp[k][i] = v;
                        denom += v;
                    }
                    if (denom <= 0.0) {
                        const double u = 1.0 / K;
                        for (int k = 0; k < K; ++k) resp[k][i] = u;
                    } else {
                        double inv = 1.0 / denom;
                        for (int k = 0; k < K; ++k) resp[k][i] *= inv;
                    }
                }

                // M-step (weighted by w_i)
                for (int k = 0; k < K; ++k) {
                    double Nk = 0.0;
                    for (size_t i = 0; i < N; ++i) Nk += w[i] * resp[k][i];
                    Nk = std::max(Nk, 1e-12);

                    double mk = 0.0;
                    for (size_t i = 0; i < N; ++i) mk += w[i] * resp[k][i] * x[i];
                    mk /= Nk;

                    double vk = 0.0;
                    for (size_t i = 0; i < N; ++i) {
                        double d = x[i] - mk;
                        vk += w[i] * resp[k][i] * d * d;
                    }
                    vk = std::max(vk / Nk, min_var);

                    comps[k].mean = mk;
                    comps[k].var  = vk;
                    comps[k].weight = Nk; // temp; normalize below
                }
                // normalize mixture weights
                double sum_pi = 0.0;
                for (int k = 0; k < K; ++k) sum_pi += comps[k].weight;
                if (sum_pi <= 0.0) {
                    for (int k = 0; k < K; ++k) comps[k].weight = 1.0 / K;
                } else {
                    for (int k = 0; k < K; ++k) comps[k].weight /= sum_pi;
                }

                // weighted log-likelihood
                double ll = 0.0;
                for (size_t i = 0; i < N; ++i) {
                    double s = 0.0;
                    for (int k = 0; k < K; ++k)
                        s += comps[k].weight * normal_pdf_1d(x[i], comps[k].mean, comps[k].var);
                    if (s <= 0) s = 1e-300;
                    ll += w[i] * std::log(s);
                }
                if (std::abs(ll - prev_ll) < tol) break;
                prev_ll = ll;
            }

            // BIC
            double ll = 0.0;
            for (size_t i = 0; i < N; ++i) {
                double s = 0.0;
                for (int k = 0; k < K; ++k)
                    s += comps[k].weight * normal_pdf_1d(x[i], comps[k].mean, comps[k].var);
                if (s <= 0) s = 1e-300;
                ll += w[i] * std::log(s);
            }
            int p = (K - 1) + 2 * K; // (weights-1) + mean + var per comp
            double bic = -2.0 * ll + p * std::log(std::max(Neff, 1.0));

            // evaluate components and mixture on grid
            std::vector<std::vector<double>> comp_y(K, std::vector<double>(N));
            std::vector<double> mix_y(N, 0.0);
            for (int k = 0; k < K; ++k) {
                for (size_t i = 0; i < N; ++i) {
                    double v = comps[k].weight * normal_pdf_1d(x[i], comps[k].mean, comps[k].var);
                    comp_y[k][i] = v;
                    mix_y[i] += v;
                }
            }

            return gmm_fit{ std::move(comps), K, bic, std::move(comp_y), std::move(mix_y) };
        };

        // Sweep K and keep the best-BIC model
        for (int K = Kmin; K <= Kmax; ++K) {
            gmm_fit cur = run_em(K);
            if (cur.bic < best.bic) best = std::move(cur);
        }
        return best;
    }

    inline gmm_sweep fit_gmm_sweep(const std::vector<double>& x, const std::vector<double>& y,
                               int Kmin, int Kmax, const std::vector<double>* init_means = nullptr,
                               int max_iters = 200, double tol = 1e-6, double min_var = 1e-6) {
        gmm_sweep out;
        out.fits.reserve(Kmax - Kmin + 1);
        double best_bic = std::numeric_limits<double>::infinity();
        for (int K = Kmin; K <= Kmax; ++K) {
            auto f = fit_gmm_weighted(x, y, K, K, max_iters, tol, min_var, init_means);
            out.fits.push_back(std::move(f));
            if (out.fits.back().bic < best_bic) {
                best_bic = out.fits.back().bic;
                out.best_idx = (int)out.fits.size() - 1;
            }
        }
        return out;
    }
};
