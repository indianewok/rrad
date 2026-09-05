#pragma once
// ============================================================================
// seqspec.hpp
//
// Translate a Lior Pachter lab "seqspec" YAML file
// (https://github.com/pachterlab/seqspec) into RAD's native read-layout CSV,
// so that `rad prep`/`rad demux` can accept a seqspec file directly via
// --layout. Detection + translation is transparent (auto-detect): when the
// resolved --layout argument looks like a seqspec YAML it is converted to an
// equivalent RAD layout CSV which is then fed, unchanged, to the existing
// ReadLayout::prep_new_layout() pipeline.
//
// What we translate
// -----------------
//   * The library molecule tree: `library_spec` (seqspec >= 0.2.0) or its
//     legacy name `assay_spec` (0.0.x/0.1.x). Each entry is a hierarchical
//     `Region`. We pick the relevant modality region, flatten its leaf
//     regions in 5'->3' order, and emit one RAD layout row per kept leaf.
//   * Barcode `onlist`s become RAD whitelist references: recognized onlists
//     (e.g. 10x 3M-february-2018) map to RAD's built-in optimized kit keys;
//     otherwise we require the literal local onlist file. The embedded R core
//     never downloads arbitrary URLs found in a seqspec document.
//
// We intentionally DROP the flowcell adapters (illumina_p5/illumina_p7) and
// the sample indices (index5/index7): these sit outside the cDNA reads RAD
// models, and dropping them reproduces RAD's curated layouts exactly. RAD
// derives the reverse-complement strand and read direction itself, so only the
// forward (5'->3') molecule is emitted here; the seqspec `sequence_spec`
// (Read) tree and `strand` are not needed for the layout.
//
// This header relies on namespaces defined in misc_utils.hpp (seq_utils,
// path_utils) which is included earlier via rad_headers.h.
// ============================================================================

#include <fkYAML/node.hpp>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <unordered_set>
#include <vector>

#include <filesystem>

#include "embedded_console.hpp"

// Keep third-party and standard-library declarations above at global scope;
// only RAD's translator itself belongs to the isolated embedded namespace.
namespace rrad_rad_core {
namespace seqspec {

namespace bfs = std::filesystem;

// ----------------------------------------------------------------------------
// small string helpers
// ----------------------------------------------------------------------------
inline std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return s;
}

inline std::string trim(const std::string& s) { return seq_utils::trim(s); }

inline bool environment_flag_enabled(const char* name) {
    const char* raw = std::getenv(name);
    if (!raw || !*raw) return false;
    const std::string value = to_lower(trim(raw));
    if (value.empty() || value == "0" || value == "false" ||
        value == "no" || value == "off") {
        return false;
    }
    // Unknown non-empty values fail closed instead of unexpectedly enabling
    // network access.
    return true;
}

inline std::string read_file_text(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("seqspec: cannot open file: " + path);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

// Strip the known seqspec YAML node tags (!Assay, !Region, !Onlist, ...) that
// appear in 0.2.x/0.3.x specs, so the parser sees plain mappings. Uses a word
// boundary so free-text values are never corrupted.
inline std::string strip_yaml_tags(const std::string& text) {
    static const std::regex tag_re(
        R"(![ ]?(Assay|Read|Region|Onlist|File|SeqKit|LibKit|SeqProtocol|LibProtocol|SeqProtocols|LibProtocols|SeqSpec)\b)");
    return std::regex_replace(text, tag_re, "");
}

// CSV-escape a single field (quote if it contains a comma, quote or newline).
inline std::string csv_escape(const std::string& field) {
    if (field.find_first_of(",\"\n\r") == std::string::npos) return field;
    std::string out = "\"";
    for (char c : field) {
        if (c == '"') out += "\"\"";
        else out += c;
    }
    out += '"';
    return out;
}

// ----------------------------------------------------------------------------
// typed-ish node accessors (seqspec scalars may parse as string/int/bool)
// ----------------------------------------------------------------------------
inline bool node_has(const fkyaml::node& n, const char* key) {
    return n.is_mapping() && n.contains(key);
}

inline std::string get_str(const fkyaml::node& n, const char* key,
                           const std::string& def = "") {
    if (!node_has(n, key)) return def;
    const fkyaml::node& v = n[key];
    if (v.is_null()) return def;
    if (v.is_string()) return v.get_value<std::string>();
    if (v.is_integer()) return std::to_string(v.get_value<std::int64_t>());
    if (v.is_boolean()) return v.get_value<bool>() ? "true" : "false";
    if (v.is_float_number()) {
        std::ostringstream ss;
        ss << v.get_value<double>();
        return ss.str();
    }
    return def;
}

inline std::optional<long long> get_int(const fkyaml::node& n, const char* key) {
    if (!node_has(n, key)) return std::nullopt;
    const fkyaml::node& v = n[key];
    if (v.is_integer()) return static_cast<long long>(v.get_value<std::int64_t>());
    if (v.is_string()) {
        try {
            return std::stoll(trim(v.get_value<std::string>()));
        } catch (...) {
            return std::nullopt;
        }
    }
    return std::nullopt;
}

// ----------------------------------------------------------------------------
// detection
// ----------------------------------------------------------------------------

// Heuristic: does this path look like a seqspec YAML (vs a RAD layout CSV)?
//   * .yaml / .yml extension  -> treat as seqspec candidate
//   * .csv / .tsv / .txt      -> not seqspec
//   * otherwise               -> sniff the first few KB for seqspec markers
inline bool is_seqspec_file(const std::string& path) {
    bfs::path p(path);
    std::string ext = to_lower(p.extension().string());
    if (ext == ".yaml" || ext == ".yml") return true;
    if (ext == ".csv" || ext == ".tsv" || ext == ".txt") return false;

    std::ifstream in(path, std::ios::binary);
    if (!in) return false;
    std::string head(8192, '\0');
    in.read(&head[0], static_cast<std::streamsize>(head.size()));
    head.resize(static_cast<size_t>(in.gcount()));
    return head.find("seqspec_version") != std::string::npos ||
           head.find("library_spec") != std::string::npos ||
           head.find("assay_spec") != std::string::npos;
}

// ----------------------------------------------------------------------------
// onlist -> whitelist resolution
// ----------------------------------------------------------------------------

// Recognize common barcode onlists and return the matching RAD kit key, or ""
// if unknown. Matched on a normalized substring of the onlist filename/url.
inline std::string onlist_to_kit(const std::string& filename_or_url) {
    bfs::path p(filename_or_url);
    std::string base = to_lower(p.filename().string());
    struct Entry { const char* needle; const char* kit; };
    // Order: more specific substrings first.
    static const Entry table[] = {
        {"3m-february-2018", "10x_3v3"},
        {"3m-3pgex",         "10x_3v4"},
        {"3m-5pgex",         "10x_5v3"},
        {"737k-august-2016", "10x_3v2"},
        {"visium-v5",        "10x_Vis_V5"},
        {"visium-v4",        "10x_Vis_V4"},
        {"visium-v3",        "10x_Vis_V3"},
        {"visium-v2",        "10x_Vis_V2"},
        {"visium-v1",        "10x_Vis_V1"},
    };
    for (const auto& e : table) {
        if (base.find(e.needle) != std::string::npos) return e.kit;
    }
    return "";
}

inline bfs::path seqspec_cache_dir() {
    if (const char* home = std::getenv("HOME"); home && *home)
        return bfs::path(home) / ".rad" / "seqspec_cache";
    return bfs::current_path() / ".rad_seqspec_cache";
}

inline std::string stable_url_key(const std::string& url) {
    // Stable FNV-1a key: this is a cache identity, not a security primitive.
    std::uint64_t value = UINT64_C(14695981039346656037);
    for (unsigned char byte : url) {
        value ^= byte;
        value *= UINT64_C(1099511628211);
    }
    std::ostringstream out;
    out << std::hex << std::setfill('0') << std::setw(16) << value;
    return out.str();
}

inline std::string safe_onlist_basename(const std::string& url,
                                        const std::string& filename) {
    std::string source = url;
    const auto fragment = source.find_first_of("?#");
    if (fragment != std::string::npos) source.resize(fragment);
    std::string base = bfs::path(source).filename().string();
    if (base.empty() || base == "/" || base == "." || base == "..") {
        base = bfs::path(filename).filename().string();
    }
    std::string safe;
    safe.reserve(base.size());
    for (unsigned char c : base) {
        if (std::isalnum(c) || c == '.' || c == '_' || c == '-') {
            safe.push_back(static_cast<char>(c));
        } else {
            safe.push_back('_');
        }
    }
    while (!safe.empty() && safe.front() == '.') safe.erase(safe.begin());
    return safe.empty() ? "onlist.txt" : safe;
}

// Download a remote onlist into the cache (best-effort, via curl/wget). Returns
// the cached path. Throws with actionable guidance on failure.
inline std::string download_remote_onlist(const std::string& url,
                                          const std::string& filename,
                                          bool verbose) {
    bfs::path cache = seqspec_cache_dir();
    std::error_code ec;
    bfs::create_directories(cache, ec);
    if (ec) {
        throw std::runtime_error("seqspec: cannot create onlist cache: " +
                                 ec.message());
    }

    const std::string base = safe_onlist_basename(url, filename);
    bfs::path target = cache / (stable_url_key(url) + "-" + base);

    if (bfs::is_regular_file(target, ec) && !ec &&
        bfs::file_size(target, ec) > 0 && !ec) {
        if (verbose)
            RAD_COUT << "[seqspec] using cached onlist: " << target.string() << "\n";
        return bfs::canonical(target).string();
    }

    if (environment_flag_enabled("RRAD_OFFLINE")) {
        throw std::runtime_error(
            "seqspec: remote onlist is unavailable because RRAD_OFFLINE is "
            "enabled. Download it separately and point the onlist at a local "
            "file, or supply --custom-whitelist/--kit.");
    }

    if (verbose)
        RAD_COUT << "[seqspec] downloading remote onlist: " << url << "\n";

    auto sh_quote = [](const std::string& s) {
        std::string q = "'";
        for (char c : s) {
            if (c == '\'') q += "'\\''";
            else q += c;
        }
        q += "'";
        return q;
    };
    static std::atomic<unsigned> download_counter{0};
    bfs::path tmp = target.string() + ".part." +
                    std::to_string(static_cast<long>(::getpid())) + "." +
                    std::to_string(download_counter.fetch_add(1));
    std::string curl = "curl -fsSL -o " + sh_quote(tmp.string()) + " " + sh_quote(url);
    int rc = std::system(curl.c_str());
    if (rc != 0) {
        std::string wget = "wget -q -O " + sh_quote(tmp.string()) + " " + sh_quote(url);
        rc = std::system(wget.c_str());
    }
    if (rc != 0 || !bfs::is_regular_file(tmp, ec) || ec ||
        bfs::file_size(tmp, ec) == 0 || ec) {
        bfs::remove(tmp, ec);
        throw std::runtime_error(
            "seqspec: failed to download remote onlist '" + url +
            "'. Pre-download it and either point the onlist at a local file or "
            "supply the whitelist via --custom-whitelist/--kit.");
    }
    bfs::rename(tmp, target, ec);
    if (ec) {
        // A concurrent process may have completed the same URL first.
        std::error_code check_error;
        if (bfs::is_regular_file(target, check_error) && !check_error &&
            bfs::file_size(target, check_error) > 0 && !check_error) {
            bfs::remove(tmp, check_error);
        } else {
            bfs::remove(tmp, check_error);
            throw std::runtime_error("seqspec: cannot finalize download: " +
                                     ec.message());
        }
    }
    return bfs::canonical(target).string();
}

// Resolve an onlist node to a RAD whitelist column value (a kit key when
// recognized, else a resolved file path). `spec_dir` is the directory of the
// seqspec file, used to resolve relative local paths.
inline std::string resolve_onlist(const fkyaml::node& onlist,
                                   const bfs::path& spec_dir, bool verbose) {
    std::string filename = get_str(onlist, "filename");
    std::string url      = get_str(onlist, "url");
    std::string urltype  = to_lower(get_str(onlist, "urltype"));
    std::string location = to_lower(get_str(onlist, "location"));

    // Explicit onlists are authoritative. Do not infer a different built-in
    // whitelist merely because a user-provided filename contains a familiar
    // 10x substring.
    bool remote = urltype == "http" || urltype == "https" || urltype == "ftp" ||
                  location == "remote" || url.rfind("http://", 0) == 0 ||
                  url.rfind("https://", 0) == 0 || url.rfind("ftp://", 0) == 0;

    if (remote) {
        if (url.empty())
            throw std::runtime_error(
                "seqspec: remote onlist '" + filename + "' has no url to download.");
        if (environment_flag_enabled("RRAD_OFFLINE")) {
            throw std::runtime_error(
                "seqspec: remote onlist cannot be used while RRAD_OFFLINE is "
                "enabled. Download it separately, update the seqspec onlist "
                "to reference that local file, or supply the whitelist with "
                "custom_whitelist/kit.");
        }
        if (environment_flag_enabled("RAD_EMBEDDED_ONLY")) {
            throw std::runtime_error(
                "seqspec: remote onlists are disabled in embedded rrad. "
                "Download the whitelist separately, update the seqspec onlist "
                "to reference that local file, or supply the whitelist with "
                "custom_whitelist/kit.");
        }
        return download_remote_onlist(url, filename, verbose);
    }

    // local: prefer `url` when it names a file; treat "./"-style url as a dir.
    std::string candidate = url;
    if (candidate.empty() || candidate == "./" || candidate == ".")
        candidate = filename;
    bfs::path p(candidate);
    if (p.is_relative()) p = spec_dir / p;
    if (bfs::exists(p) && bfs::is_directory(p)) p = p / filename;
    if (bfs::exists(p)) return bfs::canonical(p).string();

    bfs::path alt = spec_dir / filename;
    if (bfs::exists(alt)) return bfs::canonical(alt).string();

    throw std::runtime_error(
        "seqspec: local onlist not found (tried '" + p.string() + "' and '" +
        alt.string() + "'). Place the barcode file next to the spec or supply "
        "--custom-whitelist/--kit.");
}

// ----------------------------------------------------------------------------
// region -> RAD layout row mapping
// ----------------------------------------------------------------------------
struct LayoutRow {
    std::string id;
    std::string seq;
    std::string expected_length;
    std::string type;   // static | variable
    std::string cls;    // RAD global class
    std::string whitelist;
    std::string flags;
};

// region_types we deliberately omit (flowcell adapters + sample indices).
inline bool is_dropped_region_type(const std::string& rt_lower) {
    return rt_lower == "illumina_p5" || rt_lower == "illumina_p7" ||
           rt_lower == "index5" || rt_lower == "index7";
}

// Compute the RAD expected_length string ("", "N", or "min-max") for a
// variable element from its min/max len and N/X sequence run.
inline std::string variable_length_field(const fkyaml::node& region) {
    auto mn = get_int(region, "min_len");
    auto mx = get_int(region, "max_len");
    if (mn && mx && *mn > 0 && *mx > 0) {
        if (*mn == *mx) return std::to_string(*mx);
        return std::to_string(*mn) + "-" + std::to_string(*mx);
    }
    std::string seq = trim(get_str(region, "sequence"));
    if (!seq.empty()) return std::to_string(seq.size());
    if (mx && *mx > 0) return std::to_string(*mx);
    return "";
}

// Translate a single leaf region into a RAD layout row. Returns nullopt when
// the region should be skipped (dropped type / container / empty).
inline std::optional<LayoutRow> region_to_row(const fkyaml::node& region,
                                               const bfs::path& spec_dir,
                                               bool verbose,
                                               bool suppress_onlists,
                                               std::vector<std::string>& dropped) {
    std::string region_id   = trim(get_str(region, "region_id"));
    std::string region_type = trim(get_str(region, "region_type"));
    std::string seq_type     = to_lower(trim(get_str(region, "sequence_type")));
    std::string rt           = to_lower(region_type);
    std::string sequence     = trim(get_str(region, "sequence"));

    if (region_id.empty()) region_id = region_type.empty() ? "element" : region_type;

    if (is_dropped_region_type(rt)) {
        dropped.push_back(region_id + " (" + region_type + ")");
        return std::nullopt;
    }

    LayoutRow row;
    row.id = region_id;

    // ---- poly-tails: RAD synthesizes the regex from the class name ----
    if (rt == "poly_a") {
        row.id = "poly_a"; row.type = "static"; row.cls = "poly_a"; return row;
    }
    if (rt == "poly_t") {
        row.id = "poly_t"; row.type = "static"; row.cls = "poly_t"; return row;
    }

    // ---- barcodes / UMIs / onlists (variable elements) ----
    if (seq_type == "onlist" || rt == "barcode") {
        row.type = "variable";
        row.cls = "barcode";
        row.expected_length = variable_length_field(region);
        if (node_has(region, "onlist") && !region["onlist"].is_null() &&
            region["onlist"].is_mapping()) {
            if (suppress_onlists) {
                if (verbose) {
                    RAD_COUT
                        << "[seqspec] barcode '" << region_id
                        << "' onlist suppressed by an explicit demux "
                           "whitelist override.\n";
                }
            } else {
                row.whitelist =
                    resolve_onlist(region["onlist"], spec_dir, verbose);
            }
        } else if (verbose) {
            RAD_COUT << "[seqspec] barcode '" << region_id
                      << "' has no onlist; whitelist left blank.\n";
        }
        return row;
    }
    if (rt == "umi") {
        row.type = "variable"; row.cls = "umi";
        row.expected_length = variable_length_field(region);
        return row;
    }

    // ---- biological payload -> RAD "read" ----
    static const std::unordered_set<std::string> payload = {
        "cdna", "gdna", "dna", "rna", "atac", "protein", "crispr", "tag",
        "hic", "methyl", "sgrna_target"};
    if (payload.count(rt)) {
        row.type = "variable"; row.cls = "read";
        // leave length open-ended (RAD's read element is unconstrained)
        return row;
    }

    // ---- fixed sequences (adapters / primers / linkers / TSO) -> static ----
    if (seq_type == "fixed" || !sequence.empty()) {
        // sequence may legitimately be all-N for masked static regions.
        row.type = "static";
        row.seq = sequence;
        if (rt == "linker" || rt == "custom_primer") row.cls = "linker";
        else if (rt == "bead_tso") row.cls = "tso";
        else row.cls = "";  // generic adapter/primer: empty class, like RAD
        if (row.seq.empty()) {
            // a "fixed" region with no literal sequence is unusable as an anchor
            if (verbose)
                RAD_COUT << "[seqspec] skipping fixed region '" << region_id
                          << "' with empty sequence.\n";
            return std::nullopt;
        }
        return row;
    }

    // ---- random non-payload (rare) -> treat as read ----
    if (seq_type == "random") {
        row.type = "variable"; row.cls = "read";
        return row;
    }

    if (verbose)
        RAD_COUT << "[seqspec] skipping unrecognized region '" << region_id
                  << "' (type=" << region_type << ", sequence_type=" << seq_type
                  << ").\n";
    return std::nullopt;
}

// Depth-first flatten of leaf regions in document (5'->3') order.
inline void flatten_regions(const fkyaml::node& region,
                            const bfs::path& spec_dir, bool verbose,
                            bool suppress_onlists,
                            std::vector<LayoutRow>& rows,
                            std::vector<std::string>& dropped) {
    bool has_children = node_has(region, "regions") &&
                        region["regions"].is_sequence() &&
                        !region["regions"].empty();
    if (has_children) {
        for (const auto& child : region["regions"])
            flatten_regions(child, spec_dir, verbose, suppress_onlists,
                            rows, dropped);
        return;
    }
    if (auto row = region_to_row(region, spec_dir, verbose,
                                 suppress_onlists, dropped))
        rows.push_back(std::move(*row));
}

// ----------------------------------------------------------------------------
// modality selection
// ----------------------------------------------------------------------------

// Pick the top-level library region to translate. Prefer an explicit RNA
// modality; otherwise the first region. Logs the choice when >1 exists.
inline const fkyaml::node& choose_modality_region(const fkyaml::node& lib_spec,
                                                  bool verbose) {
    if (!lib_spec.is_sequence() || lib_spec.empty())
        throw std::runtime_error("seqspec: library_spec/assay_spec is empty.");

    if (lib_spec.size() == 1) return lib_spec[0];

    const fkyaml::node* chosen = nullptr;
    for (const auto& r : lib_spec) {
        std::string id = to_lower(get_str(r, "region_id"));
        std::string rt = to_lower(get_str(r, "region_type"));
        if (id == "rna" || rt == "rna") { chosen = &r; break; }
    }
    if (!chosen) chosen = &lib_spec[0];
    if (verbose) {
        RAD_COUT << "[seqspec] multiple modalities present; using region '"
                  << get_str(*chosen, "region_id") << "'. Re-run on a single-"
                  << "modality spec to select a different one.\n";
    }
    return *chosen;
}

// ----------------------------------------------------------------------------
// top-level translation
// ----------------------------------------------------------------------------

// Translate a seqspec YAML file into RAD layout-CSV text.
inline std::string translate_file_to_layout_csv(const std::string& path,
                                                bool verbose,
                                                bool suppress_onlists = false) {
    bfs::path spec_path(path);
    if (spec_path.is_relative()) spec_path = bfs::absolute(spec_path);
    bfs::path spec_dir = spec_path.parent_path();

    std::string text = strip_yaml_tags(read_file_text(spec_path.string()));

    fkyaml::node root;
    try {
        root = fkyaml::node::deserialize(text.begin(), text.end());
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("seqspec: YAML parse error in '") +
                                 path + "': " + e.what());
    }
    if (!root.is_mapping())
        throw std::runtime_error("seqspec: top-level document is not a mapping: " + path);

    // library_spec (>=0.2.0) or legacy assay_spec (0.0.x/0.1.x)
    const char* spec_key = nullptr;
    if (root.contains("library_spec")) spec_key = "library_spec";
    else if (root.contains("assay_spec")) spec_key = "assay_spec";
    else
        throw std::runtime_error(
            "seqspec: no 'library_spec' or 'assay_spec' found in " + path +
            " (is this a seqspec file?)");

    if (verbose) {
        std::string ver = get_str(root, "seqspec_version", "?");
        std::string name = get_str(root, "name", get_str(root, "assay_id", "(unnamed)"));
        RAD_COUT << "[seqspec] translating '" << name << "' (seqspec_version "
                  << ver << ", key '" << spec_key << "')\n";
    }

    const fkyaml::node& lib_spec = root[spec_key];
    const fkyaml::node& modality = choose_modality_region(lib_spec, verbose);

    std::vector<LayoutRow> rows;
    std::vector<std::string> dropped;
    flatten_regions(modality, spec_dir, verbose, suppress_onlists,
                    rows, dropped);

    if (rows.empty())
        throw std::runtime_error(
            "seqspec: no usable read-layout elements found in " + path);

    bool has_barcode = std::any_of(rows.begin(), rows.end(),
                                   [](const LayoutRow& r) { return r.cls == "barcode"; });
    bool has_read = std::any_of(rows.begin(), rows.end(),
                                [](const LayoutRow& r) { return r.cls == "read"; });

    if (verbose) {
        RAD_COUT << "[seqspec] " << rows.size() << " element(s); mode="
                  << (has_barcode ? "single-cell" : "bulk") << "\n";
        if (!dropped.empty()) {
            RAD_COUT << "[seqspec] dropped flowcell/index regions:";
            for (const auto& d : dropped) RAD_COUT << " " << d;
            RAD_COUT << "\n";
        }
        if (!has_read)
            RAD_COUT << "[seqspec] WARNING: no biological read (cdna/gdna/...) "
                         "region found.\n";
    }

    // ---- emit RAD layout CSV ----
    std::ostringstream out;
    // Title row encodes the sequencing mode for ReadLayout::get_rl_mode().
    out << (has_barcode ? "Read Layout" : "Read Layout:bulk") << ",,,,,,\n";
    out << "id,seq,expected_length,type,class,whitelist,flags\n";
    for (const auto& r : rows) {
        out << csv_escape(r.id) << "," << csv_escape(r.seq) << ","
            << csv_escape(r.expected_length) << "," << csv_escape(r.type) << ","
            << csv_escape(r.cls) << "," << csv_escape(r.whitelist) << ","
            << csv_escape(r.flags) << "\n";
    }
    return out.str();
}

// Translate `path` (a seqspec file) into a temporary RAD layout CSV and return
// its path. The caller is responsible for removing the temp file once
// prep_new_layout() has consumed it.
inline std::string convert_to_temp_layout_csv(const std::string& path,
                                              bool verbose,
                                              bool suppress_onlists = false) {
    std::string csv = translate_file_to_layout_csv(
        path, verbose, suppress_onlists);

    std::error_code ec;
    const bfs::path temp_dir = bfs::temp_directory_path(ec);
    if (ec) {
        throw std::runtime_error("seqspec: cannot resolve temp directory: " +
                                 ec.message());
    }

    std::string pattern = (temp_dir / "rad_seqspec_XXXXXX").string();
    std::vector<char> writable_pattern(pattern.begin(), pattern.end());
    writable_pattern.push_back('\0');
    int fd = ::mkstemp(writable_pattern.data());
    if (fd < 0) {
        throw std::runtime_error(
            "seqspec: cannot create exclusive temp layout: " +
            std::string(std::strerror(errno)));
    }

    const bfs::path tmp(writable_pattern.data());
    size_t offset = 0;
    while (offset < csv.size()) {
        const ssize_t written = ::write(
            fd, csv.data() + offset, csv.size() - offset);
        if (written < 0 && errno == EINTR) continue;
        if (written <= 0) {
            const std::string message = std::strerror(errno);
            ::close(fd);
            fd = -1;
            bfs::remove(tmp, ec);
            throw std::runtime_error(
                "seqspec: cannot write temp layout: " + message);
        }
        offset += static_cast<size_t>(written);
    }
    if (::close(fd) != 0) {
        const std::string message = std::strerror(errno);
        fd = -1;
        bfs::remove(tmp, ec);
        throw std::runtime_error(
            "seqspec: cannot close temp layout: " + message);
    }
    fd = -1;
    if (verbose)
        RAD_COUT << "[seqspec] wrote translated layout to " << tmp.string() << "\n";
    return tmp.string();
}

}  // namespace seqspec
}  // namespace rrad_rad_core
