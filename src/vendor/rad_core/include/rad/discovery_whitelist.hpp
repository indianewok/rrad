#pragma once


namespace rad_discovery {
namespace discovery_whitelist_detail {

struct Catalog {
    std::string path;
    std::vector<std::string> aliases;
    std::set<int> declared_lengths;
    bool builtin_bitlist = false;
};

struct LoadedCatalog {
    scan_whitelist_filter members;
    std::set<int> lengths;
};

struct UnknownBitlistLength : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// gzread also accepts uncompressed input. Unlike igzstream's EOF interface,
// checking gzerror here distinguishes truncated/corrupt gzip data from EOF.
template <typename Consumer>
static void read_catalog_lines(const std::string& path, Consumer consume, const std::function<void()>& poll = {}) {
    gzFile input = gzopen(path.c_str(), "rb");
    if (!input)
        throw std::runtime_error("Cannot open whitelist catalog: " + path);
    try {
        char buffer[65536];
        std::string pending;
        std::size_t line_number = 0;
        for (;;) {
            if (poll) poll();
            const int count = gzread(input, buffer, sizeof(buffer));
            if (count <= 0) {
                int code = Z_OK;
                const char* message = gzerror(input, &code);
                if (count < 0 || (code != Z_OK && code != Z_STREAM_END))
                    throw std::runtime_error("Cannot read whitelist catalog " +
                                             path + ": " + message);
                break;
            }
            for (int i = 0; i < count; ++i) {
                if (buffer[i] == '\n') {
                    consume(pending, ++line_number);
                    pending.clear();
                } else {
                    if (pending.size() >= 65536)
                        throw std::runtime_error("Overlong whitelist row in " + path);
                    pending.push_back(buffer[i]);
                }
            }
        }
        if (!pending.empty()) consume(pending, ++line_number);
    } catch (...) {
        gzclose(input);
        throw;
    }
    if (gzclose(input) != Z_OK)
        throw std::runtime_error("Failed to close whitelist catalog: " + path);
}

static bool decimal_token(const std::string& token) {
    std::size_t first = (!token.empty() && (token[0] == '-' || token[0] == '+')) ? 1 : 0;
    return first < token.size() &&
           std::all_of(token.begin() + first, token.end(),
                       [](unsigned char c) { return std::isdigit(c); });
}

static bool catalog_header(std::string token) {
    std::transform(token.begin(), token.end(), token.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return token == "whitelist_bcs" || token == "barcode" || token == "barcodes" ||
           token == "bitlist" || token == "barcode_bits";
}

static LoadedCatalog load_catalog(const Catalog& catalog, const Options& options) {
    LoadedCatalog loaded;
    bool content_seen = false;
    bool header_seen = false;
    std::optional<bool> bitlist;
    read_catalog_lines(catalog.path, [&](std::string line, std::size_t number) {
        if (number == 1 && line.compare(0, 3, "\xEF\xBB\xBF") == 0)
            line.erase(0, 3);
        const auto delimiter = line.find_first_of(",\t");
        std::string token = seq_utils::trim(line.substr(0, delimiter));
        if (token.empty()) return;
        if (!content_seen && !header_seen && catalog_header(token)) {
            header_seen = true;
            return;
        }
        const bool numeric = decimal_token(token);
        if (!bitlist) bitlist = numeric;
        if (*bitlist != numeric)
            throw std::runtime_error("Mixed or malformed whitelist catalog " +
                                     catalog.path + " at row " + std::to_string(number));
        int64_seq encoded;
        if (numeric) {
            // An integer's leading A bases are not represented in its value.
            // Only an explicit, unambiguous catalog/layout contract supplies L.
            if (catalog.declared_lengths.size() != 1)
                throw UnknownBitlistLength("Bitlist has no unambiguous declared barcode length: " +
                                           catalog.path);
            const auto length = static_cast<uint16_t>(*catalog.declared_lengths.begin());
            if (!parse_packed_decimal(token, length, encoded, true))
                throw std::runtime_error("Invalid packed barcode in " + catalog.path +
                                         " at row " + std::to_string(number));
        } else {
            std::transform(token.begin(), token.end(), token.begin(),
                           [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
            if (token.empty() || token.size() > 64 ||
                !std::all_of(token.begin(), token.end(), [](char c) {
                    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
                }))
                throw std::runtime_error("Invalid DNA barcode in " + catalog.path +
                                         " at row " + std::to_string(number));
            encoded = int64_seq(token);
        }
        content_seen = true;
        loaded.lengths.insert(encoded.length);
        if (encoded.length <= 31) {
            const packed_scan_barcode key =
                (uint64_t{1} << (2 * encoded.length)) | encoded.as_uint64();
            loaded.members.packed_barcodes.insert(key);
        } else {
            loaded.members.long_barcodes.insert(encoded);
        }
    }, options.boundary_poll);
    if (!content_seen || loaded.members.empty())
        throw std::runtime_error("Whitelist catalog contains no valid barcodes: " + catalog.path);
    return loaded;
}

static std::map<std::string, Catalog> registered_catalogs(const Options& options) {
    std::map<std::string, Catalog> catalogs;
    if (!options.catalogs.empty()) {
        for (const auto& source : options.catalogs) {
            const std::string path = std::filesystem::canonical(source.path).string();
            auto& catalog = catalogs[path];
            catalog.path = path;
            catalog.aliases.push_back(source.name);
            if (source.length > 0) catalog.declared_lengths.insert(source.length);
        }
        return catalogs;
    }
    const auto& defaults = whitelist_utils::default_kit_wl_paths();
    const auto& registered = whitelist_utils::list_whitelist_paths();
    auto spatial_alias = [](std::string alias) {
        std::transform(alias.begin(), alias.end(), alias.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        return alias.rfind("10x_vis_", 0) == 0 || alias.rfind("visium", 0) == 0 ||
               alias.rfind("stereo", 0) == 0 || alias.rfind("curio", 0) == 0 ||
               alias.rfind("air_space", 0) == 0 || alias.rfind("spatial", 0) == 0;
    };
    // Spatial catalogs are outside this command's current identification scope.
    // Resolve their aliases first so a second alias cannot reintroduce the same file.
    std::set<std::string> spatial_paths;
    for (const auto& entry : registered) {
        if (!spatial_alias(entry.first)) continue;
        try { spatial_paths.insert(whitelist_utils::kit_to_path(entry.first)); }
        catch (const std::exception&) {} // An unsupported catalog need not be installed.
    }
    for (const auto& entry : registered) {
        if (spatial_alias(entry.first)) continue;
        try {
            const std::string path = whitelist_utils::kit_to_path(entry.first);
            if (spatial_paths.count(path)) continue;
            auto& catalog = catalogs[path];
            catalog.path = path;
            catalog.aliases.push_back(entry.first);
            const auto builtin = defaults.find(entry.first);
            if (builtin != defaults.end() && builtin->second == entry.second &&
                builtin->second.find("_bitlist.") != std::string::npos) {
                catalog.builtin_bitlist = true;
                catalog.declared_lengths.insert(
                    entry.first == "splitseq_bc1" || entry.first == "splitseq_bc2" ? 8 : 16);
            }
        } catch (const std::filesystem::filesystem_error& error) {
            throw std::runtime_error("Cannot access registered whitelist " + entry.first +
                                     ": " + error.what());
        } catch (const std::exception& error) {
            // The resolver explicitly diagnoses missing resources. Other I/O
            // failures are errors, rather than an apparent zero-match result.
            const std::string message = error.what();
            if (message.find("RRAD_RESOURCE_MISSING:") != std::string::npos) throw;
            if (message.find("not found") == std::string::npos) throw;
            if (options.unavailable) options.unavailable->push_back(entry.first);
            RAD_CERR << "[discover] Skipping missing whitelist " << entry.first
                      << ": " << message << '\n';
        }
    }
    for (auto& entry : catalogs)
        std::sort(entry.second.aliases.begin(), entry.second.aliases.end());
    return catalogs;
}

// Custom integer catalogs can obtain their length from registered layouts.
// A range or conflicting declarations cannot safely decode a bitlist.
static void collect_layout_lengths(std::map<std::string, Catalog>& catalogs) {
    std::set<std::string> seen;
    for (const auto& entry : config_utils::list_read_layouts()) {
        try {
            const std::string path = config_utils::get_read_layout(entry.first);
            if (!seen.insert(path).second) continue;
            std::ifstream input(path);
            std::string first;
            if (!std::getline(input, first))
                throw std::runtime_error("Cannot read layout " + path);
            if (first.compare(0, 3, "\xEF\xBB\xBF") == 0) first.erase(0, 3);
            csv::CSVFormat format;
            format.header_row(first.rfind("Read Layout", 0) == 0 ? 1 : 0);
            csv::CSVReader reader(path, format);
            for (auto& row : reader) {
                if (row["class"].get<std::string>() != "barcode") continue;
                const std::string text = seq_utils::trim(row["expected_length"].get<std::string>());
                if (text.empty() || !std::all_of(text.begin(), text.end(),
                                               [](unsigned char c) { return std::isdigit(c); }))
                    continue;
                const unsigned long length = std::stoul(text);
                if (length == 0 || length > 64) continue;
                for (const auto& spec : whitelist_utils::parse_whitelist_specs(
                         row["whitelist"].get<std::string>())) {
                    if (spec.empty()) continue;
                    try {
                        const std::string catalog_path = whitelist_utils::kit_to_path(spec);
                        auto found = catalogs.find(catalog_path);
                        if (found != catalogs.end())
                            found->second.declared_lengths.insert(static_cast<int>(length));
                    } catch (const std::exception&) {
                        // Unrelated unavailable layout references do not
                        // invalidate the registered catalogs being compared.
                    }
                }
            }
        } catch (const std::exception& error) {
            RAD_CERR << "[discover] Cannot read barcode-length metadata for layout "
                      << entry.first << ": " << error.what() << '\n';
        }
    }
}

struct Scan {
    Recipe recipe;
    std::vector<std::string> ids;
    barcode_count_result observations;
};

} // namespace discovery_whitelist_detail

std::vector<CatalogMatch> identify_whitelists(const std::string& sample_path,
                                            const std::vector<Recipe>& recipes,
                                            const Options& options) {
    using namespace discovery_whitelist_detail;
    using Geometry = std::tuple<std::string, int, int, int>;
    std::map<Geometry, Scan> scans;
    for (const auto& recipe : recipes) {
        if (recipe.pattern.empty() || recipe.length < 1 || recipe.length > 64)
            throw std::invalid_argument("Whitelist discovery needs a pattern and barcode length 1..64");
        const auto captured = scan_masked_capture_length(recipe.pattern);
        if (captured > 0 && captured != static_cast<std::size_t>(recipe.length))
            throw std::invalid_argument("Masked pattern length differs from barcode length for " + recipe.id);
        const Geometry key{recipe.pattern, recipe.length, recipe.left_margin, recipe.right_margin};
        auto inserted = scans.try_emplace(key);
        if (inserted.second) inserted.first->second.recipe = recipe;
        auto& ids = inserted.first->second.ids;
        if (std::find(ids.begin(), ids.end(), recipe.id) == ids.end()) ids.push_back(recipe.id);
    }
    if (scans.empty()) return {};

    auto catalogs = registered_catalogs(options);
    collect_layout_lengths(catalogs);
    for (auto& entry : scans) {
        auto& scan = entry.second;
        // The materialized sample is already bounded by the discover command.
        // Never apply a candidate catalog here: raw evidence must be identical
        // for every catalog, independently of later cell-calling thresholds.
        scan.observations = process_fastq(
            sample_path, scan.recipe.pattern, scan.recipe.length,
            scan.recipe.left_margin, scan.recipe.right_margin,
            static_cast<int>(std::min(options.sample_reads,
                                     static_cast<std::size_t>(std::numeric_limits<int>::max()))),
            options.error_ratio, 10000, options.threads, nullptr,
            options.boundary_poll, scan.recipe.pattern.find_first_of("Nn") != std::string::npos);
        if (!scan.observations.succeeded)
            throw std::runtime_error("Whitelist discovery could not scan sample for " + scan.recipe.id);
    }

    std::vector<CatalogMatch> matches;
    for (const auto& entry : catalogs) {
        if (options.boundary_poll) options.boundary_poll();
        const auto& catalog = entry.second;
        if (catalog.builtin_bitlist && catalog.declared_lengths.size() == 1) {
            const int length = *catalog.declared_lengths.begin();
            if (std::none_of(scans.begin(), scans.end(), [&](const auto& scan) {
                    return scan.second.recipe.length == length;
                })) continue;
        }
        std::optional<LoadedCatalog> loaded;
        try {
            loaded = load_catalog(catalog, options);
        } catch (const UnknownBitlistLength& error) {
            // An explicit comparison must not turn an undecodable catalog
            // into an apparent biological no-hit result.
            if (!options.catalogs.empty()) throw;
            if (options.unavailable) {
                if (catalog.aliases.empty()) options.unavailable->push_back(catalog.path);
                else options.unavailable->insert(options.unavailable->end(),
                    catalog.aliases.begin(), catalog.aliases.end());
            }
            RAD_CERR << "[discover] Skipping whitelist: " << error.what() << '\n';
            continue;
        }
        // Only this catalog is resident. Exact lookup and catalog_size both
        // deduplicate repeated reference rows using scan-wl's compact sets.
        for (const auto& scan_entry : scans) {
            const auto& scan = scan_entry.second;
            if (!loaded->lengths.count(scan.recipe.length)) continue;
            CatalogMatch match;
            match.path = options.catalogs.empty() && !catalog.aliases.empty()
                ? catalog.aliases.front() : catalog.path;
            match.aliases = catalog.aliases;
            match.length = scan.recipe.length;
            match.sampled_reads = scan.observations.reads_processed;
            match.extracted = scan.observations.total_extractions;
            match.observed_distinct = scan.observations.counts.size();
            match.catalog_size = loaded->members.size();
            for (const auto& observed : scan.observations.counts) {
                packed_scan_barcode key = 0;
                bool found = false;
                if (pack_scan_barcode(observed.first, key)) {
                    found = loaded->members.packed_barcodes.find(key) !=
                            loaded->members.packed_barcodes.end();
                } else if (observed.first.size() > 31 && observed.first.size() <= 64) {
                    const int64_seq encoded(observed.first);
                    found = loaded->members.long_barcodes.contains(encoded);
                }
                if (found) {
                    match.reads += observed.second;
                    ++match.distinct;
                }
            }
            for (const auto& id : scan.ids) {
                match.recipe_id = id;
                matches.push_back(match);
            }
        }
    }
    std::sort(matches.begin(), matches.end(), [](const CatalogMatch& a, const CatalogMatch& b) {
        if (a.recipe_id != b.recipe_id) return a.recipe_id < b.recipe_id;
        if (a.reads != b.reads) return a.reads > b.reads;
        return a.path < b.path;
    });
    return matches;
}

} // namespace rad_discovery
