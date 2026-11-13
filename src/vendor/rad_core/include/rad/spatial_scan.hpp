#pragma once

// Included by whitelist_generator.hpp after the extraction/calling primitives.
// Spatial discovery reuses those primitives. Coordinates come from a slide
// resource or Visium HD's ordered component catalogues, never abundance sorting.

namespace rad_scan_spatial {
struct Options {
    bool enabled = false;
    bool slide_filter = false;
    std::string platform, slide, slide_map, slide_dir, bc1_axis, bc2_axis;
    int platform_flags = 0;
    bool masked_capture = false;
    int bc1_length = 0, bc2_length = 0;
    std::function<void()> boundary_poll;
};

// The callback is confined to the calling thread and restored on every exit.
// No R API is called from extraction workers.
inline thread_local const std::function<void()>* active_poll = nullptr;
inline void poll() { if (active_poll && *active_poll) (*active_poll)(); }
class InterruptScope {
    const std::function<void()>* previous_;
public:
    explicit InterruptScope(const std::function<void()>& callback) : previous_(active_poll) { if (callback) callback(); active_poll = &callback; }
    ~InterruptScope() { active_poll = previous_; }
};

class Text {
    gzFile file_ = nullptr;
    size_t operations_ = 0;
public:
    explicit Text(const std::string& path, const char* mode = "rb") {
        file_ = gzopen(path.c_str(), mode);
        if (!file_) throw std::runtime_error("Cannot open slide file: " + path);
        gzbuffer(file_, 1 << 20);
    }
    Text(const Text&) = delete;
    Text& operator=(const Text&) = delete;
    ~Text() { if (file_) gzclose(file_); }
    bool line(std::string& value) {
        if (++operations_ % 1024 == 0) poll();
        value.clear();
        char buffer[65536];
        while (gzgets(file_, buffer, sizeof(buffer))) {
            value += buffer;
            if (!value.empty() && value.back() == '\n') break;
        }
        int error = 0;
        const char* message = gzerror(file_, &error);
        if (error != Z_OK && error != Z_STREAM_END)
            throw std::runtime_error(std::string("Slide stream: ") + message);
        if (value.empty()) return false;
        if (value.back() == '\n') value.pop_back();
        if (!value.empty() && value.back() == '\r') value.pop_back();
        return true;
    }
    void write(const std::string& value) {
        if (++operations_ % 1024 == 0) poll();
        if (gzwrite(file_, value.data(), static_cast<unsigned>(value.size())) !=
            static_cast<int>(value.size()))
            throw std::runtime_error("Failed writing coordinate export");
    }
    void close() {
        auto f = file_; file_ = nullptr;
        if (f && gzclose(f) != Z_OK) throw std::runtime_error("Failed closing coordinate stream");
    }
};

inline std::string dna(std::string value) {
    value = seq_utils::trim(value);
    // Space Ranger's per-library barcode suffix is not part of the sequence.
    if (value.size() > 2 && value.substr(value.size() - 2) == "-1")
        value.resize(value.size() - 2);
    for (char& c : value) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (value.empty() || value.find_first_not_of("ACGT") != std::string::npos)
        throw std::runtime_error("Expected a nucleotide barcode, got: " + value +
            ". Spot/bin IDs are not a sequence-to-coordinate map.");
    return value;
}

inline double number(const std::string& value) {
    size_t used = 0;
    double result = std::stod(value, &used);
    if (used != value.size() || !std::isfinite(result) || std::abs(result) > 1e12)
        throw std::runtime_error("Invalid spatial coordinate: " + value);
    return result;
}

class Table {
    Text input_;
    char separator_ = '\0';
    std::map<std::string, size_t> columns_;
    std::vector<std::string> row_;
    bool first_row_pending_ = false;
public:
    explicit Table(const std::string& path) : input_(path) {
        std::string line;
        while (input_.line(line)) {
            if (line.empty() || line[0] == '#') continue;
            row_ = scan_wl_rescan_split(line, separator_);
            // Curio's native *_BeadBarcodes.txt files have no header.
            if (row_.size() == 3 && !row_[0].empty() &&
                row_[0].find_first_not_of("ACGTacgt") == std::string::npos) {
                columns_ = {{"barcode", 0}, {"x", 1}, {"y", 2}};
                first_row_pending_ = true;
                return;
            }
            for (size_t i = 0; i < row_.size(); ++i) {
                auto key = scan_wl_rescan_lower(scan_wl_rescan_trim(row_[i]));
                if (i == 0 && key.size() >= 3 && key.substr(0, 3) == "\xef\xbb\xbf")
                    key.erase(0, 3);
                if (!columns_.emplace(key, i).second)
                    throw std::runtime_error("Duplicate slide column: " + key);
            }
            return;
        }
        throw std::runtime_error("Empty slide table: " + path);
    }
    bool has(const std::string& key) const { return columns_.count(key) != 0; }
    std::string get(const std::string& key) const {
        auto col = columns_.find(key);
        if (col == columns_.end() || col->second >= row_.size())
            throw std::runtime_error("Missing slide column/value: " + key);
        return scan_wl_rescan_trim(row_[col->second]);
    }
    bool next() {
        if (first_row_pending_) { first_row_pending_ = false; return true; }
        std::string line;
        while (input_.line(line)) {
            if (line.empty() || line[0] == '#') continue;
            row_ = scan_wl_rescan_split(line, separator_);
            if (row_.size() != columns_.size())
                throw std::runtime_error("Inconsistent slide table column count");
            return true;
        }
        return false;
    }
};

struct Point { double x = 0, y = 0; };
struct SpotId { int row = 0, col = 0; std::string canonical; };
inline std::optional<SpotId> spot_id(const std::string& text) {
    if (text.rfind("s_", 0) != 0) return std::nullopt;
    static const std::regex pattern("s_002um_([0-9]+)_([0-9]+)(-1)?");
    std::smatch match;
    if (!std::regex_match(text, match, pattern))
        throw std::runtime_error("Expected native HD ID s_002um_ROW_COL[-1]: " + text);
    long long row = std::stoll(match[1]), col = std::stoll(match[2]);
    if (row > 1000000000 || col > 1000000000)
        throw std::runtime_error("HD spot ID is outside the supported coordinate range");
    std::ostringstream key;
    key << "s_002um_" << std::setfill('0') << std::setw(5) << row << '_' << std::setw(5) << col;
    return SpotId{static_cast<int>(row), static_cast<int>(col), key.str()};
}
struct Slide {
    std::unordered_map<std::string, Point> points;
    std::unordered_map<std::string, int> axis1, axis2;
    std::unordered_map<int, std::string> columns, rows;
    std::unordered_map<std::string, std::string> spot_barcodes;
    split_whitelist bc1, bc2;
    std::string source, units;
    bool axes = false;
    bool spot_labels = false;
    std::string mapping = "explicit_coordinates";
    std::string canonical_barcode(const std::string& value) const {
        const auto id = spot_id(value);
        if (!id) return dna(value);
        auto explicit_id = spot_barcodes.find(id->canonical);
        if (explicit_id != spot_barcodes.end()) return explicit_id->second;
        if (!columns.empty() && !rows.empty()) {
            auto a = columns.find(id->col), b = rows.find(id->row);
            if (a == columns.end() || b == rows.end())
                throw std::runtime_error("HD spot ID lies outside the component catalogues: " + value);
            return a->second + b->second;
        }
        return id->canonical;
    }
    uint64_t size() const {
        return axes ? static_cast<uint64_t>(axis1.size()) * axis2.size() : points.size();
    }
    std::optional<Point> lookup(const std::string& barcode) const {
        if (!axes) {
            auto it = points.find(barcode);
            return it == points.end() ? std::nullopt : std::optional<Point>(it->second);
        }
        std::optional<Point> result;
        for (int length : bc1.lengths) {
            if (barcode.size() <= static_cast<size_t>(length)) continue;
            auto a = axis1.find(barcode.substr(0, length));
            auto b = axis2.find(barcode.substr(length));
            if (a == axis1.end() || b == axis2.end()) continue;
            Point p{2.0 * a->second, 2.0 * b->second};
            if (result && (result->x != p.x || result->y != p.y)) return std::nullopt;
            result = p;
        }
        return result;
    }
};

inline void add_component(split_whitelist& wl, const std::string& barcode) {
    wl.seqs.insert(barcode);
    if (std::find(wl.lengths.begin(), wl.lengths.end(), barcode.size()) == wl.lengths.end())
        wl.lengths.push_back(static_cast<int>(barcode.size()));
}

inline std::string resolve(const Options& options) {
    namespace fs = std::filesystem;
    if (!options.slide_map.empty()) {
        if (!fs::is_regular_file(options.slide_map))
            throw std::runtime_error("Slide map does not exist: " + options.slide_map);
        return options.slide_map;
    }
    if (fs::is_regular_file(options.slide)) return options.slide;
    if (options.slide.empty() || options.slide.find_first_not_of(
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.-") != std::string::npos ||
        options.slide == "." || options.slide == "..")
        throw std::runtime_error("--slide must be an existing map file or a tile ID");
    std::vector<fs::path> roots;
    if (!options.slide_dir.empty()) roots.emplace_back(options.slide_dir);
    else {
        if (const char* dir = std::getenv("RAD_SLIDE_DIR")) roots.emplace_back(dir);
        roots.emplace_back((path_utils::find_resource_root() / "spatial").string());
    }
    for (const auto& root : roots)
        for (const auto* extension : {".tsv.gz", ".csv.gz", ".tsv", ".csv"}) {
            auto candidate = root / options.platform / (options.slide + extension);
            if (fs::is_regular_file(candidate)) return candidate.string();
        }
    throw std::runtime_error("No registered coordinate map for " + options.platform +
        " slide " + options.slide + ". Supply --slide-map FILE or register " +
        options.platform + "/" + options.slide +
        ".tsv.gz under --slide-dir (or RAD_SLIDE_DIR). A bare whitelist cannot supply coordinates.");
}

inline Slide load_slide(const Options& options, const std::string& bc1_file,
                        const std::string& bc2_file) {
    Slide slide;
    const bool hd = options.platform == "visium-hd";
    if (!options.bc1_axis.empty() || !options.bc2_axis.empty()) {
        if (!hd || options.bc1_axis.empty() || options.bc2_axis.empty() || !options.slide_map.empty())
            throw std::runtime_error("Both --bc1-axis and --bc2-axis require --visium-hd and replace --slide-map");
        slide.axes = true;
        slide.units = "um";
        slide.source = options.bc1_axis + " + " + options.bc2_axis;
        auto load_axis = [](const std::string& path, const char* column,
                            std::unordered_map<std::string, int>& axis,
                            std::unordered_map<int, std::string>& inverse, split_whitelist& wl) {
            Table table(path);
            std::set<int> coordinates;
            while (table.next()) {
                const auto barcode = dna(table.get("barcode"));
                const double value = number(table.get(column));
                if (value < 0 || value > 1e9 || value != std::floor(value))
                    throw std::runtime_error("HD axes require nonnegative native 2 um grid indices");
                if (!axis.emplace(barcode, static_cast<int>(value)).second ||
                    !coordinates.insert(static_cast<int>(value)).second)
                    throw std::runtime_error("HD axis contains duplicate barcode/coordinate");
                add_component(wl, barcode);
                inverse.emplace(static_cast<int>(value), barcode);
            }
            if (axis.empty()) throw std::runtime_error("Empty HD axis: " + path);
            std::sort(wl.lengths.begin(), wl.lengths.end());
        };
        load_axis(options.bc1_axis, "array_col", slide.axis1, slide.columns, slide.bc1);
        load_axis(options.bc2_axis, "array_row", slide.axis2, slide.rows, slide.bc2);
        return slide;
    }
    if (hd) {
        if (bc1_file.empty() != bc2_file.empty())
            throw std::runtime_error("Supply both -1 and -2 component catalogues");
        std::string a, b;
        if (!bc1_file.empty()) {
            a = whitelist_utils::kit_to_path(bc1_file);
            b = whitelist_utils::kit_to_path(bc2_file);
        } else {
            try {
                a = whitelist_utils::kit_to_path("visium_hd_bc1");
                b = whitelist_utils::kit_to_path("visium_hd_bc2");
            } catch (const std::exception&) { a.clear(); b.clear(); }
        }
        if (!a.empty() && !b.empty()) {
            slide.bc1 = load_split_whitelist_csv(a, options.bc1_length > 0 ? std::optional<uint16_t>(options.bc1_length) : std::nullopt);
            slide.bc2 = load_split_whitelist_csv(b, options.bc2_length > 0 ? std::optional<uint16_t>(options.bc2_length) : std::nullopt);
            if (slide.bc1.seqs.empty() || slide.bc2.seqs.empty())
                throw std::runtime_error("HD component catalogue is empty");
            for (size_t i = 0; i < slide.bc1.ordered.size(); ++i)
                slide.columns.emplace(static_cast<int>(i), dna(slide.bc1.ordered[i]));
            for (size_t i = 0; i < slide.bc2.ordered.size(); ++i)
                slide.rows.emplace(static_cast<int>(i), dna(slide.bc2.ordered[i]));
        }
    }
    try { slide.source = resolve(options); }
    catch (const std::exception&) {
        // Visium HD's catalogue order is its native axis order. An ID can
        // therefore name the slide without requiring an expanded position map.
        if (!hd || slide.columns.empty() || slide.rows.empty() || !options.slide_map.empty() ||
            options.slide.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.-") != std::string::npos)
            throw;
        slide.axes = true;
        slide.units = "um";
        slide.mapping = "visium_hd_catalog_order_bc1_column_bc2_row";
        slide.source = (bc1_file.empty() ? "visium_hd_bc1" : bc1_file) + std::string(" + ") +
                       (bc2_file.empty() ? "visium_hd_bc2" : bc2_file);
        for (const auto& a : slide.columns) slide.axis1.emplace(a.second, a.first);
        for (const auto& b : slide.rows) slide.axis2.emplace(b.second, b.first);
        return slide;
    }
    Table table(slide.source);
    const bool components = table.has("bc1") && table.has("bc2");
    const bool grid = table.has("array_row") && table.has("array_col");
    const bool microns = table.has("x_um") && table.has("y_um");
    const bool exported_units = table.has("coordinate_units") && table.has("x") && table.has("y");
    const bool ids_only = hd && table.has("barcode") && !table.has("x") && !grid && !microns;
    if (hd && !grid && !microns && !ids_only && !exported_units)
        throw std::runtime_error("HD maps require native 2 um array_row/array_col or x_um/y_um; x/y units cannot be guessed");
    const std::string xcol = grid ? "array_col" : microns ? "x_um" : "x";
    const std::string ycol = grid ? "array_row" : microns ? "y_um" : "y";
    if ((!table.has("barcode") && !components) || (!ids_only && (!table.has(xcol) || !table.has(ycol))))
        throw std::runtime_error("Slide table needs barcode,x,y (or x_um,y_um / array_row,array_col); HD also accepts bc1,bc2");
    slide.units = (hd || microns) ? "um" : grid ? "array" : "source";
    while (table.next()) {
        if (exported_units) {
            const auto units = table.get("coordinate_units");
            if (hd && units != "um") throw std::runtime_error("HD exported coordinates must use um units");
            if (!slide.points.empty() && units != slide.units)
                throw std::runtime_error("Mixed coordinate units in slide map");
            slide.units = units;
        }
        std::string barcode;
        if (components) {
            auto a = dna(table.get("bc1")), b = dna(table.get("bc2"));
            barcode = a + b;
            add_component(slide.bc1, a); add_component(slide.bc2, b);
            if (table.has("barcode") && !spot_id(table.get("barcode")) && dna(table.get("barcode")) != barcode)
                throw std::runtime_error("barcode disagrees with bc1+bc2 in slide map");
        } else barcode = hd ? slide.canonical_barcode(table.get("barcode")) : dna(table.get("barcode"));
        const auto id = hd && table.has("barcode") ? spot_id(table.get("barcode")) : std::nullopt;
        if (ids_only && !id) throw std::runtime_error("Coordinate-free HD map rows require native spot IDs");
        Point p = ids_only ? Point{2.0 * id->col, 2.0 * id->row} :
                            Point{number(table.get(xcol)), number(table.get(ycol))};
        if (hd && (p.x < 0 || p.y < 0 || (grid &&
                (p.x != std::floor(p.x) || p.y != std::floor(p.y)))))
            throw std::runtime_error("Invalid HD coordinate/grid index");
        if (hd && grid) { p.x *= 2; p.y *= 2; }
        if (id && (p.x != 2.0 * id->col || p.y != 2.0 * id->row))
            throw std::runtime_error("HD spot ID disagrees with its supplied coordinates");
        if (id) {
            if (components) {
                auto pair = slide.spot_barcodes.emplace(id->canonical, barcode);
                if (!pair.second && pair.first->second != barcode)
                    throw std::runtime_error("HD spot ID maps to conflicting nucleotide barcodes");
                if (!slide.columns.empty() && !slide.rows.empty() &&
                    (slide.columns.find(id->col) == slide.columns.end() ||
                     slide.rows.find(id->row) == slide.rows.end() ||
                     slide.columns.at(id->col) + slide.rows.at(id->row) != barcode))
                    throw std::runtime_error("HD spot ID components disagree with catalogue order");
            }
            slide.mapping = slide.columns.empty() ? "visium_hd_spot_id" : "visium_hd_catalog_order_bc1_column_bc2_row";
            slide.spot_labels = !components && slide.columns.empty();
        }
        auto inserted = slide.points.emplace(barcode, p);
        if (!inserted.second && (inserted.first->second.x != p.x || inserted.first->second.y != p.y))
            throw std::runtime_error("Barcode maps to conflicting spatial positions: " + barcode);
    }
    if (slide.points.empty()) throw std::runtime_error("Slide table has no barcodes");
    std::sort(slide.bc1.lengths.begin(), slide.bc1.lengths.end());
    std::sort(slide.bc2.lengths.begin(), slide.bc2.lengths.end());
    return slide;
}

inline void export_slide(const Slide& slide, const std::string& prefix) {
    Text out(prefix + "_coordinates.tsv.gz", "wb");
    out.write("barcode\tx\ty\tcoordinate_units\n");
    auto emit = [&](const std::string& barcode, Point p) {
        std::ostringstream row;
        row << std::setprecision(17) << barcode << '\t' << p.x << '\t' << p.y
            << '\t' << slide.units << '\n';
        out.write(row.str());
    };
    if (slide.axes) {
        for (const auto& a : slide.axis1) for (const auto& b : slide.axis2) {
            const auto sequence = a.first + b.first;
            const auto point = slide.lookup(sequence);
            if (!point) throw std::runtime_error("Ambiguous concatenated barcode in HD axis product");
            emit(sequence, *point);
        }
    } else for (const auto& row : slide.points) emit(row.first, row.second);
    out.close();
}

using Region = std::pair<int64_t, int64_t>;
// Eight native 2 um spots per side: 64 distinct spots in each 16 um region.
// Occupancy counts sites, not reads or aggregated 8 um bins.
inline Region region(Point p) {
    return {static_cast<int64_t>(std::floor(p.x / 16)), static_cast<int64_t>(std::floor(p.y / 16))};
}
inline size_t subbin(Point p) {
    return static_cast<size_t>((static_cast<int64_t>(std::floor(p.y / 2)) % 8) * 8 +
                               static_cast<int64_t>(std::floor(p.x / 2)) % 8);
}

inline int run(const Options& options, const std::string& input, const std::string& prefix,
               std::string adapter, int length, bool length_explicit, bool rescan,
               int left, int right, int max_reads, double error, int chunk_size, int threads,
               int umi, int offset_min, int offset_max, const std::string& bc1_file,
               const std::string& bc2_file, bool af, bool hs, bool verbose,
               whitelist_scan_stats* result_stats = nullptr, uint64_t* input_rows = nullptr,
               int* effective_length = nullptr) {
        InterruptScope interrupt_scope(options.boundary_poll);
        if (options.platform_flags != 1)
            throw std::runtime_error("--spatial requires exactly one of --visium, --visium-hd, --curio-seeker, --curio-trekker, --stereo-seq");
        if (options.slide.empty()) throw std::runtime_error("--spatial requires --slide FILE_OR_ID");
        if (prefix.empty()) throw std::runtime_error("Spatial scans require --output-prefix");
        if (threads < 0 || max_reads < 0 || left < 0 || right < 0 ||
            !std::isfinite(error) || error < 0 || error > 1 || umi < 0 || offset_min > offset_max)
            throw std::runtime_error("Invalid spatial scan parameter");
        const bool hd = options.platform == "visium-hd";
        const bool regular_visium = options.platform == "visium" && !options.slide_filter;
        if (hd && !options.slide_filter && (af || hs)) throw std::runtime_error("HD uses 32-of-64 native 2 um spots per 16 um region; --af-bcs/--hs-bcs do not apply");
        if (!hd && (!bc1_file.empty() || !bc2_file.empty()))
            throw std::runtime_error("Spatial split-barcode options require --visium-hd");
        auto slide = load_slide(options, bc1_file, bc2_file);
        // Do not let an output alias the input/map (including symlinks).
        for (const auto* suffix : {".csv", ".txt", "_rad.csv", "_rad.txt", "_coordinates.tsv.gz",
                                   "_regions.tsv", "_scan_wl.log"}) {
            const auto output = std::filesystem::weakly_canonical(prefix + suffix);
            for (const auto& source : {input, slide.source, options.bc1_axis, options.bc2_axis,
                                      bc1_file, bc2_file})
                if (!source.empty() && std::filesystem::is_regular_file(source) &&
                    output == std::filesystem::weakly_canonical(source))
                    throw std::runtime_error("Spatial output would overwrite input: " + source);
        }
        RAD_COUT << "Spatial platform: " << options.platform << "\nSlide: " << options.slide
                  << "\nMap: " << slide.source << "\nMapped barcodes: " << slide.size() << '\n';
        if (!ensure_output_directory(prefix)) throw std::runtime_error("Could not create spatial output directory");
        if (input.empty()) {
            if (rescan) throw std::runtime_error("--rescan requires --input");
            export_slide(slide, prefix);
            RAD_COUT << "Slide coordinates: " << prefix << "_coordinates.tsv.gz\n";
            return 0;
        }
        const auto started = std::chrono::steady_clock::now();
        barcode_count_result found;
        if (rescan) {
            const auto retained = read_scan_wl_rescan_counts(input, options.boundary_poll);
            if (input_rows) *input_rows = retained.input_rows;
            found = retained.counts;
            barcode_count_map canonical;
            for (const auto& item : found.counts)
                canonical[hd ? slide.canonical_barcode(item.first) : dna(item.first)] += item.second;
            found.counts = std::move(canonical);
        }
        else {
            bool masked_capture = options.masked_capture;
            if (adapter.empty()) {
                if (options.platform == "curio-seeker") { adapter = "NNNNNNNNTCTTCAGCGTTCCCGAGANNNNNN"; masked_capture = true; }
                else if (options.platform == "visium" || hd) adapter = "CTACACGACGCTCTTCCGATCT";
                else throw std::runtime_error("Supply -p with the spatial tag/CID adapter or capture pattern for " + options.platform);
            }
            const size_t masked = masked_capture ? scan_masked_capture_length(adapter) : 0;
            if (masked && (hd || left || right || (length_explicit && length != static_cast<int>(masked))))
                throw std::runtime_error("Masked capture conflicts with HD, margins, or explicit barcode length");
            if (hd) {
                if (slide.spot_labels)
                    throw std::runtime_error("FASTQ scanning needs HD component catalogues to translate spot IDs to nucleotide sequences");
                if (left || right || length_explicit)
                    throw std::runtime_error("HD uses component lengths and UMI/offset options, not -n/-l/-r");
                if (bc1_file.empty() != bc2_file.empty())
                    throw std::runtime_error("Supply both -1 and -2 component catalogues");
                if (!bc1_file.empty()) {
                    slide.bc1 = load_split_whitelist_csv(whitelist_utils::kit_to_path(bc1_file), options.bc1_length > 0 ? std::optional<uint16_t>(options.bc1_length) : std::nullopt);
                    slide.bc2 = load_split_whitelist_csv(whitelist_utils::kit_to_path(bc2_file), options.bc2_length > 0 ? std::optional<uint16_t>(options.bc2_length) : std::nullopt);
                }
                if (slide.bc1.seqs.empty() || slide.bc2.seqs.empty())
                    throw std::runtime_error("HD FASTQ scanning needs bc1/bc2 map columns, axis maps, or -1/-2 catalogues");
                const auto pairs = process_fastq_split_barcode(input, adapter, umi,
                    slide.bc1.seqs, slide.bc1.lengths, slide.bc2.seqs, slide.bc2.lengths,
                    offset_min, offset_max, max_reads, error, chunk_size, threads,
                    options.boundary_poll, masked_capture);
                found.reads_processed = pairs.reads_processed;
                found.total_extractions = pairs.total_extractions;
                for (const auto& pair : pairs.pair_counts) {
                    auto barcode = pair.first;
                    barcode.erase(barcode.find('\t'), 1);
                    found.counts[barcode] += pair.second;
                }
            } else {
                std::set<size_t> lengths;
                scan_whitelist_filter filter;
                for (const auto& p : slide.points) {
                    lengths.insert(p.first.size());
                    packed_scan_barcode key;
                    if (pack_scan_barcode(p.first, key)) filter.packed_barcodes.insert(key);
                    else filter.long_barcodes.insert(int64_seq(p.first));
                }
                if (lengths.size() != 1) throw std::runtime_error("Single-barcode spatial maps must use one barcode length");
                const int map_length = static_cast<int>(*lengths.begin());
                if (map_length < 1 || map_length > 64)
                    throw std::runtime_error("Single-barcode slide maps require 1 to 64 bases per barcode");
                if ((length_explicit && length != map_length) || (masked && static_cast<int>(masked) != map_length))
                    throw std::runtime_error("Capture length disagrees with slide barcode length");
                length = map_length;
                if (effective_length) *effective_length = length;
                found = process_fastq(input, adapter, length, left, right, max_reads, error,
                                      chunk_size, threads, options.slide_filter ? nullptr : &filter, options.boundary_poll, masked_capture);
                if (!found.succeeded) throw std::runtime_error("Spatial scan could not read input");
            }
        }
        const uint64_t extracted = found.total_extractions;
        uint64_t mapped_observations = 0;
        std::map<Region, std::bitset<64>> occupancy;
        size_t inspected = 0;
        for (auto it = found.counts.begin(); it != found.counts.end();) {
            if (++inspected % 1024 == 0) poll();
            auto p = slide.lookup(it->first);
            if (!p) {
                if (options.slide_filter) ++it;
                else it = found.counts.erase(it);
                continue;
            }
            mapped_observations += it->second;
            if (hd) occupancy[region(*p)].set(subbin(*p));
            ++it;
        }
        const uint64_t unmapped = extracted - mapped_observations;
        // Ordinary Visium retains the regular caller's all-extraction
        // normalization and thresholds. Coordinates only annotate its calls.
        if (!options.slide_filter) {
            found.total_extractions = regular_visium ? extracted : mapped_observations;
            found.filtered_to_whitelist = true;
            found.reference_whitelist_size = slide.size();
        }
        if (!mapped_observations) throw std::runtime_error("No observed barcodes match this slide; check input, map and capture pattern");

        std::unordered_set<std::string> selected;
        whitelist_scan_stats stats;
        if (!hd || options.slide_filter) {
            if (!count_perfect_matches_with_stats(found, {}, prefix + "_rad.csv", prefix + "_rad.txt",
                    verbose, af ? scan_whitelist_selection::above_floor : scan_whitelist_selection::high_specificity,
                    &stats, options.boundary_poll, !regular_visium)) throw std::runtime_error("Spatial abundance calling failed");
            Text calls(prefix + "_rad.txt");
            std::string barcode;
            while (calls.line(barcode)) selected.insert(barcode);
        } else {
            for (const auto& item : found.counts)
                if (occupancy.at(region(*slide.lookup(item.first))).count() >= 32) selected.insert(item.first);
            stats.reads_processed = found.reads_processed;
            stats.total_extractions = found.total_extractions;
            stats.unique_sequences = found.counts.size();
            stats.total_perfect_matches = found.total_extractions;
            stats.unique_perfect_matches = found.counts.size();
            stats.match_rate_percent = 100;
            stats.final_barcodes = stats.selected_barcodes = selected.size();
            stats.threshold_rule = "hd_32_of_64_2um_spots_in_16um_region";
        }
        const size_t abundance_selected = selected.size();
        const size_t abundance_high_specificity = stats.final_barcodes;
        std::unordered_set<std::string> final_calls;
        // Slide occupancy refines the common abundance calls. The historical
        // --slide route retains its original occupancy-only HD behavior.
        if (options.slide_filter) {
            auto passes_slide = [&](const std::string& barcode) {
                const auto point = slide.lookup(barcode);
                return point && (!hd || occupancy.at(region(*point)).count() >= 32);
            };
            for (auto it = selected.begin(); it != selected.end();) {
                if (!passes_slide(*it))
                    it = selected.erase(it);
                else
                    ++it;
            }
            Table abundance(prefix + "_rad.csv");
            while (abundance.next()) {
                const auto barcode = abundance.get("barcode");
                if (abundance.get("final_barcode") == "TRUE" && passes_slide(barcode))
                    final_calls.insert(barcode);
            }
            stats.final_barcodes = final_calls.size();
            stats.selected_barcodes = selected.size();
        }
        std::vector<std::string> keys;
        for (const auto& item : found.counts)
            if (slide.lookup(item.first)) keys.push_back(item.first);
        std::sort(keys.begin(), keys.end());
        std::ofstream csv(prefix + ".csv"), text(prefix + ".txt");
        csv << "barcode,count,x,y,coordinate_units,selected,final_barcode,occupied_2um_spots,region_x,region_y\n";
        Text coordinates(prefix + "_coordinates.tsv.gz", "wb");
        coordinates.write("barcode\tx\ty\tcoordinate_units\n");
        for (const auto& barcode : keys) {
            if (++inspected % 1024 == 0) poll();
            auto p = *slide.lookup(barcode);
            bool pass = selected.count(barcode);
            const bool final = options.slide_filter ? final_calls.count(barcode) != 0 : pass;
            csv << std::setprecision(17) << barcode << ',' << found.counts.at(barcode) << ','
                << p.x << ',' << p.y << ',' << slide.units << ',' << (pass ? "TRUE," : "FALSE,")
                << (final ? "TRUE," : "FALSE,");
            if (hd) {
                auto r = region(p);
                csv << occupancy.at(r).count() << ',' << r.first << ',' << r.second;
            } else csv << ",,";
            csv << '\n';
            if (pass) {
                text << barcode << '\n';
                std::ostringstream row;
                row << std::setprecision(17) << barcode << '\t' << p.x << '\t' << p.y << '\t' << slide.units << '\n';
                coordinates.write(row.str());
            }
        }
        coordinates.close(); csv.close(); text.close();
        if (!csv || !text) throw std::runtime_error("Failed writing spatial scan outputs");
        if (hd) {
            std::ofstream regions(prefix + "_regions.tsv");
            regions << "region_x\tregion_y\toccupied_2um_spots\ttotal_spots\tselected\n";
            for (const auto& r : occupancy) {
                if (++inspected % 1024 == 0) poll();
                regions << r.first.first << '\t' << r.first.second << '\t' << r.second.count()
                        << "\t64\t" << (r.second.count() >= 32 ? "TRUE" : "FALSE") << '\n';
            }
            regions.close();
            if (!regions) throw std::runtime_error("Failed writing occupancy output");
        }
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
        if (hd && !options.slide_filter) {
            std::ofstream summary(prefix + "_scan_wl.log");
            summary << "RAD spatial whitelist scan summary\ninput=" << input
                    << "\nselection=hd_occupancy\nreads_processed=" << found.reads_processed
                    << "\ntotal_extractions=" << found.total_extractions
                    << "\nunique_sequences=" << found.counts.size()
                    << "\nthreshold_rule=" << stats.threshold_rule
                    << "\nnative_spot_size_um=2\nregion_size_um=16"
                    << "\noccupied_2um_spots_required=32\nspots_per_region=64"
                    << "\nwall_time_seconds=" << elapsed << '\n';
            summary.close();
            if (!summary) throw std::runtime_error("Failed writing HD scan summary");
        } else if (!write_whitelist_scan_log(prefix + "_scan_wl.log", input, stats, elapsed)) throw std::runtime_error("Could not write spatial scan summary");
        std::ofstream log(prefix + "_scan_wl.log", std::ios::app);
        log << "spatial_platform=" << options.platform << "\nslide=" << options.slide
            << "\ncoordinate_source=" << slide.source << "\ncoordinate_units=" << slide.units
            << "\ncoordinate_mapping=" << slide.mapping
            << "\nbarcode_kind=" << (slide.spot_labels ? "hd_spot_id" : "nucleotide")
            << "\nslide_barcodes=" << slide.size() << "\nextracted_observations=" << extracted
            << "\nmapped_observations=" << mapped_observations
            << "\nunmapped_observations=" << unmapped << '\n';
        if (!options.slide_filter)
            log << "calling_profile=" << (hd ? "hd_occupancy" : regular_visium ? "rad_standard" : "rad_spatial_low_count") << '\n';
        else
            log << "abundance_selected_barcodes=" << abundance_selected
                << "\nabundance_high_specificity_barcodes=" << abundance_high_specificity
                << "\nspatial_selected_barcodes=" << selected.size()
                << "\nspatial_filter_rule=" << (hd ? "hd_32_of_64_2um_spots_in_16um_region" : "coordinate_membership") << '\n';
        log << "selected_barcodes=" << selected.size() << '\n';
        log.close();
        if (!log) throw std::runtime_error("Failed writing spatial provenance");
        RAD_COUT << "Spatial whitelist: " << selected.size() << " barcodes; " << prefix << ".txt\n";
        if (result_stats) *result_stats = stats;
        return 0;
}
} // namespace rad_scan_spatial
