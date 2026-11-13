#pragma once
// Spatial identity is supplied explicitly, never inferred from whitelist order.

namespace rad_spatial {
class Gzip {
    gzFile file_ = nullptr;
    std::function<void()> poll_;
public:
    Gzip(const std::string& path, const char* mode, std::function<void()> poll = {}) : poll_(std::move(poll)) {
        file_ = gzopen(path.c_str(), mode);
        if (!file_) throw std::runtime_error("Cannot open " + path);
        gzbuffer(file_, 1 << 20);
    }
    ~Gzip() { if (file_) gzclose(file_); }
    Gzip(const Gzip&) = delete;
    Gzip& operator=(const Gzip&) = delete;
    gzFile get() const { return file_; }
    void check() const {
        int error = Z_OK;
        const char* message = gzerror(file_, &error);
        if (error != Z_OK && error != Z_STREAM_END)
            throw std::runtime_error(std::string("gzip stream: ") + message);
    }
    void write(const std::string& text) {
        size_t offset = 0;
        while (offset < text.size()) {
            unsigned size = static_cast<unsigned>(std::min<size_t>(text.size() - offset, 1 << 20));
            if (gzwrite(file_, text.data() + offset, size) != static_cast<int>(size))
                throw std::runtime_error("Failed writing spatial output");
            offset += size;
        }
    }
    bool line(std::string& text) {
        text.clear();
        char buffer[65536];
        while (gzgets(file_, buffer, sizeof(buffer))) {
            if (poll_) poll_();
            text += buffer;
            if (!text.empty() && text.back() == '\n') break;
        }
        check();
        if (text.empty()) return false;
        if (text.back() == '\n') text.pop_back();
        if (!text.empty() && text.back() == '\r') text.pop_back();
        return true;
    }
    void close() {
        gzFile closing = file_;
        file_ = nullptr;
        if (closing && gzclose(closing) != Z_OK)
            throw std::runtime_error("Failed closing spatial output");
    }
};

std::vector<std::string> split(const std::string& text, char separator) {
    std::vector<std::string> result;
    size_t start = 0;
    while (true) {
        size_t end = text.find(separator, start);
        result.push_back(text.substr(start, end - start));
        if (end == std::string::npos) return result;
        start = end + 1;
    }
}

bool token(const std::string& text) {
    return !text.empty() && std::all_of(text.begin(), text.end(), [](unsigned char c) {
        return c > 32 && c < 127;
    });
}

int integer(const std::string& text) {
    if (text.empty() || !std::all_of(text.begin(), text.end(), [](char c) { return c >= '0' && c <= '9'; }))
        throw std::runtime_error("Expected a nonnegative integer, got: " + text);
    size_t used = 0;
    long long value = std::stoll(text, &used);
    if (used != text.size() || value > std::numeric_limits<int>::max())
        throw std::runtime_error("Coordinate or bin size out of range: " + text);
    return static_cast<int>(value);
}

class Table {
    Gzip input_;
    std::unordered_map<std::string, size_t> columns_;
    std::vector<std::string> values_;
    size_t line_ = 1;
public:
    explicit Table(const std::string& path, std::function<void()> poll = {}) : input_(path, "rb", std::move(poll)) {
        std::string text;
        if (!input_.line(text)) throw std::runtime_error("Empty table: " + path);
        auto header = split(text, '\t');
        for (size_t i = 0; i < header.size(); ++i)
            if (header[i].empty() || !columns_.emplace(header[i], i).second)
                throw std::runtime_error("Empty or duplicate table column: " + path);
    }
    void require(const std::string& column) const {
        if (!columns_.count(column)) throw std::runtime_error("Required TSV column missing: " + column);
    }
    bool next() {
        std::string text;
        if (!input_.line(text)) return false;
        ++line_;
        values_ = split(text, '\t');
        if (values_.size() != columns_.size())
            throw std::runtime_error("Wrong TSV field count at line " + std::to_string(line_));
        return true;
    }
    std::string get(const std::string& column) const {
        auto found = columns_.find(column);
        return found == columns_.end() ? "" : values_.at(found->second);
    }
};

struct Position {
    int row = 0, col = 0;
    std::string in_tissue, cell, pixel_row, pixel_col;
};

std::string spot_id(int row, int col, int bin_size = 2) {
    std::ostringstream out;
    out << "s_" << std::setfill('0') << std::setw(3) << bin_size << "um_"
        << std::setw(5) << row << '_' << std::setw(5) << col << "-1";
    return out.str();
}

Position position(const Table& table) {
    Position value;
    value.row = integer(table.get("array_row"));
    value.col = integer(table.get("array_col"));
    value.in_tissue = table.get("in_tissue");
    value.cell = table.get("cell_id");
    value.pixel_row = table.get("pxl_row_in_fullres");
    value.pixel_col = table.get("pxl_col_in_fullres");
    if (!value.in_tissue.empty() && value.in_tissue != "0" && value.in_tissue != "1")
        throw std::runtime_error("in_tissue must be 0, 1, or empty");
    if (!value.cell.empty() && !token(value.cell)) throw std::runtime_error("Invalid cell_id");
    for (const auto& pixel : {value.pixel_row, value.pixel_col}) {
        if (pixel.empty()) continue;
        size_t used = 0;
        double number = std::stod(pixel, &used);
        if (used != pixel.size() || !std::isfinite(number))
            throw std::runtime_error("Invalid pixel coordinate: " + pixel);
    }
    return value;
}

using Axis = std::unordered_map<std::string, int>;
Axis load_axis(const std::string& path, const std::string& dimension, const std::function<void()>& poll) {
    Table table(path, poll);
    table.require("barcode"); table.require(dimension);
    Axis axis;
    std::set<int> positions;
    while (table.next()) {
        std::string barcode = table.get("barcode");
        if (barcode.empty() || barcode.find_first_not_of("ACGT") != std::string::npos)
            throw std::runtime_error("Axis barcodes must contain A/C/G/T only");
        int index = integer(table.get(dimension));
        if (!axis.emplace(barcode, index).second || !positions.insert(index).second)
            throw std::runtime_error("Duplicate barcode or coordinate in axis map: " + path);
    }
    if (axis.empty()) throw std::runtime_error("Empty axis map: " + path);
    return axis;
}

struct Reference {
    Axis bc1, bc2;
    std::set<size_t> lengths;
    std::unordered_map<std::string, Position> direct, metadata;
    bool use_metadata = false;

    std::vector<Position> resolve(const std::string& barcode) const {
        if (!direct.empty()) {
            auto found = direct.find(barcode);
            return found == direct.end() ? std::vector<Position>{} : std::vector<Position>{found->second};
        }
        std::set<std::pair<int, int>> candidates;
        auto add = [&](const std::string& first, const std::string& second) {
            auto col = bc1.find(first), row = bc2.find(second);
            if (col != bc1.end() && row != bc2.end()) candidates.emplace(row->second, col->second);
        };
        auto parts = split(barcode, '-');
        if (parts.size() == 2) {
            add(parts[0], parts[1]); add(parts[1], parts[0]);
        } else if (parts.size() == 1) {
            for (size_t length : lengths) {
                if (length >= barcode.size()) continue;
                add(barcode.substr(0, length), barcode.substr(length));
                add(barcode.substr(barcode.size() - length), barcode.substr(0, barcode.size() - length));
            }
        }
        std::vector<Position> result;
        for (auto [row, col] : candidates) { Position p; p.row = row; p.col = col; result.push_back(p); }
        return result;
    }
};

void load_positions(const std::string& path, std::unordered_map<std::string, Position>& values, bool metadata, const std::function<void()>& poll) {
    Table table(path, poll);
    table.require("barcode"); table.require("array_row"); table.require("array_col");
    while (table.next()) {
        std::string key = table.get("barcode");
        Position p = position(table);
        if (!token(key)) throw std::runtime_error("Invalid mapping barcode");
        if (metadata && key != spot_id(p.row, p.col))
            throw std::runtime_error("Metadata must use native 2 um row/column spot IDs: " + key);
        if (!values.emplace(key, std::move(p)).second)
            throw std::runtime_error("Duplicate mapping barcode: " + key);
    }
    if (values.empty()) throw std::runtime_error("Empty position map: " + path);
}

struct Tags {
    std::string cb, umi, retained;
    bool invalid = false;
};
Tags tags(const std::string& comment) {
    Tags result;
    std::istringstream input(comment);
    std::set<std::string> seen;
    std::string item;
    while (input >> item) {
        if (item.rfind("CB:", 0) == 0 || item.rfind("UB:", 0) == 0) {
            std::string key = item.substr(0, 2);
            if (!seen.insert(key).second || item.substr(2, 3) != ":Z:") result.invalid = true;
            else if (key == "CB") result.cb = item.substr(5);
            else result.umi = item.substr(5);
        } else {
            // These tags belong to this export; refuse accidental reprocessing.
            if (item.rfind("XB:", 0) == 0 || item.rfind("XS:", 0) == 0 || item.rfind("XP:", 0) == 0)
                result.invalid = true;
            if (!result.retained.empty()) result.retained += '\t';
            result.retained += item;
        }
    }
    return result;
}
} // namespace rad_spatial

spatial_reformat_result reformat_spatial(const spatial_reformat_options& options) {
    using namespace rad_spatial;
        const auto& input = options.input;
        const auto& output = options.output;
        const auto& sample = options.sample;
        const auto& map_path = options.map_path;
        const auto& bc1_path = options.bc1_path;
        const auto& bc2_path = options.bc2_path;
        const auto& metadata_path = options.metadata_path;
        const auto& unit = options.unit;
        const int bin_size = options.bin_size;
        const bool tissue_only = options.tissue_only;
        if (input.empty() || output.empty() || !token(sample))
            throw std::runtime_error("--fastq, --outdir, and a whitespace-free --sample are required");
        if (sample.find(':') != std::string::npos) throw std::runtime_error("Sample ID cannot contain ':'");
        if (bin_size < 2 || bin_size % 2) throw std::runtime_error("--bin-size must be a positive multiple of 2 um");
        if (unit != "bin" && unit != "cell") throw std::runtime_error("--unit must be bin or cell");
        if (unit == "cell" && bin_size != 2)
            throw std::runtime_error("Cell assignment uses native 2 um spots; omit --bin-size with --unit cell");
        const bool axes = !bc1_path.empty() && !bc2_path.empty();
        if ((map_path.empty() && !axes) || (!map_path.empty() && (!bc1_path.empty() || !bc2_path.empty())))
            throw std::runtime_error("Supply --map OR both --bc1-axis and --bc2-axis");

        Reference reference;
        if (!map_path.empty()) load_positions(map_path, reference.direct, false, options.boundary_poll);
        else {
            reference.bc1 = load_axis(bc1_path, "array_col", options.boundary_poll);
            reference.bc2 = load_axis(bc2_path, "array_row", options.boundary_poll);
            for (const auto& entry : reference.bc1) reference.lengths.insert(entry.first.size());
        }
        if (!metadata_path.empty()) {
            reference.use_metadata = true;
            load_positions(metadata_path, reference.metadata, true, options.boundary_poll);
        }
        if (map_path.empty() && metadata_path.empty() && (tissue_only || unit == "cell"))
            throw std::runtime_error("Tissue filtering or cell assignment requires metadata");

        Gzip source(input, "rb", options.boundary_poll);
        if (!std::filesystem::create_directory(output))
            throw std::runtime_error("Output directory already exists; choose a new directory: " + output);
        const auto directory = std::filesystem::path(output);
        Gzip reads((directory / "reads.fastq.gz").string(), "wb1");
        Gzip assignments((directory / "assignments.tsv.gz").string(), "wb1");
        assignments.write("record_index\tread_id\tsample_id\tcorrected_barcode\tumi\tstatus\tspot_id\tarray_row\tarray_col\tbin_id\tbin_row\tbin_col\tcell_id\tunit_id\tsample_unit_id\tin_tissue\tpxl_row_in_fullres\tpxl_col_in_fullres\n");
        const auto start = std::chrono::steady_clock::now();
        std::map<std::string, size_t> counts;
        size_t total = 0;
        size_t emitted = 0;
        std::string header, sequence, plus, quality;
        while (source.line(header)) {
            ++total;
            if (header.empty() || header[0] != '@' || !source.line(sequence) ||
                !source.line(plus) || !source.line(quality) || plus.empty() || plus[0] != '+' ||
                sequence.size() != quality.size() ||
                sequence.find_first_of(" \t") != std::string::npos ||
                !std::all_of(quality.begin(), quality.end(), [](unsigned char c) { return c >= 33 && c <= 126; }))
                throw std::runtime_error("Expected complete four-line FASTQ at input record " + std::to_string(total));
            size_t comment_start = header.find_first_of(" \t");
            const std::string name = header.substr(1, comment_start == std::string::npos ? comment_start : comment_start - 1);
            if (!token(name)) throw std::runtime_error("Invalid FASTQ read identifier");
            Tags tag = tags(comment_start == std::string::npos ? "" : header.substr(comment_start + 1));
            std::string status = "assigned", native, bin, location, qualified;
            Position p;
            int bin_row = 0, bin_col = 0;
            bool located = false;
            if (tag.invalid) status = "invalid_tags";
            else if (tag.cb.empty()) status = "missing_barcode";
            else if (tag.umi.empty()) status = "missing_umi";
            else {
                auto candidates = reference.resolve(tag.cb);
                if (candidates.empty()) status = "unknown_barcode";
                else if (candidates.size() != 1) status = "ambiguous_barcode";
                else {
                    located = true;
                    p = candidates.front();
                    native = spot_id(p.row, p.col);
                    if (reference.use_metadata) {
                        auto found = reference.metadata.find(native);
                        if (found == reference.metadata.end()) status = "missing_metadata";
                        else p = found->second;
                    }
                    bin_row = p.row / (bin_size / 2); bin_col = p.col / (bin_size / 2);
                    bin = spot_id(bin_row, bin_col, bin_size);
                    if (status == "assigned" && tissue_only) {
                        if (p.in_tissue.empty()) status = "missing_tissue_status";
                        else if (p.in_tissue != "1") status = "off_tissue";
                    }
                    if (status == "assigned" && unit == "cell" && p.cell.empty()) status = "unassigned_cell";
                    if (status == "assigned") {
                        location = unit == "cell" ? p.cell : bin;
                        qualified = sample + ":" + location;
                    }
                }
            }
            ++counts[status];
            std::ostringstream row;
            row << total << '\t' << name << '\t' << sample << '\t' << tag.cb << '\t' << tag.umi << '\t' << status
                << '\t' << native << '\t' << (located ? std::to_string(p.row) : "")
                << '\t' << (located ? std::to_string(p.col) : "") << '\t' << bin
                << '\t' << (located ? std::to_string(bin_row) : "") << '\t' << (located ? std::to_string(bin_col) : "")
                << '\t' << p.cell << '\t' << location << '\t' << qualified << '\t' << p.in_tissue
                << '\t' << p.pixel_row << '\t' << p.pixel_col << '\n';
            assignments.write(row.str());
            if (status == "assigned") {
                std::string output_header = "@" + name + "\tCB:Z:" + location + "\tUB:Z:" + tag.umi
                    + "\tXB:Z:" + tag.cb + "\tXS:Z:" + sample + "\tXP:Z:" + native;
                if (!tag.retained.empty()) output_header += '\t' + tag.retained;
                reads.write(output_header + '\n' + sequence + '\n' + plus + '\n' + quality + '\n');
                ++emitted;
            } else if (!options.strict) {
                // Retained reads keep their original identity and tags. Their
                // audit row, not a fabricated coordinate, records the failure.
                reads.write(header + '\n' + sequence + '\n' + plus + '\n' + quality + '\n');
                ++emitted;
            }
        }
        source.check();
        reads.close(); assignments.close();
        std::ofstream summary(directory / "summary.tsv.partial");
        summary.exceptions(std::ios::failbit | std::ios::badbit);
        summary << "metric\tvalue\nstatus\tcomplete\ninput\t" << input << "\nsample_id\t" << sample
                << "\nmap\t" << map_path << "\nbc1_axis\t" << bc1_path << "\nbc2_axis\t" << bc2_path
                << "\nmetadata\t" << metadata_path << "\ncoordinate_convention\tzero_based_row_col_2um\nunit\t" << unit
                << "\nbin_size_um\t" << bin_size << "\nin_tissue_only\t" << tissue_only
                << "\nspatial_strict\t" << options.strict << "\ninput_reads\t" << total << "\nemitted_reads\t" << emitted
                << "\nassigned_reads\t" << counts["assigned"] << "\nelapsed_seconds\t"
                << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << '\n';
        for (const auto& entry : counts) summary << "status_" << entry.first << '\t' << entry.second << '\n';
        summary.close();
        std::filesystem::rename(directory / "summary.tsv.partial", directory / "summary.tsv");
        RAD_COUT << "[reformat " << (options.strict ? "--spatial-strict" : "--spatial") << "] "
                  << counts["assigned"] << '/' << total << " reads assigned; " << emitted << " emitted; audit: "
                  << (directory / "assignments.tsv.gz").string() << '\n';
        return {total, emitted, std::move(counts)};
}
