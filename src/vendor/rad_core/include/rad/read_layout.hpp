#pragma once
#include "rad_headers.h"

// Adapter edit-distance limits and the misalignment calibration sample.
namespace adapter_thresholds {

constexpr std::size_t fallback_error_num = 3;
constexpr std::size_t fallback_error_denom = 10;

/**
 * Return the shared fallback edit-distance limit for an uncalibrated adapter.
 *
 * Keep this integer-only so layout preparation and read processing cannot
 * drift because of different rounding rules or duplicated constants.
 */
constexpr int fallback_max_edit_distance(std::size_t adapter_length) {
    const auto rounded_up =
        (adapter_length * fallback_error_num +
         fallback_error_denom - 1) /
        fallback_error_denom;
    return std::max(1, static_cast<int>(rounded_up));
}

static_assert(fallback_max_edit_distance(22) == 7,
              "22 nt adapters must allow seven fallback edits");
static_assert(fallback_max_edit_distance(23) == 7,
              "23 nt adapters must allow seven fallback edits");

// misalign_lower = largest edit distance with at most this share of the
// element's chance hits (the null) at or below it.
constexpr double misalign_lower_chance_share = 0.05;

// Chance hits an element needs for a calibrated threshold; below this,
// update_read_layout uses fallback_max_edit_distance.
constexpr std::size_t min_misal_obs = 100;

// Misalignment calibration sample: the first block (rad prep: -n), then
// further blocks until every element has min_misal_obs chance
// hits, up to cal_max_reads. Whole blocks keep the sample independent
// of the thread count.
constexpr std::size_t cal_block_reads = 50000;
constexpr std::size_t cal_max_reads = 1000000;

}  // namespace adapter_thresholds

/**
 * @brief Represents a position in the ReadLayout. Each string from position strings is processed in this manner.
 * @param ref_id Reference element ID this position is based on
 * @param is_start True if this is a start position, false for stop
 * @param offset Offset from the reference position
 * @param add_flags Additional flags (e.g., "terminal_linked", "skipped")
 */
struct ParsedPosition {
    std::string ref_id; // Reference element ID this position is based on
    bool is_start = false; // True if this is a start position, false for stop
    int offset = 0; // Offset from the reference position
    std::string add_flags; // Additional flags (e.g., "terminal_linked", "skipped")

    std::string to_string() const {
        if (ref_id.empty()) {
            return "";
        }

        std::stringstream ss;
        ss << ref_id << "|"
           << (is_start ? "start" : "stop")
           << (offset >= 0 ? "+" : "") << offset;
        if (!add_flags.empty()) {
            ss << "|" << add_flags;
        }
        return ss.str();
    }

    static ParsedPosition from_string(const std::string& pos_str, bool verbose) {
        ParsedPosition pos;
        if(pos_str.empty()){
            if(verbose) RAD_COUT << "[from_string] Empty position string\n";
            return pos;
        }

        // First split by pipes to see what we're dealing with
        std::vector<std::string> parts;
        std::string temp;
        std::istringstream pipe_stream(pos_str);
        if(verbose){
            RAD_COUT << "[from_string] Parsing position string: '" << pos_str << "'\n";
        }
        if (!pipe_stream.good()) {
            throw std::runtime_error("could not parse layout position string");
        }
        while(std::getline(pipe_stream, temp, '|')) {
            if(verbose){
                RAD_COUT << "[from_string] Flags: eof =" << pipe_stream.eof()  << " fail =" << pipe_stream.fail()  << " bad =" << pipe_stream.bad() << "\n";
            }
            parts.push_back(temp);
            if(verbose){
                RAD_COUT << "[from_string] Found part: '" << temp << "'\n";
            }
        }

        if(parts.empty() || parts[0].empty()) {
            if(verbose) {
                RAD_COUT << "[from_string] Empty reference ID; treating position as unset\n";
            }
            return pos;
        }
                
        if(parts.size() >= 1) {
            pos.ref_id = parts[0];
            if(verbose) RAD_COUT << "[from_string] Set ref_id: " << pos.ref_id << "\n";      
        }

        if(parts.size() >= 2) {
            // Handle the start/stop and offset part
            std::string pos_and_offset = parts[1];
            pos.is_start = (pos_and_offset.find("start") != std::string::npos);
            std::string is_start = pos.is_start ? "yes" : "no"; 
            if(verbose) RAD_COUT << "[from_string] Is start: " << is_start << "\n";

            // Find the offset
            size_t plus_pos = pos_and_offset.find('+');
            size_t minus_pos = pos_and_offset.find('-');
            
            if(plus_pos != std::string::npos) {
                std::string offset_str = pos_and_offset.substr(plus_pos + 1);
                pos.offset = std::stoi(offset_str);
                if(verbose) RAD_COUT << "[from_string] Found positive offset: " << pos.offset << "\n";
            } else if(minus_pos != std::string::npos) {
                std::string offset_str = pos_and_offset.substr(minus_pos + 1);
                pos.offset = -std::stoi(offset_str);
                if(verbose) RAD_COUT << "[from_string] Found negative offset: " << pos.offset << "\n";
            }
        }

        if(parts.size() >= 3) {
            pos.add_flags = parts[2];
            if(verbose) RAD_COUT << "[from_string] Set flags: " << pos.add_flags << "\n";
        }

        return pos;
    }

};

/**
 * @brief Groups primary and secondary positions for an element. These are the stored breakdowns of the string for easy access.
 */
struct ReferencePositions {
    ParsedPosition primary_start;      //< Primary starting position
    ParsedPosition primary_stop;       //< Primary stopping position
    ParsedPosition secondary_start;    //< Backup starting position
    ParsedPosition secondary_stop;     //< Backup stopping position
};

/**
 * @brief These are alignment positions for start and stop info in non-variable regions.
 * @param start_stats `std::pair<double,double>` Mean and variance of start positions
 * @param stop_stats `std::pair<double,double>` Mean and variance of stop positions
 */
struct AlignmentPositions {
    std::pair<double, double> start_stats;
    std::pair<double, double> stop_stats;
};

// Multi-index container tags
struct id_tag {};               // Tag for accessing by ID
struct length_tag {};           // Tag for accessing by length
struct type_tag {};             // Tag for accessing by type
struct order_tag {};            // Tag for accessing by order
struct direction_tag {};        // Tag for accessing by direction
struct global_class_tag {};     // Tag for accessing by global class
struct dir_order_tag {};        // Tag for accessing by order...and direction

/**
 * @brief Core element structure for read processing.
 * @brief This structure is used to store information about each element in the read layout.
 * @param class_id Unique identifier for the element class
 * @param seq Raw sequence data
 * @param masked_seq Masked version of the sequence (if static)
 * @param expected_length Expected sequence length if known
 * @param type Element type (e.g., "static", "variable")
 * @param order Position in the layout
 * @param direction Orientation ("forward" or "reverse")
 * @param unidirectional Whether the element is unidirectional ("forward_only" or "reverse_only")
 * @param global_class Global classification
 * @param whitelist_path Path to whitelist file (if applicable)
 * @param misalignment_threshold Threshold for misalignment detection
 * @param aligned_positions Alignment position storage
 * @param misaligned_positions Alignment position storage
 * @param ref_pos Position information
 * 
 */
struct ReadElement {
    std::string class_id;        // Unique identifier for the element class
    std::string seq;             // Raw sequence data
    std::string masked_seq;      // Masked version of the sequence (if static)
    std::optional<int> expected_length;  // Expected sequence length if known
    std::vector<int> length_candidates;  // Allowed lengths for variable elements (e.g. 15-16)
    std::string type;            // Element type (e.g., "static", "variable")
    int order;                   // Position in the layout
    std::string direction;       // Orientation ("forward" or "reverse")
    bool unidirectional;     // Whether the element is unidirectional
    std::string global_class;    // Global classification
    std::string whitelist_path;     //path to whitelist_file (if applicable)
    std::string flags;           // Optional layout flags (e.g. joint_barcode)
    std::optional<std::tuple<int, int, int>> misalignment_threshold;  // Threshold for misalignment detection
    std::optional<AlignmentPositions> aligned_positions;  // Alignment position storage
    std::optional<AlignmentPositions> misaligned_positions;  // Alignment position storage
    std::optional<ReferencePositions> ref_pos;  // Position information

    ReadElement(
        const std::string& class_id,
        const std::string& seq,
        const std::string& masked_seq,
        const std::optional<int>& expected_length,
        const std::string& type,
        int order,
        const std::string& direction,
        const std::string& global_class,
        const std::string& whitelist_path = "",
        const bool unidirectional = false,
        const std::optional<std::tuple<int, int, int>>& misalignment_threshold = std::nullopt,
        const std::optional<AlignmentPositions>& aligned_positions = std::nullopt,
        const std::optional<AlignmentPositions>& misaligned_positions = std::nullopt,
        const std::optional<ReferencePositions>& ref_pos = std::nullopt
    ) : class_id(std::move(class_id)),
        seq(std::move(seq)),
        masked_seq(std::move(masked_seq)),
        expected_length(expected_length),
        type(std::move(type)),
        order(order),
        direction(std::move(direction)),
        unidirectional(unidirectional),
        global_class(std::move(global_class)),
        whitelist_path(std::move(whitelist_path)),
        misalignment_threshold(misalignment_threshold),
        aligned_positions(aligned_positions),
        misaligned_positions(misaligned_positions),
        ref_pos(ref_pos) {}
};

/**
 * @brief Multi-index container for read layout elements
 * @param ReadElement Core element structure
 */
typedef boost::multi_index::multi_index_container<
    ReadElement,
    boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<
            boost::multi_index::tag<id_tag>,
            boost::multi_index::member<ReadElement, std::string, &ReadElement::class_id>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<length_tag>,
            boost::multi_index::member<ReadElement, std::optional<int>, &ReadElement::expected_length>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<type_tag>,
            boost::multi_index::member<ReadElement, std::string, &ReadElement::type>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<order_tag>,
            boost::multi_index::member<ReadElement, int, &ReadElement::order>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<direction_tag>,
            boost::multi_index::member<ReadElement, std::string, &ReadElement::direction>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<global_class_tag>,
            boost::multi_index::member<ReadElement, std::string, &ReadElement::global_class>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<dir_order_tag>,
            boost::multi_index::composite_key<
                ReadElement,
                boost::multi_index::member<ReadElement, std::string, &ReadElement::direction>,
                boost::multi_index::member<ReadElement, int, &ReadElement::order>
            >
        >
    >
> Layout_Struct;

/**
 * @brief Full read layout representation, containing whitelist map, layout elements, and position map
 * @param wl_map Whitelist map for barcode correction
 * @param layout Multi-index container for read layout elements
 * @param sequencing_type Sequencing type information (bulk, single-cell)
 * @param position_map Map of reference positions
 */
class ReadLayout {
public:
    whitelist wl_map;
    Layout_Struct layout;
    std::string sequencing_type;
    std::unordered_map<std::string, ReferencePositions> position_map;

    // --concat-hmm: calibrated hmm_* parameters and the shared read-only model.
    concat_hmm::Params concat_params;
    std::shared_ptr<const concat_hmm::Model> concat_model;
    // --concat-hmm-abstain=legacy: abstained reads take the default path instead of being dropped.
    bool concat_abstain_legacy = false;

    //A bunch of access-based methods for the multi-index container
    auto& by_id() { return layout.get<id_tag>(); }
    auto& by_length() { return layout.get<length_tag>(); }
    auto& by_type() { return layout.get<type_tag>(); }
    auto& by_order() { return layout.get<order_tag>(); }
    auto& by_direction() { return layout.get<direction_tag>(); }
    auto& by_global_class() { return layout.get<global_class_tag>(); }
    auto& by_dir_order() { return layout.get<dir_order_tag>(); }

    // Const versions for read-only access
    const auto& by_id() const { return layout.get<id_tag>(); }
    const auto& by_length() const { return layout.get<length_tag>(); }
    const auto& by_type() const { return layout.get<type_tag>(); }
    const auto& by_order() const { return layout.get<order_tag>(); }
    const auto& by_direction() const { return layout.get<direction_tag>(); }
    const auto& by_global_class() const { return layout.get<global_class_tag>(); }
    const auto& by_dir_order() const { return layout.get<dir_order_tag>(); }

    // ReadLayout -> concat_hmm::layout_spec.
    concat_hmm::layout_spec to_concat_layout_spec() const {
        concat_hmm::layout_spec spec;
        spec.name = "rad_read_layout";
        spec.mode = sequencing_type;
        spec.bulk = (sequencing_type == "bulk");
        for (const auto& e : by_order()) {
            concat_hmm::layout_element x;
            x.id = e.class_id;
            x.seq = e.seq;
            x.is_static = (e.type == "static");
            x.klass = e.global_class.empty() ? e.class_id : e.global_class;
            if (x.klass == "poly_t" || x.klass == "poly_a") x.klass = "poly_tail";
            x.direction = (e.direction.rfind("reverse", 0) == 0) ? 'R' : 'F';
            x.order = e.order;
            x.length_candidates = e.length_candidates;
            if (x.length_candidates.empty() && e.expected_length && *e.expected_length > 0)
                x.length_candidates = {*e.expected_length};
            if (e.misalignment_threshold) x.max_edits = std::get<0>(*e.misalignment_threshold);
            spec.elements.push_back(std::move(x));
        }
        return spec;
    }

    // On failure the model stays unset and every read takes the default path.
    bool build_concat_model(bool verbose) {
        try {
            const concat_hmm::layout_spec spec = to_concat_layout_spec();
            const concat_hmm::Params* p = concat_params.empty() ? nullptr : &concat_params;
            auto model = std::make_shared<concat_hmm::Model>(concat_hmm::Model::build(spec, p));
            RAD_COUT << "[concat_hmm] Model built: " << model->tpl.size() << " templates, "
                      << model->A << " anchors, params="
                      << (p ? "position-map hmm_* columns" : "layout priors")
                      << (!model->param_warnings.empty()
                              ? " (invalid hmm_* cells: every read uses the existing path; remove the hmm_* columns to recalibrate)"
                              : model->guard_failed ? " (calibration guard FAILED: every read uses the existing path)" : "")
                      << "\n";
            for (const auto& w : model->param_warnings) RAD_CERR << "[concat_hmm] WARNING: " << w << "\n";
            if (p && concat_params.topology != 0 && concat_params.topology != model->topo_hash)
                RAD_COUT << "[concat_hmm] WARNING: hmm_* columns were calibrated on a different layout topology; "
                             "structural tables fall back to priors\n";
            if (verbose) RAD_COUT << model->describe();
            concat_model = std::move(model);
            return true;
        } catch (const concat_hmm::layout_spec_error& ex) {
            RAD_CERR << "[concat_hmm] WARNING: cannot build the HMM for this layout (" << ex.what()
                      << "); --concat-hmm disabled, every read uses the existing path\n";
            concat_model.reset();
            return false;
        }
    }

    static std::optional<std::string> get_optional_csv_string(const csv::CSVRow& row, const std::string& field) {
        try {
            auto f = row[field];
            if (f.is_null()) return std::nullopt;
            auto v = seq_utils::trim(f.get<std::string>());
            if (v.empty()) return std::nullopt;
            return v;
        } catch (...) {
            return std::nullopt;
        }
    }

    static std::vector<int> parse_length_candidates(const std::optional<std::string>& spec_opt) {
        std::vector<int> out;
        if (!spec_opt.has_value()) return out;

        std::string spec = seq_utils::trim(spec_opt.value());
        if (spec.empty()) return out;

        auto push_if_valid = [&](int v) {
            if (v > 0) out.push_back(v);
        };

        try {
            auto dash = spec.find('-');
            if (dash != std::string::npos) {
                std::string lhs = seq_utils::trim(spec.substr(0, dash));
                std::string rhs = seq_utils::trim(spec.substr(dash + 1));
                int a = std::stoi(lhs);
                int b = std::stoi(rhs);
                if (a > b) std::swap(a, b);
                for (int v = a; v <= b; ++v) push_if_valid(v);
            } else {
                push_if_valid(std::stoi(spec));
            }
        } catch (...) {
            return {};
        }

        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        return out;
    }

    static std::string format_length_candidates(const std::vector<int>& lengths) {
        if (lengths.empty()) return "";
        if (lengths.size() == 1) return std::to_string(lengths.front());
        return std::to_string(lengths.front()) + "-" + std::to_string(lengths.back());
    }

    // Hand-written layouts carry a "Read Layout[:mode]" title row ahead of the column
    // header; caches written before the title was persisted do not, so the header row
    // has to be detected rather than assumed.
    static int detect_header_row(const std::string& input_file) {
        std::ifstream fin(input_file);
        std::string first;
        if (!std::getline(fin, first)) return 0;
        if (first.rfind("\xEF\xBB\xBF", 0) == 0) first.erase(0, 3);
        return first.rfind("Read Layout", 0) == 0 ? 1 : 0;
    }

    // A generated cache carries the internal schema: reverse elements already derived
    // and orders already assigned. Re-deriving those collides on the unique class_id
    // index, so the rows would be dropped rather than merged.
    static bool is_generated_cache(const std::string& input_file) {
        const int header = detect_header_row(input_file);
        std::ifstream fin(input_file);
        std::string line;
        for (int i = 0; i <= header; ++i) {
            if (!std::getline(fin, line)) return false;
        }
        return line.find("class_id") != std::string::npos
            && line.find("order") != std::string::npos;
    }

    // get read layout mode
    std::string get_rl_mode(std::string input_file) const {
        std::string mode;
        std::ifstream fin(input_file);
        std::string title;
        std::getline(fin, title);           // e.g. "Read Layout:R1"
        std::stringstream ss(title);

        std::string drop;
        std::getline(ss, drop, ':');        // drop “Read Layout”
        if (!std::getline(ss, mode, ':')) { // try to read the part after the first colon
            mode.clear();
        }

        // now trim whitespace and any double-quotes
        const char* trim_chars = " \r\n\t\",";
        auto start = mode.find_first_not_of(trim_chars);
        if (start == std::string::npos) {
            mode.clear();
        } else {
            auto end = mode.find_last_not_of(trim_chars);
            mode = mode.substr(start, end - start + 1);
        }
        return mode;
    }

    // Import data from CSV
    void prep_new_layout(const std::string& input_file, bool verbose) {

        std::unordered_map<std::string, int> class_id_counts;
        std::unordered_map<std::string, int> class_id_total_counts;
        bool build_forward_only, build_reverse_only = false;
        std::vector<std::pair<ReadElement, size_t>> deferred_reverse_only;
        std::unordered_set<std::string> single_sided_ids;

        auto log_elem = [&](const ReadElement &e){
            RAD_COUT << "[element] order=" << e.order
                      << " id=" << e.class_id
                      << " type=" << e.type
                      << " dir=" << e.direction
                      << " class=" << e.global_class
                      << " seq=\"" << e.seq << "\""
                      << " exp_len=" << (e.expected_length? std::to_string(*e.expected_length) : "none")
                      << " len_cands=" << (e.length_candidates.empty() ? "none" : format_length_candidates(e.length_candidates))
                      << " flags=" << (e.flags.empty() ? "none" : e.flags)
                      << "\n";
        };
       
        std::string mode = get_rl_mode(input_file);
        //std::string sequencing_type;

        std::regex_search(mode, std::regex("bulk")) ? sequencing_type = "bulk" : "";

        // set flags
        build_forward_only = std::regex_search(mode, std::regex("R1"));
        build_reverse_only = std::regex_search(mode, std::regex("R2"));

        if(verbose){
            RAD_COUT << "[prep_new_layout] Mode: '" << mode << "'\n";
            RAD_COUT << "[prep_new_layout] Build mode: " << (build_forward_only ? "forward"  : build_reverse_only  ? "reverse" : "both")<< "\n";    
        }
        
        if (is_generated_cache(input_file)) {
            throw std::runtime_error(
                "'" + input_file + "' is a generated layout; pass the original to prep.");
        }

        // Configure CSV reader
        csv::CSVFormat format;
        format.delimiter(',').header_row(detect_header_row(input_file));
        csv::CSVReader reader(input_file, format);

        const auto& headers = reader.get_col_names();
        bool has_direction = std::find(headers.begin(), headers.end(), "direction") != headers.end();

        // First pass: Collect total counts of class_ids
        std::vector<csv::CSVRow> rows;
        for (auto& row : reader) {
            rows.push_back(row);
            std::string class_id = row["class"].is_null() || row["class"].get<std::string>().empty()
                                    ? row["id"].get<std::string>() : row["class"].get<std::string>();
            class_id_total_counts[class_id]++;
        }

        // Initialize order counter and insert seq_start
        int order_counter = 1;
        ReadElement seq_start("seq_start", "", "", 0, "static", order_counter++, "forward", "start");
        //if(!build_reverse_only) 
        layout.insert(seq_start);
        size_t row_idx = 0;

        // Process forward elements
        for (auto& row : rows) {
            std::string id = row["id"].is_null() ? "" : row["id"].get<std::string>();
            std::string seq = row["seq"].is_null() ? "" : row["seq"].get<std::string>();
            std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
            std::string type = row["type"].is_null() ? "variable" : row["type"].get<std::string>();
            std::string global_class = row["class"].is_null() ? "" : row["class"].get<std::string>();
            std::string direction = (!has_direction || row["direction"].is_null()) ? "forward" : row["direction"].get<std::string>();
            std::transform(direction.begin(), direction.end(), direction.begin(), ::tolower);
            std::string class_id = row["class"].is_null() ? id : row["class"].get<std::string>();
            std::string whitelist_path = row["whitelist"].is_null() ? "" : row["whitelist"].get<std::string>();
            std::string flags = get_optional_csv_string(row, "flags").value_or("");

            bool is_forward_only = (direction == "forward_only");
            bool is_reverse_only = (direction == "reverse_only");
            bool unidirectional = is_forward_only || is_reverse_only;
            
            std::string normalized_direction = direction;

            if (is_forward_only){
                normalized_direction = "forward";
            }
            if (is_reverse_only) {
                normalized_direction = "reverse";
            }


            // Fill empty class, class_id, or expected length
            auto length_spec = get_optional_csv_string(row, "expected_length");
            std::vector<int> length_candidates = parse_length_candidates(length_spec);
            std::optional<int> expected_length = length_candidates.empty()
                ? std::nullopt
                : std::optional<int>(length_candidates.back());
            if (global_class.empty()) {
                global_class = class_id;
            }
            if (class_id.empty()) {
                class_id = global_class;
            }

            //populate poly-tail cases
            if (class_id == "poly_t") {
                expected_length = expected_length.value_or(12);
                seq = "T{" + std::to_string(*expected_length) + ",}+";
                global_class = "poly_tail";
                if (length_candidates.empty()) length_candidates.push_back(*expected_length);
            } else if (class_id == "poly_a") {
                expected_length = expected_length.value_or(12);
                seq = "A{" + std::to_string(*expected_length) + ",}+";
                global_class = "poly_tail";
                if (length_candidates.empty()) length_candidates.push_back(*expected_length);
            }

            if (class_id_total_counts[class_id] > 1) {
                int &cnt = class_id_counts[class_id];
                ++cnt;
                std::string original = class_id;
                class_id += "_" + std::to_string(cnt);
                if (verbose) {
                    RAD_COUT << "[prep_new_layout] Duplicate class " << original << "' found (total = " 
                              << class_id_total_counts[original] << "), renaming to " << class_id << "\n";
                }
            }

            if (!seq.empty() && !expected_length.has_value()) {
                expected_length = static_cast<int>(seq.length());
                if (length_candidates.empty()) length_candidates.push_back(*expected_length);
            }
            std::string masked_seq = (type == "static" && global_class != "poly_tail") ? "MASKED" : "";

           ReadElement elem(
                            class_id, 
                            seq, 
                            masked_seq, 
                            expected_length, 
                            type, 
                            order_counter, 
                            normalized_direction, 
                            global_class, 
                            whitelist_path, 
                            unidirectional
                        );
            elem.length_candidates = std::move(length_candidates);
            elem.flags = std::move(flags);

            if (is_forward_only) {
                single_sided_ids.insert(elem.class_id);
                ++order_counter;
                layout.insert(elem);
                if (verbose) {
                    RAD_COUT << "[read_layout] Inserted element: ";
                    log_elem(elem);
                }
            } else if (is_reverse_only) {
                single_sided_ids.insert(elem.class_id);
                deferred_reverse_only.push_back({elem, row_idx});
            } else {
                ++order_counter;
                layout.insert(elem);
                if (verbose) {
                    RAD_COUT << "[read_layout] Inserted element: ";
                    log_elem(elem);
                }
            }
            ++row_idx;
        }

        // Insert seq_stop
        ReadElement seq_stop(
            "seq_stop", 
            "", 
            "", 
            0, 
            "static", 
            order_counter++, 
            "forward", 
            "stop"
        );
        layout.insert(seq_stop);
        if(build_forward_only){
            RAD_COUT << "[prep_new_layout] Forward-only (R1) layout generated.\n";
            return;
        }
        // Generate reverse complement entries
        std::vector<ReadElement> reverse_elements;
        // Determine the starting point for reverse complement orders
        int reverse_order_start = order_counter;
        const auto& ordered_index = by_order();

        for (const auto& elem : ordered_index) {
            if (single_sided_ids.count(elem.class_id) || elem.direction != "forward") {
                if (verbose) {
                    RAD_COUT << "[prep_new_layout] Skipping reverse-complement generation for "
                              << elem.class_id << " (direction = " << elem.direction << ")\n";
                }
                continue;
            }
            ReadElement reverse_elem = elem;
            reverse_elem.direction = "reverse";
            if (elem.class_id == "poly_t") {
                reverse_elem.global_class = "poly_tail";
                reverse_elem.seq = "A{" + std::to_string(elem.expected_length.value_or(12)) + ",}+";
                reverse_elem.class_id = "poly_a";
            } else if (elem.class_id == "poly_a") {
                reverse_elem.global_class = "poly_tail";
                reverse_elem.seq = "T{" + std::to_string(elem.expected_length.value_or(12)) + ",}+";
                reverse_elem.class_id = "poly_t";
            } else {
                if(!elem.unidirectional){
                    reverse_elem.seq = seq_utils::revcomp(elem.seq);
                }
                reverse_elem.class_id = "rc_" + elem.class_id;
            }
            reverse_elements.push_back(reverse_elem);
        }
        // Organize reverse complement rows
        auto start_it = std::stable_partition(reverse_elements.begin(), reverse_elements.end(),
            [](const ReadElement& elem) { 
                return elem.global_class == "start"; 
            });

        auto stop_it = std::stable_partition(start_it, reverse_elements.end(),
            [](const ReadElement& elem) { 
                return elem.global_class != "stop"; 
            });

        // Reverse the middle portion
        std::reverse(start_it, stop_it);
        // Assign order to reverse complement elements and add them to layout
        for (auto& reverse_elem : reverse_elements) {
            reverse_elem.order = reverse_order_start++;
            layout.insert(reverse_elem);
            if (verbose) {
                RAD_COUT << "[prep_new_layout] Inserted element: ";
                log_elem(reverse_elem);
            }
        }

        order_counter = reverse_order_start;
        std::sort(deferred_reverse_only.begin(), deferred_reverse_only.end(),
                  [](const auto& a, const auto& b) { return a.second < b.second; });
        for (auto elem_pair : deferred_reverse_only) {
            auto elem = elem_pair.first;
            elem.order = order_counter++;
            layout.insert(elem);
            if (verbose) {
                RAD_COUT << "[read_layout] Inserted element: ";
                log_elem(elem);
            }
        }

        if (build_reverse_only) {
            // Get the index keyed on direction:
            auto& dir_index = layout.get<direction_tag>();
            // Find the sub-range where direction == "forward"
            auto range = dir_index.equal_range("forward");
            // Erase them
            dir_index.erase(range.first, range.second);
            RAD_COUT << "[prep_new_layout] Removed all forward elements; reverse-only layout ready. New layout :\n";
            for(auto& elem : layout) {
                log_elem(elem);
            }
            return;
        }
    }

    void discover_layout(const std::string& parts_csv, const std::string& fastq_path, const std::string& output_base, 
        size_t max_reads = 25000, size_t chunk_size = 2000, int num_threads = 1, bool verbose = false) {
        layout.clear();
        position_map.clear();

        struct adapter_stats {
            std::string id;
            std::string seq;
            std::string rc_seq;
            std::vector<double> forward_starts;
            std::vector<double> forward_stops;
            std::vector<double> reverse_starts;
            std::vector<double> reverse_stops;
        };

         struct ordered_element {
            std::string id;
            std::string class_id;
            std::string seq;
            std::string direction;
            double mean_start;
            size_t dominant_hits;
            size_t alternate_hits;
        };

        auto parts = [&]() {
            std::vector<adapter_stats> parsed;
            csv::CSVFormat fmt;
            fmt.delimiter(',').quote('"').header_row(0);
            csv::CSVReader reader(parts_csv, fmt);

            for (auto& row : reader) {
                if (row["id"].is_null() || row["seq"].is_null()) {
                    continue;
                }
                std::string id = row["id"].get<std::string>();
                std::string seq = row["seq"].get<std::string>();
                if (id.empty() || seq.empty()) {
                    continue;
                }
                std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
                if (verbose) {
                    RAD_COUT << "[generate_layout_from_parts] Part: " << id
                              << " Length:  " << seq.length() << " Seq: " << seq << "\n";
                }
                parsed.push_back({id, seq, seq_utils::revcomp(seq), {}, {}, {}, {}});
            }

            if (parsed.empty()) {
                throw std::runtime_error("No parts parsed from parts list CSV");
            }
            return parsed;
        }();

        auto record_alignment = [](const std::string& query,
                                   const std::string& read_seq,
                                   std::vector<double>& starts,
                                   std::vector<double>& stops) {
            if (read_seq.empty()) return;

            static const EdlibEqualityPair kWildcardEqualities[] = {
                {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'},
                {'A', 'N'}, {'C', 'N'}, {'G', 'N'}, {'T', 'N'}
            };

            EdlibAlignResult result = edlibAlign(
                query.c_str(), query.length(),
                read_seq.c_str(), read_seq.length(),
                edlibNewAlignConfig(-1,
                                    EDLIB_MODE_HW,
                                    EDLIB_TASK_LOC,
                                    kWildcardEqualities,
                                    8)
            );
            int max_allowed_edits = static_cast<int>(query.length() * 0.1); // 10% error tolerance hardcoded

            if (result.status == EDLIB_STATUS_OK &&
                result.editDistance <= max_allowed_edits &&
                result.numLocations > 0) {
                double norm_start = (static_cast<double>(result.startLocations[0]) / read_seq.length()) * 100.0;
                double norm_stop = (static_cast<double>(result.endLocations[0]) / read_seq.length()) * 100.0;
                starts.push_back(norm_start);
                stops.push_back(norm_stop);
            }
            edlibFreeAlignResult(result);
        };

        using ChunkFunc = std::function<void(const std::vector<read_streaming::sequence>&, const std::string&)>;
        chunk_streaming<read_streaming::sequence, ChunkFunc> streamer(chunk_size);
        size_t total_reads = 0;

        ChunkFunc chunk_func = [&](const std::vector<read_streaming::sequence>& chunk, const std::string&) {
            total_reads += chunk.size();
            for (const auto& read : chunk) {
                for (auto& part : parts) {
                    record_alignment(part.seq, read.seq, part.forward_starts, part.forward_stops);
                    record_alignment(part.rc_seq, read.seq, part.reverse_starts, part.reverse_stops);
                }
            }
        };
        streamer.process_chunks(fastq_path, chunk_func, num_threads, max_reads);


        auto mean_or_zero = [](const std::vector<double>& vals) {
            if (vals.empty()) return 0.0;
            return std::accumulate(vals.begin(), vals.end(), 0.0) / static_cast<double>(vals.size());
        };

        std::vector<ordered_element> ordered_parts;

        for (const auto& part : parts) {
            size_t fwd_hits = part.forward_starts.size();
            size_t rev_hits = part.reverse_starts.size();
            if (fwd_hits == 0 && rev_hits == 0) {
                if (verbose) {
                    RAD_COUT << "[generate_layout_from_parts] Skipping part " << part.id
                              << " (no perfect alignments found)\n";
                }
                continue;
            }
            auto forw_ratio = static_cast<double>(fwd_hits) /
                               static_cast<double>(fwd_hits + rev_hits);
                               
            auto rev_ratio = static_cast<double>(rev_hits) /
                               static_cast<double>(fwd_hits + rev_hits);

            // Record a forward-oriented entry if it has any support.
            if (fwd_hits > 0) {
                ordered_parts.push_back({
                    part.id,
                    part.id,
                    part.seq,
                    (forw_ratio >= 0.95 ? std::string("forward_only") : std::string("forward")),
                    mean_or_zero(part.forward_starts),
                    fwd_hits,
                    rev_hits
                });
            }

            // Record a reverse-oriented entry if it has any support.
            if (rev_hits > 0) {
                ordered_parts.push_back({
                    part.id,
                    "rc_" + part.id,
                    part.rc_seq,
                    (rev_ratio >= 0.95 ? std::string("reverse_only") : std::string("reverse")),
                    mean_or_zero(part.reverse_starts),
                    rev_hits,
                    fwd_hits
                });
            }
        }

        auto share_same_site = [](const ordered_element& a, const ordered_element& b) {
            double pos_diff = std::abs(a.mean_start - b.mean_start);
            if (pos_diff > 1.0) return false;  // within ~1% of the read length

            size_t min_len = std::min(a.seq.size(), b.seq.size());
            return a.seq.compare(0, min_len, b.seq, 0, min_len) == 0;
        };

        auto prefer_part = [](const ordered_element& a, const ordered_element& b) {
            if ((a.dominant_hits + a.alternate_hits) != (b.dominant_hits + b.alternate_hits)) {
                return (a.dominant_hits + a.alternate_hits) > (b.dominant_hits + b.alternate_hits);
            }
            if (a.seq.size() != b.seq.size()) {
                return a.seq.size() > b.seq.size();
            }
            if (a.dominant_hits != b.dominant_hits) {
                return a.dominant_hits > b.dominant_hits;
            }
            return a.id < b.id;
        };

        std::vector<ordered_element> deduped_parts;
        for (const auto& part : ordered_parts) {
            bool merged = false;
            for (auto& kept : deduped_parts) {
                if (share_same_site(part, kept)) {
                    if (prefer_part(part, kept)) {
                        kept = part;
                    }
                    merged = true;
                    break;
                }
            }
            if (!merged) {
                deduped_parts.push_back(part);
            }
        }

        // Remove subset/overlapping elements that sit at the same site but are weaker.
        std::vector<bool> drop(deduped_parts.size(), false);
        for (size_t i = 0; i < deduped_parts.size(); ++i) {
            if (drop[i]) continue;
            for (size_t j = i + 1; j < deduped_parts.size(); ++j) {
                if (drop[j]) continue;

                const auto& a = deduped_parts[i];
                const auto& b = deduped_parts[j];
                if (std::abs(a.mean_start - b.mean_start) > 1.0) continue;

                bool a_contains_b = a.seq.find(b.seq) != std::string::npos;
                bool b_contains_a = b.seq.find(a.seq) != std::string::npos;
                if (!a_contains_b && !b_contains_a) continue;

                const auto& keep = prefer_part(a, b) ? a : b;
                const auto& drop_elem = (&keep == &a) ? b : a;

                drop[(&keep == &a) ? j : i] = true;
                if (verbose) {
                    RAD_COUT << "[generate_layout_from_parts] Dropping subset element "
                              << drop_elem.class_id << " at ~" << drop_elem.mean_start
                              << "% in favor of " << keep.class_id << "\n";
                }
            }
        }
        {
            std::vector<ordered_element> kept;
            for (size_t i = 0; i < deduped_parts.size(); ++i) {
                if (!drop[i]) kept.push_back(deduped_parts[i]);
            }
            deduped_parts.swap(kept);
        }

        if (deduped_parts.empty()) {
            throw std::runtime_error("No alignments found for any parts; cannot build layout");
        }

        const double min_hit_fraction = 0.5;      // At least 50% of sampled reads
        const size_t min_absolute_hits = 3;         // And at least 3 total hits
        std::vector<ordered_element> filtered_parts;
        for (const auto& part : deduped_parts) {
            size_t total_hits = part.dominant_hits + part.alternate_hits;
            bool keep = total_hits >= min_absolute_hits;
            if (total_reads > 0) {
                keep = keep &&
                       (static_cast<double>(total_hits) >= min_hit_fraction * static_cast<double>(total_reads));
            }

            if (!keep) {
                if (verbose) {
                    RAD_COUT << "[generate_layout_from_parts] Filtering " << part.class_id
                              << " (support=" << total_hits
                              << ", sampled_reads=" << total_reads << ")\n";
                }
                continue;
            } else {
                if (verbose) {
                    RAD_COUT << "[generate_layout_from_parts] Keeping " << part.class_id
                              << " (support=" << total_hits
                              << ", sampled_reads=" << total_reads << ")\n";
                }
            }
            filtered_parts.push_back(part);
        }

        if (filtered_parts.empty()) {
            throw std::runtime_error("All candidate parts were filtered as spurious; cannot build layout");
        }

        deduped_parts.swap(filtered_parts);

        std::unordered_map<std::string, ordered_element> base_parts_by_id;
        std::unordered_map<std::string, ordered_element> forward_only_by_id;
        std::unordered_map<std::string, ordered_element> reverse_only_by_id;

        for (const auto& part : deduped_parts) {
            if (part.direction == "forward_only") {
                auto it = forward_only_by_id.find(part.id);
                if (it == forward_only_by_id.end() || prefer_part(part, it->second)) {
                    forward_only_by_id[part.id] = part;
                }
                continue;
            }
            if (part.direction == "reverse_only") {
                auto it = reverse_only_by_id.find(part.class_id);
                if (it == reverse_only_by_id.end() || prefer_part(part, it->second)) {
                    reverse_only_by_id[part.class_id] = part;
                }
                continue;
            }
            if (part.direction != "forward") {
                continue;
            }

            auto it = base_parts_by_id.find(part.id);
            if (it == base_parts_by_id.end() || prefer_part(part, it->second)) {
                base_parts_by_id[part.id] = part;
            }
        }

        std::vector<ordered_element> base_parts;
        base_parts.reserve(base_parts_by_id.size());
        for (const auto& kv : base_parts_by_id) {
            base_parts.push_back(kv.second);
        }

        if (base_parts.empty()) {
            throw std::runtime_error("No forward-oriented parts available to seed layout");
        }

        std::sort(base_parts.begin(), base_parts.end(), [](const ordered_element& a, const ordered_element& b) {
            return a.mean_start < b.mean_start;
        });

        std::filesystem::path tmp_layout =
            std::filesystem::temp_directory_path() / "rad_parts_forward_layout.csv";
        {
            std::ofstream tmp(tmp_layout.string());
            tmp << "Read Layout\nid,seq,masked_seq,expected_length,type,class,direction,class_id,whitelist,order\n";
            int tmp_order = 1;
            for (const auto& part : base_parts) {
                std::string seq = part.seq;
                if (part.direction == "reverse") {
                    continue;
                }
                tmp << part.id << ","
                    << seq << ","
                    << "" << ","
                    << seq.length() << ","
                    << "static,"
                    << part.id << ","
                    << "forward,"
                    << part.id << ",,"
                    << tmp_order++ << "\n";
            }
        }

        prep_new_layout(tmp_layout.string(), verbose);

        auto forward_only_parts = [&]() {
            std::vector<ordered_element> parts;
            parts.reserve(forward_only_by_id.size());
            for (const auto& kv : forward_only_by_id) parts.push_back(kv.second);
            std::sort(parts.begin(), parts.end(), [](const ordered_element& a, const ordered_element& b) {
                return a.mean_start < b.mean_start;
            });
            return parts;
        }();

        auto reverse_only_parts = [&]() {
            std::vector<ordered_element> parts;
            parts.reserve(reverse_only_by_id.size());
            for (const auto& kv : reverse_only_by_id) parts.push_back(kv.second);
            std::sort(parts.begin(), parts.end(), [](const ordered_element& a, const ordered_element& b) {
                return a.mean_start > b.mean_start;
            });
            return parts;
        }();

        // Rebuild the layout to insert single-sided elements in the right spots.
        std::vector<ReadElement> ordered_elements;
        ordered_elements.reserve(layout.size());
        for (const auto& elem : by_order()) {
            ordered_elements.push_back(elem);
        }

        auto find_stop_index = [&](const std::string& dir) {
            for (size_t i = 0; i < ordered_elements.size(); ++i) {
                if (ordered_elements[i].direction == dir && ordered_elements[i].global_class == "stop") {
                    return static_cast<int>(i);
                }
            }
            return -1;
        };

        int fwd_stop_idx = find_stop_index("forward");
        int rev_stop_idx = find_stop_index("reverse");

        if (fwd_stop_idx == -1) {
            throw std::runtime_error("prep_new_layout did not produce a forward stop element");
        }

        std::vector<ReadElement> rebuilt;
        rebuilt.reserve(ordered_elements.size() + forward_only_parts.size() + reverse_only_parts.size());

        // Forward block (without the forward stop)
        for (int i = 0; i < fwd_stop_idx; ++i) {
            rebuilt.push_back(ordered_elements[i]);
        }

        for (const auto& part : forward_only_parts) {
            rebuilt.emplace_back(part.id, part.seq, "", static_cast<int>(part.seq.length()),
                                 "static", 0, "forward", part.id);
            if (verbose) {
                int total_cts = part.dominant_hits + part.alternate_hits;
                RAD_COUT << "[generate_layout_from_parts] Inserted forward_only " << part.class_id
                          << " with " << total_cts << " hits at ~" << part.mean_start << "% in the sequence\n";
            }
        }

        // Forward stop goes after forward-only inserts
        rebuilt.push_back(ordered_elements[fwd_stop_idx]);

        // Reverse block (everything after forward stop up to reverse stop if present)
        size_t reverse_start = static_cast<size_t>(fwd_stop_idx + 1);
        size_t reverse_end = (rev_stop_idx == -1) ? ordered_elements.size() : static_cast<size_t>(rev_stop_idx);
        for (size_t i = reverse_start; i < reverse_end; ++i) {
            rebuilt.push_back(ordered_elements[i]);
        }

        for (const auto& part : reverse_only_parts) {
            rebuilt.emplace_back(part.class_id, part.seq, "", static_cast<int>(part.seq.length()),
                                 "static", 0, "reverse", part.id);
            if (verbose) {
                int total_cts = part.dominant_hits + part.alternate_hits;
                RAD_COUT << "[generate_layout_from_parts] Inserted reverse_only " << part.class_id
                          << " with " << total_cts << " hits at ~" << part.mean_start << "% in the sequence\n";
            }
        }

        if (rev_stop_idx != -1) {
            rebuilt.push_back(ordered_elements[rev_stop_idx]);
            for (size_t i = static_cast<size_t>(rev_stop_idx + 1); i < ordered_elements.size(); ++i) {
                rebuilt.push_back(ordered_elements[i]);
            }
        }

        layout.clear();
        int new_order = 1;
        for (auto& elem : rebuilt) {
            elem.order = new_order++;
            layout.insert(elem);
        }

        write_to_csv(output_base, "layout");
    }

    // Export layout to CSV
    void write_to_csv(const std::string& base_path = "", const std::string& which_file = "both") const {
        std::string path_prefix = base_path;
        if (path_prefix.empty()) {
            const char* home = std::getenv("HOME");
            if (home == nullptr || *home == '\0') {
                throw std::runtime_error(
                    "Cannot use the default layout output path because HOME is not set");
            }
            path_prefix = std::string(home) + "/Desktop/";
        }

        auto checked_replace = [](const std::string& path,
                                  const std::string& payload) {
            static std::atomic<unsigned long long> counter{0};
            const std::string temporary =
                path + ".tmp." +
                std::to_string(std::chrono::steady_clock::now()
                                   .time_since_epoch()
                                   .count()) + "." +
                std::to_string(counter.fetch_add(1, std::memory_order_relaxed));
            std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
            if (!out) {
                throw std::runtime_error("Failed to open layout output: " +
                                         temporary);
            }
            out.write(payload.data(),
                      static_cast<std::streamsize>(payload.size()));
            out.flush();
            if (!out) {
                out.close();
                std::error_code remove_error;
                std::filesystem::remove(temporary, remove_error);
                throw std::runtime_error("Failed to write layout output: " +
                                         path);
            }
            out.close();
            if (out.fail()) {
                std::error_code remove_error;
                std::filesystem::remove(temporary, remove_error);
                throw std::runtime_error("Failed to close layout output: " +
                                         path);
            }
            std::error_code rename_error;
            std::filesystem::rename(temporary, path, rename_error);
            if (rename_error) {
                std::error_code remove_error;
                std::filesystem::remove(temporary, remove_error);
                throw std::runtime_error("Failed to install layout output " +
                                         path + ": " +
                                         rename_error.message());
            }
        };
        
        if(which_file == "layout" || which_file == "both") {
            // Basic layout information
            std::ostringstream layout_file;
            layout_file << "Read Layout" << (sequencing_type.empty() ? "" : ":" + sequencing_type)
                        << ",,,,,,,,,,,\n";
            layout_file << "id,seq,masked_seq,expected_length,length_candidates,flags,type,class,direction,class_id,whitelist,order\n";
            for (const auto& element : by_order()) {
                layout_file << "\"" << element.class_id << "\","
                    << "\"" << element.seq << "\","
                    << (element.masked_seq.empty() ? "" : element.masked_seq) << ","
                    << (element.expected_length.has_value() ? std::to_string(*element.expected_length) : "") << ","
                    << "\"" << format_length_candidates(element.length_candidates) << "\","
                    << "\"" << element.flags << "\","
                    << element.type << ","
                    << element.global_class << ","
                    << element.direction << ","
                    << element.class_id << ","
                    << (element.whitelist_path.empty() ? "" : element.whitelist_path) << ","
                    << element.order << "\n";
            }
            checked_replace(path_prefix + "_layout.csv", layout_file.str());
        }

        if(which_file == "pos_map" || which_file == "both") {
            // Alignment and threshold information
            std::ostringstream align_file;
            // Concat-HMM calibration (hmm_* columns) is appended only when present.
            // Readers find columns by name and ignore extra ones.
            std::map<std::string, std::map<std::string, std::string>> hmm_cells;
            std::vector<std::string> hmm_cols;
            if (!concat_params.empty()) {
                hmm_cells = concat_hmm::params_to_columns(concat_params);
                for (const auto& row : hmm_cells)
                    for (const auto& col : row.second)
                        if (std::find(hmm_cols.begin(), hmm_cols.end(), col.first) == hmm_cols.end())
                            hmm_cols.push_back(col.first);
                std::sort(hmm_cols.begin(), hmm_cols.end());
            }
            align_file << "id,primary_start,primary_stop,secondary_start,secondary_stop,"
                << "misalign_lower,misalign_mean,misalign_upper,"
                << "align_start,var_start,align_stop,var_stop,"
                << "misalign_start,mvar_start,misalign_stop,mvar_stop";
            for (const auto& col : hmm_cols) align_file << "," << col;
            align_file << "\n";
            
            for (const auto& element : by_order()) {
                align_file << element.class_id;
                
                // Position information
                auto pos_it = position_map.find(element.class_id);
                if (pos_it != position_map.end()) {
                    align_file << "," << pos_it->second.primary_start.to_string()
                            << "," << pos_it->second.primary_stop.to_string()
                            << "," << pos_it->second.secondary_start.to_string()
                            << "," << pos_it->second.secondary_stop.to_string();
                } else {
                    // Four commas for four empty fields
                    align_file << ",,,,";  
                }

                // Misalignment thresholds and positions
                if (element.misalignment_threshold) {
                    const auto& [lower, mean, upper] = *element.misalignment_threshold;
                    align_file << "," << lower << "," << mean << "," << upper;
                    
                    if (element.aligned_positions) {
                        const auto& [start_mean, start_var] = element.aligned_positions->start_stats;
                        const auto& [stop_mean, stop_var] = element.aligned_positions->stop_stats;
                        align_file << "," << start_mean << "," << start_var 
                                << "," << stop_mean << "," << stop_var;
                    } else {
                        // Four commas for missing aligned positions
                        align_file << ",,,,";  
                    }
                    
                    if (element.misaligned_positions) {
                        const auto& [mis_start_mean, mis_start_var] = element.misaligned_positions->start_stats;
                        const auto& [mis_stop_mean, mis_stop_var] = element.misaligned_positions->stop_stats;
                        align_file << "," << mis_start_mean << "," << mis_start_var 
                                << "," << mis_stop_mean << "," << mis_stop_var;
                    } else {
                        // Four commas for missing misaligned positions
                        align_file << ",,,,";
                    }
                } else {
                    // If no misalignment threshold, need commas for all missing fields
                    align_file << ",,,";  // Three commas for threshold fields
                    align_file << ",,,,";  // Four commas for aligned positions
                    align_file << ",,,,";  // Four commas for misaligned positions
                }
                if (!hmm_cols.empty()) {
                    const auto row_it = hmm_cells.find(element.class_id);
                    for (const auto& col : hmm_cols) {
                        align_file << ",";
                        if (row_it == hmm_cells.end()) continue;
                        const auto cell_it = row_it->second.find(col);
                        if (cell_it != row_it->second.end()) align_file << cell_it->second;
                    }
                }
                align_file << "\n";
            }
            // Mirrored anchors and poly-tail states need not be explicit
            // layout elements. Preserve their calibrated emissions as rows
            // with no positional fields; they do not create layout elements.
            for (const auto& row : hmm_cells) {
                if (by_id().find(row.first) != by_id().end()) continue;
                align_file << row.first << std::string(15, ',');
                for (const auto& col : hmm_cols) {
                    align_file << ',';
                    const auto cell = row.second.find(col);
                    if (cell != row.second.end()) align_file << cell->second;
                }
                align_file << '\n';
            }
            checked_replace(path_prefix + "_position_map.csv",
                            align_file.str());
        }
    }

    // Import layout from CSVs
    void import_from_csv(const std::string& base_path = "", bool verbose = false) {
        // Read main layout file
        const std::string layout_path = base_path + "_layout.csv";
        if (std::regex_search(get_rl_mode(layout_path), std::regex("bulk"))) sequencing_type = "bulk";
        csv::CSVFormat format;
        format.delimiter(',').quote('"').header_row(detect_header_row(layout_path));
        format.variable_columns(csv::VariableColumnPolicy::THROW);
        csv::CSVReader layout_reader(layout_path, format);

        //int order_counter = 1;
        for (auto& row : layout_reader) {
            auto expected_spec = get_optional_csv_string(row, "expected_length");
            auto length_spec = get_optional_csv_string(row, "length_candidates");
            auto length_candidates = parse_length_candidates(length_spec.has_value() ? length_spec : expected_spec);
            std::optional<int> expected_length = length_candidates.empty()
                ? std::nullopt
                : std::optional<int>(length_candidates.back());
            ReadElement elem(
                row["id"].get<std::string>(),
                row["seq"].get<std::string>(),
                row["masked_seq"].get<std::string>(),
                expected_length,
                row["type"].get<std::string>(),
                std::stoi(row["order"].get<std::string>()),
                row["direction"].get<std::string>(),
                row["class"].get<std::string>(),
                row["whitelist"].get<std::string>()
            );
            elem.length_candidates = std::move(length_candidates);
            elem.flags = get_optional_csv_string(row, "flags").value_or("");
            layout.insert(elem);
        }

        // Read position and alignment info
        csv::CSVFormat new_format;
        new_format.delimiter(',').quote('"').header_row(0);
        new_format.variable_columns(csv::VariableColumnPolicy::KEEP);
        csv::CSVReader pos_reader(base_path + "_position_map.csv", new_format);
        for (auto& row : pos_reader) {
            std::string id = row["id"].get<std::string>();
            if(verbose){
                RAD_COUT << "[import_from_both] Processing ID: " << id << "\n";
            }
            auto element_it = by_id().find(id);
            if (element_it == by_id().end()){
                if(verbose){
                    RAD_COUT << "[import_from_both] Skipping ID not found in layout: " << id << "\n";
                }
                continue;
            }
        // Parse position info
            if (!row["primary_start"].get<std::string>().empty()) {
            ReferencePositions pos;
            if(verbose){
                RAD_COUT << "[import_from_both] This did work!\n";
                RAD_COUT << "[import_from_both] Row data for " << id << ":\n";
                RAD_COUT << "[import_from_both] primary_start value: '" << row["primary_start"].get<std::string>() << "'\n";
            }
            try {
                    pos.primary_start = ParsedPosition::from_string(row["primary_start"].get<std::string>(), verbose);
                    pos.primary_stop = ParsedPosition::from_string(row["primary_stop"].get<std::string>(), verbose);
                    pos.secondary_start = ParsedPosition::from_string(row["secondary_start"].get<std::string>(), verbose);
                    pos.secondary_stop = ParsedPosition::from_string(row["secondary_stop"].get<std::string>(), verbose);
            
                    position_map[id] = pos;
                    if(verbose){
                        RAD_COUT << "[import_from_both] Processing element " << id << " with positions:\n"
                        << "Primary: " << pos.primary_start.to_string() << " -> " << pos.primary_stop.to_string() << "\n"
                        << "Secondary: " << pos.secondary_start.to_string() << " -> " << pos.secondary_stop.to_string() << "\n";
                    }
                    
                    by_id().modify(element_it, [pos](ReadElement& elem) {
                        elem.ref_pos = pos;
                    });
                } catch (const std::exception& e) {
                    RAD_CERR << "Error processing positions for " << id << ": " << e.what() << "\n";
                }
            } else {
                if(verbose){
                    RAD_COUT << "[import_from_both] No primary start value!\n";
                    RAD_COUT << "[import_from_both] Row data for " << id << ":\n";    
                }
        }
        // Parse misalignment and alignment stats
            if (!row["misalign_lower"].is_null()) {
                std::tuple<int, int, int> threshold{
                    std::stoi(row["misalign_lower"].get<std::string>()),
                    std::stoi(row["misalign_mean"].get<std::string>()),
                    std::stoi(row["misalign_upper"].get<std::string>())
                };

                AlignmentPositions align_pos{
                    {std::stod(row["align_start"].get<std::string>()), 
                        std::stod(row["var_start"].get<std::string>())},
                    {std::stod(row["align_stop"].get<std::string>()), 
                        std::stod(row["var_stop"].get<std::string>())}
                };

                AlignmentPositions misalign_pos{
                    {std::stod(row["misalign_start"].get<std::string>()), 
                        std::stod(row["mvar_start"].get<std::string>())},
                    {std::stod(row["misalign_stop"].get<std::string>()), 
                        std::stod(row["mvar_stop"].get<std::string>())}
                };

                by_id().modify(element_it, [&](ReadElement& elem) {
                    elem.misalignment_threshold = threshold;
                    elem.aligned_positions = align_pos;
                    elem.misaligned_positions = misalign_pos;
                });
            }
        }
    }
    
    // import a read layout .csv
    void import_read_layout(const std::string& layout_csv, bool verbose,
                            const std::string& mode_source = ""){
        using namespace csv;
        // Caches written before the mode was persisted carry no title row; fall back to
        // the layout the caller supplied, since losing bulk here drops every read.
        std::string mode = get_rl_mode(layout_csv);
        if (mode.empty() && !mode_source.empty()) mode = get_rl_mode(mode_source);
        if (std::regex_search(mode, std::regex("bulk"))) sequencing_type = "bulk";
        CSVFormat fmt;
        fmt.delimiter(',').quote('"').header_row(detect_header_row(layout_csv))
        .variable_columns(VariableColumnPolicy::THROW);
        CSVReader reader(layout_csv, fmt);
        for (auto &row : reader) {
            auto expected_spec = get_optional_csv_string(row, "expected_length");
            auto length_spec = get_optional_csv_string(row, "length_candidates");
            auto length_candidates = parse_length_candidates(length_spec.has_value() ? length_spec : expected_spec);
            std::optional<int> expected_length = length_candidates.empty()
                ? std::nullopt
                : std::optional<int>(length_candidates.back());
            ReadElement elem(
                row["id"].get<std::string>(),
                row["seq"].get<std::string>(),
                row["masked_seq"].get<std::string>(),
                expected_length,
                row["type"].get<std::string>(),
                std::stoi(row["order"].get<std::string>()),
                row["direction"].get<std::string>(),
                row["class"].get<std::string>(),
                row["whitelist"].get<std::string>()
            );
            elem.length_candidates = std::move(length_candidates);
            elem.flags = get_optional_csv_string(row, "flags").value_or("");
            layout.insert(elem);
        }
        if (verbose) {
            std::vector<ReadElement> all;
            all.reserve(layout.size());
            for (auto const &e : layout) {
                all.push_back(e);
            }
            std::vector<ReadElement> fwd, rev;
            for (auto const &e : all) {
                if (e.direction == "forward")  fwd.push_back(e);
                if (e.direction == "reverse")  rev.push_back(e);
            }
            // sort by order field
            auto cmp_order = [](auto const &a, auto const &b){
                return a.order < b.order;
            };
            std::sort(fwd.begin(), fwd.end(), cmp_order);
            std::sort(rev.begin(), rev.end(), cmp_order);

            // helper to print a list
            auto print_list = [&](auto const &vec, const char *name){
                RAD_COUT << "[read_layout]:" << name;
                for (auto const &e : vec) {
                    std::string seq = e.seq.empty() ? "" : e.seq;
                    RAD_COUT << "["
                            << e.class_id << ":"
                            << seq
                            << "]";
                }
                RAD_COUT << "\n";
            };
            print_list(fwd, "F:");
            print_list(rev, "R:");
        }
    }
    
    // import a position map .csv
    void import_position_map(const std::string& map_csv,bool verbose){
        concat_params = concat_hmm::Params{};
        concat_model.reset();
        using namespace csv;
        // set up CSV reader for the position map
        CSVFormat fmt;
        fmt.delimiter(',').quote('"').header_row(0)
            .variable_columns(VariableColumnPolicy::KEEP);
        CSVReader reader(map_csv, fmt);
        // Optional concat-HMM calibration columns (hmm_*).
        std::vector<std::string> hmm_cols;
        std::map<std::string, std::size_t> hmm_column_indices;
        const auto column_names = reader.get_col_names();
        for (std::size_t i = 0; i < column_names.size(); ++i) {
            const auto& name = column_names[i];
            if (name.rfind("hmm_", 0) == 0) {
                hmm_cols.push_back(name);
                hmm_column_indices.emplace(name, i);
            }
        }
        std::map<std::string, std::map<std::string, std::string>> hmm_cells;
        // iterate rows
        for (auto &row : reader) {
            std::string id = row["id"].get<std::string>();
            // Import parameter-only mirrored states before ignoring their
            // nonexistent position-map entry. Bounds checks also preserve
            // support for older short rows without swallowing failures.
            for (const auto& col : hmm_cols) {
                const auto index = hmm_column_indices.at(col);
                if (index >= row.size()) continue;
                auto field = row[index];
                if (field.is_null()) continue;
                std::string cell = field.get<std::string>();
                if (!cell.empty()) hmm_cells[id][col] = std::move(cell);
            }
            // find the corresponding element from the layout
            auto element_it = by_id().find(id);
            if (element_it == by_id().end()) {
                if (verbose) {
                    RAD_COUT << "[pos_map] skipping unknown id: " << id << "\n";
                }
                continue;
            }
            // only proceed if we have a non-empty primary_start
            const std::string &ps = row["primary_start"].get<std::string>();
            if (!ps.empty()) {
                // parse reference positions
                ReferencePositions rp;
                rp.primary_start   = ParsedPosition::from_string(ps, verbose);
                rp.primary_stop    = ParsedPosition::from_string(
                                            row["primary_stop"].get<std::string>(), verbose);
                rp.secondary_start = ParsedPosition::from_string(
                                            row["secondary_start"].get<std::string>(), verbose);
                rp.secondary_stop  = ParsedPosition::from_string(
                                            row["secondary_stop"].get<std::string>(), verbose);

                // store in the map and modify the element
                position_map[id] = rp;
                by_id().modify(element_it, [rp](ReadElement &e){
                    e.ref_pos = rp;
                });
            if (verbose) {
                RAD_COUT << "[pos_map] " << id
                            << " -> primary("
                            << rp.primary_start.to_string() << "-"
                            << rp.primary_stop.to_string()   << ")\n";
            }
        }

        // now, if misalignment stats are present, parse and attach
        if (!row["misalign_lower"].is_null()) {
                auto lower = std::stoi(row["misalign_lower"].get<std::string>());
                auto mean  = std::stoi(row["misalign_mean" ].get<std::string>());
                auto upper = std::stoi(row["misalign_upper"].get<std::string>());
                std::tuple<int, int, int> threshold{lower, mean, upper};
                if(verbose){
                    RAD_COUT << "[pos_map] " << id
                                << " -> misalignment threshold: "
                                << lower << "|" << mean << "|" << upper << "\n";
                }

                // alignment positions
                AlignmentPositions align_pos{
                    {
                        std::stod(row["align_start"].get<std::string>()),
                        std::stod(row["var_start"].get<std::string>())
                    },
                    {
                        std::stod(row["align_stop"].get<std::string>()),
                        std::stod(row["var_stop"].get<std::string>())
                    }
                };
                // misalignment positions
                AlignmentPositions misalign_pos{
                    {
                        std::stod(row["misalign_start"].get<std::string>()),
                        std::stod(row["mvar_start"].get<std::string>())
                    },
                    {
                        std::stod(row["misalign_stop" ].get<std::string>()),
                        std::stod(row["mvar_stop"].get<std::string>())
                    }
                };

                // attach them both in one modify call
                by_id().modify(element_it, [&](ReadElement &e){
                    e.misalignment_threshold = threshold;
                    e.aligned_positions      = align_pos;
                    e.misaligned_positions   = misalign_pos;
                });

                if (verbose) {
                    RAD_COUT << "[pos_map] " << id
                                << " -> misalignment threshold: "
                                << lower << "|" << mean << "|" << upper << "\n";
                }
            }
        }
        if (!hmm_cells.empty()) {
            // Malformed cells are retained in Params::bad_cells so model
            // construction can disable the HMM explicitly. Allocation or
            // other resource failures must reach the host unchanged.
            concat_params = concat_hmm::params_from_columns(to_concat_layout_spec(), hmm_cells);
            if (verbose)
                RAD_COUT << "[pos_map] concat-HMM calibration loaded (" << hmm_cols.size()
                          << " hmm_* columns" << (concat_params.guard_failed ? ", guard failed" : "") << ")\n";
        }
    }

    // Generate position map
    void generate_position_mapping(bool verbose = false){
        auto& dir_ordered = by_dir_order();
        auto current_direction = dir_ordered.begin();
        while (current_direction != dir_ordered.end()) {
            auto direction = current_direction->direction;
            auto end_of_direction = dir_ordered.upper_bound(
                boost::make_tuple(direction, std::numeric_limits<int>::max())
            );
            process_reads(current_direction, end_of_direction, verbose);
                auto& id_index = by_id();
                for (const auto& pos_pair : position_map) {
                    auto elem_it = id_index.find(pos_pair.first);
                    if (elem_it != id_index.end()) {
                        id_index.modify(elem_it, [positions=pos_pair.second](ReadElement& elem) {
                            elem.ref_pos = positions;
                        });
                    }
                }
            current_direction = end_of_direction;
        }
    }

    // load whitelist
    void load_wl(
        std::optional<std::string> wl_path,
        std::optional<int> mut,
        bool verbose,
        int nthreads,
        std::optional<std::string> typed_global = std::nullopt,
        std::optional<std::string> typed_true = std::nullopt) {
        using namespace std::chrono;
        auto t0 = high_resolution_clock::now();

        // Clear any old state
        wl_map.lists.clear();
        wl_map.maps.clear();

        // A joint element's axis:mask declaration is not the generic
        // global:true dual-whitelist syntax. Resolve it before catalog import.
        auto joint_group = [](const ReadElement& elem) -> std::optional<int> {
            std::string flags = elem.flags;
            std::transform(flags.begin(), flags.end(), flags.begin(), [](unsigned char c) {
                return static_cast<char>(std::tolower(c));
            });
            for (char& c : flags) {
                if (c == ',' || c == ';' || c == '|') c = ' ';
            }
            std::istringstream tokens(flags);
            std::string token;
            std::optional<int> group;
            while (tokens >> token) {
                if (token.rfind("joint_barcode", 0) != 0) continue;
                std::smatch match;
                static const std::regex grammar("^joint_barcode(?::([0-9]+)(?:_([0-9]+))?)?$");
                if (!std::regex_match(token, match, grammar) || group.has_value()) {
                    throw std::runtime_error("Invalid joint barcode flags for " + elem.class_id);
                }
                group = match[1].matched ? std::stoi(match[1].str()) : 0;
            }
            return group;
        };
        std::map<int, std::vector<const ReadElement*>> joint_groups;
        std::unordered_map<std::string, std::string> joint_axis_specs;
        std::string mask_csv_path;
        std::optional<int> masked_group;
        for (const auto& elem : by_order()) {
            if (elem.global_class != "barcode" || elem.type != "variable") continue;
            auto group = joint_group(elem);
            if (!group) continue;
            auto colon = elem.whitelist_path.find(':');
            std::string axis = elem.whitelist_path.substr(0, colon);
            joint_axis_specs[seq_utils::remove_rc(elem.class_id)] = axis;
            if (elem.direction == "forward") joint_groups[*group].push_back(&elem);
            if (colon == std::string::npos) continue;
            std::string mask = elem.whitelist_path.substr(colon + 1);
            if (axis.empty() || mask.empty() || mask.find(':') != std::string::npos) {
                throw std::runtime_error("Expected axis:mask for joint barcode " + elem.class_id);
            }
            std::string resolved = whitelist_utils::kit_to_path(mask);
            if ((masked_group && *masked_group != *group) ||
                (!mask_csv_path.empty() && mask_csv_path != resolved)) {
                throw std::runtime_error("Conflicting spatial mask declarations; only one masked joint group is supported");
            }
            masked_group = *group;
            mask_csv_path = resolved;
        }
        if (masked_group) {
            if (wl_path || typed_global || typed_true) {
                throw std::runtime_error("A whitelist override cannot be combined with a declared joint spatial mask");
            }
            if (joint_groups.size() != 1 || joint_groups[*masked_group].size() != 2) {
                throw std::runtime_error("A spatial mask requires exactly one joint group with two forward-defined barcode axes");
            }
        }

        // 1) collect static seqs for filter_bcs
        std::vector<std::string> static_seqs;
        for (auto const &elem : layout) {
            if (elem.type=="static" && !elem.seq.empty())
                static_seqs.push_back(elem.seq);
        }

        // import & alias
        for (auto const &elem : layout) {
            if (elem.global_class != "barcode" || elem.type != "variable")
                continue;

            // normalize class_id -> key
            std::string key = elem.class_id;
            if (seq_utils::is_rc(key)){
                key = seq_utils::remove_rc(key);
            }
            
            //select default length if at all, otherwise 
            size_t default_length = 16;
            if (auto lit = layout.find(elem.class_id);
                lit != layout.end() && lit->expected_length)
            {
                default_length = *lit->expected_length;
            }

            const bool typed = typed_global.has_value() || typed_true.has_value();
            std::string spec = wl_path ? wl_path.value() : elem.whitelist_path;
            if (!wl_path && joint_axis_specs.count(key)) spec = joint_axis_specs.at(key);
            if (!typed && spec.empty()) {
                throw std::invalid_argument(
                    "No whitelist was provided for barcode element " +
                    elem.class_id);
            }

            std::string path;
            if (typed) {
                // Length-prefix the role-bearing specifications so distinct
                // pairs cannot alias in the shared whitelist pool.
                const std::string global = typed_global.value_or("");
                const std::string truth = typed_true.value_or("");
                path = "typed:" + std::to_string(default_length) + ":g" +
                       std::to_string(global.size()) + ":" + global + ":t" +
                       std::to_string(truth.size()) + ":" + truth;
                spec = "global=" + (global.empty() ? "[none]" : global) +
                       ", custom=" + (truth.empty() ? "[none]" : truth);
            } else {
                path = whitelist_utils::kit_to_path(spec);
            }
            RAD_COUT << "[load_wl] [" << elem.class_id << "] [" << spec
                      << "] @ " << path << "\n";

            // import into the pool if not already
            // Isolate masked axes from other elements that reuse their catalog.
            const std::string pool_key = masked_group && joint_axis_specs.count(key)
                ? path + "#joint:" + key : path;
            auto pit = wl_map.lists.find(pool_key);
            if (pit == wl_map.lists.end()) {
                if (verbose) RAD_COUT << "[load_wl] Loading from " << path << " ...\n";
                auto t1 = high_resolution_clock::now();
                memory_utils::get_rss();
                // import & possibly generate mismatches
                whitelist::wl_entry entry = typed
                    ? wl_map.import_whitelist_typed(
                          typed_global, typed_true, verbose, default_length)
                    : wl_map.import_whitelist(spec, verbose, default_length);
                // hardcoded whitelist size limit for entry so that we don't generate too many mismatches
                memory_utils::get_rss();
                // For large true whitelists, pre-build a 3-part seed index and skip mismatch expansion.
                const size_t true_unique_size = entry.true_bcs.unique_val_size();
                if (true_unique_size > 10000) {
                    entry.build_seed_idx(verbose);
                    if (verbose) {
                        RAD_COUT << "[load_wl] Skipping generate_mismatch_barcodes for "
                                  << elem.class_id << " (deprecated path; true_unique_size="
                                  << true_unique_size << ")\n";
                    }
                    // Deprecated for large true whitelists:
                    // entry.generate_mismatch_barcodes(mut.value_or(2), verbose, nthreads);
                }
                memory_utils::get_rss();
                // build filter_bcs from static_seqs
                entry.filter_bcs.clear();
                size_t def_len = default_length;
                if (auto lit = layout.find(elem.class_id);
                    lit != layout.end() && lit->expected_length){
                    def_len = *lit->expected_length;
                }
                for (auto const &s : static_seqs) {
                    for (auto const &kmer : seq_utils::circ_kmerize(s, def_len)) {
                        int64_seq bits(kmer);
                        entry.filter_bcs.insert_bc_entry(bits);
                        auto mutations = mutation_tools::generate_lv_barcodes(bits, 2);
                        for (auto const &mut_bits : mutations) {
                            if(entry.filter_bcs.check_wl_for(mut_bits)){
                                entry.filter_bcs.insert_bc_entry(mut_bits);
                            }
                        }
                    }
                }
                auto t2 = high_resolution_clock::now();
                double ms = duration<double, std::milli>(t2 - t1).count();
                RAD_COUT << "[load_wl] loaded in " << (ms/1000.0) << " s\n";
                memory_utils::get_rss();
                pit = wl_map.lists.emplace(pool_key, std::move(entry)).first;
            }

            // alias into maps by reference
            auto &master = pit->second;
            auto [mit, inserted] = wl_map.maps.try_emplace(key, std::ref(master));
            if (inserted) {
                RAD_COUT << "[load_wl] New entry for key ='"<< key << "'\n";
                memory_utils::get_rss();

            } else if (verbose) {
                RAD_COUT << "[load_wl] Reusing entry for key ='"<< key << "'\n";
                memory_utils::get_rss();
            }

            {
                auto &E = mit->second.get();
                RAD_COUT << "[load_wl]"
                        << " ["<<elem.class_id<<":"<<key<<"]"
                        << ": global = "<<E.global_bcs.size()
                        << ", true = "  <<E.true_bcs.size()
                        << ", filter = "<<E.filter_bcs.size()
                        << "\n";
            }
        }

        // Validate the matrix strictly before the legacy loader consumes it.
        // A requested but unusable mask must never turn into an unrestricted run.
        if (masked_group) {
            const auto& members = joint_groups.at(*masked_group);
            const std::string bc1_path = whitelist_utils::kit_to_path(joint_axis_specs.at(seq_utils::remove_rc(members[0]->class_id)));
            const std::string bc2_path = whitelist_utils::kit_to_path(joint_axis_specs.at(seq_utils::remove_rc(members[1]->class_id)));
            auto axis_size = [](const std::string& path) {
                std::unordered_set<std::string> seen;
                bool first = true;
                for (const auto& raw : streaming_utils::import_text(path)) {
                    std::string line = seq_utils::trim(raw);
                    if (line.empty()) continue;
                    if (first && line == "whitelist_bcs") { first = false; continue; }
                    first = false;
                    if (line.size() > 64 || line.find_first_not_of("ACGT") != std::string::npos ||
                        !seen.insert(line).second) {
                        throw std::runtime_error("Spatial mask axes require unique DNA barcode lines: " + path);
                    }
                }
                if (seen.empty() || seen.size() > std::numeric_limits<uint16_t>::max()) {
                    throw std::runtime_error("Invalid spatial mask axis size: " + path);
                }
                return seen.size();
            };
            const size_t rows = axis_size(bc1_path), cols = axis_size(bc2_path);
            std::ifstream matrix(mask_csv_path);
            if (!matrix) throw std::runtime_error("Cannot open spatial mask: " + mask_csv_path);
            size_t row_count = 0;
            std::string line;
            while (std::getline(matrix, line)) {
                size_t col_count = 0, start = 0;
                while (true) {
                    auto comma = line.find(',', start);
                    std::string value = seq_utils::trim(line.substr(start, comma == std::string::npos ? comma : comma - start));
                    if (value != "0" && value != "1") {
                        throw std::runtime_error("spatial mask cells must be exactly 0 or 1: " + mask_csv_path);
                    }
                    ++col_count;
                    if (comma == std::string::npos) break;
                    start = comma + 1;
                }
                if (col_count != cols)
                    throw std::runtime_error("spatial mask column count does not match the BC2 axis: " + mask_csv_path);
                if (++row_count > rows)
                    throw std::runtime_error("spatial mask row count does not match the BC1 axis: " + mask_csv_path);
            }
            if (matrix.bad()) throw std::runtime_error("Could not read spatial mask: " + mask_csv_path);
            if (row_count != rows)
                throw std::runtime_error("spatial mask row count does not match the BC1 axis: " + mask_csv_path);
            auto& entry = wl_map.maps.at(seq_utils::remove_rc(members[0]->class_id)).get();
            entry.spat_wl.emplace();
            if (!entry.spat_wl->load(mask_csv_path, bc1_path, bc2_path, verbose)) {
                throw std::runtime_error("Failed to load requested spatial mask: " + mask_csv_path);
            }
            RAD_COUT << "[load_wl] Spatial mask loaded: " << entry.spat_wl->count_valid() << " valid pairs\n";
        }

        for (auto & [key, entry_ref] : wl_map.maps) {
            auto &E = entry_ref.get();
            // The true and global tables own separate storage, so removing
            // global overlaps does not invalidate this association traversal.
            E.true_bcs.for_each_association([&](const auto& barcode, const auto&) {
                E.global_bcs.remove_bc_entry(barcode);
            });

            // final report
            RAD_COUT << "[load_wl]"
                    << " [" << key << "]"
                    << ": global = "<< E.global_bcs.size()
                    << ", true = "  << E.true_bcs.size()
                    << ", filter = "<< E.filter_bcs.size()
                    << "\n";
            RAD_COUT << "[load_wl] [" << key << "] true logical associations = "
                      << E.true_bcs.association_size()
                      << ", explicit alias keys = " << E.true_bcs.physical_alias_size()
                      << "\n";
            for (auto const &kv : wl_map.maps) {
                auto const &E = kv.second.get();
                
                size_t bad_keys = E.true_bcs.validate_association_keys();
                size_t null_values = E.true_bcs.validate_association_values();
                memory_utils::get_rss();

                RAD_COUT << "Bad keys: " << bad_keys << ", Null values: " << null_values << std::endl;
            
                /*std::printf(
                    "[%s] unique_set = %.2fMiB, associations = %.2fMiB, total = %.2fMiB\n",
                        kv.first.c_str(),
                    
                        memory_utils::to_mib(bc_mem_utils::approx_unique(E.true_bcs)),
                        memory_utils::to_mib(bc_mem_utils::approx_assoc(E.true_bcs))
                );*/
                bc_mem_utils::print_memory_report(E.true_bcs, "true_bcs");
                memory_utils::get_rss();

            }
        }

        // Part tables for the exhaustive scan; true_bcs must not change after this point.
        for (auto & [pool_key, E] : wl_map.lists) {
            E.build_part_idx(verbose);
        }

        memory_utils::release_unused_heap();

        auto t3 = high_resolution_clock::now();
        double total_s = duration<double>(t3 - t0).count();
        RAD_COUT << "[load_wl] finished loading all whitelists in "
                << total_s << " s\n";
    }

    // Save whitelist to a file
    void save_wl(std::ostream &out, bool verbose, bool full = true, const std::string& whitelist_type = "true") const {
    for (auto const& [class_id, wrap] : wl_map.maps) {
        auto &entry = wrap.get();
        if (verbose) {
            RAD_COUT << "[save_wl] [" << class_id << "] type=" << whitelist_type << std::endl;
        }
        
        if (whitelist_type == "true" || whitelist_type == "both") {
            std::unordered_set<int64_seq> seen_true;
            // Allocate only for the positive-count barcodes selected below.
            entry.true_bcs.for_each_entry([&](const auto& value) {
                if (value.count.load(barcode_counts::total) > 0) {
                    seen_true.emplace(value.key());
                }
            });

            if (!seen_true.empty()) {
                entry.true_bcs.write_wl_summary(out, class_id, &seen_true);
            }
        }
        
        if ((whitelist_type == "global" || whitelist_type == "both") && full) {
            std::unordered_set<int64_seq> seen_global;
            // The catalog can be huge even when few barcodes have counts.
            entry.global_bcs.for_each_entry([&](const auto& value) {
                if (value.count.load(barcode_counts::total) > 0) {
                    seen_global.emplace(value.key());
                }
            });
            
            if (!seen_global.empty()) {
                entry.global_bcs.write_wl_summary(out, class_id, &seen_global);
            }
        }
    }
}
    
    void save_wl(const std::string &path, bool verbose, bool full = true) const {
        if (!full) {
            std::ostringstream payload;
            save_wl(payload, verbose, false, "true");

            const std::string data = payload.str();
            if (data.empty()) {
                std::error_code ec;
                std::filesystem::remove(path, ec);
                if (verbose) {
                    RAD_COUT << "[save_wl] No true whitelist entries with counts; skipping " << path << std::endl;
                }
                return;
            }

            std::ofstream out(path);
            if (!out.is_open()) {
                throw std::runtime_error(
                    "Error opening whitelist summary for writing: " + path);
            }
            out.write(data.data(), static_cast<std::streamsize>(data.size()));
            out.flush();
            if (!out) {
                throw std::runtime_error(
                    "Error writing whitelist summary: " + path);
            }
            out.close();
            if (out.fail()) {
                throw std::runtime_error(
                    "Error closing whitelist summary: " + path);
            }
        } else {
            // Extract base and extension
            std::string base = path;
            size_t dot_pos = base.find_last_of('.');
            std::string extension = "";
            if (dot_pos != std::string::npos) {
                extension = base.substr(dot_pos);
                base = base.substr(0, dot_pos);
            }
            
            auto write_nonempty_payload = [&](const std::string& out_path, const std::string& wl_type) {
                std::ostringstream payload;
                save_wl(payload, verbose, true, wl_type);

                const std::string data = payload.str();
                if (data.empty()) {
                    std::error_code ec;
                    std::filesystem::remove(out_path, ec);
                    if (verbose) {
                        RAD_COUT << "[save_wl] No " << wl_type
                                  << " whitelist entries with counts; skipping " << out_path << std::endl;
                    }
                    return false;
                }

                std::ofstream out(out_path);
                if (!out.is_open()) {
                    throw std::runtime_error(
                        "Error opening whitelist summary for writing: " +
                        out_path);
                }
                out.write(data.data(),
                          static_cast<std::streamsize>(data.size()));
                out.flush();
                if (!out) {
                    throw std::runtime_error(
                        "Error writing whitelist summary: " + out_path);
                }
                out.close();
                if (out.fail()) {
                    throw std::runtime_error(
                        "Error closing whitelist summary: " + out_path);
                }
                return true;
            };

            const std::string true_path = base + "_true" + extension;
            const std::string global_path = base + "_global" + extension;
            bool wrote_true = write_nonempty_payload(true_path, "true");
            bool wrote_global = write_nonempty_payload(global_path, "global");

            if (verbose) {
                if (wrote_true || wrote_global) {
                    RAD_COUT << "[save_wl] Written "
                              << (wrote_true ? true_path : "[no true file]")
                              << " and "
                              << (wrote_global ? global_path : "[no global file]")
                              << std::endl;
                }
            }
        }
    }

    // print full layout
void display_read_layout(bool max_verbose = false) {
    if (max_verbose) {
        auto format_position = [](const ParsedPosition& pos) {
            return pos.ref_id.empty() ? std::string("[unset]") : pos.to_string();
        };
        auto format_optional_int = [](const std::optional<int>& value) {
            return value ? std::to_string(*value) : std::string("[unset]");
        };

        RAD_COUT << "[read_layout] Complete read layout (max verbose):\n";
        for (const auto& elem : by_dir_order()) {
            RAD_COUT << "[read_layout] Element " << elem.class_id << ":\n"
                      << "  Sequence: " << (elem.seq.empty() ? "[empty]" : elem.seq) << "\n"
                      << "  Masked sequence: "
                      << (elem.masked_seq.empty() ? "[empty]" : elem.masked_seq) << "\n"
                      << "  Type: " << elem.type << "\n"
                      << "  Global class: " << elem.global_class << "\n"
                      << "  Direction: " << elem.direction << "\n"
                      << "  Order: " << elem.order << "\n"
                      << "  Expected length: "
                      << format_optional_int(elem.expected_length) << "\n"
                      << "  Length candidates: "
                      << (elem.length_candidates.empty()
                              ? "[none]"
                              : format_length_candidates(elem.length_candidates))
                      << "\n"
                      << "  Whitelist: "
                      << (elem.whitelist_path.empty() ? "[none]" : elem.whitelist_path)
                      << "\n"
                      << "  Flags: " << (elem.flags.empty() ? "[none]" : elem.flags)
                      << "\n";

            if (elem.misalignment_threshold) {
                const auto [lower, mean, upper] = *elem.misalignment_threshold;
                RAD_COUT << "  Misalignment threshold: "
                          << lower << "|" << mean << "|" << upper << "\n";
            } else {
                RAD_COUT << "  Misalignment threshold: [unset]\n";
            }

            if (elem.aligned_positions) {
                const auto& aligned = *elem.aligned_positions;
                RAD_COUT << "  Aligned start (mean|variance): "
                          << aligned.start_stats.first << "|"
                          << aligned.start_stats.second << "\n"
                          << "  Aligned stop (mean|variance): "
                          << aligned.stop_stats.first << "|"
                          << aligned.stop_stats.second << "\n";
            } else {
                RAD_COUT << "  Aligned positions: [unset]\n";
            }

            if (elem.misaligned_positions) {
                const auto& misaligned = *elem.misaligned_positions;
                RAD_COUT << "  Misaligned start (mean|variance): "
                          << misaligned.start_stats.first << "|"
                          << misaligned.start_stats.second << "\n"
                          << "  Misaligned stop (mean|variance): "
                          << misaligned.stop_stats.first << "|"
                          << misaligned.stop_stats.second << "\n";
            } else {
                RAD_COUT << "  Misaligned positions: [unset]\n";
            }

            if (elem.ref_pos) {
                const auto& positions = *elem.ref_pos;
                RAD_COUT << "  Primary Start: "
                          << format_position(positions.primary_start) << "\n"
                          << "  Secondary Start: "
                          << format_position(positions.secondary_start) << "\n"
                          << "  Primary Stop: "
                          << format_position(positions.primary_stop) << "\n"
                          << "  Secondary Stop: "
                          << format_position(positions.secondary_stop) << "\n";
            } else {
                RAD_COUT << "  Position references: [unset]\n";
            }
        }
        RAD_COUT << "\n";
        return;
    }

    RAD_COUT << "[read_layout] Complete read layout:\n";

    auto colorize = [](const std::string& text, const std::string& color_code) {
        return term_print_utils::colorize(text, color_code);
    };
    auto class_color = [&](const ReadElement& elem) -> const std::string& {
        return term_print_utils::class_color(elem.class_id, elem.seq);
    };
    auto class_color_by_id = [&](const std::string& class_id) -> const std::string& {
        return term_print_utils::class_color(class_id);
    };
    const std::string color_dir = term_print_utils::direction_label_color();

    auto fmt_double = [](double v) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(2) << v;
        return oss.str();
    };
    auto fmt_ref = [&](const ParsedPosition& pos) {
        if (pos.ref_id.empty()) {
            return std::string("[unset]");
        }
        std::ostringstream oss;
        oss << term_print_utils::colorize(pos.ref_id, class_color_by_id(pos.ref_id))
            << "|" << (pos.is_start ? "start" : "stop");
        if (pos.offset != 0) {
            oss << (pos.offset > 0 ? "+" : "") << pos.offset;
        }
        return oss.str();
    };

        auto static_segment = [&](const ReadElement& elem) {
            std::string thresh = "NA";
            if (elem.misalignment_threshold) {
                auto [min_v, mean_v, max_v] = *elem.misalignment_threshold;
                thresh = fmt_double(static_cast<double>(mean_v));
            }

            std::string start = "NA";
            std::string stop = "NA";
            if (elem.aligned_positions) {
                const auto& ap = *elem.aligned_positions;
                start = fmt_double(ap.start_stats.first);
                stop = fmt_double(ap.stop_stats.first);
            }

            std::ostringstream oss;
            if(elem.global_class == "start" || elem.global_class == "stop" || elem.global_class == "poly_tail"){
                oss << "[" << term_print_utils::colorize(elem.class_id, class_color(elem)) << "]";
            } else {
                oss << "["  << term_print_utils::colorize(elem.class_id, class_color(elem))
                    << ":" << thresh << ":" << start << "->" << stop << "]";
            }
            return oss.str();
        };

        auto variable_segment = [&](const ReadElement& elem) {
            std::string start = "start?";
            std::string stop = "stop?";
            if (elem.ref_pos) {
                const auto& rp = *elem.ref_pos;
                start = fmt_ref(rp.primary_start);
                stop = fmt_ref(rp.primary_stop);
            }
            std::ostringstream oss;
            oss << "[" << term_print_utils::colorize(elem.class_id, class_color(elem))
                << ":" << start << "->" << stop << "]";
            return oss.str();
        };

        std::map<std::string, std::string> dir_lines;
        for (const auto& elem : by_dir_order()) {
            std::string segment = (elem.type == "static") ? static_segment(elem) : variable_segment(elem);
            dir_lines[elem.direction] += segment;
        }

        for (const auto& [dir, line] : dir_lines) {
            RAD_COUT << colorize(dir, color_dir) << ": " << line << "\n";
        }
        RAD_COUT << "\n";
    }

    // Other utility methods for managing layout
    size_t size() const { return layout.size(); }
    bool empty() const { return layout.empty(); }
    void clear() { layout.clear(); }

private:

    static std::string describe_debug_element(const ReadElement* elem) {
        if (elem == nullptr) return "[none]";

        std::ostringstream oss;
        oss << elem->class_id
            << " (order=" << elem->order
            << ", type=" << elem->type
            << ", class=" << elem->global_class
            << ", direction=" << elem->direction << ")";
        return oss.str();
    }

    static std::string describe_debug_position(const ParsedPosition& position) {
        return position.ref_id.empty() ? std::string("[unset]") : position.to_string();
    }

    void print_position_context(
        const std::string& phase,
        const ReadElement& elem,
        const ReadElement* prev,
        const ReadElement* next,
        const ReadElement* prev_prev,
        const ReadElement* next_next
    ) const {
        RAD_COUT << "\n[position_map] Processing " << phase
                  << " element: " << elem.class_id << "\n"
                  << "  Element: " << describe_debug_element(&elem) << "\n"
                  << "  Previous: " << describe_debug_element(prev) << "\n"
                  << "  Next: " << describe_debug_element(next) << "\n"
                  << "  Previous-previous: "
                  << describe_debug_element(prev_prev) << "\n"
                  << "  Next-next: " << describe_debug_element(next_next) << "\n";
    }

    void print_position_assignment(
        const std::string& phase,
        const std::string& class_id,
        const ReferencePositions& positions
    ) const {
        RAD_COUT << "[position_map] " << phase << " positions for "
                  << class_id << ":\n"
                  << "  Primary Start: "
                  << describe_debug_position(positions.primary_start) << "\n"
                  << "  Secondary Start: "
                  << describe_debug_position(positions.secondary_start) << "\n"
                  << "  Primary Stop: "
                  << describe_debug_position(positions.primary_stop) << "\n"
                  << "  Secondary Stop: "
                  << describe_debug_position(positions.secondary_stop) << "\n";
    }

    void process_left_side(
        Layout_Struct::index<dir_order_tag>::type::iterator it,
        bool verbose
    ) {
        auto& ordered = by_order();
        auto next = ordered.find(it->order + 1);
        auto next_next = ordered.find(it->order + 2);

        auto prev = ordered.find(it->order - 1);
        auto prev_prev = ordered.find(it->order - 2);
        if (prev == ordered.end()) return;

        if (verbose) {
            print_position_context(
                "left",
                *it,
                &*prev,
                next != ordered.end() ? &*next : nullptr,
                prev_prev != ordered.end() ? &*prev_prev : nullptr,
                next_next != ordered.end() ? &*next_next : nullptr
            );
        }

        ReferencePositions positions;
        int offset = it->expected_length.value_or(1);

        // Always set primary positions relative to previous element
        positions.primary_start = {prev->class_id, false, 1, ""};
        positions.primary_stop = {prev->class_id, false, offset, ""};

        // Handle secondary positions based on types
        if (prev->type == "static") {
            // Check for static anchor
            if (next != ordered.end() && next->type == "static" &&
                next->expected_length.value_or(0) >= 13 &&
                next->global_class != "poly_tail") {
                    if (verbose) {
                        RAD_COUT << "  Secondary selection: next static anchor "
                                  << next->class_id << "\n";
                    }
                    positions.secondary_start = {next->class_id, true, -offset, ""};
                    positions.secondary_stop = {next->class_id, true, -1, ""};
            } else {
                if (next != ordered.end() && next->type == "variable" &&
                    next_next != ordered.end() && next_next->type == "static" &&
                    next_next->global_class != "read" &&
                    next_next->global_class != "poly_tail" &&
                    next_next->global_class != "start" &&
                    next_next->global_class != "stop" &&
                    next_next->expected_length.value_or(0) >= 13) {
                    if (verbose) {
                        RAD_COUT << "  Secondary selection: next-next static anchor "
                                  << next_next->class_id << "\n";
                    }
                    int length_offset = next->expected_length.value_or(1);
                    positions.secondary_start = {next_next->class_id, true, -length_offset - offset, ""};
                    positions.secondary_stop = {next_next->class_id, true, -length_offset - 1, ""};
                } else {
                    if (verbose) {
                        RAD_COUT << "  Secondary selection: previous static fallback "
                                  << prev->class_id << "\n";
                    }
                    positions.secondary_start = {prev->class_id, false, 1, "left_terminal_linked"};
                    positions.secondary_stop = {prev->class_id, false, offset, "left_terminal_linked"};
                }
            }
        } else if (prev->type == "variable" || prev->global_class == "poly_tail") {
            std::string flag_type = prev->type == "variable" ? "var_chained" : "poly_chained";
            positions.primary_start.add_flags = flag_type + "_start";
            positions.primary_stop.add_flags = flag_type + "_stop";
            // Prefer a reliable static anchor; use sequence boundaries only as fallbacks.
            auto prev_static = ordered.find(it->order - 2);  // Look two positions back
            if (prev_static != ordered.end() && prev_static->type == "static" && 
                prev_static->global_class != "poly_tail" &&
                prev_static->global_class != "start" &&
                prev_static->global_class != "stop" &&
                prev_static->expected_length.value_or(0) >= 13) {
                if (verbose) {
                    RAD_COUT << "  Secondary selection: previous-previous static anchor "
                              << prev_static->class_id << "\n";
                }
                int prev_length = prev->expected_length.value_or(1);
                positions.secondary_start = {prev_static->class_id, false, 1 + prev_length, ""};
                positions.secondary_stop = {prev_static->class_id, false, offset + prev_length, ""};
            } else if (next != ordered.end() && next->type == "static" &&
                       next->global_class != "poly_tail" &&
                       next->global_class != "start" &&
                       next->global_class != "stop" &&
                       next->expected_length.value_or(0) >= 13) {
                if (verbose) {
                    RAD_COUT << "  Secondary selection: next static anchor "
                              << next->class_id << "\n";
                }
                positions.secondary_start = {next->class_id, true, -offset, ""};
                positions.secondary_stop = {next->class_id, true, -1, ""};
            } else if (prev_static != ordered.end() && prev_static->type == "static" &&
                       (prev_static->global_class == "start" ||
                        prev_static->global_class == "stop")) {
                if (verbose) {
                    RAD_COUT << "  Secondary selection: previous-previous terminal fallback "
                              << prev_static->class_id << "\n";
                }
                int prev_length = prev->expected_length.value_or(1);
                positions.secondary_start = {prev_static->class_id, false, 1 + prev_length, ""};
                positions.secondary_stop = {prev_static->class_id, false, offset + prev_length, ""};
            } else if (next != ordered.end() && next->type == "static" &&
                       (next->global_class == "start" ||
                        next->global_class == "stop")) {
                if (verbose) {
                    RAD_COUT << "  Secondary selection: next terminal fallback "
                              << next->class_id << "\n";
                }
                positions.secondary_start = {next->class_id, true, -offset, ""};
                positions.secondary_stop = {next->class_id, true, -1, ""};
            } else if (auto prev_pos = position_map.find(prev->class_id); prev_pos != position_map.end()) {
                // Finally fall back to chaining
                if (verbose) {
                    RAD_COUT << "  Secondary selection: chained through previous "
                              << prev->class_id << "\n";
                }
                int prev_length = prev->expected_length.value_or(1);
                positions.secondary_start = {prev_pos->second.primary_start.ref_id, false, 1 + prev_length, "chained_start"};
                positions.secondary_stop = {prev_pos->second.primary_stop.ref_id, false, 1 + prev_length + offset, "chained_stop"};
            }
        }

        auto is_terminal = [](const ReadElement& elem) {
            return elem.global_class == "start" || elem.global_class == "stop";
        };
        auto is_reliable_static = [](const ReadElement& elem) {
            return elem.type == "static" &&
                   elem.global_class != "poly_tail" &&
                   elem.global_class != "start" &&
                   elem.global_class != "stop" &&
                   elem.expected_length.value_or(0) >= 13;
        };
        const bool primary_is_terminal_linked =
            is_terminal(*prev) ||
            (prev->type == "variable" && prev_prev != ordered.end() &&
             is_terminal(*prev_prev));
        auto& id_index = by_id();
        auto secondary = id_index.find(positions.secondary_start.ref_id);
        if (primary_is_terminal_linked &&
            positions.secondary_start.ref_id == positions.secondary_stop.ref_id &&
            secondary != id_index.end() && is_reliable_static(*secondary)) {
            if (verbose) {
                RAD_COUT << "  Primary selection: promoting static anchor "
                          << secondary->class_id
                          << " over terminal-linked chain\n";
            }
            std::swap(positions.primary_start, positions.secondary_start);
            std::swap(positions.primary_stop, positions.secondary_stop);
        }

        position_map[it->class_id] = positions;
        if (verbose) {
            print_position_assignment("Parse Left", it->class_id, positions);
        }
    }

    void process_right_side(
        Layout_Struct::index<dir_order_tag>::type::iterator it,
        bool verbose
    ) {
        auto& ordered = by_order();
        auto next = ordered.find(it->order + 1);

        auto prev = ordered.find(it->order - 1);
        auto prev_prev = ordered.find(it->order - 2);
        auto next_next = ordered.find(it->order + 2);

        if (next == ordered.end()) return;

        if (verbose) {
            print_position_context(
                "right",
                *it,
                prev != ordered.end() ? &*prev : nullptr,
                &*next,
                prev_prev != ordered.end() ? &*prev_prev : nullptr,
                next_next != ordered.end() ? &*next_next : nullptr
            );
        }

        ReferencePositions positions;
        int offset = it->expected_length.value_or(1);

        // Always set primary positions relative to next element
        positions.primary_start = {next->class_id, true, -offset, ""};
        positions.primary_stop = {next->class_id, true, -1, ""};

        // Handle secondary positions based on types
        if (next->type == "static") {
            // Check for static anchor
            if (prev != ordered.end() && prev->type == "static" &&
                prev->expected_length.value_or(0) >= 13 &&
                prev->global_class != "poly_tail") {
                    if (verbose) {
                        RAD_COUT << "  Secondary selection: previous static anchor "
                                  << prev->class_id << "\n";
                    }
                    positions.secondary_start = {prev->class_id, false, 1, ""};
                    positions.secondary_stop = {prev->class_id, false, offset, ""};
            } else {

                if (prev != ordered.end() && prev->type == "variable" &&
                    prev_prev != ordered.end() && prev_prev->type == "static" &&
                    prev_prev->global_class != "read" &&
                    prev_prev->global_class != "poly_tail" &&
                    prev_prev->global_class != "start" &&
                    prev_prev->global_class != "stop" &&
                    prev_prev->expected_length.value_or(0) >= 13) {
                    if (verbose) {
                        RAD_COUT << "  Secondary selection: previous-previous static anchor "
                                  << prev_prev->class_id << "\n";
                    }
                    int length_offset = prev->expected_length.value_or(1);
                    positions.secondary_start = {prev_prev->class_id, false, length_offset + 1, ""};
                    positions.secondary_stop = {prev_prev->class_id, false, length_offset + offset, ""};
                } else {
                    if (verbose) {
                        RAD_COUT << "  Secondary selection: next static fallback "
                                  << next->class_id << "\n";
                    }
                    positions.secondary_start = {next->class_id, true, -offset, "right_terminal_linked"};
                    positions.secondary_stop = {next->class_id, true, -1, "right_terminal_linked"};
                }
            }
        } else if (next->type == "variable" || next->global_class == "poly_tail") {
            std::string flag_type = prev->type == "variable" ? "var_chained" : "poly_chained";
            positions.primary_start.add_flags = flag_type + "_start";
            positions.primary_stop.add_flags = flag_type + "_stop";

            // Prefer a reliable static anchor; use sequence boundaries only as fallbacks.
            auto next_static = ordered.find(it->order + 2);  // Look two positions ahead
            if (next_static != ordered.end() && next_static->type == "static" && 
                next_static->global_class != "poly_tail" &&
                next_static->global_class != "start" &&
                next_static->global_class != "stop" &&
                next_static->expected_length.value_or(0) >= 13) {
                if (verbose) {
                    RAD_COUT << "  Secondary selection: next-next static anchor "
                              << next_static->class_id << "\n";
                }
                int next_length = next->expected_length.value_or(1);
                positions.secondary_start = {next_static->class_id, true, -next_length - offset, ""};
                positions.secondary_stop = {next_static->class_id, true, -next_length - 1, ""};
            } else if (prev != ordered.end() && prev->type == "static" &&
                       prev->global_class != "poly_tail" &&
                       prev->global_class != "start" &&
                       prev->global_class != "stop" &&
                       prev->expected_length.value_or(0) >= 13) {
                if (verbose) {
                    RAD_COUT << "  Secondary selection: previous static anchor "
                              << prev->class_id << "\n";
                }
                positions.secondary_start = {prev->class_id, false, 1, ""};
                positions.secondary_stop = {prev->class_id, false, offset, ""};
            } else if (next_static != ordered.end() && next_static->type == "static" &&
                       (next_static->global_class == "start" ||
                        next_static->global_class == "stop")) {
                if (verbose) {
                    RAD_COUT << "  Secondary selection: next-next terminal fallback "
                              << next_static->class_id << "\n";
                }
                int next_length = next->expected_length.value_or(1);
                positions.secondary_start = {next_static->class_id, true, -next_length - offset, ""};
                positions.secondary_stop = {next_static->class_id, true, -next_length - 1, ""};
            } else if (prev != ordered.end() && prev->type == "static" &&
                       (prev->global_class == "start" ||
                        prev->global_class == "stop")) {
                if (verbose) {
                    RAD_COUT << "  Secondary selection: previous terminal fallback "
                              << prev->class_id << "\n";
                }
                positions.secondary_start = {prev->class_id, false, 1, ""};
                positions.secondary_stop = {prev->class_id, false, offset, ""};
            } else if (auto next_pos = position_map.find(next->class_id); next_pos != position_map.end()) {
                // Finally fall back to chaining
                if (verbose) {
                    RAD_COUT << "  Secondary selection: chained through next "
                              << next->class_id << "\n";
                }
                int next_length = next->expected_length.value_or(1);
                positions.secondary_start = {next_pos->second.primary_start.ref_id, true, -offset - next_length, "chained_start"};
                positions.secondary_stop = {next_pos->second.primary_stop.ref_id, true, -next_length, "chained_stop"};
            }
        }

        auto is_terminal = [](const ReadElement& elem) {
            return elem.global_class == "start" || elem.global_class == "stop";
        };
        auto is_reliable_static = [](const ReadElement& elem) {
            return elem.type == "static" &&
                   elem.global_class != "poly_tail" &&
                   elem.global_class != "start" &&
                   elem.global_class != "stop" &&
                   elem.expected_length.value_or(0) >= 13;
        };
        const bool primary_is_terminal_linked =
            is_terminal(*next) ||
            (next->type == "variable" && next_next != ordered.end() &&
             is_terminal(*next_next));
        auto& id_index = by_id();
        auto secondary = id_index.find(positions.secondary_start.ref_id);
        if (primary_is_terminal_linked &&
            positions.secondary_start.ref_id == positions.secondary_stop.ref_id &&
            secondary != id_index.end() && is_reliable_static(*secondary)) {
            if (verbose) {
                RAD_COUT << "  Primary selection: promoting static anchor "
                          << secondary->class_id
                          << " over terminal-linked chain\n";
            }
            std::swap(positions.primary_start, positions.secondary_start);
            std::swap(positions.primary_stop, positions.secondary_stop);
        }

        position_map[it->class_id] = positions;
        if (verbose) {
            print_position_assignment("Parse Right", it->class_id, positions);
        }
    }

    void process_reads(
        Layout_Struct::index<dir_order_tag>::type::iterator start,
        Layout_Struct::index<dir_order_tag>::type::iterator end,
        bool verbose
    ) {
        auto& ordered = by_order();
        // First find the read in this direction
        auto read_it = start;
        for (; read_it != end; ++read_it) {
            if (read_it->type == "variable" && read_it->global_class == "read") {
                break;
            }
        }
        if (read_it == end) return;
        if (verbose) {
            RAD_COUT << "\n[position_map] Processing direction "
                      << read_it->direction << " around read "
                      << read_it->class_id << "\n";
        }
        // Process elements before the read with left-side anchoring
        for (auto it = start; it != read_it; ++it) {
            if (it->type == "variable") {
                process_left_side(it, verbose);
            }
        }
        // Process read positions
        auto prev = ordered.find(read_it->order - 1);
        if (prev != ordered.end()) {
            find_read_positions(*read_it, *prev, verbose);
        }
        // Process elements after the read with right-side anchoring
        auto after_read = read_it;
        ++after_read;
        for (auto it = after_read; it != end; ++it) {
            if (it->type == "variable") {
                process_right_side(it, verbose);
            }
        }
    }
    
    void find_read_positions(
        const ReadElement& read,
        const ReadElement& prev,
        bool verbose
    ) {
        auto& ordered = by_order();
        ReferencePositions positions;
        positions.primary_start = {prev.class_id, false, 1, ""};

        auto next_it = ordered.find(read.order + 1);
        auto next_next_it = ordered.find(read.order + 2);
        auto pre_prev_it = ordered.find(read.order - 2);

        if (verbose) {
            print_position_context(
                "read",
                read,
                &prev,
                next_it != ordered.end() ? &*next_it : nullptr,
                pre_prev_it != ordered.end() ? &*pre_prev_it : nullptr,
                next_next_it != ordered.end() ? &*next_next_it : nullptr
            );
        }

        if (next_it != ordered.end()) {
            positions.primary_stop = {next_it->class_id, true, -1, ""};
            if (next_next_it != ordered.end()) {
                std::string flag = next_it->type == "static" ? "truncated" : "skipped";
                int offset = next_it->expected_length.value_or(1);
                positions.secondary_stop = {
                    next_next_it->class_id,
                    true,
                    flag == "skipped" ? -(1 + offset) : -1,
                    "right_" + flag
                };
            }
        }

        if (pre_prev_it != ordered.end()) {
            std::string flag = prev.type == "static" ? "truncated" : "skipped";
            int offset = prev.expected_length.value_or(1);
            positions.secondary_start = {
                pre_prev_it->class_id,
                false,
                flag == "skipped" ? 1 + offset : 1,
                "left_" + flag
            };
        }

        if (positions.secondary_start.ref_id.empty()) {
            positions.secondary_start = positions.primary_start;
        }
        if (positions.secondary_stop.ref_id.empty()) {
            positions.secondary_stop = positions.primary_stop;
        }

        position_map[read.class_id] = positions;
        if (verbose) {
            print_position_assignment("Read", read.class_id, positions);
        }
    }

    const auto& get_pos_map() const { 
        return position_map; 
    }

};

/**
 * @brief Concat-HMM calibration (hmm_* columns). Workers get a fixed strided partition of each chunk
 * and are merged in worker order, so the result does not depend on scheduling.
 */
class concat_hmm_cal {
public:
    concat_hmm_cal(const ReadLayout& layout, int num_threads)
        : prior_(concat_hmm::Model::build(layout.to_concat_layout_spec())),
          master_(prior_),
          nthreads_(std::max(1, std::min(num_threads, 16))) {}

    void add_chunk(const std::vector<read_streaming::sequence>& chunk, size_t max_reads,
                   const std::function<void()>& boundary_poll = {}) {
        std::lock_guard<std::mutex> lock(mu_);
        const size_t room = max_reads > added_ ? max_reads - added_ : 0;
        const size_t n = std::min(chunk.size(), room);
        if (n == 0) return;
        if (boundary_poll) boundary_poll();
        const int workers = static_cast<int>(std::max<size_t>(1, std::min<size_t>(nthreads_, n / 1000)));
        std::vector<concat_hmm::Calibrator> cals;
        cals.reserve(workers);
        for (int w = 0; w < workers; ++w) cals.emplace_back(prior_);
        if (workers == 1) {
            for (size_t r = 0; r < n; ++r) {
                if (boundary_poll && r % 32 == 0) boundary_poll();
                cals[0].add_read(chunk[r].seq.data(), checked_length(chunk[r].seq.size()));
            }
        } else {
            std::atomic<bool> cancelled{false};
            std::atomic<int> remaining{workers};
            std::mutex state_mutex;
            std::condition_variable ready;
            std::exception_ptr worker_error;
            std::vector<std::thread> pool;
            pool.reserve(workers);
            auto join_all = [&]() {
                for (auto& worker : pool) if (worker.joinable()) worker.join();
            };
            try {
                for (int w = 0; w < workers; ++w) {
                    pool.emplace_back([&, w]() {
                        try {
                            for (size_t r = static_cast<size_t>(w); r < n;
                                 r += static_cast<size_t>(workers)) {
                                if (cancelled.load(std::memory_order_acquire)) break;
                                cals[w].add_read(chunk[r].seq.data(), checked_length(chunk[r].seq.size()));
                            }
                        } catch (...) {
                            std::lock_guard<std::mutex> state_lock(state_mutex);
                            if (!worker_error) worker_error = std::current_exception();
                            cancelled.store(true, std::memory_order_release);
                        }
                        remaining.fetch_sub(1, std::memory_order_release);
                        ready.notify_one();
                    });
                }
                while (remaining.load(std::memory_order_acquire) != 0) {
                    if (boundary_poll) boundary_poll();
                    std::unique_lock<std::mutex> state_lock(state_mutex);
                    ready.wait_for(state_lock, std::chrono::milliseconds(50), [&]() {
                        return remaining.load(std::memory_order_acquire) == 0;
                    });
                }
            } catch (...) {
                // Includes host cancellation and failure to create a later
                // thread: no joinable thread may escape stack unwinding.
                cancelled.store(true, std::memory_order_release);
                join_all();
                throw;
            }
            join_all();
            if (worker_error) std::rethrow_exception(worker_error);
        }
        for (const auto& cal : cals) {
            if (boundary_poll) boundary_poll();
            master_.merge(cal);
        }
        added_ += n;
    }

    // Finalizes (hard-EM + calibration guard), prints the report and stores the parameters.
    bool finish(ReadLayout& layout, const std::function<void()>& boundary_poll = {}) {
        std::lock_guard<std::mutex> lock(mu_);
        if (master_.n_reads() == 0) {
            RAD_CERR << "[concat_hmm] WARNING: no reads for calibration; layout priors will be used\n";
            return false;
        }
        const auto t0 = std::chrono::steady_clock::now();
        concat_hmm::Params params = master_.finalize(20, 0.001, boundary_poll);
        const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        RAD_COUT << "[concat_hmm] Calibration on " << master_.n_reads() << " reads (EM " << secs << " s):\n"
                  << master_.report();
        if (!master_.layout_matches())
            RAD_COUT << "[concat_hmm] WARNING: calibration guard failed (the layout does not match these reads); "
                         "the HMM is disabled and every read uses the existing path\n";
        layout.concat_params = std::move(params);
        return master_.layout_matches();
    }

    // Cached runs whose position map has no hmm_* columns: calibrate from the FASTQ head
    // without redoing the misalignment step.
    static bool calibrate_from_fastq(const std::string& fastq_path, ReadLayout& layout, int num_threads,
                                     size_t max_reads = 50000,
                                     const std::function<void()>& boundary_poll = {}) {
        // Input corruption, resource failures, and host cancellation must
        // propagate. A model/layout mismatch is handled by build_concat_model.
        if (max_reads == 0 || max_reads > static_cast<size_t>(std::numeric_limits<int64_t>::max()))
            throw std::invalid_argument("Invalid HMM calibration read count");
        std::unique_ptr<concat_hmm_cal> cal;
        try { cal = std::make_unique<concat_hmm_cal>(layout, num_threads); }
        catch (const concat_hmm::layout_spec_error& error) {
            RAD_CERR << "[concat_hmm] WARNING: cannot calibrate this layout (" << error.what() << ")\n";
            return false;
        }
        using ChunkFunc = std::function<void(const std::vector<read_streaming::sequence>&, const std::string&)>;
        ChunkFunc feed = [&](const std::vector<read_streaming::sequence>& chunk, const std::string&) {
            cal->add_chunk(chunk, max_reads, boundary_poll);
        };
        chunk_streaming<read_streaming::sequence, ChunkFunc> streamer(std::min<size_t>(max_reads, 50000));
        streamer.process_chunks(fastq_path, feed, 1, static_cast<int64_t>(max_reads), false, 0, boundary_poll);
        return cal->finish(layout, boundary_poll);
    }

private:
    static int checked_length(size_t size) {
        if (size > static_cast<size_t>(std::numeric_limits<int>::max()))
            throw std::length_error("Read is too long for concat-HMM calibration");
        return static_cast<int>(size);
    }
    concat_hmm::Model prior_;
    concat_hmm::Calibrator master_;
    int nthreads_;
    size_t added_ = 0;
    std::mutex mu_;
};

/**
 * @brief This class contains all of the setup for the misalignment threshold, as well as the masked adapters.
 * @brief contains structures for perfect matches, misalignment statistics, and position statistics.
 * @struct perfect_match
 * @struct misalignment_stats
 * @struct PositionStats
 */
class Misalignment_Setup {
public:
    inline static const std::array<EdlibEqualityPair, 8> kWildcardEqualities{{
        {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'},
        {'A', 'N'}, {'C', 'N'}, {'G', 'N'}, {'T', 'N'}
    }};
    /**
     * @brief information on perfect matches
     * @param read read sequence
     * @param direction The direction of the read ("forward" or "reverse")
     */
    struct perfect_match {
        read_streaming::sequence read;
        std::string direction;
    };

    /**
     * @brief structure for containing statistics of alignment per adapter
     * @param mean `double` mean edit distance
     * @param sum_squares  `double` sum of squares for standard deviation calculation
     * @param count  `size_t` number of reads considered
     * @param total_reads_seen `size_t` total number of reads processed
     * @param is_stable `bool` whether the statistics have stabilized
     * @param dir_counts `pair<size_t, size_t>` of counts for forward and reverse directions
     */
    struct static_alignment_stats {
        double mean;
        double sum_squares;
        size_t count;
        size_t total_reads_seen;
        bool is_stable;
        std::pair<size_t, size_t> dir_counts;
        std::vector<size_t> ed_hist;  // accrued values per edit distance (index = edit distance)

    /**
     * @brief nested structure for storing position statistics
     * @param start_positions vector of start positions
     * @param stop_positions vector of stop positions
     */
    struct PositionStats {
        std::vector<double> start_positions;
        std::vector<double> stop_positions;
        /**
         * @brief helper function to calculate mean and standard deviation
         * @param positions vector of positions
         * @return pair of mean and standard deviation
         */
        std::pair<double, double> get_stats(const std::vector<double>& positions) {
            if (positions.empty()) return {0.0, 0.0};
            
            // Calculate mean of normalized positions
            double mean = std::accumulate(positions.begin(), positions.end(), 0.0) / positions.size();
            
            // Calculate standard deviation
            double sum_squares = std::accumulate(positions.begin(), positions.end(), 0.0,
                [mean](double sum, double val) {
                    double diff = val - mean;
                    return sum + (diff * diff);
                });
                
            double sd = std::sqrt(sum_squares / (positions.size() - 1));
            return {mean, sd};
        }
        /**
         * @brief calculate mean and standard deviation for start and stop positions
         * @return AlignmentPositions containing statistics for start and stop positions
         */
        AlignmentPositions calculate() {
            return {
                get_stats(start_positions), 
                get_stats(stop_positions)
            };
        }
    } position_stats;

    static_alignment_stats() : mean(0), sum_squares(0), count(0), total_reads_seen(0), is_stable(false), dir_counts({0, 0}) {}

    //contains defaults for the misalignment stats
    /**
     * @brief update statistics with a new edit distance
     * @param new_edit_distance the new edit distance to incorporate
     * @param total_reads total number of reads processed
     * @param forward_count count of forward reads
     * @param reverse_count count of reverse reads
     * @param direction direction of the read ("forward" or "reverse")
     */
    void update_stats(int new_edit_distance, size_t total_reads, size_t forward_count, size_t reverse_count) {
            double old_mean = mean;
            count++;
            if (new_edit_distance >= 0) {
                if (ed_hist.size() <= static_cast<size_t>(new_edit_distance))
                    ed_hist.resize(static_cast<size_t>(new_edit_distance) + 1, 0);
                ed_hist[static_cast<size_t>(new_edit_distance)]++;
            }
            total_reads_seen = total_reads;
            mean = mean + (new_edit_distance - mean) / count;
            sum_squares += (new_edit_distance - old_mean) * (new_edit_distance - mean);
            double sample_ratio = static_cast<double>(count) / total_reads_seen;
            is_stable = ((std::abs(mean - old_mean) < 0.5) &&
                         (sample_ratio > 0.01)) || count >= 10000;
            dir_counts = {forward_count, reverse_count};
        }

    /**
     * @brief calculate standard deviation of edit distances
     * @return standard deviation
     */
    double get_sd() const {
            if (count < 2) return 0.0;
            return std::sqrt(sum_squares / (count - 1));
        }

    /**
     * @brief calculate confidence score based on count
     * @return confidence score
     */
    double confidence() const {
            return 1.0 - std::exp(-static_cast<double>(count) / 100.0);
        }
    };

    /**
     * @brief structure for consensus matrix
     * @param position_counts vector of maps counting bases at each position
     * @param max_edit_distance maximum edit distance for sequences in this matrix
     */
    struct consensus_matrix {
        std::vector<std::map<char, int>> position_counts;  // For each position, count of each base
        int max_edit_distance;  // Maximum edit distance for sequences in this matrix
        consensus_matrix(size_t length, int max_ed)
            : position_counts(length), max_edit_distance(max_ed) {}

        void ensure_length(size_t needed_length) {
            if (needed_length > position_counts.size()) {
                position_counts.resize(needed_length + std::max(needed_length/2, size_t(10)));
            }
        }
    };

    /**
     * @brief Chance hits of one element's misalignment edit distances. Values at or below the lowest bin
     * under the mode (v) are true copies when a higher bin lies below v; all other values are chance hits.
     */
    struct chance_hit_stats {
        size_t n_all = 0;              // accrued values
        size_t n_low = 0;              // values in the low group (true copies), removed
        size_t n_chance = 0;           // values in the chance-hit group
        int mode = -1;                 // main peak used for the split
        int split = -1;                // last edit distance of the low group; -1 when there is no low group
        double mean = 0.0;             // mean of the chance hits
        double sd = 0.0;               // sample standard deviation of the chance hits
        int lower = -1;                // largest edit distance with a cumulative chance-hit share <= chance_share
        double share_at_lower = 0.0;   // cumulative chance-hit share at lower
        double share_above_lower = 0.0;// cumulative chance-hit share at lower + 1
    };

    static chance_hit_stats chance_hits_from_histogram(const std::vector<size_t>& hist, double chance_share) {
        chance_hit_stats r;
        for (size_t c : hist) r.n_all += c;
        if (hist.empty() || r.n_all == 0) return r;
        const int last = static_cast<int>(hist.size()) - 1;
        auto argmax_from = [&](int from) {
            int m = -1;
            size_t best = 0;
            for (int d = from; d <= last; ++d)
                if (hist[d] > 0 && hist[d] >= best) { best = hist[d]; m = d; }
            return m;
        };
        int mode = argmax_from(0);
        if (mode == 0) {
            const int m1 = argmax_from(1);
            if (m1 > 0) mode = m1;
        }
        r.mode = mode;
        int v = 0;
        size_t vmin = hist[0];
        for (int d = 0; d <= mode; ++d)
            if (hist[d] <= vmin) { vmin = hist[d]; v = d; }
        size_t below_max = 0;
        for (int d = 0; d < v; ++d) below_max = std::max(below_max, hist[d]);
        const int first = (below_max > vmin) ? v + 1 : 0;   // first edit distance of the chance-hit group
        r.split = first - 1;
        double sum = 0.0;
        for (int d = 0; d <= last; ++d) {
            if (d < first) { r.n_low += hist[d]; continue; }
            r.n_chance += hist[d];
            sum += static_cast<double>(d) * static_cast<double>(hist[d]);
        }
        if (r.n_chance == 0) return r;
        r.mean = sum / static_cast<double>(r.n_chance);
        double ss = 0.0;
        for (int d = first; d <= last; ++d) {
            const double diff = static_cast<double>(d) - r.mean;
            ss += diff * diff * static_cast<double>(hist[d]);
        }
        r.sd = r.n_chance > 1 ? std::sqrt(ss / static_cast<double>(r.n_chance - 1)) : 0.0;
        auto share_at = [&](int d) {
            size_t cum = 0;
            for (int e = first; e <= std::min(d, last); ++e) cum += hist[e];
            return static_cast<double>(cum) / static_cast<double>(r.n_chance);
        };
        r.lower = -1;
        for (int d = 0; d <= last; ++d) {
            if (share_at(d) <= chance_share) r.lower = d;
            else break;
        }
        if (r.lower < 0) r.lower = 0;   // chance hits at edit distance 0 above the limit: nothing but exact copies
        r.share_at_lower = share_at(r.lower);
        r.share_above_lower = share_at(r.lower + 1);
        return r;
    }

/**
 * @brief Constructor that initializes the Misalignment_Setup with a given ReadLayout
 * @param layout ReadLayout object
 * @param misalign_lower_chance_share accepted cumulative share of chance hits at misalign_lower (default
 *        adapter_thresholds::misalign_lower_chance_share)
 */
    Misalignment_Setup(const ReadLayout& layout,
                       double misalign_lower_chance_share = adapter_thresholds::misalign_lower_chance_share)
        : misalign_lower_chance_share_(misalign_lower_chance_share) {
        if (!std::isfinite(misalign_lower_chance_share_) || misalign_lower_chance_share_ < 0 ||
            misalign_lower_chance_share_ > 1)
            throw std::invalid_argument("Misalignment chance share must be between zero and one");
        // Extract static adapters (non poly-tail, non start/stop)
        // guards against empty seq -- need to update for static indexing
        for (const auto& elem : layout.by_type()) {
            if (elem.type == "static" &&
                elem.global_class != "poly_tail" &&
                elem.global_class != "start" &&
                elem.global_class != "stop" &&
                elem.seq != ""
            ) {
                if (elem.direction == "forward") {
                    forward_adapters.push_back({elem.class_id, elem.seq});
                    init_consensus_matrices(elem.class_id, elem.seq.length(), 20);
                } else {
                    reverse_adapters.push_back({elem.class_id, elem.seq});
                    init_consensus_matrices(elem.class_id, elem.seq.length(), 20);
                }
            }
        }
    }

    /**
     * @brief Accrue the misalignment null from the FASTQ head (sample rule: adapter_thresholds) and update the layout.
     * @param min_reads size of the first block (rad prep: -n)
     * @param calibrate_concat_hmm also calibrate the concat-HMM, on the first min_reads reads only
     */
    void generate_misalignment_data(
        const std::string& fastq_path, ReadLayout& layout, int num_threads = 1,
        size_t min_reads = adapter_thresholds::cal_block_reads,
        std::function<void()> boundary_poll = {},
        bool calibrate_concat_hmm = false
    ) {
        if (min_reads == 0) throw std::invalid_argument("Calibration read count must be positive");
        if (boundary_poll) boundary_poll();
        prune_similar_reverse_adapters(layout);

        size_t forward_count = 0;
        size_t reverse_count = 0;
        size_t total_reads_processed = 0;

        //collecting statistics of both adapter stats and misalignment stats
        std::unordered_map<std::string, static_alignment_stats> adapter_stats;
        std::unordered_map<std::string, static_alignment_stats> misalignment_stats;
        auto process_block = [&](const std::vector<read_streaming::sequence>& chunk) {
            total_reads_processed += chunk.size();
            process_misalignment_chunk(
                chunk,
                forward_count,
                reverse_count,
                adapter_stats,
                misalignment_stats,
                total_reads_processed, boundary_poll
            );
            // Print progress and current misalignment thresholds.
            RAD_COUT <<
            "\r[misalignment_stats] Processed " << total_reads_processed <<
            " total reads - Forward perfect: " <<  forward_count  <<
            ", Reverse perfect: " << reverse_count << "\n";

            for (const auto& [adapter_id, stats] : adapter_stats) {
                RAD_COUT <<
                "[misalignment_stats] " << adapter_id <<
                " misalignment mean: " <<  misalignment_stats[adapter_id].mean  <<
                " (total times detected properly = " << stats.count << "," <<
                " misaligned forward: " << misalignment_stats[adapter_id].dir_counts.first <<
                ", misaligned reverse: " << misalignment_stats[adapter_id].dir_counts.second <<
                ")" <<
                (stats.is_stable ? " [STABLE]" : "") <<
                "\n";
            }
        };

        // chance hits of every element that a perfect read of the other direction feeds
        auto chance_counts = [&]() {
            std::vector<std::pair<std::string, size_t>> n;
            for (const auto* adapters : {&forward_adapters, &reverse_adapters})
                for (const auto& [id, seq] : *adapters) {
                    auto it = misalignment_stats.find(id);
                    n.emplace_back(id, it == misalignment_stats.end() ? size_t{0}
                        : chance_hits_from_histogram(it->second.ed_hist, misalign_lower_chance_share_).n_chance);
                }
            return n;
        };

        const size_t block = adapter_thresholds::cal_block_reads;
        const size_t cap = std::max(adapter_thresholds::cal_max_reads, min_reads);
        const size_t need = adapter_thresholds::min_misal_obs;
        file_streaming files(fastq_path, 4);
        read_streaming reader(files);
        size_t target = min_reads;
        size_t blocks = 0;
        bool end_of_input = false;
        std::string stop_reason;
        std::vector<std::pair<std::string, size_t>> n_chance;
        while (true) {
            // At most one block in memory; only the sequence is kept.
            while (total_reads_processed < target && !end_of_input) {
                std::vector<read_streaming::sequence> chunk;
                const size_t want = std::min(block, target - total_reads_processed);
                chunk.reserve(want);
                while (chunk.size() < want) {
                    if (boundary_poll && chunk.size() % 1024 == 0) boundary_poll();
                    auto rec = reader.next_sequence();
                    if (!rec) { end_of_input = true; break; }
                    if (!keep_perf_matches) { rec->qual = std::string(); rec->comment = std::string(); }
                    chunk.push_back(std::move(*rec));
                }
                if (!chunk.empty()) {
                    process_block(chunk);
                    ++blocks;
                }
            }
            if (reader.read_failed()) {
                std::string message = "sequence read failed (kseq code " +
                    std::to_string(reader.read_error_code()) + ") in " + reader.read_error_path();
                if (!reader.read_error_detail().empty())
                    message += " (" + reader.read_error_detail() + ")";
                throw std::runtime_error(message);
            }
            if (boundary_poll) boundary_poll();
            n_chance = chance_counts();
            if (n_chance.empty()) { stop_reason = "no static element to calibrate"; break; }
            size_t slowest = std::numeric_limits<size_t>::max();
            for (const auto& [id, n] : n_chance) slowest = std::min(slowest, n);
            if (slowest >= need) { stop_reason = "every element has enough chance hits"; break; }
            if (end_of_input) { stop_reason = "end of input"; break; }
            if (total_reads_processed >= cap) { stop_reason = "read limit reached"; break; }
            const double projected = slowest == 0 ? std::numeric_limits<double>::infinity()
                : static_cast<double>(total_reads_processed) * static_cast<double>(need) / static_cast<double>(slowest);
            if (projected > static_cast<double>(cap)) {
                stop_reason = slowest == 0 ? "an element has no chance hits"
                                           : "projected sample above the read limit";
                break;
            }
            target = std::min(total_reads_processed + block, cap);
        }
        RAD_COUT << "[misalignment_stats] Calibration sample: " << total_reads_processed << " reads in "
                  << blocks << " block(s); stop: " << stop_reason << " (need " << need
                  << " chance hits per element, read limit " << cap << "); chance hits:";
        for (const auto& [id, n] : n_chance) RAD_COUT << " " << id << "=" << n;
        RAD_COUT << "\n";

        
        // Write final results
        write_perfect_matches(keep_perf_matches);
        update_read_layout(layout, adapter_stats, misalignment_stats);
        layout.generate_position_mapping();
        // Adapter calibration can remove unobserved layout elements. Build the
        // HMM only after that final topology is known; otherwise its freshly
        // learned structural tables would be discarded as a cache mismatch.
        // A bounded second pass avoids retaining the raw sample in memory.
        if (calibrate_concat_hmm) {
            concat_hmm_cal::calibrate_from_fastq(fastq_path, layout, num_threads,
                                                 min_reads, boundary_poll);
        }
    }

private:
    double misalign_lower_chance_share_; /// accepted cumulative share of chance hits at misalign_lower
    std::unordered_map<std::string, std::vector<consensus_matrix>> consensus_matrices;
    std::vector<std::pair<std::string, std::string>> forward_adapters; /// vector of forward adapters, stored as `pair:[class_id, seq]`
    std::vector<std::pair<std::string, std::string>> reverse_adapters; /// vector of reverse adapters, stored as `pair:[class_id, seq]`
    // Debug only: keep perfect reads for write_perfect_matches (large memory).
    static constexpr bool keep_perf_matches = false;
    std::vector<perfect_match> perfect_matches; /// vector of perfect matches (only with keep_perf_matches)

    double alignment_distance_fraction(const std::string& lhs, const std::string& rhs) const {
        if (lhs.empty() && rhs.empty()) return 0.0;
        if (lhs.empty() || rhs.empty()) return 1.0;

        EdlibAlignResult result = edlibAlign(
            lhs.c_str(), lhs.length(),
            rhs.c_str(), rhs.length(),
            edlibNewAlignConfig(
                -1,
                EDLIB_MODE_HW,
                EDLIB_TASK_DISTANCE,
                kWildcardEqualities.data(),
                kWildcardEqualities.size()
            )
        );

        double distance = (result.status == EDLIB_STATUS_OK && result.editDistance >= 0)
            ? static_cast<double>(result.editDistance)
            : static_cast<double>(std::max(lhs.size(), rhs.size()));
        edlibFreeAlignResult(result);

        double denom = static_cast<double>(std::max(lhs.size(), rhs.size()));
        return denom == 0.0 ? 0.0 : distance / denom;
    }

    void prune_similar_reverse_adapters(ReadLayout& layout) {
        auto& by_id_index = layout.by_id();
        const double similarity_threshold = 0.1;
        std::vector<std::string> reverse_to_remove;

        for (const auto& [rev_id, rev_seq] : reverse_adapters) {
            double best_fraction = 1.0;
            for (const auto& [fwd_id, fwd_seq] : forward_adapters) {
                double fraction = alignment_distance_fraction(rev_seq, fwd_seq);
                best_fraction = std::min({best_fraction, fraction});
                if (best_fraction <= similarity_threshold) break;
            }
            if (best_fraction <= similarity_threshold) {
                reverse_to_remove.push_back(rev_id);
            }
        }

        for (const auto& adapter_id : reverse_to_remove) {
            if (auto it = by_id_index.find(adapter_id); it != by_id_index.end()) {
                RAD_COUT << "[prune_similar_reverse_adapters] Removing reverse adapter " << adapter_id
                          << " due to similarity with forward adapters.\n";
                by_id_index.erase(it);
            }
            consensus_matrices.erase(adapter_id);
        }

        if (!reverse_to_remove.empty()) {
            reverse_adapters.erase(
                std::remove_if(
                    reverse_adapters.begin(),
                    reverse_adapters.end(),
                    [&](const auto& entry) {
                        return std::find(reverse_to_remove.begin(), reverse_to_remove.end(), entry.first) != reverse_to_remove.end();
                    }),
                reverse_adapters.end()
            );
        }
    }

    /**
     * @brief process a chunk of reads to find perfect matches and update misalignment statistics
     * @param chunk vector of read sequences
     * @param forward_count count of perfect forward matches
     * @param reverse_count count of perfect reverse matches
     * @param adapter_stats map of adapter IDs to their alignment statistics
     * @param misalignment_stats map of adapter IDs to their misalignment statistics
     * @param total_reads total number of reads processed
     */
    void process_misalignment_chunk(
        const std::vector<read_streaming::sequence>& chunk, 
        size_t& forward_count,
        size_t& reverse_count,
        std::unordered_map<std::string, static_alignment_stats>& adapter_stats,
        std::unordered_map<std::string, static_alignment_stats>& misalignment_stats,
        size_t total_reads,
        const std::function<void()>& boundary_poll = {}
    ) {
        size_t processed = 0;
        for (const auto& read : chunk) {
            if (boundary_poll && processed++ % 128 == 0) boundary_poll();

            bool perfect_forward = find_perfect_match(read.seq, forward_adapters, adapter_stats, total_reads, "forward");
            bool perfect_reverse = find_perfect_match(read.seq, reverse_adapters, adapter_stats, total_reads, "reverse");
            
            if (perfect_forward) {
                // Update misalignment stats for REVERSE adapters
                update_misalignment_stats(
                    read.seq, reverse_adapters, adapter_stats, misalignment_stats, total_reads, forward_count, reverse_count
                );
                #pragma omp critical
                {
                    if (keep_perf_matches) perfect_matches.push_back({read, "forward"});
                    #pragma omp atomic
                    forward_count++;
                }
            }
            if (perfect_reverse) {
                // Update misalignment stats for REVERSE adapters
                update_misalignment_stats(
                    read.seq, forward_adapters, adapter_stats, misalignment_stats, total_reads, forward_count, reverse_count
                );
                #pragma omp critical
                {
                    if (keep_perf_matches) perfect_matches.push_back({read, "reverse"});
                    #pragma omp atomic
                    reverse_count++;
                }
            }
        }
    }

    /**
     * @brief update misalignment statistics for a given read sequence and adapters
     * @param read_seq the read sequence to check
     * @param adapters vector of adapter ID and sequence pairs
     * @param adapter_stats map of adapter IDs to their misalignment statistics
     * @param misalignment_stats map of adapter IDs to their misalignment position statistics
     * @param total_reads total number of reads processed
     * @param forward_counts count of perfect forward matches
     * @param reverse_counts count of perfect reverse matches
     */
    void update_misalignment_stats(const std::string& read_seq, 
                       const std::vector<std::pair<std::string, std::string>>& adapters,
                       std::unordered_map<std::string, static_alignment_stats>& adapter_stats,
                       std::unordered_map<std::string, static_alignment_stats>& misalignment_stats,
                       size_t total_reads, 
                       size_t forward_counts, 
                       size_t reverse_counts
                    ) {

        for (const auto& [id, seq] : adapters) {

            EdlibAlignResult result = edlibAlign(
                seq.c_str(), seq.length(),
                read_seq.c_str(), read_seq.length(),
                edlibNewAlignConfig(
                    -1, 
                    EDLIB_MODE_HW, 
                    EDLIB_TASK_LOC, 
                    kWildcardEqualities.data(), 
                    kWildcardEqualities.size()
                )
            );
                      
            std::string aligned_seq = collect_aligned_seq(read_seq, result);
            #pragma omp critical
            {
                //this calculates the misalignment positions for the reverse primers in forward locations, and vice versa
                misalignment_stats[id].update_stats(result.editDistance, total_reads, forward_counts, reverse_counts);  
                    double norm_misal_start = (static_cast<double>(result.startLocations[0]) / read_seq.length()) * 100.0;
                    double norm_misal_stop = (static_cast<double>(result.endLocations[0]) / read_seq.length()) * 100.0;
                    misalignment_stats[id].position_stats.start_positions.push_back(norm_misal_start);
                    misalignment_stats[id].position_stats.stop_positions.push_back(norm_misal_stop);
                update_consensus_matrices(id, aligned_seq, result.editDistance);
            }
            edlibFreeAlignResult(result);
        }
    }

    /**
     * @brief write perfect match reads to FASTQ files to the desktop
     * @param write whether to write the files (default: false)
     */
    void write_perfect_matches(bool write = false) {
        if (!write) {
            return;
        }
        const char* home = std::getenv("HOME");
        if (home == nullptr || *home == '\0') {
            throw std::runtime_error(
                "Cannot write perfect matches because HOME is not set");
        }
        const std::string desktop = std::string(home) + "/Desktop/";
        std::ofstream forward_file(desktop + "perfect_forward_reads.fastq");
        std::ofstream reverse_file(desktop + "perfect_reverse_reads.fastq");
        
        if (!forward_file.is_open() || !reverse_file.is_open()) {
            throw std::runtime_error("Failed to open output files");
        }

        for (const auto& match : perfect_matches) {
            std::ofstream& out_file = (match.direction == "forward") ? 
                                    forward_file : reverse_file;
            
            out_file << "@" << match.read.id;
            if (!match.read.comment.empty()) {
                out_file << " " << match.read.comment;
            }
            out_file << "\n"
                    << match.read.seq << "\n"
                    << "+\n"
                    << match.read.qual << "\n";
        }

        RAD_COUT << "Perfect match reads written to:\n"
                  << desktop + "perfect_forward_reads.fastq\n"
                  << desktop + "perfect_reverse_reads.fastq\n";
    }

    /**
     * @brief finds matches of adapters in a read sequence and update position statistics
     * @param read_seq the read sequence to check
     * @param adapters vector of adapter ID and sequence pairs
     * @param adapter_stats map of adapter IDs to their misalignment statistics
     * @return true if a match is found for adapters, false otherwise
     */
    bool find_perfect_match(
        const std::string& read_seq, 
        const std::vector<std::pair<std::string, std::string>>& adapters,
        std::unordered_map<std::string, static_alignment_stats>& adapter_stats,
        size_t total_reads,
        std::string direction
    ) {
        bool all_perfect = true;
        // loop over adapters
        for (const auto& [id, seq] : adapters) {
            int max_adapter_error_permissible = static_cast<int>(seq.length()*0.25); 
            // 25% error permissible, easy way to gauge w/truncants
            //there are cases like in bulk where even though the adapter hasn't been found perfectly, we should still say
            //that we can find some of this adapter within reason in the direction that we're considering it.
            //that way even if we filter the direction based on misalignment threshold, we can still have some 
            //recoverable adapters

            EdlibAlignResult result = edlibAlign(
                seq.c_str(), seq.length(),
                read_seq.c_str(), read_seq.length(),
                edlibNewAlignConfig(
                    max_adapter_error_permissible,
                    EDLIB_MODE_HW, 
                    EDLIB_TASK_LOC, 
                    kWildcardEqualities.data(), 
                    kWildcardEqualities.size())
            );
            //check status of edlib
            // this block is meant as a guard to the following
            // especially in cases where the adapters are super long, and might be truncated, or we might have gotten it wrong
            // we realy want to see if we can establish bounds for a proper adapter. this comes at the cost of
            // honestly concatenate resolution, where our boundaries for adapter presence might become a bit fuzzy. but
            // we set the perfect detected margin to be 10% or less, which means that we will only use edit distances 
            // greater than 0 when of the total reads that we've seen, 
            // the number of perfect adapter matches in total is less than 10% 

            bool edlib_status = (result.status == EDLIB_STATUS_OK);
            bool is_perfect = (result.editDistance == 0 && edlib_status);
            bool is_good_enough = (result.editDistance > 0 && edlib_status);
            if(is_good_enough){
                double perfect_detected = static_cast<double>(adapter_stats[id].count) / static_cast<double>(total_reads);
                if(perfect_detected < 0.1){
                    is_good_enough = true;
                } else {
                    is_good_enough = false;
                }
            }

            if(is_good_enough || is_perfect) {
                double normalized_start = (static_cast<double>(result.startLocations[0]) / read_seq.length()) * 100.0;
                double normalized_stop = (static_cast<double>(result.endLocations[0]) / read_seq.length()) * 100.0;
                #pragma omp critical
                {
                    adapter_stats[id].count++;
                    direction == "forward" ? adapter_stats[id].dir_counts.first++ : adapter_stats[id].dir_counts.second++;
                    adapter_stats[id].position_stats.start_positions.push_back(normalized_start);
                    adapter_stats[id].position_stats.stop_positions.push_back(normalized_stop);    
                }
            }

            edlibFreeAlignResult(result);
            // Continue through the remaining adapters even when this adapter
            // is not perfect.  Layout calibration needs observations for
            // every static element; returning here previously skipped later
            // elements such as rc_forw_primer and left their map rows blank.
            all_perfect = all_perfect && is_perfect;
        }
        return all_perfect;
    }

/**
 * @brief initialize consensus matrices for a given adapter
 * @param adapter_id the ID of the adapter
 * @param adapter_length the length of the adapter sequence
 * @param max_errors the maximum number of errors to consider
 * @brief initializes consensus matrices for edit distances from 0 to max_errors
 */
    void init_consensus_matrices(const std::string& adapter_id, size_t adapter_length, int max_errors) {
        consensus_matrices[adapter_id].clear();
        for(int ed = 0; ed <= max_errors; ed++) {
            consensus_matrices[adapter_id].emplace_back(adapter_length, ed);
        }
    }
/**
 * @brief update consensus matrices with a new aligned sequence and its edit distance
 * @param adapter_id the ID of the adapter
 * @param aligned_seq the aligned sequence
 * @param edit_distance the edit distance of the aligned sequence
 */
    void update_consensus_matrices(const std::string& adapter_id, const std::string& aligned_seq, int edit_distance) {
            // Add sequence to all matrices that accept this edit distance or higher
            for(auto& matrix : consensus_matrices[adapter_id]) {
                if(edit_distance <= matrix.max_edit_distance) {
                    matrix.ensure_length(aligned_seq.length());
                    for(size_t pos = 0; pos < aligned_seq.length(); pos++) {
                        matrix.position_counts[pos][aligned_seq[pos]]++;
                    }
                }
            }
        }
/**
 * @brief collect the aligned portion of a read sequence based on alignment result
 * @param read_seq the read sequence
 * @param result the Edlib alignment result
 * @return the aligned portion of the read sequence
 */
    std::string collect_aligned_seq(const std::string& read_seq, EdlibAlignResult& result) {
        // Extract aligned portion
        int start_pos = result.startLocations[0];
        int end_pos = result.endLocations[0] + 1;
        std::string aligned = read_seq.substr(start_pos, end_pos - start_pos);
        return aligned;
    }
/**
 * @brief generate a masked adapter sequence based on consensus matrices and misalignment statistics
 * @param adapter_id the ID of the adapter
 * @param adapter_seq the original adapter sequence
 * @param misalign_mean mean of the misalignment edit distances used for the thresholds (chance hits, or the fallback)
 * @param misalign_sd standard deviation of the same values
 * @return the masked adapter sequence
 */
    std::string generate_masked_adapter(const std::string& adapter_id, const std::string adapter_seq,
                                        double misalign_mean, double misalign_sd) {
    int threshold = static_cast<int>(std::floor(misalign_mean - misalign_sd));
    // Define max consecutive N's allowed as (threshold - 1), with a minimum of 1.
    int max_consecutive_N = (threshold > 1 ? threshold - 1 : 1);
    int max_total_N = 12;
    const auto& matrices = consensus_matrices[adapter_id];
    
    // Get matrix where max_edit_distance is >= threshold.
    auto it = std::find_if(matrices.begin(), matrices.end(),
        [threshold](const consensus_matrix& matrix) {
            return matrix.max_edit_distance >= threshold;
        });
    
    if(it == matrices.end()) return "";
    
    // Build a consensus matrix for each position.
    std::vector<std::vector<double>> consensus_matrix(it->position_counts.size(), std::vector<double>(4, 0.0));
    std::string masked_seq;
    
    for (size_t pos = 0; pos < it->position_counts.size() && pos < adapter_seq.length(); pos++) {
        int total_count = 0;
        double max_percentage = 0.0;
        char original_base = adapter_seq[pos];
    
        for (const auto& [base, count] : it->position_counts[pos]) {
            total_count += count;
        }
        // Calculate percentages and find dominant base.
        for (const auto& [base, count] : it->position_counts[pos]) {
            double percentage = (total_count > 0) ? static_cast<double>(count) / total_count * 100.0 : 0.0;
            size_t idx;
            switch(base) {
                case 'A': idx = 0; break;
                case 'C': idx = 1; break;
                case 'G': idx = 2; break;
                case 'T': idx = 3; break;
                default: continue;
            }
            consensus_matrix[pos][idx] = percentage;
            if (percentage > max_percentage) {
                max_percentage = percentage;
            }
        }
        // Decide the masked character for this position.
        // masked_seq += (max_percentage > 50.0) ? original_base : 'N';
        // if the character is 'N', enforce a maximum run length.
        char masked_char = (max_percentage > 50.0) ? original_base : 'N';
        if (masked_char == 'N') {
            // Count how many consecutive 'N's are already in masked_seq.
            int n_run = 0;
            for (int i = static_cast<int>(masked_seq.size()) - 1; i >= 0; i--) {
                if (masked_seq[i] == 'N')
                    n_run++;
                else
                    break;
            }
            // Also count the total number of 'N's in the masked sequence.
            int total_N_count = std::count(masked_seq.begin(), masked_seq.end(), 'N');
            // If we already have too many consecutive 'N's, or too many 'N's in total,
            // then use the original base instead.
            if (n_run >= max_consecutive_N || total_N_count >= max_total_N) {
                masked_char = original_base;
                //continue;
            }
        }
        masked_seq.push_back(masked_char);
    }
    return masked_seq;
}

/**
 * @brief update the ReadLayout with misalignment thresholds and masked adapters
 * @param layout the ReadLayout to update
 * @param adapter_stats map of adapter IDs to their misalignment statistics
 * @param misalignment_stats map of adapter IDs to their misalignment position statistics
 */
    void update_read_layout(
        ReadLayout& layout,
         std::unordered_map<std::string, static_alignment_stats>& adapter_stats,
         std::unordered_map<std::string, static_alignment_stats>& misalignment_stats) {
        
        auto& by_id_index = layout.by_id();
        auto& by_dir_index = layout.by_direction();
        auto color_id = [](const ReadElement& elem) {
            return term_print_utils::colorize(elem.class_id, term_print_utils::class_color(elem.class_id, elem.seq));
        };
        auto color_text = [](const std::string& text, const ReadElement& elem) {
            return term_print_utils::colorize(text, term_print_utils::class_color(elem.class_id, elem.seq));
        };

        auto has_core_static = [&](const std::string& dir) {
            auto range = by_dir_index.equal_range(dir);
            for (auto it = range.first; it != range.second; ++it) {
                if (it->type == "static" &&
                    it->global_class != "start" &&
                    it->global_class != "stop" &&
                    it->global_class != "poly_tail") {
                    return true;
                }
            }
            return false;
        };

        for(auto& [adapter_id, stats] : adapter_stats) {
            auto it = by_id_index.find(adapter_id);
            if(it != by_id_index.end()) {
                if(stats.dir_counts.first == 0 && stats.dir_counts.second == 0) {
                    RAD_COUT << "[update_read_layout] Removing adapter "
                              << color_id(*it)
                              << " due to lack of perfect matches in either direction.\n";
                    by_id_index.erase(it);
                    consensus_matrices.erase(adapter_id);
                    continue;
                }

                std::string seq = it->seq;

                AlignmentPositions align_pos = stats.position_stats.calculate();
                AlignmentPositions misalign_pos = misalignment_stats[adapter_id].position_stats.calculate();

                if(align_pos.start_stats.first == 0 && align_pos.stop_stats.first == 0){
                    RAD_COUT << "[update_read_layout] Removing element " << adapter_id << 
                    " due to insufficient alignment data.\n";
                    by_id_index.erase(it);
                    continue;
                }

                // thresholds from the chance hits only
                auto& misalignment = misalignment_stats[adapter_id];
                const chance_hit_stats chance =
                    chance_hits_from_histogram(misalignment.ed_hist, misalign_lower_chance_share_);
                const double mean_all = misalignment.mean;
                const double sd_all = misalignment.get_sd();
                double mean = chance.mean;
                double sd = chance.sd;
                int lower = chance.lower;
                const int fallback_threshold =
                    adapter_thresholds::fallback_max_edit_distance(seq.length());
                const bool use_fallback_threshold =
                    chance.n_chance < adapter_thresholds::min_misal_obs ||
                    !std::isfinite(mean) ||
                    mean < 1.0;
                if (use_fallback_threshold) {
                    mean = static_cast<double>(fallback_threshold);
                    sd = 0.0;
                    lower = fallback_threshold;
                }
                // the masking generator reads the same statistics as the thresholds
                std::string masked = generate_masked_adapter(adapter_id, seq, mean, sd);

                std::tuple<int, int, int> threshold =
                    use_fallback_threshold
                        ? std::make_tuple(
                              fallback_threshold,
                              fallback_threshold,
                              fallback_threshold)
                        : std::make_tuple(
                              lower,
                              static_cast<int>(std::floor(mean)),
                              static_cast<int>(std::round(mean + sd)));

                std::ostringstream hist_text;
                for (size_t d = 0; d < misalignment.ed_hist.size(); ++d) {
                    if (misalignment.ed_hist[d] == 0) continue;
                    if (hist_text.tellp() > 0) hist_text << " ";
                    hist_text << d << ":" << misalignment.ed_hist[d];
                }
                auto pct = [](double x) { return std::round(x * 1000.0) / 10.0; };
                std::ostringstream low_text;
                if (chance.split >= 0)
                    low_text << "ED <= " << chance.split << ": " << chance.n_low << " ("
                             << pct(static_cast<double>(chance.n_low) / static_cast<double>(std::max<size_t>(1, chance.n_all)))
                             << "%)";
                else
                    low_text << "none (no gap below the mode)";

                double misal_start = 0.0;
                double misal_stop = 0.0;
                if(misalign_pos.start_stats.first == 0.0 && misalign_pos.stop_stats.first == 0.0){
                    misal_start = std::abs(100.0 - align_pos.start_stats.first);
                    misal_stop =  std::abs(100 - align_pos.stop_stats.first);
                    misalign_pos.start_stats = {misal_start, sd};
                    misalign_pos.stop_stats = {misal_stop, sd};
                }

                RAD_COUT
                << "\n[update_read_layout] Adapter: " << color_text(adapter_id, *it) << "\n"
                << "+------------------------------+----------------------+\n"
                << "| Metric                       | Value                |\n"
                << "+------------------------------+----------------------+\n"
                << "| original seq                 | " << seq << "\n"
                << "| masked seq                   | " << masked << "\n"
                << "| misalignment observations    | " << misalignment.count << "\n"
                << "| misalignment ED histogram    | " << hist_text.str() << "\n"
                << "| true-copy group (removed)    | " << low_text.str() << "\n"
                << "| chance hits (mode ED)        | " << chance.n_chance << " (" << chance.mode << ")\n"
                << "| misalignment mean            | " << mean << "\n"
                << "| misalignment sd              | " << sd << "\n"
                << "| mean / sd of all values      | " << mean_all << " / " << sd_all << "\n"
                << "| chance share <= lower, +1    | " << pct(chance.share_at_lower) << "% / "
                << pct(chance.share_above_lower) << "% (limit " << pct(misalign_lower_chance_share_) << "%)\n"
                << "| threshold source             | "
                << (use_fallback_threshold
                        ? "fallback ceil(0.30 * adapter length)"
                        : "chance hits of the misalignment distribution")
                << "\n"
                << "| misalign lower / mean / upper| " << std::get<0>(threshold) << " / "
                << std::get<1>(threshold) << " / " << std::get<2>(threshold) << "\n"
                << "| perfect adapter count        | " << adapter_stats[adapter_id].count << "\n"
                << "| forward count                | " << adapter_stats[adapter_id].dir_counts.first << "\n"
                << "| reverse count                | " << adapter_stats[adapter_id].dir_counts.second << "\n"
                << "| aligned start pos (mean)     | " << align_pos.start_stats.first << "\n"
                << "| aligned stop pos (mean)      | " << align_pos.stop_stats.first << "\n"
                << "| misaligned start pos (mean)  | " << misalign_pos.start_stats.first << "\n"
                << "| misaligned stop pos (mean)   | " << misalign_pos.stop_stats.first << "\n"
                << "+------------------------------+----------------------+\n";

                by_id_index.modify(it, [&](ReadElement& elem) {
                    elem.masked_seq = masked;
                    elem.misalignment_threshold = threshold;
                    elem.aligned_positions = align_pos;
                    elem.misaligned_positions= misalign_pos;
                });
            }
        }

        // Remove directions that no longer contain core static adapters (excluding start/stop/poly_tail).
        std::unordered_set<std::string> present_dirs;
        for (auto it = by_dir_index.begin(); it != by_dir_index.end(); ++it) {
            present_dirs.insert(it->direction);
        }

        for (const auto& dir : present_dirs) {
            if (has_core_static(dir)) continue;
            auto range = by_dir_index.equal_range(dir);
            if (range.first == range.second) continue;
            RAD_COUT << "[update_read_layout] Removing direction " << dir
                      << " due to absence of static adapters.\n";
            while (range.first != range.second) {
                range.first = by_dir_index.erase(range.first);
            }
        }
    }
};
