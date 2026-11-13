#pragma once
#include "rad_headers.h"

/**
 * @brief Represents a processed element with alignment information
 * @param class_id ``std::string`` Unique identifier for the element
 * @param global_class ``std::string`` Global classification
 * @param edit_distance ``std::optional<int>`` Edit distance from alignment
 * @param position ``std::pair<int, int>`` Start and stop positions
 * @param type ``std::string`` Element type (variable or static)
 * @param order ``int`` Position in layout
 * @param direction ``std::string`` Orientation (forward or reverse)
 * @param element_pass ``std::optional<bool>`` Whether element passed validation--optional for variable elements
 * @param write ``std::optional<bool>`` Explicit control over whether to emit this element
 * @param seq ``std::optional<std::string>`` sequence for this element
 * @param qual ``std::optional<std::string>`` quality scores for sequencing data (if fastq)
 * @param original_seq ``std::optional<std::string>`` original (pre-corrected) sequence if available
 */
struct seq_element {
    std::string class_id;               // Unique identifier for the element
    std::string global_class;           // Global classification
    std::optional<int> edit_distance;    //Edit distance from alignment
    std::pair<int, int> position;       // Start and stop positions
    std::string type;                   // Element type
    int order;                          // Position in layout
    std::string direction;              // Orientation
    std::optional<bool> element_pass;  // Whether element passed validation
    std::optional<bool> write;         // Explicit control over whether to emit this element
    std::optional<std::string> seq;    // sequence for this element
    std::optional<std::string> qual;   // quality scores for sequencing data (if fastq)
    std::optional<std::string> original_seq;  // original (pre-corrected) sequence if available

    // Resolved during barcode correction — avoids redundant lookups in update_bc_counts.
    counter* resolved_counter = nullptr;       // pointer to the matched whitelist entry's counter
    bool     resolved_corrected = false;       // true if barcode was error-corrected (not just RC)
    bool query_complete = true;               // static-query coverage; independent of legacy acceptance

    seq_element(
        std::string class_id,
        std::string global_class,
        std::optional<int> edit_distance,
        std::pair<int, int> position,
        std::string type,
        int order,
        std::string direction,
        std::optional<bool> element_pass = std::nullopt,
        std::optional<bool> write = std::nullopt,
        std::optional<std::string> seq = std::nullopt,
        std::optional<std::string> qual = std::nullopt,
        std::optional<std::string> original_seq = std::nullopt
    ) : class_id(std::move(class_id)),
        global_class(std::move(global_class)),
        edit_distance(edit_distance),
        position(position),
        type(std::move(type)),
        order(order),
        direction(std::move(direction)),
        element_pass(element_pass),
        write(write),
        seq(std::move(seq)),
        qual(std::move(qual)),
        original_seq(std::move(original_seq)) {}
};

/**
 * @brief Represents a static processed element with alignment information
 * @param positions ``std::vector<std::pair<int, int>>`` List of start and stop positions for multiple alignments
 * @param edit_distance ``int`` Edit distance from alignment
 * @param success ``bool`` Whether alignment was successful
 * @param seq ``std::string`` Aligned sequence
 * @param cigar ``std::string`` CIGAR string representing alignment
 * @param score ``int`` Alignment score
 * @param pos ``std::pair<int, int>`` Start and stop positions of the best alignment
 */
struct static_alignments {
    std::vector<std::pair<int, int>> positions;
    int edit_distance = -1;
    bool success = false;
    bool query_complete = false;
    int query_clip_left = 0;
    int query_clip_right = 0;
    std::string seq;
    std::string cigar; 
    int score = 0;
    std::pair<int,int> pos{0, 0};
    int ref_begin = -1;   // SSW: aligned adapter offsets (clipping); -1 when not from SSW
    int ref_end = -1;
};

// Tags for the multi-index container
struct sig_id_tag {};
struct sig_global_tag {};
struct sig_ed_tag {};
struct sig_dir_tag {};
struct sig_order_tag {};
struct sig_pass_tag {};
struct sig_read_tag {};

/**
 * @brief Multi-index container for seq_element, indexed by various attributes. The intention was to make a structure that
 * can efficiently store and index processed sequencing reads to be sorted in different ways. Each seq_element should represent
 * a single element within a read, such as a barcode, UMI, or adapter, along with its alignment information and validation status.
 * Gets a little tricky when there's more than one element with the same class_id in the same read (adapters).
 *  
 * @param class_id `std::string` unique id for the element
 * @param global_class `std::string` global classification ('read', 'adapter', etc.)
 * @param edit_distance `std::optional<int>` edit distance from alignment
 * @param order `int` position in layout
 * @param direction `std::string` orientation (forward or reverse)
 * @param element_pass `std::optional<bool>` whether element passed validation
 */
typedef boost::multi_index::multi_index_container<
    seq_element,
    boost::multi_index::indexed_by<
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<sig_id_tag>,
            boost::multi_index::member<seq_element, std::string, &seq_element::class_id>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<sig_global_tag>,
            boost::multi_index::member<seq_element, std::string, &seq_element::global_class>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<sig_ed_tag>,
            boost::multi_index::member<seq_element, std::optional<int>, &seq_element::edit_distance>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<sig_order_tag>,
            boost::multi_index::member<seq_element, int, &seq_element::order>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<sig_dir_tag>,
            boost::multi_index::member<seq_element, std::string, &seq_element::direction>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<sig_pass_tag>,
            boost::multi_index::member<seq_element, std::optional<bool>, &seq_element::element_pass>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<sig_read_tag>,
            boost::multi_index::composite_key<
                seq_element,
                boost::multi_index::member<seq_element, std::string, &seq_element::direction>,
                boost::multi_index::member<seq_element, std::optional<bool>, &seq_element::element_pass>,
                boost::multi_index::member<seq_element, std::pair<int,int>, &seq_element::position>
            >
        >
    >
> SigElement;

class aligner_tools {
private:
/**
 * @brief Get the maximum number of consecutive matches from a CIGAR string
 * @param cigar CIGAR string
 * @return Maximum number of consecutive matches
 */
    int get_max_consecutive_matches(const char* cigar) {
        int max_matches = 0;
        int current_matches = 0;
        int number = 0;
        for (const char* p = cigar; *p; p++) {
            if (std::isdigit(*p)) {
                number = number * 10 + (*p - '0');
            } else {
                // For simplicity, treat 'M' and '=' as matches.
                if (*p == 'M' || *p == '=') {
                    current_matches = number;
                } else {
                    current_matches = 0;
                }
                if (current_matches > max_matches)
                    max_matches = current_matches;
                number = 0;
            }
        }
        return max_matches;
    }

/**
 * @brief Get the maximum number of consecutive matches allowing for small indels from a CIGAR string
 * @param cigar CIGAR string
 * @return Maximum number of consecutive matches with at most one indel
 */
    int get_max_consecutive_matches_with_indels(const char* cigar) {
        // First, parse the cigar string into segments.
        std::vector<std::pair<int, char>> segments;
        while (*cigar) {
            int count = 0;
            while (isdigit(*cigar)) {
                count = count * 10 + (*cigar - '0');
                ++cigar;
            }
            char op = *cigar;
            ++cigar; // skip the op
            segments.push_back({count, op});
        }
        
        int best = 0;
        // For each segment that is a match, consider it as a window start.
        for (size_t i = 0; i < segments.size(); i++) {
            // Only start at a match segment.
            if (!(segments[i].second == '=' || segments[i].second == 'M'))
                continue;
            int current = 0;
            bool used_indel = false;
            // Now extend the window from i forward.
            for (size_t j = i; j < segments.size(); j++) {
                char op = segments[j].second;
                int count = segments[j].first;
                if (op == '=' || op == 'M') {
                    current += count;
                } else if ((op == 'I' || op == 'D' || op == 'X') && count <= 2) {
                    if (!used_indel) {
                        used_indel = true;
                    } else {
                        // Already used an indel; break out.
                        break;
                    }
                } else {
                    // For any clipping (S/H) or other operation, break.
                    break;
                }
                best = std::max(best, current);
            }
        }
        return best;
    }
/**
 * @brief Compute the edit distance from a CIGAR string
 * @param cigar CIGAR string
 * @return Edit distance
 */
    int compute_edit_distance(const char* cigar) {
        int edit_distance = 0;
        int number = 0;
        while (*cigar) {
            if (std::isdigit(*cigar)) {
                number = number * 10 + (*cigar - '0');
            } else {
                // In an extended CIGAR:
                // '=' represents a match (ignored),
                // 'X' represents a mismatch,
                // 'I' represents an insertion,
                // 'D' represents a deletion.
                if (*cigar == 'X' || *cigar == 'I' || *cigar == 'D') {
                    edit_distance += number;
                }
                number = 0;
            }
            ++cigar;
        }
        return edit_distance;
    }


public:

/**
 * @brief Extract N-masked regions from an aligned sequence using CIGAR-based position mapping
 * @param result The static_alignments struct containing alignment info (seq, cigar)
 * @param query The original query/template sequence containing N's
 * @return String containing extracted N-masked region(s) joined by "_", or empty if no N's
 */
    std::string extract_n_masked_regions(const static_alignments& result, const std::string& query) {
        if (query.empty() || result.seq.empty()) {
            return "";
        }

        // Find N-region boundaries in the query
        std::vector<std::pair<size_t, size_t>> n_regions;  // start, end (exclusive)
        size_t n_start = std::string::npos;
        for (size_t i = 0; i < query.size(); i++) {
            if (query[i] == 'N' || query[i] == 'n') {
                if (n_start == std::string::npos) {
                    n_start = i;
                }
            } else {
                if (n_start != std::string::npos) {
                    n_regions.push_back({n_start, i});
                    n_start = std::string::npos;
                }
            }
        }
        // Handle trailing N's
        if (n_start != std::string::npos) {
            n_regions.push_back({n_start, query.size()});
        }

        if (n_regions.empty()) {
            return "";
        }

        // Build query-to-target position mapping using CIGAR if available
        std::vector<int> query_to_target(query.size(), -1);  // -1 = no mapping (deletion)

        if (!result.cigar.empty()) {
            // Parse CIGAR and map positions
            size_t q_pos = 0;  // position in query
            size_t t_pos = 0;  // position in target (result.seq)
            const char* cig = result.cigar.c_str();

            while (*cig) {
                int count = 0;
                while (std::isdigit(*cig)) {
                    count = count * 10 + (*cig - '0');
                    ++cig;
                }
                char op = *cig;
                ++cig;

                for (int i = 0; i < count; i++) {
                    if (op == '=' || op == 'X' || op == 'M') {
                        // Match/mismatch: both query and target advance
                        if (q_pos < query.size() && t_pos < result.seq.size()) {
                            query_to_target[q_pos] = t_pos;
                        }
                        q_pos++;
                        t_pos++;
                    } else if (op == 'I') {
                        // Insertion in query: query advances, target doesn't
                        q_pos++;
                    } else if (op == 'D') {
                        // Deletion in query: target advances, query doesn't
                        t_pos++;
                    }
                }
            }
        } else {
            // No CIGAR, assume 1:1 positional mapping
            for (size_t i = 0; i < query.size() && i < result.seq.size(); i++) {
                query_to_target[i] = i;
            }
        }

        // Extract N-masked regions using the position mapping
        std::vector<std::string> extracted;
        for (const auto& region : n_regions) {
            std::string segment;
            for (size_t q = region.first; q < region.second; q++) {
                if (q < query_to_target.size() && query_to_target[q] >= 0
                    && static_cast<size_t>(query_to_target[q]) < result.seq.size()) {
                    segment += result.seq[query_to_target[q]];
                }
            }
            if (!segment.empty()) {
                extracted.push_back(segment);
            }
        }

        // Join with "-"
        std::string output;
        for (size_t i = 0; i < extracted.size(); i++) {
            if (i > 0) output += "-";
            output += extracted[i];
        }
        return output;
    }

/**
 * @brief This is the master alignment function for static elements, using Edlib and SSW.
 * @param query ``std::string`` Query sequence
 * @param target ``std::string`` Target sequence
 * @param verbose ``bool`` Whether to print verbose output
 * @param max_edit_distance ``int`` Maximum allowed edit distance for Edlib alignment
 * @param masked_query ``std::string`` Masked query sequence (optional)
 * @param primary ``bool`` Whether this is the primary alignment attempt
 * @param expected_start ``int`` Expected start position
 * @param expected_end ``int`` Expected end position
 * @return ``static_alignments`` Struct containing alignment results for this element
 * 
 * @brief This function first attempts to align the query to the target using Edlib with specified parameters.
 * If Edlib returns exactly one candidate alignment within the expected region, that alignment is used.
 * If multiple candidates or no candidates are found, the function falls back to using the Striped Smith-Waterman (SSW)
 * algorithm for alignment. Scoring:
 * - Match: +2
 * - Mismatch: -2
 * - Gap Open: -3
 * - Gap Extend: -2
 *  The function returns a `static_alignments` struct containing the results of the alignment.
 */
    static_alignments align_static_elements(
        const std::string& query, const std::string& target, bool verbose, int max_edit_distance = -1, 
        const std::string& masked_query = "", bool primary = true, int expected_start = 1, 
        int expected_end = -1, static_alignments* additional_hits = nullptr,
        int search_lo = 1, int search_hi = -1, bool allow_ssw = true,
        bool target_at_read_start = true, bool target_at_read_end = true) {

        static_alignments alignment;
        if (additional_hits) *additional_hits = static_alignments{};
        alignment.success = false;
        alignment.edit_distance = -1;
        // Our threshold for a “good” alignment
        const int min_match_bases = 5;
        
        // If expected_end is not provided, use the target's length
        if (expected_end < 0){
                expected_end = target.size();
        }

        // Search window [search_lo, search_hi], 1-based inclusive; positions stay in target coordinates.
        // target_at_read_start / _end: false when that target edge is an HMM cut, not a physical read end.
        const bool windowed = search_lo > 1 || search_hi >= 0;
        const int win_off = windowed ? std::max(0, search_lo - 1) : 0;
        const int win_end = (search_hi < 0 || search_hi > static_cast<int>(target.size()))
            ? static_cast<int>(target.size()) : search_hi;
        const int win_len = windowed ? win_end - win_off : static_cast<int>(target.size());
        if (windowed && win_len <= 0) return alignment;

        // Custom equality rules:
        // Uppercase A, C, T, G match their lowercase forms and also the pad character 'x'.
        // The masked letter N matches uppercase A, C, T, G and 'x', but not lowercase.
        const int numEq = 13;
        EdlibEqualityPair customEqualities[numEq] = {
            {'A','a'}, {'C','c'}, {'T','t'}, {'G','g'},
            {'A','x'}, {'C','x'}, {'T','x'}, {'G','x'},
            {'N','A'}, {'N','C'}, {'N','T'}, {'N','G'},
            {'N','x'}
        };

        std::string query_to_use = query;

        // --- run Edlib to obtain candidate intervals ---
        EdlibAlignConfig config = edlibNewAlignConfig(max_edit_distance,
                                                    EDLIB_MODE_HW,
                                                    EDLIB_TASK_LOC,
                                                    customEqualities, numEq);

        EdlibAlignResult edlibResult = edlibAlign(query_to_use.c_str(), query_to_use.size(),
                                                target.c_str() + win_off, win_len,
                                                config);

        std::vector<std::pair<int,int>> candidates;
        if (edlibResult.status == EDLIB_STATUS_OK &&
            edlibResult.numLocations > 0 &&
            edlibResult.editDistance > -1) {
            std::vector<std::pair<int,int>> intervals;
            for (int i = 0; i < edlibResult.numLocations; ++i) {
                int start = edlibResult.startLocations[i] + 1 + win_off;  // Convert to 1-indexed
                int end   = edlibResult.endLocations[i] + 1 + win_off;
                intervals.push_back({start, end});
            }
            std::sort(intervals.begin(), intervals.end(),
                    [](const std::pair<int,int>& a, const std::pair<int,int>& b) {
                        return a.first < b.first;
                    });
            // Collapse overlapping intervals
            if (!intervals.empty()) {
                std::pair<int,int> current = intervals[0];
                for (size_t i = 1; i < intervals.size(); ++i) {
                    if (intervals[i].first <= current.second) {
                        current.second = std::max(current.second, intervals[i].second);
                    } else {
                        candidates.push_back(current);
                        current = intervals[i];
                    }
                }
                candidates.push_back(current);
            }
            // Retain the already allocated, actual LOC intervals rather than
            // copying them or treating an overlap union as a real alignment.
            // All share editDistance; primary/SSW selection stays unchanged.
            if (additional_hits) {
                additional_hits->positions = std::move(intervals);
                additional_hits->edit_distance = edlibResult.editDistance;
                additional_hits->success = true;
                additional_hits->query_complete = true;
            }
        }

        // --- early out: If Edlib returned exactly one candidate and it is within the expected region, use it
        if (candidates.size() == 1) {
            auto cand = candidates.front();
            if (cand.first >= expected_start && cand.second <= expected_end) {
                alignment.positions.push_back(cand);
                alignment.success = true;
                alignment.query_complete = true;
                alignment.edit_distance = edlibResult.editDistance;
                // Optionally set alignment.seq, alignment.cigar, etc.
                alignment.seq = target.substr(cand.first - 1, cand.second - cand.first + 1);
                // Return immediately without SSW
                edlibFreeAlignResult(edlibResult);
                return alignment;
            }
        }
        if (!allow_ssw) {
            if (!candidates.empty()) {
                const double centre = win_off + 1 + (win_len - 1) / 2.0;
                auto dist = [centre](const std::pair<int,int>& c) {
                    return std::fabs((c.first + c.second) / 2.0 - centre);
                };
                auto best = std::min_element(candidates.begin(), candidates.end(),
                    [&](const std::pair<int,int>& a, const std::pair<int,int>& b) { return dist(a) < dist(b); });
                alignment.positions.push_back(*best);
                alignment.success = true;
                alignment.query_complete = true;
                alignment.edit_distance = edlibResult.editDistance;
                alignment.seq = target.substr(best->first - 1, best->second - best->first + 1);
            }
            edlibFreeAlignResult(edlibResult);
            return alignment;
        }
        edlibFreeAlignResult(edlibResult);
        // --- run SSW on the entire target (for no candidate or multiple candidates) ---
        // Prepare uppercase strings for SSW.
        std::string query_ssw = query;
        std::string target_ssw = windowed ? target.substr(win_off, win_len) : target;
        std::transform(query_ssw.begin(), query_ssw.end(), query_ssw.begin(), ::toupper);
        std::transform(target_ssw.begin(), target_ssw.end(), target_ssw.begin(), ::toupper);

        int32_t maskLen = static_cast<int32_t>(query_ssw.size() / 2);
        if (maskLen < 15)
            maskLen = 15;
        //defaults of match/mismatch/gap_open/gap_extend is 2/2/3/1
        //thoughts for a shorter sequence--i'm okay with a gap being opened, but the longer the gap goes the 
        //more it should cost. so if we score gap open as a mismatch, but then gap extend as a mismatch as well?
        StripedSmithWaterman::Aligner ssw_aligner(2, 2, 3, 2);
        StripedSmithWaterman::Filter ssw_filter;
        ssw_filter.report_cigar = true;
        StripedSmithWaterman::Alignment sswAlign;
        ssw_aligner.Align(target_ssw.c_str(), query_ssw.c_str(), query_ssw.size(),
                        ssw_filter, &sswAlign, maskLen);

        if (sswAlign.query_begin >= 0 && sswAlign.query_end >= sswAlign.query_begin) {
            int ssw_start = sswAlign.query_begin + 1 + win_off; // Convert to 1-indexed
            int ssw_end = sswAlign.query_end + 1 + win_off;
            int ssw_length = ssw_end - ssw_start + 1;
            int ssw_max_matches_ind = get_max_consecutive_matches_with_indels(sswAlign.cigar_string.c_str());
            int ssw_max_matches = get_max_consecutive_matches(sswAlign.cigar_string.c_str());
            // Calculate deviation from expected region.
            int deviation = 0;
            if (ssw_start < expected_start)
                deviation += (expected_start - ssw_start);
            if (ssw_end > expected_end)
                deviation += (ssw_end - expected_end);

            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "SSW alignment: region ["
                            << ssw_start << ", " << ssw_end << "], length = "
                            << ssw_length << ", max_matches with indels = " 
                            << ssw_max_matches_ind
                            << " , max_matches = " << ssw_max_matches
                            << ", deviation = " << deviation
                            << ", cigar: " << sswAlign.cigar_string << "\n";
                    RAD_COUT << oss.str();
                }
            }
            // Accept candidate if it meets threshold. first case is within the deviation allowed. second case is for concats.
            if (
            (ssw_max_matches > min_match_bases && 
                deviation <= 100 && 
                ssw_length >= 10 && 
                ssw_max_matches_ind >= 8 && 
                ssw_max_matches >= 8
            )
            ||
            (
                //edited here for concat
                deviation > 100 && 
                ssw_length >= 10 && 
                ssw_max_matches_ind >= 10
            )) {
                alignment.positions.push_back({ssw_start, ssw_end});
                alignment.edit_distance = compute_edit_distance(sswAlign.cigar_string.c_str());
                // RAD passes the adapter as SSW's reference, not its query.
                // Read-span/edit-distance consistency cannot detect all clipping.
                alignment.query_complete = sswAlign.ref_begin == 0 &&
                    sswAlign.ref_end == static_cast<int>(query.size()) - 1;
                alignment.query_clip_left = sswAlign.ref_begin;
                alignment.query_clip_right = static_cast<int>(query.size()) - 1 - sswAlign.ref_end;
                alignment.ref_begin = sswAlign.ref_begin;
                alignment.ref_end = sswAlign.ref_end;
                // Clipped adapter bases are free only within read_end_slack of a physical read end.
                constexpr int read_end_slack = 3;
                const bool clipping_at_read_ends =
                    (sswAlign.ref_begin == 0 ||
                     (target_at_read_start && win_off + sswAlign.query_begin <= read_end_slack)) &&
                    (sswAlign.ref_end == static_cast<int>(query.size()) - 1 ||
                     (target_at_read_end &&
                      static_cast<int>(target.size()) - 1 - (win_off + sswAlign.query_end) <= read_end_slack));
                // Windowed calls: the caller applies its own clipped-hit rule.
                const bool read_end_rule = !windowed;
                alignment.success = alignment.edit_distance <= max_edit_distance &&
                    (clipping_at_read_ends || !read_end_rule);
                alignment.cigar = sswAlign.cigar_string;
                alignment.seq = target.substr(ssw_start - 1, ssw_end - ssw_start + 1);
                if (read_end_rule && !clipping_at_read_ends) {
                    // Internal hit: verify the whole adapter on the same interval, missing bases
                    // charged to the same edit budget (the interval must not grow into a barcode/UMI).
                    auto verified = edlibAlign(query_to_use.c_str(), query_to_use.size(),
                        alignment.seq.c_str(), alignment.seq.size(),
                        edlibNewAlignConfig(max_edit_distance, EDLIB_MODE_NW,
                                            EDLIB_TASK_PATH, customEqualities, numEq));
                    if (verified.status == EDLIB_STATUS_OK && verified.editDistance >= 0 &&
                        verified.alignment) {
                        char* cigar = edlibAlignmentToCigar(verified.alignment,
                            verified.alignmentLength, EDLIB_CIGAR_EXTENDED);
                        if (cigar) {
                            alignment.edit_distance = verified.editDistance;
                            alignment.query_complete = true;
                            alignment.cigar = cigar;
                            alignment.success = true;
                            free(cigar);
                        }
                    }
                    edlibFreeAlignResult(verified);
                }
            }
        }
    return alignment;
}
/**
 * @brief Find poly-base tails in a sequence
 * @param query ``std::string`` Query sequence (to determine poly-base)
 * @param sequence ``std::string`` Target sequence to search for poly-base tails
 * @param window_size ``int`` Size of the sliding window
 * @return ``static_alignments`` Struct containing positions of poly-base tails and edit distance
 * 
 * @brief This function scans the target sequence using a sliding window approach to identify regions
 * that are predominantly composed of a single base (the poly-base). It allows for small gaps (default: up to 3 bases)
 * between these regions (``min_gap``). If a window contains at least (default: 90%) of the poly-base (``min_count``), 
 * it is considered a potential poly-tail. 
 */
    static_alignments find_poly_tails(const std::string& query, const std::string& sequence, int window_size) {
        static_alignments result;
        result.success = false;
        result.edit_distance = 1;
        if (query.empty() || window_size <= 0 ||
            sequence.length() < static_cast<size_t>(window_size)) {
            return result;
        }

        // FASTQ bases are ASCII. Normalizing them directly avoids the
        // locale-aware std::toupper call in every overlapping window.
        const auto ascii_upper = [](char base) noexcept {
            const unsigned char value = static_cast<unsigned char>(base);
            return static_cast<char>(
                value >= static_cast<unsigned char>('a') &&
                value <= static_cast<unsigned char>('z')
                    ? value - static_cast<unsigned char>('a') +
                          static_cast<unsigned char>('A')
                    : value
            );
        };
        const char poly_base = ascii_upper(query[0]);
        const auto is_poly_base = [poly_base, &ascii_upper](char base) noexcept {
            return ascii_upper(base) == poly_base;
        };

        int min_count = static_cast<int>(window_size * 0.9);
        int min_gap = 3;
        const int last_window_start =
            static_cast<int>(sequence.length()) - window_size;
        int i = 0;

        // Seed the first window, then update its count in O(1) as it slides.
        int count = 0;
        for (int j = 0; j < window_size; ++j) {
            count += is_poly_base(sequence[j]) ? 1 : 0;
        }

        while (i <= last_window_start) {
            if (count >= min_count) {
                int current_start = i;
                int current_end = i + window_size - 1;
                int non_poly_count = 0;
                int last_poly_pos = current_end;
                while (current_end + 1 < static_cast<int>(sequence.length())) {
                    if (is_poly_base(sequence[current_end + 1])) {
                        if (non_poly_count <= min_gap) {
                            current_end++;
                            last_poly_pos = current_end;
                            non_poly_count = 0;
                        } else {
                            break;
                        }
                    } else {
                        current_end++;
                        non_poly_count++;
                        if (non_poly_count > min_gap) {
                            current_end = last_poly_pos;
                            break;
                        }
                    }
                }
                result.success = true;
                result.edit_distance = window_size - count;
                result.positions.emplace_back(current_start + 1, current_end + 1);
                // Move i to just after the end of current tail
                i = current_end + 1;
                // Look for next potential poly-tail after a gap
                i += min_gap;

                // A detected tail jumps over an arbitrary number of bases.
                // Re-seed the small window at the new position rather than
                // walking every skipped window.
                if (i <= last_window_start) {
                    count = 0;
                    for (int j = 0; j < window_size; ++j) {
                        count += is_poly_base(sequence[i + j]) ? 1 : 0;
                    }
                }
            } else {
                if (i == last_window_start) {
                    break;
                }
                count -= is_poly_base(sequence[i]) ? 1 : 0;
                count += is_poly_base(sequence[i + window_size]) ? 1 : 0;
                i++;
            }
        }
        return result;
    }

/**
 * @brief Longest adapter prefix ending at the read end (at_3prime) or suffix starting at the read start.
 * @return 1-based positions in target coordinates
 */
    static_alignments align_read_end_fragment(const std::string& adapter, const std::string& target,
                                              bool at_3prime, int max_keep = -1, int min_keep = 9) const {
        static_alignments out;
        const int m = static_cast<int>(adapter.size());
        const int n = static_cast<int>(target.size());
        const EdlibEqualityPair eq[13] = {
            {'A','a'}, {'C','c'}, {'T','t'}, {'G','g'},
            {'A','x'}, {'C','x'}, {'T','x'}, {'G','x'},
            {'N','A'}, {'N','C'}, {'N','T'}, {'N','G'}, {'N','x'}};
        const int t_hi = std::min({m, n, max_keep > 0 ? max_keep : m});
        for (int t = t_hi; t >= min_keep; --t) {
            const int k = std::max(1, t / 6);
            const int span = std::min(n, t + 4);
            const int off = at_3prime ? n - span : 0;
            const char* q = adapter.c_str() + (at_3prime ? 0 : m - t);
            EdlibAlignResult r = edlibAlign(q, t, target.c_str() + off, span,
                                            edlibNewAlignConfig(k, EDLIB_MODE_HW, EDLIB_TASK_LOC, eq, 13));
            int best = -1;
            if (r.status == EDLIB_STATUS_OK && r.editDistance >= 0 && r.numLocations > 0) {
                for (int i = 0; i < r.numLocations; ++i) {
                    const int s0 = off + r.startLocations[i], e0 = off + r.endLocations[i];
                    const bool anchored = at_3prime ? (e0 >= n - 3) : (s0 <= 2);
                    if (!anchored) continue;
                    if (best < 0 || (at_3prime ? e0 > off + r.endLocations[best] : s0 < off + r.startLocations[best]))
                        best = i;
                }
            }
            if (best >= 0) {
                const int s0 = off + r.startLocations[best], e0 = off + r.endLocations[best];
                out.positions.push_back({s0 + 1, e0 + 1});
                out.edit_distance = r.editDistance;
                out.success = true;
                out.query_complete = (t == m);
                out.ref_begin = at_3prime ? 0 : m - t;
                out.ref_end = at_3prime ? t - 1 : m - 1;
                out.seq = target.substr(s0, e0 - s0 + 1);
                edlibFreeAlignResult(r);
                return out;
            }
            edlibFreeAlignResult(r);
        }
        return out;
    }

/**
 * @brief Best adapter prefix (prefix=true) or suffix fragment inside the 1-based inclusive window [lo, hi],
 * scored t - 3*ED.
 * @return 1-based positions in target coordinates; query_complete is false
 */
    static_alignments align_adapter_fragment(const std::string& adapter, const std::string& target, int lo, int hi,
                                             bool prefix, int min_keep = 12, int max_keep = -1, bool strict = false) const {
        static_alignments out;
        const int m = static_cast<int>(adapter.size());
        const int n = static_cast<int>(target.size());
        const int off = std::max(0, lo - 1), end = std::min(n, hi);
        const int span = end - off;
        if (span < min_keep || m <= min_keep) return out;
        const EdlibEqualityPair eq[13] = {
            {'A','a'}, {'C','c'}, {'T','t'}, {'G','g'},
            {'A','x'}, {'C','x'}, {'T','x'}, {'G','x'},
            {'N','A'}, {'N','C'}, {'N','T'}, {'N','G'}, {'N','x'}};
        int best_sc = -(1 << 30);
        for (int t = std::min({m - 1, span, max_keep > 0 ? max_keep : m - 1}); t >= min_keep; --t) {
            const int k = strict ? (t >= 20 ? 1 : 0) : std::max(1, t / 8);
            if (t <= best_sc) break;  // no shorter t can beat the best score
            const char* q = adapter.c_str() + (prefix ? 0 : m - t);
            EdlibAlignResult r = edlibAlign(q, t, target.c_str() + off, span,
                                            edlibNewAlignConfig(k, EDLIB_MODE_HW, EDLIB_TASK_LOC, eq, 13));
            if (r.status == EDLIB_STATUS_OK && r.editDistance >= 0 && r.numLocations > 0 && t - 3 * r.editDistance > best_sc) {
                int best = 0;  // the hit at the inward edge
                for (int i = 1; i < r.numLocations; ++i)
                    if (prefix ? r.startLocations[i] < r.startLocations[best] : r.endLocations[i] > r.endLocations[best]) best = i;
                const int s0 = off + r.startLocations[best], e0 = off + r.endLocations[best];
                best_sc = t - 3 * r.editDistance;
                out = static_alignments{};
                out.positions.push_back({s0 + 1, e0 + 1});
                out.edit_distance = r.editDistance;
                out.success = true;
                out.query_complete = false;
                out.ref_begin = prefix ? 0 : m - t;
                out.ref_end = prefix ? t - 1 : m - 1;
                out.seq = target.substr(s0, e0 - s0 + 1);
            }
            edlibFreeAlignResult(r);
        }
        return out;
    }
};

/**
 * @namespace barcode_correction
 * @brief Namespace for barcode correction functions and utilities
 */
namespace barcode_correction {
    /**
     * @brief Where the whitelist lookup placed the barcode in the expanded window of `correct_barcode`.
     * Indices are 0-based in the oriented window: primer side first, UMI side last.
     */
    struct bc_position {
        enum kind_t : uint8_t { unknown = 0, window = 1, located = 2, tie = 3, clipped = 4 };
        kind_t kind = unknown;
        int end = -1;    // index of the last barcode base
        int n_end = 0;   // number of equally good ends
        int shift = 0;   // located: end minus the layout window end (positive = toward the UMI)

        void set_from_match(int match_end, int match_n_end) {
            end = match_end;
            n_end = match_n_end;
            kind = (match_n_end == 1 && match_end >= 0) ? located : tie;
        }
    };

    /// Half-width of the expanded barcode window; also the edit distance of the global-list checks.
    constexpr int expanded_window_pad = 4;

    // A correction window may contain 64 barcode bases plus eight flanks.
    // Use Edlib only for that bounded overflow; ordinary packed matching keeps
    // its existing fast path and tie rules.
    inline int barcode_distance(const int64_seq& query, const int64_seq& winner,
                                int max_dist, const std::string* expanded_query = nullptr,
                                int* match_end = nullptr, int* match_n = nullptr) {
        if (!expanded_query)
            return mutation_tools::int64_lvdist(query, winner, max_dist, match_end, match_n);
        if (!winner.is_valid() || expanded_query->find_first_not_of("ACGT") != std::string::npos)
            return -1;
        const auto target = winner.bits_to_sequence();
        auto result = edlibAlign(target.data(), static_cast<int>(target.size()),
                                 expanded_query->data(), static_cast<int>(expanded_query->size()),
                                 edlibNewAlignConfig(max_dist, EDLIB_MODE_HW, EDLIB_TASK_LOC, nullptr, 0));
        const int distance = result.status == EDLIB_STATUS_OK ? result.editDistance : -1;
        if (match_end) *match_end = distance >= 0 && result.numLocations ? result.endLocations[0] : -1;
        if (match_n) *match_n = distance >= 0 ? result.numLocations : 0;
        edlibFreeAlignResult(result);
        return distance;
    }
    inline void locate_winner(bc_position* pos, const int64_seq& query, const int64_seq& winner,
                              int max_dist, const std::string* expanded_query = nullptr) {
        if (!pos || (!expanded_query && query.length <= winner.length)) return;
        int match_end = -1, match_n = 0;
        if (barcode_distance(query, winner, max_dist, expanded_query, &match_end, &match_n) >= 0)
            pos->set_from_match(match_end, match_n);
    }

    /**
     * @brief Check if a candidate barcode passes quality checks based on counts from the whitelist
     * @param candidate ``int64_seq`` Candidate barcode sequence
     * @param wl ``const whitelist::wl_entry*`` Pointer to the whitelist entry
     * @param whitelist_source ``std::string`` Source of the whitelist ("global" or "true")
     * @param mode ``std::string`` Barcode correction mode ("defensive" or "offensive"), 
     * determines whether barcode is scanned against the global or true whitelist first and 
     * strictness of the quality checks (defensive = more strict, offensive = less strict)
     * Default quality parameter is 80% correction ratio for barcodes with >=10 total counts and >=2 raw counts
     * @param verbose ``bool`` Whether to print verbose output
     * @return ``bool`` True if the candidate passes quality checks
     */
    bool passes_quality_check(const int64_seq& candidate, const whitelist::wl_entry* wl, const std::string& whitelist_source, 
        const std::string& mode,
            bool verbose) {

            auto all_counts = wl->with_wl(whitelist_source, [&](const auto& typed_wl) {
                return typed_wl.get_all_bc_counts(candidate);
            });

            int raw_count = all_counts.load(barcode_counts::raw);
            int total_count = all_counts.load(barcode_counts::total);
            int corrected_count = all_counts.load(barcode_counts::corrected);
                    
            if(total_count < 0 & raw_count < 2){
               return false;
            }

            if (total_count == 0) {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Quality check for " << candidate.bits_to_sequence()
                            << " [" << whitelist_source << "]: total=0 (new barcode) -> PASS" << std::endl;
                        RAD_COUT << oss.str();
                    }
                }
                return true;
            }

            //reset the counter for low-set whitelist barcodes 
            //if they accrue over a certain number of reads after getting their freebies
            if(raw_count >= 2 && total_count < 0) {
                (void) wl->with_wl(whitelist_source, [&](const auto& typed_wl) {
                    typed_wl.set_bc_count(candidate, raw_count);
                });

                return true;
            }

            double correction_ratio = (static_cast<double>(corrected_count) +1) / (static_cast<double>(total_count) + 1);
            bool overall_pass = true;  // Default to pass
            int threshold;
            if(mode == "defensive") { 
                threshold = 5;
            }
            if(mode == "offensive") {
                threshold = 10;
            }
            if(total_count >= 10 & raw_count >= 2){
                overall_pass = correction_ratio <= 0.8;
            } else if (total_count >= threshold & raw_count < 2){
                (void) wl->with_wl(whitelist_source, [&](const auto& typed_wl) {
                    typed_wl.set_bc_count(candidate, -1);
                });
                //kill it if it gets 10 freebies with no associated raw counts
                return false;
            }
        return overall_pass;
    }

/**
 * @brief Tiebreaker for multiple candidate barcodes by selecting the best one based on edit distance and quality checks
 * @param query ``int64_seq`` Query barcode sequence
 * @param candidates ``const std::unordered_set<int64_seq>&`` Set of candidate barcode sequences
 * @param max_dist ``int`` Maximum allowed edit distance
 * @param verbose ``bool`` Whether to print verbose output
 * @param wl_type ``std::string`` Type of whitelist ("global" or "true")
 * @param wl ``const whitelist::wl_entry*`` Pointer to the whitelist entry
 * @return ``std::pair<std::optional<int64_seq>, std::optional<int>> Pair containing the resolved barcode (if any) and its edit distance
 */
    std::pair<std::optional<int64_seq>, std::optional<int>>
    resolve_multiple_hits_simple(const int64_seq& query, const std::unordered_set<int64_seq>& candidates,
        int max_dist, bool verbose, const std::string& wl_type, const whitelist::wl_entry* wl = nullptr,
        bc_position* pos = nullptr, const std::string* expanded_query = nullptr
    ) {
        std::map<int, std::unordered_set<int64_seq>> sorted;
        for (const auto& candidate : candidates) {
            const int distance = barcode_distance(query, candidate, max_dist, expanded_query);
            if (distance >= 0) sorted[distance].insert(candidate);
        }
        if (sorted.empty()) return {std::nullopt, std::nullopt};

        // Find candidates with raw_count > 0 at lowest edit distance
        for (const auto& [edit_dist, candidate_set] : sorted) {
            std::vector<int64_seq> valid_candidates;
            for (const auto& candidate : candidate_set) {

                int raw_count = wl->with_wl(wl_type, [&](const auto& typed_wl) {
                    return typed_wl.get_bc_count(candidate, barcode_counts::raw);
                });

                if (raw_count > 0) {
                    valid_candidates.push_back(candidate);
                }
            }
            
            if (valid_candidates.empty()) {
                continue; // Try next edit distance
            }
            
            if (valid_candidates.size() == 1) {
                if (verbose) RAD_COUT << "Winner: unique candidate at distance " << edit_dist << "\n";
                if(passes_quality_check(valid_candidates[0], wl, wl_type, "offensive", verbose)) {
                    locate_winner(pos, query, valid_candidates[0], max_dist, expanded_query);
                    return {valid_candidates[0], edit_dist};
                } else {
                    if (verbose) RAD_COUT << "Quality check failed for candidate, rejecting\n";
                    return {std::nullopt, edit_dist};
                }
            } else {
                if (verbose) RAD_COUT << "Multiple candidates at distance " << edit_dist << ", rejecting\n";
                return {std::nullopt, edit_dist};
            }
        }
        
        if (verbose) RAD_COUT << "No candidates with raw_count > 0, rejecting\n";
        return {std::nullopt, std::nullopt};
    }

/**
 * @brief Check a barcode against a whitelist and return a corrected barcode if found
 * @param bc ``int64_seq`` Barcode sequence to check
 * @param candidates ``const std::unordered_set<int64_seq>&`` Set of candidate barcode sequences
 * @param whitelist_type ``std::string`` Type of whitelist ("global" or "true")
 * @param max_dist ``int`` Maximum allowed edit distance
 * @param verbose ``bool`` Whether to print verbose output
 * @param mode ``std::string`` Barcode correction mode ("defensive" or "offensive")
 * @param wl ``const whitelist::wl_entry*`` Pointer to the whitelist entry
 * @return ``std::optional<int64_seq>`` Corrected barcode sequence if found, otherwise std::nullopt
 * 
 * @brief This function checks the provided barcode against a set of candidate barcodes derived from the specified whitelist.
 * It first identifies candidates that match the barcode within the allowed edit distance.
 * If no candidates are found, it returns std::nullopt.
 * If a single candidate is found, it calculates the Levenshtein distance to confirm the match and performs a quality check before returning the corrected barcode.
 * In cases where multiple candidates are found, the function employs an enhanced resolution strategy to identify the best match.
 * If a single best match is identified, it also undergoes a quality check before being returned.
 * If no suitable match is found or if quality checks fail, the function returns std::nullopt.
 */
    template <class Candidates>
    std::optional<int64_seq> check_against_wl(const int64_seq& bc, const Candidates& candidates,
        const std::string& whitelist_type, int max_dist, bool verbose, const std::string& mode, const whitelist::wl_entry* wl = nullptr,
        bc_position* pos = nullptr, const std::string* expanded_query = nullptr
    ) {


        // Early exit if whitelist is empty or no candidates match

        auto matched = wl->with_wl(whitelist_type, [&](const auto& typed_wl) {
            if(!typed_wl.empty()){
                return typed_wl.return_putative_correct_bcs(candidates);
            } else {
                return std::unordered_set<int64_seq>{};
            }
        });

        if (matched.empty()) {
            return std::nullopt;
        }
        if (verbose) {
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << (whitelist_type == "true" ? "MUTATION_CHECK_FOUND" : "GLOBAL_MUTATION_CHECK_FOUND")
                    << " (" << matched.size() << " candidates)\n";
                RAD_COUT << oss.str();
            }
        }
        if (matched.size() == 1) {
            // Single candidate case
            int64_seq putative_candidate = *matched.begin();
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << (whitelist_type == "true" ? "MUTATION_CHECK_CANDIDATE::" : "GLOBAL_MUTATION_CHECK_CANDIDATE::")
                        << putative_candidate.bits_to_sequence() << "\n";
                    RAD_COUT << oss.str();
                }
            }
            int res = barcode_distance(bc, putative_candidate, max_dist, expanded_query);
            if (res >= 0) {
                // Quality check
                if (passes_quality_check(putative_candidate, wl, whitelist_type, mode, verbose)) {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "LVDIST::" << res << "\n"
                                << (whitelist_type == "true" ? "MUTATION_CHECK_MATCHED" : "GLOBAL_MUTATION_CHECK_MATCHED")
                                << "\n";
                            RAD_COUT << oss.str();
                        }
                    }
                    locate_winner(pos, bc, putative_candidate, max_dist, expanded_query);
                    return putative_candidate;
                } else {
                    if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                    oss << (whitelist_type == "true" ? "MUTATION_CHECK_QUALITY_FAILED - rejecting barcode" 
                                                                    : "GLOBAL_MUTATION_CHECK_QUALITY_FAILED - rejecting barcode") << "\n";
                                RAD_COUT << oss.str();

                            }
                        }
                    return std::nullopt;
                }
            }
        } else if (matched.size() > 1) {
            // Multiple candidates case - use enhanced resolution
            auto [resolved, min_dist] = resolve_multiple_hits_simple(bc, matched, max_dist, verbose, whitelist_type, wl, pos, expanded_query);
            if (resolved.has_value()) {
                if (passes_quality_check(resolved.value(), wl, whitelist_type, mode, verbose)) {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << (whitelist_type == "true" ? "MUTATION_MULTIPLE_MATCHED_RESOLVED" 
                                                             : "GLOBAL_MUTATION_MULTIPLE_MATCHED_RESOLVED") << "\n";
                            RAD_COUT << oss.str();
                        }
                    }
                    return resolved;
                } else {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                                oss << (whitelist_type == "true" ? "MUTATION_MULTIPLE_RESOLVED_QUALITY_FAILED - rejecting barcode" 
                                                                 : "GLOBAL_MUTATION_MULTIPLE_RESOLVED_QUALITY_FAILED - rejecting barcode") << "\n";
                            RAD_COUT << oss.str();
                        }
                    }
                    return std::nullopt;
                }
            }
        }
        return std::nullopt;
    }

/**
 * @brief Exhaustively check a barcode against all entries in the whitelist
 * @param bc ``int64_seq`` Barcode sequence to check
 * @param whitelist_type ``std::string`` Type of whitelist ("global" or "true")
 * @param max_dist ``int`` Maximum allowed edit distance 
 * @param verbose ``bool`` Whether to print verbose output
 * @param mode ``std::string`` Barcode correction mode ("defensive" or "offensive")
 * @param wl ``const whitelist::wl_entry*`` Pointer to the whitelist entry
 * @return ``std::optional<int64_seq>`` Corrected barcode sequence if found, otherwise std::nullopt
 * 
 * @brief This function performs an exhaustive search of the provided barcode against all entries in the specified whitelist.
 * It calculates the Levenshtein distance between the input barcode and each whitelist entry, collecting those within the specified maximum distance.
 * If exactly one match is found, it undergoes a quality check before being returned as the corrected barcode.
 * In cases of multiple matches, the function identifies the best matches based on the smallest edit distance.
 * If a single best match is identified, it also undergoes a quality check before being returned.
 * If no matches are found or if quality checks fail, the function returns std::nullopt.
 */ 
    std::optional<int64_seq> exhaustive_check_against_wl(const int64_seq& bc, const std::string& whitelist_type,
        int max_dist,bool verbose, const std::string& mode, const whitelist::wl_entry* wl = nullptr,
        bc_position* pos = nullptr, const std::string* expanded_query = nullptr
    ) {
        // Part-table prefilter: gives the same match set as the full loop below.
        const whitelist::wl_entry::part_idx* parts = nullptr;
        if (!expanded_query && whitelist_type == "true" && wl->true_parts.covers(max_dist, wl->true_bcs.unique_val_size())) {
            parts = &wl->true_parts;
        }
        auto matches = wl->with_wl(whitelist_type, [&](const auto& typed_wl) {
            std::vector<std::pair<int64_seq, int>> match;
            if (verbose) {
            #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << (whitelist_type == "true" ? "EXHAUSTIVE_MUTATION_CHECK" : "EXHAUSTIVE_GLOBAL_MUTATION_CHECK")
                        << " (checking " << typed_wl.size() << " sequences)\n";
                    RAD_COUT << oss.str();
                }
            }

            if (parts) {
                parts->scan(bc, max_dist, match);
                return match;
            }

            auto unique_entries = typed_wl.get_unique_entries();
            for (const auto entry : unique_entries) {
                const auto barcode = entry.barcode();
                int res = barcode_distance(bc, barcode, max_dist, expanded_query);
                if (res >= 0) {
                    match.emplace_back(barcode, res);
                }
            }
            return match;
        });
        
        if (matches.empty()) {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << (whitelist_type == "true" ? "EXHAUSTIVE_MUTATION_CHECK_NO_MATCHES" 
                                                    : "EXHAUSTIVE_GLOBAL_MUTATION_CHECK_NO_MATCHES") << "\n";
                    RAD_COUT << oss.str();
                }
            }
            return std::nullopt;
        }
        
        if (verbose) {
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << (whitelist_type == "true" ? "EXHAUSTIVE_MUTATION_CHECK_FOUND" 
                                                : "EXHAUSTIVE_GLOBAL_MUTATION_CHECK_FOUND")
                    << " (" << matches.size() << " matches)\n";
                RAD_COUT << oss.str();
            }
        }
        
        if (matches.size() == 1) {
            auto [candidate, distance] = matches[0];

            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << (whitelist_type == "true" ? "EXHAUSTIVE_MUTATION_CHECK_CANDIDATE::"
                                                    : "EXHAUSTIVE_GLOBAL_MUTATION_CHECK_CANDIDATE::")
                        << candidate.bits_to_sequence() << "\n";
                    RAD_COUT << oss.str();
                }
            }

            // Quality check
            if (passes_quality_check(candidate, wl, whitelist_type, mode, verbose)) {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "LVDIST::" << distance << "\n"
                            << (whitelist_type == "true" ? "EXHAUSTIVE_MUTATION_CHECK_MATCHED"
                                                        : "EXHAUSTIVE_GLOBAL_MUTATION_CHECK_MATCHED")
                            << "\n";
                        RAD_COUT << oss.str();
                    }
                }
                locate_winner(pos, bc, candidate, max_dist, expanded_query);
                return candidate;
            } else {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << (whitelist_type == "true" ? "EXHAUSTIVE_MUTATION_CHECK_QUALITY_FAILED - rejecting barcode" 
                                                        : "EXHAUSTIVE_GLOBAL_MUTATION_CHECK_QUALITY_FAILED - rejecting barcode") << "\n";
                        RAD_COUT << oss.str();
                    }
                }
                return std::nullopt;
            }
        } else {
            // Multiple matches case - find best match(es),sort by distance (ascending)
            std::sort(matches.begin(), matches.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
            int best_distance = matches[0].second;
            // Collect all matches with the best distance
            std::vector<int64_seq> best_matches;
            for (const auto& [candidate, distance] : matches) {
                if (distance == best_distance) {
                    best_matches.push_back(candidate);
                } else {
                    break; // Since sorted, no more best matches
                }
            }
            
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << (whitelist_type == "true" ? "EXHAUSTIVE_MUTATION_MULTIPLE_FOUND" 
                                                    : "EXHAUSTIVE_GLOBAL_MUTATION_MULTIPLE_FOUND")
                        << " (" << best_matches.size() << " at distance " << best_distance << ")\n";
                    RAD_COUT << oss.str();
                }
            }
            
            if (best_matches.size() == 1) {
                // Single best match
                int64_seq best_candidate = best_matches[0];
                
            if(passes_quality_check(best_candidate, wl, whitelist_type, mode, verbose)){
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "LVDIST::" << best_distance << "\n"
                                << (whitelist_type == "true" ? "EXHAUSTIVE_MUTATION_MULTIPLE_MATCHED_RESOLVED"
                                                            : "EXHAUSTIVE_GLOBAL_MUTATION_MULTIPLE_MATCHED_RESOLVED") << "\n";
                            RAD_COUT << oss.str();
                        }
                    }
                    locate_winner(pos, bc, best_candidate, max_dist, expanded_query);
                    return best_candidate;
                } else {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << (whitelist_type == "true" ? "EXHAUSTIVE_MUTATION_MULTIPLE_RESOLVED_QUALITY_FAILED - rejecting barcode" 
                                                            : "EXHAUSTIVE_GLOBAL_MUTATION_MULTIPLE_RESOLVED_QUALITY_FAILED - rejecting barcode") << "\n";
                            RAD_COUT << oss.str();
                        }
                    }
                    return std::nullopt;
                }
            } else {
                // Multiple equally good matches - ambiguous
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << (whitelist_type == "true" ? "EXHAUSTIVE_MUTATION_MULTIPLE_AMBIGUOUS" 
                                                        : "EXHAUSTIVE_GLOBAL_MUTATION_MULTIPLE_AMBIGUOUS")
                            << " (" << best_matches.size() << " equally good matches)\n";
                        RAD_COUT << oss.str();
                    }
                }
                return std::nullopt;
            }
        }
        return std::nullopt;
    }

/**
 * @brief Perform k-mer based fuzzy search for barcode correction
 * @param original_barcode ``int64_seq`` Original barcode sequence
 * @param expanded_seq ``std::string`` Expanded sequence region to generate k-mers from
 * @param bc_len ``int`` Length of the barcode
 * @param mode ``std::string`` Barcode correction mode ("defensive" or "offensive")
 * @param wl ``const whitelist::wl_entry&`` Reference to the whitelist entry
 * @param verbose ``bool`` Whether to print verbose output
 * @param max_dist ``int`` Maximum allowed edit distance
 * @return ``std::optional<int64_seq>`` Corrected barcode sequence if found, otherwise std::nullopt
 * 
 * @brief This function implements a k-mer based fuzzy search strategy for barcode correction.
 * It generates k-mers from the provided expanded sequence region and checks them against the whitelist
 * using existing barcode correction logic. The search is performed in either "defensive" or "offensive" mode,
 * determining the order of whitelist checks. If no direct k-mer matches are found, the function attempts
 * to find matches through k-mer mutations using an exhaustive search approach.
 */
    std::optional<int64_seq> kmer_fuzzy_search(const int64_seq& original_barcode, const std::string& expanded_seq, int bc_len,
        const std::string& mode, const whitelist::wl_entry& wl, bool verbose, int max_dist,
        bc_position* pos = nullptr
    ) {
    // A direct k-mer hit sets no position here: correct_barcode finds it from the k-mer offsets.
    if (pos) *pos = bc_position{};

    if (verbose) {
        #pragma omp critical
        {
            RAD_COUT << "[kmer_fuzzy_wl_search] Original: " << original_barcode.bits_to_sequence() << std::endl;
            RAD_COUT << "[kmer_fuzzy_wl_search] Expanded region: " << expanded_seq << std::endl;
        }
    }
    
    // Step 1: Generate k-mers from expanded sequence
    std::vector<std::string> kmer_strings = seq_utils::kmerize(expanded_seq, bc_len);
    std::unordered_set<int64_seq> all_candidates;

    for (const auto& kmer_str : kmer_strings) {
        int64_seq kmer;
        kmer.sequence_to_bits(kmer_str);
        if (kmer.is_valid()) {
            all_candidates.insert(kmer);
        }
    }
    
    if (verbose) {
        #pragma omp critical
        {
            RAD_COUT << "[kmer_fuzzy_wl_search] Generated " << all_candidates.size() << " k-mers" << std::endl;
        }
    }
    
    // Step 2: Try direct k-mer matches using existing check_against_wl logic
    if (mode == "defensive") {
        // Check global first, then true
        auto global_result = check_against_wl(original_barcode, all_candidates, "global", max_dist, verbose, mode, &wl);
        if (global_result.has_value()) {
            if (verbose) {
                #pragma omp critical
                {
                    RAD_COUT << "[kmer_fuzzy_wl_search] KMER_DIRECT_HIT (global): " 
                              << global_result.value().bits_to_sequence() << std::endl;
                }
            }
            return global_result;
        }
        auto true_result = check_against_wl(original_barcode, all_candidates, "true", max_dist, verbose, mode, &wl);
        if (true_result.has_value()) {
            if (verbose) {
                #pragma omp critical
                {
                    RAD_COUT << "[kmer_fuzzy_wl_search] KMER_DIRECT_HIT (true): " 
                              << true_result.value().bits_to_sequence() << std::endl;
                }
            }
            return true_result;
        }
    } else { // offensive mode
        // Check true first, then global
        auto true_result = check_against_wl(original_barcode, all_candidates, "true", 2, verbose, mode, &wl);
        if (true_result.has_value()) {
            if (verbose) {
                #pragma omp critical
                {
                    RAD_COUT << "[kmer_fuzzy_wl_search] KMER_DIRECT_HIT (true): " 
                              << true_result.value().bits_to_sequence() << std::endl;
                }
            }
            return true_result;
        }

        auto global_result = check_against_wl(original_barcode, all_candidates, "global", max_dist, verbose, mode, &wl);
        if (global_result.has_value()) {
            if (verbose) {
                #pragma omp critical
                {
                    RAD_COUT << "[kmer_fuzzy_wl_search] KMER_DIRECT_HIT (global): " 
                              << global_result.value().bits_to_sequence() << std::endl;
                }
            }
            return global_result;
        }
    }
    
    if (verbose) {
        #pragma omp critical
        {
            RAD_COUT << "[kmer_fuzzy_wl_search] No direct k-mer hits, trying k-mer mutations..." << std::endl;
        }
    }
    constexpr size_t true_seed_threshold = 10000;
    const size_t true_unique_size = wl.true_bcs.unique_val_size();

    // Preserve the exhaustive-only path used for smaller selected-cell
    // whitelists.  For larger whitelists, use the true-barcode seed index as
    // a fast first attempt, then fall through to the same exhaustive scorer
    // below whenever the seed shortlist does not produce an accepted call.
    if (true_unique_size > true_seed_threshold && wl.true_seeds.ready) {
        auto seed_entries = wl.query_bc_seeds(original_barcode, "true");
        if (verbose) {
            #pragma omp critical
            {
                RAD_COUT << "[kmer_fuzzy_wl_search] Seed shortlist candidates: "
                          << seed_entries.size() << "\n";
            }
        }

        if (!seed_entries.empty()) {
            std::unordered_set<int64_seq> seed_candidates;
            seed_candidates.reserve(seed_entries.size());
            for (const auto e : seed_entries) {
                if (e) seed_candidates.insert(e.barcode());
            }

            auto seed_result = check_against_wl(
                original_barcode, seed_candidates, "true", 2, verbose, mode, &wl
            );
            if (seed_result.has_value()) {
                if (verbose) {
                    #pragma omp critical
                    {
                        RAD_COUT << "[kmer_fuzzy_wl_search] SEED_SHORTLIST_HIT: "
                                  << seed_result.value().bits_to_sequence() << "\n";
                    }
                }
                if (pos) pos->kind = bc_position::window;
                return seed_result;
            }
        }
    }

    // int64_seq intentionally stores at most 64 bases.  The expanded search
    // window may be wider (the barcode plus flanking padding), so retain the
    // original packed barcode for the exhaustive scorer when that happens.
    // Sliding k-mers above still examine every barcode-sized window.
    int64_seq exhaustive_query = original_barcode;
    if (expanded_seq.size() <= 64) {
        exhaustive_query.sequence_to_bits(expanded_seq);
    }
    auto true_result = exhaustive_check_against_wl(
        exhaustive_query, "true", 2, verbose, mode, &wl, pos, expanded_seq.size() > 64 ? &expanded_seq : nullptr);
    if (true_result.has_value()) {
        if (verbose) {
            #pragma omp critical
            {
                RAD_COUT << "[kmer_fuzzy_wl_search] KMER_MUTATION_HIT (true): "
                          << true_result.value().bits_to_sequence() << std::endl;
            }
        }
        return true_result;
    }
    if (verbose) {
        #pragma omp critical
        {
            RAD_COUT << "[kmer_fuzzy_wl_search] All strategies failed" << std::endl;
        }
    }
    return std::nullopt;
}
/**
 * @brief Whitelist lookup of a barcode (the body of correct_barcode)
 * @param elem ``const seq_element&`` Sequence element containing barcode information
 * @param layout ``const ReadLayout&`` Read layout containing whitelist mappings
 * @param full_read ``const read_streaming::sequence&`` Full read sequence
 * @param expanded_seq ``const std::string&`` oriented expanded window (see correct_barcode)
 * @param verbose ``bool`` Whether to print verbose output
 * @param mut_dist ``int`` Maximum allowed edit distance for mutation checks
 * @param mode ``std::string`` Barcode correction mode ("defensive" or "offensive")
 * @param pos optional out: where the lookup placed the barcode
 * @return ``std::optional<int64_seq>`` Corrected barcode sequence if found, otherwise std::nullopt
 */
    std::optional<int64_seq> correct_barcode_lookup(const seq_element& elem,  const ReadLayout& layout,
        const read_streaming::sequence& full_read, const std::string& expanded_seq, bool verbose, int mut_dist, std::string mode,
        bc_position* pos = nullptr) {
        // pick the right whitelist
        auto key = seq_utils::remove_rc(elem.class_id);
        auto wl_it = layout.wl_map.maps.find(key);
        if (wl_it == layout.wl_map.maps.end()) return std::nullopt;
        auto &wl = wl_it->second.get();
        int max_dist = expanded_window_pad;
        // extract and reverse‐complement the raw string
        // encoded barcode and reverse complement
        std::string raw = elem.seq.value();

        if (elem.direction == "reverse") {
            raw = seq_utils::revcomp(raw);
        }

        int64_seq bc, rc_bc;

        bc.sequence_to_bits(raw);
        rc_bc.sequence_to_bits(seq_utils::revcomp(raw));

        // Keep the padded region only when it fits the packed barcode
        // representation.  For a 57-64 base barcode, padding can make this
        // window wider than 64 bases; the fuzzy k-mer pass still scans that
        // full string, while packed fallback comparisons use the barcode.
        int64_seq correction_query = bc;
        if (expanded_seq.size() <= 64) {
            correction_query.sequence_to_bits(expanded_seq);
        }

        int bc_len = static_cast<int>(bc.length);
        int hp_threshold = std::max(4, static_cast<int>(bc_len * 0.4)); // 40% of barcode length, minimum 4

        // === filtering for messy barcodes ===
        bool filtered_hit = !wl.filter_bcs.empty() && (wl.filter_bcs.check_wl_for(bc) || wl.filter_bcs.check_wl_for(rc_bc));
        // A sequence can legitimately be both a barcode and an exact k-mer
        // from a static layout element. Only an identity mapping in the
        // trusted whitelist may override that collision; mutation aliases and
        // global-whitelist-only hits remain filtered.
        bool bc_true_identity =
            filtered_hit && wl.true_bcs.has_identity_mapping(bc);
        bool rc_true_identity =
            filtered_hit && wl.true_bcs.has_identity_mapping(rc_bc);
        bool exact_true_identity =
            bc_true_identity || rc_true_identity;
        bool seq_hit = !wl.filter_bcs.empty() && (seq_utils::int_kmerize(raw, 2) < 4 || mutation_tools::detect_hp(raw, hp_threshold));
        if ((filtered_hit && !exact_true_identity) || seq_hit) {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "FILTER_CHECK_FOUND\nNO_CHECK_WORKED\n";
                    RAD_COUT << oss.str();
                }
            }
            return std::nullopt;
        }
        if (exact_true_identity) {
            if (verbose) {
                #pragma omp critical
                {
                    RAD_COUT
                        << "FILTER_CHECK_EXACT_TRUE_OVERRIDE ("
                        << (bc_true_identity ? "direct" : "reverse-complement")
                        << ")\n";
                }
            }
            if (pos) pos->kind = bc_position::window;
            return bc_true_identity ? bc : rc_bc;
        }

        // === exact match in global whitelist ===
        if (!wl.global_bcs.empty() && (wl.global_bcs.check_wl_for(bc) || wl.global_bcs.check_wl_for(rc_bc))) {
            auto matched = wl.global_bcs.return_putative_correct_bcs(bc);
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "GLOBAL_CHECK_FOUND (" << matched.size() << " candidates)\n";
                    RAD_COUT << oss.str();
                }
            }
            
            if (matched.size() == 1) {
                int64_seq candidate = *matched.begin();
                if (passes_quality_check(candidate, &wl, "global", mode, verbose) || bc == candidate || rc_bc == candidate) {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "GLOBAL_CHECK_WORKED\n";
                            RAD_COUT << oss.str();
                        }
                    }
                    if (pos) pos->kind = bc_position::window;
                    return candidate;
                } else {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "GLOBAL_CHECK_QUALITY_FAILED - falling through to true whitelist\n";
                            RAD_COUT << oss.str();
                        }
                    }
                    // Fall through to true whitelist check
                }
            } else {
                // Use enhanced resolution with count-based tie breaking for global
                auto [resolved, min_dist] = resolve_multiple_hits_simple(bc, matched, max_dist, verbose, "global", &wl);
                if (resolved.has_value()) {
                    // Quality check the resolved candidate
                    if (passes_quality_check(resolved.value(), &wl, "global", mode, verbose)) {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "GLOBAL_CHECK_MULTIPLE_RESOLVED\n";
                                RAD_COUT << oss.str();
                            }
                        }
                        if (pos) pos->kind = bc_position::window;
                        return resolved;
                    } else {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "GLOBAL_CHECK_MULTIPLE_RESOLVED_QUALITY_FAILED - falling through to true whitelist\n";
                                RAD_COUT << oss.str();
                            }
                        }
                        // Fall through to true whitelist check
                    }
                } else {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "GLOBAL_CHECK_MULTIPLE_UNRESOLVED";
                            if (min_dist.has_value()) {
                                oss << " (min_dist=" << min_dist.value() << ")";
                            }
                            oss << " - falling through to true whitelist\n";
                            RAD_COUT << oss.str();
                        }
                    }
                    // Fall through to true whitelist check
                }
            }
        }
        
        // === exact match in true barcodes ===
        bool found_forw_bc = wl.true_bcs.check_wl_for(bc);
        bool found_rev_bc = wl.true_bcs.check_wl_for(rc_bc);
        if (!wl.true_bcs.empty() && (found_forw_bc || found_rev_bc)) {
            std::unordered_set<int64_seq> matched;
            if(found_forw_bc){
                matched = wl.true_bcs.return_putative_correct_bcs(bc);
            } else {
                matched = wl.true_bcs.return_putative_correct_bcs(rc_bc);
            }
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "ORIGINAL_CHECK_FOUND (" << matched.size() << " candidates)\n";
                    RAD_COUT << oss.str();
                }
            }
            
            if (matched.size() == 1) {
                int64_seq candidate = *matched.begin();
                // QC: if this fails on true whitelist, reject the barcode entirely
                if (passes_quality_check(candidate, &wl, "true", mode, verbose) || bc == candidate || rc_bc == candidate) {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "ORIGINAL_CHECK_WORKED\n";
                            RAD_COUT << oss.str();
                        }
                    }
                    if (pos) pos->kind = bc_position::window;
                    return candidate;
                } else {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "ORIGINAL_CHECK_QUALITY_FAILED - rejecting barcode\n";
                            RAD_COUT << oss.str();
                        }
                    }
                    return std::nullopt; // Fail the barcode
                }
            } else {
                // Use enhanced resolution with count-based tie breaking for true
                auto [resolved, min_dist] = resolve_multiple_hits_simple(bc, matched, 2, verbose, "true", &wl);
                if (resolved.has_value()) {
                    // Quality check the resolved candidate
                    if (passes_quality_check(resolved.value(), &wl, "true", mode, verbose)) {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "ORIGINAL_MATCH_COLLISION_RESOLVED\n";
                                RAD_COUT << oss.str();
                            }
                        }
                        if (pos) pos->kind = bc_position::window;
                        return resolved;
                    } else {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "ORIGINAL_MATCH_COLLISION_RESOLVED_QUALITY_FAILED - rejecting barcode\n";
                                RAD_COUT << oss.str();
                            }
                        }
                        return std::nullopt; // Fail the barcode
                    }
                } else {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "ORIGINAL_MATCH_COLLISION_UNRESOLVED";
                            if (min_dist.has_value()) {
                                oss << " (min_dist=" << min_dist.value() << ")";
                            }
                            oss << "\n";
                            RAD_COUT << oss.str();
                        }
                    }
                    return std::nullopt; // Fail the barcode
                }
            }
        }

        if (verbose) {
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "ORIGINAL_MATCH_NOT_FOUND\n";
                RAD_COUT << oss.str();
            }
        }
        // generate mutations
        // === k-mer fuzzy search ===

      auto kmer_fuzzy_result = kmer_fuzzy_search(bc, expanded_seq, bc_len, mode, wl, verbose, 3, pos);
      if (verbose) {
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "KMER_FUZZY_SEARCH_RESULT: " 
                    << (kmer_fuzzy_result.has_value() ? kmer_fuzzy_result->bits_to_sequence() : "NO_MATCH") 
                    << "\n";
                RAD_COUT << oss.str();
            }
        }
        if (kmer_fuzzy_result.has_value()) {
            return kmer_fuzzy_result;
        }
        
        if (pos) *pos = bc_position{};
        const auto check_true_mutations = [&]() -> std::optional<int64_seq> {
            if (wl.true_bcs.empty()) return std::nullopt;
            auto muts = mutation_tools::generate_mutated_barcodes(bc, mut_dist);
            return check_against_wl(correction_query, muts, "true", 2, verbose, mode, &wl, pos, expanded_seq.size() > 64 ? &expanded_seq : nullptr);
        };
        const auto check_global_mutations = [&]() -> std::optional<int64_seq> {
            if (wl.global_bcs.empty()) return std::nullopt;
            auto global_muts = mutation_tools::generate_lv_candidates(bc, mut_dist, wl.global_bcs.length_mask());
            return check_against_wl(correction_query, global_muts, "global", max_dist, verbose, mode, &wl, pos, expanded_seq.size() > 64 ? &expanded_seq : nullptr);
        };
        if(mode == "defensive"){
            auto global_result = check_global_mutations();
            if(global_result.has_value()){
                return(global_result);
            }

            auto true_result = check_true_mutations();
            if(true_result.has_value()){
                return(true_result);
            }
        }

        // IF OFFENSIVE: Check against true first, and then global
        
        if(mode == "offensive") {
            auto true_result = check_true_mutations();
            if (true_result.has_value()) {
                if (verbose) {
                    #pragma omp critical
                    {
                        RAD_COUT << "[kmer_fuzzy_wl_search] KMER_DIRECT_HIT (true): " 
                                << true_result.value().bits_to_sequence() << std::endl;
                    }
                }
                return true_result;
            }

            auto global_result = check_global_mutations();
            if (global_result.has_value()) {
                if (verbose) {
                    #pragma omp critical
                    {
                        RAD_COUT << "[kmer_fuzzy_wl_search] KMER_DIRECT_HIT (global): " 
                                << global_result.value().bits_to_sequence() << std::endl;
                    }
                }
                return global_result;
            }
        
        }

        // === no match ===
        if (verbose) {
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "[kmer_fuzzy_wl_search] NO_MATCH_FOUND\n";
                RAD_COUT << oss.str();
            }
        }
        return std::nullopt;
    }

/**
 * @brief Correct a barcode sequence using the provided layout and read information
 * @param elem ``const seq_element&`` Sequence element containing barcode information
 * @param layout ``const ReadLayout&`` Read layout containing whitelist mappings
 * @param full_read ``const read_streaming::sequence&`` Full read sequence
 * @param verbose ``bool`` Whether to print verbose output
 * @param mut_dist ``int`` Maximum allowed edit distance for mutation checks
 * @param mode ``std::string`` Barcode correction mode ("defensive" or "offensive")
 * @param pos optional out: where the corrected barcode lies in the oriented expanded window
 * @return ``std::optional<int64_seq>`` Corrected barcode sequence if found, otherwise std::nullopt
 */
    std::optional<int64_seq> correct_barcode(const seq_element& elem,  const ReadLayout& layout,
        const read_streaming::sequence& full_read,  bool verbose, int mut_dist, std::string mode,
        bc_position* pos = nullptr) {
        const int read_len = static_cast<int>(full_read.seq.size());
        const int win_lo = std::max(1, elem.position.first - expanded_window_pad);
        const int win_hi = std::min(read_len, elem.position.second + expanded_window_pad);
        std::string expanded_seq = seq_utils::substr_w_padding(full_read.seq, elem.position.first, elem.position.second, expanded_window_pad);
        const bool reverse = elem.direction == "reverse";
        if (reverse) {
            expanded_seq = seq_utils::revcomp(expanded_seq);
        }
        if (pos) *pos = bc_position{};

        auto result = correct_barcode_lookup(elem, layout, full_read, expanded_seq, verbose, mut_dist, mode, pos);
        if (!pos || !result.has_value()) return result;

        // Search the full string: a 64-base barcode plus padding can exceed
        // the packed representation without exceeding the supported barcode length.
        const int W = static_cast<int>(expanded_seq.size());
        const std::string winner = result->bits_to_sequence();
        int copies = 0, copy_end = -1;
        for (size_t found = expanded_seq.find(winner); found != std::string::npos;
             found = expanded_seq.find(winner, found + 1)) {
            ++copies;
            copy_end = static_cast<int>(found + winner.size() - 1);
        }
        if (copies >= 1) pos->set_from_match(copy_end, copies);
        else if (pos->kind == bc_position::located && pos->end == W - 1)
            pos->kind = bc_position::clipped;
        if (pos->kind == bc_position::located) {
            const int window_first = reverse ? (win_hi - elem.position.second) : (elem.position.first - win_lo);
            const int window_last = window_first + (elem.position.second - elem.position.first);
            pos->shift = pos->end - window_last;
        }
        return result;
    }
};

// Run counters of single-barcode position outcomes; reset by sigalign.
namespace barcode_position_stats {
inline std::atomic<size_t> moved{0}, moved_forward{0}, moved_reverse{0}, in_place{0}, window{0}, unknown{0}, tie{0},
    clipped{0}, out_of_read{0};
inline void reset() {
    for (auto* c : {&moved, &moved_forward, &moved_reverse, &in_place, &window, &unknown, &tie, &clipped, &out_of_read})
        c->store(0);
}
}  // namespace barcode_position_stats

struct sigalign_run_stats {
    size_t total_reads = 0;
    size_t reads_passing_filter = 0;
    size_t reads_demultiplexed = 0;
    size_t records_serialized = 0;
    size_t chunks_processed = 0;
    double wall_time_seconds = 0.0;
    double process_time_seconds = 0.0;
    double output_staging_time_seconds = 0.0;
    double overhead_time_seconds = 0.0;
    size_t bc_pos_moved = 0, bc_pos_moved_forward = 0, bc_pos_moved_reverse = 0, bc_pos_in_place = 0, bc_pos_window = 0,
           bc_pos_unknown = 0, bc_pos_tie = 0, bc_pos_clipped = 0, bc_pos_out_of_read = 0;
    bool concat_hmm_enabled = false;
    size_t hmm_reads = 0, hmm_single = 0, hmm_split = 0, hmm_children = 0, hmm_children_full_search = 0,
           hmm_same_molecule = 0;
    size_t hmm_legacy_guard = 0, hmm_legacy_k0 = 0, hmm_legacy_abstain = 0, hmm_legacy_artifact = 0,
           hmm_legacy_too_many = 0, hmm_full_as_missing = 0;
    size_t hmm_spacer_retry = 0, hmm_spacer_retry_ok = 0, hmm_fragment_ok = 0, hmm_read_end_ok = 0;
    bool hmm_abstain_legacy = false;
    size_t hmm_drop_abstain = 0, hmm_drop_too_many = 0, hmm_drop_d_reads = 0, hmm_drop_d_children = 0;
    size_t hmm_d_kept_left = 0, hmm_d_kept_right = 0, hmm_d_split = 0;
    size_t hmm_fold_read_start = 0, hmm_fold_read_end = 0, hmm_partner_rejected = 0, hmm_retry_from_partner = 0;
    size_t hmm_art_reads = 0, hmm_pieces = 0, hmm_trim_keep_F = 0, hmm_trim_keep_R = 0;
    size_t hmm_both_pass_dropped = 0, hmm_duplet_named = 0;
    size_t hmm_clip_hits = 0, hmm_two_unit_clip_dropped = 0;
    size_t hmm_junctions_abstained = 0, hmm_pieces_suppressed = 0, hmm_drop_unresolved = 0;
    size_t hmm_boundary_unresolved_parents = 0, hmm_boundary_unresolved_children = 0;
    size_t hmm_cut_kind[6] = {0, 0, 0, 0, 0, 0};
    size_t hmm_win[4] = {0, 0, 0, 0}, hmm_win_ok[4] = {0, 0, 0, 0};  // per concat_hmm::Status
};

#ifdef RAD_STAGE_TIMERS
namespace rad_stage_timers {
inline std::atomic<unsigned long long> read_ns{0}, hmm_ns{0}, windowed_static_ns{0}, legacy_static_ns{0}, reads{0};
inline std::atomic<unsigned long long> variable_ns{0}, filter_ns{0}, molecules{0}, molecules_passed{0};
// [0 HMM route, 1 default path][0 passed, 1 barcode correction failed, 2 other]
inline std::atomic<unsigned long long> fb_ns[2][3]{}, fb_n[2][3]{};
struct scope {
    std::atomic<unsigned long long>& acc;
    std::chrono::steady_clock::time_point t0;
    explicit scope(std::atomic<unsigned long long>& a) : acc(a), t0(std::chrono::steady_clock::now()) {}
    ~scope() {
        acc.fetch_add(static_cast<unsigned long long>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                          std::chrono::steady_clock::now() - t0).count()), std::memory_order_relaxed);
    }
};
}  // namespace rad_stage_timers
#endif

// A bounded molecule and its assigned static evidence, before any variable
// extraction or barcode-count side effects. Coordinates are parent-relative.
struct static_segment {
    std::pair<int, int> position;
    std::string direction;
    std::vector<seq_element> elements;
    int edit_distance = 0;
    bool ambiguous = false;
};

// --concat-hmm: concat_hmm::segment() names the strand and the expected static elements of each
// read or child, with a window and status for each; sigalign_static aligns only inside these windows.
struct static_check_window {
    std::string class_id;
    int lo = 1, hi = 0;  // 1-based inclusive, child coordinates
    concat_hmm::Status status = concat_hmm::Status::MISSING_EXPECTED;
    int retained_from = -1, retained_to = -1, edits = -1;
};

// Layout facts used by the windowed aligners, built once per run.
struct concat_elem_info {
    int m = 0;                   // adapter length
    bool left_outward = false;   // no payload variable before it
    bool right_outward = false;  // no payload variable after it
    int k_full = 4;
    int k_missing = 5;
    std::string poly_id;         // poly tail joined by fixed-length spacers only
    int poly_gap_lo = 0, poly_gap_hi = -1;
    bool poly_after = true;
    std::string partner_id;      // adapter joined by fixed-length spacers only
    int partner_gap_lo = 0, partner_gap_hi = -1;
    bool partner_after = true;
    bool bc_outer = false;       // opens a barcode block that partner_id closes
};

struct concat_hmm_counters {
    std::atomic<size_t> reads{0}, single{0}, split{0}, children{0}, children_full_search{0}, same_molecule{0};
    std::atomic<size_t> legacy_guard{0}, legacy_k0{0}, legacy_abstain{0}, legacy_artifact{0}, legacy_too_many{0};
    std::atomic<size_t> win[4]{}, win_ok[4]{};  // per concat_hmm::Status
    std::atomic<size_t> full_as_missing{0};
    std::atomic<size_t> spacer_retry{0}, spacer_retry_ok{0};
    std::atomic<size_t> fragment_ok{0};
    std::atomic<size_t> read_end_ok{0};
    std::atomic<size_t> drop_abstain{0}, drop_too_many{0};
    std::atomic<size_t> drop_d_reads{0}, drop_d_children{0};
    std::atomic<size_t> d_kept_left{0}, d_kept_right{0}, d_split{0};
    std::atomic<size_t> fold_read_start{0}, fold_read_end{0};
    std::atomic<size_t> partner_rejected{0};
    std::atomic<size_t> retry_from_partner{0};
    std::atomic<size_t> art_reads{0};
    std::atomic<size_t> pieces{0};
    std::atomic<size_t> trim_keep_F{0}, trim_keep_R{0};
    std::atomic<size_t> both_pass_dropped{0};
    std::atomic<size_t> clip_hits{0};
    std::atomic<size_t> two_unit_clip_dropped{0};
    std::atomic<size_t> duplet_named{0};
    std::atomic<size_t> junctions_abstained{0};
    std::atomic<size_t> pieces_suppressed{0};
    std::atomic<size_t> drop_unresolved{0};
    std::atomic<size_t> boundary_unresolved_parents{0}, boundary_unresolved_children{0};
    std::atomic<size_t> cut_kind[6]{};          // per concat_hmm CUT_* kind
};

// Barcode / UMI spacers between an outer adapter (no payload on its far side, e.g. 10x 5' forw_primer) and the
// inner element joined to it by fixed-length spacers only (10x 5' tso; a poly tail on VisiumHD).
struct concat_bc_block {
    std::string outer_id, inner_id;
    bool inner_after = true;
    int gap_lo = 0, gap_hi = -1;         // layout spacer length between them
    std::vector<std::string> barcode_ids;
};

struct concat_layout_info {
    std::vector<concat_bc_block> bc_blocks;
    std::unordered_map<std::string, concat_elem_info> elems;
    std::unordered_map<std::string, const ReadElement*> layout_by_id;
    std::vector<std::string> id_F, id_R;  // HMM spec element index -> RAD class_id per strand
    int seen_pad = 6;
    concat_hmm_counters* ctr = nullptr;
};

struct static_restriction {
    std::string direction;
    std::vector<static_check_window> windows;
    const concat_layout_info* info = nullptr;
    bool read_start = true, read_end = true;   // false at an HMM cut
};

// --concat-hmm: a piece's place in its parent read. Read ends and junction counterparts are judged
// on the parent read, because a piece edge at an HMM cut is not a physical read end.
struct static_piece_frame {
    const std::string* parent_seq = nullptr;
    int offset = 0;                           // 0-based start of the piece in the parent read
};

// A main F / R template, not a T, D or single-primer-end artifact.
inline bool concat_plain_construct(const concat_hmm::Segment& s) {
    return (s.strand == 'F' || s.strand == 'R') && !(s.flags & concat_hmm::SEG_ARTIFACT);
}

inline concat_layout_info build_concat_layout_info(const ReadLayout& layout, const concat_hmm::Model& model) {
    concat_layout_info info;
    info.seen_pad = model.opt.seen_pad;
    std::map<std::string, std::vector<const ReadElement*>> by_dir;
    for (const auto& e : layout.by_order()) by_dir[e.direction].push_back(&e);
    auto lengths = [](const ReadElement& e, int& lo, int& hi) {
        if (!e.length_candidates.empty()) {
            lo = *std::min_element(e.length_candidates.begin(), e.length_candidates.end());
            hi = *std::max_element(e.length_candidates.begin(), e.length_candidates.end());
            return true;
        }
        if (e.expected_length && *e.expected_length > 0) { lo = hi = *e.expected_length; return true; }
        return false;
    };
    for (const auto& dir_elems : by_dir) {
        const auto& v = dir_elems.second;
        for (size_t i = 0; i < v.size(); ++i) {
            const ReadElement& e = *v[i];
            if (e.type != "static" || e.seq.empty() || e.global_class == "start" ||
                e.global_class == "stop" || e.global_class == "poly_tail") continue;
            concat_elem_info x;
            x.m = static_cast<int>(e.seq.size());
            bool var_before = false, var_after = false;
            for (size_t j = 0; j < v.size(); ++j) {
                if (j == i || v[j]->type != "variable") continue;
                if (j < i) var_before = true; else var_after = true;
            }
            x.left_outward = !var_before;
            x.right_outward = !var_after;
            x.k_full = std::max(1, std::min(4, x.m * 4 / 22));
            x.k_missing = std::max(1, std::min(5, x.m * 5 / 22));
            // walk outward over fixed-length variables to the first poly tail or adapter
            auto link = [&](int step, bool want_poly) {
                int lo_sum = 0, hi_sum = 0;
                for (int j = static_cast<int>(i) + step; j >= 0 && j < static_cast<int>(v.size()); j += step) {
                    const ReadElement& o = *v[j];
                    if (o.global_class == "poly_tail") {
                        if (!want_poly) return false;
                        x.poly_id = o.class_id;
                        x.poly_gap_lo = lo_sum;
                        x.poly_gap_hi = hi_sum;
                        x.poly_after = step > 0;
                        return true;
                    }
                    if (o.type == "static" && !o.seq.empty() && o.global_class != "start" && o.global_class != "stop") {
                        if (want_poly || j == static_cast<int>(i) + step) return false;
                        x.partner_id = o.class_id;
                        x.partner_gap_lo = lo_sum;
                        x.partner_gap_hi = hi_sum;
                        x.partner_after = step > 0;
                        return true;
                    }
                    int lo = 0, hi = 0;
                    if (o.type != "variable" || o.global_class == "read" || !lengths(o, lo, hi)) return false;
                    lo_sum += lo;
                    hi_sum += hi;
                }
                return false;
            };
            if (!link(+1, true)) link(-1, true);
            if (!link(+1, false)) link(-1, false);
            info.elems[e.class_id] = x;
            info.layout_by_id[e.class_id] = &e;
            const bool adapter_partner = !x.partner_id.empty();
            const std::string& inner = adapter_partner ? x.partner_id : x.poly_id;
            const bool after = adapter_partner ? x.partner_after : x.poly_after;
            if (!inner.empty() && (after ? x.left_outward : x.right_outward)) {
                concat_bc_block b;
                b.outer_id = e.class_id;
                b.inner_id = inner;
                b.inner_after = after;
                b.gap_lo = adapter_partner ? x.partner_gap_lo : x.poly_gap_lo;
                b.gap_hi = adapter_partner ? x.partner_gap_hi : x.poly_gap_hi;
                for (int j = static_cast<int>(i) + (after ? 1 : -1); j >= 0 && j < static_cast<int>(v.size()); j += after ? 1 : -1) {
                    if (v[j]->class_id == inner) break;
                    if (v[j]->type == "variable" && v[j]->global_class == "barcode") b.barcode_ids.push_back(v[j]->class_id);
                }
                if (!b.barcode_ids.empty()) {
                    if (adapter_partner) info.elems[e.class_id].bc_outer = true;
                    info.bc_blocks.push_back(std::move(b));
                }
            }
        }
    }
    // a forward-only row maps to its reverse-complement counterpart on the R strand
    const auto& els = model.spec.elements;
    info.id_F.assign(els.size(), "");
    info.id_R.assign(els.size(), "");
    for (size_t j = 0; j < els.size(); ++j) {
        const auto& e = els[j];
        const std::string other_dir = e.direction == 'F' ? "reverse" : "forward";
        std::string other;
        for (const auto& o : layout.by_order()) {
            if (o.direction.rfind(other_dir, 0) != 0 || o.type != "static") continue;
            const bool poly = o.global_class == "poly_tail";
            if ((e.klass == "poly_tail" && poly) ||
                (!poly && !e.seq.empty() && o.seq == seq_utils::revcomp(e.seq))) { other = o.class_id; break; }
        }
        (e.direction == 'F' ? info.id_F[j] : info.id_R[j]) = e.id;
        (e.direction == 'F' ? info.id_R[j] : info.id_F[j]) = other;
    }
    return info;
}

// --concat-hmm: a barcode block in a piece after the full static search. Positions are 1-based inclusive
// in piece coordinates, -1 when absent.
struct concat_bc_state {
    const concat_bc_block* block = nullptr;
    std::string dir;
    int outer_s = -1, outer_e = -1;
    int inner_s = -1, inner_e = -1;
    bool paired = false;             // inner at the layout spacer offset from the outer (+-4 bp)
    bool outer_clipped = false;      // outer is a refused clipped hit: evidence, not an element
    bool strict() const { return outer_s > 0 && inner_s > 0 && paired; }
    bool usable() const { return inner_s > 0 && (outer_s <= 0 || paired); }
    bool bare() const { return outer_s > 0 && !paired; }
    int rank() const { return strict() ? 3 : usable() ? 2 : bare() ? 1 : 0; }
};

// What the router does with a T / single-primer-end piece.
struct concat_piece_plan {
    enum kind_t { NORMAL, TRIM } kind = NORMAL;
    std::string dir;                // strand of the only usable barcode unit; empty: none or two
    int trim_lo = 0, trim_hi = 0;   // TRIM: kept range, 0-based half-open, piece coordinates
    bool no_double_count = false;
    uint8_t clip_check = 0;         // see SigString::hmm_clip_check
    const char* note = "";
};

/**
 * @class SigString
 * @brief Container for a sequencing read and its aligned elements with multi-indexed access
 * 
 * SigString represents a single sequencing read along with all its identified and aligned
 * elements (barcodes, UMIs, adapters, etc.). It uses a SigElement multi-index container
 * to efficiently store and query elements by various criteria such as class ID, order,
 * direction, and validation status.
 * 
 * @param sig_elements Multi-index container holding SigElement objects
 * @param sequence_id Identifier for the sequencing read
 * @param sequence_length Length of the sequencing read
 * @param read_type Type of read (forward, reverse, concatenate)
 * @param additional_info Additional information about the read
 */
class SigString {
    SigElement sig_elements;
    std::string sequence_id;
    int sequence_length;
    std::string read_type;
    std::string additional_info;
    bool from_concatemer = false;
    // --concat-hmm only; the defaults are inert on every other path
    uint8_t hmm_fold_edges = 0;   // bit0 child starts at a fold-back cut, bit1 ends at one
    uint8_t hmm_fold_used = 0;
    std::string hmm_fold_dir;
    std::vector<std::string> hmm_reject_barcodes;
    struct hmm_weak_ref { std::string id; std::pair<int, int> aligned, anchored; };
    std::vector<hmm_weak_ref> hmm_weak_refs;
    std::string hmm_note;
    std::string hmm_boundary_note;
    bool hmm_no_double_count = false;
    // SSW hits refused by the read-end clipping rule in a piece's full search: evidence only, never elements.
    struct hmm_clip_hit { std::string id; std::pair<int, int> pos; bool block_adapter = false; };
    std::vector<hmm_clip_hit> hmm_clip_hits;
    uint8_t hmm_clip_check = 0;        // do not write a record whose read element holds a clip hit (bit0 fwd, bit1 rev)
    bool hmm_edge_is_fold = true;
    int hmm_parent_off = -1, hmm_parent_end = -1;  // window in the parent read, 0-based half-open (reporting only)

public:
    SigString(
        std::string id = "", 
        int length = 0, 
        std::string type = "undefined", 
        std::string info = ""
    ) : sequence_id(std::move(id)),
          sequence_length(length),
          read_type(std::move(type)),
          additional_info(std::move(info)) {}

    // Metadata accessors
    const std::string& id() const { 
        return sequence_id; 
    }
    int length() const { 
        return sequence_length; 
    }
    const std::string& type() const { 
        return read_type; 
    }
    const std::string& info() const { 
        return additional_info; 
    }

    // Container access methods
    const SigElement& elements() const { 
        return sig_elements; 
    }
    SigElement& elements() { 
        return sig_elements; 
    }

    // Index accessors
    auto& by_id() { return sig_elements.get<sig_id_tag>(); }
    auto& by_global() { return sig_elements.get<sig_global_tag>(); }
    auto& by_edit_distance() { return sig_elements.get<sig_ed_tag>(); }
    auto& by_order() { return sig_elements.get<sig_order_tag>(); }
    auto& by_direction() { return sig_elements.get<sig_dir_tag>(); }
    auto& by_pass() { return sig_elements.get<sig_pass_tag>(); }
    auto& by_read() { return sig_elements.get<sig_read_tag>(); }

    // Const versions of index accessors
    const auto& by_id() const { return sig_elements.get<sig_id_tag>(); }
    const auto& by_global() const { return sig_elements.get<sig_global_tag>(); }
    const auto& by_edit_distance() const { return sig_elements.get<sig_ed_tag>(); }
    const auto& by_order() const { return sig_elements.get<sig_order_tag>(); }
    const auto& by_direction() const { return sig_elements.get<sig_dir_tag>(); }
    const auto& by_pass() const { return sig_elements.get<sig_pass_tag>(); }
    const auto& by_read() const { return sig_elements.get<sig_read_tag>(); }

    // Element manipulation
    void add_element(const seq_element& element) {
        sig_elements.insert(element);
    }
    
    void add_element(seq_element&& element) {
        sig_elements.insert(std::move(element));
    }

    // Metadata setters
    void set_id(const std::string& id) { sequence_id = id; }
    void set_length(int length) { sequence_length = length; }
    void set_type(const std::string& type) { read_type = type; }
    void set_info(const std::string& info) { additional_info = info; }
    void mark_concatemer() { from_concatemer = true; }
    // --concat-hmm: the read element may start (end) at the child boundary
    void set_fold_edges(bool at_start, bool at_end, const std::string& dir, bool fold = true) {
        hmm_fold_edges = static_cast<uint8_t>((at_start ? 1 : 0) | (at_end ? 2 : 0));
        hmm_fold_dir = dir;
        hmm_edge_is_fold = fold;
    }
    uint8_t fold_boundary_used() const { return hmm_fold_used; }
    bool fold_edge_is_fold() const { return hmm_edge_is_fold; }
    // --concat-hmm: debug sigstrings print in parent coordinates; the mapping frame and FASTQ records do not change.
    void set_parent_frame(int off, int end) { hmm_parent_off = off; hmm_parent_end = end; }
    int parent_offset() const { return hmm_parent_off; }

    // Container operations
    size_t size() const { return sig_elements.size(); }
    bool empty() const { return sig_elements.empty(); }
    void clear() { 
        sig_elements.clear(); 
    }

    // Iterators
    auto begin() { return sig_elements.begin(); }
    auto end() { return sig_elements.end(); }
    auto begin() const { return sig_elements.begin(); }
    auto end() const { return sig_elements.end(); }

    template<typename mod> bool edit_elem(const std::string &class_id, mod m) {
      auto &idx = sig_elements.get<sig_id_tag>();
      auto it = idx.find(class_id);
      if (it == idx.end()){
        return false;
      }
      idx.modify(it, m);
      return true;
    }

private:
/**
 * @brief Generate variable elements in the SigElement container
 * @param layout_elem `pointer` to the ReadElement defining the variable element layout
 * @param static_refs multimap of static reference elements
 * @param read_seq sequencing read 
 * @param sig_elements sig_elements to look at
 * @return true if variable elements were successfully generated and embedded into the order, false otherwise
 */
    bool generate_variable_elements(
        const ReadElement* layout_elem, const std::multimap<std::string, const seq_element*>& static_refs,
                                  const std::string& read_seq, SigElement& sig_elements, bool verbose
                                ) {
        auto& id_index = sig_elements.get<sig_id_tag>();
        auto it = id_index.find(layout_elem->class_id);
        if (it == id_index.end() || !layout_elem->ref_pos){
            return false;
        }
        const auto& ref_pos = *layout_elem->ref_pos;
        std::pair<int, int> var_positions = {-1, -1};

        if (verbose) {
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "Processing variable element: " << layout_elem->class_id << std::endl;
                RAD_COUT << oss.str();
            }
        }
        return map_positions(ref_pos, static_refs, read_seq, it, id_index, var_positions, verbose);
    }

/**
 * @brief Validate variable element positions within the read length by checking boundaries
 * @param positions `std::pair<int, int>` where `std::first` is element start and `std::second` is element end
 * @param read_length read length, `size_t`
 * @return true if positions are valid, false otherwise
 */
    bool validate_var_positions(const std::pair<int, int>& positions, size_t read_length) {
        return positions.first > 0 && 
               positions.second > 0 && 
               positions.first <= static_cast<int>(read_length) &&
               positions.second <= static_cast<int>(read_length) &&
               positions.second > positions.first;
    }

/**
 * @brief --concat-hmm: child-boundary position for a read element whose references on that side are not placed.
 * @return the position next to the virtual start (stop) element, or -1
 */
    static bool hmm_ref_placed(const std::vector<const seq_element*>& ordered, const std::string& id) {
        if (id.empty()) return false;
        for (const seq_element* e : ordered)
            if (e->class_id == id && e->position.first > 0) return true;
        return false;
    }
    int hmm_fold_boundary(const std::string& ref_a, const std::string& ref_b, const std::vector<const seq_element*>& ordered,
                          const seq_element& var, bool at_start) {
        if (var.global_class != "read" || var.direction != hmm_fold_dir) return -1;
        if (hmm_ref_placed(ordered, ref_a) || hmm_ref_placed(ordered, ref_b)) return -1;
        for (const seq_element* e : ordered) {
            if (e->type != "static" || e->global_class != (at_start ? "start" : "stop")) continue;
            hmm_fold_used = static_cast<uint8_t>(hmm_fold_used | (at_start ? 1 : 2));
            return at_start ? e->position.second + 1 : e->position.first - 1;
        }
        return -1;
    }

/**
 * @brief Map variable element positions based on reference positions and static elements. 
 * Incredibly bulky and complex, but works. Future refactor necessary.
 * @param ref_pos `ReferencePositions` object containing primary and secondary reference positions
 * @param static_refs `std::multimap <std::string, const seq_element*>` of static reference elements
 * @param read_seq `std::string` sequenced read
 * @param var_it iterator to the variable element in the SigElement container
 * @param id_index index of SigElement container by sig_id_tag
 * @param var_positions modifiable `std::pair<int, int>` pair to store calculated variable element positions
 * @param verbose print verbose output
 * @return true if mapping was successful and positions were set, false otherwise
 * 
 * @brief This function maps the positions of a variable element within a sequencing read
 * based on provided reference positions and a set of static reference elements. It first
 * constructs an ordered list of static elements that share the same direction as the variable element.
 * It then attempts to determine the start position of the variable element using primary and secondary
 * reference positions. If successful, it calculates the end position based on the length of the
 * variable element's sequence. The function ensures that the calculated positions are valid
 * within the bounds of the read length before updating the provided `var_positions` pair.
 */
    bool map_positions(
        const ReferencePositions& ref_pos, 
        const std::multimap<std::string, const seq_element*>& static_refs,
        const std::string& read_seq,  
        SigElement::index<sig_id_tag>::type::iterator var_it,
        SigElement::index<sig_id_tag>::type& id_index,  
        std::pair<int, int>& var_positions,
        bool verbose
    ) {
    // Build an ordered list of static elements with the same direction as the variable element.
    std::vector<const seq_element*> ordered_elements;
    for (const auto& [key, elem] : static_refs) {
        if (elem->direction == var_it->direction) {
            ordered_elements.push_back(elem);
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Added static reference " << elem->class_id 
                              << " (order: " << elem->order 
                              << ", pos: " << elem->position.first << "-" << elem->position.second 
                              << ") to ordered list." << std::endl;
                    RAD_COUT << oss.str();
                }
            }
        }
    }
    std::sort(ordered_elements.begin(), ordered_elements.end(),
    [](const seq_element* a, const seq_element* b) {
         if (a->position.first != b->position.first)
              return a->position.first < b->position.first;
         return a->order < b->order;
    });
    
    if (verbose) {
        #pragma omp critical
        {
            std::ostringstream oss;
            oss << "Ordered static references for variable element " 
                      << var_it->class_id << ": ";
            RAD_COUT << oss.str();
            for (const auto& elem : ordered_elements) {
                std::ostringstream oss_elem;
                oss_elem << elem->class_id << "(" << elem->order << ") ";
                RAD_COUT << oss_elem.str();
            }
            RAD_COUT << std::endl;
        }
    }
    
    int start_pos = -1;
    // ----- START POSITION LOOKUP -----
    { // Primary start lookup
        auto range = static_refs.equal_range(ref_pos.primary_start.ref_id);
        const seq_element* primary_candidate = nullptr;
        for (auto it = range.first; it != range.second; ++it) {
            // Use candidate if it's present in the ordered list.
            if (std::find(ordered_elements.begin(), ordered_elements.end(), it->second) != ordered_elements.end()) {
                primary_candidate = it->second;
                break;
            }
        }
        if (primary_candidate != nullptr) {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Found " << ref_pos.primary_start.ref_id 
                              << " as the primary start reference candidate.\n"
                              << "Primary start details: is_start=" << ref_pos.primary_start.is_start 
                              << ", offset=" << ref_pos.primary_start.offset 
                              << ", static pos=(" << primary_candidate->position.first << ","
                              << primary_candidate->position.second << ")" << std::endl;
                    RAD_COUT << oss.str();
                }
            }
            auto ref_order_pos = std::find_if(ordered_elements.begin(), ordered_elements.end(),
                                              [&](const seq_element* elem) {
                                                  return elem->class_id == ref_pos.primary_start.ref_id;
                                              });
            bool is_ordered = true;
            if (ref_order_pos != ordered_elements.begin()) {
                auto prev = std::prev(ref_order_pos);
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Primary start: Previous element in ordered list is " 
                                  << (*prev)->class_id << " (pos: " 
                                  << (*prev)->position.first << "-" << (*prev)->position.second 
                                  << ", global_class: " << (*prev)->global_class << ")." << std::endl;
                        RAD_COUT << oss.str();
                    }
                }
                if ((*prev)->global_class != "poly_tail" &&
                    (*prev)->position.second > primary_candidate->position.first) {
                    is_ordered = false;
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "Primary start out-of-order: previous element's end (" 
                                      << (*prev)->position.second 
                                      << ") > candidate's beginning (" 
                                      << primary_candidate->position.first << ")." << std::endl;
                            RAD_COUT << oss.str();
                        }
                    }
                }
            } else {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Primary start candidate is the first element in the ordered list." << std::endl;
                        RAD_COUT << oss.str();
                    }
                }
            }
            if (is_ordered) {
                start_pos = ref_pos.primary_start.is_start ?
                    primary_candidate->position.first + ref_pos.primary_start.offset :
                    primary_candidate->position.second + ref_pos.primary_start.offset;
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Using primary start mapping: calculated start_pos = " << start_pos << std::endl;
                        RAD_COUT << oss.str();
                    }
                }
            }
        } else {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Primary start reference " << ref_pos.primary_start.ref_id 
                              << " not found in static_refs." << std::endl;
                    RAD_COUT << oss.str();
                }
            }
        }
    }
    
    // If primary failed, try secondary start.
    if (start_pos <= 0) {
        auto range = static_refs.equal_range(ref_pos.secondary_start.ref_id);
        const seq_element* secondary_candidate = nullptr;
        for (auto it = range.first; it != range.second; ++it) {
            if (std::find(ordered_elements.begin(), ordered_elements.end(), it->second) != ordered_elements.end()) {
                secondary_candidate = it->second;
                break;
            }
        }
        auto ref_order_pos = std::find_if(ordered_elements.begin(), ordered_elements.end(),
        [&](const seq_element* elem) {
            return elem->class_id == ref_pos.secondary_start.ref_id;
        });
        if (secondary_candidate != nullptr) {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Found " << ref_pos.secondary_start.ref_id 
                              << " as the secondary start reference candidate."
                              << "Secondary start details: is_start=" << ref_pos.secondary_start.is_start 
                              << ", offset=" << ref_pos.secondary_start.offset 
                              << ", static pos=(" << secondary_candidate->position.first << ","
                              << secondary_candidate->position.second << ")" << std::endl;
                    RAD_COUT << oss.str();
                }
            }
            bool is_ordered = true;
            if (ref_order_pos != ordered_elements.begin()) {
                auto prev = std::prev(ref_order_pos);
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Secondary start: Previous element in ordered list is " 
                                  << (*prev)->class_id << " (pos: " 
                                  << (*prev)->position.first << "-" << (*prev)->position.second 
                                  << ")." << std::endl;
                        RAD_COUT << oss.str();
                    }
                }
                if ((*prev)->position.second > secondary_candidate->position.first) {
                    is_ordered = false;
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "Secondary start out-of-order: previous element's end (" 
                                      << (*prev)->position.second 
                                      << ") > candidate's beginning (" 
                                      << secondary_candidate->position.first << ")." << std::endl;
                            RAD_COUT << oss.str();
                        }
                    }
                }
            } else {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Secondary start candidate is the first element in the ordered list." << std::endl;
                        RAD_COUT << oss.str();
                    }
                }
            }
            if (is_ordered) {
                start_pos = ref_pos.secondary_start.is_start ?
                    secondary_candidate->position.first + ref_pos.secondary_start.offset :
                    secondary_candidate->position.second + ref_pos.secondary_start.offset;
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Using secondary start mapping: calculated start_pos = " << start_pos << std::endl;
                        RAD_COUT << oss.str();
                    }
                }
            } 
        } else {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Secondary start reference " << ref_pos.secondary_start.ref_id 
                              << " not found in static_refs." << std::endl;
                    RAD_COUT << oss.str();
                }
            }
        if (!ref_pos.secondary_start.add_flags.empty() && ref_pos.secondary_start.add_flags == "left_truncated") {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Secondary start flagged as left_truncated; attempting fallback." << std::endl;
                    RAD_COUT << oss.str();
                }
            }
                // Look for the nearest static element to the right.
            auto nearest_static = std::find_if(ref_order_pos, ordered_elements.end(),
                                                [](const seq_element* elem) { return elem->type == "static"; });
                if (nearest_static != ordered_elements.end()) {
                    auto fallback_range = static_refs.equal_range((*nearest_static)->class_id);
                    const seq_element* fallback_candidate = nullptr;
                    for (auto it = fallback_range.first; it != fallback_range.second; ++it) {
                        if (it->second == *nearest_static) {
                            fallback_candidate = it->second;
                            break;
                        }
                    }
                    if (fallback_candidate != nullptr) {
                        start_pos = fallback_candidate->position.first + ref_pos.secondary_start.offset;
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "Using left terminal linked fallback: " 
                                            << fallback_candidate->class_id 
                                            << " with start_pos = " << start_pos << std::endl;
                                RAD_COUT << oss.str();
                            }
                        }
                    } else {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "Fallback static reference for secondary start not found." << std::endl;
                                RAD_COUT << oss.str();
                            }
                        }
                    }
                } else {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "No suitable fallback static element found for secondary start." << std::endl;
                            RAD_COUT << oss.str();
                        }
                    }
                }
            }
        }
    }
    
    // --concat-hmm: at a fold-back cut with no start reference, the read element starts at the child boundary.
    if (start_pos <= 0 && (hmm_fold_edges & 1))
        start_pos = hmm_fold_boundary(ref_pos.primary_start.ref_id, ref_pos.secondary_start.ref_id, ordered_elements, *var_it, true);

    // ----- STOP POSITION LOOKUP -----
    int stop_pos = -1;
    { // Primary stop lookup
        auto range = static_refs.equal_range(ref_pos.primary_stop.ref_id);
        const seq_element* primary_stop_candidate = nullptr;
        for (auto it = range.first; it != range.second; ++it) {
            if (std::find(ordered_elements.begin(), ordered_elements.end(), it->second) != ordered_elements.end()) {
                primary_stop_candidate = it->second;
                break;
            }
        }
        if (primary_stop_candidate != nullptr) {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Found " << ref_pos.primary_stop.ref_id 
                              << " as the primary stop reference candidate.\n"
                              << "Primary stop details: is_start=" << ref_pos.primary_stop.is_start 
                              << ", offset=" << ref_pos.primary_stop.offset 
                              << ", static pos=(" << primary_stop_candidate->position.first << ","
                              << primary_stop_candidate->position.second << ")" << std::endl;
                    RAD_COUT << oss.str();
                }
            }
            auto ref_order_pos = std::find_if(ordered_elements.begin(), ordered_elements.end(),
                                              [&](const seq_element* elem) {
                                                  return elem->class_id == ref_pos.primary_stop.ref_id;
                                              });
            bool is_ordered = true;
            if (ref_order_pos != ordered_elements.begin()) {
                auto prev = std::prev(ref_order_pos);
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Primary stop: Previous element in ordered list is " 
                                  << (*prev)->class_id << " (pos: " 
                                  << (*prev)->position.first << "-" << (*prev)->position.second 
                                  << ")." << std::endl;
                        RAD_COUT << oss.str();
                    }
                }
                if ((*prev)->position.second > primary_stop_candidate->position.first) {
                    is_ordered = false;
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "Primary stop out-of-order: previous element's end (" 
                                      << (*prev)->position.second 
                                      << ") > candidate's beginning (" 
                                      << primary_stop_candidate->position.first << ")." << std::endl;
                            RAD_COUT << oss.str();
                        }
                    }
                }
            } else {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Primary stop candidate is the first element in the ordered list." << std::endl;
                        RAD_COUT << oss.str();
                    }
                }
            }
            if (is_ordered) {
                stop_pos = ref_pos.primary_stop.is_start ?
                    primary_stop_candidate->position.first + ref_pos.primary_stop.offset :
                    primary_stop_candidate->position.second + ref_pos.primary_stop.offset;
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Using primary stop mapping: calculated stop_pos = " << stop_pos << std::endl;
                        RAD_COUT << oss.str();
                    }
                }
            }
        } else {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Primary stop reference " << ref_pos.primary_stop.ref_id 
                              << " not found in static_refs." << std::endl;
                    RAD_COUT << oss.str();
                }
            }
        }
    }
    
    // If primary stop mapping failed, try secondary.
    if (stop_pos <= 0) {
        if(verbose){
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "Primary stop mapping failed; attempting secondary stop mapping.\n"
                          << "Secondary stop reference ID: " << ref_pos.secondary_stop.ref_id
                          << "\nSecondary flags: " << ref_pos.secondary_stop.add_flags << std::endl;
                RAD_COUT << oss.str();
            }
        }
        auto range = static_refs.equal_range(ref_pos.secondary_stop.ref_id);
        const seq_element* secondary_stop_candidate = nullptr;
        for (auto it = range.first; it != range.second; ++it) {
            if (std::find(ordered_elements.begin(), ordered_elements.end(), it->second) != ordered_elements.end()) {
                secondary_stop_candidate = it->second;
                break;
            }
        }
        //moved this out of the verbosity loop
        auto ref_order_pos = std::find_if(
            ordered_elements.begin(), ordered_elements.end(),[&](const seq_element* elem) {
               return elem->class_id == ref_pos.secondary_stop.ref_id;
                }
        );
        if (secondary_stop_candidate != nullptr) {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Found " << ref_pos.secondary_stop.ref_id 
                              << " as the secondary stop reference candidate."
                              << "Secondary stop details: is_start=" << ref_pos.secondary_stop.is_start 
                              << ", offset=" << ref_pos.secondary_stop.offset 
                              << ", static pos=(" << secondary_stop_candidate->position.first << ","
                              << secondary_stop_candidate->position.second << ")" << std::endl;
                    RAD_COUT << oss.str();
                }
            }
            bool is_ordered = true;
            if (ref_order_pos != ordered_elements.begin()) {
                auto prev = std::prev(ref_order_pos);
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Secondary stop: Previous element in ordered list is " 
                                << (*prev)->class_id << " (pos: " 
                                << (*prev)->position.first << "-" << (*prev)->position.second 
                                << ")." << std::endl;
                        RAD_COUT << oss.str();
                    }
                }
                if ((*prev)->position.second > secondary_stop_candidate->position.first)
                    is_ordered = false;
            } else {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Secondary stop candidate is the first element in the ordered list." << std::endl;
                        RAD_COUT << oss.str();
                    }
                }
            }
            if (is_ordered) {
                stop_pos = ref_pos.secondary_stop.is_start ?
                    secondary_stop_candidate->position.first + ref_pos.secondary_stop.offset :
                    secondary_stop_candidate->position.second + ref_pos.secondary_stop.offset;
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Using secondary stop mapping: calculated stop_pos = " << stop_pos << std::endl;
                        RAD_COUT << oss.str();
                    }
                }
            }
        } else {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Secondary stop reference " << ref_pos.secondary_stop.ref_id 
                              << " not found in static_refs." << std::endl;
                    RAD_COUT << oss.str();
                }
            }
            if (!ref_pos.secondary_stop.add_flags.empty() && ref_pos.secondary_stop.add_flags == "right_truncated") {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Secondary stop flagged as right_truncated; attempting fallback." << std::endl;
                        RAD_COUT << oss.str();
                    }
                }
                // Look for the nearest static element to the left.
                auto nearest_static = std::find_if(std::make_reverse_iterator(ref_order_pos),
                                                    ordered_elements.rend(),
                                                    [](const seq_element* elem) { return elem->type == "static"; });
                if (nearest_static != ordered_elements.rend()) {
                    auto fallback_range = static_refs.equal_range((*nearest_static)->class_id);
                    const seq_element* fallback_candidate = nullptr;
                    for (auto it = fallback_range.first; it != fallback_range.second; ++it) {
                        if (it->second == *nearest_static) {
                            fallback_candidate = it->second;
                            break;
                        }
                    }
                    if (fallback_candidate != nullptr) {
                        stop_pos = fallback_candidate->position.second + ref_pos.secondary_stop.offset;
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "Using right terminal linked fallback: " 
                                            << fallback_candidate->class_id 
                                            << " with stop_pos = " << stop_pos << std::endl;
                                RAD_COUT << oss.str();
                            }
                        }
                    } else {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "Fallback static reference for secondary stop not found." << std::endl;
                                RAD_COUT << oss.str();
                            }
                        }
                    }
                } else {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "No suitable fallback static element found for secondary stop." << std::endl;
                            RAD_COUT << oss.str();
                        }
                    }
                }
            }
        }
    }    
    
    // --concat-hmm: mirror for the child end
    if (stop_pos <= 0 && (hmm_fold_edges & 2))
        stop_pos = hmm_fold_boundary(ref_pos.primary_stop.ref_id, ref_pos.secondary_stop.ref_id, ordered_elements, *var_it, false);

    var_positions = {start_pos, stop_pos};
    
    if (verbose) {
        #pragma omp critical
        {
            std::ostringstream oss;
            oss << "Final calculated positions: " << start_pos << ":" << stop_pos 
                      << " for variable element " << var_it->class_id << std::endl;
            RAD_COUT << oss.str();
        }
    }
    
    // Validate positions relative to read length.
    if (validate_var_positions(var_positions, read_seq.length())) {
        if (verbose) {
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "Positions validated for variable element " << var_it->class_id << std::endl;
                RAD_COUT << oss.str();
            }
        }
        std::string var_seq = read_seq.substr(var_positions.first - 1, 
                                              var_positions.second - var_positions.first + 1);
        id_index.modify(var_it, [&](seq_element& elem) {
            elem.position = var_positions;
            elem.seq = var_seq;
            elem.element_pass = true;
        });
        return true;
    }
    
    id_index.modify(var_it, [](seq_element& elem) {
        elem.element_pass = false;
    });
    if (verbose) {
        #pragma omp critical
        {
            std::ostringstream oss;
            oss << "Position validation failed for variable element " << var_it->class_id << std::endl;
            RAD_COUT << oss.str();
        }
    }
    return false;
}

/**
 * @brief Check if read sequence length meets the summative expected length of adapters
 * @param read_seq `std::string` sequencing read
 * @param total_expected_length total expected length of adapters, `int`
 * @return true if read length is at least 100 and meets or exceeds total expected length, false otherwise
 */
    bool read_adapter_sum_comp(const std::string& read_seq, int total_expected_length) const {
        return (read_seq.length() >= 100 && read_seq.length() >= static_cast<size_t>(total_expected_length));
    }
    
/**
 * @brief Calculate the total expected length of static adapters from the ReadLayout
 * @param layout `ReadLayout` object containing layout elements
 * @return total expected length of static adapters, `int`
 */
    // calculate the total expected length of static adapters
    int calc_total_static_len(const ReadLayout& layout) const {
        // Combined expected length of one strand's structural elements
        // (primers + poly-tail + barcode + UMI) — i.e. everything but the cDNA.
        // Summed straight from the layout's forward-orientation elements so the
        // value is a fixed property of the layout, not of a given read. Filtering
        // on elem.direction (not an "rc_" prefix) is deliberate: the reverse
        // poly-tail is named "poly_a", so a prefix test would miss it and the
        // forward "poly_t" would still be counted. Earlier this iterated the
        // read's sig_elements, which held BOTH orientations and double-counted
        // the barcode/UMI/adapters, inflating the minimum read length.
        int total = 0;
        for (const auto& elem : layout.by_order()) {
            if (elem.direction != "forward") continue;
            if (elem.expected_length) {
                total += *elem.expected_length;
            }
        }
        return total;
    }

/**
 * @brief Filter out reads whose cDNA is shorter than --min-read-length
 * @param elems elements after mask_and_trim_elements (the read element holds the written cDNA)
 * @param layout `ReadLayout` object containing layout elements
 * @param min_read_length the limit; -1 = the summed expected length of the layout's elements
 * @return true if the cDNA length meets the limit
 */
    bool filter_short_reads(const std::vector<std::reference_wrapper<const seq_element>>& elems,
                            const ReadLayout& layout, int min_read_length = -1) const {
        size_t read_len = 0;
        // find the single "read" element
        for (auto& e_ref : elems) {
            auto const& e = e_ref.get();
            if (e.global_class == "read" && e.seq.has_value()) {
                read_len = e.seq->size();
                break;
            }
        }
        // if no read (cDNA) at all -> drop
        if (read_len == 0){
            return false;
        }
        // Minimum informative cDNA length. Default (min_read_length < 0) is the
        // combined structural length (calc_total_static_len); users can override
        // via --min-read-length, including 0 to keep every read with any cDNA.
        int threshold = (min_read_length >= 0) ? min_read_length : calc_total_static_len(layout);
        // if the read is long enough, say yes; if the read is too short, say no
        return read_len >= static_cast<size_t>(threshold);
    }

/**
 * @brief Check for presence of forward direction static elements
 * @param elems vector of references to `seq_element` objects
 * @param layout `ReadLayout` object containing layout elements
 * @return false if at least one forward static element is present, true otherwise
 */
    bool filter_forward_direction_statics(
        const std::vector<std::reference_wrapper<const seq_element>>& elems, const ReadLayout& layout) const {
        bool forward_fail = false;
        int static_elements = 0;
        for (auto& e_ref : elems) {
            auto const& e = e_ref.get();
            if (
                e.direction == "forward" && 
                e.type == "static" && 
                e.global_class != "start" && 
                e.global_class != "stop" && 
                e.global_class != "poly_tail"
                ) {
                    static_elements++;
                } 
        }
        if(static_elements == 0){
            forward_fail = true;
        }
        return forward_fail;
    }

/**
 * @brief Check for presence of reverse direction static elements
 * @param elems vector of references to `seq_element` objects
 * @param layout `ReadLayout` object containing layout elements
 * @return true if at least one reverse static element is present, false otherwise
 */
    bool filter_reverse_direction_statics(
        const std::vector<std::reference_wrapper<const seq_element>>& elems, const ReadLayout& layout) const {
        bool reverse_fail = false;
        int static_elements = 0;
        for (auto& e_ref : elems) {
            auto const& e = e_ref.get();
            if (
                e.direction == "reverse" && 
                e.type == "static" && 
                e.global_class != "start" && 
                e.global_class != "stop" && 
                e.global_class != "poly_tail"
                ) {
                    static_elements++;
                } 
        }
        if(static_elements == 0){
            reverse_fail = true;
        }
        return reverse_fail;
    }

/**
 * @brief Filter static elements based on specified direction
 * @param elems vector of references to `seq_element` objects
 * @param layout `ReadLayout` object containing layout elements
 * @param direction direction string ("forward", "reverse", or other)
 * @return false if static elements are present, true if needs to be filtered
 * 
 * @note built one direction, built the other, built a wrapper, hit tab a lot, here we are
 */
    bool filter_direction_statics(
        const std::vector<std::reference_wrapper<const seq_element>>& elems, const ReadLayout& layout, 
        const std::string& direction) const {
        if(direction == "forward"){
            return filter_forward_direction_statics(elems, layout);
        } else if(direction == "reverse"){
            return filter_reverse_direction_statics(elems, layout);
        } else {
            return true;
        }
    }

/**
 * @brief Validate positions of a sequence element
 * @param e reference to `seq_element` object
 * @param direction direction string
 * @param verbose print verbose output
 * @return true if positions are valid, false otherwise
 */
    bool validate_sig_positions(seq_element& e, const std::string& direction, bool verbose){
        if (e.position.first <= 0 || e.position.second <= 0 || e.position.second <= e.position.first) {
            auto& idx = sig_elements.get<sig_id_tag>();
            idx.modify(idx.find(e.class_id), [](seq_element& x){ 
                x.element_pass = false; 
            });
            if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "  ["<< direction << "] invalid positions on ""<< e.class_id << ""\n";
                        RAD_COUT << oss.str();
                    }
                }
                return false;
            }
        return true;
    }

/**
 * @brief Filter overlapping variable elements in a set of sequence elements
 * @param elems vector of references to `seq_element` objects
 * @param verbose print verbose output
 */
    void filter_overlaps(const std::vector<std::reference_wrapper<seq_element>>& elems, bool verbose) {
        for (size_t i = 0; i < elems.size(); ++i) {
            auto &e1 = elems[i].get();
            if (e1.global_class=="start"||e1.global_class=="stop"
             || e1.type!="variable") continue;
            for (size_t j = i+1; j < elems.size(); ++j) {
                auto &e2 = elems[j].get();
                if (e2.global_class=="start"||e2.global_class=="stop"
                 || e2.type!="variable") continue;
                if (e1.position.second >= e2.position.first + 3) {
                    edit_elem(e1.class_id, [](seq_element &x){ x.element_pass = false; });
                    edit_elem(e2.class_id, [](seq_element &x){ x.element_pass = false; });
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "  Overlap detected between elements: "
                                  << e1.class_id << " and " << e2.class_id
                                  << " at positions (" << e1.position.first << "-" 
                                  << e1.position.second << ") and ("
                                  << e2.position.first << "-" 
                                  << e2.position.second << ")." << std::endl;
                            RAD_COUT << oss.str();
                        }
                    }
                }
            }
        }
    }

/**
 * @brief Group sequence elements by their direction
 * @return map of direction strings to vectors of references to `seq_element` objects
 */
    auto group_directionally(){
        std::map<std::string, std::vector<std::reference_wrapper<const seq_element>>> direction_elements;
        for (auto& elem : sig_elements)
        direction_elements[elem.direction].push_back(std::ref(elem));
        return direction_elements;
    }

/**
 * @brief Update per-barcode whitelist counters after read classification
 * @param sig  the SigString whose elements carry resolved counter pointers
 * @param layout ReadLayout (unused — counters already resolved during correction)
 * @param verbose enable verbose/debug output
 * @note Called once per read inside sigalign_filter.  Uses the counter pointer
 *       stashed on each seq_element during apply_barcode_correction — zero
 *       hash lookups, zero string operations, just atomic increments.
 */
    void update_bc_counts(
        SigString &sig,
        const ReadLayout &layout,
        bool verbose
    ) {
        const std::string& type = sig.read_type;

        for (auto const &elem : sig.elements()) {
            if (elem.global_class != "barcode" || !elem.resolved_counter) {
                continue;
            }

            auto& cnt = *elem.resolved_counter;

            if (type == "filtered" || (type != elem.direction && type != "concatenate")) {
                cnt.increment(barcode_counts::filtered);
                continue;
            }

            // Passing read — increment total + raw/corrected
            cnt.increment(barcode_counts::total);
            cnt.increment(elem.resolved_corrected
                ? barcode_counts::corrected
                : barcode_counts::raw);

            if (type == "forward" && elem.direction == "forward") {
                cnt.increment(barcode_counts::forw);
            } else if (type == "reverse" && elem.direction == "reverse") {
                cnt.increment(barcode_counts::rev);
            } else if (type == "concatenate") {
                cnt.increment(elem.direction == "forward"
                    ? barcode_counts::forw_concat
                    : barcode_counts::rev_concat);
            }
        }
    }

/** 
 * @brief Process sequence elements in a given direction: filter by length, mask overlaps, trim reads, and validate positions
 * @param direction direction string
 * @param elements vector of references to `seq_element` objects
 * @param read `read_streaming::sequence` object containing the read sequence
 * @param layout `ReadLayout` object containing layout elements
 * @param filtered_because string to append filtering reasons
 * @param verbose print verbose output
 * @return true if processing is successful, false if read is filtered
 */   
    bool process_direction_basic(
        const std::string& direction,
        std::vector<std::reference_wrapper<const seq_element>>& elements,
        const read_streaming::sequence& read,
        const ReadLayout& layout,
        std::string& filtered_because,
        bool verbose,
        int min_read_length = -1
    ) {

        if(filter_direction_statics(elements, layout, direction)){
            filtered_because += direction + "_FILTERED_NO_STATIC_ELEMENTS";
            set_info(filtered_because);
            return false;
        }

        // Sort elements by order
        std::sort(elements.begin(), elements.end(), [](const auto& a, const auto& b) {
            return a.get().order < b.get().order;
        });

        // Mask overlapping elements and trim reads
        mask_and_trim_elements(elements, read, verbose);

        // After the trim, so the length test counts the written cDNA on every route.
        if (!filter_short_reads(elements, layout, min_read_length)) {
            filtered_because += direction + "_FILTERED_READ_LENGTH";
            set_info(filtered_because);
            return false;
        }

        // Validate positions and check for overlaps
        return validate_element_positions(elements, direction, filtered_because, verbose);
    }

/**
 *  @brief Mask non-read elements in the read sequence and trim read elements
 * @param elements vector of references to `seq_element` objects
 * @param read `read_streaming::sequence` object containing the read sequence
 * @param verbose print verbose output
 */    
    void mask_and_trim_elements(std::vector<std::reference_wrapper<const seq_element>>& elements, 
            const read_streaming::sequence& read, bool verbose
    ) {

        std::string masked_read = read.seq;
        std::string masked_qual = read.is_fastq ? read.qual : "";
        
        // Step 1: Mask non-read elements
        for (const auto& elem_ref : elements) {
            const auto& elem = elem_ref.get();
            if (elem.global_class == "read" || elem.position.first <= 0 || 
                elem.position.second <= elem.position.first) {
                continue;
            }
            
            mask_element_in_sequence(masked_read, masked_qual, elem);
        }
        
        // Step 2: Extract and clean read elements
        for (const auto& elem_ref : elements) {
            const auto& elem = elem_ref.get();
            if (elem.global_class != "read" || !elem.seq) continue;
            
            extract_and_clean_read_element(elem, masked_read, masked_qual, read, verbose);
            break; // Only one read element per direction
        }
    }

/**
 * @brief Mask a non-read element in the read sequence by replacing its positions with 'N' and quality scores with the highest value.
 * @param masked_read reference to the masked read sequence string
 * @param masked_qual reference to the masked quality string
 * @param elem reference to the `seq_element` object to be masked
 */  
     void mask_element_in_sequence(std::string& masked_read, std::string& masked_qual, const seq_element& elem) {
        size_t start = elem.position.first - 1;
        size_t length = elem.position.second - elem.position.first + 1;
        
        if (start + length <= masked_read.size()) {
            std::fill(masked_read.begin() + start, masked_read.begin() + start + length, 'N');
            if (!masked_qual.empty() && masked_qual.size() == masked_read.size()) {
                std::fill(masked_qual.begin() + start, masked_qual.begin() + start + length, '\x7F');
            }
        }
    }

/**
 * @brief Extract and clean a read element from the masked read sequence.
 * @param elem reference to the `seq_element` object to be extracted and cleaned
 * @param masked_read reference to the masked read sequence string
 * @param masked_qual reference to the masked quality string
 * @param read reference to the original `read_streaming::sequence` object
 * @param verbose print verbose output
 */
    void extract_and_clean_read_element(const seq_element& elem, const std::string& masked_read, const std::string& masked_qual,
        const read_streaming::sequence& read, bool verbose
    ) {
        
        // Element coordinates are 1-based and inclusive. Validate them before
        // converting to size_t so a boundary sentinel (0 or read_length + 1)
        // cannot underflow, and allow a valid element to touch either end of
        // the read.
        if (elem.position.first < 1 ||
            elem.position.second <= elem.position.first ||
            elem.position.second > static_cast<int>(masked_read.size())) {
            if (verbose) {
                log_verbose("Skipping read element due to out-of-bounds parameters");
            }
            return;
        }

        size_t start = static_cast<size_t>(elem.position.first - 1);
        size_t length = static_cast<size_t>(
            elem.position.second - elem.position.first + 1);
        
        // Extract and clean sequence
        std::string window = masked_read.substr(start, length);
        std::string window_qual = read.is_fastq ? masked_qual.substr(start, length) : "";
        
        auto cleaned_result = remove_masked_positions(window, window_qual, read.is_fastq);
        std::string cleaned_seq = cleaned_result.first;
        std::string cleaned_qual = cleaned_result.second;
        
        // Update element
        edit_elem(elem.class_id, [cleaned_seq, cleaned_qual, is_fastq = read.is_fastq](seq_element& el) {
            *el.seq = cleaned_seq;
            if (is_fastq) {
                el.qual = cleaned_qual;
            }
        });
        
        if (verbose) {
            log_verbose("Cleaned read element: " + elem.class_id + " -> " + cleaned_seq);
        }
    }

/**
 * @brief Remove masked positions ('N') from a sequence window and its corresponding quality scores.
 * @param window sequence window strings
 * @param window_qual quality scores string
 * @param is_fastq boolean indicating if the read is in fastq format
 * @return pair of cleaned sequence and quality strings
 */
    std::pair<std::string, std::string> remove_masked_positions(const std::string& window, 
        const std::string& window_qual, bool is_fastq
    ) {
        std::string cleaned_seq, cleaned_qual;
        cleaned_seq.reserve(window.size());
        if (is_fastq) cleaned_qual.reserve(window.size());
        
        for (size_t i = 0; i < window.size(); ++i) {
            if (window[i] != 'N') {
                cleaned_seq.push_back(window[i]);
                if (is_fastq) {
                    cleaned_qual.push_back(window_qual[i]);
                }
            }
        }
        
        return {cleaned_seq, cleaned_qual};
    }

/**
 * @brief Validate positions of sequence elements and check for overlaps
 * @param elements vector of references to `seq_element` objects
 * @param direction direction string
 * @param filtered_because string to append filtering reasons
 * @param verbose print verbose output
 * @return true if positions are valid and no overlaps detected, false otherwise
 */
    bool validate_element_positions(std::vector<std::reference_wrapper<const seq_element>>& elements,
        const std::string& direction, std::string& filtered_because,
        bool verbose
    ) {

        // Check for invalid positions
        for (const auto& elem_ref : elements) {
            const auto& elem = elem_ref.get();
            if (elem.global_class == "start" || elem.global_class == "stop" || 
                elem.global_class == "poly_tail") {
                continue;
            }
            
            if (elem.position.first <= 0 || elem.position.second <= 0 || 
                elem.position.second <= elem.position.first) {
                
                mark_element_failed(elem.class_id);
                filtered_because += ":FILTERED_ELEMENT_" + elem.class_id + "_INVALID_POSITIONS";
                set_info(filtered_because);
                return false;
            }
        }
        
        // Check for overlapping variable elements
        return check_variable_element_overlaps(elements, direction, filtered_because, verbose);
    }

/**
 * @brief Check for overlapping variable elements in a set of sequence elements
 * @param elements vector of references to `seq_element` objects
 * @param direction direction string
 * @param filtered_because string to append filtering reasons
 * @param verbose print verbose output
 * @return true if no overlaps detected, false otherwise
 */
    bool check_variable_element_overlaps(std::vector<std::reference_wrapper<const seq_element>>& elements,
        const std::string& direction, std::string& filtered_because,
        bool verbose
    ) {
        
        for (size_t i = 0; i < elements.size(); i++) {
            const auto& e1 = elements[i].get();
            if (e1.global_class == "start" || e1.global_class == "stop" || e1.type != "variable") {
                continue;
            }
            
            for (size_t j = i + 1; j < elements.size(); j++) {
                const auto& e2 = elements[j].get();
                if (e2.global_class == "start" || e2.global_class == "stop" || e2.type != "variable") {
                    continue;
                }
                
                if (e1.position.second >= e2.position.first + 3) {
                    mark_element_failed(e1.class_id);
                    mark_element_failed(e2.class_id);
                    
                    if (verbose) {
                        log_verbose("Overlap detected: " + e1.class_id + " and " + e2.class_id);
                    }
                    
                    filtered_because += ":FILTERED_VARIABLE_" + e1.class_id + "_OVERLAPPING_POSITIONS";
                    set_info(filtered_because);
                    return false;
                }
            }
        }
        return true;
    }

/**
 * @brief Count valid reads in a set of sequence elements
 * @param elements vector of references to `seq_element` objects
 * @return number of valid reads
 */
    int count_valid_reads(const std::vector<std::reference_wrapper<const seq_element>>& elements
    ) {
        int count = 0;
        for (const auto& elem_ref : elements) {
            const auto& elem = elem_ref.get();
            if (elem.type == "variable" && 
                elem.element_pass && 
                elem.global_class != "barcode" && 
                elem.global_class == "read") {
                count++;
            }
        }
        return count;
    }

/**
 * @brief Parsed joint_barcode flag: group id and positional search offset range
 *
 * Flag grammar (in the `flags` column of a read layout CSV):
 *   joint_barcode          -> group 0, offset [-3, +3]  (backward compatible)
 *   joint_barcode:G        -> group G, offset [-3, +3]
 *   joint_barcode:G_O      -> group G, offset [0, O]
 *
 * Barcodes that share the same group id are validated jointly as a
 * contiguous chain.  The offset range controls how many bases of
 * positional jitter to search around the expected anchor position.
 */
    struct joint_barcode_info {
        int group = 0;
        int offset_min = -3;
        int offset_max = 3;
    };

/**
 * @brief Parse the joint_barcode flag from a layout element's flags string
 * @param flags the raw flags string from the ReadElement
 * @return parsed `joint_barcode_info` if the element carries a joint_barcode
 *         flag, or `std::nullopt` otherwise
 */
    std::optional<joint_barcode_info> parse_joint_barcode_flag(const std::string& flags) const {
        if (flags.empty()) return std::nullopt;
        std::string lowered = flags;
        std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

        // Tokenize on the same delimiters as has_layout_flag_token
        std::vector<std::string> tokens;
        std::string cur;
        for (char c : lowered) {
            if (c == ',' || c == ';' || c == '|' || c == ' ' || c == '\t' || c == '\n' || c == '\r') {
                if (!cur.empty()) { tokens.push_back(cur); cur.clear(); }
            } else {
                cur.push_back(c);
            }
        }
        if (!cur.empty()) tokens.push_back(cur);

        for (const auto& tok : tokens) {
            if (tok.rfind("joint_barcode", 0) != 0) continue;
            joint_barcode_info info;
            if (tok == "joint_barcode") return info; // bare default
            if (tok.size() > 14 && tok[13] == ':') {
                std::string suffix = tok.substr(14);
                auto underscore = suffix.find('_');
                if (underscore == std::string::npos) {
                    // joint_barcode:G
                    info.group = std::atoi(suffix.c_str());
                } else {
                    // joint_barcode:G_O
                    info.group = std::atoi(suffix.substr(0, underscore).c_str());
                    info.offset_max = std::atoi(suffix.substr(underscore + 1).c_str());
                    info.offset_min = 0;
                }
            }
            return info;
        }
        return std::nullopt;
    }

    bool has_layout_flag_token(const std::string& flags, const std::string& token) const {
        if (flags.empty() || token.empty()) return false;
        std::string lowered_flags = flags;
        std::string lowered_token = token;
        std::transform(lowered_flags.begin(), lowered_flags.end(), lowered_flags.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        std::transform(lowered_token.begin(), lowered_token.end(), lowered_token.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

        std::string cur;
        for (char c : lowered_flags) {
            if (c == ',' || c == ';' || c == '|' || c == ' ' || c == '\t' || c == '\n' || c == '\r') {
                if (cur == lowered_token) return true;
                cur.clear();
            } else {
                cur.push_back(c);
            }
        }
        return cur == lowered_token;
    }

    const ReadElement* find_layout_element(const ReadLayout& layout, const std::string& class_id) const {
        const auto& idx = layout.by_id();
        auto it = idx.find(class_id);
        if (it == idx.end()) return nullptr;
        return &(*it);
    }

    std::vector<int> effective_barcode_lengths(const ReadLayout& layout, const seq_element& elem) const {
        std::vector<int> out;
        if (const ReadElement* layout_elem = find_layout_element(layout, elem.class_id)) {
            out = layout_elem->length_candidates;
            if (out.empty() && layout_elem->expected_length.has_value()) {
                out.push_back(*layout_elem->expected_length);
            }
        }
        if (out.empty() && elem.seq.has_value() && !elem.seq->empty()) {
            out.push_back(static_cast<int>(elem.seq->size()));
        }
        out.erase(std::remove_if(out.begin(), out.end(), [](int v) { return v <= 0; }), out.end());
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        return out;
    }

    bool exact_whitelist_match(const whitelist::wl_entry& wl, const std::string& seq) const {
        int64_seq bits;
        bits.sequence_to_bits(seq);
        if (!bits.is_valid()) return false;
        // Precomputed mutation aliases are not exact barcode identities.
        return wl.true_bcs.has_identity_mapping(bits) || wl.global_bcs.check_wl_for(bits);
    }

    bool joint_barcode_pair_allowed(
        const std::vector<std::reference_wrapper<const seq_element>>& group,
        const ReadLayout& layout,
        const std::vector<int64_seq>& identities
    ) const {
        if (group.size() != identities.size()) return false;
        if (group.size() != 2) return true;

        // Pair masks use forward molecule order, including on reverse reads.
        bool reverse = group[0].get().direction == "reverse";
        size_t first = reverse
            ? (group[0].get().position.second >= group[1].get().position.second ? 0 : 1)
            : (group[0].get().position.first <= group[1].get().position.first ? 0 : 1);
        size_t second = 1 - first;
        for (const auto& ref : group) {
            auto found = layout.wl_map.maps.find(seq_utils::remove_rc(ref.get().class_id));
            if (found == layout.wl_map.maps.end()) return false;
            const auto& mask = found->second.get().spat_wl;
            if (mask.has_value() &&
                (mask->empty() || !mask->check(identities[first], identities[second]))) {
                return false;
            }
        }
        return true;
    }

/**
 * @brief Process a group of N joint barcodes that share the same joint_barcode group id
 * @param group  vector of references to the `seq_element` objects in the group, in layout order
 * @param layout `ReadLayout` containing element definitions and whitelist maps
 * @param read   `read_streaming::sequence` containing the read sequence
 * @param jb_info parsed `joint_barcode_info` (group id and search offset range)
 * @param verbose enable verbose/debug output
 * @return true if every member of the group matched its whitelist exactly, false otherwise
 * @note Members are searched as a contiguous chain: BC2 must start immediately after BC1 ends,
 *       BC3 after BC2, etc.  All length combinations are tried at each position.  The search
 *       offset is applied only to the anchor (first member) position.
 */
    bool try_process_joint_barcode_group(
        const std::vector<std::reference_wrapper<const seq_element>>& group,
        const ReadLayout& layout,
        const read_streaming::sequence& read,
        const joint_barcode_info& jb_info,
        bool verbose,
        bool* ambiguous = nullptr,
        bool* rejected_pair = nullptr
    ) {
        if (ambiguous) *ambiguous = false;
        if (rejected_pair) *rejected_pair = false;
        if (group.size() < 2) return false;

        // Validate: all must be barcodes, same direction
        const std::string& direction = group[0].get().direction;
        for (const auto& ref : group) {
            if (ref.get().global_class != "barcode") return false;
            if (ref.get().direction != direction) return false;
        }

        // Load whitelists and length candidates for each member
        struct member_info {
            const seq_element* elem;
            const whitelist::wl_entry* wl;
            std::vector<int> lengths;
            std::pair<int, int> position;
        };
        std::vector<member_info> members;
        members.reserve(group.size());

        for (const auto& ref : group) {
            const auto& elem = ref.get();
            auto wl_it = layout.wl_map.maps.find(seq_utils::remove_rc(elem.class_id));
            if (wl_it == layout.wl_map.maps.end()) return false;
            auto lens = effective_barcode_lengths(layout, elem);
            if (lens.empty()) return false;
            members.push_back({&elem, &wl_it->second.get(), std::move(lens), elem.position});
        }

        const int read_len = static_cast<int>(read.seq.size());
        if (read_len <= 0) return false;

        // Orient read and positions
        std::string oriented_read = read.seq;
        const bool reverse = (direction == "reverse");
        if (reverse) {
            oriented_read = seq_utils::revcomp(oriented_read);
            for (auto& m : members) {
                m.position = std::make_pair(
                    read_len - m.position.second + 1,
                    read_len - m.position.first + 1);
            }
        }

        // Sort members by position (leftmost first in oriented read)
        std::vector<size_t> order(members.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            return members[a].position.first < members[b].position.first;
        });

        // A fixed-length UMI directly before the chain shares its molecular
        // boundary with BC1. Primer-based extraction may have been shifted; the
        // resolved barcode placement supplies the more specific boundary.
        std::string adjacent_umi_id;
        int adjacent_umi_length = 0;
        const auto* first_layout = find_layout_element(layout, members[order[0]].elem->class_id);
        if (first_layout) {
            const auto& ordered_layout = layout.by_order();
            auto previous = ordered_layout.find(first_layout->order + (reverse ? 1 : -1));
            if (previous != ordered_layout.end() && previous->direction == direction &&
                previous->global_class == "umi" && previous->type == "variable" &&
                previous->expected_length.value_or(0) > 0 && previous->length_candidates.size() <= 1) {
                adjacent_umi_id = previous->class_id;
                adjacent_umi_length = *previous->expected_length;
            }
        }

        // Build offset search order centered on 0
        const int off_min = jb_info.offset_min;
        const int off_max = jb_info.offset_max;
        std::vector<int> offset_order;
        offset_order.reserve(static_cast<size_t>(off_max - off_min + 1));
        for (int d = 0; d <= std::max(std::abs(off_min), std::abs(off_max)); ++d) {
            if (d >= off_min && d <= off_max) offset_order.push_back(d);
            if (-d >= off_min && -d <= off_max && d != 0) offset_order.push_back(-d);
        }

        const int anchor_start = members[order[0]].position.first;
        const int oriented_len = static_cast<int>(oriented_read.size());

        // Per-member hit result
        struct member_hit {
            std::string seq;
            int start = -1;
            int end = -1;
        };

        // Recursive chain search: try all length combinations for members in order
        // cursor = current position in the oriented read (1-based)
        // depth  = which member in `order` we're trying to match
        std::vector<member_hit> hits(members.size());
        std::vector<member_hit> accepted_hits;
        std::string accepted_umi;
        bool conflicting_identity = false;
        bool mask_rejected_candidate = false;

        std::function<void(int, size_t)> search_chain = [&](int cursor, size_t depth) {
            if (conflicting_identity) return;
            if (depth == order.size()) {
                std::vector<int64_seq> identities;
                identities.reserve(hits.size());
                for (const auto& hit : hits) identities.emplace_back(hit.seq);
                if (!joint_barcode_pair_allowed(group, layout, identities)) {
                    mask_rejected_candidate = true;
                    return;
                }
                std::string candidate_umi;
                if (adjacent_umi_length > 0) {
                    const int start = hits[order[0]].start - adjacent_umi_length;
                    if (start < 1) return;
                    candidate_umi = oriented_read.substr(static_cast<size_t>(start - 1),
                                                         static_cast<size_t>(adjacent_umi_length));
                }
                if (accepted_hits.empty()) {
                    accepted_hits = hits;
                    accepted_umi = std::move(candidate_umi);
                } else {
                    if (candidate_umi != accepted_umi) conflicting_identity = true;
                    for (size_t i = 0; i < hits.size(); ++i) {
                        if (hits[i].seq != accepted_hits[i].seq) {
                            conflicting_identity = true;
                            break;
                        }
                    }
                }
                return;
            }
            size_t mi = order[depth];
            const auto& m = members[mi];
            for (int len : m.lengths) {
                int start = cursor;
                int end = start + len - 1;
                if (end > oriented_len) continue;
                std::string seq = oriented_read.substr(
                    static_cast<size_t>(start - 1), static_cast<size_t>(len));
                if (!exact_whitelist_match(*m.wl, seq)) continue;
                hits[mi] = {seq, start, end};
                search_chain(end + 1, depth + 1);
                if (conflicting_identity) return;
            }
        };

        for (int off : offset_order) {
            int start = anchor_start + off;
            if (start < 1) continue;
            search_chain(start, 0);
            if (conflicting_identity) break;
        }
        if (conflicting_identity) {
            if (ambiguous) *ambiguous = true;
            if (verbose) log_verbose("JOINT_BARCODE_AMBIGUOUS: multiple barcode or adjacent UMI identities");
            return false;
        }
        if (accepted_hits.empty()) {
            if (rejected_pair) *rejected_pair = mask_rejected_candidate;
            return false;
        }
        hits = std::move(accepted_hits);

        // Convert hits back to original coordinates and apply corrections
        auto to_original = [read_len, reverse](std::pair<int, int> p) -> std::pair<int, int> {
            if (!reverse) return p;
            return std::make_pair(read_len - p.second + 1, read_len - p.first + 1);
        };

        auto& id_index = sig_elements.get<sig_id_tag>();
        if (adjacent_umi_length > 0) {
            const int stop = hits[order[0]].start - 1;
            const auto position = to_original({stop - adjacent_umi_length + 1, stop});
            auto found = id_index.find(adjacent_umi_id);
            if (found == id_index.end() ||
                (read.is_fastq && position.second > static_cast<int>(read.qual.size()))) return false;
            id_index.modify(found, [&](seq_element& elem) {
                elem.original_seq = elem.seq;
                elem.position = position;
                elem.seq = read.seq.substr(static_cast<size_t>(position.first - 1),
                                           static_cast<size_t>(adjacent_umi_length));
                if (read.is_fastq) {
                    elem.qual = read.qual.substr(static_cast<size_t>(position.first - 1),
                                                 static_cast<size_t>(adjacent_umi_length));
                }
                elem.element_pass = true;
            });
        }
        for (size_t mi = 0; mi < members.size(); ++mi) {
            const auto& m = members[mi];
            const auto& h = hits[mi];
            auto final_pos = to_original({h.start, h.end});

            int64_seq bits;
            bits.sequence_to_bits(h.seq);
            if (!bits.is_valid()) return false;
            apply_barcode_correction(*m.elem, bits, layout, verbose);

            auto it = id_index.find(m.elem->class_id);
            if (it != id_index.end()) {
                id_index.modify(it, [&](seq_element& e) {
                    e.position = final_pos;
                });
            }
        }

        if (verbose) {
            std::ostringstream oss;
            oss << "JOINT_BARCODE_HIT group=" << jb_info.group;
            for (size_t mi = 0; mi < members.size(); ++mi) {
                size_t idx = order[mi];
                oss << " " << members[idx].elem->class_id << "=" << hits[idx].seq;
            }
            log_verbose(oss.str());
        }
        return true;
    }

/**
 * @brief process barcodes in a given direction: apply corrections and validate
 * @param direction direction string
 * @param elements vector of references to `seq_element` objects
 * @param read `read_streaming::sequence` object containing the read sequence
 * @param layout `ReadLayout` object containing layout elements
 * @param gen_mut general mutation rate for barcode correction
 * @param mode correction mode for barcode correction (offensive or defensive, which whitelist to check first)
 * @param verbose print verbose output
 * @return true if all barcodes pass, false if any barcode fails
 * @note If multiple barcodes are present, the direction passes if at least one barcode passes.
 */
    bool process_barcodes_for_direction(
        const std::string& direction, 
        std::vector<std::reference_wrapper<const seq_element>>& elements,
        const read_streaming::sequence& read, 
        const ReadLayout& layout, 
        int gen_mut,  
        const std::string& mode, 
        const std::string& joint_bc_mode,
        bool verbose
    ) {
        
        // Collect barcode elements and detect joint groups from layout flags
        struct bc_entry {
            std::string class_id;
            std::optional<joint_barcode_info> jb;
        };
        std::vector<bc_entry> barcode_entries;
        barcode_entries.reserve(elements.size());
        for (const auto& elem_ref : elements) {
            const auto& elem = elem_ref.get();
            if (elem.global_class != "barcode") continue;
            const ReadElement* le = find_layout_element(layout, elem.class_id);
            std::optional<joint_barcode_info> jb;
            if (le) jb = parse_joint_barcode_flag(le->flags);
            barcode_entries.push_back({elem.class_id, jb});
        }

        bool has_multiple_barcodes = barcode_entries.size() > 1;
        bool strict_joint_mode = (joint_bc_mode == "strict");
        bool all_barcodes_passed = true;
        bool any_barcode_passed = false;
        bool joint_identity_rejected = false;

        auto& id_index = sig_elements.get<sig_id_tag>();
        for (size_t i = 0; i < barcode_entries.size();) {
            const auto& entry = barcode_entries[i];

            // If this barcode has a joint group flag, collect the full group
            if (entry.jb.has_value()) {
                int group_id = entry.jb->group;
                std::vector<std::reference_wrapper<const seq_element>> group;
                joint_barcode_info group_info = *entry.jb;
                size_t j = i;
                while (j < barcode_entries.size() &&
                       barcode_entries[j].jb.has_value() &&
                       barcode_entries[j].jb->group == group_id) {
                    auto it = id_index.find(barcode_entries[j].class_id);
                    if (it == id_index.end()) break;
                    group.push_back(std::cref(*it));
                    ++j;
                }

                bool ambiguous = false;
                bool rejected_pair = false;
                bool exact_group = group.size() >= 2 &&
                    try_process_joint_barcode_group(group, layout, read, group_info, verbose,
                                                    &ambiguous, &rejected_pair);
                std::vector<bool> resolved(group.size(), exact_group);
                if (!exact_group && !ambiguous && !rejected_pair) {
                    // Default mode can retain independently resolved components.
                    // A failed component is not a corrected barcode and must not
                    // be written into an apparently complete molecular identity.
                    for (size_t member = 0; member < group.size(); ++member) {
                        seq_element elem = group[member].get();
                        resolved[member] = elem.seq &&
                            process_single_barcode(elem, layout, read, gen_mut, mode, verbose);
                    }
                }
                if (!ambiguous && !rejected_pair && group.size() >= 2 &&
                    std::all_of(resolved.begin(), resolved.end(), [](bool value) { return value; })) {
                    // Check the full resolved pair before output filtering: a
                    // suppressed component cannot hide a known mask violation.
                    std::vector<int64_seq> identities;
                    identities.reserve(group.size());
                    for (const auto& ref : group) identities.emplace_back(*ref.get().seq);
                    rejected_pair = !joint_barcode_pair_allowed(group, layout, identities);
                }
                const bool rejected = ambiguous || rejected_pair;
                size_t retained = 0;
                for (size_t member = 0; member < group.size(); ++member) {
                    const auto& elem = group[member].get();
                    if (rejected || !resolved[member] || !elem.write.value_or(true)) {
                        edit_elem(elem.class_id, [](seq_element& failed) {
                            failed.element_pass = false;
                            failed.write = false;
                        });
                    } else {
                        ++retained;
                    }
                }
                if (retained != group.size() || group.size() < 2) {
                    all_barcodes_passed = false;
                }
                joint_identity_rejected = joint_identity_rejected || rejected;
                any_barcode_passed = any_barcode_passed || retained > 0;
                i = std::max(i + 1, j);
                continue;
            }

            // Process as a single barcode
            auto it = id_index.find(entry.class_id);
            if (it == id_index.end()) {
                ++i;
                continue;
            }
            seq_element elem = *it;
            bool barcode_passed = process_single_barcode(elem, layout, read, gen_mut, mode, verbose, /*locate=*/true);
            if (!barcode_passed) {
                all_barcodes_passed = false;
                mark_element_failed(elem.class_id);
            } else {
                any_barcode_passed = true;
            }
            ++i;
        }

        if (joint_identity_rejected) return false;
        if (!has_multiple_barcodes) {
            return all_barcodes_passed;
        }
        if (strict_joint_mode) {
            return all_barcodes_passed;
        }
        return any_barcode_passed;
    }

/**
 * @brief Process a single barcode: apply correction and update element
 * @param elem reference to `seq_element` object representing the barcode
 * @param layout `ReadLayout` object containing layout elements
 * @param read `read_streaming::sequence` object containing the read sequence
 * @param gen_mut general mutation rate for barcode correction
 * @param mode correction mode for barcode correction (offensive or defensive, which whitelist to check
 * first)
 * @param verbose print verbose output
 * @return true if barcode correction is successful, false otherwise
 */
    bool process_single_barcode(const seq_element& elem, const ReadLayout& layout, const read_streaming::sequence& read,
        int gen_mut, const std::string& mode, bool verbose, bool locate = false
    ) {

        if (verbose) {
            log_verbose("Processing barcode: " + elem.seq.value());
        }

        barcode_correction::bc_position pos;
        auto correction_result = barcode_correction::correct_barcode(
            elem, layout, read, verbose, gen_mut, mode, locate ? &pos : nullptr
        );

        // Apply correction and update element
        if (!correction_result.has_value()) {
            return false;
        }

        apply_barcode_correction(elem, correction_result.value(), layout, verbose);
        if (locate) {
            move_barcode_position(elem, pos, layout, read, verbose);
        }
        return true;
    }

/**
 * @brief Move a single barcode to where the whitelist match found it, with the elements whose primary start and
 * stop both reference it (the UMI). The read element never moves: it was extracted and masked before correction.
 * @param elem the barcode element before correction (layout position)
 * @param pos where correct_barcode located the barcode
 */
    void move_barcode_position(const seq_element& elem, const barcode_correction::bc_position& pos,
        const ReadLayout& layout, const read_streaming::sequence& read, bool verbose
    ) {
        using bc_position = barcode_correction::bc_position;
        namespace st = barcode_position_stats;
        auto kept = [&](std::atomic<size_t>& counter, const std::string& why) {
            counter.fetch_add(1, std::memory_order_relaxed);
            if (verbose) log_verbose("Barcode position kept: " + elem.class_id + " (" + why + ")");
        };
        switch (pos.kind) {
            case bc_position::unknown: kept(st::unknown, "position unknown"); return;
            case bc_position::window: kept(st::window, "window match"); return;
            case bc_position::tie: kept(st::tie, "ambiguous: " + std::to_string(pos.n_end) + " equally good ends"); return;
            case bc_position::clipped: kept(st::clipped, "match clipped at the window edge"); return;
            case bc_position::located: break;
        }
        if (pos.shift == 0) {
            kept(st::in_place, "located in place");
            return;
        }

        const bool reverse = elem.direction == "reverse";
        // reverse strand: the UMI side is at lower read coordinates
        const int delta = reverse ? -pos.shift : pos.shift;
        const int read_len = static_cast<int>(read.seq.size());
        const std::pair<int, int> bc_pos{elem.position.first + delta, elem.position.second + delta};
        if (bc_pos.first < 1 || bc_pos.second > read_len) {
            kept(st::out_of_read, "new position outside the read");
            return;
        }

        auto& id_index = sig_elements.get<sig_id_tag>();
        auto bc_it = id_index.find(elem.class_id);
        if (bc_it == id_index.end()) return;

        struct chained_move { std::string id; std::pair<int, int> position; };
        std::vector<chained_move> chained;
        auto& type_index = layout.by_type();
        auto variable_range = type_index.equal_range("variable");
        for (auto lit = variable_range.first; lit != variable_range.second; ++lit) {
            if (!lit->ref_pos || lit->direction != elem.direction) continue;
            const auto& rp = *lit->ref_pos;
            if (rp.primary_start.ref_id != elem.class_id || rp.primary_stop.ref_id != elem.class_id) continue;
            auto it = id_index.find(lit->class_id);
            if (it == id_index.end() || it->type != "variable") continue;
            if (it->position.first < 1 || it->position.second < it->position.first) continue;  // not placed
            const int start = (rp.primary_start.is_start ? bc_pos.first : bc_pos.second) + rp.primary_start.offset;
            const int stop = (rp.primary_stop.is_start ? bc_pos.first : bc_pos.second) + rp.primary_stop.offset;
            if (start < 1 || stop > read_len || stop < start) {
                kept(st::out_of_read, "new position outside the read");
                return;
            }
            chained.push_back({lit->class_id, {start, stop}});
        }

        id_index.modify(bc_it, [&](seq_element& e) { e.position = bc_pos; });
        for (const auto& m : chained) {
            id_index.modify(id_index.find(m.id), [&](seq_element& e) {
                e.position = m.position;
                const size_t off = static_cast<size_t>(m.position.first - 1);
                const size_t len = static_cast<size_t>(m.position.second - m.position.first + 1);
                e.seq = read.seq.substr(off, len);
                if (e.qual.has_value() && read.is_fastq && read.qual.size() >= off + len) {
                    e.qual = read.qual.substr(off, len);
                }
            });
        }
        st::moved.fetch_add(1, std::memory_order_relaxed);
        (reverse ? st::moved_reverse : st::moved_forward).fetch_add(1, std::memory_order_relaxed);

        if (verbose) {
            log_verbose("Barcode position moved: " + elem.class_id + " " + std::to_string(elem.position.first) + "-" +
                        std::to_string(elem.position.second) + " -> " + std::to_string(bc_pos.first) + "-" +
                        std::to_string(bc_pos.second) + " (shift " + std::to_string(pos.shift) + ", " +
                        std::to_string(chained.size()) + " chained element(s) moved)");
        }
    }

/**
 * @brief Apply barcode correction to a sequence element and update its information
 * @param elem reference to `seq_element` object representing the barcode
 * @param corrected_barcode `int64_seq` object representing the corrected barcode
 * @param layout `ReadLayout` object containing layout elements
 * @param verbose print verbose output
 */
    void apply_barcode_correction(const seq_element& elem, const int64_seq& corrected_barcode,const ReadLayout& layout,
        bool verbose
    ) {

        auto wl_it = layout.wl_map.maps.find(seq_utils::remove_rc(elem.class_id));
        if (wl_it == layout.wl_map.maps.end()) return;
        auto& wl = wl_it->second.get();
        std::string final_bc = corrected_barcode.bits_to_sequence();

        // Single-lookup resolution: find the entry once in true, then global.
        // Stash the counter pointer on the element so update_bc_counts avoids re-lookup.
        counter* resolved = nullptr;
        bool found_in_true = false;
        bool found_in_global = false;

        resolved = wl.true_bcs.find_counter(corrected_barcode);
        if (resolved) {
            found_in_true = true;
        } else {
            resolved = wl.global_bcs.find_counter(corrected_barcode);
            if (resolved) found_in_global = true;
        }

        bool true_wl_empty = wl.true_bcs.empty();
        std::string final_wl = found_in_true ? "true" : "global";

        // Update the element — stash resolved counter and corrected flag
        auto& id_index = sig_elements.get<sig_id_tag>();
        id_index.modify(id_index.find(elem.class_id), [&](seq_element& e) {
            e.original_seq = e.seq;
            e.seq = final_bc;
            e.element_pass = true;
            e.write = found_in_true || (true_wl_empty && found_in_global);
            e.resolved_counter = resolved;
            // Check if corrected: different from original (and not just its RC)
            if (e.original_seq.has_value() && !e.original_seq->empty()) {
                const std::string& orig = e.original_seq.value();
                if (final_bc != orig) {
                    e.resolved_corrected = (final_bc != seq_utils::revcomp(orig));
                }
            }
        });

        if (verbose) {
            log_verbose("Barcode corrected: " + elem.class_id + " -> " + final_bc + " (" + final_wl + ")");
        }
    }

/**
 * @brief Determine preliminary read direction based on valid directions and pass counts
 * @param direction_valid map of direction strings to boolean indicating validity
 * @param pass_counts map of direction strings to integer counts of passing reads
 * @param verbose print verbose output
 * @return `std::pair<std::string, std::string>` containing preliminary read type and final read type
 */
    std::pair<std::string, std::string> determine_read_direction(
        const std::map<std::string, 
        bool>& direction_valid, 
        const std::map<std::string, 
        int>& pass_counts, 
        bool verbose
    ) {
        
        int forward_count = (direction_valid.count("forward") && direction_valid.at("forward")) 
                           ? pass_counts.at("forward") : 0;
        int reverse_count = (direction_valid.count("reverse") && direction_valid.at("reverse")) 
                           ? pass_counts.at("reverse") : 0;
        
        bool forward_valid = (forward_count > 0 && forward_count >= reverse_count);
        bool reverse_valid = (reverse_count > 0 && reverse_count >= forward_count);
        
        std::string preliminary_type;
        if (forward_valid && reverse_valid) {
            preliminary_type = "concatenate";
        } else if (forward_valid) {
            preliminary_type = "forward";
        } else if (reverse_valid) {
            preliminary_type = "reverse";
        } else {
            preliminary_type = "filtered";
        }
        
        if (verbose) {
            log_verbose("Preliminary read type: " + preliminary_type + 
                       " (forward: " + std::to_string(forward_count) + 
                       ", reverse: " + std::to_string(reverse_count) + ")");
        }
        
        return {preliminary_type, preliminary_type};
    }

/**
 * @brief Determine final read type based on valid directions and pass counts
 * @param direction_valid map of direction strings to boolean indicating validity
 * @param pass_counts map of direction strings to integer counts of passing reads
 * @param verbose print verbose output
 */
    void determine_final_read_type(
        const std::map<std::string, bool>& direction_valid,
        const std::map<std::string, int>& pass_counts,
        bool verbose
    ) {

        int forward_count = (direction_valid.count("forward") && direction_valid.at("forward"))
            ? pass_counts.at("forward") : 0;
        int reverse_count = (direction_valid.count("reverse") && direction_valid.at("reverse")) 
                           ? pass_counts.at("reverse") : 0;
        
        bool forward_valid = (forward_count > 0 && forward_count >= reverse_count);
        bool reverse_valid = (reverse_count > 0 && reverse_count >= forward_count);
        
        if (forward_valid && reverse_valid) {
            set_type("concatenate");
        } else if (forward_valid) {
            set_type("forward");
        } else if (reverse_valid) {
            set_type("reverse");
        } else {
            set_type("filtered");
        }
    }

/**
 * @brief Mark a sequence element as failed based on its class ID
 * @param class_id class ID string of the sequence element
 */
    inline void mark_element_failed(const std::string& class_id) {
        auto& id_index = sig_elements.get<sig_id_tag>();
        auto it = id_index.find(class_id);
        if (it != id_index.end()) {
            id_index.modify(it, [](seq_element& e) noexcept {
                e.element_pass = false;
                if (e.global_class == "barcode") {
                    e.write = false;
                }
            });
        }
    }
    
    void add_failed_barcode_to_filter(const seq_element& elem, const ReadLayout& layout, bool verbose
    ) {
        if (!elem.seq.has_value()) return;
        
        auto wl_it = layout.wl_map.maps.find(seq_utils::remove_rc(elem.class_id));
        if (wl_it == layout.wl_map.maps.end()) return;
        auto& wl = wl_it->second.get();
        int64_seq failed_bc, rc_failed_bc;
        failed_bc.sequence_to_bits(elem.seq.value());
        rc_failed_bc.sequence_to_bits(seq_utils::revcomp(elem.seq.value()));
        
        #pragma omp critical
        {
            // Check if neither the barcode nor its reverse complement are already keys
            if (!wl.filter_bcs.check_wl_for(failed_bc) && !wl.filter_bcs.check_wl_for(rc_failed_bc)) {
                if (verbose) {
                    log_verbose("Adding failed barcode to filter: " + failed_bc.bits_to_sequence());
                }
                
                // Double-check that it's still not present
                auto failed_range = wl.filter_bcs.equal_range(failed_bc);
                auto rc_failed_range = wl.filter_bcs.equal_range(rc_failed_bc);
                
                if ((failed_range.first == failed_range.second) && 
                    (rc_failed_range.first == rc_failed_range.second)) {
                    // Entry doesn't exist - insert it
                    // Insertion into filter_bcs remains disabled.
                }
            }
        }
    }
    
    int count_barcodes_in_direction(
        const std::vector<std::reference_wrapper<const seq_element>>& elements
    ) {
        int count = 0;
        for (const auto& elem_ref : elements) {
            if (elem_ref.get().global_class == "barcode") count++;
        }
        return count;
    }

    int count_all_valid_elements(
        const std::vector<std::reference_wrapper<const seq_element>>& elements
    ) {
        int count = 0;
        for (const auto& elem_ref : elements) {
            const auto& elem = elem_ref.get();
            if (elem.type == "variable" && elem.element_pass && 
                (elem.global_class == "read" || elem.global_class == "barcode")) {
                count++;
            }
        }
        return count;
    }
    
    void log_verbose(const std::string& message) const {
        #pragma omp critical
        {
            std::ostringstream oss;
            oss << message << std::endl;
            RAD_COUT << oss.str();
        }
    }

    void log_final_results(
        const std::map<std::string, bool>& direction_valid,
        const std::map<std::string, int>& pass_counts
    ) {
        #pragma omp critical
        {
            std::ostringstream oss;
            oss << "Final read type: " << read_type
                << ", forward elements: " << (pass_counts.count("forward") ? pass_counts.at("forward") : 0)
                << ", reverse elements: " << (pass_counts.count("reverse") ? pass_counts.at("reverse") : 0)
                << std::endl;
            RAD_COUT << oss.str();
        }
    }

public:
/**
 * @brief Master function for static sequence alignment and detection. 
 * @param read `read_streaming::sequence` object containing the read sequence
 * @param layout `ReadLayout` object containing layout elements
 * @param verbose print verbose output
 * 
 * @brief This function aligns static sequence elements (e.g., adapters, poly-tails)
 * within a sequencing read based on the provided layout. It handles special cases
 * for 'start' and 'stop' elements, performs poly-tail detection, and uses alignment
 * statistics to determine expected regions for other static elements. sigalign_static
 * works as a two-pass system: the first pass uses edlib to quickly find candidates within
 * expected regions, and the second pass uses the SSW aligner for the best alignment.
 * Edlib will find the majority of cases, while SSW will detect adapters that are either in 
 * misalignment regions or are truncated on the edges of reads.
 * Aligned elements in the read are masked to prevent re-alignment or overlapping.
 */
   void sigalign_static(const read_streaming::sequence &read, const ReadLayout& layout, bool verbose,
                        std::vector<seq_element>* additional_hits = nullptr,
                        const static_restriction* restrict_to = nullptr,
                        const static_piece_frame* piece = nullptr) {
    aligner_tools aligner;
    auto& type_index = layout.by_type();
    auto static_range = type_index.equal_range("static");

    if(verbose){
        #pragma omp critical
        {
            std::ostringstream oss;
            oss << "\n=== Starting static alignment for " << sequence_id << " ===" << std::endl;
            RAD_COUT << oss.str();
        }
    }
    int max_distance = -1;
    // Make a mutable copy of the read to mask aligned regions
    std::string mutable_seq = read.seq;
    std::vector<const ReadElement*> deferred_polys;
    std::vector<std::pair<const ReadElement*, static_check_window>> spacer_retries;
    size_t read_length = read.seq.length();
    // Clipped hits wait for a counterpart that may align later; no new alignment, no wider budget.
    std::vector<std::pair<const ReadElement*, static_alignments>> clipped_pending;
    // A piece edge at an HMM cut is not a physical read end.
    const bool at_read_start = !piece || piece->offset == 0;
    const bool at_read_end = !piece || !piece->parent_seq ||
        piece->offset + read_length >= piece->parent_seq->size();

    for (auto it = static_range.first; it != static_range.second; ++it) {
        // --concat-hmm: poly tails go last so a run cannot mask the edge of an adapter next to it.
        if (restrict_to && it->global_class != "start" && it->global_class != "stop") {
            if (it->direction == restrict_to->direction) {
                if (it->global_class == "poly_tail") deferred_polys.push_back(&*it);
                else align_static_in_windows(*it, read, mutable_seq, *restrict_to, aligner, verbose, &spacer_retries);
            }
            continue;
        }
        if(verbose){
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "Processing static element: " << it->class_id << std::endl;
                RAD_COUT << oss.str();
                if(it->aligned_positions){
                    std::ostringstream aligned;
                    aligned << "  Aligned positions: " 
                              << it->aligned_positions->start_stats.first << ", "
                              << it->aligned_positions->start_stats.second << std::endl;
                    RAD_COUT << aligned.str();
                } else {
                    std::ostringstream fail;
                    fail << "  No aligned positions available." << std::endl;
                    RAD_COUT << fail.str();
                }
                if(it->misaligned_positions){
                    std::ostringstream misaligned;
                    misaligned << "  Misaligned positions: " 
                              << it->misaligned_positions->start_stats.first << ", "
                              << it->misaligned_positions->start_stats.second << std::endl;
                    RAD_COUT << misaligned.str();
                } else {
                    std::ostringstream fail_again;
                    fail_again << "  No misaligned positions available." << std::endl;
                    RAD_COUT << fail_again.str();
                }
            }
        }
        // For 'start' or 'stop' types, add an element without alignment
        if (it->global_class == "start" || it->global_class == "stop") {
            add_element(seq_element(
                it->class_id,
                it->global_class,
                std::nullopt,
                // Start/stop are virtual boundaries outside the 1-based,
                // inclusive read coordinates. Position-map offsets such as
                // start+1 and stop-1 therefore resolve to the first and last
                // real bases instead of shifting terminal-derived elements.
                (it->global_class == "start")
                    ? std::make_pair(0, 0)
                    : std::make_pair(
                        static_cast<int>(read_length) + 1,
                        static_cast<int>(read_length) + 1),
                "static",
                it->order,
                it->direction,
                std::nullopt,
                std::nullopt,
                std::nullopt
            ));
            continue;
        }

        // Use the calibrated threshold when available.  Otherwise use the
        // same ceil(30%) fallback as layout preparation so an omitted map row
        // cannot silently tighten an adapter from seven edits to four.
        if (it->misalignment_threshold) {
            max_distance = std::get<0>(*it->misalignment_threshold);
        } else {
            max_distance =
                adapter_thresholds::fallback_max_edit_distance(it->seq.length());
        }

        // poly-tail solution--importantly, only will look for poly-tails if you tell it to
        if (it->global_class == "poly_tail") {
            auto result = aligner.find_poly_tails(it->seq, read.seq, 14);
            if (result.success) {
                for (const auto& positions : result.positions) {
                    add_element(seq_element(
                        it->class_id,
                        it->global_class,
                        result.edit_distance,
                        positions,
                        "static",
                        it->order,
                        it->direction,
                        true,
                        std::nullopt,
                        mutable_seq.substr(positions.first - 1,
                                             positions.second - positions.first + 1)
                    ));
                    // Mask the found region so that it is not re-aligned
                    std::fill(mutable_seq.begin() + positions.first - 1,
                              mutable_seq.begin() + positions.second, 'X');
                }
            }
            continue;
        }

        // For other static elements (non-poly_tail), use the alignment statistics
        // to define the expected region for the adapter
        int expected_start, expected_end;
        // Check if the element has aligned and misaligned positions
        if (it->aligned_positions && it->misaligned_positions) {
            const auto& [start_mean, start_var] = it->aligned_positions->start_stats;
            // Compute expected region boundaries based on the aligned stats
            expected_start = (start_mean < 50.0) ? 1 
            : static_cast<int>(std::max(0.0, ((start_mean - (start_var * 1.25)) / 100.0) * read_length));

            expected_end = (start_mean < 50.0)
                ? static_cast<int>(std::min(((start_mean + (start_var * 1.25)) / 100.0) * read_length, static_cast<double>(read_length)))
                : static_cast<int>(read_length);

        } else {
            // Default to using the entire read length if no stats are available
            expected_start = 1;
            expected_end = static_cast<int>(read_length);
        }
        if (verbose) {
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "Expected region for " << it->class_id << ": " 
                          << expected_start << " - " << expected_end << std::endl;
                RAD_COUT << oss.str();
            }
        }

        // Use the entire mutable sequence as target
        static_alignments candidates;
        auto result = aligner.align_static_elements(it->seq,
                                                    mutable_seq,
                                                    verbose,
                                                    max_distance,
                                                    it->masked_seq,
                                                    verbose,
                                                    expected_start,
                                                    expected_end,
                                                    additional_hits ? &candidates : nullptr,
                                                    1, -1, true,
                                                    at_read_start, at_read_end);
        // Forward elements are visited before their reverse counterparts.
        // Do not let an imperfect hit mask an exact opposite-strand adapter at
        // the same locus; otherwise the later read boundary falls back to a
        // barcode/UMI and can retain adapter sequence in the written cDNA.
        if (result.success && result.edit_distance > 0 && result.positions.size() == 1 &&
            exact_opposite_at_locus(*it, layout, read, result.positions.front(), max_distance)) {
            result.success = false;
            result.positions.clear();
        }
        if (candidates.edit_distance > 0) {
            candidates.positions.erase(std::remove_if(candidates.positions.begin(), candidates.positions.end(),
                [&](const auto& position) {
                    return exact_opposite_at_locus(*it, layout, read, position, max_distance);
                }), candidates.positions.end());
            if (candidates.positions.empty()) candidates.success = false;
        }
        if (!result.success && result.edit_distance >= 0 &&
            result.edit_distance <= max_distance && result.positions.size() == 1 &&
            ((result.query_clip_left > 0) != (result.query_clip_right > 0))) {
            clipped_pending.emplace_back(&*it, result);
        }
        // --concat-hmm piece: keep the refused hit as evidence only (not an element, not masked).
        if (piece && !result.success && result.edit_distance >= 0 &&
            result.edit_distance <= max_distance && result.positions.size() == 1) {
            hmm_clip_hits.push_back({it->class_id, result.positions.front(), false});
            if (verbose) log_verbose("CLIP_EVIDENCE " + it->class_id + " " +
                std::to_string(result.positions.front().first) + ":" +
                std::to_string(result.positions.front().second) + " local_edits=" +
                std::to_string(result.edit_distance) + " clip=" +
                std::to_string(result.query_clip_left) + "/" + std::to_string(result.query_clip_right));
        }
        if (additional_hits && candidates.success) {
            for (const auto& pos : candidates.positions) {
                if (pos.first < 1 || pos.second < pos.first ||
                    pos.second > static_cast<int>(read_length)) continue;
                // Alternative tracebacks at the primary locus are not extra
                // molecules. Keep legacy coordinates/metadata there.
                if (result.success && std::any_of(result.positions.begin(), result.positions.end(),
                    [&](const auto& primary_pos) {
                        return pos.first <= primary_pos.second && primary_pos.first <= pos.second;
                    })) continue;
                additional_hits->emplace_back(it->class_id, it->global_class,
                    candidates.edit_distance, pos, "static", it->order, it->direction,
                    true, std::nullopt, read.seq.substr(pos.first - 1, pos.second - pos.first + 1));
                if (verbose) log_verbose("CONCAT_CANDIDATE " + it->class_id + " " +
                    std::to_string(pos.first) + ":" + std::to_string(pos.second) +
                    " edits=" + std::to_string(candidates.edit_distance));
            }
        }
            if (result.success) {
                // Positions returned are relative to the full read
                for (const auto& pos : result.positions) {
                    int adj_start = pos.first;
                    int adj_end = pos.second;
                    if (adj_start < 1) {
                        adj_start = 1;
                    }
                    if (adj_end > static_cast<int>(read_length)) {
                        adj_end = read_length;
                    }
                    std::pair<int, int> adjusted_positions = {adj_start, adj_end};
                    std::string full_aligned_seq = mutable_seq.substr(adj_start - 1,
                                                                       adj_end - adj_start + 1);

                    // Extract N-masked regions if query contains N's
                    std::string n_extracted = aligner.extract_n_masked_regions(result, it->seq);

                    seq_element primary(
                        it->class_id,
                        it->global_class,
                        result.edit_distance,
                        adjusted_positions,
                        "static",
                        it->order,
                        it->direction,
                        true,
                        std::nullopt,
                        n_extracted.empty() ? full_aligned_seq : n_extracted,
                        std::nullopt,
                        n_extracted.empty() ? std::nullopt : std::optional<std::string>(full_aligned_seq)
                    );
                    primary.query_complete = result.query_complete;
                    add_element(std::move(primary));

                    // Mask the aligned region so that it is not re-aligned
                    std::fill(mutable_seq.begin() + adj_start - 1,
                              mutable_seq.begin() + adj_end, 'X');
                }
                // Primary region yielded a valid match; skip further processing for this element
                continue;
            }
        }
    // A clipped end may face a molecule junction: accept it only next to an exact full-length
    // reverse-complement counterpart (extraction evidence only, never a split anchor).
    for (const auto& pending : clipped_pending) {
        const auto& element = *pending.first;
        const auto& result = pending.second;
        if (!has_exact_molecular_boundary(result, element, layout, read,
                                          additional_hits, piece)) continue;
        seq_element primary(element.class_id, element.global_class,
            result.edit_distance, result.positions.front(), "static", element.order,
            element.direction, true, std::nullopt, result.seq);
        primary.query_complete = false;
        add_element(std::move(primary));
        if (piece)
            hmm_clip_hits.erase(std::remove_if(hmm_clip_hits.begin(), hmm_clip_hits.end(),
                [&](const hmm_clip_hit& h) { return h.id == element.class_id; }), hmm_clip_hits.end());
        if (verbose) log_verbose("MOLECULAR_END_CLIP " + element.class_id + " " +
            std::to_string(result.positions.front().first) + ":" +
            std::to_string(result.positions.front().second) + " local_edits=" +
            std::to_string(result.edit_distance));
    }
        for (const ReadElement* poly : deferred_polys)
            align_static_in_windows(*poly, read, mutable_seq, *restrict_to, aligner, verbose);
        for (const auto& rq : spacer_retries)
            retry_at_spacer_offset(*rq.first, rq.second, read, mutable_seq, *restrict_to, aligner, verbose);
    }

    static bool exact_opposite_at_locus(const ReadElement& element, const ReadLayout& layout,
                                        const read_streaming::sequence& read,
                                        const std::pair<int, int>& position, int max_distance) {
        if (max_distance <= 0 || element.seq.empty() ||
            element.seq.find_first_not_of("ACGT") != std::string::npos) return false;
        const std::string opposite = seq_utils::revcomp(element.seq);
        bool has_counterpart = false;
        for (const auto& other : layout.by_type()) {
            if (other.type == "static" && other.direction != element.direction &&
                other.global_class == element.global_class && other.seq == opposite) {
                has_counterpart = true;
                break;
            }
        }
        if (!has_counterpart) return false;
        const size_t first = static_cast<size_t>(std::max(0, position.first - 1 - max_distance));
        for (size_t start = read.seq.find(opposite, first); start != std::string::npos;
             start = read.seq.find(opposite, start + 1)) {
            const int left = static_cast<int>(start) + 1;
            if (left > position.first + max_distance) break;
            const int right = static_cast<int>(start + opposite.size());
            // Both boundaries must describe the same locus. An independent
            // forward repeat elsewhere in the read must still be retained.
            if (left <= position.second && right >= position.first &&
                std::abs(right - position.second) <= max_distance) return true;
        }
        return false;
    }

    bool has_exact_molecular_boundary(const static_alignments& partial,
        const ReadElement& element, const ReadLayout& layout,
        const read_streaming::sequence& read,
        const std::vector<seq_element>* additional_hits = nullptr,
        const static_piece_frame* piece = nullptr) const {
        if (partial.query_complete || partial.positions.size() != 1 ||
            ((partial.query_clip_left > 0) == (partial.query_clip_right > 0)) ||
            element.seq.find_first_not_of("ACGTacgt") != std::string::npos) return false;
        // Only a missing end that faces away from all payload can be supported.
        auto outward_edge = [&](const ReadElement& adapter) {
            bool payload_left = false, payload_right = false;
            for (const auto& item : layout.by_order()) {
                if (item.direction != adapter.direction || item.type != "variable") continue;
                payload_left = payload_left || item.order < adapter.order;
                payload_right = payload_right || item.order > adapter.order;
            }
            return payload_right && !payload_left ? -1 : payload_left && !payload_right ? 1 : 0;
        };
        const int missing_edge = partial.query_clip_left > 0 ? -1 : 1;
        if (outward_edge(element) != missing_edge) return false;
        const auto pos = partial.positions.front();
        if (pos.first < 1 || pos.second < pos.first || pos.second > static_cast<int>(read.seq.size()) ||
            read.seq.compare(pos.first - 1, pos.second - pos.first + 1, partial.seq) != 0) return false;
        const std::string counterpart = seq_utils::revcomp(element.seq);
        auto supports = [&](const seq_element& hit) {
            if (hit.type != "static" || !hit.element_pass.value_or(false) ||
                !hit.query_complete || hit.edit_distance.value_or(-1) != 0 ||
                hit.global_class != element.global_class ||
                !((element.direction == "forward" && hit.direction == "reverse") ||
                  (element.direction == "reverse" && hit.direction == "forward"))) return false;
            const auto definition = layout.by_id().find(hit.class_id);
            if (definition == layout.by_id().end() || definition->seq != counterpart ||
                outward_edge(*definition) != -missing_edge ||
                hit.position.first < 1 || hit.position.second > static_cast<int>(read.seq.size()) ||
                hit.position.second - hit.position.first + 1 != static_cast<int>(counterpart.size()) ||
                read.seq.compare(hit.position.first - 1, counterpart.size(), counterpart) != 0) return false;
            const int gap = partial.query_clip_left > 0
                ? pos.first - hit.position.second - 1
                : hit.position.first - pos.second - 1;
            return gap >= 0 && gap <= 3;
        };
        for (const auto& hit : sig_elements) if (supports(hit)) return true;
        if (additional_hits) for (const auto& hit : *additional_hits) if (supports(hit)) return true;
        // --concat-hmm piece: the counterpart may lie across the HMM cut, so test the parent read.
        if (piece && piece->parent_seq) {
            bool defined = false;
            for (const auto& item : layout.by_order()) {
                if (item.type != "static" || item.seq != counterpart ||
                    item.global_class != element.global_class ||
                    !((element.direction == "forward" && item.direction == "reverse") ||
                      (element.direction == "reverse" && item.direction == "forward")) ||
                    outward_edge(item) != -missing_edge) continue;
                defined = true;
                break;
            }
            const std::string& parent = *piece->parent_seq;
            const long len = static_cast<long>(counterpart.size());
            const long piece_lo = piece->offset;
            const long piece_hi = piece_lo + static_cast<long>(read.seq.size());
            for (int gap = 0; defined && gap <= 3; ++gap) {
                const long start = partial.query_clip_left > 0
                    ? piece_lo + (pos.first - 1) - gap - len
                    : piece_lo + pos.second + gap;
                if (start < 0 || start + len > static_cast<long>(parent.size())) continue;
                if (start >= piece_lo && start + len <= piece_hi) continue;  // already tested above
                if (parent.compare(static_cast<size_t>(start), counterpart.size(), counterpart) == 0) return true;
            }
        }
        return false;
    }

    // --concat-hmm: align one static element inside each HMM window for it, with the aligner its expected
    // status calls for. Unconfirmed barcode-side adapters are queued for retry_at_spacer_offset.
    void align_static_in_windows(const ReadElement& elem, const read_streaming::sequence& read,
                                 std::string& mutable_seq, const static_restriction& rs,
                                 aligner_tools& aligner, bool verbose,
                                 std::vector<std::pair<const ReadElement*, static_check_window>>* retry = nullptr) {
        using concat_hmm::Status;
        if (elem.global_class == "poly_tail") {
            bool named = false;
            for (const auto& w : rs.windows) named = named || w.class_id == elem.class_id;
            if (!named) return;
            auto result = aligner.find_poly_tails(elem.seq, read.seq, 14);
            if (!result.success) return;
            for (auto positions : result.positions) {
                while (positions.first <= positions.second && mutable_seq[positions.first - 1] == 'X') ++positions.first;
                while (positions.second >= positions.first && mutable_seq[positions.second - 1] == 'X') --positions.second;
                if (positions.second - positions.first + 1 < 12) continue;
                bool overlaps = false;
                for (const auto& w : rs.windows)
                    overlaps = overlaps || (w.class_id == elem.class_id && positions.first <= w.hi && w.lo <= positions.second);
                if (!overlaps) continue;
                add_element(seq_element(elem.class_id, elem.global_class, result.edit_distance, positions, "static",
                                        elem.order, elem.direction, true, std::nullopt,
                                        mutable_seq.substr(positions.first - 1, positions.second - positions.first + 1)));
                std::fill(mutable_seq.begin() + positions.first - 1, mutable_seq.begin() + positions.second, 'X');
            }
            return;
        }
        if (!rs.info || elem.seq.empty()) return;
        const auto info_it = rs.info->elems.find(elem.class_id);
        if (info_it == rs.info->elems.end()) return;
        const concat_elem_info& info = info_it->second;
        concat_hmm_counters* ctr = rs.info->ctr;
        const int read_length = static_cast<int>(read.seq.length());
        const int max_distance = elem.misalignment_threshold
            ? std::get<0>(*elem.misalignment_threshold)
            : adapter_thresholds::fallback_max_edit_distance(elem.seq.length());
        const int m = info.m;
        auto poly_support = [&](const static_alignments& hit) {
            if (info.poly_id.empty() || hit.positions.empty()) return false;
            const int as = hit.positions.front().first, ae = hit.positions.front().second;
            for (const auto& pw : rs.windows) {
                if (pw.class_id != info.poly_id || pw.status != Status::FULL) continue;
                const int gap = info.poly_after ? (pw.lo - 1 + rs.info->seen_pad) - ae
                                                : (as - 1) - (pw.hi - rs.info->seen_pad);
                if (gap >= info.poly_gap_lo - 4 && gap <= info.poly_gap_hi + 4) return true;
            }
            return false;
        };
        // At a cut (not a physical read end), the primer may reach at most 4 bp into the barcode block that a
        // FULL partner implies.
        auto clamp_to_partner = [&](static_check_window& wc) {
            if (info.partner_id.empty() || (info.partner_after ? rs.read_start : rs.read_end)) return;
            int best = -1, bd = INT_MAX;
            for (size_t i = 0; i < rs.windows.size(); ++i) {
                const auto& pw = rs.windows[i];
                if (pw.class_id != info.partner_id || pw.status != Status::FULL) continue;
                const int d = info.partner_after ? pw.lo - wc.hi : wc.lo - pw.hi;
                if (d > -m && std::abs(d) < bd) { bd = std::abs(d); best = static_cast<int>(i); }
            }
            if (best < 0) return;
            const auto& pw = rs.windows[best];
            if (info.partner_after) wc.hi = std::min(wc.hi, pw.lo + rs.info->seen_pad - 1 - info.partner_gap_lo + 4);
            else wc.lo = std::max(wc.lo, pw.hi - rs.info->seen_pad + 1 + info.partner_gap_lo - 4);
        };
        for (const auto& w0 : rs.windows) {
            if (w0.class_id != elem.class_id) continue;
            static_check_window w = w0;
            clamp_to_partner(w);
            const int st = static_cast<int>(w.status);
            if (ctr) ctr->win[st].fetch_add(1, std::memory_order_relaxed);
            if (w.hi < w.lo) {
                if (retry && w.status != Status::TRUNCATED_AT_READ_END && (!info.poly_id.empty() || !info.partner_id.empty()))
                    retry->push_back({&elem, w0});
                continue;
            }
            static_alignments hit;
            bool missing_rules = false;
            bool is_fragment = false;
            bool whole_degraded = false;
            if (w.status == Status::FULL) {
                hit = aligner.align_static_elements(elem.seq, mutable_seq, verbose, std::max(info.k_full, w.edits),
                                                    elem.masked_seq, verbose, w.lo, w.hi, nullptr, w.lo, w.hi, false);
                missing_rules = !hit.success;
            } else if (w.status == Status::PARTIAL) {
                const int rf = std::max(0, w.retained_from);
                const int rt = w.retained_to < 0 ? m : std::min(m, w.retained_to);
                whole_degraded = rf == 0 && rt == m;
                const bool faces_boundary = (rf == 0 || info.left_outward) && (rt == m || info.right_outward);
                if (rt - rf >= 10 && rt - rf < m && faces_boundary) {
                    const int a = std::max(0, rf - 2), b = std::min(m, rt + 2);
                    hit = aligner.align_static_elements(elem.seq.substr(a, b - a), mutable_seq, verbose,
                                                        std::max(1, (rt - rf) / 8), "", verbose, w.lo, w.hi,
                                                        nullptr, w.lo, w.hi, false);
                    if (hit.success) {
                        hit.query_complete = (a == 0 && b == m);
                        is_fragment = true;
                    }
                }
                if (!hit.success) {
                    const bool at3 = w.hi >= read_length - 6;
                    const bool at5 = !at3 && w.lo <= 7;
                    if ((at3 && rf == 0 && info.right_outward) || (at5 && rt == m && info.left_outward)) {
                        hit = aligner.align_read_end_fragment(elem.seq, mutable_seq, at3, rt - rf + 4);
                        is_fragment = hit.success;
                    }
                }
                missing_rules = !hit.success;
            } else if (w.status == Status::TRUNCATED_AT_READ_END) {
                const bool at3 = w.hi >= read_length - 6;
                const bool at5 = !at3 && w.lo <= 7;
                if ((at3 && info.right_outward) || (at5 && info.left_outward)) {
                    // only a retained length seen in sequence (edits >= 0) bounds the search
                    const int keep = (w.edits >= 0 && w.retained_from >= 0 && w.retained_to > w.retained_from)
                        ? w.retained_to - w.retained_from + 4 : -1;
                    hit = aligner.align_read_end_fragment(elem.seq, mutable_seq, at3, keep);
                    is_fragment = hit.success;
                }
            } else {
                missing_rules = true;
            }
            if (missing_rules) {
                const int k = std::min(max_distance, info.k_missing);
                const bool extend = max_distance > k && !info.poly_id.empty();
                hit = aligner.align_static_elements(elem.seq, mutable_seq, verbose, extend ? k + 1 : k,
                                                    elem.masked_seq, verbose, w.lo, w.hi, nullptr, w.lo, w.hi, true);
                if (hit.success && hit.edit_distance > k && !poly_support(hit)) hit.success = false;
                if (hit.success && !hit.query_complete && !hit.positions.empty()) {
                    const int len = hit.positions.front().second - hit.positions.front().first + 1;
                    const bool lclip = hit.ref_begin > 0;
                    const bool rclip = hit.ref_end >= 0 && hit.ref_end < m - 1;
                    // outside PARTIAL windows a clipped edge must face outward
                    const bool seen = w.status == Status::PARTIAL;
                    const int ed = hit.edit_distance;
                    const bool long_enough = seen ? ((len >= 12 && ed <= 1) || (len >= 16 && ed * 6 <= len))
                                                  : (len >= 16 && ed * 8 <= len);
                    if (!long_enough || (!seen && ((lclip && !info.left_outward) || (rclip && !info.right_outward))))
                        hit.success = false;
                }
                if (hit.success && w.status == Status::FULL && ctr)
                    ctr->full_as_missing.fetch_add(1, std::memory_order_relaxed);
                if (!hit.success && w.status == Status::PARTIAL && info.right_outward != info.left_outward) {
                    // the inward part is retained: prefix of a closer, suffix of an opener
                    const bool prefix = info.right_outward;
                    const int rf = std::max(0, w.retained_from);
                    const int rt = w.retained_to < 0 ? m : std::min(m, w.retained_to);
                    if (prefix ? rf == 0 : rt == m) {
                        const int t_max = std::min(m - 1, (prefix ? rt : m - rf) + 2);
                        hit = aligner.align_adapter_fragment(elem.seq, mutable_seq, w.lo, w.hi, prefix, whole_degraded ? 12 : 10, t_max);
                        is_fragment = hit.success;
                        if (hit.success && ctr) ctr->fragment_ok.fetch_add(1, std::memory_order_relaxed);
                    }
                }
                if (!hit.success && w.status != Status::FULL) {
                    const bool at3 = w.hi >= read_length - 6 && info.right_outward;
                    const bool at5 = !at3 && w.lo <= 7 && info.left_outward;
                    if (at3 || at5) {
                        hit = aligner.align_read_end_fragment(elem.seq, mutable_seq, at3, -1);
                        is_fragment = hit.success;
                    }
                }
            }
            if (!hit.success && w.status == Status::TRUNCATED_AT_READ_END && info.right_outward != info.left_outward) {
                const bool prefix = info.right_outward;
                const int rf = std::max(0, w.retained_from), rt = w.retained_to < 0 ? m : std::min(m, w.retained_to);
                const int t_max = (w.edits >= 0 && rt > rf) ? std::min(m - 1, rt - rf + 2) : m - 1;
                hit = aligner.align_adapter_fragment(elem.seq, mutable_seq, w.lo, w.hi, prefix, 10, t_max);
                is_fragment = hit.success;
                if (hit.success && ctr) ctr->fragment_ok.fetch_add(1, std::memory_order_relaxed);
            }
            // a fragment that kept the outward edge and lost the inward one
            if (!hit.success && (w.status == Status::PARTIAL || w.status == Status::TRUNCATED_AT_READ_END) &&
                info.right_outward != info.left_outward) {
                hit = aligner.align_adapter_fragment(elem.seq, mutable_seq, w.lo, w.hi, !info.right_outward, 12, m - 1, true);
                is_fragment = hit.success;
                if (hit.success && ctr) ctr->fragment_ok.fetch_add(1, std::memory_order_relaxed);
            }
            if (!hit.success || hit.positions.empty()) {
                if (retry && w.status != Status::TRUNCATED_AT_READ_END && (!info.poly_id.empty() || !info.partner_id.empty()))
                    retry->push_back({&elem, w});
                continue;
            }
            if (ctr) ctr->win_ok[st].fetch_add(1, std::memory_order_relaxed);
            const int adj_start = std::max(1, hit.positions.front().first);
            const int adj_end = std::min(read_length, hit.positions.front().second);
            if (adj_end < adj_start) continue;
            const std::string full_aligned_seq = mutable_seq.substr(adj_start - 1, adj_end - adj_start + 1);
            const std::string n_extracted = is_fragment ? std::string() : aligner.extract_n_masked_regions(hit, elem.seq);
            seq_element primary(
                elem.class_id, elem.global_class, hit.edit_distance, std::make_pair(adj_start, adj_end), "static",
                elem.order, elem.direction, true, std::nullopt,
                n_extracted.empty() ? full_aligned_seq : n_extracted, std::nullopt,
                n_extracted.empty() ? std::nullopt : std::optional<std::string>(full_aligned_seq));
            primary.query_complete = hit.query_complete;
            add_element(std::move(primary));
            std::fill(mutable_seq.begin() + adj_start - 1, mutable_seq.begin() + adj_end, 'X');
            if (verbose)
                log_verbose("CONCAT_HMM_STATIC " + elem.class_id + " " + concat_hmm::status_name(w.status) + " " +
                            std::to_string(adj_start) + ":" + std::to_string(adj_end) +
                            " edits=" + std::to_string(hit.edit_distance) +
                            (hit.query_complete ? "" : " partial"));
        }
    }

    // --concat-hmm: retry an unconfirmed barcode-side adapter at the layout spacer offset from an aligned
    // partner, then jointly with its adapter partner, then alone when the partner ran off the read end.
    // Runs after every other element of the strand.
    void retry_at_spacer_offset(const ReadElement& elem, const static_check_window& w, const read_streaming::sequence& read,
                                std::string& mutable_seq, const static_restriction& rs, aligner_tools& aligner, bool verbose) {
        if (!rs.info) return;
        const auto info_it = rs.info->elems.find(elem.class_id);
        if (info_it == rs.info->elems.end()) return;
        const concat_elem_info& info = info_it->second;
        concat_hmm_counters* ctr = rs.info->ctr;
        const int n = static_cast<int>(read.seq.length());
        const int m = info.m;
        const int k = elem.misalignment_threshold ? std::get<0>(*elem.misalignment_threshold)
                                                  : adapter_thresholds::fallback_max_edit_distance(elem.seq.length());
        const bool seen = w.status == concat_hmm::Status::PARTIAL && w.edits >= 0;
        const int tol = seen ? 8 : 4;
        const int wlo = w.lo - m, whi = w.hi + m, plo_lim = w.lo - 60, phi_lim = w.hi + 60;
        {   // already aligned here, e.g. by its partner's joint search
            auto range = by_id().equal_range(elem.class_id);
            for (auto it = range.first; it != range.second; ++it)
                if (it->position.first <= whi && wlo <= it->position.second) return;
        }
        if (ctr) ctr->spacer_retry.fetch_add(1, std::memory_order_relaxed);
        auto record = [&](const ReadElement& el, const static_alignments& hit, int hs, int he, const std::string& why) {
            const std::string full_aligned_seq = mutable_seq.substr(hs - 1, he - hs + 1);
            const std::string n_extracted = aligner.extract_n_masked_regions(hit, el.seq);
            seq_element primary(
                el.class_id, el.global_class, hit.edit_distance, std::make_pair(hs, he), "static",
                el.order, el.direction, true, std::nullopt,
                n_extracted.empty() ? full_aligned_seq : n_extracted, std::nullopt,
                n_extracted.empty() ? std::nullopt : std::optional<std::string>(full_aligned_seq));
            primary.query_complete = true;
            add_element(std::move(primary));
            std::fill(mutable_seq.begin() + hs - 1, mutable_seq.begin() + he, 'X');
            if (verbose)
                log_verbose("CONCAT_HMM_STATIC " + el.class_id + " SPACER_RETRY " + std::to_string(hs) + ":" + std::to_string(he) +
                            " edits=" + std::to_string(hit.edit_distance) + " " + why);
        };
        auto try_partner = [&](const std::string& pid, int glo, int ghi, bool after) -> bool {
            if (pid.empty()) return false;
            auto range = by_id().equal_range(pid);
            std::vector<std::pair<int, int>> partners;
            for (auto it = range.first; it != range.second; ++it) partners.push_back(it->position);
            const bool junction_side = after ? !rs.read_start : !rs.read_end;
            const int tin = pid == info.partner_id && junction_side ? std::min(tol, 4) : tol;
            for (const auto& pp : partners) {
                const int ps = pp.first, pe = pp.second;
                int lo, hi, edge_lo, edge_hi;             // edge: accepted adapter end (after) / start (before)
                if (after) { edge_hi = ps - 1 - glo + tin; edge_lo = ps - 1 - ghi - tol; lo = edge_lo - m + 1 - k; hi = edge_hi; }
                else { edge_lo = pe + 1 + glo - tin; edge_hi = pe + 1 + ghi + tol; lo = edge_lo; hi = edge_hi + m - 1 + k; }
                lo = std::max({lo, 1, plo_lim});
                hi = std::min({hi, n, phi_lim});
                if (hi - lo + 1 < m - k) continue;
                static_alignments hit = aligner.align_static_elements(elem.seq, mutable_seq, verbose, k, elem.masked_seq, verbose,
                                                                      lo, hi, nullptr, lo, hi, false);
                if (!hit.success || hit.positions.empty() || !hit.query_complete) continue;
                const int hs = std::max(1, hit.positions.front().first), he = std::min(n, hit.positions.front().second);
                if (he < hs || (after ? (he < edge_lo || he > edge_hi) : (hs < edge_lo || hs > edge_hi))) continue;
                record(elem, hit, hs, he, "partner=" + pid);
                if (ctr) ctr->spacer_retry_ok.fetch_add(1, std::memory_order_relaxed);
                if (pid == info.partner_id && info.bc_outer)
                    mark_weak_primer(elem.class_id, hs, he, ps, pe, glo, ghi, after);
                return true;
            }
            return false;
        };
        if (try_partner(info.poly_id, info.poly_gap_lo, info.poly_gap_hi, info.poly_after)) return;
        if (try_partner(info.partner_id, info.partner_gap_lo, info.partner_gap_hi, info.partner_after)) return;
        if (info.partner_id.empty() || by_id().count(info.partner_id)) return;
        auto try_joint = [&]() -> bool {
            const ReadElement* pel = nullptr;
            {
                auto found = rs.info->layout_by_id.find(info.partner_id);
                if (found != rs.info->layout_by_id.end()) pel = found->second;
            }
            if (!pel) return false;
            const auto pinfo_it = rs.info->elems.find(info.partner_id);
            if (pinfo_it == rs.info->elems.end()) return false;
            const int mp = pinfo_it->second.m;
            const int kp = pel->misalignment_threshold ? std::get<0>(*pel->misalignment_threshold)
                                                       : adapter_thresholds::fallback_max_edit_distance(pel->seq.length());
            int lo = std::max(1, wlo), hi = std::min(n, whi);
            if (hi - lo + 1 < m - k) return false;
            static_alignments he_hit = aligner.align_static_elements(elem.seq, mutable_seq, verbose, k, elem.masked_seq, verbose,
                                                                     lo, hi, nullptr, lo, hi, false);
            if (!he_hit.success || he_hit.positions.empty() || !he_hit.query_complete) return false;
            const int hs = std::max(1, he_hit.positions.front().first), he = std::min(n, he_hit.positions.front().second);
            if (he < hs) return false;
            const bool after = info.partner_after;
            const int glo = info.partner_gap_lo, ghi = info.partner_gap_hi;
            const bool junction_side = after ? !rs.read_start : !rs.read_end;
            const int tin = junction_side ? std::min(tol, 4) : tol;
            int edge_lo, edge_hi, plo, phi;  // edge: accepted partner start (after) / end (before)
            if (after) { edge_lo = he + 1 + glo - tin; edge_hi = he + 1 + ghi + tol; plo = edge_lo; phi = edge_hi + mp - 1 + kp; }
            else { edge_hi = hs - 1 - glo + tin; edge_lo = hs - 1 - ghi - tol; plo = edge_lo - mp + 1 - kp; phi = edge_hi; }
            plo = std::max(1, plo);
            phi = std::min(n, phi);
            bool named = false;
            for (const auto& pw : rs.windows)
                named = named || (pw.class_id == info.partner_id && pw.lo <= phi + mp && plo - mp <= pw.hi);
            if (!named || phi - plo + 1 < mp - kp) return false;
            static_alignments hp = aligner.align_static_elements(pel->seq, mutable_seq, verbose, kp, pel->masked_seq, verbose,
                                                                 plo, phi, nullptr, plo, phi, false);
            if (!hp.success || hp.positions.empty() || !hp.query_complete) return false;
            const int ps = std::max(1, hp.positions.front().first), pe = std::min(n, hp.positions.front().second);
            if (pe < ps || (ps <= he && hs <= pe)) return false;
            if (after ? (ps < edge_lo || ps > edge_hi) : (pe < edge_lo || pe > edge_hi)) return false;
            record(elem, he_hit, hs, he, "joint with " + info.partner_id);
            record(*pel, hp, ps, pe, "joint with " + elem.class_id);
            if (ctr) ctr->spacer_retry_ok.fetch_add(2, std::memory_order_relaxed);
            if (info.bc_outer) mark_weak_primer(elem.class_id, hs, he, ps, pe, glo, ghi, after);
            return true;
        };
        if (try_joint()) return;
        // The partner ran off the physical read end: accept the adapter alone at one edit above its anchor grade.
        const bool after = info.partner_after;
        if (after ? !rs.read_end : !rs.read_start) return;
        const auto pinfo_it = rs.info->elems.find(info.partner_id);
        if (pinfo_it == rs.info->elems.end()) return;
        const int mp = pinfo_it->second.m, glo = info.partner_gap_lo, t0 = 4;
        const int kr = std::min(k, std::max(1, std::min(4, m * 4 / 22)) + 1);
        int edge_lo, edge_hi, lo, hi;  // edge: accepted adapter end (after) / start (before)
        if (after) { edge_lo = n - mp - glo + t0 + 1; edge_hi = n - glo; lo = edge_lo - m + 1 - kr; hi = edge_hi; }
        else { edge_lo = glo + 1; edge_hi = mp + glo - t0; lo = edge_lo; hi = edge_hi + m - 1 + kr; }
        lo = std::max({lo, 1, wlo});
        hi = std::min({hi, n, whi});
        if (edge_hi < edge_lo || hi - lo + 1 < m - kr) return;
        static_alignments hit = aligner.align_static_elements(elem.seq, mutable_seq, verbose, kr, elem.masked_seq, verbose,
                                                              lo, hi, nullptr, lo, hi, false);
        if (!hit.success || hit.positions.empty() || !hit.query_complete) return;
        const int hs = std::max(1, hit.positions.front().first), he = std::min(n, hit.positions.front().second);
        if (he < hs || (after ? (he < edge_lo || he > edge_hi) : (hs < edge_lo || hs > edge_hi))) return;
        record(elem, hit, hs, he, "read-end partner=" + info.partner_id);
        if (ctr) {
            ctr->spacer_retry_ok.fetch_add(1, std::memory_order_relaxed);
            ctr->read_end_ok.fetch_add(1, std::memory_order_relaxed);
        }
    }

    // --concat-hmm: the barcode and UMI of a retry-only outer primer are placed from the partner anchor (the primer
    // laid at the layout spacer offset). Only an outer adapter: moving the inner one would move the cDNA edge.
    void mark_weak_primer(const std::string& class_id, int start, int end, int ps, int pe, int glo, int ghi, bool partner_after) {
        constexpr int min_shift = 3;
        int a = start, b = end;
        if (partner_after) { b = std::min(std::max(end, ps - 1 - ghi), ps - 1 - glo); a = b - (end - start); }
        else { a = std::min(std::max(start, pe + 1 + glo), pe + 1 + ghi); b = a + (end - start); }
        if (std::abs(a - start) < min_shift) return;
        hmm_weak_refs.push_back({class_id, {start, end}, {a, b}});
    }
    size_t weak_primer_count() const { return hmm_weak_refs.size(); }

    // Outer-adapter hits of a block: aligned elements, then refused clip hits (clipped = true, evidence only).
    template <class Fn>  // fn(position, clipped)
    void concat_hmm_each_outer_hit(const concat_bc_block& b, Fn&& fn) const {
        auto orng = by_id().equal_range(b.outer_id);
        for (auto o = orng.first; o != orng.second; ++o)
            if (o->position.first > 0) fn(o->position, false);
        for (const auto& h : hmm_clip_hits)
            if (h.id == b.outer_id) fn(h.pos, true);
    }
    size_t concat_hmm_classify_clip_hits(const concat_layout_info& info) {
        for (auto& h : hmm_clip_hits)
            for (const auto& b : info.bc_blocks) h.block_adapter = h.block_adapter || h.id == b.outer_id || h.id == b.inner_id;
        return hmm_clip_hits.size();
    }
    void mark_hmm_clip_check(uint8_t dirs) { hmm_clip_check = dirs; }
    // True when a written record's read element holds a refused clip hit, i.e. it would join two constructs.
    // Call after sigalign_filter, when the element positions are final.
    bool hmm_clip_hit_in_written_read() const {
        if (!hmm_clip_check || read_type == "filtered" || read_type == "skipped") return false;
        for (const auto& e : sig_elements) {
            if (e.global_class != "read" || e.position.first <= 0 || !e.seq.has_value() || e.seq->empty()) continue;
            if (read_type != "concatenate" && e.direction != read_type) continue;
            if (!(hmm_clip_check & (e.direction == "forward" ? 1 : 2))) continue;
            for (const auto& h : hmm_clip_hits)
                if (!h.block_adapter && h.pos.first >= e.position.first && h.pos.second <= e.position.second) return true;
        }
        return false;
    }

    // --concat-hmm: a block keeps its barcode only when its outer adapter is absent or its inner partner is at the
    // layout spacer offset (+-4 bp). Returns the number of rejected blocks.
    int concat_hmm_partner_rule(const concat_layout_info& info) {
        constexpr int tol = 4;
        int rejected = 0;
        for (const auto& b : info.bc_blocks) {
            bool outer = false, paired = false;
            concat_hmm_each_outer_hit(b, [&](const std::pair<int, int>& o, bool) {
                outer = true;
                auto irng = by_id().equal_range(b.inner_id);
                for (auto in = irng.first; in != irng.second; ++in) {
                    if (in->position.first <= 0) continue;
                    const int gap = b.inner_after ? in->position.first - o.second - 1
                                                  : o.first - in->position.second - 1;
                    paired = paired || (gap >= b.gap_lo - tol && gap <= b.gap_hi + tol);
                }
            });
            if (!outer || paired) continue;
            ++rejected;
            for (const auto& id : b.barcode_ids) hmm_reject_barcodes.push_back(id);
        }
        return rejected;
    }
    void concat_hmm_apply_barcode_rejections() {
        for (const auto& id : hmm_reject_barcodes)
            edit_elem(id, [](seq_element& x) { x.position = {-1, -1}; x.element_pass = false; });
    }

    void set_hmm_note(const std::string& n) { hmm_note = n; }
    void set_hmm_boundary_note(const std::string& n) { hmm_boundary_note = n; }
    const std::string& hmm_note_str() const { return hmm_note; }
    void mark_hmm_no_double_count() { hmm_no_double_count = true; }
    bool hmm_no_double_count_set() const { return hmm_no_double_count; }

    std::vector<concat_bc_state> concat_hmm_block_states(const concat_layout_info& info) const {
        constexpr int tol = 4;
        std::vector<concat_bc_state> out;
        for (const auto& b : info.bc_blocks) {
            concat_bc_state st;
            st.block = &b;
            auto lit = info.layout_by_id.find(b.outer_id);
            st.dir = lit != info.layout_by_id.end() ? lit->second->direction : std::string();
            auto irng = by_id().equal_range(b.inner_id);
            concat_hmm_each_outer_hit(b, [&](const std::pair<int, int>& o, bool clipped) {
                if (st.outer_s <= 0) { st.outer_s = o.first; st.outer_e = o.second; st.outer_clipped = clipped; }
                for (auto in = irng.first; in != irng.second; ++in) {
                    if (in->position.first <= 0) continue;
                    const int gap = b.inner_after ? in->position.first - o.second - 1
                                                  : o.first - in->position.second - 1;
                    if (gap >= b.gap_lo - tol && gap <= b.gap_hi + tol && !st.paired) {
                        st.paired = true;
                        st.outer_s = o.first; st.outer_e = o.second; st.outer_clipped = clipped;
                        st.inner_s = in->position.first; st.inner_e = in->position.second;
                    }
                }
            });
            if (st.inner_s <= 0)
                for (auto in = irng.first; in != irng.second; ++in)
                    if (in->position.first > 0) { st.inner_s = in->position.first; st.inner_e = in->position.second; break; }
            out.push_back(std::move(st));
        }
        return out;
    }

    // Plan a T / single-primer-end piece: TRIM at a bare outermost primer facing the one usable unit, else NORMAL.
    concat_piece_plan concat_hmm_plan_piece(const concat_layout_info& info, int len) const {
        concat_piece_plan plan;
        const auto states = concat_hmm_block_states(info);
        const concat_bc_state* F = nullptr;
        const concat_bc_state* R = nullptr;
        for (const auto& s : states) {
            const concat_bc_state*& slot = s.dir == "reverse" ? R : F;
            if (!slot || s.rank() > slot->rank()) slot = &s;
        }
        const bool u_f = F && F->usable(), u_r = R && R->usable();
        if (u_f && u_r) {
            plan.no_double_count = true;
            plan.note = "ART_BOTH";
            // check records only against a complete unit: an inner anchor alone is often a chance hit
            bool far_hit = false;
            for (const auto& h : hmm_clip_hits) far_hit = far_hit || !h.block_adapter;
            if (far_hit) plan.clip_check = static_cast<uint8_t>((R->strict() ? 1 : 0) | (F->strict() ? 2 : 0));
            return plan;
        }
        if (u_f || u_r) {
            const concat_bc_state* unit = u_f ? F : R;
            const concat_bc_state* other = u_f ? R : F;
            plan.dir = unit->dir;
            if (other && other->bare()) {
                const bool far_right = u_f && other->outer_s > unit->inner_e;
                const bool far_left = u_r && other->outer_e < unit->inner_s;
                bool outermost = far_right || far_left;
                if (outermost) {
                    for (const auto& e : sig_elements) {
                        if (e.type != "static" || e.position.first <= 0 || e.global_class == "start" ||
                            e.global_class == "stop" || e.global_class == "poly_tail") continue;
                        if (e.class_id == other->block->outer_id && e.position.first == other->outer_s) continue;
                        if ((far_right && e.position.second > other->outer_e) || (far_left && e.position.first < other->outer_s)) {
                            outermost = false;
                            break;
                        }
                    }
                }
                if (outermost) {
                    plan.kind = concat_piece_plan::TRIM;
                    plan.trim_lo = far_left ? other->outer_e : 0;
                    plan.trim_hi = far_right ? other->outer_s - 1 : len;
                    plan.note = u_f ? "P7_TRIM_F" : "P7_TRIM_R";
                    return plan;
                }
            }
            plan.note = "ART_SINGLE";
            return plan;
        }
        plan.note = "ART_NONE";
        return plan;
    }

    // Best concat_bc_state::rank on the given strand (any strand when dir is empty).
    int concat_hmm_unit_rank(const concat_layout_info& info, const std::string& dir) const {
        int r = 0;
        for (const auto& st : concat_hmm_block_states(info))
            if (dir.empty() || st.dir == dir) r = std::max(r, st.rank());
        return r;
    }

    void concat_hmm_reject_strand(const concat_layout_info& info, const std::string& dir) {
        for (const auto& b : info.bc_blocks) {
            auto lit = info.layout_by_id.find(b.outer_id);
            if (lit == info.layout_by_id.end() || lit->second->direction != dir) continue;
            for (const auto& id : b.barcode_ids) hmm_reject_barcodes.push_back(id);
        }
    }

    // Copy the parent's static alignments inside [offset0, offset0 + new_len), shifted to this molecule's
    // coordinates, plus virtual start / stop boundaries. No new alignment is run.
    void adopt_statics_from(const SigString& src, int offset0, int new_len, const ReadLayout& layout) {
        for (const auto& e : src.sig_elements) {
            if (e.type != "static" || e.global_class == "start" || e.global_class == "stop") continue;
            if (e.position.first <= 0 || e.position.first <= offset0 || e.position.second > offset0 + new_len) continue;
            seq_element c = e;
            c.position.first -= offset0;
            c.position.second -= offset0;
            add_element(std::move(c));
        }
        for (const auto& elem : layout.by_type()) {
            if (elem.type != "static" || (elem.global_class != "start" && elem.global_class != "stop")) continue;
            const int boundary = elem.global_class == "start" ? 0 : new_len + 1;
            add_element(seq_element(elem.class_id, elem.global_class, std::nullopt, {boundary, boundary},
                                    "static", elem.order, elem.direction));
        }
    }

    // Records to_fastqa_append would write: {forward, reverse}.
    std::pair<int, int> concat_hmm_pending_records() const {
        if (read_type == "skipped" || read_type == "filtered") return {0, 0};
        std::pair<int, int> out{0, 0};
        const bool bulk = additional_info == "bulk";
        for (const char* dir : {"forward", "reverse"}) {
            if (read_type != "concatenate" && read_type != dir) continue;
            bool has_bc = false, has_read = false;
            for (const auto& elem : sig_elements) {
                if (!elem.seq.has_value() || elem.direction != dir) continue;
                if (elem.global_class == "barcode") has_bc = true;
                else if (elem.global_class == "read" && !elem.seq->empty()) has_read = true;
                else if (bulk && elem.type == "static" && elem.original_seq.has_value()) has_bc = true;
            }
            if ((has_bc || bulk) && has_read) (std::string(dir) == "forward" ? out.first : out.second) += 1;
        }
        return out;
    }

    // Group occurrences, not alignment scores. This runs before variables are
    // instantiated, and never changes the legacy signature on an unresolved
    // read. In particular, an ordinary F/R duplet stays on its existing path.
    std::vector<static_segment> group_static_candidates(
        const ReadLayout& layout, const std::vector<seq_element>& additional_hits,
        bool verbose = false, int min_read_length = 0) const {
        if (additional_hits.empty()) return {};

        auto is_anchor = [&](const seq_element& elem) {
            if (elem.type != "static" || !elem.element_pass.value_or(false) ||
                elem.global_class == "start" || elem.global_class == "stop" ||
                elem.global_class == "poly_tail") return false;
            const auto found = layout.by_id().find(elem.class_id);
            if (found == layout.by_id().end() || !elem.edit_distance) return false;
            // Reject spans incompatible with a full-query alignment. This is
            // necessary but not sufficient to rule out SSW clipping; exact
            // query coverage is checked separately before a primary can root
            // a repeated orientation.
            const int span = elem.position.second - elem.position.first + 1;
            return std::abs(span - static_cast<int>(found->seq.size())) <= *elem.edit_distance;
        };
        // Most unused LOC results are an alternative to a failed/positional
        // primary, not a repeat. Avoid building grouping indices for them.
        bool repeated_anchor = false;
        for (size_t index = 0; index < additional_hits.size() && !repeated_anchor; ++index) {
            const auto& hit = additional_hits[index];
            if (!is_anchor(hit)) continue;
            auto same_id = by_id().equal_range(hit.class_id);
            auto separate = [&](const seq_element& other) {
                return other.class_id == hit.class_id && is_anchor(other) &&
                    (other.position.second < hit.position.first || hit.position.second < other.position.first);
            };
            for (auto it = same_id.first; it != same_id.second; ++it) repeated_anchor = repeated_anchor || separate(*it);
            for (size_t other = 0; other < index && !repeated_anchor; ++other) repeated_anchor = separate(additional_hits[other]);
        }
        if (!repeated_anchor) return {};
        const int read_length_floor = std::max(1, min_read_length >= 0
            ? min_read_length : calc_total_static_len(layout));

        std::map<std::string, std::vector<const seq_element*>> hits;
        for (const auto& elem : sig_elements) {
            if (is_anchor(elem)) hits[elem.class_id].push_back(&elem);
        }
        for (const auto& elem : additional_hits) {
            if (is_anchor(elem)) hits[elem.class_id].push_back(&elem);
        }
        for (auto& entry : hits) {
            auto& positions = entry.second;
            std::stable_sort(positions.begin(), positions.end(), [](const auto* a, const auto* b) {
                return a->position < b->position;
            });
            std::vector<const seq_element*> distinct;
            int cluster_end = -1;
            for (const auto* hit : positions) {
                if (distinct.empty() || hit->position.first > cluster_end) {
                    distinct.push_back(hit);
                } else if (hit->edit_distance.value_or(INT_MAX) <
                           distinct.back()->edit_distance.value_or(INT_MAX)) {
                    distinct.back() = hit;
                }
                cluster_end = std::max(cluster_end, hit->position.second);
            }
            positions = std::move(distinct);
        }

        std::vector<static_segment> segments;
        std::map<std::string, bool> incomplete_directions;
        std::map<std::string, std::vector<std::pair<std::string, std::string>>> anchor_patterns;
        for (const std::string direction : {"forward", "reverse"}) {
            std::vector<const ReadElement*> anchors;
            int first_read = INT_MAX, last_read = -1;
            for (const auto& elem : layout.by_order()) {
                if (elem.direction != direction) continue;
                if (elem.type == "variable" && elem.global_class == "read") {
                    first_read = std::min(first_read, elem.order);
                    last_read = std::max(last_read, elem.order);
                }
                if (elem.type == "static" && elem.global_class != "start" &&
                    elem.global_class != "stop" && elem.global_class != "poly_tail") {
                    anchors.push_back(&elem);
                }
            }
            // A single adapter or a homopolymer is insufficient evidence for
            // a new molecule. Require real layout anchors around the cDNA.
            if (anchors.size() < 2 || first_read == INT_MAX || first_read != last_read ||
                anchors.front()->order >= first_read || anchors.back()->order <= last_read) continue;
            bool unbounded_variables = false;
            for (const auto& elem : layout.by_order()) {
                if (elem.direction == direction && elem.type == "variable" &&
                    (elem.order < anchors.front()->order || elem.order > anchors.back()->order)) {
                    unbounded_variables = true;
                }
            }
            if (unbounded_variables) continue;
            auto& pattern = anchor_patterns[direction];
            for (const auto* anchor : anchors) {
                pattern.emplace_back(seq_utils::remove_rc(anchor->class_id),
                    direction == "reverse" ? seq_utils::revcomp(anchor->seq) : anchor->seq);
            }
            std::sort(pattern.begin(), pattern.end());

            std::vector<int> min_gaps(anchors.size(), 0);
            for (size_t index = 1; index < anchors.size(); ++index) {
                for (const auto& elem : layout.by_order()) {
                    if (elem.direction != direction || elem.type != "variable" ||
                        elem.order <= anchors[index - 1]->order || elem.order >= anchors[index]->order) continue;
                    int length = elem.expected_length.value_or(0);
                    if (!elem.length_candidates.empty()) {
                        length = *std::min_element(elem.length_candidates.begin(), elem.length_candidates.end());
                    }
                    // A necessarily too-short cDNA cannot lend support to a
                    // cut of an otherwise valid parent. This is an optimistic
                    // bound (poly-tails are not charged here); normal extraction
                    // and filtering still decide whether the child passes.
                    min_gaps[index] += std::max(elem.global_class == "read" ? read_length_floor : 0, length);
                }
            }
            const auto& heads = hits[anchors.front()->class_id];
            bool incomplete = false;
            for (size_t head = 0; head < heads.size(); ++head) {
                const int limit = head + 1 < heads.size() ? heads[head + 1]->position.first : sequence_length + 1;
                std::vector<const seq_element*> chain{heads[head]};
                bool ambiguous = false;
                for (size_t index = 1; index < anchors.size(); ++index) {
                    const seq_element* next = nullptr;
                    bool tied = false;
                    for (const auto* hit : hits[anchors[index]->class_id]) {
                        if (hit->position.first <= chain.back()->position.second + min_gaps[index] ||
                            hit->position.second >= limit) continue;
                        if (!next || *hit->edit_distance < *next->edit_distance) {
                            next = hit;
                            tied = false;
                        } else if (*hit->edit_distance == *next->edit_distance) {
                            tied = true;
                        }
                    }
                    if (!next) break;
                    chain.push_back(next);
                    ambiguous = ambiguous || tied;
                }
                if (chain.size() != anchors.size()) {
                    incomplete = true;
                    continue;
                }
                static_segment segment{{chain.front()->position.first, chain.back()->position.second}, direction, {}};
                segment.ambiguous = ambiguous;
                for (const auto* hit : chain) {
                    segment.elements.push_back(*hit);
                    segment.edit_distance += *hit->edit_distance;
                }
                for (const auto& elem : sig_elements) {
                    if (elem.global_class == "poly_tail" && elem.direction == direction &&
                        elem.position.first >= segment.position.first && elem.position.second <= segment.position.second) {
                        segment.elements.push_back(elem);
                    }
                }
                segments.push_back(std::move(segment));
            }
            incomplete_directions[direction] = incomplete;
        }
        // Let complete, stronger chains outvote incidental opposite-strand
        // matches in their interval. Equal-strength conflicting evidence is
        // unresolved, not a reason to invent a cut. Custom asymmetric layouts
        // cannot compare raw edit sums across different anchor sets.
        const bool comparable_orientations = anchor_patterns["forward"] == anchor_patterns["reverse"];
        std::sort(segments.begin(), segments.end(), [](const auto& a, const auto& b) {
            return std::tie(a.edit_distance, a.position) < std::tie(b.edit_distance, b.position);
        });
        std::vector<static_segment> selected;
        std::map<std::string, size_t> direction_counts;
        for (auto& segment : segments) {
            bool overlaps_better = false;
            for (const auto& other : selected) {
                if (segment.position.first <= other.position.second && other.position.first <= segment.position.second) {
                    if (segment.edit_distance == other.edit_distance ||
                        (segment.direction != other.direction && !comparable_orientations)) {
                        if (verbose) log_verbose("CONCAT_UNRESOLVED equal-strength/incomparable overlapping chains");
                        return {};
                    }
                    overlaps_better = true;
                }
            }
            if (overlaps_better) continue;
            if (segment.ambiguous || incomplete_directions[segment.direction]) {
                if (verbose) log_verbose("CONCAT_UNRESOLVED ambiguous/incomplete anchor chain: " + segment.direction);
                return {};
            }
            ++direction_counts[segment.direction];
            selected.push_back(std::move(segment));
        }
        if (direction_counts["forward"] < 2 && direction_counts["reverse"] < 2) return {};
        // A repeat must extend an accepted full-query primary, not consist
        // entirely of LOC alternatives whose legacy SSW validation failed.
        // Otherwise weak repeated motifs can replace a valid opposite strand.
        for (const auto& entry : direction_counts) {
            if (entry.second < 2) continue;
            bool supported = false;
            for (const auto& primary : sig_elements) {
                if (primary.direction != entry.first || !primary.query_complete || !is_anchor(primary)) continue;
                for (const auto& segment : selected) {
                    if (segment.direction != entry.first) continue;
                    for (const auto& hit : segment.elements) {
                        if (hit.class_id == primary.class_id &&
                            hit.position.first <= primary.position.second &&
                            primary.position.first <= hit.position.second) supported = true;
                    }
                }
            }
            if (!supported) {
                if (verbose) log_verbose("CONCAT_UNRESOLVED repeat without primary support: " + entry.first);
                return {};
            }
        }
        std::sort(selected.begin(), selected.end(), [](const auto& a, const auto& b) { return a.position < b.position; });
        // Clipped/one-sided primaries cannot justify a cut, but they can still
        // support legacy extraction. Do not abandon their intervals. Coverage
        // is strand-independent so incidental opposite hits inside (or across
        // touching) accepted children do not veto otherwise resolved groups.
        for (const auto& primary : sig_elements) {
            if (primary.type != "static" || !primary.element_pass.value_or(false) ||
                primary.global_class == "start" || primary.global_class == "stop" ||
                primary.global_class == "poly_tail" || layout.by_id().find(primary.class_id) == layout.by_id().end() ||
                primary.position.first < 1 || primary.position.second < primary.position.first ||
                primary.position.second > sequence_length) continue;
            int uncovered = primary.position.first;
            for (const auto& segment : selected) {
                if (segment.position.second < uncovered) continue;
                if (segment.position.first > uncovered) break;
                uncovered = segment.position.second + 1;
                if (uncovered > primary.position.second) break;
            }
            if (uncovered <= primary.position.second) {
                if (verbose) log_verbose("CONCAT_UNRESOLVED uncovered primary: " + primary.class_id);
                return {};
            }
        }
        if (verbose) {
            log_verbose("CONCAT_GROUPS " + std::to_string(selected.size()) + " (reused static alignments)");
            for (const auto& segment : selected) {
                log_verbose("CONCAT_SEGMENT " + segment.direction + " " +
                    std::to_string(segment.position.first) + ":" + std::to_string(segment.position.second) +
                    " edits=" + std::to_string(segment.edit_distance));
            }
        }
        return selected;
    }

    // Materialize only an accepted occurrence. No child adapter search is
    // performed; the only possible alignment here is a local traceback for
    // an indel-containing N-masked capture whose LOC result has no path.
    bool assign_static_segment(const static_segment& segment, const ReadLayout& layout,
                               const read_streaming::sequence& read) {
        from_concatemer = true;
        const int offset = segment.position.first - 1;
        aligner_tools aligner;
        for (auto elem : segment.elements) {
            elem.position.first -= offset;
            elem.position.second -= offset;
            const auto layout_it = layout.by_id().find(elem.class_id);
            if (layout_it != layout.by_id().end() && layout_it->seq.find('N') != std::string::npos &&
                !elem.original_seq) {
                static_alignments capture;
                capture.seq = read.seq.substr(elem.position.first - 1, elem.position.second - elem.position.first + 1);
                if (elem.edit_distance.value_or(-1) != 0) {
                    const EdlibEqualityPair equalities[] = {
                        {'A','a'}, {'C','c'}, {'T','t'}, {'G','g'},
                        {'A','x'}, {'C','x'}, {'T','x'}, {'G','x'},
                        {'N','A'}, {'N','C'}, {'N','T'}, {'N','G'}, {'N','x'}};
                    auto result = edlibAlign(layout_it->seq.c_str(), layout_it->seq.size(),
                        capture.seq.c_str(), capture.seq.size(),
                        edlibNewAlignConfig(elem.edit_distance.value_or(-1), EDLIB_MODE_NW,
                                            EDLIB_TASK_PATH, equalities, 13));
                    if (result.status == EDLIB_STATUS_OK && result.editDistance >= 0 && result.alignment) {
                        char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
                        if (cigar) { capture.cigar = cigar; free(cigar); }
                    }
                    edlibFreeAlignResult(result);
                    if (capture.cigar.empty()) return false;
                }
                elem.original_seq = capture.seq;
                elem.seq = aligner.extract_n_masked_regions(capture, layout_it->seq);
            }
            add_element(std::move(elem));
        }
        for (const auto& elem : layout.by_type()) {
            if (elem.type != "static" || (elem.global_class != "start" && elem.global_class != "stop")) continue;
            const int boundary = elem.global_class == "start" ? 0 : sequence_length + 1;
            add_element(seq_element(elem.class_id, elem.global_class, std::nullopt, {boundary, boundary},
                                    "static", elem.order, elem.direction));
        }
        return true;
    }

/**
 * @brief Master function for variable sequence alignment and detection.
 * @param read `read_streaming::sequence` object containing the read sequence
 * @param layout `ReadLayout` object containing layout elements
 * @param verbose print verbose output
 * 
 * @brief This function embeds variable sequence elements (e.g., barcodes, reads)
 * within a sequencing read based on the provided layout. It first separates variable
 * elements into read and non-read categories. Non-read variables are further partitioned
 * by direction (forward/reverse) and processed accordingly. The function generates
 * variable elements by referencing previously aligned static elements to determine
 * their positions within the read. Elements to the left of a read are processed from top-down in terms
 * of order, while those to the right are processed bottom-up.
 * Read variables are processed last to ensure accurate positioning. 
 */
    void sigalign_variable(const read_streaming::sequence &read, const ReadLayout& layout, bool verbose) {
        if(verbose){
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "\n=== Starting variable alignment for " << sequence_id << " ===" << std::endl;
                RAD_COUT << oss.str();
            }
        }

        // Split variables into read and non-read
        std::vector<const ReadElement*> non_read_vars;
        std::vector<const ReadElement*> read_vars;
        
        auto& type_index = layout.by_type();
        auto variable_range = type_index.equal_range("variable");
        
        for (auto it = variable_range.first; it != variable_range.second; ++it) {
            // Add a placeholder element into sig_elements container
            add_element(seq_element(
                it->class_id,
                it->global_class,
                std::nullopt,
                {-1, -1},
                "variable",
                it->order,
                it->direction,
                std::nullopt,
                std::nullopt,
                std::nullopt
            ));
        
            if (it->global_class == "read") {
                read_vars.push_back(&(*it));
            } else {
                non_read_vars.push_back(&(*it));
            }
        }
        
        // Partition non-read variables by direction
        std::vector<const ReadElement*> non_read_forward;
        std::vector<const ReadElement*> non_read_reverse;
        for (const auto* var : non_read_vars) {
            if (var->direction == "reverse") {
                non_read_reverse.push_back(var);
            } else {
                non_read_forward.push_back(var);
            }
        }
        
        int positioned_count = 0;

        // --concat-hmm: map_positions sees retry-only primers at their partner-anchored position (copies).
        std::vector<seq_element> hmm_anchored;
        std::vector<const seq_element*> hmm_anchored_of;
        if (!hmm_weak_refs.empty()) {
            hmm_anchored.reserve(hmm_weak_refs.size());  // no reallocation: pointers into it must stay valid
            for (const auto& elem : sig_elements) {
                if (elem.type != "static") continue;
                for (const auto& w : hmm_weak_refs) {
                    if (w.id != elem.class_id || w.aligned != elem.position || hmm_anchored.size() == hmm_weak_refs.size()) continue;
                    hmm_anchored.push_back(elem);
                    hmm_anchored.back().position = w.anchored;
                    hmm_anchored_of.push_back(&elem);
                    break;
                }
            }
        }
        auto ref_ptr = [&](const seq_element& elem) -> const seq_element* {
            for (size_t i = 0; i < hmm_anchored_of.size(); ++i)
                if (hmm_anchored_of[i] == &elem) return &hmm_anchored[i];
            return &elem;
        };
        
        // Process forward non-read variables
        std::multimap<std::string, const seq_element*> static_refs;
        for (const auto& elem : sig_elements) {
            if (elem.type == "static") {
                static_refs.insert({elem.class_id, hmm_weak_refs.empty() ? &elem : ref_ptr(elem)});
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Found static reference: " << elem.class_id 
                                  << " at position " << elem.position.first 
                                  << ":" << elem.position.second 
                                  << " with edit distance " << (elem.edit_distance ? std::to_string(elem.edit_distance.value()) : "none")
                                  << " and sequence " << (elem.seq ? elem.seq.value() : "none") 
                                  << std::endl;
                        RAD_COUT << oss.str();
                    }
                }
            }
        }
        
        // Process forward non-read variables
        for (const auto* var : non_read_forward) {
            if (generate_variable_elements(var, static_refs, read.seq, sig_elements, verbose)) {
                positioned_count++;
                // Add the newly generated variable element(s) to the static_refs multimap so that secondary positions can be calculated
                auto found = sig_elements.get<sig_id_tag>().find(var->class_id);
                if (found != sig_elements.get<sig_id_tag>().end()) {
                    static_refs.insert({var->class_id, &(*found)});
                }
            }
        }
        
        // Process reverse non-read variables in reverse order
        for (auto it = non_read_reverse.rbegin(); it != non_read_reverse.rend(); ++it) {
            const auto* var = *it;
            if (generate_variable_elements(var, static_refs, read.seq, sig_elements, verbose)) {
                positioned_count++;
                auto found = sig_elements.get<sig_id_tag>().find(var->class_id);
                if (found != sig_elements.get<sig_id_tag>().end()) {
                    static_refs.insert({var->class_id, &(*found)});
                }
            }
        }
        
        //  total ref container
        std::multimap<std::string, const seq_element*> total_refs;
        for (const auto& elem : sig_elements) {
            total_refs.insert({elem.class_id, hmm_anchored_of.empty() ? &elem : ref_ptr(elem)});
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Found reference: " << elem.class_id 
                              << " at position " << elem.position.first 
                              << ":" << elem.position.second << std::endl;
                    RAD_COUT << oss.str();
                }
            }
        }
        // Process read variables last
        for (const auto* var : read_vars) {
            if (generate_variable_elements(var, total_refs, read.seq, sig_elements, verbose)) {
                positioned_count++;
            }
        }
        if (verbose) {
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "Successfully positioned " << positioned_count << " variable elements\n"
                 << "Final sigstring elements: " << sig_elements.size() << std::endl;
                RAD_COUT << oss.str();
            }
        }
    }
   
/**
 * @brief Master function for filtering and barcode correction.
 * @param read `read_streaming::sequence` object containing the read sequence
 * @param layout `ReadLayout` object containing layout elements
 * @param gen_mut general mutation rate for barcode correction
 * @param verbose print verbose output
 * @param mode correction mode for barcode correction (offensive or defensive, which whitelist to check first)
 * 
 * @brief Master function for barcode correction and sequence-level filtering. 
 * This function processes variable sequence elements (e.g., barcodes, reads)
 * within a sequencing read based on the provided layout. It first groups variable
 * elements by direction (forward/reverse) and performs basic validation, including
 * position checks and overlap detection. Valid directions are then assessed to
 * determine the preliminary read type. Barcode correction is applied to barcode
 * elements, and the direction validity is updated accordingly. Finally, the
 * final read type is determined based on the updated direction validity. 
 * This function also includes concatenate resolution logic.
 */
    void sigalign_filter(const read_streaming::sequence &read, const ReadLayout& layout, int gen_mut, bool verbose,
        std::string mode, const std::string& joint_bc_mode, int min_read_length = -1
    ) {
        
        auto direction_elements = group_directionally();
        
        std::string filtered_because = "";
        std::map<std::string, bool> direction_valid;
        std::map<std::string, int> pass_counts;
        
        if (verbose) {
            log_verbose("Starting sigalign_filter processing");
        }
        
        // part 1: Process each direction for basic validation and masking
        for (auto& [direction, elements] : direction_elements) {
            direction_valid[direction] = process_direction_basic(
                direction, elements, read, layout, filtered_because, verbose, min_read_length
            );
            
            if (!direction_valid[direction]) {
                pass_counts[direction] = 0;
                continue;
            }
            
            // Count valid reads
            pass_counts[direction] = count_valid_reads(elements);
        }
        
        // part 2: Determine read direction based on non-barcode elements
        auto [final_direction, read_type_preliminary] = determine_read_direction(
            direction_valid, pass_counts, verbose
        );
        
        // part 3: Barcode correction (only for valid directions)
        for (auto& [direction, elements] : direction_elements) {
            if (!direction_valid[direction]){
                continue;
            }
            bool barcode_success = process_barcodes_for_direction(
                direction, elements, read, layout, gen_mut, mode, joint_bc_mode, verbose
            );
            
            if (!barcode_success) {
                direction_valid[direction] = false;
                pass_counts[direction] = 0;
                filtered_because += ":" + direction + ":BARCODE_CORRECTION_FAILED";
                set_info(filtered_because);
            } else {
                // Recount elements including successfully corrected barcodes
                pass_counts[direction] = count_all_valid_elements(elements);
            }
        }
        
        // part 4: Final read type determination
        determine_final_read_type(direction_valid, pass_counts, verbose);
        
        // part 5: Update barcode counts
        update_bc_counts(*this, layout, verbose);
        
        if (verbose) {
            log_final_results(direction_valid, pass_counts);
        }
    }
    
/**
 * @brief The main function to perform sigalign on a FASTQ file, producing aligned and demultiplexed outputs.
 * @param fastq_path `std::string` path to the input FASTQ file
 * @param layout `ReadLayout` object defining the layout of the reads
 * @param output_prefix `std::string` prefix for the output files
 * @param gen_mut `std::optional<int>` general mutation rate for barcode correction, default is 2
 * @param verbose flag to enable verbose logging
 * @param num_threads `int` number of threads to use for parallel processing
 * @param chunk_size `size_t` number of reads to process in each chunk
 * @param max_reads `size_t` maximum number of reads to process from the input file
 * @param write_debug flag to enable writing debug output files
 * @param mode `std::string` correction mode for barcode correction (offensive or defensive, which whitelist to check first)
 */
    // --concat-hmm per-read router: segments the read with the HMM, aligns each construct in its windows and
    // passes each molecule to process_molecule. Returns false when the caller must run the default path
    // (guard failed, k = 0, or an abstain with --concat-hmm-abstain=legacy); true when handled, possibly
    // without a record.
    template <class mol_fn>
    static bool concat_hmm_route(const read_streaming::sequence& read, const ReadLayout& layout,
                                 const concat_layout_info& info, mol_fn&& process_molecule, bool verbose,
                                 const std::function<void(const read_streaming::sequence&)>& preserve_unresolved = {}) {
        thread_local concat_hmm::Scratch scratch;
        thread_local concat_hmm::Result res;
        concat_hmm_counters& ctr = *info.ctr;
        const int len = static_cast<int>(read.seq.size());
        ctr.reads.fetch_add(1, std::memory_order_relaxed);
        if (layout.concat_model->guard_failed) {
            ctr.legacy_guard.fetch_add(1, std::memory_order_relaxed);
            return false;
        }
        {
#ifdef RAD_STAGE_TIMERS
            rad_stage_timers::scope t_hmm(rad_stage_timers::hmm_ns);
#endif
            concat_hmm::segment(*layout.concat_model, read.seq.data(), len, scratch, res);
        }
        // A confident repeat transition need not locate the physical junction.
        // Keep the original parent, and do not assign the unresolved interval to
        // either barcode merely by choosing an arithmetic midpoint.
        std::vector<uint8_t> boundary_blocked(static_cast<size_t>(std::max(0, res.k)), 0);
        bool has_unresolved_boundary = res.boundary_unresolved_seen;
        for (size_t j = 0; j < res.cut_flags.size() && j + 1 < boundary_blocked.size(); ++j) {
            if (!(res.cut_flags[j] & concat_hmm::BOUNDARY_UNRESOLVED)) continue;
            boundary_blocked[j] = boundary_blocked[j + 1] = 1;
            has_unresolved_boundary = true;
        }
        if (has_unresolved_boundary) {
            ctr.boundary_unresolved_parents.fetch_add(1, std::memory_order_relaxed);
            ctr.boundary_unresolved_children.fetch_add(
                static_cast<size_t>(std::count(boundary_blocked.begin(), boundary_blocked.end(), uint8_t{1})),
                std::memory_order_relaxed);
            if (preserve_unresolved) preserve_unresolved(read);
        }
        const char* legacy = nullptr;
        const char* drop = nullptr;  // written without any record
        bool junction_abstain = false;
        // Legacy fallback must not reintroduce a cut through an interval we
        // have explicitly classified as coordinate-unresolved.
        const bool abstain_legacy = layout.concat_abstain_legacy && !has_unresolved_boundary;
        if (res.flags & concat_hmm::RES_GUARD_FAILED) { legacy = "guard"; ctr.legacy_guard.fetch_add(1, std::memory_order_relaxed); }
        else if (res.k <= 0) { legacy = "k0"; ctr.legacy_k0.fetch_add(1, std::memory_order_relaxed); }
        else if (res.flags & concat_hmm::RES_TOO_MANY) {
            if (abstain_legacy) { legacy = "too_many"; ctr.legacy_too_many.fetch_add(1, std::memory_order_relaxed); }
            else { drop = "too_many"; ctr.drop_too_many.fetch_add(1, std::memory_order_relaxed); }
        } else if (res.flags & concat_hmm::RES_OPPOSITE_EVIDENCE) {
            if (abstain_legacy) { legacy = "abstain"; ctr.legacy_abstain.fetch_add(1, std::memory_order_relaxed); }
            else { drop = "abstain"; ctr.drop_abstain.fetch_add(1, std::memory_order_relaxed); }
        } else if (res.flags & concat_hmm::RES_ABSTAIN) {
            // k >= 2: abstain only on the uncertain junctions (resolved below)
            if (abstain_legacy) { legacy = "abstain"; ctr.legacy_abstain.fetch_add(1, std::memory_order_relaxed); }
            else if (res.k == 1) { drop = "abstain"; ctr.drop_abstain.fetch_add(1, std::memory_order_relaxed); }
            else junction_abstain = true;
        } else if (res.k == 1 && (res.segs[0].flags & concat_hmm::SEG_DOUBLE_BC)) {
            // D construct with no clear cDNA direction: write nothing rather than guess the unit
            drop = "d";
            ctr.drop_d_reads.fetch_add(1, std::memory_order_relaxed);
        }
        if (verbose) {
            std::ostringstream oss;
            oss << "CONCAT_HMM " << read.id << " k=" << res.k << " p_single=" << res.p_single
                << " strand=" << res.strand_call << " flags=" << res.flags << " cuts=";
            for (size_t c = 0; c < res.cuts.size(); ++c) oss << (c ? ";" : "") << res.cuts[c];
            oss << " branch=" << (legacy ? std::string("legacy:") + legacy
                                  : drop ? std::string("drop:") + drop : (res.k == 1 ? "single" : "split"));
            #pragma omp critical
            RAD_COUT << oss.str() << std::endl;
        }
        if (legacy) return false;
        if (drop) return true;

        // check regions of construct c -> windows in child coordinates
        auto restriction_for = [&](int c, int child_start, int child_end) {
            static_restriction r;
            const char strand = res.segs[c].strand;
            r.direction = strand == 'F' ? "forward" : "reverse";
            r.info = &info;
            r.read_start = child_start == 0;
            r.read_end = child_end == len;
            const int child_len = child_end - child_start;
            for (const auto& chk : res.checks) {
                if (chk.construct != c || chk.element < 0 || chk.element >= static_cast<int>(info.id_F.size())) continue;
                const std::string& id = strand == 'F' ? info.id_F[chk.element] : info.id_R[chk.element];
                if (id.empty()) continue;
                static_check_window w;
                w.class_id = id;
                w.lo = std::max(chk.start, child_start) - child_start + 1;
                w.hi = std::min(chk.end, child_end) - child_start;
                if (w.hi < w.lo) {
                    if (chk.status != concat_hmm::Status::TRUNCATED_AT_READ_END) continue;
                    w.lo = w.hi = std::max(1, std::min(child_len, w.hi));  // read-end marker
                }
                w.status = chk.status;
                w.retained_from = chk.retained_from;
                w.retained_to = chk.retained_to;
                w.edits = chk.edits;
                r.windows.push_back(std::move(w));
            }
            return r;
        };
        const bool split = res.k >= 2;
        if (!split) ctr.single.fetch_add(1, std::memory_order_relaxed);
        else {
            ctr.split.fetch_add(1, std::memory_order_relaxed);
            ctr.children.fetch_add(static_cast<size_t>(res.k), std::memory_order_relaxed);
        }
        for (uint8_t kd : res.cut_kind) if (kd < 6) ctr.cut_kind[kd].fetch_add(1, std::memory_order_relaxed);
        struct piece_t {
            int c = 0, start = 0, end = 0;
            bool plain = true, trimmed = false, suppressed = false;
            char strand = '?';
            read_streaming::sequence rs;
            SigString sig;
        };
        std::vector<piece_t> pieces;
        pieces.reserve(static_cast<size_t>(res.k));
        const bool win = res.cut_lo.size() == res.cuts.size() && res.cut_hi.size() == res.cuts.size();
        for (int c = 0; c < res.k; ++c) {
            if (boundary_blocked[static_cast<size_t>(c)]) {
                continue;
            }
            if (res.segs[c].flags & concat_hmm::SEG_SAME_MOLECULE_R) ctr.same_molecule.fetch_add(1, std::memory_order_relaxed);
            // children may overlap at a cut (cut_lo / cut_hi) so each keeps its whole barcode block
            int start = c == 0 ? 0 : (win ? res.cut_lo[c - 1] : res.cuts[c - 1]);
            int end = c + 1 == res.k ? len : (win ? res.cut_hi[c] : res.cuts[c]);
            // D kept unit: the piece ends at the removed unit's inner anchor, so that unit is never inside it
            const bool d_kept = (res.segs[c].flags & concat_hmm::SEG_D_KEPT) != 0;
            const bool d_left = d_kept && res.segs[c].strand == 'F';
            if (d_kept) { if (d_left) end = std::min(end, res.segs[c].end); else start = std::max(start, res.segs[c].start); }
            if (end <= start) continue;
            if (res.segs[c].flags & concat_hmm::SEG_DOUBLE_BC) {
                ctr.drop_d_children.fetch_add(1, std::memory_order_relaxed);
                continue;
            }
            const bool piece_of_read = split || d_kept;
            piece_t p;
            p.c = c; p.start = start; p.end = end;
            p.plain = concat_plain_construct(res.segs[c]);
            p.strand = res.segs[c].strand;
            p.rs = read_streaming::sequence{split ? read.id + "/seg" + std::to_string(c + 1) : read.id, read.comment,
                                            piece_of_read ? read.seq.substr(start, end - start) : read.seq, read.plus,
                                            read.is_fastq ? (piece_of_read ? read.qual.substr(start, end - start) : read.qual) : "", read.is_fastq};
            p.sig = SigString(p.rs.id, end - start, "undefined", layout.sequencing_type);
            if (split) p.sig.mark_concatemer();
            if (piece_of_read) p.sig.set_parent_frame(start, end);
            {
#ifdef RAD_STAGE_TIMERS
                rad_stage_timers::scope t_static(rad_stage_timers::windowed_static_ns);
#endif
                if (p.plain) {
                    // at a fold-back or strand-flip cut the read element may run to the child boundary
                    auto cut_edge = [&](int j) {
                        return j >= 0 && j < static_cast<int>(res.cut_kind.size()) &&
                               (res.cut_kind[j] == concat_hmm::CUT_FOLDBACK || res.cut_kind[j] == concat_hmm::CUT_STRAND_FLIP);
                    };
                    const bool fold_start = c > 0 && ((res.segs[c - 1].flags & concat_hmm::SEG_SAME_MOLECULE_R) || cut_edge(c - 1));
                    const bool fold_end = c + 1 < res.k && ((res.segs[c].flags & concat_hmm::SEG_SAME_MOLECULE_R) || cut_edge(c));
                    // likewise at the edge of a D kept unit, where no layout reference remains
                    if (fold_start || fold_end || d_kept)
                        p.sig.set_fold_edges(fold_start || (d_kept && !d_left), fold_end || d_left,
                                             res.segs[c].strand == 'F' ? "forward" : "reverse", !d_kept);
                    if (d_kept) {
                        p.sig.set_hmm_note(d_left ? "D_KEEP_L" : "D_KEEP_R");
                        (d_left ? ctr.d_kept_left : ctr.d_kept_right).fetch_add(1, std::memory_order_relaxed);
                    } else if (res.segs[c].flags & concat_hmm::SEG_D_SPLIT) {
                        p.sig.set_hmm_note("D_SPLIT");
                        if (res.segs[c].strand == 'F') ctr.d_split.fetch_add(1, std::memory_order_relaxed);  // count once per D
                    }
                    const static_restriction r = restriction_for(c, start, end);
                    p.sig.sigalign_static(p.rs, layout, verbose, nullptr, &r);
                } else {
                    // T / single-primer end: full static search in the piece, then the piece and partner rules
                    ctr.pieces.fetch_add(1, std::memory_order_relaxed);
                    if (split) ctr.children_full_search.fetch_add(1, std::memory_order_relaxed);
                    else ctr.art_reads.fetch_add(1, std::memory_order_relaxed);
                    const static_piece_frame frame{&read.seq, start};
                    p.sig.sigalign_static(p.rs, layout, verbose, nullptr, nullptr, &frame);
                    if (const size_t nclip = p.sig.concat_hmm_classify_clip_hits(info)) ctr.clip_hits.fetch_add(nclip, std::memory_order_relaxed);
                    const concat_piece_plan plan = p.sig.concat_hmm_plan_piece(info, end - start);
                    const std::string other = plan.dir == "forward" ? "reverse" : "forward";
                    if (plan.kind == concat_piece_plan::TRIM) {
                        const int tl = plan.trim_lo, th = plan.trim_hi;
                        read_streaming::sequence rs2{p.rs.id, read.comment, p.rs.seq.substr(tl, th - tl), p.rs.plus,
                                                     p.rs.is_fastq ? p.rs.qual.substr(tl, th - tl) : "", p.rs.is_fastq};
                        SigString sig2(rs2.id, th - tl, "undefined", layout.sequencing_type);
                        if (split) sig2.mark_concatemer();
                        sig2.set_parent_frame(start + tl, start + th);
                        sig2.adopt_statics_from(p.sig, tl, th - tl, layout);
                        sig2.set_fold_edges(tl > 0, th < end - start, plan.dir, false);
                        sig2.concat_hmm_reject_strand(info, other);
                        sig2.set_hmm_note(plan.note);
                        p.sig = std::move(sig2);
                        p.rs = std::move(rs2);
                        p.trimmed = true;
                        (plan.dir == "forward" ? ctr.trim_keep_F : ctr.trim_keep_R).fetch_add(1, std::memory_order_relaxed);
                    } else {
                        if (!plan.dir.empty()) p.sig.concat_hmm_reject_strand(info, other);
                        if (plan.no_double_count) p.sig.mark_hmm_no_double_count();
                        if (plan.clip_check) p.sig.mark_hmm_clip_check(plan.clip_check);
                        p.sig.set_hmm_note(plan.note);
                    }
                    p.strand = plan.dir == "forward" ? 'F' : plan.dir == "reverse" ? 'R' : '?';
                    if (const int rej = p.sig.concat_hmm_partner_rule(info))
                        ctr.partner_rejected.fetch_add(static_cast<size_t>(rej), std::memory_order_relaxed);
                }
            }
            if (p.sig.weak_primer_count()) ctr.retry_from_partner.fetch_add(p.sig.weak_primer_count(), std::memory_order_relaxed);
            auto one_sided_edge = [&](int j) {
                if (j < 0 || j >= static_cast<int>(res.cut_flags.size())) return false;
                const uint8_t f = res.cut_flags[static_cast<size_t>(j)];
                return !(f & concat_hmm::BOUNDARY_UNRESOLVED) &&
                       static_cast<bool>(f & concat_hmm::BOUNDARY_LEFT_BOUNDED) !=
                       static_cast<bool>(f & concat_hmm::BOUNDARY_RIGHT_BOUNDED);
            };
            if (one_sided_edge(c - 1) || one_sided_edge(c))
                p.sig.set_hmm_boundary_note("layout_edge");
            pieces.push_back(std::move(p));
        }
        // At a weak junction keep only the piece with the better barcode unit (ties: HMM confidence, then forward).
        if (split && !pieces.empty()) {
            const float tau = layout.concat_model ? layout.concat_model->opt.tau_junction : 0.9f;
            bool any_weak = false;
            std::vector<uint8_t> weak(static_cast<size_t>(res.k - 1), 0);
            for (int j = 0; j + 1 < res.k; ++j) {
                const bool low_post = j < static_cast<int>(res.cut_post.size()) && res.cut_post[j] < tau;
                const bool midpoint = j < static_cast<int>(res.cut_kind.size()) && res.cut_kind[j] == concat_hmm::CUT_MIDPOINT;
                weak[j] = static_cast<uint8_t>(low_post || midpoint || (junction_abstain && res.p_single > 0.01f));
                any_weak = any_weak || weak[j];
            }
            if (any_weak) {
                std::vector<int> rank(pieces.size(), 0);
                for (size_t i = 0; i < pieces.size(); ++i) {
                    const std::string dir = pieces[i].strand == 'F' ? "forward" : pieces[i].strand == 'R' ? "reverse" : "";
                    rank[i] = pieces[i].sig.concat_hmm_unit_rank(info, dir);
                }
                auto find = [&](int c) -> int { for (size_t i = 0; i < pieces.size(); ++i) if (pieces[i].c == c) return static_cast<int>(i); return -1; };
                for (int j = 0; j + 1 < res.k; ++j) {
                    if (!weak[j]) continue;
                    ctr.junctions_abstained.fetch_add(1, std::memory_order_relaxed);
                    const int L = find(j), R = find(j + 1);
                    if (L >= 0 && R >= 0) {
                        int keep = -1;
                        if (rank[L] != rank[R]) keep = rank[L] > rank[R] ? L : R;
                        else if (rank[L] > 0) {
                            const float cl = res.segs[j].conf, cr = res.segs[j + 1].conf;
                            keep = cl != cr ? (cl > cr ? L : R) : (pieces[L].strand == 'F' ? L : pieces[R].strand == 'F' ? R : L);
                        }
                        if (keep != L) pieces[L].suppressed = true;
                        if (keep != R) pieces[R].suppressed = true;
                    } else if (L >= 0 || R >= 0) {
                        const int x = L >= 0 ? L : R;
                        if (rank[x] == 0) pieces[x].suppressed = true;
                    }
                }
                size_t nsup = 0;
                for (auto& p : pieces) nsup += p.suppressed ? 1 : 0;
                ctr.pieces_suppressed.fetch_add(nsup, std::memory_order_relaxed);
                if (nsup == pieces.size()) {
                    ctr.drop_unresolved.fetch_add(1, std::memory_order_relaxed);
                    return true;
                }
            }
        }
        // one F + one R construct: duplet names (<id>-F-CT / <id>-R-CT); otherwise /segN
        if (split && res.k == 2 && pieces.size() == 2 &&
            ((pieces[0].strand == 'F' && pieces[1].strand == 'R') || (pieces[0].strand == 'R' && pieces[1].strand == 'F'))) {
            for (auto& p : pieces) { p.sig.set_id(read.id); p.rs.id = read.id; }
            ctr.duplet_named.fetch_add(1, std::memory_order_relaxed);
        }
        for (auto& p : pieces) {
            if (p.suppressed) continue;
            process_molecule(p.sig, p.rs);
            const uint8_t fu = p.sig.fold_boundary_used();
            if (fu && p.sig.fold_edge_is_fold()) {
                if (fu & 1) ctr.fold_read_start.fetch_add(1, std::memory_order_relaxed);
                if (fu & 2) ctr.fold_read_end.fetch_add(1, std::memory_order_relaxed);
            }
        }
        return true;
    }

    static sigalign_run_stats sigalign(
        const std::string& fastq_path, 
        const ReadLayout& layout, 
        const std::string& output_prefix, 
        std::optional<int> gen_mut, 
        bool verbose, 
        int num_threads, 
        size_t chunk_size, 
        size_t max_reads, 
        bool write_debug,
        std::string mode,
        std::string joint_bc_mode = "default",
        bool rc_umi = true,
        int min_read_length = -1,
        std::function<void()> boundary_poll = {},
        FILE* stdout_sink = nullptr,
        std::function<void(const std::string&)> output_callback = {}
    ) {
    const auto sigalign_wall_t0 = std::chrono::steady_clock::now();

    std::string file_out = path_utils::get_fastqa_type(fastq_path);
    std::string fastq_output_path = output_prefix + file_out;
    bool compress_fastq = true;
    if (stdout_sink && output_callback)
        throw std::invalid_argument("Specify only one sequence output sink");
    const bool stream_output = stdout_sink || static_cast<bool>(output_callback);
    if (stream_output) write_debug = false;

    std::unique_ptr<std::ofstream> metrics_file_ptr;
    std::unique_ptr<sigstring_writing> debug_sig_writer, debug_csv_writer, debug_fastqa_writer;

    // Primary FASTQA writer
    std::unique_ptr<sigstring_writing> primary_writer;
    if (stdout_sink)
        primary_writer = std::make_unique<sigstring_writing>(stdout_sink, sigstring_writing::format::FASTQA);
    else if (!output_callback)
        primary_writer = std::make_unique<sigstring_writing>(fastq_output_path, sigstring_writing::format::FASTQA,
                                                          compress_fastq, false, num_threads);
    std::unique_ptr<sigstring_writing> unresolved_writer;
    const std::string unresolved_path = output_prefix + "_boundary_unresolved" + file_out;
    if (layout.concat_model && !stream_output)
        unresolved_writer = std::make_unique<sigstring_writing>(unresolved_path,
            sigstring_writing::format::FASTQA, compress_fastq, false, num_threads);

    if (write_debug) {
        std::string sig_path = output_prefix + "_dbg.sig";
        std::string csv_path = output_prefix + "_dbg.csv";
        std::string fastq_debug_path = output_prefix + "_dbg" + file_out;
        std::string metrics_path = output_prefix + ".metrics.tsv";

        debug_sig_writer = std::make_unique<sigstring_writing>(sig_path, sigstring_writing::format::SIGSTRING,
            compress_fastq, false, num_threads
        );

        debug_csv_writer = std::make_unique<sigstring_writing>(csv_path, sigstring_writing::format::CSV,
            compress_fastq, false, num_threads
        );

        //adding header to the debug csv
        {
            SigString header("", 0);
            (*debug_csv_writer)(std::vector<SigString>{header});
        }

        debug_fastqa_writer = std::make_unique<sigstring_writing>(fastq_debug_path, sigstring_writing::format::FASTQA,
            compress_fastq, false, num_threads
        );
        // Honor --no-umi-rc in the debug FASTQA output too (it flows through
        // to_fastqa(), not the buffered to_fastqa_append() path).
        debug_fastqa_writer->set_rc_umi(rc_umi);

        metrics_file_ptr = std::make_unique<std::ofstream>(metrics_path);
        if (!*metrics_file_ptr) {
            throw std::runtime_error("Failed to open metrics file: " +
                                     metrics_path);
        }
        *metrics_file_ptr << "chunk_id\tseqs_in_chunk\tseqs_passed\tin_flight\tprocess_time_ms\tqueue_time_ms\ttotal_time_ms\trss_mb\n";
        if (!*metrics_file_ptr) {
            throw std::runtime_error("Failed to write metrics header: " +
                                     metrics_path);
        }
    }

    parallel_writer writer;

    int pigz_threads = (num_threads > 0 ? num_threads : 1);
    if (const char* e = std::getenv("RAD_PIGZ_THREADS")) {
        int v = std::atoi(e);
        if (v > 0) pigz_threads = v;
    }

    std::atomic<size_t> total_reads{0};
    std::atomic<size_t> total_passed{0};
    std::atomic<size_t> total_demultiplexed{0};
    std::atomic<size_t> total_records_written{0};
    std::atomic<size_t> chunk_id_ctr{0};
    std::atomic<long long> total_process_time_ms{0};
    std::atomic<long long> total_queue_time_ms{0};
    std::mutex metrics_mu;

    barcode_position_stats::reset();

    concat_hmm_counters hmm_ctr;
    std::unique_ptr<concat_layout_info> hmm_info;
    if (layout.concat_model) {
        hmm_info = std::make_unique<concat_layout_info>(build_concat_layout_info(layout, *layout.concat_model));
        hmm_info->ctr = &hmm_ctr;
    }

    auto process_chunk = [&](std::vector<read_streaming::sequence>& chunk, 
                             const std::string& path)
    {
        if (chunk.empty()){
            return;
        }
        writer.rethrow_if_failed();

        const size_t my_chunk_id = ++chunk_id_ctr;
        const auto wall_t0 = std::chrono::steady_clock::now();

        // Thread-local buffers keep alignment/serialization parallel.  Each
        // span records which input read contributed a contiguous slice so the
        // serial consolidation step can restore FASTQ input order regardless
        // of OpenMP's dynamic worker assignment.
        struct serialized_read_span {
            size_t input_index;
            size_t worker_index;
            size_t offset;
            size_t length;
        };
        const int worker_count = std::max(1, num_threads);
        std::vector<std::string> thread_buffers(worker_count);
        std::vector<std::vector<serialized_read_span>> thread_spans(
            worker_count);
        std::vector<std::string> unresolved_buffers(worker_count);
        std::vector<std::vector<serialized_read_span>> unresolved_spans(worker_count);
        std::atomic<size_t> passed_count{0};
        std::atomic<size_t> demultiplexed_count{0};
        std::atomic<size_t> records_written_count{0};

        std::vector<std::vector<SigString>> debug_by_input;
        if (write_debug) {
            debug_by_input.resize(chunk.size());
        }
        const size_t est_per_thread =
            (chunk.size() / static_cast<size_t>(worker_count) + 1) * 1600;
        for (int i = 0; i < worker_count; ++i) {
            thread_buffers[i].reserve(est_per_thread);
            thread_spans[i].reserve(
                chunk.size() / static_cast<size_t>(worker_count) + 1);
        }

        // Track per-thread processing time
        std::vector<double> thread_times(worker_count, 0.0);
        std::atomic<bool> processing_failed{false};
        std::exception_ptr processing_error;
        std::mutex processing_error_mutex;

        // === process and serialize in one pass ===
        #pragma omp parallel num_threads(worker_count)
        {
            int tid = omp_get_thread_num();
            auto thread_start = std::chrono::steady_clock::now();

            #pragma omp for schedule(dynamic) nowait
            for (size_t i = 0; i < chunk.size(); ++i) {
                if (processing_failed.load(std::memory_order_acquire) || writer.failed()) continue;
                try {
                    const size_t offset = thread_buffers[tid].size();
                    const size_t unresolved_offset = unresolved_buffers[tid].size();
                    auto process_one = [&]() {
                const auto& read = chunk[i];
#ifdef RAD_STAGE_TIMERS
                rad_stage_timers::scope t_read(rad_stage_timers::read_ns);
                rad_stage_timers::reads.fetch_add(1, std::memory_order_relaxed);
#endif

                // --concat-hmm: reads the router does not handle fall through to the default path
                if (hmm_info) {
                    bool hmm_passed = false;
                    size_t hmm_records = 0;
                    auto hmm_molecule = [&](SigString& molecule, const read_streaming::sequence& molecule_read) {
#ifdef RAD_STAGE_TIMERS
                        rad_stage_timers::molecules.fetch_add(1, std::memory_order_relaxed);
                        { rad_stage_timers::scope t_var(rad_stage_timers::variable_ns);
#endif
                        molecule.sigalign_variable(molecule_read, layout, verbose);
                        molecule.concat_hmm_apply_barcode_rejections();
#ifdef RAD_STAGE_TIMERS
                        }
                        rad_stage_timers::scope t_filt(rad_stage_timers::filter_ns);
                        const auto t_fb0 = std::chrono::steady_clock::now();
#endif
                        molecule.sigalign_filter(molecule_read, layout, gen_mut.value_or(2), verbose,
                                                 mode, joint_bc_mode, min_read_length);
#ifdef RAD_STAGE_TIMERS
                        {
                            const int oc = (molecule.type() != "filtered" && molecule.type() != "skipped") ? 0
                                         : molecule.info().find("BARCODE_CORRECTION_FAILED") != std::string::npos ? 1 : 2;
                            rad_stage_timers::fb_ns[0][oc].fetch_add(static_cast<unsigned long long>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                std::chrono::steady_clock::now() - t_fb0).count()), std::memory_order_relaxed);
                            rad_stage_timers::fb_n[0][oc].fetch_add(1, std::memory_order_relaxed);
                        }
#endif
                        // two barcode units that both pass would be written twice: never double-count
                        if (molecule.hmm_no_double_count_set() && molecule.read_type == "concatenate") {
                            molecule.set_type("skipped");
                            hmm_ctr.both_pass_dropped.fetch_add(1, std::memory_order_relaxed);
                        } else if (molecule.hmm_clip_hit_in_written_read()) {
                            molecule.set_type("skipped");
                            molecule.set_hmm_note(molecule.hmm_note_str() + "_CLIP");
                            hmm_ctr.two_unit_clip_dropped.fetch_add(1, std::memory_order_relaxed);
                        }
                        if (write_debug) debug_by_input[i].push_back(molecule);
                        if (molecule.read_type != "filtered" && molecule.read_type != "skipped") {
                            hmm_passed = true;
                            hmm_records += molecule.to_fastqa_append(thread_buffers[tid], rc_umi);
                        }
                    };
                    auto preserve_unresolved = [&](const read_streaming::sequence& parent) {
                        if (!unresolved_writer) return;
                        std::string& b = unresolved_buffers[static_cast<size_t>(tid)];
                        b += parent.is_fastq ? '@' : '>';
                        b += parent.id;
                        if (!parent.comment.empty()) { b += ' '; b += parent.comment; }
                        b += "\tHU:Z:boundary_unresolved\n";
                        b += parent.seq;
                        b += '\n';
                        if (parent.is_fastq) { b += '+'; b += parent.plus; b += '\n'; b += parent.qual; b += '\n'; }
                    };
                    if (concat_hmm_route(read, layout, *hmm_info, hmm_molecule, verbose, preserve_unresolved)) {
                        if (hmm_passed) passed_count.fetch_add(1, std::memory_order_relaxed);
                        if (hmm_records > 0) {
                            demultiplexed_count.fetch_add(1, std::memory_order_relaxed);
                            records_written_count.fetch_add(hmm_records, std::memory_order_relaxed);
                        }
                        return;
                    }
                }
                
                // process the read
                {
                    SigString sig(read.id, read.seq.length(),"undefined", layout.sequencing_type);
                    std::vector<seq_element> additional_hits;
#ifdef RAD_STAGE_TIMERS
                    std::optional<rad_stage_timers::scope> t_static;
                    t_static.emplace(rad_stage_timers::legacy_static_ns);
#endif
                    sig.sigalign_static(read, layout, verbose, &additional_hits);
                    auto segments = sig.group_static_candidates(layout, additional_hits, verbose, min_read_length);
#ifdef RAD_STAGE_TIMERS
                    t_static.reset();
#endif
                    bool parent_passed = false;
                    size_t parent_records = 0;
                    auto process_molecule = [&](SigString& molecule, const read_streaming::sequence& molecule_read) {
#ifdef RAD_STAGE_TIMERS
                        rad_stage_timers::molecules.fetch_add(1, std::memory_order_relaxed);
                        { rad_stage_timers::scope t_var(rad_stage_timers::variable_ns);
#endif
                        molecule.sigalign_variable(molecule_read, layout, verbose);
#ifdef RAD_STAGE_TIMERS
                        }
                        rad_stage_timers::scope t_filt(rad_stage_timers::filter_ns);
                        const auto t_fb0 = std::chrono::steady_clock::now();
#endif
                        molecule.sigalign_filter(molecule_read, layout, gen_mut.value_or(2), verbose,
                                                 mode, joint_bc_mode, min_read_length);
#ifdef RAD_STAGE_TIMERS
                        {
                            const int oc = (molecule.type() != "filtered" && molecule.type() != "skipped") ? 0
                                         : molecule.info().find("BARCODE_CORRECTION_FAILED") != std::string::npos ? 1 : 2;
                            rad_stage_timers::fb_ns[1][oc].fetch_add(static_cast<unsigned long long>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                std::chrono::steady_clock::now() - t_fb0).count()), std::memory_order_relaxed);
                            rad_stage_timers::fb_n[1][oc].fetch_add(1, std::memory_order_relaxed);
                        }
#endif
                        if (write_debug) debug_by_input[i].push_back(molecule);
                        if (molecule.read_type != "filtered" && molecule.read_type != "skipped") {
                            parent_passed = true;
                            parent_records += molecule.to_fastqa_append(thread_buffers[tid], rc_umi);
                        }
                    };
                    if (segments.empty()) {
                        process_molecule(sig, read);
                    } else {
                        for (size_t index = 0; index < segments.size(); ++index) {
                            const auto& segment = segments[index];
                            const int offset = segment.position.first - 1;
                            const int length = segment.position.second - offset;
                            read_streaming::sequence molecule_read{
                                read.id + "/seg" + std::to_string(index + 1), read.comment,
                                read.seq.substr(offset, length), read.plus,
                                read.is_fastq ? read.qual.substr(offset, length) : "", read.is_fastq};
                            SigString molecule(molecule_read.id, length, "undefined", layout.sequencing_type);
                            if (!molecule.assign_static_segment(segment, layout, molecule_read)) {
                                if (verbose) molecule.log_verbose("CONCAT_UNRESOLVED masked-capture traceback");
                                continue;
                            }
                            process_molecule(molecule, molecule_read);
                        }
                    }
                    // Success remains parent-scoped; emitted records count
                    // molecules. Filtering/counter side effects occur once.
                    if (parent_passed) passed_count.fetch_add(1, std::memory_order_relaxed);
                    if (parent_records > 0) {
                        demultiplexed_count.fetch_add(1, std::memory_order_relaxed);
                        records_written_count.fetch_add(parent_records, std::memory_order_relaxed);
                    }
                }  // sig destroyed here
                    };
                    process_one();
                    if (thread_buffers[tid].size() > offset)
                        thread_spans[tid].push_back({i, static_cast<size_t>(tid), offset,
                                                    thread_buffers[tid].size() - offset});
                    if (unresolved_buffers[tid].size() > unresolved_offset)
                        unresolved_spans[tid].push_back({i, static_cast<size_t>(tid), unresolved_offset,
                                                        unresolved_buffers[tid].size() - unresolved_offset});
                } catch (...) {
                    std::lock_guard<std::mutex> lock(processing_error_mutex);
                    if (!processing_error) processing_error = std::current_exception();
                    processing_failed.store(true, std::memory_order_release);
                }
            }
            
            auto thread_end = std::chrono::steady_clock::now();
            thread_times[tid] = std::chrono::duration_cast<std::chrono::milliseconds>(
                thread_end - thread_start).count();

        }

        // No C++ exception may cross an OpenMP structured block. Rethrow the
        // first worker failure now, on the serial caller, before staging any
        // partial output or updating run counters.
        if (processing_error) std::rethrow_exception(processing_error);

        std::vector<SigString> debug_sigs;
        if (write_debug) {
            debug_sigs.reserve(chunk.size());
            for (auto& sig : debug_by_input) {
                for (auto& molecule : sig) debug_sigs.push_back(std::move(molecule));
            }
        }

        writer.rethrow_if_failed();
        const auto wall_t1 = std::chrono::steady_clock::now();
        
        // Calculate actual work time (max across all threads)
        double actual_process_ms = 0.0;
        for (double t : thread_times) {
            actual_process_ms = std::max(actual_process_ms, t);
        }

        // Consolidate complete input-parent spans on the serial caller. R
        // callbacks must never be called by alignment or writer workers.
        auto consolidate = [](const auto& spans, const auto& buffers) {
            std::vector<serialized_read_span> ordered;
            size_t bytes = 0;
            for (const auto& local : spans) for (const auto& span : local) {
                ordered.push_back(span);
                bytes += span.length;
            }
            std::sort(ordered.begin(), ordered.end(), [](const auto& a, const auto& b) {
                return a.input_index < b.input_index;
            });
            std::string result;
            result.reserve(bytes);
            for (const auto& span : ordered)
                result.append(buffers[span.worker_index], span.offset, span.length);
            return result;
        };
        if (unresolved_writer)
            writer.write_raw_string(*unresolved_writer, consolidate(unresolved_spans, unresolved_buffers));
        if (passed_count > 0) {
            auto output = consolidate(thread_spans, thread_buffers);
            if (output_callback) {
                if (!output.empty()) output_callback(output);
            } else {
                writer.write_raw_string(*primary_writer, std::move(output));
            }
        }
        if (write_debug && !debug_sigs.empty()) {
            writer.write_debug(debug_sig_writer.get(), debug_csv_writer.get(),
                               debug_fastqa_writer.get(), debug_sigs);
        }

        const auto wall_t2 = std::chrono::steady_clock::now();

        // Update counters
        total_reads += chunk.size();
        total_passed += passed_count.load();
        total_demultiplexed += demultiplexed_count.load();
        total_records_written += records_written_count.load();

        const double wall_process_ms = std::chrono::duration_cast<std::chrono::milliseconds>(wall_t1 - wall_t0).count();
        const double queue_ms = std::chrono::duration_cast<std::chrono::milliseconds>(wall_t2 - wall_t1).count();
        
        // Use actual thread time for accounting
        total_process_time_ms.fetch_add(static_cast<long long>(actual_process_ms), std::memory_order_relaxed);
        total_queue_time_ms.fetch_add(static_cast<long long>(queue_ms), std::memory_order_relaxed);

        // Metrics
        if (write_debug && metrics_file_ptr) {
            std::lock_guard<std::mutex> lk(metrics_mu);
            const size_t in_flight =
                readmem_utils::inflight().chunks.load(std::memory_order_relaxed);
            const double rss_gib = readmem_utils::rss_gib();
            const double rss_mb = rss_gib < 0.0 ? -1.0 : rss_gib * 1024.0;
            *metrics_file_ptr << my_chunk_id << "\t" 
                              << chunk.size() << "\t" 
                              << passed_count.load() << "\t"
                              << in_flight << "\t"
                              << actual_process_ms << "\t" 
                              << queue_ms << "\t" 
                              << (actual_process_ms + queue_ms) << "\t"
                              << rss_mb << "\n";
            if (!*metrics_file_ptr) {
                throw std::runtime_error("Failed to write debug metrics");
            }
        }

        //print chunk stats
        {
            std::lock_guard<std::mutex> lk(metrics_mu);
            RAD_COUT << "[chunk_stats] " << my_chunk_id 
                      << ", processed=" << chunk.size()
                      << ", passed=" << passed_count << " (" 
                      << (chunk.size() ? (double)passed_count.load() * 100.0 / (double)chunk.size() : 0.0) 
                      << "%), wall=" << wall_process_ms / 1000.0 << "s"
                      << ", actual=" << actual_process_ms / 1000.0 << "s"
                      << ", queue=" << queue_ms / 1000.0 << "s\n";
            memory_utils::get_rss();

            if (my_chunk_id % 100 == 0) {
                for (const auto& kv : layout.wl_map.maps) {
                    auto& wl = kv.second.get();
                    bc_mem_utils::print_mem_snapshot(wl.true_bcs, wl.global_bcs, my_chunk_id);
                }
            }
        }

        // Cleanup
        std::vector<read_streaming::sequence>().swap(chunk);
        std::vector<std::string>().swap(thread_buffers);
        if (write_debug) {
            std::vector<SigString>().swap(debug_sigs);
        }
    };

    // Stream one chunk at a time and parallelize within that chunk.  Running
    // chunk_streaming with num_threads here created an outer OpenMP team which
    // then entered the num_threads-wide region in process_chunk above.  With
    // nested OpenMP enabled that could create num_threads^2 workers; with the
    // usual nested-disabled runtime it instead kept num_threads full chunks in
    // memory while each inner region ran serially.
    {
        chunk_streaming<read_streaming::sequence, decltype(process_chunk)>streamer(chunk_size, pigz_threads);
        const int64_t limit = (max_reads > 0 ? static_cast<int64_t>(max_reads) : -1);
        streamer.process_chunks(fastq_path, process_chunk,
                                /*chunk_workers=*/1, limit, false, 0,
                                boundary_poll);
    }

    writer.stop();
    if (primary_writer) primary_writer->finalize();
    if (unresolved_writer) unresolved_writer->finish();

    if (debug_sig_writer) debug_sig_writer->finalize();
    if (debug_csv_writer) debug_csv_writer->finalize();
    if (debug_fastqa_writer) debug_fastqa_writer->finalize();

    if (write_debug && metrics_file_ptr) {
        metrics_file_ptr->flush();
        if (!*metrics_file_ptr) {
            throw std::runtime_error("Failed to flush debug metrics");
        }
        metrics_file_ptr->close();
        if (metrics_file_ptr->fail()) {
            throw std::runtime_error("Failed to close debug metrics");
        }
    }

    // Summary
    const auto sigalign_wall_t1 = std::chrono::steady_clock::now();
    const double wall_total_s = std::chrono::duration_cast<std::chrono::milliseconds>(
        sigalign_wall_t1 - sigalign_wall_t0
    ).count() / 1000.0;

    const double process_s = total_process_time_ms.load() / 1000.0;
    const double queue_s   = total_queue_time_ms.load()   / 1000.0;
    const double accounted_s = process_s + queue_s;
    const double overhead_s = std::max(0.0, wall_total_s - accounted_s);

    RAD_COUT << "\n[sigalign] Performance Summary:\n";
    RAD_COUT << "────────────────────────────────────────────────────\n";
    RAD_COUT << "[sigalign] Wall-clock runtime: " << wall_total_s << " seconds\n";
    RAD_COUT << "[sigalign] Total chunks processed: " << chunk_id_ctr.load() << "\n";
    RAD_COUT << "[sigalign] Total reads processed: " << total_reads.load() << "\n";
    RAD_COUT << "[sigalign] Reads passing filter: " << total_passed.load() << " ("
              << (total_reads > 0 ? (total_passed.load() * 100.0 / (double)total_reads.load()) : 0.0) << "%)\n";
    RAD_COUT << "[sigalign] Reads demultiplexed: "
              << total_demultiplexed.load() << " ("
              << (total_reads > 0
                      ? (total_demultiplexed.load() * 100.0 /
                         static_cast<double>(total_reads.load()))
                      : 0.0)
              << "%)\n";
    RAD_COUT << "[sigalign] Output records serialized: "
              << total_records_written.load() << "\n";
    if (hmm_info) {
        const size_t legacy = hmm_ctr.legacy_guard + hmm_ctr.legacy_k0 + hmm_ctr.legacy_abstain +
                              hmm_ctr.legacy_artifact + hmm_ctr.legacy_too_many;
        RAD_COUT << "[sigalign] concat-HMM: reads=" << hmm_ctr.reads.load()
                  << " single=" << hmm_ctr.single.load() << " split=" << hmm_ctr.split.load()
                  << " (children=" << hmm_ctr.children.load()
                  << ", full-search children=" << hmm_ctr.children_full_search.load() << ") existing_path=" << legacy
                  << " [abstain=" << hmm_ctr.legacy_abstain.load() << " k0=" << hmm_ctr.legacy_k0.load()
                  << " artifact=" << hmm_ctr.legacy_artifact.load() << " guard=" << hmm_ctr.legacy_guard.load()
                  << " too_many=" << hmm_ctr.legacy_too_many.load() << "]"
                  << " dropped=" << (hmm_ctr.drop_abstain.load() + hmm_ctr.drop_too_many.load() + hmm_ctr.drop_d_reads.load())
                  << " [abstain=" << hmm_ctr.drop_abstain.load() << " too_many=" << hmm_ctr.drop_too_many.load()
                  << " D=" << hmm_ctr.drop_d_reads.load() << " unresolved=" << hmm_ctr.drop_unresolved.load()
                  << "] dropped_D_children=" << hmm_ctr.drop_d_children.load()
                  << " junctions_abstained=" << hmm_ctr.junctions_abstained.load()
                  << " trim_keep=" << hmm_ctr.trim_keep_F.load() + hmm_ctr.trim_keep_R.load()
                  << " duplets=" << hmm_ctr.duplet_named.load()
                  << " boundary_unresolved_saved="
                  << (unresolved_writer ? hmm_ctr.boundary_unresolved_parents.load() : 0)
                  << " abstain_policy=" << (layout.concat_abstain_legacy ? "legacy" : "drop") << "\n";
        if (stream_output)
            RAD_COUT << "[sigalign] unresolved boundaries flagged="
                      << hmm_ctr.boundary_unresolved_parents.load() << "; raw-parent output suppressed\n";
    }
#ifdef RAD_STAGE_TIMERS
    {
        const double n = std::max<double>(1.0, static_cast<double>(rad_stage_timers::reads.load()));
        RAD_COUT << "[stage_timers] reads=" << rad_stage_timers::reads.load()
                  << " read_us=" << rad_stage_timers::read_ns.load() / n / 1000.0
                  << " hmm_us=" << rad_stage_timers::hmm_ns.load() / n / 1000.0
                  << " windowed_static_us=" << rad_stage_timers::windowed_static_ns.load() / n / 1000.0
                  << " legacy_static_us=" << rad_stage_timers::legacy_static_ns.load() / n / 1000.0
                  << " static_stage_us="
                  << (rad_stage_timers::hmm_ns.load() + rad_stage_timers::windowed_static_ns.load() +
                      rad_stage_timers::legacy_static_ns.load()) / n / 1000.0
                  << " variable_us=" << rad_stage_timers::variable_ns.load() / n / 1000.0
                  << " filter_us=" << rad_stage_timers::filter_ns.load() / n / 1000.0
                  << " molecules=" << rad_stage_timers::molecules.load()
                  << " filter_split_us(path:pass/bcfail/other n)="
                  << rad_stage_timers::fb_ns[0][0].load() / n / 1000.0 << "/" << rad_stage_timers::fb_ns[0][1].load() / n / 1000.0 << "/"
                  << rad_stage_timers::fb_ns[0][2].load() / n / 1000.0 << " (" << rad_stage_timers::fb_n[0][0].load() << "/"
                  << rad_stage_timers::fb_n[0][1].load() << "/" << rad_stage_timers::fb_n[0][2].load() << ") legacy "
                  << rad_stage_timers::fb_ns[1][0].load() / n / 1000.0 << "/" << rad_stage_timers::fb_ns[1][1].load() / n / 1000.0 << "/"
                  << rad_stage_timers::fb_ns[1][2].load() / n / 1000.0 << " (" << rad_stage_timers::fb_n[1][0].load() << "/"
                  << rad_stage_timers::fb_n[1][1].load() << "/" << rad_stage_timers::fb_n[1][2].load() << ")"
                  << "\n";
    }
#endif
    RAD_COUT << "[sigalign] Timing breakdown:\n";
    RAD_COUT << "  - Accounted chunk runtime: " << accounted_s << " seconds\n";
    RAD_COUT << "  - Process time: " << process_s << " seconds\n";
    RAD_COUT << "  - Output staging time: " << queue_s << " seconds\n";
    RAD_COUT << "  - Overhead (writer drain/setup): " << overhead_s << " seconds\n";

    if (stream_output) {
        RAD_COUT << "\n[sigalign][fastq]: stdout (uncompressed)\n";
    } else if (write_debug) {
        RAD_COUT << "\n[sigalign] Output written to:\n"
                  << "[sigalign] [sigstring]: " << output_prefix << ".sig\n"
                  << "[sigalign]       [csv]: " << output_prefix << ".csv\n"
                  << "[sigalign]     [fastq]: " << fastq_output_path << (compress_fastq ? ".gz" : "") << "\n"
                  << "[sigalign]   [metrics]: " << output_prefix << ".metrics.tsv\n";
    } else {
        RAD_COUT << "\n[sigalign] Output written to:\n"
                  << "[sigalign][fastq]: " << fastq_output_path << (compress_fastq ? ".gz" : "") << "\n";
    }
    if (unresolved_writer)
        RAD_COUT << "[sigalign][unresolved boundaries]: " << unresolved_path << ".gz ("
                  << hmm_ctr.boundary_unresolved_parents.load() << " raw parent copies; unresolved intervals unassigned)\n";

    sigalign_run_stats run_stats{
        total_reads.load(),
        total_passed.load(),
        total_demultiplexed.load(),
        total_records_written.load(),
        chunk_id_ctr.load(),
        wall_total_s,
        process_s,
        queue_s,
        overhead_s
    };
    if (hmm_info) {
        run_stats.concat_hmm_enabled = true;
        run_stats.hmm_reads = hmm_ctr.reads;
        run_stats.hmm_single = hmm_ctr.single;
        run_stats.hmm_split = hmm_ctr.split;
        run_stats.hmm_children = hmm_ctr.children;
        run_stats.hmm_children_full_search = hmm_ctr.children_full_search;
        run_stats.hmm_same_molecule = hmm_ctr.same_molecule;
        run_stats.hmm_legacy_guard = hmm_ctr.legacy_guard;
        run_stats.hmm_legacy_k0 = hmm_ctr.legacy_k0;
        run_stats.hmm_legacy_abstain = hmm_ctr.legacy_abstain;
        run_stats.hmm_legacy_artifact = hmm_ctr.legacy_artifact;
        run_stats.hmm_legacy_too_many = hmm_ctr.legacy_too_many;
        run_stats.hmm_full_as_missing = hmm_ctr.full_as_missing;
        run_stats.hmm_spacer_retry = hmm_ctr.spacer_retry;
        run_stats.hmm_spacer_retry_ok = hmm_ctr.spacer_retry_ok;
        run_stats.hmm_read_end_ok = hmm_ctr.read_end_ok;
        run_stats.hmm_fragment_ok = hmm_ctr.fragment_ok;
        run_stats.hmm_abstain_legacy = layout.concat_abstain_legacy;
        run_stats.hmm_drop_abstain = hmm_ctr.drop_abstain;
        run_stats.hmm_drop_too_many = hmm_ctr.drop_too_many;
        run_stats.hmm_drop_d_reads = hmm_ctr.drop_d_reads;
        run_stats.hmm_drop_d_children = hmm_ctr.drop_d_children;
        run_stats.hmm_d_kept_left = hmm_ctr.d_kept_left;
        run_stats.hmm_d_kept_right = hmm_ctr.d_kept_right;
        run_stats.hmm_d_split = hmm_ctr.d_split;
        run_stats.hmm_fold_read_start = hmm_ctr.fold_read_start;
        run_stats.hmm_fold_read_end = hmm_ctr.fold_read_end;
        run_stats.hmm_partner_rejected = hmm_ctr.partner_rejected;
        run_stats.hmm_retry_from_partner = hmm_ctr.retry_from_partner;
        run_stats.hmm_art_reads = hmm_ctr.art_reads;
        run_stats.hmm_pieces = hmm_ctr.pieces;
        run_stats.hmm_trim_keep_F = hmm_ctr.trim_keep_F;
        run_stats.hmm_trim_keep_R = hmm_ctr.trim_keep_R;
        run_stats.hmm_both_pass_dropped = hmm_ctr.both_pass_dropped;
        run_stats.hmm_clip_hits = hmm_ctr.clip_hits;
        run_stats.hmm_two_unit_clip_dropped = hmm_ctr.two_unit_clip_dropped;
        run_stats.hmm_duplet_named = hmm_ctr.duplet_named;
        run_stats.hmm_junctions_abstained = hmm_ctr.junctions_abstained;
        run_stats.hmm_pieces_suppressed = hmm_ctr.pieces_suppressed;
        run_stats.hmm_drop_unresolved = hmm_ctr.drop_unresolved;
        run_stats.hmm_boundary_unresolved_parents = hmm_ctr.boundary_unresolved_parents;
        run_stats.hmm_boundary_unresolved_children = hmm_ctr.boundary_unresolved_children;
        for (int kd = 0; kd < 6; ++kd) run_stats.hmm_cut_kind[kd] = hmm_ctr.cut_kind[kd];
        for (int s = 0; s < 4; ++s) {
            run_stats.hmm_win[s] = hmm_ctr.win[s];
            run_stats.hmm_win_ok[s] = hmm_ctr.win_ok[s];
        }
    }
    run_stats.bc_pos_moved = barcode_position_stats::moved;
    run_stats.bc_pos_moved_forward = barcode_position_stats::moved_forward;
    run_stats.bc_pos_moved_reverse = barcode_position_stats::moved_reverse;
    run_stats.bc_pos_in_place = barcode_position_stats::in_place;
    run_stats.bc_pos_window = barcode_position_stats::window;
    run_stats.bc_pos_unknown = barcode_position_stats::unknown;
    run_stats.bc_pos_tie = barcode_position_stats::tie;
    run_stats.bc_pos_clipped = barcode_position_stats::clipped;
    run_stats.bc_pos_out_of_read = barcode_position_stats::out_of_read;
    return run_stats;
}

/**
 * @brief convert SigString to a sigstring text representation
 * @return `std::string` sigstring representation
 */
    std::string to_sigstring() const {
        std::stringstream ss;
        std::string overall = read_type; // e.g., "forward", "reverse", "concatenate", "filtered"

        // Group elements by direction
        std::map<std::string, std::vector<std::reference_wrapper<const seq_element>>> direction_elements;
        for (const auto& elem : sig_elements) {
            direction_elements[elem.direction].push_back(std::ref(elem));
        }
        
        bool first_direction = true;
        // Process each directional group
        for (const auto& [direction, elements] : direction_elements) {
            // If overall read type restricts the output, skip the other direction
            if ((overall == "forward" && direction != "forward") || (overall == "reverse" && direction != "reverse")) {
                 continue;
            }
            if (elements.empty()) {
                // Skip empty groups
                continue;
            }
            if (!first_direction) {
                ss << "\n";
            }
            first_direction = false;
            
            // Sort elements by order
            std::vector<std::reference_wrapper<const seq_element>> sorted_elements = elements;
            std::sort(sorted_elements.begin(), sorted_elements.end(),
                    [](const auto& a, const auto& b) {
                        return a.get().position.first < b.get().position.first;
                    });
            
            bool first_elem = true;
            // --concat-hmm piece: print in parent coordinates (debug output only)
            const bool parent_frame = hmm_parent_off >= 0;
            for (const auto& elem_ref : sorted_elements) {
                const auto& elem = elem_ref.get();

                if (!first_elem) {
                    ss << "|";
                }
                if(elem.position.first == -1 && elem.position.second == -1){
                    continue; // Skip elements with invalid positions
                }
                first_elem = false;
                if (parent_frame && (elem.global_class == "start" || elem.global_class == "stop")) {
                    const bool at_start = elem.global_class == "start";
                    const int p = at_start ? hmm_parent_off : hmm_parent_end;
                    ss << (elem.direction == "reverse" ? "rc_" : "") << (at_start ? "seg_start" : "seg_stop") << ":0:" << p << ":" << p;
                    continue;
                }
                const int shift = parent_frame ? hmm_parent_off : 0;
                ss << elem.class_id << ":"
                << (elem.edit_distance ? std::to_string(elem.edit_distance.value()) : "0") << ":"
                << elem.position.first + shift << ":"
                << elem.position.second + shift;
            }
            // Append final tag: direction abbreviated as F or R
            std::string dir = (direction == "forward") ? "F" : "R";
            std::string concat = (from_concatemer || read_type == "concatenate") ? ":C" : "";
            std::string combined_info = "";
            if (read_type == "filtered") {
                combined_info += ":filtered";
            }
            if (additional_info.length() > 2) {
                combined_info += ":" + additional_info;
            }
            if (!hmm_note.empty()) {
                combined_info += ":HMM=" + hmm_note;
            }
            ss << "<" << sequence_length << ":" << sequence_id << ":" << dir << concat <<  combined_info << ">";
        }
        // Return empty string if nothing was added
        return ss.str().empty() ? "" : ss.str();
    }

/**
 * @brief append fastqa representation of SigString to a buffer
 * @param buffer `std::string&` buffer to be appended to
 * @return number of FASTQ/FASTA records appended
 * 
 * @brief this function generates and writes the fastqa version of a sigstring. 
 * It handles different read types including "forward", "reverse", "concatenate", and "skipped".
 * For "concatenate", it processes both forward and reverse directions.
 * It collects barcode, UMI, and read sequences from the sig_elements, constructs the appropriate tags,
 * and appends the formatted output directly to the provided buffer. It's "to_fastqa_append" because
 * it appends directly to an existing string buffer rather than returning a new string, which I don't have to
 * go through the overhead of creating intermediate strings.
 */
    size_t to_fastqa_append(std::string& buffer, bool rc_umi = true) const {
        if (read_type == "skipped") {
            return 0;
        }
        size_t records_written = 0;
        std::vector<std::string> dirs;
        if (read_type == "concatenate") {
            dirs = {
                "forward", 
                "reverse"
            };
        } else {
            dirs = {
                read_type
            };
        }
        
        for (const auto& dir : dirs) {
            std::vector<std::string> bc_keys;
            std::unordered_map<std::string, std::string> bc_map;
            std::unordered_map<std::string, std::string> bc_dir;
            std::unordered_map<std::string, std::string> cr_map;
            std::string umi, read_seq, read_qual;
            
            // Collect elements
            for (const auto& elem : sig_elements) {
                if (!elem.seq.has_value()) {
                    continue;
                }
                if (elem.direction != dir) {
                    continue;
                }

                if (elem.global_class == "barcode") {
                    if (elem.write.has_value() && !elem.write.value()) {
                        continue; // Skip barcodes explicitly marked not to write
                    }
                    auto key = seq_utils::remove_rc(elem.class_id);
                    bool is_fwd = (dir == "forward");
                    if (bc_map.find(key) == bc_map.end()) {
                        bc_keys.push_back(key);
                        bc_map[key] = elem.seq.value();
                        bc_dir[key] = dir;
                        if (elem.original_seq.has_value()) {
                            cr_map[key] = elem.original_seq.value();
                        }
                    } else if (is_fwd && bc_dir[key] == "reverse") {
                        bc_map[key] = elem.seq.value();
                        bc_dir[key] = dir;
                        if (elem.original_seq.has_value()) {
                            cr_map[key] = elem.original_seq.value();
                        }
                    }
                    continue;
                }
                
                if (elem.global_class == "umi") {
                    if (elem.seq.has_value()) {
                        // Reverse reads are extracted on the minus strand, so the raw
                        // UMI comes out reverse-complemented relative to the molecule.
                        // Barcodes are flipped back to plus-strand during correction
                        // (correct_barcode revcomps before whitelist match), so unless
                        // we mirror that here CB:Z ends up plus-strand while UB:Z stays
                        // minus-strand. rc_umi (default on) keeps the two consistent;
                        // pass rc_umi=false to leave UMIs exactly as extracted.
                        umi = (rc_umi && dir == "reverse")
                                  ? seq_utils::revcomp(elem.seq.value())
                                  : elem.seq.value();
                    }
                    continue;
                }

                if (elem.global_class == "read") {
                    read_seq = elem.seq.value();
                    if (elem.qual.has_value()) {
                        read_qual = elem.qual.value();
                    }
                    continue;
                }
            }

            // In bulk mode with no barcodes, use N-masked regions from static elements as barcodes
            if (bc_map.empty() && additional_info == "bulk") {
                for (const auto& elem : sig_elements) {
                    if (elem.type != "static" || elem.direction != dir) {
                        continue;
                    }
                    // Static elements with original_seq have N-masked regions extracted
                    if (elem.original_seq.has_value() && elem.seq.has_value()) {
                        auto key = seq_utils::remove_rc(elem.class_id);
                        if (bc_map.find(key) == bc_map.end()) {
                            bc_keys.push_back(key);
                            bc_map[key] = elem.seq.value();  // N-masked region(s)
                            bc_dir[key] = dir;
                            cr_map[key] = elem.original_seq.value();  // Full aligned sequence
                        }
                    }
                }
            }

            if ((bc_map.empty() && additional_info != "bulk") || read_seq.empty()) {
                continue;
            }
            
            // Build barcode tags
            std::string cb_tag, cr_tag;
            for (size_t i = 0; i < bc_keys.size(); ++i) {
                const auto& key = bc_keys[i];
                if (i) {
                    cb_tag += '-';
                    if (!cr_tag.empty()) cr_tag += '-';
                }
                cb_tag += bc_map[key];
                if (cr_map.find(key) != cr_map.end()) {
                    cr_tag += cr_map[key];
                }
            }
            
            bool is_fastq = !read_qual.empty();
            bool is_concatenate = from_concatemer || (read_type == "concatenate");
            bool is_forward = (dir == "forward");
            
            // Append directly to buffer - no intermediate string
            buffer += (is_fastq ? '@' : '>');
            buffer += sequence_id;
            buffer += (is_forward ? "-F" : "-R");
            if (is_concatenate) buffer += "-CT";
            
            if (!cb_tag.empty()) {
                buffer += "\tCB:Z:";
                buffer += cb_tag;
            }
            
            if (!cr_tag.empty() && (cb_tag != cr_tag && seq_utils::revcomp(cb_tag) != cr_tag)) {
                buffer += "\tCR:Z:";
                buffer += cr_tag;
            } else {
                buffer += "\tCR:Z:";
            }
            
            if (!umi.empty()) {
                buffer += "\tUB:Z:";
                buffer += umi;
            }
            
            if (!hmm_boundary_note.empty()) {
                buffer += "\tHB:Z:";
                buffer += hmm_boundary_note;
            }
            buffer += (is_forward ? "\tTS:A:+\n" : "\tTS:A:-\n");
            buffer += read_seq;
            buffer += '\n';
            
            if (is_fastq) {
                buffer += "+\n";
                buffer += read_qual;
                buffer += '\n';
            }
            ++records_written;
        }
        return records_written;
    }

/** 
 * @brief convert SigString to a fastqa text representation
 * @return `std::string` fastqa representation
 * 
 * @brief This does all of the things that to_fastqa_append does, but instead of appending to a buffer,
 * it constructs and returns a new string containing the fastqa representation of the SigString. This is
 * an older, less efficient version of the writing process that's used in other contexts where we needed the string.
*/
    std::string to_fastqa(bool rc_umi = true) const {
        std::vector<std::string> dirs;
        if(read_type == "skipped"){
            return ""; // Skip if read type is "skipped"
        }
        if (read_type == "concatenate") {
            dirs = {
                "forward",
                "reverse"
            };
        } else {
            dirs = {
                read_type
            };
        }
        std::string all_records;
        for (const auto& dir : dirs) {
            std::vector<std::string> bc_keys;
            std::unordered_map<std::string, std::string> bc_map;
            std::unordered_map<std::string, std::string> bc_dir;
            std::unordered_map<std::string, std::string> cr_map;
            std::string umi, read_seq, read_qual;
            
            // First pass: collect all elements for this direction
            for (const auto& elem : sig_elements) {
                if (!elem.seq.has_value()) continue;
                if (elem.direction != dir) continue;
                if (elem.global_class == "barcode") {
                    if (elem.write.has_value() && !elem.write.value()) {
                        continue; // Skip barcodes explicitly marked not to write
                    }
                    auto key = seq_utils::remove_rc(elem.class_id);
                    bool is_fwd = (dir == "forward");
                    if (bc_map.find(key) == bc_map.end()) {
                        bc_keys.push_back(key);
                        bc_map[key] = elem.seq.value();
                        bc_dir[key] = dir;
                        if (elem.original_seq.has_value()) {
                            cr_map[key] = elem.original_seq.value();
                        }
                    } else if (is_fwd && bc_dir[key] == "reverse") {
                        bc_map[key] = elem.seq.value();
                        bc_dir[key] = dir;
                        if (elem.original_seq.has_value()) {
                            cr_map[key] = elem.original_seq.value();
                        }
                    }
                    continue;
                }
                if (elem.global_class == "umi") {
                    if (elem.seq.has_value()) {
                        // Match to_fastqa_append: flip minus-strand (reverse) UMIs
                        // back to plus-strand so UB:Z and CB:Z share an orientation.
                        umi = (rc_umi && dir == "reverse")
                                  ? seq_utils::revcomp(elem.seq.value())
                                  : elem.seq.value();
                    }
                    continue;
                }
                if (elem.global_class == "read") {
                    read_seq = elem.seq.value();
                    if (elem.qual.has_value()) {
                        read_qual = elem.qual.value();
                    }
                    continue;
                }
                if (elem.global_class == "poly_tail" || elem.global_class == "start" || elem.global_class == "stop") {
                    continue;
                }
            }

            // In bulk mode with no barcodes, use N-masked regions from static elements as barcodes
            if (bc_map.empty() && additional_info == "bulk") {
                for (const auto& elem : sig_elements) {
                    if (elem.type != "static" || elem.direction != dir) {
                        continue;
                    }
                    // Static elements with original_seq have N-masked regions extracted
                    if (elem.original_seq.has_value() && elem.seq.has_value()) {
                        auto key = seq_utils::remove_rc(elem.class_id);
                        if (bc_map.find(key) == bc_map.end()) {
                            bc_keys.push_back(key);
                            bc_map[key] = elem.seq.value();  // N-masked region(s)
                            bc_dir[key] = dir;
                            cr_map[key] = elem.original_seq.value();  // Full aligned sequence
                        }
                    }
                }
            }

            // Check if we have essential components - skip this direction if not
            if ((bc_map.empty() && additional_info != "bulk") || read_seq.empty()) {
                continue; // Skip this direction if either barcode or read are empty
            }

            std::string cb_tag, cr_tag;
            for (size_t i = 0; i < bc_keys.size(); ++i) {
                const auto& key = bc_keys[i];
                if (i) {
                    cb_tag += '-';
                    if (!cr_tag.empty()) cr_tag += '-';
                }
                cb_tag += bc_map[key];
                if (cr_map.find(key) != cr_map.end()) {
                    cr_tag += cr_map[key];
                }
            }
            bool is_fastq = !read_qual.empty();
            bool is_concatenate = from_concatemer || (read_type == "concatenate");
            bool is_forward = (dir == "forward");
            std::stringstream ss;
            //generating modified sequence id
            ss << (is_fastq ? '@' : '>') << sequence_id << (is_forward ? "-F" : "-R") << (is_concatenate ? "-CT" : "");
            //adding barcode tag
            if (!cb_tag.empty()) ss << "\tCB:Z:" << cb_tag;
            //adding corrected read tag for SAM
            //added a fix here so that it's left empty for reverse complement fixes as well, otherwise RCs show up and it's annoying to parse later
            if (!cr_tag.empty() && (cb_tag != cr_tag && seq_utils::revcomp(cb_tag) != cr_tag)) {
                ss << "\tCR:Z:" << cr_tag;
            } else {
                // Ensure CR tag is present even if empty
                ss << "\tCR:Z:";
            }
            //adding transcript tag for SAM
            if (!umi.empty()) ss << "\tUB:Z:" << umi;
            if (!hmm_boundary_note.empty()) ss << "\tHB:Z:" << hmm_boundary_note;
            if(is_forward){
                ss << "\tTS:A:+";
            } else {
                ss << "\tTS:A:-";
            }
            ss << "\n" << read_seq << "\n";
            if (is_fastq) {
                ss << "+\n" << read_qual << "\n";
            }
            all_records += ss.str();
        }
        return all_records;
    }

/**
 * @brief convert SigString to a CSV text representation
 * @param write_header `bool` whether to include the CSV header
 * @return `std::string` CSV representation
 * 
 * @brief This function converts the SigString object into a CSV. Great for debugging and per-element processing.
 * It includes an optional header row. It processes only variable elements with non-empty sequences.
 */
    std::string to_csv(bool write_header = false) const {
        std::stringstream ss;
        if (write_header) {
            ss << "id,elem,read_type,seq\n";
        }
        for (const auto& elem : sig_elements) {
            // Process only variable elements that have a non-empty sequence.
            if (elem.type == "variable" && elem.seq.has_value() && !elem.seq->empty()) {
                // Determine the element's direction.
                std::string read_mintype = "forward";
                std::string seq = elem.seq.value();
                std::string final_id = elem.class_id;
                if (elem.class_id.substr(0, 3) == "rc_") {
                    seq = seq_utils::revcomp(seq);
                    final_id = elem.class_id.substr(3);
                    read_mintype = "reverse";
                }
                // If the overall read type is "forward", skip reverse elements.
                if (read_type == "forward" && read_mintype != "forward")
                    continue;
                // If the overall read type is "reverse", skip forward elements.
                if (read_type == "reverse" && read_mintype != "reverse")
                    continue;
                if(from_concatemer || read_type == "concatenate"){
                    read_mintype = read_mintype + "_concatenate";
                }
                // If read_type is "concatenate", include both.
                ss << sequence_id << "," 
                << final_id << "," 
                << read_mintype << "," 
                << seq << "\n";
            }
        }
        return ss.str();
    }

};
