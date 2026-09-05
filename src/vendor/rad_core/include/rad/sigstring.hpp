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
    int edit_distance;
    bool success;
    std::string seq;
    std::string cigar; 
    int score;
    std::pair<int,int> pos; 
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
        int expected_end = -1) {

        static_alignments alignment;
        alignment.success = false;
        alignment.edit_distance = -1;
        // Our threshold for a “good” alignment
        const int min_match_bases = 5;
        
        // If expected_end is not provided, use the target's length
        if (expected_end < 0){
                expected_end = target.size();
        }

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
                                                target.c_str(), target.size(),
                                                config);

        std::vector<std::pair<int,int>> candidates;
        if (edlibResult.status == EDLIB_STATUS_OK &&
            edlibResult.numLocations > 0 &&
            edlibResult.editDistance > -1) {
            std::vector<std::pair<int,int>> intervals;
            for (int i = 0; i < edlibResult.numLocations; ++i) {
                int start = edlibResult.startLocations[i] + 1;  // Convert to 1-indexed
                int end   = edlibResult.endLocations[i] + 1;
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
        }

        // --- early out: If Edlib returned exactly one candidate and it is within the expected region, use it
        if (candidates.size() == 1) {
            auto cand = candidates.front();
            if (cand.first >= expected_start && cand.second <= expected_end) {
                alignment.positions.push_back(cand);
                alignment.success = true;
                alignment.edit_distance = edlibResult.editDistance;
                // Optionally set alignment.seq, alignment.cigar, etc.
                alignment.seq = target.substr(cand.first - 1, cand.second - cand.first + 1);
                // Return immediately without SSW
                edlibFreeAlignResult(edlibResult);
                return alignment;
            }
        }
        edlibFreeAlignResult(edlibResult);
        // --- run SSW on the entire target (for no candidate or multiple candidates) ---
        // Prepare uppercase strings for SSW.
        std::string query_ssw = query;
        std::string target_ssw = target;
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
            int ssw_start = sswAlign.query_begin + 1; // Convert to 1-indexed
            int ssw_end = sswAlign.query_end + 1;
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
                alignment.success = alignment.edit_distance < max_edit_distance + 1 ? true : false;
                alignment.cigar = sswAlign.cigar_string;
                alignment.seq = target.substr(ssw_start - 1, ssw_end - ssw_start + 1);
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
};

/**
 * @namespace barcode_correction
 * @brief Namespace for barcode correction functions and utilities
 */
namespace barcode_correction {
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
            int filtered_count = all_counts.load(barcode_counts::filtered);
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
        int max_dist, bool verbose, const std::string& wl_type, const whitelist::wl_entry* wl = nullptr
    ) {
        auto sorted = mutation_tools::int64_lvdist(query, candidates, max_dist);
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
    std::optional<int64_seq> check_against_wl(const int64_seq& bc, const std::unordered_set<int64_seq>& candidates,
        const std::string& whitelist_type, int max_dist, bool verbose, const std::string& mode, const whitelist::wl_entry* wl = nullptr
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
            int res = mutation_tools::int64_lvdist(bc, putative_candidate, max_dist);
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
            auto [resolved, min_dist] = resolve_multiple_hits_simple(bc, matched, max_dist, verbose, whitelist_type, wl);
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
        int max_dist,bool verbose, const std::string& mode, const whitelist::wl_entry* wl = nullptr
    ) {
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

            auto unique_entries = typed_wl.get_unique_entries();
            for (const auto* entry : unique_entries) {
                int res = mutation_tools::int64_lvdist(bc, entry->barcode, max_dist);
                if (res >= 0) {
                    match.emplace_back(entry->barcode, res);
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
        const std::string& mode, const whitelist::wl_entry& wl, bool verbose, int max_dist
    ) {
    
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
            for (const auto* e : seed_entries) {
                if (e) seed_candidates.insert(e->barcode);
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
                return seed_result;
            }
        }
    }

    int64_seq exp_bc;
    exp_bc.sequence_to_bits(expanded_seq);
    auto true_result = exhaustive_check_against_wl(exp_bc, "true", 2, verbose, mode, &wl);
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
 * @brief Correct a barcode sequence using the provided layout and read information
 * @param elem ``const seq_element&`` Sequence element containing barcode information
 * @param layout ``const ReadLayout&`` Read layout containing whitelist mappings
 * @param full_read ``const read_streaming::sequence&`` Full read sequence
 * @param verbose ``bool`` Whether to print verbose output
 * @param mut_dist ``int`` Maximum allowed edit distance for mutation checks
 * @param mode ``std::string`` Barcode correction mode ("defensive" or "offensive")
 * @return ``std::optional<int64_seq>`` Corrected barcode sequence if found, otherwise std::nullopt
 */
    std::optional<int64_seq> correct_barcode(const seq_element& elem,  const ReadLayout& layout, 
        const read_streaming::sequence& full_read,  bool verbose, int mut_dist, std::string mode) {
        // pick the right whitelist
        auto key = seq_utils::remove_rc(elem.class_id);
        auto wl_it = layout.wl_map.maps.find(key);
        if (wl_it == layout.wl_map.maps.end()) return std::nullopt;
        auto &wl = wl_it->second.get();
        int max_dist = 4;
        // extract and reverse‐complement the raw string
        // encoded barcode and reverse complement
        std::string raw = elem.seq.value();
        std::string expanded_seq = seq_utils::substr_w_padding(full_read.seq, elem.position.first, elem.position.second, max_dist);

        if (elem.direction == "reverse") {
            raw = seq_utils::revcomp(raw);
            expanded_seq = seq_utils::revcomp(expanded_seq);
        }

        int64_seq bc, rc_bc, exp_bc, exp_rcbc;

        bc.sequence_to_bits(raw);
        rc_bc.sequence_to_bits(seq_utils::revcomp(raw));
        exp_bc.sequence_to_bits(expanded_seq);
        exp_rcbc.sequence_to_bits(seq_utils::revcomp(expanded_seq));

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

      auto kmer_fuzzy_result = kmer_fuzzy_search(bc, expanded_seq, bc_len, mode, wl, verbose, 3);
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
        
        auto muts = mutation_tools::generate_mutated_barcodes(bc, mut_dist);
        if(mode == "defensive"){
            auto global_result = check_against_wl(exp_bc, muts, "global", max_dist, verbose, mode, &wl);
            if(global_result.has_value()){
                return(global_result);
            }

            auto true_result = check_against_wl(exp_bc, muts, "true", 2, verbose, mode, &wl);
            if(true_result.has_value()){
                return(true_result);
            }
        }

        // IF OFFENSIVE: Check against true first, and then global
        
        if(mode == "offensive") {
            auto true_result = check_against_wl(exp_bc, muts, "true", 2, verbose, mode, &wl);
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

            auto global_result = check_against_wl(exp_bc, muts, "global", max_dist, verbose, mode, &wl);
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
};

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
 * @brief Filter out reads that are shorter than the summative expected length of adapters
 * @param elems vector of references to `seq_element` objects
 * @param layout `ReadLayout` object containing layout elements
 * @return true if read length meets or exceeds total expected adapter length, false otherwise
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

        // Length filtering
        if (!filter_short_reads(elements, layout, min_read_length)) {
            filtered_because += direction + "_FILTERED_READ_LENGTH";
            set_info(filtered_because);
            return false;
        }
        
        // Sort elements by order
        std::sort(elements.begin(), elements.end(), [](const auto& a, const auto& b) {
            return a.get().order < b.get().order;
        });
        
        // Mask overlapping elements and trim reads
        mask_and_trim_elements(elements, read, verbose);
        
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
        return wl.true_bcs.check_wl_for(bits) || wl.global_bcs.check_wl_for(bits);
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
        bool verbose
    ) {
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

        std::function<bool(int, size_t)> search_chain = [&](int cursor, size_t depth) -> bool {
            if (depth == order.size()) return true; // all matched
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
                if (search_chain(end + 1, depth + 1)) return true;
            }
            return false;
        };

        bool found = false;
        for (int off : offset_order) {
            int start = anchor_start + off;
            if (start < 1) continue;
            if (search_chain(start, 0)) {
                found = true;
                break;
            }
        }
        if (!found) return false;

        // Validate the pair against the joint_spat_wl if loaded.
        // For a two-member group: check the (BC1, BC2) bit in the sparse mask.
        if (group.size() == 2) {
            size_t first_idx = order[0];
            size_t second_idx = order[1];
            auto wl_it = layout.wl_map.maps.find(
                seq_utils::remove_rc(members[first_idx].elem->class_id));
            if (wl_it != layout.wl_map.maps.end()) {
                const auto& wl_entry = wl_it->second.get();
                if (wl_entry.spat_wl.has_value() && !wl_entry.spat_wl->empty()) {
                    int64_seq bc1_bits, bc2_bits;
                    bc1_bits.sequence_to_bits(hits[first_idx].seq);
                    bc2_bits.sequence_to_bits(hits[second_idx].seq);
                    if (!wl_entry.spat_wl->check(bc1_bits, bc2_bits)) {
                        if (verbose) {
                            log_verbose("JOINT_BARCODE_PAIR_REJECTED: " +
                                hits[first_idx].seq + "+" + hits[second_idx].seq +
                                " not in valid pairs mask");
                        }
                        return false;
                    }
                }
            }
        }

        // Convert hits back to original coordinates and apply corrections
        auto to_original = [read_len, reverse](std::pair<int, int> p) -> std::pair<int, int> {
            if (!reverse) return p;
            return std::make_pair(read_len - p.second + 1, read_len - p.first + 1);
        };

        auto& id_index = sig_elements.get<sig_id_tag>();
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

                if (group.size() >= 2 &&
                    try_process_joint_barcode_group(group, layout, read, group_info, verbose)) {
                    any_barcode_passed = true;
                    i = j;
                    continue;
                }
                // Joint group failed — fall through to process individually
            }

            // Process as a single barcode
            auto it = id_index.find(entry.class_id);
            if (it == id_index.end()) {
                ++i;
                continue;
            }
            seq_element elem = *it;
            bool barcode_passed = process_single_barcode(elem, layout, read, gen_mut, mode, verbose);
            if (!barcode_passed) {
                all_barcodes_passed = false;
                mark_element_failed(elem.class_id);
            } else {
                any_barcode_passed = true;
            }
            ++i;
        }

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
        int gen_mut, const std::string& mode, bool verbose
    ) {

        if (verbose) {
            log_verbose("Processing barcode: " + elem.seq.value());
        }
        
        auto correction_result = barcode_correction::correct_barcode(
            elem, layout, read, verbose, gen_mut, mode
        );

        // Apply correction and update element
        if (!correction_result.has_value()) {
            return false;
        }

        apply_barcode_correction(elem, correction_result.value(), layout, verbose);
        return true;
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
        barcode_entry* entry = nullptr;
        bool found_in_true = false;
        bool found_in_global = false;

        entry = wl.true_bcs.find_entry(corrected_barcode);
        if (entry) {
            found_in_true = true;
        } else {
            entry = wl.global_bcs.find_entry(corrected_barcode);
            if (entry) found_in_global = true;
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
            e.resolved_counter = entry ? &entry->count : nullptr;
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
                
                barcode_entry failed_entry;
                failed_entry.barcode = failed_bc;
                failed_entry.filtered = true;
                
                // Double-check that it's still not present
                auto failed_range = wl.filter_bcs.equal_range(failed_bc);
                auto rc_failed_range = wl.filter_bcs.equal_range(rc_failed_bc);
                
                if ((failed_range.first == failed_range.second) && 
                    (rc_failed_range.first == rc_failed_range.second)) {
                    // Entry doesn't exist - insert it
                    //wl.filter_bcs.insert_bc_entry(failed_bc, failed_entry);
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
    
    void log_verbose(const std::string& message) {
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
   void sigalign_static(const read_streaming::sequence &read, const ReadLayout& layout, bool verbose) {
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
    size_t read_length = read.seq.length();

    for (auto it = static_range.first; it != static_range.second; ++it) {
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
            const auto& [mis_start, start_mvar] = it->misaligned_positions->start_stats;
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
        auto result = aligner.align_static_elements(it->seq,
                                                    mutable_seq,
                                                    verbose,
                                                    max_distance,
                                                    it->masked_seq,
                                                    verbose,
                                                    expected_start,
                                                    expected_end);
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

                    add_element(seq_element(
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
                    ));

                    // Mask the aligned region so that it is not re-aligned
                    std::fill(mutable_seq.begin() + adj_start - 1,
                              mutable_seq.begin() + adj_end, 'X');
                }
                // Primary region yielded a valid match; skip further processing for this element
                continue;
            }
        }
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
        
        // Process forward non-read variables
        std::multimap<std::string, const seq_element*> static_refs;
        for (const auto& elem : sig_elements) {
            if (elem.type == "static") {
                static_refs.insert({elem.class_id, &elem});
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
            total_refs.insert({elem.class_id, &elem});
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
        
        constexpr char qual_mask = '\x7F';
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
        std::function<void()> boundary_poll = {}
    ) {
    const auto sigalign_wall_t0 = std::chrono::steady_clock::now();

    std::string file_out = path_utils::get_fastqa_type(fastq_path);
    std::string fastq_output_path = output_prefix + file_out;
    bool compress_fastq = true;

    std::unique_ptr<std::ofstream> metrics_file_ptr;
    std::unique_ptr<sigstring_writing> debug_sig_writer, debug_csv_writer, debug_fastqa_writer;

    // Primary FASTQA writer
    sigstring_writing fastqa_writer(fastq_output_path, sigstring_writing::format::FASTQA, compress_fastq, false, num_threads);

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

    auto process_chunk = [&](std::vector<read_streaming::sequence>& chunk, 
                             const std::string& path)
    {
        if (chunk.empty()){
            return;
        }

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
        std::atomic<size_t> passed_count{0};
        std::atomic<size_t> demultiplexed_count{0};
        std::atomic<size_t> records_written_count{0};

        std::vector<std::optional<SigString>> debug_by_input;
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
                if (processing_failed.load(std::memory_order_acquire)) {
                    continue;
                }
                try {
                    const auto& read = chunk[i];
                    SigString sig(read.id, read.seq.length(),"undefined", layout.sequencing_type);
                    sig.sigalign_static(read, layout, verbose);
                    sig.sigalign_variable(read, layout, verbose);
                    sig.sigalign_filter(read, layout, gen_mut.value_or(2), verbose, mode, joint_bc_mode, min_read_length);

                    // Serialize directly to the worker-local buffer, retaining
                    // the byte range for this input read.  One read can emit
                    // multiple FASTQ records; keeping them in one span also
                    // preserves their intra-read order.
                    if (sig.read_type != "filtered" && sig.read_type != "skipped") {
                        passed_count.fetch_add(1, std::memory_order_relaxed);
                        const size_t offset = thread_buffers[tid].size();
                        const size_t records_written =
                            sig.to_fastqa_append(thread_buffers[tid], rc_umi);
                        if (records_written > 0) {
                            thread_spans[tid].push_back(serialized_read_span{
                                i,
                                static_cast<size_t>(tid),
                                offset,
                                thread_buffers[tid].size() - offset
                            });
                            demultiplexed_count.fetch_add(
                                1, std::memory_order_relaxed);
                            records_written_count.fetch_add(
                                records_written, std::memory_order_relaxed);
                        }
                    }
                    if (write_debug) {
                        debug_by_input[i].emplace(std::move(sig));
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
                if (sig.has_value()) {
                    debug_sigs.push_back(std::move(sig.value()));
                }
            }
        }

        const auto wall_t1 = std::chrono::steady_clock::now();
        
        // Calculate actual work time (max across all threads)
        double actual_process_ms = 0.0;
        for (double t : thread_times) {
            actual_process_ms = std::max(actual_process_ms, t);
        }

        // === consolidate and write ===
        if (passed_count > 0) {
            std::vector<serialized_read_span> ordered_spans;
            ordered_spans.reserve(demultiplexed_count.load(
                std::memory_order_relaxed));
            size_t total_size = 0;
            for (auto& local : thread_spans) {
                for (const auto& span : local) {
                    total_size += span.length;
                    ordered_spans.push_back(span);
                }
            }
            std::sort(
                ordered_spans.begin(), ordered_spans.end(),
                [](const serialized_read_span& lhs,
                   const serialized_read_span& rhs) {
                    return lhs.input_index < rhs.input_index;
                });

            std::string consolidated;
            consolidated.reserve(total_size);

            for (const auto& span : ordered_spans) {
                consolidated.append(
                    thread_buffers[span.worker_index],
                    span.offset,
                    span.length);
            }
            writer.write_raw_string(fastqa_writer, std::move(consolidated));
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
                    const std::string& class_id = kv.first;
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

    fastqa_writer.finalize();
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
    RAD_COUT << "[sigalign] Timing breakdown:\n";
    RAD_COUT << "  - Accounted chunk runtime: " << accounted_s << " seconds\n";
    RAD_COUT << "  - Process time: " << process_s << " seconds\n";
    RAD_COUT << "  - Output staging time: " << queue_s << " seconds\n";
    RAD_COUT << "  - Overhead (writer drain/setup): " << overhead_s << " seconds\n";

    if (write_debug) {
        RAD_COUT << "\n[sigalign] Output written to:\n"
                  << "[sigalign] [sigstring]: " << output_prefix << ".sig\n"
                  << "[sigalign]       [csv]: " << output_prefix << ".csv\n"
                  << "[sigalign]     [fastq]: " << fastq_output_path << (compress_fastq ? ".gz" : "") << "\n"
                  << "[sigalign]   [metrics]: " << output_prefix << ".metrics.tsv\n";
    } else {
        RAD_COUT << "\n[sigalign] Output written to:\n"
                  << "[sigalign][fastq]: " << fastq_output_path << (compress_fastq ? ".gz" : "") << "\n";
    }

    return sigalign_run_stats{
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
            for (const auto& elem_ref : sorted_elements) {
                const auto& elem = elem_ref.get();

                if (!first_elem) {
                    ss << "|";
                }
                if(elem.position.first == -1 && elem.position.second == -1){
                    continue; // Skip elements with invalid positions
                }
                first_elem = false;
                ss << elem.class_id << ":"
                << (elem.edit_distance ? std::to_string(elem.edit_distance.value()) : "0") << ":"
                << elem.position.first << ":"
                << elem.position.second;
            }
            // Append final tag: direction abbreviated as F or R
            std::string dir = (direction == "forward") ? "F" : "R";
            std::string concat = read_type == "concatenate" ? ":C" : "";
            std::string combined_info = "";
            if (read_type == "filtered") {
                combined_info += ":filtered";
            }
            if (additional_info.length() > 2) {
                combined_info += ":" + additional_info;
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
            bool is_concatenate = (read_type == "concatenate");
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
            bool is_concatenate = (read_type == "concatenate");
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
                if(read_type == "concatenate"){
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
