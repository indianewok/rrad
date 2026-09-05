#pragma once
#include "rad_headers.h"

/**
 * @enum barcode_counts
 * @brief enumeration for different barcode count types
 * This enum defines indices for various barcode count categories used in the `counter` class.
 * The categories include raw counts, forward and reverse counts, concatenation counts,
 * total counts, corrected counts, and filtered counts.
 */
enum barcode_counts {
    raw = 0,          // Raw count
    forw = 1,         // Forward count
    forw_concat = 2,  // Forward concatenation count
    rev = 3,          // Reverse count
    rev_concat = 4,   // Reverse concatenation count
    total = 5,         // Total count
    corrected = 6,     // Corrected count
    filtered = 7       // filtered count
};

/**
 * @struct counter
 * @brief Thread-safe counter class using atomic integers for multiple barcode count types
 * This class provides atomic operations for counting various barcode-related metrics.
 * It uses an array of 8 atomic integers to store counts for different categories defined in the
 * `barcode_counts` enum. The class supports incrementing, decrementing, loading individual counts,
 * and loading all counts as a tuple.
 */
struct counter {
    // Use an array of 8 atomic ints
    std::array<std::atomic<int>, 8> counts;

    // Default constructor: initialize all counters to 0.
    counter() : counts{{0, 0, 0, 0, 0, 0, 0, 0}} { }

    // A constructor that initializes all counters to a given value
    counter(int init) : counts{{init, init, init, init, init, init, init, init}} { }

    // Copy constructor: copy each counter's current value
    counter(const counter &other) {
        for (size_t i = 0; i < counts.size(); ++i) {
            counts[i].store(other.counts[i].load(std::memory_order_relaxed),
                              std::memory_order_relaxed);
        }
    }

    // Assignment operator: copy each counter's value
    counter& operator=(const counter &other) {
        for (size_t i = 0; i < counts.size(); ++i) {
            counts[i].store(other.counts[i].load(std::memory_order_relaxed),
                              std::memory_order_relaxed);
        }
        return *this;
    }

    // Increment the counter at the specified index
    void increment(int index) {
        counts[index].fetch_add(1, std::memory_order_relaxed);
    }

    // Overload: increment using one of the enum values for clarity
    void increment(barcode_counts index) {
        counts[static_cast<int>(index)].fetch_add(1, std::memory_order_relaxed);
    }

    // Decrement the counter at the specified slot
    void subtract(int index) {
        counts[index].fetch_add(-1, std::memory_order_relaxed);
    }
    
    // Overload: decrement using an enum value for clarity
    void subtract(barcode_counts index) {
        counts[static_cast<int>(index)].fetch_add(-1, std::memory_order_relaxed);
    }

    // Load the value of the counter at the specified index
    int load(int index) const {
        return counts[index].load(std::memory_order_relaxed);
    }

    // Overload: load using an enum value for clarity
    int load(barcode_counts index) const {
        return counts[static_cast<int>(index)].load(std::memory_order_relaxed);
    }

    // Load all counter values into a tuple
    std::tuple<int, int, int, int, int, int, int, int> load_all() const {
        return std::make_tuple(
            counts[barcode_counts::raw].load(std::memory_order_relaxed),
            counts[barcode_counts::forw].load(std::memory_order_relaxed),
            counts[barcode_counts::forw_concat].load(std::memory_order_relaxed),
            counts[barcode_counts::rev].load(std::memory_order_relaxed),
            counts[barcode_counts::rev_concat].load(std::memory_order_relaxed),
            counts[barcode_counts::total].load(std::memory_order_relaxed),
            counts[barcode_counts::corrected].load(std::memory_order_relaxed),
            counts[barcode_counts::filtered].load(std::memory_order_relaxed)
        );
    }
};

/**
 * @class int64_seq
 * @brief Class to represent barcodes with unsigned 2-bit encoding (up to 32 bases)
 * @param length `uint16_t` total number of bases in the sequence
 * @param bits `std::vector<uint64_t>` vector containing the packed sequence
 * @note Each nucleotide is encoded as follows: A=00, C=01, T=10, G=11
 */
class int64_seq {
    public:
        uint16_t length; // `uint16_t` total number of bases in the sequence
        // Barcode algorithms in this header operate on one 64-bit word. Keep
        // the representation unsigned so packing and shifting a 32-base
        // barcode never invokes signed-overflow or signed-shift behaviour.
        std::vector<uint64_t> bits;
        int64_seq() : length(0), bits{} {}

        explicit int64_seq(uint64_t raw_bits, uint16_t seq_length) : length(seq_length) {
            if (seq_length == 0 || seq_length > 32) {
                throw std::invalid_argument(
                    "barcode length must be between 1 and 32 bases");
            }
            bits.reserve(1);
            bits.push_back(raw_bits);
        }
    
        explicit int64_seq(const std::string& sequence) {
            sequence_to_bits(sequence);
        }
    
        void sequence_to_bits(const std::string& sequence) {
            if(sequence.empty()) {
                length = 0;
                bits = {0};
                return;
            }
            if (sequence.size() > 32) {
                throw std::invalid_argument(
                    "barcode length must be between 1 and 32 bases");
            }
            length = static_cast<uint16_t>(sequence.size());
            uint64_t result = 0;
            for (size_t j = 0; j < sequence.size(); ++j) {
                char c = sequence[j];
                result <<= 2;
                switch (c) {
                case 'A': break;
                case 'C': result |= 1; break;
                case 'T': result |= 2; break;
                case 'G': result |= 3; break;
                default:
                    length = 0;
                    bits.clear();
                    bits = {0};
                    return; // Invalid character, reset and exit
                }
            }
            bits = {result};
        }
    
        std::string bits_to_sequence() const {
            if (!is_valid()) return {};
            uint64_t int64_code = bits.front();
            std::string result;
            const int bases_to_decode = length;
            for (int j = 0; j < bases_to_decode; ++j) {
                switch (int64_code & 3) {
                case 0: result.insert(result.begin(), 'A'); break;
                case 1: result.insert(result.begin(), 'C'); break;
                case 2: result.insert(result.begin(), 'T'); break;
                case 3: result.insert(result.begin(), 'G'); break;
                }
                int64_code >>= 2;
            }
            return result;
        }
        
        bool operator==(const int64_seq &other) const {
            return (length == other.length) && (bits == other.bits);
        }
    
        void print_int64_seq() {
            RAD_COUT << "Sequence length: " << length << "\n";
            RAD_COUT << "Bit chunks: ";
            for (size_t i = 0; i < bits.size(); ++i) {
                RAD_COUT << bits[i] << " ";
            }
            RAD_COUT << "\n";
        }

        bool is_valid() const noexcept {
            return length > 0 && !bits.empty();
        }
    };

/**
 * @struct barcode_entry
 * @brief Struct representing a barcode entry with its sequence, count, and filtering status
 * @param barcode `int64_seq` representing the barcode sequence
 * @param count `counter` atomic counter for tracking counts
 * @param filtered `bool` flag indicating if the barcode is filtered
 */
struct barcode_entry {
    int64_seq barcode;
    mutable counter count;
    bool filtered;
    
    barcode_entry()
      : barcode(0,1), 
       count(0),
       filtered(false) 
       { }

    bool operator==(const barcode_entry &o) const noexcept {
        return barcode == o.barcode;
    }

    bool is_valid() const noexcept {
        return barcode.is_valid();
    }

};

// Hash specializations must live in the real global std namespace.  The RAD
// engine is wrapped by rrad_rad_core in the embedded build, so step outside
// that namespace for these declarations and then resume it below.
}  // namespace rrad_rad_core

namespace std {
    template <>
    struct hash<rrad_rad_core::int64_seq> {
        std::size_t operator()(const rrad_rad_core::int64_seq &seq) const {
             // Skip length hashing
            if (seq.bits.empty()) {
                // Fallback safe hash for empty sequences
                return std::hash<std::string>{}("EMPTY_SEQ");
            }
            return std::hash<uint64_t>{}(seq.bits[0]);
        }
    };

    template<>
    struct hash<rrad_rad_core::barcode_entry> {
        size_t operator()(rrad_rad_core::barcode_entry const& bc) const noexcept {
            // hash the barcode field
            return std::hash<rrad_rad_core::int64_seq>()(bc.barcode);
        }
    };
}

namespace rrad_rad_core {

namespace mutation_tools {
/** 
    * @brief Myers's bit-parallel algorithm for edit distance, followed from the edlib comments for making a better version for shorter strings and taking int64 values as input 
    * @param p:  `uint64_t` packed pattern
    * @param t:  `uint64_t` packed text
    * @param n:  number of barcode bases (<=32)
    * @param max_dist: as soon as even the best possible remaining score exceeds max_dist, return -1 as a sentinel
    * @return `int` minimum edit distance if <= max_dist, else -1
*/
    int bit_ld(uint64_t pattern, uint64_t text, int n, int max_dist) {
    if (n <= 0 || n > 32 || max_dist < 0) return -1;
    // build pattern equality vectors
    uint64_t Peq[4] = {0, 0, 0, 0};
    for (int i = 0; i < n; i++) {
        int nuc = (pattern >> (2 * i)) & 3;
        Peq[nuc] |= (uint64_t{1} << i);
    }
    
    // initialize
    uint64_t Pv = (uint64_t{1} << n) - 1;  // All 1s
    uint64_t Mv = 0;                       // All 0s
    int score = n;
    
    // Step 3: process text character
    for (int j = 0; j < n; j++) {
        int text_nuc = (text >> (2 * j)) & 3;  // Extract nucleotide
        uint64_t Eq = Peq[text_nuc];             // Get equality vector
        
        // Myers core computation
        uint64_t Xv = Eq | Mv;
        uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
        
        uint64_t Ph = Mv | ~(Xh | Pv);
        uint64_t Mh = Pv & Xh;
        
        // Update score
        if (Ph & (uint64_t{1} << (n-1))) score++;
        if (Mh & (uint64_t{1} << (n-1))) score--;
        
        // Early exit
        if (score - (n - j - 1) > max_dist) return -1;
        
        // Update for next iteration
        Ph <<= 1;
        Pv = (Mh << 1) | ~(Xv | Ph);
        Mv = Ph & Xv;
    }
    return score;
}

/**
 * @brief Myers's bit-parallel algorithm for partial edit distance (finding best match of pattern within text)
 * @param pattern: `int64_t` bitmask where bit j=1 means pattern[j] matches 1
 * @param text: `int64_t` bitmask for text only used to select p or np per column
 * @param pattern_len: `int` length of the pattern (currently <=32 bases)
 * @param text_len: `int` length of the text (currently <=32 bases)
 * @param max_dist: maximum allowed distance
 * @return `int` minimum edit distance if <= max_dist, else -1
 */
    int bit_partial_match(uint64_t pattern, uint64_t text, int pattern_len, int text_len, int max_dist) {
        if (pattern_len <= 0 || text_len <= 0 || pattern_len > 32 ||
            text_len > 32 || max_dist < 0) {
            return -1;
        }
        
        // Build pattern equality vectors
        uint64_t Peq[4] = {0, 0, 0, 0};
        for (int i = 0; i < pattern_len; i++) {
            int nuc = (pattern >> (2 * i)) & 3;
            Peq[nuc] |= (uint64_t{1} << i);
        }
        
        uint64_t pattern_mask = (uint64_t{1} << pattern_len) - 1;
        uint64_t Pv = pattern_mask;
        uint64_t Mv = 0;
        int score = pattern_len;
        // Track minimum across all positions
        int min_score = pattern_len;
        
        for (int j = 0; j < text_len; j++) {
            int text_nuc = (text >> (2 * j)) & 3;
            uint64_t Eq = Peq[text_nuc] & pattern_mask;
            
            uint64_t Xv = Eq | Mv;
            uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
            uint64_t Ph = Mv | ~(Xh | Pv);
            uint64_t Mh = Pv & Xh;
            
            Ph &= pattern_mask;
            Mh &= pattern_mask;
            
            if (Ph & (uint64_t{1} << (pattern_len - 1))) score++;
            if (Mh & (uint64_t{1} << (pattern_len - 1))) score--;
            // Key difference from full matching
            min_score = std::min(min_score, score);  
            // Perfect match found
            if (min_score == 0) return 0; 
            
            Ph <<= 1;
            Pv = ((Mh << 1) | ~(Xv | Ph)) & pattern_mask;
            Mv = Ph & Xv;
        }
        return (min_score <= max_dist) ? min_score : -1;
    }

/**
 * @brief Calculate Levenshtein distance between one query and a set of multiple strings (partial matching)
 * @param query: the query sequence as an `int64_seq`
 * @param targets: a `unordered_set<int64_seq>` of target sequences
 * @param max_dist: maximum allowed distance
 * @return `map<int, unordered_set<int64_seq>>` from edit distance to sets of sequences within that distance
 */
    std::map<int, std::unordered_set<int64_seq>> int64_lvdist(const int64_seq &query, const std::unordered_set<int64_seq> &targets, int max_dist = 4) {
        std::map<int, std::unordered_set<int64_seq>> results;
        if (query.bits.empty()){
            return results;
        }
        for (const auto &target : targets) {
            if (target.bits.empty()){
                continue;
            }
            
            int dist;
            if (query.length <= target.length) {
                // Query is shorter or equal, search query in target
                dist = bit_partial_match(query.bits[0], target.bits[0], query.length, target.length, max_dist);
            } else {
                // Target is shorter, search target in query
                dist = bit_partial_match(target.bits[0], query.bits[0], target.length, query.length, max_dist);
            }
            
            if(dist >= 0){
                results[dist].insert(target);
            }
        }
        return results;
    }

/**
 * @brief Calculate Levenshtein distance between one query and a set of barcode entries (partial matching)
 * @param query The query sequence to search for as an `int64_seq`
 * @param targets a `vector<barcode_entry>` of target sequences
 * @param max_dist Maximum allowed edit distance (default: 4)
 * @return `map<int, unordered_set<int64_seq>>` from edit distance to sets of sequences within that distance
 */
    std::map<int, std::unordered_set<int64_seq>> int64_lvdist(const int64_seq &query, const std::vector<barcode_entry> &targets, int max_dist = 4) {
        std::map<int, std::unordered_set<int64_seq>> results;
        if (query.bits.empty()){
            return results;
        }
        for (const auto &target : targets) {
            if (target.barcode.bits.empty()){
                continue;
            }
            
            int dist;
            if (query.length <= target.barcode.length) {
                // Query is shorter or equal, search query in target
                dist = bit_partial_match(query.bits[0], target.barcode.bits[0], query.length, target.barcode.length, max_dist);
            } else {
                // Target is shorter, search target in query
                dist = bit_partial_match(target.barcode.bits[0], query.bits[0], target.barcode.length, query.length, max_dist);
            }
            
            if(dist >= 0){
                results[dist].insert(target.barcode);
            }
        }
        return results;
    }

/**
 * @brief Calculate Levenshtein distance between two `int64_seq` sequences (partial matching)
 * * @param query The query sequence as an `int64_seq`
 * @param target The target sequence as an `int64_seq`
 * @param max_dist Maximum allowed edit distance
 * @return `int` minimum edit distance if <= max_dist, else -1
 */
    int int64_lvdist(const int64_seq &query, const int64_seq &target, int max_dist) {
        if (query.bits.empty() || target.bits.empty()){
            return -1;
        }
        
        if (query.length <= target.length) {
            // Query is shorter or equal, search query in target
            return bit_partial_match(query.bits[0], target.bits[0], query.length, target.length, max_dist);
        } else {
            // Target is shorter, search target in query
            return bit_partial_match(target.bits[0], query.bits[0], target.length, query.length, max_dist);
        }
    }

/**
 * @brief Generate all single-point mutations for a given `int64_seq`. 
 * @brief Works by iterating over each position in the sequence and substituting each possible nucleotide except the original and uses bitmasks to efficiently create mutated sequences.
 * @param seq The original sequence as an `int64_seq`
 * @return `unordered_set<int64_seq>` containing all unique single-point mutations of the input
 */
    std::unordered_set<int64_seq> generate_point_mutations(const int64_seq &seq) {
        std::unordered_set<int64_seq> mutations;
        if (seq.bits.empty() || seq.length == 0 || seq.length > 32){
            return mutations;
        }
        int sequence_length = seq.length;
        uint64_t original_value = seq.bits[0];
        for (int pos = 0; pos < sequence_length; ++pos) {
            uint64_t original_nuc = (original_value >> (2 * pos)) & uint64_t{3};
            for (uint64_t new_nuc = 0; new_nuc < 4; ++new_nuc) {
                if (new_nuc != original_nuc) {
                    uint64_t clear_mask = ~(uint64_t{3} << (2 * pos));
                    uint64_t mutated_value =
                        (original_value & clear_mask) |
                        (new_nuc << (2 * pos));
                    int64_seq candidate;
                    candidate.length = sequence_length;
                    candidate.bits.push_back(mutated_value);
                    mutations.insert(candidate);
                }
            }
        }
        return mutations;
    } 
/**
     * @brief Generate all single-nucleotide point mutations of a 2-bit encoded DNA sequence
     * @param original_value `int64_t` encoded DNA sequence (2 bits per nucleotide)
     * @param sequence_length Number of nucleotides in the sequence (max 32)
     * @return `unordered_set<int64_t>` of all possible single point mutations (3 × sequence_length mutations)
     * 
     * This function creates all possible single point mutations of a DNA sequence that is
     * encoded as a 64-bit integer where each nucleotide occupies 2 bits (allowing sequences
     * up to 32bp in length).
     * 
     * Bit encoding (per nucleotide):
     * - Each nucleotide is represented by 2 bits
     * - Position 0 is stored in bits [1:0], position 1 in bits [3:2], etc.
     * - Nucleotides are encoded as: A=00, C=01, T=10, G=11
     * 
     * Mutation generation process:
     * 1. For each position in the sequence (0 to sequence_length-1):
     *    - Extract the original 2-bit nucleotide at that position
     *    - Generate 3 mutations by substituting with each of the other 3 nucleotides
     * 2. Bit manipulation for each mutation:
     *    - Create a clear_mask: Inverts 0b11 shifted to target position, zeroing those 2 bits
     *    - Apply mask: (original_value & clear_mask) clears the target position
     *    - Insert new nucleotide: OR with (new_nuc << (2 * pos)) to set the new value
     * 
     * Example for a 3bp sequence (6 bits total):
     * Original: ATG = 00|10|11 (binary) = 0b001011 = 11 (decimal)
     * Position 1 (T=10): Generate mutations with A(00), C(01), G(11)
     * - Mask: ~(0b11 << 2) = ~0b1100 = ...11110011
     * - Clear: 0b001011 & 0b11110011 = 0b000011
     * - Mutate to A: 0b000011 | (0b00 << 2) = 0b000011 = AAG
     * - Mutate to C: 0b000011 | (0b01 << 2) = 0b000111 = ACG  
     * - Mutate to G: 0b000011 | (0b11 << 2) = 0b001111 = AGG
 */
    std::unordered_set<uint64_t> generate_point_mutations_raw(uint64_t original_value, int sequence_length) {
        std::unordered_set<uint64_t> mutations;
        if (sequence_length <= 0 || sequence_length > 32) {
            return mutations;
        }
        mutations.reserve(sequence_length * 3); // Pre-reserve for efficiency
        
        for (int pos = 0; pos < sequence_length; ++pos) {
            uint64_t original_nuc =
                (original_value >> (2 * pos)) & uint64_t{3};
            for (uint64_t new_nuc = 0; new_nuc < 4; ++new_nuc) {
                if (new_nuc != original_nuc) {
                    uint64_t clear_mask = ~(uint64_t{3} << (2 * pos));
                    uint64_t mutated_value =
                        (original_value & clear_mask) |
                        (new_nuc << (2 * pos));
                    mutations.insert(mutated_value); // Fast int64_t hash and insert
                }
            }
        }
        return mutations;
    }

/**
 * @brief Generate all possible mutated barcodes within a specified Hamming
 * distance.
 *
 * Chooses each set of mutated positions exactly once, in increasing position
 * order, and substitutes one of the three bases that differs from the
 * original at each chosen position. This produces the same union of Hamming
 * shells as round-by-round mutation, without generating and de-duplicating
 * the same distance-two-or-greater barcode through every mutation ordering.
 *
 * This routine operates on the first (and only) 64-bit chunk, so it supports
 * barcodes up to 32 bases.
 * @param seq The original barcode sequence as an `int64_seq`
 * @param mutation_rounds Number of mutation rounds (default: 1)
 * @return `unordered_set<int64_seq>` containing all unique mutated barcodes
 */
    std::unordered_set<int64_seq> generate_mutated_barcodes(const int64_seq &seq, int mutation_rounds = 1) {
        std::unordered_set<int64_seq> final_mutations;
        if (seq.bits.empty() || mutation_rounds <= 0 ||
            seq.length == 0 || seq.length > 32) {
            return final_mutations;
        }

        const int sequence_length = seq.length;
        const int max_distance = std::min(mutation_rounds, sequence_length);
        const uint64_t initial_raw = static_cast<uint64_t>(seq.bits[0]);

        // Exact through distance three, which covers the production paths.
        // For larger distances this remains a useful lower-bound reservation.
        size_t expected_total = static_cast<size_t>(sequence_length) * 3;
        if (max_distance >= 2) {
            expected_total +=
                static_cast<size_t>(sequence_length) *
                static_cast<size_t>(sequence_length - 1) * 9 / 2;
        }
        if (max_distance >= 3) {
            expected_total +=
                static_cast<size_t>(sequence_length) *
                static_cast<size_t>(sequence_length - 1) *
                static_cast<size_t>(sequence_length - 2) * 27 / 6;
        }
        final_mutations.reserve(expected_total);

        const auto enumerate_distance =
            [&](const auto& self, int next_position, int remaining,
                uint64_t raw_bits) -> void {
                if (remaining == 0) {
                    final_mutations.emplace(
                        raw_bits, static_cast<uint16_t>(sequence_length));
                    return;
                }

                // Leave enough positions after this one to finish the
                // requested Hamming shell. Increasing positions make every
                // set of substitutions unique by construction.
                const int last_position = sequence_length - remaining;
                for (int position = next_position;
                     position <= last_position; ++position) {
                    const int shift = 2 * position;
                    const uint64_t position_mask = uint64_t{3} << shift;
                    const uint64_t original_nucleotide =
                        (initial_raw & position_mask) >> shift;

                    for (uint64_t nucleotide = 0; nucleotide < 4;
                         ++nucleotide) {
                        if (nucleotide == original_nucleotide) {
                            continue;
                        }
                        const uint64_t mutated_bits =
                            (raw_bits & ~position_mask) |
                            (nucleotide << shift);
                        self(self, position + 1, remaining - 1,
                             mutated_bits);
                    }
                }
            };

        for (int distance = 1; distance <= max_distance; ++distance) {
            enumerate_distance(
                enumerate_distance, 0, distance, initial_raw);
        }

        return final_mutations;
    }

/**
 * @brief Generate all possible barcodes within a specified number of edit rounds (substitutions, insertions, deletions).
 * @param seq The original barcode sequence as an `int64_seq`
 * @param edit_rounds Number of edit rounds (default: 1)
 * @return `unordered_set<int64_seq>` containing all unique barcodes within the specified edit distance (excluding the original)
 * @brief Generate all DNA sequences within a specified edit distance using substitutions, deletions, and insertions
 * 
 * This function generates all possible DNA sequences that can be reached from the input sequence within the specified number of edit
 * operations. Each round explores sequences one additional edit away from the original.
*/
    std::unordered_set<int64_seq> generate_lv_barcodes(const int64_seq &seq, int edit_rounds = 1) {
        if (seq.bits.empty() || seq.length == 0 || seq.length > 32 ||
            edit_rounds <= 0) return {};
        
        std::unordered_set<int64_seq> all_mutations;
        std::unordered_set<int64_seq> current_round;
        current_round.insert(seq);
        all_mutations.insert(seq);
        
        for (int round = 0; round < edit_rounds; ++round) {
            std::unordered_set<int64_seq> next_round;
            
            for (const auto& candidate : current_round) {
                uint64_t bits = candidate.bits[0];
                int length = candidate.length;
                
                // SUBSTITUTIONS
                for (int pos = 0; pos < length; ++pos) {
                    uint64_t original_nuc =
                        (bits >> (2 * pos)) & uint64_t{3};
                    for (uint64_t new_nuc = 0; new_nuc < 4; ++new_nuc) {
                        if (new_nuc != original_nuc) {
                            uint64_t clear_mask =
                                ~(uint64_t{3} << (2 * pos));
                            uint64_t new_bits =
                                (bits & clear_mask) |
                                (new_nuc << (2 * pos));
                            int64_seq new_seq(
                                new_bits, static_cast<uint16_t>(length));
                            
                            if (all_mutations.insert(new_seq).second) {
                                next_round.insert(new_seq);
                            }
                        }
                    }
                }
                
                // DELETIONS
                if (length > 1) {
                    for (int del_pos = 0; del_pos < length; ++del_pos) {
                        uint64_t new_bits = 0;
                        int new_pos = 0;
                        
                        for (int i = 0; i < length; ++i) {
                            if (i != del_pos) {
                                uint64_t nuc =
                                    (bits >> (2 * i)) & uint64_t{3};
                                new_bits |= (nuc << (2 * new_pos));
                                new_pos++;
                            }
                        }
                        
                        int64_seq new_seq(
                            new_bits, static_cast<uint16_t>(length - 1));
                        if (all_mutations.insert(new_seq).second) {
                            next_round.insert(new_seq);
                        }
                    }
                }
                
                // INSERTIONS
                if (length < 32) {
                    for (int ins_pos = 0; ins_pos <= length; ++ins_pos) {
                        for (uint64_t new_nuc = 0; new_nuc < 4; ++new_nuc) {
                            uint64_t new_bits = 0;
                            int new_pos = 0;
                            
                            // Copy before insertion
                            for (int i = 0; i < ins_pos; ++i) {
                                uint64_t nuc =
                                    (bits >> (2 * i)) & uint64_t{3};
                                new_bits |= (nuc << (2 * new_pos));
                                new_pos++;
                            }
                            
                            // Insert new nucleotide
                            new_bits |= (new_nuc << (2 * new_pos));
                            new_pos++;
                            
                            // Copy after insertion
                            for (int i = ins_pos; i < length; ++i) {
                                uint64_t nuc =
                                    (bits >> (2 * i)) & uint64_t{3};
                                new_bits |= (nuc << (2 * new_pos));
                                new_pos++;
                            }
                            
                            int64_seq new_seq(
                                new_bits, static_cast<uint16_t>(length + 1));
                            if (all_mutations.insert(new_seq).second) {
                                next_round.insert(new_seq);
                            }
                        }
                    }
                }
            }
            
            if (next_round.empty()) break;
            current_round = std::move(next_round);
        }
        
        // Remove original sequence
        all_mutations.erase(seq);
        return all_mutations;
    }

/**
 * @brief Detect homopolymers at the start or end of a DNA sequence represented as an `int64_seq`
 * @param sequence The DNA sequence as a `string`
 * @param max_hp Maximum allowed homopolymer length
 * @return `bool` true if a homopolymer of length >= max_hp is detected, false otherwise
 */
    bool detect_hp(const std::string& sequence, int max_hp) {
        if (sequence.empty() || max_hp <= 0) {
            return false;
        }
        
        int seq_len = static_cast<int>(sequence.length());
        
        // Check leading homopolymer
        if (seq_len >= max_hp) {
            char first_base = std::toupper(sequence[0]);
            int leading_count = 1;
            
            for (int i = 1; i < seq_len && i < max_hp; ++i) {
                if (std::toupper(sequence[i]) == first_base) {
                    leading_count++;
                } else {
                    break;
                }
            }
            
            if (leading_count >= max_hp) {
                return true;
            }
        }
        
        // Check trailing homopolymer
        if (seq_len >= max_hp) {
            char last_base = std::toupper(sequence[seq_len - 1]);
            int trailing_count = 1;
            
            for (int i = seq_len - 2; i >= 0 && i >= seq_len - max_hp; --i) {
                if (std::toupper(sequence[i]) == last_base) {
                    trailing_count++;
                } else {
                    break;
                }
            }
            
            if (trailing_count >= max_hp) {
                return true;
            }
        }
        
        return false;
    }

/**
 * @brief Detect homopolymers at the start or end of a DNA sequence represented as an `int64_seq`
 * @param sequence The DNA sequence as an `int64_seq`
 * @param max_hp Maximum allowed homopolymer length
 * @return `bool` true if a homopolymer of length >= max_hp is detected, false otherwise
 */
    bool detect_hp(const int64_seq& sequence, int max_hp) {
        if (!sequence.is_valid() || max_hp <= 0) {
            return false;
        }
        
        std::string seq_str = sequence.bits_to_sequence();
        return detect_hp(seq_str, max_hp);
    }
};

   // === bc_flat_set: lean index for huge whitelists =======================
/**
* @brief memory-efficient barcode whitelist using parallel flat hash set
* @param key: `int64_seq` representing barcode sequences
* @param value: `barcode_entry` containing barcode and associated data
*/
    struct bc_flat_set {
        using key = int64_seq;
        using value = barcode_entry;

        phmap::parallel_flat_hash_set<value, phmap::Hash<value>, phmap::EqualTo<value>> s;

        bc_flat_set() = default;
    private:
        static inline const key& key_of(const key& k) noexcept { 
            return k; 
        }
        static inline const key& key_of(const value& v) noexcept { 
            return v.barcode; 
        }

    /**
     * @brief Extract the key (barcode) from either a key or value type, handling pointers as well
     * @tparam U Type which can be either key, value, or pointer to either
     * @param u The input object from which to extract the key
     * @return const reference to the extracted key
     */
    template <typename U>
    static inline const key& key_of_any(const U& u) {
        if constexpr (std::is_pointer_v<std::decay_t<U>>) {
            return key_of(*u);
        } else {
            return key_of(u);
        }
    }

    /**
     * @brief Trait to check if a type is a range (has begin() and end() methods)
     * @tparam T Type to check
     */
    template <typename T, typename = void>
    struct is_range : std::false_type {};
    template <typename T>
    struct is_range<T,
    std::void_t<
    decltype(std::begin(std::declval<const T&>())),
    decltype(std::end(std::declval<const T&>()))
    >> : std::true_type {};

    template <typename T>
    static constexpr bool is_key_or_value_v =
        std::is_same_v<std::decay_t<T>, key> ||
        std::is_same_v<std::decay_t<T>, value> ||
        (std::is_pointer_v<std::decay_t<T>> &&
        (std::is_same_v<std::remove_pointer_t<std::decay_t<T>>, key> ||
        std::is_same_v<std::remove_pointer_t<std::decay_t<T>>, value>));

    bool contains_key(const key& k) const {
        if (!k.is_valid()) return false;
        value probe; probe.barcode = k;
        return s.find(probe) != s.end();
    }

    public:
        //  Size & housekeeping (names to minimize churn)
        size_t size()             const { return s.size(); }
        size_t association_size() const { return s.size(); }
        size_t unique_val_size()  const { return s.size(); }
        bool   empty()            const { return s.empty(); }
        void   clear()                  { s.clear(); }

/**
 * @brief Check if any barcode in the input matches an entry in the whitelist
 * @tparam T Type which can be a key, value, or a range of either
 * @param x The input object to check for matches
 * @return `bool` true if any barcode matches, false otherwise
 */
        template<typename T> bool check_wl_for(const T& x) const {
            if constexpr (is_key_or_value_v<T>) {
                const key& k = key_of_any(x);
                return contains_key(k);
            } else if constexpr (is_range<T>::value) {
                for (const auto& y : x) {
                    const key& ky = key_of_any(y);
                    if (!ky.is_valid()) continue;
                    if (contains_key(ky)) return true;
                }
                return false;
            }
            return false;
        }

/**
 * @brief Return all barcodes from the input that match entries in the whitelist
 * @tparam T Type which can be a key, value, or a range of either
 * @param x The input object to check for matches
 * @return `unordered_set<key>` of matching barcodes
 */
        template<typename T> std::unordered_set<key> return_matching_barcodes(const T& x) const {
            std::unordered_set<key> out;
            if constexpr (is_key_or_value_v<T>) {
                const key& k = key_of_any(x);
                if (contains_key(k)) out.insert(k);
            } else if constexpr (is_range<T>::value) {
                for (const auto& y : x) {
                    const key& ky = key_of_any(y);
                    if (!ky.is_valid()) continue;
                    if (contains_key(ky)) out.insert(ky);
                }
            } 
            return out;
        }
/**
 * @brief Return all putative correct barcodes from the input that match entries in the whitelist
 * @tparam T Type which can be a key, value, or a range of either
 * @param x The input object to check for matches
 * @return `unordered_set<key>` of putative correct barcodes
 */
        template<typename T> std::unordered_set<key> return_putative_correct_bcs(const T& x) const {
            return return_matching_barcodes(x);
        }

/**
 * @brief Insert a barcode entry into the whitelist
 * @tparam T Type which can be a key, value, or a range of either
 * @param observed The observed barcode (not used in this implementation)
 * @param correct The correct barcode entry to insert
 */
        template<typename T> void insert_bc_entry(const T& /*observed*/, const value& correct) {
            if (!correct.is_valid()) return;
            s.lazy_emplace_l(correct,
                [&](auto& it) { /* present: no-op */ },
                [&](const auto& ctor){ ctor(correct); });
        }
/**
 * @brief Insert a barcode entry into the whitelist by constructing it from the observed barcode
 * @tparam T Type which can be a key, value, or a range of either
 * @param observed The observed barcode from which to construct the entry
 */
        template<typename T> void insert_bc_entry(const T& /*observed*/, value&& correct) {
            if (!correct.is_valid()) return;
            s.lazy_emplace_l(correct,
                [&](auto& it) { /* present */ },
                [&](const auto& ctor){ ctor(std::move(correct)); });
        }
/**
 * @brief Insert a barcode entry into the whitelist using only the observed barcode
 * @tparam T Type which can be a key, value, or a range of either
 * @param observed The observed barcode from which to create the entry
 */
        template<typename T> void insert_bc_entry(const T& observed) {
            value v{}; v.barcode = key_of(observed);
            if (!v.is_valid()) return;
            insert_bc_entry(observed, std::move(v));
        }
/**
 * @brief Remove a barcode entry from the whitelist by observed barcode
 * @tparam T Type which can be a key, value, or a range of either
 * @param observed The observed barcode to remove
 */
        // Remove by key 
        template<typename T> void remove_bc_entry(const T& observed) {
            value probe; probe.barcode = key_of(observed);
            s.erase(probe);
        }
/**
 * @brief Erase a barcode entry from the whitelist by observed barcode and return the number of elements removed
 * @tparam T Type which can be a key, value, or a range of either
 * @param observed The observed barcode to erase
 */
        template<typename T> size_t erase(const T& observed) {
            value probe; probe.barcode = key_of(observed);
            return s.erase(probe);
        }

/**
 * @brief Get the total barcode count for a given barcode
 * @tparam T Type which can be a key, value, or a range of either
 * @param x The input object for which to get the barcode count
 * @return `counter` object containing all barcode counts
 */
        template<typename T> counter get_all_bc_counts(const T& x) const {
            value probe; probe.barcode = key_of(x);
            if (!probe.is_valid()) return counter(0);
            auto it = s.find(probe);
            return (it!=s.end()) ? it->count : counter(0);
        }

/**
 * @brief Update the barcode count for a given barcode
 * @tparam T Type which can be a key, value, or a range of either
 * @param x The input object for which to update the barcode count
 * @param slot The specific barcode count slot to update (default is total)
 */
        template<typename T> void update_bc_count(const T& x, barcode_counts slot = barcode_counts::total) const {
            value probe; probe.barcode = key_of(x);
            if (!probe.is_valid()) return;
            auto it = s.find(probe);
            if (it != s.end()) const_cast<value&>(*it).count.increment(slot);
        }

/**
 * @brief Single-lookup entry accessor for direct counter manipulation
 * @tparam T key, value, or int64_seq
 * @param x the barcode to look up
 * @return mutable pointer to the entry, or nullptr if not found
 */
        template<typename T> value* find_entry(const T& x) const {
            value probe; probe.barcode = key_of(x);
            if (!probe.is_valid()) return nullptr;
            auto it = s.find(probe);
            if (it == s.end()) return nullptr;
            return const_cast<value*>(&(*it));
        }
/**
 * @brief Subtract from the barcode count for a given barcode
 * @tparam T Type which can be a key, value, or a range of either
 * @param x The input object for which to subtract from the barcode count
 * @param slot The specific barcode count slot to update (default is total)
 */
        template<typename T> void subtract_bc_count(const T& x, barcode_counts slot = barcode_counts::total) const {
            value probe; probe.barcode = key_of(x);
            if (!probe.is_valid()) return;
            auto it = s.find(probe);
            if (it != s.end()) const_cast<value&>(*it).count.subtract(slot);
        }
/**
 * @brief Get the barcode count for a given barcode
 * @tparam T Type which can be a key, value, or a range of either
 * @param x The input object for which to get the barcode count
 * @param slot The specific barcode count slot to retrieve (default is total)
 * @return `int` barcode count for the specified slot
 */
        template<typename T> int get_bc_count(const T& x, barcode_counts slot = barcode_counts::total) const {
            value probe; probe.barcode = key_of(x);
            if (!probe.is_valid()) return 0;
            auto it = s.find(probe);
            return (it != s.end()) ? it->count.load(slot) : 0;
        }
/**
 * @brief Set the barcode count for a given barcode
 * @tparam T Type which can be a key, value, or a range of either
 * @param x The input object for which to set the barcode count
 * @param new_count The new count value to set
 * @param slot The specific barcode count slot to update (default is total)
 */
        template<typename T> void set_bc_count(const T& x, int new_count, barcode_counts slot = barcode_counts::total) const {
            value probe; probe.barcode = key_of(x);
            if (!probe.is_valid()) return;
            auto it = s.find(probe);
            if (it != s.end()) const_cast<value&>(*it).count.counts[slot].store(new_count, std::memory_order_relaxed);
        }

/**
 * @brief Get all unique barcode entries in the whitelist
 * @return `vector<const value*>` of pointers to all unique barcode entries
 */
        std::vector<const value*> get_unique_entries() const {
            std::vector<const value*> result;
            result.reserve(s.size());
            for (const auto& val : s) {
                result.push_back(&val);
            }
            return result;
        }

        // ---- Equal-range compatibility (return 0/1 element range) ----
        struct value_iterator {
            const value* p = nullptr;
            bool end_ = true;
            value_iterator() = default;
            explicit value_iterator(const value* pv, bool e=false) : p(pv), end_(e) {}
            const value& operator*()  const { return *p; }
            const value* operator->() const { return p; }
            value_iterator& operator++(){ end_=true; return *this; }
            bool operator!=(const value_iterator& o) const { return end_!=o.end_ || p!=o.p; }
            bool operator==(const value_iterator& o) const { return !(*this!=o); }
        };
/**
 * @brief Get the equal range for a given barcode key
 * @param k The barcode key to search for
 * @return `pair<value_iterator,value_iterator>` representing the range of matching entries (0 or 1 element)
 */
        std::pair<value_iterator,value_iterator> equal_range(const key& k) const {
            value probe; probe.barcode = k;
            auto it = s.find(probe);
            if (it == s.end()) return { value_iterator(nullptr,true), value_iterator(nullptr,true) };
            return { value_iterator(&*it,false), value_iterator(nullptr,true) };
        }

/**
 * @brief Summarize barcode counts for all entries or a filtered set of barcodes
 * @param filter_keys Optional pointer to a set of barcode keys to filter the summary
 * @return `vector<tuple<string,int,int,int,int,int,int,int,int>>` containing barcode summaries
 */
        std::vector<std::tuple<std::string,int,int,int,int,int,int,int,int>>
        summarize_counts(const std::unordered_set<key>* filter_keys=nullptr) const {
            using row = std::tuple<std::string,int,int,int,int,int,int,int,int>;
            std::vector<row> rows;
            if (filter_keys) {
                rows.reserve(filter_keys->size());
                for (auto const& k : *filter_keys) {
                    value probe; probe.barcode = k;
                    auto it = s.find(probe);
                    int raw=0, tot=0, corr=0, fwd=0, rev=0, fwd_c=0, rev_c=0, filt=0;
                    if (it!=s.end()) {
                        auto const& c = it->count;
                        raw   = c.load(barcode_counts::raw);
                        tot   = c.load(barcode_counts::total);
                        corr  = c.load(barcode_counts::corrected);
                        fwd   = c.load(barcode_counts::forw);
                        rev   = c.load(barcode_counts::rev);
                        fwd_c = c.load(barcode_counts::forw_concat);
                        rev_c = c.load(barcode_counts::rev_concat);
                        filt  = c.load(barcode_counts::filtered);
                    }
                    rows.emplace_back(k.bits_to_sequence(),raw,tot,corr,fwd,rev,fwd_c,rev_c,filt);
                }
            } else {
                rows.reserve(s.size());
                for (auto const& be : s) {
                    auto const& k = be.barcode;
                    auto const& c = be.count;
                    rows.emplace_back(
                        k.bits_to_sequence(),
                        c.load(barcode_counts::raw),
                        c.load(barcode_counts::total),
                        c.load(barcode_counts::corrected),
                        c.load(barcode_counts::forw),
                        c.load(barcode_counts::rev),
                        c.load(barcode_counts::forw_concat),
                        c.load(barcode_counts::rev_concat),
                        c.load(barcode_counts::filtered)
                    );
                }
            }
            return rows;
        }

        // ---- Minimal debug hooks so memory report still compiles ----
        const auto& debug_unique_values() const noexcept { return s; }
        struct fake_assoc_view {
            struct node { key k; std::vector<const value*> vals; };
            std::vector<node> rows;
            size_t bucket_count() const { return rows.size(); }
            size_t subcnt()       const { return 1; }
            struct H { size_t operator()(const key& k) const {
                return std::hash<uint64_t>{}(k.bits.empty()?0:k.bits[0]); } };
            auto hash_function() const { return H{}; }
            auto begin() const { return rows.begin(); }
            auto end()   const { return rows.end(); }
        };
        fake_assoc_view debug_associations() const noexcept {
            fake_assoc_view v;
            v.rows.reserve(s.size());
            for (auto const& be : s) {
                typename fake_assoc_view::node n;
                n.k = be.barcode;
                n.vals.push_back(&be);
                v.rows.push_back(std::move(n));
            }
            return v;
        }
        size_t validate_association_keys()   const noexcept { return 0; }
        size_t validate_association_values() const noexcept { return 0; }

        template<typename K, typename V>
        void emplace(K&& k, V&& v) { insert_bc_entry(std::forward<K>(k), std::forward<V>(v)); }

        auto begin() const { return s.begin(); }
        auto end()   const { return s.end();   }

/**
 * @brief Write a summary of barcode counts to an output stream
 * @param out The output stream to write to
 * @param class_id A string identifier for the barcode class
 * @param filter_keys Optional pointer to a set of barcode keys to filter the summary
 * @param write_header Whether to write the header line (default: true)
 */
        void write_wl_summary(
            std::ostream &out,
            const std::string class_id,
            const std::unordered_set<key> *filter_keys = nullptr,
            bool write_header = true
        ) const {
            static bool header_printed = false;
            if (write_header && header_printed) header_printed = false;
            if (!header_printed) {
                out << "class_id,true_barcode,raw_count,total_count,corrected_count,forw_count,rev_count,"
                    "forw_concat_count,rev_concat_count,filtered_count\n";
                header_printed = true;
            }
            auto rows = summarize_counts(filter_keys);
            for (auto const & tpl : rows) {
                out << class_id << ','
                    << std::get<0>(tpl) << ','
                    << std::get<1>(tpl) << ','
                    << std::get<2>(tpl) << ','
                    << std::get<3>(tpl) << ','
                    << std::get<4>(tpl) << ','
                    << std::get<5>(tpl) << ','
                    << std::get<6>(tpl) << ','
                    << std::get<7>(tpl) << ','
                    << std::get<8>(tpl) << '\n';
            }
        }
    };

/**
* @brief concurrent multimap for barcode corrections
* @param key: `int64_seq` representing barcode sequences
* @param value: `barcode_entry` containing barcode and associated data
*
* @brief A thread-safe multimap structure for managing barcode corrections,
* allowing multiple values to be associated with a single key. Basically, there's a 
* stable storage of unique barcode entries, and then a concurrent map that maps
* keys to sets of pointers to those unique entries. 
*/
template<typename key, typename value> struct bc_multimap {
private:
    static inline const key& key_of(const key &k) noexcept { 
        return k; 
    }
    static inline const key& key_of(const value &v) noexcept { 
        return v.barcode; 
    }

public:
    // Stable value storage where barcode entries are stored singularly
    phmap::parallel_node_hash_set<value> unique_values;

    // pointer map of associations, split up into 2^number_of_shards in there to allow for concurrent read-writes
    //if there's read-write, each submap gets its own mutex to avoid rehash failures
    using pointer_map = phmap::parallel_node_hash_map<
    key,
    phmap::flat_hash_set<const value*>, 
    phmap::Hash<key>, 
    phmap::EqualTo<key>, 
    std::allocator<std::pair<const key, phmap::flat_hash_set<const value*>>>, 
    6,
    std::mutex
    >; 

    pointer_map associations;

    bc_multimap() {
        associations.reserve(500000);
    }

    size_t size() const { 
        return associations.size(); 
    }

    size_t association_size() const { 
        return associations.size(); 
    }

    size_t unique_val_size() const { 
        return unique_values.size(); 
    }

    bool empty() const { 
        return associations.empty(); 
    }

    void clear() {
        associations.clear();
        unique_values.clear();
    }

    void clear_associations() {
        associations.clear();
    }

    const auto& debug_unique_values() const noexcept { 
        return unique_values; 
    }

    const auto& debug_associations() const noexcept { 
        return associations; 
    }

    size_t validate_association_keys() const noexcept {
       size_t corrupted = 0;
       try {
           for (const auto& [k, value_set] : associations) {
               try {
                   // Quick access to key data
                   volatile auto len = k.length;
                   volatile auto bits_empty = k.bits.empty();
               } catch (...) {
                   corrupted++;
               }
           }
       } catch (...) {
           return SIZE_MAX; // Can't iterate at all
       }
       return corrupted;
   }

    size_t validate_association_values() const noexcept {
       size_t nulls = 0;
       try {
           for (const auto& [k, value_set] : associations) {
               for (const value* ptr : value_set) {
                   if (ptr == nullptr) {
                       nulls++;
                   }
               }
           }
       } catch (...) {
           return SIZE_MAX; // Can't iterate at all
       }
       return nulls;
   }

    // ==== Core Insertion Methods ====
    template<typename T>
    void insert_bc_entry(const T &observed, const value &correct) {
        const key &k = key_of(correct);
        if (!k.is_valid()) return;
        
        // validate observed key as well:
        const key& obs_key = key_of(observed);
        if (!obs_key.is_valid() || obs_key.bits.empty()) {
            return;
        }
        
        // add to unique_values (thread-safe)
        auto [it, inserted] = unique_values.emplace(correct);
        const value* ptr = &(*it);
        
        // Then add to associations
        associations.lazy_emplace_l(obs_key,
            [&](auto& kv_pair) { 
                kv_pair.second.insert(ptr);
            },
            [&](const auto& ctor) { 
                phmap::flat_hash_set<const value*> new_set;
                new_set.insert(ptr);
                ctor(obs_key, std::move(new_set));
            }
        );
    }

    template<typename T>
    void insert_bc_entry(const T &observed, value &&correct) {
        const key &k = key_of(correct);
        //validate key:
        if (!k.is_valid()){
            return;
        }
        // validate observed:
        const key& obs_key = key_of(observed);
        if (!obs_key.is_valid() || obs_key.bits.empty()) {
            return;
        }
        
        // Add to unique_values (thread-safe)
        auto [it, inserted] = unique_values.emplace(std::move(correct));
        const value* ptr = &(*it);
        
        // Add to associations
        associations.lazy_emplace_l(obs_key,
            [&](auto& kv_pair) {
                kv_pair.second.insert(ptr);
            },
            [&](const auto& ctor) {
                phmap::flat_hash_set<const value*> new_set;
                new_set.insert(ptr);
                ctor(obs_key, std::move(new_set));
            }
        );
    }

    template<typename T>
    void insert_bc_entry(const T &observed) {
        value v{};
        v.barcode = key_of(observed);
        if (!v.barcode.is_valid()){
            return;
        }
        insert_bc_entry(observed, std::move(v));
    }

    template<typename T, typename... Args>
    void emplace_bc_entry(const T &observed, Args&&... args) {
        value v{std::forward<Args>(args)...};
        insert_bc_entry(observed, std::move(v));
    }

    //For multimap-like interface compatibility
    template<typename K, typename V>
    void emplace(K&& k, V&& v) {
        insert_bc_entry(std::forward<K>(k), std::forward<V>(v));
    }

    // Multimap-like insert
    template<typename K, typename V>
    void insert(const std::pair<K, V>& kv) {
        insert_bc_entry(kv.first, kv.second);
    }

    // ==== Removal Methods ====
    template<typename T>
    void remove_bc_entry(const T &observed) {
        const key& obs_key = key_of(observed);
        associations.erase(obs_key);
    }

    // Multimap-like erase
    template<typename T>
    size_t erase(const T &observed) {
        const key& obs_key = key_of(observed);
        return associations.erase(obs_key);
    }

    // ==== Lookup Methods ====
    template<typename T>
    bool check_wl_for(const T &x) const {
        if constexpr (std::is_same_v<T, key> || std::is_same_v<T, value>) {
            const auto& k = key_of(x);
            if (!k.is_valid()) return false;
            return associations.find(k) != associations.end();
            
        } else {
            for (const auto &y : x) {
                const auto& ky = key_of(y);
                if (!ky.is_valid()) continue;
                if (check_wl_for(y)) return true;
            }
            return false;
        }
    }

    /**
     * @brief Check whether an observed key maps to itself as a canonical value.
     *
     * Unlike check_wl_for(), this deliberately excludes mutation aliases:
     * an association for `observed` is only an identity mapping when at least
     * one associated whitelist entry has the same barcode sequence.
     */
    bool has_identity_mapping(const key& observed) const {
        if (!observed.is_valid()) return false;

        auto it = associations.find(observed);
        if (it == associations.end()) return false;

        for (const value* candidate : it->second) {
            if (candidate != nullptr && candidate->barcode == observed) {
                return true;
            }
        }
        return false;
    }

    // Multimap-like find (returns iterator to first match)
    auto find(const key& k) const {
        return associations.find(k);
    }

    // Multimap-like count
    size_t count(const key& k) const {
        auto it = associations.find(k);
        return (it != associations.end()) ? it->second.size() : 0;
    }
    
    // Iterator for value set (const_iterator wrapper)
    // Multimap-like equal_range - returns range of iterators spanning all values for key
    struct value_iterator {
        typename phmap::flat_hash_set<const value*>::const_iterator inner_it;
        typename phmap::flat_hash_set<const value*>::const_iterator inner_end;
        
        value_iterator(typename phmap::flat_hash_set<const value*>::const_iterator it,
                    typename phmap::flat_hash_set<const value*>::const_iterator end)
            : inner_it(it), inner_end(end) {}
        
        const value& operator*() const { return **inner_it; }
        const value* operator->() const { return *inner_it; }
        
        value_iterator& operator++() { 
            ++inner_it; 
            return *this; 
        }
        
        bool operator!=(const value_iterator& other) const {
            return inner_it != other.inner_it;
        }
        
        bool operator==(const value_iterator& other) const {
            return inner_it == other.inner_it;
        }
    };

    std::pair<value_iterator, value_iterator> equal_range(const key& k) const {
        auto it = associations.find(k);
        if (it != associations.end()) {
            const auto& val_set = it->second;
            return {
                value_iterator(val_set.begin(), val_set.end()),
                value_iterator(val_set.end(), val_set.end())
            };
        }

        static const phmap::flat_hash_set<const value*> empty_set;
        return {
            value_iterator(empty_set.end(), empty_set.end()),
            value_iterator(empty_set.end(), empty_set.end())
        };
    }

    template<typename T>
    std::unordered_set<key> return_matching_barcodes(const T &x) const {
        std::unordered_set<key> out;
        if constexpr (std::is_same_v<T, key> || std::is_same_v<T, value>) {
            if (check_wl_for(x)) out.insert(key_of(x));
        } else {
            for (const auto &y : x) {
                if (check_wl_for(y)) out.insert(key_of(y));
            }
        }
        return out;
    }

    template<typename T>
    std::unordered_set<key> return_putative_correct_bcs(const T &x) const {
        std::unordered_set<key> out;
        if constexpr (std::is_same_v<T, key> || std::is_same_v<T, value>) {
            const auto& k = key_of(x);
            if (!k.is_valid()){
                return out;
            }
            auto it = associations.find(k);
            if (it != associations.end()) {
                for (const value* ptr : it->second) {
                    if (ptr) {
                        out.insert(ptr->barcode);
                    }
                }
            }
        } else {
            for (const auto &y : x) {
                auto sub = return_putative_correct_bcs(y);
                out.insert(sub.begin(), sub.end());
            }
        }
        return out;
    }

    // ==== Count Management Methods ====
    template<typename T>
    counter get_all_bc_counts(T const &x) const {
        const auto& k = key_of(x);
        if (!k.is_valid()) return counter(0);
        
        auto it = associations.find(k);
        if (it != associations.end() && !it->second.empty()) {
            return (*it->second.begin())->count;  // Returns copy of counter
        }
        return counter(0);
    }

    template<typename T>
    void update_bc_count(T const &x, barcode_counts slot = total) const {
        const auto& k = key_of(x);
        if (!k.is_valid()) return;
        auto it = associations.find(k);
        if (it != associations.end()) {
            for (const value* ptr : it->second) {
                const_cast<value*>(ptr)->count.increment(slot);
            }
        }
    }

    /// Single-lookup entry accessor for direct counter manipulation.
    /// Returns a mutable pointer to the first associated entry, or nullptr.
    template<typename T>
    value* find_entry(T const &x) const {
        const auto& k = key_of(x);
        if (!k.is_valid()) return nullptr;
        auto it = associations.find(k);
        if (it != associations.end() && !it->second.empty()) {
            return const_cast<value*>(*it->second.begin());
        }
        return nullptr;
    }

    template<typename T>
    void subtract_bc_count(T const &x, barcode_counts slot = total) const {
        const auto& k = key_of(x);
        if (!k.is_valid()) return;
        auto it = associations.find(k);
        if (it != associations.end()) {
            for (const value* ptr : it->second) {
                const_cast<value*>(ptr)->count.subtract(slot);
            }
        }
    }

    template<typename T>
    int get_bc_count(T const &x, barcode_counts slot = total) const {
        const auto& k = key_of(x);
        if (!k.is_valid()) return 0;
        auto it = associations.find(k);
        if (it != associations.end() && !it->second.empty()) {
            return (*it->second.begin())->count.load(slot);
        }
        return 0;
    }

    template<typename T>
    void set_bc_count(T const &x, int new_count, barcode_counts slot = barcode_counts::total) const {
        const auto& k = key_of(x);
        if (!k.is_valid()) return;
        
        auto it = associations.find(k);
        if (it != associations.end() && !it->second.empty()) {
            (*it->second.begin())->count.counts[slot].store(new_count, std::memory_order_relaxed);
        }
    }

    // ==== utility methods ====
    std::vector<const value*> get_unique_entries() const {
        std::vector<const value*> result;
        result.reserve(unique_values.size());
        for (const auto& val : unique_values) {
            result.push_back(&val);
        }
        return result;
    }

    // ==== Iterator Support ====
    class iterator {
    public:
        using OuterIter = typename pointer_map::iterator;
        using InnerIter = typename phmap::flat_hash_set<const value*>::iterator;
        
        OuterIter outer, outer_end;
        InnerIter inner, inner_end;
        
        iterator(OuterIter o, OuterIter oe) : outer(o), outer_end(oe) {
            if (outer != outer_end) {
                inner = outer->second.begin();
                inner_end = outer->second.end();
                advance_to_valid();
            }
        }
        
        void advance_to_valid() {
            while (outer != outer_end && inner == outer->second.end()) {
                ++outer;
                if (outer != outer_end) {
                    inner = outer->second.begin();
                    inner_end = outer->second.end();
                }
            }
        }

        iterator& operator++() {
            if (inner != inner_end) ++inner;
            advance_to_valid();
            return *this;
        }

        bool operator!=(const iterator& other) const {
            return outer != other.outer || (outer != outer_end && inner != other.inner);
        }

        std::pair<const key&, const value&> operator*() const {
            return { outer->first, **inner };
        }
        
        // Add arrow operator for convenience
        struct pair_proxy {
            const key& first;
            const value& second;
            pair_proxy(const key& k, const value& v) : first(k), second(v) {}
        };
        
        pair_proxy operator->() const {
            return pair_proxy(outer->first, **inner);
        }
    };
    
    iterator begin() { 
        return iterator(associations.begin(), associations.end()); 
    }
    iterator end() { 
        return iterator(associations.end(), associations.end()); 
    }


// ==== Summary and Output Methods ====
    std::vector<std::tuple<std::string,int,int,int,int,int,int,int,int>>
    summarize_counts(
        const std::unordered_set<key> *filter_keys = nullptr
    ) const {
        using row = std::tuple<std::string,int,int,int,int,int,int,int,int>;
        std::vector<row> rows;
        std::vector<key> keys;
        
        if (filter_keys) {
            keys.reserve(filter_keys->size());
            for (auto const &k : *filter_keys)
                keys.push_back(k);
        } else {
            std::unordered_set<key> seen;
            seen.reserve(associations.size());
            for (auto const &kv : associations)
                seen.insert(kv.first);
            keys.reserve(seen.size());
            for (auto const &k : seen)
                keys.push_back(k);
        }

        rows.reserve(keys.size());
        for (auto const &bc : keys) {
            auto it = associations.find(bc);
            int raw=0, tot=0, corr=0, fwd=0, rev=0, fwd_c=0, rev_c=0, filt=0;

            if (it != associations.end() && !it->second.empty()) {
                // Get counts from first value pointer (they should all be the same for a given key)
                const value* first_val = *it->second.begin();
                raw = first_val->count.load(barcode_counts::raw);
                tot = first_val->count.load(barcode_counts::total);
                corr = first_val->count.load(barcode_counts::corrected);
                fwd = first_val->count.load(barcode_counts::forw);
                rev = first_val->count.load(barcode_counts::rev);
                fwd_c = first_val->count.load(barcode_counts::forw_concat);
                rev_c = first_val->count.load(barcode_counts::rev_concat);
                filt = first_val->count.load(barcode_counts::filtered);
            }
            rows.emplace_back(bc.bits_to_sequence(), raw, tot, corr, fwd, rev, fwd_c, rev_c, filt);
        }
        return rows;
    }

    void write_wl_summary(
        std::ostream &out, 
        const std::string class_id, 
        const std::unordered_set<key> *filter_keys = nullptr, 
        bool write_header = true
    ) const {
        static bool header_printed = false;
        if(write_header && header_printed) header_printed = false;
        if(!header_printed) {
            out << "class_id,true_barcode,raw_count,total_count,corrected_count,forw_count,rev_count,"
                << "forw_concat_count,rev_concat_count,filtered_count\n";
            header_printed = true;
        }
        auto rows = summarize_counts(filter_keys);
        for (auto const & tpl : rows) {
            auto const & bc = std::get<0>(tpl);
            auto raw = std::get<1>(tpl);
            auto tot  = std::get<2>(tpl);
            auto corr = std::get<3>(tpl);
            auto forw = std::get<4>(tpl);
            auto rev  = std::get<5>(tpl);
            auto forw_concat = std::get<6>(tpl);
            auto rev_concat = std::get<7>(tpl);
            auto filtered = std::get<8>(tpl);
            out << class_id << ',' << bc << ','  << raw << "," << tot << ',' << corr << ',' << forw 
                << ',' << rev << ',' << forw_concat << ',' << rev_concat << ',' << filtered << '\n';
        }
    }

    void write_wl_summary(
        const std::string &path = "", 
        const std::unordered_set<key> *filter_keys = nullptr
    ) const {
        std::ofstream ofs(path);
        if (!ofs) throw std::runtime_error("Failed to open output file: " + path);
        write_wl_summary(ofs, "", filter_keys);
    }

};

/**
 * @brief Sparse bitmask for validating joint barcode pairs
 *
 * Stores a dense M x N bit matrix indexed by barcode whitelist position.
 * For Visium HD: 3350 x 3350 = ~1.4 MB.  Bit (i, j) is set iff the
 * (BC1_i, BC2_j) pair was observed in the valid-pairs whitelist.
 *
 * The sequence-to-index lookup maps are built once at load time from the
 * ordered whitelist files.  Validation is a single array access.
 */
struct spat_mask {
    std::vector<bool> bits;
    size_t rows = 0;   // BC1 axis size
    size_t cols = 0;   // BC2 axis size

    // sequence -> axis index
    phmap::flat_hash_map<int64_seq, uint16_t> bc1_to_idx;
    phmap::flat_hash_map<int64_seq, uint16_t> bc2_to_idx;

    bool empty() const { return bits.empty(); }

    /// Check whether a (BC1, BC2) pair is valid.
    bool check(const int64_seq& bc1, const int64_seq& bc2) const {
        auto it1 = bc1_to_idx.find(bc1);
        auto it2 = bc2_to_idx.find(bc2);
        if (it1 == bc1_to_idx.end() || it2 == bc2_to_idx.end()) return false;
        return bits[static_cast<size_t>(it1->second) * cols + it2->second];
    }

    /// Number of valid pairs (set bits).
    size_t count_valid() const {
        size_t n = 0;
        for (bool b : bits) n += b;
        return n;
    }

    /// Load a binary matrix CSV and the two axis whitelist files.
    /// The matrix has no header; row = BC1 index, col = BC2 index, values are 0 or 1.
    /// The whitelist files define the index ordering (one barcode per line, first line may be "whitelist_bcs").
    bool load(const std::string& matrix_path,
              const std::string& bc1_wl_path,
              const std::string& bc2_wl_path,
              bool verbose) {
        // Construct into local state so a malformed file never leaves a
        // partially usable mask behind.
        try {
            auto read_ordered_wl = [](const std::string& path) {
                std::vector<std::string> out;
                const auto lines = streaming_utils::import_text(path);
                bool first_content = true;
                for (const auto& raw_line : lines) {
                    std::string line = seq_utils::trim(raw_line);
                    if (line.empty()) continue;
                    const auto comma = line.find(',');
                    if (comma != std::string::npos) {
                        line = seq_utils::trim(line.substr(0, comma));
                    }
                    if (first_content && line == "whitelist_bcs") {
                        first_content = false;
                        continue;
                    }
                    first_content = false;
                    if (line.empty() || line.size() > 32 ||
                        !std::all_of(line.begin(), line.end(), [](char base) {
                            return base == 'A' || base == 'C' ||
                                   base == 'G' || base == 'T';
                        })) {
                        throw std::runtime_error(
                            "axis whitelist contains an invalid barcode");
                    }
                    out.push_back(std::move(line));
                }
                return out;
            };

            const auto bc1_list = read_ordered_wl(bc1_wl_path);
            const auto bc2_list = read_ordered_wl(bc2_wl_path);
            if (bc1_list.empty() || bc2_list.empty()) {
                throw std::runtime_error("axis whitelist is empty");
            }
            if (bc1_list.size() > std::numeric_limits<uint16_t>::max() ||
                bc2_list.size() > std::numeric_limits<uint16_t>::max()) {
                throw std::runtime_error(
                    "axis whitelist exceeds the supported 65535 entries");
            }
            if (bc2_list.size() >
                std::numeric_limits<size_t>::max() / bc1_list.size()) {
                throw std::runtime_error("spatial mask dimensions overflow");
            }

            phmap::flat_hash_map<int64_seq, uint16_t> new_bc1_to_idx;
            phmap::flat_hash_map<int64_seq, uint16_t> new_bc2_to_idx;
            new_bc1_to_idx.reserve(bc1_list.size());
            new_bc2_to_idx.reserve(bc2_list.size());
            for (size_t i = 0; i < bc1_list.size(); ++i) {
                int64_seq encoded(bc1_list[i]);
                if (!new_bc1_to_idx.emplace(
                        encoded, static_cast<uint16_t>(i)).second) {
                    throw std::runtime_error(
                        "BC1 axis whitelist contains a duplicate barcode");
                }
            }
            for (size_t j = 0; j < bc2_list.size(); ++j) {
                int64_seq encoded(bc2_list[j]);
                if (!new_bc2_to_idx.emplace(
                        encoded, static_cast<uint16_t>(j)).second) {
                    throw std::runtime_error(
                        "BC2 axis whitelist contains a duplicate barcode");
                }
            }

            const auto matrix_lines = streaming_utils::import_text(matrix_path);
            if (matrix_lines.size() != bc1_list.size()) {
                throw std::runtime_error(
                    "spatial mask row count does not match the BC1 axis");
            }

            std::vector<bool> new_bits(
                bc1_list.size() * bc2_list.size(), false);
            size_t valid_count = 0;
            for (size_t r = 0; r < matrix_lines.size(); ++r) {
                std::vector<std::string> cells;
                size_t start = 0;
                while (true) {
                    const auto comma = matrix_lines[r].find(',', start);
                    cells.push_back(seq_utils::trim(matrix_lines[r].substr(
                        start, comma == std::string::npos
                                   ? std::string::npos
                                   : comma - start)));
                    if (comma == std::string::npos) break;
                    start = comma + 1;
                }
                if (cells.size() != bc2_list.size()) {
                    throw std::runtime_error(
                        "spatial mask column count does not match the BC2 axis");
                }
                for (size_t c = 0; c < cells.size(); ++c) {
                    if (cells[c] != "0" && cells[c] != "1") {
                        throw std::runtime_error(
                            "spatial mask cells must be exactly 0 or 1");
                    }
                    if (cells[c] == "1") {
                        new_bits[r * bc2_list.size() + c] = true;
                        ++valid_count;
                    }
                }
            }

            rows = bc1_list.size();
            cols = bc2_list.size();
            bits = std::move(new_bits);
            bc1_to_idx = std::move(new_bc1_to_idx);
            bc2_to_idx = std::move(new_bc2_to_idx);

            if (verbose) {
                RAD_COUT << "[spat_mask] Loaded " << rows << "x" << cols
                          << " matrix (" << valid_count
                          << " valid pairs) from " << matrix_path << "\n";
            }
            return true;
        } catch (const std::exception& error) {
            bits.clear();
            rows = 0;
            cols = 0;
            bc1_to_idx.clear();
            bc2_to_idx.clear();
            throw std::runtime_error(
                std::string("spatial mask validation failed: ") +
                error.what());
        }
    }
};

/**
 * @brief Whitelist class managing true and filtered barcodes
 * @param wl_entry: Struct containing true barcodes, filtered barcodes, and global barcodes
 */
class whitelist {
    public:
        whitelist() = default;
        ~whitelist() = default;

/**
 * @brief Struct representing a whitelist entry with true barcodes, filtered barcodes, and global barcodes
 * @param true_bcs: `bc_multimap<int64_seq, barcode_entry>` for true barcodes
 * @param filter_bcs: `bc_multimap<int64_seq, barcode_entry>` for filtered barcodes
 * @param global_bcs: `bc_flat_set` for global barcodes
 * @param spat_wl: optional `spat_mask` for validating joint barcode pairs
 */
    struct wl_entry {
        bc_multimap<int64_seq, barcode_entry> true_bcs, filter_bcs;
        bc_flat_set global_bcs;
        std::optional<spat_mask> spat_wl;

/**
 * @brief Seed-based barcode candidate index
 *
 * Splits barcodes into 3 equal partitions and indexes each partition's
 * 2-bit packed value.  Any barcode within edit distance 1 of a query
 * must share at least 2 of 3 seed partitions (pigeonhole), so the
 * candidate set returned by query() is a superset of all ED-1 matches.
 *
 * Keys are packed uint32_t (supports seeds up to 16bp).  Values are
 * const barcode_entry* pointers into the owning whitelist — no copies.
 */
        struct seed_idx {
            struct key {
                uint8_t  partition = 0;
                uint32_t bits = 0;
                bool operator==(const key& o) const noexcept {
                    return partition == o.partition && bits == o.bits;
                }
            };
            struct key_hash {
                size_t operator()(const key& k) const noexcept {
                    return std::hash<uint64_t>{}(
                        (static_cast<uint64_t>(k.partition) << 32) | k.bits);
                }
            };

            std::unordered_map<key, std::vector<const barcode_entry*>, key_hash> index;
            size_t barcode_length = 0;
            bool ready = false;

            static std::array<std::pair<size_t, size_t>, 3> partitions(size_t len) {
                std::array<std::pair<size_t, size_t>, 3> windows{};
                size_t start = 0;
                const size_t base = len / 3;
                const size_t rem = len % 3;
                for (size_t i = 0; i < 3; ++i) {
                    const size_t seg_len = base + (i < rem ? 1 : 0);
                    const size_t end = start + seg_len;
                    windows[i] = {start, end};
                    start = end;
                }
                return windows;
            }

            /// Extract 2-bit packed seed from an int64_seq at [start, start+len).
            /// Single-chunk only (barcode length <= 32bp).
            static uint32_t extract_bits(const int64_seq& seq, size_t start, size_t len) {
                const size_t shift = 2 * (seq.length - start - len);
                const uint32_t mask = (1U << (2 * len)) - 1;
                return static_cast<uint32_t>((seq.bits[0] >> shift) & mask);
            }

            void build(const std::vector<const barcode_entry*>& entries, bool verbose) {
                index.clear();
                barcode_length = 0;
                ready = false;

                if (entries.empty()) return;

                const size_t bc_len = entries.front()->barcode.length;
                if (bc_len < 3) {
                    if (verbose) {
                        RAD_COUT << "[seed_idx] Skipping build: barcode length < 3\n";
                    }
                    return;
                }

                const auto windows = partitions(bc_len);
                index.reserve(entries.size() * 3);

                for (const auto* entry : entries) {
                    if (!entry || entry->barcode.length != bc_len) continue;
                    for (size_t pid = 0; pid < 3; ++pid) {
                        const auto [s, e] = windows[pid];
                        if (e <= s) continue;
                        key k{};
                        k.partition = static_cast<uint8_t>(pid);
                        k.bits = extract_bits(entry->barcode, s, e - s);
                        index[k].push_back(entry);
                    }
                }

                barcode_length = bc_len;
                ready = !index.empty();

                if (verbose) {
                    RAD_COUT << "[seed_idx] Built index: "
                              << index.size() << " keys from "
                              << entries.size() << " barcodes (len=" << bc_len
                              << ", partitions=3)\n";
                }
            }

            std::vector<const barcode_entry*> query(const int64_seq& barcode) const {
                std::vector<const barcode_entry*> candidates;
                if (!ready) return candidates;
                if (barcode.length != barcode_length) return candidates;

                const auto windows = partitions(barcode_length);
                for (size_t pid = 0; pid < 3; ++pid) {
                    const auto [s, e] = windows[pid];
                    if (e <= s) continue;
                    key k{};
                    k.partition = static_cast<uint8_t>(pid);
                    k.bits = extract_bits(barcode, s, e - s);

                    auto it = index.find(k);
                    if (it == index.end()) continue;
                    for (const auto* entry : it->second) {
                        candidates.push_back(entry);
                    }
                }
                return candidates;
            }
        };

        seed_idx true_seeds, global_seeds;

        void build_seed_idx(bool verbose) {
            auto true_entries = true_bcs.get_unique_entries();
            if (!true_entries.empty()) {
                true_seeds.build(true_entries, verbose);
            }
            auto global_entries = global_bcs.get_unique_entries();
            constexpr size_t max_global_seed_barcodes = 10000;
            if (global_entries.size() > max_global_seed_barcodes) {
                global_seeds = seed_idx{};
                if (verbose) {
                    RAD_COUT << "[seed_idx] Skipping global index: "
                              << global_entries.size()
                              << " unique barcodes exceed the 10,000-barcode cap\n";
                }
            } else if (!global_entries.empty()) {
                global_seeds.build(global_entries, verbose);
            } else {
                global_seeds = seed_idx{};
            }
        }

        bool has_seed_idx() const noexcept {
            return true_seeds.ready || global_seeds.ready;
        }

        std::vector<const barcode_entry*> query_bc_seeds(
            const int64_seq& query, const std::string& src) const {
            if (src == "true")  return true_seeds.query(query);
            if (src == "global") return global_seeds.query(query);
            auto r = true_seeds.query(query);
            auto g = global_seeds.query(query);
            r.insert(r.end(), g.begin(), g.end());
            return r;
        }

        template<class F> decltype(auto) with_wl(std::string_view src, F&& f) {
            if (src == "global") return f(global_bcs);
            return f(true_bcs);
        }
        template<class F> decltype(auto) with_wl(std::string_view src, F&& f) const {
            if (src == "global") return f(global_bcs);
            return f(true_bcs);
        }
/** 
 * @brief Generate mismatch barcodes by applying mutations and shifts to existing true barcodes
 * @param mutation_rounds Number of mutation rounds to apply 
*/
        void generate_mismatch_barcodes(int mutation_rounds, bool verbose, int nthreads = 1) {
                using namespace std::chrono;
                auto start_total = high_resolution_clock::now();
                
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "[generate_mismatch_barcodes] Generating shifted and mutated barcodes for "
                            << "mutation_rounds = " << mutation_rounds << "\n";
                        RAD_COUT << oss.str();
                    }
                }
                
                // TIMING: Collect original barcode entries
                auto start_collect = high_resolution_clock::now();

                std::vector<const barcode_entry*> originals = true_bcs.get_unique_entries();

                auto end_collect = high_resolution_clock::now();
                double collect_ms = duration<double, std::milli>(end_collect - start_collect).count();
                
                if (verbose) {
                    RAD_COUT << "[generate_mismatch_barcodes] Collected " << originals.size() 
                            << " unique originals in " << collect_ms << " ms\n";
                }
                
                // Memory tracking
                size_t initial_true_bcs_size = true_bcs.associations.size();
                size_t initial_memory_estimate = initial_true_bcs_size * 48; // Rough estimate: 48 bytes per association
                
                if (verbose) {
                    RAD_COUT << "[generate_mismatch_barcodes] Initial true_bcs size: " << initial_true_bcs_size 
                            << " associations (~" << (initial_memory_estimate / 1024 / 1024) << " MB)\n";
                }
                
                // TIMING: Generate mutations and shifts
                auto start_generation = high_resolution_clock::now();
                
                size_t total_mutations_generated = 0;
                size_t mutations_added = 0;
                size_t mutations_rejected_global = 0;
                
                // Track timing for individual operations
                double total_mutation_generation_ms = 0;
                double total_mutation_lookup_ms = 0;
                double total_mutation_insert_ms = 0;
                
                auto start_loop = high_resolution_clock::now();
                
                for (size_t i = 0; i < originals.size(); ++i) {
                    auto const *orig_be = originals[i];
                    const auto &orig_bits = orig_be->barcode;
                    
                    // Time individual barcode processing for first few
                    auto barcode_start = high_resolution_clock::now();
                                        
                    // === LEVENSHTEIN MUTATION GENERATION ===
                    auto mutation_gen_start = high_resolution_clock::now();
                    auto mutated = mutation_tools::generate_lv_barcodes(orig_bits, mutation_rounds);
                    auto mutation_gen_end = high_resolution_clock::now();
                    
                    double mutation_gen_time = duration<double, std::milli>(mutation_gen_end - mutation_gen_start).count();
                    total_mutation_generation_ms += mutation_gen_time;
                    total_mutations_generated += mutated.size();
                    
                    // === MUTATION LOOKUP & INSERT ===
                    auto mutation_lookup_start = high_resolution_clock::now();
                    for (auto const &m_bits : mutated) {
                        auto lookup_start = high_resolution_clock::now();
                        bool in_wl = global_bcs.check_wl_for(m_bits) || true_bcs.check_wl_for(m_bits);
                        auto lookup_end = high_resolution_clock::now();
                        total_mutation_lookup_ms += duration<double, std::milli>(lookup_end - lookup_start).count();
                        if (!in_wl) {
                            auto insert_start = high_resolution_clock::now();
                            true_bcs.insert_bc_entry(m_bits, *orig_be);
                            auto insert_end = high_resolution_clock::now();
                            total_mutation_insert_ms += duration<double, std::milli>(insert_end - insert_start).count();
                            mutations_added++;
                        } else {
                            mutations_rejected_global++;
                        }
                    }
                    
                    auto mutation_lookup_end = high_resolution_clock::now();
                    // Detailed timing for first 10 barcodes
                    if (verbose && i < 10) {
                        auto barcode_end = high_resolution_clock::now();
                        double total_time = duration<double, std::milli>(barcode_end - barcode_start).count();
                        RAD_COUT << "[generate_mismatch_barcodes] Barcode " << i + 1 << " (" 
                                << orig_bits.bits_to_sequence() << "): "
                                << mutated.size() << " mutations (" << mutation_gen_time << "ms gen), "
                                << "total " << total_time << "ms\n";
                    }
                    
                    // Memory usage check every 10000 barcodes
                    if (verbose && i % 1000 == 0 && i > 0) {
                        size_t current_true_bcs_size = true_bcs.associations.size();
                        size_t growth = current_true_bcs_size - initial_true_bcs_size;
                        size_t current_memory_estimate = current_true_bcs_size * 48;
                        
                        RAD_COUT << "[generate_mismatch_barcodes] After " << i << " barcodes: "
                                << current_true_bcs_size << " associations (+" << growth << "), "
                                << "~" << (current_memory_estimate / 1024 / 1024) << " MB total\n";
                    }
                }
                
                auto end_loop = high_resolution_clock::now();
                double loop_ms = duration<double, std::milli>(end_loop - start_loop).count();
                
                auto end_generation = high_resolution_clock::now();
                double generation_ms = duration<double, std::milli>(end_generation - start_generation).count();
                
                // Final memory usage
                size_t final_true_bcs_size = true_bcs.associations.size();
                size_t total_growth = final_true_bcs_size - initial_true_bcs_size;
                size_t final_memory_estimate = final_true_bcs_size * 48;
                
                // TIMING: Sanity check
                auto start_sanity = high_resolution_clock::now();
                auto end_sanity = high_resolution_clock::now();
                double sanity_ms = duration<double, std::milli>(end_sanity - start_sanity).count();
                
                auto end_total = high_resolution_clock::now();
                double total_ms = duration<double, std::milli>(end_total - start_total).count();
                
                if (verbose) {
                    RAD_COUT << "\n=== MISMATCH GENERATION TIMING REPORT ===\n";
                    RAD_COUT << "Total time: " << total_ms << " ms (" << (total_ms/1000.0) << " s)\n";
                    RAD_COUT << "  - Collection phase: " << collect_ms << " ms (" << (collect_ms/total_ms*100) << "%)\n";
                    RAD_COUT << "  - Generation loop: " << loop_ms << " ms (" << (loop_ms/total_ms*100) << "%)\n";
                    RAD_COUT << "  - Sanity check: " << sanity_ms << " ms (" << (sanity_ms/total_ms*100) << "%)\n";
                    
                    RAD_COUT << "\n=== DETAILED OPERATION BREAKDOWN ===\n";
                    RAD_COUT << "Mutation generation: " << total_mutation_generation_ms << " ms (" 
                            << (total_mutation_generation_ms/total_ms*100) << "%)\n";
                    RAD_COUT << "Mutation lookups: " << total_mutation_lookup_ms << " ms (" 
                            << (total_mutation_lookup_ms/total_ms*100) << "%)\n";
                    RAD_COUT << "Mutation insertions: " << total_mutation_insert_ms << " ms (" 
                            << (total_mutation_insert_ms/total_ms*100) << "%)\n";
                    
                    RAD_COUT << "\n=== GENERATION STATISTICS ===\n";
                    RAD_COUT << "Processed " << originals.size() << " original barcodes\n";
                    RAD_COUT << "Mutations: generated " << total_mutations_generated << ", added " << mutations_added 
                            << ", rejected " << mutations_rejected_global << "\n";
                    RAD_COUT << "Average mutations per barcode: " << (double)total_mutations_generated / originals.size() << "\n";
                    RAD_COUT << "Time per barcode: " << (loop_ms / originals.size()) << " ms\n";
                    
                    RAD_COUT << "\n=== MEMORY USAGE ===\n";
                    RAD_COUT << "Initial associations: " << initial_true_bcs_size 
                            << " (~" << (initial_memory_estimate / 1024 / 1024) << " MB)\n";
                    RAD_COUT << "Final associations: " << final_true_bcs_size 
                            << " (~" << (final_memory_estimate / 1024 / 1024) << " MB)\n";
                    RAD_COUT << "Growth: +" << total_growth << " associations (+~" 
                            << ((final_memory_estimate - initial_memory_estimate) / 1024 / 1024) << " MB)\n";
                    RAD_COUT << "Memory efficiency: " << (double)total_growth / (total_mutations_generated) * 100 
                            << "% of generated items were unique and added\n";
                }
            }
    };

    std::unordered_map<std::string, std::reference_wrapper<whitelist::wl_entry>> maps;
    std::unordered_map<std::string, whitelist::wl_entry> lists;
    // operator[] for easy insert/access: W["barcode_1"].insert_bc_entry(obs, corr);
    wl_entry & operator[](const std::string &class_id) {
        return maps.at(class_id).get();
    }

    //check all loaded class_ids in an entry
    std::vector<std::string> class_ids() const {
        std::vector<std::string> out;
        out.reserve(maps.size());
        for (auto const &kv : maps)
            out.push_back(kv.first);
        return out;
    }
/**
 * @brief Load barcodes from a specified kit name or file path
 * @param spec: Kit name or file path to load barcodes from
 * @param default_length: Default length of barcodes if not specified
 * @return Unordered set of `int64_seq` representing loaded barcodes
 */
    static std::unordered_set<int64_seq> load_barcodes(std::string const &spec, uint16_t default_length, bool verbose) {
        if (verbose) RAD_COUT << "[load_barcodes] Loading barcodes from path/kit: " << spec << "\n";
        // check if spec is a kit name or a file path
        if (verbose) {
            if(whitelist_utils::is_kit(spec)){
                RAD_COUT << "[load_barcodes] Kit name: " << spec << "\n";
                RAD_COUT << "[load_barcodes] Kit path: " << whitelist_utils::kit_to_path(spec) << "\n";
            } else {
                RAD_COUT << "[load_barcodes] File path: " << spec << "\n";
            }
        }
        std::string path = whitelist_utils::kit_to_path(spec);
        // detect bitlist vs char list
        bool isBitlist = whitelist_utils::check_if_bitlist(spec, verbose, /*N=*/10);
        if(verbose){
            std::string islist = isBitlist ? "IS" : "IS NOT";
            RAD_COUT << "[check_if_bitlist] Detected element " << islist << " a bit-configured whitelist!" << std::endl;
        }

        // Stream whitelist rows directly into the set.  The previous
        // import_text() call retained every source line while also building the
        // encoded set, creating an avoidable second full-whitelist peak.
        std::unique_ptr<std::istream> in;
        igzstream* gzip_in = nullptr;
        if (path.size() > 3 && path.substr(path.size() - 3) == ".gz") {
            auto stream = std::make_unique<igzstream>(path.c_str());
            gzip_in = stream.get();
            in = std::move(stream);
        } else {
            in = std::make_unique<std::ifstream>(path);
        }
        if (!in || !*in) {
            throw std::runtime_error("Failed to open whitelist: " + path);
        }

        std::unordered_set<int64_seq> out;
        std::string ln;
        while (std::getline(*in, ln)) {
            if (ln.empty()) continue;
            auto p = ln.find_first_of(",\t");
            std::string tok = seq_utils::trim(
                p == std::string::npos ? ln : ln.substr(0,p));
            if (tok.empty()) continue;

            int64_seq seq;
            if (isBitlist) {
                size_t consumed = 0;
                uint64_t raw_bits = 0;
                try {
                    raw_bits = std::stoull(tok, &consumed);
                } catch (...) { continue; }
                if (consumed != tok.size() || default_length == 0 ||
                    default_length > 32) {
                    continue;
                }
                if (default_length < 32) {
                    const uint64_t max_bits =
                        (uint64_t{1} << (2 * default_length)) - 1;
                    if (raw_bits > max_bits) continue;
                }
                seq.length = default_length;
                seq.bits = {raw_bits};
            } else {
                seq.sequence_to_bits(tok);
                if (seq.length == 0 || seq.bits.empty()) continue;
            }
            out.insert(seq);
        }
        if (gzip_in) {
            const int zlib_code = gzip_in->zlib_error_code();
            const std::string zlib_message = gzip_in->zlib_error_message();
            gzip_in->close();
            if (zlib_code != Z_OK && zlib_code != Z_STREAM_END) {
                throw std::runtime_error(
                    "whitelist gzip read failed: " + path +
                    " (zlib code " + std::to_string(zlib_code) +
                    (zlib_message.empty() ? ")" : ", " + zlib_message + ")"));
            }
            if (gzip_in->bad()) {
                throw std::runtime_error(
                    "whitelist gzip close failed: " + path);
            }
        } else if (in->bad()) {
            throw std::runtime_error("whitelist read failed: " + path);
        }
        return out;
    }
/**
 * @brief Populate a whitelist entry with true and global barcodes based on provided sets
 * @param out: Reference to the `wl_entry` to populate
 * @param sets: Vector of unordered sets of `int64_seq` representing barcode sets
 */
    void populate_entry(wl_entry &out, std::vector<std::unordered_set<int64_seq>> const &sets){
        if (sets.size() == 1) {
            auto const &A = sets[0];
            bool populate_global = false;
            size_t data_size = A.size();
            RAD_COUT << "[populate_entry] Single whitelist with " << data_size << " barcodes detected.\n";
            if (A.size() >= 100000) {
                populate_global = true;
            }
            RAD_COUT << "[populate_entry] Populating " << (populate_global ? "global_bcs" : "true_bcs") << "...\n";
            for (auto const &seq : A) {
                if (seq.bits.empty()) continue;
                barcode_entry be;
                be.barcode = seq;
                be.filtered = false;
                if(populate_global){
                    out.global_bcs.emplace(seq, std::move(be));
                } else {
                    out.true_bcs.emplace(seq, std::move(be));
                }
            }
            RAD_COUT << "[populate_entry] Population complete. Size of "
                      << (populate_global ? "global_bcs: " + std::to_string(out.global_bcs.size())
                                          : "true_bcs: " + std::to_string(out.true_bcs.associations.size()))
                      << "\n";
            return;
        }
        // two lists
        auto const &A = sets[0];
        auto const &B = sets[1];

        if ((A.size() != B.size())) {
            // pick smaller -> true, larger -> global
            auto const &small = (A.size() < B.size() ? A : B);
            auto const &large = (A.size() < B.size() ? B : A);

            for (auto const &seq : small) {
                if (seq.bits.empty()) continue;
                barcode_entry be;
                be.barcode = seq;
                be.filtered = false;
                out.true_bcs.emplace(seq, std::move(be));
            }
            for (auto const &seq : large) {
                //this line trims duplicate entries in global bcs
                if (!out.true_bcs.check_wl_for(seq)) {
                    if (seq.bits.empty()) continue;
                    barcode_entry be{};
                    be.barcode   = seq;
                    be.filtered  = false;
                    out.global_bcs.emplace(seq, be);
                }
            }
        }
        else {
            // fallback: merge both into true_bcs
            for (auto const &seq : A) {
                if (seq.bits.empty()) continue;
                barcode_entry be;
                be.barcode = seq;
                be.filtered = false;
                //be.flags = "";
                out.true_bcs.emplace(seq, std::move(be));
            }
            for (auto const &seq : B) {
                if (seq.bits.empty()) continue;
                barcode_entry be;
                be.barcode = seq;
                be.filtered = false;
                //be.flags = "";
                out.true_bcs.emplace(seq, std::move(be));
            }
        }
    }

    /**
     * Populate an entry without guessing roles from whitelist cardinality.
     * This path is used by the embedded R API, whose global_whitelist and
     * custom_whitelist arguments carry explicit, stable semantics.
     */
    void populate_entry_typed(
        wl_entry& out,
        const std::optional<std::unordered_set<int64_seq>>& global_set,
        const std::optional<std::unordered_set<int64_seq>>& true_set) {
        out.global_bcs.clear();
        out.true_bcs.clear();

        if (true_set) {
            for (const auto& seq : *true_set) {
                if (seq.bits.empty()) continue;
                barcode_entry entry;
                entry.barcode = seq;
                entry.filtered = false;
                out.true_bcs.emplace(seq, std::move(entry));
            }
        }
        if (global_set) {
            for (const auto& seq : *global_set) {
                if (seq.bits.empty() || out.true_bcs.check_wl_for(seq)) {
                    continue;
                }
                barcode_entry entry;
                entry.barcode = seq;
                entry.filtered = false;
                out.global_bcs.emplace(seq, std::move(entry));
            }
        }
    }

    wl_entry import_whitelist_typed(
        const std::optional<std::string>& global_spec,
        const std::optional<std::string>& true_spec,
        bool verbose,
        uint16_t default_length = 16) {
        if (!global_spec && !true_spec) {
            throw std::invalid_argument(
                "typed whitelist import requires a global or custom whitelist");
        }

        std::optional<std::unordered_set<int64_seq>> global_set;
        std::optional<std::unordered_set<int64_seq>> true_set;
        if (global_spec) {
            global_set = load_barcodes(*global_spec, default_length, verbose);
            if (global_set->empty()) {
                throw std::invalid_argument(
                    "global whitelist contains no valid barcodes: " +
                    *global_spec);
            }
        }
        if (true_spec) {
            true_set = load_barcodes(*true_spec, default_length, verbose);
            if (true_set->empty()) {
                throw std::invalid_argument(
                    "custom whitelist contains no valid barcodes: " +
                    *true_spec);
            }
        }

        wl_entry out;
        populate_entry_typed(out, global_set, true_set);
        return out;
    }
    
    //import a whitelist from file--can be true barcodes or not. 
    //cheap but easy way to import b/w both global or custom--just set
    //an arbitrary threshold for the number of true barcodes to be kept because ideally after a certain size of barcodes
    //you really just want to treat them both the same way
/**
 * @brief Import a whitelist from a specified field (kit name or file path)
 * @param field: Kit name or file path to import the whitelist from
 * @param default_length: Default length of barcodes if not specified
 * @return `wl_entry` containing the imported whitelist
 * @brief cheap but easy way to import b/w both global or custom--just set
 * an arbitrary threshold for the number of true barcodes to be kept because ideally after a certain size of barcodes
 * you really just want to treat them both the same way
 */
    wl_entry import_whitelist(std::string const &field, bool verbose, uint16_t default_length = 16) {
        if(verbose) RAD_COUT << "[import_whitelist] Importing whitelist from " << field << "\n";
        // split into 1 or 2 specs
        auto specs = whitelist_utils::parse_whitelist_specs(field);
        if (specs.empty()) {
            RAD_CERR << "[warning][import_whitelist] empty whitelist spec—returning empty entry\n";
            return {};
        } else {
            for (auto const &spec : specs) {
                if(verbose) RAD_COUT << "[import_whitelist] whitelist kit or path:" << spec << "\n";
            }
        }

        std::vector<std::unordered_set<int64_seq>> sets;
        for (auto const &spec : specs) {
            sets.push_back(load_barcodes(spec, default_length, verbose));
        }

        // assemble the wl_entry
        wl_entry out;
        out.global_bcs.clear();
        out.true_bcs.clear();
        populate_entry(out, sets);
        return out;
    }

};

/**
 * @namespace bc_mem_utils
 * @brief Namespace for memory utility functions related to barcode containers
 */
namespace bc_mem_utils {
// ---------- helpers used by both containers ----------
inline std::size_t get_int64seq_mem(const int64_seq& seq) {
    std::size_t base = sizeof(int64_seq);
    std::size_t vec_capacity = seq.bits.capacity() * sizeof(int64_t);
    return base + vec_capacity;
}

inline std::size_t get_bc_entry_mem(const barcode_entry& entry) {
    std::size_t base = sizeof(barcode_entry);
    std::size_t seq_vec = entry.barcode.bits.capacity() * sizeof(int64_t);
    return base + seq_vec;
}

// ===================================================
// =============== bc_multimap versions ==============
// ===================================================

template<typename Key, typename Value>
std::size_t approx_unique(const bc_multimap<Key, Value>& m) {
    const auto& set = m.debug_unique_values();
    if (set.empty()) return 0;

    std::size_t total_element_size = 0;
    std::size_t sample_count = 0;
    constexpr std::size_t max_samples = 100;

    for (const auto& elem : set) {
        total_element_size += get_bc_entry_mem(elem);
        if (++sample_count >= max_samples) break;
    }

    std::size_t avg_element_size = sample_count > 0
        ? total_element_size / sample_count
        : sizeof(Value);

    constexpr size_t node_overhead = 16;
    std::size_t element_memory = set.size() * (avg_element_size + node_overhead);

    std::size_t estimated_buckets = static_cast<size_t>(set.size() / 0.875);
    std::size_t bucket_memory = estimated_buckets * sizeof(void*);

    constexpr size_t num_submaps = 64;
    std::size_t submap_overhead = num_submaps * 128;

    return element_memory + bucket_memory + submap_overhead;
}

template<typename Key, typename Value>
std::size_t approx_assoc(const bc_multimap<Key, Value>& m) {
    const auto& mm = m.debug_associations();
    if (mm.empty()) return 0;

    std::size_t total_key_size = 0;
    std::size_t key_sample_count = 0;
    constexpr std::size_t max_key_samples = 100;

    for (const auto& [key, value_set] : mm) {
        if constexpr (std::is_same_v<Key, int64_seq>) {
            total_key_size += get_int64seq_mem(key);
        } else {
            total_key_size += sizeof(Key);
        }
        if (++key_sample_count >= max_key_samples) break;
    }

    std::size_t avg_key_size = key_sample_count > 0
        ? total_key_size / key_sample_count
        : sizeof(Key);

    std::size_t outer_buckets = mm.bucket_count() * sizeof(void*);

    constexpr size_t node_overhead = 16;
    std::size_t outer_nodes = mm.size() * (avg_key_size + sizeof(void*) + node_overhead);

    std::size_t inner_set_memory = 0;
    for (const auto& [key, value_set] : mm) {
        std::size_t inner_buckets = value_set.bucket_count() * sizeof(void*);
        std::size_t inner_elements = value_set.size() * sizeof(const Value*);
        constexpr std::size_t set_overhead = 64;
        inner_set_memory += inner_buckets + inner_elements + set_overhead;
    }

    constexpr size_t num_submaps = 64;
    std::size_t submap_overhead = num_submaps * 128;

    return outer_buckets + outer_nodes + inner_set_memory + submap_overhead;
}

template<typename Key, typename Value>
void print_submap_distribution(const bc_multimap<Key, Value>& m) {
    const auto& associations = m.debug_associations();
    if (associations.empty()) {
        RAD_COUT << "No associations to analyze\n";
        return;
    }
    size_t num_submaps = associations.subcnt();
    std::vector<size_t> submap_sizes(num_submaps, 0);
    size_t total_entries = 0;
    try {
        for (const auto& [k, value_set] : associations) {
            auto hash_val = associations.hash_function()(k);
            size_t submap_idx = hash_val & (num_submaps - 1);
            submap_sizes[submap_idx]++;
            total_entries++;
        }
    } catch (...) {
        RAD_COUT << "Error analyzing submap distribution\n";
        return;
    }
    RAD_COUT << "Submap distribution (" << total_entries << " total, "
              << num_submaps << " submaps):\n";

    size_t grid_cols = static_cast<size_t>(std::ceil(std::sqrt(num_submaps)));
    size_t grid_rows = static_cast<size_t>(std::ceil(static_cast<double>(num_submaps) / grid_cols));

    for (size_t row = 0; row < grid_rows; ++row) {
        for (size_t col = 0; col < grid_cols; ++col) {
            size_t idx = row * grid_cols + col;
            if (idx < num_submaps) {
                double pct = (static_cast<double>(submap_sizes[idx]) /
                             (total_entries ? total_entries : 1)) * 100.0;
                RAD_COUT << std::fixed << std::setprecision(1) << std::setw(6) << pct << "%";
            } else {
                RAD_COUT << "      ";
            }
            if (col < grid_cols - 1) RAD_COUT << " ";
        }
        RAD_COUT << "\n";
    }
    RAD_COUT << "Raw counts: [";
    for (size_t i = 0; i < num_submaps; ++i) {
        RAD_COUT << submap_sizes[i] << (i + 1 < num_submaps ? ", " : "");
    }
    RAD_COUT << "]\n";
    if (num_submaps > 1) {
        auto [min_it, max_it] = std::minmax_element(submap_sizes.begin(), submap_sizes.end());
        double min_pct = (total_entries ? (*min_it * 100.0 / total_entries) : 0.0);
        double max_pct = (total_entries ? (*max_it * 100.0 / total_entries) : 0.0);
        double balance_ratio = (*min_it > 0) ? static_cast<double>(*max_it) / (*min_it) : 0.0;
        RAD_COUT << "Balance: min=" << std::fixed << std::setprecision(1) << min_pct
                  << "%, max=" << max_pct << "%, ratio=" << std::setprecision(2) << balance_ratio << "\n";
    }
}

template<typename Key, typename Value>
void print_memory_report(const bc_multimap<Key, Value>& m, const std::string& name = "") {
    if (!name.empty()) RAD_COUT << "=== " << name << " ===\n";

    std::size_t unique_mem = approx_unique(m);
    std::size_t assoc_mem  = approx_assoc(m);
    std::size_t total_mem  = unique_mem + assoc_mem;

    RAD_COUT << "Counts:\n";
    RAD_COUT << "  Unique values:    " << m.unique_val_size() << "\n";
    RAD_COUT << "  Associations:     " << m.association_size() << "\n";

    RAD_COUT << "\nMemory breakdown:\n";
    RAD_COUT << "  Unique storage:   " << std::setw(10) << (unique_mem / 1024) << " KB  ("
              << std::fixed << std::setprecision(1) << (total_mem ? (100.0 * unique_mem / total_mem) : 0.0) << "%)\n";
    RAD_COUT << "  Association map:  " << std::setw(10) << (assoc_mem / 1024) << " KB  ("
              << (total_mem ? (100.0 * assoc_mem / total_mem) : 0.0) << "%)\n";
    RAD_COUT << "  Total:            " << std::setw(10) << (total_mem / 1024) << " KB  ("
              << std::setprecision(2) << (total_mem / 1024.0 / 1024.0) << " MB)\n";

    if (m.unique_val_size() > 0) {
        RAD_COUT << "\nPer-entry costs:\n";
        RAD_COUT << "  Unique storage:   " << (unique_mem / m.unique_val_size()) << " bytes/entry\n";
    }
    if (m.association_size() > 0) {
        RAD_COUT << "  Association:      " << (assoc_mem / m.association_size()) << " bytes/assoc\n";
    }

    RAD_COUT << "\nSample element sizes:\n";
    size_t sample = 0;
    for (const auto& elem : m.debug_unique_values()) {
        std::size_t elem_size = get_bc_entry_mem(elem);
        RAD_COUT << "  Entry " << sample << ": " << elem_size << " bytes "
                  << "(barcode vec capacity: " << elem.barcode.bits.capacity() << ")\n";
        if (++sample >= 5) break;
    }

    size_t bad_keys = m.validate_association_keys();
    size_t null_vals = m.validate_association_values();
    if (bad_keys > 0 || null_vals > 0) {
        RAD_COUT << "\n Issues: " << bad_keys << " bad keys, " << null_vals << " null values\n";
    }

    print_submap_distribution(m);
    RAD_COUT << "\n";
}

// ===================================================
// ================ bc_flat_set versions =============
// ===================================================

// More direct estimate for a set: element objects + buckets
inline std::size_t approx_unique(const bc_flat_set& fs) {
    const auto& set = fs.debug_unique_values();
    if (set.empty()) return 0;

    std::size_t total_element_size = 0;
    std::size_t sample_count = 0;
    constexpr std::size_t max_samples = 100;

    // Sample actual entry size (barcode_entry includes the int64_seq vector)
    for (const auto& elem : set) {
        total_element_size += get_bc_entry_mem(elem);
        if (++sample_count >= max_samples) break;
    }
    std::size_t avg_element_size = sample_count > 0
        ? total_element_size / sample_count
        : sizeof(barcode_entry);

    // Per-node allocator overhead (rough heuristic)
    constexpr size_t node_overhead = 16;
    std::size_t element_memory = set.size() * (avg_element_size + node_overhead);

    // Bucket array (phmap keeps LF <= ~0.875; we just estimate)
    std::size_t estimated_buckets = static_cast<size_t>(set.size() / 0.875);
    std::size_t bucket_memory = estimated_buckets * sizeof(void*);

    // Metadata overhead for parallel shards (phmap): estimate similar scale
    constexpr size_t shard_overhead = 64 * 64; // 64 shards * ~64B each
    return element_memory + bucket_memory + shard_overhead;
}

// Associations concept does not apply to a pure set.
inline std::size_t approx_assoc(const bc_flat_set&) {
    return 0;
}

// Reuse the fake view emitted by bc_flat_set::debug_associations()
inline void print_submap_distribution(const bc_flat_set& fs) {
    auto mm = fs.debug_associations();
    if (mm.rows.empty()) {
        RAD_COUT << "No associations to analyze\n";
        return;
    }
    size_t num_submaps = mm.subcnt(); // returns 1 in your fake view
    std::vector<size_t> submap_sizes(num_submaps, 0);
    size_t total_entries = 0;

    for (const auto& row : mm.rows) {
        auto hash_val = mm.hash_function()(row.k);
        size_t submap_idx = hash_val & (num_submaps - 1);
        submap_sizes[submap_idx]++;
        total_entries++;
    }

    RAD_COUT << "Submap distribution (" << total_entries << " total, "
              << num_submaps << " submaps):\n";
    size_t grid_cols = static_cast<size_t>(std::ceil(std::sqrt(num_submaps)));
    size_t grid_rows = static_cast<size_t>(std::ceil(static_cast<double>(num_submaps) / grid_cols));
    for (size_t row = 0; row < grid_rows; ++row) {
        for (size_t col = 0; col < grid_cols; ++col) {
            size_t idx = row * grid_cols + col;
            if (idx < num_submaps) {
                double pct = (static_cast<double>(submap_sizes[idx]) /
                             (total_entries ? total_entries : 1)) * 100.0;
                RAD_COUT << std::fixed << std::setprecision(1) << std::setw(6) << pct << "%";
            } else {
                RAD_COUT << "      ";
            }
            if (col < grid_cols - 1) RAD_COUT << " ";
        }
        RAD_COUT << "\n";
    }
    RAD_COUT << "Raw counts: [";
    for (size_t i = 0; i < num_submaps; ++i) {
        RAD_COUT << submap_sizes[i] << (i + 1 < num_submaps ? ", " : "");
    }
    RAD_COUT << "]\n";
    if (num_submaps > 1) {
        auto [min_it, max_it] = std::minmax_element(submap_sizes.begin(), submap_sizes.end());
        double min_pct = (total_entries ? (*min_it * 100.0 / total_entries) : 0.0);
        double max_pct = (total_entries ? (*max_it * 100.0 / total_entries) : 0.0);
        double balance_ratio = (*min_it > 0) ? static_cast<double>(*max_it) / (*min_it) : 0.0;
        RAD_COUT << "Balance: min=" << std::fixed << std::setprecision(1) << min_pct
                  << "%, max=" << max_pct << "%, ratio=" << std::setprecision(2) << balance_ratio << "\n";
    }
}

inline void print_memory_report(const bc_flat_set& fs, const std::string& name = "") {
    if (!name.empty()) RAD_COUT << "=== " << name << " ===\n";

    std::size_t unique_mem = approx_unique(fs);
    std::size_t assoc_mem  = approx_assoc(fs); // 0
    std::size_t total_mem  = unique_mem + assoc_mem;

    RAD_COUT << "Counts:\n";
    RAD_COUT << "  Unique values:    " << fs.unique_val_size() << "\n";
    RAD_COUT << "  Associations:     " << fs.association_size() << "\n";

    RAD_COUT << "\nMemory breakdown:\n";
    RAD_COUT << "  Unique storage:   " << std::setw(10) << (unique_mem / 1024) << " KB  ("
              << std::fixed << std::setprecision(1) << (total_mem ? (100.0 * unique_mem / total_mem) : 0.0) << "%)\n";
    RAD_COUT << "  Association map:  " << std::setw(10) << (assoc_mem / 1024) << " KB  ("
              << (total_mem ? (100.0 * assoc_mem / total_mem) : 0.0) << "%)\n";
    RAD_COUT << "  Total:            " << std::setw(10) << (total_mem / 1024) << " KB  ("
              << std::setprecision(2) << (total_mem / 1024.0 / 1024.0) << " MB)\n";

    if (fs.unique_val_size() > 0) {
        RAD_COUT << "\nPer-entry costs:\n";
        RAD_COUT << "  Unique storage:   " << (unique_mem / fs.unique_val_size()) << " bytes/entry\n";
    }
    if (fs.association_size() > 0) {
        RAD_COUT << "  Association:      " << (assoc_mem / fs.association_size()) << " bytes/assoc\n";
    }

    RAD_COUT << "\nSample element sizes:\n";
    size_t sample = 0;
    for (const auto& elem : fs.debug_unique_values()) {
        std::size_t elem_size = get_bc_entry_mem(elem);
        RAD_COUT << "  Entry " << sample << ": " << elem_size << " bytes "
                  << "(barcode vec capacity: " << elem.barcode.bits.capacity() << ")\n";
        if (++sample >= 5) break;
    }

    // Flat set fake view returns 0 for validations (by design)
    size_t bad_keys = fs.validate_association_keys();
    size_t null_vals = fs.validate_association_values();
    if (bad_keys > 0 || null_vals > 0) {
        RAD_COUT << "\n Issues: " << bad_keys << " bad keys, " << null_vals << " null values\n";
    }

    print_submap_distribution(fs);
    RAD_COUT << "\n";
}

// ===================================================
// ============ Mixed/symmetric snapshot API =========
// ===================================================

// Generic snapshot that works for *any* wl types that overloads support
// T must have unique_val_size(), association_size() and be accepted by approx_* overloads
template <typename WLTrue, typename WLGlobal>
inline void print_mem_snapshot(const WLTrue& true_bcs,
                               const WLGlobal& global_bcs,
                               size_t iteration)
{
    RAD_COUT << "\n--- Iteration " << iteration << " Memory Snapshot ---\n";

    std::size_t true_unique  = approx_unique(true_bcs);
    std::size_t true_assoc   = approx_assoc(true_bcs);
    std::size_t global_unique= approx_unique(global_bcs);
    std::size_t global_assoc = approx_assoc(global_bcs);
    std::size_t total        = true_unique + true_assoc + global_unique + global_assoc;

    RAD_COUT << "true_bcs:   " << std::setw(8) << true_bcs.unique_val_size() << " unique, "
              << std::setw(8) << true_bcs.association_size() << " assoc, "
              << std::setw(8) << ((true_unique + true_assoc) / 1024 / 1024) << " MB\n";
    RAD_COUT << "global_bcs: " << std::setw(8) << global_bcs.unique_val_size() << " unique, "
              << std::setw(8) << global_bcs.association_size() << " assoc, "
              << std::setw(8) << ((global_unique + global_assoc) / 1024 / 1024) << " MB\n";
    RAD_COUT << "TOTAL: " << (total / 1024 / 1024) << " MB\n";
}

} // namespace bc_mem_utils
