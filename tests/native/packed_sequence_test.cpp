#include "rad/rad_headers.h"
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <unordered_set>

using namespace rrad_rad_core;

namespace {
void require(bool value, const std::string& message) {
    if (!value) throw std::runtime_error(message);
}

std::string dna(std::mt19937_64& random, size_t length) {
    std::string result(length, 'A');
    for (char& base : result) base = "ACTG"[random() & 3];
    return result;
}

std::vector<int> string_distances(const std::string& pattern, const std::string& text, bool infix) {
    std::vector<int> previous(text.size() + 1), current(text.size() + 1);
    for (size_t j = 0; j <= text.size(); ++j) previous[j] = infix ? 0 : static_cast<int>(j);
    for (size_t i = 1; i <= pattern.size(); ++i) {
        current[0] = static_cast<int>(i);
        for (size_t j = 1; j <= text.size(); ++j) {
            current[j] = std::min({previous[j] + 1, current[j - 1] + 1,
                                  previous[j - 1] + (pattern[i - 1] != text[j - 1])});
        }
        previous.swap(current);
    }
    return previous;
}

struct partial_result { int distance, end, count; };
partial_result string_partial(const std::string& pattern, const std::string& text) {
    const auto scores = string_distances(pattern, text, true);
    const int minimum = *std::min_element(scores.begin(), scores.end());
    partial_result result{minimum, -1, 0};
    // Preserve the established API's sentinel when no column improves on
    // deleting the complete pattern before consuming any text.
    if (minimum == static_cast<int>(pattern.size())) return result;
    for (size_t j = 1; j < scores.size(); ++j) {
        if (scores[j] == minimum) {
            if (result.end < 0) result.end = static_cast<int>(j - 1);
            ++result.count;
        }
    }
    return result;
}

void check_partial(const std::string& pattern, const std::string& text) {
    const int64_seq packed_pattern(pattern), packed_text(text);
    const auto expected = string_partial(pattern, text);
    const std::string context = "partial lengths " + std::to_string(pattern.size()) + "/" +
                                std::to_string(text.size());
    for (int cap : {0, 1, 2, 3, 4, 64}) {
        int end = -99, count = -99;
        const int actual = packed_pattern.visit([&](auto p) {
            return packed_text.visit([&](auto t) {
                return mutation_tools::bit_partial_match(p, t, packed_pattern.length, packed_text.length,
                                                         cap, &end, &count);
            });
        });
        const int accepted = expected.distance <= cap ? expected.distance : -1;
        require(actual == accepted, context + ": distance differs from string DP");
        require(end == expected.end && count == expected.count,
                context + ": best end/ties differ from string DP");
        const int without_position = packed_pattern.visit([&](auto p) {
            return packed_text.visit([&](auto t) {
                return mutation_tools::bit_partial_match(p, t, packed_pattern.length, packed_text.length, cap);
            });
        });
        require(without_position == accepted, context + ": fast path differs from string DP");
        const auto wrapper_expected = pattern.size() <= text.size() ? expected : string_partial(text, pattern);
        const int wrapper = mutation_tools::int64_lvdist(packed_pattern, packed_text, cap);
        require(wrapper == (wrapper_expected.distance <= cap ? wrapper_expected.distance : -1),
                context + ": sequence wrapper lost a packed word");
    }
    require(packed_pattern.bits_to_sequence() == pattern && packed_text.bits_to_sequence() == text,
            "comparison changed a stored sequence");
    require(packed_pattern.is_wide() == (pattern.size() > 32) && packed_text.is_wide() == (text.size() > 32),
            "comparison changed initialized storage width");
}

std::unordered_set<std::string> decoded(const std::unordered_set<int64_seq>& sequences) {
    std::unordered_set<std::string> result;
    for (const auto& sequence : sequences) {
        require(sequence.is_valid(), "mutation produced an invalid sequence");
        require(sequence.is_wide() == (sequence.length > 32), "mutation retained an incorrect storage width");
        result.insert(sequence.bits_to_sequence());
    }
    require(result.size() == sequences.size(), "mutation retained duplicate string identities");
    return result;
}

std::unordered_set<std::string> string_substitutions(const std::string& sequence) {
    std::unordered_set<std::string> result;
    for (size_t i = 0; i < sequence.size(); ++i) for (char base : std::string("ACTG")) {
        if (base == sequence[i]) continue;
        auto changed = sequence;
        changed[i] = base;
        result.insert(std::move(changed));
    }
    return result;
}

std::unordered_set<std::string> string_edits(const std::string& sequence) {
    auto result = string_substitutions(sequence);
    if (sequence.size() > 1) for (size_t i = 0; i < sequence.size(); ++i)
        result.insert(sequence.substr(0, i) + sequence.substr(i + 1));
    if (sequence.size() < 64) for (size_t i = 0; i <= sequence.size(); ++i)
        for (char base : std::string("ACTG"))
            result.insert(sequence.substr(0, i) + base + sequence.substr(i));
    return result;
}

void check_encoding(std::mt19937_64& random) {
    static_assert(std::is_same_v<decltype(packed_barcode<uint64_t>{}.bits), uint64_t>);
    static_assert(std::is_same_v<decltype(packed_barcode<uint128_t>{}.bits), uint128_t>);
    static_assert(sizeof(packed_barcode<uint64_t>) < sizeof(packed_barcode<uint128_t>));
    for (size_t length : {size_t{1}, size_t{31}, size_t{32}, size_t{33}, size_t{63}, size_t{64}}) {
        const auto sequence = dna(random, length);
        const int64_seq encoded(sequence);
        require(encoded.bits_to_sequence() == sequence, "packed round trip failed");
        require(encoded.length == length && encoded.is_wide() == (length > 32), "wrong initialized width");
        require(encoded.visit([](auto word) { return sizeof(word); }) == (length > 32 ? 16 : 8),
                "dispatcher selected the wrong integer type");
        const auto before = encoded;
        const uint128_t promoted = encoded.as_uint128();
        require(encoded == before && encoded.is_wide() == (length > 32), "promotion changed original storage");
        require(int64_seq::from_word(promoted, encoded.length) == encoded, "promotion changed sequence bits");
        for (size_t start = 0; start < length; ++start) {
            const size_t count = std::min<size_t>(length - start, 32);
            require(encoded.slice(start, count).bits_to_sequence() == sequence.substr(start, count),
                    "slice across a packed-word boundary differs from string substring");
        }
        require(encoded.slice(0, length) == encoded, "full-width slice changed sequence");
    }
    const int64_seq a("A"), aa("AA");
    require(a.as_uint64() == aa.as_uint64() && a != aa, "leading A length identity was lost");
    const int64_seq default_empty, invalid("ACNT"), explicit_empty("");
    require(default_empty.empty() && !invalid.empty() && !invalid.is_valid(), "invalid sentinel changed");
    require(default_empty != invalid && explicit_empty == invalid, "empty/invalid identities changed");
    require(packed_mask<uint64_t>(0) == 0 && packed_mask<uint64_t>(32) == ~uint64_t{0}, "64-bit mask boundary failed");
    require(packed_mask<uint128_t>(0) == 0 && packed_mask<uint128_t>(64) == ~uint128_t{0}, "128-bit mask boundary failed");
    const int64_seq all_g(std::string(32, 'G'));
    require(all_g.as_uint64() == ~uint64_t{0}, "high sign bit was lost");
    require(int64_seq(int64_t{-1}, 32) == all_g, "legacy signed bitlist bits were sign-extended incorrectly");
    const auto raw_high = int64_seq::from_word(uint64_t{1} << 63, 16);
    require(raw_high.as_uint64() == (uint64_t{1} << 63), "legacy unused raw bits were masked away");
    const int64_seq high_c("C" + std::string(32, 'A')), high_g("G" + std::string(32, 'A'));
    require(static_cast<uint64_t>(high_c.as_uint128()) == static_cast<uint64_t>(high_g.as_uint128()),
            "wide-hash fixture must share its low word");
    require(std::hash<int64_seq>{}(high_c) != std::hash<int64_seq>{}(high_g), "hash ignores the upper packed word");
    require(std::hash<int64_seq>{}(int64_seq("A")) != std::hash<int64_seq>{}(int64_seq("AA")),
            "hash ignores leading-A sequence length");
    std::unordered_set<int64_seq> distinct{high_c, high_g, a, aa};
    require(distinct.size() == 4, "packed hash/equality merged distinct identities");
}

void check_matching(std::mt19937_64& random) {
    for (size_t p : {size_t{1}, size_t{31}, size_t{32}, size_t{33}, size_t{63}, size_t{64}})
        for (size_t t : {size_t{1}, size_t{31}, size_t{32}, size_t{33}, size_t{63}, size_t{64}})
            check_partial(dna(random, p), dna(random, t));
    check_partial("ACG", "ACGACG");
    check_partial("A", std::string(64, 'A'));
    check_partial(std::string(64, 'A'), std::string(64, 'G'));
    check_partial(std::string(64, 'G'), std::string(64, 'G'));
    for (int trial = 0; trial < 220; ++trial) {
        const size_t text_length = 1 + random() % 64;
        std::string text = dna(random, text_length);
        const size_t pattern_length = 1 + random() % text_length;
        const size_t start = random() % (text_length - pattern_length + 1);
        std::string pattern = text.substr(start, pattern_length);
        for (int edits = 0, count = trial % 5; edits < count; ++edits)
            pattern[random() % pattern.size()] = "ACTG"[random() & 3];
        if (trial % 3 == 0 && pattern.size() < 64)
            pattern.insert(random() % (pattern.size() + 1), 1, "ACTG"[random() & 3]);
        if (trial % 5 == 0 && pattern.size() > 1) pattern.erase(random() % pattern.size(), 1);
        check_partial(pattern, text);
    }
    // The historical failing window: a 25-base key embedded in 33 bases.
    const auto stereo = dna(random, 25);
    check_partial(stereo, "ACTG" + stereo + "TGCA");
    for (size_t length : {size_t{1}, size_t{31}, size_t{32}, size_t{33}, size_t{63}, size_t{64}}) {
        const auto left = dna(random, length);
        auto right = left;
        for (size_t edits = 0; edits < 5; ++edits) {
            if (edits) right[random() % length] = "ACTG"[random() & 3];
            const int expected = string_distances(left, right, false).back();
            const int64_seq p(left), t(right);
            for (int cap = 0; cap <= 4; ++cap) {
                const int actual = p.visit([&](auto pw) {
                    return t.visit([&](auto tw) { return mutation_tools::bit_ld(pw, tw, length, cap); });
                });
                require(actual == (expected <= cap ? expected : -1), "global distance differs from string DP");
            }
        }
    }
}

void check_mutations(std::mt19937_64& random) {
    for (size_t length : {size_t{1}, size_t{31}, size_t{32}, size_t{33}, size_t{63}, size_t{64}}) {
        const auto sequence = dna(random, length);
        const int64_seq encoded(sequence);
        const auto points = decoded(mutation_tools::generate_point_mutations(encoded));
        require(points == string_substitutions(sequence) && points.size() == length * 3,
                "packed point mutations differ from string substitutions");
        require(decoded(mutation_tools::generate_lv_barcodes(encoded, 1)) == string_edits(sequence),
                "packed edits differ from string insertion/deletion/substitution oracle");
        if (length >= 2) {
            std::unordered_set<std::string> expected;
            const auto core = sequence.substr(0, length - 1);
            for (char base : std::string("ACTG")) { expected.insert(base + core); expected.insert(core + base); }
            expected.erase(sequence);
            require(decoded(mutation_tools::generate_shifted_barcodes(encoded, 1)) == expected,
                    "width extension changed the legacy shifted candidate set");
        }
        if (length == 32 || length == 33) {
            const auto hamming = decoded(mutation_tools::generate_mutated_barcodes(encoded, 2));
            require(hamming.size() == length * 3 + length * (length - 1) * 9 / 2,
                    "Hamming shell has missing or duplicate substitutions");
            for (const auto& changed : hamming) {
                size_t distance = 0;
                require(changed.size() == length, "Hamming mutation changed length");
                for (size_t i = 0; i < length; ++i) distance += changed[i] != sequence[i];
                require(distance >= 1 && distance <= 2, "Hamming mutation changed too many bases");
            }
        }
    }
    const std::string small = "ACTG";
    std::unordered_set<std::string> expected{small}, frontier{small};
    for (int round = 0; round < 2; ++round) {
        std::unordered_set<std::string> next;
        for (const auto& sequence : frontier) for (const auto& changed : string_edits(sequence))
            if (expected.insert(changed).second) next.insert(changed);
        frontier = std::move(next);
    }
    expected.erase(small);
    require(decoded(mutation_tools::generate_lv_barcodes(int64_seq(small), 2)) == expected,
            "multiple edit rounds differ from string oracle");
    const int64_seq narrow(std::string(32, 'G'));
    std::unordered_set<int64_seq> legacy;
    for (int64_t raw : mutation_tools::generate_point_mutations_raw(int64_t{-1}, 32)) legacy.emplace(raw, 32);
    require(legacy == mutation_tools::generate_point_mutations(narrow), "legacy signed raw mutation API changed bits");
}
}

int main() {
    try {
        std::mt19937_64 random(0x7061636b6564ULL);
        check_encoding(random);
        check_matching(random);
        check_mutations(random);
        int64_seq parsed;
        require(parse_packed_decimal("18446744073709551615", 32, parsed) &&
                parsed.bits_to_sequence() == std::string(32, 'G'), "unsigned64 decimal boundary failed");
        require(parse_packed_decimal("340282366920938463463374607431768211455", 64, parsed) &&
                parsed.bits_to_sequence() == std::string(64, 'G'), "unsigned128 decimal boundary failed");
        require(parse_packed_decimal("-1", 33, parsed) &&
                parsed.bits_to_sequence() == "A" + std::string(32, 'G'), "signed legacy word was sign extended");
        const auto before = parsed;
        require(!parse_packed_decimal("340282366920938463463374607431768211456", 64, parsed) &&
                parsed == before, "decimal overflow changed a sequence or was accepted");
        require(!parse_packed_decimal("18446744073709551616", 32, parsed), "narrow decimal overflow accepted");
        require(parse_packed_decimal("  +123suffix", 16, parsed) && parsed.as_uint64() == 123,
                "legacy numeric prefix parsing changed");
        require(!parse_packed_decimal("123suffix", 16, parsed, true) &&
                !parse_packed_decimal("16", 1, parsed, true), "strict numeric validation failed");
        std::cout << "Packed 64/128-bit sequence and matching checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
