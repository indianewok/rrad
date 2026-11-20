#include "rad/rad_headers.h"
#include <set>
#include <iostream>
#include <random>

using namespace rrad_rad_core;
namespace {
void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
std::set<std::string> one_edit(const std::string& sequence) {
    std::set<std::string> result;
    for (size_t i = 0; i < sequence.size(); ++i) {
        for (const char base : std::string("ACGT")) {
            auto changed = sequence;
            changed[i] = base;
            if (changed != sequence) result.insert(changed);
        }
        if (sequence.size() > 1)
            result.insert(sequence.substr(0, i) + sequence.substr(i + 1));
    }
    if (sequence.size() < 64) for (size_t i = 0; i <= sequence.size(); ++i)
        for (const char base : std::string("ACGT"))
            result.insert(sequence.substr(0, i) + base + sequence.substr(i));
    return result;
}
std::set<std::string> reference_candidates(const std::string& sequence, int rounds,
                                          uint64_t lengths) {
    std::set<std::string> result, frontier{sequence};
    for (int round = 0; round < rounds; ++round) {
        std::set<std::string> next;
        for (const auto& value : frontier) {
            const auto variants = one_edit(value);
            next.insert(variants.begin(), variants.end());
        }
        result.insert(next.begin(), next.end());
        frontier = std::move(next);
    }
    result.erase(sequence);
    for (auto it = result.begin(); it != result.end();) {
        if (lengths && !(lengths & (uint64_t{1} << (it->size() - 1)))) it = result.erase(it);
        else ++it;
    }
    return result;
}
std::string random_dna(std::mt19937& rng, size_t length) {
    std::string value(length, 'A');
    for (auto& base : value) base = "ACGT"[rng() % 4];
    return value;
}
void compare_candidates(const std::string& text, int edits, uint64_t lengths) {
    const auto generated = mutation_tools::generate_lv_candidates(int64_seq(text), edits, lengths);
    std::set<std::string> actual;
    for (const auto& value : generated) {
        require(value.is_valid(), "Invalid generated candidate");
        actual.insert(value.bits_to_sequence());
    }
    require(actual.size() == generated.size(), "Duplicate generated candidate");
    require(actual == reference_candidates(text, edits, lengths), "Indel candidate set differs from string oracle");
}
}
int main() {
    std::mt19937 random(20261113);
    for (int length : {1, 2, 3, 6, 8}) for (int trial = 0; trial < 3; ++trial) {
        const auto text = random_dna(random, length);
        compare_candidates(text, 2, uint64_t{1} << (length - 1));
        compare_candidates(text, 2, 0);
    }
    for (int length : {31, 32, 33, 63, 64}) {
        const auto text = random_dna(random, length);
        compare_candidates(text, 1, 0);
        if (length < 64) compare_candidates(text, 1, uint64_t{1} << length);
    }
    // A deletion/insertion can be within two edits while more than two bases
    // differ positionally. The global candidate path must retain this identity.
    const std::string observed = "ACGTACGA", shifted = "CGTACGAT";
    auto candidates = mutation_tools::generate_lv_candidates(int64_seq(observed), 2, uint64_t{1} << 7);
    require(std::find(candidates.begin(), candidates.end(), int64_seq(shifted)) != candidates.end(),
            "Compensating indel barcode was not generated");

    // Full-width barcodes still have flanking sequence. Validate both identity
    // correction and a one-base positional refinement, in either orientation.
    for (int length : {57, 63, 64}) {
        const std::string truth = random_dna(random, length);
        whitelist::wl_entry entry;
        barcode_entry barcode;
        barcode.barcode = int64_seq(truth);
        entry.true_bcs.emplace(barcode.barcode, barcode);
        ReadLayout layout;
        layout.wl_map.lists.emplace("catalog", std::move(entry));
        layout.wl_map.maps.emplace("barcode", std::ref(layout.wl_map.lists.at("catalog")));
        for (bool reverse : {false, true}) for (int shift : {0, 1}) {
            const std::string oriented = "TGCA" + truth + "GTCA";
            const std::string sequence = reverse ? seq_utils::revcomp(oriented) : oriented;
            const int first = reverse ? 5 + shift : 5 - shift;
            seq_element element("barcode", "barcode", std::nullopt, {first, first + length - 1},
                                "variable", 1, reverse ? "reverse" : "forward");
            element.seq = sequence.substr(first - 1, length);
            read_streaming::sequence read{"wide", "", sequence, "", "", false};
            barcode_correction::bc_position position;
            auto corrected = barcode_correction::correct_barcode(element, layout, read, false, 2,
                                                                  "offensive", &position);
            require(corrected && *corrected == int64_seq(truth), "Wide padded correction lost barcode");
            require(position.kind == barcode_correction::bc_position::located && position.shift == shift,
                    "Wide padded correction lost positional refinement");
        }
        const std::string window = "TGCA" + truth + "GTCA";
        auto mismatch = truth;
        mismatch[23] = truth[23] == 'A' ? 'C' : 'A';
        int end = -1, count = 0;
        const int distance = barcode_correction::barcode_distance(
            int64_seq(truth), int64_seq(mismatch), 2, &window, &end, &count);
        require(distance == 1 && end == length + 3 && count == 1,
                "Wide window mismatch scoring or unique end differs from full-window alignment");
    }
    // Local adapter fragments are not complete internal adapters. They remain
    // valid at a physical read edge where the missing bases are unobserved.
    aligner_tools aligner;
    const std::string adapter = "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAA";
    const auto suffix = adapter.substr(adapter.size() - 16);
    const auto internal = aligner.align_static_elements(adapter, std::string(50, 'C') + suffix +
                                                        std::string(50, 'C'), false, 2);
    const auto at_start = aligner.align_static_elements(adapter, suffix + std::string(50, 'C'), false, 2);
    const auto artificial_edge = aligner.align_static_elements(adapter, suffix + std::string(50, 'C'),
                                                               false, 2, "", true, 1, -1,
                                                               nullptr, 1, -1, true, false, true);
    require(!internal.success, "Internal partial adapter bypassed full-length edit budget");
    require(at_start.success && !at_start.query_complete, "Physical read-start adapter clip was rejected");
    require(!artificial_edge.success, "Artificial molecule boundary was treated as a physical read end");
    const std::string linker = "CAGCACCTGATTACAAGGTGCTG";
    const std::string reverse_linker = seq_utils::revcomp(linker);
    ReadLayout opposite_layout;
    opposite_layout.layout.insert(ReadElement("linker", linker, linker, 23, "static", 1, "forward", "linker"));
    opposite_layout.layout.insert(ReadElement("rc_linker", reverse_linker, reverse_linker, 23,
                                               "static", 2, "reverse", "linker"));
    auto verify_statics = [&](const std::string& sequence, bool has_forward, int reverse_start) {
        read_streaming::sequence read{"opposite", "", sequence, "", "", false};
        SigString signature(read.id, static_cast<int>(sequence.size()));
        std::vector<seq_element> alternatives;
        signature.sigalign_static(read, opposite_layout, false, &alternatives);
        const auto& index = signature.by_id();
        const auto reverse = index.find("rc_linker");
        require(reverse != index.end() && reverse->position.first == reverse_start && reverse->edit_distance == 0,
                "An imperfect forward hit masked the exact reverse adapter");
        require((index.find("linker") != index.end()) == has_forward,
                "Opposite-strand priority removed an independent forward repeat");
        if (!has_forward)
            require(std::none_of(alternatives.begin(), alternatives.end(), [](const auto& hit) {
                return hit.class_id == "linker";
            }), "Weaker opposite hit leaked through concat alternatives");
    };
    verify_statics("TTTGA" + reverse_linker + std::string(40, 'A'), false, 6);
    auto imperfect_forward = linker;
    imperfect_forward[10] = imperfect_forward[10] == 'A' ? 'C' : 'A';
    verify_statics(linker + std::string(40, 'A') + reverse_linker, true, 64);
    verify_statics(imperfect_forward + std::string(40, 'A') + reverse_linker, true, 64);
    std::cout << "Indel candidates, padded barcodes and adapter completeness/priority passed\n";
}
