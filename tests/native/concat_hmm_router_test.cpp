// Adapted from RAD dev 6167540; package namespace and cancellation regression.
#include "rad/rad_headers.h"
#include <iostream>
#include <stdexcept>

using namespace rrad_rad_core;

namespace {
void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}

concat_hmm::layout_spec curio() {
    concat_hmm::layout_spec s;
    s.name = "router one-sided native capture test";
    auto add = [&](const char* id, const char* seq, bool fixed, const char* klass, int length) {
        concat_hmm::layout_element e;
        e.id = id; e.seq = seq; e.is_static = fixed; e.klass = klass;
        e.order = (int)s.elements.size() + 1;
        if (length) e.length_candidates = {length};
        s.elements.push_back(e);
    };
    add("barcode_1", "", false, "barcode", 8);
    add("forw_primer", "TCTTCAGCGTTCCCGAGA", true, "forw_primer", 18);
    add("barcode_2", "", false, "barcode", 6);
    add("umi", "", false, "umi", 7);
    add("poly_t", "T{8,}+", true, "poly_tail", 8);
    add("read", "", false, "read", 0);
    return s;
}

std::string dna(int count, unsigned seed) {
    std::string s;
    for (int i = 0; i < count; ++i) {
        seed = seed * 1664525u + 1013904223u;
        s += "ACGT"[(seed >> 25) & 3];
    }
    return s;
}

std::string unit(int length, unsigned seed, bool reverse) {
    std::string s = dna(8, seed) + "TCTTCAGCGTTCCCGAGA" + dna(6, seed + 1) +
                    dna(7, seed + 2) + std::string(20, 'T') + dna(length, seed + 3);
    return reverse ? concat_hmm::detail::revcomp(s) : s;
}

read_streaming::sequence read(const std::string& id, const std::string& seq) {
    std::string qual;
    for (size_t i = 0; i < seq.size(); ++i) qual += char('!' + i % 60);
    return {id, "original comment", seq, "original plus metadata", qual, true};
}

void calibration_cancellation() {
    ReadLayout layout;
    layout.prep_new_layout(std::string(RAD_SOURCE_DIR) + "/resources/read_layout/three_prime_read_layout.csv", false);
    concat_hmm_cal cal(layout, 2);
    std::vector<read_streaming::sequence> chunk(2000);
    for (auto& value : chunk) value.seq = dna(3000, 901);
    const auto calling_thread = std::this_thread::get_id();
    int polls = 0;
    bool caught = false;
    try {
        cal.add_chunk(chunk, chunk.size(), [&]() {
            require(std::this_thread::get_id() == calling_thread, "Host callback ran on a calibration worker");
            if (++polls == 2) throw std::runtime_error("cancel calibration");
        });
    } catch (const std::runtime_error& error) {
        caught = std::string(error.what()) == "cancel calibration";
    }
    require(caught && layout.concat_params.empty(), "Cancelled calibration published partial parameters");
}

void chance_calibration_threshold() {
    const auto stats = Misalignment_Setup::chance_hits_from_histogram({100, 0, 0, 10, 20, 50, 20}, 0.05);
    require(stats.n_all == 200 && stats.n_low == 100 && stats.n_chance == 100,
            "True-copy peak was not separated from chance hits");
    require(stats.mode == 5 && stats.split == 2 && stats.lower == 2,
            "Chance threshold exceeded the accepted cumulative share");
    require(stats.share_at_lower == 0 && std::abs(stats.share_above_lower - 0.1) < 1e-12,
            "Chance-threshold diagnostic shares disagree");
    const auto empty = Misalignment_Setup::chance_hits_from_histogram({}, 0.05);
    require(empty.n_chance == 0 && empty.lower == -1, "Empty calibration invented observations");
}

void mirrored_parameter_cache() {
    char directory[] = "/tmp/rrad-hmm-cache-XXXXXX";
    require(mkdtemp(directory) != nullptr, "Cannot create HMM cache fixture directory");
    const std::string prefix = std::string(directory) + "/cache";
    try {
        ReadLayout layout;
        layout.prep_new_layout(std::string(RAD_SOURCE_DIR) +
            "/resources/read_layout/nanopore_bulk_rapid_bc_read_layout.csv", false);
        layout.by_direction().erase("reverse");
        layout.generate_position_mapping();
        const auto model = concat_hmm::Model::build(layout.to_concat_layout_spec());
        layout.concat_params = model.export_params();
        bool has_mirrored_rows = false;
        for (const auto& row : layout.concat_params.cells)
            has_mirrored_rows |= layout.by_id().find(row.first) == layout.by_id().end();
        require(has_mirrored_rows, "Cache fixture has no implicit mirrored parameters");
        layout.write_to_csv(prefix, "pos_map");
        ReadLayout restored;
        restored.layout = layout.layout;
        restored.sequencing_type = layout.sequencing_type;
        restored.import_position_map(prefix + "_position_map.csv", false);
        require(restored.concat_params.cells == layout.concat_params.cells,
                "Position-map cache discarded or changed mirrored HMM emissions");
        require(restored.layout.size() == layout.layout.size(),
                "Parameter-only cache rows created layout elements");
        const auto cached_model = concat_hmm::Model::build(restored.to_concat_layout_spec(), &restored.concat_params);
        require(cached_model.param_warnings.empty(), "Roundtrip rejected cached HMM parameters");
        require(cached_model.export_params().cells == layout.concat_params.cells,
                "Cached model did not restore all mirrored HMM parameters");
    } catch (...) {
        std::filesystem::remove_all(directory);
        throw;
    }
    std::filesystem::remove_all(directory);
}

void tests() {
    concat_hmm::model_options opt;
    opt.gate = false;
    opt.foldback_split = false; // deliberately no sequence evidence for the F/R junction
    opt.tau_junction = 1.0f;    // exercises coordinate uncertainty together with posterior abstention
    auto model = std::make_shared<concat_hmm::Model>(concat_hmm::Model::build(curio(), nullptr, opt));
    const auto parent = read("unresolved_F_R", unit(355, 12, false) + dna(83, 17) + unit(517, 22, true));
    concat_hmm::Scratch scratch;
    concat_hmm::Result result;
    concat_hmm::segment(*model, parent.seq.data(), (int)parent.seq.size(), scratch, result);
    require(result.k == 2 && result.cut_flags.size() == 1 &&
            (result.cut_flags[0] & concat_hmm::BOUNDARY_UNRESOLVED), "Fixture did not yield an unresolved junction");
    require(result.flags & concat_hmm::RES_ABSTAIN, "Fixture must also exercise legacy-abstain precedence");

    // This is a router test: the callback stands for downstream barcode assignment. No whitelist or
    // variable extraction is needed to prove that blocked children never reach that callback.
    for (bool legacy : {false, true}) {
        ReadLayout layout;
        layout.concat_model = model;
        layout.concat_abstain_legacy = legacy;
        concat_hmm_counters counters;
        concat_layout_info info;
        info.ctr = &counters;
        int assigned = 0, saved = 0;
        const bool handled = SigString::concat_hmm_route(parent, layout, info,
            [&](SigString&, const read_streaming::sequence&) { ++assigned; }, false,
            [&](const read_streaming::sequence& original) {
                ++saved;
                require(original.id == parent.id && original.comment == parent.comment &&
                        original.seq == parent.seq && original.plus == parent.plus &&
                        original.qual == parent.qual && original.is_fastq,
                        "Unresolved parent was altered");
            });
        require(handled && saved == 1 && assigned == 0, "Unresolved junction leaked into barcode assignment/fallback");
        require(counters.boundary_unresolved_parents == 1 && counters.boundary_unresolved_children == 2,
                "Unresolved routing counters disagree");
        require(counters.legacy_abstain == 0, "Legacy route bypassed unresolved-boundary protection");
    }

    // Actual HMM inference on three native reverse capture units with unequal payloads and an unknown
    // scaffold. A located one-sided edge must reach normal downstream processing, not be quarantined.
    opt.tau_junction = 0.9f;
    model = std::make_shared<concat_hmm::Model>(concat_hmm::Model::build(curio(), nullptr, opt));
    const auto repeated = read("reverse_repeat", unit(275, 100, true) + dna(253, 101) +
                                               unit(410, 200, true) + dna(37, 201) + unit(620, 300, true));
    concat_hmm::segment(*model, repeated.seq.data(), (int)repeated.seq.size(), scratch, result);
    require(result.k == 3, "Reverse-repeat fixture did not produce three captures");
    for (uint8_t f : result.cut_flags)
        require((f & concat_hmm::BOUNDARY_LEFT_BOUNDED) && !(f & concat_hmm::BOUNDARY_UNRESOLVED),
                "Reverse-repeat safe edge is unresolved");
    ReadLayout layout;
    layout.concat_model = model;
    concat_hmm_counters counters;
    concat_layout_info info;
    info.ctr = &counters;
    int assigned = 0, saved = 0;
    const bool handled = SigString::concat_hmm_route(repeated, layout, info,
        [&](SigString&, const read_streaming::sequence& child) {
            ++assigned;
            require(child.seq.size() == child.qual.size(), "Child sequence/quality lengths differ");
            const size_t pos = repeated.seq.find(child.seq);
            require(pos != std::string::npos && repeated.qual.substr(pos, child.qual.size()) == child.qual,
                    "Child sequence/quality no longer maps to parent");
        }, false, [&](const read_streaming::sequence&) { ++saved; });
    require(handled && assigned == 3 && saved == 0, "One-sided repeat was withheld or sent to legacy");
    require(counters.boundary_unresolved_parents == 0, "Resolved repeat counted as unresolved");
    require(!result.boundary_unresolved_seen, "Unresolved parent status leaked into next read");

    // The only unknown junction is after the retained cap. Dropping the per-cut array must not
    // erase the obligation to preserve the raw parent or permit a legacy re-extraction.
    opt.k_cap = 2;
    model = std::make_shared<concat_hmm::Model>(concat_hmm::Model::build(curio(), nullptr, opt));
    const auto capped = read("unknown_after_cap", unit(275, 410, false) + unit(410, 510, false) +
                                                 unit(355, 610, false) + unit(517, 710, true));
    concat_hmm::segment(*model, capped.seq.data(), (int)capped.seq.size(), scratch, result);
    require(result.k == 2 && (result.flags & concat_hmm::RES_TOO_MANY), "Cap fixture did not exceed cap");
    require(result.boundary_unresolved_seen, "Cap erased unresolved-parent status");
    for (uint8_t f : result.cut_flags)
        require(!(f & concat_hmm::BOUNDARY_UNRESOLVED), "Cap fixture uncertainty must be beyond retained cuts");
    for (bool legacy : {false, true}) {
        ReadLayout capped_layout;
        capped_layout.concat_model = model;
        capped_layout.concat_abstain_legacy = legacy;
        concat_hmm_counters capped_counters;
        concat_layout_info capped_info;
        capped_info.ctr = &capped_counters;
        int written = 0, preserved = 0;
        const bool cap_handled = SigString::concat_hmm_route(capped, capped_layout, capped_info,
            [&](SigString&, const read_streaming::sequence&) { ++written; }, false,
            [&](const read_streaming::sequence& original) {
                ++preserved;
                require(original.seq == capped.seq && original.qual == capped.qual, "Capped parent changed");
            });
        require(cap_handled && preserved == 1 && written == 0, "Cap allowed unresolved parent to escape preservation");
        require(capped_counters.legacy_too_many == 0, "Cap reached legacy despite an unresolved boundary");
    }
}
} // namespace

int main() {
    try { tests(); calibration_cancellation(); chance_calibration_threshold(); mirrored_parameter_cache(); std::cout << "concat_hmm router regression tests passed\n"; return 0; }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
