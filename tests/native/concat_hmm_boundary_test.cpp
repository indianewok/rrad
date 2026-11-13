// Adapted from RAD dev 6167540; package namespace and cancellation regression.
#include "rad/concat_hmm.hpp"
#include <iostream>
#include <stdexcept>
#include <utility>

using namespace rrad_rad_core::concat_hmm;

using namespace rrad_rad_core;

namespace {
void require(bool value, const std::string& message) {
    if (!value) throw std::runtime_error(message);
}

layout_spec curio_spec() {
    layout_spec s;
    s.name = "Curio one-sided capture boundary regression";
    auto add = [&](const char* id, const char* seq, bool fixed, const char* klass, int length) {
        layout_element e;
        e.id = id; e.seq = seq; e.is_static = fixed; e.klass = klass;
        e.direction = 'F'; e.order = (int)s.elements.size() + 1;
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

std::string dna(int n, unsigned seed) {
    std::string s;
    for (int i = 0; i < n; ++i) {
        seed = 1664525u * seed + 1013904223u;
        s += "ACGT"[(seed >> 25) & 3];
    }
    return s;
}

int state(const Model& m, char strand, bool open) {
    for (int t : m.main_tpl) if (m.tpl[t].strand == strand)
        return open ? m.open_state[t] : m.close_state[t];
    throw std::runtime_error("Missing template");
}

detail::Event event(const Model& m, int st, int begin, int end) {
    detail::Event e{};
    e.start = begin; e.end = end; e.type = (uint8_t)m.st_type[st];
    e.cls = e.type < m.A ? detail::CL_S4 : detail::CL_POLY;
    e.ed = e.type < m.A ? 0 : -1;
    return e;
}

// Keep the former runtime walk as an independent oracle: in particular, an unbounded edge
// retains only the lengths encountered before its first insert, even when more elements follow.
void check_cached_edges(const Model& m, const std::string& context) {
    require(m.head_edge.size() == (size_t)m.S && m.tail_edge.size() == (size_t)m.S,
            context + ": missing cached states");
    for (int st = 0; st < m.S; ++st) for (bool head : {false, true}) {
        const auto& t = m.tpl[m.st_tpl[st]];
        const int element = m.st_elem[st];
        int expected_lo = 0, expected_hi = 0;
        bool expected_bounded = true;
        for (int i = head ? 0 : element + 1, end = head ? element : (int)t.el.size(); i < end; ++i) {
            if (t.el[i].kind == detail::EK_INSERT) { expected_bounded = false; break; }
            expected_lo += t.el[i].lmin;
            expected_hi += t.el[i].lmax;
        }
        const auto& actual = detail::template_edge_span(m, st, head);
        require(actual.bounded == expected_bounded && actual.lo == expected_lo && actual.hi == expected_hi,
                context + ": cached range differs from original walk at state " + std::to_string(st) +
                (head ? " head" : " tail"));
    }
}

struct Fixture {
    std::string seq;
    Scratch path;
    std::vector<std::pair<int, int>> inserts;
};

// The path is supplied from the known physical construct, independently of HMM detection. These tests
// isolate coordinate derivation: a correctly detected capture chain must never delete its payload.
Fixture known_path(const Model& m, const std::string& strands, const std::vector<int>& lengths,
                   const std::vector<int>& unknown_gaps = {}) {
    Fixture f;
    for (size_t i = 0; i < strands.size(); ++i) {
        if (i && !unknown_gaps.empty()) f.seq += dna(unknown_gaps[i - 1], 700u + (unsigned)i);
        const int begin = (int)f.seq.size();
        const std::string payload = dna(lengths[i], 901u + (unsigned)i);
        const std::string capture = dna(8, 100u + (unsigned)i) + "TCTTCAGCGTTCCCGAGA" +
                                    dna(6, 200u + (unsigned)i) + dna(7, 300u + (unsigned)i) + std::string(20, 'T');
        std::string unit = capture + payload;
        const bool forward = strands[i] == 'F';
        if (!forward) unit = detail::revcomp(unit);
        f.seq += unit;
        const int a = state(m, strands[i], true), b = state(m, strands[i], false);
        const int open_start = begin + (forward ? 8 : lengths[i]);
        const int open_len = forward ? 18 : 20;
        const int close_start = begin + (forward ? 39 : lengths[i] + 33);
        const int close_len = forward ? 20 : 18;
        f.path.ev.push_back(event(m, a, open_start, open_start + open_len));
        f.path.ev.push_back(event(m, b, close_start, close_start + close_len));
        f.path.path_ev.push_back((int)f.path.ev.size() - 2);
        f.path.path_ev.push_back((int)f.path.ev.size() - 1);
        f.path.path_st.push_back(a); f.path.path_st.push_back(b);
        f.path.path_cross.push_back(1); f.path.path_cross.push_back(0);
        f.inserts.push_back({begin + (forward ? 59 : 0), begin + (forward ? 59 : 0) + lengths[i]});
    }
    return f;
}

Result derive(const Model& m, Fixture& f) {
    Result r;
    detail::derive_result(m, f.seq.data(), f.path.ev.data(), (int)f.path.ev.size(), (int)f.seq.size(),
                          f.path, (int)f.inserts.size(), false, r);
    require(r.k == (int)f.inserts.size(), "Known path changed unit count");
    require(r.cut_flags.size() == r.cuts.size(), "Boundary flags do not match cuts");
    for (int i = 0; i < r.k; ++i) {
        const int lo = i ? r.cut_lo[i - 1] : 0;
        const int hi = i + 1 < r.k ? r.cut_hi[i] : (int)f.seq.size();
        require(lo <= f.inserts[i].first && hi >= f.inserts[i].second,
                "Known insert was clipped in unit " + std::to_string(i));
    }
    return r;
}

void tests() {
    model_options opts;
    opts.foldback_split = false; // ambiguous F/R gaps must remain unresolved without sequence evidence
    opts.posteriors = false;
    const Model m = Model::build(curio_spec(), nullptr, opts);
    check_cached_edges(m, "Curio main and artifact templates");
    const int ro = state(m, 'R', true), rc = state(m, 'R', false);
    require((m.st_flags[ro] & detail::SF_OPEN) && (m.st_flags[rc] & detail::SF_CLOSE), "Topology flags changed");
    require(!detail::template_edge_span(m, ro, true).bounded, "Reverse insert mistaken for fixed prefix");
    require(detail::template_edge_span(m, rc, false).lo == 8, "Reverse barcode suffix lost");
    uint8_t kind = 0, flags = 0;
    int lo = -1, hi = -1;
    const auto a = event(m, rc, 42, 60), b = event(m, ro, 822, 844);
    const int cut = detail::junction_cut(m, a, rc, b, ro, kind, INT32_MIN, INT32_MIN, &flags, &lo, &hi);
    require(cut == 68 && kind == CUT_GEOMETRY, "Real R2C2 regression still cuts at midpoint 441");
    require((flags & BOUNDARY_LEFT_BOUNDED) && (flags & BOUNDARY_RIGHT_INSERT) &&
            !(flags & BOUNDARY_UNRESOLVED), "One-sided boundary status incorrect");
    require(lo <= 68 && hi >= 68 && lo <= 319, "R2C2 supported payload start clipped");

    // The four boundary cases have distinct physical meanings even though all use the same anchors.
    const int fo = state(m, 'F', true), fc = state(m, 'F', false);
    struct BoundaryCase { int left, right, cut, lo, hi; uint8_t kind, flags; };
    const BoundaryCase cases[] = {
        {rc, fo, 441, 441, 441, CUT_GEOMETRY, BOUNDARY_LEFT_BOUNDED | BOUNDARY_RIGHT_BOUNDED},
        {rc, ro, 68, 66, 70, CUT_GEOMETRY, BOUNDARY_LEFT_BOUNDED | BOUNDARY_RIGHT_INSERT},
        {fc, fo, 814, 812, 816, CUT_GEOMETRY, BOUNDARY_LEFT_INSERT | BOUNDARY_RIGHT_BOUNDED},
        {fc, ro, 441, 60, 822, CUT_MIDPOINT,
         BOUNDARY_LEFT_INSERT | BOUNDARY_RIGHT_INSERT | BOUNDARY_UNRESOLVED}
    };
    for (const auto& expected : cases) {
        const int actual = detail::junction_cut(m, event(m, expected.left, 42, 60), expected.left,
                                               event(m, expected.right, 822, 844), expected.right,
                                               kind, INT32_MIN, INT32_MIN, &flags, &lo, &hi);
        require(actual == expected.cut && lo == expected.lo && hi == expected.hi &&
                kind == expected.kind && flags == expected.flags, "Physical edge decision table changed");
    }

    for (const std::string strands : {"F", "R", "FFF", "RRR", "RFR", "FRF"}) {
        std::vector<int> lengths;
        for (size_t i = 0; i < strands.size(); ++i) lengths.push_back(275 + (int)i * 167);
        for (int gap : {0, 37, 253}) {
            auto f = known_path(m, strands, lengths, std::vector<int>(strands.size() - 1, gap));
            const Result r = derive(m, f);
            for (size_t i = 0; i + 1 < strands.size(); ++i) {
                const bool unknown = strands[i] == 'F' && strands[i + 1] == 'R';
                require(bool(r.cut_flags[i] & BOUNDARY_UNRESOLVED) == unknown, "Unknown junction misclassified");
                if (unknown) {
                    require(r.cut_lo[i] <= f.inserts[i].first && r.cut_hi[i] >= f.inserts[i + 1].second,
                            "Unlocated boundary failed to retain complete insert interval");
                }
            }
            if (strands == "F") require(r.segs[0].end == (int)f.seq.size(), "Terminal forward insert excluded");
            if (strands == "R") require(r.segs[0].start == 0, "Terminal reverse insert excluded");
        }
    }

    auto forward = known_path(m, "FFF", {275, 410, 620}, {37, 253});
    auto reverse = known_path(m, "RRR", {620, 410, 275}, {253, 37});
    reverse.seq = detail::revcomp(forward.seq);
    const Result fr = derive(m, forward), rr = derive(m, reverse);
    const int total = (int)forward.seq.size();
    for (size_t i = 0; i < fr.cuts.size(); ++i) {
        const size_t j = fr.cuts.size() - 1 - i;
        require(rr.cut_lo[j] == total - fr.cut_hi[i] && rr.cut_hi[j] == total - fr.cut_lo[i],
                "Reverse complement changes retained boundaries");
    }

    // Missing terminal observations do not turn a read boundary into an internal insert boundary.
    for (const std::string strands : {"FFF", "RRR"}) {
        auto f = known_path(m, strands, {275, 410, 620});
        f.path.path_ev.erase(f.path.path_ev.begin());
        f.path.path_st.erase(f.path.path_st.begin());
        f.path.path_cross.erase(f.path.path_cross.begin());
        f.path.path_cross.front() = 1;
        f.path.path_ev.pop_back(); f.path.path_st.pop_back(); f.path.path_cross.pop_back();
        const Result partial = derive(m, f);
        require(partial.segs.front().start <= f.inserts.front().first &&
                partial.segs.back().end >= f.inserts.back().second, "Partial terminal observations clip payload");
    }

    // Actual inference through the gate and full decoder must retain the terminal insert too.
    for (bool gate_enabled : {false, true}) {
        model_options inference_opts = opts; inference_opts.gate = gate_enabled;
        const Model im = Model::build(curio_spec(), nullptr, inference_opts);
        check_cached_edges(im, "Gate inference model");
        Scratch inference_scratch;
        Result inference_result;
        for (const std::string strand : {"F", "R"}) {
            auto f = known_path(im, strand, {410});
            segment(im, f.seq.data(), (int)f.seq.size(), inference_scratch, inference_result);
            require(inference_result.k == 1, "Clean singleton inference changed count");
            require(inference_result.segs[0].start <= f.inserts[0].first &&
                    inference_result.segs[0].end >= f.inserts[0].second, "Singleton inference excluded its insert");
            require(inference_result.cut_flags.empty(), "Singleton retains stale boundary flags");
        }
    }

    // Length candidates are a range, not an exact midpoint offset.
    auto variable = curio_spec(); variable.elements[0].length_candidates = {6, 8, 12};
    const Model vm = Model::build(variable, nullptr, opts);
    check_cached_edges(vm, "Variable barcode lengths");
    const int vro = state(vm, 'R', true), vrc = state(vm, 'R', false);
    detail::junction_cut(vm, event(vm, vrc, 42, 60), vrc, event(vm, vro, 822, 844), vro,
                         kind, INT32_MIN, INT32_MIN, &flags, &lo, &hi);
    require(lo <= 66 && hi >= 72, "Variable barcode extent collapsed to nominal offset");

    // Two-sided Stereo-style terminal anchors still bound a non-insert junction.
    auto two_sided = curio_spec();
    layout_element terminal;
    terminal.id = "rev_primer"; terminal.seq = "CCCGCCTCTCAGTACGTCAGCAG";
    terminal.is_static = true; terminal.klass = "rev_primer"; terminal.direction = 'F'; terminal.order = 7;
    two_sided.elements.push_back(terminal);
    layout_element opening = terminal;
    opening.id = "outer_primer"; opening.seq = "CTTCCGATCTATGGCGACCTTATCAG"; opening.klass = "forw_primer"; opening.order = 0;
    two_sided.elements.insert(two_sided.elements.begin(), opening);
    const Model tm = Model::build(two_sided, nullptr, opts);
    check_cached_edges(tm, "Two-sided main, partial and artifact templates");
    const int tc = state(tm, 'F', false), to = state(tm, 'F', true);
    const int midpoint = detail::junction_cut(tm, event(tm, tc, 80, 100), tc, event(tm, to, 120, 145), to,
                                              kind, INT32_MIN, INT32_MIN, &flags, &lo, &hi);
    require(midpoint == 110 && kind == CUT_BOTH_ADAPTERS && !(flags & BOUNDARY_UNRESOLVED),
            "Two-sided terminal adapter behavior changed");

    // Noisy facing terminal-anchor coordinates may overlap by up to five bases.
    // Keep both complete anchor spans in the child windows in either orientation.
    for (char strand : {'F', 'R'}) {
        const int closing = state(tm, strand, false), opening_state = state(tm, strand, true);
        for (int overlap : {0, 1, 2, 3, 4, 5, 6, 12}) {
            const auto left = event(tm, closing, 74, 100);
            const auto right = event(tm, opening_state, 100 - overlap, 126 - overlap);
            const int actual = detail::junction_cut(tm, left, closing, right, opening_state,
                                                     kind, INT32_MIN, INT32_MIN, &flags, &lo, &hi);
            require(bool(flags & BOUNDARY_UNRESOLVED) == (overlap > 5),
                    "Facing terminal-anchor overlap tolerance differs from five bases");
            require((flags & BOUNDARY_LEFT_BOUNDED) && (flags & BOUNDARY_RIGHT_BOUNDED),
                    "Facing terminal-anchor bounds were lost");
            require(kind == CUT_BOTH_ADAPTERS && right.start <= actual && actual <= left.end,
                    "Facing terminal-anchor cut left its observed interval");
            require(lo <= right.start && hi >= left.end,
                    "Overlapping terminal-anchor windows clip a complete facing anchor");
        }
    }

    // A small overlap does not waive the independent observed-interval check.
    // Here a bounded four-base prefix places the nominal cut before both events.
    Model outside = tm;
    outside.head_edge[to] = detail::edge_span{true, 4, 4};
    outside.bc_head[to] = -1;
    outside.bc_tail[tc] = -1;
    const int clamped = detail::junction_cut(outside, event(outside, tc, 80, 100), tc,
                                             event(outside, to, 100, 126), to,
                                             kind, INT32_MIN, INT32_MIN, &flags, &lo, &hi);
    require((flags & BOUNDARY_LEFT_BOUNDED) && (flags & BOUNDARY_RIGHT_BOUNDED) &&
            (flags & BOUNDARY_UNRESOLVED) && clamped == 100,
            "Small overlap bypassed an out-of-interval boundary rejection");

    // Calibration copies and refinalizes models; exported/reloaded parameters must not change layout geometry.
    Model copy = vm;
    check_cached_edges(copy, "Copied model");
    copy.par = copy.prior;
    copy.finalize();
    check_cached_edges(copy, "Prior reset");
    const Params params = vm.export_params();
    copy.apply_params(params);
    copy.finalize();
    copy.finalize();
    check_cached_edges(copy, "Parameter reload and repeated finalize");
    Model moved = std::move(copy);
    check_cached_edges(moved, "Moved model");
    const Model restored = Model::build(variable, &params, opts);
    check_cached_edges(restored, "Parameter round trip during build");

    // Host cancellation during hard-EM must propagate on the caller thread.
    Calibrator calibrator(m);
    auto calibration_read = known_path(m, "F", {410});
    calibrator.add_read(calibration_read.seq.data(), static_cast<int>(calibration_read.seq.size()));
    bool cancelled = false;
    int polls = 0;
    try {
        calibrator.finalize(20, 0.001, [&]() {
            ++polls;
            throw std::runtime_error("host cancellation");
        });
    } catch (const std::runtime_error& error) {
        cancelled = std::string(error.what()) == "host cancellation";
    }
    require(cancelled && polls == 1, "HMM calibration swallowed caller cancellation");

    // All bundled layouts exercise short/long spacers, different anchor classes and both strand templates.
    for (const char* file : {"nanopore_bulk_rapid_bc_read_layout.csv", "splitseq_read_layout.csv",
                            "visium_hd_read_layout.csv", "three_prime_read_layout.csv",
                            "sctagger_sim_read_layout.csv", "visium_three_prime_read_layout.csv",
                            "five_prime_read_layout.csv"}) {
        const auto spec = parse_layout_csv(std::string(RAD_SOURCE_DIR) + "/resources/read_layout/" + file);
        check_cached_edges(Model::build(spec, nullptr, opts), file);
    }
}
} // namespace

int main() {
    try { tests(); std::cout << "concat_hmm boundary regression tests passed\n"; return 0; }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
