// concat_hmm.hpp -- layout-generic, self-calibrating concatemer segmenter for long reads (C++17, STL only).
// Per read: stage 0 seeds and poly runs, stage 1 gate (clean single constructs), stage 2 windowed Myers events,
// stage 3 event HSMM (Viterbi + forward-backward). Calibrated parameters live in hmm_* position-map columns.
// Thread safety: Model is immutable after build() and shareable; Scratch is per thread.
#pragma once
#define CONCAT_HMM_HAS_CUT_WINDOWS 1  // Result::cut_lo / cut_hi
#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#if defined(__aarch64__) && defined(__ARM_NEON)
#include <arm_neon.h>
#define CONCAT_HMM_NEON 1
#elif defined(__SSE4_1__)
#include <smmintrin.h>
#define CONCAT_HMM_SSE41 1
#endif

namespace rrad_rad_core {
namespace concat_hmm {

struct layout_element {
    std::string id;
    std::string seq;                     // static sequence ('T{12,}+' style for poly); empty for variable elements
    bool is_static = false;
    std::string klass;
    char direction = 'F';                // 'F' or 'R'
    int order = 0;
    std::vector<int> length_candidates;
    int max_edits = -1;                  // from misalign_lower; -1 if unknown
};
struct layout_spec {
    std::vector<layout_element> elements;
    bool bulk = false;
    std::string name;
    std::string mode;                    // text after "Read Layout:" (e.g. "bulk"), lower case
};

// Calibrated parameters as cached in RAD's position map: element id -> hmm_* column -> numbers.
// Structural tables apply only when `topology` matches the model's topology hash.
struct Params {
    std::map<std::string, std::map<std::string, std::vector<double>>> cells;
    bool guard_failed = false;
    uint64_t topology = 0;
    std::vector<std::string> bad_cells;  // "row/column" of hmm_* cells that did not parse as numbers
    bool empty() const { return cells.empty() && bad_cells.empty(); }
};

enum class Status : uint8_t { FULL, PARTIAL, TRUNCATED_AT_READ_END, MISSING_EXPECTED };
inline const char* status_name(Status s) {
    switch (s) {
        case Status::FULL: return "FULL";
        case Status::PARTIAL: return "PARTIAL";
        case Status::TRUNCATED_AT_READ_END: return "TRUNC";
        default: return "MISSING";
    }
}

struct check_region {
    int start, end;                  // 0-based half-open window (padded)
    int element;                     // index into layout_spec::elements
    char strand;                     // 'F' / 'R' strand of the construct ('T' / 'D' for artifact templates)
    int construct;                   // construct index (0-based, read order)
    Status status;
    int retained_from, retained_to;  // adapter offsets present [from, to), -1 if n/a
    int edits;                       // edit distance of the seen part, -1 if unknown
    float conf;
};
enum : uint16_t {
    SEG_PARTIAL_LEFT = 1,     // opening element of the construct not observed
    SEG_PARTIAL_RIGHT = 2,    // closing element not observed
    SEG_LOWCONF = 4,
    SEG_SAME_MOLECULE_R = 8,  // the junction to the right is a fold-back
    SEG_DOUBLE_BC = 16,       // D geometry with no clear cDNA direction (strand 'D': RAD writes nothing)
    SEG_ARTIFACT = 32,        // artifact template (barcode-less T, or a single-primer end)
    SEG_UNANCHORED_R = 64,    // cut to the right placed without a junction adapter
    SEG_EMPTY = 128,          // (near) zero-length insert
    SEG_D_KEPT = 256,         // D geometry: the unit at the cDNA 5' end is kept; the segment ends at the other unit's inner anchor
    SEG_D_SPLIT = 512         // D geometry: one half of a split at a cDNA strand change
};
struct Segment {
    int start, end;  // 0-based half-open
    char strand;     // 'F', 'R', 'T' (barcode-less artifact), 'D' (double barcode unit), '?' (no strand evidence)
    uint16_t flags;
    float conf;
};
enum : uint32_t {
    RES_ABSTAIN = 1,
    RES_TOO_MANY = 2,         // more than model_options::k_cap constructs (k capped)
    RES_GUARD_FAILED = 4,     // layout priors in use
    RES_FAST_PATH = 8,        // decided by the stage-1 gate
    RES_OPPOSITE_EVIDENCE = 16, // unexplained opposite-strand anchor -> strand '?'
    RES_RESCUED = 32            // a junction was added by the facing-anchor rescue
};
struct Result {
    int k = 0;
    std::vector<int> cuts;           // size k-1
    std::vector<Segment> segs;       // size k
    std::vector<check_region> checks;
    float p_single = 0.f;            // P(k == 1 | read)
    char strand_call = '?';          // 'F','R','M','?'
    uint32_t flags = 0;
    bool boundary_unresolved_seen = false; // survives k_cap truncation of per-junction arrays
    std::vector<float> cut_post;     // per-cut junction posterior (size k-1)
    // per-cut child windows: the left construct ends at cut_hi[i], the right one starts at cut_lo[i]
    // (cut_lo[i] <= cuts[i] <= cut_hi[i]; overlap retains uncertain sequence as well as barcode blocks)
    std::vector<int> cut_lo, cut_hi;
    std::vector<uint8_t> cut_kind;   // per cut: CUT_* value
    // Boundary geometry is independent of cut_post, which measures a junction somewhere between events.
    std::vector<uint8_t> cut_flags;  // per cut: BOUNDARY_* values
    std::vector<float> cut_aux;      // per cut: fold-back Jaccard or strand-flip evidence (strand_mu units), else 0
#ifdef CONCAT_HMM_STAGE_TIMERS
    uint64_t stage_ns[4] = {0, 0, 0, 0};  // stage0, gate, stage2 (Myers/events), stage3 (decode+output)
#endif
};
#define CONCAT_HMM_HAS_CUT_KIND 1      // Result::cut_kind
#define CONCAT_HMM_HAS_BOUNDARY_FLAGS 1 // Result::cut_flags
#define CONCAT_HMM_HAS_STRAND_TABLE 1  // calibrated strand k-mer table (hmm_strand column), strand-flip cuts
enum : uint8_t { CUT_BOTH_ADAPTERS = 0, CUT_ONE_ADAPTER = 1, CUT_GEOMETRY = 2, CUT_MIDPOINT = 3, CUT_FOLDBACK = 4, CUT_STRAND_FLIP = 5 };
enum : uint8_t {
    BOUNDARY_LEFT_BOUNDED = 1, BOUNDARY_RIGHT_BOUNDED = 2,
    BOUNDARY_LEFT_INSERT = 4, BOUNDARY_RIGHT_INSERT = 8,
    BOUNDARY_UNRESOLVED = 16, BOUNDARY_SEQUENCE_RESOLVED = 32
};

struct model_options {
    int artifact_templates = -1;   // -1 auto (when the layout has a weak outer anchor), 0 off, 1 on
    int d_template = -1;           // -1 auto, 0 off, 1 on
    bool gate = true;
    bool windowed_myers = true;    // false: whole-read Myers on the full path
    bool posteriors = true;
    bool check_regions = true;
    bool relaxed_openers = true;   // ED <= budget+2 opener scan after strong closers
    bool foldback_split = true;
    int k_cap = 256;
    int min_interior_len = -1;     // -1: layout-derived
    int min_terminal_len = 60;
    float tau_single = 0.99f;      // abstain when k == 1 and p_single < tau_single (full path only)
    float tau_junction = 0.9f;     // abstain when a junction posterior < tau_junction
    int check_pad = 25;            // bp padding of MISSING_EXPECTED windows
    int seen_pad = 6;              // bp padding of windows of seen elements
};

struct layout_spec_error : std::runtime_error {
    explicit layout_spec_error(const std::string& m) : std::runtime_error(m) {}
};

namespace detail {

constexpr float NEG = -1e30f;
constexpr float NEG_HALF = -5e29f;
constexpr int NFEAT = 32;
constexpr int FEAT_PREFIX = 16;    // 16..20: read-end prefix of a closing anchor (retained-length bins)
constexpr int FEAT_SEEDPART = 21;  // 21..25: seed-chain partial, retained-length bins
constexpr int FEAT_TRUNC5 = 26;    // opener suffix at the read start (seeds)
constexpr int FEAT_TRUNC3S = 27;   // closer running off the read end (seeds)
constexpr int NPOLYBIN = 8;
constexpr int NPARTBIN = 5;
constexpr int MAXPAT = 16;         // Myers patterns (chunks), both strands
constexpr int MAXANCH = 16;        // distinct static elements, both strands
constexpr int MAXPOLY = 4;
constexpr int MAXSTATE = 63;
constexpr int PRED_WINDOW = 10;
constexpr int MAX_OVERLAP = 12;
constexpr int INS_LO = -40;
constexpr int INS_HI = 16384;
constexpr int NINSBIN = 84;
constexpr int POLY_MIN = 12;
constexpr int PARTIAL_MIN = 9;
constexpr int JUNC_SPACER_MAX = 300;
constexpr int KGAP = 6;
constexpr int MAX_EVENTS = 1 << 15;
constexpr int PACK_MAX = 57;
constexpr int PAT_MAXLEN = 63;
constexpr int SEED_K = 11;
constexpr float SHORT_LLR_CAP = 8.0f;  // nats; LLR cap of short anchors
constexpr float FORCE_BONUS = 30.0f;   // nats; bonus of a junction forced by the facing-anchor rescue
constexpr int RESCUE_GMAX = 150;       // bp; rescue gap beyond the abutting barcode blocks, one anchor on the path
constexpr int RESCUE_GMAX_OFF = 30;    // bp; same, both anchors off the path
constexpr int RESCUE_GMIN = -5;        // bp; overlap of the two barcode blocks
constexpr int DEL_GMIN = -5;           // bp; overlap allowed between abutting barcode blocks (deleted primers)
constexpr int DEL_GMAX_CERT = 30;      // bp; filler up to which facing partner anchors certify the junction
constexpr int BC_WIN_PAD = 2;          // bp; child windows reach an anchor-derived barcode-block edge + this
constexpr int FOLD_ROOM = 25;          // fold_mirrored only
constexpr int FOLD_SPREAD = 6;         // fold_mirrored only
constexpr float WEAK_POST = 0.5f;      // posterior cap of a pairing junction without enough evidence
constexpr int FOLD_K = 10;             // fold-back decision: k-mer length
constexpr int FOLD_ARM_MIN = 40;       // bp; minimum arm length
constexpr float FOLD_J_LO = 0.06f;     // minimum arm Jaccard
constexpr float FOLD_NW_ID = 0.65f;    // banded identity of the mirrored core
constexpr int FOLD_REACH = 40;         // bp; mirrored diagonal reach to an arm end bounded by an inner anchor
constexpr int FOLD_REACH_OPEN = 80;    // bp; same, when that inner anchor was not observed
constexpr int FOLD_DIAG_MIN = 8;       // shared k-mers on the mirrored diagonal
constexpr int STRAND_K = 6;            // strand table k-mer length
constexpr int STRAND_N = 4096;
constexpr int STRAND_MARGIN_BP = 10;   // calibration: bases skipped next to the cDNA's bounding elements
constexpr int STRAND_MIN_CDNA = 60;    // calibration: minimum cDNA length counted
constexpr int FLIP_ARM_MIN = 30;       // bp; strand-flip cut: minimum arm length
constexpr float FLIP_MARGIN = 0.30f;   // strand-flip cut: minimum arm mean, as a fraction of strand_mu
constexpr float FLIP_DISTINCT_MIN = 0.5f; // distinct k-mers / k-mer positions per arm (low-complexity guard)
constexpr int FLIP_EVIDENCE_MIN = 40;  // strand-flip cut: minimum arm score sum, in strand_mu units
constexpr int STRAND_D_K = 9;          // D strand table k-mer length
constexpr int STRAND_D_N = 1 << (2 * STRAND_D_K);
constexpr double STRAND_D_MIN_KMERS = 2.0e6;  // calibration: minimum counted k-mers to write the D table
constexpr float D_EVIDENCE_MIN = 16.f; // D strand rule: minimum evidence, in strand_d_mu units

inline int poly_bin(int len) {
    static const int ub[7] = {14, 17, 21, 26, 33, 42, 55};
    for (int b = 0; b < 7; ++b)
        if (len <= ub[b]) return b;
    return 7;
}
inline int prefix_bin(int len) {
    if (len <= 10) return 0;
    if (len <= 12) return 1;
    if (len <= 15) return 2;
    if (len <= 18) return 3;
    return 4;
}
inline int seedpart_bin(int len) {
    if (len <= 12) return 0;
    if (len <= 15) return 1;
    if (len <= 18) return 2;
    if (len <= 21) return 3;
    return 4;
}
inline int ins_bin(int x) {
    if (x < 0) x = 0;
    int b = (int)std::floor(8.0 * std::log2((double)x + 16.0)) - 32;
    return std::min(std::max(b, 0), NINSBIN - 1);
}
inline double ins_bin_lo(int b) { return std::exp2((b + 32) / 8.0) - 16.0; }
inline int ins_bin_count(int b) {
    if (b == NINSBIN - 1) return 1 << 20;
    int lo = (int)std::ceil(ins_bin_lo(b)), hi = (int)std::ceil(ins_bin_lo(b + 1));
    return std::max(1, hi - lo);
}
inline char comp(char c) {
    switch (c) {
        case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A';
        default: return 'N';
    }
}
inline std::string revcomp(const std::string& s) {
    std::string o(s.rbegin(), s.rend());
    for (auto& c : o) c = comp(c);
    return o;
}
inline std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && (unsigned char)s[a] <= ' ') ++a;
    while (b > a && (unsigned char)s[b - 1] <= ' ') --b;
    return s.substr(a, b - a);
}
inline std::string lower(std::string s) {
    for (auto& c : s) c = (char)tolower((unsigned char)c);
    return s;
}
inline std::string upper(std::string s) {
    for (auto& c : s) c = (char)toupper((unsigned char)c);
    return s;
}
inline std::vector<std::string> csv_fields(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool q = false;
    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (q) {
            if (c == '"') {
                if (i + 1 < line.size() && line[i + 1] == '"') { cur += '"'; ++i; }
                else q = false;
            } else cur += c;
        } else if (c == '"') q = true;
        else if (c == ',') { out.push_back(trim(cur)); cur.clear(); }
        else if (c != '\r') cur += c;
    }
    out.push_back(trim(cur));
    return out;
}
inline std::string csv_quote(const std::string& s) {
    bool need = s.find_first_of(",\"\n") != std::string::npos;
    if (!need) return s;
    std::string o = "\"";
    for (char c : s) { if (c == '"') o += '"'; o += c; }
    return o + "\"";
}
// "15-16" -> {15,16}; "22" -> {22}; "8,10" / "8;10" / "8 10" -> {8,10}; "" (or na / nan / none / null) -> {}.
// `ok` is false on anything else (letters, negatives, open or reversed ranges, lengths above LENGTH_MAX).
constexpr int LENGTH_MAX = 4096;
inline std::vector<int> parse_length_list(const std::string& s_in, bool* ok_out = nullptr) {
    std::vector<int> out;
    bool ok = true;
    const std::string s = lower(trim(s_in));
    if (s.empty() || s == "na" || s == "nan" || s == "none" || s == "null") { if (ok_out) *ok_out = true; return out; }
    size_t i = 0;
    auto num = [&](long& v) {
        if (i >= s.size() || s[i] < '0' || s[i] > '9') return false;
        v = 0;
        while (i < s.size() && s[i] >= '0' && s[i] <= '9') { v = std::min<long>(v * 10 + (s[i] - '0'), 1L << 40); ++i; }
        return true;
    };
    auto skip_sep = [&]() { while (i < s.size() && (s[i] == ',' || s[i] == ';' || s[i] == ' ' || s[i] == '|')) ++i; };
    skip_sep();
    while (ok && i < s.size()) {
        long a = 0, b = 0;
        if (!num(a)) { ok = false; break; }
        b = a;
        if (i < s.size() && s[i] == '-') {
            ++i;
            if (!num(b) || b < a) { ok = false; break; }
        }
        if (b > LENGTH_MAX) { ok = false; break; }
        for (long x = a; x <= b; ++x) out.push_back((int)x);
        if (i < s.size() && !(s[i] == ',' || s[i] == ';' || s[i] == ' ' || s[i] == '|')) { ok = false; break; }
        skip_sep();
    }
    if (!ok) out.clear();
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    if (ok_out) *ok_out = ok;
    return out;
}
inline std::string read_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return "";
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}
inline bool is_poly_elem(const layout_element& e) {
    if (!e.is_static) return false;
    std::string k = lower(e.klass), id = lower(e.id);
    if (k == "poly_tail" || k == "poly_t" || k == "poly_a" || id == "poly_t" || id == "poly_a" || id.rfind("poly_", 0) == 0) return true;
    return e.seq.find('{') != std::string::npos;
}
inline char poly_base_of(const layout_element& e) {
    std::string id = lower(e.id), k = lower(e.klass);
    if (id == "poly_a" || k == "poly_a") return 'A';
    if (id == "poly_t" || k == "poly_t") return 'T';
    for (char c : upper(e.seq))
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T') return c;
    return 'T';
}
// Static sequence in model form: upper case, IUPAC codes -> N, other characters dropped.
inline std::string static_seq_of(const layout_element& e) {
    std::string s2;
    for (char c : upper(e.seq)) {
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N') s2 += c;
        else if (c != '\0' && std::strchr("RYSWKMBDHV", c)) s2 += 'N';
    }
    return s2;
}
inline std::string strip_flank_n(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && s[a] == 'N') ++a;
    while (b > a && s[b - 1] == 'N') --b;
    return s.substr(a, b - a);
}
inline std::string anchor_seq_of(const layout_element& e) { return strip_flank_n(static_seq_of(e)); }
inline int informative_bases(const std::string& s) {
    int c = 0;
    for (char x : s) c += x != 'N';
    return c;
}
inline bool is_insert_elem(const layout_element& e) {
    std::string k = lower(e.klass), id = lower(e.id);
    return !e.is_static && (k == "read" || id == "read" || id == "rc_read");
}
constexpr int ANCHOR_MIN_INFORMATIVE = 12;
// Throws layout_spec_error unless every static anchor has >= ANCHOR_MIN_INFORMATIVE non-N bases, spacer
// lengths lie in [0, LENGTH_MAX] and at least one anchor exists.
inline void validate_spec(const layout_spec& L) {
    int nstatic = 0;
    for (auto& e : L.elements) {
        if (!e.is_static) {
            if (!is_insert_elem(e))
                for (int x : e.length_candidates)
                    if (x < 0 || x > LENGTH_MAX)
                        throw layout_spec_error("concat_hmm: layout '" + L.name + "': element '" + e.id + "' has length candidate " + std::to_string(x) +
                                              " outside [0, " + std::to_string(LENGTH_MAX) + "]");
            continue;
        }
        if (is_poly_elem(e)) continue;
        const std::string kl = lower(e.klass);
        if (kl == "start" || kl == "stop") continue;
        const std::string a = anchor_seq_of(e);
        if (a.empty()) continue;  // static row without a sequence: not an anchor
        if (informative_bases(a) < ANCHOR_MIN_INFORMATIVE)
            throw layout_spec_error("concat_hmm: layout '" + L.name + "': static element '" + e.id + "' ('" + e.seq + "') has " +
                                  std::to_string(informative_bases(a)) + " informative bases (minimum " + std::to_string(ANCHOR_MIN_INFORMATIVE) + ")");
        ++nstatic;
    }
    if (nstatic == 0)
        throw layout_spec_error("concat_hmm: layout '" + L.name + "' has no usable static (anchor) element: the model would be empty");
}

}  // namespace detail

// Parses both RAD layout CSV formats; throws layout_spec_error when no usable static element remains.
inline layout_spec parse_layout_text(const std::string& text_in, const std::string& name = "") {
    using namespace detail;
    std::string text = text_in;
    if (text.size() >= 3 && (unsigned char)text[0] == 0xEF && (unsigned char)text[1] == 0xBB && (unsigned char)text[2] == 0xBF)
        text = text.substr(3);
    std::vector<std::string> lines;
    {
        std::stringstream ss(text);
        std::string ln;
        while (std::getline(ss, ln)) {
            if (ln.size() >= 3 && (unsigned char)ln[0] == 0xEF && (unsigned char)ln[1] == 0xBB && (unsigned char)ln[2] == 0xBF) ln = ln.substr(3);
            if (!trim(ln).empty()) lines.push_back(ln);
        }
    }
    layout_spec L;
    L.name = name;
    if (lines.size() < 2) throw layout_spec_error("concat_hmm: layout '" + name + "' has fewer than 2 non-empty lines");
    size_t hdr_line = std::string::npos;
    for (size_t li = 0; li < std::min<size_t>(lines.size(), 3); ++li) {
        auto f = csv_fields(lines[li]);
        if (!f.empty() && lower(f[0]).rfind("read layout", 0) == 0) {
            auto p = f[0].find(':');
            if (p != std::string::npos) L.mode = lower(trim(f[0].substr(p + 1)));
            continue;
        }
        bool has_id = false, has_seq = false, has_type = false;
        for (auto& x : f) {
            std::string lx = lower(x);
            has_id |= lx == "id";
            has_seq |= lx == "seq";
            has_type |= lx == "type";
        }
        if (has_id && (has_seq || has_type)) { hdr_line = li; break; }
    }
    if (hdr_line == std::string::npos) throw layout_spec_error("concat_hmm: layout '" + name + "': no header row with id/seq/type columns");
    L.bulk = L.mode == "bulk";
    auto hdr = csv_fields(lines[hdr_line]);
    std::map<std::string, int> col;
    for (size_t i = 0; i < hdr.size(); ++i) col[lower(hdr[i])] = (int)i;
    auto get = [&](const std::vector<std::string>& f, const char* k) -> std::string {
        auto it = col.find(k);
        if (it == col.end() || it->second >= (int)f.size()) return "";
        return f[it->second];
    };
    const bool has_dir = col.count("direction") > 0;
    int order = 0;
    for (size_t li = hdr_line + 1; li < lines.size(); ++li) {
        auto f = csv_fields(lines[li]);
        layout_element e;
        e.id = get(f, "id");
        if (e.id.empty()) continue;
        std::string type = lower(get(f, "type"));
        e.is_static = type == "static";
        std::string cls = get(f, "class");
        e.klass = cls.empty() ? e.id : cls;
        e.seq = get(f, "seq");
        std::string dir = has_dir ? lower(get(f, "direction")) : "forward";
        e.direction = dir.rfind("reverse", 0) == 0 ? 'R' : 'F';
        std::string ord = get(f, "order");
        e.order = ord.empty() ? ++order : std::atoi(ord.c_str());
        order = std::max(order, e.order);
        std::string lk = lower(e.klass);
        if (lk == "poly_t" || lk == "poly_a") e.klass = "poly_tail";
        // an unparsable length is an error for spacers only
        const bool spacer = !e.is_static && !is_insert_elem(e);
        bool ok1 = true, ok2 = true;
        const std::string lc = get(f, "length_candidates"), el = get(f, "expected_length");
        e.length_candidates = parse_length_list(lc, &ok1);
        if (e.length_candidates.empty()) e.length_candidates = parse_length_list(el, &ok2);
        if (spacer && (!ok1 || (e.length_candidates.empty() && !ok2)))
            throw layout_spec_error("concat_hmm: layout '" + name + "': element '" + e.id + "' has an unusable length ('" + (ok1 ? el : lc) +
                                  "'; expected e.g. 12, 15-16 or 8,10, each <= " + std::to_string(LENGTH_MAX) + ")");
        std::string ml = get(f, "misalign_lower");
        if (!ml.empty()) e.max_edits = std::atoi(ml.c_str());
        L.elements.push_back(e);
    }
    validate_spec(L);
    return L;
}
inline layout_spec parse_layout_csv(const std::string& path) {
    std::string t = detail::read_file(path);
    if (t.empty()) throw layout_spec_error("concat_hmm: cannot read layout file '" + path + "'");
    return parse_layout_text(t, path);
}

struct Model;
struct Scratch;


namespace detail {

enum elem_kind : uint8_t { EK_ANCHOR = 0, EK_POLY = 1, EK_SPACER = 2, EK_INSERT = 3 };

struct TElem {
    elem_kind kind = EK_SPACER;
    int anc = -1;            // anchor id (EK_ANCHOR) or poly id (EK_POLY)
    int lmin = 0, lmax = 0;
    std::string label;       // element id (merged spacers: "a+b")
    int spec = -1;           // layout_spec element index (first element of a merged spacer block)
    int fwd_spec = -1;       // derived reverse-complement slot: layout_spec index of its forward counterpart
    double nominal = 0, var = 0, hr = 0;  // nominal length, variance, half-range
};
struct Template {
    char strand = '?';       // F, R, T (barcode-less), D (double barcode unit)
    bool main = false;
    bool art = false;        // artifact / auxiliary template (T, A5, A5r, A3, A3r, D)
    bool dtpl = false;
    std::string name;
    std::vector<TElem> el;
    std::vector<int> obs;              // element indices of observable slots
    std::vector<double> pnom, phr;     // prefix sums of nominal / half-range over el (size el+1)
    std::vector<int> pins;             // prefix count of inserts
};
struct Anchor {
    std::string seq;         // full static sequence (flanking N stripped)
    int m = 0, k = 0, strong_ed = 0;
    uint16_t pblk[4] = {0, 0, 0, 0}, sblk[4] = {0, 0, 0, 0};  // first / last disjoint 5-mers (slot-completion pre-filter)
    int nblk = 0;
    int chunk0 = 0, nchunk = 1;
    bool closing = false, opening = false, weak = false, shrt = false;
    std::string row;         // serialization row (element id)
    int spec = -1;
};
struct Pattern {             // one Myers pattern (an anchor chunk of <= 63 bp)
    std::string used;
    int anchor = -1, off = 0, m = 0, k = 0;
    std::array<uint8_t, 65> inf{};  // inf[i]: informative (non-N) bases among used[0, i)
};
struct PWord {               // SWAR-packed word: A in [0,m_a), guard bit, B in [m_a+1, m_a+1+m_b)
    uint64_t peq[256];
    int pa = -1, pb = -1, m_a = 0, m_b = 0, t_a = 0, t_b = 0, k_a = 0, k_b = 0;
    uint64_t CG = ~0ULL, HB = 0, CS = ~0ULL, T = 0, H = 0, S0 = 0, FA = ~0ULL;
};
struct LWord {
    uint64_t peq[256];
    int pat = -1, m = 0, k = 0;
};
struct PState { uint64_t Pv, Mv, S; };
struct Clu { int32_t first, last, best_end, best; };
inline void clu_reset(Clu& c) { c.first = c.last = c.best_end = -1; c.best = 1 << 30; }
struct raw_ev { int32_t end, first, last; uint8_t pat, score; };

enum : uint8_t { CL_S4 = 0, CL_S6 = 1, CL_S8 = 2, CL_SP = 3, CL_POLY = 4 };
enum : uint8_t { EVF_TRUNC3 = 1, EVF_TRUNC5 = 2, EVF_SEEDPART = 4, EVF_RELAXED = 8, EVF_PREFIX = 16 };
struct Event {
    int32_t start, end;      // implied full-element span, clamped to the read
    uint8_t type;            // anchor id, or A + poly id
    uint8_t feat;
    uint8_t cls;             // evidence class (CL_*)
    uint8_t fl;              // EVF_*
    int16_t rlo, rhi;        // retained element offsets [rlo, rhi) (anchors), -1 if n/a
    int16_t ed;              // edit distance of the seen part (-1 unknown)
    int16_t len;             // poly length / retained length
};
struct seed_cluster {
    int32_t s, e;            // implied element span (may run past the read ends)
    uint8_t a, n, nx;        // anchor, seeds, exact seeds
    int16_t olo, ohi;        // min / max seed offset (element coordinates)
    bool valid() const { return nx >= 1 || n >= 2; }
};
struct poly_reg { int32_t s, e; uint8_t q; };  // q = poly type id

struct desc_hot {
    float base = 0;
    int32_t fixed = 0;
    int16_t tid = -1;
    int16_t slack = 0;       // certification slack around `fixed` (tight descriptors)
    uint8_t kind = 0;        // 0 tight, 1 one insert, 2 two inserts
    uint8_t valid = 0;
    uint8_t pal = 0;         // cross into an opening slot across a palindromic junction class
    uint8_t spacer_only = 0; // same-construct tight descriptor spanning only spacers
    // lengths of deletable barcode-adjacent primers skipped by a tight junction (closing side, opening side):
    // the gap may be shorter than `fixed` by dela, delb or dela + delb
    int16_t dela = 0, delb = 0;
};
struct desc_cold {
    std::vector<std::pair<int, bool>> pres;
    int jfrom = -1, jto = -1;
    std::string key;
    double sd = 0;
};
using PMap = std::map<std::string, std::vector<double>>;

enum : uint16_t {
    SF_ANCHOR = 1, SF_OPEN = 2, SF_CLOSE = 4, SF_CERT = 8, SF_WEAK = 16, SF_SHORT = 32, SF_MAIN = 64, SF_ART = 128,
    SF_DTPL = 256
};

struct Nbr {                 // tight neighbour relation used for windows / gate
    int other;               // other type (anchor id or A + poly id)
    int dir;                 // +1: other lies after this element, -1: before
    int lo, hi;              // gap range (this.end -> other.start for dir +1; other.end -> this.start for -1)
};
// Barcode block: a barcode-adjacent outer primer joined by spacers only to an inner partner anchor.
// Head rule: template start [prim] <block> [inner]; tail rule: [inner] <block> [prim] template end.
struct bc_rule { int prim, inner, block; };

// Physical prefix/suffix length range up to the first insert. An insert leaves the partial
// range intact but makes the edge unbounded; observable OPEN/CLOSE flags do not describe it.
struct edge_span {
    bool bounded = true;
    int lo = 0, hi = 0;
};

inline int pat_k_budget(int m) {
    int k = m < 12 ? 1 : m < 16 ? 2 : m < 20 ? 4 : m < 30 ? 6 : m < 45 ? 8 : 10;
    return std::min(k, m / 2);
}
inline int strong_ed_of(int m) { return std::max(0, std::min(10, (int)std::floor(m * 4.0 / 22.0 + 1e-9))); }

}  // namespace detail

struct Model {
    layout_spec spec;
    model_options opt;
    std::vector<detail::Template> tpl;
    std::vector<detail::Anchor> anch;
    std::vector<detail::Pattern> pats;
    std::vector<char> poly_base;
    std::vector<std::string> poly_row;
    int A = 0, Q = 0, NTYPE = 0, S = 0, NT = 0;
    std::vector<int> st_tpl, st_obs, st_type, st_elem;
    std::vector<uint16_t> st_flags;
    std::vector<std::vector<int>> states_by_type;
    std::vector<int> open_state, close_state;
    std::vector<int> head_fix, tail_fix;   // nominal non-insert prefix/suffix bp (-1: insert between); boundary retention uses layout ranges
    std::vector<detail::edge_span> head_edge, tail_edge;  // full layout ranges, fixed once during model construction
    // per state: bp of the barcode block between the slot and the barcode-adjacent primer at the template start /
    // end, primer counted as absent; -1 if the slot is not such an inner slot
    std::vector<int> bc_head, bc_tail;
    // per state: the slot belongs to a barcode-side unit (an inner slot or its barcode-adjacent primer)
    std::vector<uint8_t> bc_unit;
    std::vector<detail::bc_rule> head_rules, tail_rules;  // anchor partners only (poly runs never certify)
    std::vector<int> tight_dela, tight_delb;             // per tight table: deletable barcode-adjacent primer lengths
    std::vector<uint8_t> inner_role;                     // per type: 1 / 2 tail / head rule inner anchor, 4 / 8 its primer
    std::vector<std::vector<detail::Nbr>> nbr;
    std::vector<int> main_tpl;             // indices of F and R
    std::vector<uint8_t> tpl_fold;         // per template: fold-back test applies
    std::vector<std::vector<int>> tpl_mirror;  // per template, per observable slot: mirrored slot across the insert, -1 none
    std::vector<detail::PWord> words, rwords;
    std::vector<detail::LWord> lwords, single, single_rev;  // long patterns; single-pattern words per anchor (first chunk) and reversed
    std::vector<int> rword_anchor_pat;
    bool words_shared = false, rwords_shared = false;
    int Wov = 0, max_close_m = 0;
    int K = detail::SEED_K;
    uint64_t kmask = 0;
    std::vector<uint64_t> seed_bits, seed_hkey, seed_bloom;  // exact 4^K presence bits; label hash keys; 2^19-bit L1 pre-filter
    std::vector<uint16_t> seed_hval;
    uint64_t seed_hmask = 0;
    bool stage0_poly[2] = {false, false};  // detect A / T runs in stage 0
    int polyq_of_base[4] = {-1, -1, -1, -1};  // A C G T -> poly id
    bool other_poly = false;
    detail::PMap par, prior;
    std::vector<std::string> pres_names;
    std::vector<detail::desc_hot> same, cross;
    std::vector<detail::desc_cold> same_c, cross_c;
    std::vector<std::string> tight_keys;
    std::vector<int> tight_lo, tight_hi;
    std::vector<std::string> tight_row;    // serialization row of each tight table ("" = global)
    std::vector<float> emit, emit_pal;     // S * NFEAT
    std::vector<float> emit_bonus;         // S * NFEAT: uncapped - capped LLR of short anchors
    std::vector<std::vector<float>> tight;
    std::vector<float> ins[2];
    std::vector<float> lbegin, lend;
    float lp0 = 0, lp1 = 0;
    int ins1_p99 = 1 << 30;                // p99 of the (calibrated) one-insert length density
    // cDNA strand table: per STRAND_K-mer, sense vs antisense log-odds; empty when not calibrated
    std::vector<float> strand_lo;
    float strand_mu = 0.f;                 // mean per-position sense score of the calibration cDNA
    bool has_strand() const { return strand_lo.size() == (size_t)detail::STRAND_N && strand_mu > 0.f; }
    // D strand table: per STRAND_D_K-mer log-odds from the hmm_strand_d counts; empty when not calibrated
    std::vector<float> strand_d_lo;
    float strand_d_mu = 0.f;
    bool has_strand_d() const { return strand_d_lo.size() == (size_t)detail::STRAND_D_N && strand_d_mu > 0.f; }
    int zone5 = 150, zone3 = 120;
    float fast_p = 0.995f;
    int min_interior = 60, min_terminal = 60;
    bool guard_failed = false;
    std::vector<std::string> param_warnings;  // hmm_* cells rejected by apply_params
    uint64_t topo_hash = 0;
    bool finalized = false;

    static Model build(const layout_spec& spec, const Params* p = nullptr, const model_options& o = model_options());
    void finalize();
    void apply_params(const Params& p);
    Params export_params() const;
    std::string describe() const;
    std::string type_name(int ty) const {
        if (ty < A) return anch[ty].seq;
        return std::string("poly") + poly_base[ty - A];
    }
    std::string type_label(int ty) const {
        if (ty < A) return anch[ty].row;
        return poly_row[ty - A];
    }
};

struct Scratch {
    std::vector<uint64_t> seeds;            // (position << 16) | seed label
    std::vector<uint32_t> hits;
    std::vector<detail::seed_cluster> sc;
    std::vector<detail::poly_reg> poly;
    std::vector<detail::raw_ev> raw_a, raw_b, raw_x, raw_ch;
    std::vector<std::pair<int, int>> win;
    std::string cat;                        // evidence windows concatenated (separated by non-matching bytes)
    std::vector<int> cat_off;
    std::vector<detail::Event> ev, ev2;
    std::vector<float> dp;
    std::vector<int32_t> bp, cs;
    std::vector<float> fw, bw, scale;
    struct Tr { int32_t src, dst; float w; uint8_t cross, ha; };  // valid transitions of the last Viterbi
    std::vector<Tr> tr;
    std::vector<float> trm;
    std::vector<float> jpost;               // per event gap junction posterior
    std::vector<int> path_ev, path_st;
    std::vector<uint8_t> path_cross;
    std::vector<int> cidx;                  // per construct: first path index
    std::vector<uint8_t> ckind;
    std::vector<int> ca, cb, bcut;
    std::vector<int> jlo, jhi, segj;        // per path junction: child window bounds; per output segment: its right junction
    std::vector<uint8_t> jflags;            // per path junction: boundary geometry, separate from junction posterior
    std::vector<uint8_t> cweak;             // per pairing junction: bit 1 weak (abstain), 2 strong pairing, 4 / 8 partner-delimited
    std::vector<int> cgap;                  // per pairing junction: gap between the barcode-block edges
    std::vector<float> cpost, spost;
    std::vector<float> caux, saux;          // per path junction / per output segment: cut score
    std::vector<uint8_t> skind;             // per output segment: kind of the split inside its construct, 0 none
    std::vector<int> scut;                  // per output segment: the cut on its right
    std::vector<int> slot_ev;
    std::vector<std::pair<int, int>> force; // event pairs whose junction the next Viterbi forces
    std::vector<int> onp;                   // per event: index on the decoded path (-1: off the path)
    std::vector<uint32_t> fb_keys;
    std::vector<int32_t> fb_pos, fb_diag;
    std::vector<int> fb_hist;
    std::vector<uint32_t> jk_a, jk_b;       // fold-back decision: 10-mer sets of the two arms
    std::vector<std::pair<uint32_t, int32_t>> jp_a, jp_b;  // ... k-mers with positions (mirror reach)
    std::vector<int32_t> nw0, nw1;          // banded edit-distance rows (fold-back tie-break)
    std::vector<float> flip_ps;             // strand-flip cut: prefix sums of per-position strand scores
    std::vector<uint8_t> flip_seen;         // ... distinct k-mer table of one arm (low-complexity guard)
    std::string shuf;
    float best_score = 0, null_score = 0, l_z = 0, l_z1 = 0;
    Scratch() {
        ev.reserve(256); ev2.reserve(256); raw_a.reserve(64); raw_b.reserve(64);
        dp.reserve(4096); bp.reserve(4096); cs.reserve(4096);
        path_ev.reserve(64); path_st.reserve(64); path_cross.reserve(64);
        seeds.reserve(256); sc.reserve(64); poly.reserve(32); win.reserve(32);
    }
};

namespace detail {

inline int add_anchor(Model& M, const std::string& seq_in, const std::string& row, int spec) {
    const std::string seq = strip_flank_n(seq_in);
    for (int i = 0; i < (int)M.anch.size(); ++i)
        if (M.anch[i].seq == seq) return i;
    Anchor x;
    x.seq = seq;
    x.m = (int)seq.size();
    x.k = pat_k_budget(x.m);
    x.strong_ed = std::min(x.k, strong_ed_of(x.m));
    x.row = row;
    x.spec = spec;
    x.shrt = x.m < 16;
    M.anch.push_back(x);
    return (int)M.anch.size() - 1;
}
inline int add_poly(Model& M, char b, const std::string& row) {
    for (int i = 0; i < (int)M.poly_base.size(); ++i)
        if (M.poly_base[i] == b) return i;
    M.poly_base.push_back(b);
    M.poly_row.push_back(row);
    return (int)M.poly_base.size() - 1;
}
inline void set_nominal(const Model& M, TElem& e) {
    if (e.kind == EK_ANCHOR) {
        double m = M.anch[e.anc].m;
        e.nominal = m; e.var = (0.06 * m) * (0.06 * m) + 0.5; e.hr = 1;
    } else if (e.kind == EK_POLY) {
        e.nominal = 25; e.var = 15.0 * 15.0; e.hr = 15;
    } else if (e.kind == EK_SPACER) {
        e.nominal = 0.5 * (e.lmin + e.lmax);
        double h = 0.5 * (e.lmax - e.lmin);
        e.var = h * h + (0.04 * e.nominal) * (0.04 * e.nominal) + 1.0;
        e.hr = h + 1;
    } else { e.nominal = 0; e.var = 0; e.hr = 0; }
}
inline void finish_template(const Model& M, Template& T) {
    T.obs.clear();
    T.pnom.assign(T.el.size() + 1, 0.0);
    T.phr.assign(T.el.size() + 1, 0.0);
    T.pins.assign(T.el.size() + 1, 0);
    for (int i = 0; i < (int)T.el.size(); ++i) {
        set_nominal(M, T.el[i]);
        if (T.el[i].kind == EK_ANCHOR || T.el[i].kind == EK_POLY) T.obs.push_back(i);
        T.pnom[i + 1] = T.pnom[i] + T.el[i].nominal;
        T.phr[i + 1] = T.phr[i] + T.el[i].hr;
        T.pins[i + 1] = T.pins[i] + (T.el[i].kind == EK_INSERT);
    }
}

// SF_OPEN / SF_CLOSE describe observable order, not physical template ends. An insert can lie
// outside the first/last observable (for example Curio's poly-A / reverse linker). Keep that
// decoder topology unchanged and use the full layout whenever coordinates are being placed.
inline const edge_span& template_edge_span(const Model& M, int st, bool head) {
    return head ? M.head_edge[st] : M.tail_edge[st];
}
inline bool at_template_edge(const Model& M, int st, bool head) {
    const edge_span r = template_edge_span(M, st, head);
    return r.bounded && r.hi == 0;
}
inline int retained_template_start(const Model& M, int st, int event_start) {
    const edge_span r = template_edge_span(M, st, true);
    return r.bounded ? std::max(0, event_start - r.hi) : 0;
}
inline int retained_template_end(const Model& M, int st, int event_end, int read_len) {
    const edge_span r = template_edge_span(M, st, false);
    return r.bounded ? std::min(read_len, event_end + r.hi) : read_len;
}

// Build a template's element chain from layout elements (one direction, in read order).
inline Template chain_from_spec(Model& M, const layout_spec& L, const std::vector<int>& idx, char strand, const std::string& name) {
    Template T;
    T.strand = strand; T.name = name; T.main = true;
    for (int ii : idx) {
        const layout_element& e = L.elements[ii];
        std::string kl = lower(e.klass);
        if (kl == "start" || kl == "stop") continue;
        TElem t;
        t.label = e.id;
        t.spec = ii;
        if (e.is_static) {
            if (is_poly_elem(e)) {
                t.kind = EK_POLY;
                t.anc = add_poly(M, poly_base_of(e), e.id);
                t.lmin = POLY_MIN; t.lmax = 1000;
            } else {
                std::string s = static_seq_of(e);
                if (s.empty()) continue;
                t.kind = EK_ANCHOR;
                t.anc = add_anchor(M, s, e.id, ii);
                t.lmin = t.lmax = M.anch[t.anc].m;
            }
        } else if (is_insert_elem(e)) {
            t.kind = EK_INSERT; t.lmin = 0; t.lmax = 1 << 20;
        } else {
            t.kind = EK_SPACER;
            if (!e.length_candidates.empty()) { t.lmin = e.length_candidates.front(); t.lmax = e.length_candidates.back(); }
            else { t.lmin = 1; t.lmax = 30; }
            if (!T.el.empty() && T.el.back().kind == EK_SPACER) {
                T.el.back().label += "+" + e.id;
                T.el.back().lmin += t.lmin;
                T.el.back().lmax += t.lmax;
                continue;
            }
        }
        T.el.push_back(t);
    }
    return T;
}
inline TElem rc_elem(Model& M, const TElem& e0) {
    TElem t = e0;
    if (t.kind == EK_ANCHOR) {
        std::string row = e0.label.rfind("rc_", 0) == 0 ? e0.label.substr(3) : "rc_" + e0.label;
        t.anc = add_anchor(M, revcomp(M.anch[e0.anc].seq), row, -1);
    }
    if (t.kind == EK_POLY) {
        char b = comp(M.poly_base[e0.anc]);
        t.anc = add_poly(M, b, b == 'A' ? "poly_a" : b == 'T' ? "poly_t" : std::string("poly_") + b);
    }
    t.label = t.label.rfind("rc_", 0) == 0 ? t.label.substr(3) : "rc_" + t.label;
    t.fwd_spec = e0.spec >= 0 ? e0.spec : e0.fwd_spec;
    t.spec = -1;
    return t;
}
// Binds a slot to its layout element; sequences are compared in model form (anchor_seq_of).
inline void bind_spec(const Model& M, TElem& t) {
    if (t.spec >= 0 || (t.kind != EK_ANCHOR && t.kind != EK_POLY)) return;
    for (int i = 0; i < (int)M.spec.elements.size(); ++i) {
        const auto& e = M.spec.elements[i];
        if (t.kind == EK_ANCHOR && e.is_static && !is_poly_elem(e) && anchor_seq_of(e) == M.anch[t.anc].seq) { t.spec = i; return; }
        if (t.kind == EK_POLY && is_poly_elem(e) && poly_base_of(e) == M.poly_base[t.anc]) { t.spec = i; return; }
    }
    // user-format layouts carry no reverse rows: fall back to the forward counterpart
    for (int i = 0; i < (int)M.spec.elements.size(); ++i) {
        const auto& e = M.spec.elements[i];
        if (t.kind == EK_ANCHOR && e.is_static && !is_poly_elem(e) && anchor_seq_of(e) == revcomp(M.anch[t.anc].seq)) { t.spec = i; return; }
        if (t.kind == EK_POLY && is_poly_elem(e) && poly_base_of(e) == comp(M.poly_base[t.anc])) { t.spec = i; return; }
    }
    // last resort: the forward element the slot was derived from (same template position)
    if (t.fwd_spec >= 0 && t.fwd_spec < (int)M.spec.elements.size()) t.spec = t.fwd_spec;
}

inline uint64_t fnv(uint64_t h, const std::string& s) {
    for (unsigned char c : s) h = (h ^ c) * 1099511628211ULL;
    return (h ^ 0xFF) * 1099511628211ULL;
}

}  // namespace detail

namespace detail {

inline void fill_peq(uint64_t* peq, const std::string& p, int off) {
    for (int i = 0; i < (int)p.size(); ++i) {
        uint64_t bit = 1ULL << (off + i);
        if (p[i] == 'N') {
            for (const char* c = "ACGTacgtNn"; *c; ++c) peq[(unsigned char)*c] |= bit;
        } else {
            peq[(unsigned char)p[i]] |= bit;
            peq[(unsigned char)tolower(p[i])] |= bit;
        }
    }
}
inline uint8_t code2(unsigned char c) { return (uint8_t)(((c >> 1) ^ (c >> 2)) & 3); }  // A0 C1 G2 T3 (either case)
inline int informative(const std::string& s) {
    int c = 0;
    for (char x : s) c += x != 'N';
    return c;
}
// Split anchors into Myers patterns: anchors <= 63 bp are one pattern; longer ones are split at
// N runs >= 4 and into <= 63-bp chunks; chunks with < 12 informative bases are not searched.
inline void build_patterns(Model& M) {
    M.pats.clear();
    for (int a = 0; a < (int)M.anch.size(); ++a) {
        Anchor& x = M.anch[a];
        x.chunk0 = (int)M.pats.size();
        int inf = informative(x.seq);
        x.k = pat_k_budget(inf);
        x.strong_ed = std::min(x.k, strong_ed_of(inf));
        if (x.m <= PAT_MAXLEN) {
            Pattern p;
            p.used = x.seq; p.anchor = a; p.off = 0; p.m = x.m; p.k = pat_k_budget(inf);
            M.pats.push_back(p);
        } else {
            std::vector<std::pair<int, int>> reg;
            int i = 0;
            while (i < x.m) {
                int j = i;
                while (j < x.m) {
                    if (x.seq[j] == 'N') {
                        int r = j;
                        while (r < x.m && x.seq[r] == 'N') ++r;
                        if (r - j >= 4) break;
                        j = r;
                    } else ++j;
                }
                if (j > i) reg.push_back({i, j});
                while (j < x.m && x.seq[j] == 'N') ++j;
                i = j;
            }
            for (auto& r : reg) {
                int len = r.second - r.first, nch = (len + PAT_MAXLEN - 1) / PAT_MAXLEN, sz = (len + nch - 1) / nch;
                for (int c = 0; c < nch; ++c) {
                    int o = r.first + c * sz, l = std::min(sz, r.second - o);
                    std::string u = x.seq.substr(o, l);
                    if (informative(u) < 12) continue;
                    Pattern p;
                    p.used = u; p.anchor = a; p.off = o; p.m = l; p.k = pat_k_budget(informative(u));
                    M.pats.push_back(p);
                }
            }
        }
        x.nchunk = (int)M.pats.size() - x.chunk0;
        for (int pi = x.chunk0; pi < (int)M.pats.size(); ++pi) {
            Pattern& p = M.pats[pi];
            p.inf[0] = 0;
            for (int i = 0; i < p.m && i < 64; ++i) p.inf[i + 1] = (uint8_t)(p.inf[i] + (p.used[i] != 'N'));
        }
        if (x.nchunk == 0)
            throw layout_spec_error("concat_hmm: static element '" + x.row + "' yields no searchable pattern (every N-free part has < 12 informative bases)");
    }
    if ((int)M.pats.size() > MAXPAT)
        throw layout_spec_error("concat_hmm: layout needs " + std::to_string(M.pats.size()) + " Myers patterns (max " + std::to_string(MAXPAT) + ")");
}
inline void pack_words(const Model& M, const std::vector<int>& pids, const std::vector<int>& kk, std::vector<PWord>& out,
                       std::vector<int>* longs, bool& shared) {
    out.clear();
    std::vector<int> packable;
    for (size_t i = 0; i < pids.size(); ++i) {
        if (M.pats[pids[i]].m <= PACK_MAX) packable.push_back((int)i);
        else if (longs) longs->push_back(pids[i]);
    }
    std::sort(packable.begin(), packable.end(), [&](int a, int b) {
        int ma = M.pats[pids[a]].m, mb = M.pats[pids[b]].m;
        return ma != mb ? ma > mb : a < b;
    });
    std::vector<bool> used(pids.size(), false);
    for (size_t x = 0; x < packable.size(); ++x) {
        int a = packable[x];
        if (used[a]) continue;
        used[a] = true;
        int b = -1;
        for (size_t y = packable.size(); y-- > x + 1;) {
            int c = packable[y];
            if (!used[c] && M.pats[pids[a]].m + M.pats[pids[c]].m <= 56) { b = c; break; }
        }
        PWord W;
        memset(W.peq, 0, sizeof W.peq);
        W.pa = pids[a]; W.m_a = M.pats[pids[a]].m; W.k_a = kk[a]; W.t_a = W.m_a - 1;
        fill_peq(W.peq, M.pats[pids[a]].used, 0);
        if (b >= 0) {
            used[b] = true;
            W.pb = pids[b]; W.m_b = M.pats[pids[b]].m; W.k_b = kk[b];
            fill_peq(W.peq, M.pats[pids[b]].used, W.m_a + 1);
            W.t_b = W.m_a + W.m_b;
            W.CG = ~(1ULL << W.m_a);
            W.CS = ~(1ULL << (W.m_a + 1));
            W.HB = (1ULL << W.t_a) | (1ULL << W.t_b);
            W.S0 = ((uint64_t)W.m_a << W.t_a) + ((uint64_t)W.m_b << W.t_b);
            W.T = ((uint64_t)(W.k_a + 1) << W.t_a) + ((uint64_t)(W.k_b + 1) << W.t_b);
            W.H = (1ULL << (W.t_b - 1)) | (1ULL << 63);
            W.FA = (1ULL << (W.t_b - W.t_a)) - 1;
        } else {
            W.t_b = 64;
            W.CG = ~0ULL; W.CS = ~0ULL;
            W.HB = 1ULL << W.t_a;
            W.S0 = (uint64_t)W.m_a << W.t_a;
            W.T = (uint64_t)(W.k_a + 1) << W.t_a;
            W.H = 1ULL << 63;
            W.FA = ~0ULL;
        }
        out.push_back(W);
    }
    shared = !out.empty();
    for (auto& W : out)
        shared = shared && W.CG == out[0].CG && W.HB == out[0].HB && W.CS == out[0].CS && W.T == out[0].T && W.H == out[0].H;
}
inline void make_lword(const Model& M, int pat, int k, LWord& L) {
    memset(L.peq, 0, sizeof L.peq);
    L.pat = pat; L.m = M.pats[pat].m; L.k = k;
    fill_peq(L.peq, M.pats[pat].used, 0);
}
inline void build_engine(Model& M) {
    std::vector<int> all, kk, longs;
    for (int i = 0; i < (int)M.pats.size(); ++i) { all.push_back(i); kk.push_back(M.pats[i].k); }
    pack_words(M, all, kk, M.words, &longs, M.words_shared);
    M.lwords.clear();
    for (int p : longs) { LWord L; make_lword(M, p, M.pats[p].k, L); M.lwords.push_back(L); }
    // relaxed opener words (budget + 2), unchunked opening anchors only
    std::vector<int> rp, rk;
    for (int a = 0; a < M.A; ++a) {
        const Anchor& x = M.anch[a];
        if (!x.opening || x.nchunk != 1 || x.m > PACK_MAX || x.shrt) continue;
        rp.push_back(x.chunk0);
        rk.push_back(std::min(x.k + 2, 15));
    }
    std::vector<int> dummy;
    pack_words(M, rp, rk, M.rwords, nullptr, M.rwords_shared);
    M.single.assign(M.A, LWord());
    M.single_rev.assign(M.A, LWord());
    for (int a = 0; a < M.A; ++a) {
        const int pc = M.anch[a].chunk0 < (int)M.pats.size() ? M.anch[a].chunk0 : 0;
        make_lword(M, pc, M.anch[a].k, M.single[a]);
        LWord& R = M.single_rev[a];
        memset(R.peq, 0, sizeof R.peq);
        R.pat = pc; R.m = M.pats[pc].m; R.k = M.anch[a].k;
        std::string rv(M.pats[pc].used.rbegin(), M.pats[pc].used.rend());
        fill_peq(R.peq, rv, 0);
    }
    // slot-completion pre-filter (partial_in_window): a prefix (suffix) of i >= 10 bp at D <= floor(i/8)
    // contains one of the first (last) floor((m-1)/8) + 1 disjoint 5-mers of the adapter exactly
    for (int a = 0; a < M.A; ++a) {
        Anchor& x = M.anch[a];
        x.nblk = 0;
        if (informative(x.seq) != x.m || x.m < 16) continue;
        const int nb = std::min(4, (x.m - 1) / 8 + 1);
        if (nb * 5 > x.m) continue;
        auto code5 = [&](int off) { uint16_t v = 0; for (int i = 0; i < 5; ++i) v = (uint16_t)((v << 2) | code2((unsigned char)x.seq[off + i])); return v; };
        for (int b = 0; b < nb; ++b) { x.pblk[b] = code5(5 * b); x.sblk[b] = code5(x.m - 5 * (b + 1)); }
        x.nblk = nb;
    }
    int mx = 0;
    for (auto& p : M.pats) mx = std::max(mx, p.m + p.k + 3);
    M.Wov = mx;
    M.max_close_m = 0;
    for (auto& x : M.anch)
        if (x.closing) M.max_close_m = std::max(M.max_close_m, x.m);
}
inline void build_seed_index(Model& M) {
    const int K = M.K;
    M.kmask = (1ULL << (2 * K)) - 1;
    M.seed_bits.assign(((size_t)1 << (2 * K)) / 64 + 1, 0);
    std::vector<std::pair<uint64_t, uint16_t>> ent;
    std::map<uint64_t, uint16_t> exact;
    for (int pi = 0; pi < (int)M.pats.size(); ++pi) {
        const Pattern& p = M.pats[pi];
        for (int o = 0; o + K <= p.m; ++o) {
            int off = p.off + o;
            if (off > 127) break;
            uint64_t v = 0;
            bool ok = true;
            for (int i = 0; i < K; ++i) {
                char c = p.used[o + i];
                if (c == 'N') { ok = false; break; }
                v = (v << 2) | code2((unsigned char)c);
            }
            if (!ok || exact.count(v)) continue;
            uint16_t lab = (uint16_t)((p.anchor << 8) | off);
            exact[v] = lab;
            M.seed_bits[v >> 6] |= 1ULL << (v & 63);
            ent.push_back({v, lab});
        }
    }
    const size_t nexact = ent.size();
    for (size_t e = 0; e < nexact; ++e) {
        const uint64_t v = ent[e].first;
        for (int pos = 0; pos < K; ++pos) {
            const int sh = 2 * (K - 1 - pos);
            const uint64_t orig = (v >> sh) & 3;
            for (uint64_t b = 0; b < 4; ++b) {
                if (b == orig) continue;
                const uint64_t w = (v & ~(3ULL << sh)) | (b << sh);
                if ((M.seed_bits[w >> 6] >> (w & 63)) & 1) continue;
                M.seed_bits[w >> 6] |= 1ULL << (w & 63);
                ent.push_back({w, (uint16_t)(ent[e].second | 0x80)});
            }
        }
    }
    M.seed_bloom.assign((1u << 19) / 64, 0);
    for (auto& e : ent) {
        const uint64_t b = e.first & ((1u << 19) - 1);
        M.seed_bloom[b >> 6] |= 1ULL << (b & 63);
    }
    size_t hs = 1;
    while (hs < ent.size() * 4) hs <<= 1;
    M.seed_hkey.assign(hs, ~0ULL);
    M.seed_hval.assign(hs, 0);
    M.seed_hmask = hs - 1;
    for (auto& e : ent) {
        uint64_t h = ((e.first * 0x9E3779B97F4A7C15ULL) >> 32) & M.seed_hmask;
        while (M.seed_hkey[h] != ~0ULL && M.seed_hkey[h] != e.first) h = (h + 1) & M.seed_hmask;
        M.seed_hkey[h] = e.first;
        M.seed_hval[h] = e.second;
    }
}
inline uint16_t seed_label(const Model& M, uint64_t km) {
    uint64_t h = ((km * 0x9E3779B97F4A7C15ULL) >> 32) & M.seed_hmask;
    while (M.seed_hkey[h] != ~0ULL) {
        if (M.seed_hkey[h] == km) return M.seed_hval[h];
        h = (h + 1) & M.seed_hmask;
    }
    return 0xFFFF;
}

inline std::vector<double> normal_mix_pmf(int lo, int hi, double sd, double tail) {
    std::vector<double> p(hi - lo + 1);
    double s2 = 3 * sd + 3, tot = 0;
    for (int x = lo; x <= hi; ++x) {
        double v = 0.85 * std::exp(-0.5 * x * x / (sd * sd)) / sd + 0.15 * std::exp(-0.5 * x * x / (s2 * s2)) / s2;
        p[x - lo] = v; tot += v;
    }
    for (auto& v : p) v = v / tot * (1 - tail);
    if (tail > 0)
        for (int x = 0; x <= hi; ++x) p[x - lo] += tail / (hi + 1);
    for (auto& v : p) v = std::max(v, 1e-7);
    return p;
}
// Adds deletion regimes to a tight table whose junction skips barcode-adjacent primers of da / db bp:
// x = gap - fixed near -da, -db, -(da + db). Prior weights: 0.04 per single deletion, 0.02 for both.
inline void add_deletion_regimes(std::vector<double>& p, int lo, int hi, int da, int db) {
    int dd[3] = {0, 0, 0};
    double ww[3] = {0, 0, 0};
    int nr = 0;
    auto add = [&](int d, double w) {
        if (d <= 0) return;
        for (int i = 0; i < nr; ++i)
            if (dd[i] == d) { ww[i] += w; return; }
        dd[nr] = d; ww[nr] = w; ++nr;
    };
    add(da, 0.04);
    add(db, 0.04);
    if (da > 0 && db > 0) add(da + db, 0.02);
    double wsum = 0;
    for (int i = 0; i < nr; ++i) wsum += ww[i];
    for (auto& v : p) v *= (1 - wsum);
    std::vector<double> comp(p.size());
    for (int i = 0; i < nr; ++i) {
        double tot = 0;
        for (int x = lo; x <= hi; ++x) {
            const double y = x + dd[i];
            double v = 0;
            if (y >= DEL_GMIN) v = 0.7 * std::exp(-0.5 * y * y / 9.0) / 3.0 + (y >= 0 ? 0.3 * std::exp(-y / 40.0) / 40.0 : 0.0);
            comp[x - lo] = v;
            tot += v;
        }
        if (tot > 0)
            for (size_t j = 0; j < p.size(); ++j) p[j] += ww[i] * comp[j] / tot;
    }
    for (auto& v : p) v = std::max(v, 1e-7);
}
inline std::vector<double> lognormal_bins(double median, double sigma) {
    std::vector<double> p(NINSBIN);
    double mu = std::log(median), tot = 0;
    for (int b = 0; b < NINSBIN; ++b) {
        double lo = ins_bin_lo(b), hi = b + 1 < NINSBIN ? ins_bin_lo(b + 1) : 1e7;
        double mid = std::max(1.0, 0.5 * (lo + hi));
        double z = (std::log(mid) - mu) / sigma;
        double v = std::exp(-0.5 * z * z) * (hi - lo) / mid + 1e-6;
        p[b] = v; tot += v;
    }
    for (auto& v : p) v /= tot;
    return p;
}
inline std::vector<double> convolve_ins(const std::vector<double>& p1) {
    const int STEP = 8, NX = 16384 / STEP;
    std::vector<double> d(NX, 0.0), out(NINSBIN, 1e-9);
    for (int i = 0; i < NX; ++i) {
        int x = i * STEP + STEP / 2;
        int b = ins_bin(x);
        d[i] = p1[b] / ins_bin_count(b) * STEP;
    }
    for (int i = 0; i < NX; ++i) {
        if (d[i] < 1e-12) continue;
        for (int j = 0; i + j < NX; ++j) {
            if (d[j] < 1e-12) continue;
            int x = (i + j) * STEP + STEP;
            out[ins_bin(x)] += d[i] * d[j];
        }
    }
    double tot = 0;
    for (double v : out) tot += v;
    for (auto& v : out) v /= tot;
    return out;
}
inline std::string tkey(const Model& M, int st) { return M.tpl[M.st_tpl[st]].name + std::to_string(M.st_obs[st]); }

// Read-end truncation priors per retained-length bin: known 10x primers, generic otherwise.
inline void prefix_prior_shape(const std::string& seq, double* w) {
    static const double tsorc[5] = {0.45, 0.36, 0.17, 0.01, 0.01};   // rev_primer
    static const double fprc[5] = {0.64, 0.33, 0.01, 0.01, 0.01};    // rc_forw_primer
    static const double gen[5] = {0.45, 0.33, 0.18, 0.02, 0.02};
    const double* s = gen;
    if (seq.find("ACTCTGCGTTGATACCACTGCTT") != std::string::npos) s = tsorc;
    else if (seq == "AGATCGGAAGAGCGTCGTGTAG") s = fprc;
    for (int b = 0; b < NPARTBIN; ++b) w[b] = s[b];
}

inline void make_priors(Model& M) {
    PMap& pr = M.prior;
    pr.clear();
    pr["strand"] = {};
    for (int ty = 0; ty < M.NTYPE; ++ty) {
        std::vector<double> q(NFEAT, 0.0), nl(NFEAT, 0.0);
        if (ty < M.A) {
            const Anchor& x = M.anch[ty];
            static const double w[11] = {0.62, 0.13, 0.07, 0.05, 0.04, 0.03, 0.03, 0.02, 0.01, 0.01, 0.01};
            for (int d = 0; d <= x.k; ++d) q[d] = w[std::min(d, 10)];
            if (x.opening && !x.shrt) {
                if (x.k + 1 < FEAT_PREFIX) q[x.k + 1] = 0.012;
                if (x.k + 2 < FEAT_PREFIX) q[x.k + 2] = 0.008;
                q[FEAT_TRUNC5] = 0.004;
            }
            if (x.closing) {
                double sh[NPARTBIN];
                prefix_prior_shape(x.seq, sh);
                for (int b = 0; b < NPARTBIN; ++b) q[FEAT_PREFIX + b] = 0.10 * sh[b];
                q[FEAT_TRUNC3S] = 0.01;
            }
            static const double sp[NPARTBIN] = {0.15, 0.25, 0.3, 0.2, 0.1};
            for (int b = 0; b < NPARTBIN; ++b) q[FEAT_SEEDPART + b] = 0.02 * sp[b];
            const int inf = informative(x.seq);
            for (int d = 0; d <= std::min(x.k + 2, 15); ++d) {
                double c = 1;
                for (int i = 0; i < d; ++i) c = c * (inf - i) / (i + 1);
                nl[d] = std::min(5e-2, std::max(1e-9, c * std::pow(3.0, d) * std::pow(4.0, -inf) * inf * 4));
            }
            for (int b = 0; b < NPARTBIN; ++b) { nl[FEAT_PREFIX + b] = 1e-7 * std::pow(0.25, b); nl[FEAT_SEEDPART + b] = 2e-6 * std::pow(0.5, b); }
            nl[FEAT_TRUNC5] = 1e-7;
            nl[FEAT_TRUNC3S] = 1e-7;
        } else {
            static const double w[NPOLYBIN] = {0.08, 0.10, 0.14, 0.18, 0.20, 0.15, 0.10, 0.05};
            static const double n0[NPOLYBIN] = {0.5, 0.22, 0.12, 0.07, 0.04, 0.03, 0.015, 0.005};
            for (int b = 0; b < NPOLYBIN; ++b) { q[b] = w[b]; nl[b] = 3e-5 * n0[b]; }
        }
        double t = 0;
        for (double x : q) t += x;
        for (auto& x : q) x /= t;
        pr["q." + M.type_name(ty)] = q;
        if (ty < M.A) pr["qpal." + M.type_name(ty)] = q;
        pr["null." + M.type_name(ty)] = nl;
    }
    for (size_t i = 0; i < M.pres_names.size(); ++i) {
        const std::string& n = M.pres_names[i];
        double p = 0.95;
        if (n.rfind("int.", 0) == 0) p = n.back() == 'p' ? 0.85 : 0.92;
        else if (n.rfind("pal.close.", 0) == 0) p = 0.90;
        else if (n.rfind("pal.open.", 0) == 0) p = 0.89;
        else p = 0.985;
        pr["pres." + n] = {p};
    }
    const double artw = 0.5;
    for (int t = 0; t < M.NT; ++t) {
        std::vector<double> j(M.NT);
        double tot = 0;
        for (int u = 0; u < M.NT; ++u) {
            const Template& U = M.tpl[u];
            j[u] = U.main ? 1.0 : (U.strand == 'T' && !U.dtpl && U.name == "T") ? 0.05 : 0.5 * artw;
            tot += j[u];
        }
        for (auto& x : j) x /= tot;
        pr["junc." + M.tpl[t].name] = j;
    }
    std::vector<double> b(M.S, 0.0), e(M.S, 0.0);
    double tb = 0, te = 0;
    for (int st = 0; st < M.S; ++st) {
        const Template& T = M.tpl[M.st_tpl[st]];
        const double w = T.main ? 1.0 : (T.name == "T") ? 0.1 : artw;
        bool open = M.st_obs[st] == 0, close = M.st_obs[st] == (int)T.obs.size() - 1;
        b[st] = w * (open ? 20.0 : 1.0);
        e[st] = w * (close ? 10.0 : 1.5);
        tb += b[st]; te += e[st];
    }
    for (auto& x : b) x /= tb;
    for (auto& x : e) x /= te;
    pr["begin"] = b;
    pr["end"] = e;
    pr["p0"] = {0.01};
    pr["ins1"] = lognormal_bins(380, 0.7);
    pr["ins2"] = convolve_ins(pr["ins1"]);
    for (size_t t = 0; t < M.tight_keys.size(); ++t) {
        double sd = 2;
        for (int x = 0; x < M.S * M.S; ++x) {
            if (M.same[x].valid && M.same[x].kind == 0 && M.same[x].tid == (int)t) sd = M.same_c[x].sd;
            if (M.cross[x].valid && M.cross[x].kind == 0 && M.cross[x].tid == (int)t) sd = M.cross_c[x].sd;
        }
        const bool junction = M.tight_keys[t].back() == 'x';
        std::vector<double> pm = normal_mix_pmf(M.tight_lo[t], M.tight_hi[t], sd, junction ? 0.03 : 0.0);
        if (M.tight_dela[t] > 0 || M.tight_delb[t] > 0) add_deletion_regimes(pm, M.tight_lo[t], M.tight_hi[t], M.tight_dela[t], M.tight_delb[t]);
        pr[M.tight_keys[t]] = pm;
    }
    pr["gate"] = {150, 120, 0.995};
    pr["guard"] = {0};
}

inline bool palindromic(const Model& M, int t1, int t2) {
    const Template& A = M.tpl[t1];
    const Template& B = M.tpl[t2];
    if (A.obs.empty() || B.obs.empty()) return false;
    const TElem& c = A.el[A.obs.back()];
    const TElem& o = B.el[B.obs.front()];
    if (c.kind != EK_ANCHOR || o.kind != EK_ANCHOR) return false;
    return M.anch[o.anc].seq == revcomp(M.anch[c.anc].seq);
}
inline void build_descriptors(Model& M) {
    const int S = M.S;
    M.same.assign(S * S, desc_hot());
    M.cross.assign(S * S, desc_hot());
    M.same_c.assign(S * S, desc_cold());
    M.cross_c.assign(S * S, desc_cold());
    M.tight_keys.clear(); M.tight_lo.clear(); M.tight_hi.clear(); M.tight_row.clear();
    M.tight_dela.clear(); M.tight_delb.clear();
    M.pres_names.clear();
    std::map<std::string, int> pres_idx;
    auto pres = [&](const std::string& n) {
        auto it = pres_idx.find(n);
        if (it != pres_idx.end()) return it->second;
        int id = (int)M.pres_names.size();
        M.pres_names.push_back(n);
        pres_idx[n] = id;
        return id;
    };
    auto is_internal = [&](int t, int o) { return o > 0 && o < (int)M.tpl[t].obs.size() - 1; };
    auto has_close = [&](int t) { return M.tpl[t].obs.size() > 1; };
    auto int_name = [&](int t, int o) {
        const Template& T = M.tpl[t];
        return "int." + T.name + "." + std::to_string(o) + (T.el[T.obs[o]].kind == EK_POLY ? "p" : "a");
    };
    auto accumulate = [&](int t, int e_lo, int e_hi, double& fixed, double& var, double& hr, int& nins, bool& only_sp, bool& has_poly) {
        for (int x = e_lo; x < e_hi; ++x) {
            const TElem& e = M.tpl[t].el[x];
            if (e.kind == EK_INSERT) ++nins;
            else { fixed += e.nominal; var += e.var; hr += e.hr; }
            if (e.kind != EK_SPACER) only_sp = false;
            if (e.kind == EK_POLY) has_poly = true;
        }
    };
    auto finish = [&](desc_hot& d, desc_cold& c, double fixed, double var, double hr, int nins, const std::string& key,
                      bool poly_end, const std::string& row, int dela = 0, int delb = 0) {
        if (nins > 2) { d.valid = 0; return; }
        d.valid = 1;
        d.kind = (uint8_t)nins;
        d.fixed = (int32_t)std::lround(fixed);
        c.sd = std::sqrt(std::max(var, 1.0));
        c.key = key;
        d.slack = (int16_t)std::lround(hr + (poly_end ? 14 : 4));
        if (nins == 0) {
            int R = (int)std::ceil(4 * c.sd + 6);
            R = std::min(std::max(R, 12), 80);
            d.tid = (int16_t)M.tight_keys.size();
            d.dela = (int16_t)dela;
            d.delb = (int16_t)delb;
            M.tight_keys.push_back("tight." + key);
            // deleted primers: the gap may be up to dela + delb - RESCUE_GMIN bp shorter than `fixed`
            M.tight_lo.push_back(std::min(-R, -(dela + delb) + RESCUE_GMIN));
            M.tight_hi.push_back(key.back() == 'x' ? R + JUNC_SPACER_MAX : R);
            M.tight_row.push_back(row);
            M.tight_dela.push_back(dela);
            M.tight_delb.push_back(delb);
        }
    };
    // length of the barcode-adjacent primer skipped after element e at the template end (del_tail) / before e at
    // the template start (del_head) when only spacers lie between; else 0
    auto del_tail = [&](int t, int e) -> int {
        const Template& T = M.tpl[t];
        const int E = (int)T.el.size();
        if (e >= E - 2 || T.el[E - 1].kind != EK_ANCHOR) return 0;
        for (int x = e + 1; x < E - 1; ++x) if (T.el[x].kind != EK_SPACER) return 0;
        return M.anch[T.el[E - 1].anc].m;
    };
    auto del_head = [&](int t, int e) -> int {
        const Template& T = M.tpl[t];
        if (e < 2 || T.el[0].kind != EK_ANCHOR) return 0;
        for (int x = 1; x < e; ++x) if (T.el[x].kind != EK_SPACER) return 0;
        return M.anch[T.el[0].anc].m;
    };
    for (int sp = 0; sp < S; ++sp) {
        int t1 = M.st_tpl[sp], o1 = M.st_obs[sp];
        const Template& T1 = M.tpl[t1];
        for (int st = 0; st < S; ++st) {
            int t2 = M.st_tpl[st], o2 = M.st_obs[st];
            const Template& T2 = M.tpl[t2];
            const bool poly_end = M.st_type[sp] >= M.A || M.st_type[st] >= M.A;
            if (t1 == t2 && o2 > o1) {
                desc_hot& d = M.same[sp * S + st];
                desc_cold& c = M.same_c[sp * S + st];
                for (int x = o1 + 1; x < o2; ++x) c.pres.push_back({pres(int_name(t1, x)), false});
                if (is_internal(t2, o2)) c.pres.push_back({pres(int_name(t1, o2)), true});
                double fixed = 0, var = 0, hr = 0;
                int nins = 0;
                bool only_sp = true, hp = false;
                accumulate(t1, T1.obs[o1] + 1, T1.obs[o2], fixed, var, hr, nins, only_sp, hp);
                std::string row;
                if (only_sp && nins == 0 && T1.obs[o2] == T1.obs[o1] + 2) {
                    const TElem& spc = T1.el[T1.obs[o1] + 1];
                    if (spc.spec >= 0) row = M.spec.elements[spc.spec].id;
                }
                d.spacer_only = (uint8_t)(only_sp && nins == 0);
                finish(d, c, fixed, var, hr, nins, tkey(M, sp) + "-" + tkey(M, st) + "s", poly_end, row);
            }
            {
                desc_hot& d = M.cross[sp * S + st];
                desc_cold& c = M.cross_c[sp * S + st];
                const bool pal = palindromic(M, t1, t2);
                int nobs1 = (int)T1.obs.size();
                for (int x = o1 + 1; x < nobs1 - 1; ++x)
                    if (is_internal(t1, x)) c.pres.push_back({pres(int_name(t1, x)), false});
                const std::string jn = T1.name + ">" + T2.name;
                if (has_close(t1)) c.pres.push_back({pres((pal ? "pal.close." : "close.") + jn), o1 == nobs1 - 1});
                c.pres.push_back({pres((pal ? "pal.open." : "open.") + jn), o2 == 0});
                for (int x = 1; x < o2; ++x)
                    if (is_internal(t2, x)) c.pres.push_back({pres(int_name(t2, x)), false});
                if (is_internal(t2, o2)) c.pres.push_back({pres(int_name(t2, o2)), true});
                c.jfrom = t1; c.jto = t2;
                double fixed = 0, var = 4.0, hr = 0;
                int nins = 0;
                bool only_sp = false, hp = false;
                accumulate(t1, T1.obs[o1] + 1, (int)T1.el.size(), fixed, var, hr, nins, only_sp, hp);
                accumulate(t2, 0, T2.obs[o2], fixed, var, hr, nins, only_sp, hp);
                d.pal = (uint8_t)(pal && o2 == 0 && M.st_type[st] < M.A);
                // anchor slots on both sides (poly runs never certify): barcode-adjacent primers may be deleted
                const bool anchors = M.st_type[sp] < M.A && M.st_type[st] < M.A;
                const int da = anchors ? del_tail(t1, T1.obs[o1]) : 0, db = anchors ? del_head(t2, T2.obs[o2]) : 0;
                finish(d, c, fixed, var, hr, nins, tkey(M, sp) + "-" + tkey(M, st) + "x", poly_end, "", da, db);
            }
        }
    }
}

}  // namespace detail

// Builds the topology (and priors) from a layout; applies `p` when given; finalizes.
inline Model Model::build(const layout_spec& L, const Params* p, const model_options& o) {
    using namespace detail;
    validate_spec(L);
    Model M;
    M.spec = L;
    M.opt = o;
    std::vector<int> fidx, ridx;
    for (int i = 0; i < (int)L.elements.size(); ++i) (L.elements[i].direction == 'R' ? ridx : fidx).push_back(i);
    auto by_order = [&](int a, int b) { return L.elements[a].order != L.elements[b].order ? L.elements[a].order < L.elements[b].order : a < b; };
    std::stable_sort(fidx.begin(), fidx.end(), by_order);
    std::stable_sort(ridx.begin(), ridx.end(), by_order);
    Template F = chain_from_spec(M, L, fidx, 'F', "F");
    if (F.el.empty()) throw layout_spec_error("concat_hmm: layout '" + L.name + "' has no forward elements");
    Template R;
    bool mirror = true;
    if (!ridx.empty()) {
        R = chain_from_spec(M, L, ridx, 'R', "R");
        if (R.el.size() != F.el.size()) mirror = false;
        else
            for (size_t i = 0; i < F.el.size() && mirror; ++i) {
                const TElem& a = F.el[i];
                const TElem& b = R.el[F.el.size() - 1 - i];
                if (a.kind != b.kind) mirror = false;
                else if (a.kind == EK_ANCHOR && M.anch[b.anc].seq != revcomp(M.anch[a.anc].seq)) mirror = false;
            }
    } else {
        R.strand = 'R'; R.name = "R"; R.main = true;
        for (int i = (int)F.el.size() - 1; i >= 0; --i) R.el.push_back(rc_elem(M, F.el[i]));
    }
    M.tpl.push_back(F);
    M.tpl.push_back(R);
    // check capacity before any per-anchor array is indexed (auxiliary templates add no new types)
    if ((int)M.anch.size() > MAXANCH || (int)M.poly_base.size() > MAXPOLY)
        throw layout_spec_error("concat_hmm: layout has " + std::to_string(M.anch.size()) + " static elements (max " + std::to_string(MAXANCH) +
                              ") / " + std::to_string(M.poly_base.size()) + " poly types (max " + std::to_string(MAXPOLY) + "), both strands");
    if (M.anch.empty()) throw layout_spec_error("concat_hmm: layout '" + L.name + "' produced no anchor pattern");
    // weak outer anchors: an outermost anchor joined by spacers only to an inner anchor (main templates)
    std::vector<uint8_t> weak(M.anch.size(), 0);
    bool any_weak = false;
    for (int t = 0; t < 2; ++t) {
        Template& T = M.tpl[t];
        std::vector<int> obs;
        for (int i = 0; i < (int)T.el.size(); ++i)
            if (T.el[i].kind == EK_ANCHOR || T.el[i].kind == EK_POLY) obs.push_back(i);
        if (obs.size() < 2) continue;
        auto only_spacers = [&](int lo, int hi) {
            for (int x = lo + 1; x < hi; ++x) if (T.el[x].kind != EK_SPACER) return false;
            return true;
        };
        // the outer anchor must be the terminal element of the construct (a PCR primer: nothing beyond it)
        int e0 = obs[0], e1 = obs[1];
        if (e0 == 0 && T.el[e0].kind == EK_ANCHOR && T.el[e1].kind == EK_ANCHOR && only_spacers(e0, e1) && e1 > e0 + 1) {
            weak[T.el[e0].anc] = 1; any_weak = true;
        }
        int f1 = obs.back(), f0 = obs[obs.size() - 2];
        if (f1 == (int)T.el.size() - 1 && T.el[f1].kind == EK_ANCHOR && T.el[f0].kind == EK_ANCHOR && only_spacers(f0, f1) && f1 > f0 + 1) {
            weak[T.el[f1].anc] = 1; any_weak = true;
        }
    }
    const bool gen_art = mirror && (o.artifact_templates == 1 || (o.artifact_templates < 0 && any_weak));
    int nins = 0, r = -1;
    for (int i = 0; i < (int)F.el.size(); ++i)
        if (F.el[i].kind == EK_INSERT) { ++nins; r = i; }
    auto rc_of = [&](int fi) -> TElem {
        TElem t = M.tpl[1].el[F.el.size() - 1 - fi];
        return t;
    };
    if (mirror && nins == 1) {
        auto has_barcode = [&](int lo, int hi) {
            for (int i = lo; i < hi; ++i) if (F.el[i].kind == EK_SPACER) return true;
            return false;
        };
        auto anchors_only = [&](int lo, int hi) {
            for (int i = lo; i < hi; ++i) if (F.el[i].kind != EK_ANCHOR) return false;
            return hi > lo;
        };
        const int E = (int)F.el.size();
        if (!gen_art) {
            bool b3 = anchors_only(r + 1, E), b5 = anchors_only(0, r);
            if (b3 != b5) {  // barcode-less template-switch artifact T = bare side mirrored around the insert
                Template T;
                T.strand = 'T'; T.name = "T"; T.art = true;
                if (b3) {
                    for (int i = E - 1; i > r; --i) T.el.push_back(rc_of(i));
                    for (int i = r; i < E; ++i) T.el.push_back(F.el[i]);
                } else {
                    for (int i = 0; i <= r; ++i) T.el.push_back(F.el[i]);
                    for (int i = r - 1; i >= 0; --i) T.el.push_back(rc_of(i));
                }
                M.tpl.push_back(T);
            }
        } else {
            if (r > 0 && F.el[0].kind == EK_ANCHOR) {
                const bool bc = has_barcode(0, r), bare = anchors_only(0, r);
                Template A5;  // [5' chain] insert [rc(outer 5' anchor)]
                A5.art = true; A5.name = "A5"; A5.strand = bc ? 'F' : 'T';
                for (int i = 0; i <= r; ++i) A5.el.push_back(F.el[i]);
                A5.el.push_back(rc_of(0));
                M.tpl.push_back(A5);
                if (!bare) {
                    Template B;  // [outer 5' anchor] insert [rc(5' chain)]
                    B.art = true; B.name = "A5r"; B.strand = bc ? 'R' : 'T';
                    B.el.push_back(F.el[0]);
                    B.el.push_back(F.el[r]);
                    for (int i = r - 1; i >= 0; --i) B.el.push_back(rc_of(i));
                    M.tpl.push_back(B);
                }
                const bool want_d = o.d_template == 1 || (o.d_template < 0 && any_weak && weak[F.el[0].anc]);
                if (bc && want_d) {
                    Template D;  // [5' chain] insert [rc(5' chain)]: one molecule, barcode units at both ends
                    D.art = true; D.dtpl = true; D.name = "D"; D.strand = 'D';
                    for (int i = 0; i <= r; ++i) D.el.push_back(F.el[i]);
                    for (int i = r - 1; i >= 0; --i) D.el.push_back(rc_of(i));
                    M.tpl.push_back(D);
                }
            }
            if (r + 1 < E && F.el[E - 1].kind == EK_ANCHOR) {
                const bool bc = has_barcode(r + 1, E), bare = anchors_only(r + 1, E);
                Template A3;  // [rc(outer 3' anchor)] insert [3' chain]
                A3.art = true; A3.name = "A3"; A3.strand = bc ? 'R' : 'T';
                A3.el.push_back(rc_of(E - 1));
                for (int i = r; i < E; ++i) A3.el.push_back(F.el[i]);
                M.tpl.push_back(A3);
                if (!bare) {
                    Template B;  // [rc(3' chain)] insert [outer 3' anchor]
                    B.art = true; B.name = "A3r"; B.strand = bc ? 'F' : 'T';
                    for (int i = E - 1; i > r; --i) B.el.push_back(rc_of(i));
                    B.el.push_back(F.el[r]);
                    B.el.push_back(F.el[E - 1]);
                    M.tpl.push_back(B);
                }
            }
        }
    }
    M.NT = (int)M.tpl.size();
    for (auto& T : M.tpl) {
        for (auto& e : T.el) bind_spec(M, e);
        finish_template(M, T);
    }
    M.A = (int)M.anch.size();
    M.Q = (int)M.poly_base.size();
    if (M.A > MAXANCH || M.Q > MAXPOLY)
        throw layout_spec_error("concat_hmm: layout has " + std::to_string(M.A) + " static elements (max " + std::to_string(MAXANCH) +
                              ") / " + std::to_string(M.Q) + " poly types (max " + std::to_string(MAXPOLY) + ")");
    if (M.A == 0) throw layout_spec_error("concat_hmm: layout '" + L.name + "' produced no anchor pattern");
    M.NTYPE = M.A + M.Q;
    for (int a = 0; a < M.A; ++a) M.anch[a].weak = weak[a] != 0;
    for (auto& T : M.tpl) {
        if (T.obs.empty()) continue;
        if (T.el[T.obs.back()].kind == EK_ANCHOR && T.obs.size() > 1) M.anch[T.el[T.obs.back()].anc].closing = true;
        if (T.el[T.obs.front()].kind == EK_ANCHOR) M.anch[T.el[T.obs.front()].anc].opening = true;
    }
    for (int q = 0; q < M.Q; ++q) {
        char b = M.poly_base[q];
        int bi = b == 'A' ? 0 : b == 'C' ? 1 : b == 'G' ? 2 : 3;
        M.polyq_of_base[bi] = q;
        if (b == 'A') M.stage0_poly[0] = true;
        else if (b == 'T') M.stage0_poly[1] = true;
        else M.other_poly = true;
    }
    M.states_by_type.assign(M.NTYPE, {});
    M.open_state.assign(M.NT, -1);
    M.close_state.assign(M.NT, -1);
    for (int t = 0; t < M.NT; ++t) {
        const Template& T = M.tpl[t];
        if (T.main) M.main_tpl.push_back(t);
        for (int ob = 0; ob < (int)T.obs.size(); ++ob) {
            const TElem& e = T.el[T.obs[ob]];
            int ty = e.kind == EK_ANCHOR ? e.anc : M.A + e.anc;
            int st = (int)M.st_tpl.size();
            M.st_tpl.push_back(t); M.st_obs.push_back(ob); M.st_type.push_back(ty); M.st_elem.push_back(T.obs[ob]);
            uint16_t f = 0;
            if (e.kind == EK_ANCHOR) {
                f |= SF_ANCHOR;
                const Anchor& x = M.anch[e.anc];
                if (x.weak) f |= SF_WEAK;
                if (x.shrt) f |= SF_SHORT;
                if (!x.weak && !x.shrt) f |= SF_CERT;
            }
            if (ob == 0) { f |= SF_OPEN; M.open_state[t] = st; }
            if (ob == (int)T.obs.size() - 1 && T.obs.size() > 1) { f |= SF_CLOSE; M.close_state[t] = st; }
            if (T.main) f |= SF_MAIN;
            if (T.art) f |= SF_ART;
            if (T.dtpl) f |= SF_DTPL;
            M.st_flags.push_back(f);
            M.states_by_type[ty].push_back(st);
        }
    }
    M.S = (int)M.st_tpl.size();
    if (M.S > MAXSTATE) throw layout_spec_error("concat_hmm: too many observable slots (" + std::to_string(M.S) + ")");
    M.head_fix.assign(M.S, -1);
    M.tail_fix.assign(M.S, -1);
    M.head_edge.resize(M.S);
    M.tail_edge.resize(M.S);
    for (int st = 0; st < M.S; ++st) {
        const Template& T = M.tpl[M.st_tpl[st]];
        int ei = M.st_elem[st];
        if (T.pins[ei] == 0) M.head_fix[st] = (int)std::lround(T.pnom[ei]);
        if (T.pins[T.el.size()] - T.pins[ei + 1] == 0) M.tail_fix[st] = (int)std::lround(T.pnom[T.el.size()] - T.pnom[ei + 1]);
        auto edge_range = [&](int begin, int end) {
            edge_span r;
            for (int i = begin; i < end; ++i) {
                const TElem& e = T.el[i];
                if (e.kind == EK_INSERT) { r.bounded = false; break; }
                r.lo += e.lmin;
                r.hi += e.lmax;
            }
            return r;
        };
        M.head_edge[st] = edge_range(0, ei);
        M.tail_edge[st] = edge_range(ei + 1, (int)T.el.size());
    }
    M.bc_head.assign(M.S, -1);
    M.bc_tail.assign(M.S, -1);
    M.head_rules.clear();
    M.tail_rules.clear();
    for (int st = 0; st < M.S; ++st) {
        const Template& T = M.tpl[M.st_tpl[st]];
        const int ei = M.st_elem[st], E = (int)T.el.size();
        const bool poly = T.el[ei].kind == EK_POLY;
        auto block = [&](int lo, int hi) {  // spacer elements [lo, hi): summed length, -1 if anything else lies there
            if (hi <= lo) return -1;
            double s = 0;
            for (int x = lo; x < hi; ++x) {
                if (T.el[x].kind != EK_SPACER) return -1;
                s += poly ? T.el[x].lmax : T.el[x].nominal;
            }
            return (int)std::lround(s);
        };
        if (ei >= 2 && T.el[0].kind == EK_ANCHOR) M.bc_head[st] = block(1, ei);
        if (ei <= E - 3 && T.el[E - 1].kind == EK_ANCHOR) M.bc_tail[st] = block(ei + 1, E - 1);
        if (poly) continue;
        auto add_rule = [](std::vector<bc_rule>& v, bc_rule nr) {
            for (const auto& x : v)
                if (x.prim == nr.prim && x.inner == nr.inner && x.block == nr.block) return;
            v.push_back(nr);
        };
        if (M.bc_head[st] >= 0) add_rule(M.head_rules, bc_rule{T.el[0].anc, T.el[ei].anc, M.bc_head[st]});
        if (M.bc_tail[st] >= 0) add_rule(M.tail_rules, bc_rule{T.el[E - 1].anc, T.el[ei].anc, M.bc_tail[st]});
    }
    M.bc_unit.assign(M.S, 0);
    for (int st = 0; st < M.S; ++st) {
        if (M.bc_head[st] < 0 && M.bc_tail[st] < 0) continue;
        M.bc_unit[st] = 1;
        const int t = M.st_tpl[st];
        const Template& T = M.tpl[t];
        for (int ob = 0; ob < (int)T.obs.size(); ++ob) {
            const int e = T.obs[ob];
            if ((M.bc_head[st] >= 0 && e == 0) || (M.bc_tail[st] >= 0 && e == (int)T.el.size() - 1)) M.bc_unit[M.open_state[t] + ob] = 1;
        }
    }
    M.inner_role.assign(M.anch.size() + M.poly_base.size(), 0);
    if (!M.head_rules.empty() && !M.tail_rules.empty()) {
        for (const auto& tr : M.tail_rules) { M.inner_role[tr.inner] |= 1; M.inner_role[tr.prim] |= 4; }
        for (const auto& hr : M.head_rules) { M.inner_role[hr.inner] |= 2; M.inner_role[hr.prim] |= 8; }
    }
    M.tpl_fold.assign(M.NT, 0);
    for (int t = 0; t < M.NT; ++t) {
        const Template& T = M.tpl[t];
        if (!T.art || T.obs.size() < 2) continue;
        const TElem& a0 = T.el[T.obs.front()];
        const TElem& a1 = T.el[T.obs.back()];
        bool spacer = false;
        for (auto& e : T.el) spacer |= e.kind == EK_SPACER;
        // D and single-primer ends with a barcode unit: terminal anchors reverse-complementary across a barcode spacer
        if (spacer && a0.kind == EK_ANCHOR && a1.kind == EK_ANCHOR && M.anch[a1.anc].seq == revcomp(M.anch[a0.anc].seq)) M.tpl_fold[t] = 1;
    }
    M.tpl_mirror.assign(M.NT, {});
    for (int t = 0; t < M.NT; ++t) {
        if (!M.tpl_fold[t]) continue;
        const Template& T = M.tpl[t];
        const int E = (int)T.el.size(), no = (int)T.obs.size();
        M.tpl_mirror[t].assign(no, -1);
        for (int ob = 0; ob < no; ++ob) {
            const TElem& a = T.el[T.obs[ob]];
            const int em = E - 1 - T.obs[ob];
            if (a.kind != EK_ANCHOR || em == T.obs[ob]) continue;
            for (int o2 = 0; o2 < no; ++o2)
                if (T.obs[o2] == em && T.el[em].kind == EK_ANCHOR && M.anch[T.el[em].anc].seq == revcomp(M.anch[a.anc].seq)) M.tpl_mirror[t][ob] = o2;
        }
    }
    M.nbr.assign(M.NTYPE, {});
    for (int t = 0; t < M.NT; ++t) {
        const Template& T = M.tpl[t];
        for (size_t ob = 0; ob + 1 < T.obs.size(); ++ob) {
            int e1 = T.obs[ob], e2 = T.obs[ob + 1];
            if (T.pins[e2] - T.pins[e1 + 1] != 0) continue;
            double lo = 0, hi = 0;
            for (int x = e1 + 1; x < e2; ++x) {
                const TElem& e = T.el[x];
                if (e.kind == EK_POLY) { lo += 0; hi += 60; }
                else { lo += e.lmin; hi += e.lmax; }
            }
            auto ty = [&](const TElem& e) { return e.kind == EK_ANCHOR ? e.anc : M.A + e.anc; };
            int a = ty(T.el[e1]), b = ty(T.el[e2]);
            auto add = [&](int x, Nbr n) {
                for (auto& y : M.nbr[x])
                    if (y.other == n.other && y.dir == n.dir) { y.lo = std::min(y.lo, n.lo); y.hi = std::max(y.hi, n.hi); return; }
                M.nbr[x].push_back(n);
            };
            add(a, Nbr{b, +1, (int)lo, (int)hi});
            add(b, Nbr{a, -1, (int)lo, (int)hi});
        }
    }
    {
        const Template& T = M.tpl[0];
        double s5 = 0, s3 = 0;
        bool after = false;
        for (auto& e : T.el) {
            if (e.kind == EK_INSERT) { after = true; continue; }
            if (e.kind == EK_POLY) continue;
            (after ? s3 : s5) += e.kind == EK_SPACER ? e.lmin : e.nominal;
        }
        M.min_terminal = std::max(1, o.min_terminal_len);
        M.min_interior = o.min_interior_len >= 0 ? o.min_interior_len : std::max(M.min_terminal, (int)std::max(s5, s3));
    }
    build_patterns(M);
    build_engine(M);
    build_seed_index(M);
    build_descriptors(M);
    make_priors(M);
    uint64_t h = 1469598103934665603ULL;
    for (auto& T : M.tpl) {
        h = fnv(h, T.name + T.strand);
        for (auto& e : T.el) {
            std::string d = std::to_string((int)e.kind) + ":" + std::to_string(e.lmin) + "-" + std::to_string(e.lmax);
            if (e.kind == EK_ANCHOR) d += M.anch[e.anc].seq;
            if (e.kind == EK_POLY) d += M.poly_base[e.anc];
            h = fnv(h, d);
        }
    }
    M.topo_hash = h & 0x1FFFFFFFFFFFFFULL;  // exactly representable as a double
    M.par = M.prior;
    if (p && !p->empty()) M.apply_params(*p);
    M.finalize();
    return M;
}

// Parameter values -> log tables. Missing arrays fall back to priors.
inline void Model::finalize() {
    using namespace detail;
    auto get = [&](const std::string& k) -> const std::vector<double>& {
        auto it = par.find(k);
        if (it != par.end()) {
            auto ip = prior.find(k);
            if (ip == prior.end() || ip->second.size() == it->second.size()) return it->second;
        }
        return prior.at(k);
    };
    auto sl = [](double p) { return (float)std::log(std::max(p, 1e-300)); };
    emit.assign(S * NFEAT, NEG);
    emit_pal.assign(S * NFEAT, NEG);
    emit_bonus.assign(S * NFEAT, 0.f);
    for (int st = 0; st < S; ++st) {
        int ty = st_type[st];
        const auto& q = get("q." + type_name(ty));
        const auto& nl = get("null." + type_name(ty));
        const std::vector<double>* qp = ty < A ? &get("qpal." + type_name(ty)) : &q;
        const bool open = (st_flags[st] & SF_OPEN) != 0, close = (st_flags[st] & SF_CLOSE) != 0;
        for (int f = 0; f < NFEAT; ++f) {
            if (q[f] <= 0 || nl[f] <= 0) continue;
            if (ty < A) {
                // relaxed / read-start suffix features only in opening slots; read-end prefix features only in closing slots
                const Anchor& x = anch[ty];
                if (f > x.k && f < FEAT_PREFIX && !open) continue;
                if (f == FEAT_TRUNC5 && !open) continue;
                if ((f >= FEAT_PREFIX && f < FEAT_PREFIX + NPARTBIN) || f == FEAT_TRUNC3S) {
                    if (!close) continue;
                }
            }
            emit[st * NFEAT + f] = sl(q[f]) - sl(nl[f]);
            if ((*qp)[f] > 0) emit_pal[st * NFEAT + f] = sl((*qp)[f]) - sl(nl[f]);
            // short anchors: LLR capped at SHORT_LLR_CAP; the excess is kept in emit_bonus
            if (ty < A && anch[ty].shrt) {
                emit_bonus[st * NFEAT + f] = std::max(0.f, emit[st * NFEAT + f] - SHORT_LLR_CAP);
                emit[st * NFEAT + f] = std::min(emit[st * NFEAT + f], SHORT_LLR_CAP);
                emit_pal[st * NFEAT + f] = std::min(emit_pal[st * NFEAT + f], SHORT_LLR_CAP);
            }
        }
    }
    std::vector<double> pres(pres_names.size());
    for (size_t i = 0; i < pres.size(); ++i) pres[i] = get("pres." + pres_names[i])[0];
    auto fill = [&](std::vector<desc_hot>& H, std::vector<desc_cold>& C) {
        for (int x = 0; x < S * S; ++x) {
            if (!H[x].valid) continue;
            double b = 0;
            for (auto& pt : C[x].pres) {
                double p = std::min(std::max(pres[pt.first], 1e-4), 1 - 1e-4);
                b += std::log(pt.second ? p : 1 - p);
            }
            if (C[x].jfrom >= 0) b += std::log(std::max(get("junc." + tpl[C[x].jfrom].name)[C[x].jto], 1e-6));
            H[x].base = (float)b;
        }
    };
    fill(same, same_c);
    fill(cross, cross_c);
    tight.assign(tight_keys.size(), {});
    for (size_t t = 0; t < tight_keys.size(); ++t) {
        const auto& pm = get(tight_keys[t]);
        tight[t].resize(pm.size());
        for (size_t i = 0; i < pm.size(); ++i) tight[t][i] = sl(std::max(pm[i], 1e-7));
    }
    for (int w = 0; w < 2; ++w) {
        const auto& pb = get(w == 0 ? "ins1" : "ins2");
        std::vector<float>& T = ins[w];
        T.assign(INS_HI - INS_LO, NEG);
        for (int x = 0; x < INS_HI; ++x) {
            int b = ins_bin(x);
            T[x - INS_LO] = sl(std::max(pb[b], 1e-9) / ins_bin_count(b));
        }
        float v0 = T[0 - INS_LO];
        for (int x = INS_LO; x < 0; ++x) T[x - INS_LO] = v0 + 0.15f * x;
    }
    {
        const auto& pb = get("ins1");
        double c = 0;
        ins1_p99 = INS_HI;
        for (int b = 0; b < NINSBIN; ++b) {
            c += pb[b];
            if (c >= 0.99) { ins1_p99 = (int)std::lround(ins_bin_lo(b + 1)); break; }
        }
    }
    lbegin.assign(S, NEG);
    lend.assign(S, NEG);
    const auto& b = get("begin");
    const auto& e = get("end");
    for (int st = 0; st < S; ++st) { lbegin[st] = sl(std::max(b[st], 1e-6)); lend[st] = sl(std::max(e[st], 1e-6)); }
    double p0 = std::min(std::max(get("p0")[0], 1e-5), 0.5);
    lp0 = (float)std::log(p0);
    lp1 = (float)std::log(1 - p0);
    const auto& g = get("gate");
    zone5 = (int)std::lround(std::min(std::max(g[0], 60.0), 600.0));
    zone3 = (int)std::lround(std::min(std::max(g[1], 40.0), 600.0));
    fast_p = (float)std::min(std::max(g[2], 0.5), 0.99999);
    guard_failed = get("guard")[0] > 0.5;
    // cDNA strand table: [mean sense score, 4096 log-odds]; absent or malformed -> no strand-flip cuts
    strand_lo.clear();
    strand_mu = 0.f;
    {
        auto it = par.find("strand");
        if (it != par.end() && it->second.size() == (size_t)STRAND_N + 1) {
            bool ok = true;
            for (double x : it->second) ok = ok && std::isfinite(x);
            if (ok && it->second[0] > 0) {
                strand_mu = (float)it->second[0];
                strand_lo.resize((size_t)STRAND_N);
                for (int i = 0; i < STRAND_N; ++i) strand_lo[i] = (float)it->second[i + 1];
            }
        }
    }
    // D strand table: counts -> log-odds (pseudocount 1); unused when the mean sense score is <= 0.01 nats
    strand_d_lo.clear();
    strand_d_mu = 0.f;
    {
        auto it = par.find("strand_d");
        if (it != par.end() && it->second.size() == (size_t)STRAND_D_N) {
            const std::vector<double>& c = it->second;
            std::vector<float> lo((size_t)STRAND_D_N);
            double num = 0, tot = 0;
            for (int kk = 0; kk < STRAND_D_N; ++kk) {
                int rk = 0;
                for (int x = 0, v = kk; x < STRAND_D_K; ++x, v >>= 2) rk = (rk << 2) | (3 - (v & 3));
                const double l = std::log((c[(size_t)kk] + 1.0) / (c[(size_t)rk] + 1.0));
                lo[(size_t)kk] = (float)l;
                num += c[(size_t)kk] * l;
                tot += c[(size_t)kk];
            }
            if (tot > 0 && num / tot > 0.01) { strand_d_lo.swap(lo); strand_d_mu = (float)(num / tot); }
        }
    }
    finalized = true;
}

inline std::string Model::describe() const {
    using namespace detail;
    std::ostringstream os;
    os << "layout " << spec.name << (spec.mode.empty() ? "" : " mode=" + spec.mode) << "\n";
    os << "anchors (" << A << "):\n";
    for (int i = 0; i < A; ++i)
        os << "  a" << i << " " << anch[i].row << " " << anch[i].seq << " m=" << anch[i].m << " k=" << anch[i].k
           << " strong<=" << anch[i].strong_ed << (anch[i].opening ? " opening" : "") << (anch[i].closing ? " closing" : "")
           << (anch[i].weak ? " WEAK" : "") << (anch[i].shrt ? " SHORT" : "") << (anch[i].nchunk != 1 ? " chunks=" + std::to_string(anch[i].nchunk) : "")
           << "\n";
    os << "poly types:";
    for (char c : poly_base) os << " " << c;
    os << "\nMyers words: " << words.size() << " packed + " << lwords.size() << " long; relaxed opener words " << rwords.size()
       << "; Wov=" << Wov << "; seeds K=" << K << "\n";
    for (int t = 0; t < NT; ++t) {
        os << "template " << tpl[t].name << " (" << tpl[t].strand << (tpl[t].art ? ", aux" : "") << "): ";
        for (auto& e : tpl[t].el) {
            if (e.kind == EK_ANCHOR) os << "[" << e.label << ":a" << e.anc << "] ";
            else if (e.kind == EK_POLY) os << "(" << poly_base[e.anc] << ">=" << POLY_MIN << ") ";
            else if (e.kind == EK_SPACER) os << "<" << e.label << ":" << e.lmin << ".." << e.lmax << "> ";
            else os << "{" << e.label << "} ";
        }
        os << "\n";
    }
    for (const auto& r : head_rules)
        os << "barcode block (head): [" << anch[r.prim].row << "] <" << r.block << "> [" << anch[r.inner].row << "]\n";
    for (const auto& r : tail_rules)
        os << "barcode block (tail): [" << anch[r.inner].row << "] <" << r.block << "> [" << anch[r.prim].row << "]\n";
    os << "states (" << S << "):";
    for (int s = 0; s < S; ++s) os << " " << s << "=" << tkey(*this, s) << "/" << type_label(st_type[s]);
    int nsame = 0, ncross = 0;
    for (auto& d : same) nsame += d.valid;
    for (auto& d : cross) ncross += d.valid;
    os << "\ntransitions: " << nsame << " within-construct + " << ncross << " junction; tight tables " << tight_keys.size()
       << "; presence params " << pres_names.size() << "\nmin construct length interior " << min_interior << " terminal "
       << min_terminal << "; gate zones 5' " << zone5 << " 3' " << zone3 << "; fast-path p " << fast_p
       << (guard_failed ? "; CALIBRATION GUARD FAILED (layout priors)" : "") << "\n";
    return os.str();
}

namespace detail {

// s,e: union of detected 12-bp windows (sampled every 4th base). Extend to the exact union of all
// qualifying windows (<= 1 base different) touching the region, then trim to the base.
inline void flush_poly(const char* seq, int n, int s, int e, char want_uc, uint8_t q, std::vector<poly_reg>& out) {
    const char want = (char)(want_uc | 0x20);
    auto mis = [&](int j) { return (int)((seq[j] | 0x20) != want); };
    if (s > 0) {
        int lo = std::max(0, s - 11);
        int ws = s - 1;
        if (ws + 12 <= n) {
            int mm = 0;
            for (int j = ws; j < ws + 12; ++j) mm += mis(j);
            int best = s;
            for (;;) {
                if (mm <= 1) best = ws;
                if (ws == lo) break;
                mm += mis(ws - 1) - mis(ws + 11);
                --ws;
            }
            s = best;
        }
    }
    if (e < n) {
        int hi = std::min(n, e + 11);
        int we = e + 1;
        if (we - 12 >= 0) {
            int mm = 0;
            for (int j = we - 12; j < we; ++j) mm += mis(j);
            int best = e;
            for (;;) {
                if (mm <= 1) best = we;
                if (we == hi) break;
                mm += mis(we) - mis(we - 12);
                ++we;
            }
            e = best;
        }
    }
    while (s < e && (seq[s] | 0x20) != want) ++s;
    while (e > s && (seq[e - 1] | 0x20) != want) --e;
    if (e - s >= POLY_MIN) out.push_back(poly_reg{s, e, q});
}

// Byte-equality bitmask of 8 bytes against a broadcast byte (exact SWAR zero-byte test + gather).
inline uint64_t byte_eq_bits(uint64_t x, uint64_t pat) {
    const uint64_t y = x ^ pat;
    uint64_t t = ((y & 0x7F7F7F7F7F7F7F7FULL) + 0x7F7F7F7F7F7F7F7FULL) | y;
    t = ~t & 0x8080808080808080ULL;
    return ((t >> 7) * 0x0102040810204080ULL) >> 56;
}
// Exact homopolymer scan for poly bases other than A/T (rare layouts): window of 12 with <= 1 mismatch.
inline void poly_scan_other(const Model& M, const char* s, int n, std::vector<poly_reg>& out) {
    for (int q = 0; q < M.Q; ++q) {
        char B = M.poly_base[q];
        if (B == 'A' || B == 'T') continue;
        int mm = 0, run_s = -1, last_e = -1000;
        std::vector<std::pair<int, int>> wins;
        for (int i = 0; i < n; ++i) {
            mm += ((s[i] & 0xDF) != B);
            if (i >= 12) mm -= ((s[i - 12] & 0xDF) != B);
            if (i >= 11 && mm <= 1) {
                int ws = i - 11;
                if (run_s >= 0 && ws <= last_e + 5) last_e = i + 1;
                else { if (run_s >= 0) flush_poly(s, n, run_s, last_e, B, (uint8_t)q, out); run_s = ws; last_e = i + 1; }
            }
        }
        if (run_s >= 0) flush_poly(s, n, run_s, last_e, B, (uint8_t)q, out);
    }
    std::sort(out.begin(), out.end(), [](const poly_reg& x, const poly_reg& y) { return x.s < y.s; });
}

inline void stage0(const Model& M, const char* seq, int n, Scratch& S) {
    S.seeds.clear();
    S.poly.clear();
    const uint64_t* bits = M.seed_bits.data();
    const uint64_t* bloom = M.seed_bloom.data();
    const uint64_t mask = M.kmask;
    const int K = M.K;
    const int W = 12, MG = 5;
    const uint64_t F12 = 0x555555ULL;
    const bool want_a = M.stage0_poly[0], want_t = M.stage0_poly[1];
    const uint8_t qa = (uint8_t)std::max(0, M.polyq_of_base[0]), qt = (uint8_t)std::max(0, M.polyq_of_base[3]);
    uint64_t hist = 0x6666666666666666ULL;
    int a_s = -1, a_e = -1000, t_s = -1, t_e = -1000;
    const unsigned char* sp = (const unsigned char*)seq;
    auto seed_hit = [&](int i, uint64_t km) {
        if (i < K - 1) return;
        uint16_t lab = seed_label(M, km);
        if (lab != 0xFFFF) S.seeds.push_back(((uint64_t)(uint32_t)i << 16) | lab);
    };
    auto poly_hit = [&](int i, bool pa, bool pt) {
        if (i < W - 1) return;
        const int ws = i - W + 1;
        if (pa && want_a) {
            if (ws <= a_e + MG) a_e = i + 1;
            else { if (a_s >= 0) flush_poly(seq, n, a_s, a_e, 'A', qa, S.poly); a_s = ws; a_e = i + 1; }
        }
        if (pt && want_t) {
            if (ws <= t_e + MG) t_e = i + 1;
            else { if (t_s >= 0) flush_poly(seq, n, t_s, t_e, 'T', qt, S.poly); t_s = ws; t_e = i + 1; }
        }
    };
    if (S.hits.size() < (size_t)(n / 8 + 2)) S.hits.resize((size_t)(n / 8 + 2) * 2);
    uint32_t* hits = S.hits.data();
    const uint64_t polymask = (want_a || want_t) ? ~0ULL : 0ULL;
    int nh = 0, i = 0;
    for (; i + 8 <= n; i += 8) {
        uint64_t w;
        memcpy(&w, sp + i, 8);
        w = ((w >> 1) ^ (w >> 2)) & 0x0303030303030303ULL;
        w = __builtin_bswap64(w);
        w = (w | (w >> 6)) & 0x000F000F000F000FULL;
        w = (w | (w >> 12)) & 0x000000FF000000FFULL;
        w = (w | (w >> 24)) & 0xFFFFULL;
        hist = (hist << 16) | w;
        // L1-resident 2^19-bit pre-filter; the exact 4^K bitset is consulted when a recorded step is resolved
        const uint64_t BM = (1u << 19) - 1;
        const uint64_t b0 = (hist >> 14) & BM, b1 = (hist >> 12) & BM, b2 = (hist >> 10) & BM, b3 = (hist >> 8) & BM;
        const uint64_t b4 = (hist >> 6) & BM, b5 = (hist >> 4) & BM, b6 = (hist >> 2) & BM, b7 = hist & BM;
        const uint64_t any = ((bloom[b0 >> 6] >> (b0 & 63)) | (bloom[b1 >> 6] >> (b1 & 63)) | (bloom[b2 >> 6] >> (b2 & 63)) |
                              (bloom[b3 >> 6] >> (b3 & 63)) | (bloom[b4 >> 6] >> (b4 & 63)) | (bloom[b5 >> 6] >> (b5 & 63)) |
                              (bloom[b6 >> 6] >> (b6 & 63)) | (bloom[b7 >> 6] >> (b7 & 63))) & 1;
        // poly windows tested at bases i+3 and i+7: <= 1 non-A (non-T) base in 12 <=> z & (z-1) == 0.
        // Two-lane SIMD on NEON / SSE4.1; identical scalar fallback.
#if defined(CONCAT_HMM_NEON)
        uint64_t pfl;
        {
            const uint64x2_t xv = vcombine_u64(vcreate_u64((hist >> 8) & 0xFFFFFFULL), vcreate_u64(hist & 0xFFFFFFULL));
            const uint64x2_t f12 = vdupq_n_u64(F12), one = vdupq_n_u64(1);
            const uint64x2_t xs = vshrq_n_u64(xv, 1);
            const uint64x2_t z_a = vandq_u64(vorrq_u64(xv, xs), f12);
            const uint64x2_t z_t = vbicq_u64(f12, vandq_u64(xv, xs));
            const uint64x2_t e_a = vceqzq_u64(vandq_u64(z_a, vsubq_u64(z_a, one)));
            const uint64x2_t e_t = vceqzq_u64(vandq_u64(z_t, vsubq_u64(z_t, one)));
            const uint64x2_t m_a = vcombine_u64(vcreate_u64(2), vcreate_u64(8)), m_t = vcombine_u64(vcreate_u64(4), vcreate_u64(16));
            pfl = vaddvq_u64(vorrq_u64(vandq_u64(e_a, m_a), vandq_u64(e_t, m_t)));
        }
        const uint64_t fl = any | (pfl & polymask);
#elif defined(CONCAT_HMM_SSE41)
        uint64_t pfl;
        {
            const __m128i xv = _mm_set_epi64x((long long)(hist & 0xFFFFFFULL), (long long)((hist >> 8) & 0xFFFFFFULL));
            const __m128i f12 = _mm_set1_epi64x((long long)F12), one = _mm_set1_epi64x(1);
            const __m128i xs = _mm_srli_epi64(xv, 1);
            const __m128i z_a = _mm_and_si128(_mm_or_si128(xv, xs), f12);
            const __m128i z_t = _mm_andnot_si128(_mm_and_si128(xv, xs), f12);
            const __m128i z0 = _mm_setzero_si128();
            const __m128i e_a = _mm_cmpeq_epi64(_mm_and_si128(z_a, _mm_sub_epi64(z_a, one)), z0);
            const __m128i e_t = _mm_cmpeq_epi64(_mm_and_si128(z_t, _mm_sub_epi64(z_t, one)), z0);
            const __m128i v = _mm_or_si128(_mm_and_si128(e_a, _mm_set_epi64x(8, 2)), _mm_and_si128(e_t, _mm_set_epi64x(16, 4)));
            pfl = (uint64_t)_mm_cvtsi128_si64(v) | (uint64_t)_mm_extract_epi64(v, 1);
        }
        const uint64_t fl = any | (pfl & polymask);
#else
        const uint64_t x1 = (hist >> 8) & 0xFFFFFFULL, x2 = hist & 0xFFFFFFULL;
        const uint64_t z_a1 = (x1 | (x1 >> 1)) & F12, z_t1 = ~(x1 & (x1 >> 1)) & F12;
        const uint64_t z_a2 = (x2 | (x2 >> 1)) & F12, z_t2 = ~(x2 & (x2 >> 1)) & F12;
        const uint64_t pa1 = (z_a1 & (z_a1 - 1)) == 0, pt1 = (z_t1 & (z_t1 - 1)) == 0;
        const uint64_t pa2 = (z_a2 & (z_a2 - 1)) == 0, pt2 = (z_t2 & (z_t2 - 1)) == 0;
        const uint64_t fl = any | (((pa1 << 1) | (pt1 << 2) | (pa2 << 3) | (pt2 << 4)) & polymask);
#endif
        hits[nh] = ((uint32_t)(i >> 3) << 5) | (uint32_t)fl;  // step index (i is a multiple of 8): reads up to 2^30 bp
        nh += (fl != 0);
    }
    for (int h = 0; h < nh; ++h) {
        const int i0 = (int)((hits[h] >> 5) << 3);
        const uint32_t fl = hits[h] & 31;
        if (fl & 1) {
            uint64_t km = 0;
            const int st = i0 - K + 1;
            for (int j = std::max(0, st); j < i0; ++j) km = (km << 2) | code2(sp[j]);
            for (int j = 0; j < 8; ++j) {
                km = ((km << 2) | code2(sp[i0 + j])) & mask;
                if ((bits[km >> 6] >> (km & 63)) & 1) seed_hit(i0 + j, km);
            }
        }
        if (fl & 6) poly_hit(i0 + 3, (fl >> 1) & 1, (fl >> 2) & 1);
        if (fl & 24) poly_hit(i0 + 7, (fl >> 3) & 1, (fl >> 4) & 1);
    }
    for (; i < n; ++i) {
        hist = (hist << 2) | code2(sp[i]);
        const uint64_t km = hist & mask;
        if ((bits[km >> 6] >> (km & 63)) & 1) seed_hit(i, km);
        if (polymask) {
            const uint64_t x = hist & 0xFFFFFFULL;
            const uint64_t z_a = (x | (x >> 1)) & F12, z_t = ~(x & (x >> 1)) & F12;
            const bool pa = (z_a & (z_a - 1)) == 0, pt = (z_t & (z_t - 1)) == 0;
            if (pa | pt) poly_hit(i, pa, pt);
        }
    }
    if (a_s >= 0) flush_poly(seq, n, a_s, a_e, 'A', qa, S.poly);
    if (t_s >= 0) flush_poly(seq, n, t_s, t_e, 'T', qt, S.poly);
    if (M.other_poly) poly_scan_other(M, seq, n, S.poly);
    if (S.poly.size() > 1)
        std::sort(S.poly.begin(), S.poly.end(), [](const poly_reg& x, const poly_reg& y) { return x.s < y.s; });
    // seed clusters per anchor (implied element start within 4 + n of the cluster)
    S.sc.clear();
    int cur[MAXANCH];
    for (int a = 0; a < MAXANCH; ++a) cur[a] = -1;
    for (uint64_t v : S.seeds) {
        const int pos = (int)(v >> 16), a = (int)((v >> 8) & 0xFF), off = (int)(v & 0x7F);
        const uint8_t ex = (v & 0x80) ? 0 : 1;
        const int st = pos - K + 1 - off;
        const int en = st + M.anch[a].m;
        const int ci = cur[a];
        if (ci >= 0 && std::abs(st - S.sc[ci].s) <= 4 + (int)S.sc[ci].n) {
            seed_cluster& c = S.sc[ci];
            c.e = std::max(c.e, en);
            c.s = std::min(c.s, st);
            if (c.n < 255) c.n++;
            if (c.nx < 255) c.nx += ex;
            c.olo = (int16_t)std::min<int>(c.olo, off);
            c.ohi = (int16_t)std::max<int>(c.ohi, off);
        } else {
            cur[a] = (int)S.sc.size();
            S.sc.push_back(seed_cluster{st, en, (uint8_t)a, 1, ex, (int16_t)off, (int16_t)off});
        }
    }
}

inline void clu_add(Clu& c, int j, int sc, int pat, std::vector<raw_ev>& out) {
    if (c.last >= 0 && j - c.last > KGAP) {
        out.push_back(raw_ev{c.best_end, c.first, c.last, (uint8_t)pat, (uint8_t)c.best});
        clu_reset(c);
    }
    if (c.first < 0) c.first = j;
    c.last = j;
    if (sc < c.best) { c.best = sc; c.best_end = j; }
}
inline void clu_flush(Clu& c, int pat, std::vector<raw_ev>& out) {
    if (c.first >= 0) out.push_back(raw_ev{c.best_end, c.first, c.last, (uint8_t)pat, (uint8_t)c.best});
    clu_reset(c);
}
template <int NW>
struct SWords { uint64_t S[NW]; };
template <int NW>
#if defined(__GNUC__)
__attribute__((noinline))
#endif
void report_hits(const PWord* Wd, SWords<NW> sw, int j, Clu* cl, std::vector<raw_ev>& out) {
    for (int w = 0; w < NW; ++w) {
        const PWord& W = Wd[w];
        int a = (int)((sw.S[w] >> W.t_a) & W.FA);
        if (a <= W.k_a) clu_add(cl[W.pa], j, a, W.pa, out);
        if (W.pb >= 0) {
            int b = (int)(sw.S[w] >> W.t_b);
            if (b <= W.k_b) clu_add(cl[W.pb], j, b, W.pb, out);
        }
    }
}
template <int NW>
inline SWords<NW> swords(const PState* st) {
    SWords<NW> r;
    for (int w = 0; w < NW; ++w) r.S[w] = st[w].S;
    return r;
}

#define CHMM_STEP(st, w, c)                                          \
    {                                                                \
        const int wc = SH ? 0 : w;                                   \
        const uint64_t Eq = Wd[w].peq[c];                            \
        const uint64_t Xv = Eq | st.Mv;                              \
        const uint64_t sum = (Eq & st.Pv) + (st.Pv & CG[wc]);        \
        const uint64_t Xh = (sum ^ st.Pv) | Eq;                      \
        uint64_t Ph = st.Mv | ~(Xh | st.Pv);                         \
        uint64_t Mh = st.Pv & Xh;                                    \
        st.S += (Ph & HB[wc]);                                       \
        st.S -= (Mh & HB[wc]);                                       \
        Ph = (Ph << 1) & CS[wc];                                     \
        Mh = (Mh << 1) & CS[wc];                                     \
        st.Pv = Mh | ~(Xv | Ph);                                     \
        st.Mv = Ph & Xv;                                             \
    }

// Packed Myers over NW words (2 patterns each): raw events to raw_a / raw_b, final column in fin[].
template <int NW, bool SH>
inline void myers_scan(const PWord* Wd, int Wov, int NP, const char* s, int n, std::vector<raw_ev>& raw_a, std::vector<raw_ev>& raw_b, PState* fin) {
    uint64_t CG[NW], HB[NW], CS[NW], TT[NW], HH[NW];
    PState A[NW], B[NW];
    for (int w = 0; w < NW; ++w) {
        CG[w] = Wd[w].CG; HB[w] = Wd[w].HB; CS[w] = Wd[w].CS; TT[w] = Wd[w].T; HH[w] = Wd[w].H;
        A[w] = PState{~0ULL, 0, Wd[w].S0};
        B[w] = A[w];
    }
    Clu cl_a[MAXPAT], cl_b[MAXPAT];
    for (int p = 0; p < NP; ++p) { clu_reset(cl_a[p]); clu_reset(cl_b[p]); }
    raw_a.clear(); raw_b.clear();
    const unsigned char* p1 = (const unsigned char*)s;
    uint64_t NH[NW];
    for (int w = 0; w < NW; ++w) NH[w] = SH ? 0 : ~HH[w];
    const uint64_t HALL = SH ? HH[0] : ~0ULL;
    if (n < 4 * Wov) {
        for (int j = 0; j < n; ++j) {
            const unsigned c = p1[j];
            uint64_t z = ~0ULL;
            for (int w = 0; w < NW; ++w) {
                CHMM_STEP(A[w], w, c)
                z &= ((A[w].S | HH[SH ? 0 : w]) - TT[SH ? 0 : w]) | NH[w];
            }
            if (__builtin_expect((z & HALL) != HALL, 0)) report_hits<NW>(Wd, swords<NW>(A), j, cl_a, raw_a);
        }
        for (int w = 0; w < NW; ++w) fin[w] = A[w];
    } else {
        const int h = n / 2, s2 = h - Wov;
        const unsigned char* p2 = p1 + s2;
        for (int j = 0; j < h; ++j) {
            const unsigned c1 = p1[j], c2 = p2[j];
            uint64_t za = ~0ULL, zb = ~0ULL;
            for (int w = 0; w < NW; ++w) {
                CHMM_STEP(A[w], w, c1)
                CHMM_STEP(B[w], w, c2)
                za &= ((A[w].S | HH[SH ? 0 : w]) - TT[SH ? 0 : w]) | NH[w];
                zb &= ((B[w].S | HH[SH ? 0 : w]) - TT[SH ? 0 : w]) | NH[w];
            }
            if (__builtin_expect(((za & zb) & HALL) != HALL, 0)) {
                if ((za & HALL) != HALL) report_hits<NW>(Wd, swords<NW>(A), j, cl_a, raw_a);
                if ((zb & HALL) != HALL && j >= Wov) report_hits<NW>(Wd, swords<NW>(B), s2 + j, cl_b, raw_b);
            }
        }
        const int end2 = n - s2;
        for (int j = h; j < end2; ++j) {
            const unsigned c2 = p2[j];
            uint64_t zb = ~0ULL;
            for (int w = 0; w < NW; ++w) {
                CHMM_STEP(B[w], w, c2)
                zb &= ((B[w].S | HH[SH ? 0 : w]) - TT[SH ? 0 : w]) | NH[w];
            }
            if (__builtin_expect((zb & HALL) != HALL, 0)) report_hits<NW>(Wd, swords<NW>(B), s2 + j, cl_b, raw_b);
        }
        for (int w = 0; w < NW; ++w) fin[w] = B[w];
    }
    for (int p = 0; p < NP; ++p) { clu_flush(cl_a[p], p, raw_a); clu_flush(cl_b[p], p, raw_b); }
    if (!raw_b.empty()) {  // merge clusters straddling the split point
        for (int p = 0; p < NP; ++p) {
            raw_ev* a = nullptr;
            for (int i = (int)raw_a.size() - 1; i >= 0; --i)
                if (raw_a[i].pat == p) { a = &raw_a[i]; break; }
            if (!a) continue;
            for (auto& b : raw_b) {
                if (b.pat != p || b.first < 0) continue;
                if (b.first - a->last <= KGAP) {
                    if (b.score < a->score) { a->score = b.score; a->end = b.end; }
                    a->last = b.last;
                    b.first = -1;
                }
                break;
            }
        }
    }
}
#undef CHMM_STEP

inline void myers_all(const std::vector<PWord>& words, bool shared, int Wov, int NP, const char* s, int n, Scratch& W, PState* fin,
                      std::vector<raw_ev>& out) {
    const int NW = (int)words.size();
    for (int wb = 0; wb < NW; wb += 4) {
        const int nw = std::min(NW - wb, 4);
        const PWord* Wd = words.data() + wb;
        if (shared) {
            switch (nw) {
                case 1: myers_scan<1, true>(Wd, Wov, NP, s, n, W.raw_a, W.raw_b, fin + wb); break;
                case 2: myers_scan<2, true>(Wd, Wov, NP, s, n, W.raw_a, W.raw_b, fin + wb); break;
                case 3: myers_scan<3, true>(Wd, Wov, NP, s, n, W.raw_a, W.raw_b, fin + wb); break;
                default: myers_scan<4, true>(Wd, Wov, NP, s, n, W.raw_a, W.raw_b, fin + wb); break;
            }
        } else {
            switch (nw) {
                case 1: myers_scan<1, false>(Wd, Wov, NP, s, n, W.raw_a, W.raw_b, fin + wb); break;
                case 2: myers_scan<2, false>(Wd, Wov, NP, s, n, W.raw_a, W.raw_b, fin + wb); break;
                case 3: myers_scan<3, false>(Wd, Wov, NP, s, n, W.raw_a, W.raw_b, fin + wb); break;
                default: myers_scan<4, false>(Wd, Wov, NP, s, n, W.raw_a, W.raw_b, fin + wb); break;
            }
        }
        for (auto& r : W.raw_a) if (r.first >= 0) out.push_back(r);
        for (auto& r : W.raw_b) if (r.first >= 0) out.push_back(r);
    }
}
// Scalar Myers for one pattern: appends clusters (<= L.k) and returns the final column.
inline void long_scan(const LWord& L, const char* s, int n, std::vector<raw_ev>& out, uint64_t& Pv_out, uint64_t& Mv_out) {
    uint64_t Pv = ~0ULL, Mv = 0;
    int sc = L.m;
    const int sh = L.m - 1;
    Clu c;
    clu_reset(c);
    for (int j = 0; j < n; ++j) {
        uint64_t Eq = L.peq[(unsigned char)s[j]];
        uint64_t Xv = Eq | Mv;
        uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
        uint64_t Ph = Mv | ~(Xh | Pv);
        uint64_t Mh = Pv & Xh;
        sc += (int)((Ph >> sh) & 1) - (int)((Mh >> sh) & 1);
        Ph <<= 1; Mh <<= 1;
        Pv = Mh | ~(Xv | Ph);
        Mv = Ph & Xv;
        if (__builtin_expect(sc <= L.k, 0)) clu_add(c, j, sc, L.pat, out);
    }
    clu_flush(c, L.pat, out);
    Pv_out = Pv; Mv_out = Mv;
}
// Minimum semi-global edit distance of one pattern in s[0,n) and its (first) end position.
inline int best_hit(const LWord& L, const char* s, int n, int& best_end, uint64_t* Pv_out = nullptr, uint64_t* Mv_out = nullptr) {
    uint64_t Pv = ~0ULL, Mv = 0;
    int sc = L.m, best = 1 << 20;
    best_end = -1;
    const int sh = L.m - 1;
    for (int j = 0; j < n; ++j) {
        uint64_t Eq = L.peq[(unsigned char)s[j]];
        uint64_t Xv = Eq | Mv;
        uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
        uint64_t Ph = Mv | ~(Xh | Pv);
        uint64_t Mh = Pv & Xh;
        sc += (int)((Ph >> sh) & 1) - (int)((Mh >> sh) & 1);
        Ph <<= 1; Mh <<= 1;
        Pv = Mh | ~(Xv | Ph);
        Mv = Ph & Xv;
        if (sc < best) { best = sc; best_end = j; }
    }
    if (Pv_out) *Pv_out = Pv;
    if (Mv_out) *Mv_out = Mv;
    return best;
}
// Read-end prefix from the final Myers column (bits [off, off+m)): returns the retained length i maximising
// v - 3D (v = non-N bases of the prefix, v >= PARTIAL_MIN, D <= v/8), or -1; d_out gets D.
inline int prefix_from_column(const Pattern& p, uint64_t Pv, uint64_t Mv, int off, int& d_out) {
    int best_i = -1, best_d = 0;
    double best_sc = -1e9;
    for (int i = PARTIAL_MIN; i < p.m; ++i) {
        const int v = p.inf[i];
        if (v < PARTIAL_MIN) continue;
        uint64_t mask = ((1ULL << i) - 1) << off;
        int d = __builtin_popcountll(Pv & mask) - __builtin_popcountll(Mv & mask);
        if (d > v / 8) continue;
        double sc = v - 3.0 * d;
        if (sc > best_sc) { best_sc = sc; best_i = i; best_d = d; }
    }
    d_out = best_d;
    return best_i;
}

}  // namespace detail

namespace detail {

inline uint8_t anchor_class(const Anchor& x, int ed) {
    if (ed <= x.strong_ed) return CL_S4;
    if (ed <= x.k) return CL_S6;
    return CL_S8;
}
inline void push_anchor_event(std::vector<Event>& ev, int n, int start, int end, int a, int feat, uint8_t cls, uint8_t fl,
                              int rlo, int rhi, int ed) {
    if ((int)ev.size() >= MAX_EVENTS) return;
    Event e;
    e.start = std::max(0, start);
    e.end = std::min(n, std::max(end, e.start + 1));
    e.type = (uint8_t)a;
    e.feat = (uint8_t)feat;
    e.cls = cls;
    e.fl = fl;
    e.rlo = (int16_t)rlo;
    e.rhi = (int16_t)rhi;
    e.ed = (int16_t)ed;
    e.len = (int16_t)std::max(0, rhi - rlo);
    ev.push_back(e);
}
// Raw pattern-level Myers clusters (positions shifted by `base`) -> anchor events.
inline void raw_to_events(const Model& M, const std::vector<raw_ev>& raw, int base, int n, bool relaxed, std::vector<Event>& ev,
                          std::vector<raw_ev>& chunked) {
    for (const raw_ev& r : raw) {
        const Pattern& p = M.pats[r.pat];
        const Anchor& x = M.anch[p.anchor];
        if (x.nchunk != 1 || p.off != 0 || p.m != x.m) {
            raw_ev c = r;
            c.end += base;
            chunked.push_back(c);
            continue;
        }
        const int end = r.end + base + 1;
        const int ed = r.score;
        if (relaxed && ed <= x.k) continue;  // already reported by the normal scan
        push_anchor_event(ev, n, end - x.m, end, p.anchor, ed, anchor_class(x, ed), relaxed ? EVF_RELAXED : 0, 0, x.m, ed);
    }
}
// Chunked anchors: chunk hits with a consistent start form one event; a missing chunk costs its budget + 1.
inline void merge_chunks(const Model& M, std::vector<raw_ev>& ch, int n, std::vector<Event>& ev) {
    if (ch.empty()) return;
    auto istart = [&](const raw_ev& r) { const Pattern& p = M.pats[r.pat]; return r.end + 1 - p.m - p.off; };
    std::sort(ch.begin(), ch.end(), [&](const raw_ev& a, const raw_ev& b) {
        int aa = M.pats[a.pat].anchor, bb = M.pats[b.pat].anchor;
        return aa != bb ? aa < bb : istart(a) < istart(b);
    });
    size_t i = 0;
    while (i < ch.size()) {
        const int a = M.pats[ch[i].pat].anchor;
        const Anchor& x = M.anch[a];
        const int s0 = istart(ch[i]);
        int best[MAXPAT];
        for (int c = 0; c < MAXPAT; ++c) best[c] = 1 << 20;
        size_t j = i;
        while (j < ch.size() && M.pats[ch[j].pat].anchor == a && istart(ch[j]) - s0 <= 8) {
            best[ch[j].pat] = std::min(best[ch[j].pat], (int)ch[j].score);
            ++j;
        }
        int ed = 0, rlo = 1 << 20, rhi = -1, found = 0;
        for (int c = x.chunk0; c < x.chunk0 + x.nchunk; ++c) {
            if (best[c] < (1 << 20)) {
                ed += best[c];
                rlo = std::min(rlo, M.pats[c].off);
                rhi = std::max(rhi, M.pats[c].off + M.pats[c].m);
                ++found;
            } else ed += M.pats[c].k + 1;
        }
        if (found > 0) {
            if (ed <= x.k) push_anchor_event(ev, n, s0, s0 + x.m, a, ed, anchor_class(x, ed), 0, 0, x.m, ed);
            else push_anchor_event(ev, n, s0, s0 + x.m, a, FEAT_SEEDPART + seedpart_bin(rhi - rlo), CL_SP, EVF_SEEDPART, rlo, rhi, -1);
        }
        i = j;
    }
}
// Retained offsets [rlo, rhi) and mismatches of an adapter whose seed-implied span crosses a read end.
inline void seed_truncation(const Model& M, const seed_cluster& c, const char* seq, int n, int& rlo, int& rhi, int& ed) {
    const Anchor& x = M.anch[c.a];
    if (c.e > n) { rlo = 0; rhi = std::min(x.m, std::max(0, n - c.s)); }
    else { rlo = std::min(x.m, std::max(0, -c.s)); rhi = x.m; }
    ed = 0;
    for (int o = rlo; o < rhi; ++o) {
        const int p = c.s + o;
        if (p < 0 || p >= n) continue;
        const char a = x.seq[o], b = (char)(seq[p] & 0xDF);
        ed += a != 'N' && a != b;
    }
}
inline void sort_events(std::vector<Event>& ev) {
    std::sort(ev.begin(), ev.end(), [](const Event& a, const Event& b) {
        return a.start != b.start ? a.start < b.start : a.end != b.end ? a.end < b.end : a.type < b.type;
    });
}
// Same-anchor events whose ends are <= 12 bp apart are one occurrence: keep the better one.
inline void dedupe_events(const Model& M, std::vector<Event>& ev) {
    size_t w = 0;
    bool replaced = false;
    for (size_t i = 0; i < ev.size(); ++i) {
        Event& e = ev[i];
        bool dup = false;
        if (e.type < M.A) {
            for (size_t k = w; k-- > 0 && w - k <= 8;) {
                Event& f = ev[k];
                if (f.type != e.type) continue;
                if (std::abs(f.end - e.end) <= 12 || (e.start < f.end && f.start < e.end && (e.fl & (EVF_SEEDPART | EVF_TRUNC3 | EVF_TRUNC5)))) {
                    auto rank = [](const Event& x) { return x.cls == CL_SP ? 100 : x.ed; };
                    if (rank(e) < rank(f)) { f = e; replaced = true; }
                    dup = true;
                    break;
                }
            }
        }
        if (!dup) ev[w++] = e;
    }
    ev.resize(w);
    if (replaced) sort_events(ev);
}

// Drops primer events laid >= 8 bp over the barcode block of a facing partner anchor (anchor grade, ED <= 1).
inline void drop_phantom_primers(const Model& M, std::vector<Event>& ev) {
    if (M.head_rules.empty() || M.tail_rules.empty()) return;
    const size_t n = ev.size();
    size_t i0 = 0;
    while (i0 < n && !(ev[i0].type < M.A && (M.inner_role[ev[i0].type] & 12) && ev[i0].cls <= CL_S8)) ++i0;
    if (i0 == n) return;
    size_t w = i0;
    for (size_t i = i0; i < n; ++i) {
        const Event e = ev[i];
        bool phantom = false;
        if (e.type < M.A && (M.inner_role[e.type] & 12) && e.cls <= CL_S8) {
            for (const bc_rule& tr : M.tail_rules) {  // closing primer over the barcode block of an opening unit after it
                if (tr.prim != e.type) continue;
                for (const bc_rule& hr : M.head_rules)
                    for (size_t j = i + 1; j < n && ev[j].start < e.end + hr.block && !phantom; ++j) {
                        const Event& t = ev[j];
                        if (t.type != hr.inner || t.cls != CL_S4 || t.ed < 0 || t.ed > 1 || t.start < e.end - 4) continue;
                        phantom = std::min(e.end, t.start) - std::max(e.start, t.start - hr.block) >= 8;
                    }
            }
            for (const bc_rule& hr : M.head_rules) {  // opening primer over the barcode block of a closing unit before it
                if (phantom || hr.prim != e.type) continue;
                for (const bc_rule& tr : M.tail_rules)
                    for (size_t j = w; j-- > 0 && !phantom;) {
                        const Event& t = ev[j];
                        if (t.start < e.start - tr.block - 80) break;
                        if (t.type != tr.inner || t.cls != CL_S4 || t.ed < 0 || t.ed > 1 || t.end > e.start + 4) continue;
                        phantom = std::min(e.end, t.end + tr.block) - std::max(e.start, t.end) >= 8;
                    }
            }
        }
        if (!phantom) ev[w++] = e;
    }
    ev.resize(w);
}

// Stage 2: events for the decoder; whole = scan the whole read, else only the stage-0 evidence windows.
inline void extract_events(const Model& M, const char* seq, int n, Scratch& W, bool whole) {
    std::vector<Event>& ev = W.ev;
    ev.clear();
    if (n <= 0) return;
    std::vector<raw_ev>& raw = W.raw_x;
    std::vector<raw_ev>& chunked = W.raw_ch;
    raw.clear(); chunked.clear();
    PState fin[MAXPAT];
    PState fin_end[MAXPAT];
    uint64_t l_pv[MAXPAT], l_mv[MAXPAT];
    bool have_end = false;
    auto scan = [&](int lo, int hi) {
        raw.clear();
        myers_all(M.words, M.words_shared, M.Wov, (int)M.pats.size(), seq + lo, hi - lo, W, fin, raw);
        raw_to_events(M, raw, lo, n, false, ev, chunked);
        if (hi == n) {
            for (size_t w = 0; w < M.words.size(); ++w) fin_end[w] = fin[w];
            have_end = true;
        }
        for (size_t l = 0; l < M.lwords.size(); ++l) {
            raw.clear();
            uint64_t Pv, Mv;
            long_scan(M.lwords[l], seq + lo, hi - lo, raw, Pv, Mv);
            raw_to_events(M, raw, lo, n, false, ev, chunked);
            if (hi == n) { l_pv[l] = Pv; l_mv[l] = Mv; }
        }
    };
    if (whole) {
        scan(0, n);
    } else {
        auto& Wn = W.win;
        Wn.clear();
        auto add = [&](int lo, int hi) {
            lo = std::max(0, lo);
            hi = std::min(n, hi);
            if (hi > lo) Wn.push_back({lo, hi});
        };
        bool open5 = false, close3 = false;
        for (const seed_cluster& c : W.sc) {
            if (!c.valid()) continue;
            const Anchor& x = M.anch[c.a];
            const int pad = 8 + x.k;
            // exact full-length seed cluster: emit it directly and scan only towards its junction partner
            const bool perfect = c.s >= 0 && c.e <= n && x.nchunk == 1 && c.olo == 0 && c.ohi == x.m - M.K &&
                                 (int)c.nx >= x.m - M.K + 1 && informative(x.seq) == x.m && x.m <= PAT_MAXLEN;
            if (x.opening && c.s <= M.zone5) open5 = true;
            if (x.closing && c.e >= n - M.zone3) close3 = true;
            if (perfect) {
                push_anchor_event(ev, n, c.s, c.e, c.a, 0, CL_S4, 0, 0, x.m, 0);
                if (x.opening) add(c.s - 40 - pad, c.s + 4);
                if (x.closing) add(c.e - 4, c.e + 40 + pad);
            } else {
                int lo = c.s - pad, hi = c.e + pad;
                if (x.opening) lo -= 40;
                if (x.closing) hi += 40;
                add(lo, hi);
            }
            for (const Nbr& nb : M.nbr[c.a]) {
                if (nb.other >= M.A) continue;
                const int mb = M.anch[nb.other].m, pb = 8 + M.anch[nb.other].k;
                if (nb.dir > 0) add(c.e + nb.lo - pb, c.e + nb.hi + mb + pb);
                else add(c.s - nb.hi - mb - pb, c.s - nb.lo + pb);
            }
        }
        for (const poly_reg& p : W.poly) {
            const int ty = M.A + p.q;
            for (const Nbr& nb : M.nbr[ty]) {
                if (nb.other >= M.A) continue;
                const Anchor& x = M.anch[nb.other];
                const int mb = x.m, pb = 10 + x.k;
                int lo, hi;
                if (nb.dir > 0) { lo = p.e + nb.lo - pb - 12; hi = p.e + nb.hi + mb + pb + 6; }
                else { lo = p.s - nb.hi - mb - pb - 6; hi = p.s - nb.lo + pb + 12; }
                if (x.opening) lo -= 30;
                if (x.closing) hi += 30;
                add(lo, hi);
            }
        }
        if (!open5) add(0, M.zone5 + 40);
        if (!close3) add(n - M.zone3 - 40, n);
        add(n - (M.max_close_m + 12), n);
        std::sort(Wn.begin(), Wn.end());
        int m = 0;
        const int SEP = M.Wov;  // a separator of >= m + k non-matching bytes decouples two windows
        for (size_t i = 0; i < Wn.size(); ++i) {
            if (m && Wn[i].first <= Wn[m - 1].second + SEP) Wn[m - 1].second = std::max(Wn[m - 1].second, Wn[i].second);
            else Wn[m++] = Wn[i];
        }
        Wn.resize(m);
        if (m == 1) scan(Wn[0].first, Wn[0].second);
        else if (m > 1) {
            std::string& C = W.cat;
            C.clear();
            W.cat_off.clear();
            for (int i = 0; i < m; ++i) {
                if (i) C.append((size_t)SEP, '\0');
                W.cat_off.push_back((int)C.size());
                C.append(seq + Wn[i].first, (size_t)(Wn[i].second - Wn[i].first));
            }
            const int nc = (int)C.size();
            auto map_pos = [&](int p) {
                int lo = 0, hi = m - 1;
                while (lo < hi) { const int mid = (lo + hi + 1) / 2; if (W.cat_off[mid] <= p) lo = mid; else hi = mid - 1; }
                const int len = Wn[lo].second - Wn[lo].first;
                const int q = std::min(p - W.cat_off[lo], len - 1);
                return Wn[lo].first + q;
            };
            raw.clear();
            myers_all(M.words, M.words_shared, M.Wov, (int)M.pats.size(), C.data(), nc, W, fin, raw);
            for (auto& r : raw) { r.end = map_pos(r.end); r.first = map_pos(r.first); r.last = map_pos(r.last); }
            raw_to_events(M, raw, 0, n, false, ev, chunked);
            if (Wn[m - 1].second == n) {
                for (size_t w = 0; w < M.words.size(); ++w) fin_end[w] = fin[w];
                have_end = true;
            }
            for (size_t l = 0; l < M.lwords.size(); ++l) {
                raw.clear();
                uint64_t Pv, Mv;
                long_scan(M.lwords[l], C.data(), nc, raw, Pv, Mv);
                for (auto& r : raw) { r.end = map_pos(r.end); r.first = map_pos(r.first); r.last = map_pos(r.last); }
                raw_to_events(M, raw, 0, n, false, ev, chunked);
                if (Wn[m - 1].second == n) { l_pv[l] = Pv; l_mv[l] = Mv; }
            }
        }
    }
    merge_chunks(M, chunked, n, ev);
    // read-end prefixes of closing anchors (final Myers column at the read end)
    if (have_end && n >= PARTIAL_MIN + 4) {
        auto try_prefix = [&](int pat, uint64_t Pv, uint64_t Mv, int off) {
            const Pattern& p = M.pats[pat];
            const Anchor& x = M.anch[p.anchor];
            if (!x.closing || p.off != 0) return;
            for (auto& e : ev)
                if (e.type == p.anchor && e.cls != CL_SP && e.end >= n - x.m / 2) return;
            int d = 0;
            int L = prefix_from_column(p, Pv, Mv, off, d);
            if (L < 0) return;
            push_anchor_event(ev, n, n - L - d, n - L - d + x.m, p.anchor, FEAT_PREFIX + prefix_bin(p.inf[L]), CL_SP, EVF_TRUNC3 | EVF_PREFIX, 0, L, d);
        };
        for (size_t w = 0; w < M.words.size(); ++w) {
            const PWord& Wd = M.words[w];
            try_prefix(Wd.pa, fin_end[w].Pv, fin_end[w].Mv, 0);
            if (Wd.pb >= 0) try_prefix(Wd.pb, fin_end[w].Pv, fin_end[w].Mv, Wd.m_a + 1);
        }
        for (size_t l = 0; l < M.lwords.size(); ++l) try_prefix(M.lwords[l].pat, l_pv[l], l_mv[l], 0);
    }
    // relaxed opener scan right after strong closers (junction spacing)
    if (M.opt.relaxed_openers && !M.rwords.empty()) {
        const size_t ne = ev.size();
        int last_hi = -1;
        for (size_t i = 0; i < ne; ++i) {
            const Event e = ev[i];
            if (e.type >= M.A || e.cls != CL_S4 || !M.anch[e.type].closing || M.anch[e.type].weak) continue;
            if (n - e.end < M.min_terminal) continue;  // no room for another construct
            int lo = std::max(e.end - 10, last_hi), hi = std::min(n, e.end + 48);
            if (hi - lo < 16) continue;
            bool have_open = false;
            for (size_t k2 = 0; k2 < ne; ++k2) {
                const Event& f = ev[k2];
                if (f.type < M.A && M.anch[f.type].opening && f.cls <= CL_S6 && f.start >= e.end - 12 && f.start <= e.end + 12) { have_open = true; break; }
            }
            if (have_open) continue;
            raw.clear();
            myers_all(M.rwords, M.rwords_shared, M.Wov, (int)M.pats.size(), seq + lo, hi - lo, W, fin, raw);
            raw_to_events(M, raw, lo, n, true, ev, chunked);
            last_hi = hi;
        }
    }
    // seed-only events: adapters running off a read end; internal partial chains
    for (const seed_cluster& c : W.sc) {
        if (!c.valid()) continue;
        const Anchor& x = M.anch[c.a];
        if (x.nchunk != 1) continue;
        bool covered = false;
        for (const Event& e : ev)
            if (e.type == c.a && e.start < c.e + 6 && c.s < e.end + 6) { covered = true; break; }
        if (covered) continue;
        if (c.e > n && x.closing && c.nx >= 1) {
            int rlo, rhi, ed;
            seed_truncation(M, c, seq, n, rlo, rhi, ed);
            push_anchor_event(ev, n, c.s, c.e, c.a, FEAT_TRUNC3S, CL_SP, EVF_TRUNC3, rlo, rhi, ed);
        } else if (c.s < 0 && x.opening && c.nx >= 1) {
            int rlo, rhi, ed;
            seed_truncation(M, c, seq, n, rlo, rhi, ed);
            push_anchor_event(ev, n, c.s, c.e, c.a, FEAT_TRUNC5, CL_SP, EVF_TRUNC5, rlo, rhi, ed);
        } else if (c.s >= 0 && c.e <= n && c.nx >= 1 && c.n >= 2) {
            const int rlo = c.olo, rhi = std::min(x.m, (int)c.ohi + M.K);
            push_anchor_event(ev, n, c.s, c.e, c.a, FEAT_SEEDPART + seedpart_bin(rhi - rlo), CL_SP, EVF_SEEDPART, rlo, rhi, -1);
        }
    }
    for (const poly_reg& p : W.poly) {
        if ((int)ev.size() >= MAX_EVENTS) break;
        Event e;
        e.start = p.s; e.end = p.e;
        e.type = (uint8_t)(M.A + p.q);
        e.feat = (uint8_t)poly_bin(p.e - p.s);
        e.cls = CL_POLY;
        e.fl = 0;
        e.rlo = e.rhi = -1;
        e.ed = -1;
        e.len = (int16_t)std::min(p.e - p.s, 30000);
        ev.push_back(e);
    }
    sort_events(ev);
    dedupe_events(M, ev);
    drop_phantom_primers(M, ev);
}

}  // namespace detail

namespace detail {

// Expected window [ws, we) of slot o from the nearest seen slot without an insert between, else from [cs, ce).
// es_out: expected element start, INT32_MIN if unplaced. cs_slack / ce_slack widen it inward.
inline void expected_window(const Model& M, int t, int o, const int* ostart, const int* oend, int cs, int ce, int pad, int& ws, int& we,
                            int& spread, int* es_out = nullptr, int cs_slack = 0, int ce_slack = 0) {
    if (es_out) *es_out = INT32_MIN;
    const Template& T = M.tpl[t];
    const int nobs = (int)T.obs.size();
    const int eo = T.obs[o];
    const int len = (int)std::lround(T.el[eo].nominal);
    const int E = (int)T.el.size();
    int es = INT32_MIN;
    double sp = 0;
    for (int ob = o - 1; ob >= 0; --ob) {
        const int e1 = T.obs[ob];
        if (T.pins[eo] - T.pins[e1 + 1] != 0) break;
        if (ostart[ob] >= 0) { es = oend[ob] + (int)std::lround(T.pnom[eo] - T.pnom[e1 + 1]); sp = T.phr[eo] - T.phr[e1 + 1]; break; }
    }
    if (es == INT32_MIN)
        for (int oa = o + 1; oa < nobs; ++oa) {
            const int e2 = T.obs[oa];
            if (T.pins[e2] - T.pins[eo + 1] != 0) break;
            if (ostart[oa] >= 0) { es = ostart[oa] - (int)std::lround(T.pnom[e2] - T.pnom[eo + 1]) - len; sp = T.phr[e2] - T.phr[eo + 1]; break; }
        }
    int lslack = 0, rslack = 0;
    if (es == INT32_MIN) {
        if (T.pins[eo] == 0) { es = cs + (int)std::lround(T.pnom[eo]); sp = T.phr[eo]; rslack = cs_slack; }
        else if (T.pins[E] - T.pins[eo + 1] == 0) { es = ce - (int)std::lround(T.pnom[E] - T.pnom[eo + 1]) - len; sp = T.phr[E] - T.phr[eo + 1]; lslack = ce_slack; }
        else { ws = cs; we = ce; spread = 1 << 20; return; }
    }
    spread = (int)std::lround(sp);
    ws = es - spread - pad - lslack;
    we = es + len + spread + pad + rslack;
    if (es_out) *es_out = es;
}

inline int check_element(const Model&, const TElem& e) { return e.spec; }

inline void add_check(Result& R, const Model& M, const TElem& el, char strand, int construct, Status st, int ws, int we, int rf, int rt,
                      int ed, float conf, int n) {
    check_region c;
    c.start = std::max(0, ws);
    c.end = std::min(n, we);
    if (c.end <= c.start) {
        if (st != Status::TRUNCATED_AT_READ_END) return;
        c.start = std::min(c.start, n);
        c.end = c.start;
    }
    c.element = check_element(M, el);
    // check_region::element indexes layout_spec::elements; bind_spec() binds every anchor / poly slot
    if (c.element < 0 || c.element >= (int)M.spec.elements.size()) {
        throw std::runtime_error("concat_hmm: invalid layout element in check region");
    }
    c.strand = strand;
    c.construct = construct;
    c.status = st;
    c.retained_from = rf;
    c.retained_to = rt;
    c.edits = ed;
    c.conf = conf;
    R.checks.push_back(c);
}

// Slot completion: finds a retained prefix / suffix (>= 10 bp, D <= i/8) of anchor a inside [ws, we).
// Returns its adapter offsets [rlo, rhi), read span [ps, pe) and edits.
inline bool partial_in_window(const Model& M, int a, const char* s, int n, int ws, int we, int& ps, int& pe, int& rlo, int& rhi, int& ed) {
    const Anchor& x = M.anch[a];
    if (x.nchunk != 1 || x.m < 16 || x.m > PAT_MAXLEN || informative(x.seq) != x.m) return false;
    ws = std::max(0, ws); we = std::min(n, we);
    if (we - ws < 10 || we - ws > 400) return false;
    // pigeonhole pre-filter (exact): no adapter 5-mer block in the window -> no prefix / suffix piece
    bool want[2] = {x.nblk == 0, x.nblk == 0};
    if (x.nblk > 0) {
        uint32_t h = 0;
        for (int j = ws; j < we && !(want[0] && want[1]); ++j) {
            h = ((h << 2) | code2((unsigned char)s[j])) & 1023u;
            if (j - ws < 4) continue;
            for (int b = 0; b < x.nblk; ++b) { want[0] |= h == x.pblk[b]; want[1] |= h == x.sblk[b]; }
        }
        if (!want[0] && !want[1]) return false;
    }
    double best = -1e9;
    bool found = false;
    for (int dir = 0; dir < 2; ++dir) {
        if (!want[dir]) continue;
        const LWord& L = dir == 0 ? M.single[a] : M.single_rev[a];
        uint64_t Pv = ~0ULL, Mv = 0;
        for (int t = 0; t < we - ws; ++t) {
            const int j = dir == 0 ? ws + t : we - 1 - t;
            const uint64_t Eq = L.peq[(unsigned char)s[j]];
            const uint64_t Xv = Eq | Mv;
            const uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
            uint64_t Ph = Mv | ~(Xh | Pv);
            uint64_t Mh = Pv & Xh;
            Ph <<= 1; Mh <<= 1;
            Pv = Mh | ~(Xv | Ph);
            Mv = Ph & Xv;
            if (t < 9) continue;
            // cheap column filter: a prefix of >= 10 bp passing D <= i/8 needs D(10) <= 2 or D(16) <= 2
            {
                const uint64_t m10 = (1ULL << 10) - 1, m16 = (1ULL << 16) - 1;
                const int d10 = __builtin_popcountll(Pv & m10) - __builtin_popcountll(Mv & m10);
                const int d16 = __builtin_popcountll(Pv & m16) - __builtin_popcountll(Mv & m16);
                if (d10 > 2 && d16 > 2) continue;
            }
            for (int i = x.m - 1; i >= 10; --i) {
                const uint64_t mask = (1ULL << i) - 1;
                const int d = __builtin_popcountll(Pv & mask) - __builtin_popcountll(Mv & mask);
                if (d > i / 8) continue;
                const double sc = i - 3.0 * d;
                if (sc > best) {
                    best = sc; found = true; ed = d;
                    if (dir == 0) { rlo = 0; rhi = i; pe = j + 1; ps = std::max(0, pe - i - d); }
                    else { rlo = x.m - i; rhi = x.m; ps = j; pe = std::min(n, ps + i + d); }
                }
                break;
            }
        }
    }
    return found;
}

// Check region of an expected slot with no event: PARTIAL / TRUNCATED_AT_READ_END from slot completion, else
// TRUNCATED_AT_READ_END when the expected extent runs past a read end, else MISSING_EXPECTED.
inline void missing_slot_check(Result& R, const Model& M, const TElem& el, char strand, int cidx, const char* seq, int n, int ws, int we,
                               int es, float conf_partial, float conf_geom) {
    if (ws >= n || we <= 0) return;
    const int pad = M.opt.seen_pad;
    if (el.kind == EK_ANCHOR) {
        int ps, pe, rlo, rhi, ed;
        if (partial_in_window(M, el.anc, seq, n, ws, we, ps, pe, rlo, rhi, ed)) {
            const int m = M.anch[el.anc].m;
            const bool at3 = pe >= n - 6 && rlo == 0, at5 = ps <= 6 && rhi == m;
            add_check(R, M, el, strand, cidx, (at3 || at5) ? Status::TRUNCATED_AT_READ_END : Status::PARTIAL, ps - rlo - pad,
                      pe + (m - rhi) + pad, rlo, rhi, ed, conf_partial, n);
            return;
        }
    }
    const bool anchor = el.kind == EK_ANCHOR;
    const int len = anchor ? M.anch[el.anc].m : (int)std::lround(el.nominal);
    Status st = Status::MISSING_EXPECTED;
    int rf = -1, rt = -1;
    if (es != INT32_MIN) {
        if (es + len > n) {
            st = Status::TRUNCATED_AT_READ_END;
            if (anchor) { rf = 0; rt = std::min(std::max(n - es, 0), len); }
        } else if (es < 0) {
            st = Status::TRUNCATED_AT_READ_END;
            if (anchor) { rf = std::min(std::max(-es, 0), len); rt = len; }
        }
    }
    add_check(R, M, el, strand, cidx, st, ws, we, rf, rt, -1, conf_geom, n);
}

// True when an outward-truncated fragment explains a degraded full-length hit better (i - 3D > m - 3 ed + 2);
// returns the fragment's PARTIAL window and offsets.
inline bool degraded_as_fragment(const Model& M, int a, const char* seq, int n, int s, int e, int ed, int& ws, int& we, int& rlo, int& rhi,
                                 int& fed) {
    const Anchor& x = M.anch[a];
    if (ed < 0 || (!x.closing && !x.opening) || (x.closing && x.opening)) return false;
    const int pad = M.opt.seen_pad;
    int ps, pe, lo, hi, d;
    if (!partial_in_window(M, a, seq, n, s - pad, e + pad, ps, pe, lo, hi, d)) return false;
    if (x.closing ? lo != 0 : hi != x.m) return false;  // the lost part must face outward
    if (hi - lo < 12 || (hi - lo) - 3 * d <= x.m - 3 * ed + 2) return false;
    ws = ps - lo - pad;
    we = pe + (x.m - hi) + pad;
    rlo = lo; rhi = hi; fed = d;
    return true;
}

// Stage 1: decides a clean single construct (k = 1) from seeds; false sends the read to the full path.
inline bool gate(const Model& M, const char* seq, int n, Scratch& S, Result& R) {
    int cnt[MAXANCH], first[MAXANCH];
    for (int a = 0; a < M.A; ++a) { cnt[a] = 0; first[a] = -1; }
    int nvalid = 0;
    for (size_t i = 0; i < S.sc.size(); ++i) {
        const seed_cluster& c = S.sc[i];
        if (!c.valid()) continue;
        if (cnt[c.a]++ == 0) first[c.a] = (int)i;
        ++nvalid;
    }
    if (nvalid == 0) return false;
    for (int t : M.main_tpl) {
        const Template& T = M.tpl[t];
        const int nobs = (int)T.obs.size();
        if (nobs > 16) continue;
        int slot_of[MAXANCH];
        for (int a = 0; a < M.A; ++a) slot_of[a] = -1;
        bool dup = false;
        for (int o = 0; o < nobs; ++o) {
            const TElem& e = T.el[T.obs[o]];
            if (e.kind != EK_ANCHOR) continue;
            if (slot_of[e.anc] >= 0) dup = true;
            slot_of[e.anc] = o;
        }
        if (dup) continue;
        bool ok = true;
        for (int a = 0; a < M.A && ok; ++a)
            if (cnt[a] > 0 && (slot_of[a] < 0 || cnt[a] > 1)) ok = false;
        if (!ok) continue;
        int cl[16], ostart[16], oend[16];
        const int st0 = M.open_state[t] - M.st_obs[M.open_state[t]];
        int maxnx = 0, nel = 0;
        bool cert = false, weak_seen = false;
        for (int o = 0; o < nobs; ++o) {
            cl[o] = -1; ostart[o] = oend[o] = -1;
            const TElem& e = T.el[T.obs[o]];
            if (e.kind != EK_ANCHOR || cnt[e.anc] == 0) continue;
            const seed_cluster& c = S.sc[first[e.anc]];
            const int st = st0 + o;
            if (M.head_fix[st] >= 0 && c.s > M.zone5 + M.head_fix[st]) { ok = false; break; }
            if (M.tail_fix[st] >= 0 && c.e < n - M.zone3 - M.tail_fix[st]) { ok = false; break; }
            cl[o] = first[e.anc];
            ostart[o] = c.s; oend[o] = c.e;
            maxnx = std::max(maxnx, (int)c.nx);
            ++nel;
            if (M.st_flags[st] & SF_CERT) cert = true;
            if (M.st_flags[st] & SF_WEAK) weak_seen = true;
        }
        if (!ok) continue;
        // order and tight spacing of seen anchors
        int prev = -1;
        for (int o = 0; o < nobs && ok; ++o) {
            if (cl[o] < 0) continue;
            if (prev >= 0) {
                if (oend[prev] > ostart[o] + 12) { ok = false; break; }
                const int e1 = T.obs[prev], e2 = T.obs[o];
                if (T.pins[e2] - T.pins[e1 + 1] == 0) {
                    int lo = 0, hi = 0;
                    for (int x = e1 + 1; x < e2; ++x) {
                        if (T.el[x].kind == EK_POLY) { hi += 60; continue; }
                        lo += T.el[x].lmin; hi += T.el[x].lmax;
                    }
                    const int g = ostart[o] - oend[prev];
                    if (g < lo - 8 || g > hi + 8) { ok = false; break; }
                }
            }
            prev = o;
        }
        if (!ok) continue;
        // every poly run must be explained by a poly slot of this template
        int poly_slot_seen[16];
        for (int o = 0; o < nobs; ++o) poly_slot_seen[o] = -1;
        for (size_t pi = 0; pi < S.poly.size() && ok; ++pi) {
            const poly_reg& p = S.poly[pi];
            bool expl = false;
            for (int o = 0; o < nobs && !expl; ++o) {
                const int eo = T.obs[o];
                const TElem& e = T.el[eo];
                if (e.kind != EK_POLY || e.anc != p.q) continue;
                const int st = st0 + o;
                bool decided = false;
                for (int ob = o - 1; ob >= 0; --ob) {
                    const int e1 = T.obs[ob];
                    if (T.pins[eo] - T.pins[e1 + 1] != 0) break;
                    if (cl[ob] >= 0) {
                        int lo = 0, hi = 0;
                        for (int x = e1 + 1; x < eo; ++x) { lo += T.el[x].lmin; hi += T.el[x].kind == EK_POLY ? 60 : T.el[x].lmax; }
                        const int g = p.s - oend[ob];
                        expl = g >= lo - 25 && g <= hi + 40;
                        decided = true;
                        break;
                    }
                }
                if (!decided)
                    for (int oa = o + 1; oa < nobs; ++oa) {
                        const int e2 = T.obs[oa];
                        if (T.pins[e2] - T.pins[eo + 1] != 0) break;
                        if (cl[oa] >= 0) {
                            int lo = 0, hi = 0;
                            for (int x = eo + 1; x < e2; ++x) { lo += T.el[x].lmin; hi += T.el[x].kind == EK_POLY ? 60 : T.el[x].lmax; }
                            const int g = ostart[oa] - p.e;
                            expl = g >= lo - 25 && g <= hi + 40;
                            decided = true;
                            break;
                        }
                    }
                if (!decided) {
                    if (M.head_fix[st] >= 0) expl = p.s <= M.zone5 + M.head_fix[st] + 40;
                    else if (M.tail_fix[st] >= 0) expl = p.e >= n - M.zone3 - M.tail_fix[st] - 40;
                }
                if (expl) {
                    if (poly_slot_seen[o] < 0) ++nel;
                    else if (S.poly[poly_slot_seen[o]].e - S.poly[poly_slot_seen[o]].s >= p.e - p.s) continue;
                    poly_slot_seen[o] = (int)pi;
                }
            }
            if (!expl) ok = false;
        }
        if (!ok) continue;
        if (nel < 2 || maxnx < 3) continue;
        if (!cert && !(weak_seen && nel >= 2)) continue;  // a weak outer anchor needs a second element
        // refine seen anchors with a tiny single-pattern Myers (edit distance + exact position)
        int oed[16], ced[16];
        for (int o = 0; o < nobs; ++o) {
            oed[o] = ced[o] = -1;
            if (cl[o] < 0) continue;
            const seed_cluster& c = S.sc[cl[o]];
            const Anchor& x = M.anch[c.a];
            if (c.s < 0 || c.e > n) continue;
            if (c.nx >= x.m - M.K + 1 && c.olo == 0 && c.ohi == x.m - M.K && x.nchunk == 1 && informative(x.seq) == x.m) {
                oed[o] = 0;
                continue;
            }
            const int lo = std::max(0, c.s - 6 - x.k), hi = std::min(n, c.e + 6 + x.k);
            int be = -1;
            int ed = best_hit(M.single[c.a], seq + lo, hi - lo, be);
            if (be >= 0 && x.nchunk == 1 && M.pats[x.chunk0].m == x.m) {
                oed[o] = ed;
                oend[o] = lo + be + 1;
                ostart[o] = std::max(0, oend[o] - x.m);
            } else if (be >= 0) ced[o] = ed;  // chunked anchor: edit distance of its first chunk
        }
        // certification on the refined edit distances (the full path's rules): seeds alone never make a construct
        {
            auto grade = [&](int o, bool strong) {
                const Anchor& x = M.anch[T.el[T.obs[o]].anc];
                if (oed[o] >= 0) return oed[o] <= (strong ? x.strong_ed : x.k);
                if (ced[o] >= 0) {
                    const int inf0 = M.pats[x.chunk0].inf[std::min(M.pats[x.chunk0].m, 64)];
                    return ced[o] <= (strong ? std::min(M.pats[x.chunk0].k, strong_ed_of(inf0)) : M.pats[x.chunk0].k);
                }
                return false;
            };
            bool cert2 = false;
            for (int o = 0; o < nobs && !cert2; ++o) {
                if (cl[o] < 0) continue;
                const uint16_t f = M.st_flags[st0 + o];
                if ((f & SF_CERT) && grade(o, true)) cert2 = true;
                else if ((f & SF_WEAK) && nel >= 2 && grade(o, true)) cert2 = true;
            }
            for (int o1 = 0; o1 < nobs && !cert2; ++o1) {
                const bool a1 = cl[o1] >= 0, p1 = poly_slot_seen[o1] >= 0;
                if (!a1 && !p1) continue;
                const int e1 = T.obs[o1];
                int lo = 0, hi = 0;
                for (int o2 = o1 + 1; o2 < nobs && !cert2; ++o2) {
                    const int e2 = T.obs[o2];
                    if (T.pins[e2] - T.pins[e1 + 1] != 0) break;
                    for (int x2 = (o2 == o1 + 1 ? e1 + 1 : T.obs[o2 - 1] + 1); x2 < e2; ++x2) {
                        const TElem& te = T.el[x2];
                        if (te.kind == EK_POLY) hi += 60;
                        else { lo += te.lmin; hi += te.lmax; }
                    }
                    const bool a2 = cl[o2] >= 0, p2 = poly_slot_seen[o2] >= 0;
                    if (a2 || p2) {
                        const bool anchor_ok = (a1 && grade(o1, false)) || (a2 && grade(o2, false));
                        if (anchor_ok) {
                            const int end1 = a1 ? oend[o1] : S.poly[poly_slot_seen[o1]].e;
                            const int start2 = a2 ? ostart[o2] : S.poly[poly_slot_seen[o2]].s;
                            const int g = start2 - end1, sl = (p1 || p2) ? 15 : 8;
                            if (g >= lo - sl && g <= hi + sl) cert2 = true;
                        }
                    }
                    // the next pair's gap also spans element o2
                    {
                        const TElem& te = T.el[e2];
                        if (te.kind == EK_POLY) hi += 60;
                        else { lo += te.lmin; hi += te.lmax; }
                    }
                }
            }
            if (!cert2) continue;
        }
        R.k = 1;
        R.cuts.clear(); R.segs.clear(); R.checks.clear(); R.cut_post.clear();
        R.cut_lo.clear(); R.cut_hi.clear(); R.cut_kind.clear(); R.cut_flags.clear(); R.cut_aux.clear();
        R.boundary_unresolved_seen = false;
        R.p_single = M.fast_p;
        R.strand_call = T.strand;
        R.flags |= RES_FAST_PATH;
        Segment sg{0, n, T.strand, 0, M.fast_p};
        int o_open = -1, o_close = -1;
        for (int o = 0; o < nobs; ++o) if (ostart[o] >= 0 || poly_slot_seen[o] >= 0) { if (o_open < 0) o_open = o; o_close = o; }
        auto pstart = [&](int o) { return cl[o] >= 0 ? ostart[o] : S.poly[poly_slot_seen[o]].s; };
        auto pend = [&](int o) { return cl[o] >= 0 ? oend[o] : S.poly[poly_slot_seen[o]].e; };
        if (o_open != 0) sg.flags |= SEG_PARTIAL_LEFT;
        sg.start = retained_template_start(M, st0 + o_open, pstart(o_open));
        int cs_slack = 0, ce_slack = 0;  // construct bounds defaulted to the read ends (see expected_window)
        if (!template_edge_span(M, st0 + o_open, true).bounded) cs_slack = M.zone5;
        if (o_close != nobs - 1) sg.flags |= SEG_PARTIAL_RIGHT;
        sg.end = retained_template_end(M, st0 + o_close, pend(o_close), n);
        if (!template_edge_span(M, st0 + o_close, false).bounded) ce_slack = M.zone3;
        if (sg.end <= sg.start) { sg.start = 0; sg.end = n; }
        if (M.tpl[t].art) sg.flags |= SEG_ARTIFACT;
        R.segs.push_back(sg);
        if (M.opt.check_regions) {
            for (int o = 0; o < nobs; ++o) {
                const TElem& el = T.el[T.obs[o]];
                if (el.kind == EK_POLY) {
                    if (poly_slot_seen[o] >= 0) {
                        const poly_reg& p = S.poly[poly_slot_seen[o]];
                        add_check(R, M, el, T.strand, 0, Status::FULL, p.s - M.opt.seen_pad, p.e + M.opt.seen_pad, -1, -1, -1, M.fast_p, n);
                    } else {
                        int ws, we, spr, es;
                        int ps2[16], pe2[16];
                        for (int x = 0; x < nobs; ++x) { ps2[x] = (cl[x] >= 0) ? ostart[x] : (poly_slot_seen[x] >= 0 ? S.poly[poly_slot_seen[x]].s : -1); pe2[x] = (ps2[x] >= 0) ? pend(x) : -1; }
                        expected_window(M, t, o, ps2, pe2, sg.start, sg.end, M.opt.check_pad, ws, we, spr, &es, cs_slack, ce_slack);
                        missing_slot_check(R, M, el, T.strand, 0, seq, n, ws, we, es, 0.5f * M.fast_p, 0.5f);
                    }
                    continue;
                }
                const Anchor& x = M.anch[el.anc];
                if (cl[o] >= 0) {
                    const seed_cluster& c = S.sc[cl[o]];
                    if (c.e > n || c.s < 0) {
                        int rlo, rhi, ed;
                        seed_truncation(M, c, seq, n, rlo, rhi, ed);
                        if (c.e > n) add_check(R, M, el, T.strand, 0, Status::TRUNCATED_AT_READ_END, c.s - M.opt.seen_pad, n, rlo, rhi, ed, M.fast_p, n);
                        else add_check(R, M, el, T.strand, 0, Status::TRUNCATED_AT_READ_END, 0, c.e + M.opt.seen_pad, rlo, rhi, ed, M.fast_p, n);
                    } else if (oed[o] >= 0 && oed[o] <= x.strong_ed) {
                        add_check(R, M, el, T.strand, 0, Status::FULL, ostart[o] - M.opt.seen_pad, oend[o] + M.opt.seen_pad, 0, x.m, oed[o], M.fast_p, n);
                    } else if (oed[o] >= 0 && oed[o] <= x.k) {
                        int fws, fwe, flo, fhi, fed;
                        if (degraded_as_fragment(M, el.anc, seq, n, ostart[o], oend[o], oed[o], fws, fwe, flo, fhi, fed))
                            add_check(R, M, el, T.strand, 0, Status::PARTIAL, fws, fwe, flo, fhi, fed, M.fast_p, n);
                        else
                            add_check(R, M, el, T.strand, 0, Status::PARTIAL, ostart[o] - M.opt.seen_pad, oend[o] + M.opt.seen_pad, 0, x.m, oed[o], M.fast_p, n);
                    } else {
                        add_check(R, M, el, T.strand, 0, Status::PARTIAL, c.s - M.opt.seen_pad, c.e + M.opt.seen_pad, c.olo, std::min(x.m, c.ohi + M.K), oed[o], M.fast_p, n);
                    }
                    continue;
                }
                int ps2[16], pe2[16];
                for (int y = 0; y < nobs; ++y) { ps2[y] = (cl[y] >= 0) ? ostart[y] : (poly_slot_seen[y] >= 0 ? S.poly[poly_slot_seen[y]].s : -1); pe2[y] = (ps2[y] >= 0) ? pend(y) : -1; }
                int ws, we, spr, es;
                expected_window(M, t, o, ps2, pe2, sg.start, sg.end, M.opt.check_pad, ws, we, spr, &es, cs_slack, ce_slack);
                if (ws >= n || we <= 0) continue;
                // a closer expected at / beyond the read end: test a truncated prefix at the read end
                if (x.closing && o == nobs - 1 && we >= n - 4 && x.nchunk == 1) {
                    const int lo = std::max(0, n - x.m - 8);
                    int be;
                    uint64_t Pv, Mv;
                    best_hit(M.single[el.anc], seq + lo, n - lo, be, &Pv, &Mv);
                    int d = 0;
                    const int L = prefix_from_column(M.pats[x.chunk0], Pv, Mv, 0, d);
                    if (L >= PARTIAL_MIN) {
                        add_check(R, M, el, T.strand, 0, Status::TRUNCATED_AT_READ_END, n - L - d - M.opt.seen_pad, n, 0, L, d, M.fast_p, n);
                        if (sg.flags & SEG_PARTIAL_RIGHT) { R.segs[0].end = n; }
                        continue;
                    }
                }
                missing_slot_check(R, M, el, T.strand, 0, seq, n, ws, we, es, 0.5f * M.fast_p, 0.5f);
            }
        }
        return true;
    }
    return false;
}

}  // namespace detail

namespace detail {

inline float gap_score(const Model& M, const desc_hot& d, int gap) {
    const int x = gap - d.fixed;
    if (d.kind == 0) {
        const int lo = M.tight_lo[d.tid];
        if (x < lo || x > M.tight_hi[d.tid]) return NEG;
        return M.tight[d.tid][x - lo];
    }
    const std::vector<float>& T = M.ins[d.kind - 1];
    const int xi = x - INS_LO;
    if (xi < 0) return NEG;
    if (xi >= (int)T.size()) return T.back() - 0.002f * (xi - (int)T.size());
    return T[xi];
}
// log density (per bp) of an insert of length x under the one- (kind 1) / two-insert (kind 2) length tables
inline float gap_score_ins(const Model& M, int kind, int x) {
    const std::vector<float>& T = M.ins[kind - 1];
    const int xi = x - INS_LO;
    if (xi < 0) return NEG;
    return xi < (int)T.size() ? T[xi] : T.back();
}
inline float fexp(float x) {  // fast exp, rel. error < 2e-6 on [-87, 88]
    if (!(x >= -87.f)) return 0.f;  // also NaN (never reaches the int conversion below)
    if (x > 88.f) x = 88.f;
    const float t = x * 1.4426950408889634f;
    const float fi = std::floor(t);
    const float f = t - fi;
    float p = 1.0f + f * (0.69314718f + f * (0.24022651f + f * (0.05550411f + f * (0.00961813f + f * 0.00133336f))));
    int32_t ii = (int32_t)fi;
    uint32_t b;
    memcpy(&b, &p, 4);
    b += (uint32_t)ii << 23;
    memcpy(&p, &b, 4);
    return p;
}
inline bool cert_alone(const Model& M, int st, const Event& e) { return (M.st_flags[st] & SF_CERT) && e.cls == CL_S4; }

// Same-construct transition (j,sp) -> (i,st): weight and the certification it adds.
inline bool tr_same(const Model& M, const Event& ej, int sp, const Event& ei, int st, int gap, int L, float& w, uint8_t& hadd) {
    const desc_hot& d = M.same[sp * M.S + st];
    if (!d.valid || ei.cls == CL_S8) return false;
    const float g = gap_score(M, d, gap);
    if (g <= NEG_HALF) return false;
    w = g + d.base;
    uint8_t h = cert_alone(M, st, ei) ? 1 : 0;
    if (!h) {
        const bool aj = ej.cls != CL_POLY, ai = ei.cls != CL_POLY;
        if (aj || ai) {
            if (d.kind == 0) {
                h = std::abs(gap - d.fixed) <= d.slack;
            } else {
                // across an insert: a weak or short anchor needs an anchor-grade partner or a poly tail
                auto weakish = [&](const Event& e, int s) {
                    if (e.cls == CL_S4 && (M.st_flags[s] & SF_WEAK)) return true;
                    if ((M.st_flags[s] & SF_SHORT) && e.cls == CL_S4 && e.ed >= 0 && e.ed <= 1) return true;
                    // a degraded (ED5-6) long anchor at the layout-implied distance from the physical read start / end
                    if (e.cls == CL_S6 && !(M.st_flags[s] & (SF_SHORT | SF_WEAK))) {
                        if (M.head_fix[s] >= 0) { const int s0 = e.start - M.head_fix[s]; if (s0 >= -4 && s0 <= M.zone5) return true; }
                        if (M.tail_fix[s] >= 0) { const int r = L - (e.end + M.tail_fix[s]); if (r >= -4 && r <= M.zone3) return true; }
                    }
                    return false;
                };
                const bool pj = ej.cls == CL_S4 || ej.cls == CL_POLY, pi = ei.cls == CL_S4 || ei.cls == CL_POLY;
                h = (weakish(ej, sp) && pi) || (weakish(ei, st) && pj);
            }
        }
    }
    // a short anchor tightly paired with another anchor: full LLR
    if (d.kind == 0 && ((M.st_flags[st] | M.st_flags[sp]) & SF_SHORT) && ej.cls != CL_POLY && ei.cls != CL_POLY &&
        std::abs(gap - d.fixed) <= d.slack) {
        if (M.st_flags[st] & SF_SHORT) w += M.emit_bonus[st * NFEAT + ei.feat];
        if (M.st_flags[sp] & SF_SHORT) w += M.emit_bonus[sp * NFEAT + ej.feat];
    }
    hadd = h;
    return true;
}
// Junction transition (j,sp) -> new construct at (i,st).
inline bool tr_cross(const Model& M, const Event& ej, int sp, const Event& ei, int st, int gap, float& w, uint8_t& hadd, uint8_t& h0ok) {
    const desc_hot& d = M.cross[sp * M.S + st];
    if (!d.valid) return false;
    const bool strong_closer = (M.st_flags[sp] & SF_CLOSE) && ej.cls == CL_S4 && !(M.st_flags[sp] & SF_WEAK);
    const bool junc = strong_closer && (M.st_flags[st] & SF_OPEN) && gap >= -10 && gap <= 10;
    // an abutting junction pair (closer ED <= budget, opener anchor-grade) also certifies the construct it closes
    h0ok = (M.st_flags[sp] & SF_CLOSE) && (ej.cls == CL_S4 || ej.cls == CL_S6) && (M.st_flags[st] & SF_OPEN) && ei.cls == CL_S4 &&
           !(M.st_flags[st] & (SF_SHORT | SF_WEAK)) && !(M.st_flags[sp] & SF_SHORT) && gap >= -8 && gap <= 8;
    if (ei.cls == CL_S8 && !junc) return false;
    const float g = gap_score(M, d, gap);
    if (g <= NEG_HALF) return false;
    w = g + d.base;
    if (d.pal) {
        const float e0 = M.emit[st * NFEAT + ei.feat], e1 = M.emit_pal[st * NFEAT + ei.feat];
        if (e1 > NEG_HALF && e0 > NEG_HALF) w += e1 - e0;
    }
    bool h = cert_alone(M, st, ei) || (junc && !(M.st_flags[st] & SF_SHORT) && (ei.cls == CL_S4 || ei.cls == CL_S6));
    if (!h && d.kind == 0 && (ej.cls == CL_S4 || ej.cls == CL_S6) && (M.st_flags[sp] & SF_ANCHOR) && (M.st_flags[st] & SF_ANCHOR)) {
        // junction geometry: the two constructs' anchors joined by fixed-length elements only
        const int x = gap - d.fixed;
        const bool geo = x >= -d.slack && x <= d.slack + 60;
        // the same with barcode-adjacent primers deleted: the blocks abut within [DEL_GMIN, DEL_GMAX_CERT] bp
        bool geo_del = false;
        if (!geo && (d.dela | d.delb)) {
            const int dl[3] = {d.dela, d.delb, d.dela + d.delb};
            for (int q = 0; q < 3 && !geo_del; ++q) geo_del = dl[q] > 0 && x + dl[q] >= DEL_GMIN && x + dl[q] <= DEL_GMAX_CERT;
        }
        if (geo || geo_del) {
            if (M.st_flags[st] & SF_SHORT) {
                h = ei.cls == CL_S4 && ei.ed >= 0 && ei.ed <= 1;
                if (h) w += M.emit_bonus[st * NFEAT + ei.feat];
            } else h = ei.cls == CL_S4 || ei.cls == CL_S6;
            // two facing short partner anchors (ED <= 1) also certify the closed construct; both get the full LLR
            if (h && geo_del && (M.st_flags[st] & SF_SHORT) && (M.st_flags[sp] & SF_SHORT) && ej.cls == CL_S4 && ej.ed >= 0 && ej.ed <= 1) {
                h0ok = 1;
                w += M.emit_bonus[sp * NFEAT + ej.feat];
            }
        }
    }
    if (!h && d.kind == 0 && strong_closer && ei.cls == CL_POLY && std::abs(gap - d.fixed) <= d.slack) {
        // one-sided junction: a FULL closer, the next construct's opener lost, its poly tail at the barcode offset
        h = true;
    }
    hadd = h;
    return true;
}
// An exact short anchor at its layout distance from the read start / end certifies a terminal construct.
inline bool begin_geom(const Model& M, int st, const Event& e) {
    if (!(M.st_flags[st] & SF_SHORT) || e.cls != CL_S4 || e.ed != 0 || M.head_fix[st] < 0) return false;
    const int s0 = e.start - M.head_fix[st];
    return s0 >= -4 && s0 <= M.zone5;
}
inline bool end_geom(const Model& M, int st, const Event& e, int L) {
    if (!(M.st_flags[st] & SF_SHORT) || e.cls != CL_S4 || e.ed != 0 || M.tail_fix[st] < 0) return false;
    const int r = L - (e.end + M.tail_fix[st]);
    return r >= -4 && r <= M.zone3;
}
inline int head_start(const Model& M, int st, const Event& e, int lo) {
    return M.head_fix[st] >= 0 ? std::max(lo, e.start - M.head_fix[st]) : lo;
}

// Viterbi over cells (event, state, h = construct certified); junctions leave, and the read ends in, certified
// cells only. Returns k (0 = no construct) and fills the path.
inline int viterbi(const Model& M, const Event* ev, int n, int L, Scratch& W) {
    const int S = M.S;
    W.path_ev.clear(); W.path_st.clear(); W.path_cross.clear();
    W.null_score = M.lp0;
    if (n == 0) { W.best_score = M.lp0; return 0; }
    const size_t NN = (size_t)n * S * 2;
    if (W.dp.size() < NN) { W.dp.resize(NN); W.bp.resize(NN); W.cs.resize(NN); }
    std::fill(W.dp.begin(), W.dp.begin() + NN, NEG);
    float* dp = W.dp.data();
    int32_t* bp = W.bp.data();
    int32_t* cs = W.cs.data();
    W.tr.clear();
    const float log_l = std::log((float)std::max(L, 100));
    float best = NEG;
    int bi = -1, bs = -1, bh = 1;
    const bool any_force = !W.force.empty();
    for (int i = 0; i < n; ++i) {
        const Event& ei = ev[i];
        for (int st : M.states_by_type[ei.type]) {
            const float e = M.emit[st * NFEAT + ei.feat];
            if (e <= NEG_HALF) continue;
            float cur[2] = {NEG, NEG};
            int32_t cb[2] = {-1, -1}, cc[2] = {0, 0};
            if (ei.cls != CL_S8) {
                const bool bg = begin_geom(M, st, ei);
                const int h = (cert_alone(M, st, ei) || bg) ? 1 : 0;
                cur[h] = M.lbegin[st] - log_l + (bg ? M.emit_bonus[st * NFEAT + ei.feat] : 0.f);
                cb[h] = -1;
                cc[h] = (head_start(M, st, ei, 0) << 2) | 1;
            }
            const int j0 = std::max(0, i - PRED_WINDOW);
            for (int j = i - 1; j >= j0; --j) {
                const Event& ej = ev[j];
                const int gap = ei.start - ej.end;
                if (gap < -MAX_OVERLAP) continue;
                for (int sp : M.states_by_type[ej.type]) {
                    const size_t pc = ((size_t)j * S + sp) * 2;
                    const float p0 = dp[pc], p1 = dp[pc + 1];
                    if (p0 <= NEG_HALF && p1 <= NEG_HALF) continue;
                    float w;
                    uint8_t ha;
                    if (tr_same(M, ej, sp, ei, st, gap, L, w, ha)) {
                        W.tr.push_back(Scratch::Tr{(int32_t)(j * S + sp), (int32_t)(i * S + st), w, 0, ha});
                        for (int hp = 0; hp < 2; ++hp) {
                            const float pv = hp ? p1 : p0;
                            if (pv <= NEG_HALF) continue;
                            const int h = hp | ha;
                            const float v = pv + w;
                            if (v > cur[h]) {
                                cur[h] = v;
                                cb[h] = (j << 8) | (sp << 2) | (hp << 1);
                                // an anchor-grade element of the construct clears the "no own anchor" bit
                                cc[h] = (ei.cls == CL_S4 || ei.cls == CL_S6) ? (cs[pc + hp] & ~2) : cs[pc + hp];
                            }
                        }
                    }
                    uint8_t h0;
                    if (tr_cross(M, ej, sp, ei, st, gap, w, ha, h0)) {
                        if (any_force)
                            for (const auto& fp : W.force)
                                if (fp.first == j && fp.second == i) { w += FORCE_BONUS; ha = 1; h0 = 1; }
                        bool rec = false;
                        for (int hp = h0 ? 0 : 1; hp < 2; ++hp) {
                            const float pv = hp ? p1 : p0;
                            if (pv <= NEG_HALF) continue;
                            const int c = cs[pc + hp];
                            const int start = c >> 2;
                            const int endc = M.tail_fix[sp] >= 0 ? ej.end + M.tail_fix[sp] : ei.start;
                            if (endc - start < ((c & 1) ? M.min_terminal : M.min_interior)) continue;
                            if ((c & 2) && ei.start - ej.end < M.min_terminal && (ej.cls == CL_POLY || ej.cls == CL_S8)) continue;
                            if (!rec) { W.tr.push_back(Scratch::Tr{(int32_t)(j * S + sp), (int32_t)(i * S + st), w, (uint8_t)(1 | (h0 ? 2 : 0)), ha}); rec = true; }
                            const float v = pv + w;
                            if (v > cur[ha]) {
                                cur[ha] = v;
                                cb[ha] = (j << 8) | (sp << 2) | (hp << 1) | 1;
                                cc[ha] = (head_start(M, st, ei, ej.end) << 2) | ((ei.cls == CL_S4 || ei.cls == CL_S6) ? 0 : 2);
                            }
                        }
                    }
                }
            }
            const size_t cc0 = ((size_t)i * S + st) * 2;
            for (int h = 0; h < 2; ++h) {
                if (cur[h] <= NEG_HALF) continue;
                dp[cc0 + h] = cur[h] + e;
                bp[cc0 + h] = cb[h];
                cs[cc0 + h] = cc[h];
            }
            for (int h = end_geom(M, st, ei, L) ? 0 : 1; h < 2; ++h) {
                if (cur[h] <= NEG_HALF) continue;
                const int start = cc[h] >> 2;
                const int endc = (M.st_flags[st] & SF_CLOSE) ? ei.end : (M.tail_fix[st] >= 0 ? std::min(L, ei.end + M.tail_fix[st]) : L);
                // a junction-opened construct without an anchor-grade element of its own must reach into cDNA
                const bool stub = (cc[h] & 2) && (ei.cls == CL_POLY || ei.cls == CL_S8) && L - ei.end < M.min_terminal;
                if (endc - start >= M.min_terminal && !stub) {
                    const float fin = dp[cc0 + h] + M.lend[st] + M.lp1 + (end_geom(M, st, ei, L) ? M.emit_bonus[st * NFEAT + ei.feat] : 0.f);
                    if (fin > best) { best = fin; bi = i; bs = st; bh = h; }
                }
            }
        }
    }
    if (bi < 0 || best <= M.lp0) { W.best_score = M.lp0; return 0; }
    W.best_score = best;
    int i = bi, st = bs, h = bh;
    while (true) {
        const int32_t b = bp[((size_t)i * S + st) * 2 + h];
        W.path_ev.push_back(i);
        W.path_st.push_back(st);
        if (b < 0) { W.path_cross.push_back(1); break; }
        W.path_cross.push_back((uint8_t)(b & 1));
        i = b >> 8; st = (b >> 2) & 63; h = (b >> 1) & 1;
    }
    std::reverse(W.path_ev.begin(), W.path_ev.end());
    std::reverse(W.path_st.begin(), W.path_st.end());
    std::reverse(W.path_cross.begin(), W.path_cross.end());
    int k = 0;
    for (auto c : W.path_cross) k += c;
    return k;
}

// Scaled forward-backward over the same graph with a count layer (k = 1 / k >= 2): fills W.l_z, W.l_z1
// and W.jpost[x] = P(a junction lies between event x and x+1). Requires viterbi() to have run.
inline void forward_backward(const Model& M, const Event* ev, int n, int L, Scratch& W, bool backward) {
    const int S = M.S;
    W.l_z = M.lp0; W.l_z1 = NEG;
    W.jpost.assign(n + 1, 0.f);
    if (n == 0) return;
    const size_t NC = (size_t)n * S;
    if (W.fw.size() < NC * 4) { W.fw.resize(NC * 4); W.bw.resize(NC * 4); }
    if (W.scale.size() < NC) W.scale.resize(NC);
    float* F = W.fw.data();
    float* B = W.bw.data();
    float* scale = W.scale.data();
    const float* dp = W.dp.data();
    const int32_t* cs = W.cs.data();
    std::fill(F, F + NC * 4, 0.f);
    for (size_t c = 0; c < NC; ++c) scale[c] = std::max(dp[c * 2], dp[c * 2 + 1]);
    const float log_l = std::log((float)std::max(L, 100));
    const float V = W.best_score;
    for (int i = 0; i < n; ++i) {
        const Event& ei = ev[i];
        if (ei.cls == CL_S8) continue;
        for (int st : M.states_by_type[ei.type]) {
            const size_t c = (size_t)i * S + st;
            const float e = M.emit[st * NFEAT + ei.feat];
            if (e <= NEG_HALF || scale[c] <= NEG_HALF) continue;
            const bool bg = begin_geom(M, st, ei);
            const int h = (cert_alone(M, st, ei) || bg) ? 1 : 0;
            F[c * 4 + h * 2] = fexp(M.lbegin[st] - log_l + e + (bg ? M.emit_bonus[st * NFEAT + ei.feat] : 0.f) - scale[c]);
        }
    }
    // transitions (grouped by target, targets ascending): m = exp(w + e_dst + scale_src - scale_dst)
    const size_t NT = W.tr.size();
    if (W.trm.size() < NT) W.trm.resize(NT);
    float* TM = W.trm.data();
    for (size_t x = 0; x < NT; ++x) {
        const Scratch::Tr& t = W.tr[x];
        const int i = t.dst / S, st = t.dst % S;
        const float e = M.emit[st * NFEAT + ev[i].feat];
        const float m = fexp(t.w + e + scale[t.src] - scale[t.dst]);
        TM[x] = m;
        const float* fp = &F[(size_t)t.src * 4];
        float* f = &F[(size_t)t.dst * 4];
        if (!t.cross) {
            for (int hp = 0; hp < 2; ++hp) {
                const int h = hp | t.ha;
                f[h * 2 + 0] += fp[hp * 2 + 0] * m;
                f[h * 2 + 1] += fp[hp * 2 + 1] * m;
            }
        } else f[t.ha * 2 + 1] += (fp[2] + fp[3] + ((t.cross & 2) ? fp[0] + fp[1] : 0.f)) * m;
    }
    auto endw = [&](int i, int st, int h) -> float {
        const Event& ei = ev[i];
        const int c = cs[((size_t)i * S + st) * 2 + h];
        const int start = c >> 2;
        const int endc = (M.st_flags[st] & SF_CLOSE) ? ei.end : (M.tail_fix[st] >= 0 ? std::min(L, ei.end + M.tail_fix[st]) : L);
        if (endc - start < M.min_terminal) return NEG;
        if ((c & 2) && (ei.cls == CL_POLY || ei.cls == CL_S8) && L - ei.end < M.min_terminal) return NEG;
        return M.lend[st] + M.lp1 + (end_geom(M, st, ei, L) ? M.emit_bonus[st * NFEAT + ei.feat] : 0.f);
    };
    double zsum = 0, z1sum = 0;
    for (int i = 0; i < n; ++i)
        for (int st : M.states_by_type[ev[i].type]) {
            const size_t c = (size_t)i * S + st;
            if (scale[c] <= NEG_HALF) continue;
            const bool eg = end_geom(M, st, ev[i], L);
            for (int h = eg ? 0 : 1; h < 2; ++h) {
                const float ew = endw(i, st, h);
                if (ew <= NEG_HALF) continue;
                const double m = std::exp(std::min(700.0, (double)(scale[c] + ew - V)));
                zsum += (F[c * 4 + h * 2] + F[c * 4 + h * 2 + 1]) * m;
                z1sum += F[c * 4 + h * 2] * m;
            }
        }
    const double z0 = std::exp(std::min(700.0, (double)(M.lp0 - V)));
    W.l_z = V + (float)std::log(std::max(zsum + z0, 1e-300));
    W.l_z1 = z1sum > 0 ? V + (float)std::log(z1sum) : NEG;
    if (!backward || !M.opt.posteriors) return;
    // backward: B(c,h,l) = exp(b - (V - scale_c)); END terms, then transitions in reverse order (push to source)
    std::fill(B, B + NC * 4, 0.f);
    for (int i = 0; i < n; ++i)
        for (int st : M.states_by_type[ev[i].type]) {
            const size_t c = (size_t)i * S + st;
            if (scale[c] <= NEG_HALF) continue;
            const bool eg = end_geom(M, st, ev[i], L);
            for (int h = eg ? 0 : 1; h < 2; ++h) {
                const float ew = endw(i, st, h);
                if (ew <= NEG_HALF) continue;
                const float m = fexp(ew - V + scale[c]);
                B[c * 4 + h * 2] = m;
                B[c * 4 + h * 2 + 1] = m;
            }
        }
    const double zf = std::exp(std::max(-700.0, std::min(700.0, (double)(V - W.l_z))));
    for (size_t x = NT; x-- > 0;) {
        const Scratch::Tr& t = W.tr[x];
        const float m = TM[x];
        const float* b2 = &B[(size_t)t.dst * 4];
        float* b = &B[(size_t)t.src * 4];
        if (!t.cross) {
            for (int hp = 0; hp < 2; ++hp) {
                const int h = hp | t.ha;
                b[hp * 2 + 0] += b2[h * 2 + 0] * m;
                b[hp * 2 + 1] += b2[h * 2 + 1] * m;
            }
        } else {
            const float v = b2[t.ha * 2 + 1] * m;
            b[2] += v;
            b[3] += v;
            if (t.cross & 2) { b[0] += v; b[1] += v; }
            const float* fs = &F[(size_t)t.src * 4];
            const double mass = (double)(fs[2] + fs[3] + ((t.cross & 2) ? fs[0] + fs[1] : 0.f)) * v * zf;
            if (mass > 0) {
                const int i = t.src / S, i2 = t.dst / S;
                W.jpost[i] += (float)mass;
                W.jpost[i2] -= (float)mass;
            }
        }
    }
    float acc = 0;
    for (int x = 0; x < n; ++x) { acc += W.jpost[x]; W.jpost[x] = std::min(1.f, std::max(0.f, acc)); }
}

}  // namespace detail

namespace detail {

// Palindrome centre of s[a, b) from reverse-complement k-mer votes (x + y + k = 2 * centre). True when the votes
// span a long arm; `relaxed` also accepts a dense shorter span.
inline bool foldback_centre(const char* s, int a, int b, Scratch& W, int& centre, float& frac, bool relaxed = false) {
    const int k = 10;
    const int len = b - a;
    centre = (a + b) / 2;
    frac = 0;
    if (len < 80) return false;
    const int step = 2;
    size_t hs = 64;
    while (hs < (size_t)(2 * len / step + 8)) hs <<= 1;
    W.fb_keys.assign(hs, 0);
    W.fb_pos.resize(hs);
    const uint32_t hm = (uint32_t)hs - 1, km = (1u << (2 * k)) - 1;
    auto lowc = [&](uint32_t v) {
        int c[4] = {0, 0, 0, 0};
        for (int i = 0; i < k; ++i) { c[v & 3]++; v >>= 2; }
        return std::max(std::max(c[0], c[1]), std::max(c[2], c[3])) >= 7;
    };
    uint32_t f = 0, r = 0;
    int valid = 0;
    for (int i = a; i < b; ++i) {
        const unsigned char ch = (unsigned char)s[i];
        f = ((f << 2) | code2(ch)) & km;
        valid = (ch | 0x20) == 'n' ? 0 : valid + 1;
        if (valid >= k && ((i - a) % step) == 0 && !lowc(f)) {
            uint32_t h = ((f * 2654435761u) >> 8) & hm;
            while (W.fb_keys[h] != 0) h = (h + 1) & hm;
            W.fb_keys[h] = f + 1;
            W.fb_pos[h] = i - k + 1;
        }
    }
    // reverse-complement lookups (complement of code c is 3 - c with A0 C1 G2 T3)
    W.fb_diag.clear();
    W.fb_hist.assign(len / 4 + 2, 0);
    valid = 0;
    for (int i = a; i < b; ++i) {
        const unsigned char ch = (unsigned char)s[i];
        r = (r >> 2) | ((3u - code2(ch)) << (2 * k - 2));
        valid = (ch | 0x20) == 'n' ? 0 : valid + 1;
        if (valid < k) continue;
        const int y = i - k + 1;
        uint32_t h = ((r * 2654435761u) >> 8) & hm;
        while (W.fb_keys[h] != 0) {
            if (W.fb_keys[h] == r + 1) {
                const int x = W.fb_pos[h];
                if (x + k <= y) {
                    const int c2 = (x + y + k) / 2;
                    const int bin = (c2 - a) / 4;
                    if (bin >= 0 && bin < (int)W.fb_hist.size()) {
                        W.fb_hist[bin]++;
                        W.fb_diag.push_back(c2);
                        W.fb_diag.push_back(x);
                    }
                }
            }
            h = (h + 1) & hm;
        }
    }
    if (W.fb_diag.size() < 16) return false;
    int bestc = 0, bestb = -1;
    const int nb = (int)W.fb_hist.size();
    for (int x = 0; x < nb; ++x) {
        int c = 0;
        for (int d = -2; d <= 2; ++d)
            if (x + d >= 0 && x + d < nb) c += W.fb_hist[x + d];
        if (c > bestc) { bestc = c; bestb = x; }
    }
    if (bestb < 0) return false;
    const int ctr0 = a + bestb * 4 + 2;
    int cnt = 0, xmin = 1 << 30, xmax = -1;
    long csum = 0;
    for (size_t q = 0; q < W.fb_diag.size(); q += 2) {
        if (std::abs(W.fb_diag[q] - ctr0) > 10) continue;
        ++cnt;
        csum += W.fb_diag[q];
        xmin = std::min(xmin, W.fb_diag[q + 1]);
        xmax = std::max(xmax, W.fb_diag[q + 1]);
    }
    if (cnt == 0) return false;
    const int ctr = (int)(csum / cnt);
    const int arm = std::min(ctr - a, b - ctr);
    const int span = xmax - xmin + k;
    frac = (float)span / std::max(1, arm);
    centre = ctr;
    if (cnt >= 20 && arm >= 50 && span >= std::max(60, (int)(0.4 * arm)) && cnt >= 0.17 * span) return true;
    return relaxed && cnt >= 16 && arm >= 50 && span >= 60 && cnt >= 0.22 * span;
}

// Jaccard index of the FOLD_K-mer sets of the left arm s[a, c) and the reverse complement of the right arm s[c, b).
inline float arm_jaccard(const char* s, int a, int c, int b, Scratch& W) {
    const int k = FOLD_K;
    if (c - a < FOLD_ARM_MIN || b - c < FOLD_ARM_MIN) return 0.f;
    const uint32_t km = (1u << (2 * k)) - 1;
    auto collect = [&](std::vector<uint32_t>& out, int lo, int hi, bool rc) {
        out.clear();
        uint32_t v = 0;
        int valid = 0;
        if (!rc) {
            for (int i = lo; i < hi; ++i) {
                const unsigned char ch = (unsigned char)s[i];
                v = ((v << 2) | code2(ch)) & km;
                valid = (ch | 0x20) == 'n' ? 0 : valid + 1;
                if (valid >= k) out.push_back(v);
            }
        } else {  // reverse complement, read from the right end: complement of code x is 3 - x
            for (int i = hi - 1; i >= lo; --i) {
                const unsigned char ch = (unsigned char)s[i];
                v = ((v << 2) | (3u - code2(ch))) & km;
                valid = (ch | 0x20) == 'n' ? 0 : valid + 1;
                if (valid >= k) out.push_back(v);
            }
        }
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
    };
    collect(W.jk_a, a, c, false);
    collect(W.jk_b, c, b, true);
    if (W.jk_a.empty() || W.jk_b.empty()) return 0.f;
    size_t i = 0, j = 0, inter = 0;
    while (i < W.jk_a.size() && j < W.jk_b.size()) {
        if (W.jk_a[i] == W.jk_b[j]) { ++inter; ++i; ++j; }
        else if (W.jk_a[i] < W.jk_b[j]) ++i;
        else ++j;
    }
    const size_t uni = W.jk_a.size() + W.jk_b.size() - inter;
    return uni ? (float)inter / (float)uni : 0.f;
}

// 1 - banded edit distance / longer length of s[l0, l1) against rc(s[r0, r1)); band = |length difference| + 32.
inline float seg_identity(const char* s, int l0, int l1, int r0, int r1, Scratch& W) {
    const int n = l1 - l0, m = r1 - r0;
    if (n <= 0 || m <= 0) return 0.f;
    const int band = std::abs(n - m) + 32;
    const int INF = 1 << 28;
    std::vector<int32_t>& prev = W.nw0;
    std::vector<int32_t>& cur = W.nw1;
    prev.assign((size_t)m + 1, INF);
    cur.assign((size_t)m + 1, INF);
    for (int j = 0; j <= std::min(m, band); ++j) prev[j] = j;
    auto rcb = [&](int j) -> int { return 3 - code2((unsigned char)s[r1 - j]); };  // j-th base (1-based) of rc(s[r0, r1))
    for (int i = 1; i <= n; ++i) {
        const int jlo = std::max(1, i - band), jhi = std::min(m, i + band);
        std::fill(cur.begin(), cur.end(), INF);
        if (i - band <= 0) cur[0] = i;
        const int ci = code2((unsigned char)s[l0 + i - 1]);
        for (int j = jlo; j <= jhi; ++j) {
            int v = prev[j - 1] + (ci == rcb(j) ? 0 : 1);
            v = std::min(v, prev[j] + 1);
            v = std::min(v, cur[j - 1] + 1);
            cur[j] = v;
        }
        std::swap(prev, cur);
    }
    const int ed = prev[m] >= INF ? std::max(n, m) : prev[m];
    return 1.f - (float)ed / (float)std::max(n, m);
}

// Dominant anti-diagonal of the k-mers shared by the left arm and rc(right arm), positions counted from the arms'
// unit-side ends: gets its first / last positions and k-mer count `non`; false when non < FOLD_DIAG_MIN.
inline bool arm_mirror_reach(const char* s, int a, int c, int b, Scratch& W, int& min_i, int& min_j, int& max_i, int& max_j, int& non) {
    const int k = FOLD_K;
    min_i = min_j = 1 << 30; max_i = max_j = -1; non = 0;
    const int n_l = c - a, n_r = b - c;
    if (n_l < k || n_r < k) return false;
    const uint32_t km = (1u << (2 * k)) - 1;
    auto collect = [&](std::vector<std::pair<uint32_t, int32_t>>& out, int lo, int hi, bool rc) {
        out.clear();
        uint32_t v = 0;
        int valid = 0, pos = 0;
        if (!rc) {
            for (int i = lo; i < hi; ++i, ++pos) {
                const unsigned char ch = (unsigned char)s[i];
                v = ((v << 2) | code2(ch)) & km;
                valid = (ch | 0x20) == 'n' ? 0 : valid + 1;
                if (valid >= k) out.emplace_back(v, pos - k + 1);
            }
        } else {
            for (int i = hi - 1; i >= lo; --i, ++pos) {
                const unsigned char ch = (unsigned char)s[i];
                v = ((v << 2) | (3u - code2(ch))) & km;
                valid = (ch | 0x20) == 'n' ? 0 : valid + 1;
                if (valid >= k) out.emplace_back(v, pos - k + 1);
            }
        }
        std::sort(out.begin(), out.end());
    };
    collect(W.jp_a, a, c, false);
    collect(W.jp_b, c, b, true);
    W.fb_diag.clear();  // (i, j) pairs of shared k-mers
    size_t i = 0, j = 0;
    while (i < W.jp_a.size() && j < W.jp_b.size()) {
        if (W.jp_a[i].first < W.jp_b[j].first) { ++i; continue; }
        if (W.jp_a[i].first > W.jp_b[j].first) { ++j; continue; }
        size_t i2 = i, j2 = j;
        while (i2 < W.jp_a.size() && W.jp_a[i2].first == W.jp_a[i].first) ++i2;
        while (j2 < W.jp_b.size() && W.jp_b[j2].first == W.jp_b[j].first) ++j2;
        if (i2 - i <= 4 && j2 - j <= 4)
            for (size_t x = i; x < i2; ++x)
                for (size_t y = j; y < j2; ++y) { W.fb_diag.push_back(W.jp_a[x].second); W.fb_diag.push_back(W.jp_b[y].second); }
        i = i2; j = j2;
    }
    if (W.fb_diag.size() < 2 * (size_t)FOLD_DIAG_MIN) return false;
    W.fb_hist.assign((size_t)n_l + n_r + 1, 0);
    for (size_t q = 0; q < W.fb_diag.size(); q += 2) W.fb_hist[W.fb_diag[q + 1] - W.fb_diag[q] + n_l]++;
    int best = -1, bestd = 0;
    for (int d = 0; d <= n_l + n_r; ++d) {
        int cnt = 0;
        for (int e = -8; e <= 8; ++e) if (d + e >= 0 && d + e <= n_l + n_r) cnt += W.fb_hist[d + e];
        if (cnt > best) { best = cnt; bestd = d - n_l; }
    }
    for (size_t q = 0; q < W.fb_diag.size(); q += 2) {
        const int pi = W.fb_diag[q], pj = W.fb_diag[q + 1];
        if (std::abs(pj - pi - bestd) > 8) continue;
        ++non;
        min_i = std::min(min_i, pi); max_i = std::max(max_i, pi + k);
        min_j = std::min(min_j, pj); max_j = std::max(max_j, pj + k);
    }
    return non >= FOLD_DIAG_MIN;
}

// Fold-back decision for the insert s[a, b): arm Jaccard >= FOLD_J_LO, the mirrored diagonal reaching both arms'
// unit-side ends within reach_l / reach_r bp, and banded identity >= FOLD_NW_ID. Gets the centre, J and identity.
inline bool fold_decision(const char* s, int a, int b, Scratch& W, int& centre, float& J, float& ident, int reach_l = FOLD_REACH,
                          int reach_r = FOLD_REACH) {
    J = 0.f; ident = 0.f;
    centre = (a + b) / 2;
    if (b - a < 2 * FOLD_ARM_MIN) return false;
    float fr;
    foldback_centre(s, a, b, W, centre, fr, true);
    centre = std::min(std::max(centre, a + FOLD_ARM_MIN), b - FOLD_ARM_MIN);
    J = arm_jaccard(s, a, centre, b, W);
    if (J < FOLD_J_LO) return false;
    int mi, mj, xi, xj, non;
    if (!arm_mirror_reach(s, a, centre, b, W, mi, mj, xi, xj, non)) return false;
    if (mi > reach_l || mj > reach_r) return false;
    const int l0 = a + mi, l1 = std::min(centre, a + xi), r1 = b - mj, r0 = std::max(centre, b - xj);
    ident = seg_identity(s, l0, l1, r0, r1, W);
    return ident >= FOLD_NW_ID;
}


// Strand score of the k-mer at s[i] (0 for k-mers with N).
inline float strand_score_at(const Model& M, const char* s, int i) {
    uint32_t v = 0;
    for (int x = 0; x < STRAND_K; ++x) {
        const unsigned char ch = (unsigned char)s[i + x];
        if ((ch | 0x20) == 'n') return 0.f;
        v = (v << 2) | code2(ch);
    }
    return M.strand_lo[v];
}

// Strand-flip cut in s[lo, hi) (sense cDNA left, antisense right): cut = argmax of left sum - right sum. True when
// both arms pass min_margin, FLIP_EVIDENCE_MIN and the low-complexity guard; margin = min(left, -right mean) / strand_mu.
inline bool strand_flip_cut(const Model& M, const char* s, int lo, int hi, Scratch& W, int& cut, float& margin,
                            float min_margin = FLIP_MARGIN, int arm_min = FLIP_ARM_MIN) {
    cut = (lo + hi) / 2;
    margin = 0.f;
    if (!M.has_strand()) return false;
    const int n = hi - lo;
    if (n < 2 * arm_min + STRAND_K) return false;
    std::vector<float>& ps = W.flip_ps;
    ps.assign((size_t)n + 1, 0.f);
    for (int i = 0; i < n; ++i) {
        const float v = i + STRAND_K <= n ? strand_score_at(M, s, lo + i) : 0.f;
        ps[i + 1] = ps[i] + v;
    }
    const float total = ps[n];
    int best = -1;
    float bestv = -1e30f;
    for (int c = arm_min; c <= n - arm_min; ++c) {
        const float v = 2.f * ps[c] - total;
        if (v > bestv) { bestv = v; best = c; }
    }
    if (best < 0) return false;
    const float left = ps[best] / (float)best, right = (total - ps[best]) / (float)(n - best);
    const float t = min_margin * M.strand_mu;
    cut = lo + best;
    margin = std::min(left, -right) / std::max(M.strand_mu, 1e-6f);
    if (!(left >= t && right <= -t)) return false;
    if (std::min(ps[best], -(total - ps[best])) < (float)FLIP_EVIDENCE_MIN * M.strand_mu) return false;
    auto distinct_frac = [&](int a, int b) {
        std::vector<uint8_t>& seen = W.flip_seen;
        seen.assign((size_t)STRAND_N, 0);
        uint32_t v = 0;
        int valid = 0, pos = 0, d = 0;
        for (int i = a; i < b; ++i) {
            const unsigned char ch = (unsigned char)s[i];
            v = ((v << 2) | code2(ch)) & (uint32_t)(STRAND_N - 1);
            valid = (ch | 0x20) == 'n' ? 0 : valid + 1;
            if (valid >= STRAND_K) { ++pos; if (!seen[v]) { seen[v] = 1; ++d; } }
        }
        return pos ? (float)d / (float)pos : 0.f;
    };
    return distinct_frac(lo, cut) >= FLIP_DISTINCT_MIN && distinct_frac(cut, hi) >= FLIP_DISTINCT_MIN;
}

// D strand score of the k-mer at s[i] (0 for k-mers with N).
inline float strand_d_score_at(const Model& M, const char* s, int i) {
    uint32_t v = 0;
    for (int x = 0; x < STRAND_D_K; ++x) {
        const unsigned char ch = (unsigned char)s[i + x];
        if ((ch | 0x20) == 'n') return 0.f;
        v = (v << 2) | code2(ch);
    }
    return M.strand_d_lo[v];
}

// Low-complexity guard: distinct STRAND_K-mers / k-mer positions of s[a, b).
inline float distinct_kmer_frac(const char* s, int a, int b, Scratch& W) {
    std::vector<uint8_t>& seen = W.flip_seen;
    seen.assign((size_t)STRAND_N, 0);
    uint32_t v = 0;
    int valid = 0, pos = 0, d = 0;
    for (int i = a; i < b; ++i) {
        const unsigned char ch = (unsigned char)s[i];
        v = ((v << 2) | code2(ch)) & (uint32_t)(STRAND_N - 1);
        valid = (ch | 0x20) == 'n' ? 0 : valid + 1;
        if (valid >= STRAND_K) { ++pos; if (!seen[v]) { seen[v] = 1; ++d; } }
    }
    return pos ? (float)d / (float)pos : 0.f;
}

// D strand rule on the cDNA s[lo, hi) of a D geometry, in strand_d_mu units: T = total score; S / S' = weaker-arm
// evidence at the best sense-antisense / antisense-sense cut. DRULE_UNCLEAR if S' >= D_EVIDENCE_MIN; DRULE_TWO (cut
// set) if S >= D_EVIDENCE_MIN; DRULE_ONE_LEFT / _RIGHT if S <= 0 and |T| >= D_EVIDENCE_MIN; else DRULE_UNCLEAR.
enum : int { DRULE_UNCLEAR = 0, DRULE_ONE_LEFT = 1, DRULE_ONE_RIGHT = 2, DRULE_TWO = 3 };
inline int d_strand_rule(const Model& M, const char* s, int lo, int hi, Scratch& W, int& cut, float& evidence) {
    cut = -1;
    evidence = 0.f;
    if (!M.has_strand_d()) return DRULE_UNCLEAR;
    const int n = hi - lo;
    if (n < STRAND_D_K) return DRULE_UNCLEAR;
    std::vector<float>& ps = W.flip_ps;
    ps.assign((size_t)n + 1, 0.f);
    for (int i = 0; i < n; ++i) ps[i + 1] = ps[i] + (i + STRAND_D_K <= n ? strand_d_score_at(M, s, lo + i) : 0.f);
    const float mu = std::max(M.strand_d_mu, 1e-6f), total = ps[n];
    float S = -1e30f, Sp = -1e30f;
    int best = -1, bestp = -1;
    if (n >= 2 * FLIP_ARM_MIN + STRAND_D_K) {
        float bestv = -1e30f, bestpv = 1e30f;
        for (int c = FLIP_ARM_MIN; c <= n - FLIP_ARM_MIN; ++c) {
            const float v = 2.f * ps[c] - total;
            if (v > bestv) { bestv = v; best = c; }
            if (v < bestpv) { bestpv = v; bestp = c; }
        }
        S = std::min(ps[best], -(total - ps[best])) / mu;
        Sp = std::min(-ps[bestp], total - ps[bestp]) / mu;
    }
    if (Sp >= D_EVIDENCE_MIN) return DRULE_UNCLEAR;  // head to head, or a third arm
    if (S >= D_EVIDENCE_MIN && distinct_kmer_frac(s, lo, lo + best, W) >= FLIP_DISTINCT_MIN &&
        distinct_kmer_frac(s, lo + best, hi, W) >= FLIP_DISTINCT_MIN) {
        cut = lo + best;
        evidence = S;
        return DRULE_TWO;
    }
    const float T = total / mu;
    if (S <= 0.f && std::fabs(T) >= D_EVIDENCE_MIN && distinct_kmer_frac(s, lo, hi, W) >= FLIP_DISTINCT_MIN) {
        evidence = std::fabs(T);
        return T > 0.f ? DRULE_ONE_LEFT : DRULE_ONE_RIGHT;
    }
    return DRULE_UNCLEAR;
}

// Cut between the construct ending at slot sa (event A) and the next one starting at slot sb (event B); kind gets
// the CUT_* value. Never enters a known barcode block. Lx / Rx: offpath_edges results (INT32_MIN: none).
inline int junction_cut(const Model& M, const Event& A, int sa, const Event& B, int sb, uint8_t& kind, int Lx = INT32_MIN,
                        int Rx = INT32_MIN, uint8_t* boundary_flags = nullptr,
                        int* retained_lo = nullptr, int* retained_hi = nullptr) {
    const edge_span tail = template_edge_span(M, sa, false), head = template_edge_span(M, sb, true);
    const bool cl = (M.st_flags[sa] & SF_CLOSE) && tail.bounded && tail.hi == 0;
    const bool op = (M.st_flags[sb] & SF_OPEN) && head.bounded && head.hi == 0;
    const int lo = std::min(A.end, B.start), hi = std::max(A.end, B.start);
    const int bt = M.bc_tail[sa], bh = M.bc_head[sb];
    bool hl = tail.bounded, hr = head.bounded;
    int ll = A.end + tail.lo, lh = A.end + tail.hi;
    int rl = B.start - head.hi, rh = B.start - head.lo;
    // Existing inner-anchor/off-path barcode edges are usable when a terminal primer was missed.
    if (bt >= 0) { hl = true; ll = lh = A.end + bt; }
    else if (!hl && Lx != INT32_MIN) { hl = true; ll = lh = Lx; }
    if (bh >= 0) { hr = true; rl = rh = B.start - bh; }
    else if (!hr && Rx != INT32_MIN) { hr = true; rl = rh = Rx; }
    uint8_t flags = (hl ? BOUNDARY_LEFT_BOUNDED : 0) | (hr ? BOUNDARY_RIGHT_BOUNDED : 0) |
                    (tail.bounded ? 0 : BOUNDARY_LEFT_INSERT) | (head.bounded ? 0 : BOUNDARY_RIGHT_INSERT);
    int cut, wlo, whi;
    const int lp = cl ? 0 : BC_WIN_PAD, rp = op ? 0 : BC_WIN_PAD;
    switch ((hl ? 1 : 0) | (hr ? 2 : 0)) {
    case 3:  // Both edges: midpoint in the interval bounded by physical layout edges.
        cut = (lh + rl) / 2;
        kind = cl && op ? CUT_BOTH_ADAPTERS : CUT_GEOMETRY;
        wlo = std::min(cut, rl - rp);
        whi = std::max(cut, lh + lp);
        // Honor the decoder's microhomology allowance at physical edges too.
        // Retained windows keep both blocks intact; larger conflicts remain unresolved.
        if (rh - ll < DEL_GMIN) flags |= BOUNDARY_UNRESOLVED;
        break;
    case 1:  // Left edge only: retain the unknown sequence on its right.
        cut = (ll + lh) / 2;
        kind = cl ? CUT_ONE_ADAPTER : CUT_GEOMETRY;
        wlo = ll - lp; whi = lh + lp;
        break;
    case 2:  // Right edge only: retain the unknown sequence on its left.
        cut = (rl + rh) / 2;
        kind = op ? CUT_ONE_ADAPTER : CUT_GEOMETRY;
        wlo = rl - rp; whi = rh + rp;
        break;
    default: // Neither edge: a midpoint is bookkeeping only, not an assigned boundary.
        cut = (A.end + B.start) / 2;
        kind = CUT_MIDPOINT;
        flags |= BOUNDARY_UNRESOLVED;
        wlo = lo; whi = hi;
        break;
    }
    if (cut < lo || cut > hi) flags |= BOUNDARY_UNRESOLVED;
    if (flags & BOUNDARY_UNRESOLVED) { wlo = std::min(wlo, lo); whi = std::max(whi, hi); }
    cut = std::min(std::max(cut, lo), hi);
    if (boundary_flags) *boundary_flags = flags;
    if (retained_lo) *retained_lo = std::min(wlo, cut);
    if (retained_hi) *retained_hi = std::max(whi, cut);
    return cut;
}

// Barcode-block edges Lx / Rx for a junction whose facing slots give none, from inner barcode-side events the decode
// left off the path, consistent with the other side's edge. lanc / ranc: the edge came from an anchor.
inline void offpath_edges(const Model& M, const Event* ev, int nev, const Scratch& W, int sa, const Event& A, int sb, const Event& B,
                          int& Lx, int& Rx, bool* lanc = nullptr, bool* ranc = nullptr) {
    Lx = Rx = INT32_MIN;
    if (lanc) *lanc = false;
    if (ranc) *ranc = false;
    const bool left_edge = (M.st_flags[sa] & SF_CLOSE) && at_template_edge(M, sa, false);
    const bool right_edge = (M.st_flags[sb] & SF_OPEN) && at_template_edge(M, sb, true);
    const bool need_l = M.bc_tail[sa] < 0 && !left_edge, need_r = M.bc_head[sb] < 0 && !right_edge;
    if (!need_l && !need_r) return;
    auto inner_of = [&](int t, bool head, int& blk) -> int {  // the template's barcode-side inner slot type (-1: none)
        const int s0 = M.open_state[t], no = (int)M.tpl[t].obs.size();
        for (int o = 0; o < no; ++o) {
            const int b = head ? M.bc_head[s0 + o] : M.bc_tail[s0 + o];
            if (b >= 0) { blk = b; return M.st_type[s0 + o]; }
        }
        return -1;
    };
    auto on_path = [&](int i) {
        for (int x = 0; x < (int)W.path_ev.size(); ++x)
            if (W.path_ev[x] == i) return true;
        return false;
    };
    auto usable = [&](const Event& e, int a) { return e.type == a && (a < M.A ? e.cls == CL_S4 : e.cls == CL_POLY); };
    const int glo = std::min(A.end, B.start) - 4, ghi = std::max(A.end, B.start) + 4;
    const int ta = M.tail_fix[sa], hb = M.head_fix[sb];
    // the other side's edge: an observed closing (opening) primer bounds its own barcode block at its start (end)
    const int Lref = left_edge ? A.start : M.bc_tail[sa] >= 0 ? A.end + M.bc_tail[sa] : ta >= 0 ? A.end + ta : INT32_MIN;
    const int Rref = right_edge ? B.end : M.bc_head[sb] >= 0 ? B.start - M.bc_head[sb] : hb >= 0 ? B.start - hb : INT32_MIN;
    if (need_l) {
        int blk = 0;
        const int a = inner_of(M.st_tpl[sa], false, blk);
        const int dmax = a >= M.A ? RESCUE_GMAX_OFF : RESCUE_GMAX, dmin = a >= M.A ? -16 : -6;
        if (a >= 0 && Rref != INT32_MIN)
            for (int i = 0; i < nev; ++i) {
                const Event& e = ev[i];
                if (e.start > ghi) break;
                if (!usable(e, a) || e.start < glo - 64 || e.end > ghi) continue;
                const int d = Rref - (e.end + blk);
                if (d >= dmin && d <= dmax && !on_path(i)) Lx = std::max(Lx, e.end + blk);
            }
        if (lanc) *lanc = Lx != INT32_MIN && a < M.A;
    }
    if (need_r) {
        int blk = 0;
        const int a = inner_of(M.st_tpl[sb], true, blk);
        const int dmax = a >= M.A ? RESCUE_GMAX_OFF : RESCUE_GMAX, dmin = a >= M.A ? -16 : -6;
        const int Lr = Lx != INT32_MIN ? Lx : Lref;
        if (a >= 0 && Lr != INT32_MIN)
            for (int i = 0; i < nev; ++i) {
                const Event& e = ev[i];
                if (e.start > ghi) break;
                if (!usable(e, a) || e.start < glo || e.end > ghi + 64) continue;
                const int d = (e.start - blk) - Lr;
                if (d >= dmin && d <= dmax && !on_path(i)) { Rx = e.start - blk; break; }
            }
        if (ranc) *ranc = Rx != INT32_MIN && a < M.A;
    }
}

// An exact 11-mer seed of the event's own anchor overlaps it (stage-0 seed clusters).
inline bool exact_seeded(const Scratch& W, const Event& e) {
    for (const seed_cluster& c : W.sc)
        if (c.a == e.type && c.nx >= 1 && c.s < e.end && e.start < c.e) return true;
    return false;
}
// True when construct path steps [x0, x1) hold cDNA-side evidence (a slot outside Model::bc_unit).
inline bool cdna_side_anchor(const Model& M, const Scratch& W, const Event* ev, int x0, int x1) {
    for (int x = x0; x < x1; ++x) {
        if (M.bc_unit[W.path_st[x]]) continue;
        const Event& e = ev[W.path_ev[x]];
        if (e.cls == CL_POLY || e.cls == CL_S4 || (e.cls == CL_SP && (e.fl & (EVF_TRUNC3 | EVF_TRUNC5)))) return true;
    }
    return false;
}
// unused
inline bool fold_mirrored(const Model& M, const Event* ev, int nev, int t, const int* ostart, const int* oend, const int* oev, int nobs,
                          int fold, int lo, int hi, bool first, bool last, int L) {
    if (t >= (int)M.tpl_mirror.size() || M.tpl_mirror[t].empty()) return false;
    const Template& T = M.tpl[t];
    const std::vector<int>& mir = M.tpl_mirror[t];
    auto seen = [&](int o) { return ostart[o] >= 0 && oev[o] >= 0; };
    int c2 = INT32_MIN, best_el = -1;
    for (int o = 0; o < nobs; ++o) {
        const int mo = mir[o];
        if (mo < 0 || mo >= nobs || T.pins[T.obs[o]] != 0 || !seen(o) || !seen(mo)) continue;
        if (T.obs[o] > best_el) { best_el = T.obs[o]; c2 = oend[o] + ostart[mo]; }
    }
    if (c2 == INT32_MIN) return false;
    int ddmin = INT32_MAX, ddmax = INT32_MIN, mind = INT32_MAX, arm = INT32_MAX;
    for (int o = 0; o < nobs; ++o) {
        const int mo = mir[o];
        if (mo < 0 || mo >= nobs || T.pins[T.obs[o]] != 0) continue;
        const bool sh = seen(o), st = seen(mo);
        if (sh && st) {
            const int d1 = fold - oend[o], d2 = ostart[mo] - fold, dd = d1 - d2;
            ddmin = std::min(ddmin, dd); ddmax = std::max(ddmax, dd);
            mind = std::min(mind, std::abs(dd));
            arm = std::min(arm, std::min(d1, d2));
        } else if (sh) {
            const int a = T.el[T.obs[mo]].anc, s0 = c2 - oend[o], e0 = s0 + M.anch[a].m;
            if (e0 <= (last ? L - FOLD_ROOM : hi)) return false;
            for (int i = 0; i < nev; ++i)
                if (ev[i].type == a && ev[i].cls <= CL_S6 && ev[i].start >= fold && ev[i].start < (last ? L : hi)) return false;
        } else if (st) {
            const int a = T.el[T.obs[o]].anc, e0 = c2 - ostart[mo], s0 = e0 - M.anch[a].m;
            if (s0 >= (first ? FOLD_ROOM : lo)) return false;
            for (int i = 0; i < nev; ++i)
                if (ev[i].type == a && ev[i].cls <= CL_S6 && ev[i].end <= fold && ev[i].end > (first ? 0 : lo)) return false;
        }
    }
    if (ddmax - ddmin > FOLD_SPREAD) return false;
    return mind <= std::max(12, arm / 8);
}

// Check region of an observed slot (anchor or poly) of a construct: FULL / PARTIAL / TRUNCATED from its evidence class.
inline void seen_check(Result& R, const Model& M, const TElem& elm, char strand, int cidx, const char* seq, int L, const Event& e, float conf) {
    const int ws = e.start - M.opt.seen_pad, we = e.end + M.opt.seen_pad;
    if (elm.kind == EK_POLY) { add_check(R, M, elm, strand, cidx, Status::FULL, ws, we, -1, -1, -1, conf, L); return; }
    const Anchor& x = M.anch[elm.anc];
    if (e.cls == CL_S4) add_check(R, M, elm, strand, cidx, Status::FULL, ws, we, 0, x.m, e.ed, conf, L);
    else if (e.cls == CL_S6 || e.cls == CL_S8) {
        int fws, fwe, flo, fhi, fed;
        if (degraded_as_fragment(M, elm.anc, seq, L, e.start, e.end, e.ed, fws, fwe, flo, fhi, fed))
            add_check(R, M, elm, strand, cidx, Status::PARTIAL, fws, fwe, flo, fhi, fed, conf, L);
        else
            add_check(R, M, elm, strand, cidx, Status::PARTIAL, ws, we, 0, x.m, e.ed, conf, L);
    }
    else if (e.fl & (EVF_TRUNC3 | EVF_TRUNC5)) add_check(R, M, elm, strand, cidx, Status::TRUNCATED_AT_READ_END, ws, (e.fl & EVF_TRUNC3) ? L : we, e.rlo, e.rhi, e.ed, conf, L);
    else add_check(R, M, elm, strand, cidx, Status::PARTIAL, ws, we, e.rlo, e.rhi, e.ed, conf, L);
}

// Check regions of a split half [lo, hi), described by main template tmain: seen checks for observed slots, and
// MISSING_EXPECTED windows only for slots tied to an observed one without an insert.
inline void emit_half_checks(Result& R, const Model& M, const char* seq, int L, int tmain, int cidx, int lo, int hi, const Event* ev,
                             const Template& T, const int* ostart, const int* oend, const int* oev, int nobs, float conf) {
    const Template& H = M.tpl[tmain];
    const int hn = std::min((int)H.obs.size(), 64);
    int hs[64], he[64], hv[64];
    for (int o = 0; o < hn; ++o) {
        hs[o] = he[o] = hv[o] = -1;
        const TElem& x = H.el[H.obs[o]];
        for (int q = 0; q < nobs; ++q) {
            if (oev[q] < 0 || ostart[q] < lo || oend[q] > hi) continue;
            const TElem& y = T.el[T.obs[q]];
            if (y.kind != x.kind || y.anc != x.anc) continue;
            hs[o] = ostart[q]; he[o] = oend[q]; hv[o] = oev[q];
            break;
        }
    }
    for (int o = 0; o < hn; ++o) {
        const TElem& elm = H.el[H.obs[o]];
        if (hv[o] >= 0) { seen_check(R, M, elm, H.strand, cidx, seq, L, ev[hv[o]], conf); continue; }
        // expected position from a tight neighbour (no insert between), as in expected_window
        const int eo = H.obs[o];
        const int len = (int)std::lround(elm.nominal);
        int es = INT32_MIN;
        double sp = 0;
        for (int ob = o - 1; ob >= 0 && es == INT32_MIN; --ob) {
            const int e1 = H.obs[ob];
            if (H.pins[eo] - H.pins[e1 + 1] != 0) break;
            if (hs[ob] >= 0) { es = he[ob] + (int)std::lround(H.pnom[eo] - H.pnom[e1 + 1]); sp = H.phr[eo] - H.phr[e1 + 1]; }
        }
        for (int oa = o + 1; oa < hn && es == INT32_MIN; ++oa) {
            const int e2 = H.obs[oa];
            if (H.pins[e2] - H.pins[eo + 1] != 0) break;
            if (hs[oa] >= 0) { es = hs[oa] - (int)std::lround(H.pnom[e2] - H.pnom[eo + 1]) - len; sp = H.phr[e2] - H.phr[eo + 1]; }
        }
        if (es == INT32_MIN) continue;
        const int spread = (int)std::lround(sp);
        const int ws = std::max(lo, es - spread - M.opt.check_pad), we = std::min(hi, es + len + spread + M.opt.check_pad);
        if (we <= ws) continue;
        missing_slot_check(R, M, elm, H.strand, cidx, seq, L, ws, we, es, 0.5f * conf, 0.5f * conf);
    }
}

// Builds the Result (cuts, segments, check regions, strand call, abstain flag) from the decoded path.
inline void derive_result(const Model& M, const char* seq, const Event* ev, int nev, int L, Scratch& W, int k, bool fb_done, Result& R) {
    R.k = 0;
    R.boundary_unresolved_seen = false;
    R.cuts.clear(); R.segs.clear(); R.checks.clear(); R.cut_post.clear(); R.cut_lo.clear(); R.cut_hi.clear(); R.cut_kind.clear(); R.cut_flags.clear(); R.cut_aux.clear();
    R.p_single = 0.f;
    if (k == 0) { R.strand_call = '?'; return; }
    if (fb_done) R.p_single = (float)std::exp(std::min(0.f, W.l_z1 - W.l_z));
    else R.p_single = k == 1 ? 1.f : 0.f;
    const int np = (int)W.path_ev.size();
    std::vector<int>& cstart = W.cidx;
    cstart.clear();
    for (int x = 0; x < np; ++x)
        if (W.path_cross[x]) cstart.push_back(x);
    const int nc = (int)cstart.size();
    cstart.push_back(np);
    int prev_cut = 0;
    std::vector<int>& cuts = R.cuts;
    std::vector<uint8_t>& ckind = W.ckind;
    std::vector<int>& ca = W.ca;
    std::vector<int>& cb = W.cb;
    ckind.clear(); ca.clear(); cb.clear();
    W.jlo.clear(); W.jhi.clear(); W.jflags.clear(); W.cweak.clear(); W.segj.clear(); W.cgap.clear(); W.caux.clear(); W.saux.clear(); W.skind.clear(); W.scut.clear();
    for (int c = 0; c + 1 < nc; ++c) {
        const int xa = cstart[c + 1] - 1, xb = cstart[c + 1];
        const int sa = W.path_st[xa], sb = W.path_st[xb];
        const Event& A = ev[W.path_ev[xa]];
        const Event& B = ev[W.path_ev[xb]];
        int Lx = INT32_MIN, Rx = INT32_MIN;
        bool lanc = false, ranc = false;
        offpath_edges(M, ev, nev, W, sa, A, sb, B, Lx, Rx, &lanc, &ranc);
        uint8_t kind, boundary_flags;
        int boundary_lo, boundary_hi;
        float aux = 0.f;
        int cut = junction_cut(M, A, sa, B, sb, kind, Lx, Rx, &boundary_flags, &boundary_lo, &boundary_hi);
        const int bt = M.bc_tail[sa], bh = M.bc_head[sb];
        // a junction made only by a pairing of partner anchors needs a strong pairing or a cDNA-side anchor in both
        // constructs; otherwise its posterior is capped at WEAK_POST
        uint8_t weak = 0;
        int gap_pair = INT32_MIN;
        {
            const bool il = bt >= 0 && A.type < M.A, ir = bh >= 0 && B.type < M.A;
            auto primer_cert = [&](const Event& e) { return e.type < M.A && e.cls == CL_S4 && !M.anch[e.type].shrt; };
            if ((il && ir) || (il && !primer_cert(B)) || (ir && !primer_cert(A))) {
                const int g = (B.start - (ir ? bh : 0)) - (A.end + (il ? bt : 0));
                auto strong = [&](const Event& e) {
                    return e.type < M.A && e.cls == CL_S4 && (!M.anch[e.type].shrt || (e.ed >= 0 && e.ed <= 1) || exact_seeded(W, e));
                };
                const bool pair = g <= DEL_GMAX_CERT && strong(A) && strong(B);
                const bool ok = pair || (cdna_side_anchor(M, W, ev, cstart[c], cstart[c + 1]) && cdna_side_anchor(M, W, ev, cstart[c + 1], cstart[c + 2]));
                weak = (uint8_t)((ok ? 0 : 1) | (pair ? 2 : 0) | (il ? 4 : 0) | (ir ? 8 : 0));
                gap_pair = g;
            }
        }
        W.cweak.push_back(weak);
        W.cgap.push_back(gap_pair);
        if (kind == 3) {
            // no located junction: fold-back > strand flip (F then R) > duration-MAP midpoint
            const int ja = std::min(A.end, B.start), jb = std::max(A.end, B.start);
            int ctr;
            float fj, fid;
            const int rl = (A.type < M.A && (M.inner_role[A.type] & 3)) ? FOLD_REACH : FOLD_REACH_OPEN;
            const int rr = (B.type < M.A && (M.inner_role[B.type] & 3)) ? FOLD_REACH : FOLD_REACH_OPEN;
            if (M.opt.foldback_split && fold_decision(seq, ja, jb, W, ctr, fj, fid, rl, rr)) {
                cut = std::min(std::max(ctr, ja), jb);
                kind = CUT_FOLDBACK;
                aux = fj;
            } else if (M.has_strand() && M.tpl[M.st_tpl[sa]].strand == 'F' && M.tpl[M.st_tpl[sb]].strand == 'R') {
                int fc;
                float mg;
                if (strand_flip_cut(M, seq, ja, jb, W, fc, mg)) { cut = fc; kind = CUT_STRAND_FLIP; aux = mg; }
            }
            if (kind == CUT_FOLDBACK || kind == CUT_STRAND_FLIP) {
                boundary_flags = (uint8_t)((boundary_flags & ~BOUNDARY_UNRESOLVED) | BOUNDARY_SEQUENCE_RESOLVED);
                boundary_lo = boundary_hi = cut;
            }
        }
        cut = std::max(cut, prev_cut + 1);
        // child windows: an anchor-derived barcode block stays whole, but never past the facing construct's events
        const int Le = (bt >= 0 && A.type < M.A) ? A.end + bt : lanc ? Lx : INT32_MIN;
        const int Re = (bh >= 0 && B.type < M.A) ? B.start - bh : ranc ? Rx : INT32_MIN;
        int whi = std::max(cut, boundary_hi), wlo = std::min(cut, boundary_lo);
        if (Le != INT32_MIN) whi = std::max(whi, std::min(std::max(cut, Le + BC_WIN_PAD), std::max(cut, B.start)));
        if (Re != INT32_MIN) wlo = std::min(wlo, std::max(std::min(cut, Re - BC_WIN_PAD), std::min(cut, A.end)));
        whi = std::min(whi, L); wlo = std::max(wlo, 0);
        W.jlo.push_back(wlo);
        W.jhi.push_back(whi);
        W.jflags.push_back(boundary_flags);
        R.boundary_unresolved_seen = R.boundary_unresolved_seen || (boundary_flags & BOUNDARY_UNRESOLVED);
        cuts.push_back(cut);
        ckind.push_back(kind);
        W.caux.push_back(aux);
        ca.push_back(W.path_ev[xa]);
        cb.push_back(W.path_ev[xb]);
        prev_cut = cut;
    }
    std::vector<float>& cpost = W.cpost;
    cpost.clear();
    for (int c = 0; c + 1 < nc; ++c) {
        float p = 1.f;
        if (fb_done && M.opt.posteriors && k >= 2) {
            p = 0.f;
            for (int x = ca[c]; x < cb[c]; ++x) p = std::max(p, W.jpost[x]);
        }
        if (W.cweak[c] & 1) p = std::min(p, WEAK_POST);
        cpost.push_back(p);
    }
    std::vector<int>& bcut = W.bcut;  // cuts between path constructs (segments use these; fold-back cuts are added later)
    bcut.assign(cuts.begin(), cuts.end());
    std::vector<float>& spost = W.spost;  // posterior of the cut to the right of each output segment
    spost.clear();
    bool any_f = false, any_r = false, any_o = false, any_u = false;
    int ostart[64], oend[64], oev[64];
    for (int c = 0; c < nc; ++c) {
        const int x0 = cstart[c], x1 = cstart[c + 1];
        const int t = M.st_tpl[W.path_st[x0]];
        const Template& T = M.tpl[t];
        const int nobs = std::min((int)T.obs.size(), 64);
        for (int o = 0; o < nobs; ++o) { ostart[o] = oend[o] = -1; oev[o] = -1; }
        for (int x = x0; x < x1; ++x) {
            const int o = M.st_obs[W.path_st[x]];
            const Event& e = ev[W.path_ev[x]];
            if (o < nobs) { ostart[o] = e.start; oend[o] = e.end; oev[o] = W.path_ev[x]; }
        }
        const int sfirst = W.path_st[x0], slast = W.path_st[x1 - 1];
        const Event& ef = ev[W.path_ev[x0]];
        const Event& el = ev[W.path_ev[x1 - 1]];
        const bool has_open = (M.st_flags[sfirst] & SF_OPEN) != 0;
        const bool has_close = (M.st_flags[slast] & SF_CLOSE) != 0;
        Segment sg;
        sg.strand = T.strand;
        sg.flags = (uint16_t)((has_open ? 0 : SEG_PARTIAL_LEFT) | (has_close ? 0 : SEG_PARTIAL_RIGHT) | (T.art && !T.dtpl ? SEG_ARTIFACT : 0));
        int cs_slack = 0, ce_slack = 0;
        if (c > 0) sg.start = bcut[c - 1];
        else {
            // read-end bound at the barcode-block edge when the barcode-adjacent primer was not observed
            sg.start = retained_template_start(M, sfirst, ef.start);
            if (!template_edge_span(M, sfirst, true).bounded) cs_slack = M.zone5;
        }
        if (c + 1 < nc) {
            sg.end = bcut[c];
            if (ckind[c] >= 2) sg.flags |= SEG_UNANCHORED_R;
            if (ckind[c] == 3) sg.flags |= SEG_LOWCONF;
            if (ckind[c] == 4) sg.flags |= SEG_SAME_MOLECULE_R;
        } else {
            sg.end = retained_template_end(M, slast, el.end, L);
            if (!template_edge_span(M, slast, false).bounded) ce_slack = M.zone3;
        }
        sg.end = std::max(sg.end, sg.start + 1);
        float conf = R.p_single;
        if (nc > 1) {
            conf = 1.f;
            if (c > 0) conf = std::min(conf, cpost[c - 1]);
            if (c + 1 < nc) conf = std::min(conf, cpost[c]);
        }
        sg.conf = conf;
        if (conf < M.opt.tau_junction) sg.flags |= SEG_LOWCONF;
        {
            const int E = (int)T.el.size();
            int ins_el = -1;
            for (int x = 0; x < E; ++x) if (T.el[x].kind == EK_INSERT) { ins_el = x; break; }
            if (ins_el >= 0) {
                bool seen5 = false, seen3 = false;
                for (int o = 0; o < nobs; ++o) {
                    if (ostart[o] < 0) continue;
                    if (T.obs[o] < ins_el) seen5 = true; else seen3 = true;
                }
                const double l5 = seen5 ? T.pnom[ins_el] : 0, l3 = seen3 ? T.pnom[E] - T.pnom[ins_el + 1] : 0;
                if (seen5 && seen3 && (sg.end - sg.start) - l5 - l3 < 30) sg.flags |= SEG_EMPTY;
            }
        }
        // a single-primer end without its barcode-side partner anchor has no strand evidence
        if (T.art && !T.dtpl && (T.strand == 'F' || T.strand == 'R')) {
            bool partner = false;
            for (int o = 0; o < nobs && !partner; ++o) {
                const int s2 = M.open_state[t] + o;
                partner = ostart[o] >= 0 && (M.bc_head[s2] >= 0 || M.bc_tail[s2] >= 0);
            }
            if (!partner) sg.strand = '?';
        }
        // fold-back -> two constructs; otherwise, for a D construct, the D strand rule
        int fold_cut = -1;
        float fold_aux = 0.f, flip_conf = -1.f;  // flip_conf >= 0: cut posterior of a D strand rule split
        uint8_t in_kind = CUT_FOLDBACK;
        int d_keep = DRULE_UNCLEAR;
        const int cut_right = sg.end;    // stays at the junction even when a kept-left D construct trims sg.end
        if (M.tpl_fold[t]) {
            int a = sg.start, b = sg.end;
            int rl = FOLD_REACH_OPEN, rr = FOLD_REACH_OPEN;  // an arm bounded by an observed inner anchor is closed
            for (int o = 0; o < nobs; ++o) {
                if (ostart[o] < 0) continue;
                const int eo = T.obs[o];
                const bool inner = T.el[eo].kind == EK_ANCHOR && (M.inner_role[T.el[eo].anc] & 3);
                if (T.pins[eo] == 0) { if (oend[o] >= a) { a = oend[o]; rl = inner ? FOLD_REACH : FOLD_REACH_OPEN; } }
                else if (ostart[o] <= b) { b = ostart[o]; rr = inner ? FOLD_REACH : FOLD_REACH_OPEN; }
            }
            int ctr;
            float fj, fid;
            if (M.opt.foldback_split && b - a >= 2 * FOLD_ARM_MIN && fold_decision(seq, a, b, W, ctr, fj, fid, rl, rr)) { fold_cut = ctr; fold_aux = fj; }
            else if (T.dtpl) {
                // D strand rule on the cDNA between the units' inner anchors (layout offset when one was not observed)
                int lo = a, hi = b;
                int ins_el = -1;
                for (int x = 0; x < (int)T.el.size(); ++x) if (T.el[x].kind == EK_INSERT) { ins_el = x; break; }
                if (ins_el >= 0) {
                    int o_l = -1, o_r = -1;
                    for (int o = 0; o < nobs; ++o) {
                        if (ostart[o] < 0) continue;
                        if (T.obs[o] < ins_el) o_l = o; else if (o_r < 0) o_r = o;
                    }
                    if (o_l >= 0) lo = oend[o_l] + (int)std::lround(T.pnom[ins_el] - T.pnom[T.obs[o_l] + 1]);
                    if (o_r >= 0) hi = ostart[o_r] - (int)std::lround(T.pnom[T.obs[o_r]] - T.pnom[ins_el + 1]);
                    lo = std::max(a, std::min(lo, L));
                    hi = std::min(b, std::max(hi, 0));
                }
                int dc = -1;
                float dev = 0.f;
                const int rule = hi > lo ? d_strand_rule(M, seq, lo, hi, W, dc, dev) : (int)DRULE_UNCLEAR;
                if (rule == DRULE_TWO) {
                    fold_cut = dc; fold_aux = dev; in_kind = CUT_STRAND_FLIP;
                    flip_conf = 1.f;
                } else if (rule == DRULE_ONE_LEFT || rule == DRULE_ONE_RIGHT) {
                    d_keep = rule;
                    if (rule == DRULE_ONE_LEFT) { sg.strand = 'F'; sg.end = std::max(hi, sg.start + 1); sg.flags |= SEG_PARTIAL_RIGHT; }
                    else { sg.strand = 'R'; sg.start = std::min(lo, sg.end - 1); sg.flags |= SEG_PARTIAL_LEFT; }
                    sg.flags = (uint16_t)((sg.flags & ~SEG_LOWCONF) | SEG_D_KEPT);
                } else {
                    sg.flags |= SEG_DOUBLE_BC;
                    // insert longer than one cDNA explains (two inserts more likely): uncertain
                    const int X = b - a;
                    if (X > M.ins1_p99 && gap_score_ins(M, 2, X) > gap_score_ins(M, 1, X)) sg.flags |= SEG_LOWCONF;
                }
            }
        }
        const int ci = (int)R.segs.size();
        if (fold_cut > sg.start && fold_cut < sg.end) {
            // the halves are plain F and R constructs; the split's confidence is the construct's own
            Segment s1 = sg, s2 = sg;
            const uint16_t dsplit = (T.dtpl && in_kind == CUT_STRAND_FLIP) ? SEG_D_SPLIT : 0;
            s1.end = fold_cut; s1.strand = 'F';
            s1.flags = (uint16_t)((sg.flags & ~(SEG_PARTIAL_RIGHT | SEG_SAME_MOLECULE_R | SEG_UNANCHORED_R | SEG_ARTIFACT | SEG_DOUBLE_BC)) | SEG_PARTIAL_RIGHT |
                                  (in_kind == CUT_FOLDBACK ? SEG_SAME_MOLECULE_R : SEG_UNANCHORED_R) | dsplit);
            s2.start = fold_cut; s2.strand = 'R'; s2.flags = (uint16_t)((sg.flags & ~(SEG_PARTIAL_LEFT | SEG_ARTIFACT | SEG_DOUBLE_BC)) | SEG_PARTIAL_LEFT | dsplit);
            // a D strand rule split: the halves take the rule's decision as confidence
            const float fpost = flip_conf >= 0.f ? flip_conf : conf;
            if (flip_conf >= 0.f) {
                s1.conf = s2.conf = fpost;
                s1.flags = (uint16_t)(s1.flags & ~SEG_LOWCONF);
                s2.flags = (uint16_t)(s2.flags & ~SEG_LOWCONF);
            }
            R.segs.push_back(s1);
            R.segs.push_back(s2);
            spost.push_back(fpost);
            spost.push_back(c + 1 < nc ? cpost[c] : 1.f);
            W.segj.push_back(-1);
            W.segj.push_back(c + 1 < nc ? c : -1);
            W.saux.push_back(fold_aux);
            W.saux.push_back(c + 1 < nc ? W.caux[c] : 0.f);
            W.skind.push_back(in_kind);
            W.skind.push_back(0);
            W.scut.push_back(fold_cut);
            W.scut.push_back(cut_right);
            any_f = true; any_r = true;
        } else {
            R.segs.push_back(sg);
            spost.push_back(c + 1 < nc ? cpost[c] : 1.f);
            W.segj.push_back(c + 1 < nc ? c : -1);
            W.saux.push_back(c + 1 < nc ? W.caux[c] : 0.f);
            W.skind.push_back(0);
            W.scut.push_back(cut_right);
            if (sg.strand == 'F') any_f = true;
            else if (sg.strand == 'R') any_r = true;
            else any_o = true;
            if (sg.strand == '?') any_u = true;
        }
        if (!M.opt.check_regions) continue;
        if (fold_cut > sg.start && fold_cut < sg.end) {
            // fold-back halves: the main templates describe them (F left of the fold, R right of it)
            int t_f = -1, t_r = -1;
            for (int tm : M.main_tpl) (M.tpl[tm].strand == 'F' ? t_f : t_r) = tm;
            if (t_f >= 0) emit_half_checks(R, M, seq, L, t_f, ci, sg.start, fold_cut, ev, T, ostart, oend, oev, nobs, conf);
            if (t_r >= 0) emit_half_checks(R, M, seq, L, t_r, ci + 1, fold_cut, sg.end, ev, T, ostart, oend, oev, nobs, conf);
            continue;
        }
        if (d_keep != DRULE_UNCLEAR) {
            // kept D unit: the main template of its strand describes the kept span
            int tm = -1;
            for (int x : M.main_tpl) if (M.tpl[x].strand == sg.strand) tm = x;
            if (tm >= 0) emit_half_checks(R, M, seq, L, tm, ci, sg.start, sg.end, ev, T, ostart, oend, oev, nobs, conf);
            continue;
        }
        const int left_kind = c > 0 ? ckind[c - 1] : -1, right_kind = c + 1 < nc ? ckind[c] : -1;
        // junction-facing slots next to a midpoint cut get the whole inter-construct interval; next to a strand-flip
        // cut, 40 bp on each side of it
        for (int o = 0; o < nobs; ++o) {
            const TElem& elm = T.el[T.obs[o]];
            if (oev[o] >= 0) { seen_check(R, M, elm, T.strand, ci, seq, L, ev[oev[o]], conf); continue; }
            int ws, we, spr, es;
            expected_window(M, t, o, ostart, oend, sg.start, sg.end, M.opt.check_pad, ws, we, spr, &es, cs_slack, ce_slack);
            if (o == nobs - 1 && c + 1 < nc) {
                if (right_kind == CUT_MIDPOINT) { ws = std::min(ws, ev[ca[c]].end - M.opt.check_pad); we = std::max(we, ev[cb[c]].start + M.opt.check_pad); }
                else if (right_kind == CUT_STRAND_FLIP) { ws = std::min(ws, bcut[c] - 40); we = std::max(we, bcut[c] + 40); }
            }
            if (o == 0 && c > 0) {
                if (left_kind == CUT_MIDPOINT) { ws = std::min(ws, ev[ca[c - 1]].end - M.opt.check_pad); we = std::max(we, ev[cb[c - 1]].start + M.opt.check_pad); }
                else if (left_kind == CUT_STRAND_FLIP) { ws = std::min(ws, bcut[c - 1] - 40); we = std::max(we, bcut[c - 1] + 40); }
            }
            missing_slot_check(R, M, elm, T.strand, ci, seq, L, ws, we, es, 0.5f * conf, 0.5f * conf);
        }
    }
    R.k = (int)R.segs.size();
    R.cuts.clear();
    R.cut_post.clear();
    R.cut_kind.clear();
    R.cut_flags.clear();
    for (int i = 0; i + 1 < R.k; ++i) {
        // the cut is the construct's junction (W.scut), which equals the segment end except for a kept-left D construct
        const int cut = i < (int)W.scut.size() ? W.scut[i] : R.segs[i].end, j = W.segj[i];
        R.cuts.push_back(cut);
        R.cut_post.push_back(spost[i]);
        R.cut_lo.push_back(j >= 0 ? std::min(W.jlo[j], cut) : cut);
        R.cut_hi.push_back(j >= 0 ? std::max(W.jhi[j], cut) : cut);
        R.cut_kind.push_back(j >= 0 ? ckind[j] : (i < (int)W.skind.size() && W.skind[i] ? W.skind[i] : (uint8_t)CUT_FOLDBACK));
        R.cut_flags.push_back(j >= 0 ? W.jflags[j] : (uint8_t)BOUNDARY_SEQUENCE_RESOLVED);
        R.cut_aux.push_back(i < (int)W.saux.size() ? W.saux[i] : 0.f);
    }
    if (R.k > M.opt.k_cap) {
        R.flags |= RES_TOO_MANY;
        const int kc = std::max(1, M.opt.k_cap);
        R.segs[kc - 1].end = R.segs.back().end;
        R.segs.resize(kc);
        R.cuts.resize(kc - 1);
        R.cut_post.resize(kc - 1);
        R.cut_lo.resize(kc - 1);
        R.cut_hi.resize(kc - 1);
        R.cut_kind.resize(kc - 1);
        R.cut_flags.resize(kc - 1);
        R.cut_aux.resize(kc - 1);
        R.k = kc;
        auto it = std::remove_if(R.checks.begin(), R.checks.end(), [&](const check_region& c) { return c.construct >= kc; });
        R.checks.erase(it, R.checks.end());
    }
    if (any_u && R.k == 1) R.strand_call = '?';
    else if (any_o || (any_f && any_r)) R.strand_call = 'M';
    else R.strand_call = any_f ? 'F' : any_r ? 'R' : '?';
    if (R.strand_call == 'F' || R.strand_call == 'R') {
        const int tmain = R.strand_call == 'F' ? M.main_tpl[0] : M.main_tpl[1];
        bool in_t[MAXANCH];
        for (int a = 0; a < M.A; ++a) in_t[a] = false;
        for (int o : M.tpl[tmain].obs)
            if (M.tpl[tmain].el[o].kind == EK_ANCHOR) in_t[M.tpl[tmain].el[o].anc] = true;
        int pi = 0;
        for (int i = 0; i < nev; ++i) {
            while (pi < np && W.path_ev[pi] < i) ++pi;
            if (pi < np && W.path_ev[pi] == i) continue;
            const Event& e = ev[i];
            if (e.type < M.A && e.cls == CL_S4 && !in_t[e.type] && !M.anch[e.type].shrt) {
                R.strand_call = '?';
                R.flags |= RES_OPPOSITE_EVIDENCE;
                break;
            }
        }
    }
    bool abst = false;
    // a kept D unit does not abstain on p_single: the strand table decided it
    if (R.k == 1) abst = R.p_single < M.opt.tau_single && !(R.segs[0].flags & SEG_D_KEPT);
    else if (R.k >= 2) {
        // nc >= 2: the decode must exclude k = 1; a fold-back split of one construct: that construct must be confident
        const bool flip_split = nc == 1 && R.k == 2 && !R.cut_kind.empty() && R.cut_kind[0] == CUT_STRAND_FLIP;
        abst = nc >= 2 ? R.p_single > 1.f - M.opt.tau_single : (flip_split ? false : R.p_single < M.opt.tau_single);
        for (float p : R.cut_post) abst = abst || p < M.opt.tau_junction;
    }
    if (R.flags & RES_OPPOSITE_EVIDENCE) abst = true;
    if (R.k == 1 && (R.segs[0].flags & SEG_LOWCONF)) abst = true;
    if (abst) R.flags |= RES_ABSTAIN;
}

}  // namespace detail

namespace detail {
// Rescans the gap of each unanchored junction with the full Myers engine. Returns true when events were added.
inline bool rescan_unanchored(const Model& M, const char* seq, int n, Scratch& W) {
    if (!M.opt.windowed_myers) return false;
    const int np = (int)W.path_ev.size();
    bool added = false;
    const size_t ne0 = W.ev.size();
    for (int x = 1; x < np; ++x) {
        if (!W.path_cross[x]) continue;
        const int sa = W.path_st[x - 1], sb = W.path_st[x];
        if ((M.st_flags[sa] & SF_CLOSE) || (M.st_flags[sb] & SF_OPEN)) continue;
        const int lo0 = std::max(0, W.ev[W.path_ev[x - 1]].end - 8), hi0 = std::min(n, W.ev[W.path_ev[x]].start + 8);
        if (hi0 - lo0 < 20) continue;
        int lo = lo0;
        for (const auto& w : W.win) {
            if (w.second <= lo || w.first >= hi0) continue;
            if (w.first > lo) {
                W.raw_x.clear();
                PState fin[MAXPAT];
                myers_all(M.words, M.words_shared, M.Wov, (int)M.pats.size(), seq + lo, w.first - lo, W, fin, W.raw_x);
                raw_to_events(M, W.raw_x, lo, n, false, W.ev, W.raw_ch);
            }
            lo = std::max(lo, w.second);
        }
        if (hi0 > lo) {
            W.raw_x.clear();
            PState fin[MAXPAT];
            myers_all(M.words, M.words_shared, M.Wov, (int)M.pats.size(), seq + lo, hi0 - lo, W, fin, W.raw_x);
            raw_to_events(M, W.raw_x, lo, n, false, W.ev, W.raw_ch);
        }
    }
    if (W.ev.size() > ne0) {
        sort_events(W.ev);
        dedupe_events(M, W.ev);
        drop_phantom_primers(M, W.ev);
        added = true;
    }
    return added;
}
}  // namespace detail

namespace detail {
// Scans the facing window of each anchor-grade inner barcode-side event for a seedless partner anchor (fused
// barcode-side junction). Returns true when events were added (the caller re-decodes).
inline bool facing_partner_scan(const Model& M, const char* seq, int n, Scratch& W) {
    if (M.head_rules.empty() || M.tail_rules.empty()) return false;
    std::vector<Event>& ev = W.ev;
    const size_t ne = ev.size();
    {
        bool any = false;
        for (size_t i = 0; i < ne && !any; ++i) any = ev[i].type < M.A && M.inner_role[ev[i].type] && ev[i].cls == CL_S4;
        if (!any) return false;
    }
    const int np = (int)W.path_ev.size();
    auto is_last = [&](int x) { return x + 1 == np || W.path_cross[x + 1] != 0; };
    auto covered = [&](int a, int lo, int hi) {
        for (const Event& f : ev)
            if (f.type == a && f.start < hi && lo < f.end) return true;
        return false;
    };
    bool added = false;
    auto scan_partner = [&](int a, int lo, int hi, int max_ed) {
        lo = std::max(0, lo);
        hi = std::min(n, hi);
        const Anchor& x = M.anch[a];
        if (hi - lo < x.m || x.nchunk != 1 || covered(a, lo, hi)) return;
        int be = -1;
        const int ed = best_hit(M.single[a], seq + lo, hi - lo, be);
        if (be >= 0 && ed <= std::min(max_ed, x.strong_ed)) {
            push_anchor_event(ev, n, lo + be + 1 - x.m, lo + be + 1, a, ed, CL_S4, 0, 0, x.m, ed);
            added = true;
        }
    };
    for (size_t i = 0; i < ne; ++i) {
        const Event e = ev[i];
        if (e.type >= M.A || e.cls != CL_S4 || !M.inner_role[e.type]) continue;
        int x = -1;
        for (int y = 0; y < np && x < 0 && W.path_ev[y] <= (int)i; ++y)
            if (W.path_ev[y] == (int)i) x = y;
        for (const bc_rule& tr : M.tail_rules)
            for (const bc_rule& hr : M.head_rules) {
                const int B = tr.block + hr.block;
                // the facing construct must fit min_terminal inside the read
                const int mo = M.anch[hr.inner].m, mc = M.anch[tr.inner].m;
                if (tr.inner == e.type) {
                    int lo = e.end + B + RESCUE_GMIN - 2;
                    int hi = std::min(e.end + B + (x < 0 && np ? RESCUE_GMAX_OFF : RESCUE_GMAX) + mo + 2, n - M.min_terminal + hr.block + mo);
                    bool want = x < 0 || is_last(x);
                    if (!want && is_last(x + 1)) {
                        // closed by its primer P: the next unit's partner lies beyond P.end - 6
                        const Event& p = ev[W.path_ev[x + 1]];
                        if (p.type == tr.prim) {
                            lo = std::max(lo, p.end - 6);
                            // the next construct needs evidence of its own beyond the window
                            want = ev[ne - 1].end > hi + 10;
                        }
                    }
                    if (want) scan_partner(hr.inner, lo, hi, M.anch[hr.inner].strong_ed);
                }
                if (hr.inner == e.type) {
                    int lo = std::max(e.start - B - (x < 0 && np ? RESCUE_GMAX_OFF : RESCUE_GMAX) - mc - 2, M.min_terminal - tr.block - mc);
                    int hi = e.start - B - RESCUE_GMIN + 2;
                    bool want = x < 0 || W.path_cross[x];
                    if (!want && x > 0 && W.path_cross[x - 1]) {
                        // opened by its primer P: the previous unit's partner ends before P.start + 6
                        const Event& p = ev[W.path_ev[x - 1]];
                        if (p.type == hr.prim) {
                            hi = std::min(hi, p.start + 6);
                            want = ev[0].start < lo - 10;  // the previous construct needs evidence of its own
                        }
                    }
                    if (want) scan_partner(tr.inner, lo, hi, M.anch[tr.inner].strong_ed);
                }
            }
    }
    if (added) {
        sort_events(ev);
        dedupe_events(M, ev);
        drop_phantom_primers(M, ev);
    }
    return added;
}

// Facing-anchor rescue: queues in Scratch::force each pair of facing anchor-grade partner anchors at the
// fused-junction geometry that marks a junction the decode left out. Returns true when any was queued.
inline bool facing_rescue(const Model& M, const Event* ev, int n, Scratch& W) {
    W.force.clear();
    if (M.head_rules.empty() || M.tail_rules.empty() || n < 2) return false;
    {   // quick reject: no tail-rule inner anchor followed in range by a head-rule one
        int last_c = INT32_MIN / 2;
        bool cand = false;
        int reach = RESCUE_GMAX, bt = 0, bh = 0;
        for (const auto& tr : M.tail_rules) bt = std::max(bt, tr.block);
        for (const auto& hr : M.head_rules) bh = std::max(bh, hr.block);
        reach += bt + bh;
        for (int i = 0; i < n && !cand; ++i) {
            if (ev[i].type >= M.A || ev[i].cls != CL_S4) continue;
            const uint8_t r = M.inner_role[ev[i].type];
            if (r & 1) last_c = ev[i].end;
            if ((r & 2) && ev[i].start <= last_c + reach) cand = true;
        }
        if (!cand) return false;
    }
    const int np = (int)W.path_ev.size();
    W.onp.assign((size_t)n, -1);
    for (int x = 0; x < np; ++x) W.onp[W.path_ev[x]] = x;
    W.bcut.clear();  // cuts of the decoded path (derive_result rebuilds this buffer)
    for (int x = 1; x < np; ++x)
        if (W.path_cross[x]) {
            uint8_t kd;
            W.bcut.push_back(junction_cut(M, ev[W.path_ev[x - 1]], W.path_st[x - 1], ev[W.path_ev[x]], W.path_st[x], kd));
        }
    auto is_last = [&](int x) { return x + 1 == np || W.path_cross[x + 1] != 0; };
    auto overlap = [](const Event& p, int lo, int hi) { return std::min(p.end, hi) - std::max(p.start, lo); };
    for (int ia = 0; ia < n; ++ia) {
        const Event& a = ev[ia];
        if (a.type >= M.A || a.cls != CL_S4) continue;
        for (const bc_rule& tr : M.tail_rules) {
            if (tr.inner != a.type) continue;
            for (const bc_rule& hr : M.head_rules) {
                const int B = tr.block + hr.block;
                for (int ib = ia + 1; ib < n && ev[ib].start <= a.end + B + RESCUE_GMAX; ++ib) {
                    const Event& b = ev[ib];
                    if (b.type != hr.inner || b.cls != CL_S4) continue;
                    const int g = b.start - a.end - B;
                    if (g < RESCUE_GMIN || g > RESCUE_GMAX) continue;
                    const int xa = W.onp[ia], xb = W.onp[ib];
                    if (xa >= 0 && xb >= 0) continue;
                    if (xa < 0 && xb < 0 && g > RESCUE_GMAX_OFF) continue;
                    bool ok = true;
                    if (xa >= 0 && !is_last(xa)) {
                        // look through a closing barcode-adjacent primer laid over the facing barcode block
                        const int xp = xa + 1;
                        ok = is_last(xp) && ev[W.path_ev[xp]].type == tr.prim && overlap(ev[W.path_ev[xp]], b.start - hr.block, b.start) >= 8;
                    } else if (xb >= 0 && !W.path_cross[xb]) {
                        const int xp = xb - 1;
                        ok = xp >= 0 && W.path_cross[xp] && ev[W.path_ev[xp]].type == hr.prim && overlap(ev[W.path_ev[xp]], a.end, a.end + tr.block) >= 8;
                    } else {
                        // no other path event between the two anchors
                        for (int x = 0; x < np && ok; ++x) {
                            const Event& p = ev[W.path_ev[x]];
                            if (xb >= 0) ok = !(p.end > a.end && p.end <= b.start && W.path_ev[x] != ib);
                            else ok = !(p.start >= a.end && p.start < b.start && W.path_ev[x] != ia);
                        }
                    }
                    for (int c : W.bcut) ok = ok && !(c >= a.end && c <= b.start);
                    if (!ok) continue;
                    bool dup = false;
                    for (const auto& fp : W.force) dup = dup || fp.first == ia || fp.second == ib;
                    if (!dup) W.force.push_back({ia, ib});
                }
            }
        }
    }
    return !W.force.empty();
}
}  // namespace detail

// Segments one read. Model is const and shareable; Scratch and Result are per thread.
inline void segment(const Model& M, const char* seq, int len, Scratch& S, Result& R) {
    using namespace detail;
#ifdef CONCAT_HMM_STAGE_TIMERS
    using clk = std::chrono::steady_clock;
    auto t0 = clk::now();
    auto ns = [](clk::time_point a, clk::time_point b) { return (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count(); };
#endif
    R.flags = M.guard_failed ? (uint32_t)RES_GUARD_FAILED : 0u;
    R.k = 0;
    R.boundary_unresolved_seen = false;
    R.cuts.clear(); R.segs.clear(); R.checks.clear(); R.cut_post.clear(); R.cut_lo.clear(); R.cut_hi.clear(); R.cut_kind.clear(); R.cut_flags.clear(); R.cut_aux.clear();
    R.p_single = 0.f;
    R.strand_call = '?';
    if (len <= 0 || !seq) return;
    stage0(M, seq, len, S);
#ifdef CONCAT_HMM_STAGE_TIMERS
    auto t1 = clk::now();
    R.stage_ns[0] = ns(t0, t1);
#endif
    if (M.opt.gate && gate(M, seq, len, S, R)) {
#ifdef CONCAT_HMM_STAGE_TIMERS
        R.stage_ns[1] = ns(t1, clk::now());
        R.stage_ns[2] = R.stage_ns[3] = 0;
#endif
        return;
    }
#ifdef CONCAT_HMM_STAGE_TIMERS
    auto t2 = clk::now();
    R.stage_ns[1] = ns(t1, t2);
#endif
    extract_events(M, seq, len, S, !M.opt.windowed_myers);
#ifdef CONCAT_HMM_STAGE_TIMERS
    auto t3 = clk::now();
    R.stage_ns[2] = ns(t2, t3);
#endif
    int nev = (int)S.ev.size();
    int k = viterbi(M, S.ev.data(), nev, len, S);
    if (k >= 2 && rescan_unanchored(M, seq, len, S)) {
        nev = (int)S.ev.size();
        k = viterbi(M, S.ev.data(), nev, len, S);
    }
    // reads beyond the construct cap get no second look
    if (k <= M.opt.k_cap && facing_partner_scan(M, seq, len, S)) {
        nev = (int)S.ev.size();
        k = viterbi(M, S.ev.data(), nev, len, S);
    }
    bool rescued = false;
    if (k >= 1 && k <= M.opt.k_cap && facing_rescue(M, S.ev.data(), nev, S)) {
        k = viterbi(M, S.ev.data(), nev, len, S);
        for (size_t x = 1; x < S.path_ev.size() && !rescued; ++x)
            if (S.path_cross[x])
                for (const auto& fp : S.force) rescued = rescued || (fp.first == S.path_ev[x - 1] && fp.second == S.path_ev[x]);
        S.force.clear();
    }
    bool fb = false;
    if (k > 0) { forward_backward(M, S.ev.data(), nev, len, S, k >= 2); fb = true; }
    derive_result(M, seq, S.ev.data(), nev, len, S, k, fb, R);
    if (rescued) R.flags |= RES_RESCUED;
#ifdef CONCAT_HMM_STAGE_TIMERS
    R.stage_ns[3] = ns(t3, clk::now());
#endif
}

namespace detail {
inline std::string row_or(const Model& M, const char* klass, char dir, const char* canon) {
    for (auto& e : M.spec.elements)
        if (lower(e.klass) == klass && e.direction == dir) return e.id;
    // A one-strand layout still has a mirrored decoder. Store its reverse
    // tables on an existing row, in their distinct hmm_* columns, so the
    // position-map writer does not silently discard synthetic rc_* rows.
    for (auto& e : M.spec.elements)
        if (lower(e.klass) == klass) return e.id;
    return canon;
}
inline std::string insert_row(const Model& M, char dir, const char* canon) {
    for (auto& e : M.spec.elements)
        if (is_insert_elem(e) && e.direction == dir) return e.id;
    for (auto& e : M.spec.elements)
        if (is_insert_elem(e)) return e.id;
    return canon;
}
// Packs the tight tables of `row` that differ from the prior (max |log ratio| > 0.05) as [count, then per table:
// id, lo, ncore, ncore pmf values, ntail, ntail 10-bp bin means]. Canonical: pack(unpack(pack(v))) == pack(v).
inline void pack_tables(const Model& M, const PMap& par, const std::string& row, std::vector<double>& out) {
    out.clear();
    out.push_back(0);
    std::vector<double> tail, expanded;
    for (size_t t = 0; t < M.tight_keys.size(); ++t) {
        if (M.tight_row[t] != row) continue;
        auto it = par.find(M.tight_keys[t]);
        if (it == par.end()) continue;
        const std::vector<double>& v = it->second;
        const std::vector<double>& pv = M.prior.at(M.tight_keys[t]);
        if (v.size() != pv.size()) continue;
        const int lo = M.tight_lo[t];
        const int ncore = std::min((int)v.size(), 2 * (-lo) + 1);
        const int rest = (int)v.size() - ncore;
        const int ntail = (rest + 9) / 10;
        tail.assign((size_t)ntail, 0.0);
        for (int b = 0; b < ntail; ++b) {
            const int i0 = ncore + b * 10, i1 = std::min((int)v.size(), i0 + 10);
            bool equal = true;
            double sum = 0;
            for (int i = i0; i < i1; ++i) { sum += v[i]; equal = equal && v[i] == v[i0]; }
            tail[b] = equal ? v[i0] : sum / (i1 - i0);
        }
        expanded.assign(v.begin(), v.begin() + ncore);
        for (int i = ncore; i < (int)v.size(); ++i) expanded.push_back(tail[(i - ncore) / 10]);
        double dev = 0;
        for (size_t i = 0; i < expanded.size(); ++i)
            dev = std::max(dev, std::fabs(std::log(std::max(expanded[i], 1e-12) / std::max(pv[i], 1e-12))));
        if (dev < 0.05) continue;
        out.push_back((double)t);
        out.push_back(lo);
        out.push_back(ncore);
        out.insert(out.end(), v.begin(), v.begin() + ncore);
        out.push_back(ntail);
        out.insert(out.end(), tail.begin(), tail.end());
        out[0] += 1;
    }
}
inline bool cell_int(double x, int lo, int hi, int& out) {
    if (!std::isfinite(x) || x < lo - 0.5 || x > hi + 0.5) return false;
    out = (int)std::lround(x);
    return true;
}
// A valid pmf (finite, >= 0, > 0 where the prior is); renormalised only when |sum - 1| > 1e-3 (else bit-exact).
inline bool sane_pmf(std::vector<double>& v, const std::vector<double>& prior) {
    if (v.size() != prior.size()) return false;
    double tot = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        if (!std::isfinite(v[i]) || v[i] < 0 || (prior[i] > 0 && v[i] <= 0)) return false;
        tot += v[i];
    }
    if (!(tot > 0) || !std::isfinite(tot)) return false;
    if (std::fabs(tot - 1.0) > 1e-3)
        for (auto& x : v) x /= tot;
    return true;
}
// Null rates: finite, in [0, 1], and >= 1e-12 where the emission prior is positive.
inline bool sane_rates(const std::vector<double>& v, const std::vector<double>& q_prior) {
    if (v.size() != q_prior.size()) return false;
    for (size_t i = 0; i < v.size(); ++i)
        if (!std::isfinite(v[i]) || v[i] < 0 || v[i] > 1 || (q_prior[i] > 0 && v[i] < 1e-12)) return false;
    return true;
}
// False when the packed cell is malformed; valid tables of this row are stored in `par`.
inline bool unpack_tables(const Model& M, PMap& par, const std::string& row, const std::vector<double>& in) {
    if (in.empty()) return true;
    int ntab = 0;
    if (!cell_int(in[0], 0, 1 << 20, ntab)) return false;
    size_t pos = 1;
    bool ok = true;
    for (int k = 0; k < ntab; ++k) {
        int t = 0, lo = 0, ncore = 0, ntail = 0;
        if (pos + 3 > in.size()) return false;
        if (!cell_int(in[pos], -1, 1 << 20, t) || !cell_int(in[pos + 1], -(1 << 20), 1 << 20, lo) || !cell_int(in[pos + 2], 0, 1 << 20, ncore)) return false;
        pos += 3;
        if (pos + (size_t)ncore + 1 > in.size()) return false;
        const size_t core_at = pos;
        pos += (size_t)ncore;
        if (!cell_int(in[pos], 0, 1 << 20, ntail)) return false;
        pos += 1;
        if (pos + (size_t)ntail > in.size()) return false;
        const size_t tail_at = pos;
        pos += (size_t)ntail;
        if (t < 0 || t >= (int)M.tight_keys.size() || M.tight_row[t] != row || lo != M.tight_lo[t]) continue;
        const int n = M.tight_hi[t] - M.tight_lo[t] + 1;
        if (ncore > n || ncore + ntail * 10 < n) { ok = false; continue; }
        std::vector<double> v((size_t)n);
        for (int i = 0; i < n; ++i) v[i] = i < ncore ? in[core_at + i] : in[tail_at + (size_t)((i - ncore) / 10)];
        if (!sane_pmf(v, M.prior.at(M.tight_keys[t]))) { ok = false; continue; }
        par[M.tight_keys[t]] = v;
    }
    return ok;
}
}  // namespace detail

inline Params Model::export_params() const {
    using namespace detail;
    Params P;
    auto get = [&](const std::string& k) -> const std::vector<double>& {
        auto it = par.find(k);
        if (it != par.end()) return it->second;
        return prior.at(k);
    };
    const std::string r_start = row_or(*this, "start", 'F', "seq_start"), r_stop = row_or(*this, "stop", 'F', "seq_stop");
    const std::string r_rstart = row_or(*this, "start", 'R', "rc_seq_start"), r_rstop = row_or(*this, "stop", 'R', "rc_seq_stop");
    for (int a = 0; a < A; ++a) {
        auto& row = P.cells[anch[a].row];
        row["hmm_q"] = get("q." + anch[a].seq);
        row["hmm_qpal"] = get("qpal." + anch[a].seq);
        row["hmm_null"] = get("null." + anch[a].seq);
    }
    for (int q = 0; q < Q; ++q) {
        auto& row = P.cells[poly_row[q]];
        row["hmm_q"] = get("q.poly" + std::string(1, poly_base[q]));
        row["hmm_null"] = get("null.poly" + std::string(1, poly_base[q]));
    }
    P.cells[insert_row(*this, 'F', "read")]["hmm_cdna"] = get("ins1");
    P.cells[insert_row(*this, 'R', "rc_read")]["hmm_cdna2"] = get("ins2");
    {   // cDNA strand table, only when calibrated: [mean sense score, 4096 log-odds]
        auto it = par.find("strand");
        if (it != par.end() && it->second.size() == (size_t)STRAND_N + 1) P.cells[insert_row(*this, 'F', "read")]["hmm_strand"] = it->second;
        auto itd = par.find("strand_d");
        if (itd != par.end() && itd->second.size() == (size_t)STRAND_D_N) P.cells[insert_row(*this, 'F', "read")]["hmm_strand_d"] = itd->second;
    }
    P.cells[r_start]["hmm_begin"] = get("begin");
    P.cells[r_stop]["hmm_end"] = get("end");
    {
        std::vector<double> g = {get("p0")[0]};
        for (double x : get("gate")) g.push_back(x);
        P.cells[r_start]["hmm_glob"] = g;
        P.cells[r_start]["hmm_meta"] = {2.0, (double)topo_hash, guard_failed ? 1.0 : 0.0};
    }
    {
        std::vector<double> pv;
        for (auto& n : pres_names) pv.push_back(get("pres." + n)[0]);
        P.cells[r_stop]["hmm_pres"] = pv;
    }
    {
        std::vector<double> jv;
        for (int t = 0; t < NT; ++t)
            for (double x : get("junc." + tpl[t].name)) jv.push_back(x);
        P.cells[r_rstart]["hmm_junc"] = jv;
    }
    {
        std::vector<std::string> rows;
        for (auto& r : tight_row)
            if (std::find(rows.begin(), rows.end(), r) == rows.end()) rows.push_back(r);
        for (auto& r : rows) {
            std::vector<double> v;
            pack_tables(*this, par, r, v);
            if (r.empty()) P.cells[r_rstop]["hmm_tight"] = v;
            else P.cells[r]["hmm_spacer"] = v;
        }
    }
    P.guard_failed = guard_failed;
    P.topology = topo_hash;
    return P;
}

// Applies cached parameters. An invalid cell keeps the layout prior, is listed in param_warnings and marks the
// model guard_failed.
inline void Model::apply_params(const Params& P) {
    using namespace detail;
    param_warnings.clear();
    auto reject = [&](const std::string& row, const char* col, const char* why) {
        param_warnings.push_back("hmm_* cell " + row + "/" + col + " rejected (" + why + "); layout prior used");
    };
    auto cell = [&](const std::string& row, const char* col) -> const std::vector<double>* {
        auto r = P.cells.find(row);
        if (r == P.cells.end()) return nullptr;
        auto c = r->second.find(col);
        if (c == r->second.end() || c->second.empty()) return nullptr;
        return &c->second;
    };
    // a size mismatch is a different layout version (ignored)
    auto setv = [&](const std::string& key, const std::string& row, const char* col, bool rates) {
        const std::vector<double>* v = cell(row, col);
        if (!v) return;
        auto it = prior.find(key);
        if (it == prior.end() || it->second.size() != v->size()) return;
        std::vector<double> x = *v;
        // null rates are checked against the emission prior of the same type ("null.X" -> "q.X")
        const bool ok = rates ? sane_rates(x, prior.at("q." + key.substr(5))) : sane_pmf(x, it->second);
        if (!ok) { reject(row, col, rates ? "rates must be finite, in [0, 1], positive where the prior is" : "not a valid pmf"); return; }
        par[key] = x;
    };
    const std::string r_start = row_or(*this, "start", 'F', "seq_start"), r_stop = row_or(*this, "stop", 'F', "seq_stop");
    const std::string r_rstart = row_or(*this, "start", 'R', "rc_seq_start"), r_rstop = row_or(*this, "stop", 'R', "rc_seq_stop");
    for (const auto& b : P.bad_cells) param_warnings.push_back("hmm_* cell " + b + " rejected (not a list of numbers); layout prior used");
    uint64_t th = P.topology;
    bool gf = P.guard_failed;
    if (auto m = cell(r_start, "hmm_meta")) {
        bool ok = true;
        for (double x : *m) ok = ok && std::isfinite(x);
        if (m->size() >= 2) ok = ok && (*m)[1] >= 0 && (*m)[1] < 9007199254740992.0;
        if (!ok) { reject(r_start, "hmm_meta", "non-finite or out-of-range value"); th = 0; }
        else {
            if (m->size() >= 2) th = (uint64_t)std::llround((*m)[1]);
            if (m->size() >= 3) gf = gf || (*m)[2] > 0.5;
        }
    }
    for (int a = 0; a < A; ++a) {
        setv("q." + anch[a].seq, anch[a].row, "hmm_q", false);
        setv("qpal." + anch[a].seq, anch[a].row, "hmm_qpal", false);
        setv("null." + anch[a].seq, anch[a].row, "hmm_null", true);
    }
    for (int q = 0; q < Q; ++q) {
        setv("q.poly" + std::string(1, poly_base[q]), poly_row[q], "hmm_q", false);
        setv("null.poly" + std::string(1, poly_base[q]), poly_row[q], "hmm_null", true);
    }
    setv("ins1", insert_row(*this, 'F', "read"), "hmm_cdna", false);
    setv("ins2", insert_row(*this, 'R', "rc_read"), "hmm_cdna2", false);
    if (auto sv = cell(insert_row(*this, 'F', "read"), "hmm_strand")) {
        bool ok = sv->size() == (size_t)STRAND_N + 1 && (*sv)[0] > 0;
        for (double x : *sv) ok = ok && std::isfinite(x);
        if (!ok) reject(insert_row(*this, 'F', "read"), "hmm_strand", "must hold a positive mean and 4096 finite log-odds");
        else par["strand"] = *sv;
    }
    if (auto sv = cell(insert_row(*this, 'F', "read"), "hmm_strand_d")) {
        bool ok = sv->size() == (size_t)STRAND_D_N;
        for (double x : *sv) ok = ok && std::isfinite(x) && x >= 0;
        if (!ok) reject(insert_row(*this, 'F', "read"), "hmm_strand_d", "must hold 262144 finite non-negative counts");
        else par["strand_d"] = *sv;
    }
    if (auto g = cell(r_start, "hmm_glob")) {
        const std::vector<double>& G = *g;
        bool ok = true;
        for (double x : G) ok = ok && std::isfinite(x);
        if (ok && G.size() >= 1) ok = G[0] > 0 && G[0] < 1;                       // p0
        if (ok && G.size() >= 4) ok = G[1] >= 0 && G[2] >= 0 && G[3] >= 0 && G[3] <= 1;  // zone5, zone3, fast_p
        if (!ok) reject(r_start, "hmm_glob", "p0 must lie in (0, 1), zones >= 0, fast-path p in [0, 1], all finite");
        else {
            if (G.size() >= 1) par["p0"] = {G[0]};
            if (G.size() >= 4) par["gate"] = {G[1], G[2], G[3]};
            if (G.size() >= 5) par["guard"] = {G[4]};
        }
    }
    if (th == topo_hash) {
        setv("begin", r_start, "hmm_begin", false);
        setv("end", r_stop, "hmm_end", false);
        if (auto pv = cell(r_stop, "hmm_pres"))
            if (pv->size() == pres_names.size()) {
                bool ok = true;
                for (double x : *pv) ok = ok && std::isfinite(x) && x >= 0 && x <= 1;
                if (!ok) reject(r_stop, "hmm_pres", "presence probabilities must be finite and in [0, 1]");
                else
                    for (size_t i = 0; i < pres_names.size(); ++i) par["pres." + pres_names[i]] = {(*pv)[i]};
            }
        if (auto jv = cell(r_rstart, "hmm_junc"))
            if ((int)jv->size() == NT * NT) {
                bool ok = true;
                std::vector<std::vector<double>> rows((size_t)NT);
                for (int t = 0; t < NT && ok; ++t) {
                    rows[t].assign(jv->begin() + t * NT, jv->begin() + (t + 1) * NT);
                    ok = sane_pmf(rows[t], prior.at("junc." + tpl[t].name));
                }
                if (!ok) reject(r_rstart, "hmm_junc", "not a valid pmf per template");
                else
                    for (int t = 0; t < NT; ++t) par["junc." + tpl[t].name] = rows[t];
            }
        if (auto tv = cell(r_rstop, "hmm_tight"))
            if (!unpack_tables(*this, par, "", *tv)) reject(r_rstop, "hmm_tight", "malformed or invalid duration table");
        for (size_t t = 0; t < tight_row.size(); ++t) {
            const std::string& r = tight_row[t];
            if (r.empty() || std::find(tight_row.begin(), tight_row.begin() + t, r) != tight_row.begin() + t) continue;
            if (auto sv = cell(r, "hmm_spacer"))
                if (!unpack_tables(*this, par, r, *sv)) reject(r, "hmm_spacer", "malformed or invalid duration table");
        }
    }
    if (gf || !param_warnings.empty()) par["guard"] = {1};
}

// element id -> {hmm_* column -> ';'-joined numbers}, each in its shortest exact round-trip form (lossless).
inline std::map<std::string, std::map<std::string, std::string>> params_to_columns(const Params& P) {
    std::map<std::string, std::map<std::string, std::string>> out;
    char buf[64];
    for (auto& r : P.cells)
        for (auto& c : r.second) {
            std::string s;
            for (size_t i = 0; i < c.second.size(); ++i) {
                const double x = c.second[i];
                if (std::isfinite(x) && x == std::floor(x) && std::fabs(x) < 9007199254740992.0) snprintf(buf, sizeof buf, "%.0f", x);
                else {
                    for (int prec = 15; prec <= 17; ++prec) {
                        snprintf(buf, sizeof buf, "%.*g", prec, x);
                        if (!std::isfinite(x) || std::strtod(buf, nullptr) == x) break;
                    }
                }
                if (i) s += ';';
                s += buf;
            }
            out[r.first][c.first] = s;
        }
    return out;
}
// Parses hmm_* cells as numbers only; Model::apply_params validates them.
inline Params params_from_columns(const layout_spec& spec, const std::map<std::string, std::map<std::string, std::string>>& cols) {
    (void)spec;
    Params P;
    for (auto& r : cols)
        for (auto& c : r.second) {
            if (c.first.rfind("hmm_", 0) != 0 || c.second.empty()) continue;
            std::vector<double> v;
            const char* p = c.second.c_str();
            bool ok = true;
            while (*p) {
                char* q = nullptr;
                double x = std::strtod(p, &q);
                if (q == p) { ok = false; break; }
                v.push_back(x);
                p = q;
                while (*p == ';' || *p == ' ') ++p;
            }
            if (!ok || v.empty()) { P.bad_cells.push_back(r.first + "/" + c.first); continue; }
            P.cells[r.first][c.first] = v;
            if (c.first == "hmm_meta") {
                if (v.size() >= 2 && std::isfinite(v[1]) && v[1] >= 0 && v[1] < 9007199254740992.0) P.topology = (uint64_t)std::llround(v[1]);
                if (v.size() >= 3) P.guard_failed = !(v[2] <= 0.5);  // NaN counts as failed
            }
        }
    return P;
}

// Reads a RAD position map CSV into id -> {column -> cell} (all columns).
inline std::map<std::string, std::map<std::string, std::string>> read_position_map(const std::string& path, std::vector<std::string>* header = nullptr,
                                                                                     std::vector<std::string>* row_order = nullptr) {
    using namespace detail;
    std::map<std::string, std::map<std::string, std::string>> out;
    std::string t = read_file(path);
    if (t.empty()) return out;
    std::stringstream ss(t);
    std::string ln;
    std::vector<std::string> hdr;
    while (std::getline(ss, ln)) {
        if (!ln.empty() && ln.back() == '\r') ln.pop_back();
        if (trim(ln).empty()) continue;
        auto f = csv_fields(ln);
        if (hdr.empty()) { hdr = f; continue; }
        if (f.empty() || f[0].empty()) continue;
        auto& row = out[f[0]];
        for (size_t i = 1; i < hdr.size() && i < f.size(); ++i) row[hdr[i]] = f[i];
        if (row_order) row_order->push_back(f[0]);
    }
    if (header) *header = hdr;
    return out;
}

// Calibration: hard-EM (Viterbi training), Dirichlet-smoothed toward the layout priors.
namespace detail {
struct Stats {
    std::vector<std::vector<double>> q, qpal, noise, opp, tight;
    std::vector<double> ins1, ins2, pyes, pno, begin, end;
    std::vector<std::vector<double>> junc;
    double n0 = 0, nreads = 0, opp_bases = 0;
    std::vector<int> z5, z3;
    double fast_n = 0, fast_ok = 0;
    std::vector<double> tbeg, topen, tend, tclose;  // per template: constructs starting / ending, at open / close slot
    void init(const Model& M) {
        q.assign(M.NTYPE, std::vector<double>(NFEAT, 0.0));
        qpal = q; noise = q; opp = q;
        tight.clear();
        for (size_t t = 0; t < M.tight_keys.size(); ++t) tight.push_back(std::vector<double>(M.tight_hi[t] - M.tight_lo[t] + 1, 0.0));
        ins1.assign(NINSBIN, 0.0); ins2.assign(NINSBIN, 0.0);
        pyes.assign(M.pres_names.size(), 0.0); pno = pyes;
        begin.assign(M.S, 0.0); end = begin;
        junc.assign(M.NT, std::vector<double>(M.NT, 0.0));
        n0 = nreads = opp_bases = 0;
        z5.clear(); z3.clear();
        fast_n = fast_ok = 0;
        tbeg.assign(M.NT, 0); topen = tbeg; tend = tbeg; tclose = tbeg;
    }
};
inline std::vector<double> mix_norm(const std::vector<double>& cnt, const std::vector<double>& prior, double alpha, double floor_v = 1e-7) {
    std::vector<double> out(cnt.size());
    double tp = 0;
    for (size_t i = 0; i < cnt.size(); ++i) tp += prior[i];
    double tot = 0;
    for (size_t i = 0; i < cnt.size(); ++i) {
        out[i] = cnt[i] + alpha * (tp > 0 ? prior[i] / tp : 0);
        tot += out[i];
    }
    for (size_t i = 0; i < cnt.size(); ++i) out[i] = prior[i] > 0 || cnt[i] > 0 ? std::max(out[i] / tot, floor_v) : 0.0;
    return out;
}
inline uint64_t path_sig(const Scratch& W) {
    uint64_t h = 1469598103934665603ULL;
    for (size_t x = 0; x < W.path_ev.size(); ++x) h = (h ^ (uint64_t)(W.path_ev[x] * 131 + W.path_st[x] * 7 + W.path_cross[x])) * 1099511628211ULL;
    return h;
}
inline int quant(std::vector<int>& v, double q, int dflt) {
    if (v.size() < 50) return dflt;
    std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, (size_t)(q * v.size()))];
}
}  // namespace detail

class Calibrator {
public:
    explicit Calibrator(const Model& prior) : M_(std::make_shared<Model>(prior)), rng_(0x5eed1234ULL) {
        M_->par = M_->prior;
        M_->par["guard"] = {0};
        M_->finalize();
        shuf_.assign(M_->NTYPE, std::vector<double>(detail::NFEAT, 0.0));
        strand_cnt_.assign(detail::STRAND_N, 0.0);
        strand_cnt_d_.assign(detail::STRAND_D_N, 0u);
    }
    // Accumulates one read (call on every calibration read).
    void add_read(const char* seq, int len) {
        using namespace detail;
        if (len <= 0 || !seq) return;
        const Model& M = *M_;
        stage0(M, seq, len, W_);
        Result r;
        const bool g = gate(M, seq, len, W_, r);
        gate_.push_back(g ? (uint8_t)r.strand_call : 0);
        if (g && r.k == 1 && (r.strand_call == 'F' || r.strand_call == 'R')) count_strand(seq, len, r);
        extract_events(M, seq, len, W_, true);
        ev_.insert(ev_.end(), W_.ev.begin(), W_.ev.end());
        off_.push_back((uint32_t)ev_.size());
        len_.push_back(len);
        bases_ += len;
        // base-shuffled null (whole read, normal + relaxed patterns)
        W_.shuf.assign(seq, seq + len);
        {   // per-read seed: the null does not depend on thread assignment or merge order
            uint64_t h = 1469598103934665603ULL;
            for (int i = 0; i < len; ++i) h = (h ^ (unsigned char)seq[i]) * 1099511628211ULL;
            rng_.seed(h ^ 0x5eed1234ULL);
        }
        for (int i = len - 1; i > 0; --i) std::swap(W_.shuf[i], W_.shuf[rng_() % (uint64_t)(i + 1)]);
        stage0(M, W_.shuf.data(), len, W_);
        extract_events(M, W_.shuf.data(), len, W_, true);
        for (auto& e : W_.ev) shuf_[e.type][e.feat] += 1;
        if (!M.rwords.empty()) {
            W_.raw_x.clear();
            PState fin[MAXPAT];
            myers_all(M.rwords, M.rwords_shared, M.Wov, (int)M.pats.size(), W_.shuf.data(), len, W_, fin, W_.raw_x);
            for (auto& rr : W_.raw_x) {
                const Pattern& p = M.pats[rr.pat];
                if (rr.score > M.anch[p.anchor].k && rr.score < FEAT_PREFIX) shuf_[p.anchor][rr.score] += 1;
            }
        }
        null_bases_ += len;
        ++nreads_;
    }
    void merge(const Calibrator& o) {
        const uint32_t base = (uint32_t)ev_.size();
        ev_.insert(ev_.end(), o.ev_.begin(), o.ev_.end());
        for (uint32_t x : o.off_) off_.push_back(x + base);
        len_.insert(len_.end(), o.len_.begin(), o.len_.end());
        gate_.insert(gate_.end(), o.gate_.begin(), o.gate_.end());
        for (size_t t = 0; t < shuf_.size() && t < o.shuf_.size(); ++t)
            for (size_t f = 0; f < shuf_[t].size(); ++f) shuf_[t][f] += o.shuf_[t][f];
        bases_ += o.bases_;
        null_bases_ += o.null_bases_;
        nreads_ += o.nreads_;
        for (size_t i = 0; i < strand_cnt_.size() && i < o.strand_cnt_.size(); ++i) strand_cnt_[i] += o.strand_cnt_[i];
        for (size_t i = 0; i < strand_cnt_d_.size() && i < o.strand_cnt_d_.size(); ++i) strand_cnt_d_[i] += o.strand_cnt_d_[i];
        strand_reads_ += o.strand_reads_;
        strand_bases_ += o.strand_bases_;
    }
    // Runs hard-EM and the calibration guard; returns the parameters to cache (layout priors if the guard fails).
    // The caller invokes this on its main thread. Never pass an R callback
    // into add_read(), which may run in a calibration worker.
    Params finalize(int max_iter = 20, double tol = 0.001,
                    const std::function<void()>& boundary_poll = {});
    bool layout_matches() const { return matches_; }
    std::string report() const { return report_; }
    const Model& model() const { return *M_; }
    size_t n_reads() const { return nreads_; }
    size_t strand_reads() const { return strand_reads_; }

private:
    void accumulate(const detail::Event* ev, int n, int L, int k, detail::Stats& st, uint8_t gate_strand);
    void update(const detail::Stats& st);
    // Counts the sense-strand k-mers of the cDNA of a gate-called single-strand construct. Integer counts: the merge
    // order of workers cannot change the table.
    void count_strand(const char* seq, int len, const Result& r) {
        using namespace detail;
        const Model& M = *M_;
        const char st = r.strand_call;
        // bounds: the observed non-poly elements nearest the insert (a poly tail only when its side has nothing else)
        int head = -1, tail = len + 1, head_poly = -1, tail_poly = len + 1;
        for (const auto& c : r.checks) {
            if (c.status == Status::MISSING_EXPECTED || c.strand != st || c.element < 0 || c.element >= (int)M.spec.elements.size()) continue;
            const layout_element& e = M.spec.elements[c.element];
            int ins_order = INT32_MIN;
            for (const auto& x : M.spec.elements)
                if (x.direction == e.direction && is_insert_elem(x)) { ins_order = x.order; break; }
            if (ins_order == INT32_MIN || e.order == ins_order) continue;
            int side = e.order < ins_order ? -1 : 1;
            if (e.direction != c.strand) side = -side;  // reverse-complement slot derived from a forward row: read order reversed
            const bool poly = lower(e.klass) == "poly_tail";
            if (side < 0) (poly ? head_poly : head) = std::max(poly ? head_poly : head, poly ? c.start + M.opt.seen_pad : c.end - M.opt.seen_pad);
            else (poly ? tail_poly : tail) = std::min(poly ? tail_poly : tail, poly ? c.end - M.opt.seen_pad : c.start + M.opt.seen_pad);
        }
        if (head < 0) head = head_poly;
        if (tail > len) tail = tail_poly;
        if (head < 0 || tail > len) return;
        const int a = std::max(0, head + STRAND_MARGIN_BP), b = std::min(len, tail - STRAND_MARGIN_BP);
        if (b - a < STRAND_MIN_CDNA) return;
        const uint32_t km = (uint32_t)STRAND_N - 1, kmd = (uint32_t)STRAND_D_N - 1;
        uint32_t v = 0, vd = 0;
        int valid = 0;
        if (st == 'F') {
            for (int i = a; i < b; ++i) {
                const unsigned char ch = (unsigned char)seq[i];
                const uint32_t c2 = code2(ch);
                v = ((v << 2) | c2) & km;
                vd = ((vd << 2) | c2) & kmd;
                valid = (ch | 0x20) == 'n' ? 0 : valid + 1;
                if (valid >= STRAND_K) strand_cnt_[v] += 1.0;
                if (valid >= STRAND_D_K) ++strand_cnt_d_[vd];
            }
        } else {
            for (int i = b - 1; i >= a; --i) {
                const unsigned char ch = (unsigned char)seq[i];
                const uint32_t c2 = 3u - code2(ch);
                v = ((v << 2) | c2) & km;
                vd = ((vd << 2) | c2) & kmd;
                valid = (ch | 0x20) == 'n' ? 0 : valid + 1;
                if (valid >= STRAND_K) strand_cnt_[v] += 1.0;
                if (valid >= STRAND_D_K) ++strand_cnt_d_[vd];
            }
        }
        ++strand_reads_;
        strand_bases_ += b - a;
    }
    std::vector<double> strand_cnt_;
    std::vector<uint32_t> strand_cnt_d_;
    size_t strand_reads_ = 0;
    double strand_bases_ = 0;
    std::shared_ptr<Model> M_;
    Scratch W_;
    std::vector<detail::Event> ev_;
    std::vector<uint32_t> off_;
    std::vector<int32_t> len_;
    std::vector<uint8_t> gate_;
    std::vector<std::vector<double>> shuf_;
    double bases_ = 0, null_bases_ = 0;
    size_t nreads_ = 0;
    std::mt19937_64 rng_;
    bool matches_ = false;
    std::string report_;
};

inline void Calibrator::accumulate(const detail::Event* ev, int n, int L, int k, detail::Stats& st, uint8_t gate_strand) {
    using namespace detail;
    const Model& M = *M_;
    const Scratch& W = W_;
    st.nreads += 1;
    const int np = (int)W.path_ev.size();
    if (gate_strand) {
        st.fast_n += 1;
        bool ok = k == 1 && np > 0 && M.tpl[M.st_tpl[W.path_st[0]]].strand == (char)gate_strand;
        st.fast_ok += ok;
    }
    if (np == 0) {
        st.n0 += 1;
        for (int i = 0; i < n; ++i) st.noise[ev[i].type][ev[i].feat] += 1;
        return;
    }
    int pi = 0;
    // opposite-strand real null: reads decoded as one main-template construct
    const int t0 = M.st_tpl[W.path_st[0]];
    const bool single_main = k == 1 && M.tpl[t0].main;
    bool in_t[MAXANCH];
    for (int a = 0; a < M.A; ++a) in_t[a] = false;
    if (single_main)
        for (int o : M.tpl[t0].obs)
            if (M.tpl[t0].el[o].kind == EK_ANCHOR) in_t[M.tpl[t0].el[o].anc] = true;
    if (single_main) st.opp_bases += L;
    for (int i = 0; i < n; ++i) {
        while (pi < np && W.path_ev[pi] < i) ++pi;
        if (pi < np && W.path_ev[pi] == i) continue;
        st.noise[ev[i].type][ev[i].feat] += 1;
        if (single_main && ev[i].type < M.A && !in_t[ev[i].type]) st.opp[ev[i].type][ev[i].feat] += 1;
    }
    for (int x = 0; x < np; ++x) {
        const Event& e = ev[W.path_ev[x]];
        bool pal = false;
        if (x > 0 && W.path_cross[x]) pal = M.cross[W.path_st[x - 1] * M.S + W.path_st[x]].pal != 0;
        (pal ? st.qpal : st.q)[e.type][e.feat] += 1;
    }
    st.begin[W.path_st[0]] += 1;
    st.end[W.path_st[np - 1]] += 1;
    // per-template construct starts / ends (guard statistics) and gate zones
    for (int x = 0; x < np; ++x) {
        const int s = W.path_st[x];
        const int t = M.st_tpl[s];
        if (W.path_cross[x]) { st.tbeg[t] += 1; if (M.st_flags[s] & SF_OPEN) st.topen[t] += 1; }
        if (x + 1 == np || W.path_cross[x + 1]) { st.tend[t] += 1; if (M.st_flags[s] & SF_CLOSE) st.tclose[t] += 1; }
    }
    if (single_main) {
        const int s0 = W.path_st[0], s1 = W.path_st[np - 1];
        const Event& e0 = ev[W.path_ev[0]];
        const Event& e1 = ev[W.path_ev[np - 1]];
        if (M.head_fix[s0] >= 0) st.z5.push_back(std::max(0, e0.start - M.head_fix[s0]));
        if (M.tail_fix[s1] >= 0) st.z3.push_back(std::max(0, L - (e1.end + M.tail_fix[s1])));
    }
    for (int x = 1; x < np; ++x) {
        const int sp = W.path_st[x - 1], s2 = W.path_st[x];
        const bool cr = W.path_cross[x];
        const desc_hot& d = cr ? M.cross[sp * M.S + s2] : M.same[sp * M.S + s2];
        const desc_cold& c = cr ? M.cross_c[sp * M.S + s2] : M.same_c[sp * M.S + s2];
        for (auto& pt : c.pres) (pt.second ? st.pyes : st.pno)[pt.first] += 1;
        const int gap = ev[W.path_ev[x]].start - ev[W.path_ev[x - 1]].end;
        const int xg = gap - d.fixed;
        if (d.kind == 0) {
            const int lo = M.tight_lo[d.tid];
            if (xg >= lo && xg <= M.tight_hi[d.tid]) st.tight[d.tid][xg - lo] += 1;
        } else if (d.kind == 1) st.ins1[ins_bin(xg)] += 1;
        else st.ins2[ins_bin(xg)] += 1;
        if (cr) st.junc[c.jfrom][c.jto] += 1;
    }
}

inline void Calibrator::update(const detail::Stats& st) {
    using namespace detail;
    Model& M = *M_;
    PMap& par = M.par;
    const PMap& pr = M.prior;
    for (int ty = 0; ty < M.NTYPE; ++ty) {
        const std::string nm = M.type_name(ty);
        par["q." + nm] = mix_norm(st.q[ty], pr.at("q." + nm), 20.0);
        if (ty < M.A) par["qpal." + nm] = mix_norm(st.qpal[ty], par["q." + nm], 50.0);
        const auto& qq = pr.at("q." + nm);
        std::vector<double> nl(NFEAT, 0.0);
        for (int f = 0; f < NFEAT; ++f) {
            if (qq[f] <= 0) continue;
            if (ty < M.A) {
                nl[f] = null_bases_ > 0 ? (shuf_[ty][f] + 0.5) / null_bases_ : pr.at("null." + nm)[f];
                // real-read opposite-strand null (single-construct reads) for degraded features only
                if (f > M.anch[ty].strong_ed && f < FEAT_PREFIX && st.opp_bases > 1e5)
                    nl[f] = std::max(nl[f], (st.opp[ty][f] + 0.5) / st.opp_bases);
            } else nl[f] = (st.noise[ty][f] + 0.5) / std::max(bases_, 1.0);
            nl[f] = std::max(nl[f], 1e-10);
        }
        par["null." + nm] = nl;
    }
    for (size_t t = 0; t < M.tight_keys.size(); ++t) {
        const auto& c = st.tight[t];
        std::vector<double> sm(c.size(), 0.0);
        for (size_t i = 0; i < c.size(); ++i) {
            sm[i] += 0.5 * c[i];
            if (i > 0) sm[i - 1] += 0.25 * c[i];
            if (i + 1 < c.size()) sm[i + 1] += 0.25 * c[i];
        }
        par[M.tight_keys[t]] = mix_norm(sm, pr.at(M.tight_keys[t]), 20.0, 1e-7);
    }
    {
        std::vector<double> sm(NINSBIN, 0.0);
        for (int b = 0; b < NINSBIN; ++b) {
            sm[b] += 0.6 * st.ins1[b];
            if (b > 0) sm[b - 1] += 0.2 * st.ins1[b];
            if (b + 1 < NINSBIN) sm[b + 1] += 0.2 * st.ins1[b];
        }
        par["ins1"] = mix_norm(sm, pr.at("ins1"), 50.0, 1e-8);
        std::vector<double> p2 = convolve_ins(par["ins1"]);
        par["ins2"] = mix_norm(st.ins2, p2, 30.0, 1e-9);
    }
    for (size_t i = 0; i < M.pres_names.size(); ++i) {
        const double p0 = pr.at("pres." + M.pres_names[i])[0];
        par["pres." + M.pres_names[i]] = {std::min(std::max((st.pyes[i] + 10 * p0) / (st.pyes[i] + st.pno[i] + 10), 0.001), 0.999)};
    }
    for (int t = 0; t < M.NT; ++t) par["junc." + M.tpl[t].name] = mix_norm(st.junc[t], pr.at("junc." + M.tpl[t].name), 20.0, 1e-5);
    par["begin"] = mix_norm(st.begin, pr.at("begin"), 30.0, 1e-6);
    par["end"] = mix_norm(st.end, pr.at("end"), 30.0, 1e-6);
    par["p0"] = {std::max((st.n0 + 1.0) / (st.nreads + 100.0), 1e-4)};
    std::vector<int> z5 = st.z5, z3 = st.z3;
    const double fp = st.fast_n >= 100 ? (st.fast_ok + 1) / (st.fast_n + 2) : pr.at("gate")[2];
    par["gate"] = {(double)detail::quant(z5, 0.995, 120) + 30, (double)detail::quant(z3, 0.995, 90) + 30, fp};
}

inline Params Calibrator::finalize(int max_iter, double tol,
                                  const std::function<void()>& boundary_poll) {
    using namespace detail;
    using clk = std::chrono::steady_clock;
    auto T0 = clk::now();
    Model& M = *M_;
    std::ostringstream rep;
    char b[512];
    M.par = M.prior;
    M.par["guard"] = {0};
    for (int ty = 0; ty < M.A; ++ty) {
        std::vector<double> nl = M.prior.at("null." + M.type_name(ty));
        const auto& qq = M.prior.at("q." + M.type_name(ty));
        for (int f = 0; f < NFEAT; ++f)
            if (qq[f] > 0 && null_bases_ > 0) nl[f] = std::max((shuf_[ty][f] + 0.5) / null_bases_, 1e-10);
        M.par["null." + M.type_name(ty)] = nl;
    }
    M.finalize();
    const size_t N = len_.size();
    snprintf(b, sizeof b, "calibration: %zu reads, %.0f bases (mean %.0f), %.2f events/read; shuffled null %.0f bases\n", N, bases_,
             N ? bases_ / N : 0.0, N ? (double)ev_.size() / N : 0.0, null_bases_);
    rep << b;
    std::vector<uint64_t> sig(N, 0);
    Stats st;
    bool converged = false;
    int it = 0;
    std::array<double, 6> kd{};
    for (it = 1; it <= max_iter; ++it) {
        if (boundary_poll) boundary_poll();
        auto ta = clk::now();
        st.init(M);
        size_t changed = 0;
        double ll = 0;
        kd.fill(0);
        for (size_t r = 0; r < N; ++r) {
            if (boundary_poll && r % 128 == 0) boundary_poll();
            const uint32_t o0 = r ? off_[r - 1] : 0, o1 = off_[r];
            const Event* ev = ev_.data() + o0;
            const int n = (int)(o1 - o0);
            const int k = viterbi(M, ev, n, len_[r], W_);
            ll += W_.best_score;
            const uint64_t s = path_sig(W_);
            if (s != sig[r]) ++changed;
            sig[r] = s;
            kd[std::min(k, 5)] += 1;
            accumulate(ev, n, len_[r], k, st, gate_[r]);
        }
        update(st);
        M.finalize();
        const double fr = N ? (double)changed / N : 0;
        snprintf(b, sizeof b, "  iter %2d  ll/read %8.3f  changed %.4f  k0 %.4f k1 %.4f k2 %.4f k3 %.4f k4 %.4f k5+ %.4f  (%.3fs)\n", it,
                 N ? ll / N : 0.0, fr, kd[0] / std::max<size_t>(N, 1), kd[1] / std::max<size_t>(N, 1), kd[2] / std::max<size_t>(N, 1),
                 kd[3] / std::max<size_t>(N, 1), kd[4] / std::max<size_t>(N, 1), kd[5] / std::max<size_t>(N, 1),
                 std::chrono::duration<double>(clk::now() - ta).count());
        rep << b;
        if (it > 1 && fr < tol) { converged = true; break; }
    }
    const double k0 = N ? kd[0] / N : 1.0;
    // reads carrying at least one exact (ED 0) copy of a long static element of the main templates
    size_t n_exact = 0;
    {
        std::vector<uint8_t> main_long(M.A, 0);
        for (int t : M.main_tpl)
            for (int o : M.tpl[t].obs)
                if (M.tpl[t].el[o].kind == EK_ANCHOR && !M.anch[M.tpl[t].el[o].anc].shrt) main_long[M.tpl[t].el[o].anc] = 1;
        for (size_t r = 0; r < N; ++r) {
            if (boundary_poll && r % 128 == 0) boundary_poll();
            const uint32_t o0 = r ? off_[r - 1] : 0, o1 = off_[r];
            for (uint32_t x = o0; x < o1; ++x)
                if (ev_[x].type < M.A && main_long[ev_[x].type] && ev_[x].feat == 0) { ++n_exact; break; }
        }
    }
    const double fexact = N ? (double)n_exact / N : 0.0;
    snprintf(b, sizeof b, "  reads with an exact long static element: %.1f%%; reads without a certified construct: %.1f%%\n", 100 * fexact, 100 * k0);
    rep << b;
    std::string why;
    if (N < 200) why += "too few calibration reads; ";
    if (fexact < 0.30) {
        snprintf(b, sizeof b, "only %.1f%% of reads carry an exact copy of a layout adapter (static sequences do not match the library); ", 100 * fexact);
        why += b;
    }
    if (k0 > 0.5) {
        snprintf(b, sizeof b, "%.1f%% of reads carry no certified construct (layout anchors absent); ", 100 * k0);
        why += b;
    }
    for (int t : M.main_tpl) {
        if (st.tbeg[t] < 200 || st.tend[t] < 200) continue;
        const double po = st.topen[t] / st.tbeg[t], pc = st.tclose[t] / st.tend[t];
        snprintf(b, sizeof b, "  template %s: %.0f constructs; opening slot observed at construct start %.3f, closing slot at construct end %.3f\n",
                 M.tpl[t].name.c_str(), st.tbeg[t], po, pc);
        rep << b;
        if (po < 0.30 || pc < 0.30) {
            snprintf(b, sizeof b, "template %s rarely shows its %s element (%.2f); ", M.tpl[t].name.c_str(), po < 0.30 ? "opening" : "closing", std::min(po, pc));
            why += b;
        }
    }
    matches_ = why.empty();
    snprintf(b, sizeof b, "EM: %d iteration(s), %s; fast-path agreement %.4f (n=%.0f); gate zones 5' %d / 3' %d\n", std::min(it, max_iter),
             converged ? "converged" : "NOT converged", st.fast_n > 0 ? st.fast_ok / st.fast_n : 0.0, st.fast_n,
             (int)M.par["gate"][0], (int)M.par["gate"][1]);
    rep << b;
    if (!converged && matches_) rep << "  note: EM did not converge within the iteration cap\n";
    for (int a = 0; a < M.A; ++a) {
        const auto& q = M.par["q." + M.anch[a].seq];
        double pre = 0, sp = 0;
        for (int x = 0; x < NPARTBIN; ++x) { pre += q[FEAT_PREFIX + x]; sp += q[FEAT_SEEDPART + x]; }
        snprintf(b, sizeof b, "  %-16s P(ED0..%d)=", M.anch[a].row.c_str(), std::min(M.anch[a].k, 6));
        rep << b;
        for (int d = 0; d <= std::min(M.anch[a].k, 6); ++d) { snprintf(b, sizeof b, "%s%.3f", d ? "," : "", q[d]); rep << b; }
        snprintf(b, sizeof b, "  read-end prefix %.3f  partial %.3f\n", pre, sp);
        rep << b;
    }
    // cDNA strand table: [mean sense score, log-odds]; written only when the guard passes and the sample is large enough
    {
        double tot = 0;
        for (double x : strand_cnt_) tot += x;
        if (matches_ && strand_reads_ >= 500 && tot >= 2e5) {
            std::vector<double> tab((size_t)STRAND_N + 1, 0.0);
            double num = 0;
            for (int kk = 0; kk < STRAND_N; ++kk) {
                int rk = 0;
                for (int x = 0, v = kk; x < STRAND_K; ++x, v >>= 2) rk = (rk << 2) | (3 - (v & 3));
                const double lo = std::log((strand_cnt_[kk] + 1.0) / (strand_cnt_[rk] + 1.0));
                tab[kk + 1] = lo;
                num += strand_cnt_[kk] * lo;
            }
            const double mu = num / tot;
            snprintf(b, sizeof b, "strand table: %zu single-strand reads, %.0f cDNA bases, mean sense score %.4f nats/position%s\n",
                     strand_reads_, strand_bases_, mu, mu > 0.01 ? "" : " (too weak: no strand-flip cuts)");
            rep << b;
            if (mu > 0.01) { tab[0] = mu; M.par["strand"] = tab; }
        } else {
            snprintf(b, sizeof b, "strand table: not written (%zu single-strand reads, %.0f k-mers)\n", strand_reads_, tot);
            rep << b;
        }
    }
    // D strand table: cached as counts; written only when the guard passes and STRAND_D_MIN_KMERS k-mers were counted
    {
        double tot = 0;
        for (uint32_t x : strand_cnt_d_) tot += x;
        if (matches_ && strand_reads_ >= 500 && tot >= STRAND_D_MIN_KMERS) {
            std::vector<double> tab((size_t)STRAND_D_N);
            double num = 0;
            for (int kk = 0; kk < STRAND_D_N; ++kk) {
                if (boundary_poll && kk % 4096 == 0) boundary_poll();
                int rk = 0;
                for (int x = 0, v = kk; x < STRAND_D_K; ++x, v >>= 2) rk = (rk << 2) | (3 - (v & 3));
                const double lo = std::log((strand_cnt_d_[(size_t)kk] + 1.0) / (strand_cnt_d_[(size_t)rk] + 1.0));
                tab[(size_t)kk] = strand_cnt_d_[(size_t)kk];
                num += strand_cnt_d_[(size_t)kk] * lo;
            }
            const double mu = num / tot;
            snprintf(b, sizeof b, "D strand table: %d-mers, %.0f k-mers, mean sense score %.4f nats/position%s\n", STRAND_D_K, tot, mu,
                     mu > 0.01 ? "" : " (too weak: D geometries get no clear direction)");
            rep << b;
            if (mu > 0.01) M.par["strand_d"] = tab;
        } else {
            snprintf(b, sizeof b, "D strand table: not written (%s; %zu single-strand reads, %.0f k-mers; needs 500 reads and %.0f k-mers)\n",
                     !matches_ ? "calibration guard failed" : "too few reads or k-mers", strand_reads_, tot, STRAND_D_MIN_KMERS);
            rep << b;
        }
    }
    Params out;
    if (!matches_) {
        rep << "CALIBRATION GUARD FAILED: " << why << "-> falling back to layout priors\n";
        M.par = M.prior;
        M.par["guard"] = {1};
        M.finalize();
    } else {
        rep << "calibration guard: PASS (layout matches the library)\n";
        M.finalize();
    }
    out = M.export_params();
    out.guard_failed = !matches_;
    snprintf(b, sizeof b, "calibration time (EM + guard): %.3f s\n", std::chrono::duration<double>(clk::now() - T0).count());
    rep << b;
    report_ = rep.str();
    if (boundary_poll) boundary_poll();
    return out;
}

// Debug: events and decoded path of one read.
inline std::string debug_read(const Model& M, const char* seq, int len, Scratch& S) {
    using namespace detail;
    std::ostringstream os;
    stage0(M, seq, len, S);
    os << "len " << len << "  seed clusters:";
    for (auto& c : S.sc)
        if (c.valid()) os << " " << M.anch[c.a].row << "[" << c.s << "," << c.e << ") n" << (int)c.n << "/x" << (int)c.nx;
    os << "\npolys:";
    for (auto& p : S.poly) os << " " << M.poly_base[p.q] << "[" << p.s << "," << p.e << ")";
    Result R;
    const bool g = gate(M, seq, len, S, R);
    os << "\ngate: " << (g ? "FAST" : "full") << "\n";
    extract_events(M, seq, len, S, !M.opt.windowed_myers);
    for (size_t i = 0; i < S.ev.size(); ++i) {
        const Event& e = S.ev[i];
        os << "  ev" << i << " " << M.type_label(e.type) << " [" << e.start << "," << e.end << ") feat " << (int)e.feat << " cls " << (int)e.cls
           << " ed " << e.ed << " fl " << (int)e.fl << "\n";
    }
    int k = viterbi(M, S.ev.data(), (int)S.ev.size(), len, S);
    auto path = [&](const char* what) {
        os << what << " k=" << k << " score " << S.best_score << " path:";
        for (size_t x = 0; x < S.path_ev.size(); ++x)
            os << (S.path_cross[x] ? " |" : " ") << tkey(M, S.path_st[x]) << "@ev" << S.path_ev[x];
        os << "\n";
    };
    path("viterbi");
    const size_t ne0 = S.ev.size();
    if (facing_partner_scan(M, seq, len, S)) {
        os << "facing partner scan: " << S.ev.size() - ne0 << " event(s) added\n";
        for (size_t i = 0; i < S.ev.size(); ++i) {
            const Event& e = S.ev[i];
            os << "  ev" << i << " " << M.type_label(e.type) << " [" << e.start << "," << e.end << ") ed " << e.ed << "\n";
        }
        k = viterbi(M, S.ev.data(), (int)S.ev.size(), len, S);
        path("re-decode");
    }
    if (k >= 1 && facing_rescue(M, S.ev.data(), (int)S.ev.size(), S)) {
        os << "facing-anchor rescue: forced junction(s)";
        for (const auto& fp : S.force) os << " ev" << fp.first << "|ev" << fp.second;
        os << "\n";
        k = viterbi(M, S.ev.data(), (int)S.ev.size(), len, S);
        S.force.clear();
        path("rescued");
    }
    if (k > 0) {
        forward_backward(M, S.ev.data(), (int)S.ev.size(), len, S, k >= 2);
        derive_result(M, seq, S.ev.data(), (int)S.ev.size(), len, S, k, true, R);
        for (size_t c = 0; c < S.cweak.size(); ++c) {
            os << "junction " << c << ": child windows to " << S.jhi[c] << " / from " << S.jlo[c];
            if (S.cgap[c] != INT32_MIN)
                os << "; pairing junction (partner-delimited" << ((S.cweak[c] & 4) ? " left" : "") << ((S.cweak[c] & 8) ? " right" : "")
                   << ", g " << S.cgap[c] << (S.cweak[c] & 2 ? ", strong pairing" : "")
                   << ((S.cweak[c] & 1) ? "): no strong pairing, cDNA-side anchors not on both sides (posterior capped, abstain)" : "): ok");
            os << "\n";
        }
        static const char* kind_name[6] = {"both adapters", "one adapter", "layout geometry", "duration-MAP midpoint", "fold-back centre", "strand flip"};
        for (size_t i = 0; i < R.cuts.size(); ++i)
            os << "cut " << R.cuts[i] << ": " << (i < R.cut_kind.size() && R.cut_kind[i] < 6 ? kind_name[R.cut_kind[i]] : "?")
               << (i < R.cut_aux.size() && R.cut_aux[i] > 0 ? (i < R.cut_kind.size() && R.cut_kind[i] == CUT_FOLDBACK ? " (arm Jaccard " : " (margin ") + std::to_string(R.cut_aux[i]) + ")" : "")
               << ", posterior " << (i < R.cut_post.size() ? R.cut_post[i] : 0.f)
               << ", boundary flags " << (i < R.cut_flags.size() ? (int)R.cut_flags[i] : 0) << "\n";
        os << "strand table: " << (M.has_strand() ? "present" : "absent") << "\n";
        os << "D strand table (" << STRAND_D_K << "-mers, D strand rule): " << (M.has_strand_d() ? "present" : "absent") << "\n";
    }
    return os.str();
}

}  // namespace concat_hmm

} // namespace rrad_rad_core
