// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <zlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// --------- helpers -----------------------------------------------------------

struct Block { int64_t s, e; }; // inclusive genomic coordinates

inline void parse_cigar_blocks(int64_t ref_start1b,
                               std::string_view cig,
                               std::vector<Block>& out) {
  out.clear();
  int64_t rpos = ref_start1b; // 1-based
  int64_t bs = -1, be = -1;
  auto flush = [&](){ if (bs >= 0) { out.push_back({bs, be}); bs = be = -1; } };
  
  uint32_t n = 0;
  for (char c : cig) {
    if (c >= '0' && c <= '9') { n = n*10 + (c - '0'); continue; }
    if (n == 0) { continue; }
    switch (c) {
    case 'M': case '=': case 'X':
      if (bs < 0) bs = rpos;
      be = rpos + n - 1;
      rpos += n;
      break;
    case 'D':
      // deletion: ref advances, keep block open
      rpos += n;
      break;
    case 'N':
      // splice: close block, advance ref
      flush();
      rpos += n;
      break;
    case 'I': case 'S': case 'H': case 'P':
      // read-only or padding; no ref advance
      break;
    default:
      break;
    }
    n = 0;
  }
  flush();
}

inline int32_t bin_of(int64_t pos1b, int64_t region_start1b, int32_t bin_bp) {
  if (pos1b < region_start1b) pos1b = region_start1b;
  int64_t zero = pos1b - region_start1b; // 0-based within region
  if (zero < 0) zero = 0;
  return static_cast<int32_t>(zero / bin_bp);
}

inline uint64_t pack_key(uint32_t i, uint32_t j) {
  return (static_cast<uint64_t>(i) << 32) | static_cast<uint64_t>(j);
}

struct Pair { uint64_t key; uint32_t cnt; };

// RLE reduce a sorted vector of keys into (key,count)
static inline void reduce_sorted_keys(const std::vector<uint64_t>& keys,
                                      std::vector<Pair>& out) {
  out.clear();
  if (keys.empty()) return;
  uint64_t cur = keys[0];
  uint32_t c = 1;
  for (size_t k = 1; k < keys.size(); ++k) {
    if (keys[k] == cur) { ++c; }
    else { out.push_back({cur, c}); cur = keys[k]; c = 1; }
  }
  out.push_back({cur, c});
}

// Merge two sorted (key,count) arrays into dst (also sorted)
static inline void merge_counts(const std::vector<Pair>& a,
                                const std::vector<Pair>& b,
                                std::vector<Pair>& dst) {
  dst.clear();
  dst.reserve(a.size() + b.size());
  size_t i = 0, j = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i].key < b[j].key) dst.push_back(a[i++]);
    else if (b[j].key < a[i].key) dst.push_back(b[j++]);
    else { // equal
      dst.push_back({a[i].key, a[i].cnt + b[j].cnt});
      ++i; ++j;
    }
  }
  while (i < a.size()) dst.push_back(a[i++]);
  while (j < b.size()) dst.push_back(b[j++]);
}

// Lightweight PAF tokenizer: we only need fields 5,6,8,9,12 and optional "cg:Z:"
struct PafLite {
  std::string_view tname;
  int64_t tstart = -1; // 0-based
  int64_t tend   = -1; // 0-based end
  char strand = '+';   // '+' or '-'
  int mapq = 0;
  std::string_view cg; // CIGAR if present (from cg:Z:)
};

// Return false if parse failed or cg tag missing
inline bool parse_paf_line(std::string& line, PafLite& out) {
  // PAF: qname qlen qstart qend strand tname tlen tstart tend nmatch blen mapq [tags...]
  // We want field5,6,8,9,12 and cg:Z:...
  int col = 0;
  size_t i = 0, n = line.size();
  out = PafLite(); // reset
  
  auto next_field = [&](size_t& p) -> std::pair<size_t,size_t> {
    if (p >= n) return {n,n};
    size_t s = p;
    while (p < n && line[p] != '\t' && line[p] != '\n' && line[p] != '\r') ++p;
    size_t e = p;
    if (p < n && line[p] == '\t') ++p;
    return {s, e};
  };
  
  size_t p = 0;
  std::vector<std::pair<size_t,size_t>> fields;
  fields.reserve(16);
  while (p < n) {
    auto [s,e] = next_field(p);
    fields.push_back({s,e});
    if (p >= n) break;
  }
  if (fields.size() < 12) return false;
  
  // Strand (5)
  if (fields[4].second <= fields[4].first) return false;
  out.strand = line[fields[4].first];
  
  // tname (6)
  out.tname = std::string_view(line.data() + fields[5].first,
                               fields[5].second - fields[5].first);
  
  // tstart (8), tend (9)
  out.tstart = std::strtoll(line.c_str() + fields[7].first, nullptr, 10);
  out.tend   = std::strtoll(line.c_str() + fields[8].first, nullptr, 10);
  
  // mapq (12)
  out.mapq = std::atoi(line.c_str() + fields[11].first);
  
  // find cg:Z:
  for (size_t k = 12; k < fields.size(); ++k) {
    size_t s = fields[k].first, e = fields[k].second;
    if (e - s >= 5 &&
        line[s+0]=='c' && line[s+1]=='g' && line[s+2]==':' &&
        line[s+3]=='Z' && line[s+4]==':') {
      out.cg = std::string_view(line.data() + s + 5, e - (s + 5));
      break;
    }
  }
  // Require cg (CIGAR) to detect splices; if absent, skip
  if (out.cg.empty()) return false;
  return true;
}

// ---- line readers for plain and gzip files ---------------------------------
struct ILineReader {
  virtual bool getline(std::string& out) = 0;
  virtual ~ILineReader() {}
};

struct IfsReader : ILineReader {
  std::ifstream f;
  explicit IfsReader(const std::string& path) : f(path.c_str()) {}
  bool getline(std::string& out) override {
    if (!std::getline(f, out)) return false;
    if (!out.empty() && out.back() == '\r') out.pop_back();
    return true;
  }
  bool ok() const { return f.good(); }
};

struct GzReader : ILineReader {
  gzFile f = nullptr;
  std::vector<char> buf;
  explicit GzReader(const std::string& path, size_t buf_sz = 1<<16)
    : f(gzopen(path.c_str(), "rb")), buf(buf_sz, 0) {}
  bool getline(std::string& out) override {
    out.clear();
    if (!f) return false;
    for (;;) {
      char* r = gzgets(f, buf.data(), (int)buf.size());
      if (!r) {
        // EOF or error; return any partial line collected
        return !out.empty();
      }
      size_t len = std::strlen(r);
      if (len == 0) continue;
      out.append(r, len);
      if (!out.empty() && (out.back() == '\n' || out.back() == '\r')) {
        while (!out.empty() && (out.back() == '\n' || out.back() == '\r')) out.pop_back();
        return true;
      }
      // if buffer filled with no newline, loop to append next chunk
      if (len < buf.size() - 1) return true; // safety: short read without newline
    }
  }
  bool ok() const { return f != nullptr; }
  ~GzReader() { if (f) gzclose(f); }
};

static inline bool ends_with(const std::string& s, const std::string& suf) {
  if (s.size() < suf.size()) return false;
  return std::equal(suf.rbegin(), suf.rend(), s.rbegin());
}


// --------- core --------------------------------------------------------------

/*
 * splice_contacts_from_paf
 *
 * @param paf_path       path to PAF (plain text; gunzip first if .gz)
 * @param target         process only this target/contig (e.g. "chr14"); "" = accept all
 * @param region_start   1-based inclusive (default 1)
 * @param region_end     1-based inclusive; 0 = infer from max bin seen
 * @param bin_bp         bin size in bp (default 250)
 * @param min_mapq       minimum MAPQ (default 0)
 * @param chunk_lines    lines per chunk before aggregation (default 200000)
 * @param symmetric      return dsCMatrix (upper-tri stored) if TRUE;
 *                       else return full dgCMatrix (symmetrized)
 * @param threads        number of OMP threads (<=0 => use OMP default)
 *
 * Returns: S4 sparse matrix from Matrix package.
 */
// [[Rcpp::export]]
SEXP splice_contacts_from_paf(std::string paf_path,
                              std::string target = "",
                              double region_start = 1.0,
                              double region_end = 0.0,
                              int bin_bp = 250,
                              int min_mapq = 0,
                              int chunk_lines = 200000,
                              bool symmetric = true,
                              int threads = 0,
                              bool verbose = false) {
  
#ifdef _OPENMP
  if (threads > 0) omp_set_num_threads(threads);
#else
  if (threads > 0 && verbose) Rcpp::Rcout << "OpenMP not available; running single-threaded.\n";
#endif

  if (bin_bp <= 0) stop("bin_bp must be > 0");
  if (region_start < 1.0) region_start = 1.0;
  
  const int64_t R1 = static_cast<int64_t>(region_start);
  int64_t R2 = (region_end <= 0.0) ? 0 : static_cast<int64_t>(region_end);
  
  std::ifstream fin(paf_path.c_str());
  if (!fin) stop("Failed to open PAF: %s", paf_path);
  
  std::vector<std::string> lines;
  lines.reserve(chunk_lines);
  
  std::vector<Pair> globalAgg;        // sorted (key,count) accumulator
  std::vector<Pair> mergedBuf;        // tmp buffer for merges
  uint32_t max_bin_seen = 0;
  
  auto process_chunk = [&](std::vector<std::string>& chunk){
    if (chunk.empty()) return;
    
    // Thread-local edge keys
    std::vector<std::vector<uint64_t>> tls_keys;
#ifdef _OPENMP
    int nt = 1;
#pragma omp parallel
{ #pragma omp single nowait
{ nt = omp_get_num_threads(); }
}
#else
int nt = 1;
#endif
tls_keys.resize(nt);

#ifdef _OPENMP
#pragma omp parallel
#endif
{
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  auto& keys = tls_keys[tid];
  keys.reserve(chunk.size() / (nt*2) + 1024);
  std::vector<Block> blocks;
  blocks.reserve(8);
  
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for (int i = 0; i < (int)chunk.size(); ++i) {
    PafLite pl;
    if (!parse_paf_line(chunk[i], pl)) continue;
    if (min_mapq > 0 && pl.mapq < min_mapq) continue;
    if (!target.empty() && pl.tname != target) continue;
    
    // Convert PAF tstart to 1-based for our block parser
    int64_t aln_start1b = pl.tstart + 1;
    
    parse_cigar_blocks(aln_start1b, pl.cg, blocks);
    if (blocks.size() < 2) continue;
    
    const bool neg = (pl.strand == '-');
    
    for (size_t k = 0; k + 1 < blocks.size(); ++k) {
      const auto& A = blocks[k];
      const auto& B = blocks[k+1];
      
      // Keep endpoints within region if region end is set
      if (R2 > 0) {
        if (A.e   < R1 || A.e   > R2) continue;
        if (B.s   < R1 || B.s   > R2) continue;
      }
      
      uint32_t bi = static_cast<uint32_t>(bin_of(A.e, R1, bin_bp));
      uint32_t bj = static_cast<uint32_t>(bin_of(B.s, R1, bin_bp));
      if (bi > bj) std::swap(bi, bj); // upper triangle
      
      keys.push_back(pack_key(bi, bj));
    }
  }
} // end parallel

// Merge TLS keys → single vector
size_t total = 0;
for (auto& v : tls_keys) total += v.size();
std::vector<uint64_t> all;
all.reserve(total);
for (auto& v : tls_keys) {
  all.insert(all.end(), v.begin(), v.end());
  v.clear(); v.shrink_to_fit();
}
if (all.empty()) { chunk.clear(); return; }

// Update max bin
for (auto k : all) {
  uint32_t j = static_cast<uint32_t>(k & 0xffffffffu);
  if (j > max_bin_seen) max_bin_seen = j;
}

// Sort & reduce within chunk
std::sort(all.begin(), all.end());
std::vector<Pair> localAgg;
localAgg.reserve(all.size() / 2 + 1);
reduce_sorted_keys(all, localAgg);

// Merge into global
if (globalAgg.empty()) {
  globalAgg.swap(localAgg);
} else {
  merge_counts(globalAgg, localAgg, mergedBuf);
  globalAgg.swap(mergedBuf);
  localAgg.clear(); mergedBuf.clear();
}

all.clear(); all.shrink_to_fit();
chunk.clear();
  };
  
  // Streaming read
  std::string line;
  while (std::getline(fin, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    lines.emplace_back(std::move(line));
    if ((int)lines.size() >= chunk_lines) process_chunk(lines);
  }
  if (!lines.empty()) process_chunk(lines);
  
  fin.close();
  
  if (globalAgg.empty()) {
    stop("No edges produced. Did your PAF contain cg:Z: CIGARs and match the target/region?");
  }
  
  // Determine matrix size (bins)
  uint32_t n_bins;
  if (R2 > 0) {
    int64_t region_len = (R2 - R1 + 1);
    n_bins = static_cast<uint32_t>((region_len + bin_bp - 1) / bin_bp);
  } else {
    n_bins = max_bin_seen + 1;
  }
  
  // Convert (key,count) → column-compressed format
  // We want by-column order: sort by (j,i)
  struct Entry { uint32_t j, i; double x; };
  std::vector<Entry> es; es.reserve(globalAgg.size());
  
  for (const auto& pr : globalAgg) {
    uint32_t i = static_cast<uint32_t>(pr.key >> 32);
    uint32_t j = static_cast<uint32_t>(pr.key & 0xffffffffu);
    es.push_back({j, i, static_cast<double>(pr.cnt)});
  }
  globalAgg.clear(); globalAgg.shrink_to_fit();
  
  std::sort(es.begin(), es.end(), [](const Entry& a, const Entry& b){
    if (a.j != b.j) return a.j < b.j;
    return a.i < b.i;
  });
  
  // Build dsCMatrix (upper triangle only) or full dgCMatrix (symmetrized)
  if (symmetric) {
    // dsCMatrix (upper triangle stored; uplo='U')
    // Construct slots: i (row idx), p (column ptr), x (values)
    std::vector<int>    i; i.reserve(es.size());
    std::vector<int>    p(n_bins + 1, 0);
    std::vector<double> x; x.reserve(es.size());
    
    uint32_t cur_j = 0;
    size_t k = 0;
    for (uint32_t j = 0; j < n_bins; ++j) {
      p[j] = static_cast<int>(i.size());
      while (k < es.size() && es[k].j == j) {
        i.push_back(static_cast<int>(es[k].i));
        x.push_back(es[k].x);
        ++k;
      }
    }
    p[n_bins] = static_cast<int>(i.size());
    
    // Create S4 dsCMatrix
    Environment Matrix("package:Matrix");
    S4 mat("dsCMatrix");
    mat.slot("i")    = IntegerVector(i.begin(), i.end());
    mat.slot("p")    = IntegerVector(p.begin(), p.end());
    mat.slot("x")    = NumericVector(x.begin(), x.end());
    mat.slot("Dim")  = IntegerVector::create(n_bins, n_bins);
    mat.slot("Dimnames") = List::create(R_NilValue, R_NilValue);
    mat.slot("uplo") = CharacterVector::create("U");
    mat.slot("factors") = R_NilValue;
    
    return mat;
  } else {
    // Build full symmetric dgCMatrix (store both (i,j) and (j,i), diag once)
    // First count nnz per column for CSC sizes
    // We'll produce entries for (j,i) from es AND mirrored (i,j) where i!=j.
    std::vector<int> col_nnz(n_bins, 0);
    // count originals
    for (const auto& e : es) col_nnz[e.j]++;
    // count mirrors
    for (const auto& e : es) if (e.i != e.j) col_nnz[e.i]++;
    
    std::vector<int> p(n_bins + 1, 0);
    for (uint32_t j = 0; j < n_bins; ++j) p[j+1] = p[j] + col_nnz[j];
    const int nnz = p[n_bins];
    
    std::vector<int>    i(nnz);
    std::vector<double> x(nnz);
    
    // Fill columns
    std::vector<int> cursor = p; // per-column write cursor
    // originals
    for (const auto& e : es) {
      int idx = cursor[e.j]++;
      i[idx] = static_cast<int>(e.i);
      x[idx] = e.x;
    }
    // mirrors
    for (const auto& e : es) if (e.i != e.j) {
      int idx = cursor[e.i]++;
      i[idx] = static_cast<int>(e.j);
      x[idx] = e.x;
    }
    
    Environment Matrix("package:Matrix");
    S4 mat("dgCMatrix");
    mat.slot("i")    = IntegerVector(i.begin(), i.end());
    mat.slot("p")    = IntegerVector(p.begin(), p.end());
    mat.slot("x")    = NumericVector(x.begin(), x.end());
    mat.slot("Dim")  = IntegerVector::create(n_bins, n_bins);
    mat.slot("Dimnames") = List::create(R_NilValue, R_NilValue);
    mat.slot("factors") = R_NilValue;
    
    return mat;
  }
}


// fetch a column by name, tolerant to integer/double
static inline NumericVector get_num_col(DataFrame df, const std::string& nm) {
  SEXP col = df[nm];
  if (Rf_isNull(col)) stop("Column '%s' not found", nm.c_str());
  if (TYPEOF(col) == INTSXP) return as<NumericVector>(IntegerVector(col));
  if (TYPEOF(col) == REALSXP) return NumericVector(col);
  stop("Column '%s' must be integer or numeric", nm.c_str());
  return NumericVector();
}
static inline CharacterVector get_chr_col(DataFrame df, const std::string& nm) {
  SEXP col = df[nm];
  if (Rf_isNull(col)) stop("Column '%s' not found", nm.c_str());
  if (TYPEOF(col) == STRSXP) return CharacterVector(col);
  stop("Column '%s' must be character", nm.c_str());
  return CharacterVector();
}

// strip leading "cg:Z:" if present
static inline std::string_view normalize_cg(std::string_view cg) {
  if (cg.size() >= 5 && cg.substr(0,5) == std::string_view("cg:Z:"))
    return cg.substr(5);
  return cg;
}

// overlap within the donor's anchor bin (bin containing A.e)
inline int anchor_overlap_bp_donor(const Block& a, int64_t R1, int32_t bin_bp) {
  uint32_t bi = static_cast<uint32_t>(bin_of(a.e, R1, bin_bp));
  int64_t bs = R1 + static_cast<int64_t>(bi) * bin_bp;
  int64_t be = bs + bin_bp - 1;
  int64_t ov = std::min<int64_t>(a.e, be) - std::max<int64_t>(a.s, bs) + 1;
  return (int)std::max<int64_t>(0, ov);
}
// overlap within the acceptor's anchor bin (bin containing B.s)
inline int anchor_overlap_bp_acceptor(const Block& b, int64_t R1, int32_t bin_bp) {
  uint32_t bj = static_cast<uint32_t>(bin_of(b.s, R1, bin_bp));
  int64_t bs = R1 + static_cast<int64_t>(bj) * bin_bp;
  int64_t be = bs + bin_bp - 1;
  int64_t ov = std::min<int64_t>(b.e, be) - std::max<int64_t>(b.s, bs) + 1;
  return (int)std::max<int64_t>(0, ov);
}


// [[Rcpp::export]]
SEXP splice_contacts_from_paf_dt(
    DataFrame df,
    std::string target = "",
    double region_start = 1.0,
    double region_end   = 0.0,
    int     bin_bp      = 250,
    int     min_mapq    = 0,
    int     min_anchor_bp = 1,   // <-- NEW
    int     chunk_rows  = 1000000,
    bool    symmetric   = true,
    int     threads     = 0,
    bool    verbose     = false,
    std::string col_tname = "tname",
    std::string col_tstart= "tstart",
    std::string col_tend  = "tend",
    std::string col_cg    = "cg",
    std::string col_mapq  = "mapq",
    std::string col_strand= "strand"
) {
#ifdef _OPENMP
  if (threads > 0) omp_set_num_threads(threads);
#else
  if (threads > 0 && verbose) Rcpp::Rcout << "OpenMP not available; running single-threaded.\n";
#endif
  
  if (bin_bp <= 0) stop("bin_bp must be > 0");
  if (region_start < 1.0) region_start = 1.0;
  
  const int64_t R1 = static_cast<int64_t>(region_start);
  int64_t R2 = (region_end <= 0.0) ? 0 : static_cast<int64_t>(region_end);
  
  // columns
  CharacterVector v_tname = get_chr_col(df, col_tname);
  NumericVector   v_tstart= get_num_col(df, col_tstart);
  // tend is not required but try to fetch if present
  bool have_tend = !Rf_isNull(df[col_tend]);
  NumericVector   v_tend  = have_tend ? get_num_col(df, col_tend) : NumericVector(df.nrows(), NA_REAL);
  CharacterVector v_cg    = get_chr_col(df, col_cg);
  
  bool have_mapq = !Rf_isNull(df[col_mapq]);
  NumericVector   v_mapq  = have_mapq ? get_num_col(df, col_mapq) : NumericVector(df.nrows(), 60.0);
  
  const int64_t n = df.nrows();
  if (n == 0) stop("Empty data.frame");
  
  // global accumulator (sorted key,count)
  std::vector<Pair> globalAgg;
  std::vector<Pair> mergeBuf;
  uint32_t max_bin_seen = 0;
  
  auto process_range = [&](int64_t a, int64_t b){
    if (a >= b) return;
    // thread-local key buffers
    std::vector<std::vector<uint64_t>> tls_keys;
#ifdef _OPENMP
    int nt = 1;
#pragma omp parallel
    { #pragma omp single nowait { nt = omp_get_num_threads(); } }
#else
    int nt = 1;
#endif
    tls_keys.resize(nt);
    
#ifdef _OPENMP
#pragma omp parallel
#endif
{
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  auto& keys = tls_keys[tid];
  keys.reserve((b - a) / (nt*2) + 1024);
  std::vector<Block> blocks; blocks.reserve(8);
  
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for (int64_t i = a; i < b; ++i) {
    // filter by target if provided
    if (!target.empty()) {
      SEXP s = v_tname[i];
      if (s == NA_STRING) continue;
      if (std::string_view(CHAR(s)) != std::string_view(target)) continue;
    }
    // mapq filter
    if (have_mapq && min_mapq > 0) {
      double mq = v_mapq[i];
      if (!Rcpp::NumericVector::is_na(mq) && mq < min_mapq) continue;
    }
    
    // cg
    SEXP scg = v_cg[i];
    if (scg == NA_STRING) continue;
    std::string_view cgsv = normalize_cg(std::string_view(CHAR(scg)));
    
    // tstart (PAF 0-based) -> 1-based for our block parser
    double dts = v_tstart[i];
    if (Rcpp::NumericVector::is_na(dts)) continue;
    int64_t tstart1b = static_cast<int64_t>(dts) + 1;
    
    // build blocks
    parse_cigar_blocks(tstart1b, cgsv, blocks);
    if (blocks.size() < 2) continue;
    for (size_t k = 0; k + 1 < blocks.size(); ++k) {
      const auto& A = blocks[k];
      const auto& B = blocks[k+1];
      
      if (R2 > 0) {
        if (A.e < R1 || A.e > R2) continue;
        if (B.s < R1 || B.s > R2) continue;
      }
      
      // --- NEW: require minimum anchor overlap in their anchor bins ---
      if (min_anchor_bp > 1) {
        int l_ov = anchor_overlap_bp_donor(A, R1, bin_bp);
        int r_ov = anchor_overlap_bp_acceptor(B, R1, bin_bp);
        if (l_ov < min_anchor_bp || r_ov < min_anchor_bp) continue;
      }
      
      uint32_t bi = static_cast<uint32_t>(bin_of(A.e, R1, bin_bp));
      uint32_t bj = static_cast<uint32_t>(bin_of(B.s, R1, bin_bp));
      if (bi > bj) std::swap(bi, bj); // keep upper triangle
      keys.push_back(pack_key(bi, bj));
    }
  }
} // end parallel

// merge TLS -> one vector
size_t total = 0;
for (auto& v : tls_keys) total += v.size();
std::vector<uint64_t> all;
all.reserve(total);
for (auto& v : tls_keys) {
  all.insert(all.end(), v.begin(), v.end());
  v.clear();
  v.shrink_to_fit();
}
if (all.empty()) return;

// track max bin for dimension inference
for (auto k : all) {
  uint32_t j = static_cast<uint32_t>(k & 0xffffffffu);
  if (j > max_bin_seen) max_bin_seen = j;
}

// sort+reduce chunk
std::sort(all.begin(), all.end());
std::vector<Pair> localAgg;
localAgg.reserve(all.size() / 2 + 1);
reduce_sorted_keys(all, localAgg);

if (globalAgg.empty()) {
  globalAgg.swap(localAgg);
} else {
  merge_counts(globalAgg, localAgg, mergeBuf);
  globalAgg.swap(mergeBuf);
  localAgg.clear(); mergeBuf.clear();
}
  };
  
  // chunked pass to bound memory
  if (chunk_rows <= 0) chunk_rows = 1000000;
  for (int64_t a = 0; a < n; a += chunk_rows) {
    int64_t b = std::min<int64_t>(n, a + chunk_rows);
    if (verbose) Rcpp::Rcout << "Processing rows [" << a << "," << b << ")\n";
    process_range(a, b);
  }
  
  if (globalAgg.empty()) stop("No edges produced (check cg/mapq/target filters).");
  
  // matrix dimension (bins)
  uint32_t n_bins;
  if (R2 > 0) {
    int64_t region_len = (R2 - R1 + 1);
    n_bins = static_cast<uint32_t>((region_len + bin_bp - 1) / bin_bp);
  } else {
    n_bins = max_bin_seen + 1;
  }
  
  // convert to CSC (column = j)
  struct Entry { uint32_t j, i; double x; };
  std::vector<Entry> es; es.reserve(globalAgg.size());
  for (const auto& pr : globalAgg) {
    uint32_t i = static_cast<uint32_t>(pr.key >> 32);
    uint32_t j = static_cast<uint32_t>(pr.key & 0xffffffffu);
    es.push_back({j, i, static_cast<double>(pr.cnt)});
  }
  globalAgg.clear(); globalAgg.shrink_to_fit();
  
  std::sort(es.begin(), es.end(), [](const Entry& a, const Entry& b){
    if (a.j != b.j) return a.j < b.j;
    return a.i < b.i;
  });
  
  if (symmetric) {
    // dsCMatrix (upper triangle, uplo='U')
    std::vector<int>    i; i.reserve(es.size());
    std::vector<int>    p(n_bins + 1, 0);
    std::vector<double> x; x.reserve(es.size());
    
    size_t k = 0;
    for (uint32_t j = 0; j < n_bins; ++j) {
      p[j] = static_cast<int>(i.size());
      while (k < es.size() && es[k].j == j) {
        i.push_back(static_cast<int>(es[k].i));
        x.push_back(es[k].x);
        ++k;
      }
    }
    p[n_bins] = static_cast<int>(i.size());
    
    S4 mat("dsCMatrix");
    mat.slot("i")    = IntegerVector(i.begin(), i.end());
    mat.slot("p")    = IntegerVector(p.begin(), p.end());
    mat.slot("x")    = NumericVector(x.begin(), x.end());
    mat.slot("Dim")  = IntegerVector::create(n_bins, n_bins);
    mat.slot("Dimnames") = List::create(R_NilValue, R_NilValue);
    mat.slot("uplo") = CharacterVector::create("U");
    mat.slot("factors") = Rcpp::List::create();
    
    return mat;
  } else {
    // full symmetric dgCMatrix (store both (i,j) and (j,i))
    std::vector<int> col_nnz(n_bins, 0);
    for (const auto& e : es) col_nnz[e.j]++;
    for (const auto& e : es) if (e.i != e.j) col_nnz[e.i]++;
    
    std::vector<int> p(n_bins + 1, 0);
    for (uint32_t j = 0; j < n_bins; ++j) p[j+1] = p[j] + col_nnz[j];
    const int nnz = p[n_bins];
    
    std::vector<int>    i(nnz);
    std::vector<double> x(nnz);
    std::vector<int> cursor = p;
    
    for (const auto& e : es) {
      int idx = cursor[e.j]++; i[idx] = static_cast<int>(e.i); x[idx] = e.x;
    }
    for (const auto& e : es) if (e.i != e.j) {
      int idx = cursor[e.i]++; i[idx] = static_cast<int>(e.j); x[idx] = e.x;
    }
    
    S4 mat("dgCMatrix");
    mat.slot("i")    = IntegerVector(i.begin(), i.end());
    mat.slot("p")    = IntegerVector(p.begin(), p.end());
    mat.slot("x")    = NumericVector(x.begin(), x.end());
    mat.slot("Dim")  = IntegerVector::create(n_bins, n_bins);
    mat.slot("Dimnames") = List::create(R_NilValue, R_NilValue);
    mat.slot("factors") = Rcpp::List::create();
    
    return mat;
  }
}