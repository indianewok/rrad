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
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <zlib.h>   // <-- gzip streaming

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// --------------------- helpers shared with *_dt --------------------------------

struct Block { int64_t s, e; }; // inclusive genomic coords

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
      rpos += n; // ref advances, keep block open
      break;
    case 'N':
      flush();   // splice: close block, advance ref
      rpos += n;
      break;
    case 'I': case 'S': case 'H': case 'P':
      break;     // read-only or padding; no ref advance
    default: break;
    }
    n = 0;
  }
  flush();
}

inline int32_t bin_of(int64_t pos1b, int64_t region_start1b, int32_t bin_bp) {
  if (pos1b < region_start1b) pos1b = region_start1b;
  int64_t zero = pos1b - region_start1b;
  if (zero < 0) zero = 0;
  return static_cast<int32_t>(zero / bin_bp);
}

inline uint64_t pack_key(uint32_t i, uint32_t j) {
  return (static_cast<uint64_t>(i) << 32) | static_cast<uint64_t>(j);
}

struct Pair { uint64_t key; uint32_t cnt; };

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

static inline void merge_counts(const std::vector<Pair>& a,
                                const std::vector<Pair>& b,
                                std::vector<Pair>& dst) {
  dst.clear();
  dst.reserve(a.size() + b.size());
  size_t i = 0, j = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i].key < b[j].key) dst.push_back(a[i++]);
    else if (b[j].key < a[i].key) dst.push_back(b[j++]);
    else { dst.push_back({a[i].key, a[i].cnt + b[j].cnt}); ++i; ++j; }
  }
  while (i < a.size()) dst.push_back(a[i++]);
  while (j < b.size()) dst.push_back(b[j++]);
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

// ------------------------- PAF parsing ----------------------------------------

struct PafLite {
  std::string_view tname;
  int64_t tstart = -1; // 0-based
  int64_t tend   = -1; // 0-based end
  char strand = '+';   // '+' or '-'
  int mapq = 0;
  std::string_view cg; // optional cg:Z: CIGAR
};

inline bool parse_paf_line(std::string& line, PafLite& out) {
  // PAF fields: qname qlen qstart qend strand tname tlen tstart tend nmatch blen mapq [tags...]
  size_t n = line.size();
  size_t p = 0;
  auto next_field = [&](size_t& pos)->std::pair<size_t,size_t>{
    if (pos >= n) return {n,n};
    size_t s = pos;
    while (pos < n && line[pos] != '\t' && line[pos] != '\n' && line[pos] != '\r') ++pos;
    size_t e = pos;
    if (pos < n && line[pos] == '\t') ++pos;
    return {s,e};
  };
  
  std::vector<std::pair<size_t,size_t>> f;
  f.reserve(16);
  while (p < n) {
    auto se = next_field(p);
    if (se.first == se.second && se.first == n) break;
    f.push_back(se);
  }
  if (f.size() < 12) return false;
  
  out = PafLite();
  // strand (5th, 0-based idx 4)
  if (f[4].second <= f[4].first) return false;
  out.strand = line[f[4].first];
  
  // tname (6th)
  out.tname = std::string_view(line.data() + f[5].first, f[5].second - f[5].first);
  
  // tstart (8th), tend (9th)
  out.tstart = std::strtoll(line.c_str() + f[7].first, nullptr, 10);
  out.tend   = std::strtoll(line.c_str() + f[8].first, nullptr, 10);
  
  // mapq (12th)
  out.mapq = std::atoi(line.c_str() + f[11].first);
  
  // find cg:Z:
  for (size_t k = 12; k < f.size(); ++k) {
    size_t s = f[k].first, e = f[k].second;
    if (e > s + 5 &&
        line[s+0]=='c' && line[s+1]=='g' && line[s+2]==':' &&
        line[s+3]=='Z' && line[s+4]==':') {
      out.cg = std::string_view(line.data() + s + 5, e - (s + 5));
      break;
    }
  }
  if (out.cg.empty()) return false; // need CIGAR with 'N' to detect splices
  return true;
}

// --------------------- gzip/plain line readers --------------------------------

struct ILineReader { virtual bool getline(std::string& out)=0; virtual ~ILineReader(){} };

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
      if (!r) return !out.empty();       // EOF or error; return any partial
      size_t len = std::strlen(r);
      if (len) out.append(r, len);
      if (!out.empty() && (out.back()=='\n' || out.back()=='\r')) {
        while (!out.empty() && (out.back()=='\n' || out.back()=='\r')) out.pop_back();
        return true;
      }
      if (len < buf.size()-1) return true; // short read without newline
    }
  }
  bool ok() const { return f != nullptr; }
  ~GzReader(){ if (f) gzclose(f); }
};

static inline bool ends_with(const std::string& s, const std::string& suf) {
  return s.size() >= suf.size() && std::equal(suf.rbegin(), suf.rend(), s.rbegin());
}

// --------------------------- MAIN: file streamer ------------------------------

/*
 * splice_contacts_from_paf
 *
 * Stream a .paf or .paf.gz, count splice junction bin pairs, and return a sparse matrix.
 *
 * @param path     path to PAF (.paf or .paf.gz)
 * @param target       contig to keep ("" = all)
 * @param region_start 1-based inclusive (default 1)
 * @param region_end   1-based inclusive; 0 = infer from max bin seen
 * @param bin_bp       bin size in bp (default 250)
 * @param min_mapq     minimum MAPQ (default 0)
 * @param chunk_lines  lines per chunk before aggregation (default 200000)
 * @param symmetric    TRUE -> dsCMatrix (upper triangle stored, uplo='U'); FALSE -> full dgCMatrix
 * @param threads      OpenMP threads (<=0 use OMP default)
 * @param verbose      print chunk progress
 */
// [[Rcpp::export]]
SEXP stream_paf2mat(std::string path,
                              std::string target = "",
                              double region_start = 1.0,
                              double region_end   = 0.0,
                              int     bin_bp      = 250,
                              int     min_mapq    = 0,
                              int     chunk_lines = 200000,
                              bool    symmetric   = true,
                              int     threads     = 0,
                              bool    verbose     = false) {
#ifdef _OPENMP
  if (threads > 0) omp_set_num_threads(threads);
#else
  if (threads > 0 && verbose) Rcpp::Rcout << "OpenMP not available; running single-threaded.\n";
#endif
  if (bin_bp <= 0) stop("bin_bp must be > 0");
  if (region_start < 1.0) region_start = 1.0;
  
  const int64_t R1 = static_cast<int64_t>(region_start);
  int64_t R2 = (region_end <= 0.0) ? 0 : static_cast<int64_t>(region_end);
  
  // open reader
  std::unique_ptr<ILineReader> reader;
  if (ends_with(path, ".gz")) {
    auto gz = std::make_unique<GzReader>(path);
    if (!gz->ok()) stop("Failed to open gzip PAF: " + path);
    reader = std::move(gz);
  } else {
    auto ifs = std::make_unique<IfsReader>(path);
    if (!ifs->ok()) stop("Failed to open PAF: " + path);
    reader = std::move(ifs);
  }
  
  std::vector<Pair> globalAgg;
  std::vector<Pair> mergedBuf;
  uint32_t max_bin_seen = 0;
  
  auto process_chunk = [&](std::vector<std::string>& chunk){
    if (chunk.empty()) return;
    
    // per-thread key buffers
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
  keys.reserve(chunk.size() / (nt*2) + 1024);
  std::vector<Block> blocks; blocks.reserve(8);
  
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for (int i = 0; i < (int)chunk.size(); ++i) {
    PafLite pl;
    if (!parse_paf_line(chunk[i], pl)) continue;
    if (min_mapq > 0 && pl.mapq < min_mapq) continue;
    if (!target.empty() && pl.tname != target) continue;
    
    int64_t tstart1b = pl.tstart + 1; // PAF is 0-based
    parse_cigar_blocks(tstart1b, pl.cg, blocks);
    if (blocks.size() < 2) continue;
    
    for (size_t k = 0; k + 1 < blocks.size(); ++k) {
      const auto& A = blocks[k];
      const auto& B = blocks[k+1];
      
      if (R2 > 0) {
        if (A.e < R1 || A.e > R2) continue;
        if (B.s < R1 || B.s > R2) continue;
      }
      
      uint32_t bi = static_cast<uint32_t>(bin_of(A.e, R1, bin_bp));
      uint32_t bj = static_cast<uint32_t>(bin_of(B.s, R1, bin_bp));
      if (bi > bj) std::swap(bi, bj); // keep upper triangle
      keys.push_back(pack_key(bi, bj));
    }
  }
} // end parallel

// merge tls → single vector
size_t total = 0; for (auto& v : tls_keys) total += v.size();
std::vector<uint64_t> all; all.reserve(total);
for (auto& v : tls_keys) { all.insert(all.end(), v.begin(), v.end()); v.clear(); v.shrink_to_fit(); }
if (all.empty()) { chunk.clear(); return; }

// track max bin
for (auto k : all) { uint32_t j = static_cast<uint32_t>(k & 0xffffffffu); if (j > max_bin_seen) max_bin_seen = j; }

// sort + reduce within chunk
std::sort(all.begin(), all.end());
std::vector<Pair> localAgg; localAgg.reserve(all.size() / 2 + 1);
reduce_sorted_keys(all, localAgg);

// merge into global
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
  
  // streaming read loop
  std::vector<std::string> lines; lines.reserve(chunk_lines);
  std::string line;
  while (reader->getline(line)) {
    lines.emplace_back(std::move(line));
    if ((int)lines.size() >= chunk_lines) {
      if (verbose) Rcout << "Processing " << lines.size() << " lines...\n";
      process_chunk(lines);
    }
  }
  if (!lines.empty()) {
    if (verbose) Rcout << "Processing " << lines.size() << " lines...\n";
    process_chunk(lines);
  }
  
  if (globalAgg.empty()) stop("No edges produced (check cg/mapq/target filters).");
  
  // decide matrix dimension
  uint32_t n_bins = (R2 > 0)
    ? static_cast<uint32_t>(((R2 - R1 + 1) + bin_bp - 1) / bin_bp)
      : (max_bin_seen + 1);
  
  // convert to column-sorted entries (col=j)
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
    // dsCMatrix (upper triangle stored; uplo='U')
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
    mat.slot("factors") = Rcpp::List::create();   // must be a list (not NULL)
    return mat;
    
  } else {
    // full symmetric dgCMatrix
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

// [[Rcpp::export]]
SEXP stream_paf2mat_anch(std::string path,
                    std::string target = "",
                    double region_start = 1.0,
                    double region_end   = 0.0,
                    int     bin_bp      = 250,
                    int     min_mapq    = 0,
                    int     min_anchor_bp = 1,   // <-- NEW
                    int     chunk_lines = 200000,
                    bool    symmetric   = true,
                    int     threads     = 0,
                    bool    verbose     = false) {
#ifdef _OPENMP
  if (threads > 0) omp_set_num_threads(threads);
#else
  if (threads > 0 && verbose) Rcpp::Rcout << "OpenMP not available; running single-threaded.\n";
#endif
  if (bin_bp <= 0) stop("bin_bp must be > 0");
  if (region_start < 1.0) region_start = 1.0;
  
  const int64_t R1 = static_cast<int64_t>(region_start);
  int64_t R2 = (region_end <= 0.0) ? 0 : static_cast<int64_t>(region_end);
  
  // open reader (ILineReader/GzReader/IfsReader assumed defined elsewhere)
  std::unique_ptr<ILineReader> reader;
  if (ends_with(path, ".gz")) {
    auto gz = std::make_unique<GzReader>(path);
    if (!gz->ok()) stop("Failed to open gzip PAF: " + path);
    reader = std::move(gz);
  } else {
    auto ifs = std::make_unique<IfsReader>(path);
    if (!ifs->ok()) stop("Failed to open PAF: " + path);
    reader = std::move(ifs);
  }
  
  std::vector<Pair> globalAgg, mergedBuf;
  uint32_t max_bin_seen = 0;
  
  auto process_chunk = [&](std::vector<std::string>& chunk){
    if (chunk.empty()) return;
    
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
  keys.reserve(chunk.size() / (nt*2) + 1024);
  std::vector<Block> blocks; blocks.reserve(8);
  
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for (int i = 0; i < (int)chunk.size(); ++i) {
    PafLite pl;
    if (!parse_paf_line(chunk[i], pl)) continue;
    if (min_mapq > 0 && pl.mapq < min_mapq) continue;
    if (!target.empty() && pl.tname != target) continue;
    
    int64_t tstart1b = pl.tstart + 1; // PAF is 0-based
    parse_cigar_blocks(tstart1b, pl.cg, blocks);
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

// merge tls → single vector
size_t total = 0; for (auto& v : tls_keys) total += v.size();
std::vector<uint64_t> all; all.reserve(total);
for (auto& v : tls_keys) { all.insert(all.end(), v.begin(), v.end()); v.clear(); v.shrink_to_fit(); }
if (all.empty()) { chunk.clear(); return; }

// track max bin
for (auto k : all) {
  uint32_t j = static_cast<uint32_t>(k & 0xffffffffu);
  if (j > max_bin_seen) max_bin_seen = j;
}

// sort + reduce within chunk
std::sort(all.begin(), all.end());
std::vector<Pair> localAgg; localAgg.reserve(all.size() / 2 + 1);
reduce_sorted_keys(all, localAgg);

// merge into global
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
  
  // streaming read loop
  std::vector<std::string> lines; lines.reserve(chunk_lines);
  std::string line;
  while (reader->getline(line)) {
    lines.emplace_back(std::move(line));
    if ((int)lines.size() >= chunk_lines) {
      if (verbose) Rcout << "Processing " << lines.size() << " lines...\n";
      process_chunk(lines);
    }
  }
  if (!lines.empty()) {
    if (verbose) Rcout << "Processing " << lines.size() << " lines...\n";
    process_chunk(lines);
  }
  
  if (globalAgg.empty()) stop("No edges produced (check cg/mapq/target filters).");
  
  // decide matrix dimension
  uint32_t n_bins = (R2 > 0)
    ? static_cast<uint32_t>(((R2 - R1 + 1) + bin_bp - 1) / bin_bp)
      : (max_bin_seen + 1);
  
  // convert to column-sorted entries (col=j)
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
    // dsCMatrix (upper triangle stored; uplo='U')
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
    mat.slot("i")        = IntegerVector(i.begin(), i.end());
    mat.slot("p")        = IntegerVector(p.begin(), p.end());
    mat.slot("x")        = NumericVector(x.begin(), x.end());
    mat.slot("Dim")      = IntegerVector::create(n_bins, n_bins);
    mat.slot("Dimnames") = List::create(R_NilValue, R_NilValue);
    mat.slot("uplo")     = CharacterVector::create("U");
    mat.slot("factors")  = Rcpp::List::create();   // must be a list (not NULL)
    return mat;
    
  } else {
    // full symmetric dgCMatrix
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
    mat.slot("i")        = IntegerVector(i.begin(), i.end());
    mat.slot("p")        = IntegerVector(p.begin(), p.end());
    mat.slot("x")        = NumericVector(x.begin(), x.end());
    mat.slot("Dim")      = IntegerVector::create(n_bins, n_bins);
    mat.slot("Dimnames") = List::create(R_NilValue, R_NilValue);
    mat.slot("factors")  = Rcpp::List::create();
    return mat;
  }
}

