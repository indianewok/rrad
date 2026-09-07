# Vendored RAD core

This directory contains the C++ headers required by rrad's embedded RAD
demultiplexing, whitelist-discovery, and reformatting backends.

- Upstream: <https://github.com/indianewok/rad>
- RAD version: `1.0.2`
- Source commit: `24d5ce47e172222c96cd7dd19e94758d32aa5009`
- Copied path: upstream `include/`
- Adapted paths: the complete machine-readable list is in
  `RAD_VENDOR.json`. It covers modified headers under `include/`, the new
  portable OpenMP compatibility header, and the `reformat` section of
  upstream `src/main.cpp` exposed as `include/rad/reformat.hpp`.
- Upstream RAD license: MIT; see `LICENSE`
- Bundled dependency licenses: see the installed
  `THIRD_PARTY_NOTICES.md` and `licenses/` files

The small matching runtime data from upstream `resources/` is installed from
`inst/rad/resources/`. The three large 3M whitelist payloads are deliberately
left out of the package tarball and catalogued instead in
`inst/rad/resources/REMOTE_RESOURCES.tsv`. Each catalog row fixes the upstream
commit URL, exact byte count, and SHA-256 digest. rrad downloads a required
payload only on first use and reuses the verified copy from its versioned user
cache. `rad_download_whitelists()` can download those files ahead of an
offline or batch run, and `rad_resource_path()` finds one saved whitelist.
The cache root can be set with `options(rrad.cache_dir = ...)`; offline
operation can be enforced with `options(rrad.offline = TRUE)` or
`RRAD_OFFLINE=true`. The
ignored, local-only `resources/wl/visium_hd_coordinates.csv.gz` file is
intentionally absent.

Embedding adaptations are deliberately limited to the vendored copy:

1. RAD declarations are isolated in `rrad_rad_core` so upstream global names
   cannot collide with other native libraries loaded into R.
2. Boost.Filesystem is replaced with C++17 `std::filesystem`.
3. Boost.Iostreams gzip streams are replaced with the vendored gzstream/zlib
   implementation.
4. Boost.Lockfree's writer queue is replaced with a bounded standard-library
   queue with worker-error propagation.
5. The upstream `std::hash` specializations are qualified for the isolated
   namespace.
6. When `RAD_EMBEDDED_ONLY` is set, the vendored configuration loaders ignore
   the standalone CLI's `~/.rad/layout_overrides.tsv` and
   `~/.rad/whitelist_overrides.tsv`. The bridge sets and restores this flag
   around every demultiplexing run so bundled resource resolution is
   reproducible and agrees with the reported resource digest.
7. The upstream seqspec translator is isolated as
   `rrad_rad_core::seqspec`, uses `std::filesystem`, and is included outside
   the wrapper that protects third-party global namespaces. The bridge invokes
   it only on layout-cache misses and removes its temporary CSV on success or
   error.
8. Chunked sequence streaming distinguishes clean EOF from kseq parse errors
   and gzip integrity failures, and throws for malformed/truncated records.
   Existing files with unsupported extensions are rejected explicitly;
   regular `fopen` streams and `popen`/pigz streams have distinct close paths.
9. Debug output is additive: passed reads are always written to the primary
   demultiplexed FASTQ, while requested debug artifacts are written
   independently. Debug metrics emit the eight fields declared by their
   header.
10. Explicit embedded-API whitelist arguments use typed loading: kit/global
    inputs always populate the global catalog and custom/auto-generated inputs
    always populate the true catalog, independent of list size.
11. Chunk acquisition is serial, per-read OpenMP exceptions are captured and
    rethrown on the caller, and optional host interrupt callbacks run only at
    serial chunk boundaries (including auto-whitelist scans).
12. Gzip whitelist reads are integrity checked. Output streams have explicit,
    checked finalization after the asynchronous writer drains; layout and
    position-map cache writes use checked atomic replacement; whitelist and
    metrics write/close failures propagate.
13. Remote seqspec onlists use a URL-derived stable cache key plus a sanitized
    basename and an atomically renamed temporary download. The resolved cache
    path is retained in the generated layout artifact for auditability.
14. Cached calibration state is reused only when its layout and position-map
    files are both present. Explicit global and custom whitelist rows are
    whitespace-trimmed; packed-integer values are bounded by the declared
    barcode length; and files with no valid barcode are rejected instead of
    silently producing an empty catalog.
15. Automatic whitelist detection retains a reference whitelist declared by
    the layout as the global correction catalog for the subsequent demultiplex
    pass, matching the standalone RAD pipeline.
16. Before any output is opened, every derived artifact is checked against the
    FASTQ, layout, explicit and layout-derived whitelist inputs, reused caches,
    and imported companion position maps by normalized path and file identity,
    so neither direct paths nor hard/symbolic links can be overwritten. The
    full gzip integrity pass over the primary output polls R interrupts at
    bounded intervals and closes its stream during cancellation.
17. Whitelist discovery exposes the single-barcode, statistical-rescan, and
    split-barcode engines to the bridge. It validates packed barcode lengths,
    treats `N` as a wildcard in exact adapters, accepts case-insensitive FASTX
    and gzip suffixes, reports malformed or corrupt input, captures worker
    exceptions, and polls host interrupts only on serial boundaries, including
    statistical rescans and barcode selection. Split observations are kept as
    bounded aggregate pair counts, and sorted barcode statistics make artifacts
    deterministic across chunk schedules. Every scan artifact is staged and
    committed as one checked transaction, and returned configuration records
    inferred SPLiT-seq barcode lengths. Fuzzy matches retain Edlib's complete
    adapter interval so insertions and deletions cannot shift the barcode;
    forward and reverse extraction require the full mirrored margin window.
18. The CLI's reformat implementation is adapted into a typed, reusable API.
    Record transformation can use OpenMP while ordered writes remain serial;
    split writers are LRU-bounded and staged transactionally, unsafe barcode
    filenames are rejected, FASTA split files use a FASTA suffix, and collapsed
    identifiers are parsed only when explicitly enabled. Gzip sequence and
    coordinate inputs are integrity checked. Legal zero-length records and
    FASTQ plus-line payloads are retained, coordinate conversion replaces only
    its owned tags while preserving unrelated metadata, and conflicting header
    transformations are rejected. Aggregate, in-place, and split outputs are
    installed atomically only after successful validation; the no-overwrite
    path uses a no-clobber filesystem operation.
19. Reformat can emit FASTA directly from FASTQ in the same streamed pass as
    header collapse, coordinate conversion, or barcode splitting. Aggregate
    and per-barcode FASTA outputs retain the existing transactional publication
    and gzip validation guarantees, while unchanged headers preserve complete
    RAD SAM-style tag comments for downstream tools.
20. Reformat can produce canonical pRESTO identifiers for IgBLAST as
    `QNAME|BARCODE=BC_OR_CB|UMI=UB`. It prefers a complete nonempty spatial
    `BC:Z:` value, falls back to `CB:Z:`, recognizes fields independent of their
    order, and leaves split-file selection keyed by `CB:Z:`. It omits absent
    annotations and rejects whitespace plus pRESTO's reserved field,
    assignment, and list delimiters before transactional output publication.

The bridge resolves the installed `rrad` resource bundle through
`system.file("rad/resources", package = "rrad")`. For an unchanged built-in
kit mapping whose file is not bundled, it may also read the corresponding
relative path beneath the read-only `RRAD_RESOURCE_CACHE` overlay prepared by
the R wrapper. The overlay mirrors the packaged resource layout (for example,
`wl/3M-february-2018-3v3.txt_bitlist.csv.gz`). Explicit and custom whitelist
paths never use this fallback. There is no cwd/source-tree fallback, and the R
layer links only individually verified payloads into a private per-call
overlay after validating each manifest byte count and SHA-256 digest. Other
files in the persistent cache remain invisible to native resolution. Existing
layout and position-map caches are ignored unless the caller explicitly
enables cache reuse.

Package configuration link-probes C++17 `std::filesystem` and adds
`-lstdc++fs` only for older toolchains that require it. The default macOS
build is deliberately serial because CRAN's Apple-clang/OpenMP runtime pair is
not portable; the public result distinguishes requested and effective thread
counts. Native execution restores the caller's prior OpenMP thread setting.

`embedded_source_digest` identifies the final patched `include/` tree. It is
SHA-256 over a canonical UTF-8 manifest: enumerate every file beneath
`include/`, sort POSIX relative paths bytewise (`LC_ALL=C`), and append one
line per file as `<lowercase file SHA-256><two spaces><relative path><LF>`.
Hash that complete manifest with SHA-256. The provenance files outside
`include/`, including `RAD_VENDOR.json` itself, are intentionally excluded to
avoid self-reference.

`resource_tree_sha256` uses the same sorted per-file manifest algorithm, with
paths relative to `inst/rad/resources/`. That tree includes
`REMOTE_RESOURCES.tsv`; its pinned payload digests and sizes therefore make the
resource digest commit to the identities of both bundled and on-demand data
without putting the large payload bytes in the package.

The bridge compiles the vendored SSW C/C++ implementation into
`src/rad_core_bridge.cpp`. A dedicated package translation unit compiles RAD's
vendored Edlib implementation. Unused upstream sciplot, kdepp, and ancillary
Parallel Hashmap headers are intentionally omitted from the embedded include
closure.

`RAD_VENDOR.json` records machine-readable source and resource provenance,
including the remote asset set and verification policy.
