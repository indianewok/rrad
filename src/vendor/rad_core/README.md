# Vendored RAD core

This directory contains the C++ headers required by rrad's embedded RAD
demultiplexing, concatemer HMM, layout/whitelist discovery, spatial assignment,
and reformatting backends.

- Upstream: <https://github.com/indianewok/rad>
- Upstream declared version: `1.0.2` (development snapshot)
- Source commit: `616754007313ecb6fb875beb3a3fecfd3d83a267` (`dev`)
- Resource commit: `24d5ce47e172222c96cd7dd19e94758d32aa5009`; resource
  bytes are unchanged by this source update
- Copied path: upstream `include/`
- Adapted paths: the complete machine-readable list is in
  `RAD_VENDOR.json`. It covers modified headers under `include/`, the new
  portable OpenMP compatibility header, the `reformat` section of upstream
  `src/main.cpp`, and typed adaptations of `src/discovery_*.cpp` and
  `src/spatial.cpp`. The standalone command parser and process-wide stdout
  redirection are not part of the R API.
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
21. Packed barcode catalogs support up to 64 bases using narrow or wide
    storage. Search regions wider than the packed representation remain
    strings; wide alignment uses Edlib instead of overflowing fixed-width
    alignment state. Native unsigned 128-bit support is checked at configure.
22. The concatemer HMM is included in the package namespace with standard and
    SIMD headers outside that namespace. Its internal index assertion is a
    recoverable exception. Calibration workers capture failures and join
    safely on cancellation; host callbacks run only on the calling thread.
    Adapter calibration completes first, then a bounded second input pass
    calibrates the HMM on the final layout topology. Reverse-model parameters
    use existing rows when a layout has only one explicit strand; implicit
    mirrored-anchor and poly-tail states use parameter-only rows without
    positional fields. Cache import preserves those tables without creating
    layout elements.
23. Adaptive misalignment calibration retains strict FASTQ/gzip error
    reporting and polls host interrupts while reading and processing blocks.
    Its chance-hit histogram controls calibrated thresholds; insufficient
    observations retain the documented fallback threshold.
24. Layout/whitelist discovery and explicit-map spatial export are typed
    in-process APIs with input validation, cancellation and checked staged
    publication. Spatial identities come from the supplied slide map or axis
    tables rather than inferred whitelist order.
25. Optional demultiplexed streaming uses caller-owned R connections through
    the bridge on the main thread. It does not redirect process stdout or
    replace process signal handlers. Report artifacts remain staged; primary,
    debug and unresolved-read sequence files are suppressed in stream mode.
26. The umbrella header sets Boost's `BOOST_NO_AUTO_PTR` configuration before
    its includes. RAD uses modern smart pointers; C++17 builds therefore omit
    obsolete `std::auto_ptr` compatibility overloads without disabling compiler
    warnings.

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
`-lstdc++fs` only for older toolchains that require it. It also checks zlib and
POSIX threads, which the asynchronous streaming writer needs even in serial
builds, and unsigned 128-bit arithmetic for 64-base packed barcodes. Supported
64-bit GCC and Clang toolchains provide this type; an unsupported compiler
fails configuration with a specific diagnostic. OpenMP is optional and is enabled only after compiling and loading a
probe with R's selected C++17 compiler. The default macOS build is deliberately
serial because Apple's compiler does not include an OpenMP runtime. Set
`RRAD_OPENMP=no` before installation to force a serial build on any platform,
or `RRAD_OPENMP=yes` to require a working OpenMP toolchain. The default is
`RRAD_OPENMP=auto`. The public result distinguishes requested and effective
thread counts. Native execution restores the caller's prior OpenMP thread
setting.

The layout parser reports stream failures as exceptions instead of process
assertions. Unused diagnostic variables and lambdas have been removed, and
member initialization follows declaration order so compiler warnings remain
useful without being suppressed. These changes are recorded in the vendoring
manifest along with the updated source digest.

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
