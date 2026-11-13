# rrad 1.3.0

* Update the embedded RAD engine to development commit
  `616754007313ecb6fb875beb3a3fecfd3d83a267`, retaining rrad's validated input,
  captured diagnostics, cancellation, and transactional output protections.
* Support barcode catalogs up to 64 bases with compact storage and indexed
  candidate lookup. Short barcodes continue to use narrow storage.
* Calibrate adapter thresholds from the observed chance-hit distribution,
  expanding the sample when needed. Refine ordinary concatemer grouping,
  adapter clipping, UMI placement, and minimum cDNA length handling.
* Add opt-in concatemer HMM segmentation to `rad_demux()`, with cached model
  parameters, configurable abstention, conservative boundary handling, and
  preservation of unresolved parents in file-output runs.
* Add `rad_discover()` for layout and whitelist discovery and
  `rad_spatial_reformat()` for explicit slide-map assignment to spots, bins,
  or cells with tissue filtering and per-read audit tables.
* Extend whitelist scanning with masked barcode capture, spatial abundance
  selection, coordinate catalogs, and occupancy reporting.
* Allow demultiplexed output to stream to an R connection while report
  artifacts retain checked publication.
* Keep exact reverse-strand anchors from being masked by weaker forward
  matches, preventing retained adapter sequence in extracted cDNA. Reject
  ambiguous gaps at masked barcode boundaries and preserve mirrored HMM
  parameters through calibration-cache reuse.
* Add native regression suites and Linux address/undefined-behavior checks.
  Expand CI to exercise serial and explicit OpenMP builds on both Apple
  Silicon and Intel macOS alongside Linux R release, devel, and oldrel.

# rrad 1.2.2

* Probe native compiler, zlib, thread, and optional OpenMP support during
  installation, and check serial and parallel Linux builds plus Apple Silicon
  and Intel macOS builds in CI.
* Avoid changing R's random-number state in deterministic native routines.
* Preserve literal leading `@` characters in FASTQ identifiers and separate
  appended records when the existing file has no final newline. Refuse append
  operations on truncated or corrupt gzip files without changing their bytes.
* Reject fractional or invalid scalar arguments instead of silently coercing
  them, and allow interruption of long packed-read operations.
* Replace the remaining native assertion with a recoverable error and remove
  compiler warnings without suppressing diagnostics.
* Add a runnable offline vignette covering compact read storage, streaming,
  whitelist discovery, demultiplexing, and pRESTO FASTA export.
* Expand installation and platform documentation, clarify storage and return
  value semantics, and add executable help examples.
* Document the source layout, native provenance maintenance, and release checks.

# rrad 1.2.1

* Provide an embedded, version-pinned RAD engine for demultiplexing, whitelist
  discovery, and tagged-read reformatting without a separate RAD installation.
* Support FASTA conversion and pRESTO header output for downstream IgBLAST use.
* Provide compact read storage, streaming FASTQ input, checked output writes,
  and integrity-verified downloads for large built-in barcode whitelists.
