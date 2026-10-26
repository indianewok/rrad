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
