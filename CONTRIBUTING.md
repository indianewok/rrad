# Maintaining rrad

`rrad` is the RAD companion package used by SIRUP. Keep its public `rad_*` API
compatible when updating the embedded C++ code, and record user-visible changes
in `NEWS.md`.

## Source layout

- `R/packed_reads.R` and `src/packed_reads.cpp`: compact read tables, FASTQ
  streaming, quality summaries, and persistence.
- `R/rad_core.R`, `R/rad_scan_wl.R`, and `R/rad_reformat.R`: public native
  workflow interfaces and argument validation.
- `R/rad_resources.R`: verified resource downloads and cache management.
- `src/rad_core_bridge.cpp` and `src/rad_scan_wl_bridge.hpp`: the R/native
  boundary, result contracts, host interrupts, and native execution context.
- `src/vendor/rad_core/`: pinned RAD source, embedding adaptations, and
  provenance. `inst/rad/` holds the matching resource manifest and small data.
- `man/`: maintained Rd documentation; `vignettes/`: executable user guide.
- `tests/testthat/`: offline regression tests, including malformed inputs,
  failed writes, native interfaces, and provenance checks.

The Rd files and `NAMESPACE` are maintained explicitly. The source's roxygen
comments are useful API notes but are not a complete regeneration specification;
do not run `roxygen2::roxygenise()` blindly over the package. Update the R
signature, matching Rd usage and arguments, and relevant examples together.

Likewise, inspect changes from `Rcpp::compileAttributes()` carefully. The
`RcppExports` files contain package-specific native registrations and bridge
wrappers as well as generated packed-read bindings; replacing them wholesale
can remove the embedded RAD entry points.

## Checks before release

Install the package's `Depends`, `Imports`, `LinkingTo`, and `Suggests`
dependencies, plus Pandoc to render the HTML vignette and a LaTeX toolchain to
build the PDF reference manual. From the package checkout's parent directory:

```sh
R CMD build rrad
R CMD check --as-cran rrad_<version>.tar.gz
```

Replace `<version>` with the version in `DESCRIPTION`. If LaTeX is unavailable,
`--no-manual` permits a partial local check, but report that limitation and run
the full check in an environment that can build the manual before submission.
Do not skip vignettes, examples, suggested dependencies, or tests for the final
release check. Tests and examples use local synthetic fixtures and must not
download runtime data. Keep examples and tests to at most two native threads.

Run the committed GitHub Actions matrix on Linux and macOS, inspect all results,
and retain the logs for the exact release commit. The Linux release job with
OpenMP builds the PDF manual; all matrix jobs build and check the vignette.
Include R release, development, and old-release checks where practical.
Platform support must be based on
completed installation and checks, not solely on the presence of a workflow.
Review the current [CRAN repository policy](https://cran.r-project.org/web/packages/policies.html)
and [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html)
before an actual submission. CRAN acceptance is a separate review.

## Updating the RAD snapshot

Record the upstream version, commit, copied files, and local adaptations in
`src/vendor/rad_core/RAD_VENDOR.json` and its README. Keep the source and resource
digests synchronized with `inst/rad/RAD_RESOURCES.json` and the native metadata;
the vendored README specifies the canonical hashing algorithm. The provenance
tests verify these identities. Review new or changed dependencies against
`inst/THIRD_PARTY_NOTICES.md` and preserve their complete license texts under
`inst/licenses/`.

For remote whitelists, pin an immutable HTTPS URL, exact byte count, and SHA-256
digest in `inst/rad/resources/REMOTE_RESOURCES.tsv`. Avoid adding large payloads
to the package archive. Exercise offline mode, corrupted-cache recovery, and
failure cleanup without depending on a live network service.
