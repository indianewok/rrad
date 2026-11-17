# Maintaining rrad

`rrad` is the RAD companion package used by SIRUP. Keep the public `rad_*` API
compatible and record user-visible changes in `NEWS.md`.

## Source layout

- `R/packed_reads.R`, `src/packed_reads.cpp`: compact read tables, streaming,
  quality summaries, and persistence.
- `R/rad_core.R`, `R/rad_scan_wl.R`, `R/rad_reformat.R`, `R/rad_discovery.R`:
  native workflow interfaces.
- `R/rad_resources.R`: verified resource downloads and cache management.
- `src/rad_core_bridge.cpp`, `src/rad_scan_wl_bridge.hpp`,
  `src/rad_discovery_bridge.hpp`: R/native bindings, interrupts, and output
  publication.
- `src/vendor/rad_core/`: pinned RAD source and embedding adaptations.
  `inst/rad/` contains the resource manifest and small reference data.
- `man/`, `vignettes/`: API documentation and the executable workflow guide.
- `tests/testthat/`, `tests/native/`: R workflows and native algorithm,
  memory, and I/O regressions.

Rd files and `NAMESPACE` are maintained explicitly; the roxygen comments are
not a complete regeneration specification. `RcppExports` also contains manual
RAD registrations. Update matching R signatures, Rd entries, and native
registrations together, and review any generator output before replacing them.

## Validation

Use local synthetic fixtures, keep tests offline, and limit native test runs
to two threads. Add tests for observable behavior or a reproduced failure;
source spelling, exact diagnostic prose, and incidental object layout are not
API contracts.

From the checkout's parent directory, with package dependencies, Pandoc, and
LaTeX installed:

```sh
R CMD build rrad
R CMD check --as-cran rrad_<version>.tar.gz
```

Use the version in `DESCRIPTION`. `--no-manual` allows a local check without
LaTeX; the release check also builds the PDF manual.

Normal CI runs four package configurations: Linux R release/devel and Mac
ARM/Intel release. Linux also runs the native suites with address and
undefined-behavior sanitizers. Before a release, run the workflow manually with
`extended` enabled. It adds Linux oldrel, serial and ARM builds, Mac OpenMP,
PDF manual validation, and native suites on each release configuration.
Inspect results for the release commit before submission.

The standalone native runner uses R's C++17 compiler and installed BH headers:

```sh
bash tests/native/run-tests.sh
RRAD_NATIVE_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
  bash tests/native/run-tests.sh
```

It builds in a temporary directory. Compiler/runtime support is required for
sanitizers or `RRAD_NATIVE_FLAGS="-fopenmp"`; use `RRAD_NATIVE_LIBS` for link-only
options. The disk-full fixture uses Linux's `/dev/full`; other output failure
and cancellation checks run on both systems.

## Updating RAD

Record upstream commits, adapted paths, and patches in
`src/vendor/rad_core/RAD_VENDOR.json`. Its [README](src/vendor/rad_core/README.md)
describes embedding and the canonical source/resource hashes. Update native
metadata and `inst/rad/RAD_RESOURCES.json` when their corresponding inputs change.
Source and resource commits may differ if an update leaves resource bytes
unchanged.

Preserve dependency licenses in `inst/licenses/` and update
`inst/THIRD_PARTY_NOTICES.md` for dependency changes. Large whitelists belong in
`inst/rad/resources/REMOTE_RESOURCES.tsv`, with an immutable HTTPS URL, byte
count, and SHA-256, rather than in the package archive.

For submission requirements, consult the [CRAN policy](https://cran.r-project.org/web/packages/policies.html)
and [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html).
