# rrad

`rrad` is the RAD companion R package for SIRUP. It embeds the
[RAD](https://github.com/indianewok/rad) C++ engines for demultiplexing,
barcode whitelist discovery, and tagged-read reformatting. Native code is
compiled during installation; no separate `rad` executable is required.

The package also provides streaming FASTQ readers and compact read tables for
identifiers, sequences, and Phred qualities. It is intended for Linux and macOS.

## Installation

Install the published GitHub release:

```r
install.packages("remotes")
remotes::install_github("indianewok/rrad@v1.2.2")
```

For the current development version:

```r
install.packages(c("remotes", "knitr", "rmarkdown"))
remotes::install_github("indianewok/rrad", build_vignettes = TRUE)
```

Building the vignette also requires Pandoc (available with RStudio or as a
standalone installation). Omit `build_vignettes = TRUE` to install without
building the guide.

Source installation requires R >= 4.0.0, a C++17 compiler, GNU make, zlib
development headers, and POSIX threads. R dependencies are installed by `remotes`.

- On macOS, install R and Apple's Command Line Tools (`xcode-select --install`).
  The default build uses one native worker and does not require a separate
  OpenMP runtime.
- On Debian or Ubuntu, install R's development tools and zlib headers, for
  example `sudo apt-get install r-base-dev zlib1g-dev`. Linux builds use R's
  configured OpenMP support when it passes the compiler and runtime probes.

To force a serial installation, set `Sys.setenv(RRAD_OPENMP = "no")` before
installing. `RRAD_OPENMP=yes` requires a working OpenMP toolchain and fails
configuration if it is unavailable; it is intended for explicitly configured
compilers. Unset the variable to return to the platform default.

Check the installed build with `rad_core_info()$build`. A request for more than
one worker on a serial build produces a warning and runs with one worker;
results record both requested and effective thread counts. Windows is outside
the package's current supported scope (`OS_type: unix`).

## Start with a small read table

This example runs without external data or network access:

```r
library(rrad)

reads <- rad_pack_reads(
  seq_id = c("read-1", "read-2"),
  seq = c("ACGTN", "TGCAT"),
  qual = c("IIIII", "!!!!!")
)
rad_qmean(reads)
selected <- reads[reads$mean_q >= 20, ]
rad_materialize(selected)
rad_validate(selected, deep = TRUE)
```

Subsetting retains the shared native storage. Use `rad_compact()` after a large
filter to rebuild storage for just the retained reads, `rad_rbind_reads()` to
combine independently read chunks, and `rad_save()` / `rad_load()` to save and
restore tables across R sessions. Plain `saveRDS()` cannot preserve the native
backing pointers. See `?rad_pack_reads` for packed-column and join semantics.

## Native RAD workflow

For your own FASTQ file and a matching read layout:

```r
library(rrad)

run <- rad_demux(
  layout = "five_prime",
  fastq = "reads.fq.gz",
  out_dir = "run",
  output = "sample",
  threads = 1L
)

igblast_input <- rad_reformat(
  run,
  output_fastq = "run/sample-igblast.fasta",
  to_fasta = TRUE,
  header_format = "presto"
)
igblast_input$artifacts$fasta
```

Choose the layout and barcode kit for the actual library chemistry; a layout
may also be a custom CSV or seqspec YAML. `rad_scan_wl()` discovers a whitelist,
and `rad_demux(auto_whitelist = TRUE)` performs discovery before demultiplexing.
The scan's selected whitelist is available in `scan$artifacts$whitelist` for
use with `rad_demux(custom_whitelist = ...)`.

The [workflow vignette](vignettes/rrad-workflow.Rmd) runs a complete synthetic
example, including whitelist scanning, demultiplexing, pRESTO FASTA export, and
bounded-memory streaming. After installing with vignettes, open it with
`vignette("rrad-workflow", package = "rrad")`.

## Offline operation and reproducibility

Small built-in layouts and whitelists are bundled. The large built-in 10x
whitelists are downloaded only when requested, checked against a pinned byte
count and SHA-256 digest, and cached for reuse. Prepare them before offline work:

```r
rad_download_whitelists(c("10x_3v3", "10x_3v4", "10x_5v3"))
options(rrad.offline = TRUE)
```

Use `options(rrad.cache_dir = "/path/to/cache")` or `RRAD_CACHE_DIR` to select
the cache location. `RRAD_OFFLINE=true` also disables downloads.
`rad_resource_path("10x_3v3", download = FALSE)` checks for an available cached
whitelist without downloading it. `rad_clear_resource_cache(confirm = TRUE)`
removes the package-managed resource cache. See `?rad_resource_path` for details.

Keep each run's returned `core`, `config`, and `artifacts` fields with your
analysis record. `rad_core_info()` reports the upstream RAD commit, patched
source digest, and resource digest. These distinguish the embedded code from
an independently installed RAD executable. Package installation, loading, and
the included examples, vignette, and tests do not download runtime resources.

## Development and licensing

See [CONTRIBUTING.md](CONTRIBUTING.md) for package checks and the source layout.
The CI workflow checks Linux and macOS; review the completed checks before
release.

`rrad` is licensed under GPL-3 or later. Vendored components retain their own
licenses and attribution; see [third-party notices](inst/THIRD_PARTY_NOTICES.md)
and [RAD provenance](src/vendor/rad_core/README.md).
