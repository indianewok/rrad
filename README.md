# rrad

`rrad` is a standalone R package that embeds RAD's demultiplexing,
whitelist-scanning, and tagged-read reformatting engines. Its native code is
compiled when the R package is installed; users do not need a separate `rad`
executable.

The package also provides streaming FASTQ readers and compact vector classes
for read identifiers, nucleotide sequences, and Phred qualities.

## Installation

`rrad` is public and can be installed directly from its pinned GitHub release;
no GitHub credential or separate RAD installation is required.

```r
install.packages("remotes")
remotes::install_github("indianewok/rrad@v1.2.1")
```

## Native RAD workflow

```r
library(rrad)

rad_core_info()

run <- rad_demux(
  layout = "five_prime",
  fastq = "reads.fq.gz",
  out_dir = "run",
  output = "sample",
  threads = 8
)

igblast_fasta <- rad_reformat(
  run,
  output_fastq = "run/sample-igblast.fasta",
  to_fasta = TRUE,
  header_format = "presto"
)
```

`rad_scan_wl()` exposes native whitelist discovery. `rad_stream()` and the
packed read classes support custom R-native demultiplexing and analysis
workflows.

Large built-in 10x whitelists are downloaded only when first requested. Each
download is checked against a pinned byte count and SHA-256 digest and then
reused from R's per-user cache. Use `rad_download_whitelists()` to prepare a
machine before working offline.
