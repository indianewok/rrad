\name{rad-core}
\alias{rad-core}
\alias{rad_core_info}
\alias{rad_demux}
\title{Embedded RAD demultiplexing core}
\description{
Inspect and run the version-pinned RAD core compiled into the R package. These
functions do not discover or launch an external \command{rad} executable.
}
\usage{
rad_core_info()

rad_demux(
  layout,
  fastq,
  out_dir = ".",
  output = "demux",
  kit = NULL,
  global_whitelist = NULL,
  custom_whitelist = NULL,
  threads = 1L,
  chunk_size = 5000L,
  max_reads = NULL,
  auto_whitelist = FALSE,
  scan_adapter = NULL,
  scan_barcode_length = NULL,
  scan_max_error = 0.3,
  scan_max_reads = NULL,
  scan_chunk_size = NULL,
  scan_threads = NULL,
  scan_selection = c("high_specificity", "above_floor"),
  whitelist_mutation = NULL,
  generated_mutation = NULL,
  bc_correction_mode = "offensive",
  joint_bc_mode = "default",
  rc_umi = TRUE,
  min_read_length = NULL,
  write_debug = FALSE,
  verbose = FALSE,
  reuse_cache = FALSE
)
}
\arguments{
  \item{layout}{A built-in RAD layout key or path to a layout CSV or seqspec.}
  \item{fastq}{Path to an input FASTQ or gzip-compressed FASTQ.}
  \item{out_dir}{Output directory, created recursively when necessary.}
  \item{output}{Output filename prefix within \code{out_dir}.}
  \item{kit}{Optional built-in global-whitelist kit key. This is an
  alternative to \code{global_whitelist}.}
  \item{global_whitelist}{Optional path to a global whitelist. Supplying this
  together with \code{kit} is an error.}
  \item{custom_whitelist}{Optional path to a custom whitelist.}
  \item{threads}{Positive number of requested worker threads. Builds without
  OpenMP, including the default macOS build, warn and use one effective
  thread. The result records both values.}
  \item{chunk_size}{Positive number of reads per processing chunk.}
  \item{max_reads}{Optional positive read limit. \code{NULL} processes the
  complete input.}
  \item{auto_whitelist}{Whether to run RAD's barcode-whitelist scan before
  demultiplexing. This is an alternative to \code{custom_whitelist}.}
  \item{scan_adapter}{Optional adapter sequence for automatic whitelist
  scanning. By default RAD derives the nearest upstream static element from
  the read layout.}
  \item{scan_barcode_length}{Optional positive barcode length for automatic
  scanning. By default RAD derives it from the layout's barcode element.}
  \item{scan_max_error}{Non-negative maximum adapter edit-distance ratio used
  during automatic whitelist scanning.}
  \item{scan_max_reads}{Optional positive scan-specific read limit. The
  demultiplexing \code{max_reads} value is used when this is \code{NULL}.}
  \item{scan_chunk_size}{Optional positive scan-specific chunk size. The
  demultiplexing \code{chunk_size} value is used when this is \code{NULL}.}
  \item{scan_threads}{Optional positive scan-specific requested thread count.
  The demultiplexing \code{threads} value is used when this is \code{NULL}.}
  \item{scan_selection}{Automatic-whitelist selection policy: RAD's
  \code{"high_specificity"} calls or every \code{"above_floor"} barcode.}
  \item{whitelist_mutation,generated_mutation}{Optional non-negative mutation
  distances.}
  \item{bc_correction_mode}{Barcode correction policy: \code{"offensive"} or
  \code{"defensive"}.}
  \item{joint_bc_mode}{Multi-barcode policy: \code{"default"} or
  \code{"strict"}.}
  \item{rc_umi}{Whether to reverse-complement UMIs extracted from
  reverse-oriented reads.}
  \item{min_read_length}{Optional non-negative minimum cDNA length.}
  \item{write_debug}{Whether to write RAD debug artifacts in addition to the
  primary demultiplexed FASTQ.}
  \item{verbose}{Whether to print native progress captured during the
  completed run.}
  \item{reuse_cache}{Whether to reuse an existing layout and position-map
  cache for the selected output prefix. The default, \code{FALSE}, rebuilds
  both files to prevent stale state from an earlier run being selected
  implicitly.}
}
\details{
\code{rad_core_info()} is the runtime compatibility handshake. It reports the
core API version, exact upstream RAD source commit, patched embedded-source
digest, resource schema and digest, features, and build details.
\code{rad_demux()} verifies this handshake before every run and records it in
the returned result.

The implementation executes in the package process through Rcpp. Users do not
need to install the RAD command-line program or modify \code{PATH}. Small
layouts and whitelist tables ship with rrad. A required 3M whitelist is
downloaded automatically on first use from a version-pinned RAD source URL,
verified against its packaged byte count and SHA-256 digest, and then reused
on later runs. Custom whitelist files are used directly and are never copied.
Use \code{\link{rad_download_whitelists}()} to download
the built-in whitelists before an offline run,
\code{\link{rad_resource_path}()} to find one saved whitelist file, and
\code{\link{rad_clear_resource_cache}()} to remove downloaded files. The
cache root can be overridden with
\code{options(rrad.cache_dir = "/path")} or \code{RRAD_CACHE_DIR}. The
resource helper functions also accept a per-call \code{cache_dir} argument.
Set \code{options(rrad.offline = TRUE)} or \code{RRAD_OFFLINE=true} to
prohibit network access.

Set \code{reuse_cache = TRUE} only when intentionally continuing a prior run
with the same input layout and compatible calibration data.
}
\value{
\code{rad_core_info()} returns a \code{rad_core_info} list.
\code{rad_demux()} returns a \code{rad_demux_result} list with fields
\code{status}, \code{backend}, \code{core}, \code{input}, \code{config},
\code{artifacts}, \code{stats}, \code{scan}, and \code{log}. \code{scan} is
\code{NULL} unless \code{auto_whitelist = TRUE}, when it is a complete
\code{rad_scan_wl_result}; \code{config$scan_effective} records its
layout-derived and explicit effective settings.
}
\examples{
rad_core_info()

\dontrun{
run <- rad_demux(
  layout = "three_prime",
  fastq = "reads.fastq.gz",
  out_dir = "rad-output"
)
run$artifacts$fastq
}
}
