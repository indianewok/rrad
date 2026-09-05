\name{rad-scan-wl}
\alias{rad-scan-wl}
\alias{rad_scan_wl}
\title{Discover a Barcode Whitelist with the Embedded RAD Core}
\usage{
rad_scan_wl(
  input = NULL,
  output_prefix = NULL,
  batch_csv = NULL,
  adapter_seq = NULL,
  barcode_length = 16L,
  left_margin = 0L,
  right_margin = 0L,
  max_reads = NULL,
  max_error = 0.3,
  whitelist = NULL,
  threads = 1L,
  chunk_size = 10000L,
  selection = c("high_specificity", "above_floor"),
  rescan = FALSE,
  bc1_whitelist = NULL,
  bc2_whitelist = NULL,
  bc1_length = NULL,
  bc2_length = NULL,
  umi_length = 9L,
  offset_min = 0L,
  offset_max = 3L,
  verbose = FALSE
)
}
\arguments{
  \item{input}{Input FASTQ, FASTA, or gzip-compressed FASTQ/FASTA path. With
  \code{rescan = TRUE}, a CSV, TSV, or whitespace-delimited barcode/count
  table instead. Supply exactly one of \code{input} and \code{batch_csv}.}

  \item{output_prefix}{Prefix for generated files. The default is
  \code{"barcodes"} for sequence scans and \code{<input>_rescan} for
  statistical rescans.}

  \item{batch_csv}{Optional two-column CSV containing an input FASTQ/FASTA path
  and output prefix per row. A conventional header is accepted. Every row is
  validated before the first native scan begins.}

  \item{adapter_seq}{Adapter or primer immediately upstream of the barcode.
  Required for sequence scans and unused for statistical rescans.}

  \item{barcode_length}{Positive number of barcode bases to extract.}

  \item{left_margin, right_margin}{Non-negative RAD extraction offsets. For a
  forward adapter, \code{left_margin} skips bases before extraction and
  \code{right_margin} extends its end; reverse matches use the corresponding
  upstream window.}

  \item{max_reads}{Optional positive maximum number of reads to inspect;
  \code{NULL} scans the entire input.}

  \item{max_error}{Non-negative maximum adapter edit-distance ratio.}

  \item{whitelist}{Optional built-in whitelist key or whitelist file used to
  restrict and validate single-barcode observations.}

  \item{threads}{Positive requested worker-thread count. A build without
  OpenMP warns and uses one effective thread.}

  \item{chunk_size}{Positive number of reads per processing chunk.}

  \item{selection}{Which barcode set is written to the selected text file:
  RAD's \code{"high_specificity"} final calls or every
  \code{"above_floor"} barcode.}

  \item{rescan}{Re-run only RAD's statistical selector on barcode/count data;
  no sequence records are scanned.}

  \item{bc1_whitelist, bc2_whitelist}{Optional built-in keys or files for a
  two-part barcode scan. Supply both or neither. Two-part mode is an
  alternative to \code{whitelist} and is unavailable for rescans.}

  \item{bc1_length, bc2_length}{Optional packed-whitelist barcode lengths,
  each from 1 through 32. They are needed for custom numeric bit-lists, which
  do not encode sequence length; bundled SPLiT-seq keys infer 8.}

  \item{umi_length}{Non-negative UMI length between the adapter and BC1 in
  two-part mode.}

  \item{offset_min, offset_max}{Inclusive BC1 search-offset range. The minimum
  must not exceed the maximum.}

  \item{verbose}{Whether to print native progress captured during each
  completed run.}
}
\value{
A \code{rad_scan_wl_result} containing the embedded-core provenance, the
effective native configuration (including inferred split-barcode lengths),
artifact paths, scan statistics, and captured log output. Batch mode returns a
\code{rad_scan_wl_batch_result}; its \code{runs}
element contains one \code{rad_scan_wl_result} per validated CSV row.
}
\description{
Scans FASTQ or FASTA reads for barcodes following an adapter, optionally restricts
observations to a reference whitelist, and runs RAD's barcode caller without
launching an external executable. Two-part barcode scans, statistical rescans,
and batch input use the same packaged native engine.
}
\details{
A sequence scan writes \code{<prefix>.txt}, \code{<prefix>.csv}, and
\code{<prefix>_scan_wl.log}. The text file contains the selected barcode set;
the CSV retains per-barcode counts and call annotations. A two-part scan also
writes valid-pair counts and a spatial mask.

The selected whitelist can be passed directly to
\code{rad_demux(custom_whitelist = result$artifacts$whitelist)}. No standalone
RAD installation or PATH change is involved. Small reference whitelists ship
with rrad; when a built-in 3M reference is requested, rrad downloads it once,
checks its packaged byte count and SHA-256 digest, and reuses the verified user
copy on later calls. Use \code{\link{rad_download_whitelists}()} to download
the built-in whitelists before an offline or batch run; set
\code{options(rrad.offline = TRUE)} or \code{RRAD_OFFLINE=true} to prohibit
network access.
}
\examples{
\dontrun{
scan <- rad_scan_wl(
  input = "reads.fastq.gz",
  output_prefix = "run/cells",
  adapter_seq = "CTACACGACGCTCTTCCGATCT",
  barcode_length = 16L,
  whitelist = "10x_3v3"
)

demux <- rad_demux(
  layout = "three_prime",
  fastq = "reads.fastq.gz",
  out_dir = "run/demux",
  custom_whitelist = scan$artifacts$whitelist
)
}
}
