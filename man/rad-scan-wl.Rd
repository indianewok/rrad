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
  verbose = FALSE,
  capture_mode = c("after_adapter", "masked"),
  spatial = FALSE,
  spatial_slide = NULL,
  spatial_platform = NULL,
  slide_map = NULL,
  slide_dir = NULL,
  bc1_axis = NULL,
  bc2_axis = NULL,
  spatial_selection = c("abundance", "legacy")
)
}
\arguments{
  \item{input}{Input FASTQ, FASTA, or gzip-compressed FASTQ/FASTA path. With
  \code{rescan = TRUE}, a CSV, TSV, or whitespace-delimited barcode/count
  table instead. Supply exactly one of \code{input} and \code{batch_csv},
  or omit both to export a slide coordinate catalog.}

  \item{output_prefix}{Prefix for generated files. The default is
  \code{"barcodes"} for sequence scans and \code{<input>_rescan} for
  statistical rescans.}

  \item{batch_csv}{Optional two-column CSV containing an input FASTQ/FASTA path
  and output prefix per row. A conventional header is accepted. Every row is
  validated before the first native scan begins.}

  \item{adapter_seq}{Adapter or primer immediately upstream of the barcode.
  Required for sequence scans and unused for statistical rescans. Slide scans
  provide defaults for Visium, Visium HD, and Curio Seeker.}

  \item{barcode_length}{Number of barcode bases to extract, from 1 through 64.
  Masked capture infers this from the pattern when the argument is omitted.}

  \item{left_margin, right_margin}{Non-negative RAD extraction offsets. For a
  forward adapter, \code{left_margin} skips bases before extraction and
  \code{right_margin} extends its end; reverse matches use the corresponding
  upstream window.}

  \item{max_reads}{Optional positive maximum number of reads to inspect;
  \code{NULL} scans the entire input.}

  \item{max_error}{Maximum adapter edit-distance ratio from 0 to 1.}

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
  alternative to \code{whitelist}. For HD slide rescans these catalogs can
  translate native spot IDs; other rescans do not use split catalogs.}

  \item{bc1_length, bc2_length}{Optional packed-whitelist barcode lengths,
  each from 1 through 64. They are needed for custom numeric bit-lists, which
  do not encode sequence length; bundled SPLiT-seq keys infer 8.}

  \item{umi_length}{Non-negative UMI length between the adapter and BC1 in
  two-part mode.}

  \item{offset_min, offset_max}{Inclusive BC1 search-offset range. The minimum
  must not exceed the maximum.}

  \item{verbose}{Whether to print native progress captured during each
  completed run.}

  \item{capture_mode}{With \code{"after_adapter"}, N retains its legacy
  wildcard meaning and the barcode follows the adapter. With \code{"masked"},
  uppercase N positions are concatenated into the barcode and lowercase n
  matches a base without capturing it. Margins must be zero. Split scans use
  lowercase n only; component catalogs define their barcodes.}

  \item{spatial}{Use the spatial low-count abundance profile. Supported for
  single, split, batch, and rescan input without coordinate data.}

  \item{spatial_slide}{Slide map path or registered slide ID. Supply
  \code{spatial_platform} as well. Omit input to export the coordinate catalog.
  Slide scans do not accept batch input.}

  \item{spatial_platform}{One of \code{"visium"}, \code{"visium-hd"},
  \code{"curio-seeker"}, \code{"curio-trekker"}, or \code{"stereo-seq"}.}

  \item{slide_map}{Explicit sequence-to-coordinate CSV/TSV map, optionally gzip.}

  \item{slide_dir}{Registry directory containing
  \code{<platform>/<slide>.tsv.gz} (or .csv.gz, .tsv, .csv). Defaults to
  \env{RAD_SLIDE_DIR} or packaged resources.}

  \item{bc1_axis, bc2_axis}{Alternative HD axis tables with \code{barcode} and
  \code{array_col} or \code{array_row} columns, in native 2-micrometre indices.}

  \item{spatial_selection}{\code{"abundance"} applies the low-count abundance
  caller before slide membership and HD occupancy filtering. \code{"legacy"}
  retains the historical slide route: ordinary Visium uses the standard caller,
  HD uses occupancy alone, and other platforms use the low-count caller.}

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

Masked patterns can capture separate blocks, for example
\code{"NNNNACGTnnNN"} captures six bases. Reads are searched in both
orientations and reverse hits are returned in adapter orientation. Only A, C,
G, and T can be captured. An alignment gap within or next to a captured block
is rejected. In masked mode, the edit-distance budget uses the number of
literal pattern positions, excluding both kinds of wildcard.

A slide map contains nucleotide \code{barcode}, \code{x}, and \code{y}
columns; native Curio three-column bead files without a header are also read.
HD maps must supply native \code{array_row}/\code{array_col},
\code{x_um}/\code{y_um}, or exported coordinates with
\code{coordinate_units}. Separate \code{bc1}/\code{bc2} columns and native
\code{s_002um_ROW_COL} identifiers are supported. HD axis tables or ordered
component catalogs avoid expanding their Cartesian product in memory. A bare
whitelist cannot supply physical coordinates for other platforms.

HD selection requires at least 32 distinct native 2-micrometre spots in each
aligned 16-micrometre region (64 possible spots). Repeated observations of one
spot do not increase occupancy. The abundance route retains pre-filter calls
in \code{<prefix>_rad.csv} and \code{<prefix>_rad.txt}; selected slide
coordinates are in \code{<prefix>_coordinates.tsv.gz}, and HD occupancy is in
\code{<prefix>_regions.tsv}. With no input, only the complete coordinate catalog
and a log are written. Native spot-ID rescans require a map or catalogs that
resolve the IDs. Coordinate files are loaded locally; rrad does not download
slide maps.

All files for one run are staged privately and published together after a
successful scan. A failed scan preserves the previous published files. A batch
is validated before scanning, but its completed runs are separate transactions.

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
local({
  work <- tempfile("rrad-scan-example-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- file.path(work, "reads.fasta")
  writeLines(c(">read-1", "ACGTAACC", ">read-2", "ACGTTTGA"), input)
  scan <- rad_scan_wl(
    input, file.path(work, "cells"), adapter_seq = "ACGT",
    barcode_length = 4L, max_error = 0, selection = "above_floor"
  )
  readLines(scan$artifacts$whitelist)
})

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
