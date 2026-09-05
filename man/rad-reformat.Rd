\name{rad-reformat}
\alias{rad-reformat}
\alias{rad_reformat}
\title{Reformat or Split RAD-Tagged Sequence Files with the Embedded Core}
\usage{
rad_reformat(
  fastq,
  out_dir = NULL,
  split_bc = FALSE,
  reformat_header = FALSE,
  collapsed_input = FALSE,
  delimiter = "_",
  coordinate = NULL,
  bin_size = 2L,
  output_fastq = NULL,
  threads = 1L,
  chunk_size = 5000L,
  verbose = FALSE
)
}
\arguments{
  \item{fastq}{Input FASTQ, FASTA, or gzip-compressed FASTQ/FASTA path.}

  \item{out_dir}{Output directory for barcode-specific files. Required when
  \code{split_bc = TRUE} and unused otherwise.}

  \item{split_bc}{Whether to write one gzip file per observed \code{CB:Z}
  barcode, using \code{.fq.gz} for FASTQ input and \code{.fa.gz} for FASTA
  input.}

  \item{reformat_header}{Whether to collapse the read name, cell barcode, and
  UMI to \code{QNAME<delimiter>CB<delimiter>UB}. For RAD CLI parity, the full
  input comment is removed, including fields other than the \code{CB} and
  \code{UB} tags. This cannot be combined with \code{coordinate}.}

  \item{collapsed_input}{Whether to parse tagless input names as already
  collapsed \code{QNAME<delimiter>CB<delimiter>UB} values. The default is
  \code{FALSE}, so ordinary read names containing delimiters cannot be
  mistaken for cell barcodes.}

  \item{delimiter}{A single-byte delimiter used for collapsed headers and for
  recognizing already-collapsed headers when \code{collapsed_input = TRUE}.}

  \item{coordinate}{Optional coordinate mode. \code{NULL} disables coordinate
  rewriting; \code{"vizHD-v1"} maps bundled Visium HD BC1/BC2 barcodes to
  spatial bin identifiers. Coordinate mode preserves unrelated comment fields,
  is mutually exclusive with \code{reformat_header}, and counts as a requested
  action.}

  \item{bin_size}{Positive spatial bin size in microns.}

  \item{output_fastq}{Optional aggregate output path for a non-splitting run.
  When omitted, a header or coordinate rewrite replaces \code{fastq}
  atomically in place. It cannot be supplied with \code{split_bc = TRUE}.}

  \item{threads}{Positive requested worker-thread count. A build without
  OpenMP warns and uses one effective thread.}

  \item{chunk_size}{Positive number of records per streaming chunk.}

  \item{verbose}{Whether to print native progress captured during the
  completed run.}
}
\value{
A \code{rad_reformat_result} containing embedded-core provenance, normalized
configuration, exact output paths, processing statistics, and captured log
output.
}
\description{
Streams a sequence file through RAD's native post-processing engine. It can
collapse barcode tags into the read name, split records by cell barcode, or
translate Visium HD barcode pairs into spatial coordinate headers without
launching an external RAD or compression executable.
}
\details{
At least one of \code{split_bc}, \code{reformat_header}, or \code{coordinate}
must request work. Split-only runs leave the input untouched. If splitting and
header reformatting are both requested, records are rewritten before they are
sent to barcode-specific files.

A non-splitting run with no \code{output_fastq} performs a checked temporary
write followed by atomic replacement. The original input is retained if
parsing or writing fails. The exact aggregate and split-file paths are returned
in \code{result$artifacts}. Legal zero-length FASTQ/FASTA records and FASTQ
plus-line identifiers or comments are preserved.
}
\examples{
\dontrun{
# Preserve the demultiplexed input and write one gzip per cell barcode.
split <- rad_reformat(
  "run/demux.fq.gz",
  out_dir = "run/by-cell",
  split_bc = TRUE,
  reformat_header = TRUE
)

# Write collapsed headers to a separate aggregate FASTQ.
rewritten <- rad_reformat(
  "run/demux.fq.gz",
  reformat_header = TRUE,
  output_fastq = "run/demux-collapsed.fq.gz"
)
}
}
