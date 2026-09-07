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
  verbose = FALSE,
  to_fasta = FALSE,
  header_format = NULL
)
}
\arguments{
  \item{fastq}{Input FASTQ, FASTA, or gzip-compressed FASTQ/FASTA path, or a
  \code{rad_demux_result} whose aggregate \code{artifacts$fastq} should be
  consumed.}

  \item{out_dir}{Output directory for barcode-specific files. Required when
  \code{split_bc = TRUE} and unused otherwise.}

  \item{split_bc}{Whether to write one gzip file per observed \code{CB:Z}
  barcode, using \code{.fq.gz} for FASTQ input and \code{.fa.gz} for FASTA
  input.}

  \item{reformat_header}{Whether to collapse the read name, cell barcode, and
  UMI to \code{QNAME<delimiter>CB<delimiter>UB}. For RAD CLI parity, the full
  input comment is removed, including fields other than the \code{CB} and
  \code{UB} tags. This is the legacy equivalent of
  \code{header_format = "collapsed"} and cannot be combined with a different
  explicit \code{header_format} or with \code{coordinate}.}

  \item{collapsed_input}{Whether to parse tagless input names as already
  collapsed \code{QNAME<delimiter>CB<delimiter>UB} values. The default is
  \code{FALSE}, so ordinary read names containing delimiters cannot be
  mistaken for cell barcodes.}

  \item{delimiter}{A single-byte delimiter used for collapsed headers and for
  recognizing already-collapsed headers when \code{collapsed_input = TRUE}.}

  \item{coordinate}{Optional coordinate mode. \code{NULL} disables coordinate
  rewriting; \code{"vizHD-v1"} maps bundled Visium HD BC1/BC2 barcodes to
  spatial bin identifiers. Coordinate mode preserves unrelated comment fields,
  is mutually exclusive with transforming header formats, and counts as a
  requested action.}

  \item{bin_size}{Positive spatial bin size in microns.}

  \item{output_fastq}{Optional aggregate output path for a non-splitting run.
  When omitted, a header or coordinate rewrite replaces \code{fastq}
  atomically in place. Despite the legacy argument name, this is also the
  FASTA output path when \code{to_fasta = TRUE}. It cannot be supplied with
  \code{split_bc = TRUE}.}

  \item{threads}{Positive requested worker-thread count. A build without
  OpenMP warns and uses one effective thread.}

  \item{chunk_size}{Positive number of records per streaming chunk.}

  \item{verbose}{Whether to print native progress captured during the
  completed run.}

  \item{to_fasta}{Whether to convert FASTQ records to FASTA by omitting plus
  and quality lines, like \code{seqkit fq2fa}. The read name and complete
  comment, including RAD's SAM-style tags, are preserved unless another
  requested header transformation changes them. Aggregate conversion requires
  an explicit output ending in \code{.fa}, \code{.fasta}, \code{.fa.gz}, or
  \code{.fasta.gz}; split conversion writes \code{.fa.gz} files.}

  \item{header_format}{Optional header transformation. \code{NULL} preserves
  legacy behavior: \code{reformat_header = TRUE} selects
  \code{"collapsed"}; otherwise \code{"preserve"} is used.
  \code{"preserve"} keeps the QNAME and full comment, \code{"collapsed"}
  emits the positional legacy form, and \code{"presto"} emits the
  whitespace-free pRESTO annotation form
  \code{QNAME|BARCODE=BC_OR_CB|UMI=UB}. \code{BARCODE} uses a complete,
  nonempty SAM-style \code{BC:Z:} spatial barcode when present and otherwise
  falls back to \code{CB:Z:}; \code{UMI} uses \code{UB:Z:}. Fields may occur in
  any order, and missing fields are omitted. QNAME and emitted barcode/UMI
  values containing whitespace or the reserved delimiters \code{|}, \code{=},
  or \code{,} are rejected.}
}
\value{
A \code{rad_reformat_result} containing embedded-core provenance, normalized
configuration, exact output paths, processing statistics, and captured log
output.
}
\description{
Streams a sequence file through RAD's native post-processing engine. It can
collapse barcode tags into the read name, emit canonical pRESTO annotations
for IgBLAST, convert FASTQ to FASTA, split records by cell barcode, or translate
Visium HD barcode pairs into spatial coordinate headers without launching an
external RAD or compression executable.
}
\details{
At least one of \code{split_bc}, a transforming \code{header_format},
\code{coordinate}, or \code{to_fasta} must request work. Split-only runs leave
the input untouched. If splitting and header reformatting are both requested,
records are rewritten before they are sent to barcode-specific files. FASTA
conversion composes with header rewriting. With header formatting disabled,
the complete RAD-tagged comment is retained on the FASTA header.

The pRESTO format always writes QNAME first, followed by \code{BARCODE} and
\code{UMI} annotations in that canonical order when their source tags are
present. Other SAM-style fields are intentionally not copied into the
whitespace-free identifier. Barcode-specific splitting remains keyed by
\code{CB:Z:} even when a distinct \code{BC:Z:} value is projected into the
pRESTO \code{BARCODE} annotation. Missing barcode annotations remain absent so
downstream code can fail or filter explicitly.

A non-splitting run with no \code{output_fastq} performs a checked temporary
write followed by atomic replacement. The original input is retained if
parsing or writing fails. The exact aggregate and split-file paths are returned
in \code{result$artifacts}. Legal zero-length FASTQ/FASTA records and FASTQ
plus-line identifiers or comments are preserved unless FASTA output
intentionally removes FASTQ-only fields.
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

# Convert to a canonical pRESTO header in an IgBLAST FASTA query file.
igblast_input <- rad_reformat(
  "run/demux.fq.gz",
  header_format = "presto",
  output_fastq = "run/igblast-input.fa",
  to_fasta = TRUE
)
}
}
