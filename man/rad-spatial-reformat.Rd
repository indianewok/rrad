\name{rad_spatial_reformat}
\alias{rad_spatial_reformat}
\title{Export spatial assignments using explicit coordinate maps}
\description{
Assign demultiplexed FASTQ reads to native spots, spatial bins, or cells using
an explicit coordinate map. Retain an assignment audit for every input read.
}
\usage{
rad_spatial_reformat(fastq, out_dir, sample, map = NULL,
  bc1_axis = NULL, bc2_axis = NULL, metadata = NULL,
  unit = c("bin", "cell"), bin_size = 2L,
  tissue_only = FALSE, strict = FALSE)
}
\arguments{
\item{fastq}{Input four-line FASTQ, optionally compressed with gzip.}
\item{out_dir}{New output directory. Existing directories are refused.}
\item{sample}{Nonempty sample identifier containing printable ASCII characters,
without whitespace or a colon.}
\item{map}{TSV mapping corrected barcodes to native coordinates. Required
columns are \code{barcode}, \code{array_row}, and \code{array_col}.
Optional columns are \code{in_tissue}, \code{cell_id},
\code{pxl_row_in_fullres}, and \code{pxl_col_in_fullres}.}
\item{bc1_axis, bc2_axis}{Alternative to \code{map}: two TSVs with
\code{barcode,array_col} and \code{barcode,array_row} columns, respectively.
Both are required. Delimited or concatenated barcode pairs are supported.}
\item{metadata}{Optional TSV with the same position columns, keyed by native
spot identifiers such as \code{s_002um_00003_00007-1}. This supplies tissue,
cell, or pixel data for either map-based or axis-based assignment. A matching
metadata row replaces that spot's annotations; missing values do not inherit
values from \code{map}. All TSV inputs may be gzipped.}
\item{unit}{Use a spatial bin or a cell identifier for the output cell barcode.}
\item{bin_size}{Even number of microns, at least 2. Coordinates are zero-based
rows and columns on the native 2-micron grid. Cell assignment requires 2.}
\item{tissue_only}{Require \code{in_tissue = 1}. A missing tissue annotation
does not pass this filter.}
\item{strict}{Exclude reads without a valid assignment. The default retains
them with their original identity and records the reason in the audit.}
}
\details{
Input headers must have complete \code{CB:Z:} and \code{UB:Z:} fields for an
assignment. Assigned records receive the selected unit as \code{CB:Z:};
\code{XB:Z:} preserves the original barcode, \code{XS:Z:} identifies the
sample, and \code{XP:Z:} identifies the native spot. Read identifiers,
sequences, qualities, plus-line payloads, and unrelated tags are preserved.
Duplicate or conflicting reserved tags are reported as invalid assignments.
This function does not deduplicate UMIs.

Coordinates come from explicit maps; they are not inferred from whitelist
order. This is separate from the historical \code{vizHD-v1} coordinate mode
of \code{\link{rad_reformat}}. Processing is serial and interruptible.
Outputs are staged and validated before publication; failed runs remove their
own partial outputs. Input files and existing output directories are preserved.
}
\value{
A \code{rad_spatial_result} list with \code{success}, \code{records_read},
\code{records_written}, counts by assignment status, core provenance, and
output paths. Files are \code{reads.fastq.gz}, \code{assignments.tsv.gz},
and \code{summary.tsv}. The audit includes every read, including excluded reads.
}
\examples{
work <- tempfile("rad-spatial-")
dir.create(work)
reads <- file.path(work, "reads.fq")
map <- file.path(work, "map.tsv")
writeLines(c("@read1 CB:Z:ACGT UB:Z:TGCA", "ACGT", "+", "IIII"), reads)
writeLines(c("barcode\tarray_row\tarray_col", "ACGT\t3\t7"), map)
result <- rad_spatial_reformat(reads, file.path(work, "output"),
                               "sample1", map = map)
result$counts
unlink(work, recursive = TRUE)
}
\seealso{\code{\link{rad_reformat}}, \code{\link{rad_scan_wl}}}
