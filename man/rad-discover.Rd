\name{rad_discover}
\alias{rad_discover}
\title{Identify reference barcode catalogs and infer read layouts}
\description{
Compare raw barcode observations with reference catalogs, or propose a read
layout from known static sequence parts. Layout inference is experimental.
}
\usage{
rad_discover(input, out_dir = ".", output = "discovery",
  mode = c("auto", "whitelist", "layout", "both"),
  pattern = NULL, barcode_length = NULL, parts = NULL,
  catalogs = NULL, catalog_lengths = NULL,
  sample_reads = 10000L, sample_bases = 50000000L,
  max_error = 0.15, threads = 1L, overwrite = FALSE)
}
\arguments{
\item{input}{Existing FASTA/FASTQ file, optionally compressed with gzip.}
\item{out_dir}{Output directory. It is created when absent.}
\item{output}{Filename prefix within \code{out_dir}.}
\item{mode}{\code{"auto"} selects whitelist identification when \code{pattern}
is supplied and both outputs otherwise. \code{"whitelist"}, \code{"layout"},
and \code{"both"} request the corresponding outputs explicitly.}
\item{pattern}{Known flanking primer or masked capture pattern, for whitelist
mode only. Uppercase \code{N} captures a barcode base. Lowercase \code{n}
matches a base without capturing it. Other bases must be A, C, G, or T.
At least one literal anchor base is required.}
\item{barcode_length}{Barcode length from 1 to 64. Required for an ordinary
primer. A masked pattern infers the length from its uppercase \code{N} count;
an explicitly supplied value must agree.}
\item{parts}{Optional CSV of additional known static parts with \code{id} and
\code{seq} columns. Ordinary parts must contain at least eight A, C, G, or T
bases; lowercase bases are converted to uppercase. Shorter or noncanonical
rows are skipped. Their row order does not prescribe read order.}
\item{catalogs}{Optional character vector of reference file paths, with
optional names. DNA and packed decimal barcode files are supported. When
\code{NULL}, compare the registered nonspatial catalogs from the pinned
resource bundle; missing downloadable catalogs use the verified resource
cache. Explicit paths define a reproducible, restricted comparison.}
\item{catalog_lengths}{One integer from 0 to 64 per explicit catalog.
Zero infers length from DNA strings. Packed decimal catalogs require their
actual length because leading A bases cannot be recovered from the integer.
An explicit numeric catalog without a known length is an error; an automatic
registered catalog without a known length is reported as unavailable.}
\item{sample_reads, sample_bases}{Positive limits on the prefix sample.
Only complete records that fit both limits are included.}
\item{max_error}{Maximum adapter alignment error ratio from 0 to 1.}
\item{threads}{Positive worker count. Serial builds use one effective worker.}
\item{overwrite}{Allow replacement of this prefix's previous discovery
artifacts. The default refuses existing outputs.}
}
\details{
Each extracted barcode is compared independently with each eligible catalog.
The catalog with the greatest exact matching read count wins; repeated reads
each count. Equal maxima are retained as ties. Distinct barcode counts and
catalog sizes do not break ties. A matching catalog does not establish a
unique chemistry or identify cells. Spatial catalog identification is outside
the scope of this function. Use \code{\link{rad_scan_wl}} to call observed
barcodes after selecting a reference.

Layout inference uses known parts, their orientation, order, and spacing.
It does not learn novel adapter sequences. Unresolved structures are written
to a partial-layout file and must not be used as resolved demultiplexing
recipes. Sampling takes the first records, rather than a random sample.

Outputs are staged, checked, and published together. Existing inputs are
protected from overwrite. No partial result is published after an execution
error. Tests and examples use explicit small local catalogs and require no
network access.
}
\value{
A \code{rad_discovery_result} list. \code{success} describes execution;
\code{complete}, \code{layout_complete}, and \code{whitelist_complete}
describe inference. Ties, no hits, and partial layouts are valid results with
\code{complete = FALSE}. The list includes sample counts, requested/effective
threads, core provenance, notes, unavailable catalog names, output paths, and
the generated \code{whitelists}, \code{whitelist_match}, and \code{parts}
tables when requested. \code{files} uses artifact names such as
\code{whitelist_match.csv} and \code{layout.partial.csv}.
}
\examples{
work <- tempfile("rad-discover-")
dir.create(work)
primer <- "CTACACGACGCTCTTCCGATCT"
barcode <- "ACGTTGCA"
reads <- file.path(work, "reads.fasta")
catalog <- file.path(work, "reference.txt")
writeLines(c(">read1", paste0(primer, barcode)), reads)
writeLines(barcode, catalog)
result <- rad_discover(reads, work, pattern = primer,
  barcode_length = 8L, catalogs = catalog, max_error = 0)
result$whitelist_match
unlink(work, recursive = TRUE)
}
\seealso{\code{\link{rad_demux}}, \code{\link{rad_scan_wl}}}
