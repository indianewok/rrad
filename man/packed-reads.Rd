\name{packed-reads}
\alias{packed-reads}
\alias{rad_pack_reads}
\alias{as_rad_reads}
\alias{rad_stream}
\alias{rad_stream_next}
\alias{rad_stream_done}
\alias{rad_stream_close}
\alias{rad_read_fastq}
\alias{rad_materialize}
\alias{rad_write_fastq}
\alias{rad_storage_stats}
\alias{rad_compact}
\alias{rad_rbind_reads}
\alias{rad_save}
\alias{rad_load}
\alias{rad_validate}
\alias{rad_qual_values}
\alias{rad_qmean}
\alias{rad_qmin}
\alias{rad_qfraction}
\alias{as.character.rad_id}
\alias{as.character.rad_seq}
\alias{as.character.rad_phred}
\alias{as.integer.rad_phred}
\alias{format.rad_id}
\alias{format.rad_seq}
\alias{format.rad_phred}
\alias{print.rad_id}
\alias{print.rad_seq}
\alias{print.rad_phred}
\alias{is.na.rad_id}
\alias{is.na.rad_seq}
\alias{is.na.rad_phred}
\alias{anyNA.rad_id}
\alias{anyNA.rad_seq}
\alias{anyNA.rad_phred}
\alias{Ops.rad_id}
\alias{duplicated.rad_id}
\alias{anyDuplicated.rad_id}
\alias{unique.rad_id}
\alias{xtfrm.rad_id}
\alias{sort.rad_id}
\alias{[.rad_reads}
\alias{[<-.rad_reads}
\alias{[[<-.rad_reads}
\alias{$<-.rad_reads}
\alias{[<-.rad_id}
\alias{[[<-.rad_id}
\alias{[<-.rad_seq}
\alias{[[<-.rad_seq}
\alias{[<-.rad_phred}
\alias{[[<-.rad_phred}
\alias{c.rad_id}
\alias{c.rad_seq}
\alias{c.rad_phred}
\alias{rep.rad_id}
\alias{rep.rad_seq}
\alias{rep.rad_phred}
\title{Compact handle-backed read storage}
\description{
Store read identifiers, nucleotide sequences, and Phred qualities in shared
packed C++ buffers while exposing one atomic handle per row of a real
\code{data.table}. Identifiers use lossless UUID or front coding, sequences
use a two-bit A/C/T/G core with an exception lane, and qualities use lossless
seven-bit numeric Phred values.
}
\usage{
rad_pack_reads(seq_id, seq, qual = NULL, seq_comment = NULL,
  phred_offset = 33L)

as_rad_reads(x, id_col = if ("seq_id" \%in\% names(x)) "seq_id" else "id",
  seq_col = "seq", qual_col = "qual", comment_col = NULL,
  phred_offset = 33L)

rad_stream(path, phred_offset = 33L, keep_comment = TRUE)
rad_stream_next(stream, n = 50000L)
rad_stream_done(stream)
rad_stream_close(stream)

rad_read_fastq(path, max_reads = .Machine$integer.max,
  phred_offset = 33L, keep_comment = TRUE)
rad_write_fastq(x, path, phred_offset = 33L, append = FALSE)
rad_materialize(x, phred_offset = 33L)

rad_storage_stats(x)
rad_compact(x)
rad_rbind_reads(..., use.names = TRUE, fill = TRUE, phred_offset = 33L)
rad_save(x, file, compress = "xz")
rad_load(file)
rad_validate(x, deep = FALSE)

rad_qual_values(x)
rad_qmean(x)
rad_qmin(x)
rad_qfraction(x, threshold = 20L)
}
\arguments{
  \item{seq_id}{Character vector of primary read identifiers.}
  \item{seq}{Character vector of nucleotide sequences.}
  \item{qual}{Optional character vector of FASTQ quality strings.}
  \item{seq_comment}{Optional character vector of header comments.}
  \item{phred_offset}{Integer quality-character offset; normally 33. It is
  never inferred automatically.}
  \item{x}{A conventional read table, a \code{rad_reads} table, or a packed
  handle column as appropriate.}
  \item{id_col,seq_col,qual_col,comment_col}{Column names used when converting
  a conventional table.}
  \item{path,file}{Input or output path. Gzip FASTQ is selected by a
  \file{.gz} suffix.}
  \item{keep_comment}{Whether to retain FASTQ comments separately.}
  \item{stream}{A \code{rad_stream} iterator.}
  \item{n,max_reads}{Maximum number of records to read.}
  \item{append}{Append to an existing FASTQ output.}
  \item{use.names,fill}{Arguments used when combining tables.}
  \item{compress}{Compression accepted by \code{saveRDS}.}
  \item{deep}{Also decode data and validate record lengths and cached quality
  summaries.}
  \item{threshold}{Phred threshold for the returned fraction.}
  \item{...}{Additional packed tables or method arguments.}
}
\details{
The \code{seq_id}, \code{seq}, and \code{qual} columns are eight-byte
\code{integer64} subclasses. Their values are descriptors, not packed payloads;
the table columns retain a shared external pointer that owns the buffers.
Subsetting and reordering therefore copy handles only.

Common QC summaries are ordinary numeric columns computed during ingestion.
Use \code{rad_materialize()} only when downstream software requires character
vectors. Use \code{rad_rbind_reads()} rather than binding tables backed by
different stores directly. Use \code{rad_save()} rather than \code{saveRDS()}
because an external pointer is not portable across R sessions.

Canonical lower-case UUID identifiers occupy 16 payload bytes. Other IDs are
front-coded with restart points. ID handles use a stable 63-bit fingerprint so
duplicate IDs group together; construction rejects any within-store fingerprint
collision. IDs and comments supplied through the R API are normalized to UTF-8
before hashing and storage.

Fingerprint equality permits grouping within a store without decoding. Native
R comparisons decode IDs and are exact. A \code{data.table} join between
independent stores compares 63-bit fingerprints directly, so it is
probabilistic rather than collision-verified; materialize both keys when exact
cross-store equality is required. A handle-only index table can participate in
such a join, but a full \code{rad_reads[i]} result that would mix an ID column
from one store with sequence and quality columns from another is rejected.
An outer join that carries unmatched packed IDs from multiple independent
stores cannot retain both backing pointers in one atomic vector. Materialize
the join key with \code{as.character()} before such a union join. Likewise,
\code{data.table::setorder()} sees the numeric fingerprint; use base
\code{sort()} or order the decoded key when lexical order is required.

Packed columns can be subset, reordered, repeated, and concatenated while
retaining their store. Concatenation is limited to the same packed field in the
same backing store; use \code{rad_rbind_reads()} to combine stores. Direct
replacement is intentionally rejected because changing a handle without
repacking could detach a row from its sequence or quality payload. The usual
\code{DT[, column := value]}, \code{DT$column <- value}, and
\code{DT[[column]] <- value} forms reject assignments to packed core columns;
do not use the lower-level \code{data.table::set()} escape hatch to replace
them. Every rrad operation validates alignment before materializing or writing,
so a table altered through that low-level escape hatch fails explicitly rather
than silently reassociating read data. Ordinary annotation columns remain
freely mutable. A projection or grouped result that omits packed core columns
is returned as an ordinary \code{data.table}, not a malformed
\code{rad_reads} object. Save,
compaction, and safe row binding validate alignment before rebuilding a store
and reject additional columns containing nested native references.

Non-appending and appending FASTQ writes, as well as \code{rad_save()}, are
assembled in a private sibling staging directory and installed at the requested
path only after the complete write succeeds. Append therefore copies the
existing destination before extending it; this favors failure safety over
in-place append speed.

The functions whose names begin with \code{packed_} are registered native
bridges used by the public R interface and are not intended for direct use.
}
\value{
\code{rad_pack_reads()}, \code{rad_read_fastq()}, and
\code{rad_stream_next()} return a \code{rad_reads} object inheriting from
\code{data.table}. Quality functions return numeric summaries or a list of
integer vectors. Stream state functions return logical values. Other functions
return the object or path described above.
}
\seealso{\code{data.table::data.table}}
