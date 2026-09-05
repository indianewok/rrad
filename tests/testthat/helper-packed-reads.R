quality_string <- function(scores, offset = 33L) {
  if (!length(scores)) return("")
  intToUtf8(as.integer(scores) + as.integer(offset))
}

boundary_sequence <- function(n) {
  n <- as.integer(n)
  if (n == 0L) return("")
  bases <- rep(c("A", "C", "T", "G"), length.out = n)
  exceptional <- intersect(c(1L, 7L, 31L, 32L, 63L, 64L, 65L, n), seq_len(n))
  replacements <- rep(c("N", "a", "R", "-"), length.out = length(exceptional))
  bases[exceptional] <- replacements
  paste0(bases, collapse = "")
}

boundary_quality <- function(n, offset = 33L) {
  n <- as.integer(n)
  quality_string(if (n) (seq_len(n) - 1L) %% 94L else integer(), offset)
}

write_text_bytes <- function(path, text) {
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(charToRaw(text), connection)
  invisible(path)
}

fastq_text <- function(ids, sequences, qualities, comments = NULL) {
  stopifnot(length(ids) == length(sequences), length(ids) == length(qualities))
  if (is.null(comments)) comments <- rep(NA_character_, length(ids))
  stopifnot(length(ids) == length(comments))
  records <- Map(
    function(id, sequence, quality, comment) {
      header <- paste0("@", id, if (is.na(comment)) "" else paste0(" ", comment))
      paste0(header, "\n", sequence, "\n+\n", quality, "\n")
    },
    ids, sequences, qualities, comments,
    USE.NAMES = FALSE
  )
  paste0(records, collapse = "")
}

materialized_payload <- function(x, phred_offset = 33L) {
  out <- rad_materialize(x, phred_offset = phred_offset)
  wanted <- intersect(c("seq_id", "seq_comment", "seq", "qual"), names(out))
  as.data.frame(out[, wanted, with = FALSE], stringsAsFactors = FALSE)
}
