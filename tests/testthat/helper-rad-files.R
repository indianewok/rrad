rad_test_fastq <- function(path, sequences, comments = rep("", length(sequences)),
                            plus = rep("+", length(sequences))) {
  records <- unlist(lapply(seq_along(sequences), function(i) c(
    paste0("@read", i, if (nzchar(comments[[i]])) paste0(" ", comments[[i]]) else ""),
    sequences[[i]], plus[[i]], strrep("I", nchar(sequences[[i]])))), use.names = FALSE)
  writeLines(records, path)
  invisible(path)
}

rad_test_lines <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  readLines(con, warn = FALSE)
}

rad_test_fasta <- function(path, sequences) {
  writeLines(as.vector(rbind(paste0(">read-", seq_along(sequences)), sequences)), path)
}
