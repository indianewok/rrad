test_that("FASTQ streams return independent packed chunks", {
  ids <- c("read1", "read2", "read3", "read4", "read5")
  comments <- c("run=alpha ch=1", NA_character_, "third read", "fourth", NA_character_)
  sequences <- c("ACGT", "NN", "tgcA", "G-RY", "A")
  qualities <- c("!5I~", "II", "?@AB", "#$%&", "~")
  path <- tempfile(fileext = ".fastq")
  on.exit(unlink(path), add = TRUE)
  write_text_bytes(path, fastq_text(ids, sequences, qualities, comments))

  stream <- rad_stream(path, keep_comment = TRUE)
  on.exit(rad_stream_close(stream), add = TRUE)
  expect_false(rad_stream_done(stream))

  first <- rad_stream_next(stream, n = 2L)
  expect_identical(materialized_payload(first), data.frame(
    seq_id = ids[1:2], seq_comment = comments[1:2], seq = sequences[1:2],
    qual = qualities[1:2], stringsAsFactors = FALSE
  ))
  expect_false(rad_stream_done(stream))

  second <- rad_stream_next(stream, n = 10L)
  expect_identical(materialized_payload(second), data.frame(
    seq_id = ids[3:5], seq_comment = comments[3:5], seq = sequences[3:5],
    qual = qualities[3:5], stringsAsFactors = FALSE
  ))
  expect_true(rad_stream_done(stream))
  expect_null(rad_stream_next(stream, n = 1L))

  rm(first)
  invisible(gc())
  expect_identical(as.character(second$seq_id), ids[3:5])
})

test_that("streams can omit comments and enforce close semantics", {
  path <- tempfile(fileext = ".fq")
  on.exit(unlink(path), add = TRUE)
  write_text_bytes(path, fastq_text("one", "AC", "!!", "discard me"))

  stream <- rad_stream(path, keep_comment = FALSE)
  chunk <- rad_stream_next(stream, n = 1L)
  expect_false("seq_comment" %in% names(chunk))
  expect_identical(as.character(chunk$seq_id), "one")
  rad_stream_close(stream)
  rad_stream_close(stream)
  expect_error(rad_stream_next(stream, n = 1L), "FASTQ stream is closed")
})

test_that("plain and gzip FASTQ writers round-trip reordered rows", {
  ids <- c("plain-1", "plain-2", "plain-3")
  comments <- c("alpha=1", NA_character_, "gamma=3")
  sequences <- c("ACGTN", "tgcA", "G-RY")
  qualities <- c("!5I~~", "?@AB", "#$%&")
  packed <- rad_pack_reads(ids, sequences, qualities, comments)
  order <- c(3L, 1L, 2L)
  reordered <- packed[order]

  for (extension in c(".fastq", ".fastq.gz")) {
    path <- tempfile(fileext = extension)
    on.exit(unlink(path), add = TRUE)
    expect_identical(
      invisible(rad_write_fastq(reordered, path)),
      normalizePath(path, mustWork = TRUE)
    )
    reread <- rad_read_fastq(path, keep_comment = TRUE)
    expect_identical(materialized_payload(reread), data.frame(
      seq_id = ids[order], seq_comment = comments[order], seq = sequences[order],
      qual = qualities[order], stringsAsFactors = FALSE
    ))
  }
})

test_that("plain FASTQ writer produces standard four-line records", {
  packed <- rad_pack_reads(
    c("r1", "r2"), c("AC", "TG"), c("!~", "5I"),
    c("comment one", NA_character_)
  )
  path <- tempfile(fileext = ".fastq")
  on.exit(unlink(path), add = TRUE)
  rad_write_fastq(packed, path)

  expect_identical(readLines(path, warn = FALSE),
                   c("@r1 comment one", "AC", "+", "!~",
                     "@r2", "TG", "+", "5I"))
})

test_that("FASTQ append is staged and preserves complete records", {
  first <- rad_pack_reads("first", "AC", "!!")
  second <- rad_pack_reads("second", "TG", "II")

  for (extension in c(".fastq", ".fastq.gz")) {
    path <- tempfile(fileext = extension)
    on.exit(unlink(path), add = TRUE)
    rad_write_fastq(first, path)
    rad_write_fastq(second, path, append = TRUE)

    restored <- rad_read_fastq(path)
    expect_identical(as.character(restored$seq_id), c("first", "second"))
    expect_identical(as.character(restored$seq), c("AC", "TG"))
    expect_identical(as.character(restored$qual), c("!!", "II"))
  }
})

test_that("rad_read_fastq honors max_reads and handles an empty file", {
  path <- tempfile(fileext = ".fastq")
  empty <- tempfile(fileext = ".fastq")
  on.exit(unlink(c(path, empty)), add = TRUE)
  write_text_bytes(path, fastq_text(c("a", "b", "c"), c("A", "C", "G"),
                                    c("!", "5", "I")))
  write_text_bytes(empty, "")

  limited <- rad_read_fastq(path, max_reads = 2L)
  expect_identical(as.character(limited$seq_id), c("a", "b"))
  expect_identical(as.character(limited$seq), c("A", "C"))

  no_reads <- rad_read_fastq(empty)
  expect_s3_class(no_reads, "rad_reads")
  expect_equal(nrow(no_reads), 0L)
  expect_true("seq_comment" %in% names(no_reads))
})

test_that("malformed FASTQ and FASTA input fail explicitly", {
  truncated <- tempfile(fileext = ".fastq")
  fasta <- tempfile(fileext = ".fasta")
  bad_phred <- tempfile(fileext = ".fastq")
  on.exit(unlink(c(truncated, fasta, bad_phred)), add = TRUE)
  write_text_bytes(truncated, "@bad\nACGT\n+\nIII\n")
  write_text_bytes(fasta, ">not-fastq\nACGT\n")
  write_text_bytes(bad_phred, "@bad-score\nA\n+\n!\n")

  expect_error(rad_read_fastq(truncated), "truncated|length-mismatched")
  expect_error(rad_read_fastq(fasta), "FASTA record")
  expect_error(rad_read_fastq(bad_phred, phred_offset = 64L), "outside Q0-Q93")
})

test_that("gzip FASTQ streams reject missing or corrupt integrity trailers", {
  path <- tempfile(fileext = ".fastq.gz")
  corrupt <- tempfile(fileext = ".fastq.gz")
  on.exit(unlink(c(path, corrupt)), add = TRUE)
  connection <- gzfile(path, open = "wb")
  writeLines(c("@read", "ACGT", "+", "IIII"), connection)
  close(connection)

  bytes <- readBin(path, what = "raw", n = file.info(path)$size)
  expect_gt(length(bytes), 8L)
  corrupt_bytes <- bytes
  corrupt_bytes[length(corrupt_bytes) - 7L] <- as.raw(bitwXor(
    as.integer(corrupt_bytes[length(corrupt_bytes) - 7L]), 1L
  ))
  writeBin(corrupt_bytes, corrupt)
  expect_error(rad_read_fastq(corrupt), "error while reading FASTQ stream")

  writeBin(bytes[seq_len(length(bytes) - 8L)], path)

  expect_error(
    rad_read_fastq(path),
    "gzip stream was complete|unexpected end of file"
  )
})

test_that("a stream cannot resume after a destructive parsing failure", {
  path <- tempfile(fileext = ".fastq")
  on.exit(unlink(path), add = TRUE)
  writeLines(c(
    "@good", "A", "+", "!",
    "@bad", "AC", "+", "!",
    "@would-be-skipped", "G", "+", "I"
  ), path)

  stream <- rad_stream(path)
  on.exit(rad_stream_close(stream), add = TRUE)
  expect_error(
    rad_stream_next(stream, n = 2L),
    "truncated|length-mismatched"
  )
  expect_error(
    rad_stream_next(stream, n = 1L),
    "cannot be resumed after a previous read failure"
  )
})
