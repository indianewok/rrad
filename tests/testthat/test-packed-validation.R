test_that("rad_pack_reads validates vector types and record lengths", {
  expect_error(rad_pack_reads(1, "A", "!"), "seq_id must be a character")
  expect_error(rad_pack_reads("r", 1, "!"), "seq must be a character")
  expect_error(rad_pack_reads("r", "A", 1), "qual must be NULL or a character")
  expect_error(rad_pack_reads("r", "A", "!", 1),
               "seq_comment must be NULL or a character")

  expect_error(rad_pack_reads(c("r1", "r2"), "A", "!"),
               "ids and sequences must have equal length")
  expect_error(rad_pack_reads("r", "A", c("!", "!")),
               "qualities and sequences must have equal length")
  expect_error(rad_pack_reads("r", "A", "!", c("one", "two")),
               "comments and sequences must have equal length")
  expect_error(rad_pack_reads("r", NA_character_, NA_character_),
               "sequences cannot contain NA")
  expect_error(rad_pack_reads("r", "AC", "!"),
               "sequence and quality lengths differ")
})

test_that("Phred ranges and offsets are validated losslessly", {
  expect_error(rad_pack_reads("r", "A", "!", phred_offset = -1L),
               "phred_offset must be between")
  expect_error(rad_pack_reads("r", "A", "!", phred_offset = 127L),
               "phred_offset must be between")
  expect_error(rad_pack_reads("r", "A", " ", phred_offset = 33L),
               "outside Q0-Q93")

  too_high <- rawToChar(as.raw(127L))
  expect_error(rad_pack_reads("r", "A", too_high, phred_offset = 33L),
               "outside the supported FASTQ ASCII range")

  non_fastq_byte <- rawToChar(as.raw(200L))
  expect_error(rad_pack_reads("r", "A", non_fastq_byte, phred_offset = 126L),
               "outside the supported FASTQ ASCII range")

  highest <- rad_pack_reads("r", "A", "~")
  expect_identical(rad_qual_values(highest)[[1L]], 93L)
  expect_error(as.character(highest$qual, phred_offset = 34L),
               "cannot be represented")
  expect_error(rad_qfraction(highest, -1L), "threshold must be between")
  expect_error(rad_qfraction(highest, 94L), "threshold must be between")
})

test_that("public helpers reject incompatible objects and arguments", {
  expect_error(as_rad_reads(list(seq_id = "r", seq = "A", qual = "!")),
               "x must be a data.frame or data.table")
  expect_error(as_rad_reads(data.frame(seq_id = "r", seq = "A")),
               "missing columns: qual")
  expect_error(rad_materialize(data.frame()), "x must be a rad_reads")
  expect_error(rad_write_fastq(data.frame(), tempfile()), "x must be a rad_reads")
  expect_error(rad_storage_stats(data.frame()), "x must be rad_reads")
  expect_error(rad_qual_values(1:3), "x must be rad_phred")

  packed <- rad_pack_reads(c("r1", "r2"), c("A", "C"), c("!", "5"))
  expect_error(as.integer(packed$qual), "requires one read")
  expect_error(rad_stream_next(list(), 1L), "stream must be a rad_stream")
  expect_error(rad_stream_done(list()), "stream must be a rad_stream")
  expect_error(rad_stream_close(list()), "stream must be a rad_stream")

  path <- tempfile(fileext = ".fastq")
  on.exit(unlink(path), add = TRUE)
  write_text_bytes(path, fastq_text("r", "A", "!"))
  stream <- rad_stream(path)
  on.exit(rad_stream_close(stream), add = TRUE)
  expect_error(rad_stream_next(stream, 0L), "positive integer")
  expect_error(rad_stream_next(stream, NA_integer_), "positive integer")
  expect_error(rad_stream_next(stream, -1L), "positive integer")
})

test_that("missing and mismatched handles cannot be written as FASTQ", {
  missing_quality <- rad_pack_reads("r", "A", NA_character_)
  path <- tempfile(fileext = ".fastq")
  on.exit(unlink(path), add = TRUE)
  expect_error(rad_write_fastq(missing_quality, path), "missing qualities")

  packed <- rad_pack_reads(c("r1", "r2"), c("A", "CC"), c("!", "55"))
  mismatched <- data.table::copy(packed)
  data.table::set(mismatched, j = "qual", value = packed$qual[c(2L, 1L)])
  expect_error(rad_write_fastq(mismatched, path),
               "quality and sequence handles are not aligned")

  missing_store <- packed$seq
  attr(missing_store, ".rad_store") <- NULL
  expect_error(as.character(missing_store), "lost its backing RAD store")

  missing_sequence <- packed$seq
  expect_error(
    missing_sequence[1L] <- bit64::as.integer64(NA),
    "packed RAD handle columns are immutable"
  )
})

test_that("native external pointers are type checked before casting", {
  packed <- rad_pack_reads("r", "A", "!")
  path <- tempfile(fileext = ".fastq")
  on.exit(unlink(path), add = TRUE)
  write_text_bytes(path, fastq_text("r", "A", "!"))
  stream <- rad_stream(path)
  on.exit(rad_stream_close(stream), add = TRUE)

  expect_error(
    packed_fastq_stream_state_cpp(attr(packed$seq, ".rad_store")),
    "not a RAD FASTQ stream"
  )
  fake_handles <- bit64::as.integer64(1L)
  attr(fake_handles, ".rad_store") <- stream$pointer
  expect_error(packed_store_stats_cpp(fake_handles),
               "not a RAD packed-read store")

  expect_identical(as.character(packed$seq), "A")
  expect_s3_class(rad_stream_next(stream, 1L), "rad_reads")
})

test_that("FASTQ writing rejects fields containing line breaks", {
  path <- tempfile(fileext = ".fastq")
  on.exit(unlink(path), add = TRUE)

  bad_id <- rad_pack_reads("id\ncontinued", "A", "!")
  expect_error(rad_write_fastq(bad_id, path), "cannot contain carriage returns or newlines")

  bad_sequence <- rad_pack_reads("id", "A\n", "!!")
  expect_error(rad_write_fastq(bad_sequence, path),
               "cannot contain carriage returns or newlines")

  bad_comment <- rad_pack_reads("id", "A", "!", "comment\rcontinued")
  expect_error(rad_write_fastq(bad_comment, path),
               "cannot contain carriage returns or newlines")
})

test_that("FASTQ validation preserves an existing destination", {
  sentinel <- charToRaw("existing destination\n")
  bad_id <- rad_pack_reads("id\ncontinued", "A", "!")
  bad_output_phred <- rad_pack_reads("id", "A", "~")

  for (extension in c(".fastq", ".fastq.gz")) {
    path <- tempfile(fileext = extension)
    on.exit(unlink(path), add = TRUE)

    writeBin(sentinel, path)
    expect_error(
      rad_write_fastq(bad_id, path),
      "cannot contain carriage returns or newlines"
    )
    expect_identical(readBin(path, what = "raw", n = length(sentinel)), sentinel)

    writeBin(sentinel, path)
    expect_error(
      rad_write_fastq(bad_output_phred, path, phred_offset = 34L),
      "cannot be represented"
    )
    expect_identical(readBin(path, what = "raw", n = length(sentinel)), sentinel)

    writeBin(sentinel, path)
    expect_error(
      rad_write_fastq(rad_pack_reads("id", "A", "!"), path, append = "yes"),
      "append must be TRUE or FALSE"
    )
    expect_identical(readBin(path, what = "raw", n = length(sentinel)), sentinel)

    writeBin(sentinel, path)
    expect_error(
      rad_write_fastq(bad_id, path, append = TRUE),
      "cannot contain carriage returns or newlines"
    )
    expect_identical(readBin(path, what = "raw", n = length(sentinel)), sentinel)
  }
})
