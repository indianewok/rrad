test_that("IDs, sequences, and seven-bit Phred round-trip at packing boundaries", {
  lengths <- c(0L, 1L, 2L, 7L, 8L, 9L, 31L, 32L, 33L,
               63L, 64L, 65L, 127L, 128L, 129L)
  ids <- sprintf("instrument:run:flowcell:read:%04d", seq_along(lengths))
  sequences <- vapply(lengths, boundary_sequence, character(1))
  qualities <- vapply(lengths, boundary_quality, character(1))

  packed <- rad_pack_reads(ids, sequences, qualities)
  materialized <- rad_materialize(packed)

  expect_s3_class(packed, "rad_reads")
  expect_true(data.table::is.data.table(packed))
  expect_s3_class(packed$seq_id, "rad_id")
  expect_s3_class(packed$seq, "rad_seq")
  expect_s3_class(packed$qual, "rad_phred")
  expect_identical(materialized$seq_id, ids)
  expect_identical(materialized$seq, sequences)
  expect_identical(materialized$qual, qualities)
  expect_equal(materialized$length, as.numeric(lengths))

  expected_scores <- lapply(lengths, function(n) {
    if (n) (seq_len(n) - 1L) %% 94L else integer()
  })
  expect_identical(rad_qual_values(packed), expected_scores)
  expect_identical(as.integer(packed$qual[[1L]]), integer())
  expect_identical(as.integer(packed$qual[[length(lengths)]]), expected_scores[[length(lengths)]])

  stats <- rad_storage_stats(packed)
  total_bases <- sum(lengths)
  expect_equal(stats$sequence_bases, total_bases)
  expect_equal(stats$quality_scores, total_bases)
  expect_equal(stats$sequence_packed_bytes, ceiling(total_bases * 2 / 8))
  expect_equal(stats$quality_packed_bytes, ceiling(total_bases * 7 / 8))
  expect_lt(stats$quality_packed_bytes, stats$quality_raw_bytes)
  expect_gt(stats$sequence_exception_bytes, 0)
})

test_that("the UUID and front-coded ID codecs are lossless", {
  uuids <- c(
    "00000000-0000-0000-0000-000000000000",
    "01234567-89ab-cdef-0123-456789abcdef",
    "ffffffff-ffff-ffff-ffff-ffffffffffff"
  )
  uuid_reads <- rad_pack_reads(uuids, rep("A", length(uuids)), rep("!", length(uuids)))
  uuid_stats <- rad_storage_stats(uuid_reads)

  expect_identical(as.character(uuid_reads$seq_id), uuids)
  expect_identical(uuid_stats$id_codec, "uuid16")
  expect_equal(uuid_stats$id_packed_bytes, 16 * length(uuids))
  expect_lt(uuid_stats$id_packed_bytes, uuid_stats$id_raw_bytes)

  front_ids <- sprintf("movie/flowcell/shared-prefix/%06d/0_%05d", 1:130, 10000 + 1:130)
  front_reads <- rad_pack_reads(front_ids, rep("AC", 130), rep("!!", 130))
  front_stats <- rad_storage_stats(front_reads)

  expect_identical(as.character(front_reads$seq_id), front_ids)
  expect_identical(as.character(front_reads$seq_id[130:1]), rev(front_ids))
  expect_identical(front_stats$id_codec, "front")
  expect_lt(front_stats$id_packed_bytes, front_stats$id_raw_bytes)

  mixed_case <- c(uuids[[1L]], toupper(uuids[[2L]]))
  mixed_reads <- rad_pack_reads(mixed_case, c("A", "C"), c("!", "!"))
  expect_identical(as.character(mixed_reads$seq_id), mixed_case)
  expect_identical(rad_storage_stats(mixed_reads)$id_codec, "front")
})

test_that("missing and empty qualities remain distinct", {
  packed <- rad_pack_reads(
    seq_id = c("missing", "empty", NA_character_, "present"),
    seq = c("AC", "", "G", "ACT"),
    qual = c(NA_character_, "", NA_character_, quality_string(c(0L, 20L, 93L))),
    seq_comment = c(NA_character_, "empty read", "missing id", NA_character_)
  )
  materialized <- rad_materialize(packed)

  expect_identical(materialized$seq_id, c("missing", "empty", NA_character_, "present"))
  expect_identical(materialized$seq_comment,
                   c(NA_character_, "empty read", "missing id", NA_character_))
  expect_identical(materialized$qual,
                   c(NA_character_, "", NA_character_, quality_string(c(0L, 20L, 93L))))
  expect_identical(rad_qual_values(packed),
                   list(NULL, integer(), NULL, c(0L, 20L, 93L)))
  expect_equal(rad_qmean(packed), c(NA, NA, NA, 113 / 3))
  expect_identical(rad_qmin(packed), c(NA_integer_, NA_integer_, NA_integer_, 0L))
  expect_equal(rad_qfraction(packed, 20L), c(NA, NA, NA, 2 / 3))
  expect_true(is.na(packed$qual[[1L]]))
  expect_false(is.na(packed$qual[[2L]]))

  no_quality <- rad_pack_reads(c("a", "b"), c("A", "CC"), qual = NULL)
  expect_identical(rad_materialize(no_quality)$qual, c(NA_character_, NA_character_))
  expect_identical(rad_qual_values(no_quality), list(NULL, NULL))
})

test_that("quality summaries are computed on numeric Phred values", {
  scores <- list(c(0L, 20L, 40L, 93L), c(19L, 20L, 21L), c(93L), integer())
  sequences <- vapply(scores, function(x) paste(rep("A", length(x)), collapse = ""), character(1))
  qualities <- vapply(scores, quality_string, character(1))
  packed <- rad_pack_reads(paste0("q", seq_along(scores)), sequences, qualities)

  expect_equal(rad_qmean(packed), c(38.25, 20, 93, NA))
  expect_identical(rad_qmin(packed), c(0L, 19L, 93L, NA_integer_))
  expect_equal(rad_qfraction(packed, 20L), c(0.75, 2 / 3, 1, NA))
  expect_equal(rad_qfraction(packed, 0L), c(1, 1, 1, NA))
  expect_equal(rad_qfraction(packed, 93L), c(0.25, 0, 1, NA))
  expect_equal(packed$mean_q, rad_qmean(packed))
  expect_identical(packed$min_q, rad_qmin(packed))
  expect_equal(packed$frac_q20, rad_qfraction(packed, 20L))
})

test_that("alternate Phred offsets preserve scores and can be selected on output", {
  scores <- c(0L, 1L, 20L, 40L, 62L)
  original <- quality_string(scores, offset = 64L)
  packed <- rad_pack_reads("offset64", strrep("A", length(scores)), original,
                           phred_offset = 64L)

  expect_identical(rad_qual_values(packed)[[1L]], scores)
  expect_identical(as.character(packed$qual, phred_offset = 64L), original)
  expect_identical(as.character(packed$qual, phred_offset = 33L), quality_string(scores))
  expect_identical(rad_materialize(packed, phred_offset = 64L)$qual, original)
})

test_that("empty packed tables retain the full schema", {
  packed <- rad_pack_reads(character(), character(), character(), character())
  expect_equal(nrow(packed), 0L)
  expect_identical(names(packed),
                   c("read_key", "seq_id", "seq_comment", "seq", "qual", "length",
                     "mean_q", "min_q", "frac_q20"))
  expect_identical(as.character(packed$seq_id), character())
  expect_identical(as.character(packed$seq), character())
  expect_identical(as.character(packed$qual), character())
  expect_equal(rad_storage_stats(packed)$reads, 0)
})
