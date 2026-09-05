test_that("packed handle missingness is exact for every column type", {
  packed <- rad_pack_reads(
    seq_id = c("present", NA_character_, ""),
    seq = c("A", "C", ""),
    qual = c("!", NA_character_, ""),
    seq_comment = c(NA_character_, "present comment", "")
  )

  expect_identical(is.na(packed$seq_id), c(FALSE, TRUE, FALSE))
  expect_identical(is.na(packed$seq_comment), c(TRUE, FALSE, FALSE))
  expect_identical(is.na(packed$seq), c(FALSE, FALSE, FALSE))
  expect_identical(is.na(packed$qual), c(FALSE, TRUE, FALSE))
  expect_true(anyNA(packed$seq_id))
  expect_true(anyNA(packed$seq_comment))
  expect_false(anyNA(packed$seq))
  expect_true(anyNA(packed$qual))
  expect_true(anyNA(packed$qual, recursive = TRUE))

  reordered <- packed[c(3L, 1L, 2L)]
  expect_identical(is.na(reordered$seq_id), c(FALSE, FALSE, TRUE))
  expect_identical(is.na(reordered$qual), c(FALSE, FALSE, TRUE))
  expect_identical(as.character(reordered[is.na(seq_id)]$seq), "C")
  expect_identical(as.character(reordered[is.na(qual)]$seq_id), NA_character_)
})

test_that("raw rbindlist cannot silently decode handles from different stores", {
  first <- rad_pack_reads(c("a1", "a2"), c("A", "CC"), c("!", "55"))
  second <- rad_pack_reads(c("b1", "b2"), c("GGG", "TTTT"), c("III", "~~~~"))

  combined <- data.table::rbindlist(list(first, second), use.names = TRUE)
  expect_s3_class(combined$seq_id, "rad_id")
  expect_identical(as.character(combined$seq_id[1:2]), c("a1", "a2"))
  expect_error(as.character(combined$seq_id),
               "packed string handle is absent from its backing RAD store")
  expect_error(as.character(combined$seq),
               "handles from different RAD stores.*rad_rbind_reads")
  expect_error(rad_qual_values(combined$qual),
               "handles from different RAD stores.*rad_rbind_reads")
  class(combined) <- c("rad_reads", "data.table", "data.frame")
  expect_error(rad_validate(combined, deep = TRUE),
               "handles from different RAD stores.*rad_rbind_reads")
})

test_that("rad_rbind_reads safely repacks stores and preserves extra columns", {
  first <- rad_pack_reads(
    c("a1", "a2"), c("A", "CC"), c("!", NA_character_),
    c("first", NA_character_)
  )
  first[, batch := "A"]
  first[, first_only := c(10L, 11L)]

  second <- rad_pack_reads(c("b1", "b2"), c("GGG", "TTTT"), c("III", "~~~~"))
  second[, batch := "B"]
  second[, second_only := c(TRUE, FALSE)]

  combined <- rad_rbind_reads(first, second)
  expect_true(isTRUE(rad_validate(combined, deep = TRUE)))
  expect_identical(materialized_payload(combined), data.frame(
    seq_id = c("a1", "a2", "b1", "b2"),
    seq_comment = c("first", NA_character_, NA_character_, NA_character_),
    seq = c("A", "CC", "GGG", "TTTT"),
    qual = c("!", NA_character_, "III", "~~~~"),
    stringsAsFactors = FALSE
  ))
  expect_identical(combined$batch, c("A", "A", "B", "B"))
  expect_identical(combined$first_only, c(10L, 11L, NA_integer_, NA_integer_))
  expect_identical(combined$second_only, c(NA, NA, TRUE, FALSE))
  expect_identical(bit64::as.character.integer64(combined$read_key), as.character(1:4))
  expect_equal(rad_storage_stats(combined)$reads, 4)
  expect_false(rad_storage_stats(combined)$store_id %in%
                 c(rad_storage_stats(first)$store_id, rad_storage_stats(second)$store_id))

  from_list <- rad_rbind_reads(list(first, second))
  expect_identical(materialized_payload(from_list), materialized_payload(combined))
  expect_true(isTRUE(rad_validate(from_list, deep = TRUE)))
})

test_that("rad_rbind_reads handles empty input and rejects non-packed input", {
  empty <- rad_rbind_reads()
  expect_s3_class(empty, "rad_reads")
  expect_equal(nrow(empty), 0L)
  expect_true(isTRUE(rad_validate(empty, deep = TRUE)))

  packed <- rad_pack_reads("r", "A", "!")
  expect_error(rad_rbind_reads(packed, data.frame(seq_id = "x")),
               "every input must be a rad_reads")
})

test_that("rad_compact releases filtered payload while retaining table data", {
  ids <- sprintf("compact/shared-prefix/%04d", 1:20)
  sequences <- vapply(1:20, function(i) strrep(c("A", "N", "C", "R")[[i %% 4L + 1L]], i),
                      character(1))
  scores <- lapply(1:20, function(i) rep((i * 3L) %% 63L, i))
  qualities <- vapply(scores, quality_string, character(1), offset = 64L)
  original <- rad_pack_reads(ids, sequences, qualities, phred_offset = 64L)
  original[, source_row := seq_len(.N)]
  original_stats <- rad_storage_stats(original)

  kept_rows <- c(19L, 2L, 11L)
  filtered <- original[kept_rows]
  expect_equal(rad_storage_stats(filtered)$reads, 20)
  compacted <- rad_compact(filtered)
  compacted_stats <- rad_storage_stats(compacted)

  expect_true(isTRUE(rad_validate(compacted, deep = TRUE)))
  expect_identical(materialized_payload(compacted, phred_offset = 64L),
                   materialized_payload(filtered, phred_offset = 64L))
  expect_identical(compacted$source_row, kept_rows)
  expect_equal(compacted_stats$reads, length(kept_rows))
  expect_equal(compacted_stats$phred_offset, 64)
  expect_false(compacted_stats$store_id == original_stats$store_id)
  expect_lt(compacted_stats$sequence_bases, original_stats$sequence_bases)
  expect_lt(compacted_stats$quality_scores, original_stats$quality_scores)
})

test_that("rad_save and rad_load rebuild a valid independent store", {
  original <- rad_pack_reads(
    c("save1", NA_character_, "save3"),
    c("ACGTN", "", "G-RY"),
    c(quality_string(c(0L, 20L, 40L, 62L, 1L), 64L), "", NA_character_),
    c("first comment", "empty", NA_character_),
    phred_offset = 64L
  )
  original[, sample := factor(c("alpha", "beta", "alpha"))]
  original[, selected := c(TRUE, FALSE, NA)]
  path <- tempfile(fileext = ".rrad.rds")
  on.exit(unlink(path), add = TRUE)

  expect_identical(invisible(rad_save(original, path)), normalizePath(path, mustWork = TRUE))
  restored <- rad_load(path)

  expect_true(isTRUE(rad_validate(restored, deep = TRUE)))
  expect_identical(materialized_payload(restored, phred_offset = 64L),
                   materialized_payload(original, phred_offset = 64L))
  expect_identical(restored$sample, original$sample)
  expect_identical(restored$selected, original$selected)
  expect_equal(rad_storage_stats(restored)$phred_offset, 64)
  expect_false(rad_storage_stats(restored)$store_id == rad_storage_stats(original)$store_id)

  rm(original)
  invisible(gc())
  expect_identical(as.character(restored$seq), c("ACGTN", "", "G-RY"))
  expect_identical(is.na(restored$qual), c(FALSE, FALSE, TRUE))
})

test_that("rad_save preserves an existing destination on serialization failure", {
  packed <- rad_pack_reads("save", "ACGT", "IIII")
  path <- tempfile(fileext = ".rrad.rds")
  on.exit(unlink(path), add = TRUE)
  sentinel <- charToRaw("existing-valid-content")
  writeBin(sentinel, path)

  expect_error(rad_save(packed, path, compress = "not-a-codec"))
  expect_identical(
    readBin(path, what = "raw", n = length(sentinel)),
    sentinel
  )
  expect_length(
    list.files(dirname(path), pattern = "^\\.rrad-stage-"),
    0L
  )
})

test_that("rad_load rejects non-RAD and unsupported serialized payloads", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  saveRDS(data.frame(seq_id = "r"), path)
  expect_error(rad_load(path), "not a supported RAD read store")

  saveRDS(list(format = "rrad_reads", version = 2L, data = data.frame()), path)
  expect_error(rad_load(path), "not a supported RAD read store")

  expect_error(rad_save(data.frame(), path), "x must be a rad_reads")
})

test_that("rad_save rejects additional packed handle columns", {
  packed <- rad_pack_reads(
    c("s1", "s2"), c("AC", "GT"), c("!!", "II"),
    c("one", "two")
  )
  packed[, sequence_alias := seq]
  packed[, id_alias := seq_id]
  packed[, quality_alias := qual]
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  expect_error(
    rad_save(packed, path),
    "additional columns cannot be saved.*sequence_alias.*id_alias.*quality_alias"
  )
  expect_false(file.exists(path))

  nested <- rad_pack_reads(c("n1", "n2"), c("A", "C"), c("!", "5"))
  nested[, packed_list := list(list(seq[1L]), list(seq[2L]))]
  expect_error(rad_save(nested, path),
               "additional columns cannot be saved.*packed_list")
  expect_false(file.exists(path))

  captured <- rad_pack_reads(c("c1", "c2"), c("A", "C"), c("!", "5"))
  captured_handle <- captured$seq[1L]
  captured[, callback := list(
    function() as.character(captured_handle),
    function() as.character(captured_handle)
  )]
  expect_error(rad_save(captured, path),
               "additional columns cannot be saved.*callback")
  expect_false(file.exists(path))
})

test_that("repacking operations reject misaligned input instead of legitimizing it", {
  packed <- rad_pack_reads(
    c("r1", "r2", "r3"), c("A", "CC", "GGG"), c("!", "55", "III")
  )
  damaged <- data.table::copy(packed)
  data.table::set(damaged, j = "seq", value = packed$seq[c(3L, 2L, 1L)])
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  expect_error(rad_compact(damaged), "seq_id and sequence handles are not aligned")
  expect_error(rad_rbind_reads(damaged),
               "seq_id and sequence handles are not aligned")
  expect_error(rad_save(damaged, path), "seq_id and sequence handles are not aligned")
  expect_false(file.exists(path))
})

test_that("deep validation checks lengths and aligned sequence-quality handles", {
  packed <- rad_pack_reads(
    c("v1", "v2", "v3"), c("A", "CC", "GGG"), c("!", "55", NA_character_)
  )
  expect_true(isTRUE(rad_validate(packed)))
  expect_true(isTRUE(rad_validate(packed, deep = TRUE)))

  bad_length <- data.table::copy(packed)
  data.table::set(bad_length, i = 1L, j = "length", value = 999)
  expect_true(isTRUE(rad_validate(bad_length)))
  expect_error(rad_validate(bad_length, deep = TRUE),
               "stored sequence lengths do not match")

  bad_alignment <- data.table::copy(packed)
  data.table::set(bad_alignment, j = "qual",
                  value = packed$qual[c(2L, 1L, 3L)])
  expect_error(rad_validate(bad_alignment, deep = TRUE),
               "quality and sequence handles are not aligned")

  other <- rad_pack_reads(c("o1", "o2", "o3"), c("A", "CC", "GGG"), c("!", "55", "III"))
  different_store <- data.table::copy(packed)
  data.table::set(different_store, j = "seq", value = other$seq)
  expect_error(rad_validate(different_store), "do not share one backing store")
})

test_that("deep validation detects tampered cached quality summaries", {
  packed <- rad_pack_reads(
    c("q1", "q2", "q3", "q4"),
    c("AC", "GGG", "T", ""),
    c(quality_string(c(0L, 40L)), quality_string(c(19L, 20L, 21L)),
      NA_character_, "")
  )
  expect_true(isTRUE(rad_validate(packed, deep = TRUE)))

  bad_mean <- data.table::copy(packed)
  data.table::set(bad_mean, i = 1L, j = "mean_q",
                  value = bad_mean$mean_q[1L] + 1)
  expect_true(isTRUE(rad_validate(bad_mean)))
  expect_error(rad_validate(bad_mean, deep = TRUE),
               "stored mean_q values do not match packed qualities")

  bad_min <- data.table::copy(packed)
  data.table::set(bad_min, i = 2L, j = "min_q",
                  value = bad_min$min_q[2L] + 1L)
  expect_true(isTRUE(rad_validate(bad_min)))
  expect_error(rad_validate(bad_min, deep = TRUE),
               "stored min_q values do not match packed qualities")

  fractional_min <- data.table::copy(packed)
  fractional_value <- as.numeric(fractional_min$min_q)
  fractional_value[2L] <- packed$min_q[2L] + 0.5
  data.table::set(fractional_min, j = "min_q", value = fractional_value)
  expect_error(rad_validate(fractional_min, deep = TRUE),
               "stored min_q values do not match packed qualities")

  bad_fraction <- data.table::copy(packed)
  data.table::set(bad_fraction, i = 1L, j = "frac_q20", value = 0.25)
  expect_true(isTRUE(rad_validate(bad_fraction)))
  expect_error(rad_validate(bad_fraction, deep = TRUE),
               "stored frac_q20 values do not match packed qualities")

  missing_mean <- data.table::copy(packed)
  data.table::set(missing_mean, i = 3L, j = "mean_q", value = 0)
  expect_error(rad_validate(missing_mean, deep = TRUE),
               "stored mean_q values do not match packed qualities")

  empty_min <- data.table::copy(packed)
  data.table::set(empty_min, i = 4L, j = "min_q", value = 0L)
  expect_error(rad_validate(empty_min, deep = TRUE),
               "stored min_q values do not match packed qualities")
})

test_that("rad_validate rejects damaged schemas and handle classes", {
  expect_error(rad_validate(data.table::data.table()),
               "x must be a rad_reads data.table")

  packed <- rad_pack_reads("r", "A", "!")
  missing_core <- data.table::copy(packed)
  data.table::set(missing_core, j = "qual", value = NULL)
  expect_error(rad_validate(missing_core), "core columns are missing")

  bad_class <- data.table::copy(packed)
  bad_sequence <- bad_class$seq
  class(bad_sequence) <- "integer64"
  data.table::set(bad_class, j = "seq", value = bad_sequence)
  expect_error(rad_validate(bad_class), "handle column classes are invalid")

  missing_store <- data.table::copy(packed)
  missing_store_sequence <- missing_store$seq
  attr(missing_store_sequence, ".rad_store") <- NULL
  data.table::set(missing_store, j = "seq", value = missing_store_sequence)
  expect_error(rad_validate(missing_store), "do not share one backing store")
})
