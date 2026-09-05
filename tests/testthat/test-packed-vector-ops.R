test_that("c combines packed columns from the same store and field", {
  packed <- rad_pack_reads(
    c("id-1", "id-2", "id-3"),
    c("A", "CC", "GGG"),
    c("!", NA_character_, "III"),
    c("comment-1", NA_character_, "comment-3")
  )

  ids <- c(packed$seq_id[c(3L, 1L)], packed$seq_id[2L], packed$seq_id[NA_integer_])
  sequences <- c(packed$seq[c(3L, 1L)], packed$seq[2L])
  qualities <- c(packed$qual[c(3L, 1L)], packed$qual[2L])
  comments <- c(packed$seq_comment[c(3L, 1L)], packed$seq_comment[2L])

  expect_s3_class(ids, "rad_id")
  expect_s3_class(sequences, "rad_seq")
  expect_s3_class(qualities, "rad_phred")
  expect_s3_class(comments, "rad_id")
  expect_identical(attr(ids, ".rad_store"), attr(packed$seq_id, ".rad_store"))
  expect_identical(attr(sequences, ".rad_store"), attr(packed$seq, ".rad_store"))
  expect_identical(attr(qualities, ".rad_store"), attr(packed$qual, ".rad_store"))
  expect_identical(attr(ids, ".rad_field"), "id")
  expect_identical(attr(comments, ".rad_field"), "comment")

  expect_identical(as.character(ids),
                   c("id-3", "id-1", "id-2", NA_character_))
  expect_identical(as.character(sequences), c("GGG", "A", "CC"))
  expect_identical(as.character(qualities), c("III", "!", NA_character_))
  expect_identical(rad_qual_values(qualities), list(c(40L, 40L, 40L), 0L, NULL))
  expect_identical(as.character(comments),
                   c("comment-3", "comment-1", NA_character_))
})

test_that("c rejects cross-store, cross-class, and cross-field handles", {
  first <- rad_pack_reads(
    c("same", "first"), c("A", "CC"), c("!", "55"),
    c("first comment", "second comment")
  )
  second <- rad_pack_reads(
    c("same", "second"), c("G", "TT"), c("I", "~~"),
    c("first comment", "other comment")
  )

  expect_error(c(first$seq_id, second$seq_id),
               "only be combined from the same store and field")
  expect_error(c(first$seq, second$seq),
               "only be combined from the same store and field")
  expect_error(c(first$qual, second$qual),
               "only be combined from the same store and field")
  expect_error(c(first$seq_id, first$seq_comment),
               "only be combined from the same store and field")
  expect_error(c(first$seq_id, first$seq),
               "only be combined from the same store and field")
  expect_error(c(first$qual, bit64::as.integer64(1L)),
               "only be combined from the same store and field")
})

test_that("rep preserves packed column classes, stores, fields, and decoding", {
  packed <- rad_pack_reads(
    c("id-1", "id-2", "id-3"),
    c("A", "CC", "GGG"),
    c("!", NA_character_, "III"),
    c("comment-1", NA_character_, "comment-3")
  )

  ids <- rep(packed$seq_id[1:2], times = 2L)
  sequences <- rep(packed$seq[3:1], length.out = 5L)
  qualities <- rep(packed$qual[1:2], each = 2L)
  comments <- rep(packed$seq_comment[2:3], times = 2L)

  expect_s3_class(ids, "rad_id")
  expect_s3_class(sequences, "rad_seq")
  expect_s3_class(qualities, "rad_phred")
  expect_s3_class(comments, "rad_id")
  expect_identical(attr(ids, ".rad_store"), attr(packed$seq_id, ".rad_store"))
  expect_identical(attr(sequences, ".rad_store"), attr(packed$seq, ".rad_store"))
  expect_identical(attr(qualities, ".rad_store"), attr(packed$qual, ".rad_store"))
  expect_identical(attr(ids, ".rad_field"), "id")
  expect_identical(attr(comments, ".rad_field"), "comment")

  expect_identical(as.character(ids), c("id-1", "id-2", "id-1", "id-2"))
  expect_identical(as.character(sequences), c("GGG", "CC", "A", "GGG", "CC"))
  expect_identical(as.character(qualities),
                   c("!", "!", NA_character_, NA_character_))
  expect_identical(rad_qual_values(qualities), list(0L, 0L, NULL, NULL))
  expect_identical(as.character(comments),
                   c(NA_character_, "comment-3", NA_character_, "comment-3"))
})
