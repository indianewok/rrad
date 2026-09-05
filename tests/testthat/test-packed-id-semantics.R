test_that("duplicate compressed IDs group by their decoded identity", {
  ids <- c("duplicate", "solo", "duplicate", NA_character_,
           NA_character_, "", "", "last")
  packed <- rad_pack_reads(
    ids,
    rep("A", length(ids)),
    rep("!", length(ids))
  )

  grouped <- packed[, .(records = .N, rows = paste(.I, collapse = ",")), by = seq_id]
  expect_identical(as.character(grouped$seq_id),
                   c("duplicate", "solo", NA_character_, "", "last"))
  expect_identical(grouped$records, c(2L, 1L, 2L, 2L, 1L))
  expect_identical(grouped$rows, c("1,3", "2", "4,5", "6,7", "8"))
  expect_s3_class(grouped$seq_id, "rad_id")

  decoded_groups <- as.character(grouped$seq_id)
  expect_identical(grouped$records[is.na(decoded_groups)], 2L)
  expect_identical(grouped$records[match("", decoded_groups)], 2L)
  expect_identical(grouped$records[match("duplicate", decoded_groups)], 2L)
})

test_that("duplicated, anyDuplicated, and unique use ID text semantics", {
  ids <- c("duplicate", "solo", "duplicate", NA_character_,
           NA_character_, "", "", "last")
  packed <- rad_pack_reads(ids, rep("A", length(ids)), rep("!", length(ids)))
  handles <- packed$seq_id

  expect_identical(duplicated(handles), duplicated(ids))
  expect_identical(duplicated(handles, fromLast = TRUE), duplicated(ids, fromLast = TRUE))
  expect_identical(anyDuplicated(handles), anyDuplicated(ids))
  expect_identical(anyDuplicated(handles, fromLast = TRUE), anyDuplicated(ids, fromLast = TRUE))

  first_unique <- unique(handles)
  last_unique <- unique(handles, fromLast = TRUE)
  expect_s3_class(first_unique, "rad_id")
  expect_s3_class(last_unique, "rad_id")
  expect_identical(as.character(first_unique), unique(ids))
  expect_identical(as.character(last_unique), unique(ids, fromLast = TRUE))
})

test_that("compressed ID comparisons follow ordinary character semantics", {
  ids <- c("b", "a", "b", NA_character_, "")
  packed <- rad_pack_reads(ids, rep("A", length(ids)), rep("!", length(ids)))
  handles <- packed$seq_id

  expect_identical(handles == "b", ids == "b")
  expect_identical("b" == handles, "b" == ids)
  expect_identical(handles != "a", ids != "a")
  expect_identical(handles < "b", ids < "b")
  expect_identical(handles == handles[c(3L, 2L, 1L, 4L, 5L)],
                   ids == ids[c(3L, 2L, 1L, 4L, 5L)])
  expect_error(handles + 1, "not defined for compressed read IDs")
})

test_that("base sorting of compressed IDs is lexical and retains the store", {
  ids <- c("beta", "", "alpha", NA_character_, "aa", "a")
  handles <- rad_pack_reads(ids, rep("A", length(ids)), rep("!", length(ids)))$seq_id

  ascending <- sort(handles, na.last = TRUE)
  descending <- sort(handles, decreasing = TRUE, na.last = FALSE)
  expect_s3_class(ascending, "rad_id")
  expect_s3_class(descending, "rad_id")
  expect_identical(as.character(ascending), sort(ids, na.last = TRUE))
  expect_identical(as.character(descending),
                   sort(ids, decreasing = TRUE, na.last = FALSE))
  expect_identical(attr(ascending, ".rad_store"), attr(handles, ".rad_store"))
})

test_that("equivalent R string encodings receive the same compressed ID", {
  utf8 <- enc2utf8("\u00e9")
  latin1 <- iconv(utf8, from = "UTF-8", to = "latin1")
  Encoding(latin1) <- "latin1"
  expect_true(utf8 == latin1)

  packed <- rad_pack_reads(
    c(utf8, latin1), c("A", "C"), c("!", "5"), c(latin1, utf8)
  )
  expect_identical(unclass(packed$seq_id[1L]), unclass(packed$seq_id[2L]))
  expect_identical(unclass(packed$seq_comment[1L]),
                   unclass(packed$seq_comment[2L]))
  expect_identical(as.character(packed$seq_id), c(utf8, utf8))
  expect_identical(as.character(packed$seq_comment), c(utf8, utf8))
  expect_identical(packed[, .N, by = seq_id]$N, 2L)
})

test_that("the same ID compares and joins across independent stores", {
  left <- rad_pack_reads(
    c("shared", "left-only", "also-shared", NA_character_),
    rep("A", 4L), rep("!", 4L)
  )
  right <- rad_pack_reads(
    c("also-shared", "shared", "right-only", NA_character_),
    rep("C", 4L), rep("5", 4L)
  )

  expect_identical(
    left$seq_id == right$seq_id[c(2L, 3L, 1L, 4L)],
    c(TRUE, FALSE, TRUE, NA)
  )
  expect_identical(left$seq_id[c(1L, 3L)] == right$seq_id[c(2L, 1L)],
                   c(TRUE, TRUE))

  left_index <- data.table::data.table(
    seq_id = left$seq_id[1:3],
    left_value = c(10L, 20L, 30L)
  )
  right_index <- data.table::data.table(
    seq_id = right$seq_id[1:3],
    right_value = c("also", "shared", "right")
  )
  joined <- left_index[right_index, on = "seq_id", nomatch = NULL]

  expect_identical(as.character(joined$seq_id), c("also-shared", "shared"))
  expect_identical(joined$left_value, c(30L, 10L))
  expect_identical(joined$right_value, c("also", "shared"))
})
