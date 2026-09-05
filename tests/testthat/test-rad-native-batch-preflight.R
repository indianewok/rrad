test_that("batch scan protects shared whitelist inputs before any run", {
  work <- tempfile("rad-native-batch-whitelist-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  input_one <- file.path(work, "one.fastq")
  input_two <- file.path(work, "two.fastq")
  writeLines(c("@one", "ACGTAACC", "+", "IIIIIIII"), input_one)
  writeLines(c("@two", "ACGTTTGA", "+", "IIIIIIII"), input_two)

  whitelist <- file.path(work, "protected.txt")
  writeLines(c("AACC", "TTGA"), whitelist)
  first_prefix <- file.path(work, "first")
  colliding_prefix <- sub("\\.txt$", "", whitelist)
  batch <- file.path(work, "batch.csv")
  utils::write.table(
    data.frame(
      input = c(input_one, input_two),
      output_prefix = c(first_prefix, colliding_prefix)
    ),
    batch,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE
  )

  expect_error(
    rad_scan_wl(
      batch_csv = batch,
      adapter_seq = "ACGT",
      barcode_length = 4L,
      whitelist = whitelist,
      max_error = 0
    ),
    "reference whitelist"
  )

  expect_identical(readLines(whitelist), c("AACC", "TTGA"))
  expect_false(any(file.exists(c(
    paste0(first_prefix, ".csv"),
    paste0(first_prefix, ".txt"),
    paste0(first_prefix, "_scan_wl.log")
  ))))
})

test_that("batch scan detects hard-link aliases before any run", {
  work <- tempfile("rad-native-batch-hardlink-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  input_one <- file.path(work, "one.fastq")
  input_two <- file.path(work, "two.fastq")
  writeLines(c("@one", "ACGTAACC", "+", "IIIIIIII"), input_one)
  writeLines(c("@two", "ACGTTTGA", "+", "IIIIIIII"), input_two)

  whitelist <- file.path(work, "protected.txt")
  hard_link <- file.path(work, "alias.txt")
  writeLines(c("AACC", "TTGA"), whitelist)
  if (!isTRUE(file.link(whitelist, hard_link))) {
    skip("the test filesystem does not support hard links")
  }

  first_prefix <- file.path(work, "first")
  colliding_prefix <- sub("\\.txt$", "", hard_link)
  batch <- file.path(work, "batch.csv")
  utils::write.table(
    data.frame(
      input = c(input_one, input_two),
      output_prefix = c(first_prefix, colliding_prefix)
    ),
    batch,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE
  )

  expect_error(
    rad_scan_wl(
      batch_csv = batch,
      adapter_seq = "ACGT",
      barcode_length = 4L,
      whitelist = whitelist,
      max_error = 0
    ),
    "existing file alias"
  )

  expect_identical(readLines(whitelist), c("AACC", "TTGA"))
  expect_false(any(file.exists(c(
    paste0(first_prefix, ".csv"),
    paste0(first_prefix, ".txt"),
    paste0(first_prefix, "_scan_wl.log")
  ))))
})
