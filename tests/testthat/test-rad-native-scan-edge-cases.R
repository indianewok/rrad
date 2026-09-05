scan_edge_write_fastq <- function(path, sequences) {
  ids <- paste0("scan-edge-", seq_along(sequences))
  qualities <- vapply(
    sequences,
    function(sequence) paste(rep("I", nchar(sequence)), collapse = ""),
    character(1)
  )
  writeLines(as.vector(rbind(paste0("@", ids), sequences, "+", qualities)),
             path)
}

scan_edge_read_calls <- function(result) {
  readLines(result$artifacts$whitelist, warn = FALSE)
}

scan_edge_revcomp <- function(sequence) {
  paste0(rev(strsplit(chartr("ACGT", "TGCA", sequence), "")[[1L]]),
         collapse = "")
}

test_that("fuzzy scan extraction follows insertion and deletion endpoints", {
  work <- tempfile("rad-native-scan-indels-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "reads.fastq")

  # AGTCCGTA with one inserted A spans nine read bases; with one C deleted it
  # spans seven. In both cases the barcode starts at the native match end.
  scan_edge_write_fastq(
    input,
    c("AGTCACGTAGACA", "AGTCGTATTGC")
  )

  result <- rad_scan_wl(
    input,
    file.path(work, "cells"),
    adapter_seq = "AGTCCGTA",
    barcode_length = 4L,
    max_error = 0.13,
    threads = 1L,
    selection = "above_floor"
  )

  expect_equal(result$stats$reads_processed, 2)
  expect_equal(result$stats$total_extractions, 2)
  expect_setequal(scan_edge_read_calls(result), c("GACA", "TTGC"))
})

test_that("scan rejects truncated forward and reverse barcode windows", {
  work <- tempfile("rad-native-scan-truncated-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "reads.fastq")
  adapter <- "AGTCCGTA"

  scan_edge_write_fastq(
    input,
    c(
      paste0(adapter, "GAC"),
      scan_edge_revcomp(paste0(adapter, "TGC"))
    )
  )

  result <- rad_scan_wl(
    input,
    file.path(work, "cells"),
    adapter_seq = adapter,
    barcode_length = 4L,
    max_error = 0,
    threads = 1L,
    selection = "above_floor"
  )

  expect_equal(result$stats$reads_processed, 2)
  expect_equal(result$stats$total_extractions, 0)
  expect_length(scan_edge_read_calls(result), 0)
})

test_that("reverse scan margins mirror the documented forward window", {
  work <- tempfile("rad-native-scan-reverse-margins-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "reads.fastq")
  adapter <- "AGTCCGTA"

  # In adapter orientation: skip TT (left margin), extract GACA (barcode),
  # then include CG (right-margin extension). Reverse-complement the complete
  # molecule so the native reverse-strand branch must recover GACACG.
  oriented <- paste0(adapter, "TT", "GACA", "CG")
  scan_edge_write_fastq(input, scan_edge_revcomp(oriented))

  result <- rad_scan_wl(
    input,
    file.path(work, "cells"),
    adapter_seq = adapter,
    barcode_length = 4L,
    left_margin = 2L,
    right_margin = 2L,
    max_error = 0,
    threads = 1L,
    selection = "above_floor"
  )

  expect_equal(result$stats$total_extractions, 1)
  expect_identical(scan_edge_read_calls(result), "GACACG")
})

test_that("scan publication lock cannot alias a whitelist input", {
  input <- test_path("fixtures", "rad-core", "tiny.fastq")
  work <- tempfile("rad-native-scan-lock-input-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  prefix <- file.path(work, "cells")
  lock_path <- file.path(work, ".cells.rrad-commit.lock")
  writeLines(c("AACC", "TTGA"), lock_path)
  before <- readBin(lock_path, "raw", n = file.info(lock_path)$size)

  expect_error(
    rad_scan_wl(
      input,
      prefix,
      adapter_seq = "ACGTACGT",
      barcode_length = 4L,
      whitelist = lock_path,
      threads = 1L
    ),
    "scan publication lock would overwrite reference whitelist"
  )
  expect_identical(
    readBin(lock_path, "raw", n = file.info(lock_path)$size),
    before
  )
})
