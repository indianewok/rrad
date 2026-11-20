test_that("spatial export preserves sequence records and records assignment status", {
  work <- tempfile("spatial-export-"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- rad_test_fastq(file.path(work, "input.fq"), c("ACGT", "", "TGCA"),
    c("CB:Z:AAAACCCC UB:Z:TGCA ZZ:Z:keep", "CB:Z:AAAACCCC UB:Z:AAAA", "CB:Z:TTTTTTTT UB:Z:AAAA"),
    c("+instrument payload", "+empty", "+third"))
  map <- file.path(work, "map.tsv")
  writeLines(c("barcode\tarray_row\tarray_col\tin_tissue\tcell_id",
               "AAAACCCC\t5\t9\t1\tcellA"), map)
  output <- file.path(work, "result")
  result <- rad_spatial_reformat(input, output, "sample1", map = map, bin_size = 4L)
  expect_equal(result$records_read, 3)
  expect_equal(result$records_written, 3)
  expect_equal(result$counts$assigned, 2)
  lines <- rad_test_lines(result$files[["reads.fastq.gz"]])
  expect_match(lines[[1]], "CB:Z:s_004um_00002_00004-1", fixed = TRUE)
  expect_match(lines[[1]], "ZZ:Z:keep", fixed = TRUE)
  expect_identical(lines[c(3, 7, 11)], c("+instrument payload", "+empty", "+third"))
  expect_identical(lines[c(2, 6, 10)], c("ACGT", "", "TGCA"))
  audit <- utils::read.delim(gzfile(result$files[["assignments.tsv.gz"]]))
  expect_identical(audit$status, c("assigned", "assigned", "unknown_barcode"))
  strict <- rad_spatial_reformat(input, file.path(work, "strict"), "sample1", map = map,
                                 unit = "cell", tissue_only = TRUE, strict = TRUE)
  expect_equal(strict$records_written, 2)
  expect_match(rad_test_lines(strict$files[["reads.fastq.gz"]])[[1]], "CB:Z:cellA", fixed = TRUE)
})


test_that("spatial export validates maps and rolls back failed writes", {
  work <- tempfile("spatial-failure-"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- rad_test_fastq(file.path(work, "input.fq"), "ACGT", "CB:Z:AAAACCCC UB:Z:TTTT")
  x <- file.path(work, "x.tsv"); y <- file.path(work, "y.tsv")
  writeLines(c("barcode\tarray_col", "AAAA\t7"), x)
  writeLines(c("barcode\tarray_row", "CCCC\t3"), y)
  result <- rad_spatial_reformat(input, file.path(work, "axes"), "sample", bc1_axis = x, bc2_axis = y)
  expect_equal(result$counts$assigned, 1)
  expect_error(rad_spatial_reformat(input, file.path(work, "axes"), "sample", bc1_axis = x, bc2_axis = y), "new directory")
  writeLines(c("barcode\tarray_col", "AAAA\t7", "AAAA\t8"), x)
  target <- file.path(work, "invalid")
  expect_error(rad_spatial_reformat(input, target, "sample", bc1_axis = x, bc2_axis = y), "Duplicate")
  expect_false(dir.exists(target))
  writeLines(c("barcode\tarray_col", "AAAA\t7"), x)
  writeLines(c("@bad", "AAAA", "+", "II"), input)
  expect_error(rad_spatial_reformat(input, target, "sample", bc1_axis = x, bc2_axis = y), "complete four-line")
  expect_false(dir.exists(target))
  expect_false(any(grepl("rrad-demux-stage", list.files(work, all.files = TRUE))))
})


test_that("discovery and spatial export protect publication lock aliases", {
  work <- tempfile("rad-publication-input-"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- rad_test_fastq(file.path(work, ".result.rrad-demux-commit.lock"),
                           "ACGT", "CB:Z:AAAACCCC UB:Z:TTTT")
  Sys.chmod(input, "0644")
  mode <- file.info(input)$mode
  original <- readBin(input, "raw", n = 100000L)
  expect_error(rad_discover(input, work, "result", mode = "layout"),
               "publication lock|collid|same file")
  map <- file.path(work, "map.tsv")
  writeLines(c("barcode\tarray_row\tarray_col", "AAAACCCC\t0\t0"), map)
  expect_error(rad_spatial_reformat(input, file.path(work, "result"), "sample", map = map),
               "publication lock|collid|same file")
  expect_identical(readBin(input, "raw", n = 100000L), original)
  expect_identical(file.info(input)$mode, mode)
  expect_false(dir.exists(file.path(work, "result")))
})
