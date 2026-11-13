discovery_fastq <- function(path, sequences, comments = rep("", length(sequences)),
                            plus = rep("+", length(sequences))) {
  records <- unlist(lapply(seq_along(sequences), function(i) c(
    paste0("@read", i, if (nzchar(comments[[i]])) paste0(" ", comments[[i]]) else ""),
    sequences[[i]], plus[[i]], strrep("I", nchar(sequences[[i]])))), use.names = FALSE)
  writeLines(records, path)
  invisible(path)
}

discovery_lines <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  readLines(con, warn = FALSE)
}

test_that("catalog discovery counts independent exact reads and retains ties", {
  work <- tempfile("discovery-"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  primer <- "CTACACGACGCTCTTCCGATCT"
  input <- discovery_fastq(file.path(work, "input.fq"),
                           paste0(primer, c("ACGTACGT", "ACGTACGT", "TTGGCCAA")))
  a <- file.path(work, "a.txt"); b <- file.path(work, "b.txt")
  writeLines(c("ACGTACGT", "ACGTACGT"), a)
  writeLines("TTGGCCAA", b)
  run <- function(output, ...) rad_discover(input, work, output,
    pattern = primer, barcode_length = 8L, catalogs = c(a = a, b = b),
    max_error = 0, ...)
  result <- run("match")
  expect_true(result$complete)
  expect_identical(result$whitelist_match$status, "matched")
  expect_equal(result$whitelist_match$matched_reads, 2)
  expect_equal(result$whitelists$catalog_size, c(1, 1))
  expect_equal(result$whitelists$matched_reads, c(2, 1))
  writeLines("ACGTACGT", b)
  tied <- run("tie")
  expect_false(tied$complete)
  expect_identical(tied$whitelist_match$status, c("tie", "tie"))
  writeLines("GGGGAAAA", a); writeLines("CCCCAAAA", b)
  absent <- run("absent")
  expect_false(absent$complete)
  expect_identical(absent$whitelist_match$status, "no_hits")
  expect_equal(absent$whitelists$matched_reads, c(0, 0))
})

test_that("catalog discovery supports masked captures, long barcodes and sampling limits", {
  work <- tempfile("discovery-capture-"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  linker <- "TCTTCAGCGTTCCCGAGA"
  input <- discovery_fastq(file.path(work, "input.fq"),
    rep(paste0("ACGT", linker, "TGCA"), 3))
  catalog <- file.path(work, "catalog.txt"); writeLines("ACGTTGCA", catalog)
  result <- rad_discover(input, work, "masked", pattern = paste0("NNNN", linker, "NNNN"),
                        catalogs = catalog, max_error = 0, sample_reads = 2L)
  expect_equal(result$sample_reads, 2)
  expect_equal(result$whitelist_match$matched_reads, 2)
  primer <- "CTACACGACGCTCTTCCGATCT"
  barcode <- paste0(strrep("ACGT", 15), "TTGA")
  discovery_fastq(input, rep(paste0(primer, barcode), 3))
  writeLines(barcode, catalog)
  long <- rad_discover(input, work, "long", pattern = primer,
                       barcode_length = 64L, catalogs = catalog, max_error = 0,
                       sample_bases = 2L * nchar(paste0(primer, barcode)))
  expect_equal(long$sample_reads, 2)
  expect_equal(long$whitelist_match$matched_reads, 2)
})

test_that("discovery handles incomplete layouts and protects outputs", {
  work <- tempfile("discovery-layout-"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- discovery_fastq(file.path(work, "input.fq"), rep("CTACACGACGCTCTTCCGATCTACGTACGT", 4))
  result <- rad_discover(input, work, "partial", mode = "layout", max_error = 0)
  expect_true(result$success)
  expect_false(result$layout_complete)
  expect_true(file.exists(result$files[["layout.partial.csv"]]))
  old <- readBin(result$files[["layout.partial.csv"]], "raw", n = 100000L)
  expect_error(rad_discover(input, work, "partial", mode = "layout"), "already exists")
  expect_identical(readBin(result$files[["layout.partial.csv"]], "raw", n = 100000L), old)
  expect_error(rad_discover(input, work, "invalid", pattern = "NNN", catalogs = input), "literal anchor")
  expect_false(any(grepl("^invalid", list.files(work))))
  malformed <- file.path(work, "bad.fq")
  writeLines(c("@bad", "AAAA", "+", "II"), malformed)
  expect_error(rad_discover(malformed, work, "bad", mode = "layout"), "malformed|parser")
  expect_false(file.exists(file.path(work, "bad_discover.log")))
})

test_that("spatial export preserves sequence records and records assignment status", {
  work <- tempfile("spatial-export-"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- discovery_fastq(file.path(work, "input.fq"), c("ACGT", "", "TGCA"),
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
  lines <- discovery_lines(result$files[["reads.fastq.gz"]])
  expect_match(lines[[1]], "CB:Z:s_004um_00002_00004-1", fixed = TRUE)
  expect_match(lines[[1]], "ZZ:Z:keep", fixed = TRUE)
  expect_identical(lines[c(3, 7, 11)], c("+instrument payload", "+empty", "+third"))
  expect_identical(lines[c(2, 6, 10)], c("ACGT", "", "TGCA"))
  audit <- utils::read.delim(gzfile(result$files[["assignments.tsv.gz"]]))
  expect_identical(audit$status, c("assigned", "assigned", "unknown_barcode"))
  strict <- rad_spatial_reformat(input, file.path(work, "strict"), "sample1", map = map,
                                 unit = "cell", tissue_only = TRUE, strict = TRUE)
  expect_equal(strict$records_written, 2)
  expect_match(discovery_lines(strict$files[["reads.fastq.gz"]])[[1]], "CB:Z:cellA", fixed = TRUE)
})

test_that("spatial export validates maps and rolls back failed writes", {
  work <- tempfile("spatial-failure-"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- discovery_fastq(file.path(work, "input.fq"), "ACGT", "CB:Z:AAAACCCC UB:Z:TTTT")
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

test_that("demux connection output matches ordered file output", {
  work <- tempfile("demux-stream-"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- test_path("fixtures", "rad-core", "nanopore-smoke.fastq")
  file <- rad_demux("nanopore_rapid_bc", input, work, "file", min_read_length = 0L)
  stream <- rawConnection(raw(), "wb")
  on.exit(close(stream), add = TRUE)
  result <- rad_demux("nanopore_rapid_bc", input, work, "stream", min_read_length = 0L,
                      output_connection = stream, write_debug = TRUE)
  expect_true(result$stats$streamed)
  expect_identical(result$artifacts$fastq, "")
  expect_equal(result$stats$records_serialized, file$stats$records_serialized)
  expected <- paste0(paste(discovery_lines(file$artifacts$fastq), collapse = "\n"), "\n")
  expect_identical(rawToChar(rawConnectionValue(stream)), expected)
  expect_false(file.exists(file.path(work, "stream.fq.gz")))
  expect_false(file.exists(file.path(work, "stream_dbg.fq.gz")))
  text <- character()
  text_stream <- textConnection("text", "w", local = TRUE)
  on.exit(close(text_stream), add = TRUE)
  rad_demux("nanopore_rapid_bc", input, work, "text", min_read_length = 0L,
            output_connection = text_stream)
  expect_identical(text, discovery_lines(file$artifacts$fastq))
  saved <- readBin(file$artifacts$fastq, "raw", n = 100000L)
  rad_demux("nanopore_rapid_bc", input, work, "file", min_read_length = 0L,
            output_connection = stream)
  expect_identical(readBin(file$artifacts$fastq, "raw", n = 100000L), saved)
})

test_that("a failed connection callback leaves existing reports intact", {
  work <- tempfile("demux-stream-failure-"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- test_path("fixtures", "rad-core", "nanopore-smoke.fastq")
  prior <- rad_demux("nanopore_rapid_bc", input, work, "same", min_read_length = 0L)
  paths <- list.files(work, full.names = TRUE)
  paths <- paths[!dir.exists(paths)]
  before <- lapply(paths, readBin, what = "raw", n = 1000000L)
  native <- get("rad_demux_cpp", envir = asNamespace("rrad"))
  calls <- 0L
  local_mocked_bindings(rad_demux_cpp = function(options) {
    options$output_callback <- function(text) {
      calls <<- calls + 1L
      stop("deliberate destination failure")
    }
    native(options)
  }, .package = "rrad")
  connection <- rawConnection(raw(), "wb")
  on.exit(close(connection), add = TRUE)
  expect_error(rad_demux("nanopore_rapid_bc", input, work, "same",
                        min_read_length = 0L, chunk_size = 1L,
                        output_connection = connection), "deliberate destination failure")
  expect_equal(calls, 1L)
  expect_identical(lapply(paths, readBin, what = "raw", n = 1000000L), before)
  expect_false(any(grepl("rrad-demux-stage", list.files(work, all.files = TRUE))))
  expect_true(isOpen(connection, "w"))
})

test_that("discovery and spatial export protect publication lock aliases", {
  work <- tempfile("rad-publication-input-"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- discovery_fastq(file.path(work, ".result.rrad-demux-commit.lock"),
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
