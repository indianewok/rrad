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
  expected <- paste0(paste(rad_test_lines(file$artifacts$fastq), collapse = "\n"), "\n")
  expect_identical(rawToChar(rawConnectionValue(stream)), expected)
  expect_false(file.exists(file.path(work, "stream.fq.gz")))
  expect_false(file.exists(file.path(work, "stream_dbg.fq.gz")))
  text <- character()
  text_stream <- textConnection("text", "w", local = TRUE)
  on.exit(close(text_stream), add = TRUE)
  rad_demux("nanopore_rapid_bc", input, work, "text", min_read_length = 0L,
            output_connection = text_stream)
  expect_identical(text, rad_test_lines(file$artifacts$fastq))
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
