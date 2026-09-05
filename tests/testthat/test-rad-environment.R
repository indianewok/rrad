test_that("native demultiplexing is safe when HOME is unset", {
  rscript <- file.path(R.home("bin"), "Rscript")
  library_path <- normalizePath(.libPaths()[[1L]], mustWork = TRUE)
  fastq <- normalizePath(
    test_path("fixtures", "rad-core", "nanopore-smoke.fastq"),
    mustWork = TRUE
  )
  out_dir <- tempfile("rrad-no-home-")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  quote_r <- function(value) encodeString(value, quote = "\"")
  code <- paste0(
    "Sys.unsetenv('HOME'); ",
    "library(rrad, lib.loc = ", quote_r(library_path), "); ",
    "result <- rad_demux(",
    "layout = 'nanopore_rapid_bc', ",
    "fastq = ", quote_r(fastq), ", ",
    "out_dir = ", quote_r(out_dir), ", ",
    "output = 'smoke', threads = 1L, chunk_size = 1L, ",
    "max_reads = 1L, min_read_length = 0L); ",
    "stopifnot(inherits(result, 'rad_demux_result'), ",
    "identical(result$status, 'success'))"
  )

  output <- suppressWarnings(system2(
    rscript,
    c("--vanilla", "-e", shQuote(code)),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status", exact = TRUE)
  if (is.null(status)) status <- 0L

  expect_identical(
    as.integer(status),
    0L,
    info = paste(output, collapse = "\n")
  )
})
