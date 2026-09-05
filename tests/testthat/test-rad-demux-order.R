test_that("multithreaded demux preserves input FASTQ order", {
  template <- readLines(
    test_path("fixtures", "rad-core", "nanopore-smoke.fastq"),
    warn = FALSE
  )
  ids <- sprintf("ordered-%03d", seq_len(24L))
  records <- unlist(lapply(ids, function(id) {
    record <- template
    record[[1L]] <- paste0("@", id)
    record
  }), use.names = FALSE)

  work <- tempfile("rad-core-parallel-order-")
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  dir.create(work)
  fastq <- file.path(work, "ordered.fastq")
  writeLines(records, fastq)

  serial <- rad_demux(
    layout = "nanopore_rapid_bc",
    fastq = fastq,
    out_dir = file.path(work, "serial"),
    output = "ordered",
    threads = 1L,
    chunk_size = length(ids),
    min_read_length = 0L
  )
  parallel_call <- function() {
    rad_demux(
      layout = "nanopore_rapid_bc",
      fastq = fastq,
      out_dir = file.path(work, "parallel"),
      output = "ordered",
      threads = 2L,
      chunk_size = length(ids),
      min_read_length = 0L
    )
  }
  if (isTRUE(rad_core_info()$build$openmp)) {
    parallel <- parallel_call()
  } else {
    parallel <- NULL
    expect_warning(
      parallel <- parallel_call(),
      "no OpenMP support"
    )
  }

  serial_output <- readLines(gzfile(serial$artifacts$fastq), warn = FALSE)
  parallel_output <- readLines(gzfile(parallel$artifacts$fastq), warn = FALSE)
  output_headers <- parallel_output[seq.int(1L, length(parallel_output), 4L)]
  output_ids <- sub("-[FR].*$", "", sub("^@", "", output_headers))

  expect_identical(output_ids, ids)
  expect_identical(parallel_output, serial_output)
})
