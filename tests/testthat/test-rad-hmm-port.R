hmm_port_read <- function(path) {
  connection <- gzfile(path, "rt")
  on.exit(close(connection), add = TRUE)
  readLines(connection, warn = FALSE)
}

test_that("HMM stays opt-in and preserves calibrated columns on cache reuse", {
  work <- tempfile("rrad-hmm-cache-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- test_path("fixtures", "rad-core", "nanopore-smoke.fastq")
  ordinary <- rad_demux("nanopore_rapid_bc", input, file.path(work, "ordinary"),
                        threads = 1L, chunk_size = 1L, min_read_length = 0L)
  run <- function() rad_demux(
    "nanopore_rapid_bc", input, file.path(work, "hmm"),
    threads = 1L, chunk_size = 1L, min_read_length = 0L,
    concat_hmm = TRUE, reuse_cache = TRUE
  )
  first <- run()
  first_output <- hmm_port_read(first$artifacts$fastq)
  parameters <- utils::read.csv(first$artifacts$position_map,
                                check.names = FALSE, stringsAsFactors = FALSE)
  hmm_columns <- startsWith(names(parameters), "hmm_")
  expect_true(any(hmm_columns))
  expect_false(ordinary$stats$concat_hmm)
  expect_true(first$stats$concat_hmm)
  # One calibration read must not claim a reliable fitted HMM. Its guarded
  # fallback still retains the ordinary extractor's output.
  expect_match(first$log, "calibration guard failed", ignore.case = TRUE)
  expect_identical(first_output, hmm_port_read(ordinary$artifacts$fastq))
  expect_equal(first$stats$hmm_boundary_unresolved_parents, 0)
  expect_length(hmm_port_read(first$artifacts$boundary_unresolved), 0L)
  expect_false(grepl("different layout topology|has no layout element", first$log))

  cached <- run()
  cached_parameters <- utils::read.csv(cached$artifacts$position_map,
                                       check.names = FALSE, stringsAsFactors = FALSE)
  expect_identical(cached_parameters[, hmm_columns, drop = FALSE],
                   parameters[, hmm_columns, drop = FALSE])
  expect_identical(hmm_port_read(cached$artifacts$fastq), first_output)
  expect_false(grepl("Calibration on", cached$log, fixed = TRUE))
  expect_false(grepl("different layout topology|has no layout element", cached$log))
})

test_that("a fitted HMM runs through the R demux interface", {
  work <- tempfile("rrad-hmm-fitted-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  template <- readLines(test_path("fixtures", "rad-core", "nanopore-smoke.fastq"))
  cdna <- paste(rep("GACTCGACTAGCTACGATCGTAGCATGCTAGCTGACTAGCA", 10L),
                collapse = "")
  single <- sub("GATTACAGATTACA", cdna, template[[2L]], fixed = TRUE)
  # The first 200 complete molecules satisfy the HMM's calibration minimum.
  # A tandem parent follows them to exercise the non-single routing path.
  sequences <- c(rep(single, 200L), paste0(single, single))
  input <- file.path(work, "reads.fastq")
  records <- unlist(lapply(seq_along(sequences), function(i) {
    sequence <- sequences[[i]]
    c(paste0("@molecule", i), sequence, "+",
      paste(rep("I", nchar(sequence)), collapse = ""))
  }), use.names = FALSE)
  writeLines(records, input)
  result <- rad_demux(
    "nanopore_rapid_bc", input, file.path(work, "out"),
    threads = 1L, chunk_size = 50L, min_read_length = 0L,
    concat_hmm = TRUE, calibration_reads = 200L
  )
  expect_true(result$stats$concat_hmm)
  expect_equal(result$stats$total_reads, 201)
  expect_match(result$log, "Calibration on 200 reads", fixed = TRUE)
  expect_false(grepl("CALIBRATION GUARD FAILED", result$log, fixed = TRUE))
  expect_false(grepl("different layout topology|has no layout element", result$log))
  expect_match(result$log, "guard=0", fixed = TRUE)
  output <- hmm_port_read(result$artifacts$fastq)
  expect_equal(length(output), result$stats$records_serialized * 4)
  expect_gte(result$stats$records_serialized, 200)
  headers <- output[seq.int(1L, length(output), 4L)]
  expect_identical(sub("-F.*$", "", headers[seq_len(200L)]),
                   paste0("@molecule", seq_len(200L)))
  expect_equal(nchar(output[seq.int(2L, length(output), 4L)]),
               nchar(output[seq.int(4L, length(output), 4L)]))
})
