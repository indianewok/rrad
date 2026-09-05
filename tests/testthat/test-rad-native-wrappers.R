wrapper_core_info <- function(features = c("demux", "scan-wl", "reformat"),
                              openmp = TRUE,
                              source_digest = paste(rep("c", 64L), collapse = "")) {
  list(
    core_version = "1.0.2",
    core_api_version = 1L,
    source_commit = paste(rep("a", 40L), collapse = ""),
    embedded_source_digest = source_digest,
    resource_schema = 1L,
    resource_digest = paste(rep("b", 64L), collapse = ""),
    features = features,
    build = list(
      cxx_standard = 17L,
      openmp = openmp,
      resources_available = TRUE,
      resource_dir = "/test/rrad/rad/resources"
    ),
    vendored = TRUE
  )
}

wrapper_scan_result <- function(options, core = wrapper_core_info()) {
  option_or <- function(value, default) {
    if (is.null(value)) default else value
  }
  requested_threads <- option_or(options$threads, 1L)
  effective_threads <- if (isTRUE(core$build$openmp)) requested_threads else 1L
  list(
    success = TRUE,
    backend = "embedded",
    core = core,
    input = options$input,
    output_prefix = options$output_prefix,
    config = list(
      input = options$input,
      output_prefix = options$output_prefix,
      adapter_seq = option_or(options$adapter_seq, ""),
      barcode_length = option_or(options$barcode_length, 16L),
      left_margin = option_or(options$left_margin, 0L),
      right_margin = option_or(options$right_margin, 0L),
      max_reads = option_or(options$max_reads, NULL),
      max_error = option_or(options$max_error, 0.3),
      whitelist = option_or(options$whitelist, NULL),
      threads = effective_threads,
      threads_requested = requested_threads,
      threads_effective = effective_threads,
      chunk_size = option_or(options$chunk_size, 10000L),
      selection = option_or(options$selection, "high_specificity"),
      rescan = option_or(options$rescan, FALSE),
      bc1_whitelist = option_or(options$bc1_whitelist, NULL),
      bc2_whitelist = option_or(options$bc2_whitelist, NULL),
      bc1_length = option_or(options$bc1_length, NULL),
      bc2_length = option_or(options$bc2_length, NULL),
      umi_length = option_or(options$umi_length, 9L),
      offset_min = option_or(options$offset_min, 0L),
      offset_max = option_or(options$offset_max, 3L),
      verbose = option_or(options$verbose, FALSE)
    ),
    files = list(
      whitelist = paste0(options$output_prefix, ".txt"),
      statistics = paste0(options$output_prefix, ".csv"),
      scan_log = paste0(options$output_prefix, "_scan_wl.log")
    ),
    stats = list(reads_processed = 10, selected_barcodes = 2),
    log = "mock scan"
  )
}

wrapper_reformat_result <- function(options, core = wrapper_core_info()) {
  split <- isTRUE(options$split_by_barcode)
  output <- if (split) {
    options$split_output_dir
  } else if (!is.null(options$output_path)) {
    options$output_path
  } else {
    options$input_path
  }
  list(
    success = TRUE,
    backend = "embedded",
    core = core,
    input = options$input_path,
    output = output,
    files = list(
      fastq = if (split) "" else output,
      split_fastq = if (split) {
        stats::setNames(file.path(output, "AACC.fq.gz"), "AACC")
      } else {
        character()
      }
    ),
    stats = list(records_processed = 2, records_written = 2),
    log = "mock reformat"
  )
}

test_that("rad_scan_wl passes a normalized typed configuration", {
  work <- tempfile("rrad-scan-wrapper-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "reads.fastq")
  whitelist <- file.path(work, "reference.txt")
  writeLines(c("@r", "ACGT", "+", "!!!!"), input)
  writeLines("AACC", whitelist)
  captured <- NULL
  core <- wrapper_core_info()

  local_mocked_bindings(
    rad_core_info = function() core,
    .rad_scan_wl_cpp = function(options) {
      captured <<- options
      wrapper_scan_result(options, core)
    },
    .package = "rrad"
  )

  result <- rad_scan_wl(
    input = input,
    output_prefix = file.path(work, "nested", "cells"),
    adapter_seq = "ACGTACGT",
    barcode_length = 4L,
    left_margin = 1L,
    right_margin = 2L,
    max_reads = 100L,
    max_error = 0.25,
    whitelist = whitelist,
    threads = 3L,
    chunk_size = 64L,
    selection = "ABOVE_FLOOR",
    verbose = TRUE
  )

  expect_s3_class(result, "rad_scan_wl_result")
  expect_identical(result$status, "success")
  expect_identical(result$backend, "embedded")
  expect_s3_class(result$core, "rad_core_info")
  expect_identical(captured$input, normalizePath(input))
  expect_identical(captured$whitelist, normalizePath(whitelist))
  expect_identical(captured$selection, "above_floor")
  expect_identical(captured$threads, 3L)
  expect_identical(captured$chunk_size, 64L)
  expect_identical(captured$max_error, 0.25)
  expect_true(grepl("nested/cells$", captured$output_prefix))
  expect_identical(result$config$threads_requested, 3L)
  expect_identical(result$config$threads_effective, 3L)
})

test_that("rad_scan_wl supports typed two-part whitelist options", {
  work <- tempfile("rrad-split-scan-wrapper-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "reads.fastq")
  bc1 <- file.path(work, "bc1.txt")
  bc2 <- file.path(work, "bc2.txt")
  writeLines(c("@r", "ACGT", "+", "!!!!"), input)
  writeLines("295", bc1)
  writeLines("456", bc2)
  captured <- NULL
  core <- wrapper_core_info()

  local_mocked_bindings(
    rad_core_info = function() core,
    .rad_scan_wl_cpp = function(options) {
      captured <<- options
      wrapper_scan_result(options, core)
    },
    .package = "rrad"
  )

  result <- rad_scan_wl(
    input,
    file.path(work, "split"),
    adapter_seq = "ACGT",
    bc1_whitelist = bc1,
    bc2_whitelist = bc2,
    bc1_length = 8L,
    bc2_length = 8L,
    umi_length = 10L,
    offset_min = 0L,
    offset_max = 2L
  )

  expect_s3_class(result, "rad_scan_wl_result")
  expect_identical(captured$bc1_length, 8L)
  expect_identical(captured$bc2_length, 8L)
  expect_identical(captured$offset_min, 0L)
  expect_identical(captured$offset_max, 2L)
})

test_that("rad_scan_wl derives a rescan prefix and omits nullable options", {
  work <- tempfile("rrad-rescan-wrapper-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "counts.csv.gz")
  writeLines("mock", input)
  captured <- NULL
  core <- wrapper_core_info()

  local_mocked_bindings(
    rad_core_info = function() core,
    .rad_scan_wl_cpp = function(options) {
      captured <<- options
      wrapper_scan_result(options, core)
    },
    .package = "rrad"
  )

  result <- rad_scan_wl(input = input, rescan = TRUE)

  expect_s3_class(result, "rad_scan_wl_result")
  expect_identical(captured$rescan, TRUE)
  expect_false("adapter_seq" %in% names(captured))
  expect_false("max_reads" %in% names(captured))
  expect_identical(captured$output_prefix,
                   file.path(normalizePath(work), "counts_rescan"))
})

test_that("rad_scan_wl batch mode validates every row before native work", {
  work <- tempfile("rrad-scan-batch-wrapper-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input_one <- file.path(work, "one.fastq")
  input_two <- file.path(work, "two.fastq")
  writeLines(c("@r1", "ACGT", "+", "!!!!"), input_one)
  writeLines(c("@r2", "TGCA", "+", "!!!!"), input_two)
  prefix_one <- file.path(work, "out", "one")
  prefix_two <- file.path(work, "out", "two")
  batch <- file.path(work, "batch.csv")
  utils::write.table(
    data.frame(input = c(input_one, input_two),
               output_prefix = c(prefix_one, prefix_two)),
    batch,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE
  )
  calls <- list()
  core <- wrapper_core_info()

  local_mocked_bindings(
    rad_core_info = function() core,
    .rad_scan_wl_cpp = function(options) {
      calls[[length(calls) + 1L]] <<- options
      wrapper_scan_result(options, core)
    },
    .package = "rrad"
  )

  result <- rad_scan_wl(batch_csv = batch, adapter_seq = "ACGT")
  expect_s3_class(result, "rad_scan_wl_batch_result")
  expect_length(result$runs, 2L)
  expect_true(all(vapply(
    result$runs, inherits, logical(1), what = "rad_scan_wl_result"
  )))
  expect_length(calls, 2L)
  expect_identical(calls[[2L]]$input, normalizePath(input_two))
  expect_identical(
    result$config,
    result$runs[[1L]]$config[
      setdiff(names(result$runs[[1L]]$config), c("input", "output_prefix"))
    ]
  )

  bad_batch <- file.path(work, "bad-batch.csv")
  utils::write.table(
    data.frame(input = c(input_one, file.path(work, "missing.fastq")),
               output_prefix = c(prefix_one, prefix_two)),
    bad_batch,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    quote = TRUE
  )
  calls <- list()
  expect_error(
    rad_scan_wl(batch_csv = bad_batch, adapter_seq = "ACGT"),
    "batch input row 2 does not identify an existing file"
  )
  expect_length(calls, 0L)

  duplicate_batch <- file.path(work, "duplicate-batch.csv")
  utils::write.table(
    data.frame(input = c(input_one, input_two),
               output_prefix = c(prefix_one, prefix_one)),
    duplicate_batch,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    quote = TRUE
  )
  expect_error(
    rad_scan_wl(batch_csv = duplicate_batch, adapter_seq = "ACGT"),
    "output prefixes must be unique"
  )
  expect_length(calls, 0L)

  case_collision_batch <- file.path(work, "case-collision-batch.csv")
  utils::write.table(
    data.frame(
      input = c(input_one, input_two),
      output_prefix = c(file.path(work, "Sample"),
                        file.path(work, "sample"))
    ),
    case_collision_batch,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    quote = TRUE
  )
  expect_error(
    rad_scan_wl(batch_csv = case_collision_batch, adapter_seq = "ACGT"),
    "colliding output artifacts"
  )
  expect_length(calls, 0L)

  self_overwrite_batch <- file.path(work, "self.csv")
  writeLines(paste(input_one, file.path(work, "self"), sep = ","),
             self_overwrite_batch)
  expect_error(
    rad_scan_wl(batch_csv = self_overwrite_batch, adapter_seq = "ACGT"),
    "would overwrite an input or batch_csv"
  )
  expect_length(calls, 0L)

  split_collision_batch <- file.path(work, "split-collision.csv")
  utils::write.table(
    data.frame(
      input = c(input_one, input_two),
      output_prefix = c(
        file.path(work, "sample"),
        file.path(work, "sample_valid_pairs")
      )
    ),
    split_collision_batch,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    quote = TRUE
  )
  expect_error(
    rad_scan_wl(
      batch_csv = split_collision_batch,
      adapter_seq = "ACGT",
      bc1_whitelist = "splitseq_bc1",
      bc2_whitelist = "splitseq_bc2"
    ),
    "colliding output artifacts"
  )
  expect_length(calls, 0L)
})

test_that("rad_scan_wl rejects invalid and contradictory configurations", {
  work <- tempfile("rrad-scan-validation-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "reads.fastq")
  batch <- file.path(work, "batch.csv")
  writeLines(c("@r", "ACGT", "+", "!!!!"), input)
  writeLines(paste(input, file.path(work, "out"), sep = ","), batch)

  expect_error(rad_scan_wl(), "exactly one of input and batch_csv")
  expect_error(
    rad_scan_wl(input, batch_csv = batch, adapter_seq = "ACGT"),
    "exactly one of input and batch_csv"
  )
  expect_error(rad_scan_wl(input), "adapter_seq is required")
  unsupported <- file.path(work, "reads.txt")
  writeLines("not a supported sequence suffix", unsupported)
  expect_error(
    rad_scan_wl(unsupported, adapter_seq = "ACGT"),
    "must have a .fastq, .fq, .fasta, or .fa extension"
  )
  expect_error(rad_scan_wl(input, adapter_seq = "ACGT", barcode_length = 0),
               "barcode_length")
  expect_error(
    rad_scan_wl(input, adapter_seq = "ACGT", barcode_length = 33L),
    "barcode_length must be between 1 and 32"
  )
  expect_error(rad_scan_wl(input, adapter_seq = "ACGT", left_margin = -1),
               "left_margin")
  expect_error(rad_scan_wl(input, adapter_seq = "ACGT", max_reads = 0),
               "max_reads")
  expect_error(rad_scan_wl(input, adapter_seq = "ACGT", max_error = -0.1),
               "max_error")
  expect_error(rad_scan_wl(input, adapter_seq = "ACGT", max_error = 1.01),
               "max_error.*between 0 and 1")
  expect_error(rad_scan_wl(input, adapter_seq = "ACGT", whitelist = ":"),
               "one whitelist file or built-in kit key")
  expect_error(rad_scan_wl(input, adapter_seq = "ACGT", threads = 0),
               "threads")
  expect_error(rad_scan_wl(input, adapter_seq = "ACGT", chunk_size = 0),
               "chunk_size")
  expect_error(rad_scan_wl(input, adapter_seq = "ACGT", selection = "all"),
               "one of")
  expect_error(
    rad_scan_wl(
      input, file.path(work, "nested", ".."), adapter_seq = "ACGT"
    ),
    "filename prefix"
  )
  expect_error(
    rad_scan_wl(input, adapter_seq = "ACGT", bc1_whitelist = "splitseq_bc1"),
    "must be supplied together"
  )
  expect_error(
    rad_scan_wl(
      input, adapter_seq = "ACGT", whitelist = "10x_3v3",
      bc1_whitelist = "splitseq_bc1", bc2_whitelist = "splitseq_bc2"
    ),
    "are alternatives"
  )
  expect_error(
    rad_scan_wl(
      input, rescan = TRUE,
      bc1_whitelist = "splitseq_bc1", bc2_whitelist = "splitseq_bc2"
    ),
    "rescan does not support"
  )
  expect_error(
    rad_scan_wl(input, rescan = TRUE, whitelist = "10x_3v3"),
    "rescan does not support a reference whitelist"
  )
  expect_error(
    rad_scan_wl(input, adapter_seq = "ACGT", bc1_length = 8L),
    "require two-part"
  )
  expect_error(
    rad_scan_wl(
      input, adapter_seq = "ACGT",
      bc1_whitelist = "splitseq_bc1", bc2_whitelist = "splitseq_bc2",
      bc1_length = 33L
    ),
    "bc1_length must be between 1 and 32"
  )
  expect_error(
    rad_scan_wl(input, adapter_seq = "ACGT", offset_min = 2L,
                offset_max = 1L),
    "offset_min must not exceed"
  )
  expect_error(
    rad_scan_wl(input, adapter_seq = "ACGT", offset_min = -1L),
    "offset_min"
  )
  expect_error(
    rad_scan_wl(batch_csv = batch, adapter_seq = "ACGT", rescan = TRUE),
    "rescan does not support batch_csv"
  )
  expect_error(
    rad_scan_wl(batch_csv = batch, output_prefix = "unused",
                adapter_seq = "ACGT"),
    "output_prefix must be NULL"
  )
})

test_that("rad_scan_wl enforces feature and result handshakes", {
  work <- tempfile("rrad-scan-handshake-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "reads.fastq")
  writeLines(c("@r", "ACGT", "+", "!!!!"), input)

  local_mocked_bindings(
    rad_core_info = function() wrapper_core_info(features = "demux"),
    .rad_scan_wl_cpp = function(options) stop("must not run"),
    .package = "rrad"
  )
  expect_error(
    rad_scan_wl(input, adapter_seq = "ACGT"),
    "does not advertise the scan-wl feature"
  )

  core <- wrapper_core_info()
  local_mocked_bindings(
    rad_core_info = function() core,
    .rad_scan_wl_cpp = function(options) 1L,
    .package = "rrad"
  )
  expect_error(
    rad_scan_wl(input, adapter_seq = "ACGT"),
    "invalid whitelist-scan result"
  )

  changed <- wrapper_core_info(
    source_digest = paste(rep("d", 64L), collapse = "")
  )
  local_mocked_bindings(
    rad_core_info = function() core,
    .rad_scan_wl_cpp = function(options) wrapper_scan_result(options, changed),
    .package = "rrad"
  )
  expect_error(
    rad_scan_wl(input, adapter_seq = "ACGT"),
    "identity changed during the whitelist scan"
  )
})

test_that("scan and reformat record effective serial thread counts", {
  work <- tempfile("rrad-wrapper-serial-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "reads.fastq")
  writeLines(c("@r", "ACGT", "+", "!!!!"), input)
  core <- wrapper_core_info(openmp = FALSE)
  scan_options <- NULL
  reformat_options <- NULL

  local_mocked_bindings(
    rad_core_info = function() core,
    .rad_scan_wl_cpp = function(options) {
      scan_options <<- options
      wrapper_scan_result(options, core)
    },
    .rad_reformat_cpp = function(options) {
      reformat_options <<- options
      wrapper_reformat_result(options, core)
    },
    .package = "rrad"
  )

  expect_warning(
    scan <- rad_scan_wl(input, adapter_seq = "ACGT", threads = 4L),
    "no OpenMP support"
  )
  expect_identical(scan_options$threads, 4L)
  expect_identical(scan$config$threads_effective, 1L)

  output <- file.path(work, "rewritten.fastq")
  expect_warning(
    reformatted <- rad_reformat(
      input, reformat_header = TRUE, output_fastq = output, threads = 4L
    ),
    "no OpenMP support"
  )
  expect_identical(reformat_options$threads, 4L)
  expect_identical(reformatted$config$threads_effective, 1L)
})

test_that("rad_reformat delegates split and rewrite options", {
  work <- tempfile("rrad-reformat-wrapper-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "demux.fq.gz")
  writeLines("mock", input)
  captured <- NULL
  core <- wrapper_core_info()

  local_mocked_bindings(
    rad_core_info = function() core,
    .rad_reformat_cpp = function(options) {
      captured <<- options
      wrapper_reformat_result(options, core)
    },
    .package = "rrad"
  )

  result <- rad_reformat(
    input,
    out_dir = file.path(work, "by-cell"),
    split_bc = TRUE,
    reformat_header = TRUE,
    collapsed_input = TRUE,
    delimiter = "-",
    threads = 3L,
    chunk_size = 128L,
    verbose = TRUE
  )

  expect_s3_class(result, "rad_reformat_result")
  expect_identical(result$status, "success")
  expect_s3_class(result$core, "rad_core_info")
  expect_identical(captured$input_path, normalizePath(input))
  expect_identical(captured$split_by_barcode, TRUE)
  expect_identical(captured$reformat_header, TRUE)
  expect_identical(captured$parse_collapsed_id, TRUE)
  expect_identical(captured$delimiter, "-")
  expect_identical(captured$threads, 3L)
  expect_identical(captured$chunk_size, 128L)
  expect_false("output_path" %in% names(captured))
  expect_identical(result$output, normalizePath(file.path(work, "by-cell")))
})

test_that("rad_reformat delegates coordinate and aggregate output paths", {
  work <- tempfile("rrad-coordinate-wrapper-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "demux.fq.gz")
  output <- file.path(work, "nested", "spatial.fq.gz")
  writeLines("mock", input)
  captured <- NULL
  core <- wrapper_core_info()

  local_mocked_bindings(
    rad_core_info = function() core,
    .rad_reformat_cpp = function(options) {
      captured <<- options
      wrapper_reformat_result(options, core)
    },
    .package = "rrad"
  )

  result <- rad_reformat(
    input,
    coordinate = "vizHD-v1",
    bin_size = 8L,
    output_fastq = output
  )

  expect_s3_class(result, "rad_reformat_result")
  expect_identical(captured$coordinate_mode, "vizHD-v1")
  expect_identical(captured$bin_size_um, 8L)
  expect_identical(
    captured$output_path,
    file.path(normalizePath(dirname(output), mustWork = TRUE), basename(output))
  )
  expect_identical(result$output, captured$output_path)
})

test_that("rad_reformat rejects invalid and contradictory configurations", {
  work <- tempfile("rrad-reformat-validation-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "demux.fastq")
  writeLines(c("@r", "ACGT", "+", "!!!!"), input)

  expect_error(rad_reformat(input), "request at least one")
  expect_error(rad_reformat(input, split_bc = TRUE), "out_dir is required")
  expect_error(
    rad_reformat(
      input, out_dir = file.path(work, "split"), split_bc = TRUE,
      output_fastq = file.path(work, "aggregate.fastq")
    ),
    "output_fastq cannot be supplied"
  )
  expect_error(
    rad_reformat(input, out_dir = work, reformat_header = TRUE),
    "out_dir is only used"
  )
  expect_error(rad_reformat(input, reformat_header = TRUE, delimiter = "::"),
               "delimiter")
  expect_error(rad_reformat(input, reformat_header = TRUE, delimiter = "\n"),
               "delimiter cannot")
  expect_error(rad_reformat(input, reformat_header = TRUE, delimiter = " "),
               "delimiter cannot")
  expect_error(rad_reformat(input, coordinate = "vizHD-v1", bin_size = 0),
               "bin_size")
  expect_error(rad_reformat(input, reformat_header = TRUE, threads = 0),
               "threads")
  expect_error(rad_reformat(input, reformat_header = TRUE, chunk_size = 0),
               "chunk_size")
  expect_error(rad_reformat(input, split_bc = 1, out_dir = work),
               "split_bc must be TRUE or FALSE")
})

test_that("rad_reformat enforces feature and result handshakes", {
  work <- tempfile("rrad-reformat-handshake-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "demux.fastq")
  output <- file.path(work, "out.fastq")
  writeLines(c("@r", "ACGT", "+", "!!!!"), input)

  local_mocked_bindings(
    rad_core_info = function() wrapper_core_info(features = "demux"),
    .rad_reformat_cpp = function(options) stop("must not run"),
    .package = "rrad"
  )
  expect_error(
    rad_reformat(input, reformat_header = TRUE, output_fastq = output),
    "does not advertise the reformat feature"
  )

  core <- wrapper_core_info()
  local_mocked_bindings(
    rad_core_info = function() core,
    .rad_reformat_cpp = function(options) list(success = TRUE),
    .package = "rrad"
  )
  expect_error(
    rad_reformat(input, reformat_header = TRUE, output_fastq = output),
    "invalid reformat result"
  )

  changed <- wrapper_core_info(
    source_digest = paste(rep("d", 64L), collapse = "")
  )
  local_mocked_bindings(
    rad_core_info = function() core,
    .rad_reformat_cpp = function(options) {
      wrapper_reformat_result(options, changed)
    },
    .package = "rrad"
  )
  expect_error(
    rad_reformat(input, reformat_header = TRUE, output_fastq = output),
    "identity changed during the reformat run"
  )
})
