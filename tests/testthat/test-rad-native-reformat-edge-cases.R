reformat_edge_read_lines <- function(path) {
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(connection), add = TRUE)
  readLines(connection, warn = FALSE)
}

reformat_edge_first_barcode <- function(path) {
  lines <- reformat_edge_read_lines(path)
  if (tolower(trimws(lines[[1L]])) == "whitelist_bcs") {
    lines[[2L]]
  } else {
    sub(",.*", "", lines[[1L]])
  }
}

test_that("native in-place reformat preserves empty FASTQ and FASTA records", {
  work <- tempfile("rad-native-empty-reformat-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  fastq <- file.path(work, "empty.fastq")
  writeLines(c(
    "@empty-fastq CB:Z:AACC UB:Z:TGCA",
    "",
    "+empty-fastq retained-plus-comment",
    ""
  ), fastq)

  fastq_result <- rad_reformat(
    fastq, reformat_header = TRUE, threads = 1L, chunk_size = 1L
  )

  expect_identical(
    reformat_edge_read_lines(fastq),
    c("@empty-fastq_AACC_TGCA", "",
      "+empty-fastq retained-plus-comment", "")
  )
  expect_equal(fastq_result$stats$records_read, 1)
  expect_equal(fastq_result$stats$records_written, 1)
  expect_equal(fastq_result$stats$empty_sequences_skipped, 0)

  fasta <- file.path(work, "empty.fasta")
  writeLines(c(">empty-fasta CB:Z:AACC", ""), fasta)

  fasta_result <- rad_reformat(
    fasta, reformat_header = TRUE, threads = 1L, chunk_size = 1L
  )

  expect_identical(
    reformat_edge_read_lines(fasta),
    c(">empty-fasta_AACC", "")
  )
  expect_equal(fasta_result$stats$records_read, 1)
  expect_equal(fasta_result$stats$records_written, 1)
  expect_equal(fasta_result$stats$empty_sequences_skipped, 0)
})

test_that("native split-only reformat preserves the FASTQ plus line", {
  work <- tempfile("rad-native-plus-reformat-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "tagged.fastq")
  writeLines(c(
    "@read-one CB:Z:AACC UB:Z:TGCA XX:Z:keep",
    "ACGT",
    "+read-one retained plus-line comment",
    "IIII"
  ), input)

  result <- rad_reformat(
    input,
    out_dir = file.path(work, "split"),
    split_bc = TRUE,
    threads = 1L,
    chunk_size = 1L
  )
  output <- reformat_edge_read_lines(
    result$artifacts$split_fastq[["AACC"]]
  )

  expect_identical(output[[3L]], "+read-one retained plus-line comment")
  expect_identical(output[[4L]], "IIII")
})

test_that("reformat recognizes only complete SAM-style barcode tag fields", {
  work <- tempfile("rad-native-exact-tag-")
  dir.create(work)
  input <- file.path(work, "reads.fastq")
  split_dir <- file.path(work, "split")
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  writeLines(c(
    "@false-positive XX:Z:prefixCB:Z:wrong UB:Z:AACC",
    "ACGT", "+", "IIII",
    "@real XX:Z:prefixCB:Z:wrong CB:Z:RIGHT UB:Z:AACC",
    "TGCA", "+", "IIII"
  ), input)

  result <- rad_reformat(
    input, out_dir = split_dir, split_bc = TRUE,
    threads = 1L, chunk_size = 1L
  )

  expect_identical(unname(names(result$artifacts$split_fastq)), "RIGHT")
  expect_equal(result$stats$records_read, 2)
  expect_equal(result$stats$records_written, 1)
  expect_equal(result$stats$records_skipped_missing_cb, 1)
})

test_that("split reformat rejects barcode filenames that collide by case", {
  work <- tempfile("rad-native-case-collision-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  probe <- file.path(work, "case-probe")
  writeLines("probe", probe)
  case_insensitive <- file.exists(file.path(work, "CASE-PROBE"))
  unlink(probe)
  skip_if_not(case_insensitive, "filesystem distinguishes filename case")

  input <- file.path(work, "reads.fastq")
  split_dir <- file.path(work, "split")
  writeLines(c(
    "@lower CB:Z:abc UB:Z:AACC", "A", "+", "!",
    "@upper CB:Z:ABC UB:Z:AACC", "C", "+", "5"
  ), input)

  expect_error(
    rad_reformat(input, out_dir = split_dir, split_bc = TRUE,
                 threads = 1L, chunk_size = 2L),
    "same split-output filename"
  )
  expect_length(list.files(split_dir, pattern = "\\.fq\\.gz$"), 0L)
})

test_that("split publication rolls back allocation failures after installs", {
  work <- tempfile("rad-native-split-allocation-rollback-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "reads.fastq")
  split_dir <- file.path(work, "split")
  writeLines(c(
    "@first CB:Z:AACC UB:Z:AAAA", "ACGT", "+", "IIII",
    "@second CB:Z:CCGG UB:Z:CCCC", "TGCA", "+", "IIII"
  ), input)

  failpoint <- "RRAD_TEST_REFORMAT_FAIL_AFTER_SPLIT_INSTALL"
  previous <- Sys.getenv(failpoint, unset = NA_character_)
  on.exit({
    if (is.na(previous)) {
      Sys.unsetenv(failpoint)
    } else {
      do.call(Sys.setenv, setNames(list(previous), failpoint))
    }
  }, add = TRUE)
  for (fail_after in c("1", "2")) {
    do.call(Sys.setenv, setNames(list(fail_after), failpoint))
    expect_error(
      rad_reformat(
        input, out_dir = split_dir, split_bc = TRUE,
        threads = 1L, chunk_size = 2L
      )
    )
    expect_length(list.files(split_dir, pattern = "\\.fq\\.gz$"), 0L)
    expect_length(list.files(
      split_dir, pattern = "^\\.rrad-reformat-stage-", all.files = TRUE
    ), 0L)
  }

  Sys.unsetenv(failpoint)
  result <- rad_reformat(
    input, out_dir = split_dir, split_bc = TRUE,
    threads = 1L, chunk_size = 2L
  )
  expect_setequal(names(result$artifacts$split_fastq), c("AACC", "CCGG"))
})

test_that("R result construction completes before reformat publication", {
  work <- tempfile("rad-native-result-allocation-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "reads.fastq")
  writeLines(c(
    "@first CB:Z:AACC UB:Z:AAAA", "ACGT", "+", "IIII",
    "@second CB:Z:CCGG UB:Z:CCCC", "TGCA", "+", "IIII"
  ), input)
  original <- readBin(input, "raw", n = file.info(input)$size)

  failpoint <- "RRAD_TEST_REFORMAT_FAIL_RESULT_CONSTRUCTION"
  previous <- Sys.getenv(failpoint, unset = NA_character_)
  on.exit({
    if (is.na(previous)) {
      Sys.unsetenv(failpoint)
    } else {
      do.call(Sys.setenv, setNames(list(previous), failpoint))
    }
  }, add = TRUE)
  do.call(Sys.setenv, setNames(list("1"), failpoint))

  split_dir <- file.path(work, "split")
  expect_error(rad_reformat(
    input, out_dir = split_dir, split_bc = TRUE,
    threads = 1L, chunk_size = 2L
  ))
  expect_length(list.files(split_dir, pattern = "\\.fq\\.gz$"), 0L)
  expect_length(list.files(
    split_dir, pattern = "^\\.rrad-reformat-stage-", all.files = TRUE
  ), 0L)

  expect_error(rad_reformat(
    input, reformat_header = TRUE, threads = 1L, chunk_size = 2L
  ))
  expect_identical(
    readBin(input, "raw", n = file.info(input)$size),
    original
  )
  expect_length(list.files(
    work, pattern = "^\\.rrad-reformat-stage-", all.files = TRUE
  ), 0L)
})

test_that("header collapse and coordinate mapping are mutually exclusive", {
  work <- tempfile("rad-native-coordinate-conflict-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "tagged.fastq")
  writeLines(c(
    "@read-one CB:Z:AACC UB:Z:TGCA", "ACGT", "+", "IIII"
  ), input)

  expect_error(
    rad_reformat(
      input,
      reformat_header = TRUE,
      coordinate = "vizHD-v1",
      output_fastq = file.path(work, "public.fastq")
    ),
    "reformat_header cannot be combined with coordinate"
  )

  expect_error(
    rrad:::.rad_reformat_cpp(list(
      input_path = normalizePath(input),
      reformat_header = TRUE,
      coordinate_mode = "vizHD-v1",
      threads = 1L,
      chunk_size = 1L
    )),
    "reformat_header cannot be combined with coordinate_mode"
  )
})

test_that("coordinate reformat preserves unrelated tags and FASTQ plus data", {
  resources <- system.file("rad/resources", package = "rrad", mustWork = TRUE)
  x_barcode <- reformat_edge_first_barcode(file.path(
    resources, "wl", "visium_hd_bc1.csv.gz"
  ))
  y_barcode <- reformat_edge_first_barcode(file.path(
    resources, "wl", "visium_hd_bc2.csv.gz"
  ))

  work <- tempfile("rad-native-coordinate-metadata-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "spatial.fastq")
  output <- file.path(work, "spatial-out.fastq")
  writeLines(c(
    paste0(
      "@spot BC:Z:old CB:Z:", x_barcode, "-", y_barcode,
      " UB:Z:AACC RG:Z:old XX:Z:keep FI:i:7"
    ),
    "ACGT",
    "+spot retained-coordinate-plus-comment",
    "IIII"
  ), input)

  result <- rad_reformat(
    input,
    coordinate = "vizHD-v1",
    bin_size = 2L,
    output_fastq = output,
    threads = 1L,
    chunk_size = 1L
  )
  lines <- reformat_edge_read_lines(output)
  fields <- strsplit(lines[[1L]], "\t", fixed = TRUE)[[1L]]

  expect_equal(result$stats$coordinate_mapped, 1)
  expect_identical(fields[[1L]], "@spot_00000_00000")
  expect_true("BC:Z:s_002um_00000_00000-1" %in% fields)
  expect_true(paste0("CB:Z:", x_barcode, "-", y_barcode) %in% fields)
  expect_true("UB:Z:AACC" %in% fields)
  expect_true("RG:Z:s002um" %in% fields)
  expect_true(all(c("XX:Z:keep", "FI:i:7") %in% fields))
  expect_false(any(c("BC:Z:old", "RG:Z:old") %in% fields))
  expect_identical(lines[[3L]], "+spot retained-coordinate-plus-comment")
})
