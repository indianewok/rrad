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
  local_mocked_bindings(
    .rad_scan_wl_cpp = function(options) stop("unexpected native execution"),
    .package = "rrad"
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
  expect_error(
    rad_scan_wl(batch_csv = bad_batch, adapter_seq = "ACGT"),
    "batch input row 2 does not identify an existing file"
  )

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

  self_overwrite_batch <- file.path(work, "self.csv")
  writeLines(paste(input_one, file.path(work, "self"), sep = ","),
             self_overwrite_batch)
  expect_error(
    rad_scan_wl(batch_csv = self_overwrite_batch, adapter_seq = "ACGT"),
    "would overwrite an input or batch_csv"
  )

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
    rad_scan_wl(input, adapter_seq = "ACGT", barcode_length = 65L),
    "barcode_length must be between 1 and 64"
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
      bc1_length = 65L
    ),
    "bc1_length must be between 1 and 64"
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
  expect_error(rad_reformat(input, to_fasta = 1),
               "to_fasta must be TRUE or FALSE")
  expect_error(rad_reformat(input, to_fasta = TRUE),
               "output_fastq is required")
  expect_error(
    rad_reformat(
      input,
      output_fastq = file.path(work, "wrong.fastq"),
      to_fasta = TRUE
    ),
    "must end in .fa"
  )
  malformed_demux <- structure(
    list(artifacts = list()),
    class = c("rad_demux_result", "list")
  )
  expect_error(
    rad_reformat(
      malformed_demux,
      output_fastq = file.path(work, "out.fa"),
      to_fasta = TRUE
    ),
    "does not contain one aggregate FASTQ"
  )
  expect_error(
    rad_reformat(input, header_format = "unknown"),
    "header_format"
  )
  expect_error(
    rad_reformat(
      input,
      reformat_header = TRUE,
      header_format = "presto"
    ),
    "conflicts with header_format"
  )
  expect_error(
    rad_reformat(
      input,
      header_format = "presto",
      collapsed_input = TRUE
    ),
    "requires actual SAM tags"
  )
})

test_that("scan and reformat reject unsupported core features", {
  input <- test_path("fixtures", "rad-core", "tiny.fastq")
  output <- tempfile("rad-feature-", fileext = ".fa")
  core <- rad_core_info()
  core$features <- "demux"
  local_mocked_bindings(
    rad_core_info = function() core,
    .rad_scan_wl_cpp = function(options) stop("unexpected native execution"),
    .rad_reformat_cpp = function(options) stop("unexpected native execution"),
    .package = "rrad"
  )
  expect_error(rad_scan_wl(input, adapter_seq = "ACGT"),
               "does not advertise the scan-wl feature")
  expect_error(rad_reformat(input, reformat_header = TRUE),
               "does not advertise the reformat feature")
  core$features <- c("demux", "reformat")
  expect_error(rad_reformat(input, to_fasta = TRUE, output_fastq = output),
               "does not advertise FASTA reformatting")
  expect_error(rad_reformat(input, header_format = "presto"),
               "does not advertise pRESTO header formatting")
})
