native_read_lines <- function(path) {
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(connection), add = TRUE)
  readLines(connection, warn = FALSE)
}

native_write_fastq <- function(path, ids, sequences, comments = NULL) {
  if (is.null(comments)) comments <- rep("", length(ids))
  stopifnot(length(ids) == length(sequences),
            length(ids) == length(comments))
  headers <- paste0("@", ids, ifelse(nzchar(comments), " ", ""), comments)
  qualities <- vapply(
    sequences,
    function(sequence) paste(rep("I", nchar(sequence)), collapse = ""),
    character(1)
  )
  writeLines(as.vector(rbind(headers, sequences, "+", qualities)), path)
}

test_that("native whitelist scan writes selected calls and statistics", {
  input <- test_path("fixtures", "rad-core", "tiny.fastq")
  work <- tempfile("rad-native-scan-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  result <- rad_scan_wl(
    input = input,
    output_prefix = file.path(work, "cells"),
    adapter_seq = "ACGTACGT",
    barcode_length = 4L,
    max_error = 0,
    threads = 1L,
    chunk_size = 1L,
    selection = "above_floor"
  )

  expect_s3_class(result, "rad_scan_wl_result")
  expect_identical(result$status, "success")
  expect_true(all(c("scan-wl", "reformat") %in% result$core$features))
  expect_true(all(file.exists(unlist(result$artifacts[c(
    "csv", "whitelist", "scan_log"
  )]))))
  expect_setequal(native_read_lines(result$artifacts$whitelist),
                  c("AACC", "TTGA"))
  expect_equal(result$stats$reads_processed, 2)
  expect_equal(result$stats$total_extractions, 2)
  expect_equal(result$stats$selected_barcodes, 2)
  expect_identical(result$stats$selection, "above_floor")
})

test_that("native whitelist scan accepts FASTA input", {
  work <- tempfile("rad-native-fasta-scan-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "reads.fasta")
  writeLines(c(">read-one", "ACGTAACC", ">read-two", "ACGTTTGA"), input)

  result <- rad_scan_wl(
    input, file.path(work, "cells"), adapter_seq = "ACGT",
    barcode_length = 4L, max_error = 0, selection = "above_floor"
  )

  expect_equal(result$stats$reads_processed, 2)
  expect_setequal(native_read_lines(result$artifacts$whitelist),
                  c("AACC", "TTGA"))
})

test_that("native scan verbose output reports committed artifact paths", {
  input <- test_path("fixtures", "rad-core", "tiny.fastq")
  work <- tempfile("rad-native-scan-verbose-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  prefix <- file.path(work, "cells")

  console <- capture.output(result <- rad_scan_wl(
    input, prefix, adapter_seq = "ACGTACGT", barcode_length = 4L,
    max_error = 0, selection = "above_floor", verbose = TRUE
  ))

  expected <- paste0("Results written to: ", result$artifacts$statistics)
  expect_match(result$log, expected,
               fixed = TRUE)
  expect_match(paste(console, collapse = "\n"),
               expected, fixed = TRUE)
  expect_false(grepl(".rrad-tmp-", result$log, fixed = TRUE))
})

test_that("native whitelist rescan accepts a barcode count table", {
  work <- tempfile("rad-native-rescan-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  counts <- file.path(work, "counts.csv")
  writeLines(c(
    "barcode,count", "AAAA,100", "CCCC,80", "GGGG,2", "TTTT,1"
  ), counts)

  result <- rad_scan_wl(
    counts,
    file.path(work, "recalled"),
    rescan = TRUE,
    selection = "high_specificity"
  )

  expect_identical(result$stats$input_rows, 4)
  expect_equal(result$stats$total_extractions, 183)
  expect_setequal(native_read_lines(result$artifacts$whitelist),
                  c("AAAA", "CCCC"))
})

test_that("native two-part whitelist scan emits pairs and spatial mask", {
  work <- tempfile("rad-native-split-scan-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  bc1 <- file.path(work, "bc1.txt")
  bc2 <- file.path(work, "bc2.txt")
  input <- file.path(work, "reads.fastq")
  writeLines(c("barcode", "AACC", "CCGG"), bc1)
  writeLines(c("barcode", "TTGA", "GGTT"), bc2)
  native_write_fastq(
    input,
    c("r1", "r1-duplicate", "r2", "r3"),
    c("ACGTGGAACCTTGA", "ACGTGGAACCTTGA",
      "ACGTTTCCGGGGTT", "ACGTCCAAAATTTT")
  )

  result <- rad_scan_wl(
    input,
    file.path(work, "split"),
    adapter_seq = "ACGT",
    bc1_whitelist = bc1,
    bc2_whitelist = bc2,
    umi_length = 2L,
    offset_min = 0L,
    offset_max = 0L,
    max_error = 0,
    selection = "above_floor"
  )

  expect_equal(result$stats$reads_processed, 4)
  expect_equal(result$stats$total_extractions, 3)
  expect_true(file.exists(result$artifacts$valid_pairs))
  expect_true(file.exists(result$artifacts$spatial_mask))
  pairs <- utils::read.csv(result$artifacts$valid_pairs)
  expect_setequal(pairs$pair, c("AACCTTGA", "CCGGGGTT"))
  expect_equal(pairs$count[match(c("AACCTTGA", "CCGGGGTT"), pairs$pair)],
               c(2, 1))
  expect_identical(native_read_lines(result$artifacts$spatial_mask),
                   c("1,0", "0,1"))

  expect_error(
    rad_scan_wl(
      input, file.path(work, "mismatched-length"), adapter_seq = "ACGT",
      bc1_whitelist = bc1, bc2_whitelist = bc2,
      bc1_length = 8L, bc2_length = 4L, umi_length = 2L,
      offset_min = 0L, offset_max = 0L, max_error = 0
    ),
    "bc1_length.*does not match"
  )
})

test_that("native scan rejects malformed FASTQ without publishing outputs", {
  work <- tempfile("rad-native-bad-scan-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "bad.fastq")
  prefix <- file.path(work, "cells")
  writeLines(c("@good", "ACGTAACC", "+", "IIIIIIII",
               "@truncated", "ACGT", "+"), input)

  expect_error(
    rad_scan_wl(input, prefix, adapter_seq = "ACGT", barcode_length = 4L),
    "sequence read failed|could not read"
  )
  expect_false(any(file.exists(c(
    paste0(prefix, ".csv"), paste0(prefix, ".txt"),
    paste0(prefix, "_scan_wl.log")
  ))))
})

test_that("native scan accepts uppercase FASTQ and gzip suffixes", {
  source <- test_path("fixtures", "rad-core", "tiny.fastq")
  work <- tempfile("rad-native-uppercase-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  plain <- file.path(work, "READS.FASTQ")
  compressed <- file.path(work, "READS.FASTQ.GZ")
  expect_true(file.copy(source, plain))
  output <- gzfile(compressed, open = "wt")
  writeLines(readLines(source, warn = FALSE), output)
  close(output)

  for (input in c(plain, compressed)) {
    prefix <- file.path(work, paste0("cells-", basename(input)))
    result <- rad_scan_wl(
      input, prefix, adapter_seq = "ACGTACGT", barcode_length = 4L,
      max_error = 0, selection = "above_floor"
    )
    expect_equal(result$stats$reads_processed, 2)
    expect_setequal(native_read_lines(result$artifacts$whitelist),
                    c("AACC", "TTGA"))
  }
})

test_that("exact scan adapters honor N as a canonical-base wildcard", {
  work <- tempfile("rad-native-wildcard-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "reads.fastq")
  native_write_fastq(input, "read", "AACGAACC")

  result <- rad_scan_wl(
    input, file.path(work, "cells"), adapter_seq = "AACN",
    barcode_length = 4L, max_error = 0, selection = "above_floor"
  )
  expect_identical(native_read_lines(result$artifacts$whitelist), "AACC")
})

test_that("corrupt gzip rescans fail without publishing artifacts", {
  work <- tempfile("rad-native-corrupt-rescan-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "counts.csv.gz")
  prefix <- file.path(work, "rescanned")
  output <- gzfile(input, open = "wt")
  writeLines(c("barcode,count", "AAAA,100", "CCCC,50"), output)
  close(output)
  payload <- readBin(input, what = "raw", n = file.info(input)$size)
  writeBin(payload[seq_len(length(payload) - 5L)], input)

  expect_error(
    rad_scan_wl(input, prefix, rescan = TRUE),
    "gzip|compressed|rescan"
  )
  expect_false(any(file.exists(c(
    paste0(prefix, ".csv"), paste0(prefix, ".txt"),
    paste0(prefix, "_scan_wl.log")
  ))))
})

test_that("split scanning bounds very large search offsets to each read", {
  work <- tempfile("rad-native-offset-bound-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "reads.fastq")
  bc1 <- file.path(work, "bc1.txt")
  bc2 <- file.path(work, "bc2.txt")
  writeLines("AACC", bc1)
  writeLines("TTGA", bc2)
  native_write_fastq(input, "read", "ACGTGGAAAA")

  result <- rad_scan_wl(
    input, file.path(work, "split"), adapter_seq = "ACGT",
    bc1_whitelist = bc1, bc2_whitelist = bc2, umi_length = 2L,
    offset_min = 0L, offset_max = .Machine$integer.max,
    max_error = 0, selection = "above_floor"
  )
  expect_equal(result$stats$reads_processed, 1)
  expect_equal(result$stats$total_extractions, 0)
})

test_that("bundled SPLiT-seq packed bitlists infer their eight-base length", {
  decode_packed <- function(value, width = 8L) {
    alphabet <- c("A", "C", "T", "G")
    codes <- integer(width)
    for (index in rev(seq_len(width))) {
      codes[[index]] <- value %% 4
      value <- floor(value / 4)
    }
    paste0(alphabet[codes + 1L], collapse = "")
  }
  resources <- system.file("rad/resources", package = "rrad",
                           mustWork = TRUE)
  first_value <- function(name) {
    connection <- gzfile(file.path(resources, "wl", name), open = "rt")
    on.exit(close(connection), add = TRUE)
    as.numeric(readLines(connection, n = 1L, warn = FALSE)[[1L]])
  }
  bc1 <- decode_packed(first_value("splitseq_bc1_bitlist.csv.gz"))
  bc2 <- decode_packed(first_value("splitseq_bc2_bitlist.csv.gz"))
  work <- tempfile("rad-native-splitseq-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "reads.fastq")
  native_write_fastq(input, "read", paste0("ACGTGG", bc1, bc2))

  result <- rad_scan_wl(
    input, file.path(work, "splitseq"), adapter_seq = "ACGT",
    bc1_whitelist = "splitseq_bc1", bc2_whitelist = "splitseq_bc2",
    umi_length = 2L, offset_min = 0L, offset_max = 0L,
    max_error = 0, selection = "above_floor"
  )
  pairs <- utils::read.csv(result$artifacts$valid_pairs)
  expect_identical(pairs$bc1, bc1)
  expect_identical(pairs$bc2, bc2)
  expect_equal(pairs$count, 1)
  expect_identical(result$config$bc1_length, 8L)
  expect_identical(result$config$bc2_length, 8L)
})

test_that("scan output preflight preserves earlier artifacts on obstruction", {
  input <- test_path("fixtures", "rad-core", "tiny.fastq")
  work <- tempfile("rad-native-scan-transaction-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  prefix <- file.path(work, "cells")
  csv <- paste0(prefix, ".csv")
  whitelist <- paste0(prefix, ".txt")
  writeLines("preexisting statistics", csv)
  dir.create(whitelist)

  expect_error(
    rad_scan_wl(
      input, prefix, adapter_seq = "ACGTACGT", barcode_length = 4L,
      max_error = 0
    ),
    "directory"
  )
  expect_identical(readLines(csv), "preexisting statistics")
  expect_true(dir.exists(whitelist))
  expect_false(file.exists(paste0(prefix, "_scan_wl.log")))
})

test_that("scan publication waits on one interprocess transaction lock", {
  python <- Sys.which("python3")
  skip_if(!nzchar(python), "requires python3 for an external fcntl lock")

  input <- test_path("fixtures", "rad-core", "tiny.fastq")
  work <- tempfile("rad-native-scan-lock-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  prefix <- file.path(work, "cells")
  targets <- c(
    paste0(prefix, ".csv"),
    paste0(prefix, ".txt"),
    paste0(prefix, "_scan_wl.log"),
    paste0(prefix, "_valid_pairs.csv"),
    paste0(prefix, "_spat_mask.csv")
  )
  sentinels <- c(
    "old csv", "old whitelist", "old log", "old pairs", "old mask"
  )
  Map(writeLines, sentinels, targets)

  ready <- file.path(work, "holder-ready")
  observed <- file.path(work, "holder-observed")
  stage_observed <- file.path(work, "holder-stage-observed")
  holder_log <- file.path(work, "holder.log")
  lock_path <- file.path(work, ".cells.rrad-commit.lock")
  holder <- file.path(work, "hold_lock.py")
  writeLines(c(
    "import fcntl, pathlib, stat, sys, time",
    "lock_path, ready, observed, stage_observed = sys.argv[1:5]",
    "targets = [pathlib.Path(value) for value in sys.argv[5:]]",
    paste0(
      "sentinels = [b'old csv\\n', b'old whitelist\\n', b'old log\\n', ",
      "b'old pairs\\n', b'old mask\\n']"
    ),
    "with open(lock_path, 'a+b') as lock:",
    "    fcntl.lockf(lock, fcntl.LOCK_EX)",
    "    pathlib.Path(ready).write_text('ready\\n')",
    "    deadline = time.monotonic() + 3.0",
    "    intact = True",
    "    stage_report = None",
    "    while time.monotonic() < deadline:",
    "        intact = all(path.exists() and path.read_bytes() == expected",
    "                     for path, expected in zip(targets, sentinels))",
    "        if not intact:",
    "            break",
    "        parent = targets[0].parent",
    paste0(
      "        stages = [path for path in parent.iterdir() ",
      "if path.name.startswith('.rrad-scan-stage-') and path.is_dir()]"
    ),
    paste0(
      "        legacy = [path for path in parent.iterdir() ",
      "if '.rrad-tmp-' in path.name]"
    ),
    "        for stage in stages:",
    "            children = sorted(path.name for path in stage.iterdir()",
    "                              if path.is_file())",
    "            if len(children) >= 3:",
    "                stage_report = '\\n'.join([",
    "                    'stage_count=' + str(len(stages)),",
    "                    'name=' + stage.name,",
    "                    'mode=' + format(stat.S_IMODE(stage.stat().st_mode), 'o'),",
    "                    'artifacts=' + ','.join(children),",
    "                    'legacy_siblings=' + str(len(legacy))])",
    "        time.sleep(0.01)",
    "    pathlib.Path(observed).write_text('intact\\n' if intact else 'changed\\n')",
    "    if stage_report is None:",
    paste0(
      "        legacy = [path for path in targets[0].parent.iterdir() ",
      "if '.rrad-tmp-' in path.name]"
    ),
    "        stage_report = '\\n'.join([",
    "            'stage_count=0',",
    "            'legacy_siblings=' + str(len(legacy))])",
    "    pathlib.Path(stage_observed).write_text(stage_report + '\\n')"
  ), holder)

  quote_arg <- function(value) shQuote(value)
  launch <- suppressWarnings(system2(
    python,
    vapply(c(holder, lock_path, ready, observed, stage_observed, targets),
           quote_arg, character(1)),
    stdout = holder_log,
    stderr = holder_log,
    wait = FALSE
  ))
  expect_identical(launch, 0L)
  deadline <- Sys.time() + 10
  while (!file.exists(ready) && Sys.time() < deadline) Sys.sleep(0.02)
  expect_true(
    file.exists(ready),
    info = paste(readLines(holder_log, warn = FALSE), collapse = "\n")
  )

  started <- proc.time()[["elapsed"]]
  result <- rad_scan_wl(
    input, prefix, adapter_seq = "ACGTACGT", barcode_length = 4L,
    max_error = 0, selection = "above_floor"
  )
  elapsed <- proc.time()[["elapsed"]] - started

  expect_true(file.exists(observed))
  expect_identical(readLines(observed, warn = FALSE), "intact")
  expect_true(file.exists(stage_observed))
  stage_report <- readLines(stage_observed, warn = FALSE)
  expect_true("stage_count=1" %in% stage_report)
  expect_true("mode=700" %in% stage_report)
  expect_true("legacy_siblings=0" %in% stage_report)
  expect_true(paste0(
    "artifacts=cells.csv,cells.txt,cells_scan_wl.log"
  ) %in% stage_report)
  stage_name <- sub("^name=", "", grep("^name=", stage_report, value = TRUE))
  expect_match(stage_name, "^\\.rrad-scan-stage-cells-[^/]{6}$")
  expect_true(elapsed >= 2.5, info = paste("elapsed:", elapsed))
  expect_identical(as.character(file.info(lock_path)$mode), "600")
  expect_setequal(native_read_lines(result$artifacts$whitelist),
                  c("AACC", "TTGA"))
  expect_false(any(vapply(
    Map(readLines, targets[seq_len(3L)]),
    function(lines) any(grepl("^old ", lines)),
    logical(1)
  )))
  expect_false(any(file.exists(targets[4:5])))
  expect_length(
    list.files(work, pattern = "^\\.rrad-scan-stage-", all.files = TRUE),
    0L
  )
})

test_that("native reformat collapses headers and splits safely", {
  input <- test_path("fixtures", "rad-core", "tagged.fastq")
  work <- tempfile("rad-native-reformat-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  output <- file.path(work, "collapsed.fastq.gz")

  aggregate <- rad_reformat(
    input,
    reformat_header = TRUE,
    output_fastq = output,
    threads = 1L,
    chunk_size = 1L
  )
  aggregate_lines <- native_read_lines(output)
  expect_s3_class(aggregate, "rad_reformat_result")
  expect_identical(aggregate_lines[c(1L, 5L, 9L)],
                   c("@read-one_AACC_TGCA", "@read-two_TTGA_CCCC",
                     "@read-three"))
  expect_equal(aggregate$stats$records_read, 3)
  expect_equal(aggregate$stats$records_written, 3)

  split_dir <- file.path(work, "by-cell")
  split <- rad_reformat(
    input,
    out_dir = split_dir,
    split_bc = TRUE,
    reformat_header = TRUE,
    threads = 1L,
    chunk_size = 1L
  )
  expect_identical(names(split$artifacts$split_fastq), c("AACC", "TTGA"))
  expect_true(all(file.exists(split$artifacts$split_fastq)))
  expect_equal(split$stats$records_written, 2)
  expect_equal(split$stats$records_skipped_missing_cb, 1)
  expect_identical(
    native_read_lines(split$artifacts$split_fastq[["AACC"]])[[1L]],
    "@read-one_AACC_TGCA"
  )
})

test_that("native split reformat bounds writers and preserves repeated cells", {
  work <- tempfile("rad-native-writer-lru-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "many.fastq")
  barcodes <- c(sprintf("BC%03d", 1:66), "BC001")
  native_write_fastq(
    input,
    paste0("r", seq_along(barcodes)),
    rep("ACGT", length(barcodes)),
    paste0("CB:Z:", barcodes, " UB:Z:AACC")
  )

  result <- rad_reformat(
    input,
    out_dir = file.path(work, "split"),
    split_bc = TRUE,
    threads = 1L,
    chunk_size = 67L
  )

  expect_equal(result$stats$split_barcodes, 66)
  expect_gt(result$stats$gzip_members_written,
            result$stats$split_barcodes)
  expect_length(native_read_lines(
    result$artifacts$split_fastq[["BC001"]]
  ), 8L)
})

test_that("native coordinate reformat uses packaged Visium HD axes", {
  resources <- system.file("rad/resources", package = "rrad",
                           mustWork = TRUE)
  first_barcode <- function(path) {
    lines <- native_read_lines(path)
    if (tolower(trimws(lines[[1L]])) == "whitelist_bcs") {
      lines[[2L]]
    } else {
      sub(",.*", "", lines[[1L]])
    }
  }
  x_barcode <- first_barcode(file.path(
    resources, "wl", "visium_hd_bc1.csv.gz"
  ))
  y_barcode <- first_barcode(file.path(
    resources, "wl", "visium_hd_bc2.csv.gz"
  ))
  work <- tempfile("rad-native-coordinate-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "spatial.fastq")
  output <- file.path(work, "spatial-out.fastq.gz")
  native_write_fastq(
    input, c("spot", "unmapped_AACC_TGCA"), c("ACGT", "TGCA"),
    c(paste0("CB:Z:", x_barcode, "-", y_barcode, " UB:Z:AACC"), "")
  )

  result <- rad_reformat(
    input,
    coordinate = "vizHD-v1",
    collapsed_input = TRUE,
    bin_size = 2L,
    output_fastq = output,
    threads = 1L
  )
  headers <- native_read_lines(output)[c(1L, 5L)]

  expect_equal(result$stats$coordinate_mapped, 1)
  expect_equal(result$stats$coordinate_unmapped, 1)
  expect_match(headers[[1L]], "^@spot_00000_00000\\t")
  expect_match(headers[[1L]], "BC:Z:s_002um_00000_00000-1", fixed = TRUE)
  expect_match(headers[[1L]], "RG:Z:s002um", fixed = TRUE)
  expect_identical(headers[[2L]], "@unmapped_AACC_TGCA")
})

test_that("native reformat failure preserves in-place input and staging", {
  work <- tempfile("rad-native-reformat-rollback-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "bad.fastq")
  writeLines(c(
    "@good CB:Z:AACC UB:Z:TGCA", "ACGT", "+", "IIII",
    "@truncated", "ACGT", "+"
  ), input)
  before <- readBin(input, what = "raw", n = file.info(input)$size)

  expect_error(
    rad_reformat(input, reformat_header = TRUE, threads = 1L,
                 chunk_size = 1L),
    "sequence read failed"
  )
  after <- readBin(input, what = "raw", n = file.info(input)$size)
  expect_identical(after, before)
  expect_length(list.files(
    work, pattern = "^\\.rrad-reformat-stage-", all.files = TRUE
  ), 0L)
})

test_that("native reformat commits a successful in-place rewrite", {
  source <- test_path("fixtures", "rad-core", "tagged.fastq")
  work <- tempfile("rad-native-in-place-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "tagged.fastq")
  expect_true(file.copy(source, input))

  result <- rad_reformat(input, reformat_header = TRUE)
  headers <- native_read_lines(input)[c(1L, 5L, 9L)]
  expect_identical(result$output, normalizePath(input))
  expect_identical(headers,
                   c("@read-one_AACC_TGCA", "@read-two_TTGA_CCCC",
                     "@read-three"))
})

test_that("collapsed IDs are opt-in and explicit tags remain authoritative", {
  work <- tempfile("rad-native-collapsed-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  ordinary <- file.path(work, "ordinary.fastq")
  native_write_fastq(ordinary, "sample_with_underscores", "ACGT")
  ordinary_result <- rad_reformat(
    ordinary, out_dir = file.path(work, "ordinary-split"), split_bc = TRUE
  )
  expect_equal(ordinary_result$stats$records_skipped_missing_cb, 1)
  expect_length(ordinary_result$artifacts$split_fastq, 0L)

  collapsed <- file.path(work, "collapsed.fastq")
  collapsed_id <- "sample_read_AACC_TGCA"
  native_write_fastq(collapsed, collapsed_id, "ACGT")
  collapsed_result <- rad_reformat(
    collapsed, out_dir = file.path(work, "collapsed-split"),
    split_bc = TRUE, collapsed_input = TRUE
  )
  expect_identical(names(collapsed_result$artifacts$split_fastq), "AACC")
  expect_identical(
    native_read_lines(collapsed_result$artifacts$split_fastq[["AACC"]])[[1L]],
    paste0("@", collapsed_id)
  )

  tagged <- file.path(work, "tagged.fastq")
  native_write_fastq(
    tagged, collapsed_id, "ACGT", "CB:Z:CCCC UB:Z:GGGG"
  )
  tagged_result <- rad_reformat(
    tagged, out_dir = file.path(work, "tagged-split"),
    split_bc = TRUE, collapsed_input = TRUE
  )
  expect_identical(names(tagged_result$artifacts$split_fastq), "CCCC")
  expect_identical(
    native_read_lines(tagged_result$artifacts$split_fastq[["CCCC"]])[[1L]],
    paste0("@", collapsed_id, "\tCB:Z:CCCC UB:Z:GGGG")
  )
})

test_that("native FASTA splitting preserves records and case-folded suffixes", {
  work <- tempfile("rad-native-fasta-split-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "READS.FASTA.GZ")
  connection <- gzfile(input, open = "wt")
  writeLines(c(">fa-read CB:Z:AACC UB:Z:TGCA", "ACGT"), connection)
  close(connection)

  result <- rad_reformat(
    input, out_dir = file.path(work, "split"), split_bc = TRUE
  )
  path <- result$artifacts$split_fastq[["AACC"]]
  expect_match(path, "\\.fa\\.gz$")
  expect_identical(native_read_lines(path),
                   c(">fa-read\tCB:Z:AACC UB:Z:TGCA", "ACGT"))
})

test_that("native split reformat rejects unsafe barcode filenames", {
  work <- tempfile("rad-native-reformat-unsafe-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "unsafe.fastq")
  split_dir <- file.path(work, "split")
  native_write_fastq(
    input, "read", "ACGT", "CB:Z:../escape UB:Z:AACC"
  )

  expect_error(
    rad_reformat(input, out_dir = split_dir, split_bc = TRUE,
                 threads = 1L),
    "unsafe for a split-output filename"
  )
  expect_false(file.exists(file.path(work, "escape.fq.gz")))
  expect_length(list.files(split_dir, all.files = TRUE,
                           no.. = TRUE), 0L)
})
