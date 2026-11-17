test_that("rad_core_info exposes a strict embedded-core handshake", {
  info <- rad_core_info()

  expect_s3_class(info, "rad_core_info")
  expect_identical(info$core_api_version, 1L)
  expect_true(info$vendored)
  expect_true("demux" %in% info$features)
  expect_true("seqspec" %in% info$features)
  expect_true("scan-wl" %in% info$features)
  expect_true("reformat" %in% info$features)
  expect_type(info$build, "list")
  expect_true(nzchar(info$core_version))
  expect_true(nzchar(info$source_commit))
  expect_match(info$embedded_source_digest, "^[[:xdigit:]]{64}$")
})

test_that("rad_demux runs a complete one-read embedded demultiplex", {
  fastq <- test_path("fixtures", "rad-core", "nanopore-smoke.fastq")
  out_dir <- tempfile("rad-core-native-")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  dir.create(out_dir)

  # A matching output prefix must not implicitly select stale native state.
  writeLines("poisoned stale layout", file.path(out_dir, "smoke_layout.csv"))
  writeLines("poisoned stale map",
             file.path(out_dir, "smoke_position_map.csv"))

  result <- rad_demux(
    layout = "nanopore_rapid_bc",
    fastq = fastq,
    out_dir = out_dir,
    output = "smoke",
    threads = 1L,
    chunk_size = 1L,
    max_reads = 1L,
    min_read_length = 0L
  )

  expect_s3_class(result, "rad_demux_result")
  expect_identical(result$status, "success")
  expect_identical(result$backend, "embedded")
  expect_equal(result$stats$total_reads, 1)
  expect_equal(result$stats$reads_demultiplexed, 1)
  expect_equal(result$stats$records_serialized, 1)
  expect_true(file.exists(result$artifacts$fastq))
  expect_true(file.exists(result$artifacts$layout))
  expect_true(file.exists(result$artifacts$position_map))
  expect_true(file.exists(result$artifacts$demux_log))
  expect_identical(length(readLines(gzfile(result$artifacts$fastq), n = 4L)),
                   4L)
  expect_false(any(grepl("poisoned stale", readLines(result$artifacts$layout),
                         fixed = TRUE)))
})

test_that("an all-filtered run returns a valid empty primary gzip", {
  fastq <- test_path("fixtures", "rad-core", "nanopore-smoke.fastq")
  out_dir <- tempfile("rad-core-all-filtered-")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  result <- rad_demux(
    layout = "nanopore_rapid_bc",
    fastq = fastq,
    out_dir = out_dir,
    output = "filtered",
    threads = 1L,
    chunk_size = 1L,
    min_read_length = 1000000L
  )

  expect_equal(result$stats$records_serialized, 0)
  expect_true(file.exists(result$artifacts$fastq))
  expect_length(readLines(gzfile(result$artifacts$fastq)), 0L)
})

test_that("rad_demux translates seqspec YAML on the cache-miss path", {
  fastq <- test_path("fixtures", "rad-core", "nanopore-smoke.fastq")
  seqspec <- test_path("fixtures", "rad-core", "nanopore-smoke.yaml")
  out_dir <- tempfile("rad-core-seqspec-")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  temp_root <- dirname(tempdir())
  temp_pattern <- "^rad_seqspec_[[:alnum:]]{6}$"
  before <- list.files(temp_root, pattern = temp_pattern, full.names = TRUE)

  result <- rad_demux(
    layout = seqspec,
    fastq = fastq,
    out_dir = out_dir,
    output = "yaml",
    threads = 1L,
    chunk_size = 1L,
    max_reads = 1L,
    min_read_length = 0L
  )

  after <- list.files(temp_root, pattern = temp_pattern, full.names = TRUE)
  layout_lines <- readLines(result$artifacts$layout)
  expect_s3_class(result, "rad_demux_result")
  expect_equal(result$stats$total_reads, 1)
  expect_equal(result$stats$reads_demultiplexed, 1)
  expect_true(any(grepl("linker_1", layout_lines, fixed = TRUE)))
  expect_true(any(grepl("linker_3", layout_lines, fixed = TRUE)))
  expect_setequal(after, before)
})

test_that("seqspec uses explicit local onlists without filename substitution", {
  work <- tempfile("rad-core-seqspec-local-onlist-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  onlist <- file.path(work, "my-3m-february-2018-custom.txt")
  spec <- file.path(work, "local-onlist.yaml")
  writeLines("AACC", onlist)
  writeLines(c(
    "seqspec_version: 0.3.0",
    "library_spec:",
    "  - region_id: RNA",
    "    region_type: RNA",
    "    regions:",
    "      - region_id: anchor",
    "        region_type: adapter",
    "        sequence_type: fixed",
    "        sequence: ACGTACGT",
    "      - region_id: barcode",
    "        region_type: barcode",
    "        sequence_type: onlist",
    "        min_len: 4",
    "        max_len: 4",
    "        onlist:",
    paste0("          filename: ", basename(onlist)),
    "          location: local",
    "          urltype: file",
    "          url: ./",
    "      - region_id: poly_t",
    "        region_type: poly_t",
    "        sequence_type: fixed",
    "        sequence: TTTTTTTT",
    "      - region_id: read",
    "        region_type: cdna",
    "        sequence_type: random"
  ), spec)

  old <- options(rrad.offline = TRUE)
  on.exit(options(old), add = TRUE)
  result <- rad_demux(
    spec, test_path("fixtures", "rad-core", "tiny.fastq"),
    file.path(work, "out"), output = "local", threads = 1L,
    chunk_size = 1L, max_reads = 1L, min_read_length = 0L
  )

  expect_identical(result$status, "success")
  expect_equal(result$stats$loaded_true_barcodes, 1)
  expect_true(any(grepl(
    normalizePath(onlist), readLines(result$artifacts$layout), fixed = TRUE
  )))
})

test_that("embedded seqspec rejects remote onlists without network access", {
  work <- tempfile("rad-core-seqspec-remote-onlist-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  spec <- file.path(work, "remote-onlist.yaml")
  writeLines(c(
    "seqspec_version: 0.3.0",
    "library_spec:",
    "  - region_id: RNA",
    "    region_type: RNA",
    "    regions:",
    "      - region_id: anchor",
    "        region_type: adapter",
    "        sequence_type: fixed",
    "        sequence: ACGTACGT",
    "      - region_id: barcode",
    "        region_type: barcode",
    "        sequence_type: onlist",
    "        min_len: 4",
    "        max_len: 4",
    "        onlist:",
    "          filename: remote-custom.txt",
    "          location: remote",
    "          urltype: https",
    "          url: https://example.invalid/remote-custom.txt",
    "      - region_id: read",
    "        region_type: cdna",
    "        sequence_type: random"
  ), spec)
  fastq <- test_path("fixtures", "rad-core", "tiny.fastq")

  old_env <- Sys.getenv("RRAD_OFFLINE", unset = NA_character_)
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("RRAD_OFFLINE") else {
      Sys.setenv(RRAD_OFFLINE = old_env)
    }
  }, add = TRUE)
  Sys.setenv(RRAD_OFFLINE = "true")
  expect_error(
    rad_demux(spec, fastq, file.path(work, "offline"), threads = 1L),
    "RRAD_OFFLINE is enabled"
  )

  Sys.unsetenv("RRAD_OFFLINE")
  expect_error(
    rad_demux(spec, fastq, file.path(work, "embedded"), threads = 1L),
    "remote onlists are disabled in embedded rrad"
  )

  # An explicit R-side whitelist is authoritative and must prevent seqspec
  # metadata from initiating (or requiring) a network request.
  local_whitelist <- test_path(
    "fixtures", "rad-core", "tiny-whitelist.txt"
  )
  Sys.setenv(RRAD_OFFLINE = "true")
  custom <- rad_demux(
    spec, fastq, file.path(work, "custom-override"),
    custom_whitelist = local_whitelist, threads = 1L,
    chunk_size = 1L, max_reads = 1L, min_read_length = 0L
  )
  expect_identical(custom$status, "success")
  expect_equal(custom$stats$loaded_true_barcodes, 2)

  global <- rad_demux(
    spec, fastq, file.path(work, "global-override"),
    global_whitelist = local_whitelist, threads = 1L,
    chunk_size = 1L, max_reads = 1L, min_read_length = 0L
  )
  expect_identical(global$status, "success")
  expect_equal(global$stats$loaded_global_barcodes, 2)
})

test_that("padded barcode correction works through 64 bases and rejects 65", {
  work <- tempfile("rad-core-64-base-barcode-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  longest_barcode <- paste0(
    "GACCTAGTTCAGGATCCGTAACGTTCGAGTGC",
    "TAGCATCGAGTCCTAGACTGACGATCTGATCG"
  )
  whitelist <- file.path(work, "barcodes.txt")
  layout <- file.path(work, "layout.csv")
  fastq <- file.path(work, "reads.fastq")
  for (barcode_length in c(24L, 25L, 32L, 33L, 57L, 63L, 64L)) {
    barcode <- substr(longest_barcode, 1L, barcode_length)
    writeLines(barcode, whitelist)
    writeLines(c(
      "Read Layout,,,,,",
      "id,seq,expected_length,type,class,whitelist",
      paste0(
        "anchor_1,AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAA,,",
        "static,linker,"
      ),
      paste0("barcode,,", barcode_length,
             ",variable,barcode,", whitelist),
      "anchor_2,CAGCACCTGATTACAAGGTGCTG,,static,linker,",
      "read,,,variable,read,"
    ), layout)
    sequence <- paste0(
      "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAA", barcode,
      "CAGCACCTGATTACAAGGTGCTG", "GATTACAGATTACA"
    )
    writeLines(c(
      "@read", sequence, "+",
      paste(rep("I", nchar(sequence)), collapse = "")
    ), fastq)

    result <- rad_demux(
      layout, fastq, file.path(work, paste0("out-", barcode_length)),
      output = paste0("packed", barcode_length), threads = 1L,
      chunk_size = 1L, min_read_length = 0L,
      whitelist_mutation = 0L, generated_mutation = 0L
    )
    expect_identical(result$status, "success")
    expect_equal(result$stats$loaded_true_barcodes, 1)
    expect_equal(result$stats$reads_demultiplexed, 1)
    expect_equal(result$stats$records_serialized, 1)
    output_header <- readLines(gzfile(result$artifacts$fastq), n = 1L)
    expect_match(output_header, paste0("CB:Z:", barcode), fixed = TRUE)
  }

  writeLines(paste0(longest_barcode, "G"), whitelist)
  layout_lines <- readLines(layout)
  layout_lines[4] <- paste0("barcode,,65,variable,barcode,", whitelist)
  writeLines(layout_lines, layout)
  expect_error(
    rad_demux(
      layout, fastq, file.path(work, "too-long"), output = "packed65",
      threads = 1L, whitelist_mutation = 0L, generated_mutation = 0L
    ),
    "64 bases"
  )
})

test_that("malformed compound whitelist specs fail safely", {
  work <- tempfile("rad-core-malformed-whitelist-spec-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  layout <- file.path(work, "layout.csv")
  lines <- readLines(test_path("fixtures", "rad-core", "tiny-layout.csv"))
  lines[4] <- "barcode,,4,variable,barcode,:"
  writeLines(lines, layout)

  expect_error(
    rad_demux(
      layout, test_path("fixtures", "rad-core", "tiny.fastq"),
      file.path(work, "out"), threads = 1L
    ),
    "Malformed whitelist specification"
  )
})

test_that("spatial masks require exact binary dimensions", {
  work <- tempfile("rad-core-spatial-mask-validation-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  axis1 <- file.path(work, "axis1.txt")
  axis2 <- file.path(work, "axis2.txt")
  mask <- file.path(work, "mask.csv")
  layout <- file.path(work, "layout.csv")
  writeLines(c("AACC", "CCGG"), axis1)
  writeLines(c("TTGA", "GGTT"), axis2)
  writeLines(c(
    "Read Layout,,,,,,",
    "id,seq,expected_length,type,class,whitelist,flags",
    "anchor,ACGTACGT,,static,linker,,",
    paste0("barcode_1,,4,variable,barcode,", axis1, ":", mask,
           ",joint_barcode"),
    paste0("barcode_2,,4,variable,barcode,", axis2, ":", mask,
           ",joint_barcode"),
    "read,,,variable,read,,"
  ), layout)
  fastq <- test_path("fixtures", "rad-core", "tiny.fastq")
  run <- function(name) {
    rad_demux(
      layout, fastq, file.path(work, name), output = name, threads = 1L,
      whitelist_mutation = 0L, generated_mutation = 0L
    )
  }

  writeLines("1,0", mask)
  expect_error(run("short-rows"), "row count does not match")

  writeLines(c("1", "0"), mask)
  expect_error(run("short-columns"), "column count does not match")

  writeLines(c("1,2", "0,1"), mask)
  expect_error(run("non-binary"), "cells must be exactly 0 or 1")
})

test_that("spatial axis gzip corruption fails closed", {
  work <- tempfile("rad-core-spatial-axis-gzip-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  axis1 <- file.path(work, "axis1.txt.gz")
  axis2 <- file.path(work, "axis2.txt")
  mask <- file.path(work, "mask.csv")
  layout <- file.path(work, "layout.csv")
  con <- gzfile(axis1, "wb")
  writeLines(c("AACC", "CCGG"), con)
  close(con)
  bytes <- readBin(axis1, "raw", n = file.info(axis1)$size)
  bytes[length(bytes)] <- as.raw(bitwXor(as.integer(bytes[length(bytes)]),
                                         255L))
  writeBin(bytes, axis1)
  writeLines(c("TTGA", "GGTT"), axis2)
  writeLines(c("1,0", "0,1"), mask)
  writeLines(c(
    "Read Layout,,,,,,",
    "id,seq,expected_length,type,class,whitelist,flags",
    "anchor,ACGTACGT,,static,linker,,",
    paste0("barcode_1,,4,variable,barcode,", axis1, ":", mask,
           ",joint_barcode"),
    paste0("barcode_2,,4,variable,barcode,", axis2, ":", mask,
           ",joint_barcode"),
    "read,,,variable,read,,"
  ), layout)

  expect_error(
    rad_demux(
      layout, test_path("fixtures", "rad-core", "tiny.fastq"),
      file.path(work, "out"), threads = 1L,
      whitelist_mutation = 0L, generated_mutation = 0L
    ),
    "gzip"
  )
})

test_that("auto-whitelist preserves a layout-derived global catalog", {
  fixture_dir <- test_path("fixtures", "rad-core")
  fastq <- file.path(fixture_dir, "tiny.fastq")
  layout <- tempfile("rad-core-auto-layout-", fileext = ".csv")
  reference <- tempfile("rad-core-auto-reference-", fileext = ".txt")
  out_dir <- tempfile("rad-core-auto-reference-")
  on.exit(
    unlink(c(layout, reference, out_dir), recursive = TRUE),
    add = TRUE
  )

  writeLines(c("AACC", "TTGA", "CCCC"), reference)
  layout_lines <- readLines(file.path(fixture_dir, "tiny-layout.csv"))
  layout_lines[4] <- paste0(
    "barcode,,4,variable,barcode,", normalizePath(reference)
  )
  writeLines(layout_lines, layout)

  result <- rad_demux(
    layout, fastq, out_dir, output = "auto-reference",
    auto_whitelist = TRUE, threads = 1L, chunk_size = 1L,
    max_reads = 2L, min_read_length = 0L
  )

  expect_equal(result$stats$loaded_true_barcodes, 2)
  expect_equal(result$stats$loaded_global_barcodes, 1)
  expect_s3_class(result$scan, "rad_scan_wl_result")
  expect_identical(result$scan$config$adapter_seq, "ACGTACGT")
  expect_identical(result$scan$config$barcode_length, 4L)
  expect_identical(result$scan$config$whitelist, normalizePath(reference))
  expect_true(all(file.exists(unlist(result$scan$artifacts[c(
    "csv", "whitelist", "scan_log"
  )]))))
  expect_identical(result$config$scan_effective, result$scan$config)
})

test_that("layout and position caches are reused only as a complete pair", {
  fastq <- test_path("fixtures", "rad-core", "nanopore-smoke.fastq")
  out_dir <- tempfile("rad-core-cache-pair-")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  dir.create(out_dir)

  writeLines("poisoned orphan position map", file.path(
    out_dir, "orphan-position_position_map.csv"
  ))
  position_result <- rad_demux(
    layout = "nanopore_rapid_bc", fastq = fastq, out_dir = out_dir,
    output = "orphan-position", threads = 1L, chunk_size = 1L,
    max_reads = 1L, min_read_length = 0L, reuse_cache = TRUE
  )
  expect_identical(position_result$status, "success")
  expect_false(any(grepl(
    "poisoned orphan",
    readLines(position_result$artifacts$position_map, warn = FALSE),
    fixed = TRUE
  )))

  writeLines("poisoned orphan layout", file.path(
    out_dir, "orphan-layout_layout.csv"
  ))
  layout_result <- rad_demux(
    layout = "nanopore_rapid_bc", fastq = fastq, out_dir = out_dir,
    output = "orphan-layout", threads = 1L, chunk_size = 1L,
    max_reads = 1L, min_read_length = 0L, reuse_cache = TRUE
  )
  expect_identical(layout_result$status, "success")
  expect_false(any(grepl(
    "poisoned orphan",
    readLines(layout_result$artifacts$layout, warn = FALSE),
    fixed = TRUE
  )))
})

test_that("debug output is additive to the primary demultiplexed FASTQ", {
  fastq <- test_path("fixtures", "rad-core", "nanopore-smoke.fastq")
  plain_dir <- tempfile("rad-core-plain-")
  debug_dir <- tempfile("rad-core-debug-")
  on.exit(unlink(c(plain_dir, debug_dir), recursive = TRUE), add = TRUE)

  args <- list(
    layout = "nanopore_rapid_bc",
    fastq = fastq,
    output = "compare",
    threads = 1L,
    chunk_size = 1L,
    max_reads = 1L,
    min_read_length = 0L
  )
  plain <- do.call(rad_demux, c(args, list(out_dir = plain_dir)))
  debug <- do.call(rad_demux, c(
    args, list(out_dir = debug_dir, write_debug = TRUE)
  ))

  expect_true(file.exists(debug$artifacts$fastq))
  expect_identical(
    readLines(gzfile(debug$artifacts$fastq)),
    readLines(gzfile(plain$artifacts$fastq))
  )
  expect_true(all(file.exists(unlist(debug$artifacts[c(
    "debug_sig", "debug_csv", "debug_fastq", "metrics"
  )]))))
  metrics <- read.delim(debug$artifacts$metrics, check.names = FALSE)
  expect_named(metrics, c(
    "chunk_id", "seqs_in_chunk", "seqs_passed", "in_flight",
    "process_time_ms", "queue_time_ms", "total_time_ms", "rss_mb"
  ))
  expect_equal(ncol(metrics), 8L)
  expect_true(all(metrics$in_flight >= 1))
  expect_equal(
    metrics$total_time_ms,
    metrics$process_time_ms + metrics$queue_time_ms
  )
})

test_that("native demux rejects non-bundled resource trees", {
  fastq <- test_path("fixtures", "rad-core", "nanopore-smoke.fastq")
  foreign <- tempfile("rad-core-foreign-resources-")
  on.exit(unlink(foreign, recursive = TRUE), add = TRUE)
  dir.create(file.path(foreign, "read_layout"), recursive = TRUE)
  dir.create(file.path(foreign, "wl"), recursive = TRUE)

  expect_error(
    rad_demux_cpp(list(
      layout = "nanopore_rapid_bc",
      fastq = fastq,
      resource_dir = foreign
    )),
    "canonical bundled RAD resource directory"
  )
  expect_error(
    rad_demux_cpp(list(
      layout = "nanopore_rapid_bc", fastq = fastq,
      kit = "10x_3v3", global_whitelist = fastq
    )),
    "kit.*global_whitelist.*alternatives"
  )
  expect_error(
    rad_demux_cpp(list(
      layout = "nanopore_rapid_bc", fastq = fastq,
      auto_whitelist = TRUE, custom_whitelist = fastq
    )),
    "auto_whitelist.*custom_whitelist.*alternatives"
  )
})

test_that("native demux rejects malformed and unsupported sequence input", {
  good <- readLines(test_path(
    "fixtures", "rad-core", "nanopore-smoke.fastq"
  ))
  seqspec <- test_path("fixtures", "rad-core", "nanopore-smoke.yaml")
  bad_fastq <- tempfile("rad-core-truncated-", fileext = ".fastq")
  unsupported <- tempfile("rad-core-unsupported-", fileext = ".txt")
  out_dir <- tempfile("rad-core-invalid-input-")
  on.exit(unlink(c(bad_fastq, unsupported, out_dir), recursive = TRUE),
          add = TRUE)
  writeLines(c(good, "@truncated", "ACGT", "+"), bad_fastq)
  writeLines(good, unsupported)

  expect_error(
    rad_demux(
      layout = "nanopore_rapid_bc",
      fastq = bad_fastq,
      out_dir = out_dir,
      threads = 1L,
      chunk_size = 2L,
      min_read_length = 0L
    ),
    "sequence read failed \\(kseq code -2\\)"
  )
  temp_pattern <- "^rad_seqspec_[[:alnum:]]{6}$"
  before_error <- list.files(dirname(tempdir()), pattern = temp_pattern,
                             full.names = TRUE)
  expect_error(
    rad_demux(
      layout = seqspec,
      fastq = unsupported,
      out_dir = out_dir,
      output = "unsupported",
      threads = 1L,
      min_read_length = 0L
    ),
    "unsupported sequence input"
  )
  after_error <- list.files(dirname(tempdir()), pattern = temp_pattern,
                            full.names = TRUE)
  expect_setequal(after_error, before_error)
})

test_that("compressed input and whitelist corruption fail closed", {
  good <- readLines(test_path(
    "fixtures", "rad-core", "nanopore-smoke.fastq"
  ))
  layout <- test_path("fixtures", "rad-core", "tiny-layout.csv")
  tiny_fastq <- test_path("fixtures", "rad-core", "tiny.fastq")
  out_dir <- tempfile("rad-core-corrupt-gzip-")
  truncated <- tempfile("rad-core-body-", fileext = ".fastq.gz")
  bad_crc <- tempfile("rad-core-crc-", fileext = ".fastq.gz")
  bad_whitelist <- tempfile("rad-core-wl-", fileext = ".txt.gz")
  on.exit(
    unlink(c(out_dir, truncated, bad_crc, bad_whitelist), recursive = TRUE),
    add = TRUE
  )

  write_gzip <- function(path, lines) {
    connection <- gzfile(path, open = "wb")
    on.exit(close(connection), add = TRUE)
    writeLines(lines, connection)
  }
  write_gzip(truncated, good)
  bytes <- readBin(truncated, what = "raw", n = file.info(truncated)$size)
  writeBin(bytes[seq_len(length(bytes) - 8L)], truncated)

  write_gzip(bad_crc, good)
  bytes <- readBin(bad_crc, what = "raw", n = file.info(bad_crc)$size)
  bytes[length(bytes) - 7L] <- as.raw(bitwXor(
    as.integer(bytes[length(bytes) - 7L]), 255L
  ))
  writeBin(bytes, bad_crc)

  for (input in c(truncated, bad_crc)) {
    expect_error(
      rad_demux(
        layout = "nanopore_rapid_bc", fastq = input,
        out_dir = out_dir, output = basename(input), threads = 1L,
        chunk_size = 1L, min_read_length = 0L
      ),
      "sequence read failed.*gzip|gzip.*sequence read failed"
    )
  }

  write_gzip(bad_whitelist, c("AACC", "TTGA"))
  bytes <- readBin(
    bad_whitelist, what = "raw", n = file.info(bad_whitelist)$size
  )
  writeBin(bytes[seq_len(length(bytes) - 8L)], bad_whitelist)
  expect_error(
    rad_demux(
      layout = layout, fastq = tiny_fastq, out_dir = out_dir,
      output = "bad-whitelist", custom_whitelist = bad_whitelist,
      threads = 1L, chunk_size = 1L, min_read_length = 0L
    ),
    "whitelist.*gzip|gzip.*whitelist|gzip text read failed"
  )
  expect_error(
    rad_demux(
      layout = layout, fastq = tiny_fastq, out_dir = out_dir,
      output = "bad-scan-whitelist",
      global_whitelist = bad_whitelist, auto_whitelist = TRUE,
      threads = 1L, chunk_size = 1L, min_read_length = 0L
    ),
    "scan whitelist.*gzip|gzip.*scan whitelist|gzip text read failed"
  )
})

test_that("derived outputs cannot overwrite any explicit input artifact", {
  fixture_dir <- test_path("fixtures", "rad-core")
  fastq_lines <- readLines(file.path(fixture_dir, "tiny.fastq"))
  plain_fastq <- file.path(fixture_dir, "tiny.fastq")
  work <- tempfile("rad-core-collision-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  write_fastq_gzip <- function(path) {
    connection <- gzfile(path, open = "wb")
    on.exit(close(connection), add = TRUE)
    writeLines(fastq_lines, connection)
  }

  primary_input <- file.path(work, "reads.fq.gz")
  write_fastq_gzip(primary_input)
  primary_hash <- unname(tools::md5sum(primary_input))
  expect_error(
    rad_demux(
      "nanopore_rapid_bc", primary_input, work, output = "reads",
      threads = 1L, min_read_length = 0L
    ),
    "primary FASTQ would overwrite FASTQ input"
  )
  expect_identical(unname(tools::md5sum(primary_input)), primary_hash)

  linked_input <- file.path(work, "linked-source.fq.gz")
  linked_output <- file.path(work, "linked-run.fq.gz")
  write_fastq_gzip(linked_input)
  expect_true(file.link(linked_input, linked_output))
  linked_hash <- unname(tools::md5sum(linked_input))
  expect_error(
    rad_demux(
      "nanopore_rapid_bc", linked_input, work, output = "linked-run",
      threads = 1L, min_read_length = 0L
    ),
    "primary FASTQ would overwrite FASTQ input"
  )
  expect_identical(unname(tools::md5sum(linked_input)), linked_hash)

  debug_input <- file.path(work, "run_dbg.fq.gz")
  write_fastq_gzip(debug_input)
  debug_hash <- unname(tools::md5sum(debug_input))
  expect_error(
    rad_demux(
      "nanopore_rapid_bc", debug_input, work, output = "run",
      threads = 1L, min_read_length = 0L, write_debug = TRUE
    ),
    "debug FASTQ would overwrite FASTQ input"
  )
  expect_identical(unname(tools::md5sum(debug_input)), debug_hash)

  layout_input <- file.path(work, "layout-run_layout.csv")
  file.copy(file.path(fixture_dir, "tiny-layout.csv"), layout_input)
  layout_hash <- unname(tools::md5sum(layout_input))
  expect_error(
    rad_demux(
      layout_input, plain_fastq, work, output = "layout-run",
      custom_whitelist = file.path(fixture_dir, "tiny-whitelist.txt"),
      threads = 1L, min_read_length = 0L
    ),
    "layout cache would overwrite layout input"
  )
  expect_identical(unname(tools::md5sum(layout_input)), layout_hash)

  scan_reference <- file.path(work, "scan-run_scanwl.txt")
  writeLines(c("AACC", "TTGA", "CCCC"), scan_reference)
  reference_hash <- unname(tools::md5sum(scan_reference))
  expect_error(
    rad_demux(
      file.path(fixture_dir, "tiny-layout.csv"), plain_fastq, work,
      output = "scan-run", global_whitelist = scan_reference,
      auto_whitelist = TRUE, threads = 1L, min_read_length = 0L
    ),
    "auto-whitelist list would overwrite global whitelist input"
  )
  expect_identical(unname(tools::md5sum(scan_reference)), reference_hash)

  layout_reference <- file.path(work, "layout-run_scanwl.txt")
  writeLines(c("AACC", "TTGA", "CCCC"), layout_reference)
  derived_layout <- file.path(work, "derived-layout.csv")
  layout_lines <- readLines(file.path(fixture_dir, "tiny-layout.csv"))
  layout_lines[4] <- paste0(
    "barcode,,4,variable,barcode,", normalizePath(layout_reference)
  )
  writeLines(layout_lines, derived_layout)
  layout_reference_hash <- unname(tools::md5sum(layout_reference))
  expect_error(
    rad_demux(
      derived_layout, plain_fastq, work, output = "layout-run",
      auto_whitelist = TRUE, threads = 1L, min_read_length = 0L
    ),
    "auto-whitelist list would overwrite layout whitelist input"
  )
  expect_identical(
    unname(tools::md5sum(layout_reference)), layout_reference_hash
  )

  joint_mask <- file.path(work, "joint-run_whitelist_true.csv")
  writeLines("mask-input-must-survive", joint_mask)
  joint_layout <- file.path(work, "joint-layout.csv")
  writeLines(c(
    "Read Layout,,,,,,",
    "id,seq,expected_length,type,class,whitelist,flags",
    "anchor,ACGTACGT,,static,linker,,",
    paste0(
      "barcode,,4,variable,barcode,",
      normalizePath(file.path(fixture_dir, "tiny-whitelist.txt")), ":",
      normalizePath(joint_mask), ",joint_barcode"
    ),
    "poly_t,TTTTTTTT,,static,poly_t,,",
    "read,,,variable,read,,"
  ), joint_layout)
  joint_mask_hash <- unname(tools::md5sum(joint_mask))
  expect_error(
    rad_demux(
      joint_layout, plain_fastq, work, output = "joint-run",
      global_whitelist = file.path(fixture_dir, "tiny-whitelist.txt"),
      threads = 1L, min_read_length = 0L
    ),
    "true-whitelist summary would overwrite layout whitelist input"
  )
  expect_identical(unname(tools::md5sum(joint_mask)), joint_mask_hash)
})

test_that("output open failures cannot be reported as success", {
  fastq <- test_path("fixtures", "rad-core", "nanopore-smoke.fastq")
  out_dir <- tempfile("rad-core-output-failure-")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  dir.create(out_dir)
  dir.create(file.path(out_dir, "blocked.fq.gz"))

  expect_error(
    rad_demux(
      layout = "nanopore_rapid_bc", fastq = fastq,
      out_dir = out_dir, output = "blocked", threads = 1L,
      chunk_size = 1L, min_read_length = 0L
    ),
    "Failed to open file|output.*failed"
  )

  dir.create(file.path(out_dir, "layout-blocked_layout.csv"))
  expect_error(
    rad_demux(
      layout = "nanopore_rapid_bc", fastq = fastq,
      out_dir = out_dir, output = "layout-blocked", threads = 1L,
      chunk_size = 1L, min_read_length = 0L
    ),
    "Failed to install layout output"
  )

  if (file.exists("/dev/full")) {
    expect_true(file.symlink(
      "/dev/full", file.path(out_dir, "full.fq.gz")
    ))
    expect_error(
      rad_demux(
        layout = "nanopore_rapid_bc", fastq = fastq,
        out_dir = out_dir, output = "full", threads = 1L,
        chunk_size = 1L, min_read_length = 0L
      ),
      "gzip close failed|output (write|flush) failed"
    )
  }
  if (file.exists("/dev/null")) {
    expect_true(file.symlink(
      "/dev/null", file.path(out_dir, "null.fq.gz")
    ))
    expect_error(
      rad_demux(
        layout = "nanopore_rapid_bc", fastq = fastq,
        out_dir = out_dir, output = "null", threads = 1L,
        chunk_size = 1L, min_read_length = 0L
      ),
      "primary demultiplexed gzip is missing or empty"
    )
  }
})

test_that("an empty barcode whitelist fails safely in a subprocess", {
  layout <- normalizePath(test_path(
    "fixtures", "rad-core", "tiny-layout.csv"
  ))
  fastq <- normalizePath(test_path("fixtures", "rad-core", "tiny.fastq"))
  out_dir <- tempfile("rad-core-empty-wl-subprocess-")
  script <- tempfile("rad-core-empty-wl-", fileext = ".R")
  output <- tempfile("rad-core-empty-wl-", fileext = ".log")
  on.exit(unlink(c(out_dir, script, output), recursive = TRUE), add = TRUE)

  r_string <- function(value) paste(deparse(value), collapse = "")
  writeLines(c(
    "library(rrad)",
    sprintf("layout <- %s", r_string(layout)),
    sprintf("fastq <- %s", r_string(fastq)),
    sprintf("out_dir <- %s", r_string(out_dir)),
    paste0(
      "tryCatch({ rad_demux(layout, fastq, out_dir, threads=1L); ",
      "quit(status=2L) }, error=function(e) { ",
      "cat(conditionMessage(e)); quit(status=0L) })"
    )
  ), script)
  status <- system2(
    file.path(R.home("bin"), "Rscript"), script,
    stdout = output, stderr = output
  )

  expect_identical(status, 0L)
  expect_match(paste(readLines(output, warn = FALSE), collapse = "\n"),
               "No whitelist was provided")
})

test_that("OpenMP sequence failures return through R instead of terminating", {
  skip_if_not(isTRUE(rad_core_info()$build$openmp))
  good <- readLines(test_path(
    "fixtures", "rad-core", "nanopore-smoke.fastq"
  ))
  bad_fastq <- tempfile("rad-core-openmp-bad-", fileext = ".fastq")
  out_dir <- tempfile("rad-core-openmp-subprocess-")
  script <- tempfile("rad-core-openmp-", fileext = ".R")
  output <- tempfile("rad-core-openmp-", fileext = ".log")
  on.exit(
    unlink(c(bad_fastq, out_dir, script, output), recursive = TRUE),
    add = TRUE
  )
  writeLines(c(good, "@truncated", "ACGT", "+"), bad_fastq)

  r_string <- function(value) paste(deparse(normalizePath(
    value, mustWork = FALSE
  )), collapse = "")
  writeLines(c(
    "library(rrad)",
    sprintf("fastq <- %s", r_string(bad_fastq)),
    sprintf("out_dir <- %s", r_string(out_dir)),
    paste0(
      "tryCatch({ rad_demux('nanopore_rapid_bc', fastq, out_dir, ",
      "threads=2L, chunk_size=1L, min_read_length=0L); quit(status=2L) ",
      "}, error=function(e) { cat(conditionMessage(e)); quit(status=0L) })"
    )
  ), script)
  status <- system2(
    file.path(R.home("bin"), "Rscript"), script,
    stdout = output, stderr = output
  )

  expect_identical(status, 0L)
  expect_match(paste(readLines(output, warn = FALSE), collapse = "\n"),
               "sequence read failed")
})

test_that("explicit whitelist arguments preserve their declared roles", {
  skip_on_cran()
  layout <- test_path("fixtures", "rad-core", "tiny-layout.csv")
  fastq <- test_path("fixtures", "rad-core", "tiny.fastq")
  work <- tempfile("rad-core-typed-whitelists-")
  small_global <- tempfile("rad-core-global-", fileext = ".txt")
  equal_global <- tempfile("rad-core-equal-global-", fileext = ".txt")
  equal_custom <- tempfile("rad-core-equal-custom-", fileext = ".txt")
  large_custom <- tempfile("rad-core-large-custom-", fileext = ".txt")
  on.exit(
    unlink(c(work, small_global, equal_global, equal_custom, large_custom),
           recursive = TRUE),
    add = TRUE
  )
  writeLines(c("AACC", "TTGA"), small_global)
  writeLines(c("GGGG", "CCCC"), equal_global)
  writeLines(c("AACC", "TTGA"), equal_custom)

  encode_base4 <- function(value, width = 9L) {
    alphabet <- c("A", "C", "G", "T")
    digits <- character(width)
    for (position in rev(seq_len(width))) {
      digits[position] <- alphabet[(value %% 4L) + 1L]
      value <- value %/% 4L
    }
    paste0(digits, collapse = "")
  }
  large_values <- vapply(
    0:99999, encode_base4, character(1), USE.NAMES = FALSE
  )
  writeLines(c("AACC", large_values), large_custom)

  run <- function(name, ...) {
    rad_demux(
      layout = layout, fastq = fastq, out_dir = work, output = name,
      threads = 1L, chunk_size = 1L, max_reads = 1L,
      whitelist_mutation = 0L, generated_mutation = 0L,
      min_read_length = 0L, ...
    )
  }

  global <- run("small-global", global_whitelist = small_global)
  expect_equal(global$stats$loaded_global_barcodes, 2)
  expect_equal(global$stats$loaded_true_barcodes, 0)

  custom <- run("large-custom", custom_whitelist = large_custom)
  expect_equal(custom$stats$loaded_true_barcodes, 100001)
  expect_equal(custom$stats$loaded_global_barcodes, 0)

  paired <- run(
    "equal-pair", global_whitelist = equal_global,
    custom_whitelist = equal_custom
  )
  expect_equal(paired$stats$loaded_global_barcodes, 2)
  expect_equal(paired$stats$loaded_true_barcodes, 2)
})

test_that("explicit whitelists must contain valid barcodes", {
  layout <- test_path("fixtures", "rad-core", "tiny-layout.csv")
  fastq <- test_path("fixtures", "rad-core", "tiny.fastq")
  out_dir <- tempfile("rad-core-empty-typed-")
  empty <- tempfile("rad-core-empty-", fileext = ".txt")
  invalid <- tempfile("rad-core-invalid-", fileext = ".txt")
  on.exit(unlink(c(out_dir, empty, invalid), recursive = TRUE), add = TRUE)
  file.create(empty)
  writeLines(c("not-a-barcode", "123-not-dna"), invalid)

  expect_error(
    rad_demux(
      layout, fastq, out_dir, output = "empty-global",
      global_whitelist = empty, threads = 1L
    ),
    "global whitelist contains no valid barcodes"
  )
  expect_error(
    rad_demux(
      layout, fastq, out_dir, output = "invalid-custom",
      custom_whitelist = invalid, threads = 1L
    ),
    "custom whitelist contains no valid barcodes"
  )
})

test_that("explicit whitelist parsing trims rows and bounds packed values", {
  layout <- test_path("fixtures", "rad-core", "tiny-layout.csv")
  fastq <- test_path("fixtures", "rad-core", "tiny.fastq")
  out_dir <- tempfile("rad-core-strict-whitelist-")
  padded <- tempfile("rad-core-padded-", fileext = ".txt")
  out_of_range <- tempfile("rad-core-range-", fileext = ".txt")
  on.exit(
    unlink(c(out_dir, padded, out_of_range), recursive = TRUE),
    add = TRUE
  )

  con <- file(padded, open = "wb")
  writeBin(charToRaw("  AACC  \r\n"), con)
  close(con)
  writeLines("256", out_of_range)

  padded_result <- rad_demux(
    layout, fastq, out_dir, output = "padded",
    custom_whitelist = padded, threads = 1L,
    chunk_size = 1L, min_read_length = 0L
  )
  expect_equal(padded_result$stats$loaded_true_barcodes, 1)
  expect_identical(padded_result$status, "success")

  expect_error(
    rad_demux(
      layout, fastq, out_dir, output = "range",
      custom_whitelist = out_of_range, threads = 1L
    ),
    "custom whitelist contains no valid barcodes"
  )
})

test_that("rad_core_info rejects incompatible or incomplete metadata", {
  original <- rad_core_info_cpp()
  incompatible <- original
  incompatible$core_api_version <- 2L
  local_mocked_bindings(
    rad_core_info_cpp = function() incompatible,
    .package = "rrad"
  )
  expect_error(rad_core_info(), "expects 1 but core reports 2")

  missing_feature <- original
  missing_feature$features <- "file-fastq"
  local_mocked_bindings(
    rad_core_info_cpp = function() missing_feature,
    .package = "rrad"
  )
  expect_error(rad_core_info(), "does not advertise the demux feature")

  incomplete <- original
  incomplete$source_commit <- NULL
  local_mocked_bindings(
    rad_core_info_cpp = function() incomplete,
    .package = "rrad"
  )
  expect_error(rad_core_info(), "invalid metadata.*source_commit")

  invalid_digest <- original
  invalid_digest$embedded_source_digest <- "not-a-digest"
  local_mocked_bindings(
    rad_core_info_cpp = function() invalid_digest,
    .package = "rrad"
  )
  expect_error(rad_core_info(), "invalid embedded_source_digest")

  unavailable <- original
  unavailable$build$resources_available <- FALSE
  local_mocked_bindings(
    rad_core_info_cpp = function() unavailable,
    .package = "rrad"
  )
  expect_error(rad_core_info(), "bundled resources are unavailable")
})

test_that("rad_demux validates inputs before native execution", {
  fastq <- test_path("fixtures", "rad-core", "tiny.fastq")
  layout <- test_path("fixtures", "rad-core", "tiny-layout.csv")
  out_dir <- tempfile("rad-core-validation-")

  local_mocked_bindings(
    rad_demux_cpp = function(options) stop("unexpected native execution"),
    .package = "rrad"
  )

  expect_error(rad_demux(layout, tempfile(), out_dir),
               "fastq does not identify an existing file")
  expect_error(rad_demux(dirname(layout), fastq, out_dir),
               "layout must identify a layout file")
  expect_error(rad_demux(file.path(tempdir(), "missing-layout.csv"),
                         fastq, out_dir),
               "layout file does not exist")
  expect_error(rad_demux(layout, fastq, out_dir, output = "../escape"),
               "filename prefix within out_dir")
  expect_error(rad_demux(layout, fastq, out_dir, threads = 0), "threads")
  expect_error(rad_demux(layout, fastq, out_dir, chunk_size = 1.5),
               "chunk_size")
  expect_error(rad_demux(layout, fastq, out_dir, max_reads = 0), "max_reads")
  expect_error(rad_demux(layout, fastq, out_dir, auto_whitelist = 1),
               "auto_whitelist must be TRUE or FALSE")
  expect_error(rad_demux(layout, fastq, out_dir, kit = ":"),
               "kit must be one built-in kit key")
  expect_error(
    rad_demux(
      layout, fastq, out_dir, auto_whitelist = TRUE,
      scan_barcode_length = 65L
    ),
    "scan_barcode_length must be between 1 and 64"
  )
  expect_error(
    rad_demux(
      layout, fastq, out_dir, auto_whitelist = TRUE,
      scan_max_error = 1.01
    ),
    "scan_max_error must be between 0 and 1"
  )
  expect_error(rad_demux(layout, fastq, out_dir, whitelist_mutation = -1),
               "whitelist_mutation")
  expect_error(rad_demux(layout, fastq, out_dir,
                         bc_correction_mode = "optimistic"),
               "offensive.*defensive")
  expect_error(rad_demux(layout, fastq, out_dir, joint_bc_mode = "loose"),
               "default.*strict")
  expect_error(rad_demux(layout, fastq, out_dir, rc_umi = NA), "rc_umi")
  expect_error(rad_demux(layout, fastq, out_dir, min_read_length = -1),
               "min_read_length")
  expect_error(rad_demux(layout, fastq, out_dir, write_debug = NA),
               "write_debug")
  expect_error(rad_demux(layout, fastq, out_dir, reuse_cache = NA),
               "reuse_cache")
  expect_error(
    rad_demux(
      layout, fastq, out_dir, kit = "kit",
      global_whitelist = test_path(
        "fixtures", "rad-core", "tiny-whitelist.txt"
      )
    ),
    "kit and global_whitelist are alternatives"
  )
  expect_error(
    rad_demux(
      layout, fastq, out_dir, auto_whitelist = TRUE,
      custom_whitelist = test_path(
        "fixtures", "rad-core", "tiny-whitelist.txt"
      )
    ),
    "auto_whitelist and custom_whitelist are alternatives"
  )
})

test_that("automatic whitelist scanning requires core support", {
  fastq <- test_path("fixtures", "rad-core", "tiny.fastq")
  layout <- test_path("fixtures", "rad-core", "tiny-layout.csv")
  core <- rad_core_info_cpp()
  core$features <- setdiff(core$features, "auto-whitelist")
  local_mocked_bindings(
    rad_core_info_cpp = function() core,
    rad_demux_cpp = function(options) stop("unexpected native execution"),
    .package = "rrad"
  )
  expect_error(
    rad_demux(layout, fastq, tempfile("rad-core-auto-"), auto_whitelist = TRUE),
    "does not advertise the auto-whitelist feature"
  )
})

test_that("native automatic scan limits reject narrowing overflow", {
  fastq <- test_path("fixtures", "rad-core", "tiny.fastq")
  layout <- test_path("fixtures", "rad-core", "tiny-layout.csv")
  base <- list(
    layout = layout,
    fastq = fastq,
    output_dir = tempfile("rad-core-native-limit-"),
    output_prefix = "limit",
    auto_whitelist = TRUE
  )

  expect_error(
    rad_demux_cpp(c(base, list(scan_max_reads = 2147483648))),
    "scan_max_reads"
  )
  expect_error(
    rad_demux_cpp(c(base, list(scan_chunk_size = 2147483648))),
    "scan_chunk_size"
  )
})
