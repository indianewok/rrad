test_that("demux failures preserve the complete previously published artifact set", {
  fixture_dir <- test_path("fixtures", "rad-core")
  layout <- file.path(fixture_dir, "tiny-layout.csv")
  whitelist <- file.path(fixture_dir, "tiny-whitelist.txt")
  input <- file.path(fixture_dir, "tiny.fastq")
  work <- tempfile("rrad-demux-transaction-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  first <- rad_demux(
    layout = layout,
    fastq = input,
    out_dir = work,
    output = "stable",
    custom_whitelist = whitelist,
    threads = 1L,
    chunk_size = 1L,
    min_read_length = 0L,
    write_debug = TRUE
  )
  published <- unlist(first$artifacts, use.names = FALSE)
  published <- published[nzchar(published) & file.exists(published)]
  expect_true(length(published) >= 7L)
  before <- unname(tools::md5sum(published))

  malformed <- file.path(work, "malformed.fastq")
  writeLines(c(readLines(input), "@truncated", "ACGT", "+"), malformed)
  expect_error(
    rad_demux(
      layout = layout,
      fastq = malformed,
      out_dir = work,
      output = "stable",
      custom_whitelist = whitelist,
      threads = 1L,
      chunk_size = 1L,
      min_read_length = 0L
    ),
    "sequence read failed"
  )

  expect_true(all(file.exists(published)))
  expect_identical(unname(tools::md5sum(published)), before)
  expect_length(
    list.files(work, pattern = "^\\.rrad-demux-stage-", all.files = TRUE),
    0L
  )
})

test_that("successful demux publication removes stale optional artifacts", {
  fixture_dir <- test_path("fixtures", "rad-core")
  layout <- file.path(fixture_dir, "tiny-layout.csv")
  whitelist <- file.path(fixture_dir, "tiny-whitelist.txt")
  input <- file.path(fixture_dir, "tiny.fastq")
  work <- tempfile("rrad-demux-artifact-set-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  debug <- rad_demux(
    layout, input, work,
    output = "replace",
    custom_whitelist = whitelist,
    threads = 1L,
    chunk_size = 1L,
    min_read_length = 0L,
    write_debug = TRUE
  )
  debug_paths <- unlist(debug$artifacts[c(
    "debug_sig", "debug_csv", "debug_fastq", "metrics"
  )], use.names = FALSE)
  expect_true(all(file.exists(debug_paths)))

  ordinary <- rad_demux(
    layout, input, work,
    output = "replace",
    custom_whitelist = whitelist,
    threads = 1L,
    chunk_size = 1L,
    min_read_length = 0L,
    write_debug = FALSE
  )
  expect_true(file.exists(ordinary$artifacts$fastq))
  expect_false(any(file.exists(debug_paths)))
  expect_length(
    list.files(work, pattern = "^\\.rrad-demux-stage-", all.files = TRUE),
    0L
  )
})

test_that("demux publication waits on one interprocess transaction lock", {
  python <- Sys.which("python3")
  skip_if(!nzchar(python), "requires python3 for an external fcntl lock")

  fixture_dir <- test_path("fixtures", "rad-core")
  layout <- file.path(fixture_dir, "tiny-layout.csv")
  whitelist <- file.path(fixture_dir, "tiny-whitelist.txt")
  input <- file.path(fixture_dir, "tiny.fastq")
  work <- tempfile("rrad-demux-lock-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  prefix <- file.path(work, "stable")
  targets <- paste0(prefix, c(
    ".fq.gz", "_layout.csv", "_position_map.csv",
    "_whitelist_true.csv", "_whitelist_global.csv", "_demux.log",
    "_scanwl.csv", "_scanwl.txt", "_scan_wl.log", "_dbg.sig.gz",
    "_dbg.csv.gz", "_dbg.fq.gz", ".metrics.tsv"
  ))
  for (index in seq_along(targets)) {
    writeLines(paste("old artifact", index), targets[[index]])
  }

  ready <- file.path(work, "holder-ready")
  observed <- file.path(work, "holder-observed")
  holder_log <- file.path(work, "holder.log")
  lock_path <- file.path(work, ".stable.rrad-demux-commit.lock")
  holder <- file.path(work, "hold_lock.py")
  writeLines(c(
    "import fcntl, pathlib, sys, time",
    "lock_path, ready, observed = sys.argv[1:4]",
    "targets = [pathlib.Path(value) for value in sys.argv[4:]]",
    "before = [path.read_bytes() for path in targets]",
    "with open(lock_path, 'a+b') as lock:",
    "    fcntl.lockf(lock, fcntl.LOCK_EX)",
    "    pathlib.Path(ready).write_text('ready\\n')",
    "    deadline = time.monotonic() + 3.0",
    "    intact = True",
    "    while time.monotonic() < deadline:",
    "        intact = all(path.exists() and path.read_bytes() == expected",
    "                     for path, expected in zip(targets, before))",
    "        if not intact:",
    "            break",
    "        time.sleep(0.01)",
    "    pathlib.Path(observed).write_text(",
    "        'intact\\n' if intact else 'changed\\n')"
  ), holder)

  launch <- suppressWarnings(system2(
    python,
    vapply(c(holder, lock_path, ready, observed, targets),
           shQuote, character(1)),
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
  result <- rad_demux(
    layout, input, work, output = "stable",
    custom_whitelist = whitelist, threads = 1L, chunk_size = 1L,
    min_read_length = 0L, write_debug = FALSE
  )
  elapsed <- proc.time()[["elapsed"]] - started

  expect_identical(result$status, "success")
  expect_true(file.exists(observed))
  expect_identical(readLines(observed, warn = FALSE), "intact")
  expect_true(elapsed >= 2.5, info = paste("elapsed:", elapsed))
  expect_identical(as.character(file.info(lock_path)$mode), "600")
})

test_that("reuse_cache reads one locked snapshot of both cache artifacts", {
  python <- Sys.which("python3")
  skip_if(!nzchar(python), "requires python3 for an external fcntl lock")

  fixture_dir <- test_path("fixtures", "rad-core")
  layout <- file.path(fixture_dir, "tiny-layout.csv")
  whitelist <- file.path(fixture_dir, "tiny-whitelist.txt")
  input <- file.path(fixture_dir, "tiny.fastq")
  work <- tempfile("rrad-demux-cache-snapshot-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  first <- rad_demux(
    layout, input, work, output = "stable",
    custom_whitelist = whitelist, threads = 1L, chunk_size = 1L,
    min_read_length = 0L
  )
  expect_identical(first$stats$loaded_true_barcodes, 2)

  prefix <- file.path(work, "stable")
  layout_cache <- paste0(prefix, "_layout.csv")
  position_cache <- paste0(prefix, "_position_map.csv")
  cache_before <- unname(tools::md5sum(c(layout_cache, position_cache)))
  lock_path <- file.path(work, ".stable.rrad-demux-commit.lock")
  ready <- file.path(work, "writer-ready")
  writer_log <- file.path(work, "writer.log")
  writer <- file.path(work, "publish_torn_cache.py")
  writeLines(c(
    "import fcntl, pathlib, sys, time",
    "lock_path, layout_path, ready = sys.argv[1:4]",
    "layout = pathlib.Path(layout_path)",
    "original = layout.read_bytes()",
    "with open(lock_path, 'a+b') as lock:",
    "    fcntl.lockf(lock, fcntl.LOCK_EX)",
    "    layout.write_text('poisoned torn layout\\n')",
    "    pathlib.Path(ready).write_text('ready\\n')",
    "    time.sleep(1.5)",
    "    layout.write_bytes(original)"
  ), writer)

  launch <- suppressWarnings(system2(
    python,
    vapply(c(writer, lock_path, layout_cache, ready),
           shQuote, character(1)),
    stdout = writer_log,
    stderr = writer_log,
    wait = FALSE
  ))
  expect_identical(launch, 0L)
  deadline <- Sys.time() + 10
  while (!file.exists(ready) && Sys.time() < deadline) Sys.sleep(0.02)
  expect_true(
    file.exists(ready),
    info = paste(readLines(writer_log, warn = FALSE), collapse = "\n")
  )

  started <- proc.time()[["elapsed"]]
  reused <- rad_demux(
    layout, input, work, output = "stable",
    custom_whitelist = whitelist, threads = 1L, chunk_size = 1L,
    min_read_length = 0L, reuse_cache = TRUE
  )
  elapsed <- proc.time()[["elapsed"]] - started

  expect_identical(reused$status, "success")
  expect_identical(reused$stats$loaded_true_barcodes, 2)
  expect_true(elapsed >= 1.2, info = paste("elapsed:", elapsed))
  expect_identical(
    unname(tools::md5sum(c(layout_cache, position_cache))),
    cache_before
  )
  expect_length(
    list.files(work, pattern = "^\\.rrad-demux-stage-", all.files = TRUE),
    0L
  )
})

test_that("native scan error ratios are bounded before entering Edlib", {
  fixture_dir <- test_path("fixtures", "rad-core")
  scan_cpp <- getFromNamespace(".rad_scan_wl_cpp", "rrad")
  demux_cpp <- getFromNamespace("rad_demux_cpp", "rrad")

  expect_error(
    scan_cpp(list(
      input = file.path(fixture_dir, "tiny.fastq"),
      output_prefix = tempfile("rrad-scan-bound-"),
      adapter_seq = "ACGT",
      max_error = 1.01
    )),
    "`max_error` must be between 0 and 1"
  )
  expect_error(
    scan_cpp(list(
      input = file.path(fixture_dir, "tiny.fastq"),
      output_prefix = tempfile("rrad-scan-finite-bound-"),
      adapter_seq = "ACGT",
      max_error = NaN
    )),
    "`max_error` must be finite"
  )
  expect_error(
    scan_cpp(list(
      input = file.path(fixture_dir, "tiny.fastq"),
      output_prefix = tempfile("rrad-scan-length-bound-"),
      adapter_seq = "ACGT",
      barcode_length = 33L
    )),
    "`barcode_length` must be between 1 and 32"
  )
  expect_error(
    demux_cpp(list(
      layout = file.path(fixture_dir, "tiny-layout.csv"),
      fastq = file.path(fixture_dir, "tiny.fastq"),
      output_dir = tempfile("rrad-demux-bound-"),
      output_prefix = "bounded",
      auto_whitelist = TRUE,
      scan_max_error = 1.01
    )),
    "`scan_max_error` must be between 0 and 1"
  )
  expect_error(
    demux_cpp(list(
      layout = file.path(fixture_dir, "tiny-layout.csv"),
      fastq = file.path(fixture_dir, "tiny.fastq"),
      output_dir = tempfile("rrad-demux-length-bound-"),
      output_prefix = "bounded",
      auto_whitelist = TRUE,
      scan_barcode_length = 33L
    )),
    "`scan_barcode_length` must be between 1 and 32"
  )
  expect_error(
    demux_cpp(list(
      layout = file.path(fixture_dir, "tiny-layout.csv"),
      fastq = file.path(fixture_dir, "tiny.fastq"),
      output_dir = tempfile("rrad-demux-int64-bound-"),
      output_prefix = "bounded",
      auto_whitelist = TRUE,
      scan_max_reads = 2^63
    )),
    "`scan_max_reads` must be -1 or a non-negative whole number"
  )
})
