test_that("complete layout discovery hands an offline local catalog to demux", {
  work <- tempfile("rrad-complete-layout-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  old_offline <- Sys.getenv("RRAD_OFFLINE", unset = NA_character_)
  Sys.setenv(RRAD_OFFLINE = "true")
  on.exit(if (is.na(old_offline)) Sys.unsetenv("RRAD_OFFLINE") else
    Sys.setenv(RRAD_OFFLINE = old_offline), add = TRUE)

  primer <- "CTACACGACGCTCTTCCGATCT"
  reverse_primer <- "ACTCTGCGTTGATACCACTGCTT"
  barcodes <- c("GACCTAGTTCAGGATC", "CAGTTCGATGACCTGA")
  umi <- "AACCGTACGTCG"
  reverse_complement <- function(sequence) paste(
    rev(strsplit(chartr("ACGT", "TGCA", sequence), "", fixed = TRUE)[[1L]]),
    collapse = ""
  )
  # Both boundaries and a bounded barcode/UMI gap establish the whole scaffold.
  # Varying insert lengths distinguish the insert from a fixed spacer.
  insert_template <- strrep("GACTCGACTAGCTACGATCGTAGCATGCTAGCTGACTAGCA", 4L)
  inserts <- vapply(c(60L, 80L, 100L, 120L), function(length) {
    substr(insert_template, 1L, length)
  }, character(1))
  forward <- paste0(primer, rep(barcodes, 2L), umi, strrep("T", 20L),
                    inserts, reverse_primer)
  sequences <- c(forward, vapply(forward, reverse_complement, character(1)))
  input <- file.path(work, "reads.fastq")
  catalog <- file.path(work, "local-catalog.txt")
  writeLines(barcodes, catalog)
  records <- unlist(lapply(seq_along(sequences), function(i) {
    c(paste0("@read", i), sequences[[i]], "+", strrep("I", nchar(sequences[[i]])))
  }), use.names = FALSE)
  writeLines(records, input)

  inferred <- rad_discover(
    input, file.path(work, "discovery"), "complete", mode = "both",
    catalogs = c(local = catalog), max_error = 0, threads = 1L
  )
  expect_true(inferred$success)
  expect_true(inferred$complete)
  expect_true(inferred$layout_complete)
  expect_true(inferred$whitelist_complete)
  expect_equal(inferred$sample_reads, length(sequences))
  expect_length(inferred$unavailable_catalogs, 0L)
  expect_null(inferred$files[["layout.partial.csv"]])
  expect_identical(inferred$whitelist_match$status, "matched")
  expect_equal(inferred$whitelist_match$matched_reads, length(sequences))
  expect_equal(inferred$whitelist_match$matched_distinct, length(barcodes))
  expect_identical(inferred$whitelist_match$catalog, normalizePath(catalog))
  primer_support <- inferred$parts[inferred$parts$sequence == primer, , drop = FALSE]
  expect_equal(primer_support$forward, 4)
  expect_equal(primer_support$reverse, 4)

  layout_path <- inferred$files[["layout.csv"]]
  expect_true(file.exists(layout_path))
  layout <- utils::read.csv(layout_path, skip = 1L, stringsAsFactors = FALSE)
  expect_equal(layout$expected_length[layout$class == "barcode"], 16)
  expect_equal(layout$expected_length[layout$class == "umi"], 12)
  expect_identical(layout$whitelist[layout$class == "barcode"], normalizePath(catalog))
  expect_true(any(layout$class == "read"))

  # The generated layout carries its explicit local catalog. It must work
  # without a further whitelist override or a remote resource lookup.
  demultiplexed <- rad_demux(
    layout_path, input, file.path(work, "inferred-demux"), threads = 1L,
    chunk_size = 4L, min_read_length = 0L,
    whitelist_mutation = 0L, generated_mutation = 0L
  )
  canonical_layout <- system.file(
    "rad", "resources", "read_layout", "three_prime_read_layout.csv",
    package = "rrad"
  )
  reference <- rad_demux(
    canonical_layout, input, file.path(work, "reference-demux"),
    custom_whitelist = catalog, threads = 1L, chunk_size = 4L,
    min_read_length = 0L, whitelist_mutation = 0L, generated_mutation = 0L
  )
  read_output <- function(path) {
    connection <- gzfile(path, "rt")
    on.exit(close(connection), add = TRUE)
    readLines(connection, warn = FALSE)
  }
  output <- read_output(demultiplexed$artifacts$fastq)
  expect_identical(output, read_output(reference$artifacts$fastq))
  expect_equal(demultiplexed$stats$loaded_true_barcodes, length(barcodes))
  expect_equal(demultiplexed$stats$records_serialized, length(sequences))
  headers <- output[seq.int(1L, length(output), 4L)]
  expected_barcodes <- rep(barcodes, 4L)
  for (i in seq_along(headers)) {
    expect_match(headers[[i]], paste0("CB:Z:", expected_barcodes[[i]]), fixed = TRUE)
    expect_match(headers[[i]], paste0("UB:Z:", umi), fixed = TRUE)
    expect_match(headers[[i]], if (i <= 4L) "TS:A:+" else "TS:A:-", fixed = TRUE)
  }
  expect_equal(nchar(output[seq.int(2L, length(output), 4L)]),
               nchar(output[seq.int(4L, length(output), 4L)]))
})
