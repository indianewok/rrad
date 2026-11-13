scan_port_fasta <- function(path, sequences) {
  writeLines(as.vector(rbind(paste0(">read-", seq_along(sequences)), sequences)), path)
}

scan_port_reverse <- function(sequence) {
  paste(rev(strsplit(chartr("ACGT", "TGCA", sequence), "")[[1L]]), collapse = "")
}

scan_port_lines <- function(path) {
  connection <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection))
  readLines(connection, warn = FALSE)
}

test_that("masked captures concatenate blocks and preserve ordinary N behavior", {
  work <- tempfile("rrad-masked-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- file.path(work, "reads.fa")
  sequence <- "ACAGTCCGTATTGA"
  scan_port_fasta(input, c(sequence, scan_port_reverse(sequence),
                          "ANAGTCCGTATTGA", "ACAGTNCGTATTGA"))
  result <- rad_scan_wl(input, file.path(work, "masked"),
    adapter_seq = "NNAGTCCGTAnnNN", capture_mode = "masked",
    max_error = 0, selection = "above_floor", chunk_size = 1L)
  expect_identical(scan_port_lines(result$artifacts$whitelist), "ACGA")
  expect_equal(result$stats$total_extractions, 2)
  expect_equal(result$config$barcode_length, 4L)
  expect_identical(result$config$capture_mode, "masked")

  reference <- file.path(work, "reference.txt")
  writeLines("ACGA", reference)
  packed <- rad_scan_wl(input, file.path(work, "packed"),
    adapter_seq = "NNAGTCCGTAnnNN", capture_mode = "masked", whitelist = reference,
    max_error = 0, selection = "above_floor", chunk_size = 1L)
  expect_identical(scan_port_lines(packed$artifacts$whitelist), "ACGA")
  expect_equal(packed$stats$total_perfect_matches, 2)

  scan_port_fasta(input, "AGTCGACA")
  ordinary <- rad_scan_wl(input, file.path(work, "ordinary"),
    adapter_seq = "AGNC", barcode_length = 4L, max_error = 0,
    selection = "above_floor")
  masked <- rad_scan_wl(input, file.path(work, "one-mask"),
    adapter_seq = "AGNC", capture_mode = "masked", max_error = 0,
    selection = "above_floor")
  expect_identical(scan_port_lines(ordinary$artifacts$whitelist), "GACA")
  expect_identical(scan_port_lines(masked$artifacts$whitelist), "T")
  expect_error(rad_scan_wl(input, file.path(work, "bad"),
    adapter_seq = "AGNC", capture_mode = "masked", barcode_length = 4L),
    "disagrees")
  expect_error(rad_scan_wl(input, file.path(work, "bad"),
    adapter_seq = "AGNC", capture_mode = "masked", left_margin = 1L),
    "margins")
})

test_that("masked fuzzy capture follows literal-region indels on both strands", {
  work <- tempfile("rrad-masked-indel-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- file.path(work, "reads.fa")
  sequences <- c("ACAGTCACGTAGA", "ACAGTCGTAGA")
  scan_port_fasta(input, c(sequences, vapply(sequences, scan_port_reverse, "")))
  result <- rad_scan_wl(input, file.path(work, "calls"),
    adapter_seq = "NNAGTCCGTANN", capture_mode = "masked",
    max_error = 0.13, selection = "above_floor", chunk_size = 1L)
  expect_identical(scan_port_lines(result$artifacts$whitelist), "ACGA")
  expect_equal(result$stats$total_extractions, 4)
})

test_that("masked capture rejects gaps at captured block boundaries", {
  work <- tempfile("rrad-masked-boundary-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- file.path(work, "reads.fa")
  sequences <- c("AGTCCGTAAACTGACCTAG", "AGTCCGTAAAACCTGACCTAG")
  scan_port_fasta(input, c(sequences, vapply(sequences, scan_port_reverse, "")))
  result <- rad_scan_wl(input, file.path(work, "calls"),
    adapter_seq = "AGTCCGTANNNNTGACCTAG", capture_mode = "masked",
    max_error = 0.07, selection = "above_floor")
  expect_equal(result$stats$total_extractions, 0)
  expect_length(scan_port_lines(result$artifacts$whitelist), 0L)
})

test_that("64-base scan references retain high packed bits and sequence length", {
  work <- tempfile("rrad-long-scan-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- file.path(work, "reads.fa")
  reference <- file.path(work, "reference.txt")
  barcode <- paste(rep("G", 64), collapse = "")
  scan_port_fasta(input, c(paste0("AGTCCGTA", barcode),
                          scan_port_reverse(paste0("AGTCCGTA", barcode))))
  writeLines("340282366920938463463374607431768211455", reference)
  result <- rad_scan_wl(input, file.path(work, "calls"),
    adapter_seq = "AGTCCGTA", barcode_length = 64L, whitelist = reference,
    max_error = 0, selection = "above_floor")
  expect_identical(scan_port_lines(result$artifacts$whitelist), barcode)
  expect_equal(result$stats$total_perfect_matches, 2)
  writeLines(barcode, reference)
  dna <- rad_scan_wl(input, file.path(work, "dna"),
    adapter_seq = "AGTCCGTA", barcode_length = 64L, whitelist = reference,
    max_error = 0, selection = "above_floor")
  expect_identical(scan_port_lines(dna$artifacts$whitelist), barcode)
  masked <- rad_scan_wl(input, file.path(work, "masked-long"),
    adapter_seq = paste0("AGTCCGTA", paste(rep("N", 64), collapse = "")),
    capture_mode = "masked", whitelist = reference,
    max_error = 0, selection = "above_floor")
  expect_identical(scan_port_lines(masked$artifacts$whitelist), barcode)
  expect_equal(masked$config$barcode_length, 64L)
})

test_that("split masked patterns use lowercase n without capturing it", {
  work <- tempfile("rrad-split-mask-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- file.path(work, "reads.fa")
  first <- file.path(work, "bc1.txt")
  second <- file.path(work, "bc2.txt")
  writeLines("GACA", first)
  writeLines("TTGC", second)
  molecule <- "AGTCCGTATGACATTGC"
  scan_port_fasta(input, c(molecule, scan_port_reverse(molecule)))
  result <- rad_scan_wl(input, file.path(work, "calls"),
    adapter_seq = "AGTCCGTAn", capture_mode = "masked",
    bc1_whitelist = first, bc2_whitelist = second, umi_length = 0L,
    offset_min = 0L, offset_max = 0L, max_error = 0,
    selection = "above_floor")
  expect_identical(scan_port_lines(result$artifacts$whitelist), "GACATTGC")
  expect_equal(result$stats$total_extractions, 2)
})

test_that("spatial low-count rescan retains the one-observation floor", {
  work <- tempfile("rrad-spatial-floor-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- file.path(work, "counts.csv")
  writeLines(c("barcode,count", "AAAA,1000000", "CCCC,1"), input)
  ordinary <- rad_scan_wl(input, file.path(work, "ordinary"),
    rescan = TRUE, selection = "above_floor")
  spatial <- rad_scan_wl(input, file.path(work, "spatial"),
    rescan = TRUE, spatial = TRUE, selection = "above_floor")
  expect_identical(scan_port_lines(ordinary$artifacts$whitelist), "AAAA")
  expect_setequal(scan_port_lines(spatial$artifacts$whitelist), c("AAAA", "CCCC"))
  expect_lt(spatial$stats$floor, ordinary$stats$floor)
})

test_that("slide membership refines abundance calls and exports coordinates", {
  work <- tempfile("rrad-spatial-map-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- file.path(work, "counts.csv")
  map <- file.path(work, "slide.tsv")
  writeLines(c("barcode,count", "AAAA,10", "CCCC,5", "GGGG,5"), input)
  writeLines(c("barcode\tx\ty", "AAAA\t0\t0", "CCCC\t2\t0"), map)
  result <- rad_scan_wl(input, file.path(work, "calls"), rescan = TRUE,
    spatial_slide = map, spatial_platform = "curio-seeker", selection = "above_floor")
  expect_setequal(scan_port_lines(result$artifacts$whitelist), c("AAAA", "CCCC"))
  expect_setequal(scan_port_lines(result$artifacts$abundance_whitelist), c("AAAA", "CCCC", "GGGG"))
  expect_length(scan_port_lines(result$artifacts$coordinates), 3L)
  expect_equal(result$stats$selected_barcodes, 2)
  expect_equal(result$stats$input_rows, 3)

  exported <- rad_scan_wl(output_prefix = file.path(work, "catalog"),
    spatial_slide = map, spatial_platform = "curio-seeker")
  expect_length(scan_port_lines(exported$artifacts$coordinates), 3L)
  expect_identical(exported$input, "")
  expect_identical(exported$artifacts$whitelist, "")
  expect_error(rad_scan_wl(input, file.path(work, "bad"), rescan = TRUE,
    spatial_slide = map), "requires spatial_platform")
})

test_that("HD occupancy counts 32 distinct sites and supports native IDs", {
  work <- tempfile("rrad-hd-occupancy-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  first <- file.path(work, "columns.tsv")
  second <- file.path(work, "rows.tsv")
  bc1 <- c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT")
  bc2 <- c("GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT")
  writeLines(c("barcode\tarray_col", paste(bc1, 0:7, sep = "\t")), first)
  writeLines(c("barcode\tarray_row", paste(bc2, 0:7, sep = "\t")), second)
  input <- file.path(work, "counts.csv")
  ids <- sprintf("s_002um_%05d_%05d", rep(0:3, each = 8), rep(0:7, 4))
  run <- function(prefix) rad_scan_wl(input, file.path(work, prefix),
    rescan = TRUE, spatial_slide = "synthetic", spatial_platform = "visium-hd",
    bc1_axis = first, bc2_axis = second, spatial_selection = "legacy")
  writeLines(c("barcode,count", paste(ids, 100L, sep = ",")), input)
  selected <- run("selected")
  expect_length(scan_port_lines(selected$artifacts$whitelist), 32L)
  regions <- utils::read.delim(selected$artifacts$regions)
  expect_equal(regions$occupied_2um_spots, 32)
  expect_true(regions$selected)
  writeLines(c("barcode,count", paste(ids[-1], 1000L, sep = ",")), input)
  rejected <- run("rejected")
  expect_length(scan_port_lines(rejected$artifacts$whitelist), 0L)
  expect_equal(utils::read.delim(rejected$artifacts$regions)$occupied_2um_spots, 31)
})

test_that("spatial failures preserve published files and reject coordinate aliases", {
  work <- tempfile("rrad-spatial-atomic-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- file.path(work, "counts.csv")
  prefix <- file.path(work, "calls")
  map <- file.path(work, "slide.tsv.gz")
  writeLines(c("barcode,count", "AAAA,1"), input)
  connection <- gzfile(map, "wt")
  writeLines(c("barcode\tx\ty", "AAAA\t0\t0"), connection)
  close(connection)
  result <- rad_scan_wl(input, prefix, rescan = TRUE,
    spatial_slide = map, spatial_platform = "curio-seeker", selection = "above_floor")
  files <- unlist(result$artifacts, use.names = FALSE)
  files <- files[nzchar(files)]
  before <- tools::md5sum(files)
  bytes <- readBin(map, "raw", n = file.info(map)$size)
  writeBin(bytes[seq_len(length(bytes) - 8L)], map)
  expect_error(rad_scan_wl(input, prefix, rescan = TRUE,
    spatial_slide = map, spatial_platform = "curio-seeker"), "Slide stream|gzip|slide")
  expect_identical(tools::md5sum(files), before)
  alias <- paste0(prefix, "_coordinates.tsv.gz")
  expect_error(rad_scan_wl(input, prefix, rescan = TRUE,
    spatial_slide = "synthetic", slide_map = alias,
    spatial_platform = "curio-seeker"), "overwrite|alias|collid")
})


test_that("masked low-count batches validate every artifact before running", {
  work <- tempfile("rrad-mask-batch-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- file.path(work, "reads.fa")
  scan_port_fasta(input, "ACAGTCCGTATTGA")
  batch <- file.path(work, "batch.csv")
  prefixes <- file.path(work, c("one", "two"))
  utils::write.table(data.frame(input = rep(input, 2), output_prefix = prefixes),
    batch, sep = ",", row.names = FALSE, col.names = FALSE, quote = TRUE)
  result <- rad_scan_wl(batch_csv = batch, adapter_seq = "NNAGTCCGTAnnNN",
    capture_mode = "masked", spatial = TRUE, max_error = 0, selection = "above_floor")
  expect_length(result$runs, 2L)
  expect_true(all(vapply(result$runs, function(run) isTRUE(run$config$spatial), logical(1))))
  expect_true(all(vapply(result$runs, function(run) identical(
    scan_port_lines(run$artifacts$whitelist), "ACGA"), logical(1))))
  colliding <- file.path(work, c("collision", "collision_rad"))
  utils::write.table(data.frame(input = rep(input, 2), output_prefix = colliding),
    batch, sep = ",", row.names = FALSE, col.names = FALSE, quote = TRUE)
  expect_error(rad_scan_wl(batch_csv = batch, adapter_seq = "NNAGTCCGTAnnNN",
    capture_mode = "masked", max_error = 0), "colliding output artifacts")
  expect_false(file.exists(paste0(colliding[[1L]], ".csv")))
})

test_that("discovery lowercase-n primers match directionally without capturing", {
  work <- tempfile("rrad-primer-wildcard-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- file.path(work, "reads.fa")
  catalog <- file.path(work, "catalog.txt")
  scan_port_fasta(input, c("AGTACGTAGACA", "AGNACGTAGACA", "AGTNCGTAGACA"))
  writeLines("GACA", catalog)
  result <- rad_discover(input, work, pattern = "AGTnCGTA",
    barcode_length = 4L, catalogs = catalog, max_error = 0)
  expect_true(result$complete)
  expect_equal(result$whitelist_match$matched_reads, 1)
  expect_equal(result$whitelist_match$extracted_reads, 1)
})


test_that("explicit numeric discovery catalogs require a declared length", {
  work <- tempfile("rrad-discover-numeric-length-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  input <- file.path(work, "reads.fa")
  catalog <- file.path(work, "catalog.txt")
  scan_port_fasta(input, "AGTCCGTAGACA")
  writeLines("209", catalog)
  expect_error(rad_discover(input, work, pattern = "AGTCCGTA",
    barcode_length = 4L, catalogs = catalog, max_error = 0),
    "declared barcode length")
  expect_false(file.exists(file.path(work, "discovery_whitelist_match.csv")))
})
