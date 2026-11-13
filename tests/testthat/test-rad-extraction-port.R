port_reverse_complement <- function(sequence) {
  paste(rev(strsplit(chartr("ACGT", "TGCA", sequence), "", fixed = TRUE)[[1L]]),
        collapse = "")
}

port_write_fastq <- function(path, sequences) {
  lines <- unlist(lapply(seq_along(sequences), function(i) {
    sequence <- sequences[[i]]
    c(paste0("@", names(sequences)[[i]]), sequence, "+",
      paste(rep("I", nchar(sequence)), collapse = ""))
  }), use.names = FALSE)
  writeLines(lines, path)
}

port_read_fastq <- function(path) {
  connection <- gzfile(path, "rt")
  on.exit(close(connection), add = TRUE)
  readLines(connection, warn = FALSE)
}

test_that("barcode positioning moves a chained UMI in both orientations", {
  work <- tempfile("rrad-umi-position-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  anchor <- "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAA"
  anchor2 <- "CAGCACCTGATTACAAGGTGCTG"
  barcode <- "GACCTAGTTCAGGATC"
  umi <- "CGTATGCA"
  cdna <- "GACTCGACTAGCTACGATCGTAGCATGCTAGCTGACTAGCA"
  whitelist <- file.path(work, "barcodes.txt")
  layout <- file.path(work, "layout.csv")
  fastq <- file.path(work, "reads.fastq")
  writeLines(barcode, whitelist)
  writeLines(c(
    "Read Layout,,,,,", "id,seq,expected_length,type,class,whitelist",
    paste0("anchor,", anchor, ",,static,linker,"),
    paste0("barcode,,16,variable,barcode,", whitelist),
    "umi,,8,variable,umi,",
    paste0("anchor2,", anchor2, ",,static,linker,"),
    "read,,,variable,read,"
  ), layout)
  # This inserted base moves the real barcode and its adjacent UMI one base
  # beyond the layout prediction. Sequence and quality must move together.
  sequence <- paste0(anchor, "C", barcode, umi, anchor2, cdna)
  port_write_fastq(fastq, c(forward = sequence,
                          reverse = port_reverse_complement(sequence)))
  result <- rad_demux(layout, fastq, file.path(work, "out"), threads = 1L,
                      chunk_size = 2L, min_read_length = 0L,
                      whitelist_mutation = 0L, generated_mutation = 0L)
  records <- port_read_fastq(result$artifacts$fastq)
  headers <- records[c(1L, 5L)]
  expect_equal(result$stats$records_serialized, 2)
  expect_true(all(grepl(paste0("CB:Z:", barcode), headers, fixed = TRUE)))
  expect_true(all(grepl(paste0("UB:Z:", umi), headers, fixed = TRUE)))
  expect_match(headers[[1L]], "TS:A:+", fixed = TRUE)
  expect_match(headers[[2L]], "TS:A:-", fixed = TRUE)
  expect_equal(nchar(records[c(2L, 6L)]), nchar(records[c(4L, 8L)]))
  expect_identical(records[[2L]], cdna)
  expect_identical(records[[6L]], port_reverse_complement(cdna))
  counts <- utils::read.csv(result$artifacts$whitelist_true)
  expect_named(counts, c("class_id", "true_barcode", "raw_count", "total_count",
                        "corrected_count", "forw_count", "rev_count",
                        "forw_concat_count", "rev_concat_count", "filtered_count"))
  counted <- counts[counts$true_barcode == barcode, , drop = FALSE]
  expect_true(any(counted$forw_count == 1L & counted$rev_count == 1L))
})

test_that("minimum read length counts cleaned cDNA", {
  work <- tempfile("rrad-cleaned-length-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  anchor <- "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAA"
  anchor2 <- "CAGCACCTGATTACAAGGTGCTG"
  barcode <- "GACCTAGTTCAGGATC"
  cdna <- "GACTCGACTAGCTACGATCGTAGCATGCTAGCTGACTAGCA"
  whitelist <- file.path(work, "barcodes.txt")
  layout <- file.path(work, "layout.csv")
  fastq <- file.path(work, "reads.fastq")
  writeLines(barcode, whitelist)
  writeLines(c(
    "Read Layout,,,,,", "id,seq,expected_length,type,class,whitelist",
    paste0("anchor,", anchor, ",,static,linker,"),
    paste0("barcode,,16,variable,barcode,", whitelist),
    "read,,,variable,read,", paste0("anchor2,", anchor2, ",,static,linker,")
  ), layout)
  sequence <- paste0(anchor, barcode, substr(cdna, 1L, 20L),
                     paste(rep("N", 20L), collapse = ""),
                     substr(cdna, 21L, nchar(cdna)), anchor2)
  port_write_fastq(fastq, c(trim = sequence))
  run <- function(minimum) rad_demux(
    layout, fastq, file.path(work, paste0("min-", minimum)),
    threads = 1L, chunk_size = 1L, min_read_length = minimum
  )
  kept <- run(nchar(cdna))
  removed <- run(nchar(cdna) + 1L)
  expect_equal(kept$stats$records_serialized, 1)
  expect_identical(port_read_fastq(kept$artifacts$fastq)[[2L]], cdna)
  expect_equal(removed$stats$reads_passing_filter, 0)
  expect_length(port_read_fastq(removed$artifacts$fastq), 0L)
})

test_that("joint barcode mask violations cannot fall back to independent calls", {
  work <- tempfile("rrad-joint-pair-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  anchor <- "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAA"
  anchor2 <- "CAGCACCTGATTACAAGGTGCTG"
  axis1 <- c("AACCAGTA", "TGGCAATC")
  axis2 <- c("CGTAGGCT", "GCTCCTGA")
  paths <- file.path(work, c("axis1.txt", "axis2.txt", "mask.csv"))
  writeLines(axis1, paths[[1L]])
  writeLines(axis2, paths[[2L]])
  writeLines(c("1,0", "0,1"), paths[[3L]])
  layout <- file.path(work, "layout.csv")
  fastq <- file.path(work, "reads.fastq")
  writeLines(c(
    "Read Layout,,,,,,", "id,seq,expected_length,type,class,whitelist,flags",
    paste0("anchor,", anchor, ",,static,linker,,"),
    paste0("barcode_1,,8,variable,barcode,", paths[[1L]], ":", paths[[3L]], ",joint_barcode"),
    paste0("barcode_2,,8,variable,barcode,", paths[[2L]], ":", paths[[3L]], ",joint_barcode"),
    "umi,,8,variable,umi,,", paste0("anchor2,", anchor2, ",,static,linker,,"),
    "read,,,variable,read,,"
  ), layout)
  allowed <- paste0(anchor, axis1[[1L]], axis2[[1L]], "CGTATGCA", anchor2,
                    "GACTCGACTAGCTACGATCGTAGCATGCTAGCTGACTAGCA")
  forbidden <- sub(axis2[[1L]], axis2[[2L]], allowed, fixed = TRUE)
  port_write_fastq(fastq, c(allowed_f = allowed, forbidden_f = forbidden,
                          allowed_r = port_reverse_complement(allowed),
                          forbidden_r = port_reverse_complement(forbidden)))
  result <- rad_demux(layout, fastq, file.path(work, "out"), threads = 1L,
                      chunk_size = 4L, min_read_length = 0L,
                      whitelist_mutation = 0L, generated_mutation = 0L)
  output <- port_read_fastq(result$artifacts$fastq)
  headers <- output[seq.int(1L, length(output), 4L)]
  expect_equal(result$stats$records_serialized, 2)
  expect_true(all(grepl("@allowed_", headers, fixed = TRUE)))
  expect_false(any(grepl("forbidden", output, fixed = TRUE)))
})
