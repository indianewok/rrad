test_that("native demultiplexing is safe when HOME is unset", {
  rscript <- file.path(R.home("bin"), "Rscript")
  library_path <- normalizePath(.libPaths()[[1L]], mustWork = TRUE)
  fastq <- normalizePath(
    test_path("fixtures", "rad-core", "nanopore-smoke.fastq"),
    mustWork = TRUE
  )
  out_dir <- tempfile("rrad-no-home-")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  quote_r <- function(value) encodeString(value, quote = "\"")
  code <- paste0(
    "Sys.unsetenv('HOME'); ",
    "library(rrad, lib.loc = ", quote_r(library_path), "); ",
    "result <- rad_demux(",
    "layout = 'nanopore_rapid_bc', ",
    "fastq = ", quote_r(fastq), ", ",
    "out_dir = ", quote_r(out_dir), ", ",
    "output = 'smoke', threads = 1L, chunk_size = 1L, ",
    "max_reads = 1L, min_read_length = 0L); ",
    "stopifnot(inherits(result, 'rad_demux_result'), ",
    "identical(result$status, 'success'))"
  )

  output <- suppressWarnings(system2(
    rscript,
    c("--vanilla", "-e", shQuote(code)),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status", exact = TRUE)
  if (is.null(status)) status <- 0L

  expect_identical(
    as.integer(status),
    0L,
    info = paste(output, collapse = "\n")
  )
})

test_that("deterministic read operations preserve the user's RNG state", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  work <- tempfile("rrad-rng-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  exercise_reads <- function() {
    reads <- rad_pack_reads(c("r1", "r2"), c("ACGT", "TTAA"),
                            c("IIII", "HHHH"))
    rad_validate(reads, deep = TRUE)
    rad_materialize(reads)
    rad_qmean(reads$qual)
    rad_write_fastq(reads, file.path(work, "reads.fastq"))
    rad_read_fastq(file.path(work, "reads.fastq"))
    rad_core_info()
    invisible(NULL)
  }

  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  exercise_reads()
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

  set.seed(123)
  before <- get(".Random.seed", envir = .GlobalEnv)
  exercise_reads()
  expect_identical(get(".Random.seed", envir = .GlobalEnv), before)
})
