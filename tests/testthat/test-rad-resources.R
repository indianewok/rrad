resource_test_setup <- function(contents = c(
                                  three = "three-prime-resource\n",
                                  five = "five-prime-resource\n"
                                )) {
  root <- tempfile("rrad-resource-test-")
  dir.create(root)
  sources <- file.path(root, paste0(names(contents), ".gz"))
  Map(writeBin, lapply(contents, charToRaw), sources)
  hashes <- vapply(
    sources, digest::digest, character(1),
    algo = "sha256", file = TRUE, serialize = FALSE
  )
  sizes <- file.info(sources)$size
  urls <- paste0("https://example.invalid/", basename(sources))
  manifest <- file.path(root, "REMOTE_RESOURCES.tsv")
  table <- data.frame(
    key = c("10x_3v3", "10x_5v3"),
    path = paste0("wl/", basename(sources)),
    asset_set = c("rad-test", "rad-test"),
    url = urls,
    sha256 = unname(hashes),
    bytes = as.character(sizes),
    kits = c("10x_3v3;10x_3v3.1;10x_3HTv3.1", "10x_5v3"),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    table, manifest, sep = "\t", row.names = FALSE, quote = FALSE
  )
  source_by_url <- stats::setNames(sources, urls)
  downloads <- 0L
  downloader <- function(url, destfile, mode, quiet) {
    downloads <<- downloads + 1L
    source <- source_by_url[[url]]
    if (is.null(source) || !file.copy(source, destfile, overwrite = TRUE)) {
      stop("mock download failed")
    }
    0L
  }
  list(
    root = root,
    manifest = manifest,
    cache = file.path(root, "cache"),
    sources = sources,
    contents = contents,
    downloader = downloader,
    downloads = function() downloads
  )
}

with_resource_test_options <- function(setup, code, offline = FALSE) {
  old_options <- options(
    rrad.resource_manifest = setup$manifest,
    rrad.cache_dir = setup$cache,
    rrad.download.file = setup$downloader,
    rrad.offline = offline
  )
  old_env <- Sys.getenv("RRAD_OFFLINE", unset = NA_character_)
  old_overlay <- Sys.getenv("RRAD_RESOURCE_CACHE", unset = NA_character_)
  old_cache <- Sys.getenv("RRAD_CACHE_DIR", unset = NA_character_)
  on.exit({
    options(old_options)
    if (is.na(old_env)) Sys.unsetenv("RRAD_OFFLINE") else {
      Sys.setenv(RRAD_OFFLINE = old_env)
    }
    if (is.na(old_overlay)) Sys.unsetenv("RRAD_RESOURCE_CACHE") else {
      Sys.setenv(RRAD_RESOURCE_CACHE = old_overlay)
    }
    if (is.na(old_cache)) Sys.unsetenv("RRAD_CACHE_DIR") else {
      Sys.setenv(RRAD_CACHE_DIR = old_cache)
    }
  }, add = TRUE)
  force(code)
}

test_that("resource keys, aliases, basenames, and paths resolve identically", {
  setup <- resource_test_setup()
  on.exit(unlink(setup$root, recursive = TRUE), add = TRUE)

  with_resource_test_options(setup, {
    canonical <- rrad:::rad_resource_path("10x_3v3", quiet = TRUE)
    expect_identical(
      rrad:::rad_resource_path("10x_3v3.1", quiet = TRUE), canonical
    )
    expect_identical(
      rrad:::rad_resource_path(basename(setup$sources[[1L]]), quiet = TRUE),
      canonical
    )
    expect_identical(
      rrad:::rad_resource_path(
        paste0("wl/", basename(setup$sources[[1L]])), quiet = TRUE
      ),
      canonical
    )
    expect_identical(setup$downloads(), 1L)
    expect_match(
      canonical,
      "/rad-resources/rad-test/wl/three\\.gz$"
    )
  })
})

test_that("cache downloads are pinned and reused without network access", {
  setup <- resource_test_setup()
  on.exit(unlink(setup$root, recursive = TRUE), add = TRUE)

  with_resource_test_options(setup, {
    paths <- rrad:::rad_download_whitelists(
      c("10x_3v3", "10x_5v3"), quiet = TRUE
    )
    expect_named(paths, c("10x_3v3", "10x_5v3"))
    expect_true(all(file.exists(paths)))
    expect_identical(setup$downloads(), 2L)

    options(rrad.offline = TRUE)
    expect_identical(
      rrad:::rad_resource_path("10x_3v3", quiet = TRUE),
      paths[["10x_3v3"]]
    )
    expect_identical(setup$downloads(), 2L)
  })
})

test_that("the user-facing whitelist downloader delegates to the cache", {
  setup <- resource_test_setup()
  on.exit(unlink(setup$root, recursive = TRUE), add = TRUE)

  with_resource_test_options(setup, {
    expect_message(
      paths <- rad_download_whitelists("10x_3v3"),
      paste0(
        "Downloading the RAD barcode whitelist for '10x_3v3'.*",
        "reuse it for future runs"
      )
    )
    expect_named(paths, "10x_3v3")
    expect_true(file.exists(paths[[1L]]))
    expect_identical(setup$downloads(), 1L)

    expect_silent(
      reused <- rad_download_whitelists("10x_3v3")
    )
    expect_identical(unname(reused), unname(paths))
    expect_identical(setup$downloads(), 1L)
  })
})

test_that("missing resources can be queried without downloading", {
  setup <- resource_test_setup()
  on.exit(unlink(setup$root, recursive = TRUE), add = TRUE)

  with_resource_test_options(setup, {
    expect_identical(
      rrad:::rad_resource_path("10x_3v3", download = FALSE),
      NA_character_
    )
    expect_identical(setup$downloads(), 0L)
    options(rrad.offline = TRUE)
    expect_error(
      rrad:::rad_resource_path("10x_3v3", quiet = TRUE),
      "offline mode is enabled"
    )
    expect_identical(setup$downloads(), 0L)
  })
})

test_that("failed verification never publishes a partial resource", {
  setup <- resource_test_setup()
  on.exit(unlink(setup$root, recursive = TRUE), add = TRUE)
  source_size <- file.info(setup$sources[[1L]])$size
  writeBin(rep(as.raw(0L), source_size), setup$sources[[1L]])

  with_resource_test_options(setup, {
    expect_error(
      rrad:::rad_resource_path("10x_3v3", quiet = TRUE),
      "SHA-256 mismatch"
    )
    expected <- file.path(
      setup$cache, "rad-resources", "rad-test", "wl", "three.gz"
    )
    expect_false(file.exists(expected))
    expect_length(
      list.files(dirname(expected), pattern = "[.]part-", all.files = TRUE),
      0L
    )
  })
})

test_that("resource publication never replaces a cache directory", {
  setup <- resource_test_setup()
  on.exit(unlink(setup$root, recursive = TRUE), add = TRUE)
  target <- file.path(
    setup$cache, "rad-resources", "rad-test", "wl", "three.gz"
  )
  dir.create(target, recursive = TRUE)
  sentinel <- file.path(target, "keep")
  writeLines("user data", sentinel)

  with_resource_test_options(setup, {
    expect_error(
      rrad:::rad_resource_path("10x_3v3", quiet = TRUE),
      "cache target is a directory"
    )
    expect_identical(readLines(sentinel), "user data")
  })
})

test_that("a failed force refresh preserves a valid cached resource", {
  setup <- resource_test_setup()
  on.exit(unlink(setup$root, recursive = TRUE), add = TRUE)

  with_resource_test_options(setup, {
    cached <- rrad:::rad_resource_path("10x_3v3", quiet = TRUE)
    original <- readBin(cached, "raw", n = file.info(cached)$size)
    options(rrad.download.file = function(...) stop("network unavailable"))
    expect_error(
      rrad:::rad_download_whitelists(
        "10x_3v3", force = TRUE, quiet = TRUE
      ),
      "network unavailable"
    )
    expect_identical(
      readBin(cached, "raw", n = file.info(cached)$size),
      original
    )
    expect_identical(
      digest::digest(file = cached, algo = "sha256", serialize = FALSE),
      digest::digest(original, algo = "sha256", serialize = FALSE)
    )
  })
})

test_that("cache removal requires an explicit confirmation", {
  setup <- resource_test_setup()
  on.exit(unlink(setup$root, recursive = TRUE), add = TRUE)

  with_resource_test_options(setup, {
    cached <- rrad:::rad_resource_path("10x_3v3", quiet = TRUE)
    expect_error(
      rrad:::rad_clear_resource_cache(),
      "refusing to remove"
    )
    expect_true(file.exists(cached))
    expect_true(rrad:::rad_clear_resource_cache(confirm = TRUE))
    expect_false(file.exists(cached))
    expect_false(rrad:::rad_clear_resource_cache(confirm = TRUE))
  })
})

test_that("native missing markers download, attach, and retry once", {
  setup <- resource_test_setup()
  on.exit(unlink(setup$root, recursive = TRUE), add = TRUE)

  with_resource_test_options(setup, {
    Sys.setenv(RRAD_RESOURCE_CACHE = "previous-value")
    tampered_sibling <- file.path(
      setup$cache, "rad-resources", "rad-test", "wl", "five.gz"
    )
    dir.create(dirname(tampered_sibling), recursive = TRUE)
    writeLines("tampered sibling", tampered_sibling)
    native_options <- list(kit = "implicit-from-layout", threads = 3L)
    observed <- list()
    native_call <- function(options) {
      overlay <- Sys.getenv("RRAD_RESOURCE_CACHE", unset = "")
      observed[[length(observed) + 1L]] <<- list(
        options = options,
        overlay = overlay,
        sibling_visible = file.exists(file.path(overlay, "wl", "five.gz"))
      )
      whitelist <- file.path(overlay, "wl", "three.gz")
      if (!file.exists(whitelist)) {
        stop(
          paste0(
            "embedded RAD demultiplexing failed: ",
            "RRAD_RESOURCE_MISSING:wl/three.gz"
          ),
          call. = FALSE
        )
      }
      list(status = "success", whitelist = whitelist)
    }

    result <- rrad:::.rad_call_with_resource_retry(
      function() native_call(native_options), quiet = TRUE
    )
    expect_identical(result$status, "success")
    expect_false(file.exists(result$whitelist))
    expect_true(file.exists(rrad:::rad_resource_path(
      "10x_3v3", download = FALSE
    )))
    expect_length(observed, 2L)
    expect_identical(observed[[1L]]$options, native_options)
    expect_identical(observed[[2L]]$options, native_options)
    expect_identical(observed[[1L]]$overlay, "")
    expect_match(
      observed[[2L]]$overlay,
      "/rrad-native-overlay-"
    )
    expect_false(observed[[2L]]$sibling_visible)
    expect_true(file.exists(tampered_sibling))
    expect_false(dir.exists(observed[[2L]]$overlay))
    expect_identical(setup$downloads(), 1L)
    expect_identical(Sys.getenv("RRAD_RESOURCE_CACHE"), "previous-value")
  })
})

test_that("native overlays own verified bytes outside the shared cache", {
  setup <- resource_test_setup()
  on.exit(unlink(setup$root, recursive = TRUE), add = TRUE)

  with_resource_test_options(setup, {
    observations <- list()
    native_call <- function() {
      overlay <- Sys.getenv("RRAD_RESOURCE_CACHE", unset = "")
      whitelist <- file.path(overlay, "wl", "three.gz")
      if (!file.exists(whitelist)) {
        stop("RRAD_RESOURCE_MISSING:wl/three.gz", call. = FALSE)
      }

      cached <- rrad:::.rad_resource_target(
        rrad:::.rad_resource_match("10x_3v3")
      )
      original <- readBin(whitelist, "raw", n = file.info(whitelist)$size)
      writeBin(rep(as.raw(0x5a), length(original)), cached)
      after_cache_mutation <- readBin(
        whitelist, "raw", n = file.info(whitelist)$size
      )
      cache_cleared <- rrad:::rad_clear_resource_cache(confirm = TRUE)
      after_cache_clear <- if (file.exists(whitelist)) {
        readBin(whitelist, "raw", n = file.info(whitelist)$size)
      } else {
        raw()
      }
      observations <<- list(
        overlay = overlay,
        overlay_mode = as.character(file.info(overlay)$mode),
        resource_dir_mode = as.character(file.info(dirname(whitelist))$mode),
        whitelist_is_symlink = nzchar(Sys.readlink(whitelist)),
        mutation_preserved = identical(after_cache_mutation, original),
        clear_preserved = identical(after_cache_clear, original),
        cache_cleared = cache_cleared
      )
      list(status = "success", whitelist = whitelist)
    }

    result <- rrad:::.rad_call_with_resource_retry(
      native_call, quiet = TRUE
    )
    expect_identical(result$status, "success")
    expect_false(startsWith(
      observations$overlay,
      paste0(normalizePath(
        file.path(setup$cache, "rad-resources"),
        winslash = "/", mustWork = FALSE
      ), "/")
    ))
    expect_identical(observations$overlay_mode, "700")
    expect_identical(observations$resource_dir_mode, "700")
    expect_false(observations$whitelist_is_symlink)
    expect_true(observations$mutation_preserved)
    expect_true(observations$cache_cleared)
    expect_true(observations$clear_preserved)
    expect_false(dir.exists(observations$overlay))
  })
})

test_that("native retry is bounded and only catches exact missing markers", {
  setup <- resource_test_setup()
  on.exit(unlink(setup$root, recursive = TRUE), add = TRUE)

  with_resource_test_options(setup, {
    calls <- 0L
    always_missing <- function() {
      calls <<- calls + 1L
      stop("RRAD_RESOURCE_MISSING:wl/three.gz", call. = FALSE)
    }
    expect_error(
      rrad:::.rad_call_with_resource_retry(always_missing, quiet = TRUE),
      "still reports missing resource"
    )
    expect_identical(calls, 2L)

    calls <- 0L
    near_match <- function() {
      calls <<- calls + 1L
      stop("prefix RRAD_RESOURCE_MISSING:wl/three.gz", call. = FALSE)
    }
    expect_error(
      rrad:::.rad_call_with_resource_retry(near_match, quiet = TRUE),
      "prefix RRAD_RESOURCE_MISSING"
    )
    expect_identical(calls, 1L)

    multiple <- function() {
      stop(
        paste0(
          "RRAD_RESOURCE_MISSING:wl/three.gz: ",
          "RRAD_RESOURCE_MISSING:wl/three.gz"
        ),
        call. = FALSE
      )
    }
    expect_error(
      rrad:::.rad_call_with_resource_retry(multiple, quiet = TRUE),
      "RRAD_RESOURCE_MISSING"
    )
  })
})

test_that("RRAD_CACHE_DIR supplies the cache base", {
  setup <- resource_test_setup()
  on.exit(unlink(setup$root, recursive = TRUE), add = TRUE)

  with_resource_test_options(setup, {
    options(rrad.cache_dir = NULL)
    environment_cache <- file.path(setup$root, "environment-cache")
    Sys.setenv(RRAD_CACHE_DIR = environment_cache)
    on.exit(Sys.unsetenv("RRAD_CACHE_DIR"), add = TRUE)
    path <- rrad:::rad_resource_path("10x_3v3", quiet = TRUE)
    expect_match(path, "/environment-cache/rad-resources/rad-test/")
  })
})

test_that("RRAD_OFFLINE prohibits downloads when no option overrides it", {
  setup <- resource_test_setup()
  on.exit(unlink(setup$root, recursive = TRUE), add = TRUE)

  with_resource_test_options(setup, {
    options(rrad.offline = NULL)
    Sys.setenv(RRAD_OFFLINE = "yes")
    expect_error(
      rrad:::rad_resource_path("10x_3v3", quiet = TRUE),
      "offline mode is enabled"
    )
    expect_identical(setup$downloads(), 0L)

    Sys.setenv(RRAD_OFFLINE = "sometimes")
    expect_error(
      rrad:::rad_resource_path("10x_3v3", quiet = TRUE),
      "RRAD_OFFLINE must be"
    )
  })
})

test_that("native scan downloads an implicit resource and reuses it offline", {
  core <- tryCatch(rad_core_info(), error = function(error) error)
  if (inherits(core, "error")) {
    skip("requires an installed package resource tree")
  }

  work <- tempfile("rrad-native-resource-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  source <- file.path(work, "source.gz")
  connection <- gzfile(source, open = "wt")
  writeLines("AAAAAAAAAAAAAAAA", connection)
  close(connection)

  manifest <- file.path(work, "REMOTE_RESOURCES.tsv")
  url <- paste0(
    "https://example.invalid/",
    "3M-february-2018-3v3.txt_bitlist.csv.gz"
  )
  table <- data.frame(
    key = "10x_3v3",
    path = "wl/3M-february-2018-3v3.txt_bitlist.csv.gz",
    asset_set = "native-probe",
    url = url,
    sha256 = digest::digest(
      file = source, algo = "sha256", serialize = FALSE
    ),
    bytes = as.character(file.info(source)$size),
    kits = "10x_3v3;10x_3v3.1;10x_3HTv3.1",
    stringsAsFactors = FALSE
  )
  utils::write.table(
    table, manifest, sep = "\t", row.names = FALSE, quote = FALSE
  )

  sequence <- paste0("ACGT", "AAAAAAAAAAAAAAAA")
  input <- file.path(work, "reads.fastq")
  writeLines(c(
    "@r1", sequence, "+", strrep("I", nchar(sequence)),
    "@r2", sequence, "+", strrep("I", nchar(sequence))
  ), input)

  downloads <- 0L
  downloader <- function(url, destfile, mode, quiet) {
    downloads <<- downloads + 1L
    if (!file.copy(source, destfile, overwrite = TRUE)) {
      stop("mock native-resource download failed")
    }
    0L
  }
  setup <- list(
    manifest = manifest,
    cache = file.path(work, "cache"),
    downloader = downloader
  )

  with_resource_test_options(setup, {
    first <- rad_scan_wl(
      input = input,
      output_prefix = file.path(work, "cells"),
      adapter_seq = "ACGT",
      barcode_length = 16L,
      max_error = 0,
      whitelist = "10x_3v3",
      selection = "above_floor",
      chunk_size = 1L
    )
    expect_identical(first$status, "success")
    expect_identical(readLines(first$artifacts$whitelist),
                     "AAAAAAAAAAAAAAAA")
    expect_identical(downloads, 1L)
    expect_true(file.exists(rad_resource_path(
      "10x_3v3", download = FALSE
    )))

    options(rrad.offline = TRUE)
    second <- rad_scan_wl(
      input = input,
      output_prefix = file.path(work, "cells-offline"),
      adapter_seq = "ACGT",
      barcode_length = 16L,
      max_error = 0,
      whitelist = "10x_3v3",
      selection = "above_floor",
      chunk_size = 1L
    )
    expect_identical(second$status, "success")
    expect_identical(downloads, 1L)
  })
})

test_that("manifest safety checks reject traversal and unpinned input", {
  setup <- resource_test_setup()
  on.exit(unlink(setup$root, recursive = TRUE), add = TRUE)

  manifest <- utils::read.delim(
    setup$manifest, sep = "\t", colClasses = "character",
    check.names = FALSE
  )
  manifest$path[[1L]] <- "../outside.gz"
  manifest$url[[2L]] <- "http://example.invalid/five.gz"
  utils::write.table(
    manifest, setup$manifest, sep = "\t", row.names = FALSE, quote = FALSE
  )
  with_resource_test_options(setup, {
    expect_error(rrad:::.rad_resource_manifest(), "unsafe RAD resource path")
  })
})
