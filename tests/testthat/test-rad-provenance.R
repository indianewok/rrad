rad_provenance_command <- function(command, args = character(),
                                   env = character()) {
  output <- suppressWarnings(system2(
    command, args, stdout = TRUE, stderr = TRUE, env = env
  ))
  status <- attr(output, "status", exact = TRUE)
  if (!is.null(status) && status != 0L) {
    stop(paste(output, collapse = "\n"), call. = FALSE)
  }
  output
}

rad_canonical_tree_sha256 <- function(root, commands) {
  old_dir <- setwd(root)
  on.exit(setwd(old_dir), add = TRUE)

  paths <- rad_provenance_command(
    commands[["find"]], c(".", "-type", "f", "-print"), env = "LC_ALL=C"
  )
  if (length(paths) == 0L || any(!startsWith(paths, "./"))) {
    stop("canonical digest input must be a non-empty relative file tree",
         call. = FALSE)
  }
  paths <- substring(paths, 3L)

  path_manifest <- tempfile("rad-provenance-paths-")
  hash_manifest <- tempfile("rad-provenance-hashes-")
  on.exit(unlink(c(path_manifest, hash_manifest)), add = TRUE)
  writeLines(enc2utf8(paths), path_manifest, useBytes = TRUE)
  paths <- rad_provenance_command(
    commands[["sort"]], shQuote(path_manifest), env = "LC_ALL=C"
  )

  hash_file <- function(path) {
    hash_command <- commands[["hash"]]
    args <- if (identical(basename(hash_command), "shasum")) {
      c("-a", "256", shQuote(path))
    } else {
      shQuote(path)
    }
    output <- rad_provenance_command(hash_command, args)
    if (length(output) != 1L ||
        !grepl("^[0-9a-f]{64}[ *]", output, perl = TRUE)) {
      stop("unexpected SHA-256 output while building canonical manifest",
           call. = FALSE)
    }
    substring(output, 1L, 64L)
  }

  hashes <- vapply(paths, function(path) {
    hash_file(path)
  }, character(1L), USE.NAMES = FALSE)

  manifest <- paste0(hashes, "  ", enc2utf8(paths), "\n", collapse = "")
  writeBin(charToRaw(manifest), hash_manifest)
  hash_file(hash_manifest)
}

rad_json_sha256 <- function(path, field) {
  pattern <- sprintf(
    '^[[:space:]]*"%s"[[:space:]]*:[[:space:]]*"([0-9a-f]{64})"',
    field
  )
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  matches <- regmatches(lines, regexec(pattern, lines, perl = TRUE))
  values <- vapply(matches, function(match) {
    if (length(match) == 2L) match[[2L]] else ""
  }, character(1L))
  values <- values[nzchar(values)]
  if (length(values) != 1L) {
    stop(sprintf("could not read one `%s` value from RAD_VENDOR.json", field),
         call. = FALSE)
  }
  values[[1L]]
}

rad_bridge_sha256 <- function(path, constant) {
  source <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"),
                  collapse = "\n")
  pattern <- sprintf(
    '\\b%s[[:space:]]*=[[:space:]]*"([0-9a-f]{64})"[[:space:]]*;',
    constant
  )
  match <- regmatches(source, regexec(pattern, source, perl = TRUE))[[1L]]
  if (length(match) != 2L) {
    stop(sprintf("could not read `%s` from rad_core_bridge.cpp", constant),
         call. = FALSE)
  }
  match[[2L]]
}

test_that("RAD source and resource provenance digests stay synchronized", {
  root <- normalizePath(test_path("..", ".."), mustWork = FALSE)
  include_dir <- file.path(root, "src", "vendor", "rad_core", "include")
  resource_dir <- file.path(root, "inst", "rad", "resources")
  vendor_json <- file.path(root, "src", "vendor", "rad_core",
                           "RAD_VENDOR.json")
  bridge <- file.path(root, "src", "rad_core_bridge.cpp")

  source_inputs <- c(include_dir, resource_dir, vendor_json, bridge)
  if (!all(file.exists(source_inputs))) {
    skip("RAD provenance source trees are unavailable in this test context")
  }

  commands <- Sys.which(c("find", "sort", "shasum", "sha256sum"))
  if (!nzchar(commands[["find"]]) || !nzchar(commands[["sort"]]) ||
      (!nzchar(commands[["shasum"]]) &&
       !nzchar(commands[["sha256sum"]]))) {
    skip("RAD provenance check requires find, sort, and a SHA-256 utility")
  }
  commands[["hash"]] <- if (nzchar(commands[["shasum"]])) {
    commands[["shasum"]]
  } else {
    commands[["sha256sum"]]
  }

  source_digest <- rad_canonical_tree_sha256(include_dir, commands)
  resource_digest <- rad_canonical_tree_sha256(resource_dir, commands)
  vendor_source <- rad_json_sha256(vendor_json, "embedded_source_digest")
  vendor_resource <- rad_json_sha256(vendor_json, "resource_tree_sha256")
  bridge_source <- rad_bridge_sha256(bridge, "kEmbeddedSourceDigest")
  bridge_resource <- rad_bridge_sha256(bridge, "kResourceDigest")

  expect_identical(vendor_source, source_digest)
  expect_identical(bridge_source, source_digest)
  expect_identical(vendor_resource, resource_digest)
  expect_identical(bridge_resource, resource_digest)
})
