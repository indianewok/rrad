rad_canonical_tree_sha256 <- function(root) {
  paths <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  paths <- sort(paths, method = "radix")
  hashes <- vapply(file.path(root, paths), digest::digest, character(1),
                   algo = "sha256", file = TRUE, serialize = FALSE,
                   USE.NAMES = FALSE)
  manifest <- paste0(hashes, "  ", enc2utf8(paths), "\n", collapse = "")
  digest::digest(manifest, algo = "sha256", serialize = FALSE)
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

  source_digest <- rad_canonical_tree_sha256(include_dir)
  resource_digest <- rad_canonical_tree_sha256(resource_dir)
  vendor_source <- rad_json_sha256(vendor_json, "embedded_source_digest")
  vendor_resource <- rad_json_sha256(vendor_json, "resource_tree_sha256")
  bridge_source <- rad_bridge_sha256(bridge, "kEmbeddedSourceDigest")
  bridge_resource <- rad_bridge_sha256(bridge, "kResourceDigest")

  expect_identical(vendor_source, source_digest)
  expect_identical(bridge_source, source_digest)
  expect_identical(vendor_resource, resource_digest)
  expect_identical(bridge_resource, resource_digest)
})
