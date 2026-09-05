# Downloadable RAD resources -------------------------------------------------

.rad_resource_manifest_file <- function() {
  override <- getOption("rrad.resource_manifest", NULL)
  if (!is.null(override)) {
    if (!is.character(override) || length(override) != 1L ||
        is.na(override) || !nzchar(override)) {
      stop("option 'rrad.resource_manifest' must be one file path",
           call. = FALSE)
    }
    path <- path.expand(override)
  } else {
    path <- system.file(
      "rad", "resources", "REMOTE_RESOURCES.tsv",
      package = "rrad"
    )
  }

  if (!nzchar(path) || !file.exists(path) || dir.exists(path)) {
    stop("the installed RAD remote-resource manifest is unavailable",
         call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.rad_resource_manifest <- function() {
  path <- .rad_resource_manifest_file()
  manifest <- tryCatch(
    utils::read.delim(
      path,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      colClasses = "character",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      na.strings = character()
    ),
    error = function(error) {
      stop(
        "could not parse the RAD remote-resource manifest: ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )

  required <- c("key", "path", "asset_set", "url", "sha256", "bytes", "kits")
  missing <- setdiff(required, names(manifest))
  if (length(missing)) {
    stop(
      "the RAD remote-resource manifest is missing fields: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  manifest <- manifest[, required, drop = FALSE]
  if (!nrow(manifest)) {
    stop("the RAD remote-resource manifest is empty", call. = FALSE)
  }

  scalar_fields <- c("key", "path", "asset_set", "url", "sha256", "bytes")
  for (field in scalar_fields) {
    manifest[[field]] <- trimws(manifest[[field]])
    invalid <- is.na(manifest[[field]]) | !nzchar(manifest[[field]])
    if (any(invalid)) {
      stop(
        "the RAD remote-resource manifest has an empty '", field,
        "' value on row ", which(invalid)[[1L]],
        call. = FALSE
      )
    }
  }
  manifest$kits[is.na(manifest$kits)] <- ""
  manifest$kits <- trimws(manifest$kits)

  invalid_key <- !grepl("^[[:alnum:]][[:alnum:]_.-]*$", manifest$key)
  if (any(invalid_key)) {
    stop("invalid RAD resource key on row ", which(invalid_key)[[1L]],
         call. = FALSE)
  }
  invalid_asset <- !grepl(
    "^[[:alnum:]][[:alnum:]_.-]*$", manifest$asset_set
  )
  if (any(invalid_asset)) {
    stop("invalid RAD resource asset_set on row ",
         which(invalid_asset)[[1L]], call. = FALSE)
  }
  path_parts <- strsplit(manifest$path, "/", fixed = TRUE)
  invalid_path <- grepl("\\\\", manifest$path) |
    grepl("^/|/$|//", manifest$path) |
    vapply(
      path_parts,
      function(parts) any(!nzchar(parts) | parts %in% c(".", "..")),
      logical(1)
    )
  if (any(invalid_path)) {
    stop("unsafe RAD resource path on row ", which(invalid_path)[[1L]],
         call. = FALSE)
  }
  invalid_url <- !grepl("^https://[^[:space:]]+$", manifest$url)
  if (any(invalid_url)) {
    stop("RAD resource URLs must use HTTPS (row ",
         which(invalid_url)[[1L]], ")", call. = FALSE)
  }
  manifest$sha256 <- tolower(manifest$sha256)
  invalid_sha <- !grepl("^[[:xdigit:]]{64}$", manifest$sha256)
  if (any(invalid_sha)) {
    stop("invalid RAD resource SHA-256 on row ",
         which(invalid_sha)[[1L]], call. = FALSE)
  }
  invalid_bytes <- !grepl("^[1-9][0-9]*$", manifest$bytes)
  manifest$bytes_value <- suppressWarnings(as.numeric(manifest$bytes))
  invalid_bytes <- invalid_bytes | !is.finite(manifest$bytes_value) |
    manifest$bytes_value > 2^53
  if (any(invalid_bytes)) {
    stop("invalid RAD resource byte count on row ",
         which(invalid_bytes)[[1L]], call. = FALSE)
  }

  if (anyDuplicated(tolower(manifest$key))) {
    stop("RAD resource keys must be unique", call. = FALSE)
  }
  target_id <- tolower(paste(manifest$asset_set, manifest$path, sep = "/"))
  if (anyDuplicated(target_id)) {
    stop("RAD remote resources must have unique cache targets",
         call. = FALSE)
  }

  aliases <- lapply(
    manifest$kits,
    function(value) {
      if (!nzchar(value)) return(character())
      values <- trimws(strsplit(value, ";", fixed = TRUE)[[1L]])
      values[nzchar(values)]
    }
  )
  alias_owner <- unlist(
    Map(function(values, row) stats::setNames(rep(row, length(values)), values),
        aliases, seq_len(nrow(manifest))),
    use.names = TRUE
  )
  if (length(alias_owner) && anyDuplicated(tolower(names(alias_owner)))) {
    stop("RAD resource kit aliases must be unique", call. = FALSE)
  }
  attr(manifest, "aliases") <- aliases
  manifest
}

.rad_resource_match <- function(key, manifest = .rad_resource_manifest()) {
  if (!is.character(key) || length(key) != 1L || is.na(key) ||
      !nzchar(trimws(key))) {
    stop("key must be one non-empty character value", call. = FALSE)
  }
  key <- trimws(key)
  aliases <- attr(manifest, "aliases")
  if (is.null(aliases)) {
    aliases <- lapply(
      manifest$kits,
      function(value) {
        if (!nzchar(value)) character() else {
          trimws(strsplit(value, ";", fixed = TRUE)[[1L]])
        }
      }
    )
  }
  needle <- tolower(key)
  matched <- which(
    tolower(manifest$key) == needle |
      tolower(manifest$path) == needle |
      tolower(basename(manifest$path)) == needle |
      vapply(aliases, function(values) needle %in% tolower(values), logical(1))
  )
  if (!length(matched)) {
    stop(
      "unknown RAD whitelist or kit '", key, "'; available names: ",
      paste(manifest$key, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(matched) > 1L) {
    stop("ambiguous RAD whitelist or kit '", key, "'", call. = FALSE)
  }
  manifest[matched, , drop = FALSE]
}

.rad_resource_cache_base <- function(cache_dir = NULL) {
  if (is.null(cache_dir)) {
    cache_dir <- getOption("rrad.cache_dir", NULL)
  }
  if (is.null(cache_dir)) {
    environment_cache <- Sys.getenv("RRAD_CACHE_DIR", unset = "")
    if (nzchar(environment_cache)) cache_dir <- environment_cache
  }
  if (is.null(cache_dir)) {
    cache_dir <- tools::R_user_dir("rrad", which = "cache")
  }
  if (!is.character(cache_dir) || length(cache_dir) != 1L ||
      is.na(cache_dir) || !nzchar(cache_dir)) {
    stop("cache_dir must be one non-empty directory path", call. = FALSE)
  }
  cache_dir <- normalizePath(
    path.expand(cache_dir), winslash = "/", mustWork = FALSE
  )
  if (cache_dir %in% c("/", "//") || grepl("^[A-Za-z]:/$", cache_dir)) {
    stop("cache_dir must not be a filesystem root", call. = FALSE)
  }
  if (file.exists(cache_dir) && !dir.exists(cache_dir)) {
    stop("cache_dir identifies a file: ", cache_dir, call. = FALSE)
  }
  cache_dir
}

.rad_resource_cache_root <- function(cache_dir = NULL) {
  file.path(.rad_resource_cache_base(cache_dir), "rad-resources")
}

.rad_resource_target <- function(row, cache_dir = NULL) {
  do.call(
    file.path,
    as.list(c(
      .rad_resource_cache_root(cache_dir),
      row$asset_set[[1L]],
      strsplit(row$path[[1L]], "/", fixed = TRUE)[[1L]]
    ))
  )
}

.rad_resource_bundled_path <- function(row) {
  parts <- strsplit(row$path[[1L]], "/", fixed = TRUE)[[1L]]
  path <- do.call(
    system.file,
    c(as.list(c("rad", "resources", parts)), list(package = "rrad"))
  )
  if (nzchar(path) && file.exists(path) && !dir.exists(path)) path else ""
}

.rad_resource_integrity <- function(path, row) {
  if (!file.exists(path) || dir.exists(path)) {
    return(list(ok = FALSE, reason = "file is missing"))
  }
  size <- file.info(path)$size
  expected_size <- row$bytes_value[[1L]]
  if (is.na(size) || size != expected_size) {
    return(list(
      ok = FALSE,
      reason = paste0(
        "expected ", format(expected_size, scientific = FALSE),
        " bytes but found ", format(size, scientific = FALSE)
      )
    ))
  }
  actual <- tryCatch(
    digest::digest(file = path, algo = "sha256", serialize = FALSE),
    error = function(error) error
  )
  if (inherits(actual, "error")) {
    return(list(
      ok = FALSE,
      reason = paste0("could not calculate SHA-256: ", conditionMessage(actual))
    ))
  }
  expected <- row$sha256[[1L]]
  if (!identical(tolower(actual), expected)) {
    return(list(
      ok = FALSE,
      reason = paste0(
        "SHA-256 mismatch (expected ", expected,
        " but found ", tolower(actual), ")"
      )
    ))
  }
  list(ok = TRUE, reason = NULL)
}

.rad_resource_offline <- function() {
  option <- getOption("rrad.offline", NULL)
  if (!is.null(option)) {
    if (!is.logical(option) || length(option) != 1L || is.na(option)) {
      stop("option 'rrad.offline' must be TRUE or FALSE", call. = FALSE)
    }
    return(isTRUE(option))
  }

  value <- trimws(tolower(Sys.getenv("RRAD_OFFLINE", unset = "")))
  if (!nzchar(value) || value %in% c("0", "false", "no", "off")) {
    return(FALSE)
  }
  if (value %in% c("1", "true", "yes", "on")) {
    return(TRUE)
  }
  stop(
    "RRAD_OFFLINE must be one of true/false, yes/no, on/off, or 1/0",
    call. = FALSE
  )
}

.rad_resource_download <- function(url, destfile, quiet) {
  downloader <- getOption("rrad.download.file", NULL)
  if (is.null(downloader)) downloader <- utils::download.file
  if (!is.function(downloader)) {
    stop("option 'rrad.download.file' must be a function", call. = FALSE)
  }
  status <- downloader(
    url = url,
    destfile = destfile,
    mode = "wb",
    quiet = quiet
  )
  if (!is.null(status) && length(status) && !is.na(status[[1L]]) &&
      status[[1L]] != 0L) {
    stop("download failed with status ", status[[1L]], call. = FALSE)
  }
  invisible(destfile)
}

.rad_resource_publish <- function(temp, target, row) {
  # A concurrent process may have published the same immutable object while
  # this process was downloading. The pinned hash makes either copy valid.
  current <- .rad_resource_integrity(target, row)
  if (isTRUE(current$ok)) {
    unlink(temp)
    return(normalizePath(target, winslash = "/", mustWork = TRUE))
  }

  if (dir.exists(target)) {
    stop(
      "the RAD resource cache target is a directory and will not be replaced: ",
      target,
      call. = FALSE
    )
  }

  backup <- ""
  if (file.exists(target)) {
    backup <- tempfile(
      pattern = paste0(".", basename(target), ".invalid-"),
      tmpdir = dirname(target)
    )
    if (!file.rename(target, backup)) {
      stop("could not replace the invalid cached RAD resource: ", target,
           call. = FALSE)
    }
  }
  published <- file.rename(temp, target)
  if (!published) {
    if (nzchar(backup) && file.exists(backup) && !file.exists(target)) {
      file.rename(backup, target)
    }
    stop("could not save the downloaded RAD whitelist: ", target,
         call. = FALSE)
  }
  if (nzchar(backup) && file.exists(backup)) unlink(backup)
  normalizePath(target, winslash = "/", mustWork = TRUE)
}

.rad_resource_fetch <- function(row, cache_dir = NULL, force = FALSE,
                                quiet = FALSE) {
  target <- .rad_resource_target(row, cache_dir)
  current <- .rad_resource_integrity(target, row)
  if (isTRUE(current$ok) && !force) {
    return(normalizePath(target, winslash = "/", mustWork = TRUE))
  }
  if (.rad_resource_offline()) {
    stop(
      "The RAD barcode whitelist for '", row$key[[1L]],
      "' is not available because offline mode is enabled. ",
      "Local copy status: ", current$reason, ". Expected a verified file at ",
      target, ". ",
      "Run rrad::rad_download_whitelists(\"", row$key[[1L]],
      "\") while online, or supply an explicit whitelist path.",
      call. = FALSE
    )
  }

  target_dir <- dirname(target)
  if (!dir.exists(target_dir)) {
    created <- dir.create(
      target_dir, recursive = TRUE, showWarnings = FALSE
    )
    # Another process may have created the shared cache between the check and
    # dir.create(). Treat that race as success once the directory exists.
    if (!created && !dir.exists(target_dir)) {
      stop("could not create the directory for downloaded RAD whitelists: ",
           target_dir,
           call. = FALSE)
    }
  }
  temp <- tempfile(
    pattern = paste0(".", basename(target), ".part-"),
    tmpdir = target_dir
  )
  on.exit(if (file.exists(temp)) unlink(temp), add = TRUE)

  if (!quiet) {
    mib <- row$bytes_value[[1L]] / 1024^2
    message(
      "Downloading the RAD barcode whitelist for '", row$key[[1L]], "' (",
      format(round(mib, 1L), nsmall = 1L), " MiB). ",
      "rrad will verify it and reuse it for future runs."
    )
  }
  tryCatch(
    .rad_resource_download(row$url[[1L]], temp, quiet = quiet),
    error = function(error) {
      stop(
        "could not download RAD whitelist '", row$key[[1L]], "': ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )
  downloaded <- .rad_resource_integrity(temp, row)
  if (!isTRUE(downloaded$ok)) {
    stop(
      "downloaded RAD whitelist '", row$key[[1L]],
      "' failed integrity verification: ", downloaded$reason,
      call. = FALSE
    )
  }

  # A force refresh of an already-valid immutable object proves that the
  # remote bytes still match the pin; replacing the identical local object is
  # unnecessary and would make a failed refresh capable of damaging it.
  current <- .rad_resource_integrity(target, row)
  if (isTRUE(current$ok)) {
    unlink(temp)
    return(normalizePath(target, winslash = "/", mustWork = TRUE))
  }
  .rad_resource_publish(temp, target, row)
}

#' Download RAD barcode whitelists on first use
#'
#' Most users do not need to call these functions. [rad_demux()] and
#' [rad_scan_wl()] automatically download a large built-in barcode whitelist
#' the first time a kit needs it, then reuse the saved copy on later runs.
#' The download is written to R's standard per-user cache and checked against
#' the exact byte count and SHA-256 digest shipped with `rrad`. Package
#' installation, loading, examples, and tests never require network access.
#'
#' @param key A built-in whitelist name, RAD kit alias, or whitelist filename.
#' @param download Download and save a missing whitelist. If `FALSE`, an
#'   unavailable whitelist returns `NA_character_` without accessing the
#'   network.
#' @param cache_dir Optional cache base directory. The default is
#'   `getOption("rrad.cache_dir")`, then `RRAD_CACHE_DIR`, when set,
#'   otherwise `tools::R_user_dir("rrad", "cache")`.
#' @param quiet Suppress download progress messages.
#' @param kits Optional character vector of whitelist names or kit aliases.
#'   Downloadable names are `"10x_3v3"`, `"10x_3v4"`, and `"10x_5v3"`;
#'   `"10x_3v3.1"` and `"10x_3HTv3.1"` are accepted aliases for
#'   `"10x_3v3"`. `NULL` selects every downloadable whitelist.
#' @param force Re-download and verify whitelists even when valid saved copies
#'   already exist. A failed refresh leaves the valid saved copy untouched.
#' @param confirm Must be `TRUE` to remove downloaded whitelist files.
#'
#' @return `rad_resource_path()` returns one normalized file path, or
#'   `NA_character_` when `download = FALSE` and the whitelist is unavailable.
#'   `rad_download_whitelists()` returns a named character vector of paths.
#'   `rad_clear_resource_cache()` invisibly returns whether a cache was
#'   removed.
#'
#' @details Set `options(rrad.offline = TRUE)` or `RRAD_OFFLINE=true` to
#'   prohibit downloads. Already downloaded whitelists remain available. Set
#'   `options(rrad.cache_dir = "/shared/path")` or `RRAD_CACHE_DIR` to use a
#'   shared or project cache, which can be useful on computing clusters.
#'
#' @examples
#' # Locate a saved whitelist without making a network request.
#' rad_resource_path("10x_3v3", download = FALSE)
#'
#' \dontrun{
#' # Normal use is automatic: rad_demux() and rad_scan_wl() download a
#' # required large built-in whitelist once and reuse it on later calls.
#' whitelist <- rad_resource_path("10x_3v3")
#' rad_download_whitelists(c("10x_3v4", "10x_5v3"))
#'
#' # Optionally download all large whitelists before working offline.
#' rad_download_whitelists()
#' }
#' @name rad-resources
NULL

#' @rdname rad-resources
rad_resource_path <- function(key, download = TRUE, cache_dir = NULL,
                              quiet = FALSE) {
  if (!is.logical(download) || length(download) != 1L || is.na(download)) {
    stop("download must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet)) {
    stop("quiet must be TRUE or FALSE", call. = FALSE)
  }
  row <- .rad_resource_match(key)

  bundled <- .rad_resource_bundled_path(row)
  if (nzchar(bundled)) {
    integrity <- .rad_resource_integrity(bundled, row)
    if (!isTRUE(integrity$ok)) {
      stop(
        "bundled RAD whitelist '", row$key[[1L]],
        "' failed integrity verification: ", integrity$reason,
        call. = FALSE
      )
    }
    return(normalizePath(bundled, winslash = "/", mustWork = TRUE))
  }

  target <- .rad_resource_target(row, cache_dir)
  integrity <- .rad_resource_integrity(target, row)
  if (isTRUE(integrity$ok)) {
    return(normalizePath(target, winslash = "/", mustWork = TRUE))
  }
  if (!download) return(NA_character_)
  .rad_resource_fetch(row, cache_dir = cache_dir, quiet = quiet)
}

#' @rdname rad-resources
rad_download_whitelists <- function(kits = NULL, force = FALSE, quiet = FALSE,
                                    cache_dir = NULL) {
  if (!is.logical(force) || length(force) != 1L || is.na(force)) {
    stop("force must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet)) {
    stop("quiet must be TRUE or FALSE", call. = FALSE)
  }
  manifest <- .rad_resource_manifest()
  if (is.null(kits)) {
    rows <- seq_len(nrow(manifest))
  } else {
    if (!is.character(kits) || anyNA(kits) || any(!nzchar(trimws(kits)))) {
      stop("kits must be NULL or a character vector of whitelist names",
           call. = FALSE)
    }
    rows <- unique(vapply(
      kits,
      function(key) {
        match(.rad_resource_match(key, manifest)$key[[1L]], manifest$key)
      },
      integer(1)
    ))
  }
  if (!length(rows)) return(stats::setNames(character(), character()))

  paths <- vapply(
    rows,
    function(index) {
      row <- manifest[index, , drop = FALSE]
      bundled <- .rad_resource_bundled_path(row)
      if (nzchar(bundled)) {
        integrity <- .rad_resource_integrity(bundled, row)
        if (!isTRUE(integrity$ok)) {
          stop(
            "bundled RAD whitelist '", row$key[[1L]],
            "' failed integrity verification: ", integrity$reason,
            call. = FALSE
          )
        }
        return(normalizePath(bundled, winslash = "/", mustWork = TRUE))
      }
      .rad_resource_fetch(
        row, cache_dir = cache_dir, force = force, quiet = quiet
      )
    },
    character(1)
  )
  stats::setNames(unname(paths), manifest$key[rows])
}

#' @rdname rad-resources
rad_clear_resource_cache <- function(confirm = FALSE, cache_dir = NULL) {
  if (!is.logical(confirm) || length(confirm) != 1L || is.na(confirm)) {
    stop("confirm must be TRUE or FALSE", call. = FALSE)
  }
  if (!confirm) {
    stop(
      "refusing to remove downloaded RAD whitelists; call ",
      "rad_clear_resource_cache(confirm = TRUE) to confirm",
      call. = FALSE
    )
  }
  root <- .rad_resource_cache_root(cache_dir)
  if (!file.exists(root)) return(invisible(FALSE))
  if (!dir.exists(root)) {
    stop("the RAD resource cache path is not a directory: ", root,
         call. = FALSE)
  }
  status <- unlink(root, recursive = TRUE, force = FALSE)
  if (!identical(status, 0L) || file.exists(root)) {
    stop("could not completely remove the RAD resource cache: ", root,
         call. = FALSE)
  }
  invisible(TRUE)
}

.rad_call_with_resource_retry <- function(call, cache_dir = NULL,
                                          quiet = FALSE) {
  if (!is.function(call)) {
    stop("call must be a function with no required arguments", call. = FALSE)
  }
  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet)) {
    stop("quiet must be TRUE or FALSE", call. = FALSE)
  }
  manifest <- .rad_resource_manifest()
  previous <- Sys.getenv("RRAD_RESOURCE_CACHE", unset = NA_character_)
  overlay <- NULL
  on.exit({
    if (is.na(previous)) {
      Sys.unsetenv("RRAD_RESOURCE_CACHE")
    } else {
      Sys.setenv(RRAD_RESOURCE_CACHE = previous)
    }
  }, add = TRUE)
  on.exit({
    if (!is.null(overlay) && dir.exists(overlay)) {
      unlink(overlay, recursive = TRUE, force = FALSE)
    }
  }, add = TRUE)
  # Never let an ambient overlay bypass the manifest pins. Only attach files
  # after this helper has verified their exact contents.
  Sys.unsetenv("RRAD_RESOURCE_CACHE")

  fetched <- character()
  repeat {
    result <- tryCatch(call(), error = function(error) error)
    if (!inherits(result, "error")) return(result)

    message <- conditionMessage(result)
    markers <- gregexpr(
      "RRAD_RESOURCE_MISSING:", message, fixed = TRUE
    )[[1L]]
    marker_count <- if (identical(markers[[1L]], -1L)) 0L else length(markers)
    if (marker_count != 1L) stop(result)
    match <- regexec(
      "(^|: )RRAD_RESOURCE_MISSING:(wl/[[:alnum:]_.-]+)$",
      message,
      perl = TRUE
    )
    groups <- regmatches(message, match)[[1L]]
    if (length(groups) != 3L) stop(result)

    resource_path <- groups[[3L]]
    row <- tryCatch(
      .rad_resource_match(resource_path, manifest),
      error = function(error) NULL
    )
    # An unrecognized marker belongs to a different embedded-core/resource
    # contract. Preserve the original native error instead of guessing.
    if (is.null(row)) stop(result)
    if (resource_path %in% fetched) {
      stop(
        "the embedded RAD core still reports missing resource '",
        resource_path, "' after it was cached and attached",
        call. = FALSE
      )
    }

    cached <- rad_resource_path(
      row$key[[1L]], download = TRUE, cache_dir = cache_dir, quiet = quiet
    )

    # Attach only the file just verified, not its cache directory. A private
    # per-call overlay containing an independent byte copy prevents an
    # unrelated cache entry, concurrent cache mutation, or cache removal from
    # changing resources that the native call is actively consuming.
    if (is.null(overlay)) {
      cache_root <- normalizePath(
        .rad_resource_cache_root(cache_dir), winslash = "/", mustWork = FALSE
      )
      overlay <- tempfile(
        pattern = "rrad-native-overlay-", tmpdir = tempdir()
      )
      normalized_overlay <- normalizePath(
        overlay, winslash = "/", mustWork = FALSE
      )
      inside_cache <- identical(normalized_overlay, cache_root) ||
        startsWith(normalized_overlay, paste0(cache_root, "/"))
      if (inside_cache) {
        # R's session temp directory can be customized. Keep the overlay next
        # to, rather than inside, the subtree removed by cache clearing.
        overlay <- tempfile(
          pattern = ".rrad-native-overlay-", tmpdir = dirname(cache_root)
        )
      }
      created <- dir.create(
        overlay, recursive = FALSE, showWarnings = FALSE, mode = "0700"
      )
      if (!created) {
        stop("could not create the native RAD resource overlay", call. = FALSE)
      }
      Sys.chmod(overlay, mode = "0700", use_umask = FALSE)
      overlay_mode <- file.info(overlay)$mode
      if (is.na(overlay_mode) ||
          bitwAnd(as.integer(overlay_mode), 511L) != 448L) {
        stop("could not secure the native RAD resource overlay", call. = FALSE)
      }
    }
    relative <- strsplit(resource_path, "/", fixed = TRUE)[[1L]]
    overlay_target <- do.call(file.path, as.list(c(overlay, relative)))
    overlay_dir <- dirname(overlay_target)
    if (!dir.exists(overlay_dir)) {
      created <- dir.create(
        overlay_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700"
      )
      if (!created && !dir.exists(overlay_dir)) {
        stop("could not create the native RAD resource overlay directory",
             call. = FALSE)
      }
      Sys.chmod(overlay_dir, mode = "0700", use_umask = FALSE)
      overlay_mode <- file.info(overlay_dir)$mode
      if (is.na(overlay_mode) ||
          bitwAnd(as.integer(overlay_mode), 511L) != 448L) {
        stop("could not secure the native RAD resource overlay directory",
             call. = FALSE)
      }
    }
    copied <- file.copy(
      cached, overlay_target, overwrite = FALSE,
      copy.mode = FALSE, copy.date = FALSE
    )
    if (!copied) {
      stop(
        "could not copy verified RAD resource '", row$key[[1L]],
        "' into the private native overlay",
        call. = FALSE
      )
    }
    attached <- .rad_resource_integrity(overlay_target, row)
    if (!isTRUE(attached$ok)) {
      unlink(overlay_target)
      stop(
        "could not verify attached RAD resource '", row$key[[1L]], "': ",
        attached$reason,
        call. = FALSE
      )
    }
    fetched <- c(fetched, resource_path)
    Sys.setenv(RRAD_RESOURCE_CACHE = overlay)
  }
}
