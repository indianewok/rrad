# Embedded RAD core -------------------------------------------------------

.rad_core_api_version <- 1L

.rad_scalar_character <- function(x, name, allow_null = FALSE) {
  if (allow_null && is.null(x)) return(NULL)
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(name, " must be one non-empty character value", call. = FALSE)
  }
  enc2utf8(x)
}

.rad_scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE", call. = FALSE)
  }
  x
}

.rad_scalar_integer <- function(x, name, minimum = 0L, allow_null = FALSE) {
  if (allow_null && is.null(x)) return(NULL)
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != trunc(x) || x < minimum || x > .Machine$integer.max) {
    stop(
      name, " must be one whole number between ", minimum,
      " and ", .Machine$integer.max,
      call. = FALSE
    )
  }
  as.integer(x)
}

.rad_scalar_double <- function(x, name, minimum = -Inf, allow_null = FALSE) {
  if (allow_null && is.null(x)) return(NULL)
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < minimum) {
    stop(name, " must be one finite number >= ", minimum, call. = FALSE)
  }
  as.double(x)
}

.rad_existing_file <- function(x, name) {
  x <- .rad_scalar_character(x, name)
  if (!file.exists(x) || dir.exists(x)) {
    stop(name, " does not identify an existing file: ", x, call. = FALSE)
  }
  normalizePath(x, mustWork = TRUE)
}

.rad_optional_existing_file <- function(x, name) {
  if (is.null(x)) return(NULL)
  .rad_existing_file(x, name)
}

.rad_layout_reference <- function(x) {
  x <- .rad_scalar_character(x, "layout")
  if (file.exists(x)) {
    if (dir.exists(x)) {
      stop("layout must identify a layout file or a built-in layout key",
           call. = FALSE)
    }
    return(normalizePath(x, mustWork = TRUE))
  }
  if (grepl("[/\\\\]", x) || grepl("\\.(csv|ya?ml)$", x, ignore.case = TRUE)) {
    stop("layout file does not exist: ", x, call. = FALSE)
  }
  x
}

.rad_output_prefix <- function(x) {
  x <- .rad_scalar_character(x, "output")
  if (x %in% c(".", "..") || grepl("[/\\\\]", x)) {
    stop("output must be a filename prefix within out_dir", call. = FALSE)
  }
  x
}

.rad_output_directory <- function(x) {
  x <- .rad_scalar_character(x, "out_dir")
  if (file.exists(x) && !dir.exists(x)) {
    stop("out_dir identifies a file, not a directory: ", x, call. = FALSE)
  }
  if (!dir.exists(x) && !dir.create(x, recursive = TRUE, showWarnings = FALSE)) {
    stop("could not create out_dir: ", x, call. = FALSE)
  }
  normalizePath(x, mustWork = TRUE)
}

.rad_drop_null <- function(x) {
  x[!vapply(x, is.null, logical(1))]
}

.rad_validate_core_info <- function(info) {
  required <- c(
    "core_version", "core_api_version", "source_commit",
    "embedded_source_digest",
    "resource_schema", "resource_digest", "features", "build", "vendored"
  )
  if (!is.list(info) || is.null(names(info)) ||
      length(missing <- setdiff(required, names(info)))) {
    stop(
      "the embedded RAD core returned invalid metadata",
      if (exists("missing", inherits = FALSE) && length(missing)) {
        paste0("; missing: ", paste(missing, collapse = ", "))
      } else {
        ""
      },
      call. = FALSE
    )
  }
  if (!is.numeric(info$core_api_version) ||
      length(info$core_api_version) != 1L ||
      is.na(info$core_api_version) ||
      !is.finite(info$core_api_version) ||
      info$core_api_version != trunc(info$core_api_version) ||
      info$core_api_version < 1L ||
      info$core_api_version > .Machine$integer.max) {
    stop("the embedded RAD core returned an invalid core_api_version",
         call. = FALSE)
  }
  if (as.integer(info$core_api_version) != .rad_core_api_version) {
    stop(
      "incompatible embedded RAD core API: package expects ",
      .rad_core_api_version, " but core reports ", info$core_api_version,
      call. = FALSE
    )
  }
  for (field in c("core_version", "source_commit")) {
    if (!is.character(info[[field]]) || length(info[[field]]) != 1L ||
        is.na(info[[field]]) || !nzchar(info[[field]])) {
      stop("the embedded RAD core returned an invalid ", field, call. = FALSE)
    }
  }
  if (!grepl("^[[:xdigit:]]{40}$", info$source_commit)) {
    stop("the embedded RAD core returned an invalid source_commit",
         call. = FALSE)
  }
  if (!is.character(info$embedded_source_digest) ||
      length(info$embedded_source_digest) != 1L ||
      is.na(info$embedded_source_digest) ||
      !grepl("^[[:xdigit:]]{64}$", info$embedded_source_digest)) {
    stop("the embedded RAD core returned an invalid embedded_source_digest",
         call. = FALSE)
  }
  if (!is.numeric(info$resource_schema) || length(info$resource_schema) != 1L ||
      is.na(info$resource_schema) || !is.finite(info$resource_schema) ||
      info$resource_schema != trunc(info$resource_schema) ||
      info$resource_schema < 1L || info$resource_schema > .Machine$integer.max) {
    stop("the embedded RAD core returned an invalid resource_schema",
         call. = FALSE)
  }
  if (!is.character(info$resource_digest) ||
      length(info$resource_digest) != 1L || is.na(info$resource_digest) ||
      !grepl("^[[:xdigit:]]{64}$", info$resource_digest)) {
    stop("the embedded RAD core returned an invalid resource_digest",
         call. = FALSE)
  }
  if (!is.character(info$features) || anyNA(info$features)) {
    stop("the embedded RAD core returned an invalid feature list",
         call. = FALSE)
  }
  if (!("demux" %in% info$features)) {
    stop("the embedded RAD core does not advertise the demux feature",
         call. = FALSE)
  }
  if (!is.list(info$build)) {
    stop("the embedded RAD core returned invalid build metadata", call. = FALSE)
  }
  if (!is.logical(info$build$resources_available) ||
      length(info$build$resources_available) != 1L ||
      is.na(info$build$resources_available) ||
      !isTRUE(info$build$resources_available) ||
      !is.character(info$build$resource_dir) ||
      length(info$build$resource_dir) != 1L ||
      is.na(info$build$resource_dir) || !nzchar(info$build$resource_dir)) {
    stop("the embedded RAD core's bundled resources are unavailable",
         call. = FALSE)
  }
  if (!is.logical(info$vendored) || length(info$vendored) != 1L ||
      is.na(info$vendored) || !isTRUE(info$vendored)) {
    stop("the embedded RAD core is not identified as vendored", call. = FALSE)
  }
  info$core_api_version <- as.integer(info$core_api_version)
  info$resource_schema <- as.integer(info$resource_schema)
  info$features <- unique(info$features)
  structure(info, class = c("rad_core_info", "list"))
}

#' Inspect the embedded RAD demultiplexing core
#'
#' Returns the compile-time identity and compatibility information for the
#' copy of RAD packaged with this R package. No external executable or PATH
#' lookup is involved.
#'
#' @return A `rad_core_info` list containing the core version, core API
#'   version, source commit, resource schema and digest, capabilities, and
#'   build details.
#' @name rad-core
#' @rdname rad-core
rad_core_info <- function() {
  .rad_validate_core_info(rad_core_info_cpp())
}

#' Demultiplex reads with the embedded RAD core
#'
#' Runs RAD in-process through its packaged C++ core. Inputs and output
#' locations are checked in R before native execution. The returned object is
#' structured and records the exact embedded core identity used for the run.
#'
#' @param layout A built-in RAD layout key or path to a layout CSV/seqspec.
#' @param fastq Path to the input FASTQ or gzip-compressed FASTQ.
#' @param out_dir Output directory. It is created recursively when needed.
#' @param output Output filename prefix within `out_dir`.
#' @param kit Optional built-in global-whitelist kit key. This is an
#'   alternative to `global_whitelist`; supplying both is an error.
#' @param global_whitelist Optional path to a global whitelist.
#' @param custom_whitelist Optional path to a custom whitelist.
#' @param threads Positive number of worker threads.
#'   Builds without OpenMP (including the default macOS build) run serially;
#'   requesting more than one thread warns and records an effective value of
#'   one in the result.
#' @param chunk_size Positive number of reads per processing chunk.
#' @param max_reads Optional positive maximum number of reads to process;
#'   `NULL` processes the complete input.
#' @param auto_whitelist Run RAD's barcode-whitelist scan before demultiplexing.
#'   This is an alternative to `custom_whitelist`; supplying both is an error.
#' @param scan_adapter Optional adapter used by automatic whitelist scanning.
#'   `NULL` derives the nearest upstream static element from `layout`.
#' @param scan_barcode_length Optional barcode length from 1 to 32 for
#'   automatic scanning. `NULL` derives it from `layout`.
#' @param scan_max_error Adapter edit-distance ratio from 0 to 1 used during
#'   automatic scanning.
#' @param scan_max_reads Optional positive scan-specific read cap. `NULL`
#'   inherits `max_reads`.
#' @param scan_chunk_size Optional positive scan-specific chunk size. `NULL`
#'   inherits `chunk_size`.
#' @param scan_threads Optional positive scan-specific thread count. `NULL`
#'   inherits `threads`.
#' @param scan_selection Barcode-selection policy for automatic scanning:
#'   `"high_specificity"` or `"above_floor"`.
#' @param whitelist_mutation,generated_mutation Optional non-negative mutation
#'   distances.
#' @param bc_correction_mode Barcode correction policy, either `"offensive"`
#'   or `"defensive"`.
#' @param joint_bc_mode Multi-barcode policy, either `"default"` or
#'   `"strict"`.
#' @param rc_umi Reverse-complement UMIs extracted from reverse-oriented reads.
#' @param min_read_length Optional non-negative minimum cDNA length.
#' @param write_debug Write RAD debug artifacts in addition to the primary
#'   demultiplexed FASTQ.
#' @param verbose Print native progress captured during the completed run.
#' @param reuse_cache Reuse an existing `output` layout and position-map cache.
#'   The default, `FALSE`, rebuilds these files so an output prefix cannot
#'   silently select stale layout state from an earlier run.
#' @return A `rad_demux_result` list with `status`, `backend`, `core`, `input`,
#'   `config`, `artifacts`, `stats`, `scan`, and `log` fields. `scan` is `NULL`
#'   unless `auto_whitelist = TRUE`, when it is a complete
#'   `rad_scan_wl_result`; `config$scan_effective` records its layout-derived
#'   and explicit effective settings.
#' @details Small layouts and whitelist tables ship with `rrad`. A required
#'   3M whitelist is downloaded on first use from its commit-pinned RAD URL,
#'   verified by byte count and SHA-256, and reused on later runs. Use
#'   [rad_download_whitelists()] to download whitelists before offline runs.
#'   Explicit whitelist paths are used directly and are never copied.
#' @rdname rad-core
rad_demux <- function(layout,
                      fastq,
                      out_dir = ".",
                      output = "demux",
                      kit = NULL,
                      global_whitelist = NULL,
                      custom_whitelist = NULL,
                      threads = 1L,
                      chunk_size = 5000L,
                      max_reads = NULL,
                      auto_whitelist = FALSE,
                      scan_adapter = NULL,
                      scan_barcode_length = NULL,
                      scan_max_error = 0.3,
                      scan_max_reads = NULL,
                      scan_chunk_size = NULL,
                      scan_threads = NULL,
                      scan_selection = c("high_specificity", "above_floor"),
                      whitelist_mutation = NULL,
                      generated_mutation = NULL,
                      bc_correction_mode = "offensive",
                      joint_bc_mode = "default",
                      rc_umi = TRUE,
                      min_read_length = NULL,
                      write_debug = FALSE,
                      verbose = FALSE,
                      reuse_cache = FALSE) {
  layout <- .rad_layout_reference(layout)
  fastq <- .rad_existing_file(fastq, "fastq")
  out_dir <- .rad_scalar_character(out_dir, "out_dir")
  output <- .rad_output_prefix(output)
  kit <- .rad_scalar_character(kit, "kit", allow_null = TRUE)
  if (!is.null(kit) && grepl(":", kit, fixed = TRUE)) {
    stop("kit must be one built-in kit key without ':'", call. = FALSE)
  }
  global_whitelist <- .rad_optional_existing_file(
    global_whitelist, "global_whitelist"
  )
  custom_whitelist <- .rad_optional_existing_file(
    custom_whitelist, "custom_whitelist"
  )
  if (!is.null(kit) && !is.null(global_whitelist)) {
    stop("kit and global_whitelist are alternatives; supply only one",
         call. = FALSE)
  }
  threads <- .rad_scalar_integer(threads, "threads", minimum = 1L)
  chunk_size <- .rad_scalar_integer(chunk_size, "chunk_size", minimum = 1L)
  max_reads <- .rad_scalar_integer(
    max_reads, "max_reads", minimum = 1L, allow_null = TRUE
  )
  auto_whitelist <- .rad_scalar_logical(auto_whitelist, "auto_whitelist")
  if (auto_whitelist && !is.null(custom_whitelist)) {
    stop("auto_whitelist and custom_whitelist are alternatives; supply only one",
         call. = FALSE)
  }
  scan_adapter <- .rad_scalar_character(
    scan_adapter, "scan_adapter", allow_null = TRUE
  )
  scan_barcode_length <- .rad_scalar_integer(
    scan_barcode_length, "scan_barcode_length", minimum = 1L,
    allow_null = TRUE
  )
  if (!is.null(scan_barcode_length) && scan_barcode_length > 32L) {
    stop("scan_barcode_length must be between 1 and 32", call. = FALSE)
  }
  scan_max_error <- .rad_scalar_double(
    scan_max_error, "scan_max_error", minimum = 0
  )
  if (scan_max_error > 1) {
    stop("scan_max_error must be between 0 and 1", call. = FALSE)
  }
  scan_max_reads <- .rad_scalar_integer(
    scan_max_reads, "scan_max_reads", minimum = 1L, allow_null = TRUE
  )
  scan_chunk_size <- .rad_scalar_integer(
    scan_chunk_size, "scan_chunk_size", minimum = 1L, allow_null = TRUE
  )
  scan_threads <- .rad_scalar_integer(
    scan_threads, "scan_threads", minimum = 1L, allow_null = TRUE
  )
  scan_selection <- match.arg(tolower(scan_selection),
                              c("high_specificity", "above_floor"))
  scan_options_requested <-
    !is.null(scan_adapter) || !is.null(scan_barcode_length) ||
    !identical(scan_max_error, 0.3) || !is.null(scan_max_reads) ||
    !is.null(scan_chunk_size) || !is.null(scan_threads) ||
    !identical(scan_selection, "high_specificity")
  if (!auto_whitelist && scan_options_requested) {
    stop("scan_* options require auto_whitelist = TRUE", call. = FALSE)
  }
  whitelist_mutation <- .rad_scalar_integer(
    whitelist_mutation, "whitelist_mutation", minimum = 0L,
    allow_null = TRUE
  )
  generated_mutation <- .rad_scalar_integer(
    generated_mutation, "generated_mutation", minimum = 0L,
    allow_null = TRUE
  )
  bc_correction_mode <- .rad_scalar_character(
    bc_correction_mode, "bc_correction_mode"
  )
  bc_correction_mode <- match.arg(
    tolower(bc_correction_mode), c("offensive", "defensive")
  )
  joint_bc_mode <- .rad_scalar_character(joint_bc_mode, "joint_bc_mode")
  joint_bc_mode <- match.arg(tolower(joint_bc_mode), c("default", "strict"))
  rc_umi <- .rad_scalar_logical(rc_umi, "rc_umi")
  min_read_length <- .rad_scalar_integer(
    min_read_length, "min_read_length", minimum = 0L, allow_null = TRUE
  )
  write_debug <- .rad_scalar_logical(write_debug, "write_debug")
  verbose <- .rad_scalar_logical(verbose, "verbose")
  reuse_cache <- .rad_scalar_logical(reuse_cache, "reuse_cache")
  out_dir <- .rad_output_directory(out_dir)

  # Query immediately before execution. This is intentionally a runtime
  # handshake rather than an assumption based on the R package version.
  core <- rad_core_info()
  if (auto_whitelist && !("auto-whitelist" %in% core$features)) {
    stop(
      "the embedded RAD core does not advertise the auto-whitelist feature",
      call. = FALSE
    )
  }
  threads_requested <- threads
  threads_effective <- if (isTRUE(core$build$openmp)) threads else 1L
  if (threads_requested > 1L && threads_effective == 1L) {
    warning(
      "this rrad build has no OpenMP support; using one effective thread",
      call. = FALSE
    )
  }
  scan_threads_requested <- if (is.null(scan_threads)) {
    threads_requested
  } else {
    scan_threads
  }
  scan_threads_effective <- if (isTRUE(core$build$openmp)) {
    scan_threads_requested
  } else {
    1L
  }
  if (auto_whitelist && scan_threads_requested > 1L &&
      scan_threads_effective == 1L && threads_requested == 1L) {
    warning(
      "this rrad build has no OpenMP support; using one effective scan thread",
      call. = FALSE
    )
  }

  native_options <- .rad_drop_null(list(
    layout = layout,
    fastq = fastq,
    output_dir = out_dir,
    output_prefix = output,
    kit = kit,
    global_whitelist = global_whitelist,
    custom_whitelist = custom_whitelist,
    threads = threads_requested,
    chunk_size = chunk_size,
    max_reads = max_reads,
    auto_whitelist = auto_whitelist,
    whitelist_mutation = whitelist_mutation,
    generated_mutation = generated_mutation,
    bc_correction_mode = bc_correction_mode,
    joint_bc_mode = joint_bc_mode,
    rc_umi = rc_umi,
    min_read_length = min_read_length,
    write_debug = write_debug,
    verbose = verbose,
    reuse_cache = reuse_cache
  ))
  if (auto_whitelist) {
    native_options <- c(native_options, .rad_drop_null(list(
      scan_adapter = scan_adapter,
      scan_barcode_length = scan_barcode_length,
      scan_max_error = scan_max_error,
      scan_max_reads = scan_max_reads,
      scan_chunk_size = scan_chunk_size,
      scan_threads = scan_threads,
      scan_selection = scan_selection
    )))
  }

  native <- .rad_call_with_resource_retry(
    function() rad_demux_cpp(native_options)
  )
  required_result <- c(
    "success", "backend", "core", "input", "layout", "output_prefix",
    "files", "stats", "scan", "log"
  )
  if (!is.list(native) || is.null(names(native)) ||
      length(setdiff(required_result, names(native)))) {
    stop("the embedded RAD core returned an invalid demultiplexing result",
         call. = FALSE)
  }
  if (!isTRUE(native$success)) {
    stop("the embedded RAD core reported an unsuccessful run", call. = FALSE)
  }
  if (!is.character(native$backend) || length(native$backend) != 1L ||
      is.na(native$backend) || !identical(native$backend, "embedded")) {
    stop("the embedded RAD core returned an invalid backend", call. = FALSE)
  }
  if (!is.character(native$input) || length(native$input) != 1L ||
      is.na(native$input) || !nzchar(native$input) ||
      !is.character(native$layout) || length(native$layout) != 1L ||
      is.na(native$layout) || !nzchar(native$layout)) {
    stop("the embedded RAD core returned invalid input metadata", call. = FALSE)
  }
  if (!is.character(native$output_prefix) ||
      length(native$output_prefix) != 1L || is.na(native$output_prefix) ||
      !nzchar(native$output_prefix)) {
    stop("the embedded RAD core returned an invalid output prefix",
         call. = FALSE)
  }
  if (!is.list(native$files) || is.null(names(native$files)) ||
      !is.list(native$stats) || is.null(names(native$stats)) ||
      !is.character(native$log) || length(native$log) != 1L ||
      is.na(native$log)) {
    stop("the embedded RAD core returned invalid artifacts or statistics",
         call. = FALSE)
  }

  config <- list(
    out_dir = out_dir,
    output = output,
    kit = kit,
    global_whitelist = global_whitelist,
    custom_whitelist = custom_whitelist,
    threads = threads_effective,
    threads_requested = threads_requested,
    threads_effective = threads_effective,
    chunk_size = chunk_size,
    max_reads = max_reads,
    auto_whitelist = auto_whitelist,
    scan_adapter = scan_adapter,
    scan_barcode_length = scan_barcode_length,
    scan_max_error = scan_max_error,
    scan_max_reads = if (is.null(scan_max_reads)) max_reads else scan_max_reads,
    scan_chunk_size = if (is.null(scan_chunk_size)) {
      chunk_size
    } else {
      scan_chunk_size
    },
    scan_threads = scan_threads_effective,
    scan_threads_requested = scan_threads_requested,
    scan_threads_effective = scan_threads_effective,
    scan_selection = scan_selection,
    whitelist_mutation = whitelist_mutation,
    generated_mutation = generated_mutation,
    bc_correction_mode = bc_correction_mode,
    joint_bc_mode = joint_bc_mode,
    rc_umi = rc_umi,
    min_read_length = min_read_length,
    write_debug = write_debug,
    verbose = verbose,
    reuse_cache = reuse_cache
  )

  native_core <- .rad_validate_core_info(native$core)
  handshake_fields <- c(
    "core_version", "core_api_version", "source_commit",
    "embedded_source_digest", "resource_schema", "resource_digest"
  )
  if (!identical(unclass(core[handshake_fields]),
                 unclass(native_core[handshake_fields]))) {
    stop("the embedded RAD core identity changed during the run", call. = FALSE)
  }
  if (auto_whitelist) {
    if (!is.list(native$scan) || is.null(names(native$scan)) ||
        !is.list(native$scan$config) || is.null(names(native$scan$config))) {
      stop(
        "the embedded RAD core returned invalid auto-whitelist metadata",
        call. = FALSE
      )
    }
    scan_result <- .rad_scan_result_from_native(native$scan, core)
    config$scan_effective <- scan_result$config
  } else {
    if (!is.null(native$scan)) {
      stop(
        "the embedded RAD core returned unexpected auto-whitelist metadata",
        call. = FALSE
      )
    }
    scan_result <- NULL
  }
  result <- list(
    status = "success",
    backend = native$backend,
    core = native_core,
    input = list(
      fastq = native$input,
      layout = native$layout
    ),
    config = config,
    artifacts = native$files,
    stats = native$stats,
    scan = scan_result,
    log = native$log
  )
  result <- structure(result, class = c("rad_demux_result", "list"))
  if (verbose && nzchar(native$log)) {
    cat(native$log)
    if (!grepl("\n$", native$log)) cat("\n")
  }
  result
}
