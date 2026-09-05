# Native RAD whitelist scanning ----------------------------------------------

.rad_scan_scalar_double <- function(x, name, minimum = -Inf,
                                    maximum = Inf) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < minimum || x > maximum) {
    stop(
      name, " must be one finite number between ", minimum, " and ",
      maximum,
      call. = FALSE
    )
  }
  as.double(x)
}

.rad_scan_whitelist_reference <- function(x, name) {
  if (is.null(x)) return(NULL)
  x <- .rad_scalar_character(x, name)
  if (grepl(":", x, fixed = TRUE)) {
    stop(
      name,
      " must be one whitelist file or built-in kit key; ':' specifications ",
      "are only valid inside RAD read-layout files",
      call. = FALSE
    )
  }
  if (file.exists(x)) {
    if (dir.exists(x)) {
      stop(name, " must identify a whitelist file or built-in kit key",
           call. = FALSE)
    }
    return(normalizePath(x, mustWork = TRUE))
  }
  if (grepl("[/\\\\]", x) ||
      grepl("\\.(csv|tsv|txt|gz)$", x, ignore.case = TRUE)) {
    stop(name, " file does not exist: ", x, call. = FALSE)
  }
  x
}

.rad_scan_sequence_file <- function(x, name) {
  path <- .rad_existing_file(x, name)
  if (!grepl("\\.(fastq|fq|fasta|fa)(\\.gz)?$", path,
             ignore.case = TRUE)) {
    stop(
      name,
      " must have a .fastq, .fq, .fasta, or .fa extension, optionally .gz",
      call. = FALSE
    )
  }
  path
}

.rad_scan_output_prefix <- function(x) {
  x <- .rad_scalar_character(x, "output_prefix")
  expanded <- path.expand(x)
  if (basename(expanded) %in% c(".", "..") ||
      grepl("[/\\\\]$", expanded)) {
    stop("output_prefix must end in a filename prefix", call. = FALSE)
  }
  parent <- dirname(expanded)
  if (file.exists(parent) && !dir.exists(parent)) {
    stop("the output_prefix parent identifies a file: ", parent,
         call. = FALSE)
  }
  if (!dir.exists(parent) &&
      !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("could not create the output_prefix directory: ", parent,
         call. = FALSE)
  }
  file.path(normalizePath(parent, mustWork = TRUE), basename(expanded))
}

.rad_scan_default_rescan_prefix <- function(input) {
  prefix <- sub("\\.gz$", "", input, ignore.case = TRUE)
  extension <- tools::file_ext(prefix)
  if (nzchar(extension)) {
    prefix <- substr(prefix, 1L, nchar(prefix) - nchar(extension) - 1L)
  }
  paste0(prefix, "_rescan")
}

.rad_scan_parse_batch <- function(batch_csv, split_mode = FALSE,
                                  protected_paths = character()) {
  batch <- tryCatch(
    utils::read.csv(
      batch_csv,
      header = FALSE,
      stringsAsFactors = FALSE,
      colClasses = "character",
      check.names = FALSE,
      strip.white = TRUE,
      blank.lines.skip = TRUE,
      comment.char = "",
      na.strings = character()
    ),
    error = function(error) {
      stop("could not parse batch_csv: ", conditionMessage(error),
           call. = FALSE)
    }
  )
  if (ncol(batch) != 2L || nrow(batch) == 0L) {
    stop("batch_csv must contain exactly two non-empty columns",
         call. = FALSE)
  }

  first <- tolower(trimws(unlist(batch[1L, ], use.names = FALSE)))
  header_input <- c("input", "fastq", "input_fastq", "fastq_path")
  header_output <- c("output", "prefix", "output_prefix")
  if (first[[1L]] %in% header_input && first[[2L]] %in% header_output) {
    batch <- batch[-1L, , drop = FALSE]
  }
  if (nrow(batch) == 0L) {
    stop("batch_csv contains a header but no runs", call. = FALSE)
  }

  input <- trimws(batch[[1L]])
  output_prefix <- trimws(batch[[2L]])
  invalid <- is.na(input) | !nzchar(input) |
    is.na(output_prefix) | !nzchar(output_prefix)
  if (any(invalid)) {
    stop(
      "batch_csv has an empty input or output prefix on data row ",
      which(invalid)[[1L]],
      call. = FALSE
    )
  }

  # Validate the complete batch before the first native run so a bad later row
  # cannot leave an avoidable partially completed batch.
  input <- vapply(
    seq_along(input),
    function(index) {
      .rad_scan_sequence_file(
        input[[index]], paste0("batch input row ", index)
      )
    },
    character(1)
  )
  output_prefix <- vapply(
    output_prefix,
    .rad_scan_output_prefix,
    character(1)
  )
  duplicated_prefix <- duplicated(output_prefix) |
    duplicated(output_prefix, fromLast = TRUE)
  if (any(duplicated_prefix)) {
    stop("batch_csv output prefixes must be unique", call. = FALSE)
  }

  suffixes <- c(".csv", ".txt", "_scan_wl.log")
  if (split_mode) {
    suffixes <- c(suffixes, "_valid_pairs.csv", "_spat_mask.csv")
  }
  artifacts <- unlist(lapply(
    output_prefix,
    function(prefix) paste0(prefix, suffixes)
  ), use.names = FALSE)
  existing_directories <- dir.exists(artifacts)
  if (any(existing_directories)) {
    stop(
      "a batch output artifact identifies a directory: ",
      artifacts[[which(existing_directories)[[1L]]]],
      call. = FALSE
    )
  }
  comparable_artifacts <- vapply(
    artifacts,
    function(path) {
      if (file.exists(path)) normalizePath(path, mustWork = TRUE) else path
    },
    character(1)
  )
  # Use a conservative case-folded key so a batch is portable to the default
  # case-insensitive macOS filesystem and cannot overwrite an earlier row.
  artifact_keys <- tolower(comparable_artifacts)
  duplicated_artifact <- duplicated(artifact_keys) |
    duplicated(artifact_keys, fromLast = TRUE)
  if (any(duplicated_artifact)) {
    stop(
      "batch_csv run prefixes produce colliding output artifacts: ",
      artifacts[[which(duplicated_artifact)[[1L]]]],
      call. = FALSE
    )
  }

  protected_paths <- unlist(protected_paths, use.names = FALSE)
  if (length(protected_paths)) {
    protected_paths <- protected_paths[
      !is.na(protected_paths) & nzchar(protected_paths) &
        file.exists(protected_paths)
    ]
  } else {
    protected_paths <- character()
  }
  protected_inputs <- unique(c(batch_csv, input, protected_paths))
  comparable_inputs <- vapply(
    protected_inputs, normalizePath, character(1), mustWork = TRUE
  )
  collision <- artifact_keys %in% tolower(comparable_inputs)
  if (any(collision)) {
    stop(
      paste0(
        "a batch output artifact would overwrite an input or batch_csv",
        "/reference whitelist: "
      ),
      artifacts[[which(collision)[[1L]]]],
      call. = FALSE
    )
  }
  if (.rad_any_file_equivalent_cpp(artifacts, protected_inputs)) {
    stop(
      paste0(
        "a batch output artifact would overwrite an input or batch_csv",
        "/reference whitelist through an existing file alias"
      ),
      call. = FALSE
    )
  }

  data.frame(
    input = unname(input),
    output_prefix = unname(output_prefix),
    stringsAsFactors = FALSE
  )
}

.rad_scan_require_feature <- function(core) {
  if (!is.character(core$features) || anyNA(core$features) ||
      !("scan-wl" %in% core$features)) {
    stop("the embedded RAD core does not advertise the scan-wl feature",
         call. = FALSE)
  }
}

.rad_scan_validate_native_config <- function(config) {
  required <- c(
    "input", "output_prefix", "adapter_seq", "barcode_length",
    "left_margin", "right_margin", "max_reads", "max_error", "whitelist",
    "threads", "threads_requested", "threads_effective", "chunk_size",
    "selection", "rescan", "bc1_whitelist", "bc2_whitelist",
    "bc1_length", "bc2_length", "umi_length", "offset_min", "offset_max",
    "verbose"
  )
  if (!is.list(config) || is.null(names(config)) ||
      length(setdiff(required, names(config)))) {
    stop("the embedded RAD core returned an invalid scan configuration",
         call. = FALSE)
  }
  if (!is.character(config$input) || length(config$input) != 1L ||
      is.na(config$input) || !nzchar(config$input) ||
      !is.character(config$output_prefix) ||
      length(config$output_prefix) != 1L || is.na(config$output_prefix) ||
      !nzchar(config$output_prefix) ||
      !is.character(config$adapter_seq) || length(config$adapter_seq) != 1L ||
      is.na(config$adapter_seq) ||
      !is.character(config$selection) || length(config$selection) != 1L ||
      is.na(config$selection) ||
      !(config$selection %in% c("high_specificity", "above_floor")) ||
      !is.logical(config$rescan) || length(config$rescan) != 1L ||
      is.na(config$rescan) ||
      !is.logical(config$verbose) || length(config$verbose) != 1L ||
      is.na(config$verbose)) {
    stop("the embedded RAD core returned malformed scan configuration values",
         call. = FALSE)
  }
  config
}

.rad_scan_result_from_native <- function(native, core) {
  required <- c(
    "success", "backend", "core", "input", "output_prefix",
    "config", "files", "stats", "log"
  )
  if (!is.list(native) || is.null(names(native)) ||
      length(setdiff(required, names(native)))) {
    stop("the embedded RAD core returned an invalid whitelist-scan result",
         call. = FALSE)
  }
  if (!isTRUE(native$success)) {
    stop("the embedded RAD core reported an unsuccessful whitelist scan",
         call. = FALSE)
  }
  if (!is.character(native$backend) || length(native$backend) != 1L ||
      is.na(native$backend) || !identical(native$backend, "embedded")) {
    stop("the embedded RAD core returned an invalid scan backend",
         call. = FALSE)
  }
  if (!is.character(native$input) || length(native$input) != 1L ||
      is.na(native$input) || !nzchar(native$input) ||
      !is.character(native$output_prefix) ||
      length(native$output_prefix) != 1L || is.na(native$output_prefix) ||
      !nzchar(native$output_prefix)) {
    stop("the embedded RAD core returned invalid scan path metadata",
         call. = FALSE)
  }
  if (!is.list(native$files) || is.null(names(native$files)) ||
      !is.list(native$stats) || is.null(names(native$stats)) ||
      !is.character(native$log) || length(native$log) != 1L ||
      is.na(native$log)) {
    stop("the embedded RAD core returned invalid scan artifacts or statistics",
         call. = FALSE)
  }

  effective_config <- .rad_scan_validate_native_config(native$config)
  if (!identical(effective_config$input, native$input) ||
      !identical(effective_config$output_prefix, native$output_prefix)) {
    stop("the embedded RAD core returned inconsistent scan configuration paths",
         call. = FALSE)
  }

  native_core <- .rad_validate_core_info(native$core)
  handshake_fields <- c(
    "core_version", "core_api_version", "source_commit",
    "embedded_source_digest", "resource_schema", "resource_digest"
  )
  if (!identical(unclass(core[handshake_fields]),
                 unclass(native_core[handshake_fields]))) {
    stop("the embedded RAD core identity changed during the whitelist scan",
         call. = FALSE)
  }

  result <- structure(
    list(
      status = "success",
      backend = native$backend,
      core = native_core,
      input = native$input,
      output_prefix = native$output_prefix,
      config = effective_config,
      artifacts = native$files,
      stats = native$stats,
      log = native$log
    ),
    class = c("rad_scan_wl_result", "list")
  )
  if (isTRUE(effective_config$verbose) && nzchar(native$log)) {
    cat(native$log)
    if (!grepl("\n$", native$log)) cat("\n")
  }
  result
}

.rad_scan_run <- function(options, core) {
  native <- .rad_call_with_resource_retry(
    function() .rad_scan_wl_cpp(options)
  )
  .rad_scan_result_from_native(native, core)
}

#' Discover a barcode whitelist with the embedded RAD core
#'
#' Scans FASTQ or FASTA reads for barcodes following an adapter, optionally restricts
#' observations to a reference whitelist, and writes RAD's selected whitelist,
#' per-barcode statistics, and scan summary without launching an external RAD
#' executable. Two-part BC1/BC2 scans, statistical rescans, and batch input are
#' supported by the same native engine.
#'
#' Small reference whitelists ship with `rrad`. A requested built-in 3M
#' whitelist is downloaded once, verified by byte count and SHA-256, and
#' reused on later runs. Use [rad_download_whitelists()] to download the
#' whitelists ahead of an offline or batch run.
#'
#' @param input Input FASTQ, FASTA, or gzip-compressed FASTQ/FASTA path. With
#'   `rescan = TRUE`, a CSV, TSV, or whitespace-delimited barcode/count table
#'   instead. Supply exactly one of `input` and `batch_csv`.
#' @param output_prefix Prefix for generated files. The default is `"barcodes"`
#'   for sequence scans and `<input>_rescan` for statistical rescans.
#' @param batch_csv Optional two-column CSV containing an input FASTQ/FASTA
#'   path and output prefix per row. A conventional header is accepted. Every
#'   row is validated before the first native scan begins.
#' @param adapter_seq Adapter or primer immediately upstream of the barcode.
#'   Required for sequence scans and unused for statistical rescans.
#' @param barcode_length Number of barcode bases to extract, from 1 to 32.
#' @param left_margin,right_margin Non-negative RAD extraction offsets. For a
#'   forward adapter, `left_margin` skips bases before extraction and
#'   `right_margin` extends its end; reverse matches use the corresponding
#'   upstream window.
#' @param max_reads Optional positive maximum number of reads to inspect;
#'   `NULL` scans the entire input.
#' @param max_error Maximum adapter edit-distance ratio from 0 to 1.
#' @param whitelist Optional built-in whitelist key or whitelist file used to
#'   restrict and validate single-barcode observations.
#' @param threads Positive requested worker-thread count. A build without
#'   OpenMP warns and uses one effective thread.
#' @param chunk_size Positive number of reads per processing chunk.
#' @param selection Which barcode set is written to the selected `.txt` file:
#'   RAD's high-specificity final calls or every above-floor barcode.
#' @param rescan Re-run only RAD's statistical selector on barcode/count data;
#'   no sequence records are scanned.
#' @param bc1_whitelist,bc2_whitelist Optional built-in keys or files for a
#'   two-part barcode scan. Supply both or neither. Two-part mode is an
#'   alternative to `whitelist` and is unavailable for rescans.
#' @param bc1_length,bc2_length Optional packed-whitelist barcode lengths,
#'   each from 1 through 32. They are needed for custom numeric bit-lists,
#'   which do not encode sequence length; bundled SPLiT-seq keys infer 8.
#' @param umi_length Non-negative UMI length between the adapter and BC1 in
#'   two-part mode.
#' @param offset_min,offset_max Inclusive BC1 search-offset range. The minimum
#'   must not exceed the maximum.
#' @param verbose Print the native progress captured during each completed run.
#' @return A `rad_scan_wl_result` containing provenance, the effective native
#'   configuration (including inferred split-barcode lengths), artifact paths,
#'   scan statistics, and captured log output. Batch mode
#'   returns a `rad_scan_wl_batch_result` whose `runs` element contains one
#'   such result per validated CSV row.
#' @name rad-scan-wl
#' @rdname rad-scan-wl
rad_scan_wl <- function(input = NULL,
                        output_prefix = NULL,
                        batch_csv = NULL,
                        adapter_seq = NULL,
                        barcode_length = 16L,
                        left_margin = 0L,
                        right_margin = 0L,
                        max_reads = NULL,
                        max_error = 0.3,
                        whitelist = NULL,
                        threads = 1L,
                        chunk_size = 10000L,
                        selection = c("high_specificity", "above_floor"),
                        rescan = FALSE,
                        bc1_whitelist = NULL,
                        bc2_whitelist = NULL,
                        bc1_length = NULL,
                        bc2_length = NULL,
                        umi_length = 9L,
                        offset_min = 0L,
                        offset_max = 3L,
                        verbose = FALSE) {
  using_input <- !is.null(input)
  using_batch <- !is.null(batch_csv)
  if (using_input == using_batch) {
    stop("supply exactly one of input and batch_csv", call. = FALSE)
  }

  rescan <- .rad_scalar_logical(rescan, "rescan")
  verbose <- .rad_scalar_logical(verbose, "verbose")
  if (rescan && using_batch) {
    stop("rescan does not support batch_csv", call. = FALSE)
  }

  if (!is.null(adapter_seq)) {
    adapter_seq <- .rad_scalar_character(adapter_seq, "adapter_seq")
  }
  if (!rescan && is.null(adapter_seq)) {
    stop("adapter_seq is required for sequence scans", call. = FALSE)
  }
  barcode_length <- .rad_scalar_integer(
    barcode_length, "barcode_length", minimum = 1L
  )
  if (barcode_length > 32L) {
    stop("barcode_length must be between 1 and 32", call. = FALSE)
  }
  left_margin <- .rad_scalar_integer(
    left_margin, "left_margin", minimum = 0L
  )
  right_margin <- .rad_scalar_integer(
    right_margin, "right_margin", minimum = 0L
  )
  max_reads <- .rad_scalar_integer(
    max_reads, "max_reads", minimum = 1L, allow_null = TRUE
  )
  max_error <- .rad_scan_scalar_double(
    max_error, "max_error", minimum = 0, maximum = 1
  )
  whitelist <- .rad_scan_whitelist_reference(whitelist, "whitelist")
  threads <- .rad_scalar_integer(threads, "threads", minimum = 1L)
  chunk_size <- .rad_scalar_integer(
    chunk_size, "chunk_size", minimum = 1L
  )
  if (!is.character(selection) || !length(selection) || anyNA(selection)) {
    stop("selection must be high_specificity or above_floor", call. = FALSE)
  }
  selection <- match.arg(
    tolower(selection), c("high_specificity", "above_floor")
  )
  bc1_whitelist <- .rad_scan_whitelist_reference(
    bc1_whitelist, "bc1_whitelist"
  )
  bc2_whitelist <- .rad_scan_whitelist_reference(
    bc2_whitelist, "bc2_whitelist"
  )
  split_mode <- !is.null(bc1_whitelist) || !is.null(bc2_whitelist)
  if (xor(is.null(bc1_whitelist), is.null(bc2_whitelist))) {
    stop("bc1_whitelist and bc2_whitelist must be supplied together",
         call. = FALSE)
  }
  if (split_mode && !is.null(whitelist)) {
    stop("whitelist and two-part barcode whitelists are alternatives",
         call. = FALSE)
  }
  if (rescan && split_mode) {
    stop("rescan does not support two-part barcode whitelists",
         call. = FALSE)
  }
  if (rescan && !is.null(whitelist)) {
    stop("rescan does not support a reference whitelist", call. = FALSE)
  }

  bc1_length <- .rad_scalar_integer(
    bc1_length, "bc1_length", minimum = 1L, allow_null = TRUE
  )
  bc2_length <- .rad_scalar_integer(
    bc2_length, "bc2_length", minimum = 1L, allow_null = TRUE
  )
  if (!is.null(bc1_length) && bc1_length > 32L) {
    stop("bc1_length must be between 1 and 32", call. = FALSE)
  }
  if (!is.null(bc2_length) && bc2_length > 32L) {
    stop("bc2_length must be between 1 and 32", call. = FALSE)
  }
  if (!split_mode && (!is.null(bc1_length) || !is.null(bc2_length))) {
    stop("bc1_length and bc2_length require two-part barcode whitelists",
         call. = FALSE)
  }
  umi_length <- .rad_scalar_integer(
    umi_length, "umi_length", minimum = 0L
  )
  offset_min <- .rad_scalar_integer(
    offset_min, "offset_min", minimum = 0L
  )
  offset_max <- .rad_scalar_integer(
    offset_max, "offset_max", minimum = 0L
  )
  if (offset_min > offset_max) {
    stop("offset_min must not exceed offset_max", call. = FALSE)
  }

  if (using_input) {
    input <- if (rescan) {
      .rad_existing_file(input, "input")
    } else {
      .rad_scan_sequence_file(input, "input")
    }
    if (is.null(output_prefix)) {
      output_prefix <- if (rescan) {
        .rad_scan_default_rescan_prefix(input)
      } else {
        "barcodes"
      }
    }
    output_prefix <- .rad_scan_output_prefix(output_prefix)
    batch <- NULL
    batch_csv <- NULL
  } else {
    if (!is.null(output_prefix)) {
      stop("output_prefix must be NULL when batch_csv supplies run prefixes",
           call. = FALSE)
    }
    batch_csv <- .rad_existing_file(batch_csv, "batch_csv")
    batch <- .rad_scan_parse_batch(
      batch_csv,
      split_mode = split_mode,
      protected_paths = c(whitelist, bc1_whitelist, bc2_whitelist)
    )
  }

  core <- rad_core_info()
  .rad_scan_require_feature(core)
  threads_requested <- threads
  threads_effective <- if (isTRUE(core$build$openmp)) threads else 1L
  if (threads_requested > 1L && threads_effective == 1L) {
    warning(
      "this rrad build has no OpenMP support; using one effective thread",
      call. = FALSE
    )
  }

  run_one <- function(run_input, run_prefix) {
    native_options <- .rad_drop_null(list(
      input = run_input,
      output_prefix = run_prefix,
      adapter_seq = adapter_seq,
      barcode_length = barcode_length,
      left_margin = left_margin,
      right_margin = right_margin,
      max_reads = max_reads,
      max_error = max_error,
      whitelist = whitelist,
      threads = threads_requested,
      chunk_size = chunk_size,
      selection = selection,
      rescan = rescan,
      bc1_whitelist = bc1_whitelist,
      bc2_whitelist = bc2_whitelist,
      bc1_length = bc1_length,
      bc2_length = bc2_length,
      umi_length = umi_length,
      offset_min = offset_min,
      offset_max = offset_max,
      verbose = verbose
    ))
    .rad_scan_run(native_options, core)
  }

  if (!using_batch) {
    return(run_one(input, output_prefix))
  }

  runs <- lapply(
    seq_len(nrow(batch)),
    function(index) {
      run_one(batch$input[[index]], batch$output_prefix[[index]])
    }
  )
  names(runs) <- batch$output_prefix
  effective_config <- runs[[1L]]$config[
    setdiff(names(runs[[1L]]$config), c("input", "output_prefix"))
  ]
  structure(
    list(
      status = "success",
      backend = "embedded",
      core = core,
      batch_csv = batch_csv,
      config = effective_config,
      runs = runs
    ),
    class = c("rad_scan_wl_batch_result", "list")
  )
}
