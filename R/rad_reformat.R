# Native RAD output reformatting ---------------------------------------------

.rad_reformat_input <- function(x) {
  if (inherits(x, "rad_demux_result")) {
    candidate <- x$artifacts$fastq
    if (!is.character(candidate) || length(candidate) != 1L ||
        is.na(candidate) || !nzchar(candidate)) {
      stop(
        "the rad_demux_result does not contain one aggregate FASTQ artifact",
        call. = FALSE
      )
    }
    x <- candidate
  }
  .rad_existing_file(x, "fastq")
}

.rad_reformat_output_file <- function(x) {
  x <- .rad_scalar_character(x, "output_fastq")
  expanded <- path.expand(x)
  if (grepl("[/\\\\]$", expanded)) {
    stop("output_fastq must end in a filename", call. = FALSE)
  }
  if (file.exists(expanded) && dir.exists(expanded)) {
    stop("output_fastq identifies a directory: ", expanded, call. = FALSE)
  }
  parent <- dirname(expanded)
  if (file.exists(parent) && !dir.exists(parent)) {
    stop("the output_fastq parent identifies a file: ", parent,
         call. = FALSE)
  }
  if (!dir.exists(parent) &&
      !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("could not create the output_fastq directory: ", parent,
         call. = FALSE)
  }
  file.path(normalizePath(parent, mustWork = TRUE), basename(expanded))
}

.rad_reformat_require_fasta_suffix <- function(x) {
  if (!grepl("\\.(fa|fasta)(\\.gz)?$", x, ignore.case = TRUE)) {
    stop(
      "output_fastq must end in .fa, .fasta, .fa.gz, or .fasta.gz ",
      "when to_fasta is TRUE",
      call. = FALSE
    )
  }
  x
}

.rad_reformat_delimiter <- function(x) {
  x <- .rad_scalar_character(x, "delimiter")
  if (nchar(x, type = "bytes") != 1L) {
    stop("delimiter must be exactly one single-byte character", call. = FALSE)
  }
  if (grepl("^[[:space:]]$", x)) {
    stop("delimiter cannot be whitespace", call. = FALSE)
  }
  x
}

.rad_reformat_header_format <- function(header_format, reformat_header) {
  if (is.null(header_format)) {
    return(if (reformat_header) "collapsed" else "preserve")
  }
  header_format <- .rad_scalar_character(header_format, "header_format")
  header_format <- tolower(header_format)
  choices <- c("preserve", "collapsed", "presto")
  if (!(header_format %in% choices)) {
    stop(
      "header_format must be one of preserve, collapsed, or presto",
      call. = FALSE
    )
  }
  if (reformat_header && !identical(header_format, "collapsed")) {
    stop(
      "reformat_header = TRUE conflicts with header_format = \"",
      header_format, "\"; use header_format = \"collapsed\" or omit it",
      call. = FALSE
    )
  }
  header_format
}

.rad_reformat_require_feature <- function(core, to_fasta = FALSE,
                                          header_format = "preserve") {
  if (!is.character(core$features) || anyNA(core$features) ||
      !("reformat" %in% core$features)) {
    stop("the embedded RAD core does not advertise the reformat feature",
         call. = FALSE)
  }
  if (to_fasta && !("reformat-fasta" %in% core$features)) {
    stop(
      "the embedded RAD core does not advertise FASTA reformatting",
      call. = FALSE
    )
  }
  if (identical(header_format, "presto") &&
      !("reformat-presto-header" %in% core$features)) {
    stop(
      "the embedded RAD core does not advertise pRESTO header formatting",
      call. = FALSE
    )
  }
}

.rad_reformat_validate_result <- function(native, config, core) {
  required <- c(
    "success", "backend", "core", "input", "output", "files", "stats",
    "log"
  )
  if (!is.list(native) || is.null(names(native)) ||
      length(setdiff(required, names(native)))) {
    stop("the embedded RAD core returned an invalid reformat result",
         call. = FALSE)
  }
  if (!isTRUE(native$success)) {
    stop("the embedded RAD core reported an unsuccessful reformat run",
         call. = FALSE)
  }
  if (!is.character(native$backend) || length(native$backend) != 1L ||
      is.na(native$backend) || !identical(native$backend, "embedded")) {
    stop("the embedded RAD core returned an invalid reformat backend",
         call. = FALSE)
  }
  if (!is.character(native$input) || length(native$input) != 1L ||
      is.na(native$input) || !nzchar(native$input) ||
      !is.character(native$output) || length(native$output) != 1L ||
      is.na(native$output) || !nzchar(native$output)) {
    stop("the embedded RAD core returned invalid reformat path metadata",
         call. = FALSE)
  }
  if (!is.list(native$files) || is.null(names(native$files)) ||
      !is.list(native$stats) || is.null(names(native$stats)) ||
      !is.character(native$log) || length(native$log) != 1L ||
      is.na(native$log)) {
    stop(
      "the embedded RAD core returned invalid reformat artifacts or statistics",
      call. = FALSE
    )
  }

  native_core <- .rad_validate_core_info(native$core)
  handshake_fields <- c(
    "core_version", "core_api_version", "source_commit",
    "embedded_source_digest", "resource_schema", "resource_digest"
  )
  if (!identical(unclass(core[handshake_fields]),
                 unclass(native_core[handshake_fields]))) {
    stop("the embedded RAD core identity changed during the reformat run",
         call. = FALSE)
  }

  result <- structure(
    list(
      status = "success",
      backend = native$backend,
      core = native_core,
      input = native$input,
      output = native$output,
      config = config,
      artifacts = native$files,
      stats = native$stats,
      log = native$log
    ),
    class = c("rad_reformat_result", "list")
  )
  if (isTRUE(config$verbose) && nzchar(native$log)) {
    cat(native$log)
    if (!grepl("\n$", native$log)) cat("\n")
  }
  result
}

#' Reformat or split RAD-tagged sequence files with the embedded core
#'
#' Streams a FASTQ or FASTA file through RAD's native post-processing engine.
#' It can collapse `CB:Z` and `UB:Z` tags into the read name, split records into
#' barcode-specific gzip files, emit canonical pRESTO annotations for IgBLAST,
#' convert FASTQ to FASTA, or translate Visium HD barcode pairs into spatial
#' coordinate headers. No external RAD or compression executable is launched.
#'
#' @param fastq Input FASTQ, FASTA, or gzip-compressed FASTQ/FASTA path, or a
#'   `rad_demux_result` whose aggregate `artifacts$fastq` should be consumed.
#' @param out_dir Output directory for barcode-specific files. Required when
#'   `split_bc = TRUE` and unused otherwise.
#' @param split_bc Write one gzip file per observed `CB:Z` barcode, using
#'   `.fq.gz` for FASTQ input and `.fa.gz` for FASTA input.
#' @param reformat_header Collapse the read name, cell barcode, and UMI to
#'   `QNAME<delimiter>CB<delimiter>UB`. For RAD CLI parity, the complete input
#'   comment is removed, including fields other than the `CB` and `UB` tags.
#'   This is the legacy equivalent of `header_format = "collapsed"` and cannot
#'   be combined with a different explicit `header_format` or with `coordinate`.
#' @param collapsed_input Parse tagless input names as already-collapsed
#'   `QNAME<delimiter>CB<delimiter>UB` values. The default is `FALSE` so normal
#'   read names containing delimiters cannot be mistaken for cell barcodes.
#' @param delimiter A single-byte delimiter used for collapsed headers and for
#'   recognizing already-collapsed headers when `collapsed_input = TRUE`.
#' @param coordinate Optional coordinate mode. `NULL` disables coordinate
#'   rewriting; `"vizHD-v1"` maps bundled Visium HD BC1/BC2 barcodes to spatial
#'   bin identifiers. Coordinate mode preserves unrelated comment fields, is
#'   mutually exclusive with transforming header formats, and counts as a
#'   requested action.
#' @param bin_size Positive spatial bin size in microns.
#' @param output_fastq Optional aggregate output path for a non-splitting run.
#'   When omitted, a header or coordinate rewrite replaces `fastq` atomically
#'   in place. Despite the legacy argument name, this is also the FASTA output
#'   path when `to_fasta = TRUE`. It cannot be supplied with
#'   `split_bc = TRUE`.
#' @param threads Positive requested worker-thread count. A build without
#'   OpenMP warns and uses one effective thread.
#' @param chunk_size Positive number of records per streaming chunk.
#' @param verbose Print native progress captured during the completed run.
#' @param to_fasta Convert FASTQ input records to FASTA by omitting the plus
#'   and quality lines, like `seqkit fq2fa`. The read name and complete comment,
#'   including RAD's SAM-style tags, are preserved unless another requested
#'   header transformation changes them. Aggregate conversion requires an
#'   explicit `output_fastq` path; split conversion writes `.fa.gz` files.
#' @param header_format Optional header transformation. `NULL` preserves legacy
#'   behavior: `reformat_header = TRUE` selects `"collapsed"`; otherwise
#'   `"preserve"` is used. `"preserve"` keeps the QNAME and full comment,
#'   `"collapsed"` emits the positional legacy form, and `"presto"` emits the
#'   whitespace-free pRESTO annotation form
#'   `QNAME|BARCODE=BC_OR_CB|UMI=UB`. `BARCODE` uses a complete, nonempty
#'   SAM-style `BC:Z:` spatial barcode when present and otherwise falls back to
#'   `CB:Z:`; `UMI` uses `UB:Z:`. Fields may occur in any order, and missing
#'   fields are omitted. QNAME and emitted barcode/UMI values containing
#'   whitespace or the reserved pRESTO delimiters `|`, `=`, or `,` are rejected.
#' @return A `rad_reformat_result` containing provenance, normalized config,
#'   exact output paths, processing statistics, and captured log output.
#' @details
#' At least one of `split_bc`, a transforming `header_format`, `coordinate`, or
#' `to_fasta` must request work. Split-only runs leave the input untouched. A
#' non-splitting run with no `output_fastq` performs a checked temporary write
#' followed by atomic replacement; the original is retained if parsing or
#' writing fails. Legal zero-length FASTQ/FASTA records and FASTQ plus-line
#' identifiers or comments are preserved unless FASTA output intentionally
#' removes FASTQ-only fields.
#'
#' The pRESTO format always writes QNAME first, followed by `BARCODE` and `UMI`
#' annotations in that canonical order when their source tags are present.
#' Other SAM-style fields are intentionally not copied into the whitespace-free
#' identifier. Barcode-specific splitting remains keyed by `CB:Z:` even when a
#' distinct `BC:Z:` value is projected into the pRESTO `BARCODE` annotation.
#' This makes the result directly usable as a stable IgBLAST query identifier
#' while leaving downstream code to reject reads without a barcode.
#' @name rad-reformat
#' @rdname rad-reformat
rad_reformat <- function(fastq,
                         out_dir = NULL,
                         split_bc = FALSE,
                         reformat_header = FALSE,
                         collapsed_input = FALSE,
                         delimiter = "_",
                         coordinate = NULL,
                         bin_size = 2L,
                         output_fastq = NULL,
                         threads = 1L,
                         chunk_size = 5000L,
                         verbose = FALSE,
                         to_fasta = FALSE,
                         header_format = NULL) {
  fastq <- .rad_reformat_input(fastq)
  split_bc <- .rad_scalar_logical(split_bc, "split_bc")
  reformat_header <- .rad_scalar_logical(
    reformat_header, "reformat_header"
  )
  collapsed_input <- .rad_scalar_logical(
    collapsed_input, "collapsed_input"
  )
  verbose <- .rad_scalar_logical(verbose, "verbose")
  to_fasta <- .rad_scalar_logical(to_fasta, "to_fasta")
  header_format <- .rad_reformat_header_format(
    header_format, reformat_header
  )
  reformat_header <- identical(header_format, "collapsed")
  delimiter <- .rad_reformat_delimiter(delimiter)
  if (!is.null(coordinate)) {
    coordinate <- .rad_scalar_character(coordinate, "coordinate")
  }
  if (!identical(header_format, "preserve") && !is.null(coordinate)) {
    if (identical(header_format, "collapsed")) {
      stop("reformat_header cannot be combined with coordinate", call. = FALSE)
    }
    stop("header_format cannot be combined with coordinate", call. = FALSE)
  }
  if (identical(header_format, "presto") && collapsed_input) {
    stop(
      "header_format = \"presto\" requires actual SAM tags and cannot ",
      "be combined with collapsed_input = TRUE",
      call. = FALSE
    )
  }
  bin_size <- .rad_scalar_integer(bin_size, "bin_size", minimum = 1L)
  threads <- .rad_scalar_integer(threads, "threads", minimum = 1L)
  chunk_size <- .rad_scalar_integer(
    chunk_size, "chunk_size", minimum = 1L
  )

  if (!split_bc && identical(header_format, "preserve") &&
      is.null(coordinate) && !to_fasta) {
    stop(
      paste0(
        "request at least one of split_bc, a transforming header_format, ",
        "coordinate, or to_fasta"
      ),
      call. = FALSE
    )
  }
  if (split_bc) {
    if (is.null(out_dir)) {
      stop("out_dir is required when split_bc is TRUE", call. = FALSE)
    }
    if (!is.null(output_fastq)) {
      stop("output_fastq cannot be supplied when split_bc is TRUE",
           call. = FALSE)
    }
    out_dir <- .rad_output_directory(out_dir)
  } else {
    if (!is.null(out_dir)) {
      stop("out_dir is only used when split_bc is TRUE", call. = FALSE)
    }
    output_fastq <- if (is.null(output_fastq)) {
      NULL
    } else {
      .rad_reformat_output_file(output_fastq)
    }
    if (to_fasta && is.null(output_fastq)) {
      stop(
        "output_fastq is required for aggregate FASTA conversion",
        call. = FALSE
      )
    }
    if (to_fasta) {
      output_fastq <- .rad_reformat_require_fasta_suffix(output_fastq)
    }
  }

  core <- rad_core_info()
  .rad_reformat_require_feature(
    core, to_fasta = to_fasta, header_format = header_format
  )
  threads_requested <- threads
  threads_effective <- if (isTRUE(core$build$openmp)) threads else 1L
  if (threads_requested > 1L && threads_effective == 1L) {
    warning(
      "this rrad build has no OpenMP support; using one effective thread",
      call. = FALSE
    )
  }

  config <- list(
    fastq = fastq,
    out_dir = out_dir,
    split_bc = split_bc,
    reformat_header = reformat_header,
    header_format = header_format,
    collapsed_input = collapsed_input,
    delimiter = delimiter,
    coordinate = coordinate,
    bin_size = bin_size,
    output_fastq = output_fastq,
    threads = threads_effective,
    threads_requested = threads_requested,
    threads_effective = threads_effective,
    chunk_size = chunk_size,
    verbose = verbose,
    to_fasta = to_fasta
  )
  native_options <- .rad_drop_null(list(
    input_path = fastq,
    split_output_dir = out_dir,
    split_by_barcode = split_bc,
    reformat_header = reformat_header,
    presto_header = identical(header_format, "presto"),
    parse_collapsed_id = collapsed_input,
    delimiter = delimiter,
    coordinate_mode = coordinate,
    bin_size_um = bin_size,
    output_path = output_fastq,
    threads = threads_requested,
    chunk_size = chunk_size,
    verbose = verbose,
    to_fasta = to_fasta
  ))
  native <- .rad_reformat_cpp(native_options)
  .rad_reformat_validate_result(native, config, core)
}
