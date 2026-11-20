#' Identify reference catalogs or infer a read layout
#'
#' Uses a bounded prefix sample. Catalog selection counts exact read matches
#' independently for each reference and retains ties. Layout inference uses
#' known static parts and can return an incomplete structure.
#' @export
rad_discover <- function(input, out_dir = ".", output = "discovery",
                         mode = c("auto", "whitelist", "layout", "both"),
                         pattern = NULL, barcode_length = NULL, parts = NULL,
                         catalogs = NULL, catalog_lengths = NULL,
                         sample_reads = 10000L, sample_bases = 50000000L,
                         max_error = 0.15, threads = 1L, overwrite = FALSE) {
  input <- .rad_existing_file(input, "input")
  out_dir <- .rad_output_directory(out_dir)
  output <- .rad_output_prefix(output)
  mode <- match.arg(mode)
  pattern <- .rad_scalar_character(pattern, "pattern", allow_null = TRUE)
  if (mode == "auto") mode <- if (is.null(pattern)) "both" else "whitelist"
  parts <- .rad_optional_existing_file(parts, "parts")
  barcode_length <- .rad_scalar_integer(barcode_length, "barcode_length",
                                        minimum = 1L, allow_null = TRUE)
  if (!is.null(barcode_length) && barcode_length > 64L)
    stop("barcode_length must be between 1 and 64", call. = FALSE)
  if (mode != "whitelist" && !is.null(pattern))
    stop("pattern is supported only in whitelist mode", call. = FALSE)
  if (mode == "whitelist" && is.null(pattern))
    stop("whitelist mode requires pattern", call. = FALSE)
  if (!is.null(catalogs)) {
    if (!is.character(catalogs) || !length(catalogs) || anyNA(catalogs) ||
        any(!nzchar(catalogs)))
      stop("catalogs must be a nonempty character vector of file paths", call. = FALSE)
    catalog_names <- names(catalogs)
    if (is.null(catalog_names)) catalog_names <- basename(catalogs)
    if (anyNA(catalog_names) || any(!nzchar(catalog_names)))
      stop("catalog names must be nonempty", call. = FALSE)
    catalogs <- vapply(catalogs, .rad_existing_file, character(1), name = "catalogs")
    if (is.null(catalog_lengths)) catalog_lengths <- rep(0L, length(catalogs))
    if (!is.numeric(catalog_lengths) || length(catalog_lengths) != length(catalogs) ||
        anyNA(catalog_lengths) || any(!is.finite(catalog_lengths)) ||
        any(catalog_lengths != floor(catalog_lengths)) ||
        any(catalog_lengths < 0 | catalog_lengths > 64))
      stop("catalog_lengths must give one integer from 0 to 64 per catalog", call. = FALSE)
    catalog_lengths <- as.integer(catalog_lengths)
  } else {
    if (!is.null(catalog_lengths)) stop("catalog_lengths requires catalogs", call. = FALSE)
    catalog_names <- NULL
  }
  sample_reads <- .rad_scalar_integer(sample_reads, "sample_reads", minimum = 1L)
  sample_bases <- .rad_scalar_integer(sample_bases, "sample_bases", minimum = 1L)
  max_error <- .rad_scalar_double(max_error, "max_error", minimum = 0)
  if (max_error > 1) stop("max_error must be between 0 and 1", call. = FALSE)
  threads <- .rad_scalar_integer(threads, "threads", minimum = 1L)
  overwrite <- .rad_scalar_logical(overwrite, "overwrite")
  options <- .rad_drop_null(list(
    input = input, output_prefix = file.path(out_dir, output), pattern = pattern,
    barcode_length = barcode_length, parts = parts, catalogs = catalogs,
    catalog_names = catalog_names, catalog_lengths = catalog_lengths,
    sample_reads = sample_reads, sample_bases = sample_bases,
    max_error = max_error, threads = threads, overwrite = overwrite,
    layout = mode %in% c("layout", "both"),
    whitelist = mode %in% c("whitelist", "both")))
  value <- .rad_call_with_resource_retry(function() .rad_discover_cpp(options))
  value$config <- list(mode = mode, pattern = pattern, barcode_length = barcode_length,
                       sampling = "prefix", sample_reads = sample_reads,
                       sample_bases = sample_bases, catalogs = catalogs)
  # A no-hit, tied, or partial result is successful execution with incomplete
  # inference. Keep that distinction visible instead of raising an error.
  for (entry in c("whitelists.csv", "whitelist_match.csv", "parts.csv")) {
    if (!is.null(value$files[[entry]]))
      value[[sub("\\.csv$", "", entry)]] <- utils::read.csv(
        value$files[[entry]], stringsAsFactors = FALSE, check.names = FALSE)
  }
  class(value) <- c("rad_discovery_result", "list")
  value
}

#' Export demultiplexed reads with explicit spatial coordinates
#' @export
rad_spatial_reformat <- function(fastq, out_dir, sample, map = NULL,
                                 bc1_axis = NULL, bc2_axis = NULL,
                                 metadata = NULL, unit = c("bin", "cell"),
                                 bin_size = 2L, tissue_only = FALSE,
                                 strict = FALSE) {
  fastq <- .rad_existing_file(fastq, "fastq")
  out_dir <- .rad_scalar_character(out_dir, "out_dir")
  if (!nzchar(out_dir) || file.exists(out_dir) || dir.exists(out_dir))
    stop("out_dir must name a new directory", call. = FALSE)
  sample <- .rad_scalar_character(sample, "sample")
  if (!nzchar(sample) || grepl("[[:space:]:]", sample))
    stop("sample must be nonempty and contain no whitespace or colon", call. = FALSE)
  map <- .rad_optional_existing_file(map, "map")
  bc1_axis <- .rad_optional_existing_file(bc1_axis, "bc1_axis")
  bc2_axis <- .rad_optional_existing_file(bc2_axis, "bc2_axis")
  metadata <- .rad_optional_existing_file(metadata, "metadata")
  if (is.null(map)) {
    if (is.null(bc1_axis) || is.null(bc2_axis))
      stop("supply map or both bc1_axis and bc2_axis", call. = FALSE)
  } else if (!is.null(bc1_axis) || !is.null(bc2_axis)) {
    stop("map and barcode axes cannot be combined", call. = FALSE)
  }
  unit <- match.arg(unit)
  bin_size <- .rad_scalar_integer(bin_size, "bin_size", minimum = 2L)
  if (bin_size %% 2L) stop("bin_size must be an even number of microns", call. = FALSE)
  if (unit == "cell" && bin_size != 2L)
    stop("cell assignment requires bin_size = 2", call. = FALSE)
  tissue_only <- .rad_scalar_logical(tissue_only, "tissue_only")
  strict <- .rad_scalar_logical(strict, "strict")
  value <- .rad_spatial_reformat_cpp(.rad_drop_null(list(
    input = fastq, output_dir = out_dir, sample = sample, map = map,
    bc1_axis = bc1_axis, bc2_axis = bc2_axis, metadata = metadata,
    unit = unit, bin_size = bin_size, tissue_only = tissue_only, strict = strict)))
  class(value) <- c("rad_spatial_result", "list")
  value
}
