# Handle-backed read storage -------------------------------------------------

.as_rad_reads <- function(x) {
  data.table::setDT(x)
  data.table::setattr(x, "class", c("rad_reads", "data.table", "data.frame"))
  x
}

.rad_call_name <- function(expression) {
  if (!is.call(expression)) return(NA_character_)
  head <- expression[[1L]]
  if (is.symbol(head)) return(as.character(head))
  if (is.call(head) && as.character(head[[1L]]) %in% c("::", ":::")) {
    return(as.character(head[[3L]]))
  }
  NA_character_
}

.rad_one_byref_assignment <- function(expression, environment) {
  arguments <- as.list(expression)[-1L]
  argument_names <- names(arguments)
  named <- if (is.null(argument_names)) {
    rep(FALSE, length(arguments))
  } else {
    nzchar(argument_names)
  }
  if (length(arguments) && all(named)) return(argument_names)
  if (length(arguments) != 2L || any(named)) return(NA_character_)

  lhs <- arguments[[1L]]
  if (is.symbol(lhs)) return(as.character(lhs))
  if (is.character(lhs) && !anyNA(lhs) && all(nzchar(lhs))) return(lhs)

  # Indirect code (for example eval() or parse()) cannot be rewritten before
  # data.table executes it. Never evaluate a computed target merely to inspect
  # it: a stateful expression could return a different column on its next run.
  NA_character_
}

.rad_freeze_one_byref_assignment <- function(expression, environment) {
  arguments <- as.list(expression)[-1L]
  argument_names <- names(arguments)
  named <- if (is.null(argument_names)) {
    rep(FALSE, length(arguments))
  } else {
    nzchar(argument_names)
  }
  if (length(arguments) && all(named)) {
    return(list(expression = expression, targets = argument_names))
  }
  if (length(arguments) != 2L || any(named)) {
    return(list(expression = expression, targets = NA_character_))
  }

  lhs <- arguments[[1L]]
  if (is.symbol(lhs)) {
    return(list(expression = expression, targets = as.character(lhs)))
  }
  if (is.character(lhs) && !anyNA(lhs) && all(nzchar(lhs))) {
    return(list(expression = expression, targets = lhs))
  }

  # data.table's supported programmatic assignment spelling is `(columns) :=`.
  # Resolve a simple column-vector binding exactly once, then put that literal
  # value into the expression that data.table will execute. Arbitrary calls are
  # rejected without evaluation, so their result cannot change after preflight.
  if (is.call(lhs) && identical(.rad_call_name(lhs), "(") &&
      length(lhs) == 2L && is.symbol(lhs[[2L]])) {
    binding <- as.character(lhs[[2L]])
    targets <- tryCatch(
      get(binding, envir = environment, inherits = TRUE),
      error = function(...) NULL
    )
    if (is.character(targets) && !anyNA(targets) && all(nzchar(targets))) {
      expression[[2L]] <- targets
      return(list(expression = expression, targets = targets))
    }
  }

  list(expression = expression, targets = NA_character_)
}

.rad_freeze_byref_assignments <- function(expression, environment, depth = 0L) {
  if (depth > 16L) {
    return(list(expression = expression, targets = NA_character_))
  }
  if (!is.call(expression)) {
    return(list(expression = expression, targets = character()))
  }

  name <- .rad_call_name(expression)
  if (name %in% c(":=", "let")) {
    return(.rad_freeze_one_byref_assignment(expression, environment))
  }
  if (name %in% c("set", "setattr", "setnames")) {
    return(list(expression = expression, targets = NA_character_))
  }
  if (name %in% c("eval", "evalq", "eval.parent", "parse")) {
    return(list(
      expression = expression,
      targets = .rad_byref_assignment_targets(expression, environment, depth)
    ))
  }

  targets <- character()
  if (length(expression) >= 2L) {
    for (position in seq.int(2L, length(expression))) {
      frozen <- .rad_freeze_byref_assignments(
        expression[[position]], environment, depth + 1L
      )
      expression[[position]] <- frozen$expression
      targets <- c(targets, frozen$targets)
    }
  }
  list(expression = expression, targets = unique(targets))
}

.rad_targets_from_code_value <- function(value, environment, depth) {
  if (is.call(value)) {
    return(.rad_byref_assignment_targets(value, environment, depth + 1L))
  }
  if (is.expression(value)) {
    nested <- lapply(as.list(value), function(element) {
      .rad_targets_from_code_value(element, environment, depth + 1L)
    })
    return(unique(unlist(nested, use.names = FALSE)))
  }
  character()
}

.rad_static_parse_value <- function(expression, environment) {
  arguments <- as.list(expression)[-1L]
  argument_names <- names(arguments)
  text_position <- if (!is.null(argument_names) && "text" %in% argument_names) {
    match("text", argument_names)
  } else if (length(arguments)) {
    1L
  } else {
    NA_integer_
  }
  if (is.na(text_position)) return(list(known = FALSE, value = NULL))

  text_expression <- arguments[[text_position]]
  text <- if (is.character(text_expression)) {
    text_expression
  } else if (is.symbol(text_expression) &&
             exists(as.character(text_expression), envir = environment,
                    inherits = TRUE)) {
    tryCatch(get(as.character(text_expression), envir = environment,
                 inherits = TRUE), error = function(...) NULL)
  } else {
    NULL
  }
  if (!is.character(text)) return(list(known = FALSE, value = NULL))
  parsed <- tryCatch(parse(text = text, keep.source = FALSE),
                     error = function(...) NULL)
  list(known = !is.null(parsed), value = parsed)
}

.rad_static_eval_value <- function(expression, environment) {
  if (is.symbol(expression)) {
    name <- as.character(expression)
    if (!exists(name, envir = environment, inherits = TRUE)) {
      # A symbol resolved only inside data.table's column environment evaluates
      # to column data, not to another expression that can perform assignment.
      return(list(known = TRUE, value = NULL))
    }
    value <- tryCatch(get(name, envir = environment, inherits = TRUE),
                      error = function(...) NULL)
    return(list(known = !is.null(value), value = value))
  }
  if (!is.call(expression)) return(list(known = TRUE, value = expression))

  name <- .rad_call_name(expression)
  if (identical(name, "quote") && length(expression) == 2L) {
    return(list(known = TRUE, value = expression[[2L]]))
  }
  if (identical(name, "parse")) return(.rad_static_parse_value(expression, environment))
  list(known = FALSE, value = NULL)
}

.rad_byref_assignment_targets <- function(expression, environment, depth = 0L) {
  if (depth > 16L) return(NA_character_)
  if (!is.call(expression)) return(character())
  name <- .rad_call_name(expression)
  if (name %in% c(":=", "let")) {
    return(.rad_one_byref_assignment(expression, environment))
  }
  if (name %in% c("set", "setattr", "setnames")) {
    # These low-level data.table mutators intentionally bypass S3 replacement
    # methods. They remain an explicitly documented escape hatch when called
    # directly, but must not hide inside j where an error after evaluation
    # would leave the caller's table modified.
    return(NA_character_)
  }
  if (name %in% c("eval", "evalq", "eval.parent")) {
    if (length(expression) < 2L) return(NA_character_)
    resolved <- .rad_static_eval_value(expression[[2L]], environment)
    if (!isTRUE(resolved$known)) return(NA_character_)
    return(.rad_targets_from_code_value(resolved$value, environment, depth + 1L))
  }
  if (identical(name, "parse")) {
    resolved <- .rad_static_parse_value(expression, environment)
    if (!isTRUE(resolved$known)) return(NA_character_)
    return(.rad_targets_from_code_value(resolved$value, environment, depth + 1L))
  }

  nested <- lapply(as.list(expression)[-1L], function(argument) {
    .rad_byref_assignment_targets(argument, environment, depth + 1L)
  })
  unique(unlist(nested, use.names = FALSE))
}

.rad_plain_data_table <- function(x) {
  data.table::setattr(x, "class", setdiff(class(x), "rad_reads"))
  data.table::setattr(x, ".rad_store", NULL)
  x
}

.rad_finalize_table_result <- function(x) {
  if (!inherits(x, "rad_reads") || !data.table::is.data.table(x)) return(x)
  required <- c(
    "read_key", "seq_id", "seq", "qual", "length",
    "mean_q", "min_q", "frac_q20"
  )
  if (!all(required %in% names(x))) return(.rad_plain_data_table(x))

  handle_classes <- c(
    inherits(x$seq_id, "rad_id"),
    inherits(x$seq, "rad_seq"),
    inherits(x$qual, "rad_phred")
  )
  if (!all(handle_classes)) return(.rad_plain_data_table(x))

  # A data.table join can source its shared join column from i. When i is an
  # independently packed table that leaves a superficially complete result
  # whose ID handles no longer describe the sequence/quality records. Refuse
  # that result immediately instead of allowing a stale rad_reads object to
  # escape and fail only when it is eventually materialized.
  rad_validate(x)
  x
}

.rad_core_state_matches <- function(current, snapshot) {
  if (!inherits(current, "rad_reads") || !data.table::is.data.table(current)) {
    return(FALSE)
  }
  current_core <- intersect(.rad_core_columns, names(current))
  snapshot_core <- intersect(.rad_core_columns, names(snapshot))
  if (!identical(current_core, snapshot_core) || nrow(current) != nrow(snapshot)) {
    return(FALSE)
  }
  if (!isTRUE(tryCatch(rad_validate(current), error = function(...) FALSE))) {
    return(FALSE)
  }

  current_order <- order(current$read_key, na.last = TRUE)
  snapshot_order <- order(snapshot$read_key, na.last = TRUE)
  all(vapply(snapshot_core, function(column) {
    identical(current[[column]][current_order],
              snapshot[[column]][snapshot_order])
  }, logical(1)))
}

.rad_restore_table <- function(current, snapshot) {
  # Restore through data.table's by-reference primitives so the caller's object,
  # not merely this local binding, is repaired after a failing j expression.
  data.table::setattr(current, "class", c("data.table", "data.frame"))
  extras <- setdiff(names(current), names(snapshot))
  if (length(extras)) data.table::set(current, j = extras, value = NULL)
  for (column in names(snapshot)) {
    data.table::set(current, j = column, value = snapshot[[column]])
  }
  data.table::setcolorder(current, names(snapshot))

  structural <- c("names", ".internal.selfref")
  current_attributes <- setdiff(names(attributes(current)), structural)
  snapshot_attributes <- setdiff(names(attributes(snapshot)), structural)
  for (attribute in setdiff(current_attributes, snapshot_attributes)) {
    data.table::setattr(current, attribute, NULL)
  }
  for (attribute in snapshot_attributes) {
    data.table::setattr(
      current, attribute, attr(snapshot, attribute, exact = TRUE)
    )
  }
  invisible(current)
}

`[.rad_reads` <- function(x, i, j, ..., env = NULL) {
  snapshot <- NULL
  if (!missing(j)) {
    snapshot <- data.table::copy(x)
    expression <- substitute(j)
    if (!is.null(env)) {
      expression <- tryCatch(
        do.call(data.table::substitute2, list(expression, env)),
        error = function(...) NULL
      )
      if (is.null(expression)) {
        stop(
          "cannot verify an env-substituted operation on rad_reads; ",
          "materialize or repack the table first",
          call. = FALSE
        )
      }
    }
    frozen <- .rad_freeze_byref_assignments(expression, parent.frame())
    expression <- frozen$expression
    targets <- frozen$targets
    if (anyNA(targets)) {
      stop(
        "cannot verify a dynamic by-reference assignment on rad_reads; ",
        "materialize or repack the table first",
        call. = FALSE
      )
    }
    protected <- intersect(targets, .rad_core_columns)
    if (length(protected)) {
      stop(
        "packed RAD core columns are immutable; materialize or repack before ",
        "replacing: ", paste(protected, collapse = ", "),
        call. = FALSE
      )
    }
  }
  call <- match.call()
  if (!missing(j)) call$j <- expression
  if (!is.null(env)) call$env <- NULL
  call[[1L]] <- utils::getS3method("[", "data.table")
  result <- tryCatch(
    eval(call, envir = parent.frame()),
    error = function(condition) {
      if (!is.null(snapshot)) .rad_restore_table(x, snapshot)
      stop(condition)
    }
  )
  if (!is.null(snapshot) && !.rad_core_state_matches(x, snapshot)) {
    .rad_restore_table(x, snapshot)
    stop(
      "a data.table j expression attempted to mutate packed RAD core columns; ",
      "all changes were rolled back",
      call. = FALSE
    )
  }
  tryCatch(
    .rad_finalize_table_result(result),
    error = function(condition) {
      if (!is.null(snapshot)) .rad_restore_table(x, snapshot)
      stop(condition)
    }
  )
}

.rad_stop_core_replacement <- function(targets) {
  protected <- intersect(targets, .rad_core_columns)
  if (length(protected)) {
    stop(
      "packed RAD core columns are immutable; materialize or repack before ",
      "replacing: ", paste(protected, collapse = ", "),
      call. = FALSE
    )
  }
}

`$<-.rad_reads` <- function(x, name, value) {
  .rad_stop_core_replacement(name)
  out <- data.table::copy(x)
  data.table::set(out, j = name, value = value)
  .rad_finalize_table_result(out)
}

`[[<-.rad_reads` <- function(x, i, ..., value) {
  if (length(i) != 1L || is.na(i)) {
    stop("rad_reads column replacement requires one column", call. = FALSE)
  }
  if (is.character(i)) {
    target <- i
  } else if (is.numeric(i) && i == trunc(i) && i >= 1L && i <= ncol(x)) {
    target <- names(x)[[as.integer(i)]]
  } else {
    stop("rad_reads column replacement requires an existing position or a name",
         call. = FALSE)
  }
  .rad_stop_core_replacement(target)
  out <- data.table::copy(x)
  data.table::set(out, j = target, value = value)
  .rad_finalize_table_result(out)
}

`[<-.rad_reads` <- function(x, i, j, value) {
  if (missing(j)) {
    targets <- names(x)
  } else {
    targets <- tryCatch(eval(substitute(j), envir = parent.frame()),
                        error = function(...) NULL)
    if (is.numeric(targets) && all(is.finite(targets)) &&
        all(targets == trunc(targets)) && all(targets >= 1L) &&
        all(targets <= ncol(x))) {
      targets <- names(x)[as.integer(targets)]
    }
    if (!is.character(targets) || anyNA(targets) || any(!nzchar(targets))) {
      stop(
        "cannot verify a dynamic replacement on rad_reads; materialize or ",
        "repack the table first",
        call. = FALSE
      )
    }
  }
  .rad_stop_core_replacement(targets)
  call <- match.call()
  if (!missing(j)) call$j <- targets
  call[[1L]] <- utils::getS3method("[<-", "data.table")
  .rad_finalize_table_result(eval(call, envir = parent.frame()))
}

.rad_core_columns <- c(
  "read_key", "seq_id", "seq_comment", "seq", "qual",
  "length", "mean_q", "min_q", "frac_q20"
)

.rad_has_native_reference <- function(value) {
  if (inherits(value, c("rad_id", "rad_seq", "rad_phred"))) return(TRUE)
  if (typeof(value) %in% c("externalptr", "weakref", "environment")) return(TRUE)
  if (is.function(value)) return(TRUE)
  if (is.list(value) &&
      any(vapply(value, .rad_has_native_reference, logical(1)))) return(TRUE)
  attrs <- attributes(value)
  !is.null(attrs) && any(vapply(attrs, .rad_has_native_reference, logical(1)))
}

.rad_assert_safe_extras <- function(x, action) {
  extras <- setdiff(names(x), .rad_core_columns)
  unsafe <- extras[vapply(extras, function(column) {
    .rad_has_native_reference(x[[column]])
  }, logical(1))]
  if (length(unsafe)) {
    stop(
      "additional columns cannot be ", action,
      " because they contain packed handles or native references: ",
      paste(unsafe, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rad_handle_subset <- function(x, i) {
  if (missing(i)) return(x)
  old_class <- class(x)
  store <- attr(x, ".rad_store", exact = TRUE)
  field <- attr(x, ".rad_field", exact = TRUE)
  # Delegate insertion of missing/out-of-bounds values to bit64 so it writes
  # INT64_MIN, not an IEEE NA_real_ bit pattern that would look like a handle.
  class(x) <- "integer64"
  out <- x[i]
  class(out) <- old_class
  attr(out, ".rad_store") <- store
  if (!is.null(field)) attr(out, ".rad_field") <- field
  out
}

.rad_handle_c <- function(..., recursive = FALSE) {
  values <- list(...)
  prototype <- values[[1L]]
  handle_class <- class(prototype)[[1L]]
  store <- attr(prototype, ".rad_store", exact = TRUE)
  field <- attr(prototype, ".rad_field", exact = TRUE)

  compatible <- vapply(values, function(value) {
    inherits(value, handle_class) &&
      identical(attr(value, ".rad_store", exact = TRUE), store) &&
      identical(attr(value, ".rad_field", exact = TRUE), field)
  }, logical(1))
  if (is.null(store) || !all(compatible)) {
    stop("packed RAD handles can only be combined from the same store and field",
         call. = FALSE)
  }

  integer64_values <- lapply(values, function(value) {
    class(value) <- "integer64"
    value
  })
  out <- do.call(bit64::c.integer64,
                 c(integer64_values, list(recursive = recursive)))
  class(out) <- class(prototype)
  attr(out, ".rad_store") <- store
  if (!is.null(field)) attr(out, ".rad_field") <- field
  out
}

.rad_handle_rep <- function(x, ...) {
  old_class <- class(x)
  store <- attr(x, ".rad_store", exact = TRUE)
  field <- attr(x, ".rad_field", exact = TRUE)
  class(x) <- "integer64"
  out <- rep(x, ...)
  class(out) <- old_class
  attr(out, ".rad_store") <- store
  if (!is.null(field)) attr(out, ".rad_field") <- field
  out
}

#' Pack reads into a compact data.table
#'
#' Creates a genuine `data.table` whose `seq_id`, `seq`, and `qual` columns
#' contain 64-bit handles into a shared C++ store. IDs are front-coded (or
#' stored as 16-byte canonical UUIDs), sequences use a two-bit A/C/T/G core
#' with a lossless exception lane, and qualities use seven-bit numeric Phred.
#'
#' @param seq_id Character vector of read identifiers.
#' @param seq Character vector of sequences.
#' @param qual Optional FASTQ quality strings. `NA` represents missing quality.
#' @param seq_comment Optional FASTQ header comments, stored separately.
#' @param phred_offset Integer FASTQ offset. No automatic detection is done.
#' @return A `rad_reads` object inheriting from `data.table`.
rad_pack_reads <- function(seq_id, seq, qual = NULL, seq_comment = NULL,
                           phred_offset = 33L) {
  if (!is.character(seq_id)) stop("seq_id must be a character vector", call. = FALSE)
  if (!is.character(seq)) stop("seq must be a character vector", call. = FALSE)
  if (!is.null(qual) && !is.character(qual)) {
    stop("qual must be NULL or a character vector", call. = FALSE)
  }
  if (!is.null(seq_comment) && !is.character(seq_comment)) {
    stop("seq_comment must be NULL or a character vector", call. = FALSE)
  }
  seq_id <- enc2utf8(seq_id)
  if (!is.null(seq_comment)) seq_comment <- enc2utf8(seq_comment)
  out <- packed_reads_build_cpp(
    ids = seq_id,
    sequences = seq,
    qualities_sexp = qual,
    comments_sexp = seq_comment,
    phred_offset = as.integer(phred_offset)
  )
  .as_rad_reads(out)
}

#' Convert a conventional read table into compact RAD storage
#'
#' @param x A data.frame or data.table.
#' @param id_col,seq_col,qual_col Column names. `qual_col = NULL` creates
#'   records without quality.
#' @param comment_col Optional header-comment column name.
#' @param phred_offset Integer FASTQ offset.
#' @return A compact `rad_reads` data.table.
as_rad_reads <- function(x, id_col = if ("seq_id" %in% names(x)) "seq_id" else "id",
                         seq_col = "seq", qual_col = "qual", comment_col = NULL,
                         phred_offset = 33L) {
  if (!is.data.frame(x)) stop("x must be a data.frame or data.table", call. = FALSE)
  needed <- c(id_col, seq_col, qual_col, comment_col)
  needed <- needed[!is.null(needed)]
  missing_columns <- setdiff(needed, names(x))
  if (length(missing_columns)) {
    stop("missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  rad_pack_reads(
    seq_id = x[[id_col]],
    seq = x[[seq_col]],
    qual = if (is.null(qual_col)) NULL else x[[qual_col]],
    seq_comment = if (is.null(comment_col)) NULL else x[[comment_col]],
    phred_offset = phred_offset
  )
}

#' Open a streaming FASTQ iterator
#'
#' The parser keeps only its current kseq record. Each call to
#' [rad_stream_next()] returns a new independently owned packed chunk.
#'
#' @param path FASTQ or gzip-compressed FASTQ path.
#' @param phred_offset Integer FASTQ offset.
#' @param keep_comment Preserve the comment following the primary read ID in a
#'   separate `seq_comment` column.
#' @return A `rad_stream` iterator.
rad_stream <- function(path, phred_offset = 33L, keep_comment = TRUE) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("path must be one non-empty filename", call. = FALSE)
  }
  path <- normalizePath(path, mustWork = TRUE)
  pointer <- packed_fastq_stream_open_cpp(
    path = path,
    phred_offset = as.integer(phred_offset),
    keep_comment = isTRUE(keep_comment)
  )
  structure(
    list(pointer = pointer, path = path, phred_offset = as.integer(phred_offset)),
    class = "rad_stream"
  )
}

#' Read the next packed FASTQ chunk
#'
#' @param stream A `rad_stream` returned by [rad_stream()].
#' @param n Maximum records in the chunk.
#' @return A `rad_reads` data.table, or `NULL` at end of file.
rad_stream_next <- function(stream, n = 50000L) {
  if (!inherits(stream, "rad_stream")) stop("stream must be a rad_stream", call. = FALSE)
  if (length(n) != 1L || is.na(n) || n <= 0 || n > .Machine$integer.max) {
    stop("n must be one positive integer", call. = FALSE)
  }
  out <- packed_fastq_stream_next_cpp(stream$pointer, as.integer(n))
  if (is.null(out)) return(NULL)
  .as_rad_reads(out)
}

#' Test whether a FASTQ stream has reached EOF
rad_stream_done <- function(stream) {
  if (!inherits(stream, "rad_stream")) stop("stream must be a rad_stream", call. = FALSE)
  isTRUE(unname(packed_fastq_stream_state_cpp(stream$pointer)[["done"]]))
}

#' Close a FASTQ stream
rad_stream_close <- function(stream) {
  if (!inherits(stream, "rad_stream")) stop("stream must be a rad_stream", call. = FALSE)
  state <- packed_fastq_stream_state_cpp(stream$pointer)
  if (!isTRUE(unname(state[["closed"]]))) packed_fastq_stream_close_cpp(stream$pointer)
  invisible(stream)
}

#' Read a FASTQ directly into packed storage
#'
#' @param path FASTQ or gzip-compressed FASTQ path.
#' @param max_reads Maximum number of records.
#' @param phred_offset Integer FASTQ offset.
#' @param keep_comment Preserve header comments separately.
#' @return A compact `rad_reads` data.table, or an empty packed table for an
#'   empty file.
rad_read_fastq <- function(path, max_reads = .Machine$integer.max,
                           phred_offset = 33L, keep_comment = TRUE) {
  stream <- rad_stream(path, phred_offset = phred_offset, keep_comment = keep_comment)
  on.exit(rad_stream_close(stream), add = TRUE)
  out <- rad_stream_next(stream, n = max_reads)
  if (!is.null(out)) return(out)
  rad_pack_reads(
    character(), character(), character(),
    if (isTRUE(keep_comment)) character() else NULL,
    phred_offset = phred_offset
  )
}

#' Materialize packed read columns
#'
#' This is the explicit escape hatch for code that requires ordinary character
#' vectors. It copies only the selected rows represented by `x`.
rad_materialize <- function(x, phred_offset = 33L) {
  if (!inherits(x, "rad_reads")) stop("x must be a rad_reads table", call. = FALSE)
  rad_validate(x)
  sequence <- as.character(x$seq)
  quality <- as.character(x$qual, phred_offset = phred_offset)
  .rad_validate_cached_summaries(x, sequence, quality)
  out <- data.table::copy(x)
  data.table::setattr(out, "class", c("data.table", "data.frame"))
  data.table::set(out, j = "seq_id", value = as.character(x$seq_id))
  if ("seq_comment" %in% names(x)) {
    data.table::set(out, j = "seq_comment", value = as.character(x$seq_comment))
  }
  data.table::set(out, j = "seq", value = sequence)
  data.table::set(out, j = "qual", value = quality)
  attr(out, ".rad_store") <- NULL
  out
}

.rad_private_staging_path <- function(target, filename = "payload") {
  parent <- dirname(target)
  for (attempt in seq_len(32L)) {
    staging_directory <- tempfile(pattern = ".rrad-stage-", tmpdir = parent)
    if (dir.create(staging_directory, mode = "0700", showWarnings = FALSE)) {
      return(file.path(staging_directory, filename))
    }
  }
  stop("could not create a private staging directory beside: ", target,
       call. = FALSE)
}

#' Write packed reads directly to FASTQ
rad_write_fastq <- function(x, path, phred_offset = 33L, append = FALSE) {
  if (!inherits(x, "rad_reads")) stop("x must be a rad_reads table", call. = FALSE)
  if (!all(c("seq_id", "seq", "qual") %in% names(x))) {
    stop("x must contain seq_id, seq, and qual", call. = FALSE)
  }
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path)) {
    stop("path must be one non-empty filename", call. = FALSE)
  }
  if (!is.numeric(phred_offset) || length(phred_offset) != 1L ||
      is.na(phred_offset) || !is.finite(phred_offset) ||
      phred_offset != trunc(phred_offset) || phred_offset < 0L ||
      phred_offset > 126L) {
    stop("phred_offset must be one integer between 0 and 126", call. = FALSE)
  }
  if (!is.logical(append) || length(append) != 1L || is.na(append)) {
    stop("append must be TRUE or FALSE", call. = FALSE)
  }
  rad_validate(x)
  if (anyNA(x$qual)) stop("FASTQ output cannot contain missing qualities", call. = FALSE)
  target <- path.expand(path)
  parent <- dirname(target)
  if (!dir.exists(parent)) {
    stop("FASTQ output directory does not exist: ", parent, call. = FALSE)
  }
  staging <- .rad_private_staging_path(
    target,
    if (grepl("\\.gz$", target, ignore.case = TRUE)) "payload.gz" else "payload"
  )
  on.exit(unlink(dirname(staging), recursive = TRUE), add = TRUE)
  if (append && file.exists(target) && !file.copy(target, staging)) {
    stop("could not stage existing FASTQ output: ", target, call. = FALSE)
  }
  comments <- if ("seq_comment" %in% names(x)) x$seq_comment else NULL
  packed_write_fastq_cpp(
    id_handles = x$seq_id,
    sequence_handles = x$seq,
    quality_handles = x$qual,
    comment_handles = comments,
    path = staging,
    phred_offset = as.integer(phred_offset),
    append = append
  )
  if (!file.rename(staging, target)) {
    stop("could not atomically install FASTQ output: ", target, call. = FALSE)
  }
  invisible(normalizePath(target, mustWork = TRUE))
}

#' Inspect packed-store payload and index sizes
rad_storage_stats <- function(x) {
  handles <- if (inherits(x, "rad_reads")) x$seq else x
  if (!inherits(handles, c("rad_id", "rad_seq", "rad_phred"))) {
    stop("x must be rad_reads or a packed RAD handle column", call. = FALSE)
  }
  packed_store_stats_cpp(handles)
}

.rad_repack_materialized <- function(x, phred_offset = 33L) {
  .rad_assert_safe_extras(x, "repacked")
  comment <- if ("seq_comment" %in% names(x)) x[["seq_comment"]] else NULL
  packed <- rad_pack_reads(
    seq_id = x[["seq_id"]],
    seq = x[["seq"]],
    qual = x[["qual"]],
    seq_comment = comment,
    phred_offset = phred_offset
  )
  extras <- setdiff(names(x), .rad_core_columns)
  for (column in extras) data.table::set(packed, j = column, value = x[[column]])
  packed
}

#' Repack a filtered RAD table
#'
#' The ordinary subset operation is cheap because it copies handles. After a
#' large filter, `rad_compact()` rewrites only surviving records into a new
#' store so discarded payload can be released.
rad_compact <- function(x) {
  if (!inherits(x, "rad_reads")) stop("x must be a rad_reads table", call. = FALSE)
  .rad_assert_safe_extras(x, "repacked")
  rad_validate(x, deep = TRUE)
  phred_offset <- rad_storage_stats(x)$phred_offset
  .rad_repack_materialized(rad_materialize(x, phred_offset), phred_offset)
}

#' Safely combine packed read tables
#'
#' Raw `data.table::rbindlist()` must not be used across different stores,
#' because a handle is store-local. This function decodes and repacks the
#' records into one new store while preserving non-core columns.
rad_rbind_reads <- function(..., use.names = TRUE, fill = TRUE,
                            phred_offset = 33L) {
  inputs <- list(...)
  if (length(inputs) == 1L && is.list(inputs[[1L]]) &&
      !inherits(inputs[[1L]], "rad_reads")) {
    inputs <- inputs[[1L]]
  }
  if (!length(inputs)) {
    return(rad_pack_reads(character(), character(), character(),
                          phred_offset = phred_offset))
  }
  if (!all(vapply(inputs, inherits, logical(1), what = "rad_reads"))) {
    stop("every input must be a rad_reads table", call. = FALSE)
  }
  for (input in inputs) {
    .rad_assert_safe_extras(input, "combined")
    rad_validate(input, deep = TRUE)
  }
  materialized <- lapply(inputs, rad_materialize, phred_offset = phred_offset)
  combined <- data.table::rbindlist(materialized, use.names = use.names, fill = fill)
  .rad_repack_materialized(combined, phred_offset)
}

#' Save packed reads safely
#'
#' External pointers are not portable through `saveRDS()`. This method writes
#' a versioned, compressed, materialized payload and [rad_load()] recreates the
#' packed store on load. A native block-file format can replace this bridge
#' without changing the table API.
rad_save <- function(x, file, compress = "xz") {
  if (!inherits(x, "rad_reads")) stop("x must be a rad_reads table", call. = FALSE)
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("file must be one non-empty filename", call. = FALSE)
  }
  .rad_assert_safe_extras(x, "saved")
  rad_validate(x, deep = TRUE)
  phred_offset <- rad_storage_stats(x)$phred_offset
  payload <- list(
    format = "rrad_reads",
    version = 1L,
    phred_offset = phred_offset,
    data = rad_materialize(x, phred_offset = phred_offset)
  )
  target <- path.expand(file)
  parent <- dirname(target)
  if (!dir.exists(parent)) {
    stop("rad_save output directory does not exist: ", parent, call. = FALSE)
  }
  staging <- .rad_private_staging_path(target)
  on.exit(unlink(dirname(staging), recursive = TRUE), add = TRUE)
  saveRDS(payload, file = staging, compress = compress)
  if (!file.rename(staging, target)) {
    stop("could not atomically install RAD read store: ", target, call. = FALSE)
  }
  invisible(normalizePath(target, mustWork = TRUE))
}

#' Load reads written by [rad_save()]
rad_load <- function(file) {
  payload <- readRDS(file)
  if (!is.list(payload) || !identical(payload$format, "rrad_reads") ||
      !identical(payload$version, 1L) || !is.data.frame(payload$data)) {
    stop("file is not a supported RAD read store", call. = FALSE)
  }
  .rad_repack_materialized(payload$data, payload$phred_offset)
}

.rad_validate_cached_summaries <- function(x, sequence, quality) {
  expected <- nchar(sequence, type = "bytes", allowNA = FALSE, keepNA = FALSE)
  observed <- nchar(quality, type = "bytes", allowNA = TRUE, keepNA = TRUE)
  if (!identical(as.numeric(expected), as.numeric(x$length))) {
    stop("stored sequence lengths do not match the length column", call. = FALSE)
  }
  mismatch <- !is.na(observed) & observed != expected
  if (any(mismatch)) stop("stored sequence and quality lengths differ", call. = FALSE)

  same_numeric <- function(left, right, tolerance = sqrt(.Machine$double.eps)) {
    both_missing <- is.na(left) & is.na(right)
    both_present <- !is.na(left) & !is.na(right)
    all(both_missing | (both_present & abs(left - right) <= tolerance))
  }
  if (!same_numeric(x$mean_q, rad_qmean(x$qual))) {
    stop("stored mean_q values do not match packed qualities", call. = FALSE)
  }
  if (!same_numeric(as.numeric(x$min_q), as.numeric(rad_qmin(x$qual)),
                    tolerance = 0)) {
    stop("stored min_q values do not match packed qualities", call. = FALSE)
  }
  if (!same_numeric(x$frac_q20, rad_qfraction(x$qual, 20L))) {
    stop("stored frac_q20 values do not match packed qualities", call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate packed handles and sequence/quality invariants
rad_validate <- function(x, deep = FALSE) {
  if (!inherits(x, "rad_reads") || !data.table::is.data.table(x)) {
    stop("x must be a rad_reads data.table", call. = FALSE)
  }
  required <- c("read_key", "seq_id", "seq", "qual", "length",
                "mean_q", "min_q", "frac_q20")
  if (!all(required %in% names(x))) stop("rad_reads core columns are missing", call. = FALSE)
  if (!inherits(x$seq_id, "rad_id") || !inherits(x$seq, "rad_seq") ||
      !inherits(x$qual, "rad_phred")) {
    stop("rad_reads handle column classes are invalid", call. = FALSE)
  }
  if (!bit64::is.integer64(x$read_key)) {
    stop("rad_reads read_key column must use integer64 storage", call. = FALSE)
  }
  if (!all(vapply(c("length", "mean_q", "min_q", "frac_q20"),
                  function(column) is.numeric(x[[column]]), logical(1)))) {
    stop("rad_reads length and quality summary columns must be numeric", call. = FALSE)
  }
  stores <- list(
    attr(x$seq_id, ".rad_store", exact = TRUE),
    attr(x$seq, ".rad_store", exact = TRUE),
    attr(x$qual, ".rad_store", exact = TRUE)
  )
  if (any(vapply(stores, is.null, logical(1))) ||
      !identical(stores[[1L]], stores[[2L]]) || !identical(stores[[1L]], stores[[3L]])) {
    stop("rad_reads columns do not share one backing store", call. = FALSE)
  }
  comments <- if ("seq_comment" %in% names(x)) x$seq_comment else NULL
  packed_validate_handles_cpp(x$seq_id, x$seq, x$qual, comments)
  sequence_keys <- x$seq
  class(sequence_keys) <- "integer64"
  expected_keys <- sequence_keys %% bit64::as.integer64(2^32)
  if (length(x$read_key) != length(expected_keys) || anyNA(x$read_key) ||
      any(x$read_key != expected_keys)) {
    stop("read_key and sequence handles are not aligned", call. = FALSE)
  }
  if (isTRUE(deep)) {
    sequence <- as.character(x$seq)
    quality <- as.character(x$qual)
    .rad_validate_cached_summaries(x, sequence, quality)
  }
  invisible(TRUE)
}

#' Return numeric quality scores
#'
#' @return A list with one integer vector per read. Missing qualities are NULL.
rad_qual_values <- function(x) {
  if (inherits(x, "rad_reads")) x <- x$qual
  if (!inherits(x, "rad_phred")) stop("x must be rad_phred", call. = FALSE)
  packed_quality_values_cpp(x)
}

#' Compute mean Phred for packed qualities
rad_qmean <- function(x) {
  if (inherits(x, "rad_reads")) x <- x$qual
  if (!inherits(x, "rad_phred")) stop("x must be rad_phred", call. = FALSE)
  packed_quality_mean_cpp(x)
}

#' Compute minimum Phred for packed qualities
rad_qmin <- function(x) {
  if (inherits(x, "rad_reads")) x <- x$qual
  if (!inherits(x, "rad_phred")) stop("x must be rad_phred", call. = FALSE)
  packed_quality_min_cpp(x)
}

#' Compute the fraction of scores at or above a threshold
rad_qfraction <- function(x, threshold = 20L) {
  if (inherits(x, "rad_reads")) x <- x$qual
  if (!inherits(x, "rad_phred")) stop("x must be rad_phred", call. = FALSE)
  packed_quality_fraction_cpp(x, as.integer(threshold))
}

.rad_handle_is_na <- function(x) packed_handle_is_na_cpp(x)

`[.rad_id` <- function(x, i, ...) .rad_handle_subset(x, i)
`[[.rad_id` <- function(x, i, ...) .rad_handle_subset(x, i)
`[.rad_seq` <- function(x, i, ...) .rad_handle_subset(x, i)
`[[.rad_seq` <- function(x, i, ...) .rad_handle_subset(x, i)
`[.rad_phred` <- function(x, i, ...) .rad_handle_subset(x, i)
`[[.rad_phred` <- function(x, i, ...) .rad_handle_subset(x, i)

c.rad_id <- .rad_handle_c
c.rad_seq <- .rad_handle_c
c.rad_phred <- .rad_handle_c
rep.rad_id <- .rad_handle_rep
rep.rad_seq <- .rad_handle_rep
rep.rad_phred <- .rad_handle_rep

.rad_handle_replace <- function(x, i, ..., value) {
  stop("packed RAD handle columns are immutable; repack the affected reads instead",
       call. = FALSE)
}
`[<-.rad_id` <- .rad_handle_replace
`[[<-.rad_id` <- .rad_handle_replace
`[<-.rad_seq` <- .rad_handle_replace
`[[<-.rad_seq` <- .rad_handle_replace
`[<-.rad_phred` <- .rad_handle_replace
`[[<-.rad_phred` <- .rad_handle_replace

as.character.rad_id <- function(x, ...) packed_id_decode_cpp(x)
as.character.rad_seq <- function(x, ...) packed_sequence_decode_cpp(x)
as.character.rad_phred <- function(x, ..., phred_offset = 33L) {
  packed_quality_decode_cpp(x, as.integer(phred_offset))
}

as.integer.rad_phred <- function(x, ...) {
  if (length(x) != 1L) {
    stop("as.integer.rad_phred() requires one read; use rad_qual_values() for a column",
         call. = FALSE)
  }
  packed_quality_values_cpp(x)[[1L]]
}


is.na.rad_id <- .rad_handle_is_na
is.na.rad_seq <- .rad_handle_is_na
is.na.rad_phred <- .rad_handle_is_na
anyNA.rad_id <- function(x, recursive = FALSE) any(.rad_handle_is_na(x))
anyNA.rad_seq <- function(x, recursive = FALSE) any(.rad_handle_is_na(x))
anyNA.rad_phred <- function(x, recursive = FALSE) any(.rad_handle_is_na(x))

Ops.rad_id <- function(e1, e2) {
  comparisons <- c("==", "!=", "<", "<=", ">", ">=")
  if (!.Generic %in% comparisons) {
    stop(.Generic, " is not defined for compressed read IDs", call. = FALSE)
  }
  if (inherits(e1, "rad_id")) e1 <- as.character(e1)
  if (!missing(e2) && inherits(e2, "rad_id")) e2 <- as.character(e2)
  do.call(.Generic, list(e1, e2), envir = baseenv())
}

duplicated.rad_id <- function(x, incomparables = FALSE, fromLast = FALSE, nmax = NA, ...) {
  duplicated(as.character(x), incomparables = incomparables,
             fromLast = fromLast, nmax = nmax, ...)
}

anyDuplicated.rad_id <- function(x, incomparables = FALSE, fromLast = FALSE, ...) {
  anyDuplicated(as.character(x), incomparables = incomparables,
                fromLast = fromLast, ...)
}

unique.rad_id <- function(x, incomparables = FALSE, fromLast = FALSE, ...) {
  x[!duplicated.rad_id(x, incomparables = incomparables,
                       fromLast = fromLast, ...)]
}

xtfrm.rad_id <- function(x) xtfrm(as.character(x))

sort.rad_id <- function(x, decreasing = FALSE, ...) {
  x[order(as.character(x), decreasing = decreasing, ...)]
}

format.rad_id <- function(x, ...) as.character(x)
format.rad_seq <- function(x, ..., width = 24L) {
  packed_sequence_format_cpp(x, as.integer(width))
}
format.rad_phred <- function(x, ...) packed_quality_format_cpp(x)

print.rad_id <- function(x, ...) {
  print(format(x, ...), quote = FALSE)
  invisible(x)
}
print.rad_seq <- function(x, ...) {
  print(format(x, ...), quote = FALSE)
  invisible(x)
}
print.rad_phred <- function(x, ...) {
  print(format(x, ...), quote = FALSE)
  invisible(x)
}
