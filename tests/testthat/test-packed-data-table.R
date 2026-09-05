test_that("data.table subset and reorder operations preserve packed handles", {
  ids <- paste0("read-", 1:6)
  sequences <- c("A", "CC", "TGN", "acgt", "GATTACA", "N-RY")
  scores <- list(10L, c(20L, 21L), c(30L, 31L, 32L),
                 c(40L, 41L, 42L, 43L), rep(50L, 7), c(0L, 1L, 2L, 3L))
  qualities <- vapply(scores, quality_string, character(1))
  packed <- rad_pack_reads(ids, sequences, qualities)

  selected <- packed[c(6L, 2L, 2L, 4L, 1L)]
  expected <- c(6L, 2L, 2L, 4L, 1L)
  expect_s3_class(selected, "rad_reads")
  expect_s3_class(selected$seq_id, "rad_id")
  expect_s3_class(selected$seq, "rad_seq")
  expect_s3_class(selected$qual, "rad_phred")
  expect_identical(as.character(selected$seq_id), ids[expected])
  expect_identical(as.character(selected$seq), sequences[expected])
  expect_identical(as.character(selected$qual), qualities[expected])
  expect_identical(rad_qual_values(selected), scores[expected])

  selected[, requested_order := c(4L, 1L, 5L, 3L, 2L)]
  data.table::setorder(selected, requested_order)
  expected <- expected[order(c(4L, 1L, 5L, 3L, 2L))]
  expect_identical(as.character(selected$seq_id), ids[expected])
  expect_identical(as.character(selected$seq), sequences[expected])
  expect_identical(as.character(selected$qual), qualities[expected])
  expect_s3_class(selected$qual, "rad_phred")
})

test_that("data.table assignment cannot silently replace packed core columns", {
  packed <- rad_pack_reads(c("a", "b"), c("A", "CC"), c("!", "55"))

  expect_error(
    packed[, seq := rev(seq)],
    "packed RAD core columns are immutable"
  )
  expect_identical(as.character(packed$seq), c("A", "CC"))

  protected <- c("seq_id", "qual")
  expect_error(
    packed[, (protected) := list(rev(seq_id), rev(qual))],
    "seq_id, qual"
  )
  expect_identical(as.character(packed$seq_id), c("a", "b"))

  extra <- "programmatic_annotation"
  packed[, (extra) := c("x", "y")]
  expect_identical(packed$programmatic_annotation, c("x", "y"))

  evaluations <- 0L
  changing_target <- function() {
    evaluations <<- evaluations + 1L
    if (evaluations == 1L) "another_annotation" else "seq"
  }
  expect_error(
    packed[, (changing_target()) := "unsafe"],
    "cannot verify a dynamic by-reference assignment"
  )
  expect_identical(evaluations, 0L)
  expect_identical(as.character(packed$seq), c("A", "CC"))
  expect_false("another_annotation" %in% names(packed))

  before_dynamic_set <- data.table::copy(packed)
  expect_error(
    packed[, {
      do.call(
        data.table::set,
        list(x = packed, j = "seq", value = rev(packed$seq))
      )
      NULL
    }],
    "all changes were rolled back"
  )
  expect_true(isTRUE(rad_validate(packed, deep = TRUE)))
  expect_identical(rad_materialize(packed), rad_materialize(before_dynamic_set))

  set_alias <- data.table::set
  expect_error(
    packed[, {
      set_alias(packed, j = "qual", value = rev(packed$qual))
      NULL
    }],
    "all changes were rolled back"
  )
  expect_true(isTRUE(rad_validate(packed, deep = TRUE)))
  expect_identical(rad_materialize(packed), rad_materialize(before_dynamic_set))

  expect_error(
    packed[, { seq := rev(seq) }],
    "packed RAD core columns are immutable"
  )
  expect_error(
    packed[, `:=`(seq = rev(seq))],
    "packed RAD core columns are immutable"
  )
  expect_error(
    packed[, J, env = list(J = quote(seq := rev(seq)))],
    "packed RAD core columns are immutable"
  )
  expect_error(
    packed[, eval(parse(text = "seq := rev(seq)"))],
    "packed RAD core columns are immutable"
  )
  indirect <- quote(seq := rev(seq))
  expect_error(
    packed[, eval(indirect)],
    "packed RAD core columns are immutable"
  )
  expect_error(
    packed[, eval(parse(
      text = "data.table::set(packed, j = 'seq', value = rev(packed$seq))"
    ))],
    "cannot verify a dynamic by-reference assignment"
  )
  expect_identical(as.character(packed$seq), c("A", "CC"))

  packed[, annotation := c("first", "second")]
  expect_identical(packed$annotation, c("first", "second"))
  packed[, `:=`(named_annotation = c(1L, 2L))]
  expect_identical(packed$named_annotation, c(1L, 2L))

  expect_error(packed$seq <- rev(packed$seq),
               "packed RAD core columns are immutable")
  expect_error(packed[["qual"]] <- rev(packed[["qual"]]),
               "packed RAD core columns are immutable")
  expect_error(packed[, "length"] <- list(rev(packed$length)),
               "packed RAD core columns are immutable")
  packed$annotation <- c("updated-first", "updated-second")
  packed[["named_annotation"]] <- c(3L, 4L)
  expect_identical(packed$annotation, c("updated-first", "updated-second"))
  expect_identical(packed$named_annotation, c(3L, 4L))

  corrupted <- data.table::copy(packed)
  data.table::set(corrupted, j = "seq", value = rev(corrupted$seq))
  expect_error(rad_materialize(corrupted), "not aligned")

  corrupted_summary <- data.table::copy(packed)
  data.table::set(corrupted_summary, j = "length", value = c(9, 9))
  expect_error(rad_materialize(corrupted_summary),
               "stored sequence lengths do not match")

  corrupted_key <- data.table::copy(packed)
  data.table::set(corrupted_key, j = "read_key",
                  value = rev(corrupted_key$read_key))
  expect_error(rad_validate(corrupted_key),
               "read_key and sequence handles are not aligned")
})

test_that("data.table projections and joins cannot masquerade as packed read tables", {
  packed <- rad_pack_reads(
    c("a", "b", "c"), c("A", "CC", "GGG"), c("!", "55", "III")
  )

  projected <- packed[, .(seq_id, seq)]
  expect_s3_class(projected, "data.table")
  expect_false(inherits(projected, "rad_reads"))
  expect_identical(as.character(projected$seq_id), c("a", "b", "c"))
  projected_columns <- quote(.(seq_id, seq))
  evaluated_projection <- packed[, eval(projected_columns)]
  expect_s3_class(evaluated_projection, "data.table")
  expect_false(inherits(evaluated_projection, "rad_reads"))
  grouped <- packed[, .N, by = seq_id]
  expect_s3_class(grouped, "data.table")
  expect_false(inherits(grouped, "rad_reads"))

  other <- rad_pack_reads(c("c", "a"), c("T", "G"), c("~", "I"))
  index <- data.table::data.table(seq_id = other$seq_id, label = c("C", "A"))
  expect_error(
    packed[index, on = "seq_id", nomatch = NULL],
    "do not share one backing store"
  )
  expect_true(isTRUE(rad_validate(packed)))

  merged <- merge(packed, index, by = "seq_id")
  expect_s3_class(merged, "rad_reads")
  expect_true(isTRUE(rad_validate(merged, deep = TRUE)))
  expect_identical(as.character(merged$seq_id), c("a", "c"))

  expect_error(merge(packed, other, by = "seq_id"),
               "do not share one backing store")

  keyed <- data.table::copy(packed)
  data.table::setkey(keyed, seq_id)
  expect_true(isTRUE(rad_validate(keyed, deep = TRUE)))
  expect_identical(
    setNames(as.character(keyed$seq), as.character(keyed$seq_id))[c("a", "b", "c")],
    c(a = "A", b = "CC", c = "GGG")
  )
})

test_that("copied and subset tables keep their stores alive through garbage collection", {
  survivor <- local({
    original <- rad_pack_reads(
      paste0("gc/", 1:4),
      c("ACGT", "NNNN", "tgcA", "G-RY"),
      vapply(list(0:3, 20:23, 40:43, 90:93), quality_string, character(1))
    )
    copied <- data.table::copy(original[c(4L, 1L, 3L)])
    copied[, ordinary_column := c("d", "a", "c")]
    copied
  })

  invisible(gc())
  expect_identical(as.character(survivor$seq_id), c("gc/4", "gc/1", "gc/3"))
  expect_identical(as.character(survivor$seq), c("G-RY", "ACGT", "tgcA"))
  expect_identical(rad_qual_values(survivor), list(90:93, 0:3, 40:43))
  expect_identical(survivor$ordinary_column, c("d", "a", "c"))

  column_only <- local({
    temporary <- rad_pack_reads(c("x", "y"), c("A", "N"), c("!", "~"))
    temporary$seq[c(2L, 1L)]
  })
  invisible(gc())
  expect_identical(as.character(column_only), c("N", "A"))
})

test_that("as_rad_reads converts conventional tables without modifying them", {
  source <- data.frame(
    id = c("r1", "r2"),
    comment = c("first", NA_character_),
    seq = c("AC", "TG"),
    quality = c("!!", "I~"),
    stringsAsFactors = FALSE
  )
  original <- source
  packed <- as_rad_reads(source, id_col = "id", seq_col = "seq",
                         qual_col = "quality", comment_col = "comment")

  expect_identical(source, original)
  expect_identical(materialized_payload(packed), data.frame(
    seq_id = source$id,
    seq_comment = source$comment,
    seq = source$seq,
    qual = source$quality,
    stringsAsFactors = FALSE
  ))

  no_quality <- as_rad_reads(source, id_col = "id", seq_col = "seq",
                             qual_col = NULL, comment_col = "comment")
  expect_identical(rad_materialize(no_quality)$qual, c(NA_character_, NA_character_))
})

test_that("packed columns encode missing integer subset positions safely", {
  packed <- rad_pack_reads(
    c("id-1", "id-2", "id-3"),
    c("A", "CC", "GGG"),
    c("!", "55", "III")
  )
  index <- c(1L, NA_integer_, 3L)

  ids <- packed$seq_id[index]
  sequences <- packed$seq[index]
  qualities <- packed$qual[index]

  expect_s3_class(ids, "rad_id")
  expect_s3_class(sequences, "rad_seq")
  expect_s3_class(qualities, "rad_phred")
  expect_identical(is.na(ids), c(FALSE, TRUE, FALSE))
  expect_identical(is.na(sequences), c(FALSE, TRUE, FALSE))
  expect_identical(is.na(qualities), c(FALSE, TRUE, FALSE))
  expect_true(anyNA(ids))
  expect_true(anyNA(sequences))
  expect_true(anyNA(qualities))
  expect_identical(as.character(ids), c("id-1", NA_character_, "id-3"))
  expect_identical(as.character(qualities), c("!", NA_character_, "III"))
  expect_identical(rad_qual_values(qualities), list(0L, NULL, c(40L, 40L, 40L)))
  expect_equal(rad_qmean(qualities), c(0, NA, 40))
  expect_identical(rad_qmin(qualities), c(0L, NA_integer_, 40L))
  expect_equal(rad_qfraction(qualities, 20L), c(0, NA, 1))
  expect_error(as.character(sequences), "missing packed handle is not permitted")
  expect_identical(as.character(sequences[!is.na(sequences)]), c("A", "GGG"))
})

test_that("packed columns encode logical NA subset positions safely", {
  packed <- rad_pack_reads(
    c("id-1", "id-2", "id-3"),
    c("A", "CC", "GGG"),
    c("!", "55", "III")
  )
  index <- c(TRUE, NA, FALSE)

  ids <- packed$seq_id[index]
  sequences <- packed$seq[index]
  qualities <- packed$qual[index]

  expect_identical(is.na(ids), c(FALSE, TRUE))
  expect_identical(is.na(sequences), c(FALSE, TRUE))
  expect_identical(is.na(qualities), c(FALSE, TRUE))
  expect_identical(as.character(ids), c("id-1", NA_character_))
  expect_identical(as.character(qualities), c("!", NA_character_))
  expect_identical(rad_qual_values(qualities), list(0L, NULL))
  expect_error(as.character(sequences), "missing packed handle is not permitted")
  expect_identical(as.character(sequences[!is.na(sequences)]), "A")
})

test_that("packed columns encode out-of-bounds subset positions safely", {
  packed <- rad_pack_reads(
    c("id-1", "id-2", "id-3"),
    c("A", "CC", "GGG"),
    c("!", "55", "III")
  )
  index <- c(2L, 99L)

  ids <- packed$seq_id[index]
  sequences <- packed$seq[index]
  qualities <- packed$qual[index]

  expect_identical(is.na(ids), c(FALSE, TRUE))
  expect_identical(is.na(sequences), c(FALSE, TRUE))
  expect_identical(is.na(qualities), c(FALSE, TRUE))
  expect_identical(as.character(ids), c("id-2", NA_character_))
  expect_identical(as.character(qualities), c("55", NA_character_))
  expect_identical(rad_qual_values(qualities), list(c(20L, 20L), NULL))
  expect_error(as.character(sequences), "missing packed handle is not permitted")
  expect_identical(as.character(sequences[!is.na(sequences)]), "CC")
})

test_that("scalar packed views preserve inserted missing handles", {
  packed <- rad_pack_reads(c("id-1", "id-2"), c("A", "CC"), c("!", "55"))

  for (index in list(NA_integer_, 99L)) {
    id <- packed$seq_id[[index]]
    sequence <- packed$seq[[index]]
    quality <- packed$qual[[index]]

    expect_length(id, 1L)
    expect_length(sequence, 1L)
    expect_length(quality, 1L)
    expect_true(is.na(id))
    expect_true(is.na(sequence))
    expect_true(is.na(quality))
    expect_identical(as.character(id), NA_character_)
    expect_identical(as.character(quality), NA_character_)
    expect_identical(rad_qual_values(quality), list(NULL))
    expect_error(as.character(sequence), "missing packed handle is not permitted")
  }
})
