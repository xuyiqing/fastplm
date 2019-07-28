create.indicators <- function(inds) {
  CHECK.INPUT(inds, "inds", "matrix")

  effect.names <-
    if (is.null(colnames(inds)))
      sapply(SEQ(1, ncol(inds)), function(col) paste("effect", col, sep = "."))
    else
      colnames(inds)

  factors     <- lapply(SEQ(1, ncol(inds)), function(col) factor(inds[, col]))
  levels      <- lapply(factors, levels)
  level.sizes <- sapply(levels, length)
  inds        <- sapply(factors, as.numeric)
  uid         <- sample.int(1000000, 1)

  r <- mget(c("effect.names", "levels", "level.sizes", "inds", "uid"), environment())
  class(r) <- "indicators"
  r
}

create.subindicators <- function(sub.inds, model = NULL, inds = NULL, fe = NULL) {
  CheckComponents <- NULL
  inds <- CHECK.INIT.INPUT(inds, "inds", function() model$inds, "indicators")
  fe   <- CHECK.INIT.INPUT(fe, "fe", function() model$fe, "fixed.effects")
  CHECK.INPUT(sub.inds, "sub.inds", "matrix")

  ASSERT.MATRIX.DIM(sub.inds, "sub.inds", ncol(inds$inds), is.width = TRUE)

  factored <- sapply(SEQ(1, ncol(sub.inds)),
    function(col) as.numeric(factor(sub.inds[, col], levels = inds$levels[[col]])))

  for (i in which(is.na(factored))) {
    col <- (i - 1) %/% nrow(sub.inds) + 1
    row <- i - (col - 1) * nrow(sub.inds)

    warning(WARN.INDS.unknown.indicator(row, col, sub.inds))
    factored[row, ] <- NA
  }

  cross.components <- CheckComponents(fe$ptr, factored)
  invalid.rows <- sapply(cross.components, function(err) {
    warning(WARN.INDS.cross.component(err, sub.inds, inds$levels))
    err$row
  })

  if (length(invalid.rows) > 0)
    factored[invalid.rows, ] <- NA

  sub.inds <- list(inds = factored, parent = inds$uid)
  class(sub.inds) <- "sub.indicators"
  sub.inds
}

WARN.INDS.unknown.indicator <- function(row, col, inds) {
  col.str <- paste("column", valid.or.default(colnames(inds)[col], col))
  row.str <- paste("row",    valid.or.default(rownames(inds)[row], row))

  sprintf("Unknown indicator \"%s\" in [%s, %s]. Row ignored.",
    paste(inds[row, col]), col.str, row.str, row.str)
}

WARN.INDS.cross.component <- function(error, inds, levels) {
  group.x <- group.y <- value.x <- value.y <- NULL
  list2env(error, environment())

  row.str <- paste("row", valid.or.default(rownames(inds)[row], row))

  group.x.str <- paste("column", valid.or.default(colnames(inds)[group.x], group.x))
  group.y.str <- paste("column", valid.or.default(colnames(inds)[group.y], group.y))

  value.x.str <- paste((levels[[group.x]])[value.x])
  value.y.str <- paste((levels[[group.y]])[value.y])

  sprintf("In [%s], indicator \"%s\" in [%s] and indicator \"%s\" in [%s] are from disconnected components; the predicted value for %s will be NA.",
    row.str, value.x.str, group.x.str, value.y.str, group.y.str, row.str)
}
