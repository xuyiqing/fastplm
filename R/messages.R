valid.or.default <- function (x, def) {
  if (is.null(x) || is.na(x)) def else x
}

warning.predict.fixed.effects.with.unknown.indicators <- function(row, col, inds) {
  col.str <- paste("column", valid.or.default(colnames(inds)[col], col))
  row.str <- paste("row",    valid.or.default(rownames(inds)[row], row))

  sprintf("[%s, %s]: indicator \"%s\" has not occured before and cannot be estimated; the predicted value for %s will be NA.",
    col.str, row.str, as.character(inds[row, col]), row.str)
}
