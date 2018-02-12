valid.or.default <- function (x, def) {
  if (is.null(x) || is.na(x)) def else x
}

warning.predict.fixed.effects.with.unknown.indicators <- function(row, col, inds) {
  col.str <- paste("column", valid.or.default(colnames(inds)[col], col))
  row.str <- paste("row",    valid.or.default(rownames(inds)[row], row))

  sprintf("Indicator \"%s\" in [%s, %s] has not occured before and cannot be estimated; the predicted value for %s will be NA.",
    as.character(inds[row, col]), col.str, row.str, row.str)
}

warning.contain.multiple.components <- function() {
  sprintf("Given indicators contain multiple disconnected components. Prediction of entries containing fixed effect values in different components is impossible.")
}

warning.cross.component <- function(error, inds, group.levels) {
  list2env(error, environment())

  row.str <- paste("row", valid.or.default(rownames(inds)[row], row))

  group.x.str <- paste("column", valid.or.default(colnames(inds)[group.x], group.x))
  group.y.str <- paste("column", valid.or.default(colnames(inds)[group.y], group.y))

  value.x.str <- paste((group.levels[[group.x]])[value.x])
  value.y.str <- paste((group.levels[[group.y]])[value.y])

  sprintf("In [%s], indicator \"%s\" in [%s] and indicator \"%s\" in [%s] are from disconnected components; the predicted value for %s will be NA.",
    row.str, value.x.str, group.x.str, value.y.str, group.y.str, row.str)
}
