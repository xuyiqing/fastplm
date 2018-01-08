solveFE <- function(rawData, rawFixedEffects, coreNum = 1, estimateFE = FALSE) {
  solve.fixed.effects(rawData, rawFixedEffects, coreNum)
}

solve.fixed.effects <- function(data, inds, core.num = 1) {
  if (!is.null(inds)) {
    N <- ncol(inds)

    factors <- lapply(1 : N, function(col) factor(inds[, col]))
    group.levels <- lapply(1 : N, function(col) levels(factors[[col]]))
    factored.inds <- sapply(1 : N, function(col) as.numeric(factors[[col]]))
  }

  result <- internalSolveFE(data, factored.inds, core.num)

  result$.group.levels <- group.levels

  names(result$FEcoefs) <-
    if (is.null(colnames(inds)))
      sapply(1 : ncol(inds), function(col) paste("effect", col, sep = ""))
    else
      colnames(inds)

  result
}

predictFE <- function(model, newX, FEValues = NULL, grandMean = 0) {
  predict.fixed.effects(model, newX, FEValues)
}

predict.fixed.effects <- function(model, x, inds = NULL) {
  if (is.null(x) && is.null(inds))
    stop("x and inds cannot both be NULL.")

  if (is.null(x))
    x <- matrix(data = NA, nrow = ncol(inds), ncol = 0)

  y <- x %*% model$coefficients + model$intercept

  if (is.null(inds))
    return(y)

  factored.inds <- sapply(1 : ncol(inds),
    function(col) as.numeric(factor(inds[, col], levels = model$.group.levels[[col]])))

  sapply(which(is.na(factored.inds)), function(i) {
    col <- (i - 1) %/% nrow(inds) + 1
    row <- i - (col - 1) * nrow(inds)
    warning(warning.predict.fixed.effects.with.unknown.indicators(row, col, inds))
  })

  with.effects <- function(row) sum(sapply(1 : ncol(inds),
    function(col) (model$FEcoefs[[col]])[factored.inds[row, col]]))

  y <- y + sapply(1 : nrow(inds), with.effects)

  return(y)
}
