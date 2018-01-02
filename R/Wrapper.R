solveFE <- function(rawData, rawFixedEffects, coreNum = 1, estimateFE = FALSE) {
  groupLevels <- list()
  if (!is.null(rawFixedEffects)) {
    height <- dim(rawFixedEffects)[1]
    width <- dim(rawFixedEffects)[2]

    newFixedEffects <- matrix(0, height, width)

    for (i in 1 : width) {
      theFactor <- factor(rawFixedEffects[, i])
      newFixedEffects[, i] <- as.numeric(theFactor)
      groupLevels[[i]] <- levels(theFactor)
    }

    rawFixedEffects <- newFixedEffects
  }

  result <- internalSolveFE(rawData, rawFixedEffects, coreNum, estimateFE)
  
  if (estimateFE) {
    width <- dim(rawFixedEffects)[2]
    effectNames <- c()
    i <- 1
    while (i <= width) {
      if (is.null(colnames(rawFixedEffects)[i]))
        effectNames <- c(effectNames, paste("effect", i, sep = ""))
      else
        effectNames <- c(effectNames, colnames(rawFixedEffects)[i])
      i <- i + 1
    }
    names(result$FEcoefs) <- effectNames

    result$.group.levels <- groupLevels
  }
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
