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
  if (is.null(newX)) {
    if (is.null(FEValues))
      return("newX and FEValues cannot be both NULL.")
    newX <- matrix(data = NA, nrow = dim(FEValues)[1], ncol = 0)
  }
  
  if (dim(newX)[2] != dim(model$coefficients)[1])
    return("incompatible width for input X")
    
  y <- newX %*% model$coefficients + model$intercept + grandMean
  
  if (!is.null(FEValues)) {
    height <- dim(FEValues)[1]
    width <- dim(FEValues)[2]
    
    for (i in 1 : width)
      FEValues[, i] <- as.numeric(factor(FEValues[, i], levels = model$.group.levels[[i]]))

    for (i in 1 : width) {
      for (j in 1 : height) {
        y[j] <- y[j] + (model$FEcoefs[[i]])[FEValues[j, i]]
      }
    }
  }
  
  return(y)
}
