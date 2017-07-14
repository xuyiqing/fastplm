makeGPData <- function(paramCount = 4,
                       timeCount = 10, indivCount = 8,
                       timeParamCount = 3, indivParamCount = 2, seed = 57) {
  
  obsCount <- timeCount * indivCount
  
  errors <- matrix(rnorm(obsCount, 2), obsCount, 1)
  observations <- matrix(rnorm(obsCount * paramCount, 6), obsCount, paramCount)
  params <- matrix(rnorm(paramCount, 5), paramCount, 1)
  
  toisData <- matrix(rnorm(timeCount * indivParamCount, 4), timeCount, indivParamCount)
  iotsData <- matrix(rnorm(indivCount * timeParamCount, 4), indivCount, timeParamCount)
  
  toisParams <- matrix(rnorm(indivParamCount * indivCount, 4), indivParamCount, indivCount)
  iotsParams <- matrix(rnorm(timeParamCount * timeCount, 4), timeParamCount, timeCount)
  
  sum <- errors + observations %*% params
    + matrix(toisData %*% toisParams, obsCount, 1)
    + matrix(t(iotsData %*% iotsParams), obsCount, 1)

  list(paramCount = paramCount, timeCount = timeCount, indivCount = indivCount,
       timeParamCount = timeParamCount, indivParamCount = indivParamCount,
       sum = sum, errors = errors, observations = observations, params = params,
       toisData = toisData, iotsData = iotsData, toisParams = toisParams, iotsParams = iotsParams)
}

knockOutGPData <- function(data, missingCellCount) {
  obsCount <- data$timeCount * data$indivCount
  missingCellIds <- sample(1:obsCount, missingCellCount)
  sum <- data$sum
  
  for (id in missingCellIds)
    sum[id] <- NaN
  
  data$knockedOutSum <- data$sum[-missingCellIds,]
  data$knockedOutObservations <- data$observations[-missingCellIds,]
  data$sum <- sum
  
  data
}

writeGPData <- function(data, path) {
  printMatrix <- function(m) {
    write.table(format(m, justify="right"), row.names=F, col.names=F, quote=F)
  }
  
  sink(path)
  
  cat(data$paramCount, data$timeCount, data$indivCount, "\n")
  printMatrix(matrix(data$sum, data$timeCount, data$indivCount))
  cat("\n")
  for (i in 1:data$paramCount) {
    mat <- data$observations[,i]
    printMatrix(matrix(mat, data$timeCount, data$indivCount))
    cat("\n")
  }
  
  cat(data$indivParamCount, "\n")
  printMatrix(data$toisData)
  cat("\n")
  
  cat(data$timeParamCount, "\n")
  printMatrix(data$iotsData)
  cat("\n")
  
  sink()
}

data <- makeGPData()
data_ <- knockOutGPData(data, 10)

