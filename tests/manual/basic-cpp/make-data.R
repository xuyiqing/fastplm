source("../../testthat/make-small-data.R", local = TRUE)

N <- 64000
LEVEL <- 64000
list2env(make.small.data(), environment())

data <- cbind(y, x)
write.table(data, file = "data.csv", sep=",", row.names = FALSE, col.names = FALSE)
write.table(inds, file = "inds.csv", sep=",", row.names = FALSE, col.names = FALSE)
