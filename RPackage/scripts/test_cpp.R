set.seed(123)
n <- 1000
nlvl <- 20
x1 <- rnorm(n, 3)
x2 <- rnorm(n, 3)
e <- rnorm(n, 1) # error

## generate 3 group indicators, each of 20 levels
gp <- matrix(sample(1:nlvl, n*3, replace = TRUE), n, 3)
colnames(gp) <- c("gp1","gp2","gp3")

## generate group effect, stored in a nlvl*3 matrix
gp.coef <- matrix(runif(nlvl * 3),nlvl,3)

## assign group eff to each observation
gp.eff <- matrix(NA, n, 3)
for (i in 1:n) {
  for (j in 1:3) {
    gp.eff[i,j] <- gp.coef[gp[i,j],j]
  } 
}
colnames(gp.eff) <- c("gp_eff1","gp_eff2","gp_eff3")

## outcome
y <- 5 + 1 * x1 + 3 * x2 + gp.eff[,1] + gp.eff[,2] + gp.eff[,3] + e 

sink("/Users/selveskii/Desktop/data.txt")
cat(1000, 3, "\n")
cat(y, "\n")
cat(x1, "\n")
cat(x2, "\n")
cat("\n\n")
cat(1000, 3, "\n")
cat(gp[,1], "\n")
cat(gp[,2], "\n")
cat(gp[,3], "\n")
sink()
