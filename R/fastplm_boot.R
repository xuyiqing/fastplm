fastplm.boot <- function(seed, 
                         y, 
                         x, 
                         ind, 
                         sfe.index, 
                         cfe.index, 
                         cluster,
                         parallel = FALSE,
                         wild = FALSE, 
                         jackknife = FALSE,
                         refinement = FALSE, 
                         pos = NULL,
                         nboots, 
                         core.num) {  

    if (!is.null(seed)) {
        set.seed(seed)
    }
    p <- dim(x)[2]

    if (jackknife == 1) {
        wild <- 0
    }

    ## se <- ifelse(refinement == 1, 1, 0) 
    ## model.cl <- ifelse(refinement == 1, cluster, NULL)
    model <- fastplm.core(y = y, x = x, ind = ind, 
                          sfe.index = sfe.index, cfe.index = cfe.index, 
                          se = 1, cl = cluster, core.num = core.num)

    ## function to get two-sided p-values
    get.pvalue <- function(vec) {
        if (NaN%in%vec|NA%in%vec) {
            nan.pos <- is.nan(vec)
            na.pos <- is.na(vec)
            pos <- c(which(nan.pos),which(na.pos))
            vec.a <- vec[-pos]
            a <- sum(vec.a >= 0)/(length(vec)-sum(nan.pos|na.pos)) * 2
            b <- sum(vec.a <= 0)/(length(vec)-sum(nan.pos|na.pos)) * 2  
        } else {
            a <- sum(vec >= 0)/length(vec) * 2
            b <- sum(vec <= 0)/length(vec) * 2  
        }
        return(min(as.numeric(min(a, b)),1))
    }

    ## wild bootstrap
    if (wild == TRUE) {        
        
        raw.cluster <- level.cluster <- level.count <- NULL
        if (!is.null(cluster)) {
            ## one-way cluster 
            raw.cluster <- as.numeric(as.factor(cluster))
            level.cluster <- unique(raw.cluster)
            level.count <- table(raw.cluster)
        }

        if (refinement == 0) {
            ## wild bootstrap
            boot.coef <- matrix(NA, p, nboots)
            res <- model$residuals
            y.fit <- y - res 
            core.num <- ifelse(parallel == TRUE, 1, core.num)

            one.boot <- function(num = NULL) {
                
                if (!is.null(cluster)) {
                    ## wild boot by cluster 
                    res.p.level <- sample(c(-1, 1), length(level.cluster), replace = TRUE)
                    res.p <- c()
                    for (j in 1:length(level.cluster)) {
                        res.p <- c(res.p, rep(res.p.level[j], level.count[j]))
                    }
                    res.p <- as.matrix(res.p)
                } else {
                    res.p <- as.matrix(sample(c(-1, 1), dim(y)[1], replace = TRUE))
                }

                y.boot <- y.fit + res.p * res 
                boot.model <- try(fastplm.core(y = y.boot, x = x, ind = ind, 
                                               sfe.index = sfe.index, cfe.index = cfe.index, 
                                               se = 0, cl = NULL, core.num = core.num))

                if ('try-error' %in% class(boot.model)) {
                    return(rep(NA, p))
                } else {
                    return(c(boot.model$coefficients))
                }
            }

            if (parallel == TRUE) {
                boot.out <- foreach(j = 1:nboots, 
                                    .inorder = FALSE,
                                    .export = c("fastplm.core"),
                                    .packages = c("fastplm")
                                    ) %dopar% {
                                        return(one.boot())
                                    }
                for (j in 1:nboots) { 
                    boot.coef[, j] <- c(boot.out[[j]])
                }
            } else {
                for (i in 1:nboots) {
                    boot.coef[, i] <- one.boot()
                }
            }

            P_x <- as.matrix(apply(boot.coef, 1, get.pvalue))
            stderror <- as.matrix(apply(boot.coef, 1, sd, na.rm = TRUE))
            CI <- t(apply(boot.coef, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
            ## rewrite uncertainty estimates
            est.coefficients <- cbind(model$coefficients, stderror, P_x, CI)
            colnames(est.coefficients) <- c("Coef", "Std. Error", "P Value", "CI_lower", "CI_upper")
            model$est.coefficients <- est.coefficients
        
        } else {
            ## bootstrap refinement 
            ## obtain residuals under H0
            if (p == 1) {
                model.res <- fastplm.core(y = y, x = NULL, ind = ind, 
                          sfe.index = sfe.index, cfe.index = cfe.index, 
                          se = 0, cl = NULL, core.num = core.num)
            } else {
                model.res <- fastplm.core(y = y, x = as.matrix(x[, -pos]), 
                          ind = ind, sfe.index = sfe.index, cfe.index = cfe.index, 
                          se = 0, cl = NULL, core.num = core.num)
            }
            res <- model.res$residuals 
            y.fit <- y - res
            boot.wald <- rep(NA, nboots)
            ## bootstrap refinement: wald
            core.num <- ifelse(parallel == TRUE, 1, core.num)

            one.boot <- function(num = NULL) {
                ## wild boot by cluster 
                res.p.level <- sample(c(-1, 1), length(level.cluster), replace = TRUE)
                res.p <- c()
                for (j in 1:length(level.cluster)) {
                    res.p <- c(res.p, rep(res.p.level[j], level.count[j]))
                }
                res.p <- as.matrix(res.p)

                y.boot <- y.fit + res.p * res 
                boot.model <- try(fastplm.core(y = y.boot, x = x, ind = ind, 
                                               sfe.index = sfe.index, cfe.index = cfe.index, 
                                               se = 1, cl = cluster, core.num = core.num))
                if ('try-error' %in% class(boot.model)) {
                    return(NA)
                } else {
                    return(boot.model$est.coefficients[pos, "Z value"])
                }
            }

            if (parallel == TRUE) {
                boot.out <- foreach(j = 1:nboots, 
                                    .inorder = FALSE,
                                    .export = c("fastplm.core"),
                                    .packages = c("fastplm")
                                    ) %dopar% {
                                        return(one.boot())
                                    }
                boot.wald <- c(unlist(boot.out)) 
            } else {
                for (i in 1:nboots) {
                    boot.wald[i] <- one.boot()
                }
            }

            wald.all <- model$est.coefficients[pos, "Z value"] 
            refine.p <- sum(abs(boot.wald) > abs(wald.all), na.rm = 1) / (nboots - sum(is.na(boot.wald))) 
            wald.refinement <- matrix(c(wald.all, refine.p), 1, 2)
            colnames(wald.refinement) <- c("Z value", "P value")
            model$refinement <- list(wald.refinement = wald.refinement, boot.wald = boot.wald)

        }
    } else {    
        
        core.num <- ifelse(parallel == TRUE, 1, core.num)

        ## non-parametric bootstrap
        if (jackknife ==1 || is.null(cluster)) {
            ## non-parametric bootstrap
            jack.pos <- 1
            if (jackknife == 1) {
                nboots <- length(unique(ind[,1]))
                jack.pos <- as.numeric(as.factor(ind[,1]))
            }
            boot.coef <- matrix(NA, p, nboots)

            one.boot <- function(num = NULL) {
                if (is.null(num)) {
                    boot.id <- sample(1:dim(y)[1], dim(y)[1], replace = TRUE)
                } else {
                    boot.id <- (1:dim(y)[1])[which(jack.pos != num)]
                }
                

                boot.model <- try(fastplm.core(y = as.matrix(y[boot.id,]), 
                                               x = as.matrix(x[boot.id,]), 
                                               ind = as.matrix(ind[boot.id,]), 
                                               sfe.index = sfe.index, cfe.index = cfe.index, 
                                               se = 0, cl = NULL, core.num = core.num))
                if ('try-error' %in% class(boot.model)) {
                    return(rep(NA, p))
                } else {
                    return(c(boot.model$coefficients))
                }
            }

            boot.seq <- NULL
            if (jackknife == 1) {
                boot.seq <- unique(jack.pos)
            }
            
            if (parallel == TRUE) {
                boot.out <- foreach(j = 1:nboots, 
                                    .inorder = FALSE,
                                    .export = c("fastplm.core"),
                                    .packages = c("fastplm")
                                    ) %dopar% {
                                        return(one.boot(boot.seq[j]))
                                    }
                for (j in 1:nboots) { 
                    boot.coef[, j] <- c(boot.out[[j]])
                }
            } else {
                for (i in 1:nboots) {
                    boot.coef[, i] <- one.boot(boot.seq[i])
                }
            }

            if (jackknife == 0) {
                P_x <- as.matrix(apply(boot.coef, 1, get.pvalue))
                stderror <- as.matrix(apply(boot.coef, 1, sd, na.rm = TRUE))
                CI <- t(apply(boot.coef, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                ## rewrite uncertainty estimates
                est.coefficients <- cbind(model$coefficients, stderror, P_x, CI)
            } else {
                beta.j <- jackknifed(model$coefficients, boot.coef, 0.05)
                est.coefficients <- cbind(model$coefficients, beta.j$se, beta.j$P, beta.j$CI.l, beta.j$CI.u)

            }
            colnames(est.coefficients) <- c("Coef", "Std. Error", "P Value", "CI_lower", "CI_upper")
            model$est.coefficients <- est.coefficients
        
        } else {
            
            if (refinement == 0) {
                
                boot.coef <- matrix(NA, p, nboots)
                
                if (dim(cluster)[2] == 1) {
                    
                    ## one-way cluster 
                    raw.cluster <- as.numeric(as.factor(cluster))
                    level.cluster <- unique(raw.cluster)
                    ## level.count <- table(raw.cluster)
                    split.id <- split(1:dim(y)[1], raw.cluster)

                    one.boot <- function(num = NULL) {
                        level.id <- sample(level.cluster, length(level.cluster), replace = TRUE)
                        boot.id <- c()
                        for (j in 1:length(level.id)) {
                            boot.id <- c(boot.id, split.id[[level.id[j]]])
                        }
                        boot.model <- try(fastplm.core (y = as.matrix(y[boot.id,]), 
                                                        x = as.matrix(x[boot.id,]), 
                                                        ind = as.matrix(ind[boot.id,]), 
                                                        sfe.index = sfe.index, cfe.index = cfe.index, 
                                                        se = 0, cl = NULL, core.num = core.num))
                        if ('try-error' %in% class(boot.model)) {
                            return(rep(NA, p))
                        } else {
                            return(c(boot.model$coefficients))
                        }
                    }

                    if (parallel == TRUE) {
                        boot.out <- foreach(j = 1:nboots, 
                                            .inorder = FALSE,
                                            .export = c("fastplm.core"),
                                            .packages = c("fastplm")
                                            ) %dopar% {
                                                return(one.boot())
                                            }
                        for (j in 1:nboots) { 
                            boot.coef[, j] <- c(boot.out[[j]])
                        }
                    } else {
                        for (i in 1:nboots) {
                            boot.coef[, i] <- one.boot()
                        }
                    }

                    P_x <- as.matrix(apply(boot.coef, 1, get.pvalue))
                    stderror <- as.matrix(apply(boot.coef, 1, sd, na.rm = TRUE))
                    CI <- t(apply(boot.coef, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                    ## rewrite uncertainty estimates
                    est.coefficients <- cbind(model$coefficients, stderror, P_x, CI)
                    colnames(est.coefficients) <- c("Coef", "Std. Error", "P Value", "CI_lower", "CI_upper")
                    model$est.coefficients <- est.coefficients
                }
                
                else if (dim(cluster)[2] == 2) {
                    ## two-way cluster
                    boot.coef1 <- boot.coef2 <- boot.coef3 <- matrix(NA, p, nboots)
                    ## level 1
                    raw.cluster1 <- as.numeric(as.factor(cluster[,1]))
                    level.cluster1 <- unique(raw.cluster1)
                    ## level 2
                    raw.cluster2 <- as.numeric(as.factor(cluster[,2]))
                    level.cluster2 <- unique(raw.cluster2)
                    ## intersection
                    raw.cluster3 <- as.numeric(as.factor(paste(cluster[,1], "-:-", cluster[,2], sep = "")))
                    level.cluster3 <- unique(raw.cluster3)

                    split.id1 <- split(1:dim(y)[1], raw.cluster1)
                    split.id2 <- split(1:dim(y)[1], raw.cluster2)
                    split.id3 <- split(1:dim(y)[1], raw.cluster3)

                    one.boot <- function(num = NULL) {
                        ## level 1
                        level.id1 <- sample(level.cluster1, length(level.cluster1), replace = TRUE)
                        boot.id1 <- c()
                        for (j in 1:length(level.id1)) {
                            boot.id1 <- c(boot.id1, split.id1[[level.id1[j]]])
                        }

                        ## level 2
                        level.id2 <- sample(level.cluster2, length(level.cluster2), replace = TRUE)
                        boot.id2 <- c()
                        for (j in 1:length(level.id2)) {
                            boot.id2 <- c(boot.id2, split.id2[[level.id2[j]]])
                        }

                        ## intersection
                        level.id3 <- sample(level.cluster3, length(level.cluster3), replace = TRUE)
                        boot.id3 <- c()
                        for (j in 1:length(level.id3)) {
                            boot.id3 <- c(boot.id3, split.id3[[level.id3[j]]])
                        }

                        ## level 1
                        boot.model1 <- try(fastplm.core (y = as.matrix(y[boot.id1,]), 
                                                         x = as.matrix(x[boot.id1,]), 
                                                         ind = as.matrix(ind[boot.id1,]), 
                                                         sfe.index = sfe.index, cfe.index = cfe.index, 
                                                         se = 0, cl = NULL, core.num = core.num))
                        if ('try-error' %in% class(boot.model1)) {
                            oneboot.coef1 <- rep(NA, p)
                        } else {
                            oneboot.coef1 <- c(boot.model1$coefficients)
                        }

                        ## level 2
                        boot.model2 <- try(fastplm.core (y = as.matrix(y[boot.id2,]), 
                                                         x = as.matrix(x[boot.id2,]), 
                                                         ind = as.matrix(ind[boot.id2,]), 
                                                         sfe.index = sfe.index, cfe.index = cfe.index, 
                                                         se = 0, cl = NULL, core.num = core.num))
                        if ('try-error' %in% class(boot.model2)) {
                            oneboot.coef2 <- rep(NA, p)
                        } else {
                            oneboot.coef2 <- c(boot.model2$coefficients)
                        }

                        ## intersection
                        boot.model3 <- try(fastplm.core (y = as.matrix(y[boot.id3,]), 
                                                         x = as.matrix(x[boot.id3,]), 
                                                         ind = as.matrix(ind[boot.id3,]), 
                                                         sfe.index = sfe.index, cfe.index = cfe.index, 
                                                         se = 0, cl = NULL, core.num = core.num))
                        if ('try-error' %in% class(boot.model3)) {
                            oneboot.coef3 <- rep(NA, p)
                        } else {
                            oneboot.coef3 <- c(boot.model3$coefficients)
                        }

                        return(list(oneboot.coef1 = oneboot.coef1, 
                                    oneboot.coef2 = oneboot.coef2,
                                    oneboot.coef3 = oneboot.coef3))

                    }

                    if (parallel == TRUE) {
                        boot.out <- foreach(j = 1:nboots, 
                                            .inorder = FALSE,
                                            .export = c("fastplm.core"),
                                            .packages = c("fastplm")
                                            ) %dopar% {
                                                return(one.boot())
                                            }
                        for (j in 1:nboots) { 
                            boot.coef1[, j] <- c(boot.out[[j]]$oneboot.coef1)
                            boot.coef2[, j] <- c(boot.out[[j]]$oneboot.coef2)
                            boot.coef3[, j] <- c(boot.out[[j]]$oneboot.coef3)
                        }
                    } else {
                        for (i in 1:nboots) {
                            boot.sub <- one.boot()
                            boot.coef1[, i] <- c(boot.sub$oneboot.coef1)
                            boot.coef2[, i] <- c(boot.sub$oneboot.coef2)
                            boot.coef3[, i] <- c(boot.sub$oneboot.coef3)
                        }
                    }
                    
                    stderror1 <- as.matrix(apply(boot.coef1, 1, sd, na.rm = TRUE))
                    stderror2 <- as.matrix(apply(boot.coef2, 1, sd, na.rm = TRUE))
                    stderror3 <- as.matrix(apply(boot.coef3, 1, sd, na.rm = TRUE))
                    stderror <- stderror1 + stderror2 - stderror3
                    CI <- cbind(c(model$coefficients) - qnorm(0.975)*c(stderror), c(model$coefficients) + qnorm(0.975)*c(stderror))
                    Zx <- c(model$coefficients)/c(stderror)
                    P_z <- 2 * min(1 - pnorm(Zx), pnorm(Zx))
                    ## rewrite uncertainty estimates
                    est.coefficients <- cbind(model$coefficients, stderror, Zx, P_z, CI)
                    colnames(est.coefficients) <- c("Coef", "Std. Error", "Z value", "Pr(>|Z|)", "CI_lower", "CI_upper")
                    model$est.coefficients <- est.coefficients
                }
            } else {
                
                ## bootstrap refinement 
                boot.wald <- matrix(NA, p, nboots)
                ## one-way cluster 
                raw.cluster <- as.numeric(as.factor(cluster))
                level.cluster <- unique(raw.cluster)
                level.count <- table(raw.cluster)
                split.id <- split(1:dim(y)[1], raw.cluster)

                one.boot <- function(num = NULL) {
                    level.id <- sample(level.cluster, length(level.cluster), replace = TRUE)
                    boot.id <- c()
                    boot.cluster <- c()
                    for (j in 1:length(level.id)) {
                        boot.id <- c(boot.id, split.id[[level.id[j]]])
                        boot.cluster <- c(boot.cluster, rep(j, level.count[level.id[j]]))
                    }
                    boot.cluster <- as.matrix(boot.cluster)

                    boot.model <- try(fastplm.core (y = as.matrix(y[boot.id,]), 
                                                    x = as.matrix(x[boot.id,]), 
                                                    ind = as.matrix(ind[boot.id,]), 
                                                    sfe.index = sfe.index, cfe.index = cfe.index, 
                                                    se = 1, cl = boot.cluster, core.num = core.num))
                    if ('try-error' %in% class(boot.model)) {
                        return(rep(NA, p))
                    } else {
                        return(c(boot.model$coefficients - model$coefficients)/c(boot.model$est.coefficients[, "Std. Error"]))
                    }
                }

                if (parallel == TRUE) {
                    boot.out <- foreach(j = 1:nboots, 
                                        .inorder = FALSE,
                                        .export = c("fastplm.core"),
                                        .packages = c("fastplm")
                                        ) %dopar% {
                                            return(one.boot())
                                        }
                    for (j in 1:nboots) { 
                        boot.wald[, j] <- c(boot.out[[j]])
                    }
                } else {
                    for (i in 1:nboots) {
                        boot.wald[, i] <- one.boot()
                    }
                }

                wald.all <- model$est.coefficients[, "Z value"] 
                wald.CI <- t(apply(boot.wald, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                ## colnames(wald.CI) <- c("CI_lower", "CI_upper") 
                refine.p <- c()
                for (i in 1:length(wald.all)) {
                    refine.p <- c(refine.p, sum(abs(boot.wald[i]) > abs(wald.all[i]), na.rm = 1) / (nboots - sum(is.na(boot.wald[i,]))))
                }
                wald.refinement <- cbind(wald.all, wald.CI, refine.p)
                colnames(wald.refinement) <- c("Z value", "Lower 2.5 percentile", "Upper 97.5 percentile", "P value")
                model$refinement <- list(wald.refinement = wald.refinement, boot.wald = boot.wald)

            }
        }
        

    }

    return(model)
}



## jackknife se
jackknifed <- function(x,  ## ols estimates
                       y,
                       alpha) { ## sub-sample ols estimates) 

    p <- length(x)
    N <- dim(y)[2]  ## sample size
    if (N == 1) {
        y <- t(y)
        N <- dim(y)[2]
    }

    X <- matrix(rep(c(x), N), p, N) * N
    Y <- X - y * (N - 1)

    Yvar <- apply(Y, 1, var, na.rm = TRUE)
    vn <- N - apply(is.na(y), 1, sum) 

    Ysd <- sqrt(Yvar/vn)  ## jackknife se

    CI.l <- Ysd * qnorm(alpha/2) + c(x)
    CI.u <- Ysd * qnorm(1 - alpha/2) + c(x)

    ## wald test
    P <- NULL
    for (i in 1:p) {
        subz <- pnorm(c(x)[i]/Ysd[i])
        P <- c(P, 2 * min(1 - subz, subz))
    }

    ## P <- 2 * min(1 - pnorm(c(x)/Ysd), pnorm(c(x)/Ysd))

    out <- list(se = Ysd, CI.l = CI.l, CI.u = CI.u, P = P)

    return(out)
    
}