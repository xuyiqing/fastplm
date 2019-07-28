## core function
fastplm.core <- function(y = NULL, ## outcome vector
                         x = NULL, ## covariate matrix
                         ind = NULL, ## indicator matrix
                         se = 0,
                         robust = FALSE,
                         cl = NULL,
                         sfe.index = NULL, ## index for simple fe , a vector of integer
                         cfe.index = NULL, ## index for complex fe, a list of double, e.g. c(effect, influence)
                         core.num = 1)  {
  
    ## fe <- NULL
    ## inds <- create.indicators(ind)
    ## CHECK.INPUT(inds, "inds", "indicators")

    ## data <- CHECK.INIT.INPUT(data, "data, x, y",
    ##   init = function() cbind(y, x))
    ## fe  <- CHECK.INIT.INPUT(fe, "fe", arg.class = "fixed.effects",
    ##   init = function() create.fixed.effects(inds))

    ## ----- data ------ ## 
    ## if (is.null(data)) {
    data <- cbind(y, x)
    ## }

    ## create indicators
    inds <- create.indicators(ind) ## an indicator class object

    ## ------ fe ------- ##
    fe <- cfes <- sub.cfe.index <- num.cfe.index <- NULL

    ## create complex fe
    gcfe <- function(i) {
        sub.cfe.index <- cfe.index[[i]]
        inf.weight <- as.numeric(inds$levels[[sub.cfe.index[2]]])
        cfe.sub <- create.complex.effect(inds, sub.cfe.index[1], sub.cfe.index[2], t(as.matrix(inf.weight)))
        return(cfe.sub)
    }

    ## complex fe
    if (!is.null(cfe.index)) {    
        ## length
        cfe.length <- length(cfe.index)
        cfes <- lapply(1:cfe.length, gcfe)
    }

    ## create an fe class object
    fe <- create.fixed.effects(inds, sfes = sfe.index, cfes = cfes)

    ## ---------- solve the model ------------ ##
    model <- SolveFixedEffects(data, fe$ptr, core.num)
    model <- name.fe.model(model, inds, fe)
    model$inds <- inds
    model$fe <- fe

    ## ------------ inference ---------------- ##
    if (se == 1) {
        
        dx <- model$demeaned$x
        res <- model$residuals
        gtot <- 0 ## total loss of degree of freedom

        ## simple fe levels: 1 redundant parameter for each fe 
        if (is.null(sfe.index)) {
            gtot <- length(unlist(model$inds$levels)) - length(model$inds$levels)
        } else {
            ## gtot <- 0
            for (i in 1:length(sfe.index)) {
                gtot <- gtot + length(model$inds$levels[[sfe.index[i]]]) - 1 ## 1 reference level
            }
        }

        ## we test nest among the first two sets of fixed effects
        adj.dof <- 0
        need.adj.dof <- FALSE
        if (is.null(sfe.index)) {
            if (dim(ind)[2] >= 2) {
                need.adj.dof <- TRUE
            }
        } else {
            if (length(sfe.index) >= 2) {
                need.adj.dof <- TRUE
            }
        }

        if (need.adj.dof) {
            if (is.null(sfe.index)) {
                pos1 <- 1
                pos2 <- 2
            } else {
                pos1 <- sfe.index[1]
                pos2 <- sfe.index[2]
            }
            level1 <- ind[, pos1]
            level2 <- ind[, pos2]
            c.level1 <- length(inds$levels[[pos1]])
            c.level2 <- length(inds$levels[[pos2]])
            
            if (c.level1 > c.level2) {
                fe.level <- table(unlist(tapply(level1, level2, unique)))
            } else {
                fe.level <- table(unlist(tapply(level2, level1, unique)))
            }
            if (max(fe.level) == 1) {
                adj.dof <- min(c.level1, c.level2) - 1
            }
        }
        
        gtot <- gtot - adj.dof

        ## complex fe levels: 1 redundant parameter for each effect variable: cfe (effect, influence)
        if (!is.null(cfe.index)) {
            ## for (i in 1:length(cfe.index)) {
            ##     if (is.null(sfe.index)) { ## in this case, simple fe contains redundant parameters 
            ##         gtot <- gtot + length(model$inds$levels[[cfe.index[[i]][1]]]) - 1
            ##     } else {
            ##         if (cfe.index[[i]][2]%in%sfe.index) {
            ##             gtot <- gtot + length(model$inds$levels[[cfe.index[[i]][1]]]) - 1
            ##         } else {
            ##             gtot <- gtot + length(model$inds$levels[[cfe.index[[i]][1]]]) 
            ##         }
            ##     }
            ## }
            sum.cfe.coef <- NULL 
            for (i in 1:length(cfe.index)) {
                sum.cfe.coef <- sum(c(model$cfe.coefs[[i]])) 
                if (sum.cfe.coef <= 1e-7) {
                    gtot <- gtot + length(model$inds$levels[[cfe.index[[i]][1]]]) - 1
                } else {
                    gtot <- gtot + length(model$inds$levels[[cfe.index[[i]][1]]])
                }
            }
        }

        ## intercept
        gtot <- gtot + 1

        ## (x^{\prime}x)^{-1}
        invx <- solve(t(dx)%*%dx)

        ## degree of freedom
        df <- dim(x)[1] - gtot - dim(x)[2]
        ## estimated sigma hat
        sig2 <- sum(t(res)%*%res)/df
        RMSE <- sqrt(sig2)

        R2 <- 1 - sum(t(res)%*%res)/sum((data[, 1] - mean(data[, 1]))^2)
        Adj_R2 <- 1 - (1 - R2) * (dim(x)[1] - 1) / df

        projR2 <- 1 - sum(t(res)%*%res)/sum((c(model$demeaned$y) - mean(model$demeaned$y))^2)
        projAdj_R2 <- 1 - (1 - projR2) * (dim(x)[1] - 1) / df
      
        if (is.null(cl)) { ## ols se
        
            if (robust == FALSE) {
                stderror <- as.matrix(sqrt(sig2*diag(invx)))
            } else {
                meat <- dx * matrix(rep(c(res), dim(x)[2]), dim(x)[1], dim(x)[2])
                meat <- t(meat) %*% meat
                vcov <- invx %*% meat %*% invx 
                vcov <- dim(x)[1] / df * vcov

                stderror <- c()
                for (i in 1:dim(x)[2]) {
                    stderror <- c(stderror, sqrt(vcov[i, i]))
                }
                stderror <- as.matrix(stderror) 

            }
            ## stderror <- as.matrix(sqrt(sig2*diag(invx)))
            Tx <- c(model$coefficients)/c(stderror)
            P_t <- 2 * min(1 - pt(Tx, df), pt(Tx, df))
            CI <- cbind(c(model$coefficients) - qnorm(0.975)*c(stderror), c(model$coefficients) + qnorm(0.975)*c(stderror))

            model$df <- df

            ## F test
            if (robust == FALSE) {
                ## full model
                fF <- R2/(1 - R2) * df/(dim(x)[1] - df - 1)
                ## projection model
                projF <- projR2/(1 - projR2) * df/dim(x)[2]

                F_statistic <- cbind(fF, dim(x)[1] - df - 1, df, 1 - pf(fF, dim(x)[1] - df - 1, df))
                proj_F_statistic <- cbind(projF, dim(x)[2], df, 1 - pf(fF, dim(x)[2], df))

                colnames(F_statistic) <- c("F_statistic", "df1", "df2", "P value")
                colnames(proj_F_statistic) <- c("F_statistic", "df1", "df2", "P value")

                model$F_statistic <- F_statistic
                model$proj_F_statistic <- proj_F_statistic
            } else {
                ## robust wald test 
                fF <- t(model$coefficients) %*% solve(vcov) %*% model$coefficients
                fF <- c(fF) / dim(x)[2]

                F_statistic <- cbind(fF, dim(x)[2], df, 1 - pf(fF, dim(x)[2], df))
                colnames(F_statistic) <- c("F_statistic", "df1", "df2", "P value")

                model$F_statistic <- F_statistic

            }

            est.coefficients <- cbind(model$coefficients, stderror, Tx, P_t, CI)
            colnames(est.coefficients) <- c("Coef", "Std. Error", "t value", "Pr(>|t|)", "CI_lower", "CI_upper")


        } else { ## clustered robust se: wald test

            if (dim(cl)[2] == 1) {
                stderror <- cluster.se(cl = cl, x = x, res = res, dx =dx, 
                                       sfe.index = sfe.index, 
                                       cfe.index = cfe.index,
                                       ind = ind, df = df, invx = invx, model = model)
            } 
            else if (dim(cl)[2] == 2) {
                
                stderror1 <- cluster.se(cl = as.matrix(cl[,1]), x = x, res = res, 
                                        dx =dx, sfe.index = sfe.index, 
                                        cfe.index = cfe.index, ind = ind,
                                        df = df, invx = invx, model = model)

                stderror2 <- cluster.se(cl = as.matrix(cl[,2]), x = x, res = res, 
                                        dx =dx, sfe.index = sfe.index, 
                                        cfe.index = cfe.index, ind = ind,
                                        df = df, invx = invx, model = model)

                cl_inter <- as.numeric(as.factor(paste(cl[,1], "-:-", cl[,2], sep = "")))

                stderror3 <- cluster.se(cl = as.matrix(cl_inter), x = x, res = res, 
                                        dx =dx, sfe.index = sfe.index, 
                                        cfe.index = cfe.index, ind = ind,
                                        df = df, invx = invx, model = model)

                stderror <- sqrt(stderror1^2 + stderror2^2 - stderror3^2)
            }

            ## wald test 
            Zx <- c(model$coefficients)/c(stderror)
            P_z <- 2 * min(1 - pnorm(Zx), pnorm(Zx))
            CI <- cbind(c(model$coefficients) - qnorm(0.975)*c(stderror), c(model$coefficients) + qnorm(0.975)*c(stderror))
            est.coefficients <- cbind(model$coefficients, stderror, Zx, P_z, CI)
            colnames(est.coefficients) <- c("Coef", "Std. Error", "Z value", "Pr(>|Z|)", "CI_lower", "CI_upper")
        }      

        model$R2 <- R2
        model$Adj_R2 <- Adj_R2
        model$projR2 <- projR2
        model$projAdj_R2 <- projAdj_R2
        model$RMSE <- RMSE

        model$est.coefficients <- est.coefficients
      
    }

    ## class(model) <- "feModel"
    ## class(model) <- "fastplm"
    model
}

## subfunction for clustered se
cluster.se <- function(cl, 
                       x,
                       res,
                       dx,
                       sfe.index,
                       cfe.index,
                       ind,
                       df,
                       invx,
                       model) {

    ## replicate degree of freedom 
    df.cl <- df
    
    raw.cl <- as.numeric(as.factor(cl))
    level.cl <- unique(raw.cl)
    ## raw.cl <- cl
    ## level.cl <- unique(cl)
        
    meat <- matrix(0, dim(x)[2], dim(x)[2])
    for (i in 1:length(level.cl)) {
        sub.id <- which(raw.cl == level.cl[i])
        if (length(sub.id) == 1) {
            hmeat <- res[sub.id] * dx[sub.id,]
            meat <- meat + hmeat^2
        } else {
            hmeat <- t(as.matrix(res[sub.id,])) %*% as.matrix(dx[sub.id,])
            meat <- meat + t(hmeat) %*% hmeat
        }
    }

    ## check cluster for the first 2 levels of fixed effects
    cluster.adj <- 0
    contain.cl <- NULL
    if (is.null(sfe.index)) {
        contain.cl <- apply(ind, 2, function(vec) sum(vec == cl)) == dim(x)[1]
    } else {
        contain.cl <- apply(as.matrix(ind[, sfe.index]), 2, function(vec) sum(vec == cl)) == dim(x)[1]
    }
    ## adjust q formula for finite sample
    if (max(contain.cl) == 1) {
        df.cl <- df.cl + length(level.cl) - 1
    }

    if (length(contain.cl) >= 2) {
        if (max(contain.cl[c(1,2)]) == 1) {
            cluster.adj <- 1 
        }
    }

    ## check within cluster fe
    if (is.null(sfe.index)) {
        for (i in 1:dim(ind)[2]) {
            ## if (sum(c(cl)==c(ind[,i])) != length(cl)) {
                ## fe.level <- unlist(tapply(c(cl), ind[,i], unique))
                ## if (length(unique(fe.level)) == length(fe.level)) {
                ##   df <- df + length(fe.level) - 1
                ## }
                fe.level <- as.numeric(table(unlist(tapply(ind[,i], c(cl), unique))))
                if (length(fe.level) != length(level.cl)) {
                    if (sum(fe.level == 1) == length(fe.level)) {
                        if ( (i <= 2) && (cluster.adj == 1)) {
                            df.cl <- df.cl + length(fe.level) - length(level.cl)
                        } else {
                            df.cl <- df.cl + length(fe.level) - 1
                        }
                    }
                }   
            ## }
        }
    } else {
        for (i in 1:length(sfe.index)) {
            ## if (sum(c(cl)==c(ind[,sfe.index[i]])) != length(cl)) {
                ## fe.level <- unlist(tapply(c(cl), ind[,sfe.index[i]], unique))
                ## if (length(unique(fe.level)) == length(fe.level)) {
                ##   df <- df + length(fe.level) - 1
                ## }
                fe.level <- unlist(tapply(ind[,sfe.index[i]], c(cl), unique))
                if (length(fe.level) != length(level.cl)) {
                    if (sum(fe.level == 1) == length(fe.level)) {
                        if ( (i <= 2) && (cluster.adj == 1)) {
                            df.cl <- df.cl + length(fe.level) - length(level.cl)
                        } else {
                            df.cl <- df.cl + length(fe.level) - 1
                        }
                    }
                }
            ## }
        }
    }

    if (!is.null(cfe.index)) {
        for (i in 1:length(cfe.index)) {
            ## if (sum(c(cl)==c(ind[,cfe.index[[i]][1]])) != length(cl)) {
                ## fe.level <- unlist(tapply(c(cl), ind[,sfe.index[i]], unique))
                ## if (length(unique(fe.level)) == length(fe.level)) {
                ##   df <- df + length(fe.level) - 1
                ## }
                fe.level <- unlist(tapply(ind[,cfe.index[[i]][1]]), c(cl), unique)
                if (length(fe.level) != length(level.cl)) {
                    if (sum(fe.level == 1) == length(fe.level)) {
                        df.cl <- df.cl + length(fe.level) - 1
                    }
                }
            ## }
        }
    }
    q <- length(level.cl)/(length(level.cl) - 1) * (dim(x)[1] - 1)/(df.cl + 1)
    ## q <- length(level.cl)/(length(level.cl)-1)
    stderror <- sqrt(q) * as.matrix(sqrt(diag(invx %*% meat %*% invx)))
    #Tx <- c(model$coefficients)/c(stderror)
    #P_t <- 2 * (1 - pt(Tx, length(level.cl) - 1))
    #CI <- cbind(c(model$coefficients) - qnorm(0.975)*c(stderror), c(model$coefficients) + qnorm(0.975)*c(stderror))
    return(stderror)
}


