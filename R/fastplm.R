## generic function
fastplm <- function(formula = NULL, ## 
                    data = NULL, ## a matrix or a dataframe
                    index = NULL, ## index name    
                    y = NULL, ## outcome vector
                    x = NULL, ## covariate matrix
                    ind = NULL, ## indicator matrix
                    sfe = NULL, ## index for simple fe , a vector of integer
                    cfe = NULL, ## index for complex fe, a list of double, e.g. c(effect, influence)
                    PCA = TRUE,
                    se = 1, ## uncertainty estimates for covariates
                    vce = "robust", ## standard, robust, clustered, bootstrap
                    cluster = NULL, ## cluster name
                    wild = FALSE,
                    refinement = FALSE,
                    test_x = NULL,
                    parallel = FALSE, 
                    nboots = 200, ## bootstrap number
                    seed = NULL,
                    core.num = 1) {
  
    UseMethod("fastplm")
}

## formula
fastplm.formula <- function(formula = NULL, ## 
                            data = NULL, ## a matrix or a dataframe
                            index = NULL, ## index name    
                            y = NULL, ## outcome vector
                            x = NULL, ## covariate matrix
                            ind = NULL, ## indicator matrix
                            sfe = NULL, ## index for simple fe , a vector of integer
                            cfe = NULL, ## index for complex fe, a list of double, e.g. c(effect, influence)
                            PCA = TRUE,
                            se = 1, ## uncertainty estimates for covariates
                            vce = "robust", ## standard, robust, clustered
                            cluster = NULL, ## cluster name
                            wild = FALSE,
                            refinement = FALSE,
                            test_x = NULL,
                            parallel = FALSE,
                            nboots = 200, ## bootstrap number
                            seed = NULL,
                            core.num = 1) {
  
    ## receive a data frame with variable name
    ## if (!is.null(formula)) {
    ## indicator matrix
    varnames <- all.vars(formula)
    if (!is.null(cluster)) {
        if (class(cluster) != "character") {
            stop("\"cluster\" should be a character value.\n")
        }
    }
    data.name <- names(data)
    all.name <- unique(c(varnames,index,cluster))
    for (i in 1:length(all.name)) {
        if (!all.name[i] %in% data.name) {
            stop(paste("variable ", all.name[i], " is not in the dataset.\n", sep = ""))
        }
    }
    data <- data[,unique(c(varnames,index,cluster))]
    if (sum(is.na(data)) > 0) {
        data <- na.omit(data)
    }
    y <- as.matrix(data[, varnames[1]])
    x <- NULL
    if (length(varnames) >= 2) {
        x <- as.matrix(data[, varnames[2:length(varnames)]])
    }
    ind <- as.matrix(data[, index])

    ## save cluster name
    cluster.level <- cluster
    if (!is.null(cluster)) {
        cluster <- as.matrix(data[, cluster])
    }
    ## data <- as.matrix(data[, varnames])
    ## colnames(data) <- varnames
    ## } else {
        ## receive a matrix 
    ##     if (!is.null(data)) {
    ##         if (sum(is.na(data)) > 0) {
    ##             data <- na.omit(data)
    ##         }
    ##         p <- dim(data)[2]
    ##         y <- as.matrix(data[, 1])
    ##         x <- NULL
    ##         if (p >= 2) {
    ##             x <- as.matrix(data[, 2:p])
    ##         } 
    ##     }
    ## }

    out <- fastplm.default(formula = formula, data = data,
                           index = index, y = y, 
                           x = x, ind = ind, 
                           sfe = sfe, cfe = cfe,
                           PCA = PCA, se = se, vce = vce,
                           cluster = cluster, 
                           wild = wild, refinement = refinement, 
                           test_x = test_x, 
                           parallel = parallel, nboots = nboots, 
                           seed = seed, core.num = core.num)

    out <- c(out, list(cluster.level = cluster.level))
    class(out) <- "fastplm"

    return(out)

}

## default function
fastplm.default <- function(formula = NULL, ## 
                            data = NULL, ## a matrix or a dataframe
                            index = NULL, ## index name    
                            y = NULL, ## outcome vector
                            x = NULL, ## covariate matrix
                            ind = NULL, ## indicator matrix
                            sfe = NULL, ## index for simple fe , a vector of integer
                            cfe = NULL, ## index for complex fe, a list of double, e.g. c(effect, influence)
                            PCA = TRUE,
                            se = 1, ## uncertainty estimates for covariates
                            vce = "robust", ## standard, robust, clustered
                            cluster = NULL, ## cluster in default is a column vector
                            wild = FALSE,
                            refinement = FALSE,
                            test_x = NULL,  
                            parallel = FALSE,
                            nboots = 200, ## bootstrap number
                            seed = NULL,
                            core.num = 1)  {
  
    ## fe <- NULL
    ## inds <- create.indicators(ind)
    ## CHECK.INPUT(inds, "inds", "indicators")

    ## data <- CHECK.INIT.INPUT(data, "data, x, y",
    ##   init = function() cbind(y, x))
    ## fe  <- CHECK.INIT.INPUT(fe, "fe", arg.class = "fixed.effects",
    ##   init = function() create.fixed.effects(inds))

    ## ----- data ------ ## 
    ## if (is.null(y)) {
    ##     y <- as.matrix(data[, 1])
    ##     p <- dim(data)[2] - 1
    ##     if (p > 0) {
    ##         x <- as.matrix(data[, 2:(p+1)])
    ##     } else {
    ##         x <- NULL
    ##     }
    ## } else {
    ##     p <- ifelse(is.null(x), 0, dim(x)[2])
    ## }

    ## receive a matrix 
    if (is.null(formula) && !is.null(data)) {
        if (sum(is.na(data)) > 0) {
            data <- na.omit(data)
        }
        p <- dim(data)[2]
        y <- as.matrix(data[, 1])
        x <- NULL
        if (p >= 2) {
            x <- as.matrix(data[, 2:p])
        } 
    }

    ## number of covariates
    p <- ifelse(is.null(x), 0, dim(x)[2])
    #if (sum(is.na(y)) > 0) {
    #    stop("Missing values in dependent variable.\n")
    #}
    #if (p > 0) {
    #    if (sum(is.na(x)) > 0) {
    #        stop("Missing values in covariates.\n")
    #    }
    #}

    ## index 
    if (is.null(ind)) {
        stop("No indicators.\n")
    } else {
        if (class(ind) != "matrix") {
            stop("\"ind\" must be a matrix.\n")
        }
    }

    robust <- FALSE
    if (vce == "cl") {
        vce <- "clustered"
    }
    if (vce == "boot") {
        vce <- "bootstrap"
    }

    bootstrap <- 0
    if (se == 1) {
        if (p == 0) {
            cat("No covariates.\n")
            refinement <- 0
            se <- 0
            bootstrap <- 0
        }
        if (!vce %in% c("standard", "clustered", "robust", "bootstrap")) {
            stop("Choose \" vce \" from c(\" standard \", \" robust \", \" clustered \", \" bootstrap \").\n")
        } else {
            if (vce == "standard") {
                cat("Please consider clustering the standard errors or using block bootstraps.\n")
            } 
            else if (vce == "robust") {
                robust <- 1
            }
            else if (vce == "bootstrap") {
                bootstrap <- 1
            }
        }
    } else {
        refinement <- 0
    }

    if (refinement == 1) {
        bootstrap <- 1
    }

    n.cluster <- NULL
    if (!is.null(cluster)) {
        if (dim(cluster)[2] > 2) {
            stop("The number of clustered variables should not be larger than 2.\n")
        }
        else if (dim(cluster)[2] == 2) {
            if (wild == TRUE) {
                stop("For cluster wild bootstrap, please only specify one cluster variable.\n")
            }
            if (refinement == TRUE) {
                stop("For cluster bootstrap refinement, please only specify one cluster variable.\n")

            }
            n.cluster <- c(length(unique(cluster[,1])), length(unique(cluster[,2]))) 
        }
        else {
            n.cluster <- c(length(unique(c(cluster))))
        }
        if (vce %in% c("standard", "robust")) {
            vce <- "clustered"
        }
    } else {
        if (refinement == TRUE) {
            cat("No cluster variable. Cluster bootstrap refinement will not be performed.\n")
            refinement <- 0
        }
    }

    if (is.null(core.num) == FALSE) {
        if (core.num <= 0) {
            stop("\"core.num\" option misspecified. Try, for example, core.num = 1.")
        }
    } else {
        core.num <- detectCores()
    }

    pos <- NULL
    if (refinement == 1 && wild == TRUE) {
        if (is.null(test_x)) {
            stop("For wild bootstrap refinement, variables of interest should be specified.\n")
        }
        if (length(test_x) >= 2) {
            stop("Please specify only one variable of interest.\n")
        }
        if (class(test_x) == "character") {
            pos <- which(colnames(data)[2:(p+1)] == test_x)
        } else {
            pos <- test_x
        }
    }

    ## ------ fe ------- ## 
    num.sfe.index <- num.cfe.index <- NULL
    sfe.index <- sfe
    cfe.index <- cfe

    ## 1.simple fe index
    if (!is.null(sfe.index)) {
        if (class(index) == "character" && class(sfe.index) == "character") {
            num.sfe.index <- sapply(1:length(sfe.index), function(i) which(index == sfe.index[i]))
            sfe.index <- num.sfe.index
        }
    }

    ## 2.complex fe index
    raw_pc.ref <- NULL ## raw influence matrix
    if (!is.null(cfe.index)) {    
        ## length
        cfe.length <- length(cfe.index)
        raw.influence.index <- rep(NA, cfe.length)

        ## position in index
        if (class(index) == "character") {
            for (i in 1:cfe.length) {
                if (class(cfe.index[[i]]) == "character") {
                    num.cfe.index <- sapply(1:2, function(j) which(index == cfe.index[[i]][j]) )
                    cfe.index[[i]] <- num.cfe.index
                }
            }
        }
        
        ## pca for correlated factors
        if (PCA == TRUE) {
            raw.influence.index <- unique(sapply(1:cfe.length, function(i)cfe.index[[i]][2]))
            raw.influence <- as.matrix(ind[,raw.influence.index])

            pc.influence <- prcomp(raw.influence)$x ## pca
            ## qr.influence <- qr(raw.influence)
            ## pc.influence <- qr.Q(qr.influence)

            ## reference for raw and principal com influence
            get.ref <- function(mat) {return(as.matrix(mat[!duplicated(mat),],,2))}
            raw_pc.ref <- lapply(1:length(raw.influence.index), function(i)get.ref(cbind(raw.influence[,i], pc.influence[,i])))

            ## change ind
            ind[,raw.influence.index] <- pc.influence
        }
    }

    ## variance type 
    variance.type <- NULL
    refinement.type <- NULL
    if (se == 1) {
        if (bootstrap == 0) {
            if (is.null(cluster)) {
                if (vce == "standard") {
                    variance.type <- "Standard"
                } 
                else if (vce == "robust") {
                    variance.type <- "Robust"
                } 
                
            } else {
                variance.type <- "Clustered"
            }
        } else { ## bootstrap
            if (is.null(cluster)) {
                if (wild == 0) {
                    variance.type <- "Non-parametric Bootstrap"
                } else {
                    variance.type <- "Wild Bootstrap"
                }
            } else {
                if (wild == 0) {
                    if (refinement == 0) {
                        variance.type <- "Cluster Bootstrap"
                    } else {
                        variance.type <- "Clustered"
                        refinement.type <- "Percentile-t Bootstrap Refinement"
                    }
                } else {
                    if (refinement == 0) {
                        variance.type <- "Wild Cluster Bootstrap"
                    } else {
                        variance.type <- "Clustered"
                        refinement.type <- "Wild Cluster Bootstrap Refinement"
                    }
                }
            }
        }
    }

    if (bootstrap == 0) {
      
        model <- fastplm.core(y = y, x = x, ind = ind, 
                              sfe.index = sfe.index, cfe.index = cfe.index, 
                              se = se, robust = robust,
                              cl = cluster, core.num = core.num)
    } else {
      
        if (parallel == TRUE) {
            para.clusters <- makeCluster(core.num)
            registerDoParallel(para.clusters)
            cat("Parallel computing ...\n")
        }
        model <- fastplm.boot(seed = seed, y = y, x = x, ind = ind, 
                              sfe.index = sfe.index, cfe.index = cfe.index, 
                              cluster = cluster, parallel = parallel,
                              wild = wild, 
                              refinement = refinement, pos = pos,  
                              nboots = nboots, core.num = core.num)
        if (parallel == TRUE) {
            stopCluster(para.clusters)
        }
    }

    ## complex fixed effect
    model$cfe.index <- cfe.index ## a numeric list
    model$raw_pc.ref <- raw_pc.ref ## a ref list for complex fe, used for prediction
    model$PCA <- PCA

    if (p > 0 & !is.null(colnames(data))) {
        rownames(model$coefficients) <- colnames(data)[2:(p+1)]
        if (se == 1) {
            rownames(model$est.coefficients) <- colnames(data)[2:(p+1)]
        }
        if (!is.null(model$refinement)) {
            if(wild == FALSE) {
                rownames(model$refinement$wald.refinement) <- colnames(data)[2:(p+1)]
            } else {
                rownames(model$refinement$wald.refinement) <- colnames(data)[2:(p+1)][pos]
            }
        }
    }

    model <- c(model, list(call = match.call(),
                           n.cluster = n.cluster,
                           variance.type = variance.type,
                           refinement.type = refinement.type))
    class(model) <- "fastplm"
    model
}

name.fe.model <- function(model, inds, fe) {
  ## CHECK.INPUT(model, "model", "feModel")
  
  coefs <- model$sfe.coefs
  names(coefs) <- fe$sfe.names
  for (i in SEQ(1, length(fe$sfes)))
    rownames(coefs[[i]]) <- inds$levels[[fe$sfes[i]]]
  model$sfe.coefs <- coefs

  coefs <- model$cfe.coefs
  names(coefs) <- fe$cfe.names
  for (i in SEQ(1, length(fe$cfes)))
    rownames(coefs[[i]]) <- inds$levels[[fe$cfe.effs[i]]]
  model$cfe.coefs <- coefs

  model
}

## ------- method -------- ##
## predict
predict.fastplm <- function(object, data = NULL, x = NULL, ind, ...) {
    model <- object
    CHECK.INPUT(model, "model", "fastplm")
    if (!is.null(data)) { ## receive a dataframe
        coef <- model$coefficients
        if (is.null(coef)) {
            x <- NULL
        } else {
            x.name <- rownames(coef)
            if (is.null(x.name)) {
                stop("No covariate names in the model.")
            } else {
                if (sum(x.name %in% names(data)) < length(x.name)) {
                    stop("The dataset doesn\'t contain some covariates in the model.")
                } else {
                    x <- as.matrix(data[, x.name])
                }
            }
        }

        ind.name <- model$inds$effect.names
        if (sum(ind.name %in% names(data)) < length(ind.name)) {
            stop("The dataset doesn\'t contain some indicators in the model.")
        } else {
            ind <- as.matrix(data[, ind.name])
        }
    } else { ## receive matrix

        if (!is.null(x)) {
            CHECK.INPUT(x, "x", "matrix")
            coef <- model$coefficients 
            if (is.null(coef)) {
                stop("There is not any covariate in the model.")
            } else {
                if (dim(coef)[1] != dim(x)[2]) {
                    stop("The number of covariates should be the same as that in the model.")
                }
            }
        }

        if (dim(ind)[2] != length(model$inds$effect.names)) {
            stop("The number of indicators should be the same as that in the model.")
        }

    }
    
    if (!is.null(model$cfe.index) && model$PCA == TRUE) {
        sub.raw.influence.index <- unique(sapply(1:length(model$cfe.index), function(i) model$cfe.index[[i]][2]))
        sub.raw.influence <- as.matrix(ind[,sub.raw.influence.index])
        sub.pc.influence <- matrix(NA, dim(sub.raw.influence)[1], dim(sub.raw.influence)[2])
        raw_pc.ref <- model$raw_pc.ref
        for (i in 1:dim(sub.raw.influence)[2]) {
            sub.ref <- raw_pc.ref[[i]]
            for (j in 1:dim(sub.ref)[1]) {
                sub.pc.influence[which(sub.raw.influence[,i] == sub.ref[j,1]),i] <- sub.ref[j,2]
            }
        }
        if (sum(is.na(c(sub.pc.influence))) > 0) {
            stop("There may be some unspecific levels within some fixed effects in the model.\n")
        }
        ind[,sub.raw.influence.index] <- sub.pc.influence
    }

    inds <- create.subindicators(sub.inds = ind, model = model)
    ## CHECK.INPUT(inds, "inds", "sub.indicators")

    if (!is.null(x)) {
        ASSERT.MATRIX.DIM(x, "x", length(model$coefficients), is.width = TRUE)
    }

    if (inds$parent != model$inds$uid)
        stop(ERR.FE.predict.invalid.inds())
    
    if (!is.null(x)) {
        y <- x %*% model$coefficients + model$intercept
    } else {
        y <- as.matrix(rep(model$intercept, dim(ind)[1]))
    }
    fe <- model$fe

    for (col in fe$sfes) {
        effs <- model$sfe.coefs[[col]]
        sum  <- sapply(SEQ(1, dim(ind)[1]), function(i) effs[inds$inds[i, col]])
        y    <- y + sum
    }

    for (i in SEQ(1, length(fe$cfe.effs))) {
        eff.col <- fe$cfe.effs[i]
        inf.col <- fe$cfe.infs[i]
        effs    <- model$cfe.coefs[[i]]
        weights <- fe$weights[[i]]

        sum <- sapply(SEQ(1, dim(ind)[1]), function(i)
            effs[inds$inds[i, eff.col], ] %*% weights[, inds$inds[i, inf.col]])
        y   <- y + sum
    }
    y
}

## Print
summary.fastplm <- function(object,  
                            ...) {
    
    ## cat("Call:\n")
    ## print(object$call, digits = 4)

    if (is.null(object$coefficients)) {
        cat("\n(No coefficients)\n")
    } else {
        if (is.null(object$est.coefficients)) {
            cat("\n")
            coefficients <- object$coefficients
            colnames(coefficients) <- "Coef"
            print(round(coefficients, 3))
        } else {
            cat("\n")
            print(round(object$est.coefficients, 3))
            cat("\n---\n\n")
            if (is.null(object$F_statistic)) {
                cat("Residual standard error: ", sprintf("%.3f",object$RMSE), sep = "")
            } else {
                cat("Residual standard error: ", sprintf("%.3f",object$RMSE), " on ", sprintf("%.f", object$df), " degrees of freedom", sep = "")
            }
            if (!is.null(object$R2)) {
                cat("\nMultiple R-squared(full model): ", sprintf("%.3f", object$R2), "\tAdjusted R-squared: ", sprintf("%.3f", object$Adj_R2), sep = "")
            }
            if (!is.null(object$projR2)) {
                cat("\nMultiple R-squared(proj model): ", sprintf("%.3f", object$projR2), "\tAdjusted R-squared: ", sprintf("%.3f", object$projAdj_R2), sep = "")
            }
            if (!is.null(object$F_statistic)) {
                if (object$variance.type == "Robust") {
                    cat("\nWald F-statistic: ", sprintf("%.3f", object$F_statistic[1,1]), " on ", sprintf("%.f", object$F_statistic[1,2]), " and ", sprintf("%.f", object$F_statistic[1,3]), " DF, p-value: < ", sprintf("%.3f", object$F_statistic[1,4]), sep = "")
                } else {
                    cat("\nF-statistic(full model): ", sprintf("%.3f", object$F_statistic[1,1]), " on ", sprintf("%.f", object$F_statistic[1,2]), " and ", sprintf("%.f", object$F_statistic[1,3]), " DF, p-value: < ", sprintf("%.3f", object$F_statistic[1,4]), sep = "")
                }
                
            }
            if (!is.null(object$proj_F_statistic)) {
                cat("\nF-statistic(proj model): ", sprintf("%.3f", object$proj_F_statistic[1,1]), " on ", sprintf("%.f", object$proj_F_statistic[1,2]), " and ", sprintf("%.f", object$proj_F_statistic[1,3]), " DF, p-value: < ", sprintf("%.3f", object$proj_F_statistic[1,4]), sep = "")
            }

            cat("\n\n---\n\n") 
            cat("Variance Type: ", object$variance.type)
            cat("\n")
            if (!is.null(object$n.cluster)) {
                clname <- object$cluster.level[1]
                if (is.null(clname)) {
                    if(length(object$n.cluster) == 1) {
                        clname <- "Cluster variable"
                    } else {
                        clname <- "Cluster variable 1"
                    }
                    cat("Number of clusters: ", object$n.cluster[1], " in ", clname, ".", sep = "")
                } else {
                    cat("Number of clusters: ", object$n.cluster[1], " in variable \"", clname, "\".", sep = "")
                }
                if (length(object$n.cluster) > 1) {
                    clname <- object$cluster.level[2]
                    if (is.null(clname)) {
                        clname <- "Cluster variable 2"
                        cat("\nNumber of clusters: ", object$n.cluster[2], " in ", clname, ".", sep = "")
                    } else {
                        cat("\nNumber of clusters: ", object$n.cluster[2], " in variable \"", clname, "\".", sep = "")
                    }
                }
            }
            if ((length(object$fe$sfe.names) >= 3 || length(object$fe$cfe.names) >= 2) && dim(object$est.coefficients)[2] > 5) { ## rule out bootstrap
                cat("\n* Degree of freedom may be too conservative due to redundant parameters.")
            }
        }
    }

    if (!is.null(object$refinement$wald.refinement)) {
        cat("\n\n---\n\n")
        ## cat("\nBootstrap refinement:\n")
        cat(object$refinement.type)
        cat("\n\n")
        print(round(object$refinement$wald.refinement, 3))
    }  
}

ERR.FE.predict.invalid.inds <- function() {
    sprintf("The given indicators do not belong to the indicators in the model.")
}
