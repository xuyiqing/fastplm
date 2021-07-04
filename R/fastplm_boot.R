fastplm.boot <- function(seed, 
                         y, 
                         x, 
                         ind, 
                         sfe.index, 
                         cfe.index, 
                         cluster,
                         robust = FALSE,
                         parallel = FALSE,
                         wild = FALSE, 
                         jackknife = FALSE,
                         refinement = FALSE, 
                         pos = NULL,
                         test.value = NULL,
                         nboots = 200, 
                         bootcluster = NULL, #2-way clusters in wild (refinement)
                         core.num) {  

    if (!is.null(seed)) {
        set.seed(seed)
    }
    p <- dim(x)[2]

    if (jackknife == 1) {
        wild <- 0
    }

    if(parallel==TRUE){
        core.num <- 1
        requireNamespace("doParallel")
		## require(iterators)
		maxcores <- detectCores()
		cores <- min(maxcores, 4)
		pcl <- future::makeClusterPSOCK(cores)
		doParallel::registerDoParallel(pcl)
	    cat("Parallel computing with", cores,"cores...\n")
        use.fun <- c("SEQ","name.fe.model","ivfastplm.core","get.dof","iv.cluster.se","fastplm.core","cluster.se",
                     "solveiv","solvecpp")
    }

    ## se <- ifelse(refinement == 1, 1, 0) 
    ## model.cl <- ifelse(refinement == 1, cluster, NULL)
    model <- fastplm.core(y = y, x = x, ind = ind, robust = robust,
                          sfe.index = sfe.index, cfe.index = cfe.index, 
                          se = 1, cl = cluster, core.num = core.num)

    #save df.use & df,cl.use for wild bootstrap(only)
    df.use <- model$df.original
    df.cl.use <- model$df.cl

    ## function to get two-sided p-values
    get.pvalue <- function(vec,to_test=0) {
        if (NaN%in%vec|NA%in%vec) {
            nan.pos <- is.nan(vec)
            na.pos <- is.na(vec)
            pos <- c(which(nan.pos),which(na.pos))
            vec.a <- vec[-pos]
            a <- sum(vec.a >= to_test)/(length(vec)-sum(nan.pos|na.pos)) * 2
            b <- sum(vec.a <= to_test)/(length(vec)-sum(nan.pos|na.pos)) * 2  
        } else {
            a <- sum(vec >= to_test)/length(vec) * 2
            b <- sum(vec <= to_test)/length(vec) * 2  
        }
        return(min(as.numeric(min(a, b)),1))
    }

    ## wild bootstrap
    if (wild == TRUE) {        
        
        raw.cluster <- level.cluster <- level.count <- NULL
        if (!is.null(cluster)) {
            if(dim(cluster)[2]==1){
                ## one-way cluster 
                raw.cluster <- as.numeric(as.factor(cluster))
                level.cluster <- unique(raw.cluster)
                level.count <- table(raw.cluster)
            }
            if(dim(cluster)[2]==2){
                if(is.null(bootcluster)){
                    level.cluster.1 <- unique(as.numeric(as.factor(cluster[,1])))
                    level.cluster.2 <- unique(as.numeric(as.factor(cluster[,2])))
                    if(length(level.cluster.1)<=length(level.cluster.2)){
                        raw.cluster <- as.numeric(as.factor(cluster[,1]))
                        bootcluster <- colnames(cluster)[1]
                    }else{
                        raw.cluster <- as.numeric(as.factor(cluster[,2]))
                        bootcluster <- colnames(cluster)[2]
                    }
                    level.cluster <- unique(raw.cluster)
                    level.count <- table(raw.cluster)
                }else{
                    raw.cluster <- as.numeric(as.factor(cluster[,bootcluster]))
                    level.cluster <- unique(raw.cluster)
                    level.count <- table(raw.cluster)
                }
                cat(paste0("Mutiple Clustered Wild Bootstrap: (Boot)Clustered at ",bootcluster," level.\n"))
            }
        }

        if (refinement == 0) {
            ## wild bootstrap
            cat("Wild Bootstrap without Percentile-t Refinement.\n")
            boot.coef <- matrix(NA, p, nboots)
            res <- model$residuals
            y.fit <- y - res 
            core.num <- ifelse(parallel == TRUE, 1, core.num)

            one.boot <- function(num = NULL) {  
                if (!is.null(cluster)) {
                    ## wild boot by cluster 
                    res.p <- rep(0,dim(cluster)[1])
                    for (j in 1:length(level.cluster)) {
                        temp.level <- sample(c(-1,1), 1)
                        res.p[which(raw.cluster==level.cluster[j])] <- temp.level
                    }
                    res.p <- as.matrix(res.p)
                } else {
                    res.p <- as.matrix(sample(c(-1, 1), dim(y)[1], replace = TRUE))
                }

                y.boot <- y.fit + res.p * res 
                boot.model <- try(fastplm.core(y = y.boot, x = x, ind = ind, 
                                               sfe.index = sfe.index, cfe.index = cfe.index, 
                                               se = 0, cl = NULL, core.num = core.num,
                                               df.use = df.use,
                                               df.cl.use = df.cl.use))

                if ('try-error' %in% class(boot.model)) {
                    return(rep(NA, p))
                } else {
                    return(c(boot.model$coefficients))
                }
            }

            if (parallel == TRUE) {
                boot.out <- foreach(j = 1:nboots, 
                                    .inorder = FALSE,
                                    .export = use.fun,
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
                    if (i%%50==0) cat(i) else cat(".")
                }
            }
            
            P_x <- as.matrix(apply(boot.coef, 1, get.pvalue))
            stderror <- as.matrix(apply(boot.coef, 1, sd, na.rm = TRUE))
            CI <- t(apply(boot.coef, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
            ## rewrite uncertainty estimates
            est.coefficients <- cbind(model$coefficients, stderror,model$coefficients/stderror, P_x, CI)
            colnames(est.coefficients) <- c("Coef", "Std. Error","Z Value", "P Value", "CI_lower", "CI_upper")
            rownames(est.coefficients) <- colnames(x)
            model$est.coefficients <- est.coefficients
            vcov <- as.matrix(cov(t(boot.coef),use="na.or.complete"))
            colnames(vcov) <- colnames(x)
            rownames(vcov) <- colnames(x)
            model$vcov <- vcov
        
        } else {
            if(is.null(pos)==TRUE){
                cat("Unrestricted Wild Bootstrap with Percentile-t Refinement.\n")
                boot.coef <- matrix(NA, p, nboots)
                res <- model$residuals 
                y.fit <- y - res
                boot.wald <- rep(NA, nboots)

                ## bootstrap refinement: wald
                core.num <- ifelse(parallel == TRUE, 1, core.num)

                one.boot <- function(num = NULL) {
                    if (!is.null(cluster)) {
                        ## wild boot by cluster 
                        res.p <- rep(0,dim(cluster)[1])
                        for (j in 1:length(level.cluster)) {
                            temp.level <- sample(c(-1,1), 1)
                            res.p[which(raw.cluster==level.cluster[j])] <- temp.level
                        }
                        res.p <- as.matrix(res.p)
                    } else {
                        res.p <- as.matrix(sample(c(-1, 1), dim(y)[1], replace = TRUE))
                    }
                

                    y.boot <- y.fit + res.p * res 
                    boot.model <- fastplm.core(y = y.boot, x = x, ind = ind, robust = robust,
                                                    sfe.index = sfe.index, cfe.index = cfe.index, 
                                                    se = 1, cl = cluster, core.num = core.num,
                                                    df.use = df.use,
                                                    df.cl.use = df.cl.use)
                    if ('try-error' %in% class(boot.model)) {
                        return(NA)
                    } else {
                        boot.use <- (boot.model$est.coefficients[, "Coef"]-model$est.coefficients[, "Coef"])/boot.model$est.coefficients[, "Std. Error"]
                        return(c(boot.use))
                    }
                }

                if (parallel == TRUE) {
                    boot.out <- foreach(j = 1:nboots, 
                                        .inorder = FALSE,
                                        .export = use.fun,
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
                        if (i%%50==0) cat(i) else cat(".")
                    }
                }

                wald.refinement <- matrix(NA,nrow=0,ncol=6)
                null.totest <- model$est.coefficients[,'t value']
                for(k in 1:length(null.totest)){
                    p.refine <- get.pvalue(boot.coef[k,],to_test=null.totest[k])
                    CI.refine <- quantile(boot.coef[k,], probs = c(0.025,0.975), na.rm = TRUE)*model$est.coefficients[k,'Std. Error'] + model$est.coefficients[k,'Coef']
                    refine.sub <- c(model$est.coefficients[k,'Coef'],
                                    model$est.coefficients[k,'Std. Error'],
                                    model$est.coefficients[k,'t value'],
                                    p.refine,
                                    CI.refine[1],
                                    CI.refine[2])
                    wald.refinement <- rbind(wald.refinement,refine.sub)
                }
                colnames(wald.refinement) <- c("Coef","Std Error","t value",
                                                "Refined P value","Refined CI_lower",
                                                "Refined CI_upper")
                rownames(wald.refinement) <- colnames(x)

                model$refinement <- list(wald.refinement = wald.refinement, boot.wald = boot.coef)
            }

            if(is.null(pos)==FALSE){

                cat("Restricted Wild Bootstrap with Percentile-t Refinement.\n")
                if(is.null(test.value)==TRUE){
                    test.value <- 0
                }
                
                if (p == 1) {
                    #under restricted conditions
                    y.restricted <- y - test.value*x
                    model.res <- fastplm.core(y = y.restricted, x = NULL, ind = ind, 
                                              sfe.index = sfe.index, cfe.index = cfe.index, 
                                              se = 0, cl = NULL, core.num = core.num)
                } else {
                    #under restricted conditions
                    y.restricted <- y - test.value*x[,pos]
                    pos.index <- which(colnames(x)==pos)
                    model.res <- fastplm.core(y = y.restricted, x = as.matrix(x[, -pos.index]), 
                                              ind = ind, sfe.index = sfe.index, cfe.index = cfe.index, 
                                              se = 0, cl = NULL, core.num = core.num)
                }

                res <- model.res$residuals 
                y.fit <- y - res #(y.fit=y.restricted-res+test.value*x=y-res)
                boot.wald <- rep(NA, nboots)

                one.boot <- function(num = NULL) {
                    if (!is.null(cluster)) {
                        ## wild boot by cluster 
                        res.p <- rep(0,dim(cluster)[1])
                        for (j in 1:length(level.cluster)) {
                            temp.level <- sample(c(-1,1), 1)
                            res.p[which(raw.cluster==level.cluster[j])] <- temp.level
                        }
                        res.p <- as.matrix(res.p)
                    } else {
                        res.p <- as.matrix(sample(c(-1, 1), dim(y)[1], replace = TRUE))
                    }

                    y.boot <- y.fit + res.p * res 
                    boot.model <- try(fastplm.core(y = y.boot, x = x, ind = ind, robust = robust,
                                                   sfe.index = sfe.index, cfe.index = cfe.index, 
                                                   se = 1, cl = cluster, core.num = core.num,
                                                   df.use = df.use,
                                                   df.cl.use = df.cl.use))

                    if ('try-error' %in% class(boot.model)) {
                        return(NA)
                    } else {
                        boot.use <- (boot.model$est.coefficients[pos, "Coef"]-test.value)/boot.model$est.coefficients[pos, "Std. Error"]
                        return(c(boot.use))
                    }
                }

                if (parallel == TRUE) {
                    boot.out <- foreach(j = 1:nboots, 
                                        .inorder = FALSE,
                                        .export = use.fun,
                                        .packages = c("fastplm")
                                        ) %dopar% {
                                            return(one.boot())
                                        }
                    boot.wald <- c(unlist(boot.out)) 
                } else {
                    for (i in 1:nboots) {
                        boot.wald[i] <- one.boot()
                        if (i%%50==0) cat(i) else cat(".")
                    }
                }

                p.refine <- get.pvalue(boot.wald,to_test=(model$est.coefficients[pos, "Coef"]-test.value)/model$est.coefficients[pos, "Std. Error"])
                CI.refine <- quantile(boot.wald, probs = c(0.025,0.975), na.rm = TRUE)*model$est.coefficients[pos,'Std. Error'] + model$est.coefficients[pos, "Coef"]

                refine.sub <- c(model$est.coefficients[pos,'Coef'],
                                model$est.coefficients[pos,'Std. Error'],
                                model$est.coefficients[pos,'t value'],
                                p.refine,
                                CI.refine[1],
                                CI.refine[2])

                names(refine.sub) <- c("Coef","Std Error","t value",
                                       "Refined P value","Refined CI_lower",
                                       "Refined CI_upper")
                model$refinement <- list(wald.refinement = refine.sub, boot.wald = boot.wald)
            }
        }
    } 
    else {    
        core.num <- ifelse(parallel == TRUE, 1, core.num)
        ## non-parametric bootstrap

        if (jackknife ==1 || is.null(cluster)) { #jackknife or non-clustered bootstrap
            jack.pos <- 1
            if (jackknife == 1) {
                nboots <- length(unique(ind[,1]))
                jack.pos <- as.numeric(as.factor(ind[,1]))
            }
            boot.coef <- matrix(NA, p, nboots)
            
            if (refinement == 0){ ## non-parametric bootstrap

                if (jackknife == 1){
                    cat("Jackknife without Percentile-t Refinement.\n")
                }else{
                    cat("Pairs Bootstrap without Percentile-t Refinement.\n")
                }

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
                                        .export = use.fun,
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
                        if (i%%50==0) cat(i) else cat(".")
                    }
                }

                if (jackknife == 0) {
                    P_x <- as.matrix(apply(boot.coef, 1, get.pvalue))
                    stderror <- as.matrix(apply(boot.coef, 1, sd, na.rm = TRUE))
                    CI <- t(apply(boot.coef, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                    ## rewrite uncertainty estimates
                    est.coefficients <- cbind(model$coefficients, stderror,model$coefficients/stderror, P_x, CI)
                    vcov <- as.matrix(cov(t(boot.coef),use="na.or.complete"))
                } else {
                    beta.j <- jackknifed(model$coefficients, boot.coef, 0.05)
                    est.coefficients <- cbind(model$coefficients, beta.j$se,model$coefficients/beta.j$se, beta.j$P, beta.j$CI.l, beta.j$CI.u)
                    vcov <- beta.j$Yvcov
                }
                colnames(est.coefficients) <- c("Coef", "Std. Error","Z Value", "P Value", "CI_lower", "CI_upper")
                rownames(est.coefficients) <- colnames(x)
                model$est.coefficients <- est.coefficients
                colnames(vcov) <- colnames(x)
                rownames(vcov) <- colnames(x)
                model$vcov <- vcov
            }else{
                if (jackknife == 1){
                    stop("Can't refine a Jackknife P-Value,\n")
                    #cat("Jackknife with Percentile-t Refinement.\n")
                }else{
                    cat("Pairs Bootstrap with Percentile-t Refinement.\n")
                }

                
                    one.boot <- function(num = NULL) {
                        boot.id <- sample(1:dim(y)[1], dim(y)[1], replace = TRUE)
                        boot.model <- fastplm.core(y = as.matrix(y[boot.id,]), 
                                                       x = as.matrix(x[boot.id,]), 
                                                       ind = as.matrix(ind[boot.id,]), 
                                                       robust = robust,
                                                       sfe.index = sfe.index, cfe.index = cfe.index, 
                                                       se = 1, cl = NULL, core.num = core.num)
                
                        if ('try-error' %in% class(boot.model)) {
                            return(rep(NA, p))
                        } else {
                            boot.use <- (boot.model$est.coefficients[, "Coef"]-model$est.coefficients[, "Coef"])/boot.model$est.coefficients[, "Std. Error"]
                            return(c(boot.use))
                        }
                    }

                    if (parallel == TRUE) {
                        boot.out <- foreach(j = 1:nboots, 
                                            .inorder = FALSE,
                                            .export = use.fun,
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
                            if (i%%50==0) cat(i) else cat(".")
                        }
                    }

                    wald.refinement <- matrix(NA,nrow=0,ncol=6)
                    null.totest <- model$est.coefficients[,'t value']
                    for(k in 1:length(null.totest)){
                        p.refine <- get.pvalue(boot.coef[k,],to_test=null.totest[k])
                        CI.refine <- quantile(boot.coef[k,], probs = c(0.025,0.975), na.rm = TRUE)*model$est.coefficients[k,'Std. Error'] + model$est.coefficients[k,'Coef']
                        refine.sub <- c(model$est.coefficients[k,'Coef'],
                                        model$est.coefficients[k,'Std. Error'],
                                        model$est.coefficients[k,'t value'],
                                        p.refine,
                                        CI.refine[1],
                                        CI.refine[2])
                        wald.refinement <- rbind(wald.refinement,refine.sub)
                    }
                    colnames(wald.refinement) <- c("Coef","Std Error","t value",
                                                   "Refined P value","Refined CI_lower",
                                                   "Refined CI_upper")
                    rownames(wald.refinement) <- colnames(x)
                    model$refinement <- list(wald.refinement = wald.refinement, boot.wald = boot.coef)
                
            }
        } else { #clustered-bootstrap
            
            if (refinement == 0) {
                cat("Pairs Bootstrap without Percentile-t Refinement.\n")
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
                                            .export = use.fun,
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
                            if (i%%50==0) cat(i) else cat(".")
                        }
                    }

                    P_x <- as.matrix(apply(boot.coef, 1, get.pvalue))
                    stderror <- as.matrix(apply(boot.coef, 1, sd, na.rm = TRUE))
                    CI <- t(apply(boot.coef, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                    ## rewrite uncertainty estimates
                    est.coefficients <- cbind(model$coefficients, stderror,model$coefficients/stderror, P_x, CI)
                    colnames(est.coefficients) <- c("Coef", "Std. Error","Z Value", "P Value", "CI_lower", "CI_upper")
                    rownames(est.coefficients) <- colnames(x)
                    model$est.coefficients <- est.coefficients
                    vcov <- as.matrix(cov(t(boot.coef),use="na.or.complete"))
                    colnames(vcov) <- colnames(x)
                    rownames(vcov) <- colnames(x)
                    model$vcov <- vcov
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
                                            .export = use.fun,
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
                            if (i%%50==0) cat(i) else cat(".")
                        }
                    }
                    
                    stderror1 <- as.matrix(apply(boot.coef1, 1, sd, na.rm = TRUE))
                    stderror2 <- as.matrix(apply(boot.coef2, 1, sd, na.rm = TRUE))
                    stderror3 <- as.matrix(apply(boot.coef3, 1, sd, na.rm = TRUE))
                    stderror <- stderror1 + stderror2 - stderror3

                    vcov1 <- as.matrix(cov(t(boot.coef1),use="na.or.complete"))
                    vcov2 <- as.matrix(cov(t(boot.coef2),use="na.or.complete"))
                    vcov3 <- as.matrix(cov(t(boot.coef3),use="na.or.complete"))
                    vcov <- vcov1 + vcov2 - vcov3
                    colnames(vcov) <- colnames(x)
                    rownames(vcov) <- colnames(x)
                    model$vcov <- vcov


                    CI <- cbind(c(model$coefficients) - qnorm(0.975)*c(stderror), c(model$coefficients) + qnorm(0.975)*c(stderror))
                    Zx <- c(model$coefficients)/c(stderror)
                    P_z <- 2 * min(1 - pnorm(Zx), pnorm(Zx))
                    ## rewrite uncertainty estimates
                    est.coefficients <- cbind(model$coefficients, stderror, Zx, P_z, CI)
                    colnames(est.coefficients) <- c("Coef", "Std. Error", "Z value", "Pr(>|Z|)", "CI_lower", "CI_upper")
                    rownames(est.coefficients) <- colnames(x)
                    model$est.coefficients <- est.coefficients
                }

            } else {
                cat("Pairs Bootstrap with Percentile-t Refinement.\n")
                ## bootstrap refinement 
                boot.coef <- matrix(NA, p, nboots)

                if(dim(cluster)[2]==1){
                    ## one-way cluster 
                    raw.cluster <- as.numeric(as.factor(cluster))
                    level.cluster <- unique(raw.cluster)
                    level.count <- table(raw.cluster)
                }
                if(dim(cluster)[2]==2){
                    stop("For pairs bootstrap with refinement, please only specify one cluster variable.\n")
                    ##cl_inter <- as.numeric(as.factor(paste(cluster[,1], "-:-", cluster[,2], sep = "")))
                    #if cl_inter is unique; then stderror3 is robust standard error
                    ##cl_inter.count <- table(cl_inter)
                    ##if(max(cl_inter.count)==1){
                       ##cat("Warning: Pairs Bootstrap with Percentile-t Refinement may be significantly Slower than the Wild Bootstrap.\n") 
                    ##}
                    if(is.null(bootcluster)){
                        level.cluster.1 <- unique(as.numeric(as.factor(cluster[,1])))
                        level.cluster.2 <- unique(as.numeric(as.factor(cluster[,2])))
                        if(length(level.cluster.1)<=length(level.cluster.2)){
                            raw.cluster <- as.numeric(as.factor(cluster[,1]))
                            bootcluster <- colnames(cluster)[1]
                        }else{
                            raw.cluster <- as.numeric(as.factor(cluster[,2]))
                            bootcluster <- colnames(cluster)[2]
                        }
                        level.cluster <- unique(raw.cluster)
                        level.count <- table(raw.cluster)
                    }else{
                        raw.cluster <- as.numeric(as.factor(cluster[,bootcluster]))
                        level.cluster <- unique(raw.cluster)
                        level.count <- table(raw.cluster)
                    }
                    cat(paste0("Mutiple Clustered Pairs Bootstrap with Refinement: (Boot)Clustered at ",bootcluster," level.\n"))
                }

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
                                                    robust = robust,
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
                                        .export = use.fun,
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
                        if (i%%50==0) cat(i) else cat(".")
                    }
                }
            
                wald.refinement <- matrix(NA,nrow=0,ncol=6)
                null.totest <- model$est.coefficients[,'t value']
                for(k in 1:length(null.totest)){
                    p.refine <- get.pvalue(boot.coef[k,],to_test=null.totest[k])
                    CI.refine <- quantile(boot.coef[k,], probs = c(0.025,0.975), na.rm = TRUE)*model$est.coefficients[k,'Std. Error'] + model$est.coefficients[k,'Coef']
                    refine.sub <- c(model$est.coefficients[k,'Coef'],
                                    model$est.coefficients[k,'Std. Error'],
                                    model$est.coefficients[k,'t value'],
                                    p.refine,
                                    CI.refine[1],
                                    CI.refine[2])
                    wald.refinement <- rbind(wald.refinement,refine.sub)
                }
                colnames(wald.refinement) <- c("Coef","Std Error","t value",
                                               "Refined P value","Refined CI_lower",
                                               "Refined CI_upper")
                rownames(wald.refinement) <- colnames(x)

                model$refinement <- list(wald.refinement = wald.refinement, boot.wald = boot.coef)

            }
        }
        

    }



    if(parallel==TRUE){
        suppressWarnings(stopCluster(pcl))
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

    Yvcov <- cov(t(Y),use="na.or.complete")
    Yvcov <- Yvcov/vn 

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

    out <- list(se = Ysd, CI.l = CI.l, CI.u = CI.u, P = P, Yvcov=Yvcov)

    return(out)
    
}