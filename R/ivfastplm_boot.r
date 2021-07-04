#iv fastplm boot
ivfastplm.boot <- function(seed, 
                           y = NULL, 
                           x = NULL, 
                           z = NULL,
                           ind = NULL, 
                           sfe.index = NULL, 
                           cfe.index = NULL, 
                           cluster = NULL,
                           robust = FALSE,
                           parallel = FALSE,
                           wild = FALSE, 
                           jackknife = FALSE,
                           refinement = FALSE, 
                           pos = NULL,
                           test.value = NULL,
                           nboots = 200, 
                           bootcluster = NULL, #2-way clusters in wild (refinement)
                           core.num = 1){
    p <- dim(x)[2]
    if(is.null(z)){
        stop("No exogenous variables.\n")
    }

    if(is.null(colnames(x))){
        colnames(x) <- paste0("x.",c(1:dim(x)[2]))
    }

    if(is.null(colnames(z))){
        colnames(z) <- paste0("z.",c(1:dim(z)[2]))
    }

    if(is.null(colnames(y))){
        colnames(y) <- 'y'
    }

    if (jackknife == 1) {
        wild <- 0
    }

    # get 
    model.iv <- ivfastplm.core(y = y, x = x, z = z, 
                               ind = ind, se = 1, 
                               robust = robust,
                               sfe.index = sfe.index, 
                               cfe.index = cfe.index, 
                               cl = cluster, 
                               core.num = core.num, 
                               need.fe = TRUE,
                               iv.test = FALSE)

    dy <- model.iv$demeaned$dy
    dx <- model.iv$demeaned$dx
    dz <- model.iv$demeaned$dz

    y.fit <- fitted.y <- model.iv$fitted.values
    res <- residual.y <- model.iv$residuals

    iv.name <- model.iv$names
    en.var.x <- iv.name$en.var.x
    ex.var.x <- iv.name$ex.var.x
    include.iv <- iv.name$include.iv
    exclude.iv <- iv.name$exclude.iv
    
    iv.matrix <- model.iv$matrix
    inv.z <- iv.matrix$inv.z # inv(z'z)
    inv.xPzx <- iv.matrix$inv.xPzx

    #save df.original & df.cl for wild bootstrap(only)
    df.use <- df.original <- model.iv$df.original
    df.cl.use <- df.cl <- model.iv$df.cl

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

    if(parallel==TRUE){
        core.num <- 1
        requireNamespace("doParallel")
		## require(iterators)
		maxcores <- detectCores()
		cores <- min(maxcores, 4)
		pcl<-makeCluster(cores)  
		doParallel::registerDoParallel(pcl)
	    cat("Parallel computing with", cores,"cores...\n")
        use.fun <- c("SEQ","name.fe.model","ivfastplm.core","get.dof","iv.cluster.se","fastplm.core")
    }

    if(wild==TRUE){
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

        if(length(en.var.x)>0){
           en.var.x.matrix <- matrix(dx[,en.var.x],ncol=length(en.var.x))
            sub.lm <- simpleols(Y=en.var.x.matrix,X=dz)
            en.var.x.coef <- sub.lm[['coef']]
            colnames(en.var.x.coef) <- en.var.x
            en.var.x.hat <- sub.lm[['Y_hat']]
            colnames(en.var.x.hat) <- en.var.x
            res.en.x <- en.var.x.res <- sub.lm[["u_hat"]]
            en.var.x.fitted <- matrix(x[,en.var.x],ncol=length(en.var.x)) - en.var.x.res
            colnames(en.var.x.fitted) <- colnames(en.var.x.res) <- colnames(en.var.x.hat) <- en.var.x 
        }
        ####

        if (refinement == 0) {
                ## wild bootstrap
                cat("Wild Bootstrap without Percentile-t Refinement.\n")
                boot.coef <- matrix(NA, p, nboots)
                
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
                    
                    if(length(en.var.x)>0){
                       x.boot.en <- matrix(en.var.x.fitted + sweep(en.var.x.res, MARGIN=1, res.p, `*`),ncol=length(en.var.x)) 
                       colnames(x.boot.en) <- en.var.x
                    }else{
                        x.boot.en <- matrix(NA,ncol=0,nrow=dim(dy)[1])
                    }
                    
                    if(length(ex.var.x)>0){
                        x.boot.ex <- matrix(x[,ex.var.x],ncol=length(ex.var.x))
                        colnames(x.boot.ex) <- ex.var.x
                    }else{
                        x.boot.ex <- matrix(NA,ncol=0,nrow=dim(dy)[1])
                    }

                    x.boot <- cbind(x.boot.en,x.boot.ex)
                    x.boot <- x.boot[,colnames(x)]
                    #only need coefficients, don't need se
                    boot.model.iv <- try(ivfastplm.core(y = y.boot, x = x.boot,z=z, ind = ind, 
                                                   sfe.index = sfe.index, cfe.index = cfe.index, 
                                                   se = 0, cl = NULL, core.num = core.num,
                                                   df.use = df.use, df.cl.use = df.cl.use, need.fe = FALSE, iv.test = FALSE
                                                   ))

                    if ('try-error' %in% class(boot.model.iv)) {
                        return(rep(NA, p))
                    } else {
                        return(c(boot.model.iv$coefficients))
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
                est.coefficients <- cbind(model.iv$coefficients, stderror,model.iv$coefficients/stderror, P_x, CI)
                colnames(est.coefficients) <- c("Coef", "Std. Error","Z Value", "P Value", "CI_lower", "CI_upper")
                rownames(est.coefficients) <- colnames(dx)
                model.iv$est.coefficients <- est.coefficients

                vcov <- as.matrix(cov(t(boot.coef),use="na.or.complete"))
                colnames(vcov) <- colnames(dx)
                rownames(vcov) <- colnames(dx)
                model.iv$vcov <- vcov
        }

        if(refinement==1){
            if(is.null(pos)==TRUE){
                cat("Unrestricted Wild Bootstrap with Percentile-t Refinement.\n")
                boot.coef <- matrix(NA, p, nboots)
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
                   
                    if(length(en.var.x)>0){
                       x.boot.en <- matrix(en.var.x.fitted + sweep(en.var.x.res, MARGIN=1, res.p, `*`),ncol=length(en.var.x)) 
                       colnames(x.boot.en) <- en.var.x
                    }else{
                        x.boot.en <- matrix(NA,ncol=0,nrow=dim(dy)[1])
                    }
                    
                    if(length(ex.var.x)>0){
                        x.boot.ex <- matrix(x[,ex.var.x],ncol=length(ex.var.x))
                        colnames(x.boot.ex) <- ex.var.x
                    }else{
                        x.boot.ex <- matrix(NA,ncol=0,nrow=dim(dy)[1])
                    } 
                    x.boot <- cbind(x.boot.en,x.boot.ex)
                    x.boot <- x.boot[,colnames(x)]

                    boot.model <- try(ivfastplm.core(y = y.boot, x = x.boot, z = z, ind = ind, robust = robust,
                                                sfe.index = sfe.index, cfe.index = cfe.index, 
                                                se = 1, cl = cluster, core.num = core.num,
                                                df.use = df.use,
                                                df.cl.use = df.cl.use))
                    if ('try-error' %in% class(boot.model)) {
                        return(NA)
                    } else {
                        boot.use <- (boot.model$est.coefficients[, "Coef"]-model.iv$est.coefficients[, "Coef"])/boot.model$est.coefficients[, "Std. Error"]
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
                null.totest <- model.iv$est.coefficients[,'t value']
                for(k in 1:length(null.totest)){
                    p.refine <- get.pvalue(boot.coef[k,],to_test=null.totest[k])
                    CI.refine <- quantile(boot.coef[k,], probs = c(0.025,0.975), na.rm = TRUE)*model.iv$est.coefficients[k,'Std. Error'] + model.iv$est.coefficients[k,'Coef']
                    refine.sub <- c(model.iv$est.coefficients[k,'Coef'],
                                    model.iv$est.coefficients[k,'Std. Error'],
                                    model.iv$est.coefficients[k,'t value'],
                                    p.refine,
                                    CI.refine[1],
                                    CI.refine[2])
                    wald.refinement <- rbind(wald.refinement,refine.sub)
                }
                colnames(wald.refinement) <- c("Coef","Std Error","t value",
                                            "Refined P value","Refined CI_lower",
                                            "Refined CI_upper")
                rownames(wald.refinement) <- colnames(dx)
                model.iv$refinement <- list(wald.refinement = wald.refinement, boot.wald = boot.coef)
            }

            if(is.null(pos)==FALSE){

                if(is.null(test.value)==TRUE){
                    test.value <- 0
                }

                if(pos %in% en.var.x){
                    cat("Restricted Wild Efficient Residual Bootstrap with Percentile-t Refinement.\n")
                    if(length(en.var.x)==1){
                        dy.restricted <- dy - test.value*dx[,en.var.x]
                        dx1 <- matrix(dx[,ex.var.x],ncol=length(ex.var.x))
                        inv.dx1 <- solve(t(dx1)%*%dx1)
                        reg1.coef <- inv.dx1%*%t(dx1)%*%dy.restricted
                        u.hat1 <- dy.restricted - dx1%*%reg1.coef
                       
                        dy2 <- dx[,en.var.x]
                        dx2 <- cbind(dz,u.hat1)
                        inv.dx2 <- solve(t(dx2)%*%dx2)
                        reg2.coef <- inv.dx2%*%t(dx2)%*%dy2
                        u.hat2 <- dy2 - dx2%*%reg2.coef + reg2.coef[dim(dz)[2]+1,1]*u.hat1

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

                            y.boot <- y-u.hat1+res.p * u.hat1
                            en.x.boot <- x[,en.var.x] - u.hat2 + res.p * u.hat2

                            if(length(ex.var.x)>0){
                                x.boot.ex <- matrix(x[,ex.var.x],ncol=length(ex.var.x))
                                colnames(x.boot.ex) <- ex.var.x
                            }else{
                                x.boot.ex <- matrix(NA,ncol=0,nrow=dim(dy)[1])
                            } 

                            x.boot <- cbind(en.x.boot,x.boot.ex)
                            colnames(x.boot) <- c(en.var.x,ex.var.x) 
                            x.boot <- x.boot[,colnames(x)] 

                            boot.model <- try(ivfastplm.core(y = y.boot, x = x.boot, z = z, ind = ind, robust = robust,
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
                    }

                    if(length(en.var.x)>1){
                        dy.restricted <- dy - test.value*dx[,pos]
                        en.dx <- dx[,en.var.x]
                        pos.index <- which(colnames(en.dx)==pos)
                        en.dx.other <- matrix(en.dx[,-pos.index],ncol=dim(en.dx)[2]-1)
                        
                        #a simple 2sls
                        dx1 <- cbind(en.dx.other,dx[,ex.var.x])
                        inv.dx1 <- solve(t(dx1)%*%dx1)
                        dz1 <- dz
                        inv.dz1 <- inv.z
                        Pdz1 <- dz1%*%inv.dz1%*%t(dz1)
                        inv.dxPzdx1 <- solve(t(dx1)%*%Pdz1%*%dx1)
                        reg1.coef <- matrix(inv.dxPzdx1%*%t(dx1)%*%Pdz1%*%dy.restricted,ncol=1)
                        u.hat1 <- matrix(dy.restricted-dx1%*%reg1.coef,ncol=1)

                        dy2 <- dx[,en.var.x]
                        dx2 <- cbind(dz,u.hat1)
                        inv.dx2 <- solve(t(dx2)%*%dx2)
                        reg2.coef <- inv.dx2%*%t(dx2)%*%dy2

                        reg2.coef.slice <- matrix(reg2.coef[(dim(dz)[2]+1),],ncol=dim(reg2.coef)[2])
                        u.hat2 <- dy2 - dx2%*%reg2.coef + kronecker(reg2.coef.slice,u.hat1)

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

                            y.boot <- y-u.hat1+res.p * u.hat1
                            en.x.boot <- matrix(x[,en.var.x] - u.hat2 + sweep(u.hat2, MARGIN=1, res.p, `*`),ncol=length(en.var.x))
                            if(length(ex.var.x)>0){
                                x.boot.ex <- matrix(x[,ex.var.x],ncol=length(ex.var.x))
                                colnames(x.boot.ex) <- ex.var.x
                            }else{
                                x.boot.ex <- matrix(NA,ncol=0,nrow=dim(dy)[1])
                            } 

                            x.boot <- cbind(en.x.boot,x.boot.ex)
                            colnames(x.boot) <- c(en.var.x,ex.var.x) 
                            x.boot <- x.boot[,colnames(x)] 

                            boot.model <- try(ivfastplm.core(y = y.boot, x = x.boot, z = z, ind = ind, robust = robust,
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
                    }

                    boot.wald <- rep(NA, nboots)
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

                    p.refine <- get.pvalue(boot.wald,to_test=(model.iv$est.coefficients[pos, "Coef"]-test.value)/model.iv$est.coefficients[pos, "Std. Error"])
                    CI.refine <- quantile(boot.wald, probs = c(0.025,0.975), na.rm = TRUE)*model.iv$est.coefficients[pos,'Std. Error'] + model.iv$est.coefficients[pos, "Coef"]

                    refine.sub <- c(model.iv$est.coefficients[pos,'Coef'],
                                    model.iv$est.coefficients[pos,'Std. Error'],
                                    (model.iv$est.coefficients[pos, "Coef"]-test.value)/model.iv$est.coefficients[pos, "Std. Error"],
                                    p.refine,
                                    #paste0(pos,"=",test.value),
                                    CI.refine[1],
                                    CI.refine[2])
                                    
                    refine.sub <- matrix(refine.sub,nrow=1)
                    #refine.sub <- as.data.frame(refine.sub)
                    rownames(refine.sub) <- pos
                    colnames(refine.sub) <- c("Coef","Std Error","t value(H0)",
                                              "Refined P value","Refined CI_lower","Refined CI_upper")
                    model.iv$refinement <- list(wald.refinement = refine.sub, boot.wald = boot.wald)
                }

                if(pos %in% ex.var.x){
                    # wild restricted residual bootstrap
                    cat("Wild Restricted Bootstrap with Percentile-t Refinement.\n")
                    boot.wald <- rep(NA, nboots)

                    dy.restricted <- dy - test.value*dx[,pos]
                    pos.index <- which(colnames(dx)==pos)
                    dx.other <- matrix(dx[,-pos.index],ncol=dim(dx)[2]-1)

                    #a simple 2sls
                    dx1 <- dx.other
                    inv.dx1 <- solve(t(dx1)%*%dx1)
                    dz1 <- dz
                    inv.dz1 <- inv.z
                    Pdz1 <- dz1%*%inv.dz1%*%t(dz1)
                    inv.dxPzdx1 <- solve(t(dx1)%*%Pdz1%*%dx1)
                    reg1.coef <- matrix(inv.dxPzdx1%*%t(dx1)%*%Pdz1%*%dy.restricted,ncol=1)
                    u.hat1 <- matrix(dy.restricted-dx1%*%reg1.coef,ncol=1)

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
                        y.boot <- y - u.hat1 + res.p * u.hat1
                        x.boot.en <- matrix(en.var.x.fitted + sweep(en.var.x.res, MARGIN=1, res.p, `*`),ncol=length(en.var.x))
                        colnames(x.boot.en) <- en.var.x

                        if(length(ex.var.x)>0){
                            x.boot.ex <- matrix(x[,ex.var.x],ncol=length(ex.var.x))
                            colnames(x.boot.ex) <- ex.var.x
                        }else{
                            x.boot.ex <- matrix(NA,ncol=0,nrow=dim(dy)[1])
                        } 

                        x.boot <- cbind(x.boot.en,x.boot.ex)
                        colnames(x.boot) <- c(en.var.x,ex.var.x) 
                        x.boot <- x.boot[,colnames(x)] 

                        boot.model <- try(ivfastplm.core(y = y.boot, x = x.boot, z = z, ind = ind, robust = robust,
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
                        for (j in 1:nboots) { 
                            boot.wald[j] <- c(unlist(boot.out[[j]])) 
                        }
                    } else {
                        for (i in 1:nboots) {
                            boot.wald[i] <- one.boot()
                            if (i%%50==0) cat(i) else cat(".")
                        }
                    }
                
                    p.refine <- get.pvalue(boot.wald,to_test=(model.iv$est.coefficients[pos, "Coef"]-test.value)/model.iv$est.coefficients[pos, "Std. Error"])
                    CI.refine <- quantile(boot.wald, probs = c(0.025,0.975), na.rm = TRUE)*model.iv$est.coefficients[pos,'Std. Error'] + model.iv$est.coefficients[pos, "Coef"]

                   refine.sub <- c(model.iv$est.coefficients[pos,'Coef'],
                                    model.iv$est.coefficients[pos,'Std. Error'],
                                    (model.iv$est.coefficients[pos, "Coef"]-test.value)/model.iv$est.coefficients[pos, "Std. Error"],
                                    p.refine,
                                    #paste0(pos,"=",test.value),
                                    CI.refine[1],
                                    CI.refine[2])
                                    
                    refine.sub <- matrix(refine.sub,nrow=1)
                    #refine.sub <- as.data.frame(refine.sub)
                    rownames(refine.sub) <- pos
                    colnames(refine.sub) <- c("Coef","Std Error","t value(H0)",
                                              "Refined P value","Refined CI_lower","Refined CI_upper")
                    model.iv$refinement <- list(wald.refinement = refine.sub, boot.wald = boot.wald)
                }
            }
        }
    }
    else{
        
        if (jackknife ==1 || is.null(cluster)) {    
                jack.pos <- 1
                if (jackknife == 1) {
                    nboots <- length(unique(ind[,1]))
                    jack.pos <- as.numeric(as.factor(ind[,1]))
                }
                boot.coef <- matrix(NA, p, nboots)
        
        if(refinement==0){
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
                
                boot.model <- try(ivfastplm.core(y = as.matrix(y[boot.id,]), 
                                                 x = as.matrix(x[boot.id,]),
                                                 z = as.matrix(z[boot.id,]), 
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
                est.coefficients <- cbind(model.iv$coefficients, stderror,model.iv$coefficients/stderror, P_x, CI)
                vcov <- as.matrix(cov(t(boot.coef),use="na.or.complete"))
            } else {
                beta.j <- jackknifed(model.iv$coefficients, boot.coef, 0.05)
                est.coefficients <- cbind(model.iv$coefficients, beta.j$se,model.iv$coefficients/beta.j$se, beta.j$P, beta.j$CI.l, beta.j$CI.u)
                vcov <- beta.j$Yvcov
            }
            colnames(est.coefficients) <- c("Coef", "Std. Error","Z Value", "P Value", "CI_lower", "CI_upper")
            rownames(est.coefficients) <- colnames(dx)
            model.iv$est.coefficients <- est.coefficients

            colnames(vcov) <- colnames(dx)
            rownames(vcov) <- colnames(dx)
            model.iv$vcov <- vcov
        }

        if(refinement==1){
            if (jackknife == 1){
                stop("Can't refine a Jackknife P-Value,\n")
                #cat("Jackknife with Percentile-t Refinement.\n")
            }else{
                cat("Pairs Bootstrap with Percentile-t Refinement.\n")
            }

            one.boot <- function(num = NULL) {
                boot.id <- sample(1:dim(y)[1], dim(y)[1], replace = TRUE)
                boot.model <- try(ivfastplm.core(y = as.matrix(y[boot.id,]), 
                                               x = as.matrix(x[boot.id,]), 
                                               z = as.matrix(z[boot.id,]),
                                               ind = as.matrix(ind[boot.id,]), 
                                               robust = robust,
                                               sfe.index = sfe.index, cfe.index = cfe.index, 
                                               se = 1, cl = NULL, core.num = core.num))
                
                if ('try-error' %in% class(boot.model)) {
                    return(rep(NA, p))
                } else {
                    boot.use <- (boot.model$est.coefficients[, "Coef"]-model.iv$est.coefficients[, "Coef"])/boot.model$est.coefficients[, "Std. Error"]
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
            null.totest <- model.iv$est.coefficients[,'t value']
            for(k in 1:length(null.totest)){
                p.refine <- get.pvalue(boot.coef[k,],to_test=null.totest[k])
                CI.refine <- quantile(boot.coef[k,], probs = c(0.025,0.975), na.rm = TRUE)*model.iv$est.coefficients[k,'Std. Error'] + model.iv$est.coefficients[k,'Coef']
                refine.sub <- c(model.iv$est.coefficients[k,'Coef'],
                                model.iv$est.coefficients[k,'Std. Error'],
                                model.iv$est.coefficients[k,'t value'],
                                p.refine,
                                CI.refine[1],
                                CI.refine[2])
                wald.refinement <- rbind(wald.refinement,refine.sub)
            }
            colnames(wald.refinement) <- c("Coef","Std Error","t value",
                                               "Refined P value","Refined CI_lower",
                                               "Refined CI_upper")
            rownames(wald.refinement) <- colnames(dx)
            model.iv$refinement <- list(wald.refinement = wald.refinement, boot.wald = boot.coef)
        }
        }
        else{ #clustered
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
                        boot.model <- try(ivfastplm.core(y = as.matrix(y[boot.id,]), 
                                                         x = as.matrix(x[boot.id,]), 
                                                         z = as.matrix(z[boot.id,]),
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
                    est.coefficients <- cbind(model.iv$coefficients, stderror,model.iv$coefficients/stderror, P_x, CI)
                    colnames(est.coefficients) <- c("Coef", "Std. Error","Z Value", "P Value", "CI_lower", "CI_upper")
                    rownames(est.coefficients) <- colnames(dx)
                    model.iv$est.coefficients <- est.coefficients

                    vcov <- as.matrix(cov(t(boot.coef),use="na.or.complete"))
                    colnames(vcov) <- colnames(dx)
                    rownames(vcov) <- colnames(dx)
                    model.iv$vcov <- vcov
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
                        boot.model1 <- try(ivfastplm.core(y = as.matrix(y[boot.id1,]), 
                                                         x = as.matrix(x[boot.id1,]),
                                                         z = as.matrix(z[boot.id1,]), 
                                                         ind = as.matrix(ind[boot.id1,]), 
                                                         sfe.index = sfe.index, cfe.index = cfe.index, 
                                                         se = 0, cl = NULL, core.num = core.num))
                        if ('try-error' %in% class(boot.model1)) {
                            oneboot.coef1 <- rep(NA, p)
                        } else {
                            oneboot.coef1 <- c(boot.model1$coefficients)
                        }

                        ## level 2
                        boot.model2 <- try(ivfastplm.core (y = as.matrix(y[boot.id2,]), 
                                                         x = as.matrix(x[boot.id2,]), 
                                                         z = as.matrix(z[boot.id2,]),
                                                         ind = as.matrix(ind[boot.id2,]), 
                                                         sfe.index = sfe.index, cfe.index = cfe.index, 
                                                         se = 0, cl = NULL, core.num = core.num))
                        if ('try-error' %in% class(boot.model2)) {
                            oneboot.coef2 <- rep(NA, p)
                        } else {
                            oneboot.coef2 <- c(boot.model2$coefficients)
                        }

                        ## intersection
                        boot.model3 <- try(ivfastplm.core (y = as.matrix(y[boot.id3,]), 
                                                         x = as.matrix(x[boot.id3,]), 
                                                         z = as.matrix(z[boot.id3,]),
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
                    colnames(vcov) <- colnames(dx)
                    rownames(vcov) <- colnames(dx)
                    model.iv$vcov <- vcov

                    CI <- cbind(c(model.iv$coefficients) - qnorm(0.975)*c(stderror), c(model.iv$coefficients) + qnorm(0.975)*c(stderror))
                    Zx <- c(model.iv$coefficients)/c(stderror)
                    P_z <- 2 * min(1 - pnorm(Zx), pnorm(Zx))
                    ## rewrite uncertainty estimates
                    est.coefficients <- cbind(model.iv$coefficients, stderror, Zx, P_z, CI)
                    colnames(est.coefficients) <- c("Coef", "Std. Error", "Z value", "Pr(>|Z|)", "CI_lower", "CI_upper")
                    rownames(est.coefficients) <- colnames(dx)
                    model.iv$est.coefficients <- est.coefficients
                }
            }

            if(refinement==1){
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

                    boot.model <- try(ivfastplm.core (y = as.matrix(y[boot.id,]), 
                                                    x = as.matrix(x[boot.id,]), 
                                                    z = as.matrix(z[boot.id,]),
                                                    ind = as.matrix(ind[boot.id,]), robust = robust,
                                                    sfe.index = sfe.index, cfe.index = cfe.index, 
                                                    se = 1, cl = boot.cluster, core.num = core.num))
                    if ('try-error' %in% class(boot.model)) {
                        return(rep(NA, p))
                    } else {
                        return(c(boot.model$coefficients - model.iv$coefficients)/c(boot.model$est.coefficients[, "Std. Error"]))
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
                null.totest <- model.iv$est.coefficients[,'t value']
                for(k in 1:length(null.totest)){
                    p.refine <- get.pvalue(boot.coef[k,],to_test=null.totest[k])
                    CI.refine <- quantile(boot.coef[k,], probs = c(0.025,0.975), na.rm = TRUE)*model.iv$est.coefficients[k,'Std. Error'] + model.iv$est.coefficients[k,'Coef']
                    refine.sub <- c(model.iv$est.coefficients[k,'Coef'],
                                    model.iv$est.coefficients[k,'Std. Error'],
                                    model.iv$est.coefficients[k,'t value'],
                                    p.refine,
                                    CI.refine[1],
                                    CI.refine[2])
                    wald.refinement <- rbind(wald.refinement,refine.sub)
                }
                colnames(wald.refinement) <- c("Coef","Std Error","t value",
                                               "Refined P value","Refined CI_lower",
                                               "Refined CI_upper")
                rownames(wald.refinement) <- colnames(dx)
                model.iv$refinement <- list(wald.refinement = wald.refinement, boot.wald = boot.coef)
            } 
        }
    }
    if(parallel==TRUE){
        suppressWarnings(stopCluster(pcl))
    }

    




    return(model.iv)
   
}