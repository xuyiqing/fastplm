ivfastplm.core <- function(y = NULL, ## outcome vector
                         x = NULL, ## regressor matrix
                         z = NULL, ## exogenous variables matrix
                         ind = NULL, ## indicator matrix
                         se = 0,
                         robust = FALSE,
                         cl = NULL,
                         sfe.index = NULL, ## index for simple fe , a vector of integer
                         cfe.index = NULL, ## index for complex fe, a list of double, e.g. c(effect, influence)
                         core.num = 1,
                         df.use = NULL,
                         df.cl.use = NULL,
                         need.fe = FALSE,
                         iv.test = FALSE,
                         first = TRUE,
                         orthog = NULL,
                         endog = NULL) {

    if(is.null(z)){
        stop("No exogenous variables.\n")
    }
    
    if(is.null(x)){
        stop("No regressors.\n")
    }

    data.reduce <- cbind(y, z)
    data <- cbind(y, x)

    ## create indicators
    inds <- create.indicators(ind) ## an indicator class object

    ## ------ fe ------- ##
    fe <- cfes <- sub.cfe.index <- num.cfe.index <- NULL

    ## create complex fe
    gcfe <- function(i) {
        sub.cfe.index <- cfe.index[[i]]
        inf.weight <- as.numeric(inds$levels[[sub.cfe.index[2]]])
        cfe.sub <- create.complex.effect(inds, sub.cfe.index[1], 
                                         sub.cfe.index[2], t(as.matrix(inf.weight)))
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
    

    ## demean all variables and then recover the fixed effects in ivregression using the coefficients
    model <- SolveFixedEffects(data, fe$ptr, core.num)
    model <- name.fe.model(model, inds, fe)
    model$inds <- inds
    model$fe <- fe

    model.reduce <- SolveFixedEffects(data.reduce, fe$ptr, core.num)
    model.reduce <- name.fe.model(model.reduce, inds, fe)
    model.reduce$inds <- inds
    model.reduce$fe <- fe

    dy <- model$demeaned$y
    dz <- model.reduce$demeaned$x
    dx <- model$demeaned$x
    colnames(dy) <- colnames(y)
    colnames(dz) <- colnames(z)
    colnames(dx) <- colnames(x)

    #inv.x <- solve(t(dx)%*%dx)
    #inv.z <- solve(t(dz)%*%dz)
    #Pz <- dz%*%inv.z%*%t(dz)
    #inv.xPzx <- solve(t(dx)%*%Pz%*%dx)
    #inv.z.zx <- inv.z%*%t(dz)%*%dx
    #iv.coef <- matrix(inv.xPzx%*%t(dx)%*%Pz%*%dy,ncol=1)
    #u.hat <- matrix(dy-dx%*%iv.coef,ncol=1)

    iv.solve <- try(solveiv(Y=dy,X=dx,Z=dz),silent = T)
    if ('try-error' %in% class(iv.solve)) {
        stop("Exact collinearity in the regression.\n")
    }
    inv.x <- iv.solve[["invxx"]]
    inv.z <- iv.solve[["invzz"]]
    inv.xPzx <- iv.solve[["invXPzX"]]
    iv.coef <- iv.solve[["coef"]]
    inv.z.zx <- iv.solve[["invzzx"]]
    rownames(iv.coef) <- colnames(x)
    u.hat <- iv.solve[["residual"]]
    ZY <- iv.solve[["ZY"]]; #t(Z)*y
    XZ <- iv.solve[["XZ"]]; #t(X)*Z

    if(sum(is.na(inv.x))>0){
        stop("Exact collinearity in the regression.")
    }
    if(sum(is.na(inv.z))>0){
        stop("Exact collinearity in the regression.")
    }

    colnames(inv.xPzx) <- rownames(inv.xPzx) <- colnames(inv.x) <- rownames(inv.x) <- colnames(dx)
    colnames(inv.z) <- rownames(inv.z) <- colnames(dz)
    rownames(XZ) <- colnames(dx)
    colnames(XZ) <- colnames(dz)
    
    colnames(ZY) <- colnames(dy)
    rownames(ZY) <- colnames(dz)
    
    model.iv <- list(coefficients=iv.coef,res=u.hat)
    model.iv$fe <- fe
    if(se == 1){
        if(is.null(df.use)){
            dof.output <- get.dof(model=model,
                                  x=x,
                                  cl=NULL,
                                  ind=ind,
                                  sfe.index=sfe.index,
                                  cfe.index=cfe.index)
            gtot <- dof.output$gtot    
            model.iv$dof <- dof.output$dof
            ## degree of freedom
            df <- dim(dy)[1] - gtot
        }else{
            df <- df.use
        }
        model.iv$df.original <- df
        model$df.original <- df
        
        if(is.null(cl)){

            sigma2 <- t(u.hat)%*%u.hat/(df)
            RMSE <- sqrt(sigma2)
            R2 <- 1 - sum(t(u.hat)%*%u.hat)/sum((data[, 1] - mean(data[, 1]))^2)
            Adj_R2 <- 1 - (1 - R2) * (dim(x)[1] - 1) / df
            projR2 <- 1 - sum(t(u.hat)%*%u.hat)/sum((c(model$demeaned$y) - mean(model$demeaned$y))^2)
            projAdj_R2 <- 1 - (1 - projR2) * (dim(x)[1] - 1) / df

            if (robust == FALSE) {
                vcov <- as.matrix(c(t(u.hat)%*%u.hat/(df))*inv.xPzx)
                stderror <- as.matrix(sqrt(diag(vcov)))
            } else{
                meat <- dz * matrix(rep(c(u.hat), dim(dz)[2]), dim(dz)[1], dim(dz)[2])
                meat <- t(meat) %*% meat
                vcov <- inv.xPzx%*%(t(inv.z.zx)%*%meat%*%inv.z.zx)%*%inv.xPzx
                vcov <- dim(dx)[1] / df*vcov
                stderror <- as.matrix(sqrt(diag(vcov)))
            }

            Tx <- c(model.iv$coefficients)/c(stderror)

            P_t <- NULL
            for (i in 1:length(Tx)) {
                subt <- pt(Tx[i], df)
                P_t <- c(P_t, 2 * min(1 - subt, subt))
            } 
            CI <- cbind(c(model.iv$coefficients) - qt(0.975,df=df)*c(stderror), c(model.iv$coefficients) + qt(0.975,df=df)*c(stderror))
            model.iv$df.residual <- df

            ## F test
            if (robust == FALSE) {
                ## full model
                ##fF <- R2/(1 - R2) * df/(dim(x)[1] - df - 1)
                ## projection model
                projF <- t(model.iv$coefficients) %*% solve(vcov) %*% model.iv$coefficients/dim(x)[2]

                ##F_statistic <- cbind(fF, dim(x)[1] - df - 1, df, 1 - pf(fF, dim(x)[1] - df - 1, df))
                proj_F_statistic <- cbind(projF, dim(x)[2], df, 1 - pf(projF, dim(x)[2], df))

                ##colnames(F_statistic) <- c("F_statistic", "df1", "df2", "P value")
                colnames(proj_F_statistic) <- c("F_statistic", "df1", "df2", "P value")

                ##model.iv$F_statistic <- F_statistic
                model.iv$proj_F_statistic <- proj_F_statistic
            } else {
                ## robust wald test 
                fF <- t(model.iv$coefficients) %*% solve(vcov) %*% model.iv$coefficients
                fF <- c(fF) / dim(x)[2]
                F_statistic <- cbind(fF, dim(x)[2], df, 1 - pf(fF, dim(x)[2], df))
                colnames(F_statistic) <- c("F_statistic", "df1", "df2", "P value")
                model.iv$proj_F_statistic <- F_statistic
            }

            est.coefficients <- cbind(model.iv$coefficients, stderror, Tx, P_t, CI)
            colnames(est.coefficients) <- c("Coef", "Std. Error", "t value", "Pr(>|t|)", "CI_lower", "CI_upper")
            rownames(est.coefficients) <- colnames(x)
            model.iv$est.coefficients <- est.coefficients
        }else{ 
            #clusted fe
            model.iv$cl <- cl
            if (dim(cl)[2] == 1) {
                stderror <- iv.cluster.se(cl = cl, x = x, u.hat = u.hat, dz=dz,
                                          sfe.index = sfe.index, 
                                          cfe.index = cfe.index,
                                          ind = ind, df = df, df.cl = df.cl.use,
                                          model = model,
                                          inv.z.zx = inv.z.zx,
                                          inv.xPzx = inv.xPzx
                                          )
                num.cl <- stderror$num.cl
                vcov <- stderror$vcov.cl
                df.cl <- stderror$df.cl
                dof <- stderror$dof
                stderror <- stderror$stderror
                model.iv$dof <- dof
            }
            else if (dim(cl)[2] == 2) {

                stderror1 <- iv.cluster.se(cl = as.matrix(cl[,1]), x = x, u.hat = u.hat, dz=dz,
                                           sfe.index = sfe.index, 
                                           cfe.index = cfe.index,
                                           ind = ind, df = df, df.cl = df.cl.use[1],
                                           model = model,
                                           inv.z.zx = inv.z.zx,
                                           inv.xPzx = inv.xPzx)
                num.cl1 <- stderror1$num.cl
                vcov1 <- stderror1$vcov.cl
                df.cl1 <- stderror1$df.cl
                dof1 <- stderror1$dof
                stderror1 <- stderror1$stderror

                stderror2 <-  iv.cluster.se(cl = as.matrix(cl[,2]), x = x, u.hat = u.hat, dz=dz,
                                            sfe.index = sfe.index, 
                                            cfe.index = cfe.index,
                                            ind = ind, df = df, df.cl = df.cl.use[2],
                                            model = model,
                                            inv.z.zx = inv.z.zx,
                                            inv.xPzx = inv.xPzx)
                num.cl2 <- stderror2$num.cl
                vcov2 <- stderror2$vcov.cl
                df.cl2 <- stderror2$df.cl
                dof2 <- stderror2$dof
                stderror2 <- stderror2$stderror

                cl_inter <- as.numeric(as.factor(paste(cl[,1], "-:-", cl[,2], sep = "")))

                #if cl_inter is unique; then stderror3 is robust standard error
                cl_inter.count <- table(cl_inter)
                if(max(cl_inter.count)==1){
                    meat3 <- dz * matrix(rep(c(u.hat), dim(dz)[2]), dim(dz)[1], dim(dz)[2])
                    meat3 <- as.matrix(t(meat3) %*% meat3)
                    vcov3 <- inv.xPzx%*%(t(inv.z.zx)%*%meat3%*%inv.z.zx)%*%inv.xPzx
                    vcov3 <- dim(dx)[1]/model.iv$df.original*vcov3
                    stderror3 <- c()
                    for (i in 1:dim(x)[2]) {
                        stderror3 <- c(stderror3, sqrt(vcov3[i, i]))
                    }
                    stderror3 <- as.matrix(stderror3)
                    num.cl3 <- length(cl_inter.count)
                    df.cl3 <- model.iv$df.original
                    dof3 <- model.iv$dof
                }else{
                    stderror3 <- iv.cluster.se(cl = as.matrix(cl_inter), x = x, u.hat = u.hat, dz=dz,
                                            sfe.index = sfe.index, 
                                            cfe.index = cfe.index,
                                            ind = ind, df = df, df.cl = df.cl.use[3],
                                            model = model,
                                            inv.z.zx = inv.z.zx,
                                            inv.xPzx = inv.xPzx)
                    num.cl3 <- stderror3$num.cl
                    vcov3 <- stderror3$vcov.cl
                    df.cl3 <- stderror3$df.cl
                    dof3 <- stderror3$dof
                    stderror3 <- stderror3$stderror
                }

                df.cl <- c(df.cl1,df.cl2,df.cl3)
                num.cl <- min(num.cl1,num.cl2,num.cl3)
                dof.output <- list()
                dof.names <- names(dof1)
                for(sub.dof.name in dof.names){
                    min.index <- which.min(c(dof1[[sub.dof.name]][1]-dof1[[sub.dof.name]][2],
                                           dof2[[sub.dof.name]][1]-dof2[[sub.dof.name]][2],
                                           dof3[[sub.dof.name]][1]-dof3[[sub.dof.name]][2]))
                    if(min.index==1){
                        dof.output[[sub.dof.name]] <- dof1[[sub.dof.name]]
                    }
                    if(min.index==2){
                        dof.output[[sub.dof.name]] <- dof2[[sub.dof.name]]
                    }
                    if(min.index==3){
                        dof.output[[sub.dof.name]] <- dof3[[sub.dof.name]]
                    }
                }

                model.iv$dof <- dof.output
                vcov <- vcov1+vcov2-vcov3
                stderror <- sqrt(stderror1^2 + stderror2^2 - stderror3^2)
            }

            df <- min(num.cl-1,df)
            Tx <- c(model.iv$coefficients)/c(stderror)
            P_t <- NULL
            for (i in 1:length(Tx)) {
                subt <- pt(Tx[i], df)
                P_t <- c(P_t, 2 * min(1 - subt, subt))
            } 

            CI <- cbind(c(model.iv$coefficients) - qt(0.975,df)*c(stderror), c(model.iv$coefficients) + qt(0.975,df)*c(stderror))
            est.coefficients <- cbind(model.iv$coefficients, stderror, Tx, P_t, CI)
            colnames(est.coefficients) <- c("Coef", "Std. Error", "t value", "Pr(>|t|)", "CI_lower", "CI_upper")
            rownames(est.coefficients) <- colnames(x)
            model.iv$est.coefficients <- est.coefficients

            #Wald Test
            
            fF <- try(t(model.iv$coefficients) %*% solve(vcov) %*% model.iv$coefficients,silent = T)
            if ('try-error' %in% class(fF)) {
                fF <- NULL
                F_statistic <- NULL
                warning("Number of clusters insufficient to calculate robust covariance matrix.\n")
            }
            else{
                fF <- c(fF) / dim(x)[2]
                F_statistic <- cbind(fF, dim(x)[2], df, 1 - pf(fF, dim(x)[2], df))
                colnames(F_statistic) <- c("F_statistic", "df1", "df2", "P value")
            }
            
            
            model.iv$proj_F_statistic <- F_statistic
            model.iv$df.residual <- df
            model.iv$df.cl <- df.cl

            df.cl.max <- max(df.cl)
            sigma2 <- t(u.hat)%*%u.hat/df.cl.max
            RMSE <- sqrt(sigma2)
            R2 <- 1 - sum(t(u.hat)%*%u.hat)/sum((data[, 1] - mean(data[, 1]))^2)
            Adj_R2 <- 1 - (1 - R2) * (dim(x)[1] - 1) / df.cl.max
            projR2 <- 1 - sum(t(u.hat)%*%u.hat)/sum((c(model$demeaned$y) - mean(model$demeaned$y))^2)
            projAdj_R2 <- 1 - (1 - projR2) * (dim(x)[1] - 1) / df.cl.max
        }
        colnames(vcov) <- rownames(vcov) <- colnames(x)
        model.iv$vcov <- vcov
        model.iv$RMSE <- RMSE
        model.iv$R2 <- R2
        model.iv$Adj_R2 <- Adj_R2
        model.iv$projR2 <- projR2
        model.iv$projAdj_R2 <- projAdj_R2
    }

    if(need.fe==TRUE){ #recover the fixed effects
        model.y <- SolveFixedEffects(matrix(y,ncol=1), fe$ptr, core.num)
        model.y <- name.fe.model(model.y, inds, fe)
        sfe.x.list <- list()
        cfe.x.list <- list()
        intercept.x <- c()

        if(identical(model.y$sfe.coefs,character(0))==FALSE){
            for(sub.sfe.name in names(model.y$sfe.coefs)){
                sfe.x.list[[sub.sfe.name]] <- model.y$sfe.coefs[[sub.sfe.name]]
            }
        }

        if(identical(model.y$cfe.coefs,character(0))==FALSE){
            for(sub.cfe.name in names(model.y$cfe.coefs)){
                cfe.x.list[[sub.cfe.name]] <- model.y$cfe.coefs[[sub.cfe.name]]
            }
        }

        intercept.x <- model.y$intercept

        if(dim(x)[2]!=0){
            for(i in 1:dim(x)[2]){
                model.x <- SolveFixedEffects(matrix(x[,i],ncol=1), fe$ptr, core.num)
                model.x <- name.fe.model(model.x, inds, fe)

                if(identical(model.y$sfe.coefs,character(0))==FALSE){
                    for(sub.sfe.name in names(model.y$sfe.coefs)){
                        sfe.x.list[[sub.sfe.name]] <- sfe.x.list[[sub.sfe.name]]-iv.coef[i,]*model.x$sfe.coefs[[sub.sfe.name]]
                    }
                }

                if(identical(model.y$cfe.coefs,character(0))==FALSE){
                    for(sub.cfe.name in names(model.y$cfe.coefs)){
                        cfe.x.list[[sub.cfe.name]] <- cfe.x.list[[sub.cfe.name]]-iv.coef[i,]*model.x$cfe.coefs[[sub.cfe.name]]
                    }
                }
                intercept.x <- intercept.x-iv.coef[i,]*model.x$intercept
            }
        }

        
        if(identical(cfe.x.list,list())){
            cfe.x.list <- character(0)
        }
        model.iv$sfe.coefs <- sfe.x.list
        model.iv$cfe.coefs <- cfe.x.list
        model.iv$intercept <- intercept.x
    }
    
    model.iv$y <- y
    model.iv$x <- x
    model.iv$z <- z 
    model.iv$residuals <- u.hat
    model.iv$fitted.values <- y - u.hat
    model.iv$ols.model <- model
    model.iv$reduce.model <- model.reduce
    model.iv$demeaned <- list(dx=dx,dy=dy,dz=dz)
    model.iv$inds <- inds 
    model.iv$matrix <- list(inv.x=inv.x,inv.z=inv.z,inv.z.zx=inv.z.zx,inv.xPzx=inv.xPzx,XZ=XZ,ZY=ZY)

    #all endogenous variables, extracted from x
    en.var.x <- c()
    ex.var.x <- c()
    #from z
    include.iv <- c()
    exclude.iv <- c()
    for(sub.dx in colnames(x)){
        get.ex <- matrix(sapply(as.data.frame(dz), identical, dx[,sub.dx]),nrow = 1)
        colnames(get.ex) <- colnames(z)
        if(any(get.ex)==TRUE){
            include.iv <- c(include.iv,colnames(get.ex)[which(get.ex==TRUE)])
            ex.var.x <- c(ex.var.x,sub.dx)
        }
    }
    include.iv <- unique(include.iv)
    ex.var.x <- unique(ex.var.x)
    exclude.iv <- colnames(dz)[which(!colnames(dz)%in%include.iv)]
    en.var.x <- colnames(dx)[which(!colnames(dx)%in%ex.var.x)]
    model.iv$names <- list(en.var.x=en.var.x,ex.var.x=ex.var.x,include.iv=include.iv,exclude.iv=exclude.iv)

    if(iv.test == TRUE){
        iv.test.results <- iv.test(model.iv, robust=robust, cl=cl, first = first, orthog=orthog, endog=endog)
        model.iv$tests <- iv.test.results
    }

    class(model.iv) <- "ivfastplm"
    return(model.iv)

}

