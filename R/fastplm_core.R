## core function
fastplm.core <- function(y = NULL, ## outcome vector
                         x = NULL, ## covariate matrix
                         ind = NULL, ## indicator matrix
                         se = 0,
                         robust = FALSE,
                         cl = NULL,
                         sfe.index = NULL, ## index for simple fe , a vector of integer
                         cfe.index = NULL, ## index for complex fe, a list of double, e.g. c(effect, influence)
                         core.num = 1,
                         df.use = NULL,
                         df.cl.use = NULL)  {
  
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

    ## ---------- solve the model ------------ ##
    model <- SolveFixedEffects(data, fe$ptr, core.num)
    model <- name.fe.model(model, inds, fe)
    model$inds <- inds
    model$fe <- fe
    ##model$coefficients <- c(model$intercept,model$coefficients)

    ## ------------ inference ---------------- ##
    if (se == 1) {
        
        dx <- model$demeaned$x
        #dx <- cbind(rep(1,dim(dx)[1]),dx)

        res <- model$residuals
        if(is.null(df.use)){
            dof.output <- get.dof(model=model,
                                  x=x,
                                  cl=NULL,
                                  ind=ind,
                                  sfe.index=sfe.index,
                                  cfe.index=cfe.index)
        
            gtot <- dof.output$gtot    
            model$dof <- dof.output$dof
            ## degree of freedom
            df <- dim(y)[1] - gtot
        }else{
            df <- df.use
        }
        
        model$df.original <- df

        ## (x^{\prime}x)^{-1}
        invx <- solvecpp(dx)

        if (is.null(cl)) { ## ols se
            ## estimated sigma hat
            sig2 <- sum(t(res)%*%res)/df
            RMSE <- sqrt(sig2)
            R2 <- 1 - sum(t(res)%*%res)/sum((data[, 1] - mean(data[, 1]))^2)
            Adj_R2 <- 1 - (1 - R2) * (dim(x)[1] - 1) / df
            projR2 <- 1 - sum(t(res)%*%res)/sum((c(model$demeaned$y) - mean(model$demeaned$y))^2)
            projAdj_R2 <- 1 - (1 - projR2) * (dim(x)[1] - 1) / df

            if (robust == FALSE) {
                stderror <- as.matrix(sqrt(sig2*diag(invx)))
                vcov <- as.matrix(c(sig2)*invx)
            } else {
                meat <- dx * matrix(rep(c(res), dim(dx)[2]), dim(dx)[1], dim(dx)[2])
                meat <- t(meat) %*% meat
                vcov <- invx %*% meat %*% invx 
                vcov <- dim(dx)[1] / df * vcov
                stderror <- c()
                for (i in 1:dim(x)[2]) {
                    stderror <- c(stderror, sqrt(vcov[i, i]))
                }
                stderror <- as.matrix(stderror) 
            }



            ## stderror <- as.matrix(sqrt(sig2*diag(invx)))
            Tx <- c(model$coefficients)/c(stderror)

            P_t <- NULL
            for (i in 1:length(Tx)) {
                subt <- pt(Tx[i], df)
                P_t <- c(P_t, 2 * min(1 - subt, subt))
            } 
        
            CI <- cbind(c(model$coefficients) - qt(0.975,df=df)*c(stderror), c(model$coefficients) + qt(0.975,df=df)*c(stderror))

            model$df.residual <- df

            ## F test
            if (robust == FALSE) {
                ## full model
                fF <- R2/(1 - R2) * df/(dim(x)[1] - df - 1)
                ## projection model
                projF <- projR2/(1 - projR2) * df/dim(x)[2]

                F_statistic <- cbind(fF, dim(x)[1] - df - 1, df, 1 - pf(fF, dim(x)[1] - df - 1, df))
                proj_F_statistic <- cbind(projF, dim(x)[2], df, 1 - pf(projF, dim(x)[2], df))

                colnames(F_statistic) <- c("F_statistic", "df1", "df2", "P value")
                colnames(proj_F_statistic) <- c("F_statistic", "df1", "df2", "P value")

                model$F_statistic <- F_statistic
                model$proj_F_statistic <- proj_F_statistic
            } else {
                ## robust wald test 
                fF <- try(t(model$coefficients) %*% solve(vcov) %*% model$coefficients,silent = T)
                if ('try-error' %in% class(fF)) {
                    warning("Exact collinearity in the regression.\n")
                    F_statistic <- NULL
                }
                else{
                    fF <- c(fF) / dim(x)[2]
                    F_statistic <- cbind(fF, dim(x)[2], df, 1 - pf(fF, dim(x)[2], df))
                    colnames(F_statistic) <- c("F_statistic", "df1", "df2", "P value")
                }
                model$proj_F_statistic <- F_statistic
            }

            est.coefficients <- cbind(model$coefficients, stderror, Tx, P_t, CI)
            colnames(est.coefficients) <- c("Coef", "Std. Error", "t value", "Pr(>|t|)", "CI_lower", "CI_upper")

        } else { ## clustered robust se: wald test
            if (dim(cl)[2] == 1) {
                stderror <- cluster.se(cl = cl, x = x, res = res, dx =dx, 
                                       sfe.index = sfe.index, 
                                       cfe.index = cfe.index,
                                       ind = ind, df = df, invx = invx, model = model,
                                       df.cl = df.cl.use)
                num.cl <- stderror$num.cl
                vcov <- stderror$vcov.cl
                df.cl <- stderror$df.cl
                dof <- stderror$dof
                stderror <- stderror$stderror
                model$dof <- dof
            } 
            else if (dim(cl)[2] == 2) {

                stderror1 <- cluster.se(cl = as.matrix(cl[,1]), x = x, res = res, 
                                        dx =dx, sfe.index = sfe.index, 
                                        cfe.index = cfe.index, ind = ind,
                                        df = df, invx = invx, model = model,
                                        df.cl = df.cl.use[1])
                num.cl1 <- stderror1$num.cl
                vcov1 <- stderror1$vcov.cl
                df.cl1 <- stderror1$df.cl
                dof1 <- stderror1$dof
                stderror1 <- stderror1$stderror

                stderror2 <- cluster.se(cl = as.matrix(cl[,2]), x = x, res = res, 
                                        dx =dx, sfe.index = sfe.index, 
                                        cfe.index = cfe.index, ind = ind,
                                        df = df, invx = invx, model = model,
                                        df.cl = df.cl.use[2])
                num.cl2 <- stderror2$num.cl
                vcov2 <- stderror2$vcov.cl
                df.cl2 <- stderror2$df.cl
                dof2 <- stderror2$dof
                stderror2 <- stderror2$stderror

                cl_inter <- as.numeric(as.factor(paste(cl[,1], "-:-", cl[,2], sep = "")))

                #if cl_inter is unique; then stderror3 is robust standard error
                cl_inter.count <- table(cl_inter)
                if(max(cl_inter.count)==1){
                    meat3 <- dx * matrix(rep(c(res), dim(dx)[2]), dim(dx)[1], dim(dx)[2])
                    meat3 <- as.matrix(t(meat3) %*% meat3)
                    vcov3 <- invx %*% meat3 %*% invx 
                    vcov3 <- dim(dx)[1]/model$df.original*vcov3
                    stderror3 <- c()
                    for (i in 1:dim(x)[2]) {
                        stderror3 <- c(stderror3, sqrt(vcov3[i, i]))
                    }
                    stderror3 <- as.matrix(stderror3)
                    num.cl3 <- length(cl_inter.count)
                    df.cl3 <- model$df.original
                    dof3 <- model$dof
                }else{
                    stderror3 <- cluster.se(cl = as.matrix(cl_inter), x = x, res = res, 
                                        dx =dx, sfe.index = sfe.index, 
                                        cfe.index = cfe.index, ind = ind,
                                        df = df, invx = invx, model = model,
                                        df.cl = df.cl.use[3])
                    num.cl3 <- stderror3$num.cl
                    vcov3 <- stderror3$vcov.cl
                    df.cl3 <- stderror3$df.cl
                    dof3 <- stderror3$dof
                    stderror3 <- stderror3$stderror
                }

                #dof.output <- list()
                #dof.output[[colnames(cl)[1]]]=dof1
                #dof.output[[colnames(cl)[2]]]=dof2
                #dof.output[[paste0(colnames(cl)[1],"*",colnames(cl)[2])]]=dof3
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

                model$dof <- dof.output
                vcov <- vcov1+vcov2-vcov3
                stderror <- sqrt(stderror1^2 + stderror2^2 - stderror3^2)
            }

            ##
            df <- min(num.cl-1,df)
            Tx <- c(model$coefficients)/c(stderror)
            P_t <- NULL
            for (i in 1:length(Tx)) {
                subt <- pt(Tx[i], df)
                P_t <- c(P_t, 2 * min(1 - subt, subt))
            } 

            CI <- cbind(c(model$coefficients) - qt(0.975,df)*c(stderror), c(model$coefficients) + qt(0.975,df)*c(stderror))
            est.coefficients <- cbind(model$coefficients, stderror, Tx, P_t, CI)
            colnames(est.coefficients) <- c("Coef", "Std. Error", "t value", "Pr(>|t|)", "CI_lower", "CI_upper")
        
            #Wald Test
            fF <- try(t(model$coefficients) %*% solve(vcov) %*% model$coefficients,silent = T)
            if ('try-error' %in% class(fF)) {
                fF <- NULL
                F_statistic <- NULL
                warning("Number of clusters insufficient to calculate robust covariance matrix or exact collinearity in the regression.\n")
            }else{
                fF <- c(fF) / dim(x)[2]
                F_statistic <- cbind(fF, dim(x)[2], df, 1 - pf(fF, dim(x)[2], df))
                colnames(F_statistic) <- c("F_statistic", "df1", "df2", "P value")
            }
            
            model$proj_F_statistic <- F_statistic
            model$df.residual <- df
            model$df.cl <- df.cl

            df.cl.max <- max(df.cl)
            sig2 <- sum(t(res)%*%res)/df.cl.max
            RMSE <- sqrt(sig2)
            R2 <- 1 - sum(t(res)%*%res)/sum((data[, 1] - mean(data[, 1]))^2)
            Adj_R2 <- 1 - (1 - R2) * (dim(x)[1] - 1) / df.cl.max
            projR2 <- 1 - sum(t(res)%*%res)/sum((c(model$demeaned$y) - mean(model$demeaned$y))^2)
            projAdj_R2 <- 1 - (1 - projR2) * (dim(x)[1] - 1) / df.cl.max
        }
        colnames(vcov) <- colnames(x)
        rownames(vcov) <- colnames(x)
        model$vcov <- vcov      
        model$R2 <- R2
        model$Adj_R2 <- Adj_R2
        model$projR2 <- projR2
        model$projAdj_R2 <- projAdj_R2
        model$RMSE <- RMSE
        rownames(est.coefficients) <- colnames(x)
        model$est.coefficients <- est.coefficients
      
    }

    ## class(model) <- "feModel"
    ## class(model) <- "fastplm"
    model
}




