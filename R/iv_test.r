iv.test <- function(model.iv, robust=TRUE, cl=NULL, first = TRUE, orthog=NULL, endog=NULL){
    #if(class(model.iv)!="ivfastplm"){
    #    stop("Must be an ivfastplm object.")
    #}

    dy <- model.iv$demeaned$dy
    dx <- model.iv$demeaned$dx
    dz <- model.iv$demeaned$dz
    iv.name <- model.iv$names
    
    y <- model.iv$y
    x <- model.iv$x
    z <- model.iv$z

    u.hat <- model.iv$residuals
    iv.matrix <- model.iv$matrix

    en.var.x <- iv.name$en.var.x
    ex.var.x <- iv.name$ex.var.x
    include.iv <- iv.name$include.iv
    exclude.iv <- iv.name$exclude.iv

    # first stage
    if(length(en.var.x)==1){
        multiple.endo <- 0
    }else{
        multiple.endo <- 1
    }
    pos.z.test <- which(colnames(z)%in%exclude.iv)
    
    # regress all endogenous x on all excluded variables
    inv.z <- iv.matrix$inv.z # inv(z'z)
    inv.xPzx <- iv.matrix$inv.xPzx
    

    if(is.null(cl)){
        sub.df.residual <- model.iv$df.residual + dim(x)[2] - dim(z)[2]
        sub.df <- model.iv$df.original + dim(x)[2] - dim(z)[2]
        sub.df.cl <- sub.df
    }else{
        sub.df.residual <- model.iv$df.residual
        sub.df <- model.iv$df.original
        sub.df.cl <- model.iv$df.cl + dim(x)[2] - dim(z)[2]
    }


    first.stage <- list()
    if(length(en.var.x)>0 && first==TRUE){ #first stage information   
        
        en.var.x.matrix <- matrix(dx[,en.var.x],ncol=length(en.var.x))
        #en.var.x.coef <- inv.z%*%t(dz)%*%en.var.x.matrix
        #colnames(en.var.x.coef) <- en.var.x
        #en.var.x.hat <- dz%*%en.var.x.coef
        #colnames(en.var.x.hat) <- en.var.x

        sub.lm <- simpleols(Y=en.var.x.matrix,X=dz)
        en.var.x.coef <- sub.lm[['coef']]
        colnames(en.var.x.coef) <- en.var.x
        en.var.x.hat <- sub.lm[['Y_hat']]
        colnames(en.var.x.hat) <- en.var.x

        # F test
        for(i in 1:length(en.var.x)){
            target.en.var <- en.var.x[i]
            other.en.var <- en.var.x[-i]
            
            sub.dx <- matrix(dx[,target.en.var],ncol=1)
            sub.coef <- matrix(en.var.x.coef[,target.en.var],ncol=1)
            rownames(sub.coef) <- colnames(z)
            sub.res <- matrix(sub.dx-en.var.x.hat[,target.en.var],ncol=1)
            
            sub.df.cl <- sub.df
            sub.sigma2 <- t(sub.res)%*%sub.res/(sub.df)
            sub.RMSE <- sqrt(sub.sigma2)
            if(is.null(cl)){
                if (robust == FALSE) {
                    sub.vcov <- as.matrix(c(t(sub.res)%*%sub.res/sub.df)*inv.z)
                } else{
                    meat <- dz * matrix(rep(c(sub.res), dim(dz)[2]), dim(dz)[1], dim(dz)[2])
                    meat <- t(meat) %*% meat
                    vcov <- inv.z%*%meat%*%inv.z
                    sub.vcov <- dim(dz)[1]/sub.df*vcov
                }
            }
            else{
                cl <- model.iv$cl
                sub.df.cl <- model.iv$df.cl + dim(x)[2] - dim(z)[2]
                if (dim(cl)[2] == 1) {
                    stderror <- cluster.se(cl = cl, x = z, res = sub.res, dx =dz, 
                                        sfe.index = NULL, 
                                        cfe.index = NULL,
                                        ind = NULL, df = NULL, invx = inv.z,
                                        df.cl = sub.df.cl)
                    sub.vcov <- stderror$vcov.cl
                }
                if (dim(cl)[2] == 2) {
                    stderror1 <- cluster.se(cl = as.matrix(cl[,1]), x = z, res = sub.res, dx =dz, 
                                        sfe.index = NULL, 
                                        cfe.index = NULL,
                                        ind = NULL, df = NULL, invx = inv.z,
                                        df.cl = sub.df.cl[1])
                    vcov1 <- stderror1$vcov.cl

                    stderror2 <- cluster.se(cl = as.matrix(cl[,2]), x = z, res = sub.res, dx =dz, 
                                        sfe.index = NULL, 
                                        cfe.index = NULL,
                                        ind = NULL, df = NULL, invx = inv.z,
                                        df.cl = sub.df.cl[2])
                    vcov2 <- stderror2$vcov.cl

                    cl_inter <- as.numeric(as.factor(paste(cl[,1], "-:-", cl[,2], sep = "")))

                    #if cl_inter is unique; then stderror3 is robust standard error
                    cl_inter.count <- table(cl_inter)
                    if(max(cl_inter.count)==1){
                        meat3 <- dz * matrix(rep(c(sub.res), dim(dz)[2]), dim(dz)[1], dim(dz)[2])
                        meat3 <- as.matrix(t(meat3) %*% meat3)
                        vcov3 <- inv.z %*% meat3 %*% inv.z 
                        vcov3 <- dim(dz)[1]/sub.df*vcov3
                    }else{
                        stderror3 <- cluster.se(cl = as.matrix(cl_inter), x = z, res = sub.res, dx =dz, 
                                                sfe.index = NULL, 
                                                cfe.index = NULL,
                                                ind = NULL, df = NULL, invx = inv.z,
                                                df.cl = sub.df.cl[3])
                        vcov3 <- stderror3$vcov.cl
                    }
                    sub.vcov <- vcov1+vcov2-vcov3
                }
            }

            colnames(sub.vcov) <- rownames(sub.vcov) <- colnames(z)
            sub.use.coef <- matrix(sub.coef[exclude.iv,],ncol=1)
            sub.use.vcov <- as.matrix(sub.vcov[exclude.iv,exclude.iv])
            partial.F <- t(sub.use.coef)%*%solve(sub.use.vcov)%*%sub.use.coef/length(exclude.iv)

            partial.F.statistic <- cbind(partial.F, length(exclude.iv), sub.df.residual, 1 - pf(partial.F, length(exclude.iv), sub.df.residual))
            colnames(partial.F.statistic) <- c("partial F statistics", "df1", "df2", "P value")
            sub.model <- list(coefficients=sub.coef,vcov=sub.vcov,df.original=sub.df,df.cl=sub.df.cl,df.residual=sub.df.residual,
                            partial.F=partial.F.statistic,RMSE=sub.RMSE)


            if(multiple.endo==0){
                AP.F.statistic <- SW.F.statistic <- partial.F.statistic
                sub.model$AP.F.statistic <- AP.F.statistic
                sub.model$SW.F.statistic <- SW.F.statistic
            }

            if(multiple.endo==1){
                #AP F.statistics
                target.dx.hat <- matrix(en.var.x.hat[,target.en.var],ncol=1)
                target.dx <- matrix(dx[,target.en.var],ncol=1)
                
                if(!is.null(include.iv)){
                    other.dx.hat <- cbind(matrix(en.var.x.hat[,other.en.var],ncol=length(other.en.var)),
                                          matrix(dx[,include.iv],ncol=length(include.iv)))

                    other.dx <- cbind(matrix(dx[,other.en.var],ncol=length(other.en.var)),
                                      matrix(dx[,include.iv],ncol=length(include.iv)))
                }else{
                    other.dx.hat <- matrix(en.var.x.hat[,other.en.var],ncol=length(other.en.var))
                    other.dx <- matrix(dx[,other.en.var],ncol=length(other.en.var))
                }
            
                #inv.other.dx.hat <- solve(t(other.dx.hat)%*%other.dx.hat)
                #other.coef.hat <- matrix(inv.other.dx.hat%*%t(other.dx.hat)%*%target.dx,ncol=1)
                #e1 <- target.dx - other.dx.hat%*%other.coef.hat

                sub.lm.APSW <- simpleols(Y=target.dx,X=other.dx.hat)
                other.coef.hat <- sub.lm.APSW[['coef']]
                e1 <- sub.lm.APSW[['u_hat']]

                #second stage
                AP.sub.coef <- matrix(inv.z%*%t(dz)%*%e1,ncol=1)
                rownames(AP.sub.coef) <- colnames(z)
                AP.sub.res <- matrix(e1-dz%*%AP.sub.coef,ncol=1)
                AP.sub.sigma2 <- t(AP.sub.res)%*%AP.sub.res/(sub.df)
                AP.sub.RMSE <- sqrt(AP.sub.sigma2)
                if(is.null(cl)){
                    if (robust == FALSE) {
                        AP.sub.vcov <- as.matrix(c(t(AP.sub.res)%*%AP.sub.res/sub.df)*inv.z)
                    } else{
                        AP.meat <- dz * matrix(rep(c(AP.sub.res), dim(dz)[2]), dim(dz)[1], dim(dz)[2])
                        AP.meat <- t(AP.meat) %*% AP.meat
                        AP.vcov <- inv.z%*%AP.meat%*%inv.z
                        AP.sub.vcov <- dim(dz)[1]/sub.df*AP.vcov
                    }
                }
                else{
                    if (dim(cl)[2] == 1) {
                        stderror <- cluster.se(cl = cl, x = z, res = AP.sub.res, dx =dz, 
                                        sfe.index = NULL, 
                                        cfe.index = NULL,
                                        ind = NULL, df = NULL, invx = inv.z,
                                        df.cl = sub.df.cl)
                        AP.sub.vcov <- stderror$vcov.cl
                    }
                    if (dim(cl)[2] == 2) {
                        stderror1 <- cluster.se(cl = as.matrix(cl[,1]), x = z, res = AP.sub.res, dx =dz, 
                                        sfe.index = NULL, 
                                        cfe.index = NULL,
                                        ind = NULL, df = NULL, invx = inv.z,
                                        df.cl = sub.df.cl[1])
                        AP.vcov1 <- stderror1$vcov.cl

                        stderror2 <- cluster.se(cl = as.matrix(cl[,2]), x = z, res = AP.sub.res, dx =dz, 
                                        sfe.index = NULL, 
                                        cfe.index = NULL,
                                        ind = NULL, df = NULL, invx = inv.z,
                                        df.cl = sub.df.cl[2])
                        AP.vcov2 <- stderror2$vcov.cl

                        if(max(cl_inter.count)==1){
                            meat3 <- dz * matrix(rep(c(AP.sub.res), dim(dz)[2]), dim(dz)[1], dim(dz)[2])
                            meat3 <- as.matrix(t(meat3) %*% meat3)
                            AP.vcov3 <- inv.z %*% meat3 %*% inv.z 
                            AP.vcov3 <- dim(dz)[1]/sub.df*vcov3
                        }else{
                            stderror3 <- cluster.se(cl = as.matrix(cl_inter), x = z, res = AP.sub.res, dx =dz, 
                                                sfe.index = NULL, 
                                                cfe.index = NULL,
                                                ind = NULL, df = NULL, invx = inv.z,
                                                df.cl = sub.df.cl[3])
                            AP.vcov3 <- stderror3$vcov.cl
                        }
                        AP.sub.vcov <- AP.vcov1+AP.vcov2-AP.vcov3
                    }
                }

                colnames(AP.sub.vcov) <- rownames(AP.sub.vcov) <- colnames(z)
                AP.sub.use.coef <- matrix(AP.sub.coef[exclude.iv,],ncol=1)
                AP.sub.use.vcov <- as.matrix(AP.sub.vcov[exclude.iv,exclude.iv])
                Fdf1 <- (length(exclude.iv)-length(en.var.x)+1)
                AP.partial.F <- t(AP.sub.use.coef)%*%solve(AP.sub.use.vcov)%*%AP.sub.use.coef/Fdf1
                AP.partial.F.statistic <- cbind(AP.partial.F, Fdf1, sub.df.residual, 1 - pf(AP.partial.F, Fdf1, sub.df.residual))
                colnames(AP.partial.F.statistic) <- c("Angrist-Pischke F Value", "df1", "df2", "P value")
                sub.model$AP.F.statistic <- AP.partial.F.statistic

                #second stage
                e2 <- target.dx - other.dx%*%other.coef.hat
                SW.sub.coef <- matrix(inv.z%*%t(dz)%*%e2,ncol=1)
                rownames(SW.sub.coef) <- colnames(z)
                SW.sub.res <- matrix(e2-dz%*%SW.sub.coef,ncol=1)
                SW.sub.sigma2 <- t(SW.sub.res)%*%SW.sub.res/(sub.df)
                SW.sub.RMSE <- sqrt(SW.sub.sigma2)
                if(is.null(cl)){
                    if (robust == FALSE) {
                        SW.sub.vcov <- as.matrix(c(t(SW.sub.res)%*%SW.sub.res/sub.df)*inv.z)
                    } else{
                        SW.meat <- dz * matrix(rep(c(SW.sub.res), dim(dz)[2]), dim(dz)[1], dim(dz)[2])
                        SW.meat <- t(SW.meat) %*% SW.meat
                        SW.vcov <- inv.z%*%SW.meat%*%inv.z
                        SW.sub.vcov <- dim(dz)[1]/sub.df*SW.vcov
                    }
                }
                else{
                    if (dim(cl)[2] == 1) {
                        stderror <- cluster.se(cl = cl, x = z, res = SW.sub.res, dx =dz, 
                                        sfe.index = NULL, 
                                        cfe.index = NULL,
                                        ind = NULL, df = NULL, invx = inv.z,
                                        df.cl = sub.df.cl)
                        SW.sub.vcov <- stderror$vcov.cl
                    }
                    if (dim(cl)[2] == 2) {
                        stderror1 <- cluster.se(cl = as.matrix(cl[,1]), x = z, res = SW.sub.res, dx =dz, 
                                        sfe.index = NULL, 
                                        cfe.index = NULL,
                                        ind = NULL, df = NULL, invx = inv.z,
                                        df.cl = sub.df.cl[1])
                        SW.vcov1 <- stderror1$vcov.cl

                        stderror2 <- cluster.se(cl = as.matrix(cl[,2]), x = z, res = SW.sub.res, dx =dz, 
                                        sfe.index = NULL, 
                                        cfe.index = NULL,
                                        ind = NULL, df = NULL, invx = inv.z,
                                        df.cl = sub.df.cl[2])
                        SW.vcov2 <- stderror2$vcov.cl

                        if(max(cl_inter.count)==1){
                            meat3 <- dz * matrix(rep(c(SW.sub.res), dim(dz)[2]), dim(dz)[1], dim(dz)[2])
                            meat3 <- as.matrix(t(meat3) %*% meat3)
                            SW.vcov3 <- inv.z %*% meat3 %*% inv.z 
                            SW.vcov3 <- dim(dz)[1]/sub.df*vcov3
                        }else{
                            stderror3 <- cluster.se(cl = as.matrix(cl_inter), x = z, res = SW.sub.res, dx =dz, 
                                                sfe.index = NULL, 
                                                cfe.index = NULL,
                                                ind = NULL, df = NULL, invx = inv.z,
                                                df.cl = sub.df.cl[3])
                            SW.vcov3 <- stderror3$vcov.cl
                        }
                        SW.sub.vcov <- SW.vcov1+SW.vcov2-SW.vcov3
                    }
                }

                colnames(SW.sub.vcov) <- rownames(SW.sub.vcov) <- colnames(z)
                SW.sub.use.coef <- matrix(SW.sub.coef[exclude.iv,],ncol=1)
                SW.sub.use.vcov <- as.matrix(SW.sub.vcov[exclude.iv,exclude.iv])
                
                
                SW.partial.F <- t(SW.sub.use.coef)%*%solve(SW.sub.use.vcov)%*%SW.sub.use.coef/Fdf1
                SW.partial.F.statistic <- cbind(SW.partial.F, Fdf1, sub.df.residual, 1 - pf(SW.partial.F, Fdf1, sub.df.residual))
                colnames(SW.partial.F.statistic) <- c("Sanderson-Windmeijer F Value", "df1", "df2", "P value")
                sub.model$SW.F.statistic <- SW.partial.F.statistic
            }
            first.stage[[en.var.x[i]]] <- sub.model
        }
    }

    #Anderson-Rubin Wald test
    if(length(en.var.x)>0){
        dy <- model.iv$demeaned$dy
        AR.sub.coef <- matrix(inv.z%*%t(dz)%*%dy,ncol=1)
        rownames(AR.sub.coef) <- colnames(z)
        AR.sub.res <- matrix(dy-dz%*%AR.sub.coef,ncol=1)
        AR.sub.sigma2 <- t(AR.sub.res)%*%AR.sub.res/(sub.df)
        AR.sub.RMSE <- sqrt(AR.sub.sigma2)
        if(is.null(cl)){
            if (robust == FALSE) {
                AR.sub.vcov <- as.matrix(c(t(AR.sub.res)%*%AR.sub.res/sub.df)*inv.z)
            } else{
                    AR.meat <- dz * matrix(rep(c(AR.sub.res), dim(dz)[2]), dim(dz)[1], dim(dz)[2])
                    AR.meat <- t(AR.meat) %*% AR.meat
                    AR.vcov <- inv.z%*%AR.meat%*%inv.z
                    AR.sub.vcov <- dim(dz)[1]/sub.df*AR.vcov
                }
        }
        else{
            if (dim(cl)[2] == 1) {
                stderror <- cluster.se(cl = cl, x = z, res = AR.sub.res, dx =dz, 
                                    sfe.index = NULL, 
                                    cfe.index = NULL,
                                    ind = NULL, df = NULL, invx = inv.z,
                                    df.cl = sub.df.cl)
                 AR.sub.vcov <- stderror$vcov.cl
            }
            if (dim(cl)[2] == 2) {
                stderror1 <- cluster.se(cl = as.matrix(cl[,1]), x = z, res = AR.sub.res, dx =dz, 
                                        sfe.index = NULL, 
                                        cfe.index = NULL,
                                        ind = NULL, df = NULL, invx = inv.z,
                                        df.cl = sub.df.cl[1])
                AR.vcov1 <- stderror1$vcov.cl

                stderror2 <- cluster.se(cl = as.matrix(cl[,2]), x = z, res = AR.sub.res, dx =dz, 
                                        sfe.index = NULL, 
                                        cfe.index = NULL,
                                        ind = NULL, df = NULL, invx = inv.z,
                                        df.cl = sub.df.cl[2])
                AR.vcov2 <- stderror2$vcov.cl

                cl_inter <- as.numeric(as.factor(paste(cl[,1], "-:-", cl[,2], sep = "")))

                #if cl_inter is unique; then stderror3 is robust standard error
                cl_inter.count <- table(cl_inter)

                if(max(cl_inter.count)==1){
                            meat3 <- dz * matrix(rep(c(AR.sub.res), dim(dz)[2]), dim(dz)[1], dim(dz)[2])
                            meat3 <- as.matrix(t(meat3) %*% meat3)
                            AR.vcov3 <- inv.z %*% meat3 %*% inv.z 
                            AR.vcov3 <- dim(dz)[1]/sub.df*AR.vcov3
                }else{
                            stderror3 <- cluster.se(cl = as.matrix(cl_inter), x = z, res = AR.sub.res, dx =dz, 
                                                sfe.index = NULL, 
                                                cfe.index = NULL,
                                                ind = NULL, df = NULL, invx = inv.z,
                                                df.cl = sub.df.cl[3])
                            AR.vcov3 <- stderror3$vcov.cl
                }
                AR.sub.vcov <- AR.vcov1+AR.vcov2-AR.vcov3
            }
        }
        colnames(AR.sub.vcov) <- rownames(AR.sub.vcov) <- colnames(z)
        AR.sub.use.coef <- matrix(AR.sub.coef[exclude.iv,],ncol=1)
        AR.sub.use.vcov <- as.matrix(AR.sub.vcov[exclude.iv,exclude.iv])
        AR.partial.F <- t(AR.sub.use.coef)%*%solve(AR.sub.use.vcov)%*%AR.sub.use.coef/(length(exclude.iv))
        AR.partial.F.statistic <- cbind(AR.partial.F, length(exclude.iv), sub.df.residual, 1 - pf(AR.partial.F, length(exclude.iv), sub.df.residual))
        colnames(AR.partial.F.statistic) <- c("Anderson-Rubin F Value", "df1", "df2", "P value")
        first.stage$AR.F.value <- AR.partial.F.statistic
    }
        
    #overidentification
    #en.var.x <- iv.name$en.var.x
    #ex.var.x <- iv.name$ex.var.x
    #include.iv <- iv.name$include.iv
    #exclude.iv <- iv.name$exclude.iv

    if(is.null(cl)){
        df.cl <- df.residual <- model.iv$df.residual
    }else{
        df.cl <- model.iv$df.cl
        df.residual <- model.iv$df.residual
    }
    u.hat <- model.iv$residuals

    over.identify <- try(get.over.identify(dy=dy,
                                           dx=dx,
                                           dz=dz,
                                           cl=cl,
                                           robust = robust,
                                           df.cl=df.cl,
                                           u.hat=u.hat,
                                           inv.z=inv.z,
                                           en.var.x=en.var.x,
                                           exclude.iv=exclude.iv,
                                           orthog=orthog),silent=TRUE)
    if('try-error' %in% class(over.identify)) {
        over.identify <- NULL
        warning("Collinearity/identification problems in overidentification test.\n")
    }                                   
                                                          
    #endogenous
    #en.var.x <- iv.name$en.var.x
    #ex.var.x <- iv.name$ex.var.x
    #include.iv <- iv.name$include.iv
    #exclude.iv <- iv.name$exclude.iv
    endogenous.out <- list()
    
    if(length(en.var.x)>0){
        if(is.null(endog)){ #test all instrumented covariate
            endog <- en.var.x
        }else{
            #endog must belong to en.var.x
            for(var in endog){
                if(!var%in%en.var.x){
                    stop(paste0(var," must belong to all endogenous variables."))
                }
            }
        }
        
        #put endogenous variables in z
        dx.endog <- dx
        dy.endog <- dy
        dz.endog <- as.matrix(cbind(dz,matrix(dx[,endog],ncol=length(endog))))
        colnames(dz.endog) <- c(colnames(dz),endog)

        iv.solve.endog <- solveiv(Y=dy.endog,X=dx.endog,Z=dz.endog)
        u.hat.endog <- iv.solve.endog[["residual"]]
        inv.z.endog <- iv.solve.endog[["invzz"]]
        
        orthog.endog <- endog
        en.var.x.endog <- en.var.x[which(!en.var.x%in%endog)]
        exclude.iv.endog <- c(exclude.iv,en.var.x[which(en.var.x%in%endog)])


        endog.out <- try(get.over.identify(dy=dy.endog,
                                    dx=dx.endog,
                                    dz=dz.endog,
                                    cl=cl,
                                    robust = robust,
                                    df.cl=df.cl,
                                    u.hat=u.hat.endog,
                                    inv.z=inv.z.endog,
                                    en.var.x=en.var.x.endog,
                                    exclude.iv=exclude.iv.endog,
                                    orthog=orthog.endog),silent=TRUE) 
        
        if('try-error' %in% class(endog.out)) {
            endog.out <- NULL
            endogenous <- NULL
            warning("Collinearity/identification problems in endogenous test.\n")
        }else{
            endogenous <- endog.out$C.statistic
            names(endogenous) <- c("DWH Value","df","P Value") 
            endogenous.out[['DWH']] <- endogenous
        }     
    }


    output <- list(first.stage=first.stage,over.identify=over.identify,endogenous=endogenous.out)
    return(output)
}