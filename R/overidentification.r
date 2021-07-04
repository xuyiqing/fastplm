    get.over.identify <- function(dy,
                                  dx,
                                  dz,
                                  robust = 0,
                                  cl=NULL,
                                  df.cl=NULL,
                                  u.hat,
                                  inv.z=NULL,
                                  en.var.x,
                                  exclude.iv,
                                  orthog=NULL){
        if(is.null(orthog)==FALSE){
            if(length(orthog)>(length(exclude.iv)-length(en.var.x))){
                stop("Collinearity/identification problems in orthog conditions.\n")
            }
            for(var in orthog){
                if(!var%in%colnames(dz)){
                    stop(paste0(var," is not an excluded instrument nor an included instrument variables.\n"))
                }
            }
            new.dx <- dx
            new.dz <- dz[,colnames(dz)[which(!colnames(dz)%in%orthog)]]
            new.dz <- matrix(new.dz,ncol=dim(dz)[2]-length(orthog))
            colnames(new.dz) <- colnames(dz)[which(!colnames(dz)%in%orthog)]
        }
        over.identify <- list()
        if(is.null(cl)){
            if(robust==FALSE){ #Sargan's Statistic
                #egmm.coef <- model.iv$coefficients
            }else{ #Hansen's Test
                #meat <- dz * matrix(rep(c(u.hat), dim(dz)[2]), dim(dz)[1], dim(dz)[2])
                #meat <- t(meat) %*% meat
                #meat <- meat*dim(dz)[1]/sub.df
                #egmm.coef <- solve(t(dx)%*%dz%*%solve(meat)%*%t(dz)%*%dx)%*%t(dx)%*%dz%*%solve(meat)%*%t(dz)%*%dy
                #hansen.meat <- meat
                #colnames(hansen.meat) <- rownames(hansen.meat) <- colnames(dz)
                #u.hat.egmm <- dy - dx%*%egmm.coef
                get.gmm <- solvegmm(Y=dy,X=dx,Z=dz, u_hat=u.hat)
                hansen.meat <- meat <- get.gmm[['meat']]
                inv.hansen.meat <- get.gmm[['invmeat']]
                colnames(inv.hansen.meat) <- rownames(inv.hansen.meat) <- colnames(hansen.meat) <- rownames(hansen.meat) <- colnames(dz)
                egmm.coef <- get.gmm[['coef']]
                u.hat.gmm <- get.gmm[['u_hat']]
            }
        }else{
            if (dim(cl)[2] == 1) {
                stderror <- cluster.se(cl = cl, x = NULL, res = u.hat, dx =dz, 
                                       sfe.index = NULL, 
                                       cfe.index = NULL,
                                       ind = NULL, df = NULL, invx = inv.z,
                                       df.cl = df.cl)
                meat <- stderror$meat
            }
            if (dim(cl)[2] == 2) {
                stderror1 <- cluster.se(cl = as.matrix(cl[,1]), x = NULL, res = u.hat, dx =dz, 
                                        sfe.index = NULL, 
                                        cfe.index = NULL,
                                        ind = NULL, df = NULL, invx = inv.z,
                                        df.cl = df.cl[1])
                meat1 <- stderror1$meat

                stderror2 <- cluster.se(cl = as.matrix(cl[,2]), x = NULL, res = u.hat, dx =dz, 
                                        sfe.index = NULL, 
                                        cfe.index = NULL,
                                        ind = NULL, df = NULL, invx = inv.z,
                                        df.cl = df.cl[2])
                meat2 <- stderror2$meat

                cl_inter <- as.numeric(as.factor(paste(cl[,1], "-:-", cl[,2], sep = "")))
                #if cl_inter is unique; then stderror3 is robust standard error
                cl_inter.count <- table(cl_inter)

                if(max(cl_inter.count)==1){
                    meat3 <- dz * matrix(rep(c(u.hat), dim(dz)[2]), dim(dz)[1], dim(dz)[2])
                    meat3 <- as.matrix(t(meat3) %*% meat3)
                }else{
                    stderror3 <- cluster.se(cl = as.matrix(cl_inter), x = NULL, res = u.hat, dx =dz, 
                                            sfe.index = NULL, 
                                            cfe.index = NULL,
                                            ind = NULL, df = NULL, invx = inv.z,
                                            df.cl = df.cl[3])
                    meat3 <- stderror3$meat
                }
                meat <- meat1 + meat2 - meat3
            }
            #egmm.coef <- solve(t(dx)%*%dz%*%solve(meat)%*%t(dz)%*%dx)%*%t(dx)%*%dz%*%solve(meat)%*%t(dz)%*%dy
            #hansen.meat <- meat
            #colnames(hansen.meat) <- rownames(hansen.meat) <- colnames(dz)
            #u.hat.egmm <- dy - dx%*%egmm.coef
            meat <- as.matrix(meat)
            get.gmm <- solvegmm_meat(Y=dy,X=dx,Z=dz,meat=meat)
            hansen.meat <- meat <- get.gmm[['meat']]
            inv.hansen.meat <- get.gmm[['invmeat']]
            colnames(inv.hansen.meat) <- rownames(inv.hansen.meat) <- colnames(hansen.meat) <- rownames(hansen.meat) <- colnames(dz)
            egmm.coef <- get.gmm[['coef']]
            u.hat.gmm <- get.gmm[['u_hat']]
        }

        if(length(en.var.x)<length(exclude.iv)){
            #overidentification
            if(is.null(cl)){
                if(robust==FALSE){ #Sargan's Statistic
                    sargan <- (t(u.hat)%*%dz%*%inv.z%*%t(dz)%*%u.hat/(t(u.hat)%*%u.hat))*dim(dz)[1]
                    p.sargan <- pchisq(sargan, df=length(exclude.iv)-length(en.var.x), lower.tail=FALSE)
                    sargan.statistic <- c(sargan,length(exclude.iv)-length(en.var.x),p.sargan)
                    names(sargan.statistic) <- c("Sargan Statistic","df","P value")
                    over.identify$sargan <- sargan.statistic
                }else{ #Hansen's Test
                    J.hansen <- get.gmm[["hansen"]]
                    p.hansen <- pchisq(J.hansen, df=length(exclude.iv)-length(en.var.x), lower.tail=FALSE)
                    hansen.statistic <- c(J.hansen,length(exclude.iv)-length(en.var.x),p.hansen)
                    names(hansen.statistic) <- c("Hansen Statistic","df","P value")
                    over.identify$hansen <- hansen.statistic
                }
            }else{
                J.hansen <- get.gmm[["hansen"]]
                p.hansen <- pchisq(J.hansen, df=length(exclude.iv)-length(en.var.x), lower.tail=FALSE)
                hansen.statistic <- c(J.hansen,length(exclude.iv)-length(en.var.x),p.hansen)
                names(hansen.statistic) <- c("Hansen Statistic","df","P value")
                over.identify$hansen <- hansen.statistic
            }

            if(is.null(orthog)==FALSE){ # get C statistics
                                        # use old sigma^2 & S
                #new.inv.z <- solve(t(new.dz)%*%new.dz)
                #new.Pz <- new.dz%*%new.inv.z%*%t(new.dz)
                #new.inv.xPzx <- solve(t(new.dx)%*%new.Pz%*%new.dx)
                #new.inv.z.zx <- new.inv.z%*%t(new.dz)%*%new.dx
                #new.iv.coef <- matrix(new.inv.xPzx%*%t(new.dx)%*%new.Pz%*%dy,ncol=1)
                #rownames(new.iv.coef) <- colnames(new.dx)
                #new.u.hat <- matrix(dy-new.dx%*%new.iv.coef,ncol=1)

                iv.solve.new <- solveiv(Y=dy,X=new.dx,Z=new.dz)
                new.u.hat <- iv.solve.new[["residual"]]
                new.iv.coef <- iv.solve.new[["coef"]]
                rownames(new.iv.coef) <- colnames(new.dx)
                new.inv.z <- iv.solve.new[["invzz"]]

                if(is.null(cl)){
                    if(robust==FALSE){ #Sargan's Statistic
                        new.sargan <- (t(new.u.hat)%*%new.dz%*%new.inv.z%*%t(new.dz)%*%new.u.hat/(t(u.hat)%*%u.hat))*dim(dz)[1] #use original sigma^2 (u.hat)
                        C.value <- sargan-new.sargan
                        C.p.value <- pchisq(C.value, df=length(orthog), lower.tail=FALSE)
                        C.statistic <- c(C.value,length(orthog),C.p.value)
                        names(C.statistic) <- c("C Statistic","df","P value")
                        over.identify$C.statistic <- C.statistic
                    }else{ #Hansen's Test 
                           #Use original meat
                        new.meat <- as.matrix(hansen.meat[colnames(new.dz),colnames(new.dz)])
                        new.get.gmm <- solvegmm_meat(Y=dy,X=new.dx,Z=new.dz,meat=new.meat)
                        new.J.hansen <- new.get.gmm[['hansen']]
                        #new.egmm.coef <- solve(t(new.dx)%*%new.dz%*%solve(new.meat)%*%t(new.dz)%*%new.dx)%*%t(new.dx)%*%new.dz%*%solve(new.meat)%*%t(new.dz)%*%dy
                        #new.u.hat.egmm <- dy - new.dx%*%new.egmm.coef
                        #new.J.hansen <- t(new.u.hat.egmm)%*%new.dz%*%solve(new.meat)%*%t(new.dz)%*%new.u.hat.egmm
                        C.value <- J.hansen-new.J.hansen
                        C.p.value <- pchisq(C.value, df=length(orthog), lower.tail=FALSE)
                        C.statistic <- c(C.value,length(orthog),C.p.value)
                        names(C.statistic) <- c("C Statistic","df","P value")
                        over.identify$C.statistic <- C.statistic
                    }
                }else{
                    new.meat <- as.matrix(hansen.meat[colnames(new.dz),colnames(new.dz)])
                    new.get.gmm <- solvegmm_meat(Y=dy,X=new.dx,Z=new.dz,meat=new.meat)
                    new.J.hansen <- new.get.gmm[['hansen']]
                    #new.egmm.coef <- solve(t(new.dx)%*%new.dz%*%solve(new.meat)%*%t(new.dz)%*%new.dx)%*%t(new.dx)%*%new.dz%*%solve(new.meat)%*%t(new.dz)%*%dy
                    #new.u.hat.egmm <- dy - new.dx%*%new.egmm.coef
                    #new.J.hansen <- t(new.u.hat.egmm)%*%new.dz%*%solve(new.meat)%*%t(new.dz)%*%new.u.hat.egmm
                    C.value <- J.hansen-new.J.hansen
                    C.p.value <- pchisq(C.value, df=length(orthog), lower.tail=FALSE)
                    C.statistic <- c(C.value,length(orthog),C.p.value)
                    names(C.statistic) <- c("C Statistic","df","P value")
                    over.identify$C.statistic <- C.statistic
                }
            }
        }
        return(over.identify)
    }

    
    


