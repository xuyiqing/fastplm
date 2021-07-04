#1, degree of freedom calculation
#2, singletons detection and drop
#3, clustered standard error

get.dof <- function(model,
                    x,
                    cl=NULL,
                    ind,
                    sfe.index,
                    cfe.index){
    df.info <- list()                   
    #calculate the degree of freedom
    gtot <- 0 #total loss of degree of freedom
    if(is.null(sfe.index)){
        for (i in 1:length(model$inds$levels)) {
            cat.num <- length(model$inds$levels[[i]]) 
            redund.num <- 1
            out.cum <- c(cat.num,redund.num)
            names(out.cum) <- c("Category","Redundant")
            df.info[[i]] <- out.cum
        }
    }else{
        for (i in 1:length(sfe.index)) {
            cat.num <- length(model$inds$levels[[sfe.index[i]]]) 
            redund.num <- 1
            out.cum <- c(cat.num,redund.num)
            names(out.cum) <- c("Category","Redundant")
            df.info[[model$fe$sfe.names[i]]] <- out.cum
        }
    }

    ## we test nest among the first two sets of fixed effects, that is "firstpair" option in reghdfe
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
        c.level1 <- length(model$inds$levels[[pos1]])
        c.level2 <- length(model$inds$levels[[pos2]])
            
        if (c.level1 > c.level2) {
            fe.level <- table(unlist(tapply(level1, level2, unique)))
        } else {
            fe.level <- table(unlist(tapply(level2, level1, unique)))
        }

        if (max(fe.level) == 1) {
            if(c.level1 > c.level2){
                df.info[[model$fe$sfe.names[1]]]['Redundant'] <- c.level2
            }
            if(c.level1 <= c.level2){
                df.info[[model$fe$sfe.names[2]]]['Redundant'] <- c.level1
            }
            #adj.dof <- min(c.level1, c.level2) - 1
        }
    }

    if (!is.null(cfe.index)) {
            sum.cfe.coef <- NULL 
            for (i in 1:length(cfe.index)) {
                sum.cfe.coef <- sum(c(model$cfe.coefs[[i]])) 
                if (sum.cfe.coef <= 1e-7) {
                    cat.num <- length(model$inds$levels[[cfe.index[[i]][1]]])
                    redund.num <- 1
                    out.cum <- c(cat.num,redund.num)
                    names(out.cum) <- c("Category","Redundant")
                    df.info[[model$fe$cfe.names[i]]] <- out.cum
                } else {
                    cat.num <- length(model$inds$levels[[cfe.index[[i]][1]]])
                    redund.num <- 0
                    out.cum <- c(cat.num,redund.num)
                    names(out.cum) <- c("Category","Redundant")
                    df.info[[model$fe$cfe.names[i]]] <- out.cum
                }
            }
    }

    if(is.null(cl)==FALSE){
        raw.cl <- as.numeric(as.factor(cl))
        level.cl <- unique(raw.cl)
        
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
            correct.index <- which.max(contain.cl)
            if(is.null(sfe.index)){
                df.info[[correct.index]]["Redundant"] <- df.info[[i]]["Category"]  
            }else{
                df.info[[model$fe$sfe.names[correct.index]]]["Redundant"] <- df.info[[model$fe$sfe.names[correct.index]]]["Category"]  
            }
        }

        if (length(contain.cl) >= 2) {
            if (max(contain.cl[c(1,2)]) == 1) {
                cluster.adj <- 1 
            }
        }

        ## check fixed effects nested within clusters
        if (is.null(sfe.index)) {
            for (i in 1:dim(ind)[2]) {
                fe.level <- as.numeric(table(unlist(tapply(ind[,i], c(cl), unique))))
                #if (length(fe.level) != length(level.cl)) {
                    if (sum(fe.level == 1) == length(fe.level)) {
                        df.info[[i]]['Redundant'] <- df.info[[i]]['Category']
                    }
                #}   
            }
        } else {
            for (i in 1:length(sfe.index)) {
                fe.level <- as.numeric(table(unlist(tapply(ind[,sfe.index[i]], c(cl), unique))))
                #if (length(fe.level) != length(level.cl)) {
                    if (sum(fe.level == 1) == length(fe.level)) {
                        df.info[[model$fe$sfe.names[i]]]['Redundant'] <- df.info[[model$fe$sfe.names[i]]]['Category']
                    }
                #}
            }
        }

        if (!is.null(cfe.index)) {
            for (i in 1:length(cfe.index)) {
                fe.level <- unlist(tapply(ind[,cfe.index[[i]][1]], c(cl), unique))
                #if (length(fe.level) != length(level.cl)) {
                    if (sum(fe.level == 1) == length(fe.level)) {
                        df.info[[model$fe$cfe.names[i]]]['Redundant'] <- df.info[[model$fe$cfe.names[i]]]['Category']
                    }
                #}
            }
        }
    }

    gtot<-0
    for(element in df.info){
        gtot <- gtot + element[1]-element[2]
    }

    if(is.null(x)){
        gtot <- gtot + 1
    }else{
        gtot <- gtot + 1 + dim(x)[2]
    }
    names(gtot) <- NULL
    output <- list(gtot=gtot,dof=df.info)
    return(output)

}

have.singletons <- function(ind=NULL){ #Matrix
    for(k in 1:dim(ind)[2]){
        if(min(table(ind[,k]))==1){
            return(TRUE)
        }
    }
    return(FALSE)            
} 

drop.singletons <- function(ind=NULL){ #Matrix
    to.check <- ind
    colnames(to.check) <- NULL
    keep.all <- c()
    result <- apply(to.check, 2, sub.check)   
                      
    for(i in 1:length(result)){
        if(!is.null(result[[i]])){
            tar.all <- result[[i]][which(result[[i]][,2]==1),1]
            tar.index <- which(!to.check[,i]%in%c(tar.all))
            keep.all <- c(keep.all,tar.index)
        }
    }

    keep.all <- unique(keep.all)
    return(keep.all)                      
}

sub.check <- function(col){
    t <- aggregate(col, by=list(col), FUN=length)
    if(min(t[,2])==1){
      return(t)
    }
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
                       df.cl = NULL,
                       invx,
                       model = NULL) {

    ## replicate degree of freedom 
    df.Redundant <- 0

    raw.cl <- as.numeric(as.factor(cl))
    level.cl <- unique(raw.cl)
    #meat <- matrix(0, dim(dx)[2], dim(dx)[2])
    #get.sub.meat <- function(cl.index){
    #    sub.id <- which(raw.cl == cl.index)
    #   if(length(sub.id)>1){
    #        hmeat <- as.matrix(t(as.matrix(dx[sub.id,])) %*% as.matrix(res[sub.id,]))
    #    }
    #    if(length(sub.id)==1){
    #        hmeat <- t(as.matrix(t(as.matrix(dx[sub.id,]))*res[sub.id]))
    #    }
    #    return(hmeat%*%t(hmeat))
    #}
    #asx <- lapply(level.cl,get.sub.meat)
    #meat <- Reduce(`+`, asx)
    cal.cl <- 0
    if(is.null(df.cl)){
        cal.cl <- 1
        dof.output.cl <- get.dof(model=model,
                                 x=x,
                                 cl=cl,
                                 ind=ind,
                                 sfe.index=sfe.index,
                                 cfe.index=cfe.index)
        gtot <- dof.output.cl$gtot    
        df.cl <- length(raw.cl)-gtot
    }
    q <- length(level.cl)/(length(level.cl) - 1) * (dim(dx)[1] - 1)/(df.cl)

    
    clustercpp.output <- clustercpp(rawcl = matrix(raw.cl,ncol=1),X = dx,Res = res, q = q,invX = invx)
    vcov <- clustercpp.output[["vcov"]]
    meat <- clustercpp.output[["meat"]]
    stderror <- as.matrix(sqrt(diag(vcov)))
    
    
    #vcov <- q*invx %*% meat %*% invx
    #stderror <- sqrt(q) * as.matrix(sqrt(diag(invx %*% meat %*% invx)))

    if(cal.cl==1){
       out <- list(stderror=stderror,num.cl=length(level.cl),meat=meat,vcov.cl=vcov,df.cl=df.cl,dof=dof.output.cl$dof) 
    }
    else{
       out <- list(stderror=stderror,num.cl=length(level.cl),meat=meat,vcov.cl=vcov,df.cl=df.cl,dof=NULL) 
    }
    return(out)
}

## subfunction for iv clustered se
iv.cluster.se <- function(cl, 
                       x,
                       u.hat,
                       dz,
                       sfe.index,
                       cfe.index,
                       ind,
                       df,
                       df.cl = NULL,
                       inv.xPzx,
                       inv.z.zx,
                       model) {

    ## replicate degree of freedom 
    df.Redundant <- 0

    raw.cl <- as.numeric(as.factor(cl))
    level.cl <- unique(raw.cl)

    #meat <- matrix(0, dim(dz)[2], dim(dz)[2])
    #get.sub.meat <- function(cl.index){
    #    sub.id <- which(raw.cl == cl.index)
    #    if(length(sub.id)>1){
    #        hmeat <- as.matrix(t(as.matrix(dz[sub.id,])) %*% as.matrix(u.hat[sub.id,]))
    #    }
    #    if(length(sub.id)==1){
    #        hmeat <- t(as.matrix(t(as.matrix(dz[sub.id,]))*u.hat[sub.id]))
    #    }
    #    return(hmeat%*%t(hmeat))
    #}
    #asx <- lapply(level.cl,get.sub.meat)
    #meat <- Reduce(`+`, asx)

    cal.cl <- 0
    if(is.null(df.cl)){
        cal.cl <- 1
        dof.output.cl <- get.dof(model=model,
                                 x=x,
                                 cl=cl,
                                 ind=ind,
                                 sfe.index=sfe.index,
                                 cfe.index=cfe.index)
        gtot <- dof.output.cl$gtot    
        df.cl <- length(raw.cl)-gtot
    }
    q <- length(level.cl)/(length(level.cl) - 1) * (dim(x)[1] - 1)/(df.cl)

    clustercpp.output <- ivclustercpp(rawcl = matrix(raw.cl,ncol=1),X = dz,Res = u.hat, q = q,invxPzx=inv.xPzx, invzzx=inv.z.zx)
    vcov <- clustercpp.output[["vcov"]]
    meat <- clustercpp.output[["meat"]]
    stderror <- as.matrix(sqrt(diag(vcov)))
    

    #vcov <- inv.xPzx%*%(t(inv.z.zx)%*%meat%*%inv.z.zx)%*%inv.xPzx
    #vcov <- q*vcov
    #stderror <- as.matrix(sqrt(diag(vcov)))

    if(cal.cl==1){
       out <- list(stderror=stderror,num.cl=length(level.cl),meat=meat,vcov.cl=vcov,df.cl=df.cl,dof=dof.output.cl$dof) 
    }
    else{
       out <- list(stderror=stderror,num.cl=length(level.cl),meat=meat,vcov.cl=vcov,df.cl=df.cl,dof=NULL) 
    }
    return(out)
}
