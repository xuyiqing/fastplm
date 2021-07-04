## generic function
fastplm <- function(data =NULL,
                    formula = NULL, ## 
                    index = NULL, ## index name    
                    y = NULL, ## outcome vector
                    x = NULL, ## covariate matrix or endogenous variables in IV regression
                    z = NULL, ## all exogenous variables, if not null, use IV regression or GMM
                    ind = NULL, ## indicator matrix
                    sfe = NULL, ## index for simple fe , a vector of integer
                    cfe = NULL, ## index for complex fe, a list of double, e.g. c(effect, influence)
                    PCA = TRUE,
                    sp = NULL,
                    knots = NULL, ## b-spline
                    degree = 3, ## order of curve
                    se = 1, ## uncertainty estimates for covariates
                    vce = "robust", ## standard, robust, clustered, bootstrap, jackknife
                    cluster = NULL, ## cluster name
                    wild = TRUE,
                    refinement = FALSE,
                    test_x = NULL,
                    parallel = TRUE, 
                    nboots = 200, ## bootstrap number
                    seed = NULL,
                    core.num = 1,
                    drop.singletons = FALSE,
                    bootcluster = NULL,
                    iv.test = TRUE,
                    first = FALSE,
                    orthog = NULL,
                    endog = NULL) {
    
    cluster.level <- NULL

    ## formula and data are present
    ## omit y, x and ind
    if(!is.null(formula)){
        varnames <- all.vars(formula) #contains y,x

        if(is.null(data)){
            stop("\"data\" should be specified when \"formula\" is present.\n")
        }

        if(is.null(index)){
            stop("\"index\" should be specified when \"formula\" is present.\n")
        }

        if (!is.null(cluster)) {
            if (class(cluster)[1] != "character") {
                stop("\"cluster\" should be (a) character value(s) when \"formula\" is present..\n")
            }
        }
        if (!is.null(sp)) {
            if (class(sp)[1] != "character") {
                stop("\"sp\" should be a character value when \"formula\" is present..\n")
            }
        }
        if (!is.null(z)) {
            if (class(z)[1] != "character") {
                stop("\"z\" should be (a) character value(s) when \"formula\" is present.\n")
            }
        }

        data.name <- names(data)
        all.name <- unique(c(varnames,index,cluster,sp,z))
        for (i in 1:length(all.name)) {
            if (!all.name[i] %in% data.name) {
                stop(paste("variable ", all.name[i], " is not in the dataset.\n", sep = ""))
            }
        }
        data <- data[,unique(c(varnames,index,cluster,sp,z))]
        if (sum(is.na(data)) > 0) {
            data <- na.omit(data)
        }
        
        ## bootcluster must be one of 2 clusters
        if(length(cluster)==2){
            if(!is.null(bootcluster)){
                if(!bootcluster%in%cluster){
                    stop("bootcluster must be one of the cluster names or \"interaction\".\n")
                }
            }
        }

        y <- as.matrix(data[, varnames[1]])
        x <- NULL
        if (length(varnames) >= 2) {
            x <- as.matrix(data[, varnames[2:length(varnames)]])
            colnames(x) <- varnames[2:length(varnames)]
        }

        if(!is.null(z)){
            z.name <- z
            z <- as.matrix(data[,z.name])
            colnames(z) <- z.name
        }

        ind <- as.matrix(data[,index])
        ind <- matrix(ind ,nrow=dim(data)[1])
        colnames(ind) <- index

        if(is.null(sfe)){
            sfe <- colnames(ind)                
        }

        ## save cluster name
        cluster.level <- cluster
        if (!is.null(cluster)) {
            cluster <- as.matrix(data[, cluster])
        }
        if (!is.null(sp)) {
            sp.name <- sp
            if(sp%in%colnames(x)){
                stop("sp shouldn't be one of regressors.\n")
            }
            if(is.null(z)==FALSE){
                if(sp%in%z){
                    stop("sp shouldn't be one of exogenous variables.\n")
                }
            }
            sp <- as.matrix(data[, sp])
            colnames(sp) <- sp.name
        }
    }
    else if (is.null(formula) && !is.null(data) && is.null(y) && is.null(x) && is.null(z)) {
        if (sum(is.na(data)) > 0) {
            data <- na.omit(data)
        }
        p <- dim(data)[2]
        y <- as.matrix(data[, 1])
        x <- NULL
        if (p >= 2) {
            x <- as.matrix(data[, 2:p])
        }
        if(is.null(ind)==TRUE){
            stop("\"ind\" should be specified if \"formula\" is absent.\n")
        }

        if(is.null(colnames(ind))){
            colnames(ind) <- paste0("FE.",c(1:dim(ind)[2]))
        }

        if(is.null(sfe)){
            index <- colnames(ind)
            sfe <- colnames(ind)                
        } 
    }
    else if(is.null(y)==FALSE){

        if(is.null(ind)==TRUE){
            stop("\"ind\" should be specified if \"formula\" is absent.\n")
        }

        if(class(y)[1]!='matrix'){
            stop("y should be a matrix.\n")
        }

        if(!is.null(x)){
            if(class(x)[1]!='matrix'){
                stop("x should be a matrix.\n")
            }            
        }

        if(class(ind)[1]!='matrix'){
            stop("ind should be a matrix.\n")
        }

        if(!is.null(x)){
            if(dim(x)[1]==dim(y)[1] && dim(x)[1]==dim(ind)[1]){
                ##
            }else{
                stop("y,x and ind don't have the same number of observations.\n")
            }            
        }
        
        if(is.null(colnames(ind))){
            colnames(ind) <- paste0("FE.",c(1:dim(ind)[2]))
        }

        if(is.null(sfe)){
            index <- colnames(ind)
            sfe <- colnames(ind)                
        }

        if(is.null(x)==FALSE){
            if(is.null(colnames(x))){
                colnames(x) <- paste0("x.",c(1:dim(x)[2]))
            }
        }

        if(is.null(colnames(y))){
            y <- matrix(y,ncol=1)
            colnames(y) <- 'y'
        }

        if(is.null(z)==FALSE){
            if(class(z)[1]!='matrix'){
                stop("z should be a matrix.\n")
            }
            if(dim(z)[1]!=dim(y)[1]){
                stop("y,x,z and ind don't have the same number of observations.\n")
            }
            if(is.null(colnames(z))){
                colnames(z) <- paste0("z.",c(1:dim(z)[2]))
            }
        }

        if(is.null(data)){
            data <- cbind(y,x,z,ind)
            colnames(data) <- c(colnames(y),colnames(x),colnames(z),colnames(ind))
        }
    }


    ## number of covariates
    p.old <- p <- ifelse(is.null(x), 0, dim(x)[2])
    L.old <- L <- ifelse(is.null(z), 0, dim(z)[2])

    if(p>L && L!=0){ #
        stop("Underidentification in IV regression. The number of exogenous variables is less than the number of endogenous variables.")
    }

    if(L!=0){
        iv.reg <- 1
    }else{
        iv.reg <- 0 
    }

    #if (sum(is.na(y)) > 0) {
    #    stop("Missing values in dependent variable.\n")
    #}
    #if (p > 0) {
    #    if (sum(is.na(x)) > 0) {
    #        stop("Missing values in covariates.\n")
    #    }
    #}
    ## ------ b-spline ----- ##
    ## bsp <- 0
    sp.name <- NULL
    if (!is.null(sp)) {
        if (!class(sp)[1] %in% c("matrix", "numeric", "integer")) {
            stop("\"sp\" must be a numeric vector.\n")
        }
        if (class(sp)[1] != "matrix") {
            sp <- as.matrix(sp)
        }
        sp.name <- colnames(sp)
        ## bsp <- 1
    }

    sp.raw <- sp.matrix <- NULL
    if (!is.null(sp)) {
        get.sp <- try(bs(c(sp), knots = knots, degree = degree, 
                      intercept = FALSE), silent = TRUE)

        if ('try-error' %in% class(get.sp)) {
            stop("Cannot generate matrix for b-spline curve.\n")
        } else {
            sp.matrix <- as.matrix(get.sp)
            ## gen reference level 
            ref.pos <- !duplicated(c(sp))
            sp.ref <- c(sp)[ref.pos]
            sp.matrix.ref <- sp.matrix[ref.pos,]
            sp.raw <- list(sp.ref = sp.ref, sp.matrix.ref = sp.matrix.ref)
        }

        sp.name2 <- c()
        for (i in 1:dim(sp.matrix)[2]) {
            sp.name2 <- c(sp.name2, paste("sp", i, sep = ""))
        }
        cov.name <- c(colnames(x), sp.name2)
        x <- cbind(x, sp.matrix)
        colnames(x) <- cov.name
        p <- dim(x)[2]

        if(L>0){
           cov.name.z <- c(colnames(z),sp.name2)
            z <- cbind(z,sp.matrix)
            colnames(z) <- cov.name.z
            L <- dim(z)[2] 
        }  
    }

    ## index 
    if (is.null(ind)) {
        stop("No indicators.\n")
    } else {
        if (class(ind)[1] != "matrix") {
            stop("\"ind\" must be a matrix.\n")
        }
    }
    
    #drop singletons
    #iteratively
    num.singletons <- 0
    use.index <- as.matrix(c(1:dim(ind)[1]))
    drop.index <- NULL
    if(drop.singletons==TRUE){
        dim.old <- dim(ind)[1]
        ind.index <- as.matrix(c(1:dim(ind)[1]))
        while(have.singletons(ind=ind)==TRUE){
            keep.index <- drop.singletons(ind=ind)
            y <- y[keep.index,,drop=FALSE]

            if(!is.null(x)){
                x <- x[keep.index,,drop=FALSE]
            }

            if(!is.null(z)){
                z <- z[keep.index,,drop=FALSE] 
            }

            if(!is.null(cluster)){
                cluster <- cluster[keep.index,,drop=FALSE]
            }
            
            if(!is.null(sp)){
                sp <- sp[keep.index,,drop=FALSE]
            }
            
            if(!is.null(ind)){
                ind <- ind[keep.index,,drop=FALSE]
            }  
            ind.index <- ind.index[keep.index,,drop=FALSE]          
        }
        dim.new <- dim(y)[1]
        num.singletons <- dim.old-dim.new
        use.index <- ind.index
        drop.index <- setdiff(c(1:dim.old),ind.index)
        if(dim.new[1] < dim.old[1]){
            cat(paste0("Drop ",num.singletons," Singletons.\n"))
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
            ##cat("No covariates.\n")
            ##se <- 0
            se <- 0
            bootstrap <- 0
            refinement <- 0
        } else {
            if (p.old == 0) {
                refinement <- 0
            }
        }
        if (!vce %in% c("standard", "clustered", "robust", "bootstrap", "jackknife")) {
            stop("Choose \" vce \" from c(\" standard \", \" robust \", \" clustered \", \" jackknife \", \" bootstrap \").\n")
        } else {
            if (vce == "standard") {
                #cat("Please consider clustering the standard errors or using block bootstraps.\n")
                robust <- 0
            } 
            else if (vce == "robust") {
                robust <- 1
            }
            else if (vce %in% c("bootstrap", "jackknife")) {
                bootstrap <- 1
            }
        }
    } else {
        refinement <- 0
        iv.test <- FALSE
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
            #if (wild == TRUE) {
            #    stop("For cluster wild bootstrap, please only specify one cluster variable.\n")
            #}
            #if (refinement == TRUE) {
            #    stop("For cluster bootstrap refinement, please only specify one cluster variable.\n")
            #}
            n.cluster <- c(length(unique(cluster[,1])), length(unique(cluster[,2]))) 
        }
        else {
            n.cluster <- c(length(unique(c(cluster))))
        }
        if (vce %in% c("standard", "robust")) {
            vce <- "clustered"
        }
    } else {
        #if (refinement == TRUE) {
        #    cat("No cluster variable. Cluster bootstrap refinement will not be performed.\n")
        #    refinement <- 0
        #}
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
        #if (is.null(test_x)) {
        #    stop("For wild bootstrap refinement, variables of interest should be specified.\n")
        #}
        if (length(test_x) >= 2) {
            stop("Please specify only one variable of interest.\n")
        }
        #if (class(test_x) == "character") {
        #    pos <- which(colnames(data)[2:(p.old+1)] == test_x)
        #} else {
        pos <- test_x
        #}
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
    raw.influence.pca <- NULL
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
            colnames(raw.influence) <- NULL
            raw.influence.pca <- prcomp(raw.influence)
            pc.influence <- raw.influence.pca$x ## pca
            ## qr.influence <- qr(raw.influence)
            ## pc.influence <- qr.Q(qr.influence)
            ## reference for raw and principal com influence
            get.ref <- function(mat) {return(as.matrix(mat[!duplicated(mat),],,2))}
            raw_pc.ref <- lapply(1:length(raw.influence.index), 
                                 function(i)get.ref(cbind(raw.influence[,i], pc.influence[,i])))

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
            if (vce == "jackknife") {
                variance.type <- "Jackknife"
            } else {
                if (is.null(cluster)) {
                    if (wild == 0) {
                        variance.type <- "Pairs Bootstrap"
                        if(refinement == 1){
                            variance.type <- "Robust"
                            if(is.null(pos)==TRUE){
                                refinement.type <- "Unrestricted Percentile-t Pairs Bootstrap Refinement"
                            }else{
                                refinement.type <- "Restricted Percentile-t Pairs Bootstrap Refinement"
                            }
                        }
                    } else {
                        variance.type <- "Wild Bootstrap"
                        if(refinement == 1){
                            variance.type <- "Robust"
                            if(is.null(pos)==TRUE){
                                refinement.type <- "Unrestricted Percentile-t Wild Bootstrap Refinement"
                            }else{
                                refinement.type <- "Restricted Percentile-t Wild Bootstrap Refinement"
                            }
                        }
                    }
                } else {
                    if (wild == 0) {
                        variance.type <- "Clustered Pairs Bootstrap"
                        if(refinement == 1){
                            variance.type <- "Clustered"
                            if(is.null(pos)==TRUE){
                                refinement.type <- "Unrestricted Percentile-t Clustered Pairs Bootstrap Refinement"
                            }else{
                                refinement.type <- "Restricted Percentile-t Clustered Pairs Bootstrap Refinement"
                            }
                        }
                    } else {
                        variance.type <- "Wild Clustered Bootstrap"
                        if(refinement == 1){
                            variance.type <- "Clustered"
                            if(is.null(pos)==TRUE){
                                refinement.type <- "Unrestricted Percentile-t Clustered Wild Bootstrap Refinement"
                            }else{
                                refinement.type <- "Restricted Percentile-t Clustered Wild Bootstrap Refinement"
                            }
                        }
                    }
                }
            }        
        }
    }

    if(iv.reg==0){
        if (bootstrap == 0) {
            model <- fastplm.core(y = y, x = x, ind = ind, 
                                sfe.index = sfe.index, cfe.index = cfe.index, 
                                se = se, robust = robust,
                                cl = cluster, core.num = core.num)
        } else {
            #if (parallel == TRUE) {
            #    para.clusters <- makeCluster(core.num)
            #    registerDoParallel(para.clusters)
            #    cat("Parallel computing ...\n")
            #}
            jackknife <- ifelse(vce == "jackknife", 1, 0)
            model <- fastplm.boot(seed = seed, y = y, x = x, ind = ind, 
                                sfe.index = sfe.index, cfe.index = cfe.index, robust = robust,
                                cluster = cluster, parallel = parallel,
                                wild = wild, jackknife = jackknife,
                                refinement = refinement, pos = pos,  
                                nboots = nboots, core.num = core.num,
                                bootcluster = bootcluster)
            #if (parallel == TRUE) {
            #    stopCluster(para.clusters)
            #}
        }
    }

    if(iv.reg==1){
        if (bootstrap == 0) {
            model <- ivfastplm.core(y = y, ## outcome vector
                                    x = x, ## regressor matrix
                                    z = z, ## exogenous variables matrix
                                    ind = ind, ## indicator matrix
                                    se = se,
                                    robust = robust,
                                    cl = cluster,
                                    sfe.index = sfe.index, ## index for simple fe , a vector of integer
                                    cfe.index = cfe.index,
                                    core.num = core.num,
                                    need.fe = TRUE,
                                    iv.test = iv.test,
                                    first = first,
                                    orthog = orthog,
                                    endog = endog) 
        }else{
            jackknife <- ifelse(vce == "jackknife", 1, 0)
            model <- ivfastplm.boot(seed=seed, 
                                    y = y, 
                                    x = x, 
                                    z = z,
                                    ind = ind, 
                                    sfe.index = sfe.index, 
                                    cfe.index = cfe.index, 
                                    cluster = cluster,
                                    robust = robust,
                                    parallel = parallel,
                                    wild = wild, 
                                    jackknife = jackknife,
                                    refinement = refinement, 
                                    pos = pos,
                                    nboots = nboots, 
                                    bootcluster = bootcluster, #2-way clusters in wild (refinement)
                                    core.num = core.num)
        }
    }
    

    ## complex fixed effect
    model$cfe.index <- cfe.index ## a numeric list
    model$raw_pc.ref <- raw_pc.ref ## a ref list for complex fe, used for prediction
    model$PCA <- PCA
    model$raw.influence.pca <- raw.influence.pca
    #print(model$coefficients)
    if (p > 0 & !is.null(colnames(data))) {

        x.name <- NULL
        if (p.old > 0) {
            x.name <- colnames(data)[2:(p.old+1)]
        }
        
        sp.length <- p - p.old
        sp.name2 <- NULL
        if (sp.length > 0) {
            for (i in 1:sp.length) {
                sp.name2 <- c(sp.name2, paste("sp", i, sep = ""))
            }
        }
        cov.name <- c(x.name, sp.name2)
        rownames(model$coefficients) <- cov.name
    }

    if (!is.null(index)) {
        model$inds$effect.names <- index
    }

    model <- c(model, list(call = match.call(),
                           N = dim(y)[1],
                           n.cluster = n.cluster,
                           variance.type = variance.type,
                           refinement.type = refinement.type,
                           num.singletons = num.singletons,
                           drop.index = drop.index,
                           sp.raw = sp.raw,
                           sp.name = sp.name,
                           cluster.level = cluster.level,
                           first=first,
                           endog=endog,
                           orthog=orthog))
    class(model) <- "fastplm"
    
    return(model)

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
predict.fastplm <- function(object, data = NULL, x = NULL, 
                                    sp = NULL, ind = NULL, ...) {
    model <- object
    CHECK.INPUT(model, "model", "fastplm")

    sp.raw <- model$sp.raw
    sp.ref <- sp.matrix.ref <- NULL
    p.sp <- 0
    if (!is.null(sp.raw)) {
        sp.ref <- sp.raw$sp.ref
        sp.matrix.ref <- sp.raw$sp.matrix.ref
        p.sp <- dim(sp.matrix.ref)[2]
    }
    
    
    if (!is.null(data)) { ## receive a dataframe
        
        coef <- model$coefficients
        if (is.null(rownames(coef))) {
            stop("Cannot receive a data frame. Try matrix.\n")
        }
        p <- 0
        if (!is.null(coef)) {
            p <- length(c(coef))
        }

        x <- NULL
        if (!is.null(coef)) {

            if (p == p.sp) {
                x <- NULL
            } else {
                x.name <- rownames(coef)[1:(p - p.sp)]
                if (sum(x.name %in% names(data)) < length(x.name)) {
                    stop("The dataset doesn\'t contain some covariates in the model.")
                } else {
                    x <- as.matrix(data[, x.name])
                }
            }
            
            sp.matrix <- NULL
            sp.name <- model$sp.name
            if (!is.null(sp.name)) {
                if (!is.null(sp.raw)) {
                    sp.ind <- data[, sp.name]
                    sp.pos <- sapply(1:length(sp.ind), function(vec){return(which(sp.ref == sp.ind[vec]))})
                    sp.matrix <- sp.matrix.ref[sp.pos, ]
                }
            }
            x <- cbind(x, sp.matrix)
            if (dim(x)[2] != p) {
                stop("The number of covariates should be the same as that in the model.")
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
        }
        sp.matrix <- NULL
        if (!is.null(sp)) {
            sp <- c(sp)
            if (!is.null(sp.raw)) {
                sp.pos <- sapply(1:length(sp), function(vec){return(which(sp.ref == sp[vec]))})
                sp.matrix <- sp.matrix.ref[sp.pos, ]
            }

        }
        x <- cbind(x, sp.matrix)
        coef <- model$coefficients 
        if (!is.null(x)) {
            ASSERT.MATRIX.DIM(x, "x", length(model$coefficients), is.width = TRUE)
        }
        if (dim(ind)[2] != length(model$inds$effect.names)) {
            stop("The number of indicators should be the same as that in the model.")
        }
    }
    
    if (!is.null(model$cfe.index) && model$PCA == TRUE) {
        sub.raw.influence.index <- unique(sapply(1:length(model$cfe.index), function(i) model$cfe.index[[i]][2]))
        sub.raw.influence <- as.matrix(ind[,sub.raw.influence.index])
        #sub.pc.influence <- matrix(NA, dim(sub.raw.influence)[1], dim(sub.raw.influence)[2])
        sub.pc.influence <- predict(model$raw.influence.pca, sub.raw.influence)
        raw_pc.ref <- model$raw_pc.ref
        #for (i in 1:dim(sub.raw.influence)[2]) {
            ## doesn't allow extrapolation here
            ## correct this!
        #    sub.ref <- raw_pc.ref[[i]]
        #    for (j in 1:dim(sub.ref)[1]) {
        #        sub.pc.influence[which(sub.raw.influence[,i] == sub.ref[j,1]),i] <- sub.ref[j,2]
        #    }
        #}
        if (sum(is.na(c(sub.pc.influence))) > 0) {
            stop("There may be some unspecific levels within some fixed effects in the model.\n")
        }
        ind[,sub.raw.influence.index] <- sub.pc.influence
    }

    #inds <- create.subindicators(sub.inds = ind, model = model)
    #print(inds)
    #if (!is.null(x)) {
    #    ASSERT.MATRIX.DIM(x, "x", length(model$coefficients), is.width = TRUE)
    #}
    #if (inds$parent != model$inds$uid)
    #    stop(ERR.FE.predict.invalid.inds())

    get.eff <- function(index,effs){
        index <- as.character(index)
        if(index %in% rownames(effs)){
            out <- effs[index,]
            names(out) <- NULL
            return(out)
        }
        else{
            return(NA)
        }
    }

    if (!is.null(x)) {
        y <- x %*% model$coefficients + model$intercept
    } 
    else {
        y <- as.matrix(rep(model$intercept, dim(ind)[1]))
    }
    fe <- model$fe

    
    for (col in SEQ(1,length(fe$sfes))) {
        effs <- model$sfe.coefs[[col]]
        sum <- rep(NA,dim(y)[1])
        use.index <- which(as.character(ind[,col]) %in% rownames(effs))
        sum[use.index] <- effs[as.character(ind[use.index,col]),]
        #sum  <- sapply(ind[,col],function(x) get.eff(x,effs))
        y <- y + sum
    }

    for (i in SEQ(1, length(fe$cfe.effs))) {
        eff.col <- fe$cfe.effs[i]
        inf.col <- fe$cfe.infs[i]
        effs    <- model$cfe.coefs[[i]]
        rownames(effs) <- model$inds$levels[[eff.col]]
        #sum <- sapply(ind[,eff.col],function(x) get.eff(x,effs))
        sum <- rep(NA,dim(y)[1])
        use.index <- which(as.character(ind[,eff.col]) %in% rownames(effs))
        sum[use.index] <- effs[as.character(ind[use.index,eff.col]),]
        sum <- sum * ind[,inf.col]
        y   <- y + sum
    }

    return(y)
}

## Print
summary.fastplm <- function(object,
                            detail = FALSE,  
                            ...) {
    

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
                cat("\n* Degree of freedom may be too conservative due to redundant parameters.\n")
            }
        }
    }

    if (!is.null(object$refinement$wald.refinement)) {
        cat("\n\n---")
        cat(object$refinement.type)
        cat("\n")
        print(round(object$refinement$wald.refinement, 3))
    }  

    if(!is.null(object$dof)){
        cat("\n---\n") 
        cat("Absorbed degrees of freedom:\n")
        for(sub.dof.name in names(object$dof)){
            print.out <- object$dof[[sub.dof.name]]
            print.append <- print.out[1]-print.out[2]
            print.out <- matrix(c(print.out,print.append),nrow=1)
            rownames(print.out) <- sub.dof.name
            colnames(print.out) <- c("Categories","-Redundant","=Num.Coefs")
            print(print.out)
        }
    }

    if(!is.null(object$tests) & detail == TRUE){
        cat("\n---\n") 
        cat("Diagnostic tests for IV regression:\n")
        cat("\n")
        if(!is.null(object$tests$first.stage)){
            if((is.list(object$tests$first.stage) & length(object$tests$first.stage) != 0)){
                cat("--First Stage Regression--\n")
                print.out <- object$tests$first.stage$AR.F.value
                rownames(print.out) <- "Anderson-Rubin F Test"
                colnames(print.out) <- c("F-Value","DF1","DF2","P-Value")
                print(print.out)
                cat("\n")
                en.var.x <- object$names$en.var.x
                for(sub.en.var in en.var.x){
                    if(!is.null(object$tests$first.stage[[sub.en.var]]) ){
                        cat(paste0(sub.en.var,':\n'))
                        print.out <- rbind(object$tests$first.stage[[sub.en.var]][["partial.F"]],
                                        object$tests$first.stage[[sub.en.var]][["AP.F.statistic"]],
                                        object$tests$first.stage[[sub.en.var]][["SW.F.statistic"]])
                        rownames(print.out) <- c("Partial F Statistics","Angrist&Pischke F Statistics","Sanderson&Windmeijer F Statistics")
                        colnames(print.out) <- c("F-Value","df1","df2","P-Value")
                        print(print.out)
                        cat("\n")
                    }
                }
                cat("\n")
            }
        }

        if(!is.null(object$tests$over.identify) & length(object$tests$over.identify)>0){
            cat("--Overidentification Tests--\n")
            if(object$variance.type!='Standard'){
                print.out <- matrix(object$tests$over.identify$hansen,nrow=1)
                rownames(print.out) <- "Hansen Tests:"
                colnames(print.out) <- c("Hansen Statistics","df","P-Value")
                print(print.out) 
            }else{
                print.out <- matrix(object$tests$over.identify$sargan,nrow=1)
                rownames(print.out) <- "Sargan Tests:"
                colnames(print.out) <- c("Sargan Statistics","df","P-Value")
                print(print.out)
            }
            
            if(!is.null(object$tests$over.identify$C.statistic)){
                cat("\n")
                cat(paste0("Overidentification tests for ",paste(object$orthog,collapse=","),"\n"))
                print.out <- matrix(object$tests$over.identify$C.statistic,nrow=1)
                rownames(print.out) <- "C Tests:"
                colnames(print.out) <- c("C Statistics","df","P-Value")
                print(print.out)
            }
            cat("\n")
        }

        if(!is.null(object$tests$endogenous)){
            if((is.list(object$tests$endogenous) & length(object$tests$endogenous) != 0)){
                cat("--Endogeneity Tests--\n")
                if(is.null(object$endog)){
                    cat("Regressors Tested: All endogenous regressors.")
                }else{
                    cat(paste0("Regressors Tested: ",paste(object$endog,collapse=',')))
                }
                cat("\n")
                print.out <- matrix(object$tests$endogenous$DWH,nrow=1)
                colnames(print.out) <- c("DWH Value","df","P-Value")
                rownames(print.out) <- "Endogeneity Test"
                print(print.out)
            } 
        }
        cat("\n")
    }
}

ERR.FE.predict.invalid.inds <- function() {
    sprintf("The given indicators do not belong to the indicators in the model.")
}
