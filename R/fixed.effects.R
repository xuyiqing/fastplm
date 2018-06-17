locate.effect <- function(inds, id) {
  if (mode(id) == "numeric") {
    if (!is.null(inds[id]))
      return(as.integer(id))
  }
  else if (mode(id) == "character") {
    effect.names <- inds$effect.names
    for (i in SEQ(1, length(effect.names)))
      if (effect.names[i] == id)
        return(i)
  }

  stop(ERR.locate.effect(id))
}

ERR.FE.locate.effect <- function(id) {
  if (mode(id) == 'numeric')
    return(sprintf("Specified effect index %s out of bound", paste(id)))

  if (mode(id) == 'character')
    return(sprintf("Specified effect name %s is not found.", paste(id)))

  sprintf("Effect cannot be specified using id of mode %s.", paste(mode(id)))
}

create.complex.effect <- function(inds, eff, inf, weight) {
  CHECK.INPUT(inds, "inds", "indicators")
  CHECK.INPUT(weight, "weight", "matrix")

  eff.id <- locate.effect(inds, eff)
  inf.id <- locate.effect(inds, inf)

  ASSERT.MATRIX.DIM(weight, "weight", inds$level.sizes[inf.id], is.width = TRUE)

  r <- list(eff = eff.id, inf = inf.id, weight = weight)
  class(r) <- "complex.effect"
  r
}

create.fixed.effects <- function(inds, sfes = NULL, cfes = NULL) {
  CHECK.INPUT(inds, "inds", "indicators")
  CHECK.INPUT(cfes, "cfes", "list of cfes",
    check.type = function(cfes) {
    for (cfe in cfes)
      if (class(cfe) != "complex.effect")
        return(FALSE)
    return(TRUE)
  })

  sfes <-
    if (is.null(sfes))
      SEQ(1, length(inds$levels))
    else
      sapply(sfes, function(id) locate.effect(inds, id))

  cfe.effs <- sapply(cfes, function(x) x$eff)
  cfe.infs <- sapply(cfes, function(x) x$inf)
  weights  <- lapply(cfes, function(x) x$weight)

  # If "cfes" is empty, we need to force those variables' "types".
  if (length(cfe.effs) == 0) {
    cfe.effs <- vector("integer")
    cfe.infs <- vector("integer")
    weights  <- list()
  }

  ptr <- CreateFixedEffects(inds$level.sizes, inds$inds, sfes, cfe.effs, cfe.infs, weights)

  if (ContainMultipleComponents(ptr))
    warning(WARN.FE.multiple.components())

  sfe.names <- sapply(sfes, function(id) inds$effect.names[id])
  cfe.names <- sapply(SEQ(1, length(cfe.effs)), function(i) {
    paste(inds$effect.names[cfe.effs[i]],
          inds$effect.names[cfe.infs[i]],
          sep = ":")
  })

  fe <- mget(c("ptr", "sfes", "sfe.names", "cfe.effs", "cfe.infs", "weights", "cfe.names"), environment())
  class(fe) <- "fixed.effects"
  fe
}

WARN.FE.multiple.components <- function() {
  sprintf("The fixed effect system contains multiple disconnected components. Prediction of entries containing level indicators in different components is meaningless.")
}
