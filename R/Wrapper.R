solveFE <- function(rawData, rawFixedEffects, coreNum = 1, estimateFE = FALSE) {
  solve.fixed.effects(rawData, rawFixedEffects, coreNum)
}

name.fxied.effects <- function(model, inds) {
  names(model$FEcoefs) <-
    if (is.null(colnames(inds)))
      sapply(1 : ncol(inds), function(col) paste("effect", col, sep = "."))
    else
      colnames(inds)

  for (col in 1 : ncol(inds))
    rownames(model$FEcoefs[[col]]) <- model$group.levels[[col]]

  model
}

create.indicators <- function (inds) {
  factors      <- lapply(1 : ncol(inds), function(col) factor(inds[, col]))
  levels       <- lapply(factors, levels)
  level.sizes  <- sapply(levels, length)
  effect.names <-
    if (is.null(colnames(inds)))
      sapply(1 : ncol(inds), function(col) paste("effect", col, sep = "."))
    else
      colnames(inds)

  inds         <- sapply(factors, as.numeric)

  r <- mget(c("levels", "level.sizes", "effect.names", "inds"),
    environment())
  class(r) <- "indicators"
  r
}

locate.effect <- function(inds, id) {
  if (mode(id) == "numeric") {
    if (!is.null(inds[id]))
      return(as.integer(id))
  }
  else if (mode(id) == "character") {
    effect.names <- inds$effect.names
    for (i in 1 : length(effect.names))
      if (effect.names[i] == id)
        return(i)
  }

  stop(error.locate.effect(id))
}

create.complex.effect <- function(inds, eff, inf, weight) {
  if (class(inds) != "indicators")
    stop(error.indicators.wrong.type())

  eff.id <- locate.effect(inds, eff)
  inf.id <- locate.effect(inds, inf)

  if (!is.matrix(weight))
    stop(error.create.complex.effect.matrix.required())

  exp.ncol <- inds$level.sizes[inf.id]
  act.ncol <- ncol(weight)
  if (exp.ncol != act.ncol)
    stop(error.create.complex.effect.matrix.wrong.size(exp.ncol, act.ncol))

  r <- list(eff = eff.id, inf = inf.id, weight = weight)
  class(r) <- "complex.effect"
  r
}

create.fixed.effects <- function(inds, sfes = NULL, cfes = NULL) {
  if (class(inds) != "indicators")
    stop(error.indicators.wrong.type())

  sfes <-
    if (is.null(sfes))
      1 : length(inds$levels)
    else
      sapply(sfes, locate.effect)

  for (cfe in cfes)
    if (class(cfe) != "complex.effect")
      stop(error.cfes.wrong.type())

  cfe.effs <- vector("integer")
  cfe.infs <- vector("integer")
  weights  <- list()

  if (!is.null(cfes)) {
    cfe.effs <- sapply(cfes, function(x) x$eff)
    cfe.infs <- sapply(cfes, function(x) x$inf)
    weights  <- lapply(cfes, function(x) x$weight)
  }

  cpp.fixed.effects <- CreateFixedEffects(inds$level.sizes, inds$inds,
    sfes, cfe.effs, cfe.infs, weights)
  if (ContainMultipleComponents(cpp.fixed.effects))
    warning(warning.contain.multiple.components())
  cpp.fixed.effects
}

solve.fixed.effects <- function(data, inds, core.num = 1) {
  indicators        <- create.indicators(inds)
  cpp.fixed.effects <- create.fixed.effects(indicators)

  model <- SolveFixedEffects(data, cpp.fixed.effects, core.num)
  model$cpp.fixed.effects <- cpp.fixed.effects

  model$group.levels <- indicators$levels
  model <- name.fxied.effects(model, inds)
  model
}

predictFE <- function(model, newX, FEValues = NULL, grandMean = 0) {
  predict.fixed.effects(model, newX, FEValues)
}

predict.fixed.effects <- function(model, x, inds = NULL) {
  if (is.null(x) && is.null(inds))
    stop("x and inds cannot both be NULL.")

  if (is.null(x))
    x <- matrix(data = NA, nrow = ncol(inds), ncol = 0)

  y <- x %*% model$coefficients + model$intercept

  if (is.null(inds))
    return(y)

  factored.inds <- sapply(1 : ncol(inds),
    function(col) as.numeric(factor(inds[, col], levels = model$group.levels[[col]])))

  for (i in which(is.na(factored.inds))) {
    col <- (i - 1) %/% nrow(inds) + 1
    row <- i - (col - 1) * nrow(inds)
    warning(warning.predict.fixed.effects.with.unknown.indicators(row, col, inds))
  }

  cross.components <- CheckComponents(model$cpp.fixed.effects, factored.inds)
  rows.to.knock.out <- sapply(cross.components, function(error) {
    warning(warning.cross.component(error, inds, model$group.levels))
    error$row
  })
  if (length(rows.to.knock.out) > 0)
    factored.inds[rows.to.knock.out, ] <- NA

  with.effects <- function(row) sum(sapply(1 : ncol(inds),
    function(col) (model$FEcoefs[[col]])[factored.inds[row, col]]))

  y <- y + sapply(1 : nrow(inds), with.effects)

  return(y)
}
