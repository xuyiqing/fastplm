solve.fe.model <- function(inds, fe = NULL, data = NULL, x = NULL, y = NULL, core.num = 1) {
  CHECK.INPUT(inds, "inds", "indicators")

  data <- CHECK.INIT.INPUT(data, "data, x, y",
    init = function() cbind(y, x))
  fe  <- CHECK.INIT.INPUT(fe, "fe", arg.class = "fixed.effects",
    init = function() create.fixed.effects(inds))

  model <- SolveFixedEffects(data, fe$ptr, core.num)
  model <- name.fe.model(model, inds, fe)
  model$inds <- inds
  model$fe   <- fe
  class(model) <- "fe.model"
  model
}

name.fe.model <- function(model, inds, fe) {
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

predict.fe.model <- function(model, x, inds) {
  CHECK.INPUT(model, "model", "fe.model")
  CHECK.INPUT(x, "x", "matrix")
  CHECK.INPUT(inds, "inds", "sub.indicators")

  ASSERT.MATRIX.DIM(x, "x", length(model$coefficients), is.width = TRUE)

  if (inds$parent != model$inds$uid)
    stop(ERR.FE.predict.invalid.inds())

  y  <- x %*% model$coefficients + model$intercept
  fe <- model$fe

  for (col in fe$sfes) {
    effs <- model$sfe.coefs[[col]]
    sum  <- sapply(SEQ(1, nrow(x)), function(i) effs[inds$inds[i, col]])
    y    <- y + sum
  }

  for (i in SEQ(1, length(fe$cfe.effs))) {
    eff.col <- fe$cfe.effs[i]
    inf.col <- fe$cfe.infs[i]
    effs    <- model$cfe.coefs[[i]]
    weights <- fe$weights[[i]]

    sum <- sapply(SEQ(1, nrow(x)), function(i)
      effs[inds$inds[i, eff.col], ] %*% weights[, inds$inds[i, inf.col]])
    y   <- y + sum
  }

  y
}

ERR.FE.predict.invalid.inds <- function() {
  sprintf("The given indicators do not belong to the indicators in the model.")
}
