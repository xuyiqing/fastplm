CHECK.INIT.INPUT <- function(x, arg.name, init = NULL, arg.class = NULL, check.type = NULL) {
  if (is.null(x) && !is.null(init))
    x <- init()

  if (is.null(x))
    stop(ERR.GEN.arg.empty(arg.name))

  if (!is.null(arg.class))
    CHECK.INPUT(x, arg.name, arg.class, check.type)

  x
}

ERR.GEN.arg.empty <- function(arg.name) {
  sprintf("Argument \"%s\" is not required but not provided (= NULL).", arg.name)
}

CHECK.INPUT <- function(x, arg.name, arg.class, check.type = NULL) {
  if (!is.null(check.type)) {
    if (check.type(x))
      return()

    stop(ERR.GEN.arg.invalid(arg.name, arg.class))
  }

  if (class(x)[1] != arg.class)
    stop(ERR.GEN.arg.invalid(arg.name, arg.class, class(x)))
}

ERR.GEN.arg.invalid <- function(arg.name, expected, actual = "unknown") {
  sprintf("Argument \"%s\" has invalid class: expected \"%s\", actual \"%s\".", arg.name, expected, actual)
}

ASSERT.MATRIX.DIM <- function(mat, name, expected, is.width = TRUE) {
  if (is.width) {
    if (expected != ncol(mat))
      stop(ERR.GEN.incompatible.matrix(name, expected, ncol(mat)))
  } else {
    if (expected != nrow(mat))
      stop(ERR.GEN.incompatible.matrix(name, expected, ncol(mat), FALSE))
  }
}

ERR.GEN.incompatible.matrix <- function(name, expected, actual, is.width = TRUE) {
  dim <- if (is.width) "width" else "height"

  sprintf("Matrix %s has incompatible %s: expected %d, actual %d.", name, dim, expected, actual)
}

valid.or.default <- function(x, def) {
  if (is.null(x) || is.na(x)) def else x
}

SEQ <- function(from, to) {
  if (from <= to)
    from : to
  else
    c()
}
