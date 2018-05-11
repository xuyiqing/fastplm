# FastPLM

A package for fast fixed-effects algorithms.

## Package Structures

The package is structurally similar to most packages using Rcpp, with `src` directories containing all C/C++ codes. However, we also provide `cmake` support, with which we can build and debug a pure C++ package using whatever tool chains convenient.

## Primitives

The library provides a list of primitives to solve fixed effects problems:

* `create.indicators`
* `create.complex.effect`
* `create.fixed.effects`
* `solve.fixed.effects`

We will also expose some higher level functions simplifying some everyday usage.
Those functions will be merely a thin wrapper of those primitives,
and can be implemented in pure R without performance hits.

The following section documents those primitives. Undocumented features,
including but not limited optional parameters and properties of returned objects
should not be tampered with and are subject to change at any time.
Use them at your own risks.

- `create.indicators`:

  Creates an object describing the groups each observation belongs to.
  The result is required in other primitives.

  - parameters:

    - `inds`: A matrix where each row corresponds to the row (observation)
    in the X-Y dataset and each column represents an effect.
    The entry of row X of effect Y represents the group in Y that X belongs to.
    Groups can be specified by numbers or strings.

  - returns: An object of class `indicators` with the following properties:

    - `levels`: Distinct groups in each effect.
    It is structured as a list of the same length as the number of effects.
    Each list item is a vector of levels (distinct groups) of the corresponding
    effect. The order of levels in each list item is not guaranteed.

    - `level.sizes`: The number of occurrences of each level in the input.
    The structure and order of `level.sizes` is the same as `levels`.

    * `effect.names`: A vector of names of each effect.
    If `colnames` of the input is set, that will be used.
    Otherwise, effects will be named as `effect.1`, `effect.2`, and so on,
    in the same order as they are listed as the columns of the input.

- `create.complex.effect`

  Creates a *complex fixed effect*.
  The result is used to construct the system of fixed effects.

  A *complex fixed effect* is a generalized version of fixed effect.
  It consists a pair of two effects (I, E) interacting with each other.
  I stands for *influence*. Each level of I has an *observed* vector `weight`.
  E stands for *effect*. Each level of Ehas an *unobserved* vector `coefficient`
  to estimate.
  For each observation, the effect of some complex fixed effect (I, E)
  is the dot product of the weight of the row’s group for I
  and the coefficient of the row’s level for E.

  - parameters:

    - `inds`: The `indicators` object the complex fixed effect is about.

    - `eff`: The effect of the complex fixed effect.
    It can be specified either by effect name as in `inds$effect.names`,
    or by the index of the effect in `inds`.

    - `inf`: The influence of the complex fixed effect.
    It can be specified either by effect name as in `inds$effect.names`,
    or by the index of the effect in `inds`.

    - `weight`: A matrix of weights associated with the influence.
    Each column represents the weight of the corresponding level.
    Weight should be given in the same order as the order of levels
    for influence in `inds$levels`.

      Only matrix is accepted, even if the weight associated with each
      level is one-dimensional (i.e. a scalar). Use `as.matrix` when necessary.

  - returns: An object of class `complex.effect`.

- `create.fixed.effects`

  Creates a system of fixed effects for some model.
  The result is used to solve that model.

  - parameters:

    - `inds`: The `indicators` object the system of fixed effects is about.

    - `sfes`: A collection of simple (i.e. non-complex) fixed effects.
    Each can be specified either by effect name as in `inds$effect.names`,
    or by the index of the effect in `inds`.

      `sfes` is optional. If not specified, it will be taken to be the
      collection of all effects in `inds`.

    - `cfes`: A collection of complex fixed effects. Each should be an object
    of class `complex.effect` created by `create.complex.effect`.

      `cfes` is optional. If not specified, it will be taken as an empty list.

  - returns: An opaque object representing the system of fixed effects.

- `solve.fixed.effects`

  Solves the fixed effects model.

  - warning: The function has not been implemented.

  - parameters:

    - `data`: The X-Y dataset. The first column is the Y vector.
    The rest is the X matrix.

      `data` is optional. If not specified, then `data` will be built from
      `x` and `y` (i.e. `data <- cbind(y, x)`).
      If both `data` and `x`, `y` are present, `data` will be used used.

    - `y`: The Y vector. Optional if `data` is provided.

    - `x`: The X matrix. Optional if `data` is provided.

      X can have zero column; in such case, the function merely estimates
      these fixed effects.

    - `inds`: The `indicators` object the model is about.

    - `fe`: The system of fixed effects characterizing the model
    created by `create.fixed.effect`.

    - `core.num`: The number of cores used in demean.

      `core.num` is optional. Default is `1`.

  - returns: The object of class `fe.model` recording the solution.

    TODO: finish the specification for `fe.model`.
