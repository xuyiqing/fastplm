# FastPLM

A package for fast fixed-effects algorithms.

## Package Structures

The package is structurally similar to most packages using Rcpp, with `src` directories containing all C/C++ codes. However, we also provide `cmake` support, with which we can build and debug a pure C++ package using whatever tool chains convenient.

## Primitives

The library provides the following primitives to solve a fixed effect system:

* `create.indicators`
* `create.complex.effect`
* `create.fixed.effects`
* `solve.fe.model`

Once solved, one can use the following two primitives to predict based on a solved system:

* `create.subindicators`
* `predict.fe.model`

We also plan to expose some higher level functions simplifying some everyday usage.
Those functions will be merely a thin wrapper of those primitives,
and can be implemented in pure R without performance hits.
