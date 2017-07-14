# FastPLM

A package for fast fixed-effects algorithms.

## Package Structures

The package is structurally similar to most packages using Rcpp, with `src` directories containing all C/C++ codes. However, we also provide `cmake` support, with which we can build and debug a pure C++ package using whatever tool chains convenient.

## Usage

There are currently two functions exposed directly using Rcpp:

* `solveFE`
* `solveGP`



### SolveFE

This function solves an ordianary linear system with multi-way fixed effects. It currently takes two inputs:

* a normal data matrix with its first column being Y and the rest being X;
* a matrix specifying fixed effects; each column represents a category, and within each column it specifies which group of that category the row belongs to.

Concretely, suppose we have three individuals A, B, C and two firms X, Y. We make six observations sorted as (A, X), (B, X), (C, X), (A, Y), (B, Y), (C, Y). The fixed effect matrix would look like

| Individual | Firm |
| ---------- | ---- |
| 1          | 1    |
| 2          | 1    |
| 3          | 1    |
| 1          | 2    |
| 2          | 2    |
| 3          | 2    |

The number for group does not need to start from 1. The programs internally handle all complexities to identify different groups.

The output would including:

* coefficients or parameters or β or whatever you prefer to call it;
* intercept;
* fixed effects for each category.

Coefficients are estimated via MAP, whereas fixed effects are estimated using gradient descent.

Note that fixed effects for each category is a vector. It is identified with the given group values in an increasing order. So let us say we have three individuals labeld as A = 7, B = 2, C = 9, then the output fixed effect for the individual category would be in the order of B, A, C.

### SolveGP

This function handles the model as specified in `math/GPanel.pdf`. It only supports two way effects. For simplicity, let’s assume the two categories are time and individuals. 

The function takes five arguments:

* Y,
* X,
* time-observable individual-specific (TOIS) effects,
* individual-observable time-specific (IOTS) effects,
* a boolean value specifying whether the data is balanced (`TRUE`) or not.

A few notes on data structure:

* Logically, we assume Y to be a matrix with rows specifying time and columns specifying individuals;
* X is treated as a cube, with the third dimension meaning different variables;
* Concretely, we treat Y as a vector and X a matrix. The data is sorted in a column major. In other words, the vector would be the same as we go through the first column of Y, then the second, then the third, and so on.
* TOIS and IOTS are organized the same as specified in  `math/GPanel.pdf`.
* For unbalanced situations, knocking out missing cells in Y by marking it as `NaN` ((not a number)[https://en.wikipedia.org/wiki/NaN]).

The output is the estimated coefficients.

This function is just finished. We will polish it to add more options and outputs.

