
# selectioninterval: An R package for a selection bias sensitivity analysis

<!-- badges: start -->
<!-- badges: end -->

The goal of selectioninterval is to provide a confidence interval for a range of 
possible estimates, given a logistic model for the inverse probability weights. 
Currently, ordinary least squares and two-stage least squares are the only supported
estimators.

## Installation

To install, run

``` r
library(devtools)
install_github("matt-tudball/selectioninterval")
```

## Example

This is a basic example of a selection bias sensitivity analysis:

``` r
library(selectioninterval)
## To be done
```

## References
Tudball, M., Zhao, Q., Hughes, R., Tilling, K., & Bowden, J. An interval estimation approach to sample selection bias. *arXiv*. http://arxiv.org/abs/1906.10159 (2020).
