
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
Tudball, M., Hughes, R., Tilling, K., Bowden, J., Zhao, Q. Sample-constrained partial identification with application to selection bias. *Biometrika*. https://academic.oup.com/biomet/advance-article/doi/10.1093/biomet/asac042/6649721 (2022).
