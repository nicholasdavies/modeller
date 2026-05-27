
# modeller: interactive infectious disease modelling

<!-- badges: start -->
[![R-CMD-check](https://github.com/nicholasdavies/modeller/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nicholasdavies/modeller/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Build, run, and interactively explore compartmental infectious
disease models. Supports deterministic (ODE), stochastic (SSA), and
discrete-time (difference equation) solvers with a simple model
specification syntax. Includes an interactive Shiny app
for adjusting parameters and initial conditions, visualising dynamics, and
overlaying data. Developed for the Modelling and the Dynamics of Infectious
Diseases course at the London School of Hygiene & Tropical Medicine.

## Installation

You can install the development version of modeller like so:

``` r
install.packages("remotes")
remotes::install_github("nicholasdavies/modeller")
```
