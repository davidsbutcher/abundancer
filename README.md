Abundancer
================

Abundancer is used to estimate the isotopic abundances of <sup>12</sup>C
and <sup>14</sup>N from a peaklist of a single charge state of an intact
protein and its sequence. This is useful in cases where the abundances
of these isotopes are expected to deviate from typical values, *e.g.* in
cases where proteins are isotopically depleted.

## Installation

Install from GitHub:

``` r
remotes::install_github("davidsbutcher/abundancer")
```

## Usage

To use abundancer, see documentation for
`calculate_score_matrix_dual()`.

Arguments for `calculate_score_matrix_dual()`:

-   `MSspectrum`: Input peaklist as an `MSnbase` Spectrum1 object. If
    not proved, `mz` and `intensity` arguments must be provided.
-   `mz`: Numeric vector of peaklist mz values. Must be same length as
    `intensity`.
-   `intensity`: Numeric vector of peaklist intensity values. Must be
    same length as `mz`.
-   `sequence`: Character vector of single-letter canonical amino acids.
-   `PTMformula` Character vector containing summary difference formula
    for any PTMs on the protein in Hill notation, *e.g.* “C1H2” for
    methylation.
-   `charge`: Charge state og the protein represented by the peaklist.
-   `resolvingPower`: Resolving power at 400 m/z for the spectrum.
    Defaults to 300,000, which is typical for 21 T FT-ICR MS data.

Other functions like `calculate_score_array()` are also available, but
under heavy development and not recommended for use.

A Shiny web application for abundancer is hosted at the MagLab ICR
group’s [TDP apps site](https://tdpapps.magnet.fsu.edu). A [Docker image
of the web
application](https://hub.docker.com/repository/docker/davidsbutcher/abundancer)
is also available.
