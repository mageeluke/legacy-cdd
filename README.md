# legacy-cdd

#Introduction This repository contains the code for the analysis of the project: "The unexpected influence of legacy conspecific density dependence." The project uses data from the [Wabikon ForestGEO plot](https://www.forestgeo.si.edu/sites/north-america/wabikon), which are avialable upon request through the [Data Portal](https://www.forestgeo.si.edu/explore-data).

## File structure

Code contained in this repository is available in the source folder. The steps for reproducing the analysis are outlined below.

The full census data cannot be made freely available, but we provide a subset of 2000 stems to be used in the scripts below.

# Methods (project steps)

Working within the source code:

1.  The first step to replicate our methods is to run the [neighborhood code](https://github.com/mageeluke/legacy-cdd/blob/main/Source/Wabikon_Legacy_CDD_prep_all_species.R) on your data. We provide a subset of the Wabikon plot data for demonstration.

2.  Once the neighborhoods are calculated, save your data and run the Stan model (run_stan_cdd.R)

3.  Diagnostics are included in the Stan model script as the code to calculate average predictive comparisons for each species (Gelman & Pardoe 2007, Gustafson 2007).

4.  The code to plot the model figures is also included in the source folder. We include Stan model parameter outputs as well as the average predictive comparison code.

Additional data:

-   Source code and output are shown for grid searches for 6 and 12 neigbhorhood exponent models './source/grid searches', and "output" folders, respectively.

# Systems and R packages

All analysis were conducted in R version 4.3.1

The analyses relied on the following R packages:

-   rstan (Version 2.26.23)

-   rstanarm (Version 2.26.1)

-   lme4 (Version 1.1-34)

-   DHARMa (Version 0.4.6)

-   boot (Version 1.3-28.1)

-   purrr (Version 1.0.2)

-   data.table (Version 1.14.8)

-   sp (Version 2.1-1)
