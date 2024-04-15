# legacy-cdd

#Introduction This repository contains the code for the analysis of the project: "The unexpected influence of legacy conspecific density dependence." The project uses data from the [Wabikon ForestGEO plot](https://www.forestgeo.si.edu/sites/north-america/wabikon).

## File structure

Code contained in this repository is available in the source folder. The steps for reproducing the analysis are outlined below.

The full census data cannot be made freely available, but we provide a subset data to be used in repeating the analysis.

# Methods (project steps)

1.  The first step to replicate these methods is to run the neighborhood code on your data

2.  Once the neigbhrhoods are clacluated, the run the Stan model

3.  Diagnostics are included

4.  The post-model analysis for predicted change in survival and the resulting figures are shown in:

# Systems and R packages (dependencies)

All analysis were conducted in R version 4.3.1

The steps used in the analyses relied on the following R packages:

-   
