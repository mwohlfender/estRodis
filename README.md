# R package estRodis
R package for estimating the effective reproduction number, the dispersion parameter and the testing probability from the sequence cluster size distribution. The theoretical background and a use case of the package's functionalities are presented in (ADD LINK TO PAPER).

## (A) Overview of content of repository

At the heart of the R-package estRodis are functions to simulate the size distribution of identical sequence clusters 

## (B) General remarks

### Name convention
Whenever "model one" is mentioned in comments in the code, this refers to the standard model developed in the paper and "model two" refers to the alternative model described in the section "Sensitivity analysis" of the supplementary material.

## (C) How to install the R package estRodis

R code to install the latest stable release of the estRodis package: \
`devtools::install_github("mwohlfender/estRodis", ref = "main", force = TRUE)`

R code to install the newest development version of the estRodis package: \
`devtools::install_github("mwohlfender/estRodis@v0.0.1-zeta", ref = "main", force = TRUE)`

## (D) How to use the R package estRodis


## (E) Short description of the functionality of the R package estRodis

The R package estRodis contains nine functions. Short introductions of all functions are provided below, more detailed information on input and output values can be found in the documentation of the respective function within the R package.

### E.1 simulate cluster sizes
Simulate the size of identical sequence clusters.

### E.2 simulate transmission chain
Simulate a transmission chain.

### E.3 plot transmission chain
Create a plot of a transmission chain.

### E.4 plot transmission chain with mutations
Create a plot of a transmission chain including mutations of the viral genome.

### E.5 estimate parameters model one
Apply model one to estimate parameters from the sequence cluster size distribution.

### E.6 estimate parameters model two
Apply model two to estimate parameters from the sequence cluster size distribution.

### E.7 initialize parameters model one
Determine initial values for model one.

### E.8 initialize parameters model two
Determine initial values for model two.

### E.9 scaled beta distribution
Probability density function of a scaled version of the beta distribution. (ADD FORMULA) 
This distribution is used as prior distribution of the testing probability in model one.






