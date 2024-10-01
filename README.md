# R package estRodis
estRodis is an R package for statistical analyses in the field of infectious disease dynamics. More precisely, for estimating the effective reproduction number, the dispersion parameter and the testing probability from the sequence cluster size distribution. The theoretical background and a use case of the package's functionalities are presented in the paper "Estimating $R_e$ and overdispersion in secondary cases from the size of identical sequence clusters of SARS-CoV-2" by Emma Hodcroft et al. A preprint is available on [medRxiv](https://www.medrxiv.org/content/10.1101/2024.05.26.24307940v1).

## (A) Overview of content of R package estRodis

At the heart of the R-package estRodis are functions to simulate the size distribution of identical sequence clusters and to estimate parameters related to transmission dynamics from the sequence cluster size distribution. The simulation is based on a mathematical model of the size distribution of identical sequence clusters that takes into account the viral transmission process, the mutation of the virus and incomplete case-detection. We described the probability that a case was detected as the product of the probability that a case was confirmed by a test and the probability that the viral genome of a confirmed case is sequenced. There are six Bayesian inference models to estimate the effective reproduction number, the dispersion parameter, the probability a mutation occurs at a case before onward transmission or the yearly mutation rate of the virus and the testing probability to choose from.

### Overview of models

All models use weakly informative prior distributions for the effective reproduction number $R_{e}$ (gamma distributions with parameters (10, 10)) and the dispersion parameter $k$ (gamma distribution with parameters (5,10)). Model five is the main model and is presented in the paper "Estimating $R_e$ and overdispersion in secondary cases from the size of identical sequence clusters of SARS-CoV-2" by Emma Hodcroft et al. (see [medRxiv](https://www.medrxiv.org/content/10.1101/2024.05.26.24307940v1)). The other five models are discussed in the supplementary material.

* `model one` informative prior distribution for the yearly mutation rate (normal distribution with mean 14 and variance 0.25) and weakly informative prior distribution for the probability that a case is confirmed by a test (scaled beta distribution on the interval [0.05, 1] with parameters (1,3), see below for more details on the scaled beta distribution)
* `model two` informative prior distribution for the yearly mutation rate (normal distribution with mean 14 and variance 0.25) and constant value for the probability that a case is confirmed by a test
* `model three` constant value for the probability a mutation occurs at a case before onward transmission and weakly informative prior distribution for the probability that a case is confirmed by a test (scaled beta distribution on the interval [0.05, 1] with parameters (1,3), see below for more details on the scaled beta distribution)
* `model four` constant value for the probability a mutation occurs at a case before onward transmission and constant value for the probability that a case is confirmed by a test
* `model five` informative prior distribution for the probability a mutation occurs at a case before onward transmission (beta distribution with parameters (27,68)) and weakly informative prior distribution for the probability that a case is confirmed by a test (scaled beta distribution on the interval [0.05, 1] with parameters (1,3), see below for more details on the scaled beta distribution)
* `model six` informative prior distribution for the probability a mutation occurs at a case before onward transmission (beta distribution with parameters (27,68)) and constant value for the probability that a case is confirmed by a test

## (B) General remarks

### B.1 Use case
The code used for the statistical analysis of the paper "Estimating $R_e$ and overdispersion in secondary cases from the size of identical sequence clusters of SARS-CoV-2" by Emma Hodcroft et al. (see [medRxiv](https://www.medrxiv.org/content/10.1101/2024.05.26.24307940v1)), which is a use case of the estRodis package, can be found here: [GitHub Martin Wohlfender](https://github.com/mwohlfender/R_overdispersion_cluster_size)

### B.2 Name convention
Whenever "model five" is mentioned in comments in the code, this refers to the standard model developed in the paper. model one", "model two", "model three", "model four" and "model six" refer to the alternative models described in the section "Sensitivity analysis" of the supplementary material of the paper  "Estimating $R_e$ and overdispersion in secondary cases from the size of identical sequence clusters of SARS-CoV-2" by Emma Hodcroft et al.

## (C) How to install the R package estRodis

R code to install the latest stable release of the estRodis package: \
```devtools::install_github("mwohlfender/estRodis@v1.0.0", ref = "main", force = TRUE)```

R code to install the newest development version of the estRodis package: \
```devtools::install_github("mwohlfender/estRodis", ref = "main", force = TRUE)```

## (D) How to use the R package estRodis

In the following we present code examples how to use the main functionalities of the R package estRodis.
### (D1) Simulate data

We simulate 1000 identical sequence clusters. 
```
simulated_clusters <- estRodis_simulate_cluster_sizes_v2(n_clusters = 1000,
                                                         max_cluster_size = 2500,
                                                         R = 0.8,
                                                         k = 0.3,
                                                         mutation_proba = 0.2,
                                                         testing_proba = 0.6,
                                                         sequencing_proba = 0.4)
```
Result:

|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|__size__| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 13 | 14 | 15 | 18 | 19 | 21 | 25 | 28 | 31 |
|__frequency__| 742 | 104 | 49 | 31 | 22 | 13 | 9 | 6 | 2 | 4 | 2 | 6 | 1 | 2 | 2 | 1 | 1 | 1 | 1 | 1 |
|__percentage__| 0.742 | 0.104 | 0.049 | 0.031 | 0.022 | 0.013 | 0.009 | 0.006 | 0.002 | 0.004 | 0.002 | 0.006 | 0.001 | 0.002 | 0.002 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 |

### (D2) Apply model five

We apply model five to estimate the effective reproduction number, the dispersion parameter, the mutation probability and the testing probability from the clusters we simulated in section (D1).
```
options(mc.cores = parallelly::availableCores())

estRodis_estimate_parameters_five(clusters_size = simulated_clusters$size,
                                  clusters_freq = simulated_clusters$frequency,
                                  sequencing_proba = 0.44)
```
Result:
```
Inference for Stan model: estRodis_stan_model_estimate_parameters_five.
4 chains, each with iter=2000; warmup=1000; thin=1; 
post-warmup draws per chain=1000, total post-warmup draws=4000.

                    mean se_mean   sd     2.5%      25%      50%      75%    97.5% n_eff Rhat
R                   0.89    0.00 0.09     0.74     0.82     0.88     0.94     1.08  1144    1
k                   0.29    0.00 0.07     0.16     0.24     0.28     0.33     0.44  1905    1
mutation_proba      0.29    0.00 0.05     0.20     0.25     0.29     0.32     0.38  1559    1
testing_proba       0.65    0.01 0.23     0.20     0.49     0.68     0.84     0.98  1204    1
detection_proba     0.29    0.00 0.10     0.09     0.21     0.30     0.37     0.43  1204    1
lp__            -1165.95    0.04 1.48 -1169.70 -1166.68 -1165.62 -1164.86 -1164.08  1199    1

Samples were drawn using NUTS(diag_e) at Tue Oct  1 15:42:56 2024.
For each parameter, n_eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor on split chains (at 
convergence, Rhat=1).
```
### (D3) Apply model six

We apply model six to estimate the effective reproduction number, the dispersion parameter, the mutation probability and the testing probability from the clusters we simulated in section (D1).
```
options(mc.cores = parallelly::availableCores())

estRodis_estimate_parameters_six(clusters_size = simulated_clusters$size,
                                 clusters_freq = simulated_clusters$frequency,
                                 testing_proba = 0.55,
                                 sequencing_proba = 0.44)
```
Result:
```
Inference for Stan model: estRodis_stan_model_estimate_parameters_six.
4 chains, each with iter=2000; warmup=1000; thin=1; 
post-warmup draws per chain=1000, total post-warmup draws=4000.

                   mean se_mean   sd     2.5%      25%      50%      75%    97.5% n_eff Rhat
R                  0.90    0.00 0.06     0.79     0.86     0.90     0.94     1.04  1439    1
k                  0.28    0.00 0.06     0.18     0.23     0.27     0.31     0.42  1748    1
mutation_proba     0.29    0.00 0.04     0.20     0.26     0.29     0.32     0.38  1467    1
lp__           -1163.88    0.03 1.19 -1166.78 -1164.44 -1163.58 -1163.00 -1162.51  1345    1

Samples were drawn using NUTS(diag_e) at Tue Oct  1 15:47:19 2024.
For each parameter, n_eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor on split chains (at 
convergence, Rhat=1).
```

## (E) Short description of the functionality of the R package estRodis

The R package estRodis contains 18 functions. Short introductions of all functions are provided below, more detailed information on input and output values can be found in the documentation of the respective function within the R package.

### E.1 simulate cluster sizes
Simulate the size of identical sequence clusters implementing the mutation process via the yearly mutation rate and the mean generation interval.

### E.2 simulate cluster sizes v2
Simulate the size of identical sequence clusters implementing the mutation process directly via the probability a mutation occurs at a case before onward transmission.

### E.3 simulate transmission chain
Simulate a transmission chain.

### E.4 plot transmission chain
Create a plot of a transmission chain.

### E.5 plot transmission chain with mutations
Create a plot of a transmission chain including mutations of the viral genome.

### E.6 estimate parameters model one, two , three, four, five and six
Apply model one, two, three, four, five or six to estimate parameters from the sequence cluster size distribution.

### E.7 initialize parameters model one, two , three, four, five and six
Determine initial values for model one, two, three, four, five or six.

### E.8 scaled beta distribution
Probability density function of a scaled version of the beta distribution. 

$$ f \left( x \right) = \frac{ \left( x - 0.05 \right)^{ \alpha - 1 } \left( 1 - x \right)^{ \beta - 1 }}{ \left( 1 - 0.05 \right)^{ \alpha + \beta - 1 }} \frac{\Gamma \left( \alpha + \beta \right) }{ \Gamma \left( \alpha \right) \Gamma \left( \beta \right)} \quad \forall x  \in \left[ 0.05, 1 \right]  $$

This distribution is used as prior distribution of the testing probability in models one, thee ad five.






