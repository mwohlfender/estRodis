# R package estRodis
estRodis is an R package for statistical analyses in the field of infectious disease dynamics. More precisely, for estimating the effective reproduction number, the dispersion parameter and the testing probability from the sequence cluster size distribution. The theoretical background and a use case of the package's functionalities are presented in (ADD LINK TO PAPER).

## (A) Overview of content of R package estRodis

At the heart of the R-package estRodis are functions to simulate the size distribution of identical sequence clusters and to estimate parameters related to transmission dynamics from the sequence cluster size distribution. The simulation is based on a mathematical model of the size distribution of identical sequence clusters that takes into account the viral transmission process, the mutation of the virus and incomplete case-detection. We described the probability that a case was detected as the product of the probability that a case was confirmed by a test and the probability that the viral genome of a confirmed case is sequenced. There are two Bayesian inference models to estimate the effective reproduction number, the dispersion parameter, the yearly mutation rate and the testing probability to choose from. The first model uses a weakly informative prior distribution for the probability that a case is confirmed by a test, whereas the second model takes a fixed value for this probability as input. Both models use weakly informative prior distributions for the effective reproduction number, the dispersion parameter and the yearly mutation rate and require a fixed value for the probability that the viral genome of a confirmed case is sequenced as input. 

## (B) General remarks

### B.1 Use case
The code used for the statistical analysis of the paper (ADD LINK TO PAPER), which is a use case of the estRodis package, can be found here: (ADD LINK TO GITHUB REPO)

### B.2 Name convention
Whenever "model one" is mentioned in comments in the code, this refers to the standard model developed in the paper (prior distribution for testing probability) and "model two" refers to the alternative model described in the section "Sensitivity analysis" of the supplementary material (fixed value for testing probability).

## (C) How to install the R package estRodis

R code to install the latest stable release of the estRodis package: \
`devtools::install_github("mwohlfender/estRodis", ref = "main", force = TRUE)`

R code to install the newest development version of the estRodis package: \
`devtools::install_github("mwohlfender/estRodis@v0.0.1-zeta", ref = "main", force = TRUE)`

## (D) How to use the R package estRodis

### (D1) Simulate data and apply model one

We simulate 1000 identical sequence clusters. 
```
simulated_clusters <- estRodis_simulate_cluster_sizes(n_clusters = 1000,
                                                      max_cluster_size = 2500,
                                                      R = 0.8,
                                                      k = 0.3,
                                                      yearly_mutation_rate = 14,
                                                      mean_generation_interval = 5.2,
                                                      testing_proba = 0.6,
                                                      sequencing_proba = 0.4)
```
Result:

|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|__size__| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 19 | 20 | 22 | 24 | 26 |
|__frequency__| 725 | 111 | 59 | 31 | 16 | 13| 10 | 3 | 6 | 2 | 3 | 4 | 4 | 2 | 4 | 1 | 2 | 1 | 1 | 1 | 1 |
|__percentage__| 0.725 | 0.111 | 0.059 | 0.031 | 0.016 | 0.013 | 0.01 | 0.003 | 0.006 | 0.002 | 0.003 | 0.004 | 0.004 | 0.002 | 0.004 | 0.001 | 0.002 | 0.001 | 0.001 | 0.001 | 0.001|

Next, we apply model one to estimate the effective reproduction number, the dispersion parameter, the yearly mutation rate and the testing probability from the simulated clusters.
```
options(mc.cores = parallelly::availableCores())

estRodis_estimate_parameters_one(
  clusters_size = simulated_clusters$size,
  clusters_freq = simulated_clusters$frequency,
  sequencing_proba = 0.44)
```
Result:
```
Inference for Stan model: estRodis_stan_model_estimate_parameters_one.
4 chains, each with iter=2000; warmup=1000; thin=1; 
post-warmup draws per chain=1000, total post-warmup draws=4000.

                            mean se_mean   sd     2.5%      25%      50%      75%    97.5% n_eff Rhat
R                           0.81    0.00 0.06     0.71     0.76     0.80     0.85     0.95  1171    1
k                           0.34    0.00 0.08     0.19     0.28     0.33     0.39     0.53  1813    1
number_yearly_mutations    13.99    0.01 0.49    13.00    13.66    13.99    14.31    14.95  2254    1
testing_proba               0.60    0.01 0.24     0.15     0.41     0.62     0.81     0.99  1061    1
mutation_proba              0.18    0.00 0.01     0.17     0.18     0.18     0.18     0.19  2252    1
detection_proba             0.27    0.00 0.11     0.07     0.18     0.27     0.36     0.43  1061    1
lp__                    -1156.67    0.04 1.54 -1160.57 -1157.41 -1156.29 -1155.53 -1154.84  1200    1

Samples were drawn using NUTS(diag_e) at Tue Mar 26 11:01:36 2024.
For each parameter, n_eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor on split chains (at 
convergence, Rhat=1).
```
### (D2) Simulate data and apply model two

```
simulated_clusters <- estRodis_simulate_cluster_sizes(n_clusters = 1000,
max_cluster_size = 2500, R = 0.8, k = 0.3, yearly_mutation_rate = 14,
mean_generation_interval = 5.2, testing_proba = 0.6, sequencing_proba = 0.4)

options(mc.cores = parallelly::availableCores())

estRodis_estimate_parameters_two(
  clusters_size = simulated_clusters$size,
  clusters_freq = simulated_clusters$frequency,
  testing_proba = 0.55,
  sequencing_proba = 0.44)
```

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






