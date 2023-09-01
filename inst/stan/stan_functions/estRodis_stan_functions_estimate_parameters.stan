
functions {

  // scaled beta distribution with shape parameters alpha and beta on interval [p,q]
  real estRodis_stan_scaled_beta_lpdf(real alpha, real beta, real p, real q, real x) {

    real result = 0;

    if (x  >= p && x <= q) {

      result = ((x-p)^(alpha-1) * (q - x)^(beta-1)) / ((q - p)^(alpha+beta-1) * tgamma(alpha) * tgamma(beta) / tgamma(alpha + beta));

    }

    return(result);

  }



  int determine_scaling_factor(real x) {

    int result = 0;

    if (x >= 1000) {

      result = 1000;

    } else {

      while (result < x) result += 1;

    }

    return(result);

  }


  // random variable X: size of identical sequence cluster (1, 2, 3, ...)
  // random variable Y: number of detected cases of an identical sequence cluster (0, 1, 2, ...)
  // random variable Z: number of observed cases of an identical sequence cluster (1, 2, ...)
  // as we can not deduce the number of identical sequence clusters of which no case was detected, we need to renormalize the distribution of Y
  // to get the distribution of the size of those identical sequence clusters which we can observe (P(Z=n) = P(Y=n) / P(Y=0))

  real estRodis_stan_likelihood_log(int[] clusters_size, int[] clusters_freq, real R, real k, real mutation_proba, real detection_proba) {

    // here the number of different identical sequence cluster sizes is stored
    int n_clusters_size = num_elements(clusters_size);

    // here the upper limit up to which the probability of observing an identical sequence cluster of a certain size is calculated is stored
    int upper_limit_size_identical_sequence_clusters = max(determine_scaling_factor(2 / detection_proba) * max(clusters_size), max(clusters_size) + 500);

    // here a variable needed for storing results at the right place in a list is stored
    int index = 2;

    // here the logarithm of R * (1 -  mutation_proba) / k is stored
    real log_R_d_k = log(R * (1 - mutation_proba) / k);

    // here the logarithm of 1 + R * (1 - mutation_proba) / k is stored
    real log_one_p_R_d_k = log(1 + R * (1 - mutation_proba) / k);

    // here the logarithm of detection_proba will be stored (if detection_proba < 1)
    real log_detection_proba;

    // here the logarithm of 1 - detection_proba will be stored (if detection_proba < 1)
    real log_1_m_detection_proba;

    // here the logarithm of 1 - P(Y=0)  will be stored
    real log_proba_detect_ident_seq_tree_size_bigger_zero;

    // here the logarithms of 1!, 2!, 3!, ..., upper_limit_size_identical_sequence_clusters! will be stored (if detection_proba < 1)
    real log_sums[upper_limit_size_identical_sequence_clusters] = rep_array(0.0, upper_limit_size_identical_sequence_clusters);

    // here intermediate results needed during the computations will be stored
    real result_temp[upper_limit_size_identical_sequence_clusters] = rep_array(-10000.0, upper_limit_size_identical_sequence_clusters);

    // here the logarithms of P(X=1), P(X=2), ..., P(X=upper_limit_size_identical_sequence_clusters) will be stored
    real distribution_size_ident_seq_tree[upper_limit_size_identical_sequence_clusters] = rep_array(-10000.0, upper_limit_size_identical_sequence_clusters);

    // here the logarithms of P(Y=0) and P(Y=n) for all n contained in clusters_size will be stored
    real distribution_size_ident_seq_tree_detection_0[n_clusters_size + 1] =  rep_array(-10000.0, n_clusters_size + 1);

    // here the logarithms of P(Z=n) for all n contained in clusters_size will be stored
    real distribution_size_ident_seq_tree_detection[n_clusters_size] = rep_array(-10000.0, n_clusters_size);

    // here the likelihood will be stored
    real likelihood[n_clusters_size] = rep_array(-10000.0, n_clusters_size);

    // P(X=1), P(X=2), ..., P(X=upper_limit_size_identical_sequence_clusters) are determined
    for(ii in 1:upper_limit_size_identical_sequence_clusters) {

      // P(X = ii) is determined
      distribution_size_ident_seq_tree[ii] = lgamma(k*ii + ii - 1) - lgamma(k*ii) - lgamma(ii + 1) + (ii - 1) * log_R_d_k - (k*ii + ii - 1)*log_one_p_R_d_k;

    }

    if (detection_proba == 1.0) {

      for (ii in 1:n_clusters_size) distribution_size_ident_seq_tree_detection[ii] = distribution_size_ident_seq_tree[ii];

    } else {

      // here the logarithm of detection_proba is stored
      log_detection_proba = log(detection_proba);

      // here the logarithm of 1 - detection_proba is stored
      log_1_m_detection_proba = log(1 - detection_proba);

      // here the logarithms of 1!, 2!, 3!, ..., upper_limit_size_identical_sequence_clusters! are stored
      for (ii in 2:upper_limit_size_identical_sequence_clusters) log_sums[ii] += log_sums[ii-1] + log(ii);

      // determine the logarithm of P(Y=0)
      for (ii in 1:upper_limit_size_identical_sequence_clusters) result_temp[ii] = ii*log_1_m_detection_proba + distribution_size_ident_seq_tree[ii];

      distribution_size_ident_seq_tree_detection_0[1] = log_sum_exp(result_temp);

      // determine the logarithm of P(Y = n) for all n contained in clusters_size
      for (jj in clusters_size) {

        // vector needed to store intermediate results
        result_temp = rep_array(-10000.0, upper_limit_size_identical_sequence_clusters);

        for (ii in jj:upper_limit_size_identical_sequence_clusters) result_temp[ii] = log_sums[ii] - log_sums[jj] - log_sums[max(ii-jj, 1)] + jj * log_detection_proba + (ii-jj) * log_1_m_detection_proba + distribution_size_ident_seq_tree[ii];

        // determine the logarithm of P(Y = jj)
        distribution_size_ident_seq_tree_detection_0[index] = log_sum_exp(tail(result_temp,upper_limit_size_identical_sequence_clusters-jj+1));

        index += 1;

      }

      // determine the logarithm of P(Y > 0)
      log_proba_detect_ident_seq_tree_size_bigger_zero = log(1 - exp(distribution_size_ident_seq_tree_detection_0[1]));

      // determine the logarithm of P(Z=n) for all n contained in clusters_size
      for (ii in 1:n_clusters_size) distribution_size_ident_seq_tree_detection[ii] = distribution_size_ident_seq_tree_detection_0[ii+1] -  log_proba_detect_ident_seq_tree_size_bigger_zero;

    }

    for (ii in 1:n_clusters_size) {

      likelihood[ii] = clusters_freq[ii] * distribution_size_ident_seq_tree_detection[ii];

    }

    return(sum(likelihood));

  }

}
