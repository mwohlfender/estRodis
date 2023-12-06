
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
  
  
  
  real function_extinction(real x, real R, real k, real mutation_proba) {
    
    return(x - 1 / ((1 + ((R * (1 - mutation_proba) * (1 - x)) / k))^k));
    
  }
  
  
  
  real derivative_function_extinction(real x, real R, real k, real mutation_proba) {
    
    return(1 + k * 1 / ((1 + ((R * (1 - mutation_proba) * (1 - x)) / k))^(k+1)) * (- R / k * (1 - mutation_proba)));
    
  }
  
  
  
  real absolute_value(real x) {
    
    real result = 0;
    
    if (x >= 0) {
      
      result = x;
      
    } else {
      
      result = -x;
      
    }
    
    return(result);
    
  }
  
  
  
  real root_newton_extinction(real R, real k, real mutation_proba, int max_iter, real tol, real initial_guess) {
    
    real x_n = -1;
    real x_n_p_one = initial_guess;
    
    int index = 1;
    
    while(absolute_value(x_n_p_one - x_n) >= tol && index <= max_iter) {
      
      x_n = x_n_p_one;
      
      x_n_p_one = x_n - function_extinction(x_n, R, k, mutation_proba) / derivative_function_extinction(x_n, R, k, mutation_proba);
      
      index = index + 1;
      
    }
    
    return(x_n_p_one);
    
  }
  
  
  
  real estRodis_stan_likelihood_log_smart(int[] clusters_size, int[] clusters_freq, int upper_limit_cluster_size, real R, real k, real mutation_proba, real detection_proba, real tol_inf) {
    
    // here the number of different identical sequence cluster sizes is stored
    int n_clusters_size = num_elements(clusters_size);
    
    // maximum of largest element of `clusters_size` and `upper_limit_cluster_size`
    int max_cluster_size = max(max(clusters_size), upper_limit_cluster_size);
    
    // upper limit up to which the probability P(X = n) is calculated
    int upper_limit_size_identical_sequence_clusters = max(determine_scaling_factor(2 / detection_proba) * max_cluster_size, max_cluster_size + 500);
    
    // logarithm of R * (1 -  mutation_proba) / k
    real log_R_d_k = log(R * (1 - mutation_proba) / k);
    
    // logarithm of 1 + R * (1 - mutation_proba) / k
    real log_one_p_R_d_k = log(1 + R * (1 - mutation_proba) / k);
    
    // logarithm of detection_proba
    real log_detection_proba;
    
    // logarithm of 1 - detection_proba
    real log_1_m_detection_proba;
    
    // logarithms of 1!, 2!, 3!, ..., upper_limit_size_identical_sequence_clusters!
    real log_sums[upper_limit_size_identical_sequence_clusters] = rep_array(0.0, upper_limit_size_identical_sequence_clusters);
    
    // intermediate result needed during the computation
    real result_temp[upper_limit_size_identical_sequence_clusters] = rep_array(-10000.0, upper_limit_size_identical_sequence_clusters);
    
    // logarithm of P(X=n) for all 1 <= n <= `upper_limit_size_identical_sequence_clusters`
    real distribution_size_ident_seq_cluster[upper_limit_size_identical_sequence_clusters] = rep_array(-10000.0, upper_limit_size_identical_sequence_clusters);
    
    // logarithm of P(Y=0)
    real log_proba_detect_ident_seq_tree_size_zero = -10000.0;
    
    // logarithm of P(Y>0)
    real log_proba_detect_ident_seq_tree_size_bigger_zero = -10000.0;
    
    // logarithm of P(Z=n) for all n contained in `clusters_size` and P(Z=inf)
    real distribution_size_ident_seq_cluster_detection[max_cluster_size+1] = rep_array(-10000.0, max_cluster_size+1);
    
    // index of an element in `clusters_size`
    int index = 1;
    
    // ingredient of extinction probability of identical sequence cluster
    real extinction_proba_zero = root_newton_extinction(R, k, mutation_proba, 1000, 1e-12, 0.5);
    
    // extinction probability of identical sequence cluster
    real extinction_proba = 0;
    
    // P(Z = inf)
    real proba_inf = 0;
    
    // check whether calculation shall be done or not
    int do_calculation = 0;
    
    // lowest cluster size that is considered an infinite cluster
    int start_of_infinity = max_cluster_size;
    
    // result
    real result = 0;
    
    // determine P(X=n) for all 1 <= n <= `upper_limit_size_identical_sequence_clusters`
    for(ii in 1:upper_limit_size_identical_sequence_clusters) {
      
      // determine P(X = ii)
      distribution_size_ident_seq_cluster[ii] = lgamma(k*ii + ii - 1) - lgamma(k*ii) - lgamma(ii + 1) + (ii - 1) * log_R_d_k - (k*ii + ii - 1)*log_one_p_R_d_k;
      
    }
    
    if (detection_proba == 1.0) {
      
      // determine P(Z=n) for all 1 <= n <= `max_cluster_size`
      for (ii in 1:max_cluster_size) distribution_size_ident_seq_cluster_detection[ii] = distribution_size_ident_seq_cluster[ii];
      
      extinction_proba = extinction_proba_zero;
      
      proba_inf = 1 - extinction_proba;
      
      if (proba_inf > 0.0) {
        
        distribution_size_ident_seq_cluster_detection[max_cluster_size+1] = log(proba_inf);
        
      } else {
        
        distribution_size_ident_seq_cluster_detection[max_cluster_size+1] = -10000.0;
        
      }
      
      for (jj in 2:max_cluster_size) {
        
        if (start_of_infinity == max_cluster_size) {
          
          if ((extinction_proba - sum(exp(distribution_size_ident_seq_cluster_detection[1:(jj-1)]))) / (1 - sum(exp(distribution_size_ident_seq_cluster_detection[1:(jj-1)]))) < tol_inf) {
            
            start_of_infinity = jj;
            
          }
          
        }
        
      }
      
    } else {
      
      // logarithm of detection_proba
      log_detection_proba = log(detection_proba);
      
      // logarithm of 1 - detection_proba
      log_1_m_detection_proba = log(1 - detection_proba);
      
      // determine logarithm of n! for all 1 <= n <= `upper_limit_size_identical_sequence_clusters`
      for (ii in 2:upper_limit_size_identical_sequence_clusters) log_sums[ii] += log_sums[ii-1] + log(ii);
      
      // determine logarithm of P(Y>0)
      for (ii in 1:upper_limit_size_identical_sequence_clusters) result_temp[ii] = ii*log_1_m_detection_proba + distribution_size_ident_seq_cluster[ii];
      
      log_proba_detect_ident_seq_tree_size_zero = log_sum_exp(result_temp);
      
      log_proba_detect_ident_seq_tree_size_bigger_zero = log(1 - exp(log_proba_detect_ident_seq_tree_size_zero));
      
      extinction_proba = floor(1e12*(extinction_proba_zero - exp(log_proba_detect_ident_seq_tree_size_zero)) / exp(log_proba_detect_ident_seq_tree_size_bigger_zero))/1e12;
      
      proba_inf = 1 - extinction_proba;
      
      if (proba_inf > 0.0) {
        
        distribution_size_ident_seq_cluster_detection[max_cluster_size+1] = log(proba_inf);
        
      } else {
        
        distribution_size_ident_seq_cluster_detection[max_cluster_size+1] = -10000.0;
        
      }
      
      
      // determine logarithm of P(Z = n) for all n contained in `clusters_size`
      for (jj in 1:max_cluster_size) {
        
        if (jj == 1 ) {
          
          do_calculation = 1;
          
          if (jj == clusters_size[index]) {
            
            index = index + 1;
            
          }
          
        }
        
        if (jj >= 2) {
          
          if (R * (1 - mutation_proba) < 1) {
            
            if (jj == clusters_size[index]) {
              
              do_calculation = 1;
              
              index = index + 1;
              
            } else {
              
              do_calculation = 0;
              
            }
            
          } else {
            
            if (jj < start_of_infinity) {
              
              if ((extinction_proba - sum(exp(distribution_size_ident_seq_cluster_detection[1:(jj-1)]))) / (1 - sum(exp(distribution_size_ident_seq_cluster_detection[1:(jj-1)]))) < tol_inf) {
                
                do_calculation = 0;
                
                start_of_infinity = jj;
                
                distribution_size_ident_seq_cluster_detection[max_cluster_size+1] = log(1 - sum(exp(distribution_size_ident_seq_cluster_detection[1:max_cluster_size])));
                
              } else {
                
                do_calculation = 1;
                
              }
              
            } else {
              
              do_calculation = 0;
              
            }
            
          }
          
        }
        
        if (do_calculation == 1) {
          
          // vector needed to store intermediate results
          result_temp = rep_array(-10000.0, upper_limit_size_identical_sequence_clusters);
          
          for (ii in jj:max_cluster_size) result_temp[ii] = log_sums[ii] - log_sums[jj] - log_sums[max(ii-jj, 1)] + jj * log_detection_proba + (ii-jj) * log_1_m_detection_proba + distribution_size_ident_seq_cluster[ii];
          
          // determine logarithm of P(Z=jj)
          distribution_size_ident_seq_cluster_detection[jj] = log_sum_exp(tail(result_temp, upper_limit_size_identical_sequence_clusters-jj+1)) - log_proba_detect_ident_seq_tree_size_bigger_zero;
          
        }
        
      }
      
    }
    
    for (ii in 1:n_clusters_size) {
      
      if (clusters_size[ii] < start_of_infinity) {
        
        result = result + clusters_freq[ii] * distribution_size_ident_seq_cluster_detection[clusters_size[ii]];
        
      } else {
        
        result = result + clusters_freq[ii] * distribution_size_ident_seq_cluster_detection[max_cluster_size+1];
        
      }
      
    }
    
    // print(start_of_infinity);
    //
    // print(distribution_size_ident_seq_cluster_detection[1:125]);
    //
    // print(distribution_size_ident_seq_cluster_detection[max_cluster_size+1]);
    
    return(result);
    
  }
  
  
  
  real stan_cluster_size_proba_bigger_zero(int max_cluster_size, real R, real k, real mutation_proba, real detection_proba) {
    
    // here the upper limit up to which the probability of observing an identical sequence cluster of a certain size is calculated is stored
    int upper_limit_size_identical_sequence_clusters = max(determine_scaling_factor(2 / detection_proba) * max_cluster_size, max_cluster_size + 500);
    
    // here the logarithm of R * (1 -  mutation_proba) / k is stored
    real log_R_d_k = log(R * (1 - mutation_proba) / k);
    
    // here the logarithm of 1 + R * (1 - mutation_proba) / k is stored
    real log_one_p_R_d_k = log(1 + R * (1 - mutation_proba) / k);
    
    // here the logarithm of detection_proba will be stored (if detection_proba < 1)
    real log_detection_proba;
    
    // here the logarithm of 1 - detection_proba will be stored (if detection_proba < 1)
    real log_1_m_detection_proba;
    
    // here the logarithms of 1!, 2!, 3!, ..., upper_limit_size_identical_sequence_clusters! will be stored (if detection_proba < 1)
    real log_sums[upper_limit_size_identical_sequence_clusters] = rep_array(0.0, upper_limit_size_identical_sequence_clusters);
    
    // here intermediate results needed during the computations will be stored
    real result_temp[upper_limit_size_identical_sequence_clusters] = rep_array(-10000.0, upper_limit_size_identical_sequence_clusters);
    
    // here the logarithms of P(X=1), P(X=2), ..., P(X=upper_limit_size_identical_sequence_clusters) will be stored
    real distribution_size_ident_seq_tree[upper_limit_size_identical_sequence_clusters] = rep_array(-10000.0, upper_limit_size_identical_sequence_clusters);
    
    // here P(Y>0) will be stored
    real proba_detect_ident_seq_tree_size_bigger_zero = 1;
    
    // P(X=1), P(X=2), ..., P(X=upper_limit_size_identical_sequence_clusters) are determined
    for(ii in 1:upper_limit_size_identical_sequence_clusters) {
      
      // P(X = ii) is determined
      distribution_size_ident_seq_tree[ii] = lgamma(k*ii + ii - 1) - lgamma(k*ii) - lgamma(ii + 1) + (ii - 1) * log_R_d_k - (k*ii + ii - 1)*log_one_p_R_d_k;
      
    }
    
    if (detection_proba == 1.0) {
      
      // determine P(Y > 0)
      proba_detect_ident_seq_tree_size_bigger_zero = 1;
      
    } else {
      
      // here the logarithm of detection_proba is stored
      log_detection_proba = log(detection_proba);
      
      // here the logarithm of 1 - detection_proba is stored
      log_1_m_detection_proba = log(1 - detection_proba);
      
      // here the logarithms of 1!, 2!, 3!, ..., upper_limit_size_identical_sequence_clusters! are stored
      for (ii in 2:upper_limit_size_identical_sequence_clusters) log_sums[ii] += log_sums[ii-1] + log(ii);
      
      // determine the logarithm of P(Y=0)
      for (ii in 1:upper_limit_size_identical_sequence_clusters) result_temp[ii] = ii*log_1_m_detection_proba + distribution_size_ident_seq_tree[ii];
      
      // determine P(Y > 0)
      proba_detect_ident_seq_tree_size_bigger_zero = 1 - exp(log_sum_exp(result_temp));
      
    }
    
    return(proba_detect_ident_seq_tree_size_bigger_zero);
    
  }
  
  
  
  real stan_proba_inf(int max_cluster_size, real R, real k, real mutation_proba, real detection_proba) {
    
    // ingredient of extinction probability of identical sequence cluster
    real extinction_proba_zero = root_newton_extinction(R, k, mutation_proba, 1000, 1e-12, 0.5);
    
    // P(Y > 0)
    real proba_bigger_zero = stan_cluster_size_proba_bigger_zero(max_cluster_size, R, k, mutation_proba, detection_proba);
    
    // P(Z = inf)
    real result = 1 - ((extinction_proba_zero - (1 - proba_bigger_zero)) / proba_bigger_zero);
    
    return(result);
    
  }
  
  
  real stan_likelihood_log(int[] clusters_size, int[] clusters_freq, real R, real k, real mutation_proba, real detection_proba) {
    
    // number of different identical sequence cluster sizes
    int n_clusters_size = num_elements(clusters_size);
    
    // upper limit up to which the probability P(X = n) is calculated
    int upper_limit_size_identical_sequence_clusters = max(determine_scaling_factor(2 / detection_proba) * max(clusters_size), max(clusters_size) + 500);
    
    // logarithm of R * (1 -  mutation_proba) / k
    real log_R_d_k = log(R * (1 - mutation_proba) / k);
    
    // logarithm of 1 + R * (1 - mutation_proba) / k
    real log_one_p_R_d_k = log(1 + R * (1 - mutation_proba) / k);
    
    // logarithm of detection_proba
    real log_detection_proba = log(detection_proba);
    
    // logarithm of 1 - detection_proba
    real log_1_m_detection_proba = log(1 - detection_proba);
    
    // logarithm of probability P(Y > 0)
    real log_proba_detect_ident_seq_tree_size_bigger_zero;
    
    // logarithms of 1!, 2!, 3!, ..., upper_limit_size_identical_sequence_clusters!
    real log_sums[upper_limit_size_identical_sequence_clusters] = rep_array(0.0, upper_limit_size_identical_sequence_clusters);
    
    // intermediate result needed during the computation
    real result_temp[upper_limit_size_identical_sequence_clusters] = rep_array(-10000.0, upper_limit_size_identical_sequence_clusters);
    
    // logarithms of P(X=1), P(X=2), ..., P(X=upper_limit_size_identical_sequence_clusters)
    real distribution_size_ident_seq_tree[upper_limit_size_identical_sequence_clusters] = rep_array(-10000.0, upper_limit_size_identical_sequence_clusters);
    
    // logarithms of P(Z=n) for all n contained in `clusters_size`
    real distribution_size_ident_seq_tree_detection[n_clusters_size] = rep_array(-10000.0, n_clusters_size);
    
    // result
    real likelihood[n_clusters_size] = rep_array(-10000.0, n_clusters_size);
    
    // P(X=1), P(X=2), ..., P(X=upper_limit_size_identical_sequence_clusters) are determined
    for(ii in 1:upper_limit_size_identical_sequence_clusters) {
      
      // P(X = ii) is determined
      distribution_size_ident_seq_tree[ii] = lgamma(k*ii + ii - 1) - lgamma(k*ii) - lgamma(ii + 1) + (ii - 1) * log_R_d_k - (k*ii + ii - 1)*log_one_p_R_d_k;
      
    }
    
    if (detection_proba == 1.0) {
      
      for (ii in 1:n_clusters_size) distribution_size_ident_seq_tree_detection[ii] = distribution_size_ident_seq_tree[clusters_size[ii]];
      
    } else {
      
      // determine logarithms of 1!, 2!, 3!, ..., upper_limit_size_identical_sequence_clusters!
      for (ii in 2:upper_limit_size_identical_sequence_clusters) log_sums[ii] += log_sums[ii-1] + log(ii);
      
      // determine logarithm of P(Y=0)
      for (ii in 1:upper_limit_size_identical_sequence_clusters) result_temp[ii] = ii*log_1_m_detection_proba + distribution_size_ident_seq_tree[ii];
      
      // determine logarithm of P(Y > 0)
      log_proba_detect_ident_seq_tree_size_bigger_zero = log(1 - exp(log_sum_exp(result_temp)));
      
      // determine logarithm of P(Y = n) for all n contained in `clusters_size`
      for (jj in 1:n_clusters_size) {
        
        if (clusters_freq[jj] != 0) {
          
          // vector needed to store intermediate results
          result_temp = rep_array(-10000.0, upper_limit_size_identical_sequence_clusters);
          
          for (ii in clusters_size[jj]:upper_limit_size_identical_sequence_clusters) result_temp[ii] = log_sums[ii] - log_sums[clusters_size[jj]] - log_sums[max(ii-clusters_size[jj], 1)] + clusters_size[jj] * log_detection_proba + (ii-clusters_size[jj]) * log_1_m_detection_proba + distribution_size_ident_seq_tree[ii];
          
          // determine logarithms of P(Z=clusters_size[1]), P(Z=clusters_size[2]), ..., P(Z=clusters_size[n_clusters_size])
          distribution_size_ident_seq_tree_detection[jj] = log_sum_exp(tail(result_temp, upper_limit_size_identical_sequence_clusters-clusters_size[jj]+1)) - log_proba_detect_ident_seq_tree_size_bigger_zero;
          
        }
        
      }
      
    }
    
    for (ii in 1:n_clusters_size) {
      
      likelihood[ii] = clusters_freq[ii] * distribution_size_ident_seq_tree_detection[ii];
      
    }
    
    // print(distribution_size_ident_seq_tree_detection);
    
    return(sum(likelihood));
    
  }
  
  
  
  real estRodis_stan_likelihood_log_inf_improved(int[] clusters_size, int[] clusters_freq, int max_cluster_size, real R, real k, real mutation_proba, real detection_proba) {
    
    // here the number of different identical sequence cluster sizes is stored
    int n_clusters_size = num_elements(clusters_size);
    
    // here the sizes of the finite clusters will be stored
    int clusters_size_finite[n_clusters_size - 1] = rep_array(0, n_clusters_size - 1);
    
    // here the frequencies of the finite clusters will be stored
    int clusters_freq_finite[n_clusters_size - 1] = rep_array(0, n_clusters_size - 1);
    
    // P(Z = inf)
    real proba_inf = stan_proba_inf(max_cluster_size, R, k, mutation_proba, detection_proba);
    
    // log(P(Z = inf))
    real log_proba_inf = 0.0;
    
    // here an intermediate result will be stored
    real likelihood_finite = 0.0;
    
    // here the result will be stored
    real result = 0.0;
    
    if ((1.0 - mutation_proba) * R < 1.0) {
      
      result = stan_likelihood_log(clusters_size, clusters_freq, R, k, mutation_proba, detection_proba);
      
    } else {
      
      for (ii in 1:(n_clusters_size-1)) {
        
        clusters_size_finite[ii] = clusters_size[ii];
        
        clusters_freq_finite[ii] = clusters_freq[ii];
        
      }
      
      likelihood_finite = stan_likelihood_log(clusters_size_finite, clusters_freq_finite, R, k, mutation_proba, detection_proba);
      
      if (proba_inf > 0.0) {
        
        log_proba_inf = log(proba_inf);
        
      }
      
      result = likelihood_finite + clusters_freq[n_clusters_size] * log_proba_inf;
      
    }
    
    return(result);
    
  }
  
}
