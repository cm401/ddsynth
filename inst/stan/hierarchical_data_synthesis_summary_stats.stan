// hierarchical_data_synthesis_summary_stats.stan
functions {
  // CDF for each distribution
  real dist_cdf_fun(real x, int dist_type, real loc, real phi) {
    if (dist_type == 1) {  // lognormal
      return lognormal_cdf(x | loc, phi);
    } else if (dist_type == 2) {  // gamma
      real mean_d = exp(loc);
      real shape = phi;
      real rate = shape / mean_d;
      return gamma_cdf(x | shape, rate);
    } else if (dist_type == 3) {  // weibull
      real scale = exp(loc);
      real shape = phi;
      return weibull_cdf(x | shape, scale);
    }
    return 0;
  }
  
  // "log pdf" for each distribution
  real dist_logpdf_fun(real x, int dist_type, real loc, real phi) {
    if (dist_type == 1) {  // lognormal
      return lognormal_lpdf(x | loc, phi);
    } else if (dist_type == 2) {  // gamma
      real mean_d = exp(loc);
      real shape = phi;
      real rate = shape / mean_d;
      return gamma_lpdf(x | shape, rate);
    } else if (dist_type == 3) {  // weibull
      real scale = exp(loc);
      real shape = phi;
      return weibull_lpdf(x | shape, scale);
    }
    return 0;
  }
  
  // Log likelihood for order statistic:
  // k-th order statistic out of n observations
  real order_stat_logpdf_fun(real x, int n, int k, int dist_type, real loc, real phi) {
    real F = dist_cdf_fun(x, dist_type, loc, phi);
    real log_f = dist_logpdf_fun(x, dist_type, loc, phi);
  
    F = fmax(fmin(F, 0.99999), 0.00001);
  
    real log_dens = lchoose(n, k) + log_f + (k - 1) * log(F) + (n - k) * log1m(F);
  
    return log_dens;
  }

  // Helper function to compute gamma quantile approximation
  real gamma_quantile_approx(real p, real shape, real scale) {
    // Wilson-Hilferty approximation for gamma quantiles
    real z = inv_Phi(p);
    if (shape > 1) {
      return shape * scale * pow(1 - 1.0/(9*shape) + z/(3*sqrt(shape)), 3);
    } else {
      return shape * scale * (1 + z / sqrt(shape));   // Fallback for small shape
    }
  }
}

data {
  int<lower=1> n_datasets;              // Number of datasets
  array[n_datasets] int<lower=1> n_obs; // Sample sizes for each dataset
  array[n_datasets] int<lower=1,upper=5> summary_type; // 1=median+range, 2=median+IQR, 3=mean+sd, 4=raw freq table, 5=interval-censored freq table
  int<lower=1,upper=3> dist_type;       // 1=lognormal, 2=gamma, 3=weibull
  
  // Observed summaries - organized by dataset (used for summary_type 1, 2, 3)
  array[n_datasets] real<lower=0> obs_stat1;     // median or mean
  array[n_datasets] real<lower=0> obs_stat2;     // min, q25, or sd
  array[n_datasets] real<lower=0> obs_stat3;     // max, q75, or placeholder

  // --- Frequency table data for summary_type == 4 ---
  // All datasets' frequency tables are stored in flat arrays.
  // For dataset d, its entries occupy indices freq_start[d] .. freq_start[d] + freq_len[d] - 1.
  int<lower=0> n_freq_total;                        // Total number of (value, count) pairs across all type-4 datasets
  array[n_freq_total] real<lower=0> freq_value;     // Observed day values for type 4 (must be > 0 for continuous distributions)
  array[n_freq_total] int<lower=1>  freq_count;     // Number of individuals with that day value (types 4 and 5)
  array[n_datasets]   int<lower=0>  freq_start;     // 1-based start index into freq arrays for dataset d
  array[n_datasets]   int<lower=0>  freq_len;       // Number of distinct values for dataset d (0 if not type 4 or 5)

  // Interval bounds for summary_type == 5 (interval-censored frequency table).
  // For type 4 datasets these arrays are ignored (populate with zeros).
  // When freq_lower[i] == freq_upper[i] the contribution falls back to the log-PDF.
  array[n_freq_total] real<lower=0> freq_lower;    // lower bound of censoring interval
  array[n_freq_total] real<lower=0> freq_upper;    // upper bound of censoring interval

  // Prior hyperparameters for mu0 ~ normal(mu0_mean, mu0_sd)
  real mu0_mean;
  real<lower=0> mu0_sd;
  
  // Prior hyperparameters for log_tau ~ normal(log_tau_mean, log_tau_sd)
  real log_tau_mean;
  real<lower=0> log_tau_sd;
  
  // Prior hyperparameters for log_phi ~ normal(log_phi_mean, log_phi_sd)
  real log_phi_mean;
  real<lower=0> log_phi_sd;
}

transformed data {
  real z_q25 = -0.6745;       
  real z_q75 = 0.6745;       
  real z_q90 = 1.28155;
  real z_q95 = 1.64485;
  
  // Validate data
  for (d in 1:n_datasets) {
    if (summary_type[d] == 1) {
      if (obs_stat2[d] > obs_stat1[d] || obs_stat1[d] > obs_stat3[d]) {
        reject("For summary_type=1, must have min <= median <= max");
      }
    } else if (summary_type[d] == 2) {
      if (obs_stat2[d] > obs_stat1[d] || obs_stat1[d] > obs_stat3[d]) {
        reject("For summary_type=2, must have q25 <= median <= q75");
      }
    } else if (summary_type[d] == 4) {
      if (freq_len[d] == 0) {
        reject("For summary_type=4, freq_len must be > 0");
      }
    } else if (summary_type[d] == 5) {
      if (freq_len[d] == 0) {
        reject("For summary_type=5, freq_len must be > 0");
      }
    }
  }
}

parameters {
  real mu0;                        // Population mean (location)
  real log_tau;                    // Log of between-study SD
  real log_phi;                    // Log of distribution-specific parameter
  vector[n_datasets] loc_d_raw;    // Non-centered parameterization
}

transformed parameters {
  real<lower=0> tau = exp(log_tau);
  real<lower=0> phi = exp(log_phi);
  vector[n_datasets] loc_d = mu0 + tau * loc_d_raw;
}

model {
  mu0 ~ normal(mu0_mean, mu0_sd);
  log_tau ~ normal(log_tau_mean, log_tau_sd);
  log_phi ~ normal(log_phi_mean, log_phi_sd);
  
  loc_d_raw ~ std_normal();
  
  for (d in 1:n_datasets) {
    real loc = loc_d[d];
    int n = n_obs[d];
    
    if (summary_type[d] == 1) {  // median + range (min, max)
      int k_median = (n + 1) %/% 2;
      target += order_stat_logpdf_fun(obs_stat1[d], n, k_median, dist_type, loc, phi);
      target += order_stat_logpdf_fun(obs_stat2[d], n, 1, dist_type, loc, phi);
      target += order_stat_logpdf_fun(obs_stat3[d], n, n, dist_type, loc, phi);
    }
    
    else if (summary_type[d] == 2) {  // median + IQR (q25, q75)
      int k_median = (n + 1) %/% 2;
      target += order_stat_logpdf_fun(obs_stat1[d], n, k_median, dist_type, loc, phi);
      
      int k_q25 = (n + 1) %/% 4;
      if (k_q25 < 1) k_q25 = 1;
      target += order_stat_logpdf_fun(obs_stat2[d], n, k_q25, dist_type, loc, phi);

      int k_q75 = (3 * (n + 1)) %/% 4;
      if (k_q75 <= k_q25) k_q75 = k_q25 + 1;
      if (k_q75 > n) k_q75 = n;
      target += order_stat_logpdf_fun(obs_stat3[d], n, k_q75, dist_type, loc, phi);
    }
    
    else if (summary_type[d] == 3) {  // mean + sd
      real expected_mean;
      real expected_sd;
      real se_mean;
      real se_sd;
      
      if (dist_type == 1) {  // lognormal
        expected_mean = exp(loc + phi^2 / 2);
        real var_ = (exp(phi^2) - 1) * exp(2 * loc + phi^2);
        expected_sd = sqrt(var_);
        se_mean = expected_sd / sqrt(n);
        se_sd = expected_sd / sqrt(2 * (n - 1));
        
      } else if (dist_type == 2) {  // gamma
        real mean_d = exp(loc);
        real shape = phi;
        real scale_param = mean_d / shape;
        expected_mean = mean_d;
        expected_sd = sqrt(shape * scale_param^2);  
        se_mean = expected_sd / sqrt(n);
        se_sd = expected_sd / sqrt(2 * (n - 1));
      } else if (dist_type == 3) {  // weibull
        real scale = exp(loc);
        real shape = phi;
        expected_mean = scale * tgamma(1 + 1.0 / shape);
        real var_ = scale^2 * (tgamma(1 + 2.0 / shape) - pow(tgamma(1 + 1.0 / shape), 2));
        expected_sd = sqrt(var_);
        se_mean = expected_sd / sqrt(n);
        se_sd = expected_sd / sqrt(2 * (n - 1));
      }
      
      obs_stat1[d] ~ normal(expected_mean, se_mean);
      obs_stat2[d] ~ normal(expected_sd, se_sd);
    }
    
    else if (summary_type[d] == 4) {  // raw frequency table
      // Direct likelihood: for each distinct observed value, add count * log_pdf(value).
      // This is equivalent to fitting the distribution directly to all individual observations,
      // but using the compressed frequency-table representation.
      int s = freq_start[d];
      int len = freq_len[d];
      for (i in s:(s + len - 1)) {
        target += freq_count[i] * dist_logpdf_fun(freq_value[i], dist_type, loc, phi);
      }
    }

    else if (summary_type[d] == 5) {  // interval-censored frequency table
      // Likelihood: count * log[ F(upper) - F(lower) ] for each interval.
      // When lower == upper (point observation), falls back to count * log_pdf(value)
      // to avoid log(0).
      int s = freq_start[d];
      int len = freq_len[d];
      for (i in s:(s + len - 1)) {
        if (freq_lower[i] == freq_upper[i]) {
          target += freq_count[i] * dist_logpdf_fun(freq_lower[i], dist_type, loc, phi);
        } else {
          real cdf_u = dist_cdf_fun(freq_upper[i], dist_type, loc, phi);
          real cdf_l = dist_cdf_fun(freq_lower[i], dist_type, loc, phi);
          target += freq_count[i] * log(cdf_u - cdf_l);
        }
      }
    }
  }
}

generated quantities {
  real pred_mean;
  real pred_median;
  real pred_q25;
  real pred_q75;
  real pred_q90;
  real pred_q95;
  real pred_sd;
  vector[n_datasets * 3] log_lik;
  
  // Predicted quantities:
  // - When n_datasets < 5, tau is not identifiable from data (prior-dominated).
  //   pred_* are computed at mu0 directly to avoid tau^2/2 inflation.
  //   See: Higgins & Thompson (2002) doi:10.1002/sim.1186
  //        Gelman (2006) doi:10.1214/06-BA117A
  //        Rover et al. (2021) doi:10.1002/jrsm.1475
  // - When n_datasets >= 5, tau is identifiable; sample from Normal(mu0, tau)
  //   to include between-study heterogeneity. L=2000 for MC stability.
  {
    if (n_datasets < 5) {
      // Use mu0 directly — tau unidentifiable with fewer than 5 studies
      if (dist_type == 1) {  // lognormal
        pred_mean   = exp(mu0 + phi^2 / 2);
        pred_median = exp(mu0);
        pred_q25    = exp(mu0 - 0.6745  * phi);
        pred_q75    = exp(mu0 + 0.6745  * phi);
        pred_q90    = exp(mu0 + 1.28155 * phi);
        pred_q95    = exp(mu0 + 1.64485 * phi);
        pred_sd     = sqrt((exp(phi^2) - 1) * exp(2 * mu0 + phi^2));

      } else if (dist_type == 2) {  // gamma
        real mean_d      = exp(mu0);
        real scale_param = mean_d / phi;
        pred_mean   = mean_d;
        pred_sd     = sqrt(mean_d * scale_param);
        pred_median = gamma_quantile_approx(0.5,  phi, scale_param);
        pred_q25    = gamma_quantile_approx(0.25, phi, scale_param);
        pred_q75    = gamma_quantile_approx(0.75, phi, scale_param);
        pred_q90    = gamma_quantile_approx(0.90, phi, scale_param);
        pred_q95    = gamma_quantile_approx(0.95, phi, scale_param);

      } else if (dist_type == 3) {  // weibull
        real scale  = exp(mu0);
        pred_mean   = scale * tgamma(1 + 1.0 / phi);
        pred_median = scale * pow(log(2),         1.0 / phi);
        pred_q25    = scale * pow(log(4.0 / 3.0), 1.0 / phi);
        pred_q75    = scale * pow(log(4.0),       1.0 / phi);
        pred_q90    = scale * pow(log(10.0),      1.0 / phi);
        pred_q95    = scale * pow(log(20.0),      1.0 / phi);
        pred_sd     = sqrt(scale^2 * (tgamma(1 + 2.0/phi) - pow(tgamma(1 + 1.0/phi), 2)));
      }

    } else {
      // n_datasets >= 5: tau identifiable; include between-study heterogeneity
      // via Monte Carlo integration over Normal(mu0, tau). L=2000 for stability.
      int L = 2000;
      vector[L] means;
      vector[L] medians;
      vector[L] q25s;
      vector[L] q75s;
      vector[L] q90s;
      vector[L] q95s;
      vector[L] sds;

      for (l in 1:L) {
        real loc_sample = normal_rng(mu0, tau);

        if (dist_type == 1) {  // lognormal
          means[l]   = exp(loc_sample + phi^2 / 2);
          medians[l] = exp(loc_sample);
          q25s[l]    = exp(loc_sample - 0.6745  * phi);
          q75s[l]    = exp(loc_sample + 0.6745  * phi);
          q90s[l]    = exp(loc_sample + 1.28155 * phi);
          q95s[l]    = exp(loc_sample + 1.64485 * phi);
          sds[l]     = sqrt((exp(phi^2) - 1) * exp(2 * loc_sample + phi^2));

        } else if (dist_type == 2) {  // gamma
          real mean_d      = exp(loc_sample);
          real scale_param = mean_d / phi;
          means[l]   = mean_d;
          sds[l]     = sqrt(mean_d * scale_param);
          medians[l] = gamma_quantile_approx(0.5,  phi, scale_param);
          q25s[l]    = gamma_quantile_approx(0.25, phi, scale_param);
          q75s[l]    = gamma_quantile_approx(0.75, phi, scale_param);
          q90s[l]    = gamma_quantile_approx(0.90, phi, scale_param);
          q95s[l]    = gamma_quantile_approx(0.95, phi, scale_param);

        } else if (dist_type == 3) {  // weibull
          real scale    = exp(loc_sample);
          means[l]      = scale * tgamma(1 + 1.0 / phi);
          medians[l]    = scale * pow(log(2),         1.0 / phi);
          q25s[l]       = scale * pow(log(4.0 / 3.0), 1.0 / phi);
          q75s[l]       = scale * pow(log(4.0),       1.0 / phi);
          q90s[l]       = scale * pow(log(10.0),      1.0 / phi);
          q95s[l]       = scale * pow(log(20.0),      1.0 / phi);
          real var_weib = scale^2 * (tgamma(1 + 2.0/phi) - pow(tgamma(1 + 1.0/phi), 2));
          sds[l]        = sqrt(var_weib);
        }
      }

      pred_mean   = mean(means);
      pred_median = mean(medians);
      pred_q25    = mean(q25s);
      pred_q75    = mean(q75s);
      pred_q90    = mean(q90s);
      pred_q95    = mean(q95s);
      pred_sd     = mean(sds);
    }
  }
  
  // Log likelihood
  {
    int idx = 1;
    
    for (d in 1:n_datasets) {
      real loc = loc_d[d];
      int n = n_obs[d];
      
      if (summary_type[d] == 1) {  // median + range
        int k_median = (n + 1) %/% 2;
        log_lik[idx]     = order_stat_logpdf_fun(obs_stat1[d], n, k_median, dist_type, loc, phi);
        log_lik[idx + 1] = order_stat_logpdf_fun(obs_stat2[d], n, 1,         dist_type, loc, phi);
        log_lik[idx + 2] = order_stat_logpdf_fun(obs_stat3[d], n, n,         dist_type, loc, phi);
        
      } else if (summary_type[d] == 2) {  // median + IQR
        int k_median = (n + 1) %/% 2;
        int k_q25    = (n + 1) %/% 4;
        if (k_q25 < 1) k_q25 = 1;
        int k_q75    = (3 * (n + 1)) %/% 4;
        if (k_q75 <= k_q25) k_q75 = k_q25 + 1;
        if (k_q75 > n) k_q75 = n;
        
        log_lik[idx]     = order_stat_logpdf_fun(obs_stat1[d], n, k_median, dist_type, loc, phi);
        log_lik[idx + 1] = order_stat_logpdf_fun(obs_stat2[d], n, k_q25,    dist_type, loc, phi);
        log_lik[idx + 2] = order_stat_logpdf_fun(obs_stat3[d], n, k_q75,    dist_type, loc, phi);
        
      } else if (summary_type[d] == 3) {  // mean + sd
        real expected_mean;
        real expected_sd;
        real se_mean;
        real se_sd;
        
        if (dist_type == 1) {
          expected_mean = exp(loc + phi^2 / 2);
          real var_ = (exp(phi^2) - 1) * exp(2 * loc + phi^2);
          expected_sd = sqrt(var_);
          se_mean = expected_sd / sqrt(n);
          se_sd = expected_sd / sqrt(2 * (n - 1));
        } else if (dist_type == 2) {
          real mean_d = exp(loc);
          real shape = phi;
          real scale_param = mean_d / shape;
          expected_mean = mean_d;
          expected_sd = sqrt(mean_d * scale_param);
          se_mean = expected_sd / sqrt(n);
          se_sd = expected_sd / sqrt(2 * (n - 1));
        } else if (dist_type == 3) {
          real scale = exp(loc);
          real shape = phi;
          expected_mean = scale * tgamma(1 + 1.0 / shape);
          real var_ = scale^2 * (tgamma(1 + 2.0 / shape) - pow(tgamma(1 + 1.0 / shape), 2));
          expected_sd = sqrt(var_);
          se_mean = expected_sd / sqrt(n);
          se_sd = expected_sd / sqrt(2 * (n - 1));
        }
        
        log_lik[idx]     = normal_lpdf(obs_stat1[d] | expected_mean, expected_sd / sqrt(n));
        log_lik[idx + 1] = normal_lpdf(obs_stat2[d] | expected_sd,   expected_sd / sqrt(2 * (n - 1)));
        log_lik[idx + 2] = 0;  // placeholder
        
      } else if (summary_type[d] == 4) {  // raw frequency table
        // Sum log-likelihoods over all individuals, using the frequency table.
        // Stored as a single scalar in log_lik[idx]; slots idx+1 and idx+2 are 0 (unused).
        real ll_type4 = 0;
        int s = freq_start[d];
        int len = freq_len[d];
        for (i in s:(s + len - 1)) {
          ll_type4 += freq_count[i] * dist_logpdf_fun(freq_value[i], dist_type, loc, phi);
        }
        log_lik[idx]     = ll_type4;
        log_lik[idx + 1] = 0;  // unused
        log_lik[idx + 2] = 0;  // unused

      } else if (summary_type[d] == 5) {  // interval-censored frequency table
        real ll_type5 = 0;
        int s = freq_start[d];
        int len = freq_len[d];
        for (i in s:(s + len - 1)) {
          if (freq_lower[i] == freq_upper[i]) {
            ll_type5 += freq_count[i] * dist_logpdf_fun(freq_lower[i], dist_type, loc, phi);
          } else {
            real cdf_u = dist_cdf_fun(freq_upper[i], dist_type, loc, phi);
            real cdf_l = dist_cdf_fun(freq_lower[i], dist_type, loc, phi);
            ll_type5 += freq_count[i] * log(cdf_u - cdf_l);
          }
        }
        log_lik[idx]     = ll_type5;
        log_lik[idx + 1] = 0;  // unused
        log_lik[idx + 2] = 0;  // unused
      }
      
      idx += 3;
    }
  }
}
