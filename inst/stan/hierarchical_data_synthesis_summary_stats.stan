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
    real log_f = dist_logpdf_fun(x, dist_type, loc, phi);  // CHANGED: get log directly
  
    // Improved numerical stability with tighter bounds
    F = fmax(fmin(F, 0.99999), 0.00001);  // CHANGED: tighter bounds
  
    // Use log_f directly instead of log(exp(log_f))
    real log_dens = lchoose(n, k) + log_f + (k - 1) * log(F) + (n - k) * log1m(F);
  
    return log_dens;
  }

  // Helper function to compute gamma quantile approximation
  real gamma_quantile_approx(real p, real shape, real scale) {
    // Wilson-Hilferty approximation for gamma quantiles
    real z = inv_Phi(p);  // standard normal quantile
    if (shape > 1) {
      return shape * scale * pow(1 - 1.0/(9*shape) + z/(3*sqrt(shape)), 3);
    } else {
      // Fallback for small shape
      return shape * scale * (1 + z / sqrt(shape));
  }
}
}

data {
  int<lower=1> n_datasets;              // Number of datasets
  array[n_datasets] int<lower=1> n_obs; // Sample sizes for each dataset
  array[n_datasets] int<lower=1,upper=3> summary_type; // 1=median+range, 2=median+IQR, 3=mean+sd
  int<lower=1,upper=3> dist_type;       // 1=lognormal, 2=gamma, 3=weibull
  
  // Observed summaries - organized by dataset
  array[n_datasets] real<lower=0> obs_stat1;     // median or mean
  array[n_datasets] real<lower=0> obs_stat2;     // min, q25, or sd
  array[n_datasets] real<lower=0> obs_stat3;     // max, q75, or placeholder
  
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
  // Precompute constants for quantiles
  real z_q25 = 0.6745;      // ~Phi^(-1)(0.25)
  real z_q75 = 0.6745;      // ~Phi^(-1)(0.75)
  real z_q90 = 1.28155;     // ~Phi^(-1)(0.90)
  real z_q95 = 1.64485;     // ~Phi^(-1)(0.95)
  
  // Validate data
  for (d in 1:n_datasets) {
    if (summary_type[d] == 1) {  // median + range
      if (obs_stat2[d] > obs_stat1[d] || obs_stat1[d] > obs_stat3[d]) {
        reject("For summary_type=1, must have min <= median <= max");
      }
    } else if (summary_type[d] == 2) {  // median + IQR
      if (obs_stat2[d] > obs_stat1[d] || obs_stat1[d] > obs_stat3[d]) {
        reject("For summary_type=2, must have q25 <= median <= q75");
      }
    }
  }
}

parameters {
  real mu0;                    // Population mean (location)
  real log_tau;                // Log of between-study SD
  real log_phi;                // Log of distribution-specific parameter
  vector[n_datasets] loc_d_raw; // Non-centered parameterization
}

transformed parameters {
  real<lower=0> tau = exp(log_tau);
  real<lower=0> phi = exp(log_phi);
  vector[n_datasets] loc_d = mu0 + tau * loc_d_raw;
}

model {
  // Priors (using configured or default values)
  mu0 ~ normal(mu0_mean, mu0_sd);
  log_tau ~ normal(log_tau_mean, log_tau_sd);
  log_phi ~ normal(log_phi_mean, log_phi_sd);
  
  // Non-centered parameterization
  loc_d_raw ~ std_normal();
  
  // Likelihood for each dataset
  for (d in 1:n_datasets) {
    real loc = loc_d[d];
    int n = n_obs[d];
    
    if (summary_type[d] == 1) {  // median + range (min, max)
      // Median is the middle order statistic
      int k_median = (n + 1) %/% 2;
      target += order_stat_logpdf_fun(obs_stat1[d], n, k_median, dist_type, loc, phi);
      
      // Min is 1st order statistic
      target += order_stat_logpdf_fun(obs_stat2[d], n, 1, dist_type, loc, phi);
      
      // Max is n-th order statistic
      target += order_stat_logpdf_fun(obs_stat3[d], n, n, dist_type, loc, phi);
    }
    
    else if (summary_type[d] == 2) {  // median + IQR (q25, q75)
      // Median
      int k_median = (n + 1) %/% 2;
      target += order_stat_logpdf_fun(obs_stat1[d], n, k_median, dist_type, loc, phi);
      
      // Q25 is approximately the n/4-th order statistic (at least 1)
      int k_q25 = (n + 1) %/% 4;
      if (k_q25 < 1) k_q25 = 1;
      target += order_stat_logpdf_fun(obs_stat2[d], n, k_q25, dist_type, loc, phi);

      // Q75 is approximately the 3n/4-th order statistic (at least k_q25 + 1)
      int k_q75 = (3 * (n + 1)) %/% 4;
      if (k_q75 <= k_q25) k_q75 = k_q25 + 1;
      if (k_q75 > n) k_q75 = n;  // optional safety
      target += order_stat_logpdf_fun(obs_stat3[d], n, k_q75, dist_type, loc, phi);
    }
    
    else if (summary_type[d] == 3) {  // mean + sd
      // For mean and SD, we use moment-based likelihood (normal approximation)
      
      real expected_mean;
      real expected_sd;
      real se_mean;
      real se_sd;
      
      if (dist_type == 1) {  // lognormal
        expected_mean = exp(loc + phi^2 / 2);
        real var_ = (exp(phi^2) - 1) * exp(2 * loc + phi^2);
        expected_sd = sqrt(var_);
        se_mean = expected_sd / sqrt(n);
        se_sd = expected_sd / sqrt(2 * (n - 1));  // approximate
        
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
      
      // Likelihood for observed mean and sd
      obs_stat1[d] ~ normal(expected_mean, se_mean);
      obs_stat2[d] ~ normal(expected_sd, se_sd);
      // obs_stat3[d] is ignored (placeholder)
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
  
  {
    int L = 100;
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
        means[l] = exp(loc_sample + phi^2 / 2);
        medians[l] = exp(loc_sample);
        q25s[l] = exp(loc_sample - 0.6745 * phi);
        q75s[l] = exp(loc_sample + 0.6745 * phi);
        q90s[l] = exp(loc_sample + 1.28155 * phi);
        q95s[l] = exp(loc_sample + 1.64485 * phi);  
        sds[l] = sqrt((exp(phi^2) - 1) * exp(2 * loc_sample + phi^2));
        
      } else if (dist_type == 2) {  // gamma
        real mean_d = exp(loc_sample);
        real shape = phi;
        real scale_param = mean_d / shape;
        
        means[l] = mean_d;
        sds[l] = sqrt(mean_d * scale_param);
        
        medians[l] = gamma_quantile_approx(0.5, shape, scale_param);
        q25s[l] = gamma_quantile_approx(0.25, shape, scale_param);
        q75s[l] = gamma_quantile_approx(0.75, shape, scale_param);
        q90s[l] = gamma_quantile_approx(0.90, shape, scale_param);
        q95s[l] = gamma_quantile_approx(0.95, shape, scale_param); 
          
      } else if (dist_type == 3) {  // weibull
        real scale = exp(loc_sample);
        real shape = phi;
        
        means[l] = scale * tgamma(1 + 1.0 / shape);
        medians[l] = scale * pow(log(2), 1.0 / shape);
        q25s[l] = scale * pow(log(4.0 / 3.0), 1.0 / shape);
        q75s[l] = scale * pow(log(4.0), 1.0 / shape);
        q90s[l] = scale * pow(log(10.0), 1.0 / shape);
        q95s[l] = scale * pow(log(20.0), 1.0 / shape);
        
        real var_weib = scale^2 * (tgamma(1 + 2.0 / shape) - pow(tgamma(1 + 1.0 / shape), 2));
        sds[l] = sqrt(var_weib);
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
        if (k_q75 > n) k_q75 = n;  // optional
        
        log_lik[idx]     = order_stat_logpdf_fun(obs_stat1[d], n, k_median, dist_type, loc, phi);
        log_lik[idx + 1] = order_stat_logpdf_fun(obs_stat2[d], n, k_q25,    dist_type, loc, phi);
        log_lik[idx + 2] = order_stat_logpdf_fun(obs_stat3[d], n, k_q75,    dist_type, loc, phi);
        
      } else if (summary_type[d] == 3) {  // mean + sd
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
          expected_sd = sqrt(mean_d * scale_param);
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
        
        log_lik[idx]     = normal_lpdf(obs_stat1[d] | expected_mean, expected_sd / sqrt(n));
        log_lik[idx + 1] = normal_lpdf(obs_stat2[d] | expected_sd,   expected_sd / sqrt(2 * (n - 1)));
        log_lik[idx + 2] = 0;  // placeholder
      }
      
      idx += 3;
    }
  }
}
