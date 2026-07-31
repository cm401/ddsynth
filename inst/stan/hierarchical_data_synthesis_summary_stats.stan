// hierarchical_data_synthesis_summary_stats.stan
functions {
  // CDF for each distribution
  real dist_cdf_fun(real x, int dist_type, real loc, real phi, real kappa) {
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
    } else if (dist_type == 4) {  // burr XII: lambda=exp(loc), c=phi, k=kappa
      real u = log(x) - loc;
      return -expm1(-kappa * log1p_exp(phi * u));
    } else if (dist_type == 5) {  // generalised gamma (Prentice): mu=loc, sigma=phi, Q=kappa
      real gamma_shape = 1.0 / (kappa * kappa);
      real w = (log(x) - loc) / phi;
      return gamma_cdf(gamma_shape * exp(kappa * w) | gamma_shape, 1);
    }
    return 0;
  }

  // "log pdf" for each distribution
  real dist_logpdf_fun(real x, int dist_type, real loc, real phi, real kappa) {
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
    } else if (dist_type == 4) {  // burr XII: lambda=exp(loc), c=phi, k=kappa
      real u = log(x) - loc;
      real log_term = log1p_exp(phi * u);  // log(1 + (x/lambda)^c)
      return log(phi) + log(kappa) + (phi - 1) * u - loc - (kappa + 1) * log_term;
    } else if (dist_type == 5) {  // generalised gamma (Prentice): mu=loc, sigma=phi, Q=kappa
      real gamma_shape = 1.0 / (kappa * kappa);
      real w = (log(x) - loc) / phi;
      return log(kappa) - log(phi) - log(x)
             + gamma_shape * log(gamma_shape)
             + gamma_shape * kappa * w
             - gamma_shape * exp(kappa * w)
             - lgamma(gamma_shape);
    }
    return 0;
  }

  // Log-scale CDF for each distribution.
  // Used in order statistic and interval-censored likelihoods to avoid
  // probability-scale clipping and catastrophic cancellation.
  real dist_log_cdf_fun(real x, int dist_type, real loc, real phi, real kappa) {
    if (dist_type == 1) {  // lognormal
      return lognormal_lcdf(x | loc, phi);
    } else if (dist_type == 2) {  // gamma
      real mean_d = exp(loc);
      real rate   = phi / mean_d;
      return gamma_lcdf(x | phi, rate);
    } else if (dist_type == 3) {  // weibull
      real scale = exp(loc);
      return weibull_lcdf(x | phi, scale);
    } else if (dist_type == 4) {  // burr XII: lambda=exp(loc), c=phi, k=kappa
      real u = log(x) - loc;
      return log1m_exp(-kappa * log1p_exp(phi * u));
    } else if (dist_type == 5) {  // generalised gamma (Prentice): mu=loc, sigma=phi, Q=kappa
      real gamma_shape = 1.0 / (kappa * kappa);
      real w = (log(x) - loc) / phi;
      return gamma_lcdf(gamma_shape * exp(kappa * w) | gamma_shape, 1);
    }
    return negative_infinity();
  }

  // Log likelihood for order statistic:
  // k-th order statistic out of n observations.
  //
  // Uses log-scale CDF (_lcdf) and CCDF (_lccdf) directly rather than
  // computing the CDF on the probability scale and clipping. Clipping
  // (fmax/fmin) is not differentiable at its boundaries, causing zero
  // gradients whenever F ~ 0 or F ~ 1 — a common occurrence for tightly
  // concentrated distributions such as gamma with large shape.
  //
  // Safeguard: terms with coefficient 0 are omitted explicitly to avoid
  // 0 * (-Inf) = NaN in Stan's autodiff:
  //   k == 1 (minimum): coefficient of log_F  is (k-1) = 0 — term omitted
  //   k == n (maximum): coefficient of log_1mF is (n-k) = 0 — term omitted
  real order_stat_logpdf_fun(real x, int n, int k, int dist_type, real loc, real phi, real kappa) {
    real log_f;
    real log_F   = 0.0;  // log CDF  = log P(X <= x); default safe (unused when k==1)
    real log_1mF = 0.0;  // log CCDF = log P(X >  x); default safe (unused when k==n)

    if (dist_type == 1) {  // lognormal
      log_f   = lognormal_lpdf(x  | loc, phi);
      log_F   = lognormal_lcdf(x  | loc, phi);
      log_1mF = lognormal_lccdf(x | loc, phi);

    } else if (dist_type == 2) {  // gamma
      real mean_d = exp(loc);
      real rate   = phi / mean_d;
      log_f   = gamma_lpdf(x  | phi, rate);
      log_F   = gamma_lcdf(x  | phi, rate);
      log_1mF = gamma_lccdf(x | phi, rate);

    } else if (dist_type == 3) {  // weibull
      real scale = exp(loc);
      log_f   = weibull_lpdf(x  | phi, scale);
      log_F   = weibull_lcdf(x  | phi, scale);
      log_1mF = weibull_lccdf(x | phi, scale);

    } else if (dist_type == 4) {  // burr XII: lambda=exp(loc), c=phi, k=kappa
      real u        = log(x) - loc;
      real log_term = log1p_exp(phi * u);  // log(1 + (x/lambda)^c)
      log_f   = log(phi) + log(kappa) + (phi - 1) * u - loc - (kappa + 1) * log_term;
      log_1mF = -kappa * log_term;                    // log CCDF = log((1+(x/lam)^c)^(-k))
      log_F   = log1m_exp(log_1mF);                   // log CDF  = log(1 - CCDF)

    } else if (dist_type == 5) {  // generalised gamma (Prentice): mu=loc, sigma=phi, Q=kappa
      real gamma_shape = 1.0 / (kappa * kappa);
      real w    = (log(x) - loc) / phi;
      real arg  = gamma_shape * exp(kappa * w);  // gamma_shape * exp(Q*w)
      log_f   = log(kappa) - log(phi) - log(x)
                + gamma_shape * log(gamma_shape)
                + gamma_shape * kappa * w
                - arg
                - lgamma(gamma_shape);
      if (k > 1) log_F   = gamma_lcdf(arg  | gamma_shape, 1);
      if (k < n) log_1mF = gamma_lccdf(arg | gamma_shape, 1);
    }

    real log_dens = lchoose(n, k) + log_f;
    if (k > 1) log_dens += (k - 1) * log_F;
    if (k < n) log_dens += (n - k) * log_1mF;
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

  // Log of integral_{ex_l}^{ex_r} F_D(T - e) de — the truncation denominator
  // for type 7 (right-truncated, doubly censored).  Uses 7-point GL quadrature.
  // Degenerates to log F_D(T - ex_l) when ex_l == ex_r (point exposure).
  real trunc_denom_log(real ex_l, real ex_r, real T,
                       int dist_type, real loc, real phi, real kappa) {
    if (ex_l == ex_r) {
      return dist_log_cdf_fun(T - ex_l, dist_type, loc, phi, kappa);
    }
    array[7] real t = {-0.9491079123427585, -0.7415311855993945,
                       -0.4058451513773832,  0.0,
                        0.4058451513773832,  0.7415311855993945,
                        0.9491079123427585};
    array[7] real w = { 0.1294849661688697,  0.2797053914892767,
                        0.3818300505051189,  0.4179591836734694,
                        0.3818300505051189,  0.2797053914892767,
                        0.1294849661688697};
    real mid  = 0.5 * (ex_r + ex_l);
    real half = 0.5 * (ex_r - ex_l);
    array[7] real log_terms;
    for (k in 1:7) {
      real e_k      = mid + half * t[k];
      real max_del  = T - e_k;
      log_terms[k]  = log(w[k]) + dist_log_cdf_fun(max_del, dist_type, loc, phi, kappa);
    }
    return log_sum_exp(log_terms);
  }

  // Log of integral_{ex_l}^{ex_r} [1 - F_D(T - e)] de — right-censored contribution
  // for type 7 (onset not yet observed by T).  Uses 7-point GL quadrature.
  real right_censor_log_lik(real ex_l, real ex_r, real T,
                             int dist_type, real loc, real phi, real kappa) {
    if (ex_l == ex_r) {
      return log1m_exp(dist_log_cdf_fun(T - ex_l, dist_type, loc, phi, kappa));
    }
    array[7] real t = {-0.9491079123427585, -0.7415311855993945,
                       -0.4058451513773832,  0.0,
                        0.4058451513773832,  0.7415311855993945,
                        0.9491079123427585};
    array[7] real w = { 0.1294849661688697,  0.2797053914892767,
                        0.3818300505051189,  0.4179591836734694,
                        0.3818300505051189,  0.2797053914892767,
                        0.1294849661688697};
    real mid  = 0.5 * (ex_r + ex_l);
    real half = 0.5 * (ex_r - ex_l);
    array[7] real log_terms;
    for (k in 1:7) {
      real e_k     = mid + half * t[k];
      real max_del = T - e_k;
      log_terms[k] = log(w[k]) + log1m_exp(dist_log_cdf_fun(max_del, dist_type, loc, phi, kappa));
    }
    return log_sum_exp(log_terms);
  }

  // Double-censored log-likelihood for one observation (summary type 6).
  //
  // Computes log P(E in [ex_l, ex_r], O in [ev_l, ev_r]) proportional to theta,
  // where D = O - E ~ f_D(d; theta), assuming E uniform on [ex_l, ex_r] and
  // O uniform on [ev_l, ev_r].
  //
  // The normalising constant (ex_r - ex_l)(ev_r - ev_l) is constant in theta
  // and is dropped, consistent with how type 5 drops 1/(ex_r - ex_l).
  //
  // Four cases:
  //   Both endpoints point-observed          -> log f_D(ev_l - ex_l)
  //   Point exposure, interval event         -> log[F(ev_r - ex_l) - F(ev_l - ex_l)]
  //   Interval exposure, point event (= type 5) -> log[F(ev_l - ex_l) - F(ev_l - ex_r)]
  //   Both interval-censored (general)       -> 7-point Gauss-Legendre quadrature
  real dc_log_lik(real ex_l, real ex_r, real ev_l, real ev_r,
                  int dist_type, real loc, real phi, real kappa) {
    // Case 1: both endpoints point-observed
    if (ex_l == ex_r && ev_l == ev_r) {
      return dist_logpdf_fun(ev_l - ex_l, dist_type, loc, phi, kappa);
    }
    // Case 2: point exposure, interval event
    if (ex_l == ex_r) {
      return log_diff_exp(
        dist_log_cdf_fun(ev_r - ex_l, dist_type, loc, phi, kappa),
        dist_log_cdf_fun(ev_l - ex_l, dist_type, loc, phi, kappa)
      );
    }
    // Case 3: interval exposure, point event — reduces to type 5
    if (ev_l == ev_r) {
      return log_diff_exp(
        dist_log_cdf_fun(ev_l - ex_l, dist_type, loc, phi, kappa),
        dist_log_cdf_fun(ev_l - ex_r, dist_type, loc, phi, kappa)
      );
    }
    // Case 4: both interval-censored — 7-point Gauss-Legendre quadrature on [ex_l, ex_r].
    // Nodes and weights on [-1, 1] (Abramowitz & Stegun table 25.4).
    // The Jacobian factor (ex_r - ex_l)/2 is constant in theta and is dropped.
    array[7] real t = {-0.9491079123427585, -0.7415311855993945,
                       -0.4058451513773832,  0.0,
                        0.4058451513773832,  0.7415311855993945,
                        0.9491079123427585};
    array[7] real w = { 0.1294849661688697,  0.2797053914892767,
                        0.3818300505051189,  0.4179591836734694,
                        0.3818300505051189,  0.2797053914892767,
                        0.1294849661688697};
    real mid  = 0.5 * (ex_r + ex_l);
    real half = 0.5 * (ex_r - ex_l);
    array[7] real log_terms;
    for (k in 1:7) {
      real e_k = mid + half * t[k];
      log_terms[k] = log(w[k]) + log_diff_exp(
        dist_log_cdf_fun(ev_r - e_k, dist_type, loc, phi, kappa),
        dist_log_cdf_fun(ev_l - e_k, dist_type, loc, phi, kappa)
      );
    }
    return log_sum_exp(log_terms);
  }
}

data {
  int<lower=1> n_datasets;              // Number of datasets
  array[n_datasets] int<lower=1> n_obs; // Sample sizes for each dataset
  array[n_datasets] int<lower=1,upper=7> summary_type; // 1=median+range, 2=median+IQR, 3=mean+sd, 4=raw freq table, 5=interval-censored freq table, 6=double interval-censored freq table, 7=double interval-censored with right truncation/censoring
  int<lower=1,upper=5> dist_type;       // 1=lognormal, 2=gamma, 3=weibull, 4=burr XII, 5=gen. gamma

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

  // Event window bounds for summary_type == 6 and 7 (double interval-censored).
  // For all other summary types these arrays are populated with zeros.
  // For types 6 and 7, freq_lower / freq_upper carry the exposure window (expo_lower, expo_upper).
  array[n_freq_total] real<lower=0> event_lower;   // lower bound of event window; 0 for right-censored (type 7)
  array[n_freq_total] real<lower=0> event_upper;   // upper bound of event window; 0 for right-censored (type 7)

  // Type 7: right truncation / right censoring.
  // event_observed[i] = 1 if onset was recorded, 0 if right-censored (onset not yet seen).
  // truncation_time[d] = analysis date T in days from the reference date; ignored for other types.
  array[n_freq_total] int<lower=0, upper=1> event_observed;  // 1=onset seen, 0=right-censored
  array[n_datasets]   real<lower=0>         truncation_time; // analysis date T per dataset

  // Prior hyperparameters for mu0 ~ normal(mu0_mean, mu0_sd)
  real mu0_mean;
  real<lower=0> mu0_sd;

  // Prior hyperparameters for log_tau ~ normal(log_tau_mean, log_tau_sd)
  real log_tau_mean;
  real<lower=0> log_tau_sd;

  // Prior hyperparameters for log_phi ~ normal(log_phi_mean, log_phi_sd)
  real log_phi_mean;
  real<lower=0> log_phi_sd;

  // Prior hyperparameters for log_kappa ~ normal(log_kappa_mean, log_kappa_sd)
  // kappa = exp(log_kappa) > 0; used by dist_type 4 (Burr XII k) and 5 (GG Q).
  // For dist_type 1-3 supply wide uninformative priors (e.g. mean=0, sd=1);
  // kappa will be sampled from its prior but does not enter the likelihood.
  real log_kappa_mean;
  real<lower=0> log_kappa_sd;
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
    } else if (summary_type[d] == 6) {
      if (freq_len[d] == 0) {
        reject("For summary_type=6, freq_len must be > 0");
      }
    } else if (summary_type[d] == 7) {
      if (freq_len[d] == 0) {
        reject("For summary_type=7, freq_len must be > 0");
      }
      if (truncation_time[d] <= 0) {
        reject("For summary_type=7, truncation_time must be > 0");
      }
    }
  }
}

parameters {
  real mu0;                        // Population mean (location)
  real log_tau;                    // Log of between-study SD
  real log_phi;                    // Log of distribution-specific shape/scale parameter
  real log_kappa;                  // Log of 3rd distribution parameter (Burr XII k; GG Q)
  vector[n_datasets] loc_d_raw;    // Non-centered parameterization
}

transformed parameters {
  real<lower=0> tau   = exp(log_tau);
  real<lower=0> phi   = exp(log_phi);
  real<lower=0> kappa = exp(log_kappa);
  vector[n_datasets] loc_d = mu0 + tau * loc_d_raw;
}

model {
  mu0       ~ normal(mu0_mean, mu0_sd);
  log_tau   ~ normal(log_tau_mean, log_tau_sd);
  log_phi   ~ normal(log_phi_mean, log_phi_sd);
  log_kappa ~ normal(log_kappa_mean, log_kappa_sd);

  loc_d_raw ~ std_normal();

  for (d in 1:n_datasets) {
    real loc = loc_d[d];
    int n = n_obs[d];

    if (summary_type[d] == 1) {  // median + range (min, max)
      int k_median = (n + 1) %/% 2;
      target += order_stat_logpdf_fun(obs_stat1[d], n, k_median, dist_type, loc, phi, kappa);
      target += order_stat_logpdf_fun(obs_stat2[d], n, 1, dist_type, loc, phi, kappa);
      target += order_stat_logpdf_fun(obs_stat3[d], n, n, dist_type, loc, phi, kappa);
    }

    else if (summary_type[d] == 2) {  // median + IQR (q25, q75)
      int k_median = (n + 1) %/% 2;
      target += order_stat_logpdf_fun(obs_stat1[d], n, k_median, dist_type, loc, phi, kappa);

      int k_q25 = (n + 1) %/% 4;
      if (k_q25 < 1) k_q25 = 1;
      target += order_stat_logpdf_fun(obs_stat2[d], n, k_q25, dist_type, loc, phi, kappa);

      int k_q75 = (3 * (n + 1)) %/% 4;
      if (k_q75 <= k_q25) k_q75 = k_q25 + 1;
      if (k_q75 > n) k_q75 = n;
      target += order_stat_logpdf_fun(obs_stat3[d], n, k_q75, dist_type, loc, phi, kappa);
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

      } else if (dist_type == 4) {  // burr XII: lambda=exp(loc), c=phi, k=kappa
        // E[X^r] = lambda^r * k * B(k - r/c, 1 + r/c), requires k*c > r
        // Mean requires kappa*phi > 1; variance requires kappa*phi > 2.
        // Hard barrier: reject this region of parameter space.
        if (kappa * phi <= 2.0) {
          target += negative_infinity();
        } else {
          real lam = exp(loc);
          expected_mean = lam * kappa * exp(lbeta(kappa - 1.0/phi, 1.0 + 1.0/phi));
          real e2       = lam^2 * kappa * exp(lbeta(kappa - 2.0/phi, 1.0 + 2.0/phi));
          expected_sd   = sqrt(fabs(e2 - expected_mean^2));
          se_mean = expected_sd / sqrt(n);
          se_sd   = expected_sd / sqrt(2 * (n - 1));
          obs_stat1[d] ~ normal(expected_mean, se_mean);
          obs_stat2[d] ~ normal(expected_sd,   se_sd);
        }

      } else if (dist_type == 5) {  // generalised gamma (Prentice): mu=loc, sigma=phi, Q=kappa
        // E[T^r] = exp(r*loc) * kappa^(2*r*phi/kappa)
        //          * Gamma(gamma_shape + r*phi/kappa) / Gamma(gamma_shape)
        // where gamma_shape = 1/kappa^2
        real gamma_shape = 1.0 / (kappa * kappa);
        real log_ET  = loc + 2.0*phi/kappa * log(kappa)
                       + lgamma(gamma_shape + phi/kappa) - lgamma(gamma_shape);
        real log_ET2 = 2.0*loc + 4.0*phi/kappa * log(kappa)
                       + lgamma(gamma_shape + 2.0*phi/kappa) - lgamma(gamma_shape);
        expected_mean = exp(log_ET);
        expected_sd   = sqrt(fabs(exp(log_ET2) - expected_mean^2));
        se_mean = expected_sd / sqrt(n);
        se_sd   = expected_sd / sqrt(2 * (n - 1));
        obs_stat1[d] ~ normal(expected_mean, se_mean);
        obs_stat2[d] ~ normal(expected_sd,   se_sd);
      }

      // For dist_type 1-3 the likelihood is added below (after the if-else chain).
      // For dist_type 4-5 it is added inside their own branches above.
      if (dist_type <= 3) {
        obs_stat1[d] ~ normal(expected_mean, se_mean);
        obs_stat2[d] ~ normal(expected_sd, se_sd);
      }
    }

    else if (summary_type[d] == 4) {  // raw frequency table
      // Direct likelihood: for each distinct observed value, add count * log_pdf(value).
      // This is equivalent to fitting the distribution directly to all individual observations,
      // but using the compressed frequency-table representation.
      int s = freq_start[d];
      int len = freq_len[d];
      for (i in s:(s + len - 1)) {
        target += freq_count[i] * dist_logpdf_fun(freq_value[i], dist_type, loc, phi, kappa);
      }
    }

    else if (summary_type[d] == 5) {  // interval-censored frequency table
      // Likelihood: count * log[ F(upper) - F(lower) ] for each interval.
      // Uses log_diff_exp(log_F_upper, log_F_lower) for numerical stability —
      // avoids catastrophic cancellation when the two CDF values are close.
      // When lower == upper (point observation), falls back to count * log_pdf
      // since log_diff_exp(a, a) = -Inf.
      int s = freq_start[d];
      int len = freq_len[d];
      for (i in s:(s + len - 1)) {
        if (freq_lower[i] == freq_upper[i]) {
          target += freq_count[i] * dist_logpdf_fun(freq_lower[i], dist_type, loc, phi, kappa);
        } else {
          real log_cdf_u = dist_log_cdf_fun(freq_upper[i], dist_type, loc, phi, kappa);
          real log_cdf_l = dist_log_cdf_fun(freq_lower[i], dist_type, loc, phi, kappa);
          target += freq_count[i] * log_diff_exp(log_cdf_u, log_cdf_l);
        }
      }
    }

    else if (summary_type[d] == 6) {  // double interval-censored frequency table
      // Likelihood: count * dc_log_lik(expo_lower, expo_upper, event_lower, event_upper).
      // freq_lower / freq_upper carry the exposure window bounds for type 6.
      // Uses 7-point GL quadrature; degenerates correctly to type 5 when
      // event_lower[i] == event_upper[i] (point event observation).
      int s = freq_start[d];
      int len = freq_len[d];
      for (i in s:(s + len - 1)) {
        target += freq_count[i] * dc_log_lik(
          freq_lower[i], freq_upper[i],
          event_lower[i], event_upper[i],
          dist_type, loc, phi, kappa
        );
      }
    }

    else if (summary_type[d] == 7) {  // double interval-censored + right truncation/censoring
      // For observed individuals (event_observed[i] == 1):
      //   log L = dc_log_lik(...) - trunc_denom_log(...)
      //   i.e. numerator = integral F(O_R-e) - F(O_L-e) de over exposure window
      //        denominator = integral F(T-e) de over exposure window
      // For right-censored individuals (event_observed[i] == 0):
      //   log L = right_censor_log_lik(...)
      //   i.e. integral [1 - F(T-e)] de over exposure window
      int s = freq_start[d];
      int len = freq_len[d];
      real T = truncation_time[d];
      for (i in s:(s + len - 1)) {
        if (event_observed[i] == 1) {
          target += freq_count[i] * (
            dc_log_lik(freq_lower[i], freq_upper[i],
                       event_lower[i], event_upper[i],
                       dist_type, loc, phi, kappa)
            - trunc_denom_log(freq_lower[i], freq_upper[i], T,
                              dist_type, loc, phi, kappa)
          );
        } else {
          target += freq_count[i] * right_censor_log_lik(
            freq_lower[i], freq_upper[i], T,
            dist_type, loc, phi, kappa
          );
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
  //   pred_* are computed at mean(loc_d) directly to avoid tau^2/2 inflation.
  //   See: Higgins & Thompson (2002) doi:10.1002/sim.1186
  //        Gelman (2006) doi:10.1214/06-BA117A
  //        Rover et al. (2021) doi:10.1002/jrsm.1475
  // - When n_datasets >= 5, tau is identifiable; sample from Normal(mu0, tau)
  //   to include between-study heterogeneity. L=2000 for MC stability.
  {
    if (n_datasets < 5) {
      real loc_pred = mean(loc_d);

      if (dist_type == 1) {  // lognormal
        pred_mean   = exp(loc_pred + phi^2 / 2);
        pred_median = exp(loc_pred);
        pred_q25    = exp(loc_pred - 0.6745  * phi);
        pred_q75    = exp(loc_pred + 0.6745  * phi);
        pred_q90    = exp(loc_pred + 1.28155 * phi);
        pred_q95    = exp(loc_pred + 1.64485 * phi);
        pred_sd     = sqrt((exp(phi^2) - 1) * exp(2 * loc_pred + phi^2));

      } else if (dist_type == 2) {  // gamma
        real mean_d      = exp(loc_pred);
        real scale_param = mean_d / phi;
        pred_mean   = mean_d;
        pred_sd     = sqrt(mean_d * scale_param);
        pred_median = gamma_quantile_approx(0.5,  phi, scale_param);
        pred_q25    = gamma_quantile_approx(0.25, phi, scale_param);
        pred_q75    = gamma_quantile_approx(0.75, phi, scale_param);
        pred_q90    = gamma_quantile_approx(0.90, phi, scale_param);
        pred_q95    = gamma_quantile_approx(0.95, phi, scale_param);

      } else if (dist_type == 3) {  // weibull
        real scale  = exp(loc_pred);
        pred_mean   = scale * tgamma(1 + 1.0 / phi);
        pred_median = scale * pow(log(2),         1.0 / phi);
        pred_q25    = scale * pow(log(4.0 / 3.0), 1.0 / phi);
        pred_q75    = scale * pow(log(4.0),       1.0 / phi);
        pred_q90    = scale * pow(log(10.0),      1.0 / phi);
        pred_q95    = scale * pow(log(20.0),      1.0 / phi);
        pred_sd     = sqrt(scale^2 * (tgamma(1 + 2.0/phi) - pow(tgamma(1 + 1.0/phi), 2)));

      } else if (dist_type == 4) {  // burr XII: closed-form quantiles
        // Q(p) = lambda * ((1-p)^(-1/k) - 1)^(1/c)
        real lam    = exp(loc_pred);
        pred_median = lam * pow(pow(0.5,  -1.0/kappa) - 1.0, 1.0/phi);
        pred_q25    = lam * pow(pow(0.75, -1.0/kappa) - 1.0, 1.0/phi);
        pred_q75    = lam * pow(pow(0.25, -1.0/kappa) - 1.0, 1.0/phi);
        pred_q90    = lam * pow(pow(0.10, -1.0/kappa) - 1.0, 1.0/phi);
        pred_q95    = lam * pow(pow(0.05, -1.0/kappa) - 1.0, 1.0/phi);
        // E[X^r] = lambda^r * k * B(k-r/c, 1+r/c), requires k*c > r
        if (kappa * phi > 1.0) {
          pred_mean = lam * kappa * exp(lbeta(kappa - 1.0/phi, 1.0 + 1.0/phi));
        } else {
          pred_mean = positive_infinity();
        }
        if (kappa * phi > 2.0) {
          real e2 = lam^2 * kappa * exp(lbeta(kappa - 2.0/phi, 1.0 + 2.0/phi));
          pred_sd = sqrt(fabs(e2 - pred_mean^2));
        } else {
          pred_sd = positive_infinity();
        }

      } else if (dist_type == 5) {  // generalised gamma: Monte Carlo via gamma_rng
        // Sample T = exp(mu + sigma/Q * log(Q^2 * Y)), Y ~ Gamma(1/Q^2, 1)
        int L_gg = 200;
        vector[L_gg] gg_samples;
        real gs = 1.0 / (kappa * kappa);
        for (l in 1:L_gg) {
          real y = gamma_rng(gs, 1);
          gg_samples[l] = exp(loc_pred + phi / kappa * log(kappa * kappa * y));
        }
        gg_samples  = sort_asc(gg_samples);
        pred_mean   = mean(gg_samples);
        pred_sd     = sd(gg_samples);
        pred_median = gg_samples[100];
        pred_q25    = gg_samples[50];
        pred_q75    = gg_samples[150];
        pred_q90    = gg_samples[180];
        pred_q95    = gg_samples[190];
      }

    } else {
      // n_datasets >= 5: tau identifiable; include between-study heterogeneity
      // via Monte Carlo integration over Normal(mu0, tau). L=2000 for stability.
      int L = 2000;

      if (dist_type == 5) {
        // GG: sample one T per study draw to get the full predictive distribution
        // T = exp(mu + sigma/Q * log(Q^2 * Y)), Y ~ Gamma(1/Q^2, 1)
        int L_gg = 200;
        vector[L_gg] gg_t_samples;
        real gs = 1.0 / (kappa * kappa);
        for (l in 1:L_gg) {
          real loc_sample = normal_rng(mu0, tau);
          real y = gamma_rng(gs, 1);
          gg_t_samples[l] = exp(loc_sample + phi / kappa * log(kappa * kappa * y));
        }
        gg_t_samples = sort_asc(gg_t_samples);
        pred_mean   = mean(gg_t_samples);
        pred_sd     = sd(gg_t_samples);
        pred_median = gg_t_samples[100];
        pred_q25    = gg_t_samples[50];
        pred_q75    = gg_t_samples[150];
        pred_q90    = gg_t_samples[180];
        pred_q95    = gg_t_samples[190];

      } else {
        // dist_type 1-4: compute per-study analytic moments/quantiles, then average
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

          } else if (dist_type == 4) {  // burr XII: closed-form quantiles
            real lam_l  = exp(loc_sample);
            medians[l]  = lam_l * pow(pow(0.5,  -1.0/kappa) - 1.0, 1.0/phi);
            q25s[l]     = lam_l * pow(pow(0.75, -1.0/kappa) - 1.0, 1.0/phi);
            q75s[l]     = lam_l * pow(pow(0.25, -1.0/kappa) - 1.0, 1.0/phi);
            q90s[l]     = lam_l * pow(pow(0.10, -1.0/kappa) - 1.0, 1.0/phi);
            q95s[l]     = lam_l * pow(pow(0.05, -1.0/kappa) - 1.0, 1.0/phi);
            if (kappa * phi > 1.0) {
              means[l] = lam_l * kappa * exp(lbeta(kappa - 1.0/phi, 1.0 + 1.0/phi));
            } else {
              means[l] = positive_infinity();
            }
            if (kappa * phi > 2.0) {
              real e2_l = lam_l^2 * kappa * exp(lbeta(kappa - 2.0/phi, 1.0 + 2.0/phi));
              sds[l] = (kappa * phi > 1.0) ? sqrt(fabs(e2_l - means[l]^2)) : positive_infinity();
            } else {
              sds[l] = positive_infinity();
            }
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
  }

  // Log likelihood
  {
    int idx = 1;

    for (d in 1:n_datasets) {
      real loc = loc_d[d];
      int n = n_obs[d];

      if (summary_type[d] == 1) {  // median + range
        int k_median = (n + 1) %/% 2;
        log_lik[idx]     = order_stat_logpdf_fun(obs_stat1[d], n, k_median, dist_type, loc, phi, kappa);
        log_lik[idx + 1] = order_stat_logpdf_fun(obs_stat2[d], n, 1,        dist_type, loc, phi, kappa);
        log_lik[idx + 2] = order_stat_logpdf_fun(obs_stat3[d], n, n,        dist_type, loc, phi, kappa);

      } else if (summary_type[d] == 2) {  // median + IQR
        int k_median = (n + 1) %/% 2;
        int k_q25    = (n + 1) %/% 4;
        if (k_q25 < 1) k_q25 = 1;
        int k_q75    = (3 * (n + 1)) %/% 4;
        if (k_q75 <= k_q25) k_q75 = k_q25 + 1;
        if (k_q75 > n) k_q75 = n;

        log_lik[idx]     = order_stat_logpdf_fun(obs_stat1[d], n, k_median, dist_type, loc, phi, kappa);
        log_lik[idx + 1] = order_stat_logpdf_fun(obs_stat2[d], n, k_q25,   dist_type, loc, phi, kappa);
        log_lik[idx + 2] = order_stat_logpdf_fun(obs_stat3[d], n, k_q75,   dist_type, loc, phi, kappa);

      } else if (summary_type[d] == 3) {  // mean + sd
        real expected_mean;
        real expected_sd;

        if (dist_type == 1) {
          expected_mean = exp(loc + phi^2 / 2);
          real var_ = (exp(phi^2) - 1) * exp(2 * loc + phi^2);
          expected_sd = sqrt(var_);
        } else if (dist_type == 2) {
          real mean_d = exp(loc);
          real shape = phi;
          real scale_param = mean_d / shape;
          expected_mean = mean_d;
          expected_sd = sqrt(mean_d * scale_param);
        } else if (dist_type == 3) {
          real scale = exp(loc);
          real shape = phi;
          expected_mean = scale * tgamma(1 + 1.0 / shape);
          real var_ = scale^2 * (tgamma(1 + 2.0 / shape) - pow(tgamma(1 + 1.0 / shape), 2));
          expected_sd = sqrt(var_);
        } else if (dist_type == 4) {
          // Moments require kappa*phi > 2; model block already rejects that region,
          // so this fallback (negative_infinity) should never be reached in practice.
          if (kappa * phi > 2.0) {
            real lam = exp(loc);
            expected_mean = lam * kappa * exp(lbeta(kappa - 1.0/phi, 1.0 + 1.0/phi));
            real e2       = lam^2 * kappa * exp(lbeta(kappa - 2.0/phi, 1.0 + 2.0/phi));
            expected_sd   = sqrt(fabs(e2 - expected_mean^2));
          } else {
            log_lik[idx]     = negative_infinity();
            log_lik[idx + 1] = negative_infinity();
            log_lik[idx + 2] = 0;
          }
        } else if (dist_type == 5) {
          real gamma_shape = 1.0 / (kappa * kappa);
          real log_ET  = loc + 2.0*phi/kappa * log(kappa)
                         + lgamma(gamma_shape + phi/kappa) - lgamma(gamma_shape);
          real log_ET2 = 2.0*loc + 4.0*phi/kappa * log(kappa)
                         + lgamma(gamma_shape + 2.0*phi/kappa) - lgamma(gamma_shape);
          expected_mean = exp(log_ET);
          expected_sd   = sqrt(fabs(exp(log_ET2) - expected_mean^2));
        }

        // For dist_type 4 when moments don't exist the log_lik slots were
        // already assigned above; skip the normal_lpdf for that case.
        if (dist_type != 4 || kappa * phi > 2.0) {
          log_lik[idx]     = normal_lpdf(obs_stat1[d] | expected_mean, expected_sd / sqrt(n));
          log_lik[idx + 1] = normal_lpdf(obs_stat2[d] | expected_sd,   expected_sd / sqrt(2 * (n - 1)));
          log_lik[idx + 2] = 0;  // placeholder
        }

      } else if (summary_type[d] == 4) {  // raw frequency table
        // Sum log-likelihoods over all individuals, using the frequency table.
        // Stored as a single scalar in log_lik[idx]; slots idx+1 and idx+2 are 0 (unused).
        real ll_type4 = 0;
        int s = freq_start[d];
        int len = freq_len[d];
        for (i in s:(s + len - 1)) {
          ll_type4 += freq_count[i] * dist_logpdf_fun(freq_value[i], dist_type, loc, phi, kappa);
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
            ll_type5 += freq_count[i] * dist_logpdf_fun(freq_lower[i], dist_type, loc, phi, kappa);
          } else {
            real log_cdf_u = dist_log_cdf_fun(freq_upper[i], dist_type, loc, phi, kappa);
            real log_cdf_l = dist_log_cdf_fun(freq_lower[i], dist_type, loc, phi, kappa);
            ll_type5 += freq_count[i] * log_diff_exp(log_cdf_u, log_cdf_l);
          }
        }
        log_lik[idx]     = ll_type5;
        log_lik[idx + 1] = 0;  // unused
        log_lik[idx + 2] = 0;  // unused

      } else if (summary_type[d] == 6) {  // double interval-censored frequency table
        real ll_type6 = 0;
        int s = freq_start[d];
        int len = freq_len[d];
        for (i in s:(s + len - 1)) {
          ll_type6 += freq_count[i] * dc_log_lik(
            freq_lower[i], freq_upper[i],
            event_lower[i], event_upper[i],
            dist_type, loc, phi, kappa
          );
        }
        log_lik[idx]     = ll_type6;
        log_lik[idx + 1] = 0;  // unused
        log_lik[idx + 2] = 0;  // unused

      } else if (summary_type[d] == 7) {  // double interval-censored + right truncation/censoring
        real ll_type7 = 0;
        int s = freq_start[d];
        int len = freq_len[d];
        real T = truncation_time[d];
        for (i in s:(s + len - 1)) {
          if (event_observed[i] == 1) {
            ll_type7 += freq_count[i] * (
              dc_log_lik(freq_lower[i], freq_upper[i],
                         event_lower[i], event_upper[i],
                         dist_type, loc, phi, kappa)
              - trunc_denom_log(freq_lower[i], freq_upper[i], T,
                                dist_type, loc, phi, kappa)
            );
          } else {
            ll_type7 += freq_count[i] * right_censor_log_lik(
              freq_lower[i], freq_upper[i], T,
              dist_type, loc, phi, kappa
            );
          }
        }
        log_lik[idx]     = ll_type7;
        log_lik[idx + 1] = 0;  // unused
        log_lik[idx + 2] = 0;  // unused
      }

      idx += 3;
    }
  }
}
