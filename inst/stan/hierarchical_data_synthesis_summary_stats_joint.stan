// hierarchical_data_synthesis_summary_stats_joint.stan
//
// Identical to hierarchical_data_synthesis_summary_stats.stan except that
// summary types 1 (median+range) and 2 (median+IQR) use the joint density of
// three order statistics rather than the product of their marginals.
//
// The joint density of order statistics X_(i1) < X_(i2) < X_(i3) from n IID
// observations with CDF F and density f is:
//
//   f(x1,x2,x3) = n! / [(i1-1)! (i2-i1-1)! (i3-i2-1)! (n-i3)!]
//                 * F(x1)^{i1-1} * f(x1)
//                 * [F(x2)-F(x1)]^{i2-i1-1} * f(x2)
//                 * [F(x3)-F(x2)]^{i3-i2-1} * f(x3)
//                 * [1-F(x3)]^{n-i3}
//
// Gap terms with zero exponent are omitted to avoid 0*(-Inf) = NaN.
// Gaps are computed via log_diff_exp for numerical stability.

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

  // Log-scale CCDF (log survivor function) for each distribution.
  real dist_log_ccdf_fun(real x, int dist_type, real loc, real phi, real kappa) {
    if (dist_type == 1) {  // lognormal
      return lognormal_lccdf(x | loc, phi);
    } else if (dist_type == 2) {  // gamma
      real mean_d = exp(loc);
      real rate   = phi / mean_d;
      return gamma_lccdf(x | phi, rate);
    } else if (dist_type == 3) {  // weibull
      real scale = exp(loc);
      return weibull_lccdf(x | phi, scale);
    } else if (dist_type == 4) {  // burr XII: log CCDF = -k * log(1 + (x/lambda)^c)
      real u = log(x) - loc;
      return -kappa * log1p_exp(phi * u);
    } else if (dist_type == 5) {  // generalised gamma (Prentice)
      real gamma_shape = 1.0 / (kappa * kappa);
      real w = (log(x) - loc) / phi;
      return gamma_lccdf(gamma_shape * exp(kappa * w) | gamma_shape, 1);
    }
    return negative_infinity();
  }

  // Log likelihood for order statistic:
  // k-th order statistic out of n observations.
  real order_stat_logpdf_fun(real x, int n, int k, int dist_type, real loc, real phi, real kappa) {
    real log_f;
    real log_F   = 0.0;
    real log_1mF = 0.0;

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
      real log_term = log1p_exp(phi * u);
      log_f   = log(phi) + log(kappa) + (phi - 1) * u - loc - (kappa + 1) * log_term;
      log_1mF = -kappa * log_term;
      log_F   = log1m_exp(log_1mF);

    } else if (dist_type == 5) {  // generalised gamma (Prentice)
      real gamma_shape = 1.0 / (kappa * kappa);
      real w    = (log(x) - loc) / phi;
      real arg  = gamma_shape * exp(kappa * w);
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

  // Joint log density of three order statistics x1 < x2 < x3 with ranks
  // i1 < i2 < i3 from a sample of size n.
  //
  // log f(x1,x2,x3) = log_nc
  //   + (i1-1)*log_F(x1)                           [omitted when i1==1]
  //   + log_f(x1)
  //   + (i2-i1-1)*log[F(x2)-F(x1)]                 [omitted when i2-i1==1]
  //   + log_f(x2)
  //   + (i3-i2-1)*log[F(x3)-F(x2)]                 [omitted when i3-i2==1]
  //   + log_f(x3)
  //   + (n-i3)*log[1-F(x3)]                         [omitted when i3==n]
  //
  // where log_nc = lgamma(n+1) - lgamma(i1) - lgamma(i2-i1) - lgamma(i3-i2)
  //                             - lgamma(n-i3+1)
  real joint_3_order_stat_logpdf(
      real x1, real x2, real x3,
      int i1, int i2, int i3,
      int n,
      int dist_type, real loc, real phi, real kappa
  ) {
    real log_nc = lgamma(n + 1)
                  - lgamma(i1)
                  - lgamma(i2 - i1)
                  - lgamma(i3 - i2)
                  - lgamma(n - i3 + 1);

    real log_f1 = dist_logpdf_fun(x1, dist_type, loc, phi, kappa);
    real log_f2 = dist_logpdf_fun(x2, dist_type, loc, phi, kappa);
    real log_f3 = dist_logpdf_fun(x3, dist_type, loc, phi, kappa);

    real log_F1 = dist_log_cdf_fun(x1, dist_type, loc, phi, kappa);
    real log_F2 = dist_log_cdf_fun(x2, dist_type, loc, phi, kappa);
    real log_F3 = dist_log_cdf_fun(x3, dist_type, loc, phi, kappa);
    real log_1mF3 = dist_log_ccdf_fun(x3, dist_type, loc, phi, kappa);

    real lp = log_nc + log_f1 + log_f2 + log_f3;

    if (i1 > 1)       lp += (i1 - 1)       * log_F1;
    if (i2 - i1 > 1)  lp += (i2 - i1 - 1)  * log_diff_exp(log_F2, log_F1);
    if (i3 - i2 > 1)  lp += (i3 - i2 - 1)  * log_diff_exp(log_F3, log_F2);
    if (i3 < n)       lp += (n - i3)        * log_1mF3;

    return lp;
  }

  // Helper function to compute gamma quantile approximation
  real gamma_quantile_approx(real p, real shape, real scale) {
    real z = inv_Phi(p);
    if (shape > 1) {
      return shape * scale * pow(1 - 1.0/(9*shape) + z/(3*sqrt(shape)), 3);
    } else {
      return shape * scale * (1 + z / sqrt(shape));
    }
  }
}

data {
  int<lower=1> n_datasets;
  array[n_datasets] int<lower=1> n_obs;
  array[n_datasets] int<lower=1,upper=5> summary_type;
  int<lower=1,upper=5> dist_type;

  array[n_datasets] real<lower=0> obs_stat1;
  array[n_datasets] real<lower=0> obs_stat2;
  array[n_datasets] real<lower=0> obs_stat3;

  int<lower=0> n_freq_total;
  array[n_freq_total] real<lower=0> freq_value;
  array[n_freq_total] int<lower=1>  freq_count;
  array[n_datasets]   int<lower=0>  freq_start;
  array[n_datasets]   int<lower=0>  freq_len;

  array[n_freq_total] real<lower=0> freq_lower;
  array[n_freq_total] real<lower=0> freq_upper;

  real mu0_mean;
  real<lower=0> mu0_sd;

  real log_tau_mean;
  real<lower=0> log_tau_sd;

  real log_phi_mean;
  real<lower=0> log_phi_sd;

  real log_kappa_mean;
  real<lower=0> log_kappa_sd;
}

transformed data {
  real z_q25 = -0.6745;
  real z_q75 = 0.6745;
  real z_q90 = 1.28155;
  real z_q95 = 1.64485;

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
  real mu0;
  real log_tau;
  real log_phi;
  real log_kappa;
  vector[n_datasets] loc_d_raw;
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

    if (summary_type[d] == 1) {  // median + range: joint density of (min, median, max)
      int k_median = (n + 1) %/% 2;
      target += joint_3_order_stat_logpdf(
          obs_stat2[d], obs_stat1[d], obs_stat3[d],  // x1=min, x2=median, x3=max
          1, k_median, n,                             // i1=1, i2=k_median, i3=n
          n, dist_type, loc, phi, kappa);
    }

    else if (summary_type[d] == 2) {  // median + IQR: joint density of (q25, median, q75)
      int k_median = (n + 1) %/% 2;
      int k_q25 = (n + 1) %/% 4;
      if (k_q25 < 1) k_q25 = 1;
      int k_q75 = (3 * (n + 1)) %/% 4;
      if (k_q75 <= k_q25) k_q75 = k_q25 + 1;
      if (k_q75 > n) k_q75 = n;
      target += joint_3_order_stat_logpdf(
          obs_stat2[d], obs_stat1[d], obs_stat3[d],  // x1=q25, x2=median, x3=q75
          k_q25, k_median, k_q75,
          n, dist_type, loc, phi, kappa);
    }

    else if (summary_type[d] == 3) {  // mean + sd (unchanged from factorised model)
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

      } else if (dist_type == 4) {  // burr XII
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

      } else if (dist_type == 5) {  // generalised gamma
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

      if (dist_type <= 3) {
        obs_stat1[d] ~ normal(expected_mean, se_mean);
        obs_stat2[d] ~ normal(expected_sd, se_sd);
      }
    }

    else if (summary_type[d] == 4) {  // raw frequency table (unchanged)
      int s = freq_start[d];
      int len = freq_len[d];
      for (i in s:(s + len - 1)) {
        target += freq_count[i] * dist_logpdf_fun(freq_value[i], dist_type, loc, phi, kappa);
      }
    }

    else if (summary_type[d] == 5) {  // interval-censored frequency table (unchanged)
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
  // Joint log-likelihood per dataset (3 slots each for interface compatibility;
  // for types 1 and 2 the joint value is stored in slot idx and slots idx+1,
  // idx+2 are set to 0, since the joint density cannot be decomposed per statistic).
  vector[n_datasets * 3] log_lik;

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

      } else if (dist_type == 4) {  // burr XII
        real lam    = exp(loc_pred);
        pred_median = lam * pow(pow(0.5,  -1.0/kappa) - 1.0, 1.0/phi);
        pred_q25    = lam * pow(pow(0.75, -1.0/kappa) - 1.0, 1.0/phi);
        pred_q75    = lam * pow(pow(0.25, -1.0/kappa) - 1.0, 1.0/phi);
        pred_q90    = lam * pow(pow(0.10, -1.0/kappa) - 1.0, 1.0/phi);
        pred_q95    = lam * pow(pow(0.05, -1.0/kappa) - 1.0, 1.0/phi);
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
      int L = 2000;

      if (dist_type == 5) {
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

          } else if (dist_type == 4) {  // burr XII
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

      if (summary_type[d] == 1) {  // joint density of (min, median, max)
        int k_median = (n + 1) %/% 2;
        log_lik[idx] = joint_3_order_stat_logpdf(
            obs_stat2[d], obs_stat1[d], obs_stat3[d],
            1, k_median, n,
            n, dist_type, loc, phi, kappa);
        log_lik[idx + 1] = 0;
        log_lik[idx + 2] = 0;

      } else if (summary_type[d] == 2) {  // joint density of (q25, median, q75)
        int k_median = (n + 1) %/% 2;
        int k_q25    = (n + 1) %/% 4;
        if (k_q25 < 1) k_q25 = 1;
        int k_q75    = (3 * (n + 1)) %/% 4;
        if (k_q75 <= k_q25) k_q75 = k_q25 + 1;
        if (k_q75 > n) k_q75 = n;
        log_lik[idx] = joint_3_order_stat_logpdf(
            obs_stat2[d], obs_stat1[d], obs_stat3[d],
            k_q25, k_median, k_q75,
            n, dist_type, loc, phi, kappa);
        log_lik[idx + 1] = 0;
        log_lik[idx + 2] = 0;

      } else if (summary_type[d] == 3) {  // mean + sd (unchanged)
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

        if (dist_type != 4 || kappa * phi > 2.0) {
          log_lik[idx]     = normal_lpdf(obs_stat1[d] | expected_mean, expected_sd / sqrt(n));
          log_lik[idx + 1] = normal_lpdf(obs_stat2[d] | expected_sd,   expected_sd / sqrt(2 * (n - 1)));
          log_lik[idx + 2] = 0;
        }

      } else if (summary_type[d] == 4) {  // raw frequency table
        real ll_type4 = 0;
        int s = freq_start[d];
        int len = freq_len[d];
        for (i in s:(s + len - 1)) {
          ll_type4 += freq_count[i] * dist_logpdf_fun(freq_value[i], dist_type, loc, phi, kappa);
        }
        log_lik[idx]     = ll_type4;
        log_lik[idx + 1] = 0;
        log_lik[idx + 2] = 0;

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
        log_lik[idx + 1] = 0;
        log_lik[idx + 2] = 0;
      }

      idx += 3;
    }
  }
}
