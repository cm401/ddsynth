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
    // F_D(0) = 0 for every family here (support starts at 0), handled as an
    // exact special case rather than the general per-family formula below.
    // For Burr XII/generalised gamma (and empirically, Weibull's own
    // weibull_lcdf), that formula involves u = exp(shape * log(x/scale));
    // differentiating w.r.t. shape gives a factor of log(x/scale), which at
    // x=0 is -Inf * 0 = NaN even though the true limit is 0 (a removable
    // singularity). Verified this bypass leaves log-likelihood values
    // unchanged on non-zero x, while giving finite gradients at x=0 for all
    // 5 families.
    if (x <= 0) return negative_infinity();
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
    // 1 - F_D(0) = 1 for every family here; see dist_log_cdf_fun for why this
    // must be special-cased rather than left to the general formula.
    if (x <= 0) return 0;
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
    if (x1 >= x2 || x2 >= x3) return negative_infinity();

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

  // Rounded log-likelihood for the joint density of 3 order statistics:
  // integrates joint_3_order_stat_logpdf over a box of width `resolution`
  // around each reported (rounded) value, via 5-point-per-dimension tensor
  // Gauss-Legendre quadrature (125 nodes). Validated in R against a nested
  // adaptive-quadrature reference (Burr XII, n=15, resolution=1): relative
  // error ~5e-8, versus ~1e-4 at 3 points/dimension: 5 points gives a
  // comfortable safety margin at a still-modest node count. Lower bounds are
  // clamped just above 0 since delays are non-negative (only relevant for
  // the smallest of the three statistics, e.g. the reported min).
  //
  // Loses accuracy (~9% error in a case checked in the code review) when two
  // of the three reported statistics are within `resolution` of each other,
  // since joint_3_order_stat_logpdf's ordering constraint (x1<x2<x3) puts a
  // hard zero inside the integration box. Known limitation of this secondary
  // sensitivity model; see REVISION_TODO.md.
  real joint_3_order_stat_rounded_loglik(
      real v1, real v2, real v3,
      int i1, int i2, int i3,
      int n, int dist_type, real loc, real phi, real kappa, real resolution
  ) {
    array[5] real t = {-0.9061798459386640, -0.5384693101056831, 0.0,
                         0.5384693101056831,  0.9061798459386640};
    array[5] real w = { 0.2369268850561891,  0.4786286704993665,
                         0.5688888888888889,
                         0.4786286704993665,  0.2369268850561891};

    real half_res = resolution / 2;
    real v1_lo = fmax(v1 - half_res, 1e-6); real v1_hi = v1 + half_res;
    real v2_lo = fmax(v2 - half_res, 1e-6); real v2_hi = v2 + half_res;
    real v3_lo = fmax(v3 - half_res, 1e-6); real v3_hi = v3 + half_res;
    real mid1 = 0.5 * (v1_hi + v1_lo); real half1 = 0.5 * (v1_hi - v1_lo);
    real mid2 = 0.5 * (v2_hi + v2_lo); real half2 = 0.5 * (v2_hi - v2_lo);
    real mid3 = 0.5 * (v3_hi + v3_lo); real half3 = 0.5 * (v3_hi - v3_lo);

    array[125] real log_terms;
    int idx = 1;
    for (a in 1:5) {
      real x1 = mid1 + half1 * t[a];
      for (b in 1:5) {
        real x2 = mid2 + half2 * t[b];
        for (c in 1:5) {
          real x3 = mid3 + half3 * t[c];
          log_terms[idx] = log(w[a]) + log(w[b]) + log(w[c])
                           + joint_3_order_stat_logpdf(x1, x2, x3, i1, i2, i3, n,
                                                        dist_type, loc, phi, kappa);
          idx += 1;
        }
      }
    }
    return log_sum_exp(log_terms) + log(half1) + log(half2) + log(half3);
  }

  // Log-likelihood contribution for a raw frequency-table count that has
  // been rounded to the nearest multiple of `resolution` days: log P(v -
  // resolution/2 <= X < v + resolution/2), rather than treating v as an
  // exact continuous observation (summary_type 4's freq_value, and
  // summary_type 5's degenerate freq_lower==freq_upper point entries).
  // This is the same rounding treatment already applied to reported order
  // statistics, for consistency across every summary type where a single
  // reported value stands in for an underlying continuous, day-rounded
  // observation.
  //
  // No explicit clamp is needed on the lower bound, even when it goes
  // negative (v < resolution/2): unlike the order statistic case, this uses
  // a CDF *difference*, never the point density, so dist_log_cdf_fun's own
  // x<=0 special case (identical for any x<=0, negative or not) already
  // gives the correct result directly.
  real rounded_freq_loglik_fun(real v, int dist_type, real loc, real phi, real kappa, real resolution) {
    real half_res = resolution / 2;
    return log_diff_exp(
      dist_log_cdf_fun(v + half_res, dist_type, loc, phi, kappa),
      dist_log_cdf_fun(v - half_res, dist_type, loc, phi, kappa)
    );
  }

  // Log-density of a bivariate normal with mean (mx, my), variances (vx, vy)
  // and covariance cxy. Used for the joint mean/SD likelihood (summary_type
  // 3), which needs a nonzero covariance term for every right-skewed family
  // used here.
  //
  // Guards against a non-positive-definite or non-finite covariance (which
  // can arise from floating-point cancellation, or from the 3rd/4th central
  // moments overflowing at extreme parameter draws explored during warmup)
  // by treating it as a rejected region of parameter space (-Inf) rather
  // than propagating a NaN into target.
  real bvn_log_dens(real x, real y, real mx, real my, real vx, real vy, real cxy) {
    real det  = vx * vy - cxy * cxy;
    if (!(vx > 0) || !(vy > 0) || !(det > 0) || !(det < positive_infinity())) {
      return negative_infinity();
    }
    real dx   = x - mx;
    real dy   = y - my;
    real quad = (dx * dx * vy - 2 * dx * dy * cxy + dy * dy * vx) / det;
    return -log(2 * pi()) - 0.5 * log(det) - 0.5 * quad;
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

  // Log of integral_{ex_l}^{ex_r} F_D(T - e) de: truncation denominator
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

  // Log of integral_{ex_l}^{ex_r} [1 - F_D(T - e)] de: right-censored contribution
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
  // Identical to the version in hierarchical_data_synthesis_summary_stats.stan.
  real dc_log_lik(real ex_l, real ex_r, real ev_l, real ev_r,
                  int dist_type, real loc, real phi, real kappa, real resolution) {
    // Both endpoints point-observed: the implied delay ev_l - ex_l is itself
    // a rounded reported value, not an exact continuous observation, and is
    // treated the same way via rounded_freq_loglik_fun as summary_type 4/5's
    // degenerate case, rather than the raw point density dist_logpdf_fun
    // (singular/undefined at a delay of exactly 0 for several families).
    if (ex_l == ex_r && ev_l == ev_r) {
      return rounded_freq_loglik_fun(ev_l - ex_l, dist_type, loc, phi, kappa, resolution);
    }
    if (ex_l == ex_r) {
      return log_diff_exp(
        dist_log_cdf_fun(ev_r - ex_l, dist_type, loc, phi, kappa),
        dist_log_cdf_fun(ev_l - ex_l, dist_type, loc, phi, kappa)
      );
    }
    if (ev_l == ev_r) {
      return log_diff_exp(
        dist_log_cdf_fun(ev_l - ex_l, dist_type, loc, phi, kappa),
        dist_log_cdf_fun(ev_l - ex_r, dist_type, loc, phi, kappa)
      );
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
  int<lower=1> n_datasets;
  array[n_datasets] int<lower=1> n_obs;
  array[n_datasets] int<lower=1,upper=7> summary_type; // 1=median+range, 2=median+IQR, 3=mean+sd, 4=raw freq table, 5=interval-censored freq table, 6=double interval-censored freq table, 7=double interval-censored with right truncation/censoring
  int<lower=1,upper=5> dist_type;

  array[n_datasets] real<lower=0> obs_stat1;
  array[n_datasets] real<lower=0> obs_stat2;
  array[n_datasets] real<lower=0> obs_stat3;

  // Day-fraction rounding resolution (e.g. 1 = whole day, 1/24 = hourly); see
  // detect_resolution() in R/utils.R. Used by summary_type 1/2 (reported
  // order statistics), summary_type 4/5 (raw/degenerate frequency-table
  // counts), and summary_type 6/7's own degenerate (both-endpoints-point-
  // observed) case in dc_log_lik, all treated as rounded point values.
  // Ignored for summary_type 3 (populate with 1, unused).
  array[n_datasets] real<lower=0> resolution;

  int<lower=0> n_freq_total;
  array[n_freq_total] real<lower=0> freq_value;
  array[n_freq_total] int<lower=1>  freq_count;
  array[n_datasets]   int<lower=0>  freq_start;
  array[n_datasets]   int<lower=0>  freq_len;

  array[n_freq_total] real<lower=0> freq_lower;
  array[n_freq_total] real<lower=0> freq_upper;

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

  // Non-strict (<=), matching the factorised model: day-rounding a
  // short-incubation dataset routinely produces exact ties between
  // reported statistics (e.g. Q1 == median), which a strict check would
  // reject and abort the entire fit on. joint_3_order_stat_rounded_loglik
  // still loses accuracy near (but not at) a tie; see its docstring and
  // REVISION_TODO.md.
  for (d in 1:n_datasets) {
    if (summary_type[d] == 1) {
      if (obs_stat2[d] > obs_stat1[d] || obs_stat1[d] > obs_stat3[d]) {
        reject("For summary_type=1, must have min <= median <= max");
      }
      if (resolution[d] <= 0) {
        reject("For summary_type=1, resolution must be > 0");
      }
    } else if (summary_type[d] == 2) {
      if (obs_stat2[d] > obs_stat1[d] || obs_stat1[d] > obs_stat3[d]) {
        reject("For summary_type=2, must have q25 <= median <= q75");
      }
      if (resolution[d] <= 0) {
        reject("For summary_type=2, resolution must be > 0");
      }
    } else if (summary_type[d] == 3) {
      // A sample SD is undefined for a single observation; the mean/SD
      // likelihood also divides by (n-1), which is degenerate at n=1.
      if (n_obs[d] < 2) {
        reject("For summary_type=3, n_obs must be >= 2 (a sample SD requires at least 2 observations)");
      }
    } else if (summary_type[d] == 4) {
      if (freq_len[d] == 0) {
        reject("For summary_type=4, freq_len must be > 0");
      }
      if (resolution[d] <= 0) {
        reject("For summary_type=4, resolution must be > 0");
      }
    } else if (summary_type[d] == 5) {
      if (freq_len[d] == 0) {
        reject("For summary_type=5, freq_len must be > 0");
      }
      if (resolution[d] <= 0) {
        reject("For summary_type=5, resolution must be > 0");
      }
    } else if (summary_type[d] == 6) {
      if (freq_len[d] == 0) {
        reject("For summary_type=6, freq_len must be > 0");
      }
      if (resolution[d] <= 0) {
        reject("For summary_type=6, resolution must be > 0");
      }
    } else if (summary_type[d] == 7) {
      if (freq_len[d] == 0) {
        reject("For summary_type=7, freq_len must be > 0");
      }
      if (truncation_time[d] <= 0) {
        reject("For summary_type=7, truncation_time must be > 0");
      }
      if (resolution[d] <= 0) {
        reject("For summary_type=7, resolution must be > 0");
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

    if (summary_type[d] == 1) {  // median + range: joint density of (min, median, max), rounded
      int k_median = (n + 1) %/% 2;
      target += joint_3_order_stat_rounded_loglik(
          obs_stat2[d], obs_stat1[d], obs_stat3[d],  // x1=min, x2=median, x3=max
          1, k_median, n,                             // i1=1, i2=k_median, i3=n
          n, dist_type, loc, phi, kappa, resolution[d]);
    }

    else if (summary_type[d] == 2) {  // median + IQR: joint density of (q25, median, q75), rounded
      int k_median = (n + 1) %/% 2;
      int k_q25 = (n + 1) %/% 4;
      if (k_q25 < 1) k_q25 = 1;
      int k_q75 = (3 * (n + 1)) %/% 4;
      if (k_q75 <= k_q25) k_q75 = k_q25 + 1;
      if (k_q75 > n) k_q75 = n;
      target += joint_3_order_stat_rounded_loglik(
          obs_stat2[d], obs_stat1[d], obs_stat3[d],  // x1=q25, x2=median, x3=q75
          k_q25, k_median, k_q75,
          n, dist_type, loc, phi, kappa, resolution[d]);
    }

    else if (summary_type[d] == 3) {  // mean + sd (identical to factorised model, see there for derivation)
      real mean_d;
      real var_d;  // sigma_d^2
      real mu3;
      real mu4;
      int  moments_ok = 1;

      if (dist_type == 1) {  // lognormal
        real omega = exp(phi^2);
        mean_d = exp(loc + phi^2 / 2);
        var_d  = mean_d^2 * (omega - 1);
        mu3    = (omega + 2) * (omega - 1)^2 * mean_d^3;
        mu4    = (omega^4 + 2*omega^3 + 3*omega^2 - 3) * (omega - 1)^2 * mean_d^4;

      } else if (dist_type == 2) {  // gamma
        real shape = phi;
        mean_d = exp(loc);
        var_d  = mean_d^2 / shape;
        mu3    = 2 * mean_d^3 / shape^2;
        mu4    = (3 + 6 / shape) * var_d^2;

      } else if (dist_type == 3) {  // weibull
        real shape = phi;
        real scale = exp(loc);
        real g1 = tgamma(1 + 1.0/shape);
        real g2 = tgamma(1 + 2.0/shape);
        real g3 = tgamma(1 + 3.0/shape);
        real g4 = tgamma(1 + 4.0/shape);
        mean_d = scale * g1;
        var_d  = scale^2 * (g2 - g1^2);
        mu3    = scale^3 * (g3 - 3*g1*g2 + 2*g1^3);
        mu4    = scale^4 * (g4 - 4*g1*g3 + 6*g1^2*g2 - 3*g1^4);

      } else if (dist_type == 4) {  // burr XII: lambda=exp(loc), c=phi, k=kappa
        // 4th moment (needed for Var(S^2)) requires k*c > 4, stricter than
        // the k*c > 2 needed for the variance alone.
        if (kappa * phi <= 4.0) {
          moments_ok = 0;
          target += negative_infinity();
        } else {
          real lam = exp(loc);
          real m1 = lam   * kappa * exp(lbeta(kappa - 1.0/phi, 1 + 1.0/phi));
          real m2 = lam^2 * kappa * exp(lbeta(kappa - 2.0/phi, 1 + 2.0/phi));
          real m3 = lam^3 * kappa * exp(lbeta(kappa - 3.0/phi, 1 + 3.0/phi));
          real m4 = lam^4 * kappa * exp(lbeta(kappa - 4.0/phi, 1 + 4.0/phi));
          mean_d = m1;
          var_d  = fabs(m2 - m1^2);
          mu3    = m3 - 3*m1*m2 + 2*m1^3;
          mu4    = m4 - 4*m1*m3 + 6*m1^2*m2 - 3*m1^4;
        }

      } else if (dist_type == 5) {  // generalised gamma (Prentice): mu=loc, sigma=phi, Q=kappa
        real gamma_shape = 1.0 / (kappa * kappa);
        real m1 = exp(loc     + 2.0*phi/kappa * log(kappa) + lgamma(gamma_shape + phi/kappa)   - lgamma(gamma_shape));
        real m2 = exp(2.0*loc + 4.0*phi/kappa * log(kappa) + lgamma(gamma_shape + 2.0*phi/kappa) - lgamma(gamma_shape));
        real m3 = exp(3.0*loc + 6.0*phi/kappa * log(kappa) + lgamma(gamma_shape + 3.0*phi/kappa) - lgamma(gamma_shape));
        real m4 = exp(4.0*loc + 8.0*phi/kappa * log(kappa) + lgamma(gamma_shape + 4.0*phi/kappa) - lgamma(gamma_shape));
        mean_d = m1;
        var_d  = fabs(m2 - m1^2);
        mu3    = m3 - 3*m1*m2 + 2*m1^3;
        mu4    = m4 - 4*m1*m3 + 6*m1^2*m2 - 3*m1^4;
      }

      if (moments_ok == 1) {
        real var_mean = var_d / n;
        real var_var2 = (mu4 - ((n - 3.0) / (n - 1.0)) * var_d^2) / n;
        real cov_mv2  = mu3 / n;
        target += bvn_log_dens(obs_stat1[d], obs_stat2[d]^2, mean_d, var_d, var_mean, var_var2, cov_mv2);
      }
    }

    else if (summary_type[d] == 4) {  // raw frequency table
      // Each distinct observed value is treated as rounded to the nearest
      // multiple of `resolution` days, not as an exact continuous observation.
      int s = freq_start[d];
      int len = freq_len[d];
      for (i in s:(s + len - 1)) {
        target += freq_count[i] * rounded_freq_loglik_fun(freq_value[i], dist_type, loc, phi, kappa, resolution[d]);
      }
    }

    else if (summary_type[d] == 5) {  // interval-censored frequency table
      // When lower == upper (no reported range), the value is treated as a
      // rounded point observation, like summary_type 4.
      int s = freq_start[d];
      int len = freq_len[d];
      for (i in s:(s + len - 1)) {
        if (freq_lower[i] == freq_upper[i]) {
          target += freq_count[i] * rounded_freq_loglik_fun(freq_lower[i], dist_type, loc, phi, kappa, resolution[d]);
        } else {
          real log_cdf_u = dist_log_cdf_fun(freq_upper[i], dist_type, loc, phi, kappa);
          real log_cdf_l = dist_log_cdf_fun(freq_lower[i], dist_type, loc, phi, kappa);
          target += freq_count[i] * log_diff_exp(log_cdf_u, log_cdf_l);
        }
      }
    }

    else if (summary_type[d] == 6) {  // double interval-censored frequency table
      int s = freq_start[d];
      int len = freq_len[d];
      for (i in s:(s + len - 1)) {
        target += freq_count[i] * dc_log_lik(
          freq_lower[i], freq_upper[i],
          event_lower[i], event_upper[i],
          dist_type, loc, phi, kappa, resolution[d]
        );
      }
    }

    else if (summary_type[d] == 7) {  // double interval-censored + right truncation/censoring
      int s = freq_start[d];
      int len = freq_len[d];
      real T = truncation_time[d];
      for (i in s:(s + len - 1)) {
        if (event_observed[i] == 1) {
          target += freq_count[i] * (
            dc_log_lik(freq_lower[i], freq_upper[i],
                       event_lower[i], event_upper[i],
                       dist_type, loc, phi, kappa, resolution[d])
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

      if (summary_type[d] == 1) {  // joint density of (min, median, max), rounded
        int k_median = (n + 1) %/% 2;
        log_lik[idx] = joint_3_order_stat_rounded_loglik(
            obs_stat2[d], obs_stat1[d], obs_stat3[d],
            1, k_median, n,
            n, dist_type, loc, phi, kappa, resolution[d]);
        log_lik[idx + 1] = 0;
        log_lik[idx + 2] = 0;

      } else if (summary_type[d] == 2) {  // joint density of (q25, median, q75), rounded
        int k_median = (n + 1) %/% 2;
        int k_q25    = (n + 1) %/% 4;
        if (k_q25 < 1) k_q25 = 1;
        int k_q75    = (3 * (n + 1)) %/% 4;
        if (k_q75 <= k_q25) k_q75 = k_q25 + 1;
        if (k_q75 > n) k_q75 = n;
        log_lik[idx] = joint_3_order_stat_rounded_loglik(
            obs_stat2[d], obs_stat1[d], obs_stat3[d],
            k_q25, k_median, k_q75,
            n, dist_type, loc, phi, kappa, resolution[d]);
        log_lik[idx + 1] = 0;
        log_lik[idx + 2] = 0;

      } else if (summary_type[d] == 3) {  // mean + sd (see model block for derivation)
        real mean_d;
        real var_d;
        real mu3;
        real mu4;
        int  moments_ok = 1;

        if (dist_type == 1) {
          real omega = exp(phi^2);
          mean_d = exp(loc + phi^2 / 2);
          var_d  = mean_d^2 * (omega - 1);
          mu3    = (omega + 2) * (omega - 1)^2 * mean_d^3;
          mu4    = (omega^4 + 2*omega^3 + 3*omega^2 - 3) * (omega - 1)^2 * mean_d^4;

        } else if (dist_type == 2) {
          real shape = phi;
          mean_d = exp(loc);
          var_d  = mean_d^2 / shape;
          mu3    = 2 * mean_d^3 / shape^2;
          mu4    = (3 + 6 / shape) * var_d^2;

        } else if (dist_type == 3) {
          real shape = phi;
          real scale = exp(loc);
          real g1 = tgamma(1 + 1.0/shape);
          real g2 = tgamma(1 + 2.0/shape);
          real g3 = tgamma(1 + 3.0/shape);
          real g4 = tgamma(1 + 4.0/shape);
          mean_d = scale * g1;
          var_d  = scale^2 * (g2 - g1^2);
          mu3    = scale^3 * (g3 - 3*g1*g2 + 2*g1^3);
          mu4    = scale^4 * (g4 - 4*g1*g3 + 6*g1^2*g2 - 3*g1^4);

        } else if (dist_type == 4) {
          if (kappa * phi <= 4.0) {
            moments_ok = 0;
          } else {
            real lam = exp(loc);
            real m1 = lam   * kappa * exp(lbeta(kappa - 1.0/phi, 1 + 1.0/phi));
            real m2 = lam^2 * kappa * exp(lbeta(kappa - 2.0/phi, 1 + 2.0/phi));
            real m3 = lam^3 * kappa * exp(lbeta(kappa - 3.0/phi, 1 + 3.0/phi));
            real m4 = lam^4 * kappa * exp(lbeta(kappa - 4.0/phi, 1 + 4.0/phi));
            mean_d = m1;
            var_d  = fabs(m2 - m1^2);
            mu3    = m3 - 3*m1*m2 + 2*m1^3;
            mu4    = m4 - 4*m1*m3 + 6*m1^2*m2 - 3*m1^4;
          }

        } else if (dist_type == 5) {
          real gamma_shape = 1.0 / (kappa * kappa);
          real m1 = exp(loc     + 2.0*phi/kappa * log(kappa) + lgamma(gamma_shape + phi/kappa)   - lgamma(gamma_shape));
          real m2 = exp(2.0*loc + 4.0*phi/kappa * log(kappa) + lgamma(gamma_shape + 2.0*phi/kappa) - lgamma(gamma_shape));
          real m3 = exp(3.0*loc + 6.0*phi/kappa * log(kappa) + lgamma(gamma_shape + 3.0*phi/kappa) - lgamma(gamma_shape));
          real m4 = exp(4.0*loc + 8.0*phi/kappa * log(kappa) + lgamma(gamma_shape + 4.0*phi/kappa) - lgamma(gamma_shape));
          mean_d = m1;
          var_d  = fabs(m2 - m1^2);
          mu3    = m3 - 3*m1*m2 + 2*m1^3;
          mu4    = m4 - 4*m1*m3 + 6*m1^2*m2 - 3*m1^4;
        }

        if (moments_ok == 1) {
          real var_mean = var_d / n;
          real var_var2 = (mu4 - ((n - 3.0) / (n - 1.0)) * var_d^2) / n;
          real cov_mv2  = mu3 / n;
          log_lik[idx]  = bvn_log_dens(obs_stat1[d], obs_stat2[d]^2, mean_d, var_d, var_mean, var_var2, cov_mv2);
        } else {
          log_lik[idx]  = negative_infinity();
        }
        log_lik[idx + 1] = 0;  // unused
        log_lik[idx + 2] = 0;  // unused

      } else if (summary_type[d] == 4) {  // raw frequency table
        real ll_type4 = 0;
        int s = freq_start[d];
        int len = freq_len[d];
        for (i in s:(s + len - 1)) {
          ll_type4 += freq_count[i] * rounded_freq_loglik_fun(freq_value[i], dist_type, loc, phi, kappa, resolution[d]);
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
            ll_type5 += freq_count[i] * rounded_freq_loglik_fun(freq_lower[i], dist_type, loc, phi, kappa, resolution[d]);
          } else {
            real log_cdf_u = dist_log_cdf_fun(freq_upper[i], dist_type, loc, phi, kappa);
            real log_cdf_l = dist_log_cdf_fun(freq_lower[i], dist_type, loc, phi, kappa);
            ll_type5 += freq_count[i] * log_diff_exp(log_cdf_u, log_cdf_l);
          }
        }
        log_lik[idx]     = ll_type5;
        log_lik[idx + 1] = 0;
        log_lik[idx + 2] = 0;

      } else if (summary_type[d] == 6) {  // double interval-censored frequency table
        real ll_type6 = 0;
        int s = freq_start[d];
        int len = freq_len[d];
        for (i in s:(s + len - 1)) {
          ll_type6 += freq_count[i] * dc_log_lik(
            freq_lower[i], freq_upper[i],
            event_lower[i], event_upper[i],
            dist_type, loc, phi, kappa, resolution[d]
          );
        }
        log_lik[idx]     = ll_type6;
        log_lik[idx + 1] = 0;
        log_lik[idx + 2] = 0;

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
                         dist_type, loc, phi, kappa, resolution[d])
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
