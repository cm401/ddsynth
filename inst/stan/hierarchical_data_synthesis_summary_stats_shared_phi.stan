// hierarchical_data_synthesis_summary_stats_shared_phi.stan
//
// Pre-point-4 variant, deliberately kept alongside the hierarchical model:
// a single shared scalar `phi` per pathogen rather than per-study `phi_d`.
// Identical to hierarchical_data_synthesis_summary_stats.stan as of commit
// ead8803 (before the point-4 phi_d/omega reparameterization).
//
// Used for gamma (POINT4_LIKELIHOOD_MATHS.md Part E.5): gamma's shape
// parameter is a concentration-type parameter whose per-study Fisher
// information about phi vanishes as phi grows (exactly -1/(2*phi), the
// fastest decay of any family in this corpus, compounded by gamma
// currently sitting at the largest fitted phi of any family here), so the
// non-centered phi_d/omega hierarchy creates a genuine funnel for gamma
// rather than just needing a better prior. See Part E.4.1 for the general
// "spread-type vs. concentration-type shape parameter" intuition.
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

  // Log-density of the k-th order statistic (k out of n observations).
  // Fully log-scale: uses _lcdf/_lccdf directly and never exponentiates a
  // CDF/CCDF value, so it stays well-conditioned arbitrarily deep in either
  // tail. Terms with a zero coefficient are omitted explicitly to avoid
  // 0 * (-Inf) = NaN in Stan's autodiff (k==1: coefficient of log_F is 0;
  // k==n: coefficient of log_1mF is 0).
  real order_stat_logpdf_fun(real x, int n, int k, int dist_type, real loc, real phi, real kappa) {
    real log_f   = dist_logpdf_fun(x, dist_type, loc, phi, kappa);
    real log_F   = 0.0;
    real log_1mF = 0.0;
    if (k > 1) log_F   = dist_log_cdf_fun(x, dist_type, loc, phi, kappa);
    if (k < n) log_1mF = dist_log_ccdf_fun(x, dist_type, loc, phi, kappa);

    // Order-statistic normalising constant is n!/((k-1)!(n-k)!) = k*choose(n,k),
    // not choose(n,k) alone. See hierarchical_data_synthesis_summary_stats.stan.
    real log_dens = lchoose(n, k) + log(k) + log_f;
    if (k > 1) log_dens += (k - 1) * log_F;
    if (k < n) log_dens += (n - k) * log_1mF;
    return log_dens;
  }

  // Log-likelihood contribution for a reported order statistic that has been
  // rounded to the nearest multiple of `resolution` days (1 = whole day,
  // 1/24 = hourly, etc.): log P(v - resolution/2 <= X_(k) < v + resolution/2),
  // rather than treating v as an exact continuous observation (which discards
  // a material part of the uncertainty for short-incubation pathogens).
  // Lower bound is clamped just above 0 since delays are non-negative.
  //
  // Computed by composite Gauss-Legendre quadrature (24 panels x 7-point
  // rule = 168 nodes) directly on order_stat_logpdf_fun, entirely in log
  // space via log_sum_exp.
  //
  // This deliberately avoids computing the order-statistic CDF via Stan's
  // beta_lcdf(F_D(x) | k, n-k+1): that looks like the natural closed form
  // (F_{X(k)}(x) = I_{F_D(x)}(k, n-k+1), the regularized incomplete beta),
  // but Stan's beta_lcdf evaluates the incomplete beta on the linear
  // probability scale first and only takes its log afterwards, so it
  // silently returns -Inf whenever the true value is below ~1e-308, which
  // happens routinely for large n even when F_D(x) itself is an entirely
  // ordinary, well-represented number (e.g. n=1000, k=500, F_D=0.01 has an
  // exact log-CDF of -1618, but beta_lcdf(0.01 | 500, 501) returns -Inf).
  // A single-panel Gauss-Legendre rule over the full box has the opposite
  // failure mode at large n: the order statistic's own sampling distribution
  // can be far narrower than the box (SE ~ 1/sqrt(n)), so a low-order
  // polynomial rule cannot resolve it and does not reliably converge with
  // more points at a single panel. Composite (multi-panel) quadrature
  // resolves both problems: validated in R against a high-precision adaptive
  // reference for n = 15 to 50,000 at resolution = 1 (the widest, hence
  // hardest to resolve, box this function is used with), both near the mode
  // and deep in the tail (the exact regime that broke beta_lcdf); error is
  // below 1e-9 near the mode (where posterior mass concentrates) even at
  // n=50,000, and the residual error deep in an implausible-parameter tail
  // region is at most ~0.02 log-units, immaterial next to a log-density of
  // that magnitude (around -450 in the case tested). A narrower box (finer
  // resolution) is strictly easier to resolve, since the accuracy problem is
  // driven by the box being wide relative to the order statistic's SE.
  real order_stat_rounded_loglik_fun(real v, int n, int k, int dist_type, real loc, real phi, real kappa, real resolution) {
    int n_panels = 24;
    array[7] real t = {-0.9491079123427585, -0.7415311855993945,
                       -0.4058451513773832,  0.0,
                        0.4058451513773832,  0.7415311855993945,
                        0.9491079123427585};
    array[7] real w = { 0.1294849661688697,  0.2797053914892767,
                        0.3818300505051189,  0.4179591836734694,
                        0.3818300505051189,  0.2797053914892767,
                        0.1294849661688697};
    real half_res = resolution / 2;
    real v_lo = fmax(v - half_res, 0);
    real v_hi = v + half_res;
    real panel_width = (v_hi - v_lo) / n_panels;
    array[n_panels * 7] real log_terms;
    int idx = 1;
    for (p in 1:n_panels) {
      real a    = v_lo + (p - 1) * panel_width;
      real mid  = a + 0.5 * panel_width;
      real half = 0.5 * panel_width;
      for (j in 1:7) {
        real x = mid + half * t[j];
        log_terms[idx] = log(w[j]) + log(half)
                         + order_stat_logpdf_fun(x, n, k, dist_type, loc, phi, kappa);
        idx += 1;
      }
    }
    return log_sum_exp(log_terms);
  }

  // Log-likelihood contribution for a raw frequency-table count that has
  // been rounded to the nearest multiple of `resolution` days: log P(v -
  // resolution/2 <= X < v + resolution/2), rather than treating v as an
  // exact continuous observation (summary_type 4's freq_value, and
  // summary_type 5's degenerate freq_lower==freq_upper point entries).
  // This is the same rounding treatment already applied to reported order
  // statistics (see order_stat_rounded_loglik_fun), for consistency across
  // every summary type where a single reported value stands in for an
  // underlying continuous, day-rounded observation.
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
    // Wilson-Hilferty approximation for gamma quantiles
    real z = inv_Phi(p);
    if (shape > 1) {
      return shape * scale * pow(1 - 1.0/(9*shape) + z/(3*sqrt(shape)), 3);
    } else {
      return shape * scale * (1 + z / sqrt(shape));   // Fallback for small shape
    }
  }

  // Log of integral_{ex_l}^{ex_r} F_D(T - e) de: the truncation denominator
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
  //
  // Computes log P(E in [ex_l, ex_r], O in [ev_l, ev_r]) proportional to theta,
  // where D = O - E ~ f_D(d; theta), assuming E uniform on [ex_l, ex_r] and
  // O uniform on [ev_l, ev_r].
  //
  // The normalising constant (ex_r - ex_l)(ev_r - ev_l) is constant in theta
  // and is dropped, consistent with how type 5 drops 1/(ex_r - ex_l).
  //
  // Four cases:
  //   Both endpoints point-observed          -> rounded to `resolution`, like summary_type 4/5
  //   Point exposure, interval event         -> log[F(ev_r - ex_l) - F(ev_l - ex_l)]
  //   Interval exposure, point event (= type 5) -> log[F(ev_l - ex_l) - F(ev_l - ex_r)]
  //   Both interval-censored (general)       -> 7-point Gauss-Legendre quadrature
  real dc_log_lik(real ex_l, real ex_r, real ev_l, real ev_r,
                  int dist_type, real loc, real phi, real kappa, real resolution) {
    // 7-point Gauss-Legendre nodes/weights on [-1, 1] (Abramowitz & Stegun
    // table 25.4), shared by every case below. Every case quadratures over
    // the exposure side (its genuinely reported window in Cases 2-4, or a
    // resolution-wide window standing in for a point report in Cases 1/3),
    // treating the event side as a closed-form probability at each exposure
    // node: dist_log_cdf_fun's CDF-difference when event is a genuine
    // window (Cases 2, 4), or rounded_freq_loglik_fun's resolution-wide
    // probability when event is itself a point report (Cases 1, 3). This
    // keeps every case's Jacobian-dropping convention identical to Case 4's
    // pre-existing one: calling this function with a point widened to its
    // own resolution window gives the exact same result, node for node, as
    // calling Case 4 directly with that window as an explicit interval.
    array[7] real t = {-0.9491079123427585, -0.7415311855993945,
                       -0.4058451513773832,  0.0,
                        0.4058451513773832,  0.7415311855993945,
                        0.9491079123427585};
    array[7] real w = { 0.1294849661688697,  0.2797053914892767,
                        0.3818300505051189,  0.4179591836734694,
                        0.3818300505051189,  0.2797053914892767,
                        0.1294849661688697};
    real half_res = resolution / 2;

    // Case 1: both endpoints point-observed. Exposure and event are each
    // independently day-rounded, so the delay's implied distribution is the
    // convolution of two independent Uniform(-resolution/2, resolution/2)
    // rounding errors: a triangular kernel of half-width `resolution`, not
    // the single rounded_freq_loglik_fun evaluation used before. Integrating
    // out the exposure side analytically collapses the general double
    // integral to this 1D quadrature over rounded_freq_loglik_fun exactly,
    // no approximation in that reduction step.
    if (ex_l == ex_r && ev_l == ev_r) {
      real d = ev_l - ex_l;
      array[7] real log_terms;
      for (k in 1:7) {
        real v_k = half_res * t[k];
        log_terms[k] = log(w[k]) + rounded_freq_loglik_fun(d - v_k, dist_type, loc, phi, kappa, resolution);
      }
      return log_sum_exp(log_terms);
    }
    // Case 2: point exposure, interval event. Exposure is day-rounded; event
    // stays exactly as reported (a genuine censoring window, not itself
    // rounded further). Quadrature averages the point-exposure formula over
    // exposure's own rounding window instead of evaluating it once at the
    // exact reported value.
    if (ex_l == ex_r) {
      array[7] real log_terms;
      for (k in 1:7) {
        real e_k = ex_l + half_res * t[k];
        log_terms[k] = log(w[k]) + log_diff_exp(
          dist_log_cdf_fun(ev_r - e_k, dist_type, loc, phi, kappa),
          dist_log_cdf_fun(ev_l - e_k, dist_type, loc, phi, kappa)
        );
      }
      return log_sum_exp(log_terms);
    }
    // Case 3: interval exposure, point event. Quadrature over the genuinely
    // reported exposure window, exactly like Case 4, but with the event
    // side's closed-form CDF-difference replaced by rounded_freq_loglik_fun
    // (event is itself a point report, standing in for its own resolution
    // window) instead of the exact reported event window Case 4 uses.
    if (ev_l == ev_r) {
      real mid = 0.5 * (ex_r + ex_l);
      real half = 0.5 * (ex_r - ex_l);
      array[7] real log_terms;
      for (k in 1:7) {
        real e_k = mid + half * t[k];
        log_terms[k] = log(w[k]) + rounded_freq_loglik_fun(ev_l - e_k, dist_type, loc, phi, kappa, resolution);
      }
      return log_sum_exp(log_terms);
    }
    // Case 4: both interval-censored (quadrature on [ex_l, ex_r]).
    // The Jacobian factor (ex_r - ex_l)/2 is constant in theta and is dropped.
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

  // Day-fraction rounding resolution (e.g. 1 = whole day, 1/24 = hourly); see
  // detect_resolution() in R/utils.R. Used by summary_type 1/2 (reported
  // order statistics), summary_type 4/5 (raw/degenerate frequency-table
  // counts), and summary_type 6/7's dc_log_lik whenever exposure and/or
  // event is point-observed (Cases 1-3), all treated as rounded values.
  // Ignored for summary_type 3 (populate with 1, unused).
  array[n_datasets] real<lower=0> resolution;

  // --- Frequency table data for summary_type == 4 ---
  // All datasets' frequency tables are stored in flat arrays.
  // For dataset d, its entries occupy indices freq_start[d] .. freq_start[d] + freq_len[d] - 1.
  int<lower=0> n_freq_total;                        // Total number of (value, count) pairs across all type-4 datasets
  array[n_freq_total] real<lower=0> freq_value;     // Observed day values for type 4, rounded to `resolution`; 0 allowed
  array[n_freq_total] int<lower=1>  freq_count;     // Number of individuals with that day value (types 4 and 5)
  array[n_datasets]   int<lower=0>  freq_start;     // 1-based start index into freq arrays for dataset d
  array[n_datasets]   int<lower=0>  freq_len;       // Number of distinct values for dataset d (0 if not type 4 or 5)

  // Interval bounds for summary_type == 5 (interval-censored frequency table).
  // For type 4 datasets these arrays are ignored (populate with zeros).
  // When freq_lower[i] == freq_upper[i] (no reported range), the value is
  // treated as a rounded point observation, like summary_type 4.
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
  // T is a study-design cutoff, not a rounded observation of a random event
  // time, so it is never treated as day-rounded the way exposure/event
  // windows are (see resolution, above). Set it using the same
  // time-encoding convention as that dataset's own windows.
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

    if (summary_type[d] == 1) {  // median + range (min, max), rounded to `resolution[d]`
      int k_median = (n + 1) %/% 2;
      target += order_stat_rounded_loglik_fun(obs_stat1[d], n, k_median, dist_type, loc, phi, kappa, resolution[d]);
      target += order_stat_rounded_loglik_fun(obs_stat2[d], n, 1, dist_type, loc, phi, kappa, resolution[d]);
      target += order_stat_rounded_loglik_fun(obs_stat3[d], n, n, dist_type, loc, phi, kappa, resolution[d]);
    }

    else if (summary_type[d] == 2) {  // median + IQR (q25, q75), rounded to `resolution[d]`
      int k_median = (n + 1) %/% 2;
      target += order_stat_rounded_loglik_fun(obs_stat1[d], n, k_median, dist_type, loc, phi, kappa, resolution[d]);

      int k_q25 = (n + 1) %/% 4;
      if (k_q25 < 1) k_q25 = 1;
      target += order_stat_rounded_loglik_fun(obs_stat2[d], n, k_q25, dist_type, loc, phi, kappa, resolution[d]);

      int k_q75 = (3 * (n + 1)) %/% 4;
      if (k_q75 <= k_q25) k_q75 = k_q25 + 1;
      if (k_q75 > n) k_q75 = n;
      target += order_stat_rounded_loglik_fun(obs_stat3[d], n, k_q75, dist_type, loc, phi, kappa, resolution[d]);
    }

    else if (summary_type[d] == 3) {  // mean + sd
      // Joint likelihood for (sample mean, sample variance), using the
      // exact finite-sample moments of these two statistics rather than
      // the univariate normal-theory Var(SD) and independence assumption
      // (correct only for Gaussian data). For the unbiased sample variance
      // S^2 (divisor n-1):
      //   Var(mean)   = sigma^2 / n                          [exact]
      //   Cov(mean,S^2) = mu3 / n                             [exact]
      //   Var(S^2)    = (mu4 - ((n-3)/(n-1)) * sigma^4) / n   [exact]
      // The likelihood is built for (mean, S^2) rather than (mean, SD)
      // because these two identities hold exactly for any n, whereas
      // transforming to SD via the delta method (S = sqrt(S^2)) adds an
      // extra approximation that a Monte Carlo check showed is inaccurate
      // at realistic sample sizes for skewed families. The S^2 -> SD
      // Jacobian depends only on the observed data, not on the model
      // parameters, so it is dropped without affecting inference.
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
        // E[X^r] = lambda^r * k * B(k - r/c, 1 + r/c), requires k*c > r.
        // The 4th moment (needed for Var(S^2)) requires k*c > 4, stricter
        // than the k*c > 2 needed for the variance alone. Hard barrier:
        // reject this region of parameter space for datasets that report
        // a mean and SD.
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
        // log E[X^r] = r*loc + (2*r*phi/kappa)*log(kappa)
        //              + lgamma(gamma_shape + r*phi/kappa) - lgamma(gamma_shape)
        // where gamma_shape = 1/kappa^2. Always finite (no kc-style threshold).
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
      // multiple of `resolution` days, exactly like a reported order
      // statistic (see rounded_freq_loglik_fun), not as an exact continuous
      // observation.
      int s = freq_start[d];
      int len = freq_len[d];
      for (i in s:(s + len - 1)) {
        target += freq_count[i] * rounded_freq_loglik_fun(freq_value[i], dist_type, loc, phi, kappa, resolution[d]);
      }
    }

    else if (summary_type[d] == 5) {  // interval-censored frequency table
      // Likelihood: count * log[ F(upper) - F(lower) ] for each interval.
      // Uses log_diff_exp(log_F_upper, log_F_lower) for numerical stability:
      // avoids catastrophic cancellation when the two CDF values are close.
      // When lower == upper (point observation, no reported range), it is
      // rounded like summary_type 4 rather than falling back to a point
      // log-pdf.
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
          dist_type, loc, phi, kappa, resolution[d]
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
        log_lik[idx]     = order_stat_rounded_loglik_fun(obs_stat1[d], n, k_median, dist_type, loc, phi, kappa, resolution[d]);
        log_lik[idx + 1] = order_stat_rounded_loglik_fun(obs_stat2[d], n, 1,        dist_type, loc, phi, kappa, resolution[d]);
        log_lik[idx + 2] = order_stat_rounded_loglik_fun(obs_stat3[d], n, n,        dist_type, loc, phi, kappa, resolution[d]);

      } else if (summary_type[d] == 2) {  // median + IQR
        int k_median = (n + 1) %/% 2;
        int k_q25    = (n + 1) %/% 4;
        if (k_q25 < 1) k_q25 = 1;
        int k_q75    = (3 * (n + 1)) %/% 4;
        if (k_q75 <= k_q25) k_q75 = k_q25 + 1;
        if (k_q75 > n) k_q75 = n;

        log_lik[idx]     = order_stat_rounded_loglik_fun(obs_stat1[d], n, k_median, dist_type, loc, phi, kappa, resolution[d]);
        log_lik[idx + 1] = order_stat_rounded_loglik_fun(obs_stat2[d], n, k_q25,   dist_type, loc, phi, kappa, resolution[d]);
        log_lik[idx + 2] = order_stat_rounded_loglik_fun(obs_stat3[d], n, k_q75,   dist_type, loc, phi, kappa, resolution[d]);

      } else if (summary_type[d] == 3) {  // mean + sd (see model block for derivation)
        real mean_d;
        real var_d;  // sigma_d^2
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
        // Sum log-likelihoods over all individuals, using the frequency table.
        // Stored as a single scalar in log_lik[idx]; slots idx+1 and idx+2 are 0 (unused).
        real ll_type4 = 0;
        int s = freq_start[d];
        int len = freq_len[d];
        for (i in s:(s + len - 1)) {
          ll_type4 += freq_count[i] * rounded_freq_loglik_fun(freq_value[i], dist_type, loc, phi, kappa, resolution[d]);
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
            ll_type5 += freq_count[i] * rounded_freq_loglik_fun(freq_lower[i], dist_type, loc, phi, kappa, resolution[d]);
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
            dist_type, loc, phi, kappa, resolution[d]
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
