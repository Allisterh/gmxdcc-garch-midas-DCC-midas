# Analytic gradients (gm_score, dcc_midas_score), expected information
# matrices, and the analytic_grad option of fit_garch_midas() / fit_dcc().
# 解析梯度、期望信息矩阵与 analytic_grad 选项的测试。
# sim_data() / dcc_panel() are defined in helper-legacy.R / helper-panel.R.

# ---------------------------------------------------------------------------
# Fixtures / 夹具
# ---------------------------------------------------------------------------

# First-step data bundle and spec, built exactly as inside fit_garch_midas()
# 第一步数据束与参数规格,构造方式与 fit_garch_midas() 内部一致
ag_gm_setup <- function(returns, X, K, rv, break_dates, lag_fun) {
  prep <- gmxdcc:::gm_prepare_data(returns, X = X, K = K, rv = rv,
                                   period = "monthly", break_dates = break_dates)
  idx <- gmxdcc:::gm_par_index(length(prep$lag_mats), nrow(prep$D_mat),
                               "norm", prep$cov_names)
  data <- list(daily_ret = prep$returns, lag_mats = prep$lag_mats, K = prep$K,
               D_mat = prep$D_mat, lag_fun = lag_fun, distribution = "norm",
               idx = idx)
  spec <- gmxdcc:::gm_param_spec(
    n_cov = length(prep$lag_mats), n_regimes = nrow(prep$D_mat),
    distribution = "norm", lag_fun = lag_fun, cov_names = prep$cov_names,
    mu_start = mean(prep$returns)
  )
  list(data = data, spec = spec)
}

# Estimated ("formal") points: optimum of three short multistart runs
# 正式点:三个种子的短多起点估计值
ag_formal_points <- function(loglik_fn, data, spec, seeds = 1:3, n_starts = 2) {
  lapply(seeds, function(s) {
    opt <- gmxdcc:::run_multistart(loglik_fn, data, spec,
                                   n_starts = n_starts, seed = s)
    stats::coef(opt$est)
  })
}

# Random points inside the constraints, w2 drawn too. Slopes on the
# realized-volatility covariate (monthly RV here ranges ~6-160) are drawn in
# a narrow range: U(-1, 1) slopes give log tau around -130 and a likelihood
# of -5e17, where the finite-difference reference itself fails.
# 随机点:在约束内抽取,w2 也随机;RV 斜率取窄区间(RV 量级 6–160,
# U(-1,1) 斜率会使 log tau ≈ -130、似然 ≈ -5e17,数值差分本身失准)。
ag_random_points <- function(spec, w2_idx, w2_range, seed, n = 3, mu_idx = NULL,
                             narrow_idx = NULL) {
  set.seed(seed)
  B <- gmxdcc:::draw_starts(spec, n)
  B[, w2_idx] <- stats::runif(n * length(w2_idx), w2_range[1], w2_range[2])
  if (!is.null(mu_idx)) B[, mu_idx] <- B[, mu_idx] + stats::rnorm(n, 0, 0.05)
  if (length(narrow_idx)) {
    B[, narrow_idx] <- stats::runif(n * length(narrow_idx), -0.01, 0.01)
  }
  lapply(seq_len(n), function(i) B[i, ])
}

# Total and per-observation analytic score vs maxLik::numericGradient
# (eps = 1e-6); tolerance 1e-4 * max(1, |s|) per component.
# 总得分与逐观测得分对照数值梯度;逐分量容差 1e-4 * max(1, |s|)。
ag_check_score <- function(loglik_fn, score_fn, par, data, label) {
  S_ana <- score_fn(par, data)
  S_num <- maxLik::numericGradient(function(p) loglik_fn(p, data), t0 = par, eps = 1e-6)
  s_ana <- colSums(S_ana); s_num <- colSums(S_num)
  expect_identical(dim(S_ana), dim(S_num), label = label)
  expect_true(all(abs(s_ana - s_num) <= 1e-4 * pmax(1, abs(s_num))),
              label = sprintf("%s: total score, max |diff| = %.3g", label,
                              max(abs(s_ana - s_num))))
  expect_true(all(abs(S_ana - S_num) <= 1e-4 * pmax(1, abs(S_num))),
              label = sprintf("%s: per-observation score, max |diff| = %.3g", label,
                              max(abs(S_ana - S_num))))
}

# Second-step fixture: bivariate residuals simulated from a DCC-MIDAS process
# whose long-run correlation is driven by the macro covariate (K_c = 6).
# 第二步夹具:按 DCC-MIDAS 过程模拟二元残差,长期相关由宏观变量驱动。
ag_dcc_setup <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    d <- sim_data(); K_c <- 6
    x  <- align_sample(d$r1, K = K_c, macro = d$macro, period = "monthly")
    TT <- length(x); days <- as.Date(zoo::index(x))
    X_mat <- midas_lag_matrix(x, d$macro, K = K_c, period = "monthly")
    F_t <- gmxdcc:::midas_filter_sum(X_mat, gmxdcc:::midas_filter_weights(K_c, 3, "Beta"))
    rbar <- tanh(0.3 + 0.4 * F_t)
    set.seed(4242)
    a <- 0.05; b <- 0.9
    eta <- matrix(0, TT, 2); q <- c(1, 1, 0)       # q_11, q_22, q_12
    eta[1, ] <- stats::rnorm(2)
    for (t in 2:TT) {
      q <- (1 - a - b) * c(1, 1, rbar[t]) +
        a * c(eta[t - 1, 1]^2, eta[t - 1, 2]^2, eta[t - 1, 1] * eta[t - 1, 2]) + b * q
      rho <- q[3] / sqrt(q[1] * q[2])
      z <- stats::rnorm(2)
      eta[t, ] <- c(z[1], rho * z[1] + sqrt(1 - rho^2) * z[2])
    }
    cache <<- list(TT = TT, K_c = K_c, days = days, X_mat = X_mat, eta = eta)
    cache
  }
})

ag_dc_bundle <- function(s, break_dates = NULL, t0 = NULL, rc = FALSE) {
  D_mat <- build_break_dummies(s$days, break_dates)
  spec  <- gmxdcc:::dcc_param_spec(2, rc = rc, n_X = 1, n_regimes = nrow(D_mat),
                                   lag_fun = "Beta", X_names = "M")
  dc <- list(res = s$eta, n_assets = 2L, TT = s$TT, K_c = s$K_c, N_c = 20,
             rc = rc, X_mats = list(s$X_mat), D_mat = D_mat, lag_fun = "Beta",
             idx = spec$idx)
  if (!is.null(t0)) dc$t0 <- t0
  list(dc = dc, spec = spec)
}

ag_breaks <- as.Date(c("2014-07-01", "2017-01-01"))

# Two-asset panel from a genuine two-step DGP: bivariate DCC-MIDAS
# innovations (long-run correlation driven by the macro covariate, K_c = 6)
# scaled by asset-specific GARCH-MIDAS volatilities (K = 12). A
# well-identified surface keeps the comparison of optimizer end points
# about the gradient, not about flat directions.
# 两步 DGP 面板:二元 DCC-MIDAS 新息 × 各资产 GARCH-MIDAS 波动;
# 似然面可识别,比较停点时不受平坦方向干扰。
ag_two_step_panel <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    d   <- sim_data()
    cal <- align_sample(xts::xts(rep(0, length(d$r1)), order.by = zoo::index(d$r1)),
                        K = 12, macro = d$macro, period = "monthly")
    TT  <- length(cal)
    F12 <- gmxdcc:::midas_filter_sum(midas_lag_matrix(cal, d$macro, K = 12, period = "monthly"),
                                     gmxdcc:::midas_filter_weights(12, 3, "Beta"))
    F6  <- gmxdcc:::midas_filter_sum(midas_lag_matrix(cal, d$macro, K = 6, period = "monthly"),
                                     gmxdcc:::midas_filter_weights(6, 3, "Beta"))
    rbar <- tanh(0.3 + 0.4 * F6)
    set.seed(5151)
    a <- 0.04; b <- 0.92
    eta <- matrix(0, TT, 2); q <- c(1, 1, 0)
    eta[1, ] <- stats::rnorm(2)
    for (t in 2:TT) {
      q <- (1 - a - b) * c(1, 1, rbar[t]) +
        a * c(eta[t - 1, 1]^2, eta[t - 1, 2]^2, eta[t - 1, 1] * eta[t - 1, 2]) + b * q
      rho <- q[3] / sqrt(q[1] * q[2])
      z <- stats::rnorm(2)
      eta[t, ] <- c(z[1], rho * z[1] + sqrt(1 - rho^2) * z[2])
    }
    gm_sim <- function(e, m, th, al, be, mu) {
      tau <- exp(m + th * F12); g <- numeric(TT); g[1] <- 1
      for (t in 2:TT) g[t] <- (1 - al - be) + al * g[t - 1] * e[t - 1]^2 + be * g[t - 1]
      mu + sqrt(g * tau) * e
    }
    r1 <- gm_sim(eta[, 1], 0.2, 0.5, 0.06, 0.85, 0.05)
    r2 <- gm_sim(eta[, 2], 0.0, 0.3, 0.05, 0.88, 0.03)
    panel <- xts::xts(cbind(A1 = r1, A2 = r2), order.by = zoo::index(cal))
    cache <<- list(panel = panel, macro = d$macro)
    cache
  }
})

# ---------------------------------------------------------------------------
# 1. gm_score vs numerical gradient / 第一步得分对照数值梯度
# ---------------------------------------------------------------------------
test_that("gm_score matches the numerical gradient (Beta, Almon, X + RV)", {
  skip_on_cran()
  d <- sim_data()
  cases <- list(
    list(lab = "Beta, X, 2 breaks", X = list(M = d$macro), rv = FALSE,
         br = ag_breaks, lf = "Beta", w2r = c(1.5, 6)),
    list(lab = "Almon, X, 1 break", X = list(M = d$macro), rv = FALSE,
         br = ag_breaks[2], lf = "Almon", w2r = c(-0.05, -0.002)),
    list(lab = "Beta, RV + X, 1 break", X = list(M = d$macro), rv = TRUE,
         br = ag_breaks[1], lf = "Beta", w2r = c(1.5, 6))
  )
  for (cs in cases) {
    st <- ag_gm_setup(d$r_dgp, cs$X, K = 12, rv = cs$rv,
                      break_dates = cs$br, lag_fun = cs$lf)
    idx <- st$data$idx
    rv_theta <- if (cs$rv) idx$theta[[which(st$data$idx$names[idx$w2] == "w2_RV")]] else NULL
    pts <- c(ag_formal_points(gmxdcc:::gm_loglik, st$data, st$spec),
             ag_random_points(st$spec, idx$w2, cs$w2r, seed = 31, mu_idx = idx$mu,
                              narrow_idx = rv_theta))
    expect_length(pts, 6)
    for (i in seq_along(pts)) {
      ag_check_score(gmxdcc:::gm_loglik, gm_score, pts[[i]], st$data,
                     sprintf("%s, point %d", cs$lab, i))
    }
    expect_identical(colnames(gm_score(pts[[1]], st$data)), idx$names)
  }
})

test_that("gm_score refuses Student-t innovations", {
  d <- sim_data()
  st <- ag_gm_setup(d$r_dgp, list(M = d$macro), K = 12, rv = FALSE,
                    break_dates = NULL, lag_fun = "Beta")
  dat <- st$data
  dat$distribution <- "std"
  expect_error(gm_score(rep(0.1, 7), dat), "norm")
  expect_error(gm_expected_info(rep(0.1, 7), dat), "norm")
})

# ---------------------------------------------------------------------------
# 2. dcc_midas_score vs numerical gradient / 第二步得分对照数值梯度
# ---------------------------------------------------------------------------
test_that("dcc_midas_score matches the numerical gradient (breaks / no breaks)", {
  skip_on_cran()
  s <- ag_dcc_setup()
  for (br in list(NULL, ag_breaks)) {
    lab <- if (is.null(br)) "no breaks" else "2 breaks"
    b <- ag_dc_bundle(s, break_dates = br)
    idx <- b$dc$idx
    pts <- c(ag_formal_points(gmxdcc:::dcc_midas_loglik, b$dc, b$spec),
             ag_random_points(b$spec, unlist(idx$w2), c(1.2, 6), seed = 32))
    expect_length(pts, 6)
    for (i in seq_along(pts)) {
      ag_check_score(gmxdcc:::dcc_midas_loglik, dcc_midas_score, pts[[i]], b$dc,
                     sprintf("DCC-MIDAS %s, point %d", lab, i))
    }
    # Zero rows up to the recursion start, as in dcc_midas_loglik()
    # 递推起点之前为零行,与 dcc_midas_loglik() 一致
    S <- dcc_midas_score(pts[[1]], b$dc)
    expect_true(all(S[seq_len(s$K_c), ] == 0))
    expect_identical(colnames(S), idx$names)
  }
})

test_that("dcc_midas_score follows the burn_in recursion start (t0 field)", {
  skip_on_cran()
  s <- ag_dcc_setup()
  b <- ag_dc_bundle(s, break_dates = ag_breaks, t0 = 1L)
  par <- ag_random_points(b$spec, unlist(b$dc$idx$w2), c(1.2, 6), seed = 33, n = 1)[[1]]
  ag_check_score(gmxdcc:::dcc_midas_loglik, dcc_midas_score, par, b$dc, "t0 = 1")
  S <- dcc_midas_score(par, b$dc)
  expect_true(all(S[1, ] == 0))
  expect_true(any(S[2, ] != 0))
})

test_that("dcc_midas_score refuses unsupported configurations", {
  s <- ag_dcc_setup()
  b <- ag_dc_bundle(s, rc = TRUE)
  expect_error(dcc_midas_score(rep(0.1, b$spec$n_par), b$dc),
               "two assets, exogenous covariates only")
  b2 <- ag_dc_bundle(s)
  b2$dc$n_assets <- 3L
  expect_error(dcc_expected_info(rep(0.1, b2$spec$n_par), b2$dc),
               "two assets, exogenous covariates only")
})

# ---------------------------------------------------------------------------
# 3. Expected information: symmetric, positive definite / 期望信息对称正定
# ---------------------------------------------------------------------------
test_that("gm_expected_info and dcc_expected_info are symmetric positive definite", {
  skip_on_cran()
  d <- sim_data()
  st <- ag_gm_setup(d$r_dgp, list(M = d$macro), K = 12, rv = TRUE,
                    break_dates = ag_breaks, lag_fun = "Beta")
  for (par in ag_formal_points(gmxdcc:::gm_loglik, st$data, st$spec, seeds = 1)) {
    I1 <- gm_expected_info(par, st$data)
    expect_equal(dim(I1), rep(st$spec$n_par, 2))
    expect_true(isSymmetric(I1, tol = 1e-12))
    expect_gt(min(eigen(I1, symmetric = TRUE, only.values = TRUE)$values), 0)
  }
  s <- ag_dcc_setup()
  for (br in list(NULL, ag_breaks)) {
    b <- ag_dc_bundle(s, break_dates = br)
    for (par in ag_formal_points(gmxdcc:::dcc_midas_loglik, b$dc, b$spec, seeds = 1)) {
      I2 <- dcc_expected_info(par, b$dc)
      expect_equal(dim(I2), rep(b$spec$n_par, 2))
      expect_true(isSymmetric(I2, tol = 1e-12))
      expect_gt(min(eigen(I2, symmetric = TRUE, only.values = TRUE)$values), 0)
    }
  }
})

# ---------------------------------------------------------------------------
# 4. analytic_grad = TRUE reaches the same optimum / 解析梯度拟合到达同一最优点
# ---------------------------------------------------------------------------
# Coefficient tolerance 1e-3 * max(1, |coef|): absolute 1e-3 for
# coefficients below 1 in magnitude, relative 1e-3 above. The relative
# part matters for w2 only. The MIDAS weight parameter is the flattest
# direction of the likelihood, and BFGS (reltol 1e-8) stops some way short
# of the optimum along it with either gradient (asset 2 of the two-step
# panel: both fits stop near w2 = 12.47 and differ there by 0.003, while
# the polished optimum is 11.82 with a log-likelihood higher by only 9e-4).
# All other coefficients agree to about 1e-5.
# 系数容差 1e-3 * max(1, |coef|):|coef| < 1 时为绝对 1e-3,否则为相对 1e-3,
# 实际只影响 w2。w2 是似然最平坦的方向,两种梯度下 BFGS 都在其上提前停止
# (两步面板资产 2:两者都停在 w2 ≈ 12.47、相差 0.003,加精最优 11.82,
# 似然只高 9e-4)。
ag_coef_close <- function(e1, e0, label) {
  expect_true(all(abs(e1 - e0) <= 1e-3 * pmax(1, abs(e0))),
              label = sprintf("%s: max |diff| = %.3g (%s)", label, max(abs(e1 - e0)),
                              names(e0)[which.max(abs(e1 - e0))]))
}

test_that("fit_garch_midas(analytic_grad = TRUE) matches the default fit", {
  skip_on_cran()
  d <- sim_data()
  args <- list(d$r_dgp, X = list(M = d$macro), K = 12, rv = FALSE,
               break_dates = ag_breaks[2], distribution = "norm",
               lag_fun = "Beta", n_starts = 3, seed = 11)
  f0 <- do.call(fit_garch_midas, args)
  f1 <- do.call(fit_garch_midas, c(args, list(analytic_grad = TRUE)))
  expect_lte(abs(f1$loglik - f0$loglik), 1e-2)
  ag_coef_close(f1$maxlik$estimate, f0$maxlik$estimate, "fit_garch_midas")
  expect_true(f1$spec$analytic_grad)
  expect_false(f0$spec$analytic_grad)
  expect_error(do.call(fit_garch_midas,
                       utils::modifyList(args, list(distribution = "std",
                                                    analytic_grad = TRUE))),
               "norm")
})

test_that("fit_dcc(analytic_grad = TRUE) matches the default fit", {
  skip_on_cran()
  pn <- ag_two_step_panel()
  args <- list(pn$panel, univariate = "garch_midas",
               univariate_args = list(X = list(M = pn$macro), K = 12, rv = FALSE),
               correlation = "DCC-MIDAS", rc = FALSE,
               corr_X = list(M = pn$macro), K_c = 6,
               n_starts = 3, seed = 5)
  f0 <- suppressWarnings(do.call(fit_dcc, args))
  f1 <- suppressWarnings(do.call(fit_dcc, c(args, list(analytic_grad = TRUE))))
  for (i in 1:2) {
    expect_lte(abs(f1$u_llk[[i]] - f0$u_llk[[i]]), 1e-2)
    ag_coef_close(f1$est_univ_model_origin[[i]]$maxlik$estimate,
                  f0$est_univ_model_origin[[i]]$maxlik$estimate,
                  sprintf("fit_dcc step 1, asset %d", i))
  }
  expect_lte(abs(as.numeric(f1$c_llk) - as.numeric(f0$c_llk)), 1e-2)
  ag_coef_close(f1$maxlik$estimate, f0$maxlik$estimate, "fit_dcc step 2")
  expect_identical(f1$spec$analytic_grad_used, c(univariate = TRUE, correlation = TRUE))
})

# ---------------------------------------------------------------------------
# 5. Default unchanged / 默认行为不变
# ---------------------------------------------------------------------------
test_that("analytic_grad = FALSE is bit-identical to omitting the argument", {
  skip_on_cran()
  pn <- dcc_panel()
  args <- list(pn$panel["2014/2019"], univariate = "garch_midas",
               univariate_args = list(X = pn$macro, K = 6, rv = FALSE),
               correlation = "DCC-MIDAS", rc = FALSE,
               corr_X = list(M = pn$macro), K_c = 6,
               corr_break_dates = ag_breaks[2],
               n_starts = 2, seed = 3)
  f_omit  <- suppressWarnings(do.call(fit_dcc, args))
  f_false <- suppressWarnings(do.call(fit_dcc, c(args, list(analytic_grad = FALSE))))
  expect_identical(f_false, f_omit)
  expect_false(f_omit$spec$analytic_grad)
  expect_identical(f_omit$spec$analytic_grad_used, c(univariate = FALSE, correlation = FALSE))
})

test_that("non-DCC-MIDAS correlation falls back with one message", {
  skip_on_cran()
  pn <- dcc_panel()
  f0 <- suppressWarnings(fit_dcc(pn$panel, univariate = "sGARCH",
                                 correlation = "cDCC", distribution = "norm"))
  expect_message(
    f1 <- suppressWarnings(fit_dcc(pn$panel, univariate = "sGARCH",
                                   correlation = "cDCC", distribution = "norm",
                                   analytic_grad = TRUE)),
    "numerical gradient"
  )
  # rugarch first step + cDCC second step: nothing analytic is used, so the
  # estimates are those of the default fit. 无解析步骤,估计值与默认相同。
  expect_identical(f1$maxlik$estimate, f0$maxlik$estimate)
  expect_identical(f1$spec$analytic_grad_used, c(univariate = FALSE, correlation = FALSE))
})
