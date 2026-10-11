# burn_in option ("dccmidas" default vs "exact") and the burn_in output field.
# burn_in 选项(默认 "dccmidas" / "exact")与输出字段 burn_in 的测试。

# ---------------------------------------------------------------------------
# Shared fixture for the internal-level checks / 内部函数层面的共享夹具
# ---------------------------------------------------------------------------
bi_setup <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    d <- sim_data()
    set.seed(77)
    K_c <- 6; N_c <- 20
    x  <- align_sample(d$r1, K = K_c, macro = d$macro, period = "monthly")
    TT <- length(x); N <- 2
    L  <- chol(matrix(c(1, 0.5, 0.5, 1), 2))
    res  <- matrix(rnorm(TT * N), TT, N) %*% L
    mv_c <- midas_lag_matrix(x, d$macro, K = K_c, period = "monthly")
    D_t  <- array(0, dim = c(N, N, TT))
    for (tt in seq_len(TT)) diag(D_t[, , tt]) <- abs(rnorm(N, 1, 0.1))
    cache <<- list(TT = TT, N = N, K_c = K_c, N_c = N_c,
                   res = res, mv_c = mv_c, D_t = D_t)
    cache
  }
})

# Data bundle; t0 = NULL omits the field (fallback to K_c)
# 构造数据束;t0 = NULL 时不含该字段(回退到 K_c)
mk_dc <- function(s, rc, n_X, t0 = NULL) {
  dc <- list(res = s$res, n_assets = s$N, TT = s$TT, K_c = s$K_c, N_c = s$N_c,
             rc = rc, X_mats = if (n_X) list(s$mv_c) else list(),
             D_mat = matrix(1, 1, s$TT), lag_fun = "Beta",
             idx = gmxdcc:::dcc_par_index(s$N, rc = rc, n_X = n_X, n_regimes = 1))
  if (!is.null(t0)) dc$t0 <- t0
  dc
}

par_rc  <- c(0.03, 0.9, 2)                          # a, b, w2_RC
par_x   <- c(0.03, 0.9, 0.2, 0.3, 1.5)              # a, b, delta, theta_X, w2_X
par_rcx <- c(0.03, 0.9, 2, 0.2, 0.3, 1.5, 0.4)      # + theta_RC

# ---------------------------------------------------------------------------
# 1. Intermediate output: per-observation likelihood start / 中途输出:似然起点
# ---------------------------------------------------------------------------
test_that("per-observation likelihood is zero exactly over the burn-in", {
  s <- bi_setup()
  cases <- list(
    list(lab = "dccmidas, RC",      rc = TRUE,  nX = 0, par = par_rc, t0 = s$K_c),
    list(lab = "exact, RC",         rc = TRUE,  nX = 0, par = par_rc, t0 = s$N_c + s$K_c),
    list(lab = "dccmidas, X only",  rc = FALSE, nX = 1, par = par_x,  t0 = s$K_c),
    list(lab = "exact, X only",     rc = FALSE, nX = 1, par = par_x,  t0 = 1L),
    list(lab = "exact, RC+X",       rc = TRUE,  nX = 1, par = par_rcx, t0 = s$N_c + s$K_c)
  )
  for (cs in cases) {
    ll <- gmxdcc:::dcc_midas_loglik(cs$par, mk_dc(s, cs$rc, cs$nX, cs$t0))
    expect_length(ll, s$TT)
    expect_true(all(ll[seq_len(cs$t0)] == 0), label = cs$lab)
    expect_true(ll[cs$t0 + 1] != 0, label = cs$lab)
    expect_true(all(is.finite(ll)), label = cs$lab)
  }
})

test_that("a bundle without t0 falls back to K_c (legacy equivalence)", {
  s <- bi_setup()
  for (cfg in list(list(TRUE, 0, par_rc), list(FALSE, 1, par_x))) {
    expect_identical(
      gmxdcc:::dcc_midas_loglik(cfg[[3]], mk_dc(s, cfg[[1]], cfg[[2]])),
      gmxdcc:::dcc_midas_loglik(cfg[[3]], mk_dc(s, cfg[[1]], cfg[[2]], t0 = s$K_c))
    )
  }
})

# ---------------------------------------------------------------------------
# 2. Intermediate output: what the first estimated R_bar is built from
#    中途输出:首个估计日的长期相关由哪些已实现相关构成
# ---------------------------------------------------------------------------
test_that("exact burn-in: first R_bar uses only genuine realized correlations", {
  s <- bi_setup()
  K_c <- s$K_c; N_c <- s$N_c
  C    <- midas_rc_rcpp(N_c, s$N, s$TT, s$res)
  Rbar <- gmxdcc:::build_longrun_corr(par_rc, mk_dc(s, TRUE, 0))

  # Filter weights over (C_{t-K_c}, ..., C_t), oldest first, zero on C_t
  # 滤波权重按 (C_{t-K_c},...,C_t) 由旧到新排列,当天权重为 0
  fw <- gmxdcc:::midas_filter_weights(K_c, par_rc[3], "Beta")

  # First estimated day under "exact" (1-based) / "exact" 下首个估计日
  s_ex <- N_c + K_c + 1
  win  <- (s_ex - K_c):s_ex
  expect_equal(Rbar[1, 2, s_ex], sum(fw * C[1, 2, win]), tolerance = 1e-12)
  expect_true(all(C[1, 2, win[-length(win)]] != 0))     # all K_c lags genuine

  # Under "dccmidas" the first estimated day leans on zero-filled C_t
  # "dccmidas" 下首个估计日用到了填 0 的已实现相关
  s_dm   <- K_c + 1
  win_dm <- (s_dm - K_c):(s_dm - 1)
  expect_true(any(C[1, 2, win_dm] == 0))
})

# ---------------------------------------------------------------------------
# 3. Intermediate output: placeholder matrices / 中途输出:占位矩阵
# ---------------------------------------------------------------------------
test_that("dcc_matrices keeps placeholders exactly over the burn-in", {
  s <- bi_setup()
  S <- stats::cov(s$res)
  for (cfg in list(list(TRUE, 0, par_rc, s$K_c),
                   list(TRUE, 0, par_rc, s$N_c + s$K_c),
                   list(FALSE, 1, par_x, 1L))) {
    t0 <- cfg[[4]]
    m  <- gmxdcc:::dcc_matrices(cfg[[3]], "DCC-MIDAS", mk_dc(s, cfg[[1]], cfg[[2]], t0), s$D_t)
    for (tt in seq_len(t0)) {
      expect_identical(m$R_t[, , tt], diag(s$N))
      expect_equal(m$H_t[, , tt], S, tolerance = 1e-15)
    }
    first <- t0 + 1
    expect_true(m$R_t[1, 2, first] != 0)
    expect_equal(m$H_t[, , first],
                 s$D_t[, , first] %*% m$R_t[, , first] %*% s$D_t[, , first],
                 tolerance = 1e-12)
  }
})

# ---------------------------------------------------------------------------
# 4. The two options differ only by initialization / 两选项只差初始化
# ---------------------------------------------------------------------------
test_that("dccmidas vs exact likelihoods converge after the burn-in", {
  s <- bi_setup()
  t_ex <- s$N_c + s$K_c
  ll_dm <- gmxdcc:::dcc_midas_loglik(par_rcx, mk_dc(s, TRUE, 1, s$K_c))
  ll_ex <- gmxdcc:::dcc_midas_loglik(par_rcx, mk_dc(s, TRUE, 1, t_ex))
  d <- abs(ll_ex - ll_dm)[(t_ex + 1):s$TT]
  expect_gt(d[1], 1e-8)                       # differ right after the start
  expect_lt(max(utils::tail(d, 200)), 1e-6)   # Q recursion forgets the start
})

# ---------------------------------------------------------------------------
# 5. fit_dcc: option handling / fit_dcc:选项启动
# ---------------------------------------------------------------------------
test_that("fit_dcc validates the burn_in option", {
  skip_on_cran()
  pn <- dcc_panel()
  expect_error(fit_dcc(pn$panel, correlation = "cDCC", burn_in = "foo"))
  expect_error(
    fit_dcc(pn$panel["2018/2019"], correlation = "DCC-MIDAS", rc = TRUE,
            N_c = 600, K_c = 6, burn_in = "exact"),
    "leaves no in-sample days"
  )
})

test_that("default burn_in reproduces the previous behaviour exactly", {
  skip_on_cran()
  pn <- dcc_panel()
  a <- suppressWarnings(fit_dcc(pn$panel, correlation = "cDCC"))
  b <- suppressWarnings(fit_dcc(pn$panel, correlation = "cDCC", burn_in = "dccmidas"))
  expect_identical(coef(a), coef(b))
  expect_identical(a$H_t, b$H_t)

  m1 <- suppressWarnings(fit_dcc(pn$panel, correlation = "DCC-MIDAS", rc = FALSE,
                                 corr_X = list(M = pn$macro), K_c = 6,
                                 n_starts = 2, seed = 1))
  m2 <- suppressWarnings(fit_dcc(pn$panel, correlation = "DCC-MIDAS", rc = FALSE,
                                 corr_X = list(M = pn$macro), K_c = 6,
                                 n_starts = 2, seed = 1, burn_in = "dccmidas"))
  expect_identical(coef(m1), coef(m2))
  expect_identical(as.numeric(m1$c_llk), as.numeric(m2$c_llk))
  expect_identical(m1$R_t, m2$R_t)
  expect_identical(m1$burn_in, 6L)
})

# ---------------------------------------------------------------------------
# 6. fit_dcc end to end: every model x option / 端到端:全部模型 x 选项
# ---------------------------------------------------------------------------
test_that("burn_in field and placeholders are correct for every model", {
  skip_on_cran()
  pn <- dcc_panel()
  N_c <- 20; K_c <- 6
  grid <- list(
    list(cm = "cDCC", opt = "dccmidas", Kc = NULL, exp = 1L),
    list(cm = "cDCC", opt = "dccmidas", Kc = 5,    exp = 5L),
    list(cm = "cDCC", opt = "exact",    Kc = 5,    exp = 1L),
    list(cm = "aDCC", opt = "dccmidas", Kc = NULL, exp = 1L),
    list(cm = "aDCC", opt = "exact",    Kc = NULL, exp = 1L),
    list(cm = "DECO", opt = "dccmidas", Kc = NULL, exp = 1L),
    list(cm = "DECO", opt = "exact",    Kc = NULL, exp = 1L),
    list(cm = "X",    opt = "dccmidas", exp = K_c),
    list(cm = "X",    opt = "exact",    exp = 1L),
    list(cm = "RC",   opt = "dccmidas", exp = K_c),
    list(cm = "RC",   opt = "exact",    exp = N_c + K_c),
    list(cm = "RCXb", opt = "dccmidas", exp = K_c),
    list(cm = "RCXb", opt = "exact",    exp = N_c + K_c)
  )

  for (g in grid) {
    lab <- paste(g$cm, g$opt)
    fit <- suppressWarnings(switch(g$cm,
      X    = fit_dcc(pn$panel, correlation = "DCC-MIDAS", rc = FALSE,
                     corr_X = list(M = pn$macro), K_c = K_c,
                     n_starts = 2, seed = 1, burn_in = g$opt),
      RC   = fit_dcc(pn$panel, correlation = "DCC-MIDAS", rc = TRUE,
                     N_c = N_c, K_c = K_c, burn_in = g$opt),
      RCXb = fit_dcc(pn$panel, correlation = "DCC-MIDAS", rc = TRUE, N_c = N_c,
                     corr_X = list(M = pn$macro), K_c = K_c,
                     corr_break_dates = as.Date("2016-01-01"),
                     n_starts = 2, seed = 1, burn_in = g$opt),
      fit_dcc(pn$panel, correlation = g$cm, K_c = g$Kc, burn_in = g$opt)
    ))

    b <- fit$burn_in
    expect_identical(b, as.integer(g$exp), label = lab)
    expect_identical(fit$spec$burn_in, g$opt, label = lab)
    expect_true(all(is.finite(coef(fit))), label = lab)
    expect_true(is.finite(as.numeric(fit$c_llk)), label = lab)
    expect_identical(fit$maxlik$code, 0L, label = lab)

    # placeholders before, genuine correlations after / 之前占位,之后为真实相关
    for (tt in seq_len(b)) expect_identical(fit$R_t[, , tt], diag(2), label = lab)
    expect_true(fit$R_t[1, 2, b + 1] != 0, label = lab)
    after <- fit$R_t[, , (b + 1):dim(fit$R_t)[3], drop = FALSE]
    expect_true(all(is.finite(after)), label = lab)
    expect_true(all(abs(after) <= 1 + 1e-10), label = lab)
    expect_equal(after[1, 1, ], rep(1, dim(after)[3]), tolerance = 1e-12)
  }
})

test_that("cDCC: per-observation likelihood matches the reported c_llk", {
  skip_on_cran()
  pn <- dcc_panel()
  for (cfg in list(list("dccmidas", 5), list("exact", 5))) {
    fit <- suppressWarnings(fit_dcc(pn$panel, correlation = "cDCC",
                                    K_c = cfg[[2]], burn_in = cfg[[1]]))
    kern <- if (cfg[[1]] == "exact") NULL else cfg[[2]]
    ll <- new_dcc_loglik(coef(fit), fit$eps_t, kern)
    expect_true(all(ll[seq_len(fit$burn_in)] == 0))
    expect_equal(sum(ll), as.numeric(fit$c_llk), tolerance = 1e-6)
  }
})

test_that("out-of-sample segment is unaffected by the burn-in", {
  skip_on_cran()
  pn <- dcc_panel()
  oos <- 100
  fit <- suppressWarnings(fit_dcc(
    pn$panel, correlation = "DCC-MIDAS", rc = TRUE, N_c = 20,
    corr_X = list(M = pn$macro), K_c = 6,
    out_of_sample = oos, n_starts = 2, seed = 1, burn_in = "exact"
  ))
  expect_identical(fit$burn_in, 26L)
  expect_identical(dim(fit$R_t_oos), as.integer(c(2, 2, oos)))
  expect_true(all(is.finite(fit$R_t_oos)))
  expect_true(all(fit$R_t_oos[1, 2, ] != 0))     # no placeholder in OOS
  expect_output(summary(fit), "Burn-in: first 26")
})
