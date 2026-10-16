#' MIDAS filter weights and their derivative with respect to w2
#'
#' Returns the rolling-sum weight vector used by [build_tau()] (identical to
#' [midas_filter_weights()]: reversed, current period zero-padded, `w1`
#' fixed at 1 for Beta and 0 for Almon) together with its derivative with
#' respect to `w2`. The raw weights live on the `K + 1` grid of
#' `midas_weights(1:(K + 1), K + 1, w1, w2)`; writing them as
#' \eqn{\pi_k = w_k / \sum_{j=1}^{K+1} w_j}, the derivative is
#' \eqn{\partial \pi_k / \partial w_2 = \pi_k (c_k - \bar c)},
#' \eqn{\bar c = \sum_{j=1}^{K+1} \pi_j c_j}, with
#' \eqn{c_k = \partial \log w_k / \partial w_2 = \log(1 - k/(K+1))} (Beta) or
#' \eqn{k^2} (Almon). For Beta the `K + 1` term carries zero weight, so the
#' sum effectively runs over `1..K` (eq. (1.5) of the derivation note). For
#' Almon the `K + 1` term is part of the normalizing sum and is then dropped
#' from the filter, exactly as in [midas_filter_weights()].
#' 返回与 midas_filter_weights 完全相同格点的滤波权重及其对 w2 的导数;
#' 归一化在 K+1 格点上进行(Beta 的第 K+1 项权重为 0;Almon 的第 K+1 项
#' 进入分母后被丢弃,与包内滤波一致)。
#'
#' @param K MIDAS lag order.
#' @param w2 Second weighting parameter.
#' @param lag_fun `"Beta"` or `"Almon"`.
#'
#' @return List with `phi` and `dphi`, numeric vectors of length `K + 1` in
#'   the order expected by `roll::roll_sum()` (oldest lag first, trailing 0).
#' @keywords internal
midas_weight_derivs <- function(K, w2, lag_fun) {
  w1   <- ifelse(lag_fun == "Beta", 1, 0)
  k    <- 1:(K + 1)
  # Raw weights on the K + 1 grid, same call as midas_filter_weights()
  # 与 midas_filter_weights() 相同的调用,K+1 格点上的原始归一化权重
  pi_k <- midas_weights(k, (K + 1), w1, w2, lag_fun)
  c_k  <- if (lag_fun == "Beta") log(1 - k / (K + 1)) else k^2
  # Beta: the K + 1 term has pi = 0 and c = -Inf; exclude it from the mean
  # Beta 的第 K+1 项 pi = 0、c = -Inf,不进入加权平均
  ok    <- pi_k > 0
  c_bar <- sum(pi_k[ok] * c_k[ok])
  dpi   <- numeric(K + 1)
  dpi[ok] <- pi_k[ok] * (c_k[ok] - c_bar)
  ord <- function(v) c(rev(v)[2:(K + 1)], 0)   # 反转 + 当期补 0
  list(phi = ord(pi_k), dphi = ord(dpi))
}

#' Complete-window weighted sum of a MIDAS lag matrix
#'
#' Same call as in [build_tau()]: `roll::roll_sum()` with the given weights,
#' keeping row `K + 1` (the single complete window). Weights may be
#' negative (derivative weights).
#' 与 build_tau 相同的加权滚动和,取第 K+1 行;权重可为负(导数权重)。
#'
#' @param lag_mat `(K + 1) x TT` lag matrix.
#' @param wvec Weight vector of length `K + 1`.
#' @return Numeric vector of length `TT`.
#' @keywords internal
midas_filter_sum <- function(lag_mat, wvec) {
  K1 <- length(wvec)
  if (nrow(lag_mat) != K1) {
    stop("midas_filter_sum: lag matrix rows must equal the weight length.", call. = FALSE)
  }
  filt <- suppressWarnings(roll::roll_sum(lag_mat, c(K1), weights = wvec))
  as.numeric(filt[K1, ])
}

#' Derivative of the log long-run variance with respect to all parameters
#'
#' For \eqn{\log \tau_t = m(t) + \sum_j \theta_j(t) F_{j,t}(w_{2,j})}
#' (see [build_tau()]) returns the `TT x p` matrix of
#' \eqn{\partial \log \tau_t / \partial \psi}: columns `m_r` are the regime
#' dummies (eq. (1.4a)), `theta_j_r` the dummies times the MIDAS filter
#' (1.4b), `w2_j` the regime slope times the filter of the weight derivative
#' (1.4c); `alpha`, `beta`, `mu` columns are zero.
#' 长期成分对数对全部参数的导数(推导文 (1.4a)–(1.4c)),T × p 矩阵。
#'
#' @param param Parameter vector ordered as in [gm_par_index()].
#' @param data Data bundle as in [gm_loglik()].
#' @return `TT x p` numeric matrix.
#' @keywords internal
dlogtau <- function(param, data) {
  idx <- data$idx
  D   <- data$D_mat
  TT  <- ncol(D)
  out <- matrix(0, TT, idx$n_par)
  out[, idx$m] <- t(D)                                         # (1.4a)
  for (j in seq_along(data$lag_mats)) {
    wd <- midas_weight_derivs(data$K[j], param[idx$w2[j]], data$lag_fun)
    F_j  <- midas_filter_sum(data$lag_mats[[j]], wd$phi)
    dF_j <- midas_filter_sum(data$lag_mats[[j]], wd$dphi)
    theta_t <- as.vector(param[idx$theta[[j]]] %*% D)
    out[, idx$theta[[j]]] <- t(D) * F_j                        # (1.4b) 按列乘 F
    out[, idx$w2[j]]      <- theta_t * dF_j                    # (1.4c)
  }
  out
}

# Shared intermediate quantities of the first-step score (eqs. (1.2)-(1.7)).
# 第一步得分的共享中间量:g、tau、h、eps、d log h(v)。
gm_score_parts <- function(param, data) {
  if (!identical(data$distribution, "norm")) {
    stop("Analytic gradient (gm_score) is available for distribution = \"norm\" only; ",
         "use analytic_grad = FALSE for Student-t innovations.", call. = FALSE)
  }
  idx <- data$idx; p <- idx$n_par
  comp <- gm_components(param, data)
  g    <- as.numeric(zoo::coredata(comp$g_it))
  tau  <- as.numeric(zoo::coredata(comp$tau_d))
  h    <- g * tau
  r    <- as.numeric(zoo::coredata(data$daily_ret))
  TT   <- length(r)
  al   <- unname(param[idx$alpha]); be <- unname(param[idx$beta])
  eps  <- r - unname(param[idx$mu])
  zeta <- eps^2 / tau
  dlt  <- dlogtau(param, data)

  # Direct terms of the g recursion (1.7); row 1 is zero (g_1 = 1 fixed).
  # g 递推 (1.7) 的直接项;第 1 行为 0(g_1 = 1 为常数)。
  X <- matrix(0, TT, p)
  if (TT > 1) {
    prev <- 1:(TT - 1); cur <- 2:TT
    lam  <- c(idx$w2, idx$m, unlist(idx$theta))
    X[cur, idx$alpha] <- -1 + zeta[prev]                                  # (1.7a)
    X[cur, idx$beta]  <- -1 + g[prev]                                     # (1.7b)
    X[cur, lam]       <- -al * zeta[prev] * dlt[prev, lam, drop = FALSE]  # (1.7c)
    X[cur, idx$mu]    <- -2 * al * eps[prev] / tau[prev]                  # (1.7d)
  }
  # Column-wise recursion y_t = x_t + beta * y_{t-1} (stats::filter in C)
  # 逐列递推 y_t = x_t + beta * y_{t-1}
  dg <- apply(X, 2, function(x) as.numeric(stats::filter(x, be, method = "recursive")))
  dg <- matrix(dg, TT, p)
  v  <- dg / g + dlt                                                      # (1.3)
  list(h = h, eps = eps, v = v)
}

#' Analytic per-observation score of the GARCH-MIDAS likelihood
#'
#' Bollerslev–Wooldridge-type score of the Gaussian GARCH-MIDAS
#' quasi-likelihood computed by [gm_loglik()]:
#' \deqn{s_t = \tfrac12 (\epsilon_t^2/h_t - 1)\, v_t + (\epsilon_t/h_t)\, e_\mu,
#'   \qquad v_t = \partial \log h_t / \partial \psi,}
#' with \eqn{v_t = \partial \log g_t + \partial \log \tau_t}, the long-run part
#' in closed form and the short-run part by the recursion of Fiorentini,
#' Calzolari and Panattoni (1996) extended by the GARCH-MIDAS term
#' \eqn{-\alpha (\epsilon_{t-1}^2/\tau_{t-1})\, \partial \log \tau_{t-1}}.
#' Supports any number of covariates (including realized volatility), any
#' number of regimes, and both Beta and Almon weights; Gaussian innovations
#' only. Pass it to `maxLik` as `grad` (see the `analytic_grad` argument of
#' [fit_garch_midas()]).
#' GARCH-MIDAS 正态拟似然的解析逐观测得分;支持任意协变量个数(含 RV)、
#' 任意区间数、Beta 与 Almon 权重;只支持正态分布(std 时报错)。
#'
#' @param param Numeric parameter vector ordered as in [gm_par_index()].
#' @param data Data bundle as in [gm_loglik()] (`daily_ret`, `lag_mats`,
#'   `K`, `D_mat`, `lag_fun`, `distribution`, `idx`).
#'
#' @return `TT x p` matrix of per-observation scores, columns named by the
#'   canonical parameter names. Its column sums are the gradient of
#'   `sum(gm_loglik(param, data))`.
#' @references Bollerslev, T. and Wooldridge, J. M. (1992). Quasi-maximum
#'   likelihood estimation and inference in dynamic models with
#'   time-varying covariances. *Econometric Reviews* 11(2), 143–172.
#'
#'   Fiorentini, G., Calzolari, G. and Panattoni, L. (1996). Analytic
#'   derivatives and the computation of GARCH estimates. *Journal of
#'   Applied Econometrics* 11(4), 399–417.
#' @seealso [gm_expected_info()], [fit_garch_midas()]
#' @export
gm_score <- function(param, data) {
  P <- gm_score_parts(param, data); idx <- data$idx
  s <- 0.5 * (P$eps^2 / P$h - 1) * P$v                         # (1.2) 第一项
  s[, idx$mu] <- s[, idx$mu] + P$eps / P$h                     # (1.2) mu 项
  colnames(s) <- idx$names
  s
}

#' Expected information matrix of the GARCH-MIDAS likelihood
#'
#' Sum over observations of the conditional expected negative Hessian under
#' correct specification of the Gaussian model (eq. (2.2) of the
#' derivation note): \eqn{I_T = \sum_t [\tfrac12 v_t v_t' + h_t^{-1} e_\mu e_\mu']}.
#' Needs first derivatives only.
#' 期望信息矩阵(推导文 (2.2)),只用一阶导数。
#'
#' @inheritParams gm_score
#' @return `p x p` symmetric matrix with parameter names as dimnames.
#' @seealso [gm_score()]
#' @export
gm_expected_info <- function(param, data) {
  P <- gm_score_parts(param, data); idx <- data$idx
  I <- 0.5 * crossprod(P$v)                                    # (2.2) 第一项
  I[idx$mu, idx$mu] <- I[idx$mu, idx$mu] + sum(1 / P$h)        # + sum(1/h) e_mu e_mu'
  dimnames(I) <- list(idx$names, idx$names)
  I
}
