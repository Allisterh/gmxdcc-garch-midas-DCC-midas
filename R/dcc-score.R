# Supported configuration of the analytic DCC-MIDAS score.
# 解析得分支持的设定:二元、只含外生协变量(rc = FALSE)、至少 1 个协变量。
dcc_score_supported <- function(data) {
  isTRUE(data$n_assets == 2) && !isTRUE(data$rc) && length(data$X_mats) >= 1L
}

dcc_score_check <- function(data) {
  if (!dcc_score_supported(data)) {
    stop("analytic gradient is available for two assets, exogenous covariates only (rc = FALSE)",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Shared intermediate quantities of the second-step score (eqs. (3.1)-(3.6)).
# The recursion starts where dcc_midas_loglik() starts it:
# t0 = data$t0 %||% data$K_c (Q_{t0} = I, derivatives zero up to t0);
# the MIDAS lag grid uses data$K_c.
# 第二步得分的共享中间量;递推起点与 dcc_midas_loglik() 相同
# (t0 = data$t0 %||% data$K_c),MIDAS 滞后格点用 data$K_c。
dcc_score_parts <- function(param, data) {
  dcc_score_check(data)
  idx <- data$idx; p <- idx$n_par; TT <- data$TT; K_c <- data$K_c
  t0  <- data$t0 %||% data$K_c
  D   <- data$D_mat; R <- nrow(D)
  a <- unname(param[idx$a]); b <- unname(param[idx$b])

  # Long-run correlation from the package builder; off-diagonal element
  # 长期相关取自包内构造器 build_longrun_corr 的非对角元
  rbar <- as.numeric(build_longrun_corr(param, data)[2, 1, ])

  # d z_t / d phi (eq. (3.5)); R == 1: pair-specific delta/theta,
  # R > 1: regime-varying m/theta. w2 is pair-specific in both cases.
  # z 对参数的导数 (3.5):无断点为逐对 delta/theta,有断点为区间 m/theta。
  dz <- matrix(0, TT, p)
  if (R == 1L) dz[, idx$delta] <- 1 else dz[, idx$m] <- t(D)
  for (j in seq_along(data$X_mats)) {
    wd   <- midas_weight_derivs(K_c, param[idx$w2[[j]]], data$lag_fun)
    F_j  <- midas_filter_sum(data$X_mats[[j]], wd$phi)
    dF_j <- midas_filter_sum(data$X_mats[[j]], wd$dphi)
    if (R == 1L) {
      slope_t <- rep(unname(param[idx$theta[[j]]]), TT)
      dz[, idx$theta[[j]]] <- F_j
    } else {
      slope_t <- as.vector(param[idx$theta[[j]]] %*% D)
      dz[, idx$theta[[j]]] <- t(D) * F_j
    }
    dz[, idx$w2[[j]]] <- slope_t * dF_j
  }
  drbar <- (1 - rbar^2) * dz                                   # d tanh(z) = (1 - tanh^2) dz

  # Q recursion on t = t0+1..TT (MODEL-SPEC sec. 2), Q_{t0} = I
  # Q 递推,t = t0+1..TT,Q_{t0} = I
  eta  <- data$res
  e_ii <- eta[, 1]^2; e_mm <- eta[, 2]^2; e_im <- eta[, 1] * eta[, 2]
  tt <- (t0 + 1):TT; pv <- tt - 1
  qrec <- function(x, init) {
    out <- rep(init, TT)
    out[tt] <- as.numeric(stats::filter(x, b, method = "recursive", init = init))
    out
  }
  q_ii <- qrec((1 - a - b) + a * e_ii[pv], 1)
  q_mm <- qrec((1 - a - b) + a * e_mm[pv], 1)
  q_im <- qrec((1 - a - b) * rbar[tt] + a * e_im[pv], 0)
  rho  <- q_im / sqrt(q_ii * q_mm)

  # Direct terms of the derivative recursions (3.4a)-(3.4d)
  # 导数递推的直接项 (3.4a)–(3.4d);对角元不依赖长期参数
  lam <- setdiff(seq_len(p), c(idx$a, idx$b))
  n <- length(tt)
  Xim <- matrix(0, n, p); Xii <- matrix(0, n, p); Xmm <- matrix(0, n, p)
  Xim[, idx$a] <- -rbar[tt] + e_im[pv]                         # (3.4a)
  Xii[, idx$a] <- -1 + e_ii[pv]
  Xmm[, idx$a] <- -1 + e_mm[pv]
  Xim[, idx$b] <- -rbar[tt] + q_im[pv]                         # (3.4b)
  Xii[, idx$b] <- -1 + q_ii[pv]
  Xmm[, idx$b] <- -1 + q_mm[pv]
  Xim[, lam]   <- (1 - a - b) * drbar[tt, lam, drop = FALSE]   # (3.4c)
  drec <- function(X) {
    out <- matrix(0, TT, p)
    out[tt, ] <- apply(X, 2, function(x)
      as.numeric(stats::filter(x, b, method = "recursive", init = 0)))
    out
  }
  dq_im <- drec(Xim); dq_ii <- drec(Xii); dq_mm <- drec(Xmm)

  d <- 1 - rho^2
  drho <- dq_im / sqrt(q_ii * q_mm) - 0.5 * rho / q_ii * dq_ii -
    0.5 * rho / q_mm * dq_mm                                   # (3.3)
  dl_drho <- rho / d + (e_im * (1 + rho^2) - rho * (e_ii + e_mm)) / d^2   # (3.1)
  list(t0 = t0, tt = tt, rho = rho, drho = drho, dl_drho = dl_drho)
}

#' Analytic per-observation score of the DCC-MIDAS correlation likelihood
#'
#' Score of the bivariate DCC-MIDAS quasi-likelihood computed by
#' [dcc_midas_loglik()], with the first-step standardized residuals held
#' fixed: \eqn{s_t = (\partial \ell_t / \partial \rho_t)\, \partial \rho_t / \partial \phi},
#' where \eqn{\rho_t = q_{12,t} / \sqrt{q_{11,t} q_{22,t}}} and the derivatives
#' of the three distinct `Q_t` elements follow linear recursions (the
#' diagonal elements depend on `a`, `b` only). Available for two assets,
#' exogenous covariates only (`rc = FALSE`, at least one covariate), with or
#' without breaks: without breaks the long-run parameters are the pair
#' coefficients `delta`, `theta_<X>`; with breaks the regime coefficients
#' `m_r`, `theta_<X>_r`; `w2_<X>` in both cases. Rows up to the recursion
#' start (`data$t0`, default `data$K_c`) are zero, matching the zero rows of
#' [dcc_midas_loglik()].
#' 二元 DCC-MIDAS 相关拟似然的解析逐观测得分(第一步残差视为常数);
#' 只支持两资产、外生协变量(rc = FALSE);有无断点均可;递推起点之前的
#' 行为 0,与 dcc_midas_loglik 的零值行一致。
#'
#' @param param Parameter vector in [dcc_par_index()] order.
#' @param data Data bundle as in [build_longrun_corr()] / [dcc_midas_loglik()]
#'   (`res`, `n_assets`, `TT`, `K_c`, optional `t0`, `rc`, `X_mats`, `D_mat`,
#'   `lag_fun`, `idx`).
#'
#' @return `TT x p` matrix of per-observation scores with canonical
#'   parameter names; column sums are the gradient of
#'   `sum(dcc_midas_loglik(param, data))`.
#' @references Engle, R. F. (2002). Dynamic conditional correlation: A simple
#'   class of multivariate generalized autoregressive conditional
#'   heteroskedasticity models. *Journal of Business & Economic Statistics*
#'   20(3), 339–350.
#'
#'   Colacito, R., Engle, R. F. and Ghysels, E. (2011). A component model
#'   for dynamic correlations. *Journal of Econometrics* 164(1), 45–59.
#' @seealso [dcc_expected_info()], [fit_dcc()]
#' @export
dcc_midas_score <- function(param, data) {
  P <- dcc_score_parts(param, data)
  s <- P$dl_drho * P$drho                                      # (3.6)
  s[seq_len(P$t0), ] <- 0
  colnames(s) <- data$idx$names
  s
}

#' Expected information matrix of the DCC-MIDAS correlation likelihood
#'
#' \eqn{I_T = \sum_{t > t_0} (1 + \rho_t^2)/(1 - \rho_t^2)^2\,
#' (\partial \rho_t/\partial \phi)(\partial \rho_t/\partial \phi)'}
#' (eq. (3.7) of the derivation note; the factor is the Fisher information
#' of a bivariate standard normal correlation). Same scope as
#' [dcc_midas_score()].
#' 第二步期望信息矩阵(推导文 (3.7)),适用范围同 dcc_midas_score。
#'
#' @inheritParams dcc_midas_score
#' @return `p x p` symmetric matrix with parameter names as dimnames.
#' @seealso [dcc_midas_score()]
#' @export
dcc_expected_info <- function(param, data) {
  P <- dcc_score_parts(param, data); tt <- P$tt
  w <- (1 + P$rho[tt]^2) / (1 - P$rho[tt]^2)^2                 # (3.2)
  I <- crossprod(P$drho[tt, , drop = FALSE] * sqrt(w))         # (3.7)
  dimnames(I) <- list(data$idx$names, data$idx$names)
  I
}
