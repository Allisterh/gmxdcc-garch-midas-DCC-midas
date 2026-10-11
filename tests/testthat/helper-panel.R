# Correlated 2-asset panel with GARCH dynamics in the first column, shared by
# test-dcc-fit.R and test-burn-in.R.
# 两资产相关面板(第一列含真实 GARCH 动态),供多个测试文件共享。
dcc_panel <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    d <- sim_data()
    set.seed(505)
    r1 <- d$r_dgp
    z  <- xts::xts(rnorm(length(r1), 0, 0.8), order.by = zoo::index(r1))
    r2 <- 0.6 * r1 + z                     # correlated second asset
    panel <- xts::merge.xts(A1 = r1, A2 = r2)
    colnames(panel) <- c("A1", "A2")
    cache <<- list(panel = panel, macro = d$macro)
    cache
  }
})
