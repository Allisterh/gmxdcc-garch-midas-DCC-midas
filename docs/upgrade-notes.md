# gmxdcc 升级说明 / Upgrade notes (0.1.0 → 0.1.1)

> **TL;DR (English)**
> - Upgrading to 0.1.1 requires **no code changes**: with default arguments
>   every number is identical to 0.1.0 (verified by the test suite).
> - New: `fit_dcc(..., burn_in = "dccmidas" | "exact")` (default
>   `"dccmidas"` = old behaviour) and the output field `fit$burn_in`.
> - **Action recommended:** drop the first `fit$burn_in` days of `H_t`,
>   `R_t`, `R_t_bar` before computing conditional betas or losses — they are
>   placeholders (this was already true in 0.1.0, just undocumented).
> - Exact old code: `remotes::install_github("TLib124/gmxdcc@v0.1.0")`.

本文件写给**正在使用 gmxdcc 的其他项目**(及其维护者/助手):说明 0.1.1 改了什么、旧代码要不要改、以及如何得到与旧版完全相同的结果。

---

## 1. 先确认你装的是哪个版本

```r
packageVersion("gmxdcc")                 # 0.1.1 表示已是新版
"burn_in" %in% names(formals(gmxdcc::fit_dcc))   # TRUE 表示有新选项
```

---

## 2. 改了什么

### 0.1.1(本次)

| 变化 | 类型 | 对旧代码的影响 |
|---|---|---|
| `fit_dcc()` 新增参数 `burn_in = c("dccmidas", "exact")` | 新增,默认值保持旧行为 | **无**。不写这个参数时结果与 0.1.0 逐位相同 |
| 结果对象新增字段 `fit$burn_in`(整数)与 `fit$spec$burn_in` | 新增字段 | 无。只是多了信息 |
| `summary(fit)` 多打印一行 `Burn-in: ...` | 输出格式 | 若有代码解析 summary 的文本输出,需留意多了一行 |

### 0.1.0 发布后的修正(若你的安装早于这些修正,也会一并获得)

| 变化 | 对旧代码的影响 |
|---|---|
| `control$start_matrix` 只作用于第二步(相关性);第一步起点改由 `univariate_args$start_matrix` 注入(单个矩阵或按资产的列表) | 只影响**自己注入起点**的代码(如滚动重估热重启)。以前在 `univariate = "garch_midas"` 下注入相关步起点会报列数错误,现在正常 |
| 注入的起点若贴着约束边界,会自动微调进可行域内部 | 以前会报 `initial value in 'vmmin' is not finite`,现在正常 |

C++ 内核没有任何改动。

---

## 3. 我的旧代码需要改吗?

**不需要。** 默认 `burn_in = "dccmidas"` 与 0.1.0 完全一致(测试套件中有断言:不写参数与显式写 `"dccmidas"` 的系数、似然、`H_t`、`R_t` 逐位相同;所有对照旧研究脚本的金标准测试照常通过)。

**但建议加一步:剔除占位天。** 相关递推需要起始期,输出数组不截短,而是在开头保留占位值:

- `R_t = 单位阵`(条件相关恒为 0);
- `H_t = 标准化残差的协方差`,**不是收益尺度**。实测包内置数据上占位天的 `H_t[1,2]` 约为 0.89,而正常天约为 1e-5,差了 5 个数量级;
- 似然贡献为 0。

占位天在 0.1.0 中就存在(默认 DCC-MIDAS 为前 `K_c` 天),只是之前没有字段告诉你。**如果你以前直接用全长的 `H_t` 算条件 beta、损失函数或画图,开头那几十天混入了占位值,建议回头检查。**

```r
b    <- fit$burn_in
keep <- (b + 1):dim(fit$H_t)[3]
H    <- fit$H_t[, , keep]
beta <- H[1, 2, ] / H[2, 2, ]          # 按你的资产顺序调整下标
days <- fit$Days[keep]
```

0.1.0 下没有 `fit$burn_in` 字段时,占位天数可按下表自己算(默认约定):DCC-MIDAS 为 `K_c`;cDCC/aDCC/DECO 为 1(若传了 `K_c` 则为 `K_c`)。

---

## 4. `burn_in` 两个选项的区别

| 模型 | `"dccmidas"`(默认) | `"exact"` |
|---|---|---|
| DCC-MIDAS | 占位 `K_c` 天 | 含 RC 项:`N_c + K_c` 天;仅 X 协变量:1 天 |
| cDCC / aDCC / DECO | 1 天(传了 `K_c` 则 `K_c` 天) | 1 天 |

- `"dccmidas"`:沿用 dccmidas 包(Candila)和旧研究脚本的约定,便于和文献、旧结果对照。
- `"exact"`:更严格。含 RC 项时,RC 的日度已实现相关前 `N_c` 天填 0,默认约定下第 `K_c+1` 到 `N_c+K_c` 天的长期相关部分由这些 0 构成并进入似然;`"exact"` 把这段也跳过。仅 X 协变量时 `K_c` 的单位是**月**,日历对齐后从第 1 天起就完整,`"exact"` 只跳 1 天。

**注意:** 两个选项计入似然的天数不同,`c_llk` 不能直接拿来比较谁更好;在真实数据上两者的参数估计也可能有明显差别(结构突变模型的似然面多峰)。建议把 `"exact"` 当作稳健性检验报告,并适当加大 `n_starts`。

---

## 5. 想得到与旧版完全相同的结果

按需求从轻到重:

1. **什么都不做**:装 0.1.1,不写 `burn_in`,数值与 0.1.0 相同。
2. **显式写出来**(推荐写在正式脚本里,意图更清楚):
   ```r
   fit_dcc(..., burn_in = "dccmidas")
   ```
3. **精确锁定旧代码**(例如论文复现包要求代码一字不差):
   ```r
   remotes::install_github("TLib124/gmxdcc@v0.1.0", force = TRUE, upgrade = "never")
   packageVersion("gmxdcc")   # 应显示 0.1.0
   ```
   回到新版:
   ```r
   remotes::install_github("TLib124/gmxdcc@v0.1.1", force = TRUE, upgrade = "never")
   ```

安装前请先重启 R 会话(RStudio:Session → Restart R)。Windows 下已加载的包无法被覆盖安装,会出现"正在使用,因此不会被安装"的警告且实际没有装上。

---

## 6. 参考

- 完整说明:`docs/gmxdcc-manual.md`(中文)、`docs/gmxdcc-manual-en.md`(英文),
  第 4 节"烧入期(burn-in)与占位天"
- 变更记录:`NEWS.md`
- 帮助页:`?fit_dcc`(参数 `burn_in`、返回值说明)
