# Helmholtz 实用自适应 LOD：论文数值实验最小实施计划

> 状态日期：2026-08-12
>
> 对应论文：`../LOD_paper/helmholtz_lod_certified_amsart.tex`，重点对应 Practical simplified adaptive loop 和 Adaptive numerical experiments。
>
> 目标：在 LOD-C++ 中实现论文的核心自适应流程，并以尽可能少的重复计算完成可发表的数值验证。
>
> 结果定位：默认运行是 **practical adaptive LOD**。除非所有常数和谱量都采用严格上、下界验证，否则结果只能称为 practical/conditional，不能标成 fully certified。

### 冻结的论文口径与版本基线

本计划与论文数值章节统一采用以下口径：

1. R1 仅在 \(\kappa=16\) 做小规模校准；R2a 和 S 在 \(\kappa=16\) 做五方法
   主比较，并只对 R2a/S 的代表性三方法补 \(\kappa=8,32\)。
2. 生产表只列 practical adaptive LOD 及四个基线，不列 fully certified adaptive
   LOD。每步不重复计算 reference inf--sup、\(\widehat\delta_{\mathrm{loc}}\)、
   稳定性余量或分立理论常数。
3. ambient shadow 初始化每个 reference epoch 只执行一次；一次 \(H\)-加密及 ambient
   ratio 修复之后，状态机返回 localization check，而不是重新进入 epoch 初始化。
4. \(\theta_H=0.5\) 固定使用，完整 Dörfler 参数敏感性不属于生产实验。

论文源文件当前未纳入 Git。正式实现使用下列不可变基线：

```text
file: ../LOD_paper/helmholtz_lod_certified_amsart.tex
sha256: 03d83e0eb7128aa5ef00002c6dac110f548351e52e92e56d7adf709880854d20
frozen: 2026-08-11
```

同一份哈希记录保存在
`experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE.sha256`。论文源文件每次修改后
都必须生成新哈希；已有实验的 `run.json` 继续保留原哈希，不得用新哈希覆盖旧结果。

## 1. 总体取舍

本计划不继续扩展旧的 `fine/mixed/macro` 三分支，也不在每步计算
\(C_F,C_\pi,C_a,c_W,C_{\mathrm{sd}},C_{\mathrm{ov}}\) 和参考离散 inf--sup 常数。
数值实验只保留论文算法的两个核心量：

\[
\eta_H \quad\text{和}\quad
\Theta_{\mathrm{loc}}^{\mathrm{amb}\to\mathrm{ref}}.
\]

- \(\eta_H\) 估计固定参考空间中的误差
  \(\lVert u_h-U\rVert_\kappa\)，并驱动粗网格 \(H\)-加密。
- \(\Theta_{\mathrm{loc}}^{\mathrm{amb}\to\mathrm{ref}}\) 检查局部化 corrector
  是否充分；不满足阈值时只执行全局 \(\ell\leftarrow\ell+1\)。
- ambient 网格只是跟随 \(H\)-加密的 shadow fine mesh，用于 corrector 扰动证书；
  它既不定义误差目标，也不触发 reference/fine 加密。
- 正式主实验尽量在一个 reference epoch 内完成。reference refresh 仅作为显式、独立的
  兜底操作，不属于自适应内循环。

这样可以删除旧计划中最昂贵、但对论文核心结论并非必要的部分：逐步验证全部理论常数、
每个容差重新运行、多个 Dörfler 参数全扫描、局部 \(\ell_T\)、corrector fine-mesh
自适应，以及每步导出完整场数据。

## 2. 数值空间和误差合同

一个 reference epoch 内固定

\[
V_H\subset V_h^{\mathrm{ref}}\subset V_{\mathrm{amb}}\subset V.
\]

三张网格必须在程序和输出中使用不同名称，禁止继续把 `fine`、`audit` 和
`reference` 混为一谈。

| 论文对象 | 建议代码对象 | 是否可在内循环改变 | 用途 |
|---|---|---:|---|
| \(\mathcal T_H\) | `coarse_mesh` | 是 | LOD 粗网格和 \(H\)-标记 |
| \(\mathcal T_h^{\mathrm{ref}}\) | `reference_mesh` | 否 | corrector、参考解 \(u_h\)、真实误差目标 |
| \(\mathcal T_{h,\mathrm{amb}}\) | `ambient_mesh` | 是 | shadow 网格和局部化扰动证书 |

epoch 开始时设置

\[
\mathcal T_{h,\mathrm{amb}}:=\mathcal T_h^{\mathrm{ref}}.
\]

每次 \(H\)-加密后，只加密对应区域的 ambient 网格并做相容闭包，直到

\[
\rho_{\mathrm{amb}}
=\max_{T\in\mathcal T_H}
  \max_{K\subset T,\,K\in\mathcal T_{h,\mathrm{amb}}}
  \frac{h_K}{H_T}
\le \rho_\star.
\]

实现必须计算实际的直径比并检查上界；不能假设 `h/H` 严格等于某个常数。
默认可把 \(\rho_\star=1/4\) 作为初始值，最终值由一次小规模校准确定。

论文中

\[
\lVert Q_\infty^*-Q_{\infty,h}^*\rVert
\le C_s\rho_h^s
\]

只说明统一控制网格比可以避免 ideal-corrector 离散误差在 \(H\)-自适应中恶化；
固定网格比本身不代表该误差收敛。数值报告不得作更强的结论。

### 2.1 reference refresh

只有在粗网格继续加密会破坏 \(V_H\subset V_h^{\mathrm{ref}}\)，或者外部的连续误差
策略明确要求时，才允许结束当前 epoch：

\[
\mathcal T_h^{\mathrm{ref}}\leftarrow
\mathcal T_{h,\mathrm{amb}},\qquad
V_h^{\mathrm{ref}}\leftarrow V_{\mathrm{amb}}.
\]

随后重新计算 \(u_h\)，epoch 编号加一，并把 ambient 网格重新初始化为 reference
网格。由于误差目标已经改变，不同 epoch 的误差曲线不能无标记地连接。
\(\Theta_{\mathrm{loc}}\) 超阈值绝不能触发 reference refresh。

## 3. 最小自适应算法

### 3.1 手动设置的参数

生产运行只接受下列参数：

| 参数 | 建议初值 | 含义 |
|---|---:|---|
| `c_H` | `0.5` | schema-v3 兼容的粗分辨率诊断值；不触发 practical driver 加密 |
| `theta_loc` | `0.19985934547160381` | R1 校准后的 corrector 局部化证书阈值 |
| `C0_usr` | `0.49135072057990814` | R1 observed effectivity 的 practical 基础系数 |
| `C1_usr` | `0.059576473311033412` | R1 corrector 扰动的 practical 附加系数 |
| `theta_H` | `0.5` | \(\eta_{H,T}\) 的 Dörfler 参数 |
| `rho_star` | `0.25` | ambient/coarse 最大直径比 |
| `ell0` | `2` | 初始全局 oversampling 层数 |
| `ell_max` | `6` | 防止异常无限增加的工作上限 |
| `max_H_steps` | `10` | 单条自适应轨迹的粗加密上限 |
| `trajectory_policy` | `fixed_work_horizon`（E1 取轨迹） | `practical_indicator` 仅用于独立算法停止研究 |
| `practical_stop_tolerance` | 独立研究参数 | (U_{prac}) 的绝对停止阈值，不用于 E1 论文曲线截断 |
| plateau policy | 三点窗：几何平均 ratio `>=0.9`，波动 `<=15%` | 只读 reference-error 平台诊断，不反馈 MARK/STOP |

这些数值是启动配置，不是理论常数。R1 校准后，将同一网格族、插值算子、边界条件和
波数区间使用的 `theta_loc/C0_usr/C1_usr` 固定下来，主实验中不得逐算例调参。

停止量为

\[
U_{\mathrm{prac}}
=\bigl(C_0^{\mathrm{usr}}
      +C_1^{\mathrm{usr}}\Theta_{\mathrm{loc}}^{\mathrm{amb}\to\mathrm{ref}}\bigr)
  \eta_H.
\]

若用户常数没有严格验证，`U_prac` 只叫 practical stopping indicator，不叫
rigorous upper bound。

### 3.2 状态机

每个 reference epoch 严格执行以下顺序：

```text
输入: T_H, T_h^ref, ell, theta_loc, C0_usr, C1_usr,
      theta_H, rho_star, trajectory_policy, practical_stop_tolerance, work limits
输出: T_H, T_h,amb, ell, U, eta_H, Theta_loc, 完整迭代日志

1. T_h,amb <- T_h^ref。
2. 记录 max_T kappa*H_T 但不以先验 c_H 触发加密；同步 ambient shadow，
   直到 rho_amb <= rho_star。预渐近区只由已有 reference-error 轨迹后验识别。
3. 在固定 reference 空间计算 reference localized correctors 和 Theta_loc。
4. 若 Theta_loc > theta_loc：
      ell <- ell + 1（全局）；
      只重算受影响 correctors 与 ambient Riesz 问题；返回 3。
   禁止执行 reference/fine refinement。
5. 解 Petrov--Galerkin LOD 问题，计算 eta_{H,T}、eta_H 和 U_prac。
6. `practical_indicator` 算法研究模式若 U_prac <= practical_stop_tolerance：结束；
   E1 与校准的 `fixed_work_horizon` 模式忽略该阈值，直到预先冻结的 H-step 上限。
7. 按 theta_H 对 eta_{H,T} 做 Dörfler 标记，加密 T_H 并做相容闭包。
8. 验证 V_H subset V_h^ref；若失败，返回 ReferenceRefreshRequired，
   不得在粗加密函数内部偷偷改变 reference 网格。
9. 只更新 ambient shadow，直到 rho_amb <= rho_star；返回 3。
```

必须写状态转换测试，保证任何 `Theta_loc > theta_loc` 的下一动作都是
`IncreaseGlobalEll`，从状态枚举中删除/禁用 `RefineCorrectorFine` 分支。
论文相对 reference-energy 目标始终是事后评价量，不得作为在线停止阈值。校准模式可
报告相邻误差比、负对数改善和 DOF 归一化斜率。schema-v4 固定查看最后三个误差点：
两步几何平均比值位于 `[0.9,1]` 且窗口 `(max-min)/min<=15%` 时标记
`plateau_observed`；允许一次小幅反弹，但强下降不得误判为平台。

## 4. 当前项目差距和代码工作包

当前项目已有 NVB、三层网格雏形、LOD corrector/PG solve、局部 kernel Riesz、
Dörfler 标记、论文算例注册表和 JSON 配置框架。这些应复用。以下缺口必须先补齐，
否则旧 `bench_helmholtz_certified` 的输出不能冒充论文新算法。

### WP1：冻结 reference，ambient 单独跟随粗网格

**状态：已完成（2026-08-11）。** 新增独立 `ReferenceEpochHierarchy`，旧
`AdaptiveMeshHierarchy` 保持 legacy 行为。粗网格 refinement 先在候选对象上验证
`coarse -> reference` 嵌入，失败返回 `ReferenceRefreshRequired` 且不提交部分状态；
ambient ratio 修复同样在候选 shadow mesh 上完成后一次性提交。G0 Release 基线为
tests-only 30/30，加入本工作包后为 31/31；启用 benchmarks/smoke 的最终完整
Release 回归为 38/38。

涉及文件：

- `include/helmholtz/adaptive/hierarchy.h`
- `src/helmholtz/adaptive/hierarchy.cpp`
- `tests/test_helmholtz_adaptive.cpp`
- `tests/test_helmholtz_reference_epoch_hierarchy.cpp`

当前 `AdaptiveMeshHierarchy::refine()` 会在粗网格逼近 fine 网格时自动细化 fine
网格。这与固定 reference epoch 冲突。建议增加 reference-epoch 模式，或新增
`ReferenceEpochHierarchy`，明确提供：

```cpp
refine_coarse_preserving_reference(marked_H);
enforce_ambient_ratio(rho_star);
reference_embedding_holds();
refresh_reference_from_ambient();
```

验收：

1. 连续三次局部 \(H\)-加密后，`reference_mesh` 的节点、单元和 generation 完全不变。
2. `ambient_mesh` 只在违反比值的 coarse 父区域及其闭包内改变。
3. 每步均有 `rho_amb <= rho_star + tolerance`。
4. 如果 coarse/reference 嵌入将失效，返回结构化状态而非自动细化 reference。

实现记录：连续三次局部 (H)-加密保持 reference 几何、层级和 generation 不变；
ambient 使用实际三角形直径和 element-parent 映射检查比值，并保留未受影响区域；容量
耗尽测试验证失败路径完全事务化。显式 `refresh_reference_from_ambient()` 会推进 epoch，
随后原先被拒绝的 coarse refinement 可继续。单位方形和正式 mixed-boundary L 型网格
均通过 prolongation composition、坐标重构、(I_H P_H) right-inverse 和边界测度检查。

### WP2：把现有 kernel residual 泛化成两种用途

**状态：已完成（2026-08-11）。** 共享的局部 kernel Riesz 求解器现可按
`KernelRieszSpace` 显式选择固定 reference 或 ambient shadow；两条路径返回不同结果
类型，ambient 结果不含 `eta_H` 标记接口。历史 audit/certified 路径保留为兼容回归，
不作为新 practical driver 的默认接口。

涉及文件：

- `include/helmholtz/adaptive/kernel_residual.h`
- `src/helmholtz/adaptive/kernel_residual.cpp`
- `tests/test_helmholtz_kernel_residual.cpp`

现有代码以 audit kernel 命名。应抽象成可选择离散空间的局部 saddle-point Riesz
装配：

- `ReferenceResidualRiesz`：在 \(W_h^{\mathrm{ref}}\) 中计算
  \(\eta_{H,z}\) 和 \(\eta_H\)，这是误差估计器和标记量；
- `AmbientDefectRiesz`：在 \(W_{\mathrm{amb},z}\) 中计算
  \(\Theta_{\mathrm{loc}}\) 所需的局部 Gram 贡献。

不得把 ambient Riesz 量直接改名为 reference estimator。验收时分别检查 kernel
约束、Riesz 恒等式、局部平方和与单元分配守恒。

实现记录：`compute_reference_residual_riesz()` 仅在 `reference_mesh()`、
`reference_quasi_interpolation()` 上计算节点指标、总 η_H、守恒的单元分配和 Dörfler
集合；ambient shadow 加密后，相同 reference 输入的 η_H 保持不变。
`compute_ambient_defect_riesz()` 仅在 `ambient_mesh()` 上接收按粗基排列的 defect RHS，
返回逐 patch 的 Riesz representatives、局部 Gram 贡献及其总和 G_loc。当前 RHS 的
ambient-to-reference 回缩构造仍严格留给 WP3，本工作包没有用未回缩 residual 冒充该
defect。Saddle-point 与显式 kernel-basis 两种实现逐 patch 对照通过；kernel 约束、
Riesz 恒等式、局部平方和、Gram Hermitian/半正定性和 reference 单元分配守恒均由
`test_helmholtz_kernel_residual.cpp` 覆盖。受影响的历史 certificates/driver/backend
定向回归为 4/4 通过；启用 benchmarks/smoke 的最终 Release 全量回归为 38/38。

### WP3：实现 ambient-to-reference retraction 和新证书

**状态：已完成（2026-08-11）。** 新增独立 `reference_retraction` 模块；普通
浮点结果固定标记为 `ImplementationStudy`，不声称严格 verified certificate。旧
`Theta_total/Theta_h/delta_h/q_h` 链未接入 practical 路径。

已新增：

- `include/helmholtz/adaptive/reference_retraction.h`
- `src/helmholtz/adaptive/reference_retraction.cpp`
- `tests/test_helmholtz_reference_retraction.cpp`

按论文构造局部稳定投影 \(J_{\mathrm{ref}}:V_{\mathrm{amb}}\to
V_h^{\mathrm{ref}}\)，并实现

\[
\mathfrak R_{\mathrm{ref}}
=(I_{V_h^{\mathrm{ref}}}-P_{H,\mathrm{ref}}I_H)
J_{\mathrm{ref}}\big|_{W_{\mathrm{amb}}}.
\]

证书装配使用

\[
\mathscr D_{\mathrm{loc}}^{\mathrm{amb}\to\mathrm{ref}}(v_H)(w)
=a\!\left(\mathfrak R_{\mathrm{ref}}w,
          (I-Q_{\ell,h}^{*,\mathrm{ref}})v_H\right),
\]

而不是旧的未回缩 ambient residual。由局部 Riesz Gram 矩阵
\(G_{\mathrm{loc}}\) 和粗能量矩阵 \(M_\kappa\) 求

\[
\Theta_{\mathrm{loc}}^2
=\lambda_{\max}(G_{\mathrm{loc}},M_\kappa).
\]

首版可用稀疏 Lanczos/LOBPCG 和上一轮特征向量 warm start；仅在很小的 R1
验证网格上用稠密特征值作交叉检查。

验收：

1. 对 \(w\in W_h^{\mathrm{ref}}\)，
   \(\mathfrak R_{\mathrm{ref}}w=w\) 达到线性求解容差。
2. 任意 ambient 输入的回缩结果满足 \(I_H\mathfrak R_{\mathrm{ref}}w=0\)。
3. 增大 \(\ell\) 时证书在测试问题上下降；若用 reference ideal corrector 的
   小矩阵直接算
   \(\|Q_{\infty,h}^*-Q_{\ell,h}^*\|\)，证书给出论文要求的一侧控制。
4. 禁止退回旧 `Theta_total/Theta_h/delta_h/q_h` 链作为正式自适应输入。

实现记录：`J_ref` 使用每个自由 reference 节点选定的单参考三角形及其局部 L2 对偶
基函数装配，是 Scott--Zhang 型局部投影；Dirichlet 行显式归零。
`ReferenceRetraction` 实施
`(I_ref-P_H,ref I_H,ref) J_ref`，并检查 projector identity、reference-kernel
identity、`I_H R_ref=0` 和齐次迹。回缩 defect RHS 按论文的复数内积约定装配为
`R_ref^* A_ref^* Y_ell`，随后复用 WP2 的 `AmbientDefectRiesz` 构造 `G_loc`。
最大广义特征值使用可 warm start 的迭代求解，小维数自动与稠密解交叉检查。

G3 小矩阵检查使用 R1、κ=4（仅算法单元验证，不是 E0 的 κ=16 校准）：
`ell=1,2,3` 得到的 `Theta_loc` 依次为约
`1.26534, 0.113464, 0`；ideal reference corrector 的回缩 defect 达到机器精度。
直接 corrector 扰动为 `0.240449`，有限维定理上界为 `0.378976`，一侧方向通过；
同时得到 `C_J≈1.61411`、`C_ret≈1.66631`、`c_W≈0.838429`、
`C_sd≈0.251115`。测试还拒绝与 localized basis 不匹配的陈旧证书。证书 builder
会从当前两层网格和 PDE 系数重装 operators 并比对，且要求 ambient
子单元继承对应 reference 父单元系数；陈旧/错配 operators 同样 fail closed。WP2/WP3
及旧 certificates/driver/backend 的 Release 定向回归为 5/5 通过；启用 benchmarks/smoke
的最终 Release 全量回归为 39/39。

### WP4：新增 practical driver，不强行复用旧 certified 状态机

**状态：已完成（2026-08-11）。**

已新增：

- `include/helmholtz/adaptive/practical_driver.h`
- `src/helmholtz/adaptive/practical_driver.cpp`
- `tests/test_helmholtz_practical_driver.cpp`

状态只保留：

```text
CoarseAdmissibility
LocalizationCheck
SolveAndEstimate
RefineCoarse
Converged
ReferenceRefreshRequired
WorkLimitReached
Failed
```

旧 `certified_driver.*` 和 `numerical_backend.*` 可保留做历史回归，但不应继续
扩展成论文默认 runner。首版每次网格变化允许 full rebuild，以正确性优先；通过
等价性测试后再加入 patch cache。

执行状态（2026-08-11）：WP4/G4 已完成。独立的 `PracticalAdaptiveDriver` 已将
`ReferenceEpochHierarchy`、reference residual Riesz、ambient-to-reference
localization certificate 和 Petrov--Galerkin LOD solve 接成真实数值链；每次
$H$-变化或 global \(\ell\) 变化均执行 full rebuild。状态与上列合同一致，动作
枚举中不存在 corrector/reference fine-refinement；非退化的 `H4/reference6`
回归锁定 `Theta_loc` 超阈值后的下一动作只能为 `IncreaseGlobalEll`。reference
residual 的 Dörfler 标记会事务性加密 $H$ 并恢复 ambient ratio；固定 reference
容量不足时返回 `ReferenceRefreshRequired`，不提交部分粗网格且不隐式刷新 epoch。
工作步数、网格规模、$H$-步数和 wall time 均有结构化 `WorkLimitReached` 终止。
新增测试同时覆盖正常收敛、G4、真实 estimator 驱动的 $H$-加密和 reference
capacity exhaustion；启用 benchmarks/smoke 的 Release 全量回归为 40/40。
WP5 论文 runner、schema v2、多容差轨迹抽取和最终 CSV/JSON/VTU 输出尚未开始，
因此这里的完成只表示 practical 算法合同与状态机完成，不表示生产实验完成。

### WP5：统一论文 runner 和精简配置

**状态：已完成（2026-08-11，PALOD runner 与 v2 合同）。**

已增加可执行文件 `bench_helmholtz_adaptive_paper`，复用：

- `include/helmholtz/benchmarks/paper_cases.h`
- `src/helmholtz/benchmarks/paper_cases.cpp`
- `include/helmholtz/experiments/paper_config.h`
- `src/helmholtz/experiments/paper_config.cpp`
- `experiments/helmholtz_adaptive_paper/`

新配置版本应显式含 `reference_mesh`、`ambient_mesh`、`reference_epoch`、
`rho_star`、`theta_loc`、`C0_usr`、`C1_usr`，并将旧的 `theta_h`、`q_h`、
`refine_corrector_fine` 等字段从 practical 模式的有效合同中移除。

runner 必须支持从一条轨迹提取多个容差首次命中点。不要为
\(10^{-1},5\times10^{-2},2\times10^{-2},10^{-2}\) 各自重跑一次。

执行状态：legacy certified schema v1 保持不变；新增独立
`PracticalPaperConfig`、`schema-v2.json` 和 `output-schema-v2.json`。v2 显式冻结
`reference_mesh/ambient_mesh/reference_epoch`、`rho_star`、`theta_loc`、
`C0_usr/C1_usr`、单一 `theta_H=0.5`、solver、quadrature、容差、资源上限、
Git/build provenance 和论文 SHA-256，并严格拒绝旧 `theta_h/q_h` 等字段。
run ID 覆盖完整 canonical v2 配置。runner 在数值工作前核对版本化论文哈希，当前
允许真实 `PALOD`、`HLOD-fixed`、`SLOD`、`UFEM` 和 `AFEM` backend，任何名称都不能
由 legacy proxy 冒充。
`HLOD-fixed` 与 PALOD 共用 reference-kernel `eta_H` 和 H 标记，但固定全局 ell，且
不计算 `Theta_loc`、ambient certificate 或 ambient shadow refinement；`UFEM` 只做
一致 conforming P1 加密/求解，并把候选解延拓到同一固定 reference 空间后事后计算误差。
`SLOD` 在一致粗网格序列上固定协议常数 \(c_{\mathrm{prior}}=1\)，即
`ell=ceil(log2(kappa))`，每层重建真实 LOD；它不计算 PALOD/HLOD 的 `eta_H`、
`Theta_loc`、`U_prac` 或 ambient certificate。
`AFEM` 在当前 conforming \(P_1\) 网格上累计
\(h_T^2\|f+\kappa^2u_H\|_T^2\)、内部通量跳跃与 impedance 边界残差，使用同一
`theta_H=0.5` 做 Dörfler 标记和局部 NVB 加密，不访问 reference solution 做 MARK/STOP。

每条 PALOD 运行先独立计算一次 evaluation reference solution，driver 本身只导出
已完成 MARK/STOP 后的候选解；reference error 在 driver 返回后计算，其时间与
reference solve 一并单列并从 method time 排除。四个共同误差目标通过
`extract_practical_target_hits` 从同一 journal 提取首次命中点，未命中写
`not_reached`。每个 run 生成 `iterations.csv`、`summary.csv`、`run.json`、
`ell_history.csv` 和带 `eta_H_T`/availability 字段的 `final_mesh.vtu`；输出 claim
固定为 `implementation-study`。配置 round-trip/hash/legacy-field rejection、单轨迹
多目标提取、reference 隔离、真实 PALOD chain 与五文件产物均有回归测试；启用
benchmarks/smoke 的 Release 全量回归为 41/41。这里的 WP5 完成不表示 E1 五方法
backend 或正式生产矩阵已经完成；它们仍按后续实验阶段推进。

资源 pilot 后合同升级为 schema v3（2026-08-12）：v2 曾把绝对 (U_{prac}) 停止阈值
隐式设成最小相对 reference-error 目标 0.01，R2a 因而在相对误差约 0.0189 时提前停止。
v3 将 `practical_stop_tolerance` 与 `relative_energy_targets` 分离，并新增
`trajectory_policy=fixed_work_horizon` 的非论文校准模式。该模式继续使用同一算法标记，
但不由 reference error 或 (U_{prac}) 提前停止；reference error 只在动作已经固定后流式
评估并输出平台诊断。v2 schema、配置与已完成 pilot 保留为历史证据，不重写 run ID。

## 5. 简化后的数值实验矩阵

实验分三层推进。只有上一层通过验收才运行下一层，以避免在错误实现上消耗大规模算力。

### E0：小规模单元验证和参数校准

**状态：已完成（2026-08-11）。** 冻结参数与原始证据位于
`experiments/helmholtz_adaptive_paper/calibration/R1-k16-v1/`，可由
`experiments/helmholtz_adaptive_paper/run_e0_calibration.sh` 在 Release 构建中重建。
层级 4/6 因 κ=16 时 ambient kernel 不具 Helmholtz coercivity 被 fail closed；最终采用
固定 coarse/reference 层级 6/8，并把 ambient ratio 控制到 0.25。

使用 R1，主波数 \(\kappa=16\)。只做以下工作：

1. 在一个小的固定 \(H/h\) 层级上取 \(\ell=1,2,3,4\)，检查
   `Theta_loc` 随 \(\ell\) 的变化。
2. 在可承受的小矩阵上直接计算 ideal reference corrector 扰动，与回缩证书比较。
3. 检查 \(\eta_H\) 的 Riesz 恒等式，以及 smooth case 不产生明显虚假局部集中。
4. 依据 observed effectivity 只选定一组保守的
   `theta_loc/C0_usr/C1_usr`；不计算全部理论常数。

E0 是调试和校准，不作为主性能结果。参数一旦冻结，后续不再逐案例调优。

实现记录：\(\ell=1,2,3,4\) 的 `Theta_loc` 分别约为
1.83348、0.181690、0.0436981、0.00270569；direct corrector perturbation 均被
ambient-to-reference 一侧上界控制。局部 constraint、Riesz stationarity、能量恒等式和
单元分配残差均在约 \(10^{-15}\) 或更小，smooth R1 的局部 effectivity 分布无异常长尾。
冻结 `ell=2`、`theta_loc=0.19985934547160381`、
`C0_usr=0.49135072057990814`、`C1_usr=0.059576473311033412` 和
`rho_star=0.25`；常数来自 observed multiplier 加 1.25 安全系数，不是理论常数。
真实 practical driver smoke 记录 5 个 journal 项，严格执行
`ell=1 -> IncreaseGlobalEll -> ell=2 -> Complete`，未发生 H 加密或 reference refresh。
新增 `helmholtz_e0_R1_k16_calibration` Release gate 后，完整回归为 42/42；该 gate
独立重建 CSV/JSON 和 driver smoke，当前耗时约 150 秒。

E1 backend 准备状态（2026-08-12）：真实 `HLOD-fixed`、`SLOD`、`UFEM` 与 `AFEM` 已接入同一 practical
paper runner。HLOD-fixed 强制 `ell0==ell_max`，复用 `eta_H`/H 标记但完全跳过 PALOD
的 localization certificate 和 ambient refinement；SLOD 执行一致粗网格、固定先验 ell
的真实 LOD 轨迹；UFEM 执行一致 conforming P1 加密/求解；AFEM 用体残差、内部通量跳跃
和 impedance 边界项做 Dörfler 局部加密。四条基线路径都只在运行后将候选解与同一
reference solution 比较，都有方法专属 action、fail-closed smoke 和五文件输出；不适用
字段留空，不得写成伪零值（SLOD 仅 ell 适用；AFEM 的 `eta_H` 明确定义为标准 P1 残差
指标而非 reference-kernel Riesz 指标）。加入四条基线 smoke 后完整 Release 回归为 46/46。
这表示 E1 五方法 backend 基础设施闭合，但不表示 R2a/S 生产实验已经完成。

### E1：核心主实验

主波数统一取 \(\kappa=16\)，只保留两个真正体现自适应价值的算例：

- R2a：\(\sigma=2^{-5}\) 的局部光滑源，检验源附近和波传播方向的粗网格集中；
- S：L 型区域混合边界的 \(r^{2/3}\) 奇性，检验低正则性下的角点加密。

每个算例比较五种方法：

1. `PALOD`：本计划的 practical adaptive LOD；
2. `HLOD-fixed`：相同 \(\eta_H\) 标记，但 \(\ell\) 固定为校准值；
3. `SLOD`：准一致粗网格、冻结 \(c_{\mathrm{prior}}=1\) 的标准先验
   \(\ell=\lceil\log_2\kappa\rceil\)；
4. `UFEM`：一致加密的 conforming \(P_1\) FEM；
5. `AFEM`：标准体残差、通量跳跃和 impedance 边界项驱动的 \(P_1\) AFEM。

这 10 条轨迹是论文最重要的数据。AFEM backend 已达到共同 reference-error 和工作统计
合同；下一步必须用生产资源上限验证 R2a/S 轨迹，不能仅凭 R1 smoke 写入论文结果。

S 必须使用论文规定的 \(\Gamma_D/\Gamma_R\) 混合边界和制造解。若生产求解链仍不支持
Dirichlet 节点，S 处于 blocked 状态；全 Robin 的 L 型替代算例只能标为 surrogate，
不能写成正式 S 结果。

### E2：最小波数敏感性

不对所有配置做笛卡尔积扫描。只在 R2a 和 S 上补
\(\kappa\in\{8,32\}\)，比较：

- `PALOD`；
- `SLOD`；
- `UFEM` 或 `AFEM` 中在该算例更有代表性的一个。

加上 E1 的 \(\kappa=16\)，即可展示 \(8,16,32\) 的波数趋势。固定 \(\ell\) HLOD
和两种 FEM 不需要在所有波数重复。R2b（\(\sigma=2^{-6}\)）只在 R2a 的局部性不够
明显时补做，不是默认生产矩阵。

因此默认生产规模约为：E0 的 4 个小测试、E1 的 10 条主轨迹、E2 的 12 条补充轨迹，
总计约 26 条，而不是“所有算例 × 所有波数 × 所有方法 × 所有容差 × 所有参数”的
大规模笛卡尔积。

论文数值章节已经同步采用上述分层协议。若以后扩展实验矩阵，应先同时更新论文协议、
本计划和 manuscript hash，再启动新增计算。

## 6. 共同运行规则

### 6.1 reference 解与真实误差

- 每个 `(case, kappa, reference_epoch)` 只求一次 reference solution \(u_h\)，随后复用。
- 所有方法使用同一参考问题、数据、求积和线性求解容差。
- 正式误差为
  \[
  E_\kappa^{\mathrm{ref}}
  =\frac{\lVert u_h-U\rVert_\kappa}{\lVert u_h\rVert_\kappa},
  \qquad
  E_0^{\mathrm{ref}}
  =\frac{\lVert u_h-U\rVert_{L^2}}{\lVert u_h\rVert_{L^2}}.
  \]
- reference/exact error 只用于事后验证，不参与标记和停止。
- reference 解时间单列，不混入方法 online time；若发生 refresh，单独记录其代价。

### 6.2 一条轨迹复用多个误差目标

对每种方法运行到最小目标或工作上限，保存整条轨迹，然后事后提取首次达到

\[
10^{-1},\quad 5\times10^{-2},\quad
2\times10^{-2},\quad10^{-2}
\]

的迭代。达不到的目标记为 `not_reached`，不无限增加参考网格或 \(\ell\)。

### 6.3 资源上限

每个配置必须有：

- `ell_max`、`max_H_steps`、`max_unknowns`、`max_wall_seconds`；
- 结构化退出原因，而不是进程崩溃或静默截断；
- 开发阶段只计时一次；仅最终选入论文的曲线重复 3 次并取中位数；
- 只输出初始、两张代表性中间网格和最终网格，不在每步写大型 VTU；
- full field/checkpoint 仅在失败或显式 `--debug-artifacts` 时保存。

## 7. 最小输出合同

`iterations.csv` 每行至少包含：

```text
schema_version, case, method, kappa, run_id,
reference_epoch, iteration, action, stop_reason,
N_H, N_ref, N_amb, ell,
kappa_H_max, mu_H, rho_amb,
eta_H, Theta_loc, U_prac,
reference_energy_error, reference_L2_error,
marked_H, rebuilt_correctors, ambient_refined_elements,
time_mesh, time_corrector, time_certificate, time_solve,
time_estimator, time_total_cumulative, peak_memory_mb
```

另生成：

- `summary.csv`：每个共同误差目标的首次命中点和累计工作；
- `run.json`：完整参数、代码 git commit、论文 SHA-256、编译选项、硬件和状态；
- `final_mesh.vtu`：最终粗网格及 \(\eta_{H,T}\)；
- `ell_history.csv`：`Theta_loc`、\(\ell\) 和 global corrector rebuild 历史。

旧 `theta_total/theta_h/delta_h/q_h` 可以保留在 legacy/debug 文件中，但不得出现在
practical 主表中造成“仍有 h 分支”的误解。

## 8. 性能优化顺序

只有通过正确性 gate 后才做优化，顺序如下：

1. **reference solve 复用**：同一 case/kappa/epoch 只计算一次。
2. **ambient-defect 多右端复用**：同一 patch 的 energy/constraint saddle 矩阵只分解一次，
   批量求解全部 coarse-input 右端；必须保持 Gram、`Theta_loc` 和残差诊断不变。
3. **corrector patch cache**：只有 source element、patch、reference 子空间和插值约束
   都不变时才复用；\(\ell\) 全局增加时仅复用仍然数学等价的装配/分解。
4. **局部 H 加密复用**：只重算受 coarse patch 变化影响的 correctors 和 Riesz 问题；
   用 full rebuild 作小规模逐步对照，误差需在求解容差内。
5. **特征值 warm start**：`Theta_loc` 的最大广义特征值沿用上一轮向量；当估计值与
   `theta_loc` 已有充分间隔时提前停止迭代。
6. **稀疏 saddle-point Riesz**：生产中不显式构造整个 kernel 基；显式基只用于单元测试。
7. **输出节流**：日志写标量，场数据只保存代表性快照。

不建议为节省时间降低到单精度；Helmholtz 稳定性和证书比较统一使用 double。

当前性能实现边界（2026-08-12）：已完成跨进程 reference 解磁盘复用、单轨迹多容差事后抽取、
`Theta_loc` warm start、E0 dense hierarchy 常数去重、HLOD-fixed 对 PALOD
localization/ambient 工作的完全跳过，以及第 2 项 patch 内多右端 Riesz 复用。R2a、
κ=16 development profile（E0 诊断 `c_H`、零 H-step，不是论文结果）显示 certificate
由 10.107 秒降至 1.404 秒，method time 由 10.217 秒降至 1.509 秒，`Theta_loc`、标记和
误差在舍入范围内不变；证据见 `profiles/R2a-k16-development-v1.md`。practical driver
在 H 或 ell 改变后仍以 full rebuild 为正确性基线，尚未启用 production corrector/patch
cache，也未根据单个 development profile 调整求解器容差。生产实验不再为满足先验
`c_H` 扩大 reference 网格；完整绘制已有 reference-error 轨迹，并仅在误差随加密持续下降、
局部对数斜率已脱离明显预渐近波动的区段比较方法。若数据不支持，则标注
`pre_asymptotic`，不追加昂贵 reference refinement 来调出预期结果。第 3 项 cache 必须
与 full rebuild 逐步对照后才可启用。

runner 现按轨迹点流式计算 reference error，评估完成后立即释放 candidate，不再在 journal
中累计保存全部 reference-size 复向量；该评估时间继续从 method time 排除。磁盘 reference
cache 的 key 绑定 reference 网格、边界标签、完整装配算子、载荷、求解器格式、Git commit
和 benchmark 二进制摘要；损坏、维数或输入变化均 fail closed 为 cache miss。连续两个独立
R1 进程已验证第一条 miss、第二条 hit，算法与误差列逐字节一致。localization 最大特征值
的幂迭代在小谱隙下若未收敛，维数不超过 512 时改用带残差验收的 Hermitian 稠密 eigensolve；
更大问题仍 fail closed，避免把停滞误写成证书。

服务器资源 pilot（AMD EPYC 9554，非论文数据）比较 4/8/16 worker，三组非计时轨迹
逐字段一致。R2a wall time 为 9.32/9.17/9.24 秒，RSS 为 247/261/289 MiB；S 为
35.94/35.42/35.74 秒与 1.29/1.31/1.35 GiB，均无 swap。4 worker 相对 8 worker 只慢
约 1.5% 且占用更低，因此冻结 `PATCH_THREADS=4`；16 worker 不再使用。certificate
占 method time 约 85%--89%，后续性能工作应针对 certificate/eigensolve，不再扩大线程数。
R2a 误差 0.190/0.0958/0.0490/0.0189，每步仍改善 49%--61%，明确未平台；S 为
0.578/0.445/0.369，只有两个改善区间，不足以判平台。下一 gate 是 schema-v3 延长校准：
R2a 6 H-step、S 5 H-step，均以 4 worker 串行运行。pilot/calibration 输出不得入论文。
v2 优化合并后的 Release 全量回归为 48/48；加入 schema-v3、平台诊断与
fixed-horizon 端到端 gate 后，Release 全量回归为 49/49（2026-08-12）。

延长校准结果（提交 `f3bf6a2`）：R2a 七点误差最终为
`0.00651/0.00678/0.00605`，三点窗几何平均比约 `0.9645`、相对波动约 `12%`，
schema-v4 将其识别为经验平台；S 六个已算点最终为 `0.235/0.169/0.0856`，
几何平均比约 `0.604`，仍处于强下降区。R2a level 10→11 独立 reference audit
得到相对 reference 差 `0.09260`，是终点报告误差的 `15.29` 倍，远超冻结的 `25%`
adequacy 门槛，因此必须进入显式 epoch 1，不能把 level-10 平台当作 PALOD 收敛。
下一 gate 只运行 R2a/reference-level-11/epoch-1 六步校准与 S 的第六个 H-step；
二者仍是非论文结果。fixed horizon 达到计划步数记为 `TrajectoryComplete/success`，
只有真实时间、内存、迭代或未知量上限才记作 censored。
加入 schema-v4、非零 epoch 传播、reference adequacy 独立 gate 后，Release 全量回归为
50/50（2026-08-12）。

2026-08-12 的服务器结果与后续本地 gate 审计给出了新的阻塞结论。S epoch 2、
reference level 12 的六步校准经 `.done`、`run.json`、最终 CSV 动作、`.time`、
前后 `meminfo` 和 14/14 payload 哈希独立核验后，确认为
`success/TrajectoryComplete`。七个误差点为
`0.615676/0.496907/0.425012/0.319588/0.238550/0.149087/0.101268`；最后两个比值
为 `0.624974/0.679256`，几何平均 `0.651550`，不是平台，且严格的 `0.1` 目标尚未命中。
该进程 peak RSS 约 17.31 GiB、process swaps 为 0；但上传的 build identity 记录
`PATCH_THREADS=8`，与冻结的 4 线程协议不符，所以只保留作数值校准证据，不作为 4 线程
资源证据。新脚本对 calibration/custom 的非 4 线程值 fail closed，并生成相对路径
`SHA256SUMS`，避免服务器绝对路径导致回传后无法直接验证。

S 最终网格把约 99.16% 的 `sum(eta_H_T^2)` 集中在 `r<0.5`，远场保持粗网格；最细区域
主要落在 cut-off 过渡环 `0.25<r<0.5`，而重入角点相邻六个单元仍为 `h=0.125`。
因此当前 calibration 支持局部化，但还不能完成 G7 的角点恢复主张。正式 S 结果仍须检查
角点、两条 Dirichlet 重入边和 Robin 边界附近的网格分布。

epoch-2 的独立 12→13 audits 均失败：R2a 的相邻 reference 差为 `0.0434044`、终点分数
为 `8.51361`；S 分别为 `0.116608` 和 `1.15148`，都远大于 0.25。R2a 随后在本机继续
完成 epoch 3/4/5、reference level 13/14/15 的同一六步轨迹，终点误差分别为
`0.00627629/0.00724420/0.00587381`；相应 13→14、14→15、15→16 audit 分数仍为
`4.61608/2.93656/2.44423`。level-15 运行 peak RSS 约 8.40 GiB、wall time 约 24:49、
swap 为 0；level 16 按资源增长已接近 23 GiB WSL 的安全边界，转服务器执行。

当前下一 gate 只包含 R2a epoch 6/reference level 16 与 S epoch 3/reference level 13 的
六步校准，随后分别做 16→17 和 13→14 reference audit。两条校准在服务器串行运行，
`PATCH_THREADS=4`，输出仍是 calibration，不得进入论文图表。只有两项最终
`terminal_error_fraction<=0.25` 才允许冻结十个正式 E1 配置；当前禁止启动十配置生产 E1。
完整逐点数据、run ID、资源记录和判定见
`experiments/helmholtz_adaptive_paper/calibration/reference-epoch-v4-2026-08-12.md`。

2026-08-13 经用户明确授权，先运行了两条 S/PALOD、κ=16、两次 H 加密的探索性 E1，
仅用于验证制造解后处理和绘图合同。epoch 0/reference level 10 与 epoch 1/level 11 均为
`success/TrajectoryComplete`，4 patch 线程、零 swap；峰值 RSS 分别约 1.29/1.74 GiB。
runner 现同时输出真正的非约束粗自由度 `DoF_H`、相对 reference error、制造解绝对误差和
制造解相对误差。误差求值为单向 post-processing sink，不参与 MARK/STOP。两条短轨迹的
reference gate 未通过/未执行，结果必须标为 exploratory，不能绕过上述正式 E1 gate。
这两条旧轨迹分别重建 hierarchy，epoch 1 从初始粗网格重新开始，因此只能作为制造解
后处理校准，不能解释为论文中的连续 reference epoch。论文的 refresh 语义是
`T_h_ref <- T_h_amb`，而当前 `T_H`、累计 H-step 与当前 `ell` 均继承。

2026-08-13 随后补齐单次 driver 运行内的显式 reference refresh：配置用累计
`reference_refresh_H_steps` 指定 epoch 边界，刷新前后 fail-closed 核对粗节点数、
粗自由度、粗单元数和 coarse mesh version，evaluation sink 按每个 epoch 的实际
reference mesh 分别构造/缓存 FEM reference。S/PALOD、κ=16 的四步受控 E1 试跑在
第 2 个 H-step 后从 epoch 0 刷新到 epoch 1；边界前后均为 `N_H=288`、`DoF_H=271`，
新 reference 的 `N_ref=4209` 等于旧 ambient 的 `N_amb=4209`，epoch 1 首次求解仍在
同一 `DoF_H=271` 上。轨迹最终为 `success/TrajectoryComplete`，末端
`DoF_H=425`，reference energy error 约 `0.212844`，制造解相对 energy error 约
`0.292275`。该短轨迹仍未通过正式 reference-adequacy gate，标为受控 E1 pilot，
不作为论文正式数据。
完整证据见
`experiments/helmholtz_adaptive_paper/calibration/exploratory-e1-exact-error-2026-08-13.md`。

corrector patch cache 的 correctness scaffold 已加入：相同离散状态全命中且与 full
rebuild 在 corrector、基、粗算子和 LOD 解上逐项一致；PDE 或实际 patch system 改变均
fail closed；若增大 ell 后 patch 已被边界截断为同一个系统，则允许数学等价复用。
当前一次局部 H 加密会改变全局准插值约束，小规模测试没有可严格复用的远场 patch，
因此缓存保持 opt-in、未接入 production driver，避免无命中时增加哈希和
存储开销。局部 H 复用须先把约束依赖局部化，再重新通过同一 full-rebuild 对照。

## 9. 验收门槛

| Gate | 必须通过的检查 | 未通过时禁止 |
|---|---|---|
| G0 基线 | 现有 Helmholtz/FEM/LOD/NVB 测试通过 | 修改实验 runner |
| G1 层级 | reference 固定、ambient ratio、嵌入失败结构化返回 | 自适应轨迹 |
| G2 Riesz | reference \(\eta_H\) 恒等式和局部分配守恒 | H 标记 |
| G3 回缩证书 | retraction identity/kernel 约束；小矩阵一侧控制 | \(\ell\) 决策 |
| G4 状态机 | 超阈值后只 `ell++`；无 h-refine；全局 \(\ell\) | 主算例 |
| G5 R1 校准 | 参数冻结，smooth 网格无异常集中 | E1 |
| G6 R2a | 五方法至少完成共同中等误差目标 | E2/R2 图表 |
| G7 S | 混合边界、制造解、角点加密和 FEM 基线均正确 | S 图表 |
| G8 复现 | clean build 可由配置生成 CSV、JSON、VTU 和图表；`run.json` 含冻结论文哈希 | 写入论文 |

G3 还需特别验证 ambient 证书确实控制 reference-space quantity
\(\lVert Q_{\infty,h}^*-Q_{\ell,h}^*\rVert\)。若当前离散投影不能证明或数值验证
retraction 性质，应让 runner fail closed 并报告 `retraction_not_verified`；不得重新引入
corrector h-refinement 分支绕开问题。

## 10. 论文最终需要的图表

只生成能够支持核心结论的图表：

1. R2a、S：reference energy error 对 \(N_H\) 和累计时间，各一图；
2. R2a、S：`eta_H`、`Theta_loc`、`ell` 随自适应步变化，各一图；
3. R2a、S：PALOD 最终粗网格与 `eta_{H,T}` 分布；
4. \(\kappa=8,16,32\)：达到共同中等误差目标时的 DOF、时间、\(\ell\) 汇总表；
5. R1 小规模：direct corrector perturbation 与 ambient-to-reference certificate 对照表。

主文不需要展示所有中间网格、所有敏感性扫描或所有内部常数。额外诊断放补充材料。

论文中的核心结论应限制为：

- `Theta_loc` 超阈值时 global \(\ell\) 更新能恢复 corrector 可接受性；
- admissible 后，\(\eta_H\) 能在局部源和角点奇性处有效引导 \(H\)-加密；
- 相比准一致 standard LOD 和普通 FEM，PALOD 在给定 reference error 下改善粗自由度
  或总工作量；
- ambient shadow 始终满足统一网格比，未观察到 ideal-corrector 离散误差恶化。

若数据不支持某一项，就删去该结论，不增加新参数扫描来“调出”预期结果。

## 11. 建议执行顺序和完成定义

### 第一阶段：算法合同（必须先完成）

1. WP1：reference epoch hierarchy；
2. WP2：reference/ambient 两套 kernel Riesz；
3. WP3：retraction 与 `Theta_loc`；
4. WP4：practical driver 和状态机测试。

### 第二阶段：统一实验基础设施

5. WP5：paper runner、schema v2/v3、轨迹式容差抽取；
6. R1 小矩阵端到端 smoke 和参数冻结；
7. corrector/cache full-rebuild 等价性测试。

### 第三阶段：生产实验

8. R2a、\(\kappa=16\) 的五方法主实验；
9. S、\(\kappa=16\) 的混合边界五方法主实验；
10. 两个算例的 \(\kappa=8,32\) 最小波数补充；
11. 最终配置重复 3 次、生成汇总图表并更新论文数值章节。

以下全部满足才算完成：

- 默认 driver 中不存在 corrector/reference `h`-refinement 决策；
- \(\eta_H\) 在 reference kernel 中计算，误差目标始终是 \(u_h\)；
- ambient-to-reference certificate 通过 retraction 验证；
- reference epoch 和任何 refresh 在日志中可追踪；
- R2a 和正式 mixed-boundary S 的主实验完成；
- standard LOD 与普通 FEM 基线使用同一 reference solution；
- 一条轨迹复用多个误差目标，资源限制和失败原因完整；
- 所有论文数字均可从版本化配置和原始 CSV 重建；
- practical 与 fully certified 的表述严格区分。

## 12. 暂缓项目

以下内容不阻塞这篇论文的核心数值实验：

- 对每个理论常数做区间验证和 directed rounding；
- 每一步重新计算参考 inf--sup 常数；
- 局部 \(\ell_T\) 或 corrector fine-mesh 自适应；
- 高对比系数、三维、并行强扩展；
- R2b 和 \(\theta_H\in\{0.3,0.7\}\) 的完整敏感性；
- 所有方法在所有波数、所有容差上的笛卡尔积；
- reference-to-continuous error 的完全认证。

这些项目只有在核心图表已经稳定、且论文确实需要补充证据时再开启。
