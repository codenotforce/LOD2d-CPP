# Reference-epoch 自适应 Helmholtz LOD：数值实验实施计划

> 状态日期：2026-08-23
>
> 唯一理论口径：`../LOD_paper/helmholtz_lod_certified_amsart_revised.tex`
>
> 核心目标：实现并验证论文中的 reference-epoch 自适应机制，以最少的重复计算完成
> localized smooth、L-shaped low-regularity、epoch behavior 和 wave-number dependence
> 四类数值结果。
>
> 默认结果口径：生产运行采用 practical/conditional 模式。只有在所有离线常数和代数
> 误差均给出严格界时，才把对应结果称为 certified。

### 论文版本基线

```text
file: ../LOD_paper/helmholtz_lod_certified_amsart_revised.tex
sha256: 94b0c1469312ce006f3b76d08b30f920115d274f442e3912ca660ccf919bd3f9
paper repository base commit: 30f6c05abeafcb5b76e74cb9f23d07ee1e13b13a
frozen: 2026-08-23
```

每次论文修改后必须生成新 SHA-256。生产结果的 `run.json` 保存实际使用的论文哈希，
已有结果不得被新哈希覆盖。

## 1. 冻结的算法口径

新版论文使用三张嵌套网格

\[
  \mathcal T_H\preceq\mathcal T_h\preceq\mathcal T_c,
  \qquad
  V_H\subset V_h\subset V_c\subset V.
\]

三者职责必须在代码、配置和输出中严格分开：

| 对象 | 代码名称 | 一个 epoch 内是否变化 | 唯一职责 |
|---|---|---:|---|
| \(\mathcal T_H\) | `coarse_mesh` | 是 | LOD 粗空间；由 \(\eta_{H,T}\) 标记 |
| \(\mathcal T_h\) | `reference_mesh` | 否 | corrector、参考解 \(u_h\)、误差目标 \(\|u_h-U\|_\kappa\) |
| \(\mathcal T_c\) | `candidate_mesh` | 是 | 为下一 epoch 准备；由 candidate indicator 局部加密 |

epoch 开始时仅执行一次

\[
  \mathcal T_c:=\mathcal T_h.
\]

candidate 不是 ambient certification mesh，也不定义当前误差目标。旧计划中的
`ambient_mesh/rho_amb/reference_retraction` 口径已经失效，不能继续作为论文 runner
的正式输出。

### 1.1 一个 epoch 内的状态顺序

```text
EpochInit
  -> CorrectorCheck
       Theta_loc / delta_loc_hat 过大: global ell++ -> CorrectorCheck
  -> LodSolveAndReferenceEstimate
       compute U, eta_H,T, eta_H, practical/reference bound
  -> ProposeCoarseRefinement
  -> EnrichCandidate
       candidate indicator + conforming hierarchy closure
  -> LazyDualDecision
       未触发: CommitCoarse -> CorrectorCheck
       decrease/interval/tolerance 触发: CandidateDualCheck -> RefreshOrContinue
       structural/level-gap 触发: 直接 ReferenceRefresh（不做冗余 dual solve）
  -> Terminate / ReferenceRefresh / CommitCoarse
```

必须满足：

- corrector certificate 只触发全局 \(\ell\leftarrow\ell+1\)，不触发 reference 或
  candidate 加密；
- \(\eta_H\) 只标记当前 coarse mesh；
- candidate equilibrated indicator 只决定 candidate enrichment；
- candidate dual-Riesz 只用于 epoch-switch 检查，并且懒惰调用；
- 算法永远不为 marking 或 switching 求解 candidate Helmholtz Galerkin 问题；
- reference refresh 后误差目标改变，输出必须增加 epoch 编号，不能把两段曲线当作同一
  reference error 轨迹。

## 2. 资源节约原则

### 2.1 不在每一步重复计算的量

下列量采用一次性离线输入或手动设置：

\[
 C_{\mathrm{app}},\ C_{\mathrm{ol}}(1),\ C_{\mathrm{sd}},\
 C_{\mathrm{ov}},\ C_a,\ c_W,\ C_I,\ C_F.
\]

生产循环不计算 reference inf--sup 常数，也不逐网格重新估计这些理论常数。推荐输入：

```text
c_res_usr       coarse resolution threshold for kappa H_T
theta_loc_usr   threshold applied directly to Theta_loc
C_rel_usr       aggregate multiplier for practical reference indicator
theta_H         coarse Dörfler parameter
theta_c         candidate marking parameter
q_dual          relative-decrease trigger for a lazy dual check
m_dual          maximum coarse steps between dual checks
tau_ep          epoch-switch ratio
```

默认初值仅用于启动校准：

```text
c_res_usr = 0.5
theta_H   = 0.5
theta_c   = 0.5
q_dual    = 0.5
m_dual    = 3
tau_ep    = 0.5
ell0      = 2
ell_max   = 6
```

`theta_loc_usr` 和 `C_rel_usr` 在一个小规模 smooth 校准中冻结。若这些数值不是严格
理论界，输出列使用 `U_practical`，不得写成 rigorous `U_ref`。

### 2.2 必须避免的重复工作

1. 每个 `(case, kappa, epoch)` 只解一次 reference problem；
2. 不解 candidate Helmholtz problem；只在很小的验证网格上允许离线求 \(u_c\) 检查
   dual-gap 实现；
3. 一条轨迹运行到最小误差目标，再事后抽取多个容差首次命中点；
4. 开发阶段不做重复计时；最终进入论文的配置只重复 3 次并取中位数；
5. corrector、reference Riesz 和 candidate dual Riesz 使用 patch factorization cache；
6. candidate dual check 只在 decrease/interval/tolerance lazy trigger 时执行；结构嵌入失败或
   `reference_refresh_level_gap` 已经确定必须刷新时直接切换 epoch，避免计算不影响决策的
   candidate dual；
7. candidate flux 默认只在 coarse marked set 的固定层 active region 上重构；正式声明
   这是 marking heuristic。只在小规模 audit 中运行全局 reconstruction 检查可靠性；
8. `Theta_loc` 使用 matrix-free Lanczos/LOBPCG、上一轮向量 warm start 和阈值分离停止；
   稠密 Gram/eigensolve 仅用于小矩阵单元测试；
9. 默认每步只写标量日志，网格场数据只保存 epoch 初始、switch 前和最终快照；唯一例外是 E1 的第一轮 epoch：为绘制真实网格演化图，保存选定 checkpoint 的 coarse/candidate 网格，并单独保存一次固定 reference 网格；
10. 不做所有参数、波数、方法和容差的笛卡尔积扫描。

## 3. 当前 LOD-C++ 与论文的差距

当前分支已经具备 NVB hierarchy、固定 reference residual Riesz、LOD corrector/PG solve、
粗网格 Dörfler 标记、reference error evaluation、practical driver、统一论文 runner 和
PALOD/HLOD/FEM 基线。以下是新版论文新增且尚未形成正式生产链的部分。

### WP1：将 legacy ambient 语义迁移为 candidate epoch hierarchy

涉及：

- `include/helmholtz/adaptive/hierarchy.h`
- `src/helmholtz/adaptive/hierarchy.cpp`
- `tests/test_helmholtz_adaptive.cpp`

新增或重命名接口：

```cpp
begin_reference_epoch();                  // candidate <- reference
propose_coarse_refinement(marked_H);
enrich_candidate(marked_c);
close_candidate_over_proposed_coarse();   // T_H^+ <= T_c^+
candidate_contains_proposed_coarse();
commit_coarse_refinement();
refresh_reference_from_candidate();
```

验收：reference 在 epoch 内的节点、单元和 generation 完全不变；candidate 是 reference
的 conforming refinement；候选 coarse mesh 只有在 reference 可容纳时才在本 epoch
提交，否则必须先 refresh。

### WP2：reference corrector certificate

涉及：

- `include/helmholtz/adaptive/kernel_residual.h`
- `src/helmholtz/adaptive/kernel_residual.cpp`
- `include/helmholtz/adaptive/certificates.h`
- `src/helmholtz/adaptive/certificates.cpp`

新版论文的 defect 完全位于 reference kernel：

\[
D_{\mathrm{loc}}^h(v_H)(w)
=a\bigl(w,(I-Q_{\ell,h}^*)v_H\bigr),
\qquad w\in W_h.
\]

实现局部 Riesz、\(G_{\mathrm{loc}}\) 的 matrix-free action 和

\[
\Theta_{\mathrm{loc}}^2
=\lambda_{\max}(G_{\mathrm{loc}},M_\kappa).
\]

现有 ambient-to-reference retraction certificate 仅保留 legacy 回归，不进入新版 runner。
小网格验收应直接组装 \(Q_{\infty,h}^*-Q_{\ell,h}^*\)，检查
\(\Theta_{\mathrm{loc}}\) 的上下界和随 \(\ell\) 的下降。

### WP3：compatibility-corrected candidate flux indicator

建议新增：

- `include/helmholtz/adaptive/candidate_flux.h`
- `src/helmholtz/adaptive/candidate_flux.cpp`
- `tests/test_helmholtz_candidate_flux.cpp`

正式实现采用二维 \(\operatorname{RT}_2/P_2\) patch mixed reconstruction，约定
\(\operatorname{RT}_2(K)=[P_2(K)]^2+xP_2(K)\) 且
\(\nabla\cdot\operatorname{RT}_2(K)=P_2(K)\)。这是当前 \(P_1\) primal field、
candidate hat function 和 \(f_c\in P_1(\mathcal T_c)\) 所产生的二次
divergence/normal-trace 数据所需的最低匹配 RT 阶数。必须包括：

- \(f_c\) 的固定阶多项式投影；
- 非 Dirichlet patch 的 compatibility defect \(\chi_a^c\)；
- patchwise correction
  \(\bar c_a=\chi_a^c/|P_a^c|\)；
- \(\delta_{\mathrm{PG}}=\sum_a\bar c_a\mathbf 1_{P_a^c}\)；
- Neumann/Robin normal-trace 约束；
- \(\eta_{\mathrm{eq},K}^c\) 和 candidate Dörfler marking；
- optional active-region reconstruction。

可以另设 RT0/P0 冒烟测试来快速检查 patch 装配，但它必须先投影二次 source 和
boundary data，并把投影缺陷计入 oscillation；RT0 不得用于论文主实验、正式
candidate marking 或 estimator effectivity 结论。若改用 BDM，则对应选择
\(\operatorname{BDM}_3/P_2\)，而不是 BDM2。

验收：patch compatibility 达到线性求解容差；全局 divergence 和 boundary flux 恒等式
成立；在小网格上，全局 indicator 控制直接计算的 residual dual norm。

### WP1--WP5 实施状态（2026-08-23）

当前工作树已经执行到 WP5 的组件级验收；新版论文 runner/schema 仍属于 WP7，
因此这里不把旧 PALOD 轨迹重新解释为 reference-epoch 生产结果，也不把 practical
常数字段称为 rigorous certificate。

- **WP1 / G1 已通过**：`ReferenceEpochHierarchy` 增加
  `begin_reference_epoch/propose_coarse_refinement/enrich_candidate/
  close_candidate_over_proposed_coarse/candidate_contains_proposed_coarse/
  commit_coarse_refinement/refresh_reference_from_candidate`。coarse proposal 在提交前不改变
  当前 coarse mesh；candidate enrichment 不改变固定 reference；proposal 越过 reference 时
  返回结构化 `ReferenceRefreshRequired`，refresh 后方可提交。legacy ambient accessor 仅为旧
  runner 回归保留。
- **WP2 / G2 已通过 implementation-study 验收**：reference defect 直接由
  `A_h^* (I-Q_{ell,h}^*)Phi_H` 形成，并在 `W_h` patch 上求 Riesz；新增
  matrix-free `G_loc` action，小维数才形成 dense Gram 交叉检查。localized smooth、
  `kappa=16`、`ell=1,2,3` 得
  `Theta_loc = 2.50414, 1.41851, 4.61105e-08`，direct defect 为
  `4.55780, 3.09622, 0`，均落在输入离线常数给出的 theorem bracket 内并随 `ell` 下降。
- **WP3 / G3 已通过 implementation-study 验收**：正式路径采用 contravariant-Piola
  二维 RT2 vertex-patch reconstruction 与 P1 source projection，二次 divergence/normal
  trace 逐点约束；包括 non-Dirichlet compatibility correction、Dirichlet free trace、
  Robin trace、`delta_PG`、element indicator、candidate Dörfler marking 和 optional active
  region。RT0/P0 API 仅保留为显式 smoke 路径。`kappa=16` 小网格中，R1 的
  compatibility/divergence/normal-jump/boundary 误差分别为
  `0/1.14e-13/2.76e-14/1.17e-14`，S mixed-boundary case 为
  `0/2.84e-14/6.22e-15/1.00e-15`；两例均满足
  `eta_eq >= discrete residual dual norm`。
- **WP4 / G4 已通过小网格验收**：`KernelRieszSpace::CandidateResidual` 在 (W_c) 上
  复用局部 saddle-point Riesz 框架，生产入口只接受已有的 (U|_{V_c})，没有求解
  candidate Helmholtz problem 的接口。结果分别输出 `L_gap_c` 或
  `L_gap_practical_c`。E0 直接求 (u_c) 的隔离审计得到
  `eta_dual_c=0.961709`、`L_c=L_gap_c=0.0150267`、
  `||u_c-U||=||u_c-u_h||=1.13069`，两个下界均成立。
- **WP5 / G5 已通过状态机验收**：新增独立 `ReferenceEpochPracticalDriver`，完整覆盖
  计划中的 12 个状态。转换测试验证了 corrector 失败只做 global `ell++`、首次 lazy
  dual 可跳过、终止前强制 dual、结构性 hierarchy trigger 跳过冗余 dual 并且必须先
  refresh 再 commit、
  numerical gap 优先于 termination，以及 `H/epoch/dual/reference/candidate/wall/ell`
  上限的结构化退出。backend contract 不暴露 candidate Galerkin solve。

可重建记录见
`experiments/helmholtz_adaptive_paper/calibration/reference-epoch-wp1-wp3-2026-08-23/README.md`。

### WP4：candidate dual-Riesz gap certificate

复用 `ReferenceResidualRiesz` 的局部 saddle-point 框架，使 kernel mesh 可选择
`reference` 或 `candidate`。在 \(W_c\) 上计算

\[
L_c=\frac{\eta_{\mathrm{dual}}^c}{C_aC_{\mathrm{ov}}},
\qquad
L_{\mathrm{gap}}^c=[L_c-U_{\mathrm{ref}}]_+.
\]

practical 模式可用冻结的 aggregate multiplier，但字段名必须表明是 practical gap。
小网格上允许求解 \(u_c\) 仅作验证，检查
\(L_c\le\|u_c-U\|_\kappa\) 和
\(L_{\mathrm{gap}}^c\le\|u_c-u_h\|_\kappa\)。生产 runner 禁止求 \(u_c\)。

### WP5：reference-epoch practical driver

建议新增独立状态机，避免继续扩大旧 PALOD driver：

```text
EpochInit
CorrectorCheck
SolveEstimate
ProposeCoarse
CandidateEnrich
LazyDualDecision
CandidateDualCheck
CommitCoarse
ReferenceRefresh
Converged
WorkLimitReached
Failed
```

每个状态有唯一允许的 mesh mutation，并写转换测试。`CorrectorCheck` 失败的下一状态
只能是 global `ell++` 后再次 `CorrectorCheck`。

### WP6：singularity-aware hybrid variant

L-shaped 主实验需要：

\[
\Omega_S=N^\ell(S),\quad
\Omega_F=N^{2\ell}(S),\quad
V_H|_{\Omega_F}=V_h|_{\Omega_F}.
\]

实现要求：

- 在 \(\Omega_F\) 内 coarse 与 reference 匹配；
- 若 corrector check 导致 \(\ell\) 增大，则同步重定义 \(\Omega_S=N^\ell(S)\)、\(\Omega_F=N^{2\ell}(S)\)，并在继续 corrector solve 前恢复 matching condition；
- 对 \(T\subset\Omega_S\) 跳过零 kernel corrector；
- coarse marking 仅在 regular region，global reference certificate 保持不变；
- candidate 仍可在 singular region 加密，为下一 epoch 准备；
- 记录跳过的 corrector 数量、维数和时间。

### WP7：runner、配置和输出 schema

将 `bench_helmholtz_adaptive_paper` 升级为 reference-epoch schema。旧配置仍可读，但
新版生产配置不得继续输出 `ambient_mesh/rho_amb/reference_retraction` 作为核心字段。

### WP6--WP7 与 E0 实施状态（2026-08-23）

- **WP6 / G6 已通过组件级验收**：实现
  `classify_singular_regions/restore_hybrid_reference_matching`，按共享顶点的
  coarse-element graph 构造 `Omega_S=N^ell(S)` 与 `Omega_F=N^(2ell)(S)`。
  每次 global `ell` 改变后，先在固定 reference 所允许的范围内事务式加密 coarse，恢复
  `V_H|Omega_F=V_h|Omega_F`，再组装 corrector；reference mesh/version 保持不变。
  `Omega_S` 中的零 kernel element corrector 被跳过并记录数量与工作维数；coarse
  Dörfler marking 仅使用 regular-region indicator mass，global reference residual
  certificate 不截断；candidate 仍可在 singular region 加密。小网格测试在 `ell=2`
  得到 `|Omega_S|=54`、`|Omega_F|=138`，首轮跳过 24 个 corrector、312 个工作单元。
- **WP7 / G7 已通过 schema-v5 数值 smoke**：现有
  `bench_helmholtz_adaptive_paper` 按 `schema_version` 分派；legacy v4 路径保持可读，
  v5 使用独立 `ReferenceEpochPaperConfig`，配置与核心输出均不含
  `ambient_mesh/rho_amb/reference_retraction`。数值 backend 连接 WP1--WP6，candidate
  阶段只延拓当前 reference 上的 LOD 值并执行 RT2/P2 flux 与 candidate dual-Riesz，
  不暴露 candidate Galerkin solve。v5 输出 `iterations.csv/run.json/summary.csv/
  epoch_history.csv/mesh_manifest.csv/corrector_work.csv`；lazy 数值写 `NA`，并记录
  coarse/reference/candidate DoF、corrector rebuild/skip、各阶段时间和最终三类网格。
- **E0 / G0 已通过五项校准**：持久产物位于
  `experiments/helmholtz_adaptive_paper/calibration/reference_epoch_e0_v5/`。
  R1、`kappa=16`、`H-level=3/reference-level=5` 上完成 `ell=1,2,3,4`
  direct defect/`Theta_loc` 对照、reference Riesz 与 nodal-to-element mass conservation、
  candidate RT2/P2 恒等式、隔离的 direct `u_c` dual-gap 下界审计，并一次性冻结
  `theta_loc_usr=1.961323606684523`、`C_rel_usr=1.6358263950741063`。两者均标记为
  implementation-study empirical parameters，不宣称为严格定理常数。

## 4. 最小可发表实验矩阵

实验按 gate 分四层；上一层未通过时不启动下一层。

### E0：代数与小网格校准

使用 localized smooth case，\(\kappa=16\)，一个很小的 fixed hierarchy：

1. \(\ell=1,2,3,4\) 的 direct corrector defect 与 `Theta_loc`；
2. reference residual Riesz 恒等式和 nodal-to-element mass conservation；
3. candidate RT2/P2 compatibility/divergence/boundary-flux identities；
4. 小网格直接求 \(u_c\)，验证 dual gap 两个下界；
5. 冻结 `theta_loc_usr/C_rel_usr`，主实验不再逐案例调参。

这些是验证，不做性能结论。

### E1：localized smooth 主实验，\(\kappa=16\)

比较论文规定的四种方法：

1. reference-epoch adaptive LOD；
2. fixed-oversampling LOD；
3. uniform conforming \(P_1\) FEM；
4. adaptive conforming \(P_1\) FEM。

主图包括 exact/reference energy error 对 DOF 和累计时间，以及第一轮 reference epoch 的
真实网格演化。网格图只针对 proposed reference-epoch adaptive LOD：

- 单独绘制一次固定 reference mesh \(\mathcal T_h^{0}\)；
- 在同一组 checkpoint 上绘制 coarse mesh \(\mathcal T_H^{0,n_j}\) 的变化；
- 在相同 checkpoint 上绘制 candidate mesh \(\mathcal T_c^{0,n_j}\) 的变化；
- 最后一个 checkpoint 取第一轮 epoch 的 pre-switch 状态；若第一轮未发生 switch，则取该
  epoch 内最后一个可用状态。

为使网格演化可读，E1 的初始 coarse/reference level gap 应预先设置得足以容纳多个 coarse
commit；但不得为了作图而抑制算法本应触发的 epoch switch。网格图与 \(\ell\)、epoch、
candidate dual-check 时间线使用同一 iteration/epoch 编号。

#### E1 执行状态（2026-08-23）

E1 四条本地 implementation-study 轨迹及绘图链已经完成，但旧 PALOD 的 `H2/h10/12-step`
轨迹只保留为诊断结果，不再作为最终主曲线：epoch 0 后段的 exact error 约 0.14 平台来自
level-10 reference 误差地板，而且刷新发生在总预算末端。新的 PALOD 生产配置改为
`H2/h12`、`ell0=2`、最多 15 个 H-step/4 个 epoch，并在 prospective coarse 与 reference
的最小局部 NVB level gap 达到 4 时（每个 epoch 至少 3 次 H commit 后）提前刷新。后续语义
审查又要求 refresh 前把 candidate 加深到 gap=6、`ell` 跨 epoch 继承，并禁止开启少于 4 个
解点的新 epoch；详见本节后续的 corrected production 记录。reference error 在 epoch refresh
处继续分段，不能跨目标比较；刷新前后的 manufactured exact error 才是同一口径。

fixed-oversampling LOD 的正式配置采用 `ell=3`、reference level 15。h12 冻结对照中，
`ell=3/4` 在 H-step 11 的相对精确能量误差分别为 0.0623348/0.0620395，差异约 0.48%；
同时 wall time 为 14.18/30.70 s、峰值 RSS 为 354/972 MB。因此 `ell=3` 用于 E1 主曲线，
`ell=4` 仅保留为充分性 audit；在补做同口径对照前不把 `ell=2` 用作正式结果。

本地单次运行的资源结果为：PALOD 2:36.63、5.76 GiB；fixed LOD (`ell=3`)
3:47.04、1.44 GiB；UFEM 2:59.64、1.85 GiB；AFEM 0:19.71、0.304 GiB。这些只用于
实现验证和瓶颈定位，论文最终计时仍须按冻结硬件做三次重复并报告中位数。

已落地的性能实现包括：reference-defect Gram 局部分解缓存、patch OpenMP、block action、
四向量局部最优 block Rayleigh--Ritz（保留历史方向并复用上一轮 action）和跨 H-step
特征向量延拓；fixed-LOD 的同构 patch 分组、分解复用和
有界 cache；reference residual Riesz 的同构约束系统分组、多 RHS、SPD 分解加小型 Schur
complement。针对旧轨迹累计 75.1 s 的 candidate flux，现已增加 RT2 vertex-patch OpenMP、
按 vertex 的确定性归并，并按 element 缓存 quadrature、Piola RT2 basis、15x15 local mass
与 P2 点；candidate dual patch Riesz 也并行执行。结构/level-gap 强制刷新不再计算冗余 dual。
输出细分记录 candidate close/operator/prolongation/flux/enrich，以及 dual
operator/load/prolongation/prepare/patch-solve/reduction。跨 H-step 的 candidate patch
fingerprint/LRU 分解复用仍是后续优化，不在本次结果中提前宣称完成。

`H2/h12/gap4` 的 8-step 受控 probe 已在 12 GiB WSL 上通过：2:36.06、峰值
5,128,604 KiB、无 swap，并以计划内 `maximum_H_steps reached` 结束。旧实现曾在
`N_H=161` 的 localization 迭代 300 次后失败；新 block 路径在 `N_H=161/242` 分别用
5/3 次迭代达到 `5.18e-5/8.37e-5` 残差。gap=4 的强制 refresh 行明确记录
`SkipCandidateDual`，而 gap=7/5 的普通 interval checks 仍执行 dual。细分计时同时表明：
在 epoch 1 的 `N_H=161,N_c=20722` 步，candidate 总时间 22.68 s 中 RT2 reconstruction
为 4.47 s、candidate mesh enrichment/embedding 更新为 17.94 s；因此完成 RT2/dual
并行后，下一性能重点应转为 candidate hierarchy refinement 与 embedding 的增量更新。

candidate hierarchy 已改为增量更新：每次 NVB 加细产生的 step prolongation 直接与既有
reference-to-candidate/coarse-to-candidate 的 nodal、element 和 DG prolongation 组合，parent
map 沿 step parent 传播，固定 reference 的 quasi-interpolation 不再重建；candidate 侧只重建
自身 quasi-interpolation，并执行轻量维数/边界/parent-range 检查。连续两次局部加细已与全量
几何 embedding 重建逐项比对到 `1e-12`。在 `N_H=161,N_c=20722` 的同一轨迹点，candidate
enrichment 从 17.94 s 降为 0.227 s（约 79 倍），8-step 总 wall time 从 2:36.06 降为
2:06.03（约 19.2%），数值/状态列完全一致。新增 CSV 分阶段字段为
`time_candidate_nvb_refine`、`time_candidate_embedding_composition`、
`time_candidate_parent_map_update`、`time_candidate_quasi_interpolation` 和
`time_candidate_embedding_validation`。

正式 15-step E1 PALOD 已在同一 12 GiB WSL 上完成，run id 为
`R1_PALOD-reference-epoch_k16_r0_09fa933db6d58f2f`。轨迹包含 15 个解点和 4 个 epoch，
三次 refresh 后 exact relative energy error 分别从 0.2747 降至 0.1262、从 0.05394 降至
0.03808、从 0.03047 降至 0.02406，没有旧 h10 轨迹的 0.14 reference-error 平台，也没有
refresh 反弹。10%/5% target 已在 `N_H=161/479` 达到；固定 15-step 预算内 2%/1% target
尚未达到。运行以计划内 `WorkLimitReached: maximum_H_steps reached` 结束，wall time
9:06.04，峰值 RSS 9,000,792 KiB，无 swap。为使本地运行可完成，prepared saddle patch
完成 SparseLU 后立即释放不再使用的 local energy/constraint 副本和 row-major defect RHS；
该释放不改变 Gram action，相关 kernel residual、state-machine、hierarchy 回归均通过。

上述 `09fa...` 轨迹随后被 epoch 语义审查降级为诊断结果：它在 refresh 后把 `ell` 重置为
`ell0=2`，且新 epoch 的局部 reference/coarse gap 只有 4、2，最后还产生了一个只有单点的
epoch 3，不能用于分段收敛率。修正后的生产规则为：`ell` 跨 epoch 继承；强制 refresh 阈值
仍为 4，但 refresh 前按 proposed coarse 定向加深 candidate，使提升后的 reference gap 至少
为 6；`minimum_H_steps_per_epoch=3`，并新增
`minimum_solved_points_per_new_epoch=4`，剩余预算不足时在 refresh 前结构化结束。收敛率输出
按 epoch 分段拟合，任何 exact/reference exponent 都不得跨 refresh。

proposed coarse 事务现在缓存 proposed-to-reference/candidate 的 nodal、element、DG embedding，
candidate NVB 后用 step prolongation 增量更新，containment、level-gap、closure、refresh 和 commit
均复用缓存。修正后的 8-step probe 中，epoch 1 从 gap=6、`ell=3` 开始；累计 `time_mesh`
仅 1.45 s，`LazyDualDecision` 约 0.00056 s，最后 commit 0.076 s，相比旧轨迹的
85.3/47.2/34.7 s 已消除全局几何搜索瓶颈。无跨步 Gram LRU 的 probe 为 2:32.36、峰值
5,572,280 KiB、无 swap。

reference-defect Gram 增加了轨迹私有、有界数值指纹 LRU，并在 reference refresh 时清空。
512-entry pilot 只有 33 hits/682 misses，却把峰值推到 11,593,028 KiB，因此不得采用；32-entry
版本为 21 hits/694 misses，Gram factorization 从 60.70 s 降到 57.88 s，wall time 从
2:32.36 降到 2:28.43，但峰值增至 6,394,072 KiB。生产默认仅保留 32 项的小型缓存，输出
`gram_factor_cache_hits/misses`，不宣称大幅加速。RT2 跨 H-step factor cache 尚未实现；在
mesh/commit 瓶颈消除后，它仍是下一项独立优化，candidate dual 暂不优先。

修正后的 15-step 主实验因继承 `ell=3`，在本地 4:49 时 RSS 已达约 11.11 GB、available
仅约 0.38 GB，已在 swap/OOM 前主动中止；runner 仅在结束时写正式结果，因此该中止不形成
论文数据。完整修正轨迹必须在 366 GB 服务器上运行。

### E2：L-shaped low-regularity 主实验，\(\kappa=16\)

比较：

1. singularity-aware reference-epoch LOD；
2. standard reference-epoch LOD；
3. adaptive \(P_1\) FEM；
4. fixed-oversampling LOD 仅在成本可承受时补做。

主图：corner 附近最终网格、error-vs-DOF、error-vs-time、matching region 中省掉的
corrector work。必须使用论文的 mixed \(\Gamma_D/\Gamma_R\) 边界；全 Robin surrogate
不能写成正式 L-shaped 结果。

#### E2 启动状态（2026-08-23）

- schema-v5 已显式携带并哈希 case-S 制造解参数；E2 固定
  \(\gamma=0\)、quintic cutoff 外半径 1、smooth-wave amplitude \(B=0.05\)，即非振荡
  corner singularity 为主、保留小振幅 \(e^{i\kappa x}\) 光滑波。凹角边继续使用
  Dirichlet、外边 Robin，禁止退回历史默认 \(\gamma=1,B=0\)。
- hybrid 的 level-gap guard 改为 regular region 中的最小**正** gap。\(\Omega_F\) 内
  \(H=h\) 的零 gap 是 matching condition，不是 reference exhaustion；若 closure 真正要求
  细化零-gap 单元，严格 proposed-to-reference containment 仍会立即触发 structural refresh。
- 本地 `H3/h8/4-step` hybrid probe（8 threads）wall 约 45--48 s、峰值 RSS 约
  3.6 GiB、无 swap；exact relative energy 依次为
  `0.10921, 0.08445, 0.06691, 0.05352`。每次 corrector check 跳过 54 个 patch、约
  2134--2200 个工作单元。主耗时为 `Theta_loc/Gram` 约 15 s 和 candidate RT2 flux
  约 13 s，当前规模尚无单一异常瓶颈，故不在主实验前引入未经验证的新 cache。
- 同口径 standard probe 表明 H3 初始粗空间处于强预收敛区；证书自动把 ell 提至 7，
  四点 exact error 仍由 `4.46` 降至 `2.01`。这条轨迹作为消融对照保留，不与 hybrid
  共用固定 ell，也不用于声称低正则性最优阶。
- 服务器正式配置为 hybrid/standard `H3/h12/15-step`、gap guard 4、refresh target 6、
  新 epoch 至少保留 4 个解点；hybrid `ell<=4`，standard `ell<=10`。另配 AFEM
  `H3->level20` 和可选 fixed-LOD `H3/h16,ell=3`。完整命令和验收见
  `experiments/helmholtz_adaptive_paper/E1_E2_REFERENCE_EPOCH_SERVER_2026-08-23.md`。
- 本地 probe 只证明下降趋势和资源可行性，不能代替 E2 主实验。是否观察到稳定
  \(N^{-1/2}\) 必须由服务器深轨迹按 epoch 分段拟合后判断；不得跨 refresh 或删除
  pre-asymptotic 点来获得目标斜率。

### E3：epoch behavior

直接复用 E1 中一条具有至少两个 epoch 的轨迹，不额外重跑。绘制：

\[
\|u_h-U\|_\kappa,\quad U_{\mathrm{practical/ref}},\quad
\eta_H,\quad\eta_{\mathrm{eq}}^c,
\]

并在 lazy dual check 时补
\(\eta_{\mathrm{dual}}^c,L_{\mathrm{gap}}^c\)。图中用竖线标出 reference refresh；
每个 epoch 的 reference validation error 分段显示。

### E4：最小 wave-number study

localized smooth case 取

\[
\kappa\in\{8,16,32\}.
\]

只比较 reference-epoch adaptive LOD、fixed-oversampling LOD 和 AFEM。\(\kappa=16\)
直接复用 E1；因此只新增 6 条轨迹。UFEM 不重复整个 wave-number sweep。

默认生产量约为：E0 五个小检查、E1 四条、E2 三至四条、E4 六条，总计 18--19 条
轨迹。一个轨迹同时提供全部误差目标，不按 tolerance 重跑。

## 5. 误差目标和运行上限

所有方法事后提取首次达到

\[
10^{-1},\qquad5\times10^{-2},\qquad2\times10^{-2},\qquad10^{-2}
\]

的迭代。未达到者记 `not_reached`。每条配置必须含：

```text
ell_max
maximum_H_steps
maximum_epochs
maximum_reference_unknowns
maximum_candidate_unknowns
maximum_wall_seconds
maximum_dual_checks
```

达到任何上限都结构化退出，不自动放大网格或内存预算。

## 6. 最小输出合同

`iterations.csv`：

```text
schema_version, manuscript_sha256, code_commit,
case, method, kappa, run_id, repeat_index,
epoch, iteration, state, action, stop_reason,
N_H, N_h, N_c, ell, kappa_H_max,
eta_H, Theta_loc, U_practical,
eta_eq_c, eta_dual_c, L_gap_c, dual_check_performed,
relative_reference_energy, relative_exact_energy, relative_exact_L2,
marked_H, marked_c, rebuilt_correctors, skipped_correctors,
time_corrector, time_theta, time_lod_solve, time_reference_riesz,
time_candidate_flux, time_candidate_close, time_candidate_operator_assembly,
time_candidate_prolongation, time_candidate_flux_reconstruction,
time_candidate_flux_prepare, time_candidate_flux_patch_solve,
time_candidate_flux_merge, time_candidate_flux_audit,
candidate_flux_parallel_threads, time_candidate_enrich,
time_candidate_dual, time_candidate_dual_operator_assembly,
time_candidate_dual_load_assembly, time_candidate_dual_prolongation,
time_candidate_dual_solve, time_candidate_dual_prepare,
time_candidate_dual_patch_solve, time_candidate_dual_reduction,
candidate_dual_patch_factorizations, candidate_dual_parallel_threads,
time_mesh, time_total_cumulative,
peak_memory_mb
```

未计算的 lazy 量必须写 `NA`，不能沿用上一次数值冒充本步结果。

其他文件：

- `run.json`：完整配置、代码 commit、论文 SHA-256、编译器、硬件、退出状态；
- `summary.csv`：每个误差目标的首次命中点；
- `epoch_history.csv`：每次 refresh 前后 mesh/target/dual-gap；
- `mesh_manifest.csv`：记录 `case, epoch, iteration, stage, mesh_role, filename, N_cells, N_dofs`，
  作为网格文件与 `iterations.csv` 的唯一对应关系；
- 一般运行只保存 `epoch_start/pre_switch/final` 网格快照；
- E1 第一轮 epoch 额外保存
  `mesh_E1_e000_reference.vtu`，以及每个选定 checkpoint 的
  `mesh_E1_e000_iXXX_coarse.vtu`、`mesh_E1_e000_iXXX_candidate.vtu`；
- `corrector_work.csv`：patch dimension、rebuild/reuse/skip 和累计时间。

E1 的 reference mesh 在第一轮 epoch 内只写一次，因为它严格固定。coarse/candidate checkpoint
必须共享同一 `iteration`，并优先取 `epoch_start`、若干已 commit 的中间状态和
`pre_switch`。不要保存未提交的 prospective coarse mesh 作为正式 coarse 演化状态；若希望
诊断 \(\mathcal T_H^+\)，使用单独的 debug 文件名且不进入论文图。

## 7. 正确性与资源 Gate

| Gate | 验收内容 | 未通过时禁止 |
|---|---|---|
| G0 | 现有 NVB、FEM、LOD、reference Riesz 测试通过 | 新 driver |
| G1 | coarse/reference/candidate 全程嵌套；reference epoch 内不变 | 自适应运行 |
| G2 | corrector direct defect 与 `Theta_loc` 界、\(\ell\) 下降 | oversampling 决策 |
| G3 | RT2/P2 compatibility、divergence、boundary flux、global residual audit | candidate marking |
| G4 | candidate dual gap 小网格下界 | epoch switch |
| G5 | 状态机转换和所有 work-limit 结构化退出 | E1 |
| G6 | smooth manufactured boundary、forcing、exact error 验证 | smooth 主图 |
| G7 | L-shaped mixed boundary、singular forcing、matching region | L-shaped 主图 |
| G8 | cache 与 full rebuild 在小网格上逐步等价 | 开启 cache |
| G9 | clean build 可由配置重建 CSV、VTU 和图表 | 写论文结果 |

## 8. 优化实施顺序

1. full rebuild 版本通过 G1--G7；
2. reference solution 和 reference operators 按 epoch 缓存；
3. corrector patch factorization 按 mesh/patch/\(\ell\)/\(\kappa\) fingerprint 缓存；
4. reference Riesz 与 corrector certificate 共用 \(b_\kappa\)-Gram factorization；
5. candidate RT2/P2 patch matrix 复用，只更新 RHS；
6. active candidate region；
7. lazy candidate dual check；
8. matrix-free `Theta_loc` warm start；
9. 最后才做并行 patch solve。

任何 cache 必须先与 full rebuild 做逐迭代等价测试。不要为了减少运行时间降低为单精度。

## 9. 论文最终需要的图表

只保留下列能够直接支持论文结论的结果：

1. localized smooth：四方法 error-vs-DOF 与 error-vs-time；
2. L-shaped：三至四方法 error-vs-DOF 与 error-vs-time；
3. E1 第一轮 epoch 的真实网格演化图：固定 reference mesh，以及相同 checkpoint 下的
   coarse mesh 和 candidate mesh 演化；
4. L-shaped 最终 hybrid/matching-region 网格，用于展示 singularity-aware 结构；
5. 至少一个 epoch 时间线：reference error、candidate enrichment、lazy dual check、refresh；
6. \(\kappa=8,16,32\) 的共同目标工作量表；
7. 小网格 `Theta_loc` 与 direct corrector defect 表；
8. singularity-aware 方法省掉的 corrector solve/work 表。

开发诊断、逐步场变量、完整参数扫描和失败轨迹不进入正文。E1 第一轮 epoch 的网格
checkpoint 是为正文图预先声明的例外，不属于无目的的逐步场数据输出。

## 10. E1 第一轮 epoch 的网格演化图

这里不再生成抽象的 mesh schematic，而是直接读取 E1（localized smooth, \(\kappa=16\)）
第一轮 reference epoch 的实际 VTU 输出。绘图脚本统一使用

```text
tools/visualization/plot_reference_epoch_meshes.py
```

但脚本职责改为读取 `mesh_manifest.csv` 和对应 VTU，而不是程序化构造示意网格。

### 10.1 保存时刻与状态定义

设第一轮 epoch 为 \(e=0\)。reference mesh
\(\mathcal T_h^0\) 在整个 epoch 内固定，因此只保存和绘制一次。coarse 与 candidate
网格在同一组 checkpoint \(n_0<n_1<\cdots<n_J\) 上比较：

- \(n_0\)：epoch 初始化后的已提交状态；此时
  \(\mathcal T_c^{0,n_0}=\mathcal T_h^0\)；
- 中间 \(n_j\)：完成该次 iteration 后的已提交 coarse 状态和对应 candidate 状态；
- \(n_J\)：第一轮 epoch 的 `pre_switch` 状态，即 refresh 发生前的最后一个实际状态。

如果第一轮 epoch 的 committed coarse step 很多，正文只选择 3--4 个具有代表性的 checkpoint，
例如 `start / early / late / pre_switch`；所有原始 checkpoint 仍由 `mesh_manifest.csv` 记录。
如果第一轮 epoch 很短，则展示实际存在的全部状态，不人为延迟 switch。

### 10.2 图的布局

正文优先使用一个组合图：

```text
                         fixed reference mesh T_h^0
                   [one panel spanning the full width]

coarse mesh:       T_H^{0,n0}   T_H^{0,n1}   ...   T_H^{0,nJ}
candidate mesh:    T_c^{0,n0}   T_c^{0,n1}   ...   T_c^{0,nJ}
```

coarse 与 candidate 两行必须使用完全相同的 checkpoint、坐标范围和纵横比。reference mesh
不在每一列重复，因为其固定性本身正是 reference-epoch 方法需要展示的性质。图标题或 caption
明确写出

\[
  \mathcal T_H^{0,n_j}\preceq\mathcal T_h^0\preceq\mathcal T_c^{0,n_j}
\]

对所有展示的 \(n_j\) 成立。

localized smooth 精确解的局部中心 \((x_\star,y_\star)=(3/4,1/2)\) 可用一个小标记
叠加在各 panel 中，便于观察 coarse/candidate refinement 是否向高梯度区域集中；该标记仅用于
可视化，不参与自适应算法。

### 10.3 输出文件

脚本输出

```text
../LOD_paper/figures/E1_epoch0_mesh_evolution.pdf
../LOD_paper/figures/E1_epoch0_mesh_evolution.png
../LOD_paper/figures/E1_epoch0_mesh_evolution.svg
```

论文最终应以该真实网格演化图替换当前 `Mesh schematics` 中的示意图。若开发阶段仍保留
`epoch_mesh_hierarchy.pdf`，只能作为内部诊断图，不作为最终 numerical evidence。

### 10.4 绘图验收

绘图前自动检查：

1. 所有 checkpoint 的 epoch 均为 0；
2. reference mesh hash 在这些 checkpoint 上完全一致；
3. 每个 checkpoint 满足
   \(\mathcal T_H^{0,n_j}\preceq\mathcal T_h^0\preceq\mathcal T_c^{0,n_j}\)；
4. coarse/candidate 文件来自同一 iteration；
5. `pre_switch` 图像在 refresh 修改 reference mesh 之前输出；
6. 图中 mesh cell count 与 `iterations.csv/mesh_manifest.csv` 一致。

## 11. 完成定义

以下全部满足才算“完成论文数值实验”：

- 新 runner 不再以 ambient/retraction 口径解释 candidate mesh；
- 当前误差目标始终是固定 epoch 的 \(u_h\)；
- corrector 超阈值只导致 global `ell++`；
- candidate enrichment 不求 candidate Helmholtz 解；
- lazy dual check 和结构嵌入条件都能触发正确的 refresh；
- smooth、L-shaped、epoch、wave-number 四类结果完成；
- L-shaped 使用正式 mixed boundary；
- 一条轨迹复用多个误差目标；
- 每个数字可由配置、代码 commit 和 manuscript hash 重建；
- practical 与 rigorous certified 结果在表格中明确区分；
- 最终论文 PDF 无未定义引用、无 LaTeX 警告、图像可读。

## 12. 暂缓项目

以下项目不阻塞核心论文实验：

- 每步重新计算全部理论常数或 reference inf--sup；
- candidate Helmholtz solve；
- 所有 \(\theta_H,\theta_c,q_{\mathrm{dual}},m_{\mathrm{dual}}\) 敏感性；
- 所有方法在所有波数的笛卡尔积；
- 三维、高对比、多右端项和并行扩展；
- 局部 oversampling \(\ell_T\)；
- directed rounding 和所有代数量的严格区间验证；
- 超出正文所需的网格动画和每步高分辨率 VTU。
