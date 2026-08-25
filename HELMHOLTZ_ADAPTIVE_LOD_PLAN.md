# Reference-epoch 自适应 Helmholtz LOD：数值实验实施计划

> 状态日期：2026-08-24
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
sha256: 71f59581ea3a5e4cd65659055915715b9c86c582793db7154a5f0b3b31843ca8
paper repository working-tree base: 30f6c05
frozen: 2026-08-24
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

### 1.2 2026-08-24 修订算法的强制顺序

Algorithm 1 的实现顺序冻结为：先形成 prospective coarse；完成 candidate 的 Dörfler
标记和 refinement，再添加 hierarchy closure；随后先判 structural/reserve refresh，且该分支
不得求 candidate dual。每个 epoch 第一解初始化 lazy-dual baseline，终止前强制做一次
dual check；刷新时继承旧 epoch 的最终 `ell`。

Algorithm 2 额外冻结为：

\[
  \ell_S=\min\{j:B_{R_*}(S)\subset N^j(S)\},\qquad
  \Omega_S=N^{\ell_S}(S),\qquad
  \Omega_F=N^\ell(\Omega_S).
\]

在 `Omega_F` 上 coarse/reference 精确匹配；只省略 `Omega_S` 中的零 kernel corrector；
coarse estimator 在完整 regular region 上做 bulk marking；candidate indicator 在
`Omega_F/regular` 两侧分别做 Dörfler。prospective matching 必须在 candidate estimator
之前完成。

一般 conforming NVB 上，`Omega_F` 精确匹配与界面外所有 regular 单元统一正 reserve
缺少可实现性条件。当前实现明确标为 `implementation-erratum`：采用
`min(g_tar,max(0,d_graph(T,Omega_F)-b))`，自动选择最小 spill-free collar，并把每个
trial/closure 写入 `hybrid_reserve.csv`。这不是可以隐藏的性能细节；论文下一版应同步修正
reserve 条件。coarse marker 可从 closure-safe regular 子集选择，但最终 marked mass 仍须
达到原始完整 regular mass 的 Dörfler 下界；若达不到则回退到 full-regular 标记并允许
structural refresh。

两套生产实验按既定 implementation-study 口径保留强预收敛点，不执行论文
Eq. (coarse-admissibility) 的解析分辨率门；reference/candidate well-posedness 采用数值稀疏
分解检查。因此 Algorithm 1 的 `manuscript_conformance` 记为
`implementation-study-variant`，不得写成 `direct`。Algorithm 2 另记录
`reserve_trigger_scope=far-full-target-regular-cells`，因为 graded collar 的触发域不是论文
当前公式中的完整 prospective regular set。

最终本地门禁（2026-08-24）：WSL Release 构建后 `ctest -R helmholtz -j4` 为
`50/50` 通过；schema-v6 smoke 的 `iterations.csv/corrector_work.csv` 分别为 139/49 列且
逐行一致。服务器科学门已按新版定义取消旧的 `ell_S>=ell` 判据，并把
`RefreshReference` 的旧 epoch 行映射到真正的刷新后 `epoch+1`；同一结果目录若 commit、
二进制 SHA256 或 patch thread 数不同会拒绝复用，防止跨构建混合 `.done`。

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
7. candidate RT2/P2 flux 在完整 candidate mesh 上重构；`Omega_F/regular` 两侧只分拆同一
   全局 indicator 的 marking mass，不截断 estimator；
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
\ell_S=\min\{j:B_{R_*}(S)\subset N^j(S)\},\quad
\Omega_S=N^{\ell_S}(S),\quad
\Omega_F=N^\ell(\Omega_S)=N^{\ell_S+\ell}(S),\quad
V_H|_{\Omega_F}=V_h|_{\Omega_F}.
\]

其中物理半径 `R_*>0` 在整条轨迹开始前冻结。物理半径只决定 singular core；实际
oversampling `ell` 决定 transition/matching buffer。新版论文不再要求 `ell_S>=ell`，也不再
使用旧的 `Omega_F=N^(2 ell_S)(S)`。

实现要求：

- 在 \(\Omega_F\) 内 coarse 与 reference 匹配；
- 每次 coarse/reference/`ell` 改变后由冻结 `R_*` 重算最小 `ell_S` 和
  `Omega_F=N^ell(Omega_S)`，并在 corrector solve 前恢复 matching condition；
- 对 \(T\subset\Omega_S\) 跳过零 kernel corrector；
- coarse marking 仅在 regular region，global reference certificate 保持不变；
- candidate 在 `Omega_F` 与 regular region 分别做 bulk marking，仍可在 singular region
  加密，为下一 epoch 准备；
- 记录跳过的 corrector 数量、维数和时间。

### WP7：runner、配置和输出 schema

将 `bench_helmholtz_adaptive_paper` 升级为 reference-epoch schema。旧配置仍可读，但
新版生产配置不得继续输出 `ambient_mesh/rho_amb/reference_retraction` 作为核心字段。

### WP6--WP7 与 E0 实施状态（2026-08-24）

- **WP6 / G6 已通过组件级验收**：实现新版固定物理半径版本
  `classify_singular_regions_with_physical_radius/
  restore_hybrid_reference_matching_with_physical_radius`。按共享顶点 coarse-element
  graph 选择最小 `ell_S`，使固定圆 `B_R(S)` 完全包含于
  `Omega_S=N^ell_S(S)`，再取 `Omega_F=N^ell(Omega_S)`。每次 coarse 或 global `ell`
  改变后，先在固定 reference 所允许的范围内事务式加密 coarse，恢复
  `V_H|Omega_F=V_h|Omega_F`，再组装 corrector；reference mesh/version 保持不变。
  `Omega_S` 中的零 kernel element corrector 被跳过并记录数量与工作维数；coarse
  Dörfler marking 仅使用 regular-region indicator mass，global reference residual
  certificate 不截断；candidate 的 F/R 两侧分别达到同一 `theta_c`。真实 H3/h8 回归在
  `R_*=0.125,ell=2` 得到 `ell_S=3`、`|Omega_S|=54`、`|Omega_F|=150`。
- **WP7 / G7 已通过 schema-v6 数值 smoke**：现有
  `bench_helmholtz_adaptive_paper` 按 `schema_version` 分派；legacy v4 路径保持可读，
  v6 使用独立 `ReferenceEpochPaperConfig`，配置与核心输出均不含
  `ambient_mesh/rho_amb/reference_retraction`。数值 backend 连接 WP1--WP6，candidate
  阶段只延拓当前 reference 上的 LOD 值并执行 RT2/P2 flux 与 candidate dual-Riesz，
  不暴露 candidate Galerkin solve。v6 输出 `iterations.csv/run.json/summary.csv/
  epoch_history.csv/mesh_manifest.csv/corrector_work.csv/hybrid_reserve.csv`；lazy 数值写 `NA`，并记录
  coarse/reference/candidate DoF、corrector rebuild/skip、各阶段时间和最终三类网格。
- **E0 / G0 已通过五项校准**：持久产物位于
  `experiments/helmholtz_adaptive_paper/calibration/reference_epoch_e0_v6_revised/`。
  R1、`kappa=16`、`H-level=3/reference-level=5` 上完成 `ell=1,2,3,4`
  direct defect/`Theta_loc` 对照、reference Riesz 与 nodal-to-element mass conservation、
  candidate RT2/P2 恒等式、隔离的 direct `u_c` dual-gap 下界审计，并一次性冻结
  `theta_loc_usr=2.6311141730106273`、`C_rel_usr=1.4706194473212615`。两者均标记为
  implementation-study empirical parameters，不宣称为严格定理常数。最终生产二进制还须
  生成第七个 `06-calibration-provenance.json`，闭合 code commit、binary SHA、论文 SHA、
  compiler 和线程数来源链；缺少该文件时不得启动 main。

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

#### E1 新版论文重基线（2026-08-24）

旧 E1 制造解和旧 H2/h10--h12 PALOD 结果全部降级为历史诊断。新版 R1 固定为

\[
u=x^2(1-x)^2\sin(\pi y)\exp(-80((x-0.75)^2+(y-0.5)^2))e^{i\kappa x},
\]

top/bottom Dirichlet、left Neumann、right Robin。新版 E0 后冻结
`theta_loc_usr=2.6311141730106273,C_rel_usr=1.4706194473212615`。

本地最终代码门禁采用 H4/h12、8 threads、4 个解点：wall 13.90 s、峰值约
1.94 GiB、无 swap；exact relative energy 为
`1.8345,3.0356,7.3649,0.47126`。前三点是保留的 Helmholtz pre-asymptotic 段，第四点
开始离开该区间，不能删除前三点后伪造单调曲线。累计瓶颈为 LOD corrector 5.18 s、
localization certificate 4.21 s、candidate/RT2 2.92 s；这些阶段已经 patch 并行，精确误差
积分也已并行并从 method time 单列。

首次 target-gap=9 refresh 后的 h 约 44k unknowns；8-thread WSL 的 corrector/Gram 峰值
接近 10.2 GiB，因此本地只做 4 点门禁。服务器依次运行：

1. `H4/h12/gap 4->9/5-step` refresh factor；
2. 同口径 8-step pilot；
3. 12-step main；epoch 安全上限设为 13，避免保护上限先于 H-step 预算误停。

配置和命令见 `E1_E2_REVISED_ALGORITHMS_SERVER_2026-08-24.md`。下列 2026-08-23
记录仅保留旧轨迹审计，不再定义新版主实验。

E1 主实验额外冻结 epoch 0 网格审计：每个已求解 H-step 都在
`mesh_manifest.csv` 中记录 coarse/reference/candidate 三元组、`H_step` 与
`reference_mesh_version`。同一 epoch 的 reference 只写一份 VTU，各 H-step 的 reference
行复用该文件；绘图仍在每一列显示 reference，并用 version 与 SHA256 双重确认它没有变化。
只绘制第一个 epoch，避免后续 epoch 的重复图和存储开销。操作见
`E1_REVISED_MAIN_SERVER_2026-08-24.md`。

#### E1 历史执行状态（2026-08-23，已降级）

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

#### E1 统一 H6/小 reserve 重测（2026-08-25）

最新版论文对 Algorithm 1 和 R1 制造解没有新的数学改动，因此本轮不改变 estimator、
certificate、epoch 语义或误差口径，只重新冻结更经济且更长的数值轨迹。五种比较方法统一
从 H6 开始；新增 SLOD 作为 standard uniform LOD 补充对照。自适应方法取
`theta_H=0.1`，PALOD 取 H6/h12、`theta_c=0.3`、trigger/target gap 2/6、DirectSchur、
最多 36 个 H-step；SLOD 取 H6/h10、固定 gap 4、`ell=3`、十次同步加细至 H16/h20。

本地 PALOD trigger/target 2/6 门禁为 22.57 s、峰值 4,144,700 KiB、零 swap，形成
5 点和 4 点两个可拟合 epoch，exact-error/DoF 指数分别为 3.38 和 1.51；相比旧 H4/
target-gap9 主实验的 1:56:00 和 315,573,716 KiB，确认 reference reserve 是原有主要
资源放大因素。SLOD 五步 ell=3/4 对照的最大相对误差差为 0.137%，尾部指数为
0.754/0.755，而 ell=3 将墙钟从 24.38 s 降到 12.68 s、峰值从 1,250,572 KiB 降到
766,296 KiB，因此生产采用 ell=3。AFEM/UFEM H6 本地门禁尾部指数为 0.532/0.518。

服务器模式为 `e1-revised-h6-pilot` 和 `e1-revised-h6-main`；main 按 AFEM、UFEM、
SLOD、fixed LOD、PALOD 从短到长串行计时，避免并行污染论文 wall/RSS。完整命令、
门禁和绘图见 `E1_UNIFIED_H6_GAP6_SERVER_2026-08-25.md`。

### E2：L-shaped low-regularity 主实验，\(\kappa=16\)

比较：

1. singularity-aware reference-epoch LOD；
2. standard reference-epoch LOD；
3. adaptive \(P_1\) FEM；
4. fixed-oversampling LOD 仅在成本可承受时补做。

主图：corner 附近最终网格和 error-vs-DOF；不再绘制 cumulative wall-time 图。
matching region 中省掉的 corrector work 与分阶段时间只保留在 JSON/CSV 审计中。
必须使用论文的 mixed \(\Gamma_D/\Gamma_R\) 边界；全 Robin surrogate 不能写成正式
L-shaped 结果。

#### E2 无径向 cut-off 最终重基线（2026-08-25）

> 2026-08-25 后续审计：最新版论文 Algorithm 2 明确要求 candidate
> 在 `Omega_F` 与 regular region 分别满足同一个 `theta_c` Dörfler
> 下界。早期 cutoff-free 主轨迹误用了一个全局最小标记集；23 次 candidate
> 更新中有 10 次没有标记任何 `Omega_F` 单元，最终网格呈现角点中心偏粗、
> 外围 matching interface 楔形加细的结构，末 4/6 点指数降至约
> `0.417/0.437`。该轨迹保留为反例审计，不再作为 Algorithm 2 主结果。
>
> runner 已恢复 `Omega_F/regular` 独立 bulk marking。H6/h10、
> `theta_H=theta_c=0.2` 的 16-step 本地门禁峰值约 5.64 GiB、零 swap，
> 末 4/6/8/12 点 exact-error 指数分别约
> `0.593/0.616/0.603/0.671`，且每一步两区都达到各自 bulk 下界。
> corrected 24-step 服务器主实验已通过：最终 `N_H=10154`、exact relative
> energy error `0.0148296`；末 4/6/8/12 点指数分别为
> `0.521/0.501/0.514/0.575`，相关系数均不大于 `-0.986`，最后八点误差
> 严格下降。峰值约 36.2 GiB、零 swap、墙钟约 18 分 02 秒。最终网格不再
> 出现旧轨迹的角点粗洞和长楔形 refinement tongue，因此该结果取代
> `68a3e44` 的 global-marking 轨迹，作为 Algorithm 2 的 E2 主结果。

最新版论文将制造解冻结为

\[
u=b_\partial r^{2/3}\sin(2\theta/3)
+0.25\,C_{\rm osc}xyb_\partial
e^{-25((x+1/2)^2+(y-1/2)^2)}e^{i\kappa(x+1/2)},
\qquad
b_\partial=(1-x^2)^2(1-y^2)^2.
\]

不再使用 radial cut-off；解析 forcing 按论文的 singular/oscillatory split 计算。配置新增
`singular_solution_profile=boundary-weight-gaussian` 并纳入 canonical identity；旧
`radial-cutoff` profile 仅为历史结果兼容保留。单元测试已覆盖 mixed homogeneous boundary、
有限差分梯度、远离角点的 PDE residual、角点低正则性和“无正半径 transition annulus”。

本地网格门禁表明，AFEM 与 moving PALOD 的加细同时集中于凹角 `(0,0)` 和 Gaussian
中心 `(-0.5,0.5)`，未再出现旧 cut-off 环带。论文冻结 `R_*=0.125`；H6/h10、8-step
moving pilot 用时 19.1 s、峰值 1.30 GiB、零 swap，误差由约 0.160 降至 0.0665，
因此没有采用较小半径的必要。

四种正式比较统一从 H6 开始并取 `theta_H=0.2`：

1. AFEM：局部 level 上限 24，本地得到 38 个解点，尾部 12 点指数约 0.520；
2. moving PALOD：H6/h10、`R_*=0.125`、`ell0=2`、`theta_c=0.2`，最多 36 步；
3. standard uniform LOD (SLOD)：固定 gap 4、`ell=2`，由 ambient/coarse 资源门禁
   限制，`maximum_H_steps=10` 只是上限；
4. standard reference-epoch PALOD：H6/h12、trigger/target gap 2/6、DirectSchur、
   `ell0=2` 且跨 epoch 继承、`theta_c=0.2`，最多 24 步。

性能冻结依据：积分阶从 16/20/32 降至 12/16/24、递归深度从 8 降至 6，使 SLOD
wall time 下降约 30%，最大相对 exact-error 改变仅 `5.1e-6`；SLOD `ell=3->2`
使 wall time 再降约 25%、峰值内存降约 37%，短段尾部指数仍为 0.652；fixed-LOD
六步 probe 约加速 2 倍。`theta_c=0.3->0.2` 也降低两类 PALOD candidate work，短段
没有观测到掉阶。standard PALOD 的主要剩余瓶颈为 RT2 candidate flux，生产已经使用
patch OpenMP、summary audit、DirectSchur、较低积分阶和 `ell=2`；它在服务器序列中最后运行。

服务器模式为 `e2-cutofffree-revised-pilot` 与 `e2-cutofffree-revised-main`，顺序均为
AFEM、moving PALOD、SLOD、standard PALOD。为保持论文 wall/RSS 可解释，方法间串行，
方法内部保留 16-thread patch 并行。完整命令、门禁和绘图见
`experiments/helmholtz_adaptive_paper/E2_CUTOFF_FREE_SERVER_2026-08-25.md`。

#### E2 新版论文重基线（2026-08-24，已被 2026-08-25 无 cut-off 版本取代）

> 本节 radial-cutoff/tensor-bump 制造解与暂停状态均已失效，只作为历史审计保留，
> 不得与 2026-08-25 正式 E2 数据混用。

新版制造解为 C-infinity radial cutoff 的
`r^(2/3) sin(2 theta/3)` 奇异项，加上振幅 `B=0.05`、支撑位于
`(-0.75,-0.25)x(0.25,0.75)` 的 C-infinity tensor bump wave；`R_*=0.125`。
区域采用新版 `Omega_S=N^ell_S(S), Omega_F=N^ell(Omega_S)`，不再使用旧
`ell_S>=ell/Omega_F=N^(2ell_S)` 口径。

真实 H3/h8 回归发现论文当前 uniform regular reserve 与 exact matching 在 NVB 界面
一般不兼容；无限追赶不是性能 bug。实现采用前述 adaptive conformity collar，并以
`run.json: manuscript_conformance=implementation-erratum` 和 `hybrid_reserve.csv` 显式保留。

closure-aware full-regular Dörfler 校准结果：

- `theta_H=0.5`：完整 epoch 最多 2 点；
- `theta_H=0.3`：仍最多 2 点；
- `theta_H=0.1`：至少两个完整 epoch 各有 3 点，且不改变原始完整 regular mass 的
  Dörfler 下界。

因此 E2 hybrid 冻结 `theta_H=0.1,theta_c=0.5`。h8/8-point pilot 的 exact error 为
`0.5455,0.4839,0.4666,0.4609,0.3792,0.3436,0.3257,0.2989`；两个三点 epoch 的
log(error)--log(DOF) 斜率约 `-0.47` 与 `-1.16`。后者仍是短段暂态；最终是否稳定接近
`N^-1/2` 只能由 h12 服务器深轨迹判定。

服务器顺序冻结为 hybrid/standard refresh factor，随后 8-step pilots；通过后才运行
hybrid 15-step（epoch 安全上限 16）、standard 12-step（安全上限 13）、AFEM 和可选
fixed LOD。factor/pilot/main 的 `.done` 只在零 swap、matching/reserve、Dörfler、物理半径、
patch cap、ell 继承及相应三点 epoch 自动门禁通过后生成。下列旧状态
仅保留历史审计。

#### E2 历史启动状态（2026-08-23，已降级）

- schema-v5 已显式携带并哈希 case-S 制造解参数；E2 固定
  \(\gamma=0\)、quintic cutoff 外半径 1、smooth-wave amplitude \(B=0.05\)，即非振荡
  corner singularity 为主、保留小振幅 \(e^{i\kappa x}\) 光滑波。凹角边继续使用
  Dirichlet、外边 Robin，禁止退回历史默认 \(\gamma=1,B=0\)。
- 当时的 hybrid 区域曾不再把 (l_s) 直接等同于 corrector \(\ell\)，但仍选择满足
  (l_s\ge\ell) 且 (B_{R_*}(z)\cap\Omega\subset N^{2l_s}(z)) 的最小 (l_s)；该定义现已被
  上文新版 `Omega_S=N^ell_S(S),Omega_F=N^ell(Omega_S)` 完全取代，不得用于新版门禁。
  `iterations.csv` 与
  `corrector_work.csv` 记录实际 (l_s)、请求半径和离散保证覆盖半径。
- hybrid 的 level-gap guard 改为 regular region 中的最小**正** gap。\(\Omega_F\) 内
  \(H=h\) 的零 gap 是 matching condition，不是 reference exhaustion；若 closure 真正要求
  细化零-gap 单元，严格 proposed-to-reference containment 仍会立即触发 structural refresh。
- 更新后的固定半径 `H3/h8/4-step` hybrid probe（8 threads）wall 约 67 s、峰值 RSS
  约 3.16 GiB、无 swap；exact relative energy 依次为
  `0.10921, 0.08868, 0.06430, 0.04687`。epoch 0 在粗网格取 (l_s=2)，refresh 后
  自动取 (l_s=11)，离散保证覆盖半径约 0.254，满足 (R_*=0.25)。相比旧的
  (l_s=\ell=2) 探针，额外 wall time 来自固定圆内 coarse/reference matching，未出现
  新的内存或求解瓶颈；因此不在服务器主实验前引入未经验证的新 cache。
- 同口径 standard probe 表明 H3 初始粗空间处于强预收敛区；证书自动把 ell 提至 7，
  四点 exact error 仍由 `4.46` 降至 `2.01`。这条轨迹作为消融对照保留，不与 hybrid
  共用固定 ell，也不用于声称低正则性最优阶。
- `R_*=0.25` 的 `H3/h12` 服务器首轮在 matching 后形成病态大的 transition corrector
  patch：运行约 7 小时仍由单线程 SparseLU 占据，15 个 OpenMP worker 等待，RSS 约
  36.6 GiB。这是固定物理圆与 h12 matching 的尺度组合问题，不是内存不足。该任务已
  停止，E1 完成结果保留。`0.25` 只作为失败的压力测试记录，不再视为冻结参数。
- schema-v5 hybrid 配置增加 `hybrid_maximum_corrector_patch_fine_elements`。corrector 前
  必须输出 matching 时间、(l_s)、最大和 p95 patch fine-element count；最大值超过
  `100000` 时立即失败，禁止再次无进度运行数小时。
- 提交 `2fd2ec3` 的 h12 单步 pilot 已完成：`R_*=0.0625/0.125` 的 wall 分别为
  `13.98/14.84 s`，峰值 RSS 分别约 `2.095/2.099 GiB`，两者最大 patch 均为 3360；
  exact relative energy 为 `0.100152/0.100157`，而 skipped work units 为
  `2200/6408`。二者误差与成本基本相同，较大半径省掉更多奇异区 corrector，故正式
  E2 hybrid 冻结 `hybrid_minimum_physical_radius=0.125`。
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
marked_H, hybrid_regular_indicator_mass, hybrid_admissible_indicator_mass,
hybrid_marked_H_indicator_mass, hybrid_coarse_conformity_collar,
hybrid_coarse_marking_closure_safe, hybrid_full_regular_doerfler,
marked_c, active_correctors, rebuilt_correctors, reused_correctors,
skipped_correctors, time_corrector_check_total, time_lod_build_total,
time_lod_mesh_and_interpolation, time_lod_operators, time_corrector,
time_lod_corrected_basis, time_lod_coarse_operator,
time_lod_coarse_factorization, time_theta, time_lod_solve, time_reference_riesz,
time_candidate_enrich_total, time_candidate_close, time_candidate_operator_assembly,
time_candidate_prolongation, time_candidate_flux_reconstruction,
time_candidate_flux_prepare, time_candidate_flux_patch_solve,
time_candidate_flux_merge, time_candidate_flux_audit,
candidate_flux_parallel_threads, time_candidate_enrich,
time_candidate_dual, time_candidate_dual_operator_assembly,
time_candidate_dual_load_assembly, time_candidate_dual_prolongation,
time_candidate_dual_solve, time_candidate_dual_prepare,
time_candidate_dual_patch_solve, time_candidate_dual_reduction,
candidate_dual_patch_factorizations, candidate_dual_parallel_threads,
time_mesh, time_validation_cumulative, time_artifact_capture_cumulative,
time_method_cumulative, time_total_cumulative,
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
- `corrector_work.csv`：active/rebuild/reuse/skip 和完整 LOD build/Gram 分阶段时间；
- `hybrid_reserve.csv`：Algorithm 2 每次 collar trial/closure 的 spill、profile margin、far
  gap、mesh size 和时间；非 hybrid 运行仍生成只有 header 的文件。

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
