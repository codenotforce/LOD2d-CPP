# Helmholtz 认证自适应 LOD：论文数值实验交付计划

> 状态日期：2026-08-07
>
> 唯一终点：完成论文数值实验章节所要求的全部算例、对照方法、认证量、误差与工作量统计、图表和可复现数据。
>
> 理论与实验合同来源：[helmholtz_lod_certified_amsart.tex](../LOD_paper/helmholtz_lod_certified_amsart.tex)，重点对应 Certification-driven adaptivity 和 Adaptive numerical experiments 两节。
>
> 当前实现事实：项目已经具备单位正方形、固定 \(h,\ell\)、仅自适应 \(H\) 的 Stage-1 基线；尚不能把当前 strong-residual proxy 的结果称为 certified adaptive LOD。

本文件取代 2026-07-25 版本的旧主线。原计划中的 fine/mixed/macro 强残差、经验 reliability envelope、高对比系数、局部 \(\ell_T\) 和独立 patch 细网格不再阻塞论文实验；它们仅作为诊断或论文完成后的扩展。正式实现严格跟随论文已经冻结的

\[
\eta_H,\qquad
\Theta_{\mathrm{tot}},\qquad
\Theta_h,\qquad
\widehat\delta_{\mathrm{tot}},\qquad
\widehat\delta_h,\qquad
\widehat\delta_\ell,\qquad
q_{\mathrm{tot}},q_h,q_\ell
\]

及四步自适应顺序。

## 1. 完成定义与声明边界

### 1.1 计划完成不等于论文实验完成

只有同时满足下列条件，项目才可标记为“论文数值实验完成”：

1. R1、R2、S、K 四组实验全部完成，且正式波数为 \(\kappa\in\{8,16,32\}\)。
2. certified adaptive LOD、fixed-\(h,\ell\) 的 \(H\)-adaptive LOD、两种 standard LOD 配置、uniform \(P_1\) FEM、adaptive \(P_1\) FEM 均由同一比较驱动运行。
3. 每种方法都按共同相对能量误差目标

   \[
   10^{-1},\quad 5\cdot10^{-2},\quad 2\cdot10^{-2},\quad 10^{-2}
   \]

   取第一次达到目标的迭代；未达到时保留 censored/work-limit 状态，不能删除失败样本。
4. 所有称为 certified 的运行均使用经过验证的常数界、谱包络和 inf-sup 下界；普通 Lanczos/SVD 近似只能产生 diagnostic 或 conditional 结果。
5. 原始数据、配置、编译信息、硬件信息、日志、网格快照、聚合表和绘图脚本形成完整追溯链。
6. 能从干净构建开始，用统一入口生成论文需要的表格和图片；论文中的占位表不再含破折号。

### 1.2 结果状态必须显式区分

| 状态标签 | 含义 | 允许的论文表述 |
|---|---|---|
| implementation-study | 当前 strong-residual proxy 或未完成证书链的工程结果 | 只能称实现研究或 baseline |
| conditional | 公式实现完整，但一个或多个常数、特征值或奇异值只有非验证近似 | 条件认证，不得简称 certified |
| audit-certified | 严格包围 \(\|u_{\widehat h}-U\|_\kappa\) | 可称 audit-space certified |
| continuous-certified | audit 区间与独立连续离散误差区间均经过验证 | 可称 continuous-error certified |
| empirical-reference | 只通过两级 audit saturation 检查参考解 | 可用于误差比较，不得称连续误差证书 |

每个运行目录、每一行汇总数据和每张论文表都必须包含上述状态字段。

### 1.3 不属于当前论文完成门槛的内容

- 空间变化的 oversampling \(\ell_T\)；当前论文只允许全局 \(\ell\leftarrow\ell+1\)。
- 每个 patch 彼此独立且不共形的细网格。
- 基于“预计误差下降/预计成本”的 \(H/h/\ell\) 三路竞争决策。
- 七常数经验 reliability envelope、经验 \(k\)-权重、高对比系数和 PML。
- 自适应收缩、准最优复杂度或无污染定理。

这些内容可以保留为后续研究，但不得推迟本文件中的论文实验关键路径。

## 2. 论文实验合同

### 2.1 必做算例

| ID | 区域与数据 | 正则性 | 正式参数 | 实验目的 |
|---|---|---|---|---|
| R1 | \(\Omega=(0,1)^2\)，光滑制造解 \(u_{\mathrm{reg}}=\phi(x_1)\phi(x_2)e^{\mathrm i\kappa x_1}\)，\(\phi(t)=16t^2(1-t)^2\) | 光滑 | \(\kappa=8,16,32\) | 误差、残差恒等式、effectivity、audit 校准和无虚假网格集中 |
| R2a | 单位方形，\(L^2\)-归一化 Gaussian，\(x_0=(0.35,0.55)\)，\(\sigma=2^{-5}\) | 光滑但局部化 | \(\kappa=8,16,32\) | 普通正则性下的自适应收益 |
| R2b | 与 R2a 相同，\(\sigma=2^{-6}\) | 光滑但更局部化 | \(\kappa=8,16,32\) | 对更尖锐局部源的鲁棒性 |
| S | \(\Omega_L=(-1,1)^2\setminus([0,1]\times[-1,0])\)，混合边界和已知 \(r^{2/3}\) 奇异制造解 | \(H^{1+2/3-\varepsilon}\)，通常不属于 \(H^{1+2/3}\) | \(\kappa=8,16,32\) | 高奇性下自适应粗网格是否恢复角点 |
| K | R2a、R2b、S 的波数汇总 | 依算例而定 | 汇总 \(\kappa=8,16,32\) | 同时比较误差、稳定裕量和总工作量随波数的变化 |

R2 的源必须是

\[
f_\sigma(x)
=c_\sigma\exp\!\left(-\frac{|x-x_0|^2}{2\sigma^2}\right),
\qquad
\|f_\sigma\|_{L^2(\Omega)}=1.
\]

当前代码中的 \(\exp(-80r^2)\) 和 \(\exp(-40r^2)\) 均不等价于正式 R2，不能继续进入论文数据。

S 算例固定角度约定 \(0<\theta<3\pi/2\)，凹角两条边为 \(\Gamma_D\)，其余外边界为 \(\Gamma_R\)，并使用

\[
u_{\mathrm{sing}}(r,\theta)
=\chi(r)r^{2/3}\sin(2\theta/3)e^{\mathrm i\kappa x_1}.
\]

生产配置使用真正的 \(C^\infty\) cut-off。定义

\[
\psi(t)=
\begin{cases}
0,&t\le0,\\
\exp(-1/t),&t>0,
\end{cases}
\]

并冻结

\[
\chi(r)=
\frac{\psi(r_1-r)}
{\psi(r_1-r)+\psi(r-r_0)},
\qquad
r_0=\tfrac14,\quad r_1=\tfrac12.
\]

因此 \(\chi=1\) 于 \(r\le r_0\)，\(\chi=0\) 于 \(r\ge r_1\)，且过渡在 \((r_0,r_1)\) 内光滑。该 cut-off 的支撑不接触外侧 Robin 边界。实现时必须同时提供值、梯度、Laplacian 和在 \(r=0\) 附近稳定的求值路径。

### 2.2 必较方法

| 方法 ID | 方法 | 参数规则 | 当前状态 |
|---|---|---|---|
| CALOD | 完整 certified adaptive LOD | 按论文顺序联合更新 \(H,h,\ell\) 和 audit space | 未实现 |
| HLOD | fixed-\(h,\ell\) 的 \(H\)-adaptive LOD | 使用与 CALOD 相同的 \(\eta_H\)，冻结下述 \(h_{\mathrm{prior}},\ell_{\mathrm{prior}}\) | 当前只有 strong-residual proxy |
| SLOD-prior | quasi-uniform standard LOD | 使用同一 \(h_{\mathrm{prior}},\ell_{\mathrm{prior}}\) | 核心求解已实现，缺统一驱动 |
| SLOD-matched | quasi-uniform standard LOD | 使用同一算例 CALOD 在成功或资源停止前全部有效历史中的最大 fine resolution 和最大 \(\ell\) | 未实现自动匹配 |
| UFEM | uniform conforming \(P_1\) FEM | 一致 NVB 加细 | 核心求解已实现，缺统一驱动 |
| AFEM | adaptive conforming \(P_1\) FEM | volume、flux-jump、impedance-boundary 残差及相同 \(\theta_H\) | 未实现 |

现有 strong-residual H-only 驱动保留为 HLOD-proxy，仅用于回归、诊断和展示当前实现演化；论文最终主比较中的 HLOD 必须使用同一 audit-kernel \(\eta_H\)，从而只隔离“是否自适应 \(h,\ell\)”这一因素。

HLOD 与 SLOD-prior 的固定参数不得使用正式误差回调。令 \(h_{\mathrm{prior}}\) 为共同初始 fine hierarchy 中满足下列条件的第一个 uniform level：

\[
\kappa h_{\max}\le c_h,\qquad c_h=\tfrac18,
\]

并额外要求 R2 中 \(h_{\max}\le\sigma/4\)，S 中 \(h_{\max}\le r_0/16\)。取

\[
\ell_{\mathrm{prior}}
=\left\lceil\log_2\kappa\right\rceil .
\]

这两个参数逐 case、逐 \(\kappa\) 写入冻结的 baseline_parameters.csv。SLOD-matched 读取 CALOD 在成功返回或 work/memory/time limit 停止前全部有效历史中的最大 local fine level 和最大 \(\ell\)，并把该 fine level 作为 uniform level。CALOD 被 censored 时 SLOD-matched 仍运行并继承 censored-source 标志；只有 CALOD 在产生任何有效 corrector 参数前失败时才记为 unavailable，禁止人工补参数。

### 2.3 公平比较规则

1. 所有方法使用同一 PDE 数据、边界标签、积分规则、复数约定、线性代数容差和公共 evaluation reference；CALOD 另有算法内部 cert-audit，但不能替代公共 reference。
2. CALOD 与 HLOD 使用相同 \(\theta_H\)；AFEM 也使用同一数值，便于解释标记敏感性。
3. 主实验固定 \(\theta_H=0.5\)。敏感性实验对 R1、R2a、R2b、S 的全部 \(\kappa=8,16,32\) 使用 \(\theta_H=0.3,0.7\)；若未来只保留代表波数，必须先同步修改论文实验协议。
4. 六种配置共享同一个 case-specific boundary-fitted 初始网格；CALOD、HLOD、AFEM 从该网格局部加细，SLOD/UFEM 使用它的 uniform descendants。所有方法共享相同的 case/\(\kappa\) work、memory 和 wall-time 上限定义。
5. 比较按共同 evaluation reference 给出的相对能量误差 \(E_\kappa^{\mathrm{ref}}\) 目标，不按相同迭代数或人为相同 DOF 截取；R1、S 另外报告 exact error，但不改变达标判据。
6. 只有公共 evaluation-reference 全局求解时间单列并排除；CALOD 内部 cert-audit 网格、矩阵、分解和 Riesz 计算全部计入 CALOD。
7. 被测方法的累计时间必须包含达到目标前的所有迭代、所有失败的内部 \(h/\ell\) 认证循环以及必要的重建。
8. 相同硬件、线程数、编译选项、求解容差和缓存策略下重复正式计时，报告中位数；默认重复 5 次。
9. exact solution 和 evaluation-reference error 只用于评价，绝不能进入 MARK 或算法停止逻辑；CALOD 只能读取论文定义的 cert-audit 证书。

## 3. 当前能力与关键差距

| 能力 | 状态 | 现有证据 | 距离论文交付 |
|---|---|---|---|
| NVB、稳定粗单元 ID、固定 master fine mesh、\(V_H\subset V_h\) | 已实现 | include/helmholtz/adaptive/hierarchy.h，tests/test_helmholtz_adaptive.cpp | 需扩展为 \(V_H\subset V_h\subset V_{\widehat h}\) |
| 拟插值和自然延长，\(I_HP_{Hh}\approx I\) | 已实现 | Helmholtz model 构建检查 | 需同时覆盖 mixed boundary、非均匀 \(V_h\) 和 audit space |
| R1 制造解、standard LOD、uniform FEM | 部分实现 | manufactured benchmark 和 \(k\)-scan | 缺统一方法注册、共同目标和累计工作量 |
| fixed-\(h,\ell\) H-only 自适应 | 已实现 Stage-1 | src/helmholtz/adaptive/driver.cpp | 当前估计子不是论文 \(\eta_H\) |
| fine/mixed/macro 强残差 | 已实现诊断 | src/helmholtz/adaptive/estimator.cpp | 只保留 HLOD-proxy 和 AFEM 复用，不是 CALOD 主估计子 |
| 局部无约束能量 Riesz | 已实现诊断 | estimator.cpp 的 dual calibration | 必须改为 audit-kernel 受约束 Riesz |
| R2 Gaussian | 部分实现 | adaptive 和 \(k\)-scan 各有一种 Gaussian | \(\sigma\) 与归一化均不符合论文 |
| L 型区域与 mixed Dirichlet/Robin | 未实现 | TriMesh 有 Dirichlet 节点字段 | 当前 Helmholtz model 明确拒绝 Dirichlet 节点，是 S 的硬阻塞项 |
| \(\eta_{H,z}\)、\(\eta_H\)、局部效率诊断 | 未实现 | 无 | CALOD 和最终 HLOD 的前置条件 |
| \(G_{\mathrm{tot}},G_h,\Theta_{\mathrm{tot}},\Theta_h\) | 未实现 | 无 | \(h/\ell\) 决策的前置条件 |
| verified eigen/SVD enclosure 和 inf-sup 下界 | 未实现 | 仅有小系统近似 inf-sup | 无此项只能输出 conditional |
| 局部 corrector-\(h\) 更新和全局 \(\ell\) 更新 | 未实现 | 当前 \(h,\ell\) 全程固定 | CALOD 内循环缺失 |
| adaptive \(P_1\) FEM | 未实现 | NVB、Dörfler 和残差装配可复用 | 必做 comparator |
| matched-tolerance 比较驱动 | 未实现 | benchmark 输出格式彼此独立 | 无法生成论文主表 |
| 累计 patch dimension、峰值内存、core time | 未实现或不完整 | 有阶段 wall time，patch_cost 仅用于调度 | 必须建立统一数据协议 |

旧计划中的一个错误验收必须删除：当 \(\mathcal T_H=\mathcal T_h\) 时，\(u_h-U\)、代数残差和 kernel residual 可以趋于机器精度，但 strong residual 仍含连续 FEM 离散误差，不要求趋零。

## 4. 冻结的数学与算法合同

### 4.1 模型、边界和三层空间

论文正式实验固定 \(A=I\)、\(n=1\)、\(\beta=1\)，统一模型为

\[
-\Delta u-\kappa^2u=f
\quad\text{in }\Omega,
\]

\[
u=g_D\quad\text{on }\Gamma_D,\qquad
\nabla u\cdot\nu-\mathrm i\kappa u=g_R
\quad\text{on }\Gamma_R.
\]

R1、R2 使用纯阻抗边界；S 使用凹角两边 Dirichlet、其余边 Robin。正式三个 case 均为齐次边界数据 \(g_D=g_R=0\)。全局 FEM、LOD corrector、局部 Riesz、audit、强残差和 AFEM 必须共享同一 boundary-edge 分类，禁止某个模块把 \(\Gamma_D\) 误装成 Robin。

项目可以继续保留一般 \(A,n,\beta\) 的接口，但这些系数在本计划的生产配置中不得变化；高对比和变系数实验属于第 13 节扩展。

离散空间必须始终满足

\[
V_H\subset V_h\subset V_{\widehat h}\subset V.
\]

- \(V_H\)：自适应粗空间。
- \(V_h\)：corrector 离散空间；第一版使用一个全局共形但可局部加细的空间，不采用彼此独立的 patch 网格。
- \(V_{\widehat h}\)：独立 audit 空间，也是所有可计算 kernel Riesz 问题的 ambient space。

每次网格变化都重新验证

\[
\|I_HP_{Hh}-I\|_F\le10^{-9},
\qquad
\|I_HP_{H\widehat h}-I\|_F\le10^{-9},
\]

并记录网格、插值、边界和 corrector space 的版本号。

### 4.2 论文主估计子

在

\[
W_{\widehat h}=\ker(I_H|_{V_{\widehat h}})
\]

上先冻结具体 patch policy。设 \(p_I\) 是拟插值支撑传播的已验证粗层数，取

\[
D_z=\omega_z^{m_D},
\qquad m_D=p_I+1,
\]

其中 \(\omega_z^{m_D}\) 按共享粗顶点的层邻接扩张，并在边界处截断；若最小 \(m_D\) 不能同时包含 \(\operatorname{supp}\lambda_z\) 和 \(I_H\) 的传播支撑，则模型构建直接失败。扩大 patch 严格按论文定义

\[
D_z^+
=\operatorname{int}\!\bigcup_{T\subset\overline{D_z}}
\overline{\omega_T^\dagger}.
\]

patch adjacency、\(p_I,m_D,\omega_T^\dagger\) 的构造和边界截断规则都进入 constants-registry hash。随后在粗节点局部子空间 \(W_{\widehat h,z}\) 上求受约束 Riesz 代表元

\[
b_\kappa(\xi_z,w)=R_U(w)
\qquad\forall w\in W_{\widehat h,z},
\]

\[
\eta_{H,z}=\|\xi_z\|_\kappa,\qquad
\eta_H^2=\sum_z\eta_{H,z}^2.
\]

局部代数系统使用 kernel basis 或满行秩约束后的 saddle-point 形式

\[
\begin{bmatrix}
B_{\kappa,z}&I_{H,z}^*\\
I_{H,z}&0
\end{bmatrix}
\begin{bmatrix}\xi_z\\\boldsymbol\Lambda_z\end{bmatrix}
=
\begin{bmatrix}r_z\\0\end{bmatrix}.
\]

这里 \(\boldsymbol\Lambda_z\) 是约束乘子，与粗帽函数 \(\lambda_z\) 无关。定义粗节点的单元重数

\[
m_z:=\#\{T\in\mathcal T_H:z\in T\}.
\]

粗节点量无重复分配到单元：

\[
\eta_{H,T}^2
=\sum_{z\in\mathcal N_H\cap T}\frac{\eta_{H,z}^2}{m_z},
\qquad
\sum_T\eta_{H,T}^2=\eta_H^2.
\]

随后用 \(\eta_{H,T}\) 做 \(\theta_H\)-Dörfler 标记。R1、S 及已计算 audit 解的 R2 还要报告

\[
\mathcal I_{\mathrm{loc},z}
=\frac{\eta_{H,z}}
{\|u_{\widehat h}^{\mathrm{cert}}-U\|_{\kappa,D_z^+}},
\]

但该局部真误差不能参与标记。

### 4.3 Corrector 和稳定性证书

必须装配 coarse energy matrix \(M_\kappa\)、总 corrector Gram matrix \(G_{\mathrm{tot}}\) 和 fine-corrector Gram matrix \(G_h=\sum_TG_{h,T}\)，并计算

\[
\Theta_{\mathrm{tot}}^2
=\lambda_{\max}(G_{\mathrm{tot}},M_\kappa),
\qquad
\Theta_h^2
=\lambda_{\max}(G_h,M_\kappa).
\]

由此构造并记录

\[
\underline\delta_{\mathrm{tot}}
=\frac{\Theta_{\mathrm{tot}}}{C_aC_{\mathrm{ov}}},
\qquad
\widehat\delta_{\mathrm{tot}}
=\frac{C_{\mathrm{sd}}}{c_W}\Theta_{\mathrm{tot}},
\]

\[
\widehat\delta_h
=\frac{C_{\mathrm{ol}}(\ell)^{1/2}}{c_W}\Theta_h,
\]

\[
\underline\delta_\ell
=\max\{0,\underline\delta_{\mathrm{tot}}-\widehat\delta_h\},
\]

以及

\[
\widehat\delta_\ell
=\min\!\left\{
\widehat\delta_{\mathrm{tot}}+\widehat\delta_h,\,
C_{\mathrm{loc}}C_{\mathrm{ol}}(\ell+s)^{1/2}\beta^\ell
\right\}.
\]

只有在另行验证 polynomial patch-growth 条件后，第二项才可换成多项式 overlap 形式。随后计算

\[
q_{\mathrm{tot}}
=\frac{C_\Pi C_{\mathrm{Fort}}}
{\gamma_{V_{\widehat h}}(\kappa)}
\widehat\delta_{\mathrm{tot}},
\qquad
q_h
=\frac{C_\Pi C_{\mathrm{Fort}}}
{\gamma_{V_{\widehat h}}(\kappa)}
\widehat\delta_h,
\qquad
q_\ell
=\frac{C_\Pi C_{\mathrm{Fort}}}
{\gamma_{V_{\widehat h}}(\kappa)}
\widehat\delta_\ell.
\]

Audit-space 主区间固定为

\[
L_{\mathrm{LOD}}
=\frac{\eta_H}{C_aC_{\mathrm{ov}}},
\qquad
U_{\mathrm{LOD}}
=C_{\mathrm{sd}}\left(
c_W^{-1}
+\frac{C_\Pi C_{\mathrm{Fort}}}
{\gamma_{V_{\widehat h}}(\kappa)}
\widehat\delta_{\mathrm{tot}}
\right)\eta_H.
\]

稳定裕量固定为

\[
\widehat\varepsilon
=C_\Pi\widehat\delta_{\mathrm{tot}},
\qquad
m_{\mathrm{stab}}
=\frac{\gamma_{V_{\widehat h}}(\kappa)}
{2C_{\mathrm{Fort}}}
-C_a\widehat\varepsilon(2+\widehat\varepsilon).
\]

只有 \(m_{\mathrm{stab}}\ge0\) 且 \(q_{\mathrm{tot}}\le q_0\) 时，corrector certification 才通过。

R1、R2、S 的正式离散均使用实网格、实插值权、实系数和共轭不变的边界子空间。实现必须验证 primal/adjoint corrector 的共轭关系及两侧证书相等；若 mixed-boundary 消元或任何后续扩展破坏该门槛，则分别装配 trial 和 test 两套 \(G_{\mathrm{tot}},G_h\) 与误差界，不能继续复用单侧证书。

主特征向量 \(x_h\) 给出

\[
\eta_{h,T}^2=x_h^*G_{h,T}x_h,
\qquad
\sum_T\eta_{h,T}^2=\Theta_h^2.
\]

它标记“哪个 coarse-source corrector 的细分辨率不足”，不是细单元 residual，也不是连续解离散误差指标。若主特征值成簇，必须聚合主谱簇，不能让任意特征向量选择改变标记结果。

普通特征值迭代和 SVD 结果只属于 conditional 模式。certified 模式还需：

- \(G_{\mathrm{tot}}\)、\(G_h\) 的 Hermitian/PSD 数值验证；
- 广义特征值上包络；
- \(\gamma_{V_{\widehat h}}(\kappa)\) 的验证下界；
- 线性系统和 Riesz 求解误差的可计入上界；
- 所有解析常数的来源、方向和版本。

### 4.4 论文版四步状态机

状态机固定为下列顺序，不允许改成三种指标同时竞争：

1. **COARSE_ADMISSIBILITY**
   计算

   \[
   \mu_H=C_{\mathrm{app}}\max_T\kappa H_T.
   \]

   若 \(\mu_H>\mu_0\)，标记全部违规单元、NVB 加细和共形闭包，维护三层嵌套后重新开始。

2. **CORRECTOR_CERTIFICATION**
   计算 \(\Theta_{\mathrm{tot}},\Theta_h\)、全部 \(\delta\) 界、\(q_{\mathrm{tot}}\) 和验证 inf-sup 下界。只有稳定条件通过且 \(q_{\mathrm{tot}}\le q_0\) 才能进入第 3 步。否则：

   - 若 \(\widehat\delta_h>\tau\underline\delta_{\mathrm{tot}}\)，按 \(\eta_{h,T}\) 做 \(\theta_h\)-Dörfler 标记，并在选中 corrector patch 覆盖区统一加细共形 \(V_h\)；
   - 否则执行全局 \(\ell\leftarrow\ell+1\)；
   - 更新受影响 corrector 后重复第 2 步。

3. **COARSE_ERROR_CONTROL**
   求 Petrov-Galerkin 解，计算 \(\eta_{H,z},\eta_{H,T},\eta_H\) 和 audit-space 区间。仅当任务只请求 audit-space error 且 \(U_{\mathrm{LOD}}\le\mathrm{tol}\) 时立即成功返回；其余情况，包括 continuous-error 任务下 \(U_{\mathrm{LOD}}\le\mathrm{tol}\)，都先形成并暂存 \(H\)-Dörfler 集，再进入第 4 步。该集合只有在第 4 步判断仍需 coarse refinement 时才实际加细。

4. **AUDIT_CONTROL**
   若要求连续误差，计算独立区间

   \[
   L_{\widehat h}\le
   \|u-u_{\widehat h}\|_\kappa
   \le U_{\widehat h}.
   \]

   当 \(U_{\widehat h}>\rho_{\mathrm{aud}}U_{\mathrm{LOD}}\) 时先加细 \(V_{\widehat h}\) 并回到第 2 步；否则形成

   \[
   L_{\mathrm{true}}
   :=\max\{0,L_{\mathrm{LOD}}-U_{\widehat h},
              L_{\widehat h}-U_{\mathrm{LOD}}\},
   \qquad
   U_{\mathrm{true}}
   :=U_{\mathrm{LOD}}+U_{\widehat h}.
   \]

   只有 \(U_{\mathrm{true}}\le\mathrm{tol}\) 时 continuous-error 任务成功。若仍未达到 tolerance，再执行之前形成的 \(H\)-标记并回到第 1 步。

任何 work、memory、iteration 或线性代数上限都必须返回结构化失败状态和完整历史。

### 4.5 Reference 与停止规则

- R1 和 S 有精确解；精确误差只用于验证和图表。
- R2 没有精确解。优先使用独立外部估计子；如果尚未实现，只能采用两级 audit 解差不超过被测误差 5% 的 empirical saturation 规则。
- 必须区分两类细空间：
  - \(V_{\widehat h}^{\mathrm{cert}}\) 是 CALOD 状态机内部的 certification audit space；其网格、矩阵、Riesz 问题和分解都是 CALOD 成本，不能从 method time 中剔除；
  - \(V_{\mathrm{ref}}\) 是参数冻结后为六种配置共同构造的 evaluation reference，只用于后处理误差和 matched-target 选择。它不参与任何方法的 MARK/STOP，其全局求解时间单列并从 method time 中排除。
- 证书验证量和共同比较量必须分别定义：

  \[
  e_{\mathrm{cert}}
  :=\|u_{\widehat h}^{\mathrm{cert}}-U\|_\kappa,
  \qquad
  E_\kappa^{\mathrm{ref}}
  :=\frac{\|u_{\mathrm{ref}}-U\|_\kappa}
  {\|u_{\mathrm{ref}}\|_\kappa},
  \qquad
  E_0^{\mathrm{ref}}
  :=\frac{\|u_{\mathrm{ref}}-U\|_{L^2}}
  {\|u_{\mathrm{ref}}\|_{L^2}}.
  \]

  \(\mathcal I_U=U_{\mathrm{LOD}}/e_{\mathrm{cert}}\)、local effectivity、\(L_{\mathrm{LOD}},U_{\mathrm{LOD}}\) 都属于 cert-audit；matched-target 和六方法主图只使用 \(E_\kappa^{\mathrm{ref}}\)。R1、S 再附 exact error。
- \(L_{\widehat h},U_{\widehat h}\) 专门包围 \(\|u-u_{\widehat h}^{\mathrm{cert}}\|_\kappa\)，不能用 evaluation-reference saturation 量替换。
- 当前论文数值章节仍用同一个 \(u_{\widehat h}\) 记号同时承担两种角色。正式导出前必须同步修改论文：证书部分写 \(u_{\widehat h}^{\mathrm{cert}}\)，公平比较部分写 \(u_{\mathrm{ref}}\)；或者在配置中证明两空间实际相同。默认采用修改记号的方案。
- 论文数值协议中的默认 audit 比例解释为

  \[
  \rho_{\mathrm{aud}}=0.05,\qquad
  U_{\widehat h}\le0.05\,U_{\mathrm{LOD}},
  \]

  不是逗号分隔的两个量。
- 生产 driver 的停止依据只能是证书、预设 tolerance 和资源上限，不能读取 exact/reference error。

## 5. 实施工作包

所有工作包都以 full rebuild 为正确性金标准。增量复用只能在对应数学结果逐项与 full rebuild 一致后启用。

### WP0：冻结实验协议和数据模式

**状态：已完成（2026-08-08）；配置协议、注册表、canonical hash/run ID 和严格往返测试已落地。**

任务：

1. 新增统一的 case registry：R1、R2a、R2b、S。
2. 新增统一的 method registry：CALOD、HLOD、SLOD-prior、SLOD-matched、UFEM、AFEM；保留 HLOD-proxy 但标为 diagnostic。
3. 定义版本化配置和输出 schema，禁止各 benchmark 自行发明列名。
4. 冻结共同目标、波数、\(\theta_H\)、重复次数、容差和失败状态。
5. 为每次运行生成不可变 run ID，内容至少包含 case、method、参数哈希、Git commit、编译哈希和重复编号。

建议产物：

- experiments/helmholtz_adaptive_paper/schema-v1.json
- experiments/helmholtz_adaptive_paper/cases/*.json
- experiments/helmholtz_adaptive_paper/methods/*.json
- include/helmholtz/experiments/paper_config.h
- src/helmholtz/experiments/paper_config.cpp

验收：

- 配置往返序列化不丢字段；
- 未识别字段和不同 schema version 明确报错；
- 相同配置产生相同 canonical hash；
- 所有方法共享同一 problem definition、quadrature 和 tolerance 对象。

### WP1：统一几何、边界和论文数据

**状态：已完成（2026-08-08）；显式边标签及局部/全局 NVB 标签守恒、mixed-boundary P1 FEM/LOD/residual、L 型网格、R1/R2/S 数据对象、统一高阶 QuadraturePolicy、边界标签 VTK 输出及积分收敛门槛均已完成。**

任务：

1. 将边界从“Dirichlet 节点列表”升级为显式 boundary-edge tags，至少支持 Interior、Dirichlet、Robin。
2. mixed boundary 标签在 NVB 加细、闭包、粗细映射、patch 截取和网格输出中保持不变。
3. 修改 Helmholtz operator、load、FEM solve、LOD corrector、局部 Riesz 和 residual，使它们只在 \(\Gamma_R\) 装配阻抗项，并正确消元 \(\Gamma_D\)。
4. 加入 L 型初始 NVB 网格，验证面积为 3、凹角为原点、边界分类和参考边兼容。
5. 实现 R1 的统一数据类，复用现有制造解但由 case registry 调用。
6. 实现带 \(\sigma\) 参数和 \(L^2\) 归一化的 R2；归一化常数用解析积分或高精度独立积分校验。
7. 实现 S 的 \(\chi\)、\(u_{\mathrm{sing}}\)、梯度、Laplacian、体源和边界数据。
8. 新增统一 QuadraturePolicy。R1 使用能精确/高精度解析多项式振荡项的阶数；R2 对窄 Gaussian 自动升阶或递归细分；S 对接触凹角的单元使用 Duffy 变换或等价奇异积分，其余单元使用统一高阶规则。全局 FEM、LOD、AFEM、audit、error 和 normalization 必须调用同一 policy。

建议产物：

- include/helmholtz/boundary.h
- include/helmholtz/benchmarks/paper_cases.h
- src/helmholtz/benchmarks/paper_cases.cpp
- tests/test_helmholtz_mixed_boundary.cpp
- tests/test_helmholtz_paper_cases.cpp

验收：

- R2a、R2b 的离散高阶积分满足 \(|\|f_\sigma\|_{L^2}-1|\le10^{-10}\)；
- R1 的 homogeneous impedance 检查延续通过；
- S 的两条凹角边全部且仅被标为 Dirichlet，外边界全部且仅被标为 Robin；
- S 的制造解在抽样点满足 PDE 和混合边界，误差与求积阶一致；
- mixed-boundary 全局矩阵与手工小网格参考一致；
- NVB 多轮加细后边界测度和标签守恒；
- R2 的 load、S 的 load/energy error 在连续两次升阶或细分后的相对变化低于预登记容差；当前硬编码 7 点三角公式在通过该收敛门槛前不得用于正式数据。

### WP2：三层网格、audit space 和统一误差模块

**状态：进行中（2026-08-08）；已显式维护 \(T_H,T_h,T_{\widehat h}\)、父子映射、三组节点 prolongation、fine/audit 两套 \(I_H\) 及按 fine-element 集合驱动的 cert-audit 共形局部加细；统一误差角色接口、局部误差服务和分离计时缓存待完成。**

任务：

1. 将 AdaptiveMeshHierarchy 扩展为显式维护 \(T_H,T_h,T_{\widehat h}\) 及父子关系。
2. \(V_h\) 使用全局共形局部加细；选中 corrector patch 时加细其覆盖的 fine elements 并做闭包。
3. cert-audit \(V_{\widehat h}\) 不读取 evaluation-reference error，只由 Step 4 的 audit 证书规则加细，且始终是 \(V_h\) 的共形加细。
4. 维护 \(P_{Hh},P_{h\widehat h},P_{H\widehat h}\)、\(I_H\) 和 kernel constraint restriction。
5. 建立统一细空间服务但分开两个角色：CALOD 内部可更新的 cert_audit，以及全部方法共享、只供评价的 evaluation reference。exact、外部 estimator 和 two-level saturation 通过显式角色接口调用。
6. 提供全局能量误差、\(L^2\) 误差和 \(D_z^+\) 局部能量误差。
7. 分别缓存 cert-audit 和 evaluation-reference factorization；前者计入 CALOD certificate time，后者记为 evaluation_reference_time 并单列排除。

验收：

- 每轮三层网格共形、无悬挂点并满足嵌套；
- 三个 prolongation 的复合误差不超过 \(10^{-10}\)；
- 两个 right-inverse 检查不超过 \(10^{-9}\)；
- exact/reference 两条误差计算在 R1、S 上相互收敛；
- evaluation reference 的结果不暴露给 MARK 和 STOP 接口；cert_audit 只通过论文证书量影响 CALOD 决策；
- cert_audit 成本进入 CALOD，evaluation reference 的全局求解时间单列。

### WP3：audit-kernel 局部残差估计子 \(\eta_H\)

**状态：未实现；替换旧计划的主估计子。**

任务：

1. 构造固定阶粗节点 patch \(D_z\)、扩大 patch \(D_z^+\) 及相应 audit DOF。
2. 从 \(I_H\) 提取局部约束，删除 inactive 和线性相关约束行。
3. 同时实现 kernel-basis 参考路径和 saddle-point 生产路径。
4. 从同一个全局残差限制得到每个局部右端，禁止再用 broken residual 的边归属替代。
5. 输出 \(\eta_{H,z}\)、\(\eta_H\)、\(\eta_{H,T}\) 和 Dörfler 集。
6. R1、R2、S 输出局部效率比分布；至少保存 count、min、Q1、median、Q3、P90、max。
7. 现有 fine/mixed/macro 估计子移入 diagnostics 命名空间，避免与论文 \(\eta_H\) 混名。

建议产物：

- include/helmholtz/adaptive/kernel_residual.h
- src/helmholtz/adaptive/kernel_residual.cpp
- tests/test_helmholtz_kernel_residual.cpp

验收：

- \(I_H\xi_z\) 的相对约束残差不超过 \(10^{-10}\)；
- 局部 Riesz 方程残差和能量恒等式误差不超过 \(10^{-10}\)；
- kernel-basis 与 saddle-point 的 \(\eta_{H,z}\) 在小网格上一致；
- \(\sum_T\eta_{H,T}^2=\sum_z\eta_{H,z}^2\) 的相对误差不超过 \(10^{-12}\)；
- 显式构造全局 kernel 的小网格上验证论文的上下界方向；
- 当 \(T_H=T_h=T_{\widehat h}\) 时，kernel residual 和 \(u_{\widehat h}-U\) 同时接近机器精度；strong residual 不参加此门槛。

### WP4：Corrector、谱与稳定性证书

**状态：未实现；CALOD 的核心阻塞项。**

任务：

1. 在 audit kernel 空间构造 \(\zeta_z^{\mathrm{tot}}(\phi_i)\)，装配 \(G_{\mathrm{tot}}\)。
2. 在每个 \(D_{T,\ell}\) 构造 fine-corrector residual Riesz 元 \(z_{T,h}(\phi_i)\)，装配 \(G_{h,T}\) 和 \(G_h\)。
3. 装配 coarse \(M_\kappa\)，实现广义最大特征值和主谱簇检测。
4. 实现 \(\Theta_{\mathrm{tot}},\Theta_h\) 及全部 \(\delta\)、\(q\)、stability margin、\(L_{\mathrm{LOD}},U_{\mathrm{LOD}}\)。
5. 计算 \(\gamma_{V_{\widehat h}}(\kappa)\) 的最小奇异值近似和验证下包络。
6. 增加 LOD_ENABLE_VERIFIED_CERTIFICATES 构建选项，生产后端固定为 MPFR/MPFI 定向舍入区间算术；普通 Eigen 路径只输出 diagnostic approximation。复杂 Hermitian 问题通过实 \(2\times2\) block 表示进入 verified dense kernel。
7. 用区间 Cholesky/残差包络验证 \(M_\kappa\) 正定和广义最大特征值上界；用验证逆残差或等价区间 SVD 给出 \(\gamma_{V_{\widehat h}}\) 下界；局部线性求解使用 interval residual correction 包围 Riesz 解误差。
8. 将线性求解误差、特征对残差和舍入误差显式带入 enclosure，并把 backend、precision、rounding mode 写入 metadata。
9. 实现常数注册表，逐项记录

   \[
   C_{\mathrm{app}},C_{\mathrm{st}},C_{\mathrm{sd}},
   C_{\mathrm{ov}},C_a,C_{\mathrm{Fort}},c_W,C_\Pi,
   C_{\mathrm{loc}},C_{\mathrm{ol}}(\cdot),\beta,s
   \]

   的值、界方向、来源、推导/导入脚本、适用 mesh class、patch-policy hash 和 verified 标志。\(C_{\mathrm{app}},C_{\mathrm{st}},C_{\mathrm{sd}},C_{\mathrm{ov}}\) 由有限 NVB patch type 的解析组合界或 interval generalized eigenproblem 验证；\(c_W,C_\Pi,C_a,C_{\mathrm{Fort}}\) 按论文公式由已验证输入传播；\(C_{\mathrm{loc}},\beta,s\) 若没有严格值，则 localization 上界只使用 \(\widehat\delta_{\mathrm{tot}}+\widehat\delta_h\) 分支。
10. 在正式 real-data/mixed-boundary 离散上验证 conjugation invariance 和 primal/adjoint certificate 等价；失败时装配两侧证书并取保证稳定性的较差界。
11. 根据 dominant eigenvector 或主谱簇形成 \(\eta_{h,T}\)。

建议产物：

- include/helmholtz/adaptive/certificates.h
- src/helmholtz/adaptive/certificates.cpp
- include/solver/verified_spectrum.h
- src/solver/verified_spectrum.cpp
- tests/test_helmholtz_certificates.cpp
- tests/test_verified_spectrum.cpp

验收：

- \(G_{\mathrm{tot}},G_h\) 的 Hermitian 误差低于 \(10^{-12}\)，最小特征值在允许的 enclosure 内非负；
- 稠密小问题、显式算子范数和广义特征值结果一致；
- 人为降低 corrector fine resolution 时主要增大 \(\Theta_h\)；
- 人为减小 \(\ell\) 时 \(\underline\delta_\ell\) 或 localization 上界能暴露缺陷；
- 验证上界始终不小于高精度参考值，验证下界始终不大于高精度参考值；
- MPFR/MPFI backend 在不同 precision 下给出嵌套 enclosure，Eigen-only 构建永远不能输出 certified；
- R1、R2、S 的 primal/adjoint 共轭门槛通过，或两侧独立证书路径通过；
- 缺任一 verified 输入时状态自动降为 conditional，不能由命令行强行标成 certified。

### WP5：论文版 CALOD 与 HLOD 状态机

**状态：H-only proxy 已有；完整状态机未实现。**

任务：

1. 将当前 driver 拆成可测试的状态：COARSE_ADMISSIBILITY、CORRECTOR_CERTIFICATION、COARSE_ERROR_CONTROL、AUDIT_CONTROL、DONE、WORK_LIMIT、FAILURE。
2. 实现违反 \(\mu_H\) 的强制 coarse refinement。
3. 实现 \(\widehat\delta_h>\tau\underline\delta_{\mathrm{tot}}\) 时的 corrector-\(h\) 分支。
4. 实现 reverse branch 的全局 \(\ell\leftarrow\ell+1\)。
5. corrector 合格后才允许 \(\eta_H\) 驱动 \(H\)-refinement。
6. 实现 audit refinement 和 continuous interval 停止。
7. HLOD 复用同一 \(\eta_H\) 和 coarse 标记代码，但严格冻结第 2.2 节的 \(h_{\mathrm{prior}},\ell_{\mathrm{prior}}\)。
8. HLOD-proxy 保留旧 driver，输出命名空间和状态与 HLOD 分开。

建议产物：

- include/helmholtz/adaptive/certified_driver.h
- src/helmholtz/adaptive/certified_driver.cpp
- tests/test_helmholtz_certified_driver.cpp

验收：

- 单元测试可分别强制触发 \(H\)、\(h\)、\(\ell\)、audit 四个分支；
- 状态转移顺序与论文算法逐步一致；
- exact/reference error 无法从 driver decision context 访问；
- 每个停止原因唯一、结构化且可序列化；
- 给定固定标记序列，checkpoint 恢复与连续运行逐轮一致；
- CALOD 在无法满足稳定性或资源限制时明确失败，不能静默转为 H-only。

### WP6：统一 comparator 和 matched-tolerance driver

**状态：standard LOD 和 uniform FEM 核心存在；统一比较与 AFEM 缺失。**

任务：

1. 用同一 ProblemDefinition 封装 CALOD、HLOD、SLOD-prior、SLOD-matched、UFEM、AFEM。
2. HLOD 与 SLOD-prior 共同读取 baseline_parameters.csv 中按第 2.2 节规则生成的 \(h_{\mathrm{prior}},\ell_{\mathrm{prior}}\)，不得在看正式结果后修改。
3. SLOD-matched 等待同一 case、\(\kappa\) 的 CALOD 成功或资源停止后，读取全部有效历史的最大 local fine level 和最大 \(\ell\)；CALOD censored 不取消该对照，只附加 censored-source 标志；无任何有效参数时才输出 unavailable。
4. 实现 conforming adaptive \(P_1\) FEM，固定局部指标

   \[
   \eta_{\mathrm{AFEM},T}^2
   =h_T^2\|f+\Delta u_T+\kappa^2u_T\|_T^2
   +\frac12\sum_{e\subset\partial T\cap\mathcal E_{\mathrm{int}}}
      h_e\|[\nabla u_h\cdot\nu_e]\|_e^2
   +\sum_{e\subset\partial T\cap\Gamma_R}
      h_e\|\nabla u_h\cdot\nu-\mathrm i\kappa u_h-g_R\|_e^2.
   \]

   Dirichlet 边没有边界 residual；公共内边贡献向两侧各分一半；用同一 \(\theta_H\) 选择最小或近似最小 Dörfler 集并做 NVB 闭包。
5. 所有方法每轮都计算统一 reference error，但该值只送往 evaluator。
6. matched-target evaluator 选择第一次达到每个目标的迭代，并保留未达目标状态。

建议产物：

- benchmarks/bench_helmholtz_paper.cpp
- include/helmholtz/experiments/method_runner.h
- src/helmholtz/experiments/method_runner.cpp
- src/helmholtz/adaptive/fem_driver.cpp
- tests/test_helmholtz_adaptive_fem.cpp
- tests/test_helmholtz_matched_targets.cpp

验收：

- 六个正式配置输出同一 schema；
- 相同方法经旧独立 benchmark 和新统一 runner 得到一致解；
- AFEM 的 element、interior-edge 和 Robin contribution 分别与独立高阶求积/手工小网格参考一致，公共边唯一计数且 Dirichlet 边不含 Robin 项；
- 同一误差目标的选择结果不依赖输出顺序；
- reference 时间从 method time 中剔除；
- SLOD-matched 的 \(h,\ell\) 与对应 CALOD 历史严格匹配。

### WP7：工作量、内存、复用和结果写出

**状态：有局部阶段 wall time；论文级累计量缺失。**

任务：

1. 每轮记录 mesh、corrector assembly/solve、certificate、global solve、estimate、mark、refine、I/O 的 wall 和 core time。
2. 记录 cumulative wall time、cumulative core time、peak RSS、global unknowns、最大局部 saddle unknowns 和 peak active unknowns。
3. 对每个 rebuilt corrector 记录

   \[
   \dim W_{h_n}(D_{T,\ell_n}),
   \]

   并累计

   \[
   W_{\mathrm{patch}}^{\mathrm{cum}}
   =\sum_n\sum_{T\in\mathcal R_n}
   \dim W_{h_n}(D_{T,\ell_n}).
   \]

4. 每个 corrector 保存 source element、patch、\(h\)、\(\ell\)、primal/adjoint、mesh/interpolation version 和 rebuilt/reused 标志。
5. **论文必需路径**：提供 full_rebuild，完整记录 rebuilt corrector 和 patch work；该路径可以直接用于正式实验。
6. **可选性能路径**：在不推迟正式矩阵的前提下实现 incremental，并逐轮与 full_rebuild 比较网格、basis、解、估计子、证书和停止决策。
7. 保存 iteration CSV、local indicator CSV、metadata JSON、timing JSON、mesh/field VTK、stdout/stderr 和 checksum manifest。

验收：

- 累计量严格等于各轮量求和；
- full rebuild 时 \(W_{\mathrm{patch}}^{\mathrm{cum}}\) 可由历史独立复算；
- 若启用 incremental，其解和证书必须与 full rebuild 在预定容差内一致；未实现 incremental 不阻塞 G5–G7；
- peak memory 的采样方式在全部方法中一致；
- 中断后已有运行目录保持可解析，状态为 interrupted 而不是伪成功。

### WP8：自动化测试、smoke matrix 和回归

**状态：现有 adaptive/reliability 测试只覆盖 Stage-1。**

测试层级：

1. **代数层**：复数共轭约定、残差恒等式、kernel constraint、Riesz energy identity、Gram matrix、谱 enclosure。
2. **网格与边界层**：三层嵌套、NVB 标签继承、L 型面积、mixed boundary 测度和 right inverse。
3. **证书层**：小网格直接参考、上下界方向、状态降级、局部量求和。
4. **算法层**：四分支触发、Dörfler 最小性、停止和 work-limit、checkpoint。
5. **方法层**：R1/R2/S 上六配置小规模 smoke，输出 schema 完整。
6. **回归层**：现有 elliptic LOD、Helmholtz FEM、standard LOD、corrector 和 Stage-1 测试全部继续通过。

smoke 运行可以使用较小波数或较浅网格，但必须标记 smoke，不能混入论文原始数据。

### WP9：生产实验、聚合和论文图表

**状态：必须在 WP0–WP8 的论文必需门槛全部通过后启动；WP7 的 incremental 子项不在阻塞路径。**

任务：

1. 先运行 R1，校验 exact error、audit error、证书、无虚假网格集中和全部方法接口。
2. 再运行 R2a/R2b，检查局部源附近的 coarse concentration 和 matched-target 收益。
3. 再运行 S，检查凹角网格集中、\(2/3\) 奇性和 mixed-boundary 全链路。
4. 最后聚合 K，只有误差、稳定性和工作量同时支持时才讨论波数鲁棒性或污染抑制。
5. 正式计时重复 5 次；数值解和自适应历史只要确定性一致，聚合采用时间中位数。
6. 一键生成论文表、图片、figure-data CSV 和 LaTeX 可直接读取的表格片段。
7. 同步更新论文数值章节的误差记号：cert-audit 区间/effectivity 使用 \(u_{\widehat h}^{\mathrm{cert}}\)，六方法 matched-target 使用公共 \(u_{\mathrm{ref}}\)。

验收：

- 所有正式配置、目标和重复编号都有成功或明确 censored 状态；
- 聚合脚本不读取手工编辑后的数据；
- 图中每个点能追溯到 run ID；
- 论文主表和图由脚本重新生成后无 diff；
- 原始结果、聚合结果与论文文字中的方法名称完全一致。

## 6. 正式实验矩阵与执行顺序

### 6.1 主矩阵

主矩阵统一使用 \(\theta_H=0.5\)：

| Case | 数据参数 | \(\kappa\) | 方法配置 | 目标 |
|---|---|---|---|---|
| R1 | 精确光滑解 | 8, 16, 32 | CALOD、HLOD、SLOD-prior、SLOD-matched、UFEM、AFEM | 四个共同误差目标 |
| R2a | \(\sigma=2^{-5}\) | 8, 16, 32 | 同上 | 同上 |
| R2b | \(\sigma=2^{-6}\) | 8, 16, 32 | 同上 | 同上 |
| S | \(r^{2/3}\) 奇异解 | 8, 16, 32 | 同上 | 同上 |

SLOD-matched 不是独立预设参数；它必须在对应 CALOD 历史产生后运行。K 是 R2a、R2b、S 的聚合视图，不额外改变 PDE。

### 6.2 敏感性矩阵

对 R1、R2a、R2b、S 的每个正式波数 \(\kappa=8,16,32\) 使用

\[
\theta_H\in\{0.3,0.5,0.7\}.
\]

主文至少报告 CALOD、HLOD 和 AFEM 的敏感性；standard LOD 和 UFEM 没有 Dörfler 参数，不重复运行。若 \(\theta_h\) 也做敏感性，只进入补充材料，不与 \(\theta_H\) 混为一张图。

### 6.3 执行阶段

1. **开发 smoke**：小网格验证全部模块和所有适用的结果状态，不产生论文结论。
2. **参数冻结 preflight**：估算每个方法的内存和时间，冻结 work limits、initial levels、solver tolerances、verified constants 和 standard LOD 的 \(c\)。
3. **R1 calibration**：必须首先通过；任何证书不包含 audit error、误差计算不一致或网格虚假集中都阻止后续运行。
4. **R2 regular adaptivity**：完成两种 \(\sigma\)、三种 \(\kappa\) 和主矩阵。
5. **S singular adaptivity**：完成 mixed boundary、精确奇性和全部 comparator。
6. **K aggregation**：汇总三种波数，检查稳定裕量和工作量。
7. **正式重复计时**：在参数和代码冻结后重新运行，不复用开发计时。
8. **论文导出**：生成 figure-data、表格、图片、日志索引和数据说明。

### 6.4 R1 专项验收

- exact source 和 homogeneous impedance 条件再次自动验证；
- fine/audit FEM 误差随加细下降；
- CALOD 区间包含 audit error；
- \(\eta_H\) 的全局上下界方向正确；
- coarse mesh 不应无理由集中在单个内部点；报告 level histogram、最大/最小粗层差和网格熵类诊断；
- strong-residual proxy 与论文 \(\eta_H\) 分栏报告，不混称 effectivity。

### 6.5 R2 专项验收

- 每个 \(\sigma\) 的 \(L^2\) 归一化通过；
- audit/reference 独立于所有 adaptive marked sets；
- coarse mesh 在源附近和波场显著区域产生可解释集中；
- \(\sigma=2^{-6}\) 不得因初始 \(V_h\) 分辨率不足而被误判为 \(H\) 误差；corrector-\(h\) 和 audit gate 必须先通过；
- 报告两级 audit saturation 比例及其状态标签。

### 6.6 S 专项验收

- L 型域、凹角、角度约定和 \(\Gamma_D/\Gamma_R\) 与论文一致；
- cut-off、源项和边界数据由同一制造解对象生成；
- exact solution 的局部能量误差计算能处理角点奇性；
- uniform FEM/LOD 显示受低正则性影响的效率，adaptive FEM/LOD 的网格在凹角集中；
- 不用全 Robin 平滑源替代 S；平滑 L 型算例只能作为额外诊断。

### 6.7 K 专项验收

对 R2a、R2b、S 分别绘制或列出随 \(\kappa\) 变化的：

- \(E_\kappa^{\mathrm{ref}},E_0^{\mathrm{ref}}\)；
- 首次达标累计 wall time；
- \(W_{\mathrm{patch}}^{\mathrm{cum}}\)；
- peak memory 和 peak unknowns；
- \(\mu_H,q_h,q_\ell,q_{\mathrm{tot}},m_{\mathrm{stab}}\)；
- \(\Theta_{\mathrm{tot}},\Theta_h,\eta_H\)；
- 失败或 censored 比例。

只有这些量共同支持时，论文才可以使用 wavenumber robust 或 pollution suppression 一类表述。

## 7. 数据协议

### 7.1 每轮 iteration.csv 必须包含

**身份字段**

- schema_version、run_id、case_id、method_id、repeat_id；
- \(\kappa,\sigma,\theta_H,\theta_h,\mu_0,q_0,\tau,\rho_{\mathrm{aud}}\)；
- iteration、inner_iteration、state、action、decision_reason、status。

**空间与网格**

- \(N_H,N_h,N_{\widehat h}\)、coarse elements、fine elements、audit elements；
- min/max/local coarse level，local corrector \(h\) level 分布；
- \(H_{\max}\)、当前全局 \(\ell\)；
- marked_H、marked_h、closure_added、audit_refined；
- rebuilt_correctors、reused_correctors。

**误差与 reference**

- cert-audit 真误差 \(e_{\mathrm{cert}}=\|u_{\widehat h}^{\mathrm{cert}}-U\|_\kappa\) 及其局部 \(D_z^+\) 版本；
- 公共比较误差 \(E_\kappa^{\mathrm{ref}},E_0^{\mathrm{ref}}\)；
- exact error（若存在）；
- cert-audit level、cert-audit refinement 和内部 saturation/certificate 状态；
- evaluation-reference level、evaluation saturation ratio 和 reference_status；
- \(L_{\widehat h},U_{\widehat h},L_{\mathrm{true}},U_{\mathrm{true}}\)。

**估计子与证书**

- \(\eta_H,\Theta_{\mathrm{tot}},\Theta_h\)；
- \(\underline\delta_{\mathrm{tot}},\widehat\delta_{\mathrm{tot}}\)；
- \(\widehat\delta_h,\underline\delta_\ell,\widehat\delta_\ell\)；
- \(q_{\mathrm{tot}},q_h,q_\ell,\mu_H\)；
- verified inf-sup lower bound、stability margin；
- \(L_{\mathrm{LOD}},U_{\mathrm{LOD}}\)；
- upper effectivity \(\mathcal I_U=U_{\mathrm{LOD}}/e_{\mathrm{cert}}\)、certificate gap \(\mathcal G\)；
- algebraic、Petrov、corrector、constraint residual；
- certificate_status 和所有 enclosure residual。

对 \(\mathcal I_U,\mathcal G,\mathcal I_{\mathrm{loc},z}\) 及任何比值，只有分母严格为正时才写数值。分母非正时值写 JSON/CSV typed null，并同时写 value_status=invalid_denominator；非 LOD 方法的证书字段写 typed null 和 value_status=not_applicable。禁止用 0、infinity 或字符串 NaN 代替。

**工作量**

- mesh、operator、corrector、certificate、global solve、estimate、mark、refine、I/O wall time；
- 本轮和累计 wall/core time；
- offline/online time；
- evaluation_reference_time，且明确 excluded_from_method_time；cert-audit time 包含在 CALOD certificate/method time；
- 本轮 processed patch dimension 和 \(W_{\mathrm{patch}}^{\mathrm{cum}}\)；
- peak RSS、global unknowns、max local unknowns、peak active unknowns。

### 7.2 局部数据

local_H.csv 每个节点/单元至少记录：

- coarse node ID、\(D_z,D_z^+\) 标识；
- \(\eta_{H,z}^2,\eta_{H,T}^2\)；
- local error、local effectivity；
- kernel constraint residual、Riesz residual；
- 是否被 \(H\)-Dörfler 标记。

local_h.csv 每个 coarse-source corrector 至少记录：

- stable element ID、patch ID、patch element/DOF 数；
- current \(h\) level、\(\ell\)；
- \(\eta_{h,T}^2\) 或谱簇聚合贡献；
- rebuilt/reused、cache signature；
- 是否被 \(h\)-Dörfler 标记。

### 7.3 metadata.json 必须包含

- schema version、完整 resolved configuration 和 hash；
- Git commit、dirty flag 和变更文件列表；
- compiler、build type、编译选项、依赖版本；
- CPU、内存、操作系统、线程数、NUMA/affinity；
- quadrature rules、linear/eigen/SVD tolerances；
- constants registry 及 provenance；
- 每个可空字段的 value_status 枚举定义：valid、not_applicable、not_computed、invalid_denominator、enclosure_failed；
- random seed；若算法完全确定，也显式记录 deterministic=true；
- start/end time、exit code、stop reason；
- raw 文件 checksum。

### 7.4 推荐结果目录

results/helmholtz_adaptive_paper/

- schema/
- configs/
- raw/R1/、raw/R2a/、raw/R2b/、raw/S/
- aggregate/
- figure_data/
- tables/
- logs/
- manifest/

每个 raw run 目录只写一次，不允许后处理脚本原地改写。聚合输出写入 aggregate，论文图片写入 figures/paper，图表数据同时复制到 figure_data。

## 8. 论文必须生成的表格和图

### 8.1 表格

1. **Benchmark matrix**：R1、R2a、R2b、S、K 的区域、正则性、参数、运行状态。
2. **Matched-target accuracy/work table**：case、method、target、\(N_H\)、peak unknowns、\(E_\kappa^{\mathrm{ref}}\)、累计时间、内存、状态。
3. **Certificate table**：\(\eta_H,\Theta_{\mathrm{tot}},\Theta_h,L_{\mathrm{LOD}},U_{\mathrm{LOD}},\mathcal I_U,\mathcal G,m_{\mathrm{stab}}\)。
4. **Wavenumber table**：R2/S 的 \(\kappa\)、误差、时间、patch work、稳定裕量。
5. **Reproducibility table**：硬件、线程、编译器、容差、重复次数和 commit。

SLOD-prior 与 SLOD-matched 在最终表中必须分成两行，不能继续合并为 standard LOD 一行。

### 8.2 图片

1. R1、R2a、R2b、S 的初始/最终 coarse mesh、marked sets 和 source/solution。
2. \(E_\kappa^{\mathrm{ref}}\) 对 coarse DOF、peak unknowns、累计 wall time。
3. 达到共同目标时各方法累计时间和内存的比较。
4. \(\eta_H,\Theta_{\mathrm{tot}},\Theta_h\) 与每轮 refinement action。
5. \(\mu_H,q_h,q_\ell,q_{\mathrm{tot}},m_{\mathrm{stab}}\) 的历史。
6. upper effectivity 和 certificate gap 的历史。
7. local effectivity 的 box/violin 或分位数图。
8. corrector、certificate、global solve、estimate 等累计时间分解。
9. \(W_{\mathrm{patch}}^{\mathrm{cum}}\) 对误差。
10. \(\kappa\) 对误差、时间、memory、patch work 和稳定裕量。
11. audit level、saturation ratio 和 reference error floor。
12. \(\theta_H=0.3,0.5,0.7\) 的敏感性图。

所有图只读取 aggregate/figure_data，不直接解析终端日志。

## 9. 分阶段门槛与依赖顺序

### G0：现有基线不回退

- 当前 mesh、NVB、QI、elliptic LOD、Helmholtz FEM、corrector、standard LOD 和 Stage-1 adaptive 测试全部通过。
- 记录现有 R1、H-convergence 和 \(k\)-scan 的小规模 golden 结果。

### G1：论文问题定义可用

依赖 WP0、WP1。

**状态：已通过（2026-08-08）；R1/R2a/R2b/S 数据、mixed-boundary 全链路、L 型几何、Gaussian 归一化以及 R2/S 载荷与误差积分升阶稳定性均由回归测试覆盖。**

- R1、R2a、R2b、S 的问题对象完成；
- mixed boundary 全链路测试通过；
- L 型网格和 Gaussian 归一化通过。

未过 G1 不得开始 S 的任何性能实验。

### G2：三层空间和 \(\eta_H\) 可信

依赖 WP2、WP3。

- 三层嵌套和 right inverse 通过；
- 局部 kernel Riesz 通过 basis/saddle 双实现校验；
- 局部量求和、Dörfler 和局部效率诊断通过。

未过 G2，所有 adaptive 结果只能使用 HLOD-proxy 标签。

### G3：证书链可运行

依赖 WP4。

- \(\Theta_{\mathrm{tot}},\Theta_h,\delta,q,L/U\) 全链路通过；
- 谱和 inf-sup enclosure 通过小问题验证；
- conditional/audit-certified 状态降级逻辑通过。

未过 G3，不得生成 CALOD 论文数据。

### G4：四步状态机可控

依赖 WP5。

- \(H,h,\ell,\widehat h\) 四个分支均有确定性测试；
- work limit、checkpoint 和 full-rebuild 历史一致性通过；若启用 incremental，再额外通过复用一致性；
- exact/reference 信息隔离通过。

### G5：六配置公平比较

依赖 WP6、WP7。

- AFEM、两种 SLOD 和 matched-target driver 通过；
- 所有方法输出统一 schema；
- 累计时间、reference 排除和 patch work 守恒。

### G6：缩小版整矩阵

依赖 WP8。

- R1、R2a、R2b、S × 六配置至少各有一个 smoke；
- 所有图表脚本能从 smoke 数据完整生成；
- 不允许字符串 NaN、缺列或无解释缺失值；typed null 必须带合法 value_status/reason。

### G7：正式生产与论文交付

依赖 WP9。

- 参数、代码、schema 和绘图脚本冻结；
- 完成全部主矩阵和敏感性矩阵；
- 正式重复计时、聚合和论文导出完成。

关键路径为

\[
\text{协议}
\longrightarrow
\text{mixed boundary 与 cases}
\longrightarrow
\text{三层 audit}
\longrightarrow
\eta_H
\longrightarrow
\text{corrector/spectral certificates}
\longrightarrow
\text{四步 CALOD}
\longrightarrow
\text{comparators}
\longrightarrow
\text{正式矩阵}.
\]

强残差调参、增量性能优化和额外物理算例不能插入这条路径之前。

## 10. 参数登记表

下表区分论文已经固定的值与必须在正式运行前冻结的工程值。

| 参数 | smoke 建议值 | 正式值/冻结规则 | 状态 |
|---|---:|---|---|
| \(\kappa\) | 4 或 8 | 8, 16, 32 | 论文固定 |
| \(\sigma\) | \(2^{-5}\) | \(2^{-5},2^{-6}\) | 论文固定 |
| \(\theta_H\) | 0.5 | 主实验 0.5；敏感性 0.3、0.7 | 论文固定 |
| target | \(10^{-1}\) | \(10^{-1},5\cdot10^{-2},2\cdot10^{-2},10^{-2}\) | 论文固定 |
| \(\rho_{\mathrm{aud}}\) | 0.05 | 0.05，若论文修改则全局同步 | 暂定固定 |
| \(c_h\) | \(1/8\) | 第 2.2 节 prior-fine 规则固定为 \(1/8\) | 本计划固定 |
| standard LOD \(c\) | 1 | \(\ell_{\mathrm{prior}}=\lceil\log_2\kappa\rceil\) | 本计划固定 |
| \(\theta_h\) | 0.5 | preflight 前冻结，默认 0.5 | 待冻结 |
| \(\mu_0\) | 0.5 | 由可验证常数与初始网格规则冻结 | 待冻结 |
| \(q_0\) | 0.25 | 与稳定性条件一起冻结 | 待冻结 |
| \(\tau\) | 0.5 | h/ell 人工缺陷测试后、正式结果前冻结 | 待冻结 |
| \(\ell_0\) | \(\lceil\log_2\kappa\rceil\) | 配置固定，不按正式结果回调 | 待冻结 |
| \(H_0,h_0,\widehat h_0\) | 小层级 | 对每个 \(\kappa\) 由 admissibility、corrector 和 audit preflight 冻结 | 待冻结 |
| initial coarse mesh | case 默认网格 | CALOD/HLOD/AFEM 共享同一 mesh hash | 规则固定 |
| work/memory/time limits | 较小 | 每个 case/\(\kappa\) 对全部方法统一，preflight 后冻结 | 待冻结 |
| linear relative residual | \(10^{-10}\) | 所有方法相同；certified 模式另计 enclosure | 待冻结 |
| eigen/SVD residual | \(10^{-10}\) | certified 模式使用验证包络 | 待冻结 |
| spectral cluster gap | \(10^{-6}\) 相对值 | 小问题和扰动测试后冻结 | 待冻结 |
| support propagation \(p_I\) | 由插值实现测得并验证 | constants registry 中冻结 | 待验证 |
| node patch \(m_D\) | \(p_I+1\) | \(D_z=\omega_z^{m_D}\) | 规则固定 |
| \(D_z^+\) | 论文 union 公式 | 绑定 \(\omega_T^\dagger\) 和 patch-policy hash | 规则固定 |
| verified backend | Eigen diagnostic | MPFR/MPFI directed rounding | 本计划固定 |
| timing repeats | 1 | 5，报告中位数 | 本计划固定 |
| max iterations | 较小 | 可按方法结构设置，但必须由统一 work/time limit 截断且正式结果前冻结 | 待冻结 |
| \(r_0,r_1\) | \(1/4,1/2\) | \(1/4,1/2\)，除非制造解测试失败 | 暂定固定 |

所有“待冻结”项必须在 production manifest 中有确定值和冻结日期。查看正式结果后不得只对不利方法回调参数。

## 11. 风险与处理

### 11.1 Verified constants 暂时不可得

风险：算法公式实现完成，但 \(C_{\mathrm{sd}},C_{\mathrm{ov}}\) 或谱量没有严格 enclosure。

处理：

- 先完成 conditional 模式和全部数值恒等式；
- 状态字段强制为 conditional；
- 继续实现 verified enclosure；
- 论文若在提交前仍缺严格输入，必须把数值措辞从 fully certified 降级，不能用近似值伪装验证值。

### 11.2 Audit space 成本过高

处理：

- cert-audit factorization 与 RHS 分离并缓存，但其成本仍计入 CALOD；
- evaluation-reference time 单列排除；
- 小规模运行验证证书，大规模运行可使用矩阵-free Riesz，但必须保留 enclosure；
- 达到 memory limit 时输出 censored，不静默换成较粗 reference。

### 11.3 局部 \(h\) 更新破坏嵌套或大量重建

处理：

- 第一版只允许一个全局共形、局部加细的 \(V_h\)；
- 选中 patch 的 fine elements 取并集后一次闭包；
- 每轮验证 \(V_H\subset V_h\subset V_{\widehat h}\)；
- 先用 full rebuild 确认正确，再做 signature-based reuse。

### 11.4 主特征值成簇导致标记不稳定

处理：

- 检测相对 spectral gap；
- 对主谱簇的 \(M_\kappa\)-正交基聚合 \(G_{h,T}\) 贡献；
- 在小扰动下测试 marked set 稳定性。

### 11.5 S 的奇异源求值不稳定

处理：

- 分离 \(r<r_0\)、过渡环和 \(r>r_1\) 的解析表达；
- 不在 \(r=0\) 直接使用含负幂的浮点表达；
- 用符号/高精度离线值和多阶 quadrature 交叉验证；
- 保留 exact solution point evaluation 与 weak-form manufactured load 两条测试。

### 11.6 性能优化改变数值结果

处理：

- 论文正确性基线始终保留 full_rebuild；
- incremental 每次改动都运行逐轮 equivalence；
- 正式方法必须统一说明是否启用复用；
- 不允许某一 comparator 独享未计时的预计算。

## 12. 提交论文前完成清单

### 理论定义与实现一致性

- [ ] 使用 \(V_H\subset V_h\subset V_{\widehat h}\)。
- [ ] \(\eta_{H,z}\) 来自 audit-kernel 受约束 Riesz，不是 strong residual。
- [ ] 节点到单元分配无重复且严格守恒。
- [ ] \(\Theta_{\mathrm{tot}},\Theta_h,G_{\mathrm{tot}},G_h,M_\kappa\) 与论文定义一致。
- [ ] \(q_{\mathrm{tot}},q_h,q_\ell\) 和所有 \(\delta\) 界方向正确。
- [ ] \(\ell\) 只做全局更新，没有把局部 \(\ell_T\) 混入论文算法。
- [ ] certified 标签只在全部验证输入存在时输出。

### 算例

- [ ] R1：\(\kappa=8,16,32\)。
- [ ] R2a：\(\sigma=2^{-5}\)，\(\kappa=8,16,32\)。
- [ ] R2b：\(\sigma=2^{-6}\)，\(\kappa=8,16,32\)。
- [ ] S：mixed boundary、\(r^{2/3}\) 奇异解，\(\kappa=8,16,32\)。
- [ ] K：R2a、R2b、S 的波数聚合。
- [ ] \(\theta_H=0.3,0.5,0.7\) 敏感性完成。

### 方法

- [ ] CALOD。
- [ ] fixed-\(h,\ell\) HLOD。
- [ ] SLOD-prior。
- [ ] SLOD-matched。
- [ ] UFEM。
- [ ] AFEM。
- [ ] HLOD-proxy 只作为诊断，不混入 CALOD。

### 误差、证书和工作量

- [ ] \(e_{\mathrm{cert}},E_\kappa^{\mathrm{ref}},E_0^{\mathrm{ref}}\)。
- [ ] \(\mathcal I_U,\mathcal G,m_{\mathrm{stab}}\)。
- [ ] \(L_{\mathrm{LOD}},U_{\mathrm{LOD}},L_{\mathrm{true}},U_{\mathrm{true}}\) 及 validity status。
- [ ] local effectivity 分布。
- [ ] \(\mu_H,q_h,q_\ell,q_{\mathrm{tot}}\) 历史。
- [ ] local \(h\) levels、全局 \(\ell\)、marked sets。
- [ ] cumulative wall/core time。
- [ ] evaluation-reference time 单列；cert-audit time 计入 CALOD。
- [ ] \(W_{\mathrm{patch}}^{\mathrm{cum}}\)。
- [ ] peak memory、global/max-local/peak-active unknowns。
- [ ] offline/online 分解。

### 公平性与可复现性

- [ ] 同一 PDE、积分、容差、硬件和线程。
- [ ] 每个目标取第一次达标迭代。
- [ ] 未达目标记录 censored/work-limit。
- [ ] 正式时间重复 5 次并报告中位数。
- [ ] run config、commit、dirty flag、编译器和硬件完整。
- [ ] raw 数据只读，aggregate 可重建。
- [ ] 每张图和每行表可追溯到 run ID。
- [ ] 一条命令可从 aggregate 重新生成论文图表。
- [ ] 论文已分开 \(u_{\widehat h}^{\mathrm{cert}}\) 与 \(u_{\mathrm{ref}}\) 的记号和公式。

## 13. 非论文阻塞的保留研究

下列已有或拟议内容可以继续保存在 diagnostics、补充材料或后续计划中，但只能在 G7 关键路径不受影响时推进：

- fine/mixed/macro strong-residual 指标比较；
- 经验 reliability envelope 和留出拟合；
- 高对比 \(A,n,\beta\)；
- 局部 \(\ell_T\) 和 corrector-tail 指标；
- 每个 patch 的独立细网格；
- cost-score 联合 \(H/h/\ell\) 决策；
- PML、复杂系数和非共轭不变问题；
- 自适应收缩、准最优复杂度与更强污染结论。

项目接下来的唯一执行顺序是：先完成论文实验合同，再讨论这些扩展。
