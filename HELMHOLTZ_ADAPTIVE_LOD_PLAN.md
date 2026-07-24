# Helmholtz 自适应 LOD 数值实验计划书

> 状态日期：2026-07-25
>
> 当前阶段：第一阶段工程基线和制造解校准已完成；正在加入核空间受限的局部代数对偶估计子，生产估计子尚未冻结。
>
> 暂停点：先实现并验证核空间代数对偶估计子及伴随 corrector 证书；在此之前第二、三阶段不启动。

## 执行状态

| 工作包 | 状态 | 结论 |
|---|---|---|
| A1 固定细空间、稳定 NVB 树和嵌套检查 | 已完成 | full-rebuild 金标准可用 |
| A2 broken residual、跳量、Robin 项和代数残差恒等式 | 已完成 | 装配正确性已验证 |
| A3 `fine/mixed/macro` 三类聚合指标 | 已完成 | 作为候选保留 |
| A4a 无约束、分割残差的局部能量 Riesz 对偶校准 | 已完成 | 只作离散诊断基线 |
| A4b 核空间受限的局部代数对偶估计子 | 未开始 | 新的理论主候选，需验证稳定分解与 corrector 证书 |
| A5 可靠性包络、拟合和留出工具 | 已完成 | 尚缺足够参数数据 |
| A6 制造解连续误差校准 | 已完成 | 区分了 `u-u_LOD` 与 `u_h-u_LOD` |
| A7 冻结唯一生产指标 | 待重新评估 | A4b 完成后与强残差候选统一比较 |
| A8 增量 corrector 复用 | 未开始 | 先保持 full rebuild 金标准 |
| B L 型区域低正则实验 | 未开始 | 等待 A7 |
| C `H/h/ell` 联合自适应 | 未开始 | 等待 A7 和 B |

使用方法见 [HELMHOLTZ_ADAPTIVE_GUIDE.md](HELMHOLTZ_ADAPTIVE_GUIDE.md)，
已测结果见 [DEVELOPMENT.md](DEVELOPMENT.md)。

## 1. 目标与实施原则

本计划遵循“先建立可核验的数值估计子，再讨论理论强化和联合自适应”的顺序：

1. 在单位正方形上固定细网格 `T_h` 和局部化层数 `ell`，只自适应粗网格 `T_H`。
2. 用细网格解 `u_h` 作为全离散金标准，检验残差估计子的可靠性、有效性和标记能力。
3. 在 L 型区域上比较统一粗网格 LOD 与自适应粗网格 LOD，观察低正则性下的误差-自由度和误差-时间关系。
4. 在上述基线可信后，研究分别指示 `H`、`h` 和 `ell` 的三部分估计子。
5. 每一阶段都同时检查数值正确性、嵌套性、内存、corrector 时间和增量复用收益。

第一阶段不追求立即证明波数一致的可靠性常数。Helmholtz 问题的后验可靠性会包含离散稳定性或 inf-sup 常数，因此先固定若干波数完成全离散验证，再单独研究关于 `k` 的鲁棒加权。

## 2. 数学模型与离散空间

设 `Omega` 为单位正方形，后续替换为 L 型区域。考虑

```math
-\nabla\cdot(A\nabla u)-k^2 n u=f \quad\text{in }\Omega,
```

以及齐次阻抗 Robin 边界条件

```math
A\nabla u\cdot\nu-i k\beta u=0 \quad\text{on }\Gamma=\partial\Omega.
```

复值弱形式采用“第一变量线性、第二变量共轭线性”的约定：

```math
a_k(u,v)
=(A\nabla u,\nabla v)_\Omega
-k^2(nu,v)_\Omega
-ik(\beta u,v)_\Gamma,
\qquad
F(v)=(f,v)_\Omega.
```

定义用于误差度量和局部 Riesz 问题的正定能量范数

```math
\|v\|_E^2
=(A\nabla v,\nabla v)_\Omega
+k^2(nv,v)_\Omega.
```

令 `V_h` 为固定细网格 `T_h` 上的复值连续 P1 空间，`V_H` 为自适应粗网格 `T_H` 上的复值连续 P1 空间，并始终要求

```math
V_H\subset V_h.
```

拟插值保持当前形式

```math
I_H=E_H\Pi_H^{dg},
\qquad W_h=\ker(I_H|_{V_h}).
```

对粗单元 `T` 及其 `ell_T` 层 patch `omega_T^{ell_T}`，原 corrector 满足

```math
a_{k,\omega_T^{\ell_T}}(Q_{T,h}^{\ell_T}v,w)
=a_{k,T}(v,w)
\quad
\forall w\in W_h(\omega_T^{\ell_T}).
```

人工 patch 边界施加齐次 Dirichlet，patch 与物理边界重合部分保留 Robin 项。全局局部化 corrector 为

```math
Q_{H,h}^{\boldsymbol\ell}=\sum_{T\in\mathcal T_H}Q_{T,h}^{\ell_T}.
```

Petrov-Galerkin LOD 使用原 corrector 形成 trial basis，并用伴随 corrector 形成 test basis。记全离散解为 `u_{H,h}^{boldsymbol ell}`，细网格参考解为 `u_h`。

## 3. 嵌套网格的强制不变量

### 3.1 同一棵 NVB 加细树

粗细网格必须来自同一个带参考边标记的初始三角网格，并使用相同的 newest-vertex bisection 规则。固定的主细网格应包含允许的最大加细深度内整棵 NVB 树的叶子，因此任意自适应粗网格叶子都能表示为细网格单元的并集。

每次粗网格加细后必须验证：

```math
\forall T\in\mathcal T_H,
\qquad
\overline T=\bigcup_{\tau\in\operatorname{children}_h(T)}\overline\tau.
```

并验证粗节点 P1 基函数能由细节点 P1 基函数精确再现。禁止仅通过浮点坐标查找来推断父子关系。

### 3.2 深度约束

若固定主细网格的最大 NVB 深度为 `L_h`，则任何粗叶单元必须满足

```math
L_H(T)<L_h.
```

达到细网格深度的单元不可继续作为粗网格标记对象。若标记集合包含这类单元，应停止实验或先扩充全局细网格，不能产生不满足 `V_H subset V_h` 的网格。

### 3.3 加细闭包

使用现有 `bisect_newest_vertex(..., marked_elements)` 完成相容闭包，不允许悬挂点。实验同时记录用户标记数、闭包新增数和最终加细数，以区分估计子效果与网格相容代价。

## 4. 第一阶段：固定 `h, ell`，只自适应 `H`

第一阶段先完成“残差恒等式 -> 候选指标校准 -> 可靠性结构检验”，再启动自适应循环。自适应循环始终固定 `h,ell`；只有前置校准矩阵临时改变 `h,ell`，用于辨识各误差项。这样可以避免标记算法掩盖估计子本身的问题。

### 4.1 全离散目标误差与残差

第一阶段估计相对于固定细空间的误差

```math
e_h=u_h-u_{H,h}^{\ell},
```

而不是连续误差 `u-u_{H,h}^{ell}`。细空间残差定义为

```math
R_h(v_h)=F(v_h)-a_k(u_{H,h}^{\ell},v_h),
\qquad v_h\in V_h,
```

其矩阵形式为

```math
r_h=b_h-A_hu_{H,h}^{\ell}.
```

由于 LOD 解只在多尺度测试空间上满足 Petrov-Galerkin 正交性，`R_h` 在整个 `V_h` 上通常不为零。所有积分残差必须先与 `r_h` 做逐自由度交叉验证，再用于后验估计。

### 4.2 Broken residual、边集合与唯一归属

`u_{H,h}^{ell}` 只在细单元上是 P1，因此粗单元上的散度必须按细网格 broken 方式理解。对 `tau in T_h` 定义

```math
R_\tau=f+k^2nu_{H,h}^{\ell}
+\nabla\cdot(A\nabla u_{H,h}^{\ell}),
```

对细网格内部边和物理 Robin 边分别定义

```math
J_e=[A\nabla u_{H,h}^{\ell}\cdot\nu_e],
\qquad
B_e=A\nabla u_{H,h}^{\ell}\cdot\nu-ik\beta u_{H,h}^{\ell}.
```

若 `A` 在细单元上为常数，则单个 P1 单元中的散度项为零，但系数跳跃和梯度跳跃仍通过 `J_e` 保留。

对每个粗单元 `T` 建立以下集合：

```math
\mathcal T_h(T)=\{\tau\in\mathcal T_h:\tau\subset T\},
```

以及由细边贡献分割得到的 `E_h^int(T)` 和 `E_h^R(T)`。粗单元内部的细边必须计入；位于两个粗单元公共边上的跳跃贡献向两侧各分一半；Robin 边只属于唯一相邻粗单元。由此保证局部贡献之和严格等于全局贡献，并且与边法向选择无关。

### 4.3 三个残差指标候选

第一轮同时实现三个候选，不预先指定最终标记指标。

**A. 细尺度聚合指标**

```math
\eta_{fine}(T)^2
=\sum_{\tau\in\mathcal T_h(T)}h_\tau^2\|R_\tau\|_\tau^2
+\sum_{e\in\mathcal E_h^{int}(T)}h_e\|J_e\|_e^2
+\sum_{e\in\mathcal E_h^R(T)}h_e\|B_e\|_e^2.
```

它是标准细网格残差的粗父单元聚合，最接近对 `u_h-u_LOD` 的常规全离散估计。

**B. 混合尺度宏观指标**

```math
\eta_{mix}(T)^2
=H_T^2\sum_{\tau\in\mathcal T_h(T)}\|R_\tau\|_\tau^2
+\sum_{e\in\mathcal E_h^{int}(T)}h_e\|J_e\|_e^2
+\sum_{e\in\mathcal E_h^R(T)}h_e\|B_e\|_e^2.
```

这是原 4.0 构想的严格版本：体残差采用粗尺度，边残差保留细尺度。由于解和系数只在细网格上可直接积分，它并不省略细单元遍历；它改变的是权重和粗单元聚合方式。

**C. 完全宏观指标**

```math
\eta_{macro}(T)^2
=H_T^2\sum_{\tau\in\mathcal T_h(T)}\|R_\tau\|_\tau^2
+H_T\sum_{e\in\mathcal E_h^{int}(T)}\|J_e\|_e^2
+H_T\sum_{e\in\mathcal E_h^R(T)}\|B_e\|_e^2.
```

公共粗边从相邻两侧分别使用对应的 `H_T/2` 权重。该指标更接近由 `v-I_Hv` 的粗尺度逼近性质导出的形式，但可能随一个粗单元内细边数量增长而明显高估。

三个全局指标统一记为

```math
\eta_X^2=\sum_{T\in\mathcal T_H}\eta_X(T)^2,
\qquad X\in\{fine,mix,macro\}.
```

### 4.4 两类局部对偶残差：诊断基线与理论主候选

#### 4.4.1 已实现的分割残差 Riesz 校准

当前代码在小中型网格上构造局部 Riesz 金标准。用正定内积

```math
(z,v)_E=(A\nabla z,\nabla v)+k^2(nz,v)
```

在 `omega_T^m` 上求

```math
(z_T,v_h)_{E,\omega_T^m}=R_{h,T}(v_h)
\quad\forall v_h\in V_{h,0}(\omega_T^m),
\qquad
\eta_{dual}(T)=\|z_T\|_{E,\omega_T^m}.
```

局部残差分割必须满足 `sum_T R_{h,T}=R_h`。先比较 `m=1,2` 对排序的影响；人工 patch 边界使用齐次 Dirichlet，物理边界保留正定 Riesz 内积对应的自然边界处理。该指标只负责校准，不进入大规模生产循环。

这个构造有两个刻意保留的特征：

1. 测试空间是无约束的 `V_{h,0}(omega_T^m)`，不是 LOD 核空间；
2. 右端是人为分割后的 `R_{h,T}`，不同局部问题并不直接使用同一个全局残差。

因此它能检查残差装配、局部排序和标记行为，但目前没有从局部量到
`R_h|_{W_h}` 的稳定分解双边等价。它不能直接替代下面的生产候选。

#### 4.4.2 核空间受限的局部代数对偶估计子

本节加入 `helmholtz_lod_aposteriori_report.pdf` 中的构造。该报告不存入仓库；计划书只记录可复现的数学定义和实施要求。

令

```math
W_h=V_h\cap\ker I_H,
\qquad
W_{h,j}=\{w_h\in W_h:\operatorname{supp}w_h\subset\omega_j\},
```

其中 `omega_j` 是固定阶的重叠粗 patch，优先使用粗节点星形 patch及其固定层扩张。对当前 LOD 解 `U=u_{H,h}^{ell}`，不再分割残差，直接限制同一个全局代数残差

```math
R_U(w_h)=F(w_h)-a_k(U,w_h)
```

到每个 `W_{h,j}`。以

```math
b_k(v,w)
=(A\nabla v,\nabla w)_\Omega+k^2(nv,w)_\Omega
```

为正定 Hermitian 内积，求局部 Riesz 代表元

```math
b_k(\xi_j,w_h)=R_U(w_h)
\quad\forall w_h\in W_{h,j},
\qquad
\eta_j=\|\xi_j\|_E,
\qquad
\eta_{LOD}^2=\sum_j\eta_j^2.
```

这里 `w_h in ker I_H` 必须由核空间基或拉格朗日乘子严格施加，不能用无约束 `V_{h,0}(omega_j)` 代替。矩阵形式为

```math
\begin{bmatrix}
B_{E,j} & Q_j^*\\
Q_j & 0
\end{bmatrix}
\begin{bmatrix}\xi_j\\\lambda_j\end{bmatrix}
=
\begin{bmatrix}r_j\\0\end{bmatrix},
\qquad
B_{E,j}=K_j+k^2M_j,
\qquad
\eta_j^2=\xi_j^*B_{E,j}\xi_j,
```

其中 `Q_j` 是 `I_H` 在局部自由度上的约束矩阵，`r_j` 是全局残差向量对该 patch 自由度的限制。该局部矩阵虽然是鞍点形式，但主块是正定能量矩阵；可复用现有 DirectSchur 约束求解基础设施。

理论目标不是 `sum_j R_{h,j}=R_h`，而是验证核空间的稳定重叠分解：对任意 `w_h in W_h`，存在 `w_{h,j} in W_{h,j}` 使

```math
w_h=\sum_j w_{h,j},
\qquad
\sum_j\|w_{h,j}\|_E^2\le C_{sd}^2\|w_h\|_E^2,
```

且任意 `v_j in W_{h,j}` 满足

```math
\left\|\sum_jv_j\right\|_E^2
\le C_{ov}^2\sum_j\|v_j\|_E^2.
```

对当前 \(I_H=E_H\Pi_H^{dg}\) 和自适应 NVB 非拟一致粗网格，\(C_{sd}\) 是否保持与 \(H/h/k\) 无关仍需单独证明并做数值谱检查；报告中的结论不能在未核查这一假设时直接宣称成立。

由此得到

```math
\frac1{C_{ov}}\eta_{LOD}
\le \|R_U|_{W_h}\|_{W_h'}
\le C_{sd}\eta_{LOD}.
```

这正是它与 4.4.1 的本质区别：局部量来自“同一残差在重叠核子空间上的限制”，双边等价依赖稳定空间分解，而不是依赖残差的单元归属或半边分配。

#### 4.4.3 从残差范数到 LOD 误差及伴随 corrector 证书

该估计子首先控制相对于固定细空间的误差

```math
e_h=u_h-U.
```

若近似测试 corrector 为 `tilde C_h^*`，定义

```math
\delta_h^*
=\sup_{0\ne v_H\in V_H}
\frac{\|(C_{\infty,h}^*-\widetilde C_h^*)v_H\|_E}
{\|v_H\|_E}.
```

在核空间 coercivity、理想 Petrov-Galerkin inf-sup 稳定性和上述稳定分解成立时，目标双边界为

```math
\frac1{C_aC_{ov}}\eta_{LOD}
\le\|u_h-U\|_E
\le C_{sd}
\left(c_W^{-1}+\frac{C_\Pi\delta_h^*}{\alpha_{H,\infty}(k)}\right)
\eta_{LOD}.
```

理想测试 corrector 下 `delta_h^*=0`，可靠性常数不再携带全局 Helmholtz 稳定性增长。局部化、fine-corrector 离散和局部线性求解误差统一进入 `delta_h^*`，因此仅计算 `eta_LOD` 还不够，必须同步实现可计算证书。

对任意 `v_H`，定义伴随 corrector 方程缺陷

```math
D_{v_H}^*(w_h)
=a_k\!\left(w_h,(I-\widetilde C_h^*)v_H\right),
\qquad w_h\in W_h.
```

在每个 `W_{h,j}` 上求

```math
b_k(w_h,\zeta_j(v_H))=D_{v_H}^*(w_h)
\quad\forall w_h\in W_{h,j},
```

并定义

```math
\vartheta(v_H)^2=\sum_j\|\zeta_j(v_H)\|_E^2,
\qquad
\Theta_h^*
=\sup_{0\ne v_H\in V_H}
\frac{\vartheta(v_H)}{\|v_H\|_E}.
```

对粗基逐列组装后，

```math
(\Theta_h^*)^2=\lambda_{\max}(G^*,M_E),
```

可用稠密广义特征值、Lanczos 或幂迭代计算。目标验证关系为

```math
\frac1{C_aC_{ov}}\Theta_h^*
\le\delta_h^*
\le\frac{C_{sd}}{c_W}\Theta_h^*.
```

最终以无量纲证书

```math
q^*(k)=2C_FC_\Pi C_{st,h}(k)k\,\frac{C_{sd}}{c_W}\Theta_h^*
```

控制测试空间误差；若 `q^*(k)<=q_0`，则得到关于 `k` 一致的目标形式

```math
\frac1{C_aC_{ov}}\eta_{LOD}
\le\|u_h-U\|_E
\le C_{sd}(c_W^{-1}+q_0)\eta_{LOD}.
```

在报告采用的归一化 `c_W=1/3` 下，上一式中的证书写成 `q^*(k)=6 C_F C_Pi C_st,h(k) k C_sd Theta_h^*`。一般系数情形应先按代码的 `A,n,beta` 和能量范数重新核对常数，不能机械照搬数值 `3` 和 `6`。第一版先验证缩放关系以及 `q^*` 随 `ell`、corrector 容差的趋势，不把未知理论常数拟合成看似精确的证书。

若全局粗 Petrov-Galerkin 系统没有精确求解，还需独立记录

```math
\eta_{alg}
=\sup_{0\ne\widetilde y_H\in\widetilde Y_H}
\frac{|R_U(\widetilde y_H)|}{\|\widetilde y_H\|_E}.
```

`eta_alg` 是全局线性系统的代数求解误差，不是 `eta_LOD`，也不是 `Theta_h^*`。它应作为可靠性上界的附加项，并用于设置粗系统停止容差。

#### 4.4.4 候选选择规则

候选指标按以下标准选择：全局 effectivity 的稳定性、与核空间 `eta_j` 的 Spearman 相关性、Dörfler 标记集合重合率，以及一次估计的时间和内存。首轮在看数据前登记建议门槛：局部 Spearman 相关系数的跨样本中位数不低于 `0.8`，与核空间代数对偶指标 Dörfler 集合的能量重合率不低于 `0.7`。阈值可在试运行前调整，但不能在查看正式结果后回调。

若不同指标各有优势，生产标记与理论诊断可以保留两个不同字段，不能混用同一个名称 `eta_H`。进入第三阶段后，`eta_H` 专指本阶段冻结的生产指标 `eta_X`。

### 4.5 既有先验可靠性包络及其适用条件

本节保留原计划中的参数化先验包络，作为与 4.4.2--4.4.3 代数对偶理论路线的次级对照，不再作为唯一的生产估计子冻结依据。优先顺序是：

1. 先验证 \(\eta_{LOD}\) 与 \(R_U|_{W_h}\) 的双边等价；
2. 再验证 \(\Theta_h^*\)/\(q^*\) 对伴随 corrector 质量的证书作用；
3. 最后比较本节先验包络是否提供额外、可迁移且不过松的信息。

对选定指标 `eta_X`，在 `eta_X>0` 时定义观测可靠性比值

```math
C_{rel}^{obs}=\frac{\|u_h-u_{H,h}^{\ell}\|_E}{\eta_X}.
```

若 `eta_X=0` 且离散误差也在容差内，则记为已收敛样本；若 `eta_X=0` 但离散误差非零，则该候选估计子直接判定失败，不进入常数拟合。

拟检验下列先验结构，并将右端记作 `C_rel^pred`：

```math
C_{rel}^{obs}\le
\left(C_1+C_2\rho^\ell+C_3q^{\frac12-\varepsilon}\right)
\frac{1+C_{\ell,k}^2}{1-C_{\ell,k}^2},
```

其中

```math
C_{\ell,k}=
\frac{C_4+\left(C_5\rho^\ell+C_6q^{\frac12-\varepsilon}\right)k^{d+1}}
{C_7-\left(C_5\rho^\ell+C_6q^{\frac12-\varepsilon}\right)k^{d+1}},
\qquad d=2.
```

这里的 \(\rho\in(0,1)\) 表示 corrector 的局部化衰减率；为避免与 Robin 系数 \(\beta\) 混淆，不再沿用原草稿中的衰减符号 \(\beta\)。二维实验中使用 `d+1=3`，不再使用含义不明确的 `n+1`。理论安全版本取

```math
q=q_{max}:=\max_{T\in\mathcal T_H}\max_{\tau\subset T}\frac{h_\tau}{H_T}.
```

当细网格准一致时，`q_max` 等于 `h/H_min`。它可能很保守，但保持全局最坏情形意义。另记录仅作诊断的指标加权有效比值

```math
q_{eff}^2=
\frac{\sum_T\eta_X(T)^2q_T^2}{\sum_T\eta_X(T)^2},
\qquad
q_T=\max_{\tau\subset T}h_\tau/H_T,
```

但在没有证明前，`q_eff` 不能替代可靠性上界中的 `q_max`。

其中 `0<epsilon<1/2` 是预先固定的小量，首轮取 `epsilon=0.01`。`C_1,...,C_7` 假设不随 `k,ell,h,H` 变化，但允许依赖区域、系数椭圆性上下界、网格形状正则性和拟插值稳定常数。不同系数族或不同区域不能默认共用同一组常数，必须作为迁移验证处理。

该公式只在以下稳定性条件同时成立时解释为可靠性候选：

```math
C_7-\left(C_5\rho^\ell+C_6q^{\frac12-\varepsilon}\right)k^3>0,
\qquad
0\le C_{\ell,k}<1.
```

若条件失败，应报告“超出先验公式适用区间”，不能继续计算可能为负或发散的因子。固定 `h` 而持续加细 `H` 会使 `q_max=h/H_min` 增大；因此该上界可能在 `H` 接近 `h` 时退化，即使真实全离散误差最终趋于零。这反映的是充分条件不锐利，而不是算法误差增大。第一阶段应设 `q_max<=q_limit` 的尺度分离保护，并单独报告接近 `H=h` 的退化实验。

### 4.6 可靠性结构的数值辨识与留出验证

`C_1,...,C_7` 不能在同一小组数据上无约束拟合后再用该组数据宣称验证成功。采用以下流程：

1. **独立估计衰减率**：由相邻 oversampling 层的 corrector 差

   ```math
   D_\ell=\left(\sum_T\|Q_{T,h}^{\ell+1}-Q_{T,h}^{\ell}\|_E^2\right)^{1/2}
   ```

   在进入饱和区前拟合 \(D_\ell\approx C\rho^\ell\)，得到 `rho` 及置信区间。
2. **固定结构参数**：取 `epsilon=0.01`，二维指数固定为 `k^3`；先不同时辨识指数、`rho` 和七个常数。
3. **训练数据**：在单位正方形上改变 `k, ell, h` 和统一 `H`，用约束优化寻找满足所有训练样本的最小上包络常数组；约束所有 `C_i>=0`、分母正且 `C_{ell,k}<1-delta`。由于 `C_4,...,C_7` 同比缩放不改变 `C_{ell,k}`，若理论没有预先给定 `C_7`，数值辨识时固定 `C_7=1`，等价地只拟合各常数与 `C_7` 的比值。
4. **留出数据**：使用未参与拟合的局部源位置、自适应网格和至少一个波数检验

   ```math
   C_{rel}^{obs}\le C_{rel}^{pred}.
   ```

5. **敏感性分析**：分别去掉 `rho^ell`、`q^(1/2-epsilon)` 或 `k^3` 项，比较留出集违反次数和上界松弛度，判断每一项是否被数据支持。

由于这些常数参数可能不可辨识，首选目标不是得到唯一 `C_i`，而是检验该函数结构能否形成跨参数有效的上包络。报告训练/留出样本数、最大违反率、`max(Cobs/Cpred)`、中位松弛度，以及稳定性裕量

```math
m_{den}=C_7-\left(C_5\rho^\ell+C_6q^{\frac12-\varepsilon}\right)k^3,
\qquad
m_{stab}=1-C_{\ell,k}.
```

若常数组高度不唯一，则增加参数归一化或合并常数，而不是继续扩大训练集后只报告一个拟合优度。首轮留出门槛预先设为 `max(C_rel^obs/C_rel^pred)<=1.01`；同时报告 `median(C_rel^pred/C_rel^obs)`，若大于 `10`，即使没有违反也标记为“上界过松”，不能据此宣称估计子实用。

### 4.7 Helmholtz 波数加权候选

在标准三个指标完成校准后，再比较实验性权重

```math
\alpha_\tau=\min(h_\tau,k^{-1}),
\qquad
\alpha_e=\min(h_e^{1/2},k^{-1/2}).
```

该权重不能直接称为波数鲁棒估计子。必须同时报告离散 inf-sup 值、`kH_T`、`kh_tau`、effectivity 和第 4.5 节稳定性裕量随 `k` 的变化。

### 4.8 校准后自适应循环

只有候选指标通过 4.4--4.6 的校准后，才固定一个生产标记指标并执行：

1. **SOLVE**：构造当前 corrector、trial/test basis，求 `u_{H,h}^{ell}`。
2. **ESTIMATE**：装配唯一一次细单元和细边残差，同时导出三个候选聚合；生产循环选择其中已冻结的一种。
3. **MARK**：Dörfler 标记最小集合 `M`，满足

   ```math
   \sum_{T\in\mathcal M}\eta_X(T)^2
   \ge\theta\sum_{T\in\mathcal T_H}\eta_X(T)^2,
   \qquad 0<\theta<1.
   ```

4. **REFINE**：对 `M` 做 NVB 加细和相容闭包，更新父子映射。
5. **VERIFY**：检查嵌套性、无悬挂点、尺度分离、corrector 约束残差、Petrov 残差和可靠性公式适用条件。

建议先测试 `theta=0.3,0.5,0.7`。在该定义下，`theta` 越大通常标记越多。停止准则为 `eta_X<=tol`、达到粗自由度上限、`q_max>q_limit`、所有候选达到细网格深度，或达到最大迭代数。


## 5. 单位正方形实验矩阵

### 5.1 问题设置

- 区域：`Omega=(0,1)^2`。
- 边界：全边界齐次阻抗 Robin。
- 系数：先用 `A=n=beta=1`，再加入与当前椭圆 LOD 一致的高对比系数。
- 右端项：平滑高斯源、局部尖峰源、靠近边界的局部源。
- 波数：先 `k=4,8,16`，满足稳定性基线后再扩展。
- 细网格：固定全局 NVB 主网格，确保参考误差充分小。
- 校准子阶段：对 `k=4,8,16`，交叉改变统一 `H`、fine gap 和 `ell`，并预先划分训练集与留出集。`rho` 只由 corrector 层差独立估计。
- 自适应子阶段：冻结已选指标、固定主细网格和统一 `ell=ceil(log2(k))+c`，使 localization error 小于粗网格误差；自适应数据不得回流重新拟合可靠性常数。

### 5.2 对照组

在相同 `h,k,ell` 下比较：

1. 统一加细粗网格 LOD；
2. 残差驱动自适应粗网格 LOD；
3. 细网格 P1 参考解 `u_h`；
4. 小问题上的局部对偶残差指标；
5. 三种残差指标及可靠性先验结构的训练/留出对照。

### 5.3 输出指标

每轮输出 CSV/JSON：

- 自适应轮次、粗节点/单元/自由度、最小/最大粗层级；
- 用户标记数、闭包新增数、实际加细数；
- `||u_h-u_LOD||_E`、相对能量误差和 L2 误差；
- `eta_fine`、`eta_mix`、`eta_macro` 及各自的 effectivity；
- `eta_j`、`eta_LOD`、核空间全局对偶范数、双边比值及 patch 重叠数；
- `Theta_h_star`、`q_star`、`eta_alg` 以及所用 corrector/粗系统求解容差；
- `C_rel^obs`、`C_rel^pred`、`q_max`、`q_eff`、`m_den` 和 `m_stab`；
- 局部指标与局部真误差/对偶指标的 Spearman 排序相关系数；
- 标记集合与对偶指标理想标记集合的重合率；
- corrector、估计、标记、加细、粗求解、reference 和总时间；
- corrector 重算/复用数量、峰值 RSS、Petrov 与 corrector 残差；
- `kH_T` 的最小值、分位数和最大值。

核心图形为误差/估计子对粗自由度、误差对总时间、effectivity 对迭代步，以及粗网格和指标热图。

### 5.4 第一阶段验收标准

- 所有自适应粗网格满足 `V_H subset V_h` 且无悬挂点。
- 使用相同系数表示和积分规则时，积分残差与代数残差的相对差不超过 `1e-10`。
- 当 `T_H=T_h` 时，离散误差与残差指标同时趋近机器精度。
- 对至少三个源项，自适应曲线不劣于统一加细；局部源项上应出现明确收益。
- effectivity 不随自适应轮次无控制发散；
- 选定指标在留出集上满足预先规定的可靠性违反率，并报告上界松弛度；
- 在可显式构造全局 `W_h` 的小网格上，`eta_LOD` 与 `||R_U|W_h||` 满足 4.4.2 的双边关系；
- `Theta_h_star` 随 `ell` 增大和局部求解容差收紧而下降，并能区分 localization/corrector 缺陷与 `eta_alg`；
- 使用第 4.5 节公式的样本满足分母为正且 `C_{ell,k}<1`，失效样本被显式分类而非静默丢弃；
- 冷启动完整重建与增量更新得到相同网格、corrector、解和指标。

## 6. 第二阶段：L 型区域低正则实验

采用

```math
\Omega=(-1,1)^2\setminus([0,1]\times[-1,0])
```

并使用与该区域相容的初始 NVB 标记。主实验保持全边界齐次阻抗 Robin，使用平滑源或靠近重入角的局部源，由足够细的 `V_h` 解作为参考。几何重入角会降低解的正则性，预期自适应网格集中在原点附近。

另设一个“估计子公式验证”实验：采用已知 `r^(2/3)` 型奇异制造解，并施加由该解导出的非齐次边界数据。此实验只用于核对收敛阶和边界残差符号，不与纯齐次 Robin 主实验混为一组。

比较内容：

- 统一粗网格 LOD 与自适应粗网格 LOD 的误差-DOF、误差-时间曲线；
- 普通粗 P1 FEM、统一 LOD、自适应 LOD 的 pollution 与角点奇异误差；
- `theta`、`ell` 和 `k` 对加细区域的影响；
- 网格是否同时响应角点奇异性、源项位置和高对比系数界面。

## 7. 第三阶段：`H/h/ell` 三部分估计子

建议先用误差分解指导算法，而不是一开始声称三部分严格正交：

```math
u-u_{H,h}^{\boldsymbol\ell}
=(u-u_h)
+(u_h-u_{H,h}^{\boldsymbol\infty})
+(u_{H,h}^{\boldsymbol\infty}-u_{H,h}^{\boldsymbol\ell}).
```

三项分别对应 fine discretization、coarse/model reduction 和 localization。目标形式为

```math
\|u-u_{H,h}^{\boldsymbol\ell}\|_E
\lesssim C_{stab}(k)
\bigl(\eta_h+\eta_H+\eta_\ell+\operatorname{osc}(f,A,n)\bigr).
```

### 7.1 `H` 指标

`eta_H(T)` 指示粗空间是否需要在 `T` 附近增加自由度。基线采用第 4 节的聚合残差；进阶版本使用一次 NVB 后新增的层次细节空间 `Z_T^H subset V_h`：

```math
\eta_H(T)=\|R_h\|_{(Z_T^H,E)'}.
```

这可以通过很小的局部 Riesz 问题计算，并直接衡量“加细该粗单元可能消除多少残差”。

### 7.2 `h` 指标

`eta_h(tau)` 指示当前细空间不足。可采用经典细网格体/跳跃残差，或在候选细化产生的层次细节空间 `Z_tau^h subset V_{h+}` 上定义

```math
\eta_h(\tau)=\|R\|_{(Z_\tau^h,E)'}.
```

因为 `u_h` 本身不能估计连续误差，验证 `eta_h` 时需要更细参考空间 `V_{h+}` 或制造解。细网格一旦局部改变，所有粗网格必须仍是它的祖先网格。

### 7.3 `ell` 指标

允许每个粗单元使用不同 `ell_T`。先采用可验证的双层差分指标

```math
\eta_\ell(T)^2
=\sum_{j\in\mathcal N(T)}|U_j|^2
\left\|
(Q_{T,h}^{\ell_T+1}-Q_{T,h}^{\ell_T})\varphi_j
\right\|_{E,\omega_T^{\ell_T+1}}^2,
```

并对伴随 corrector 构造对应指标。Petrov-Galerkin 情形至少记录 primal 和 adjoint 两部分，不能只检查原 corrector。

更便宜的候选是人工 patch 边界上的 corrector 通量泄漏：

```math
\eta_{\ell,flux}(T)^2
=\sum_{e\subset\partial\omega_T^{\ell_T}\setminus\partial\Omega}
h_e\|A\nabla Q_{T,h}^{\ell_T}u_H\cdot\nu\|_{L^2(e)}^2.
```

通量指标必须先与 `ell_T+1` 差分指标做相关性验证，不能直接作为最终估计子。

### 7.4 是否允许不同 patch 使用不同 `h` 和 `ell`

- **不同 `ell_T`：可行。** 每个单元 corrector 本来就是独立局部问题，全局 corrector 仍可按单元求和。分析中需要控制 `min_T ell_T`、patch 重叠常数，并分别处理 primal/adjoint localization error。
- **同一全局自适应 `V_h` 上的有效局部 `h`：可行且优先。** 每个 patch 自动继承其覆盖区域的非均匀细网格，但所有 corrector 仍属于同一个共形空间 `V_h`。
- **每个 patch 拥有彼此独立的局部细网格：原则上可行，但不是当前架构的直接推广。** 它需要局部 `I_H`/kernel、系数投影、跨网格 lift/scatter、重叠区一致积分和误差控制。该方案放在最后，不作为前三阶段依赖。

### 7.5 联合决策

先分别验证三类指标，再按“预计误差下降/新增计算代价”排序：

```math
s_H(T)=\eta_H(T)^2/\Delta\operatorname{cost}_H(T),\quad
s_h(\tau)=\eta_h(\tau)^2/\Delta\operatorname{cost}_h(\tau),\quad
s_\ell(T)=\eta_\ell(T)^2/\Delta\operatorname{cost}_\ell(T).
```

每轮只执行收益最高且满足稳定性约束的一类或一批动作。特别地，必须先满足最低 oversampling 条件和 `kH_T` 条件，再允许成本模型削减 `ell_T`。

## 8. 模块化设计

建议新增 `include/helmholtz/adaptive/` 与 `src/helmholtz/adaptive/`，保持数学模块与实验驱动解耦。

### 8.1 数据模块

- `AdaptiveMeshHierarchy`：稳定单元 ID、父子关系、层级、粗到细覆盖、NVB 参考边及版本号。
- `AdaptiveLodConfig`：`k,h,ell,theta,tol,max_dofs`、估计子和复用策略。
- `AdaptiveLodState`：当前粗网格、模型、解、指标和缓存句柄。
- `AdaptiveIterationRecord`：误差、估计子、时间、内存、标记及复用统计。

### 8.2 算法模块

- `HelmholtzResidualAssembler`：一次装配 broken element、内部边和 Robin 边的原始贡献；
- `HelmholtzResidualEstimator`：从同一组原始贡献生成 fine/mix/macro 三种 coarse aggregation；
- `PartitionedDualResidualCalibrator`：保留 4.4.1 的无约束、分割残差 Riesz 诊断；
- `KernelDualResidualEstimator`：在 `W_{h,j}` 上求 4.4.2 的受约束局部 Riesz 问题并输出 `eta_j/eta_LOD`；
- `AdjointCorrectorCertificate`：组装 `G_star,M_E`，计算或估计 `Theta_h_star` 和 `q_star`；
- `GlobalAlgebraicResidualEstimator`：计算粗测试空间上的 `eta_alg`，与 LOD 估计子分开报告；
- `ReliabilityEnvelopeCalibrator`：离线估计 `rho`、约束拟合上包络并执行留出验证，不参与生产求解；
- `DoerflerMarker`：确定性排序、最小集合和并列值规则。
- `AdaptiveNvbRefiner`：标记、闭包、父子映射及嵌套检查。
- `AdaptiveHelmholtzLodDriver`：只负责 SOLVE-ESTIMATE-MARK-REFINE 状态机。
- `AdaptiveResultWriter`：CSV/JSON、网格快照和可恢复 checkpoint。

估计器不得直接修改网格；refiner 不得访问 PDE；driver 不得重新实现 corrector 数学。

## 9. 增量计算与缓存设计

粗网格局部加细后，不应无条件重算所有 corrector。为每个粗单元 corrector 建立依赖签名：

```text
(k, coefficient_revision, fine_mesh_revision,
 interpolation_revision, coarse_element_id,
 patch_element_ids, ell_T, primal_or_adjoint)
```

仅当 patch、局部粗基函数、拟插值限制或算子发生变化时使缓存失效。初期实现必须同时提供：

- `full_rebuild`：每轮完全重建，作为正确性金标准；
- `incremental`：复用未受影响的 corrector、局部 pattern 和 factorization；
- 两种模式逐项比较 basis、粗矩阵、解、指标和残差。

可进一步复用：

- 相同 `A,k,T_h,T_H,boldsymbol ell`、不同 `f`：复用全部 corrector、basis 和粗矩阵分解；
- 仅局部粗加细：复用 patch 依赖集合未变化的 corrector；
- patch pattern 相同：复用符号分析和多右端分解；
- reference solve：同一 `A,k,T_h` 下缓存细矩阵分解，避免每轮重分解。

缓存应有内存预算和 LRU/显式释放策略；结果文件记录 cache hit、miss、eviction 和节省时间。

## 10. 测试计划

### 10.1 单元测试

- 复值体残差、内部边通量跳跃及 Robin 边符号的手算网格测试；
- 每条内部边只计一次，法向翻转不改变指标；
- 积分残差向量与 `b_h`-A_hu 一致；
- 粗单元内部细边、公共粗边半分配和 Robin 边唯一归属满足局部求和恒等式；
- `fine/mix/macro` 三个指标由同一原始残差贡献生成，切换指标不重复装配；
- Dörfler 集合满足阈值、最小性和确定性；
- NVB 闭包无悬挂点，父子面积守恒；
- 每轮 `V_H subset V_h` 的基函数再现测试；
- `T_H=T_h` 时 corrector、离散误差和残差退化正确；
- primal/adjoint 指标的共轭和独立装配测试；
- checkpoint 恢复后下一轮结果一致。

### 10.2 集成和 golden 测试

- 小方形网格上与稠密参考残差对偶范数比较；
- 同一全局残差限制到 `W_{h,j}` 的局部右端与直接矩阵限制一致，且不依赖 broken residual 的边归属；
- 受约束 Riesz 解满足 `Q_j xi_j=0`、局部方程残差和能量恒等式 `eta_j^2=xi_j^* B_E,j xi_j`；
- 显式全局核空间小问题验证 `eta_LOD` 与 `||R_U|W_h||` 的上下界及 patch 重叠常数；
- 逐粗基组装的 `G_star` 与直接算子范数一致，广义最大特征值等于 `Theta_h_star^2`；
- 人为放松 corrector 和粗系统容差，分别验证 `Theta_h_star` 与 `eta_alg` 响应不同误差源；
- 可靠性公式在分母非正或 `C_{ell,k}>=1` 时返回明确的不可用状态；
- 固定训练/留出数据的 envelope golden test，防止拟合流程发生数据泄漏；
- 固定标记序列下的网格、指标、解和 corrector golden 数据；
- `full_rebuild` 与 `incremental` 在多轮自适应中的逐轮一致性；
- 相同模型、两个不同 RHS 确认 corrector 与 factorization 未重建；
- L 型区域初始网格和多轮 NVB 相容性测试。

## 11. 实施顺序与门槛

1. **网格前置改造**：显式自适应 `T_H`、稳定 ID、NVB 树、粗细覆盖关系和 `V_H subset V_h` 检查。
2. **原始残差装配**：一次生成 broken 体残差、内部边跳跃和 Robin 边贡献，并与代数残差逐自由度交叉验证。
3. **三个指标实现**：从同一原始贡献生成 `eta_fine`、`eta_mix` 和 `eta_macro`，验证边归属及全局求和。
4. **核空间代数对偶估计子**：保留既有分割残差 Riesz 校准，新增 `W_{h,j}` 受约束局部问题；在可显式构造全局 `W_h` 的小网格上验证双边等价。
5. **伴随 corrector 与代数求解证书**：实现 `Theta_h_star` 广义特征值和独立的 `eta_alg`，改变 `ell`、局部容差和粗系统容差验证误差源分离。
6. **强残差代理比较**：比较 `fine/mix/macro` 与 `eta_j` 的排序、Dörfler 集合、时间和内存；第 4.5 节参数化包络只作次级对照。
7. **冻结生产指标**：只有核空间双边校准和 corrector 证书通过后，才按预登记门槛选择标记指标；若无候选通过则停止。
8. **方形 H-only 自适应**：完成可复现 SOLVE-ESTIMATE-MARK-REFINE，对比统一粗加细。
9. **增量复用**：在 full rebuild 金标准通过后加入局部 corrector/cache 更新。
10. **L 型区域**：不重新拟合常数，先测试方形阶段选择的指标，再完成低正则 uniform/adaptive 对比。
11. **变量 `ell_T`**：先实现 `ell+1` 差分指标，再验证通量泄漏替代指标。
12. **全局非均匀 `V_h`**：实现 `h` 指标和保持嵌套的联合加细。
13. **三类联合决策**：加入成本模型、稳定性保护和联合 benchmark。
14. **独立 patch 细网格**：仅在前述方案显示内存或时间瓶颈时研究。

每一步必须保持现有椭圆 LOD、Helmholtz FEM、corrector 和 Petrov-Galerkin 测试通过。任何增量优化都必须先与 full rebuild 比较数值结果，再报告加速比。

## 12. 首轮建议参数

为控制开发成本，第一轮建议：

```text
domain       = unit square
coefficient  = A=n=beta=1
k            = 8
initial H    = 2 or 3 NVB levels
master h     = 9 or 10 NVB levels
ell          = ceil(log2(k)) + 1
theta        = 0.5
max steps    = 8
reference    = cached solve on V_h
mode         = two-sided Petrov-Galerkin
```

先用局部高斯源产生非均匀误差分布。该设置通过后，再加入 `k=4,16`、高对比系数和 L 型区域。

## 13. 需要重点审阅的数学问题

1. `eta_mix` 中粗单元内部细边是否继续使用 `h_e`，还是理论推导要求统一改为 `H_T`。
2. 对 `I_H=E_H Pi_H^dg` 和非拟一致自适应 NVB 网格，如何构造局部右逆并证明 `W_h=sum_j W_{h,j}` 的稳定分解常数与 `H/h/k` 无关。
3. 核空间 patch 应选粗节点星形、粗单元 patch 还是固定层扩张，怎样兼顾 `C_sd/C_ov`、标记局部性和求解成本。
4. 4.4.1 的分割残差 Riesz 量是否只保留为装配诊断，以及它与核空间 `eta_j` 在何种条件下可能等价。
5. 第 4.5 节先验式中的二维波数指数是否确为 `k^3`，以及 `C_1,...,C_7` 是否能依据理论关系减少自由参数。
6. 全局理论上界是否必须使用 `q_max=h/H_min`；若希望使用局部 `q_T`，需要怎样重新组织可靠性求和。
7. Helmholtz 离散 inf-sup 因子与 `C_{ell,k}` 的关系，以及实验性 `k` 加权应在哪个阶段引入。
8. `eta_ell` 是否需要同时包含 primal 和 adjoint corrector，及两者的组合方式。
9. `eta_H+eta_h+eta_ell` 最终以可加可靠上界为目标，还是先作为误差下降/成本决策指标。

在这些问题确认前，不把候选公式写成已证明的后验误差上界。
