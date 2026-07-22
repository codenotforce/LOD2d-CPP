# Helmholtz LOD 非拟一致 NVB 网格局部逆不等式数值验证计划书

> 状态日期：2026-07-22
>
> 当前状态：本地 WSL 与 EPYC 服务器实验均已完成；服务器深细网格达到 `L_h=16`。
>
> 适用范围：二维齐次阻抗 Robin 边界 Helmholtz 方程、复值 P1 Petrov-Galerkin LOD、全局协调但非拟一致的 NVB 粗网格。

## 执行状态

| 阶段 | 内容 | 状态 | 退出条件 |
|---|---|---|---|
| P0 | 数学口径、网格族和失败判据冻结 | 已完成 | 本计划书评审通过 |
| P1 | `C=0`、局部矩阵和复广义特征值校准 | 已完成 | CTest smoke gate 和负控制通过 |
| P2 | 确定性非拟一致 NVB 网格生成与不变量检查 | 已完成 | `Gamma=1,...,8`，相邻比始终为 `sqrt(2)` |
| P3 | Helmholtz trial/test LOD 局部逆常数实现 | 已完成 | 原/伴随空间、残差和局部谱可复核 |
| P4 | 固定细空间的 WSL 主扫描 | 已完成 | 三类位置、三种 `ell` 和三种标记比例已扫描 |
| P5 | 固定粗网格的 `h` 极限与波数扫描 | 已完成 | 本地 `L_h=11,...,14`；服务器 `L_h=14,15,16`；`k=1,2,4` |
| P6 | 结论、原始 CSV、元数据写入 `DEVELOPMENT.md` | 已完成 | 本地与服务器原始数据均位于约定结果目录 |

本计划只规定实验设计，不预设局部逆不等式一定成立。正式结果必须区分“在已测试参数内未观察到增长”和“已证明与网格无关”；数值实验只能提供前者或给出反例。

### 2026-07-21 执行结论

- 所有构造的粗空间都通过 nodal、element、DG prolongation 和 `I_HP_H=I` 四重检查；主扫描最大嵌套残差为 `7.70e-15`，因此确认 `V_H subset V_h`。
- `C=0` 使用局部 `H_T` 时恒为 `8.485281`；错误使用全局 `H_max` 时从 `8.485281` 增至 `67.882251`，负控制符合预期。
- 同单元分母的 trial/test 常数对细空间和质量正谱阈值敏感；在固定最终粗网格上随 `L_h=11,12,13,14` 为 `97.30,101.44,124.31,132.50`。当前数值证据不支持宣称该形式成立。
- 一层 patch 分母在相同 `h` 扫描中为 `11.98,12.41,12.45,12.57`，表现出清晰平台；但它从均匀网格到强非拟一致网格仍有初始增长，边界链最后三点也未通过预注册平台判据，因此只能记录为“有条件支持、尚不能判定 grading-independent”。
- fine level 15 单点在约 `11,645,504 KiB` RSS 时被 signal 9 终止；该点保留为服务器任务，不改变本地结论。
- 完整数值和解释见 `DEVELOPMENT.md` 第 28 节及 `results/helmholtz_local_inverse/`。
- 已追加最大值反馈加细：每轮取 `argmax Q_T` 或 `argmax Q_T_patch`，分别测试仅加细目标单元和连同一层顶点邻域加细。四组实验均保持固定细网格并通过 `V_H subset V_h` 检查；其中 `argmax-patch`、零邻层在 6 步后以 25 个单元达到等级范围 `3:8` 和 `Gamma=5.657`，作为反馈式强非拟一致主案例。
- 分母已进一步与 oversampling 半径匹配：`ell=3` 时主 patch 指标使用 `omega_T^3`，`argmax-patch` 也由该指标驱动；`omega_T^1` 仅作为诊断保留。在固定的等级 `3:5` 粗网格上，`L_h=11,12,13` 给出 `q_max=0.125,0.0884,0.0625`，匹配分母最大值为 `3.5289,3.5181,3.5124`。`L_h=14` 反馈点在 `11,668,680 KiB` 被 signal 9 终止，作为内存上限记录保留。

### 2026-07-22 EPYC 服务器结论

- 固定同一个等级 `3:8` 粗网格，`Gamma=4 sqrt(2)`，服务器 `L_h=14,15,16` 的 `q_max` 分别为 `0.125,0.088388,0.0625`。
- 与 `ell=3` 匹配的局部逆常数为

```math
C_{\mathrm{inv},3}(h)
=\max_{T\in\mathcal T_H}H_T
\sup_{0\ne v\in V_{H,3}^{\mathrm{ms}}}
\frac{\|\nabla v\|_{L^2(T)}}
     {\|v\|_{L^2(\omega_T^3)}}
=3.2019028413,\ 3.2017193578,\ 3.2013850546.
```

- 最后两级相对变化为 `1.044e-4`，总变化为 `1.617e-4`；因此在 `h_T/H_T<=1/16` 时数值上可记录

```math
\|\nabla v\|_{L^2(T)}
\le 3.202\,H_T^{-1}\|v\|_{L^2(\omega_T^3)}
```

  对所有已测试粗单元和离散 trial LOD 函数成立。`3.202` 是实验上界，不是严格证明的网格无关常数。
- 同单元分母从 `132.0481` 增至 `162.1273`，增长 `22.78%`，不满足 `h` 平台判据。
- `L_h=16` 的最大值反馈序列位于 `[3.1788,3.5400]`，但最大值在对称边界热点间迁移，最终仅达到等级 `3:5`、`Gamma=2`；它支持反馈路径上的有界性，强非拟一致证据仍由固定 `3:8` 网格提供。
- 服务器所有行均有 `fixed_fine_mesh=1`；最大 nesting、primal corrector、constraint 残差分别为 `1.38e-13,1.32e-13,6.70e-15`，确认 `V_H subset V_h` 和局部谱计算正确。
- `L_h=16` hscan 峰值约 `180.8 GiB`，反馈实验峰值约 `237.1 GiB`；匹配常数已经收敛，不再执行 `L_h=17`。
- 原始服务器数据位于 `results/helmholtz_local_inverse_server/`，完整解释见 `DEVELOPMENT.md` 第 28 节。

## 1. 背景和验证目标

现有 `bench_inverse_inequality` 在实值椭圆 LOD、拟一致粗网格上计算

```math
H_T\sup_{v_H\ne0}
\frac{\|\nabla(1-C)v_H\|_{L^2(T)}}
     {\|(1-C)v_H\|_{L^2(T)}}.
```

`DEVELOPMENT.md` 已完成以下校准：

1. `C=0` 粗 P1 基线对全局粗网格层数稳定；
2. 固定 `H,ell` 时必须先做嵌套 `h` 加细，避免把有限细空间常数误当成连续 corrector 常数；
3. 固定 `H,h` 时，局部逆常数对 `ell` 可能先增长后进入平台；
4. 广义特征值只能在局部质量矩阵的正子空间上求解；
5. 最大值才对应不等式常数，median、p90 和 p99 只能作诊断，不能替代最大值。

本计划把该问题扩展到 Helmholtz LOD 和非拟一致粗网格。粗网格在整个区域上保持协调、形状正则和 NVB 可达，但不再拟一致。每一步只标记当前最细叶单元中的一部分继续加细，使

```math
g_j=L_{\max,j}-L_{\min,j},
\qquad
\Gamma_j=\frac{H_{\max,j}}{H_{\min,j}}
```

随迭代 `j` 增大。核心问题是：使用每个单元自己的直径 `H_T` 缩放后，Helmholtz LOD 空间的局部逆常数是否仍与全局非拟一致性 `g_j`、`Gamma_j` 无关。

## 2. 数学模型和离散空间

### 2.1 Helmholtz 模型

在 `Omega=(0,1)^2` 上考虑

```math
-\nabla\cdot(A\nabla u)-k^2nu=f
\qquad\text{in }\Omega,
```

以及齐次阻抗 Robin 边界条件

```math
A\nabla u\cdot\nu-ik\beta u=0
\qquad\text{on }\partial\Omega.
```

采用项目现有约定：第一变量线性、第二变量共轭线性，

```math
a_k(u,v)
=(A\nabla u,\nabla v)_\Omega
-k^2(nu,v)_\Omega
-ik(\beta u,v)_{\partial\Omega}.
```

第一轮固定 `A=n=beta=1`，先隔离网格非拟一致性和 Helmholtz corrector 的影响；异质系数放在扩展实验中。

### 2.2 嵌套空间和 corrector

令 `T_H^j` 为第 `j` 个非拟一致 NVB 粗网格，`T_h` 为整个扫描期间固定的主细网格，要求

```math
V_H^j\subset V_h,
\qquad
I_H^j=E_H^j\Pi_{H,j}^{dg},
\qquad
W_h^j=\ker(I_H^j|_{V_h}).
```

对粗单元 `T` 和 `ell` 层 patch `omega_T^ell`，原 corrector 定义为

```math
a_{k,\omega_T^\ell}(Q_{T,h}^\ell v,w)
=a_{k,T}(v,w)
\qquad\forall w\in W_h^j(\omega_T^\ell).
```

伴随 corrector 按项目的第二变量约定定义为

```math
a_{k,\omega_T^\ell}(w,Q_{T,h}^{*,\ell}v)
=a_{k,T}(w,v)
\qquad\forall w\in W_h^j(\omega_T^\ell).
```

全局局部化 corrector 和两类基矩阵为

```math
Q_{H,h}^\ell=\sum_{T\in\mathcal T_H^j}Q_{T,h}^\ell,
```

```math
G_j^{trial}=(I-Q_{H,h}^\ell)P_H,
\qquad
G_j^{test}=(I-Q_{H,h}^{*,\ell})P_H.
```

当前实系数实现中 test basis 应与 trial basis 的逐项共轭一致，但实验仍分别计算两者，以免把实现特例误当成一般结论。

## 3. 需要验证的局部逆不等式

### 3.1 主判据：同单元局部逆不等式

对 `X in {trial,test}` 和任意 `w=G_j^X a`，主验证目标为

```math
\|\nabla w\|_{L^2(T)}
\le C_{inv}H_T^{-1}\|w\|_{L^2(T)}
\qquad\forall T\in\mathcal T_H^j.
```

相应的数值局部常数为

```math
Q_{T,j}^X
=H_T\sup_{a:\,\|G_j^Xa\|_T>0}
\frac{\|\nabla G_j^Xa\|_{L^2(T)}}
     {\|G_j^Xa\|_{L^2(T)}},
```

全局观测常数为

```math
C_{inv,j}^X=\max_{T\in\mathcal T_H^j}Q_{T,j}^X.
```

要检验的数值命题是：当 NVB 网格族形状正则、`kH_max` 处于 LOD 可解析区间、`h` 已解析 corrector 且 `ell` 足够时，存在与 `j`、`g_j` 和 `Gamma_j` 无关的有限常数，使

```math
\sup_j C_{inv,j}^{trial}<\infty,
\qquad
\sup_j C_{inv,j}^{test}<\infty.
```

### 3.2 广义特征值形式

在粗单元 `T` 的所有细子单元上装配几何刚度矩阵 `S_T` 和质量矩阵 `M_T`。对 `X in {trial,test}` 定义

```math
K_{T,j}^X=(G_j^X)^*S_TG_j^X,
\qquad
M_{T,j}^X=(G_j^X)^*M_TG_j^X.
```

则

```math
Q_{T,j}^X
=H_T\sqrt{\lambda_{\max}(K_{T,j}^X,M_{T,j}^X)}.
```

特征值问题只在 `M_{T,j}^X` 的正谱子空间上求解。默认丢弃

```math
\lambda_m\le10^{-12}\lambda_{m,\max}
```

的质量方向，并在小规模样本上用 `10^-10,10^-12,10^-14` 做阈值敏感性检查。

### 3.3 Patch 分母版本

同单元分母可能存在局部零迹方向。为区分真正的逆估计增长与质量子空间退化，同时报告一层邻域分母版本

```math
\widetilde Q_{T,j}^X
=H_T\sup_{a:\,\|G_j^Xa\|_{\omega_T^1}>0}
\frac{\|\nabla G_j^Xa\|_{L^2(T)}}
     {\|G_j^Xa\|_{L^2(\omega_T^1)}}.
```

该量使用 `T` 上的刚度矩阵和 `omega_T^1` 上的质量矩阵。两种定义必须分列保存，不能在看到结果后选择更平坦的一种作为“主结果”。

### 3.4 Helmholtz 能量版本

同时记录

```math
Q_{1,k,T,j}^X
=H_T\sup_a
\frac{\left(\|\nabla G_j^Xa\|_T^2
      +k^2\|G_j^Xa\|_T^2\right)^{1/2}}
     {\|G_j^Xa\|_T}.
```

在 `A=n=1` 且分子分母使用同一局部区域时，应满足核对恒等式

```math
(Q_{1,k,T,j}^X)^2=(Q_{T,j}^X)^2+(kH_T)^2.
```

该恒等式用于验证装配，不取代梯度型主判据。

### 3.5 辅助诊断量

以下量只用于定位增长来源：

```math
Q_{C,T,j}^{trial}
=H_T\sup_a
\frac{\|\nabla Q_{H,h}^\ell P_Ha\|_T}
     {\|(I-Q_{H,h}^\ell)P_Ha\|_T},
```

以及 `C=0` 粗 P1 基线。`Q_C` 不是主不等式，不得与 `Q_T` 混为一列。

## 4. 非拟一致 NVB 网格族

### 4.1 基本原则

这里的“全局非拟一致网格”指覆盖整个 `Omega` 的全局协调三角剖分，而不是每轮全局一致加细。网格必须满足：

- 从 `make_helmholtz_unit_square_mesh()` 的兼容参考边标号出发；
- 只调用现有 `bisect_newest_vertex`/`AdaptiveMeshHierarchy::refine`；
- 不允许悬挂点、重复坐标节点或手工拼接三角形；
- NVB 闭包可以附加加细，但必须单独记录闭包代价；
- 整个网格族共享一个固定主细网格 `T_h`，始终保持 `V_H^j subset V_h`。

### 4.2 确定性“最细叶继续加细”规则

初始粗网格取全局 NVB level `L_0`。第 `j` 轮计算当前最大叶层数

```math
L_{max,j}=\max_{T\in\mathcal T_H^j}L_H(T),
```

并只从集合

```math
F_j=\{T:L_H(T)=L_{max,j}\}
```

中选择用户标记单元。禁止标记较粗单元，以保证细化集中在已经最细的区域。

主网格族采用内部点 `x_*=(0.25,0.25)`，按单元重心到 `x_*` 的距离、稳定 element ID 依次排序，标记

```math
m_j=\max(1,\lceil\theta|F_j|\rceil),
\qquad \theta=0.25
```

个最细叶单元。该规则是确定性的，不使用 PDE 解或误差指标，避免把逆常数变化与自适应估计器混合。

另设两个稳健性网格族：

1. `single-chain`：每轮只标记距离 `x_*` 最近的一个最细叶单元，最大化全局尺度比；
2. `boundary-chain`：取 `x_*=(0,0.5)`，检验物理 Robin 边界和人工 patch 边界的共同影响。

### 4.3 非拟一致性指标

每轮必须记录：

```math
L_{min,j},\quad L_{max,j},\quad g_j=L_{max,j}-L_{min,j},
```

```math
H_{min,j}=\min_T H_T,\quad
H_{max,j}=\max_T H_T,\quad
\Gamma_j=H_{max,j}/H_{min,j},
```

以及共享边相邻单元的局部尺度比

```math
\Gamma_{nbr,j}
=\max_{T\sim T'}\max(H_T/H_{T'},H_{T'}/H_T).
```

NVB 预计使 `Gamma_nbr` 保持受控，而 `Gamma_j` 随深层局部加细增长。若 `Gamma_j` 没有总体增长，说明所生成的网格族没有实现本实验目标。

同时记录用户标记数、闭包附加加细数、最终真正加细数、最大/最小角和重复坐标节点数。

### 4.4 层号和几何尺度

项目中的 `H`、`h` 是 NVB sweep/层号，不是字面上的 `2^-H`、`2^-h`。所有表格以实测 `H_T=max edge length` 为准，层号仅用于重现 NVB 树。不得用层号替代几何直径计算 `Q_T` 或 `Gamma_j`。

## 5. 固定细空间和 patch 约定

### 5.1 固定主细网格

对一个完整的非拟一致扫描，预先固定最大局部加细轮数 `J_max` 和主细层 `L_h`，要求

```math
L_h\ge L_{max,J_{max}}+m_h,
```

其中第一轮取 `m_h=4`。每次粗网格改变后，通过现有 `complete_to_fine_level` 重建到同一 `L_h`，并核对规范化后的节点坐标和单元连接完全一致。

如果某轮粗单元达到 `L_h`，该轮无效；不能临时提高 `L_h` 后继续与前面数据混表。

### 5.2 局部解析比

记录每个粗单元的

```math
q_T=\max_{\tau\subset T}h_\tau/H_T,
\qquad q_{max}=\max_Tq_T.
```

主扫描固定 `T_h`，但随着局部粗网格变细，最细粗单元的 `q_T` 会增大。因此必须另做固定最终粗网格、继续减小 `h` 的扫描，确认 `C_inv` 对 `h` 已进入平台后，才能把随 `Gamma_j` 的增长归因于非拟一致粗网格。

### 5.3 Patch 层数

现有 patch 由粗网格拓扑层数定义。在非拟一致网格上，相同 `ell` 的物理直径和包含的细自由度可能差异很大。因此每个粗单元还要记录：

- patch 粗单元数和细单元数；
- `diam(omega_T^ell)/H_T`；
- patch 内最大与最小粗层数之差；
- 是否接触物理边界；
- 是否位于加细过渡带。

主扫描取固定 `ell=3`，之后用 `ell=2,3,4` 检查结论是否只是 patch 层数不足造成的。

## 6. 实现计划

### 6.1 新增 benchmark

建议新增

```text
benchmarks/bench_helmholtz_local_inverse.cpp
```

并链接 `lod2d_core`。不要直接修改实值椭圆 `bench_inverse_inequality` 的语义；新 benchmark 需要复值 Helmholtz basis、自适应粗网格、纯 Robin 边界和原/伴随 corrector，单独实现更容易审计。

建议命令行接口：

```text
--k=2
--initial-H=3
--fine-level=12
--steps=6
--ell=3
--mark=fraction|single-chain|boundary-chain
--mark-fraction=0.25
--seed-x=0.25 --seed-y=0.25
--basis=trial|test|coarse
--denominator=element|patch1
--mass-threshold=1e-12
--format=csv
--mesh-out=PATH
--check
```

一次运行应构造每个网格后同时计算 trial、test、coarse 三类结果，避免三次独立运行产生不同闭包路径。`--basis` 只用于调试单一分支。

### 6.2 复值局部谱计算

局部矩阵必须使用 `adjoint()`：

```math
K=G^*S_TG,\qquad M=G^*M_TG.
```

数值舍入后先检查

```math
\|K-K^*\|/\|K\|,\qquad \|M-M^*\|/\|M\|
```

低于容差，再用 `(K+K*)/2`、`(M+M*)/2` 进入自伴广义特征值流程。不得使用 `transpose()` 代替 `adjoint()`。

步骤为：

1. 对 `M` 做 `SelfAdjointEigenSolver`；
2. 丢弃相对阈值以下的质量特征向量；
3. 构造 `Z=U_+ Lambda_+^{-1/2}`；
4. 求 `Z^*KZ` 的最大特征值；
5. 回代最大化系数向量，检查广义特征残差；
6. 保存最大模态在细网格上的值，供热点可视化。

小网格必须用独立稠密 SVD/广义特征值计算交叉验证。

### 6.3 输出字段

CSV 至少包含：

```text
mesh_family,iteration,k,ell,basis,denominator,
coarse_nodes,coarse_elements,fine_nodes,fine_elements,
L_min,L_max,level_gap,H_min,H_max,grading_ratio,neighbor_ratio,
marked,closure_added,refined_total,q_max,
min,median,p90,p99,max,argmax_element,argmax_level,
argmax_x,argmax_y,argmax_boundary,argmax_transition,
patch_elements,patch_diameter_over_H,mass_rank,mass_condition,
petrov_residual,corrector_residual,constraint_residual,
mesh_ms,corrector_ms,inverse_ms,total_ms,peak_rss_kb
```

逐单元明细另存一个 CSV，不能只保存汇总分位数。网格输出包含稳定 element ID、parent ID、层数、三个顶点坐标和对应 `Q_T`。

### 6.4 运行脚本和目录

建议新增

```text
scripts/run_helmholtz_local_inverse.sh
results/helmholtz_local_inverse/
```

每个参数点使用独立进程，脚本保存命令、Git revision/status、编译器、线程数、主机内存、CSV、标准输出和 `/usr/bin/time -v`。支持 `RESUME=1`，但只有完整 CSV 行和成功退出码同时存在时才视为完成。

## 7. 数值实验步骤

### 步骤 1：网格生成器单元测试

使用小参数 `L_0=2,L_h=8,J_max=4`，分别生成 fraction、single-chain、boundary-chain 三个网格族。每轮验证：

- 面积总和为 1；
- 无悬挂点、无重复坐标节点、所有三角形面积为正；
- 层数与父子面积比一致；
- 稳定 element ID 和 parent ID 唯一；
- `g_j` 非减，`Gamma_j` 总体增加；
- 共享边只出现一次或两次；
- 固定细网格的规范化拓扑在各轮相同。

若这些检查失败，不进入逆不等式计算。

### 步骤 2：`C=0` 校准

在每个非拟一致网格上令 `G=P_H`，计算 element 和 patch1 两种分母的 `Q_T`。预期：

- 使用局部 `H_T` 时，`max Q_T` 对 `Gamma_j` 保持稳定；
- 使用错误的全局 `H_max` 缩放时，最细区常数随 `Gamma_j` 增长；
- 最大模态和直接函数积分一致；
- trial/test 分支在 `C=0` 时完全相同。

这一负控制能检验局部直径、粗细子单元映射和广义特征值实现。若 `C=0` 已随 `Gamma_j` 明显增长，优先排查实现，不能据此否定 LOD 不等式。

### 步骤 3：小规模 Helmholtz 正确性门槛

取 `k=2,L_0=2,L_h=8,ell=2,J_max=3`：

1. 检查 `I_HP_H=I`；
2. 检查原/伴随 corrector 属于 `ker I_H`；
3. 检查 corrector、约束和 Petrov 残差；
4. 检查物理边界保留 Robin 项、人工 patch 边界施加齐次 Dirichlet；
5. 对实系数检查 `G_test-conj(G_trial)`；
6. 将若干 `Q_T` 与独立稠密计算比较；
7. 核对 Helmholtz 能量版本恒等式。

建议门槛：结构等式误差 `<1e-10`，线性求解残差 `<1e-9`，局部谱相对差 `<1e-9`。

### 步骤 4：WSL 主非拟一致性扫描

第一轮可运行规模：

```text
k=2, L_0=3, L_h=12, ell=3,
J=0,...,6, mark=fraction, theta=0.25.
```

每轮同时输出 coarse、trial、test；同时计算 element 和 patch1 分母。重点绘制或制表：

```math
C_{inv,j}^X\ \text{vs.}\ g_j,
\qquad
C_{inv,j}^X\ \text{vs.}\ \Gamma_j,
```

并按单元层数、加细核心、过渡带、未加细远场和物理边界分组报告 max/p99。正式结论必须使用 max。

### 步骤 5：标记规则和位置稳健性

固定 `k,L_0,L_h,ell,J_max`，比较：

- `fraction, theta=0.25`；
- `single-chain`；
- `boundary-chain`；
- fraction 的 `theta=0.10,0.50`。

若只有某一种规则增长，应定位 argmax 是否始终位于闭包过渡带。边界族增长而内部族稳定时，应继续区分 Robin 边界效应与 patch 几何效应。

### 步骤 6：固定最终粗网格的 `h` 极限

选取主扫描的 `j=0,3,6` 三个粗网格，保持粗网格和 `ell` 不变，取

```math
L_h=L_{max,j}+2,\ L_{max,j}+4,\ L_{max,j}+6.
```

检查 `C_inv` 是否随 `h` 进入平台。若常数仍持续随细网格增长，则当前数据只能说明离散空间行为，不能判断连续 corrector 的局部逆不等式。

### 步骤 7：oversampling 扫描

固定一个中等和一个强非拟一致网格，取 `ell=1,2,3,4`。先识别 `ell` 平台，再比较不同 `Gamma_j`。如果增大 `ell` 消除了随 grading 的增长，说明原结果更可能是局部化误差，而不是非拟一致网格本身导致的失效。

### 步骤 8：波数扫描

使用 `k=1,2,4,8`，为每个 `k` 选择初始层数使

```math
kH_{max,0}\le1,
```

并保持相近的 `kh` 和相同的 `L_h-L_{max,J}`。报告原始 `C_inv`、`kH_T` 和能量版本常数。任何违反解析条件的点必须标记为 `underresolved`，不得并入关于网格无关性的拟合。

### 步骤 9：扩展系数实验

只有单位系数主扫描通过后才增加：

- `A` 的固定棋盘系数，对比度 `10^2,10^4`；
- `n` 和 `beta` 的有界正分片常数；
- 加细核心与系数跳跃重合/不重合两种位置。

此阶段同时报告几何梯度版本和 `A` 加权版本。不同系数族不共享同一个拟合常数。

## 8. 预注册判据和结论规则

### 8.1 支持性证据

对每个固定的 `k,ell,h` 解析族，拟合

```math
\log C_{inv,j}=b+\alpha\log\Gamma_j.
```

在至少四个不同 `g_j` 且 `Gamma_j` 增长至少 4 倍时，以下条件同时满足才写为“数值上支持 grading-independent”：

1. 最后三个已解析点的 `max C_inv` 最大值/最小值 `<=1.15`；
2. 拟合斜率 `|alpha|<=0.10`；
3. element 和 patch1 两种分母给出一致的稳定/增长判断；
4. trial 和 test 均稳定；
5. `C=0` 基线稳定；
6. 增加 `h` 和 `ell` 不改变判断；
7. 所有正确性残差通过门槛。

这些阈值是实验分类标准，不是理论常数。

### 8.2 反证信号

以下任一现象应报告为“观察到不稳定或可能反例”，而不是调参隐藏：

- 最后三个已解析点连续增长且总增长超过 25%；
- `alpha>=0.20`，并在更细 `h`、更大 `ell` 后仍存在；
- 最大化模态稳定地集中在加细过渡带，且增长与 `Gamma_j` 同步；
- 同单元和 patch1 分母均增长；
- coarse `C=0` 稳定而 trial/test LOD 不稳定。

若只在质量阈值改变时增长，结论应为“局部质量退化/谱计算不稳”；若只在未解析的 `h` 或 `kH` 下增长，结论应为“参数未解析”，不能作为反例。

### 8.3 不允许的结论

- 不能用 median 或 p99 平台代替 `max` 平台；
- 不能把有限组数据称为数学证明；
- 不能把 NVB level 直接写成几何 `2^-level`；
- 不能混合不同 `h`、`ell`、标记规则或波数后只画一条无标签曲线；
- 不能删除增长点后重新选择阈值；
- 不能把 `Q_C` 当作主局部逆常数。

## 9. 资源分级和建议参数

### 9.1 本地 WSL 校准

```text
P1: k=2, L_0=2, L_h=8,  J=0:3, ell=2
P2: k=2, L_0=3, L_h=12, J=0:6, ell=3
```

先串行或 4 线程验证确定性，再用 8 线程计时。每个网格单独进程，避免保留前一轮 factorization。

### 9.2 高内存正式扫描

```text
k=2,4,8
ell=2,3,4
mark=fraction,single-chain,boundary-chain
J=0:8
L_h-L_max,J=4 and 6
```

正式扫描前先估计每个 patch 的细自由度和内存，不直接启动最大笛卡尔积。优先完成单位系数、内部 fraction 族，再扩展其他维度。

## 10. 测试和验收清单

### 网格与嵌套

- [ ] NVB 网格协调、面积守恒、无重复坐标节点；
- [ ] `V_H^j subset V_h` 和 `I_HP_H=I`；
- [ ] 固定主细网格跨迭代不变；
- [ ] `g_j`、`Gamma_j`、`Gamma_nbr,j` 正确记录；
- [ ] 用户标记和闭包附加加细分开记录。

### 局部矩阵和谱

- [ ] `S_T`、`M_T` 的子单元装配无遗漏无重复；
- [ ] 复基矩阵使用 `adjoint()`；
- [ ] Hermitian 缺陷、质量秩和条件数已输出；
- [ ] 正质量子空间阈值敏感性通过；
- [ ] 最大特征对残差和稠密交叉验证通过；
- [ ] 能量版本恒等式通过。

### Helmholtz LOD

- [ ] 原/伴随 corrector 残差和约束残差通过；
- [ ] 物理 Robin 边界和人工 patch 边界处理正确；
- [ ] trial/test 分列输出；
- [ ] `kH_max`、`kh` 和 inf-sup/稳定性诊断已记录；
- [ ] 固定粗网格的 `h` 极限已检查。

### 可重复性

- [ ] CSV 包含汇总和逐单元明细；
- [ ] 保存网格、argmax 模态、命令、Git 状态和资源统计；
- [ ] 脚本支持中断续跑且不重复成功点；
- [ ] 关键结果写入 `DEVELOPMENT.md`，稳定用法写入指南；
- [ ] WSL 与 Windows 镜像中的源文档保持一致。

## 11. 预计文件改动

实施阶段预计涉及：

```text
benchmarks/bench_helmholtz_local_inverse.cpp
benchmarks/CMakeLists.txt
scripts/run_helmholtz_local_inverse.sh
tests/test_helmholtz_local_inverse.cpp
tests/CMakeLists.txt
HELMHOLTZ_GUIDE.md
DEVELOPMENT.md
results/helmholtz_local_inverse/
```

如果局部谱装配在 benchmark 与测试之间出现重复，再提取到

```text
include/helmholtz/local_inverse.h
src/helmholtz/local_inverse.cpp
```

第一版不应为追求抽象而提前改动公共 API。

## 12. 最终交付格式

最终实验报告至少回答：

1. `C=0` 粗 P1 基线在同一非拟一致网格族上是否稳定；
2. Helmholtz trial/test 的 `max Q_T` 是否随 `Gamma_j` 增长；
3. argmax 位于加细核心、过渡带、远场还是 Robin 边界；
4. 结论对同单元/patch 分母、`h`、`ell`、波数和标记规则是否稳健；
5. 观察到的平台是否可能只是有限细空间或局部化平台；
6. 哪些参数点满足解析条件，哪些只能作为失稳诊断；
7. 数据支持、否定或尚不能判断哪一种数学命题。

只有这些问题都有原始数据支持后，才能在 `DEVELOPMENT.md` 中写出关于 Helmholtz 非拟一致 NVB 网格局部逆不等式的结论。
