# Helmholtz-LOD 实现任务书

## 执行状态

> 状态日期：2026-07-21
>
> 七个 Helmholtz 基础阶段均已完成。Patch 求解器优化、两层 Schwarz 和
> 自适应 LOD 分别由独立计划书继续跟踪。

| 阶段 | 范围 | 状态 |
|---|---|---|
| 1 | 复值 FEM 算子与制造解收敛 | 已完成 |
| 2 | 全局 NVB 层次与复值拟插值 | 已完成 |
| 3 | 原/伴随局部矫正子及 patch 边界 | 已完成 |
| 4 | 双边及 test-only Petrov-Galerkin 模型 | 已完成 |
| 5 | 波数扫描和 FEM/LOD/reference 对照 | 已完成 |
| 6 | 矫正子并行化与可复用算子数据 | 已完成 |
| 7 | 用户指南、测试和可复现扫描脚本 | 已完成 |

当前使用方法见 [HELMHOLTZ_GUIDE.md](HELMHOLTZ_GUIDE.md)，实测验证数据
统一保存在 [DEVELOPMENT.md](DEVELOPMENT.md)。

## 1. 项目目标

在现有二维 P1 LOD C++ 项目中实现 Peterseim 型 Helmholtz Petrov-Galerkin LOD，主要数值目标是验证误差对波数 `k` 的鲁棒性，并与标准 P1 Galerkin FEM 比较其 pollution effect。

第一阶段包含：

- 复值 Helmholtz 有限元与纯齐次阻抗型 Robin 边界条件；
- 基于全局 newest-vertex bisection（NVB）的嵌套粗细网格；
- 复值空间上的拟插值 `I_H = E_H * Pi_H^dg`；
- 原矫正子、伴随矫正子及双边 Petrov-Galerkin 离散；
- 标准 FEM、LOD 和细网格参考解的波数扫描实验；
- 同一算子、不同右端项时的 corrector 缓存复用。

局部自适应和 PML 不在第一阶段实现，但数据结构和接口必须为二者保留扩展空间。

主要理论参考：Daniel Peterseim, *Eliminating the pollution effect in Helmholtz problems by local subscale correction*, Mathematics of Computation 86 (2017), 1005-1036；公开版本：<https://arxiv.org/abs/1411.7512>。

## 2. 数学模型与统一约定

### 2.1 模型问题

先取 `Omega = (0,1)^2`，并令 `A = n = beta = 1`。求复值函数

```math
u \in V := H^1(\Omega;\mathbb C)
```

满足

```math
-\nabla\cdot(A\nabla u)-k^2 n u=f \qquad\text{in }\Omega,
```

```math
A\nabla u\cdot\nu-i k\beta u=0 \qquad\text{on }\partial\Omega.
```

这里的“纯齐次 Robin”表示整个 `partial Omega` 都采用上述阻抗边界，边界右端项为零；体力 `f` 可以非零。

### 2.2 复值内积和双线性型

统一采用“第一变量线性、第二变量共轭线性”的约定：

```math
a_k(u,v)=
\int_\Omega A\nabla u\cdot\nabla\overline v\,dx
-k^2\int_\Omega n u\overline v\,dx
-ik\int_{\partial\Omega}\beta u\overline v\,ds.
```

右端项为

```math
F(v)=\int_\Omega f\overline v\,dx.
```

细网格矩阵因此为

```math
\mathbf A_h(k)=\mathbf K_h-k^2\mathbf M_h-ik\mathbf B_h,
```

其中 `B_h` 是物理边界质量矩阵。当前实系数情形下 `A_h` 是复对称矩阵，但不是 Hermitian 矩阵。Helmholtz 主系统和约束鞍点系统不得使用 LLT、CG 或 PCG。

误差采用 Helmholtz 能量范数

```math
\|v\|_{1,k}^2=\|\nabla v\|_{L^2(\Omega)}^2+k^2\|v\|_{L^2(\Omega)}^2.
```

## 3. 网格层次

### 3.1 全局 NVB

从具有兼容参考边标记的初始三角网格出发，每一级标记全部单元，再执行 NVB 和必要的 closure，使网格始终协调且无悬挂点：

```math
\mathcal T_H \subset \mathcal T_h.
```

`H` 和 `h` 均使用对应网格的最大单元直径，而不能仅使用加细层编号代替。

### 3.2 必须保留的层次信息

- coarse vertex 到 fine vertex 的嵌入关系；
- fine element 的 coarse ancestor；
- coarse element 的 fine descendants；
- 元素邻接、边邻接和 patch 层数；
- 每条 patch 边属于物理边界还是人工截断边界；
- NVB refinement tree，供后续自适应加细和缓存失效分析使用。

### 3.3 网格测试

- 每轮全局 NVB 后无悬挂点；
- 面积守恒，三角形方向一致；
- 粗细 P1 嵌入在随机点上精确一致；
- ancestor/descendant 覆盖完整且不重复；
- 重复全局 NVB 的单元数和顶点数与 golden data 一致。

## 4. 复值拟插值

定义

```math
I_H=E_H\circ\Pi_H^{dg}.
```

其中 `Pi_H^dg` 是逐粗单元的 `L2` 投影：

```math
(\Pi_H^{dg}v,q_H)_T=(v,q_H)_T
\qquad\forall q_H\in P_1(T),
```

`E_H` 再将不连续的粗 P1 函数通过顶点加权平均映射到连续空间 `V_H`。

由于网格、基函数和平均权重为实数，`I_H` 的矩阵可以继续以实稀疏矩阵保存，并复线性地作用于复向量。但是所有通用复值代码必须使用共轭转置 `adjoint()`，不能依赖 `transpose()` 在当前实矩阵上的偶然等价。

定义细尺度空间

```math
W_h=\ker I_H.
```

拟插值验收条件：

- `I_H v_H = v_H`；
- `I_H(I_H v) = I_H v`；
- 对构造出的 `w_h in W_h`，`||I_H w_h||` 接近机器精度；
- 对复数 `alpha, beta` 验证复线性；
- 实部、虚部分别计算与一次复值计算结果一致。

## 5. 局部 patch 与边界条件

对每个粗单元 `T` 构造 `ell` 层 patch `omega_T^ell`，并定义

```math
W_h(\omega_T^\ell)=
\{w_h\in W_h:\operatorname{supp}w_h
\subset\overline{\omega_T^\ell}\}.
```

patch 边界必须分成两部分：

```math
\Gamma_{\rm art}=\partial\omega_T^\ell\setminus\partial\Omega,
\qquad
\Gamma_{\rm phys}=\partial\omega_T^\ell\cap\partial\Omega.
```

- 在 `Gamma_art` 上施加齐次 Dirichlet 条件；
- 在 `Gamma_phys` 上不消去自由度，并保留原问题的 Robin 边界项；
- patch 算子只装配 `Gamma_phys` 上的边界质量矩阵；
- element corrector 右端中的 `a_{k,T}` 只包含 `partial T cap partial Omega` 上的 Robin 边积分。

这一分类应由边拓扑完成，禁止通过浮点坐标比较推断边界。

## 6. 原矫正子与伴随矫正子

### 6.1 原矫正子

对 `v_H in V_H`，局部原矫正子满足

```math
a_{k,\omega_T^\ell}(C_{T,\ell}v_H,w)
=a_{k,T}(v_H,w)
\qquad\forall w\in W_h(\omega_T^\ell).
```

定义

```math
C_\ell=\sum_{T\in\mathcal T_H}C_{T,\ell}.
```

### 6.2 伴随矫正子

伴随矫正子满足

```math
a_{k,\omega_T^\ell}(w,C^*_{T,\ell}v_H)
=a_{k,T}(w,v_H)
\qquad\forall w\in W_h(\omega_T^\ell).
```

对于当前的实系数、实拟插值和复对称矩阵，存在快捷关系

```math
C_\ell^*v=\overline{C_\ell\overline v}.
```

第一版可以只计算原矫正子并通过该关系生成检验基，但必须保留显式 adjoint API 和残差测试。将来引入一般复系数或 PML 后，不允许默认该关系仍成立。

### 6.3 求解策略

正确性基线使用复值稀疏直接法，例如 `Eigen::SparseLU`。局部问题虽然在满足 `kH` 分辨率条件的核空间上具有强制性，但原始矩阵和约束鞍点矩阵均不适合 LLT/CG/PCG。

直接法通过全部正确性测试后，再评估：

- 以 Laplace constrained operator 为预条件器的 GMRES；
- patch pattern 分组与 symbolic factorization 复用；
- 原/伴随矫正子的共轭复用；
- 同一 patch 中三个局部粗基右端的 block solve。

## 7. Petrov-Galerkin LOD

令 `P` 为粗 P1 系数到细网格 P1 系数的嵌入矩阵，定义试验基和检验基

```math
\mathbf\Psi_{\rm tr}=(I-C_\ell)\mathbf P,
\qquad
\mathbf\Psi_{\rm te}=(I-C_\ell^*)\mathbf P.
```

完整双边 Petrov-Galerkin 系统为

```math
\mathbf\Psi_{\rm te}^*\mathbf A_h
\mathbf\Psi_{\rm tr}\mathbf U_H
=\mathbf\Psi_{\rm te}^*\mathbf F_h.
```

细网格上的重构解为

```math
\mathbf u_{H,\ell,h}=\mathbf\Psi_{\rm tr}\mathbf U_H.
```

同时实现论文讨论的一边修正模式：试验空间使用 `V_H`，检验空间使用 `(I-C_ell^*)V_H`。完整双边模式作为数学正确性金标准，一边修正模式作为独立的性能和精度实验项。任何模式都不得被命名为普通 Galerkin。

## 8. 软件模块设计

在现有 `LodProblemData`、`LodOperators` 和 `LodModel` 分层上增加 Helmholtz 专用类型，避免把复值、波数和边界语义塞入椭圆代码的隐式分支。

建议职责如下：

- `HelmholtzProblemData`：网格、`k`、系数 `A/n/beta`、右端项和边界类型；
- `HelmholtzOperators`：`K`、`M`、`B`、`A(k)`、load vector 和误差范数；
- `HelmholtzCorrector`：patch、原/伴随矫正子、solver 和缓存；
- `HelmholtzLodModel`：预计算、Petrov-Galerkin 组装、求解和重构；
- benchmark 层：参数扫描、计时、CSV/JSON 输出，不重复数学装配代码。

底层网格、patch、NVB、稀疏装配和参数解析工具应与椭圆模块共享。复值标量通过模板或明确的复值类型别名传播，不能把复向量拆成两套重复实现。

## 9. 正确性验证矩阵

### 9.1 算子级测试

- 单个三角形的 `K_T`、`M_T`；
- 单条边的 Robin 质量矩阵；
- `A_h = K_h - k^2 M_h - i k B_h` 的符号；
- `x^* A_h y` 与直接数值积分一致；
- 每条物理边恰好装配一次。

### 9.2 细网格 FEM

构造满足齐次 Robin 条件的光滑复值制造解，验证固定 `k` 下预期收敛阶：

```math
\|u-u_h\|_{1,k}=O(h),
\qquad
\|u-u_h\|_{L^2}=O(h^2).
```

制造解主要验证符号、共轭和边界装配；波数鲁棒性实验另用固定体源与细网格参考解。

### 9.3 corrector 测试

对随机复值粗函数和随机核空间函数验证

```math
a(C_\ell v,w)-a(v,w)\approx0,
```

```math
a(w,C_\ell^*v)-a(w,v)\approx0,
```

以及

```math
C_\ell^*v\approx\overline{C_\ell\overline v}.
```

额外验证人工 patch 边界值为零、物理边界自由度未被误删，并比较小网格局部实现与显式稠密 `ker(I_H)` 参考实现。

### 9.4 Petrov-Galerkin 测试

- 显式基矩阵组装与 element-wise 组装一致；
- coarse residual 在检验空间中接近机器精度；
- 双边模式和一边模式分别与小网格稠密参考系统一致；
- 改变 `f` 但保持算子不变时，corrector cache 命中且结果与重算一致。

## 10. 波数鲁棒性实验

### 10.1 参数设计

第一组建议取

```math
k\in\{4,8,16,32,64\}.
```

利用 NVB 层级选择粗网格，使

```math
kH\approx0.5\quad\text{或}\quad kH\approx1,
```

并选择 corrector 细网格使

```math
kh\le0.05.
```

默认 oversampling 取

```math
\ell=\lceil\log_2 k\rceil.
```

在最大 `k` 上至少增加一级细网格，确认参考解和 fully discrete corrector 的细尺度误差已经稳定。

### 10.2 对比方法

- coarse P1 Galerkin FEM；
- 双边 Petrov-Galerkin LOD；
- 一边修正 Petrov-Galerkin LOD；
- fine P1 参考解。

### 10.3 误差和稳定性指标

```math
e_{1,k}=\frac{\|u_h-u_H\|_{1,k}}{\|u_h\|_{1,k}},
\qquad
e_{L^2}=\frac{\|u_h-u_H\|_{L^2}}{\|u_h\|_{L^2}}.
```

还应输出：

- `kH`、`kh`、每波长自由度；
- `||u_h-u_H||_(1,k) / (H ||f||)`；
- coarse scaled operator 的最小奇异值估计；
- mesh、operator、corrector、coarse solve、reference solve 和总时间；
- 峰值内存、solver 迭代次数和失败状态。

固定每波长自由度时，若标准 FEM 的相对误差随 `k` 明显增长，而 LOD 误差保持有界，即观察到预期的抗 pollution 行为。不能只比较某一个 `k` 下的绝对误差。

### 10.4 oversampling 实验

对每个 `k` 比较

```math
\ell=1,
\quad
\left\lceil\tfrac12\log_2 k\right\rceil,
\quad
\lceil\log_2 k\rceil,
\quad
\left\lceil\tfrac32\log_2 k\right\rceil.
```

该实验用于区分离散误差、pollution error 和 localization error，并检查理论上的对数 oversampling 趋势。

## 11. 缓存与不同右端项复用

corrector 只依赖算子和空间，不依赖体力 `f`。缓存键至少包含：

- coarse/fine mesh 拓扑及几何 hash；
- NVB 层次关系；
- `k`；
- `A`、`n`、`beta` 及其离散数据 hash；
- Robin/Dirichlet 边界标记；
- `ell`；
- 拟插值类型和参数；
- corrector solver 及影响结果的容差。

相同 `A / n / beta / mesh / H / h / ell / k`、不同 `f` 时，应完整复用原矫正子、伴随矫正子和 coarse operator factorization。改变 `k` 后默认不得复用 corrector；多波数降阶或插值复用属于后续独立研究任务。

## 12. 自适应扩展路线

第一阶段即保留下列数据，以免未来重构：

- coarse element residual 和 interior-edge jump；
- 每个 element corrector 的局部贡献、patch 范围和 cache key；
- localization indicator 与 fine discretization indicator；
- NVB refinement tree 和受影响 patch 的反向索引。

未来自适应循环为

```text
SOLVE -> ESTIMATE -> MARK -> NVB REFINE -> UPDATE PATCHES -> SOLVE
```

局部加细后只使几何上受影响的 corrector cache 失效。误差估计需要同时区分 coarse residual、corrector localization 和细网格离散误差，不能只套用椭圆问题的残差估计器。

## 13. PML 扩展路线

PML 会把标量 Laplace 项改写成复值、各向异性的张量扩散项，并改变质量系数。为此：

- operator 接口应接受复值张量 `A(x)` 和复值质量系数 `n(x)`；
- 边界条件由策略对象管理，不能在 corrector 中硬编码 Robin；
- 原矫正子和伴随矫正子必须允许真正独立装配；
- `C^*v = conj(C conj(v))` 只能作为经过条件检查的 fast path；
- physical domain、PML layer 和 outer boundary 必须有独立 region/boundary tags。

PML 实现应在纯 Robin 版本的波数鲁棒性验证完成后开始。

## 14. 实施顺序与验收门槛

1. **复值 FEM 基线**：完成 `K/M/B/A(k)` 和制造解收敛测试。
2. **NVB 与复值拟插值**：通过嵌套、投影和核空间测试。
3. **局部 corrector**：直接法通过 primal、adjoint 和 patch 边界测试。
4. **Petrov-Galerkin 模型**：双边模式通过小网格稠密参考测试。
5. **波数实验**：完成标准 FEM/LOD/reference 的参数扫描。
6. **性能优化**：在误差不变的前提下加入 pattern、factorization 和多右端复用。
7. **文档与复现**：README 保留快速开始；HELMHOLTZ_GUIDE 维护完整运行和服务器扫描命令。

每一阶段必须保持现有椭圆 LOD 测试通过。任何求解器优化都必须同时报告：误差差异、corrector residual、总时间和峰值内存。

## 15. 第一阶段完成标准

- 标准 Helmholtz FEM 对制造解达到预期收敛阶；
- primal/adjoint corrector 残差达到设定容差；
- patch 物理边界使用 Robin，人工边界使用齐次 Dirichlet；
- 双边 Petrov-Galerkin 与小网格稠密金标准一致；
- 固定 `kH` 时形成 `k = 4...64` 的可复现实验表格和误差曲线；
- 能清楚区分 FEM pollution、localization error 和 fine-grid error；
- 相同算子、不同 `f` 时无需重算 corrector；
- 现有椭圆问题 benchmark 和测试无回归。

实现原则是先建立可信的复值直接法金标准，再进行迭代求解和性能优化。Helmholtz 代码中最优先检查的是共轭转置、Robin 符号、物理/人工 patch 边界分类以及 primal/adjoint 空间是否真正不同。
