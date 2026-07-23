# Helmholtz LOD 连续 hp 细尺度空间与 H 收敛实验计划书

> 状态日期：2026-07-24
>
> 当前状态：HP4 已完成；HP5-HP6 已完成本地正确性与缩减规模复核，`L_h=12` 全矩阵扫描待服务器运行。
>
> 实验范围：单位正方形、拟一致全局 NVB 网格、纯齐次阻抗 Robin
> 边界、粗空间连续 P1、细尺度空间连续 Pp、双边
> Petrov-Galerkin Helmholtz LOD。

## 执行状态

| 阶段 | 内容 | 状态 | 退出条件 |
|---|---|---|---|
| HP0 | 数学口径、误差分解和实验矩阵冻结 | 完成 | 本计划书通过评审 |
| HP1 | 连续三角形 Pp 空间、自由度和积分规则 | 完成 | P1 等价、连续性和积分测试通过 |
| HP2 | `I_H=E_H Pi_H^dg` 对连续 Pp 空间的扩展 | 完成 | `I_H P_H=I` 和核约束测试通过 |
| HP3 | 复值 hp Helmholtz 算子与 patch 鞍点问题 | 完成 | 原/伴随残差和边界测试通过 |
| HP4 | hp Petrov-Galerkin LOD 高级模型 | 完成 | p=1 与现有模型一致 |
| HP5 | fine-hp、局部化和 H 收敛校准 | 完成本地规模 | 三类误差已分别输出 |
| HP6 | 主收敛扫描、原始数据和结论 | 部分完成 | 本地复核完成；服务器全矩阵待运行 |

### 0.1 第三步实施记录

2026-07-24 已按实施顺序完成 HP0-HP3，并在第三步边界停止：

- 连续三角形 `P_p` 首版支持 `p=1,2,3`，自由度由顶点、定向共享边
  和单元内部自由度组成；
- 三角形与边积分使用按次数生成的 Gauss 规则，Robin 项只装配在
  物理边界；
- 已实现 coarse P1 到 fine `P_p` 的节点注入和混合矩形式的
  `I_H=E_H Pi_H^dg`，构造阶段检查 `I_H P_H=I`；
- 单 patch 组装以 hp 自由度全局 incidence 删除人工边界自由度，
  保留与物理边界重合部分的高阶 Robin 自由度；
- 原 corrector 使用现有 `DirectSaddle` 金标准，伴随 corrector 在
  当前实约束/实 coarse basis 约定下由共轭关系得到，并显式检查
  `A^*q^*+B^Tlambda^*=\overline r`；
- 内部、普通边界和角点三类 patch 的原问题、核约束和伴随残差均
  小于 `2e-10`；
- 制造解小网格回归观测到 `p=1,2,3` 的 energy 阶约为
  `0.74,1.97,2.74`；该测试只作为实现回归，正式渐近阶结论留给
  HP5-HP6。

### 0.2 第四至第六步实施记录

2026-07-24 已实现独立的 `HelmholtzHpLodModel`。模型组装全部 element
corrector，构造双边 trial/test 基、Petrov-Galerkin 粗系统，并复用
corrector、粗分解和 fine 分解处理不同右端项。`p=1` 的注入、拟插值、
fine 算子、trial/test 基、粗算子和最终解均与现有 P1 模型在
`3e-9` 以内一致。

统一 benchmark 支持 `fem`、`ell`、固定 master fine 的 `H` 扫描和
耦合 `H-h` 扫描。本地已完成 fine-hp、`ell` 校准、固定 fine 扫描、
`L_h=10` master-depth 复核和耦合复核。完整参数和结论见
`DEVELOPMENT.md`。

本地内存和 direct-saddle 成本不适合执行预注册的全部 `L_h=12`
patch 矩阵。`scripts/run_helmholtz_hp_convergence_server.sh` 使用临时 CSV
和成功 `.done` 标记；设置 `FULL=1` 执行原计划的服务器主扫描。当前
不能把本地缩减扫描写成已通过全部预注册渐近阶门槛。

本计划只研究连续细尺度空间次数 `p` 对 corrector 离散质量和 LOD
关于 `H` 的收敛表现。它不实现局部 `p_T` 自适应，也不把粗尺度约束
空间提高到高阶。

## 1. 研究问题

当前 Helmholtz LOD 使用连续 P1 细空间。新实验将细空间替换为

```math
V_h^p
=
\left\{
v\in H^1(\Omega;\mathbb C):
v|_\tau\in\mathbb P_p(\tau)
\quad\forall\tau\in\mathcal T_h
\right\},
```

并定义

```math
W_{H,h}^p
:=
\ker\left(I_H|_{V_h^p}\right).
```

局部 corrector 在

```math
W_{H,h}^{p,\ell}(T)
:=
\left\{
w\in W_{H,h}^p:
\operatorname{supp}w
\subset\overline{N^\ell(T)}
\right\}
```

中计算。实验回答以下问题：

1. 连续 Pp 离散的原、伴随 corrector 是否满足现有复值约束和
   Petrov-Galerkin 约定？
2. 增加 `p` 是否降低 fine-space/corrector 离散误差和误差平台？
3. 当粗空间仍为 P1 时，LOD 能量误差关于 `H` 是否仍表现为一阶，
   而不是随 `p` 自动升阶？
4. `p` 是否影响 corrector 的局部化速度，从而要求不同的 `ell`？
5. 在 fine-hp 误差和局部化误差都足够小时，测得的 `H` 收敛阶是否
   稳定、可复现？

### 1.1 这里的 p 不是什么

本计划中的 `p` 只表示细空间 `V_h^p` 的多项式次数。粗空间和粗尺度
约束仍为连续 P1。因此它不同于高阶 LOD 文献中提高 coarse
constraint/QOI 次数的 `p`。

预期结论是：

```math
\boxed{
\text{提高细空间 }p
\text{ 可以降低 corrector 离散误差，但不自动提高关于 }H
\text{ 的主收敛阶。}
}
```

若未来希望获得类似 `O(H^{p+2})` 的高阶主误差，必须另行提高粗约束
空间或 coarse QOI 的次数；这不属于本实验。

## 2. Helmholtz 模型

取

```math
\Omega=(0,1)^2,
\qquad A=n=\beta=1,
```

并求 `u in H^1(Omega;C)` 满足

```math
-\Delta u-k^2u=f
\qquad\text{in }\Omega,
```

```math
\partial_nu-\mathrm i k u=0
\qquad\text{on }\partial\Omega.
```

整个物理边界都使用齐次阻抗 Robin 条件，不设置 Dirichlet 边界。

采用“第一变量线性、第二变量共轭线性”的约定：

```math
a_k(u,v)
=
(\nabla u,\nabla v)_\Omega
-k^2(u,v)_\Omega
-\mathrm i k(u,v)_{\partial\Omega}.
```

能量型误差范数取为

```math
\|v\|_{1,k}^2
:=
\|\nabla v\|_{L^2(\Omega)}^2
+k^2\|v\|_{L^2(\Omega)}^2.
```

边界阻抗项不直接加入正定误差范数，但必须进入离散算子和代数残差。

## 3. 拟一致全局 NVB 网格

从当前单位正方形初始三角剖分出发，对所有活动单元执行全局
newest-vertex bisection：

```math
\mathcal T_0
\prec\mathcal T_1
\prec\cdots
\prec\mathcal T_{L_h}.
```

每一级都要求：

- 网格协调，无悬挂点；
- 形状正则常数受控；
- 所有单元来自同一 NVB 树；
- 粗网格和细网格嵌套；
- `P_node`、`P_elem` 和 `P_dg` prolongation 恒等式成立。

主实验只使用全局加细，因此网格是拟一致的。非拟一致网格和局部
`p_T` 分布留给后续联合自适应实验。

### 3.1 网格尺度

不能用 NVB sweep 数直接代替网格尺寸。每一行必须测量

```math
H
:=
\max_{T\in\mathcal T_H}\operatorname{diam}(T),
\qquad
h
:=
\max_{\tau\in\mathcal T_h}\operatorname{diam}(\tau).
```

由于一次全局 NVB sweep 只使单元数约翻倍，建议主扫描每隔两个 sweep
取一个粗层级，使 `H` 明显下降。收敛阶始终使用实测 `H` 计算。

## 4. 连续三角形 Pp 空间

### 4.1 第一版范围

第一版支持统一次数

```math
p\in\{1,2,3\},
```

`p=4` 作为可选服务器实验。所有细单元使用同一个 `p`，因此共享边上的
迹空间一致，连续性定义没有歧义。

第一版不实现任意局部 `p_tau`。未来若允许相邻单元次数不同，必须采用
层次基和 minimum rule，不能简单拼接两个不同次数的 nodal trace。

### 4.2 自由度布局

对三角形 `tau`，连续 Pp 自由度分为：

- 3 个顶点自由度；
- 每条边 `p-1` 个内部自由度；
- 每个单元 `(p-1)(p-2)/2` 个内部自由度。

全局自由度数满足

```math
N_{h,p}
=
N_v
+(p-1)N_e
+\frac{(p-1)(p-2)}2N_T.
```

边自由度必须使用全局规范边方向。局部三角形方向改变时，边内部基函数
的局部编号必须正确置换，避免跨单元迹不连续。

### 4.3 基函数选择

正确性阶段采用三角形重心格点上的 nodal Lagrange Pp 基，限制
`p<=4`。该选择容易验证 Kronecker 性质、粗 P1 注入和边方向。

若以后研究更高 `p` 或局部 `p_T` 自适应，再切换到层次 Jacobi/Legendre
基；第一版不同时实现两套基函数。

### 4.4 数值积分

常系数算子至少满足：

- 刚度积分精确到次数 `2p-2`；
- 质量积分精确到次数 `2p`；
- Robin 边质量积分精确到次数 `2p`；
- coarse-P1/fine-Pp 混合质量积分精确到次数 `p+1`。

制造解含指数函数，精确误差积分使用更高阶规则，默认积分阶至少
`2p+8`，并用再次提高积分阶后的变化量作 quadrature 误差检查。

## 5. 粗空间注入与拟插值

粗空间保持

```math
V_H
=
\mathbb P_1(\mathcal T_H)\cap H^1(\Omega).
```

由于 `V_H subset V_h^p`，需要构造注入

```math
P_H^{h,p}:V_H\longrightarrow V_h^p.
```

对 nodal Pp 空间，注入值由粗 P1 基函数在所有 fine-Pp nodal
自由度位置的函数值给出，不能只填 fine mesh 顶点。

### 5.1 Pi_H^dg 的 Pp 扩展

在每个粗单元 `T` 上，以粗 P1 DG 基 `psi_{T,i}` 定义

```math
M_T^{\mathrm{dg}}
=
\left((\psi_{T,j},\psi_{T,i})_T\right)_{i,j=1}^3,
```

以及 coarse/fine 混合矩阵

```math
G_T
=
\left((\varphi_j^{h,p},\psi_{T,i})_T\right)_{i,j}.
```

局部 L2 投影矩阵为

```math
\Pi_{H,T}^{\mathrm{dg}}
=
\left(M_T^{\mathrm{dg}}\right)^{-1}G_T.
```

继续使用现有 averaging operator `E_H`：

```math
I_H
=
E_H\Pi_H^{\mathrm{dg}}.
```

### 5.2 必须验证的不变量

对每个 `(H,h,p)` 必须检查

```math
\|I_HP_H^{h,p}-I\|_{\max}
\le \mathrm{tol}_{\mathrm{nest}},
```

以及：

- 常数和粗 P1 函数被正确重构；
- `rank(I_H)=dim(V_H)`；
- `dim ker(I_H)=N_{h,p}-dim(V_H)`；
- patch 约束行不缺失、不重复；
- `p=1` 时矩阵与当前 P1 实现一致。

如果 `I_HP_H=I` 未通过，则不得进入 corrector 实验。

## 6. hp 细尺度 corrector

对粗单元 `T` 定义 patch

```math
\omega_T^\ell=N^\ell(T).
```

人工边界为

```math
\Gamma_T^{\mathrm{art}}
=
\partial\omega_T^\ell\setminus\partial\Omega,
```

物理边界为

```math
\Gamma_T^R
=
\partial\omega_T^\ell\cap\partial\Omega.
```

局部空间中的函数在人工边界取零；这要求删除人工边上的所有高阶边
自由度，而不只是端点。物理边界上的自由度保留，并装配 Robin 项。

### 6.1 原 corrector

定义

```math
\mathcal C_{T,h,p}^\ell v_H
\in W_{H,h}^{p,\ell}(T)
```

满足

```math
a_{k,\omega_T^\ell}
\left(\mathcal C_{T,h,p}^\ell v_H,w\right)
=
a_{k,T}(v_H,w)
\qquad
\forall w\in W_{H,h}^{p,\ell}(T).
```

### 6.2 伴随 corrector

双边 Petrov-Galerkin 模型还需要

```math
\mathcal C_{T,h,p}^{*,\ell}v_H
\in W_{H,h}^{p,\ell}(T)
```

满足

```math
a_{k,\omega_T^\ell}
\left(w,\mathcal C_{T,h,p}^{*,\ell}v_H\right)
=
a_{k,T}(w,v_H)
\qquad
\forall w\in W_{H,h}^{p,\ell}(T).
```

实现中必须由当前“第一变量线性”的矩阵约定推导共轭转置关系，不允许
通过复制 primal 代码猜测符号。

### 6.3 鞍点离散

记 patch Helmholtz 矩阵为

```math
A_\omega
=
K_\omega-k^2M_\omega-\mathrm i kR_\omega,
```

约束矩阵为 `B_omega`。primal 系统写为

```math
\begin{pmatrix}
A_\omega&B_\omega^*\\
B_\omega&0
\end{pmatrix}
\begin{pmatrix}q\\\lambda\end{pmatrix}
=
\begin{pmatrix}r_T\\0\end{pmatrix}.
```

正确性阶段统一使用现有 `DirectSaddle` 金标准：

- 不对完整鞍点系统使用 PCG；
- 不在第一阶段启用 shifted-GMRES；
- 不因局部 `A_omega` 接近奇异而切换到未经验证的 Schur 路径；
- 求解失败必须携带 `T,H,h,p,ell,k` 抛出清晰错误。

迭代求解器复用留到直接法通过全部收敛测试之后。

## 7. hp Petrov-Galerkin LOD

定义全局 corrector

```math
\mathcal C_{h,p}^\ell
=
\sum_{T\in\mathcal T_H}\mathcal C_{T,h,p}^\ell,
\qquad
\mathcal C_{h,p}^{*,\ell}
=
\sum_{T\in\mathcal T_H}\mathcal C_{T,h,p}^{*,\ell}.
```

trial/test 基矩阵为

```math
G_{\mathrm{trial}}
=
\left(I-\mathcal C_{h,p}^\ell\right)P_H^{h,p},
```

```math
G_{\mathrm{test}}
=
\left(I-\mathcal C_{h,p}^{*,\ell}\right)P_H^{h,p}.
```

粗 Petrov-Galerkin 系统为

```math
G_{\mathrm{test}}^*A_{h,p}G_{\mathrm{trial}}U
=
G_{\mathrm{test}}^*F_{h,p}.
```

重构的 LOD 解为

```math
u_{H,h,p}^\ell
=
G_{\mathrm{trial}}U.
```

主实验采用 two-sided 模式。`test-only` 只作为回归诊断，不与
two-sided 数据混合拟合收敛阶。

## 8. 先验误差口径与实验假设

### 8.1 统一 fine-hp 参考空间

先定义 fine-hp FEM 解

```math
u_h^p\in V_h^p:
\qquad
a_k(u_h^p,v_h)=F(v_h)
\quad\forall v_h\in V_h^p.
```

再在同一个 `V_h^p` 中构造离散全局 corrector 和局部 corrector。
因此相对于 `u_h^p`，局部 corrector 在 `W_{H,h}^p` 中不是“近似
求解另一个离散空间”，而是当前 fine problem 的精确离散 corrector。

在分辨率、离散 inf-sup 和局部化稳定条件下，目标误差结构为

```math
\|u-u_{H,h,p}^\ell\|_{1,k}
\lesssim
\underbrace{\|u-u_h^p\|_{1,k}}_{\text{fine-hp 误差}}
+
C_{\mathrm{stab}}(k)
\underbrace{H\|f\|}_{\text{P1 粗尺度误差}}
+
C_{\mathrm{stab}}(k)
\underbrace{\ell^{d/2}\vartheta_p^\ell\|f\|}_{\text{局部化误差}}.
```

这里采用保守的 patch-overlap 因子 `ell^(d/2)`。如果后续理论能得到
更小的 `ell^((d-1)/2)`，再单独更新，不在实验前预设。

### 8.2 fine-hp 误差

若制造解逐单元属于 `H^(1+s)`，标准 hp 逼近给出

```math
\|u-u_h^p\|_{1,k}
\lesssim
\left(\frac{h}{p+1}\right)^s
|u|_{H^{1+s}(\mathcal T_h)},
\qquad 0<s\le p.
```

对本计划中的光滑制造解，增加 `p` 应显著降低 fine-space 误差平台。

### 8.3 关于 H 的预期阶

因为粗空间和约束仍为 P1，主假设为

```math
\boxed{
\|u-u_{H,h,p}^\ell\|_{1,k}=O(H)
}
```

以及在对偶问题具有足够稳定性时

```math
\boxed{
\|u-u_{H,h,p}^\ell\|_{L^2(\Omega)}=O(H^2).
}
```

提高 `p` 的预期作用是延后误差平台、改善误差常数和 corrector 精度，
而不是把上述斜率变成 `p+2` 或 `p+3`。

### 8.4 H=h 时不能使用旧的零误差结论

当 `p>1` 时，即使粗细网格几何上相同，

```math
V_H^1\subsetneq V_h^p.
```

因此

```math
W_{H,h}^p\ne\{0\}
```

一般仍成立，不能宣称 `u_h^p-u_{H,h,p}^\infty=0`。只有粗空间本身也
等于 `V_h^p` 时，旧的终点精确性结论才适用。

## 9. 制造解

采用已在 Helmholtz 自适应实验中校准的制造解：

```math
\phi(t)=16t^2(1-t)^2,
```

```math
u(x,y)
=
\phi(x)\phi(y)e^{\mathrm i kx}.
```

因为 `phi` 和 `phi'` 在 `0,1` 处都为零，所以在整个边界上

```math
u=0,
\qquad
\partial_nu=0,
```

从而严格满足齐次 Robin 条件

```math
\partial_nu-\mathrm i k u=0.
```

对应右端为

```math
f(x,y)
=
-e^{\mathrm i kx}
\left[
\phi''(x)\phi(y)
+\phi(x)\phi''(y)
+2\mathrm i k\phi'(x)\phi(y)
\right].
```

必须同时计算：

- fine-hp FEM 对精确解的能量和 L2 误差；
- LOD 对精确解的能量和 L2 误差；
- LOD 与同一 `V_h^p` fine FEM 解之间的离散误差。

只报告 `u_h^p-u_LOD` 不足以说明对连续问题的收敛阶。

## 10. 收敛实验设计

### 10.1 HP1：P1 等价基线

先取

```text
k=4, p=1, mode=two-sided
```

在相同 `H,h,ell` 下比较新 hp 路径和当前 P1 路径：

- coarse/fine 网格和 prolongation 完全一致；
- `I_H` 最大差；
- stiffness、mass、Robin 矩阵最大差；
- primal/adjoint corrector 最大差；
- coarse operator 和最终解最大差；
- 精确误差和代数残差。

该步骤必须达到浮点舍入量级一致，之后才允许运行 `p>1`。

### 10.2 HP2：fine-hp FEM 校准

暂不计算 LOD，只测试 `V_h^p`：

```text
p = 1,2,3
L_h = 4,6,8,10,12
k = 4
```

使用实测 `h` 拟合：

```math
\operatorname{rate}_h^E
=
\frac{\log(E_{j-1}^E/E_j^E)}
{\log(h_{j-1}/h_j)},
```

```math
\operatorname{rate}_h^{L^2}
=
\frac{\log(E_{j-1}^{L^2}/E_j^{L^2})}
{\log(h_{j-1}/h_j)}.
```

对光滑制造解，期望能量阶接近 `p`、L2 阶接近 `p+1`。若该校准失败，
不得用相同实现解释 LOD 的 `H` 阶。

### 10.3 HP3：ell 对 p 的校准

固定一个中等规模案例，例如

```text
k=4, L_H=4, L_h=10
p=1,2,3
ell=1,2,3,4,5
```

计算相邻层差

```math
\delta_{\ell,p}
:=
\|u_{H,h,p}^{\ell+1}
-u_{H,h,p}^{\ell}\|_{1,k}.
```

并记录 primal/adjoint corrector 的 annulus energy 和约束残差。该实验
不预设 `vartheta_p` 与 `p` 无关，而是为每个 `p` 选择能使局部化误差
低于主误差的 `ell` 规则。

### 10.4 HP4：固定 master fine space 的 H 主扫描

主扫描对每个 `p` 固定同一个 master fine space：

```text
k = 4
p = 1,2,3
L_h = 12
L_H = 2,4,6,8
mode = two-sided
```

`ell(L_H,p)` 由 HP3 预先冻结，主扫描开始后不根据最终误差临时调整。

固定 `V_h^p` 的优点是：

- 所有 `H` 行共享同一个 fine FEM 参考解；
- fine-space 误差在该组内不变；
- 更容易识别 H 误差何时碰到 fine/localization 平台。

### 10.5 HP5：master-fine 加深复核

用至少两个 master fine level 重复公共粗层级：

```text
L_h = 10,12
```

若资源允许，再加入 `L_h=14`。只有当加深 fine space 后拟合的 H 斜率
和误差常数稳定，才可把观测阶解释为 LOD 的 H 收敛阶。

### 10.6 HP6：耦合 H-h 复核

作为补充，保持 NVB 层级差

```math
L_h-L_H=g
```

固定，例如 `g=6` 或 `g=8`，令 `H,h` 同时趋零。该实验检查固定 master
fine space 得到的结论是否依赖某个有限维参考空间。

耦合扫描不能替代主扫描，因为它会同时改变 fine-hp 误差和 H 误差。

### 10.7 可选波数复核

主结论先在 `k=4` 完成。之后可增加

```text
k=1,8
```

但每一行必须满足预先选择的分辨率条件

```math
kH\le c_{\mathrm{res}}.
```

不同 `k` 的数据不混合拟合一个 H 斜率。

## 11. 收敛阶计算与有效行

对相邻有效行定义

```math
\operatorname{rate}_H^E
=
\frac{\log(E_{j-1}^E/E_j^E)}
{\log(H_{j-1}/H_j)},
```

```math
\operatorname{rate}_H^{L^2}
=
\frac{\log(E_{j-1}^{L^2}/E_j^{L^2})}
{\log(H_{j-1}/H_j)}.
```

若两个层级实测 `H` 相同，则不计算 pairwise rate。

除逐对阶外，对最后至少三个有效层级做

```math
\log E=c+r\log H
```

的最小二乘拟合，同时报告 `r`、`R^2` 和所用层级，不允许只挑选最有利
的两点。

### 11.1 误差地板判据

一行进入 H 阶拟合前，应满足

```math
\|u-u_h^p\|_{1,k}
\le 0.1\|u-u_{H,h,p}^\ell\|_{1,k},
```

以及

```math
\delta_{\ell,p}
\le 0.1\|u-u_{H,h,p}^\ell\|_{1,k}.
```

若不满足，该行仍写入 CSV，但标记为：

- `fine_floor=1`，或
- `localization_floor=1`。

它不得用于声称 H 渐近阶。

## 12. 模块化实现

不把高阶自由度逻辑塞进现有 `TriMesh` 顶点数组。建议新增以下模块。

### 12.1 通用 hp-FEM 层

```text
include/fem/hp_tri_space.h
src/fem/hp_tri_space.cpp
```

职责：

- 连续 Pp 全局自由度编号；
- 顶点、边、内部自由度分类；
- 局部到全局映射；
- 边方向置换；
- patch 人工边界自由度识别；
- P1 到 Pp 注入。

```text
include/fem/hp_triangle_basis.h
src/fem/hp_triangle_basis.cpp
```

职责：

- 参考三角形 nodal Pp 基；
- 基函数值和参考梯度；
- 仿射映射后的物理梯度。

```text
include/fem/quadrature.h
src/fem/quadrature.cpp
```

职责：

- 三角形和边积分规则；
- 按所需精确次数选择规则；
- 制造解误差的高阶积分。

### 12.2 Helmholtz hp 层

```text
include/helmholtz/hp_operators.h
src/helmholtz/hp_operators.cpp
```

组装复值 stiffness、mass、Robin、load 和误差矩阵。

```text
include/helmholtz/hp_interpolation.h
src/helmholtz/hp_interpolation.cpp
```

构造 `P_H^(h,p)`、`Pi_H^dg`、`I_H` 和 patch 约束行。

```text
include/helmholtz/hp_corrector.h
src/helmholtz/hp_corrector.cpp
```

提取 patch hp 自由度、处理人工/物理边界、组装原/伴随鞍点问题。

```text
include/helmholtz/hp_model.h
src/helmholtz/hp_model.cpp
```

提供面向实验的高级 API：

```cpp
HelmholtzHpProblemConfig config;
config.H = coarse_level;
config.h = fine_level;
config.p = polynomial_degree;
config.ell = oversampling;
config.wavenumber = k;

auto model = HelmholtzHpLodModel::build(config);
auto lod = model.solve_source(source);
auto fine = model.solve_fine_reference(load);
```

不要在第一版直接修改现有 `HelmholtzLodModel` 的内部数据布局。先建立
独立 hp 路径，并用 `p=1` golden test 证明等价；稳定后再评估是否抽取
共享模板。

### 12.3 Benchmark

新增

```text
benchmarks/bench_helmholtz_hp_H_convergence.cpp
```

建议参数：

```text
--k=4
--p-levels=1,2,3
--h=12
--H-levels=2,4,6,8
--ell=auto|N
--ell-offset=N
--mode=two-sided|test-only
--format=table|csv
--check
```

若 `--ell=auto`，必须使用已冻结的 `ell(H,p)` 规则，而不是读取当前
误差后动态选择。

## 13. 测试计划

### 13.1 hp 基础测试

1. nodal Kronecker 性质；
2. 分片单位性和梯度和为零；
3. Pp 多项式重构；
4. 相邻三角形共享边迹连续；
5. 正、反局部边方向得到同一个全局函数；
6. 自由度计数公式；
7. 三角形/边 quadrature 多项式精确性；
8. `p=1` 与现有 P1 basis 和 assembly 一致。

### 13.2 拟插值测试

1. `I_H P_H=I`；
2. coarse P1 重构；
3. `rank(I_H)=dim(V_H)`；
4. 随机核向量满足 `I_Hw=0`；
5. p=1 与现有 `E_H Pi_H^dg` 一致；
6. 不同 NVB 层级和边方向均覆盖。

### 13.3 Helmholtz 算子测试

1. stiffness Hermitian；
2. mass、Robin mass Hermitian positive semidefinite；
3. `A=K-k^2M-ikR` 符号正确；
4. 整域 Robin 边界的所有高阶边自由度被装配；
5. 制造解载荷随 quadrature 阶提高收敛；
6. fine Pp FEM 能量/L2 阶符合预期。

### 13.4 Corrector 测试

覆盖内部、物理边界和角点 patch：

- primal equation residual；
- adjoint equation residual；
- constraint residual；
- 人工边界所有 hp 自由度为零；
- 物理边界高阶自由度未被错误删除；
- direct saddle 与显式小规模 null-space 金标准一致；
- p=1 corrector 与现有实现一致。

### 13.5 全模型测试

- p=1 hp 模型与现有 `HelmholtzLodModel` 一致；
- two-sided trial/test 共轭关系正确；
- coarse Petrov residual 达标；
- repeated RHS 不重算 corrector；
- 输入 `p<1`、`H>h`、不支持的局部次数分布时清晰失败；
- 小规模制造解 H 扫描作为 CTest smoke test。

## 14. 输出规范

每个 CSV 行至少包含：

```text
k, mode, L_H, L_h, p, ell
H, h, kH, h_over_H
coarse_nodes, coarse_elements
fine_vertices, fine_elements, hp_dofs
trial_dofs, test_dofs
fine_energy_error, fine_l2_error
lod_energy_error, lod_l2_error
lod_vs_fine_energy, lod_vs_fine_l2
energy_rate_H, l2_rate_H
fine_floor, localization_floor, used_in_fit
primal_corrector_residual
adjoint_corrector_residual
constraint_residual
petrov_residual
mesh_ms, operators_ms, correctors_ms
basis_factorization_ms, fine_solve_ms, total_ms
peak_rss_kib
```

结果目录建议为

```text
results/helmholtz_hp_H_convergence/
```

包含：

- `metadata.txt`：提交号、编译器、依赖、CPU、线程数和完整命令；
- `fine_hp_calibration.csv`；
- `ell_by_p_calibration.csv`；
- `H_convergence_fixed_fine.csv`；
- `H_convergence_coupled.csv`；
- 对应 `.log` 和 `/usr/bin/time -v` 输出。

实测结论写入 `DEVELOPMENT.md`，本计划书只更新阶段状态、验收结果和
必要的数学口径，不复制完整时间表。

## 15. 验收标准

### 15.1 正确性门槛

- p=1 新旧路径最终解相对差不超过 `1e-10`；
- `I_H P_H-I` 最大范数不超过 `1e-11`；
- primal、adjoint、constraint、Petrov 残差默认不超过 `1e-9`；
- 提高 quadrature 阶后误差相对变化不超过 `1e-8`；
- 16 项现有 CTest 全部继续通过；
- 新增 hp 测试全部通过。

容差可按条件数放宽，但必须记录原因，不能静默修改。

### 15.2 数值结论门槛

只有同时满足以下条件，才可声称观察到关于 H 的渐近阶：

1. 至少三个有效粗层级；
2. fine-floor 和 localization-floor 判据均通过；
3. 回归 `R^2>=0.98`；
4. 加深 master fine space 后斜率变化不超过 `0.15`；
5. 增加 `ell` 后斜率变化不超过 `0.15`；
6. 所有代数残差低于对应离散误差的 `1%`。

预期但不强制写入验收门槛的结果是：

```math
r_E\approx1,
\qquad
r_{L^2}\approx2.
```

如果数据不符合预期，必须优先排查 fine/localization/quadrature/代数
误差地板，而不是直接修改拟合区间。

## 16. 性能与内存

第一版以正确性为主，但必须记录高阶自由度增长。二维拟一致网格上，
单个 patch 的 hp 自由度大致随

```math
N_{\omega,p}
\sim
p^2
\left(\frac{\operatorname{diam}\omega}{h}\right)^2
```

增长；直接鞍点分解的时间和内存会增长得更快。

实施规则：

- p=1,2,3 先串行或低并发校准；
- corrector 并行时限制同时存活的 patch 分解数；
- 不在 OpenMP 外层并行中再启用稠密线性代数多线程；
- 每个 `(p,H,h,ell)` 记录峰值 RSS；
- 超过预设内存预算时减少并发，而不是改变数学参数；
- direct saddle 正确后，再单独评估复用 symbolic pattern、
  shifted-GMRES 或 block multi-RHS。

## 17. 实施顺序

### 第一步：只实现连续 Pp FEM

先完成自由度、基函数、积分、Robin 装配和制造解 fine FEM 阶验证。
此时不接入 LOD。

### 第二步：实现 coarse P1 到 fine Pp 的注入和 I_H

冻结 `I_H P_H=I` 金标准，确保核空间定义正确。

### 第三步：实现单 patch direct-saddle corrector

先测一个内部 patch、一个边界 patch和一个角点 patch，再扩展到全部
element correctors。

### 第四步：构造 hp Petrov-Galerkin LOD

先完成 p=1 新旧模型等价，再启用 p=2,3。

### 第五步：依次做三类误差校准

顺序固定为：

```text
fine-hp FEM -> ell-by-p -> fixed-master-fine H scan
```

不能跳过前两项直接解释 H 斜率。

### 第六步：复核和文档

完成 master-fine 加深、耦合 H-h 和可选波数复核，更新阶段状态、
`DEVELOPMENT.md` 和结果目录。

## 18. 明确不在本计划中的内容

- 局部 `p_T` 自适应；
- 非拟一致粗网格；
- 独立 patch 网格；
- H-h-ell-p 联合后验估计子；
- 高阶 coarse constraint/QOI；
- PML；
- 三维；
- 以 GMRES 性能替代 direct-saddle 正确性金标准。

这些方向必须在本实验确认连续 `V_h^p` corrector 和 H 收敛行为之后，
再建立独立计划。
