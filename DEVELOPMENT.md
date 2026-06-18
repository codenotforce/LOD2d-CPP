# LOD2d-C++ 开发日志

> MATLAB → C++ 迁移全记录：每次尝试、成功与失败的原因。

## 项目概览

将 MATLAB 版 LOD2d 迁移到 C++20，目标是在相同参数下运行时间显著低于 MATLAB 优化版（`parfor + RCM + +lod`）。

| 指标 | MATLAB 串行 | MATLAB 优化 (4w) | C++ (4t) | C++ (16t) |
|------|------------|-----------------|----------|-----------|
| H=4,h=8 矫正子 | 25.2 s | 2.3 s (热) / 6.1 s (温) | 5.5 s | 2.9 s |
| H=4,h=8 总时间 | 27.0 s | 3.4 s (热) / 9.3 s (温) | 6.6 s | 4.3 s |

**平台:** WSL2 Ubuntu 22.04, g++ 11.4, CMake 3.16+, Eigen 3.4, SuiteSparse 5.10
**硬件:** 16 核 (i7-13700H?), 32 GB RAM

---

## 阶段 A–G：模块迁移

### Phase A: 工程骨架 ✅

建立 CMake 工程，引入 Eigen / SuiteSparse / OpenMP。

**成功路径:** `CMakeLists.txt` + 目录结构 + 最小 `main.cpp` 验证。

---

### Phase B: 网格细化 ✅ 36/36

从 MATLAB `+lod/refineMesh` 迁移到 C++。

**尝试:**
1. 使用 `std::map<Edge, int>` 排序边 → **失败**。行优先 vs MATLAB `find(sparse)` 的列优先排序不一致。
2. 改为自定义比较器 `sort by second, then first` → **成功**。列优先排序匹配 MATLAB。
3. 三角剖分按父单元分组 → **失败**。MATLAB 按子三角类型（sub-1, sub-2, …）分组。
4. 改为按子类型分组 → **成功**。
5. `kron(speye, full12x3)` 构建 P_dg → **失败**。行排序不对。
6. 改为 `[kron(E,sub1); kron(E,sub2); …]` 垂直拼接 → **成功**。
7. P_dg 零值 triplet → **失败**。Eigen `setFromTriplets` 不跳过显式零。
8. 加 `if (val != 0.0)` 过滤 → **成功**。

**最终方法:** 列优先边排序 + 子类型三角分组 + kron 垂直拼接 P_dg + 零值过滤。36/36 金标准通过。

---

### Phase C: DG 装配 ✅ 10/10

迁移 `lod.assembleDG` — 逐单元 DG 刚度矩阵装配。

**方法:** 直接计算每单元 3×3 刚度 `Ke(i,j) = cA · (∇φ_i · ∇φ_j)`，triplet 组装到全局稀疏矩阵。

**尝试:**
- 9 个 triplet 全部加入 → **失败**。Eigen 保留显式零，nnz 多 28%。
- 加 `if (val != 0.0)` → **成功**。10/10 金标准通过。

---

### Phase D: 拟插值 ✅ 3793/3793

迁移 `lod.buildQuasiInterp` — Mhdg, B, PiHdg, EH, IH 完整组装。

**方法:** 块对角 mass 矩阵 + 三步稀疏乘 + 顶点度数平均。

**尝试:**
- `IHp` 稀疏遍历替换 `coeff()` 双层循环 → **失败**（见下文 Phase F）。回退到原始 `coeff()` 循环。

**最终方法:** `coeff()` 双层循环（IH 很小：NH×Nh ≈ 81×1089）。3793/3793 金标准通过。

---

### Phase E: 补丁构建 ✅ 6/6

迁移 `lod.buildPatches` — 顶点-单元关联 → ℓ 步邻接扩展。

**方法:** `Ivt` (vertex→element) → `Itt = spones(Ivt' * Ivt)` → ℓ 次 `patch = Itt * patch` → 二值化 + 对角线。

**尝试:**
- `coeffRef(i,i) = 1.0` 设置对角线 → **失败**。Eigen 压缩稀疏矩阵 `coeffRef` 只在已分配位置生效。
- 改为 triplet 重建 → **成功**。`setFromTriplets` 正确求和重复条目。

**最终方法:** triplet 重建（阈值为 1.0，对角线显式添加）。6/6 金标准通过。

---

### Phase F: 矫正子求解 ✅ 3/3

迁移 `lod.computeCorrector` — 每个粗单元求解局部 saddle-point 问题。

**方法:** Sph 组装 → rhsp 组装 → IHp 提取 → `Sph \ RHS` → mu 计算 → 存储。

---

### Phase G: 完整 LOD ✅ 3/3

串联所有模块，与 MATLAB 金标准对比。H=3,h=5 全流程 3/3 通过。

---

## 性能优化尝试

### ✅ 成功

| # | 优化 | 方法 | 效果 |
|---|------|------|------|
| 1 | **块对角 Shdg O(1) 取块** | 对每个 patch 单元，直接 `InnerIterator(Shdg, 3e+ci)` 读取 3 列，而非扫描全矩阵 | ~10% 矫正子加速 |
| 2 | **稀疏 Cholesky** | `Eigen::SimplicialLLT` 替代 `toDense().ldlt()` | ~100×（矫正子求解器） |
| 3 | **OpenMP 并行** | `#pragma omp parallel for schedule(dynamic)` | ~12×（16 线程） |
| 4 | **`unordered_map` DOF 查表** | O(1) 替代 O(N) 线性搜索 | 粗求解 ~100× |
| 5 | **零值过滤** | 装配时 `if (val != 0.0)` 跳过 | 减少稀疏存储 |
| 6 | **Release 模式** | `-O3 -march=native -DNDEBUG` | ~8× vs Debug |
| 7 | **RHS 维数修复** | `f_coarse = ones(NH)` 替代 `f_vec = ones(Nh)` | 修复崩溃 bug |

### ❌ 失败

| # | 尝试 | 路径 | 失败原因 |
|---|------|------|---------|
| 1 | **CHOLMOD 替代 SimplicialLLT** | `include/solver/cholmod_wrapper.h` → `Eigen→triplet→CHOLMOD→solve` | N<8000 时转换开销（~3ms）超过求解器收益（~0.5ms）。Eigen SimplicialLLT 比 CHOLMOD 快 3–50×。`cholmod_l_*` (64-bit API) segfault，改用 `cholmod_*` (32-bit API)。 |
| 2 | **IHp 稀疏遍历** | `for (InnerIterator it(IH, row); …)` + `std::lower_bound(dofph)` 替代 `coeff()` 双层循环 | `IH` 的 InnerIterator 返回的非零集与 `coeff()` 不一致——`coeff()` 对某些 (i,j) 返回非零值但 InnerIterator 跳过。导致金标准测试全部失败（max err 0.56）。 |
| 3 | **symrcm 重排 Sph** | 对每个 patch 的 Sph 做 `symrcm` 再求解 | Cholesky fill-in 增加 ~5%（矩阵本身来自网格细化，已接近最优排序）。额外 ~5ms 开销抵消求解器收益。 |
| 4 | **预计算全部稀疏子矩阵** | 串行预提取每个元素的 Shdg/C_ell 子矩阵 | 预计算耗时 ~10s，但 4 线程 parfor 并行执行同样操作仅需 ~0.1s。串行化抵消了并行收益。 |
| 5 | **预计算索引集**（轻量版）| 只预计算 dofph/dofphdg（不存矩阵值） | 节省 ~3ms/elem，但预计算 0.35s + 广播开销 → 净负收益（0.83×）。 |
| 6 | **`coeffRef(i,i) = 1.0`** | 直接设置 Eigen 稀疏矩阵对角线 | Eigen 压缩稀疏矩阵的 `coeffRef` 只在既有非零位置工作。新位置需 `insert` 或重建为 triplet。 |
| 7 | **GPU (gpuArray)** | 在 MATLAB 侧测试 GPU 加速 | MATLAB `sparse gpuArray` 不支持任意索引向量——`Shdg(dofs, dofs)`（核心瓶颈）无法在 GPU 运行。 |
| 8 | **pcg + ichol** | 对 Nph>1000 的 patch 用 pcg+ichol 替代 backslash | 多 RHS（6–30 列）时 pcg 需 K×30 次迭代，总开销远超 1 次 Cholesky 分解。所有实测规模下 backslash 均更快。 |

### ⚪ 不需要

| # | 策略 | 原因 |
|---|------|------|
| 1 | `parallel.pool.Constant` | Threads 池共享内存，无需显式常量。仅 Process 池需要。 |
| 2 | 负载均衡排序 | `parfor` / `schedule(dynamic)` 已自动处理。最大 patch 与最小 patch 差距 <2×。 |
| 3 | 融合稀疏乘 `A'*(B*A)` | MATLAB 已内部融合此表达式，C++ 同理。 |

---

## 关键技术决策

### 1. 列优先边排序

MATLAB `find(sparse(row, col, 1))` 返回**列优先**顺序的条目。C++ `std::map` 默认行优先。必须用自定义比较器 `sort by (j, i)` 匹配，否则所有下游结构（节点编号、连通性、延长矩阵）全部错位。

### 2. P_dg 行排列

MATLAB `[kron(E,sub1); kron(E,sub2); …]` 垂直拼接使 P_dg 的**行按子三角类型分组**，而非按元素分组。直接 `kron(speye(nt), full12x3)` 产生按元素分组的行序——二者等价但 nnz 位置不同，导致后续稀疏乘结果不同。

### 3. 子三角排列

MATLAB 细化后三角剖分按 `[sub1(所有元素); sub2(所有元素); sub3; sub4]` 排列。C++ 原始的 `[sub1..4(元素0); sub1..4(元素1); …]` 导致 P_elem 和 P_dg 结构错位。

### 4. CHOLMOD 不可行

对于 N<8000 的 SPD 矩阵，Eigen→triplet→CHOLMOD 的转换开销（~3ms/次）完全覆盖了求解器本身的优势（~0.5ms/次）。Eigen `SimplicialLLT` 直接操作 Eigen 稀疏格式，零转换开销。CHOLMOD 只有在 N>20000 时才可能反超，而 LOD 矫正子的 Sph 规模在 1000–50000 范围。

### 5. 构建模式

Debug 模式下 C++ 比 MATLAB 慢 10×。必须使用 `-DCMAKE_BUILD_TYPE=Release`（含 `-O3 -march=native -DNDEBUG`）才能获得可比性能。此配置下 Eigen 的维度断言被禁用——这是潜在风险，也是之前 `f_vec = ones(Nh)` 维数 bug 未被检测的原因。

---

## 待解决问题

1. **C++ 4 线程 vs MATLAB 4 线程**：矫正子 5.5s vs 2.3s（热）/ 6.1s（温）。MATLAB JIT 对稀疏操作的优化（尤其是 `Shdg(idx, idx)` 和 `cg2dgh(idx, :)')` 的内部实现）是 C++ 难以匹敌的。

2. **跨平台一致性**：WSL→Windows 文件同步导致 CRLF 警告和构建缓存问题。建议将 Windows 侧 `LOD2d_C++/` 改为只读镜像，所有编辑和构建在 WSL 内完成。

3. **单元刚度预计算**：当前每个矫正子需要从 Shdg 读取 Ke 块（InnerIterator × 3 列 × 5000 元素）。若将 Shdg 的每元素 3×3 块预存为 `vector<Matrix3d>`，可消除所有 Shdg 访问。

---

## 文件结构

```
lod2d-cpp/
├── CMakeLists.txt
├── include/
│   ├── mesh/       types.h, refine.h
│   ├── fem/        assemble_dg.h
│   ├── lod/        quasi_interp.h, patches.h, corrector.h
│   └── solver/     coarse_solve.h, cholmod_wrapper.h
├── src/            实现文件（对应 include 结构）
├── tests/          金标准测试（test_mesh, test_golden, test_dg, test_qi, test_patch, test_corr, test_full, debug_cholmod）
├── benchmarks/     bench_H4h8, bench_profile, bench_refine, data_H4h8.txt
└── README.md
```
