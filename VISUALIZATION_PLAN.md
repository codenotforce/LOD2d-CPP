# LOD2d-C++ 可视化与论文制图计划书

## 1. 计划目标

本计划以复值 Helmholtz Petrov–Galerkin LOD 为主，包括当前已经实现的第一阶段自适应粗网格流程。实值椭圆型 LOD 只承担与 `basicLOD2d.m` 对照的最小回归用途，不追求与 Helmholtz 相同的可视化功能覆盖。

实施时遵守一个优先于所有图片功能的约束：

> 可视化默认关闭时，不增加额外参考求解、校正子保留、矩阵复制、逐迭代快照、文件 I/O 或与迭代次数成正比的常驻内存。

### 1.1 当前实施状态（2026-07-27）

阶段 1–3 已完成并在 WSL Release 环境通过测试：

- 新增独立的 `lod2d_visualization_io` 流式 VTU 导出库，字段使用非拥有式只读
  `std::span`；该库不链接到 `lod2d_core`；
- 新增显式运行的 `bench_helmholtz_visualization`。只有运行此程序时才计算
  fine FEM reference、构造临时粗尺度延拓/细尺度修正并写出 VTU/JSON；
- 写出粗/细网格、精确解、fine FEM、完整 LOD 解、粗尺度延拓、细尺度修正、
  LOD-reference 误差、介质系数和可追溯 `run.json`；
- 生成 Helmholtz 实部、虚部、模、相位、解分解、误差、网格与固定截线论文图；
- 从 manufactured wave-number scan 生成固定 `kH=1` 的 FEM/LOD 污染对比图；
- 从 `results/helmholtz_manufactured/validation.csv` 生成误差–DOF 双对数图。
  横轴按 DOF 递增，拟合斜率为负，曲线从左上向右下；
- C++ 测试覆盖字段长度、连接关系和复数拆分；Python 测试覆盖 manifest/VTU
  一致性、`u_lod = u_coarse + u_fine_scale`、误差字段、负斜率和 headless 制图。

现有求解器、实值 `LodModel` 的临时校正子释放策略和原 benchmark 均未加入
可视化状态、额外参考求解或默认 I/O。阶段 4 的自适应逐轮快照与动态加细动画
尚未实施。

### 1.2 Helmholtz `H` 收敛数据接口（2026-07-27）

`bench_helmholtz_H_convergence` 已为固定细网格、global-NVB 粗网格序列实现
专用 CSV 接口。该接口服务于绝对 \(k\)-加权能量误差论文图，不改变通用
VTU/自适应可视化路线：

- `all_results.csv`/`summary.csv` 使用实测 `H_max`、实际 `coarse_nodes`
  以及 `p1_energy_abs`、`lod_energy_abs` 作收敛图；
- `H_L##_coarse_nodes.csv` 和 `H_L##_coarse_elements.csv` 保存每个粗层的
  坐标、拓扑、面积、直径、P1 值和 LOD 粗系数；
- 固定细网格 `fine_L##_nodes.csv`、`fine_L##_elements.csv` 只保存一次；
- 显式 `--export-fields` 后，`H_L##_fields.csv` 保存 fine-node 上的 LOD
  重构、粗尺度分量、校正分量、P1 延拓、精确解和点误差；
- C++ 写出的实部/虚部是权威数据，模、相位和论文配色仍由可视化项目处理；
- 收敛图必须标为 absolute energy error，不能复用污染实验的 relative
  error 图注。

默认服务器注册参数和字段解释见
`HELMHOLTZ_H_CONVERGENCE_SERVER_RUNBOOK.md`。level 21 的逐层 fine fields
体量很大，单纯绘制误差-自由度曲线时保持 `EXPORT_FIELDS=0`。

计划需要满足以下三项核心需求：

1. 绘制粗、细网格，标出误差指示子较大的单元、Dörfler 标记单元和 NVB 闭包单元，并能生成动态加细过程；
2. 绘制误差关于自由度的收敛图，尤其是相对 \(\{1,k\}\) 范数（本文统一记为 \(\|\cdot\|_{1,k}\)）关于自由度的双对数图；
3. 绘制参考解、LOD 全解、粗尺度有限元部分和细尺度修正，确保 LOD 基函数带来的微小振荡不会在插值、降采样或图片压缩中丢失。

此外，计划还给出适合论文使用的扩展图表、数据格式、程序结构、复现规范和分阶段验收标准。

## 2. 阅读项目后的现状结论

### 2.1 当前项目已经具备的数据

项目不是单一的网格测试程序，而是已经包含：

- 实值椭圆型 LOD 的粗细网格、准插值、局部校正子、LOD 基和重复右端项求解；
- 复值 Helmholtz 的粗细网格、主/伴随校正子、修正试验/检验基、LOD 解和细网格参考解；
- LEB、红细化、NVB 和带闭包的局部 NVB；
- 自适应 Helmholtz LOD 的逐步误差、估计子、标记数量、闭包增长、稳定性和计时数据；
- manufactured solution、\(H\)、\(h\)、\(\ell\)、\(k\)、\(p\)、预条件器和并行性能实验的 CSV 输出。

与可视化直接相关的核心数据结构包括：

- `TriMesh.nodes`：二维顶点坐标；
- `TriMesh.elems`：三角形顶点连接关系；
- `HelmholtzLodSolution.coarse_coefficients`：复值粗尺度系数；
- `HelmholtzLodSolution.fine_values`：复值细网格上的完整 LOD 解；
- `HelmholtzLodModel::trial_basis()`、`corrected_trial_basis()` 及主/伴随校正子；
- 自适应流程中的逐单元 `indicator_squared`、层数、稳定单元 ID 和误差历史。

实值 `LodReuseSolution.uH/uHms` 只在 MATLAB 对照测试中使用，不作为本计划第一阶段的主要数据源。

现有 Helmholtz 自适应代码中的离散能量范数已经是

\[
\|v\|_{1,k,h}
=
\left(v^\ast S_h v+k^2v^\ast M_hv\right)^{1/2}
=
\left(\|\nabla v\|_{L^2}^2+k^2\|v\|_{L^2}^2\right)^{1/2}.
\]

因此用户要求的相对 \(\{1,k\}\) 范数图不需要重新定义数值算法，只需统一命名、补充自由度字段，并把现有结果稳定地导出和绘图。

### 2.2 当前可视化相关缺口

目前项目没有统一的网格/场数据导出层，也没有 Python 绘图目录。主要缺口是：

- 自适应驱动只保存标量历史和最终网格，没有保存每次迭代的网格快照；
- 最终网格 CSV 有单元层数、坐标和指示子，但没有逐步的 `marked`、`eligible`、`closure` 和 `new_child` 标志；
- 没有按需导出选定迭代的 LOD 解、参考解、粗尺度插值部分和细尺度修正；
- 当前没有“只读借用现有 Helmholtz 数据、写完立即释放”的低开销快照接口；
- 已有 benchmark CSV 的列名和自由度口径尚未完全统一；
- 没有统一的论文配色、字号、LaTeX 标签、图片尺寸和输出格式；
- 没有自动检查“图片所用数据”和 C++ benchmark 数值是否一致。

### 2.3 与 `basicLOD2d.m` 的对应关系

`D:\code\femcode\LOD2d_MATLAB\src\basicLOD2d.m` 最后的原始可视化包含四部分：

1. 随机扩散系数 \(A\)；
2. 细网格参考解 \(u_h\)；
3. 细网格上的完整 LOD 近似 \(u^{ms,h}_{H,\ell}\)；
4. 粗网格上的有限元部分 \(u^h_{H,\ell}\)。

它还计算了非相对的能量误差、LOD 的 \(L^2\) 误差和粗尺度有限元部分的 \(L^2\) 误差。

这组实值图只作为迁移正确性的可选对照基线，不应驱动 `LodModel` 长期保留局部校正子，也不排在 Helmholtz 主线之前。Helmholtz 主线进一步增加：

- 完整 LOD 解与粗尺度插值之差；
- 点态/单元误差图；
- 基函数和校正子图；
- 自适应网格和误差指示子图；
- Helmholtz 复值解的实部、虚部、模和相位图。

### 2.4 关于现有内存优化的判断

实值 `LodModel` 构造后丢弃局部校正子的临时集合，应视为有意的内存优化，而不是需要通过扩大模型状态来“修复”的可视化缺口。本计划明确：

- 不修改实值 `LodModel` 的默认对象生命周期和常驻数据；
- 不为实值对照图增加“保留全部校正子”的配置；
- 如果以后确实需要实值基函数图，应在校正子刚计算完成、尚未释放时，通过仅针对指定单元/指定列的同步回调导出，随后继续按原逻辑释放；
- Helmholtz 可视化只读取模型当前本来就持有的 `problem`、`operators`、`trial_basis`、`correctors` 和解向量，不再复制一套同规模对象；
- 自适应每轮的模型在销毁前完成可选导出，不把历轮模型或历轮解保存在 `AdaptiveHelmholtzResult` 中。

## 3. 技术选型

### 3.1 推荐结论

采用分层方案：

| 层次 | 编程语言/工具 | 职责 |
|---|---|---|
| 数值计算层 | C++20 | 求解、误差计算、网格状态和场数据导出 |
| 数据交换层 | VTK XML `.vtu`、CSV、JSON | 保存非结构三角网格、点数据、单元数据、收敛历史和运行元数据 |
| 论文静态图 | Python 3.10 + NumPy + Pandas + Matplotlib | 二维网格、等值填色、收敛曲线、组合图，输出 PDF/SVG/PNG |
| 三维/交互/动画 | Python + PyVista（基于 VTK） | 三维表面、交互查看、网格快照、MP4/GIF |
| 人工检查 | ParaView，可选 | 检查 `.vtu` 字段、单元标记和时间序列，不作为论文自动化主流程 |

不建议在 C++ 求解器内直接使用 Qt、OpenGL 或 Matplotlib C++ 封装。原因是这会把图形界面和渲染依赖带入核心库，增加 WSL 构建、服务器运行和测试的复杂度，并且论文排版控制仍然不如 Python 灵活。

### 3.2 为什么选择 VTK/PyVista

项目使用局部细化后的非结构三角网格，VTK 的 `UnstructuredGrid` 正好表达“顶点 + 单元拓扑 + 点/单元属性”。VTK XML 的 `.vtu` 是非结构网格的标准扩展名，并支持便携编码、压缩和并行扩展；PyVista 可以由数组或文件直接构造 `UnstructuredGrid`，也能离屏截图和生成电影。

第一版在 C++ 中实现一个轻量的 `.vtu` XML 写出器，不把完整 VTK C++ 库加入求解器依赖。初期可采用 ASCII XML 便于检查；数据规模增大后再加入 appended binary/compression。CSV 继续用于收敛历史，JSON 用于保存运行元数据。

相关官方资料：

- [VTK XML 文件格式与 `.vtu`](https://docs.vtk.org/en/v9.6.1/vtk_file_formats/vtkxml_file_format.html)
- [PyVista `UnstructuredGrid`](https://docs.pyvista.org/api/core/_autosummary/pyvista.unstructuredgrid)
- [PyVista 离屏截图](https://docs.pyvista.org/api/plotting/_autosummary/pyvista.plot.html)
- [PyVista MP4 写出](https://docs.pyvista.org/api/plotting/_autosummary/pyvista.plotter.open_movie)
- [ParaView 文件时间序列](https://docs.paraview.org/en/latest/UsersGuide/dataIngestion.html#handling-temporal-file-series)

### 3.3 为什么静态论文图以 Matplotlib 为主

Matplotlib 能直接在非结构三角网格上使用 `tricontourf`，适合误差、指示子和解场的俯视图；收敛图、局部放大、参考斜率和多面板排版也更容易统一。图片可直接保存为 PDF/SVG 或高分辨率 PNG。

相关官方资料：

- [Matplotlib 非结构三角网格填色 `tricontourf`](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.tricontourf.html)
- [Matplotlib 图片输出 `savefig`](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html)
- [Matplotlib 动画写出](https://matplotlib.org/stable/api/animation_api.html)

### 3.4 WSL 与 Windows 工作流

已只读核对指定 WSL 环境：

- Ubuntu 22.04；
- Python 3.10.12；
- CMake 3.22.1；
- g++ 11.4.0。

建议保持现有原则：

1. C++ 构建、测试、数值运行和图片生成全部在 WSL 仓库完成；
2. 所有脚本只使用仓库相对路径，不在代码中写死 `D:\...` 或 `\\wsl.localhost\...`；
3. 在 WSL 生成原始数据和候选图；
4. 只把需要版本控制的源码、绘图脚本、小型基准数据和最终论文图同步到 Windows 镜像；
5. GitHub 提交和推送继续在 Windows 仓库完成；
6. 大型 `.vtu`、逐帧 PNG 和 MP4 默认不进 Git，只保留生成命令、元数据和必要的小型回归样例。

Windows 仓库 `D:\code\femcode\LOD2d_C++` 与 WSL 仓库
`/home/qcxubuntu/learning/LOD2d-C++` 中受版本控制的源文件应保持一致。每次同步计划书或后续可视化源码后，用 SHA-256 逐文件校验；只有 WSL 完成构建/测试且 Windows 镜像校验一致后，才从 Windows 提交 GitHub。两边各自已有的 `build/`、大型 `results/` 和本地未提交实验数据不要求镜像一致，也不得在同步源码时被覆盖或删除。

## 4. 建议的程序结构

建议新增如下目录和文件，具体命名可在实现时微调：

```text
include/io/
  vtk_writer.h                 # 与求解器无关的轻量 VTU 数据结构和写出接口
  visualization_snapshot.h    # 非拥有式只读视图、字段选择和迭代信息

src/io/
  vtk_writer.cpp
  visualization_snapshot.cpp

tools/visualization/
  requirements.txt             # 锁定 Python 依赖版本
  style.py                     # 论文统一样式、色图、尺寸
  io.py                        # 读取 VTU/CSV/JSON 和字段检查
  plot_mesh.py                 # 网格、层数、标记、闭包和指示子
  plot_solution.py             # 实值/复值解、误差、细尺度修正
  plot_convergence.py          # DOF-误差、效应指数、参考斜率
  plot_basis.py                # 基函数、校正子、局部 patch
  animate_adaptivity.py        # 自适应加细 MP4/GIF
  reproduce_all.py             # 一条命令生成论文所需全部图片

tests/
  test_visualization_export.cpp

tools/visualization/tests/
  test_vtu_schema.py
  test_plot_smoke.py

results/visualization/
  <case-id>/                   # 原始运行数据，默认不提交

figures/
  paper/                       # 最终选定、可复现的论文图
```

### 4.1 求解器和绘图器之间的边界

C++ 不负责颜色、相机、字体或图片格式。C++ 只保证：

- 网格拓扑正确；
- 数值字段命名和含义稳定；
- 复数被无损拆成实部与虚部；
- 误差和自由度口径明确；
- 输出包含完整复现元数据。

导出器不拥有求解器对象。它只在调用期间借用网格、向量和标量数组，写出完成后立即返回；禁止把 `HelmholtzLodModel`、稀疏矩阵或大向量复制进“可视化结果对象”。

Python 不重新求解 PDE，也不重新定义论文中的误差。Python 只负责：

- 读取 C++ 输出；
- 由实部、虚部派生模和相位并交叉验证；
- 选择色标和图形布局；
- 生成静态图、交互图和动画。

### 4.2 自适应快照接口

不应让核心自适应驱动直接依赖文件系统。建议在 benchmark/实验驱动层提供可选观察器，核心自适应代码只调用抽象回调。接口采用非拥有式只读视图，例如：

```cpp
struct AdaptiveSnapshotView {
    int iteration;
    const TriMesh &coarse_mesh;
    std::span<const int> coarse_levels;
    std::span<const std::uint64_t> element_ids;
    std::span<const double> indicator_squared;
    std::span<const int> marked_elements;
    const HelmholtzLodSolution *lod;   // 未请求解场时为 nullptr
    const ComplexVector *reference;    // 未计算参考解时为 nullptr
};

std::function<void(const AdaptiveSnapshotView &)> snapshot_callback;
```

该视图只在回调执行期间有效，不允许缓存其中的引用。基础网格快照至少暴露：

- 迭代编号；
- 加细前的粗网格、层数和稳定单元 ID；
- 当前选用的指示子；三类指示子只有在命令行显式请求比较时才全部暴露；
- `eligible` 和 Dörfler `marked` 集合；
- 加细后产生的新单元 ID；
- 由闭包引起的额外加细单元；
- 当前 LOD 解和已经存在的细网格参考解，均为可空只读指针；
- 当前误差、自由度、\(H_{\max}\)、\(q_{\max}\)、\(\ell\)、\(k\)。

回调必须是同步的。第一版不使用后台写出队列，因为异步队列为了延长数据生命周期通常必须复制整个网格和解向量，反而增加峰值内存。同步导出只在用户显式开启时增加运行时间，但不会改变默认求解路径。

回调有三个优点：

1. 可以流式写盘，不必把所有大规模解快照留在内存；
2. 单元测试可以用内存回调核对标记和父子关系，不受文件 I/O 影响。
3. 回调为空时只发生一次廉价条件判断，不构造快照数据、不计算派生场。

当前 `AdaptiveMeshHierarchy::refine(marked)` 还需要返回或填充一份
`AdaptiveRefinementReport`，使“用户标记”和“NVB 闭包追加”能够按单元 ID 区分，而不只是记录数量。

### 4.3 保持原计算效率的硬性设计规则

1. **默认关闭**：所有现有 benchmark 的默认行为、输出和数值路径保持不变；只有传入 `--snapshot-dir` 或 `--export-final` 才启用导出。
2. **不触发额外求解**：普通求解不为了画图自动计算 fine FEM reference、exact quadrature error、局部对偶指示子或 inf-sup；导出只能使用该实验原本已经计算的数据，除非用户显式请求附加诊断。
3. **字段按位选择**：使用 `ExportField` 位掩码选择 `mesh`、`indicator`、`solution`、`reference`、`basis` 等字段。未选择字段不分配、不计算、不写盘。
4. **只写权威数据**：C++ 对复值场默认只写实部和虚部；模、相位、对数值和色标在 Python 中派生，避免 C++ 再分配多个 \(O(N_h)\) 向量。
5. **流式序列化**：VTU writer 直接从 `TriMesh`、Eigen 向量和 `std::span` 逐数组写入文件，不构造完整 XML 字符串、DOM 或第二份 connectivity。
6. **不累计历轮状态**：自适应动画的数据逐轮落盘，`AdaptiveHelmholtzResult` 继续只保留标量 history 和最终轻量数据，不加入 `vector<mesh>`、`vector<solution>` 或 `vector<model>`。
7. **选列而非整矩阵复制**：基函数图从当前已有 `trial_basis`/corrector 数据中读取一个指定粗自由度或粗单元，导出后立即释放临时 \(O(N_h)\) 列向量。
8. **不改变实值模型内存策略**：实值 `LodModel` 继续丢弃临时校正子；MATLAB 对照优先只画现有 `uH/uHms`，不为了基函数图改变其常驻内存。
9. **I/O 与求解计时分离**：新增 `export_ms`，不把图片/VTU 写盘时间混入 `build_ms`、`solve_ms`、`estimate_ms`，论文性能图默认排除可视化 I/O。
10. **大数据默认不导出每轮解**：自适应默认只导出粗网格、层数、指示子和标记；解场默认只导出最终迭代。只有显式 `--snapshot-fields=solution` 才逐轮保存解。

性能验收门槛：

- 导出关闭时，峰值常驻内存不得因可视化功能出现与 \(N_h\)、粗单元数或迭代数相关的增长；
- 导出关闭时，Release benchmark 中位运行时间变化应处于测量噪声内，目标不超过 1%；
- 导出开启时，额外峰值内存应保持为 \(O(N_h+N_T)\)，且不得变成 \(O(n_{\mathrm{iter}}N_h)\)；
- writer 的临时工作区应有明确上界，优先逐块写出，不能一次复制所有 point/cell fields。

## 5. 数据格式和字段约定

### 5.1 每个运行目录

```text
results/visualization/<case-id>/
  run.json
  history.csv
  final.vtu
  adaptive/
    iteration_000_pre.vtu
    iteration_000_post.vtu
    iteration_001_pre.vtu
    ...
    adaptive.vtu.series
  basis/
    coarse_node_XXXX.vtu
  figures/
    ...
```

`case-id` 建议包含问题、关键参数和短哈希，例如：

```text
helmholtz_adaptive_k8_H2_h7_ell4_gaussian_a1b2c3d
```

### 5.2 `.vtu` 点数据

按需要写入以下字段：

| 字段 | 含义 |
|---|---|
| `u_lod_real`, `u_lod_imag` | 完整 LOD 解的实、虚部，默认权威输出 |
| `u_lod_abs`, `u_lod_phase` | Python 端派生字段；除调试外不要求 C++ 重复写出 |
| `u_reference_*` | 细网格参考解或 manufactured exact solution |
| `u_coarse_interp_*` | 未经校正的粗尺度系数延拓到细网格 |
| `u_fine_scale_*` | `u_lod - u_coarse_interp`，直接突出微振荡 |
| `error_*` | `u_reference - u_lod` |
| `basis_standard_*` | 标准粗基函数在细网格上的延拓 |
| `basis_corrector_*` | 对应局部校正子 |
| `basis_lod_*` | 修正后的 LOD 基函数 |

对于复值 Helmholtz，实部和虚部是权威数据；模和相位由 Python 端计算并核对
\(\lvert u\rvert=\sqrt{(\Re u)^2+(\Im u)^2}\)。

`reference`、`error`、`basis` 和 `fine_scale` 均为按需字段。普通自适应网格动画不写这些大向量；最终解任务才请求它们。`error` 和 `fine_scale` 优先在 Python 中由已导出的权威场相减得到，避免 C++ 额外分配同规模向量。

在 \(\lvert u\rvert\) 接近零的点，相位没有数值意义；绘图时应根据相对阈值屏蔽这些点，不能显示随机相位噪声。

### 5.3 `.vtu` 单元数据

| 字段 | 含义 |
|---|---|
| `element_id` | 跨迭代稳定的单元 ID |
| `parent_id` | 父单元 ID，无父节点时使用约定哨兵值 |
| `level` | 当前 NVB 层数 |
| `indicator`, `indicator_squared` | 当前选定的局部指示子 |
| `indicator_fine_squared` | fine 聚合指示子 |
| `indicator_mixed_squared` | mixed 聚合指示子 |
| `indicator_macro_squared` | macro 聚合指示子 |
| `eligible` | 是否允许在尺度分离约束下加细 |
| `marked` | 是否被 Dörfler 策略直接标记 |
| `closure_refined` | 是否仅因 NVB 一致性闭包而加细 |
| `new_child` | 是否在本次迭代新产生 |
| `diffusion` | 单元扩散系数 |
| `refractive_index` | Helmholtz 折射率/质量项系数 |
| `patch_membership` | 绘制指定校正子时是否属于 \(\ell\)-patch |

### 5.4 `history.csv`

保留现有自适应 CSV 字段，并增加：

- `coarse_dofs_total`；
- `coarse_dofs_free`；
- `fine_dofs_total`；
- `lod_trial_dofs`；
- `reference_kind = exact|fine_fem`；
- `energy_norm_kind = 1k|A_energy`；
- `coefficient_id` 和随机种子；
- 可选的峰值内存。

论文中“误差关于自由度”默认使用实际参与粗问题求解的自由度：

- 实值齐次 Dirichlet 问题：`coarse_dofs_free`，即去掉边界后的粗自由度；
- 当前齐次 impedance Robin Helmholtz：所有粗节点均是未知量，因此通常等于 `coarse_nodes`；
- hp 扩展：使用实际粗离散空间的维数，不能用节点数替代。

图注必须明确误差是相对精确解还是相对固定细网格参考解。二者不应画成同一含义的曲线。

### 5.5 `run.json`

至少记录：

- Git commit；
- 工作区是否 dirty；
- CMake build type；
- 编译器和版本；
- Python 及绘图库版本；
- 所有 PDE、LOD、网格、估计子和线程参数；
- 系数场种子或输入文件哈希；
- 运行时间和时区；
- 生成该目录的完整命令；
- 范数定义；
- 停止原因。

## 6. 三项核心可视化的实施方案

### 6.1 网格、标记与动态加细

#### 静态图

每次迭代至少生成一张三联图：

1. 当前粗网格，按 `level` 着色；
2. 指示子 \(\eta_T\) 或 \(\eta_T^2\) 热图；
3. 标记结果：直接标记单元用红色，闭包追加单元用橙色，新生成子单元用蓝色边界。

图中显示：

- 迭代编号；
- 粗单元数和实际粗自由度；
- Dörfler 参数 \(\theta\)；
- 标记/闭包数量；
- \(\eta\)、相对 \(\|e\|_{1,k}\)、\(q_{\max}\)。

指示子跨迭代比较时使用统一的对数色标范围，避免每一帧自动缩放造成“误差似乎没有下降”的视觉误导。

#### 动画

每一轮至少安排两个关键帧：

1. `pre`：显示指示子并高亮即将加细的单元；
2. `post`：显示新网格，区分直接标记产生的子单元和闭包产生的子单元。

输出：

- MP4：论文答辩、报告和网页使用；
- GIF：快速预览；
- `.vtu` 文件序列：ParaView/PyVista 交互检查；
- 每帧 PNG：便于定位问题和重新编码。

相机、坐标范围、画布尺寸、色标范围和字体在整个动画中固定。二维网格动画以 Matplotlib 为主，三维解随网格变化的动画再使用 PyVista。

#### 必要的代码改动

- 自适应 driver 增加逐步快照回调；
- hierarchy 返回实际闭包单元和父子 ID 报告；
- benchmark 增加 `--snapshot-dir`、`--snapshot-fields`、`--snapshot-every`；
- 保留现有 `--mesh-out` 作为轻量最终网格输出，避免破坏旧脚本。

### 6.2 相对 \(\|e\|_{1,k}\) 关于自由度

主图采用双对数坐标：

\[
x=N_{\mathrm{dof}},\qquad
y=\frac{\|u_{\mathrm{ref}}-u_{\mathrm{LOD}}\|_{1,k}}
        {\|u_{\mathrm{ref}}\|_{1,k}}.
\]

同时建议绘制：

- 相对 \(L^2\) 误差；
- uniform refinement 与 adaptive refinement 对比；
- 不同 \(\ell\) 的曲线；
- 不同 \(k\) 的曲线或分面图；
- manufactured exact error 与 fine-reference error，使用不同线型；
- fine FEM discretization floor；
- 与理论一致时再添加参考斜率，不能仅凭有限数据自动宣称理论阶。

每个点应带有或可追溯到：

- \(H,h,\ell,k,p\)；
- 粗/细自由度；
- 是否通过残差正确性阈值；
- 运行 ID。

绘图脚本需要：

- 拒绝在 log 轴上静默绘制零值、负值和 NaN；
- 可选择误差条或多随机系数样本的置信区间；
- 将拟合区间与斜率输出到单独 CSV，而不是只写在图上；
- 默认输出 PDF 和 600 dpi PNG；
- 对密集三角场图采用“场栅格化、文字和坐标轴矢量化”，控制 PDF 大小。

### 6.3 最终解与微小振荡细节

#### 复值 Helmholtz 问题

这是本计划的主线。第一批至少输出：

- \(\Re u_{\mathrm{LOD}}\)；
- \(\Im u_{\mathrm{LOD}}\)；
- \(\lvert u_{\mathrm{LOD}}\rvert\)；
- \(\arg u_{\mathrm{LOD}}\)；
- 参考解与 LOD 解的同色标对比；
- 实部或模的误差图；
- 粗尺度插值和细尺度校正图。

实部/虚部使用以零为中心的对称色标；模使用从零开始的顺序色标；相位使用循环色图并固定为 \([-\pi,\pi]\)。

其中粗尺度延拓直接用已经存在的 `coarse_to_fine * coarse_coefficients` 计算一次临时向量；完整解直接引用 `fine_values`；细尺度校正优先由 Python 用二者相减。C++ 不为这些图长期保存第二套解。

#### 实值椭圆问题（可选回归）

只有在需要核对 MATLAB 迁移结果时，才复现：

1. \(A(x)\)；
2. \(u_h\)；
3. \(u^{ms,h}_{H,\ell}\)；
4. \(P_Hu_H\) 或粗网格上的 \(u_H\)；
5. \(u^{ms,h}_{H,\ell}-P_Hu_H\)。

该任务只使用现有 `uH/uHms/P_node` 和该对照实验原本需要的参考解，不要求 `LodModel` 保存已经释放的局部校正子。

#### 防止振荡细节丢失的规则

- Helmholtz 直接使用已经求得的 `fine_values`，不重复执行 `trial_basis * coarse_coefficients`；
- 实值可选对照直接使用已经求得的 `uHms`；
- 不允许只把粗系数画在粗网格后再做平滑插值；
- 不对论文主图做网格 decimation；
- 2D 场图保留细网格三角剖分；
- 3D 表面按真实细网格顶点抬升；
- 增加一个微结构区域的局部放大框；
- 增加一条固定截线，例如 \(y=0.5\)，比较参考解、完整 LOD、粗尺度部分和细尺度修正；
- 系数场采用单元常值、不平滑显示；
- 解场可以在 P1 单元内线性插值，但不能使用会改变极值或相位的额外平滑；
- hp-P2/P3 情况不能只画顶点值，应按单元高阶基在足够密的可视化采样网格上求值后再绘图。

对于细微振荡，论文中优先使用“俯视等值填色 + 局部放大 + 截线”，三维透视图作为辅助。仅使用带光照的三维表面容易被相机、光照和纵横比掩盖细节。

## 7. 其他适合论文的可视化

下表按优先级排列。

### A. 论文主体建议必做

| 图 | 科学问题 | 建议形式 |
|---|---|---|
| 系数场/介质场 | 问题是否具有高对比、微尺度结构 | 单元常值热图 |
| 标准基、校正子、LOD 基 | LOD 如何注入细尺度信息 | 三联图及局部放大 |
| \(\ell\)-patch 支撑 | 校正子的局部化区域是什么 | 网格边界 + patch 着色 |
| 校正子衰减 | 局部化误差是否随 \(\ell\) 衰减 | 半对数图 |
| 参考解/LOD/粗尺度/细尺度修正 | 方法最终近似了什么 | 四联或五联图 |
| 点态误差或单元能量误差 | 误差集中在哪里 | 对数热图 |
| uniform vs adaptive | 自适应是否以更少 DOF 达到同误差 | work–precision 双对数图 |
| 自适应指示子与最终网格 | 网格是否追踪源、界面或波结构 | 指示子与网格并排 |
| effectivity | 估计子是否稳定 | effectivity–DOF 图 |
| fine FEM floor | LOD 误差是否已受参考离散误差限制 | 收敛图上的水平/独立曲线 |

### B. Helmholtz 论文建议

| 图 | 用途 |
|---|---|
| 实部、虚部、模、相位 | 完整描述复值波场 |
| 固定截线的振幅和相位 | 观察相位滞后、数值色散和局部振荡 |
| 误差关于 \(kH\) | 展示分辨率条件 |
| 固定 DOF 下误差关于 \(k\) | 展示 pollution/高波数困难 |
| inf-sup 常数关于 \(H,k,\ell\) | 展示离散稳定性 |
| 主/伴随校正基对比 | 解释 Petrov–Galerkin 的双侧修正 |
| 不同 \(\ell\) 的误差、残差和耗时 | 说明局部化精度—代价权衡 |
| P1/P2/P3 work–precision | 展示 hp 路径的收益和成本 |

### C. 自适应估计子研究建议

| 图 | 用途 |
|---|---|
| fine/mixed/macro 三类局部指示子并排 | 比较不同聚合策略 |
| 强残差指示子 vs 局部对偶指示子散点图 | 检查相关性 |
| Spearman 相关系数随迭代变化 | 判断排序稳定性 |
| Dörfler 集合能量重合率随迭代变化 | 检查标记一致性 |
| \(\eta/\|e\|\) 的 effectivity envelope | 支持可靠性/效率讨论 |
| 标记数与闭包新增数 | 衡量网格一致性闭包开销 |
| \(q_{\max}\)、\(q_{\mathrm{effective}}\) 随迭代 | 检查尺度分离约束 |
| exact error、fine-reference error、estimator 同图 | 避免把连续误差和离散误差混淆 |

### D. 实现与性能论文/附录可选

| 图 | 用途 |
|---|---|
| corrector/coarse/reference/estimate 分阶段耗时 | 找出主要成本 |
| 线程数—加速比与并行效率 | OpenMP 扩展性 |
| 峰值内存关于 DOF | 大规模可行性 |
| patch 局部 DOF 数分布 | 解释负载不均衡 |
| 粗矩阵、校正矩阵和基矩阵稀疏图 | 展示局部化产生的稀疏结构 |
| GMRES 残差历史 | 比较预条件器 |
| 误差—总时间而非只对 DOF | 给出实际 work–precision |

“稀疏矩阵图”和“性能图”应放在实现/求解器论文或附录中；如果论文重点是 LOD 数值分析，不宜挤占主体篇幅。

## 8. 统一的论文制图规范

### 8.1 色图

- 实部、虚部和有正负的误差：感知均匀的发散色图，零点居中；
- 模、系数、层数、误差指示子：感知均匀的顺序色图；
- 相位：循环色图；
- `marked`、`closure`、`new_child`：固定离散颜色；
- 避免彩虹色图；
- 同一比较组必须共享色标范围。

### 8.2 输出

- 曲线、网格边线、文字：PDF/SVG；
- 大型三角场：600 dpi PNG，或 PDF 中只栅格化场；
- 动画：H.264 MP4，另给低分辨率 GIF 预览；
- 保留原始帧，便于无损重编码；
- 字号、线宽和成图物理尺寸由 `style.py` 统一控制；
- 图中变量用 LaTeX 数学记号，图题保持简短，参数放图注或角标。

### 8.3 可复现性

每张最终图片应能由一条命令生成，例如：

```bash
python3 tools/visualization/reproduce_all.py \
  --manifest results/visualization/<case-id>/run.json \
  --output figures/paper
```

脚本不读取未记录的全局配置，不依赖当前工作目录之外的数据，不通过手工拖动色条或相机生成最终论文图。

## 9. 分阶段实施计划

### 阶段 0：确定科学口径（0.5–1 天）

交付：

- 固定 \(\|\cdot\|_{1,k}\)、\(A\)-能量范数、相对误差和 DOF 定义；
- 明确每张图使用 exact solution 还是 fine FEM reference；
- 确定第一批 Helmholtz 复现实验的 \(k,H,h,\ell,p\)、源项和网格参数；
- 确定论文主图的尺寸、字体和输出格式。

验收：

- 定义写入代码注释、CSV 字段说明和 `run.json`；
- 不再使用含义不清的 `energy_error` 图例。

### 阶段 1：通用导出层（2–3 天）

交付：

- 轻量 `.vtu` 写出器；
- CSV history 和 JSON manifest；
- Helmholtz 复值 point data 与自适应 cell data；
- 非拥有式只读视图、字段位掩码和流式 writer；
- 一个小网格回归样例；
- C++/Python schema 测试和“导出关闭”性能基线。

验收：

- PyVista 和 ParaView 均能读取；
- 顶点数、单元数、连接关系、面积和字段长度与 C++ 一致；
- 复数拆分/重组误差在机器精度量级；
- 默认关闭导出时不改变现有 benchmark 输出、不触发额外求解或大对象分配；
- Release 基线的关闭态耗时变化目标不超过 1%，峰值内存无规模相关增长。

### 阶段 2：Helmholtz 静态网格与最终解（1.5–2 天）

交付：

- 网格图、系数图；
- Helmholtz 实部、虚部、模和相位；
- 参考解、完整 LOD 解、粗尺度延拓、细尺度修正和误差图；
- 局部放大与固定截线图。

验收：

- 图中的 `fine_values` 与当前求解结果逐点一致；
- Python 派生的模、相位、误差和细尺度修正满足代数恒等式；
- 局部放大能看见 LOD 校正引入的微小振荡；
- 最终解导出结束后没有可视化数据继续滞留在 model 之外。

### 阶段 3：收敛与论文曲线（1–2 天）

交付：

- 相对 \(\|e\|_{1,k}\)–DOF；
- 相对 \(L^2\)–DOF；
- uniform/adaptive、不同 \(\ell\)、不同 \(k\) 的绘图入口；
- 拟合斜率和数据筛选报告。

验收：

- 每个图点能回查到 CSV 行和 run ID；
- Python 从导出字段重算的范数与 C++ CSV 相对差建议小于 \(10^{-10}\)；
- exact 和 fine-reference 曲线有清楚区分。

### 阶段 4：自适应快照和动画（2–3 天）

交付：

- 自适应 snapshot callback；
- refinement report；
- `pre/post` `.vtu` 序列；
- 网格、标记、闭包和指示子静态图；
- MP4/GIF 动画。

验收：

- 快照数与迭代历史一致；
- 直接标记集合满足当前 Dörfler 规则；
- `marked` 与 `closure_refined` 可分别复核；
- 新子单元父 ID 可追溯；
- 所有动画帧使用相同坐标和色标；
- 默认快照只包含粗网格和单元字段，逐轮解场必须显式开启；
- 峰值内存不随已完成迭代数线性增长。

### 阶段 5：基函数、校正子和复值解（2–3 天）

交付：

- 标准粗基、主/伴随校正子和修正基；
- patch 边界图；
- \(\ell\) 局部化衰减图；
- Helmholtz 实部、虚部、模、相位；
- 微结构局部放大和截线图。

实现注意：

- Helmholtz 模型已经暴露基和校正子数据；
- 只从 Helmholtz 模型现有数据中抽取用户指定的一个粗自由度/粗单元，不复制整套稀疏矩阵；
- 不修改实值 `LodModel` 的校正子释放策略；
- hp-P2/P3 另加单元内采样器，不能复用 P1 顶点绘图。

验收：

- `basis_lod = basis_standard - basis_corrector` 在容差内成立；
- `u_fine_scale = u_lod - u_coarse_interp` 在容差内成立；
- 局部 patch 外的截断校正子符合实现支持范围；
- 相位零模区域已屏蔽。

### 阶段 6：可选 MATLAB 回归与论文定稿（1–2 天）

交付：

- `reproduce_all.py`；
- 固定依赖版本；
- 图片清单和图注参数表；
- README 中的可视化运行说明；
- 最终 Helmholtz 论文图；
- 如论文或回归确有需要，再用现有 `uH/uHms` 复现 MATLAB 四联图，不增加实值模型常驻状态。

验收：

- 在干净 WSL 环境中能从已有数据一键重绘；
- 无 GUI 环境下静态图生成成功；
- 最终图不依赖手工 ParaView 状态；
- Windows 镜像中只提交计划内文件。

预计总工作量约 10–15 个有效开发日。若第一阶段只做“Helmholtz 网格 + 收敛图 + 最终 P1 解”，可先在约 5–8 日内形成可用于论文初稿的最小版本；基函数、完整自适应动画、实值 MATLAB 对照和 hp 高阶场采样放在第二批。

## 10. 测试与科学正确性要求

### 10.1 数据测试

- 每个三角形顶点索引合法且面积为正；
- 导出前后节点和单元数量完全一致；
- 稳定单元 ID 在未加细单元上跨迭代不变；
- 父子关系形成无环树；
- point data 长度等于节点/DOF 数；
- cell data 长度等于单元数；
- 不允许静默写出 NaN/Inf，确有科学含义的 NaN 必须只出现在表格并写明原因。

### 10.2 数值测试

- C++ 与 Python 对 \(L^2\)、\(\|\cdot\|_{1,k}\) 的结果一致；
- manufactured exact error 与现有 benchmark 一致；
- 细网格参考误差与 exact error 明确分开；
- LOD 重构满足 `fine_values = trial_basis * coarse_coefficients`；
- 实值重构满足 `uHms = G * uH`；
- 指示子总量等于局部平方量求和后的平方根；
- Dörfler 标记达到要求的能量比例。

### 10.3 图形测试

- 小样例能在 headless WSL 中生成；
- 图片非空且具有预期分辨率；
- 所有比较面板共享指定色标；
- 相位图使用循环色图并屏蔽零模区域；
- 动画相机不漂移；
- 不以逐像素完全一致作为主要回归标准，避免渲染库小版本变化导致脆弱测试；测试数据、范围、对象数量和输出文件更重要。

### 10.4 性能与内存回归

选择一个小型 correctness case 和一个代表性 Release case，分别测量：

1. 修改前/关闭导出的 wall time 与峰值 RSS；
2. 修改后/关闭导出的 wall time 与峰值 RSS；
3. 只导出最终解；
4. 每轮只导出网格与指示子；
5. 显式开启逐轮解场导出。

计时至少重复 5 次，比较中位数；峰值内存使用 `/usr/bin/time -v` 或等价工具记录。关闭导出是必须通过的回归门，开启导出的时间只单独报告，不混入求解器性能结论。若关闭态出现超过 1% 的稳定回退或随问题规模增长的额外内存，先修复架构，不继续增加图片功能。

## 11. 风险与处理

| 风险 | 处理方式 |
|---|---|
| 细网格 VTU 和逐帧图片很大 | 默认只保存粗网格字段和最终解；提供 `--snapshot-every` 与字段位掩码；大文件不进 Git；后期使用二进制压缩 |
| Python/PyVista 在无显示服务器环境失败 | 2D 主流程用 Matplotlib `Agg`；PyVista 使用 off-screen；必要时安装 Mesa/EGL 或在 WSLg 下检查 |
| 不同迭代自动色标导致错误视觉比较 | 动画预扫描全局范围或由 manifest 固定范围 |
| 相位在零模区域出现噪声 | 按相对模阈值 mask |
| 粗节点数被错误当作所有问题的 DOF | 显式导出 total/free/trial DOF |
| exact error 与 fine-reference error 混用 | 字段名、线型和图注强制区分 |
| 高阶解只画顶点导致细节错误 | hp 单独实现单元内采样 |
| 保存所有校正子导致内存增加 | 不新增保存逻辑；只借用 Helmholtz 已有数据导出指定列；实值模型维持释放策略 |
| 异步 I/O 队列复制大向量 | 第一版使用同步回调和非拥有式视图，不建立后台队列 |
| 为图片额外计算 reference/diagnostics | 默认禁止；只有用户显式请求附加诊断时执行，并单独计时 |
| Windows 与 WSL 镜像不同步 | 计划书和后续源码同步后做 SHA-256 校验；输出记录 commit/dirty 状态 |
| 图很好看但无法证明数值正确 | 所有绘图字段与 C++ 数值回归测试绑定 |

## 12. 推荐的第一批论文图

为了尽快形成论文可用结果，建议按以下顺序制作：

1. Helmholtz 介质系数场和初始粗/细网格；
2. Helmholtz 参考解、完整 LOD 解、粗尺度延拓和细尺度修正；
3. Helmholtz 解的实部、虚部、模和相位；
4. 细尺度修正的局部放大和固定截线；
5. 相对 \(\|e\|_{1,k}\) 与相对 \(L^2\) 误差关于实际粗自由度；
6. uniform 与 adaptive 的误差–DOF 对比；
7. 自适应某三次代表性迭代的“指示子 + 标记 + 网格”静态组合图；
8. 最终自适应网格；
9. 标准粗基、主/伴随校正子和修正 LOD 基；
10. effectivity 和 \(q_{\max}\) 随迭代变化；
11. \(\ell\) 对误差、校正子衰减和耗时的影响；
12. 附录中的计时分解和并行加速图；
13. 需要迁移回归时再加入实值 MATLAB 四联图。

这组图能够同时回答：问题有多难、LOD 如何工作、细尺度信息在哪里、方法是否收敛、自适应是否有效、Helmholtz 计算是否稳定，以及实现成本如何。

## 13. 最终决策

本项目的推荐可视化路线是：

> C++20 保持为唯一数值真值来源；Helmholtz 可视化只借用求解过程中已经存在的数据并按需流式导出；用轻量 VTU/CSV/JSON 交换数据；在 WSL 中用 Python/Matplotlib 生成论文静态图，用 PyVista/VTK 生成三维图和动画；ParaView 只用于检查；Windows 镜像只承担版本控制和 GitHub 推送。

实施优先级为：

1. 先建立导出关闭态的时间和峰值内存基线；
2. 统一 Helmholtz 数据、范数和自由度口径；
3. 完成非拥有式、字段按需、流式导出层；
4. 完成 Helmholtz 最终解、细尺度修正和误差–DOF 图；
5. 完成逐步自适应网格快照和动画；
6. 完成主/伴随基函数、校正子和论文一键复现；
7. 最后按需完成不改变实值模型内存策略的 MATLAB 对照。

该路线能满足当前三项要求，同时把“默认无额外求解、无大对象复制、无历轮状态累计、实值模型不保留临时校正子”作为验收门，从而不牺牲原有计算效率，并为后续自适应、hp、Helmholtz pollution 和性能论文保留扩展空间。
