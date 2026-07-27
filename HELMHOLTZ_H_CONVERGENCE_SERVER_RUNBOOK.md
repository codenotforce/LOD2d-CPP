# Helmholtz LOD 绝对能量范数收敛实验：EPYC 服务器操作手册

## 1. 实验目标

在固定的全局 NVB 细网格 \(V_h\) 上逐次全局加细 \(V_H\)，比较

\[
E_{\mathrm{LOD}}(H)=
\left(\|\nabla(u-u_{H,h,\ell})\|_{L^2}^2+
k^2\|u-u_{H,h,\ell}\|_{L^2}^2\right)^{1/2}
\]

和同一粗网格上的 P1 Galerkin 解

\[
E_{\mathrm{P1}}(H)=
\left(\|\nabla(u-u_H)\|_{L^2}^2+
k^2\|u-u_H\|_{L^2}^2\right)^{1/2}.
\]

这里是绝对能量误差，不除以 \(\|u\|_{1,k}\)。相邻点收敛率为

\[
p_j=\frac{\log(E(H_{j-1})/E(H_j))}
{\log(H_{j-1}/H_j)}.
\]

程序把 \(u_H\) 延拓到同一个 \(V_h\) 后，与 LOD 使用相同的细网格积分。
每个点还检查 \(P_{hH}x_H=x_h\)、\(P_{hH}y_H=y_h\)，保证
\(V_H\subset V_h\)。

## 2. 默认正式参数

默认正式实验是：

```text
k = 32
L_h = 21
L_H = 10 11 12 13 14 15
ell = 5
two-sided Petrov-Galerkin
```

实际几何量约为：

```text
h_max = 2^-10
kh = 1/32
kH = sqrt(2), 1, 1/sqrt(2), 1/2, 1/(2sqrt(2)), 1/4
```

共有六个收敛点；最细的粗网格仍满足 \(h/H=1/8\)，前面的点具有更强的
尺度分离。不要把 NVB 层号直接当作 \(2^{-L}\)；绘图使用 CSV 中实测的
`H_max` 或 `coarse_nodes`。

## 3. 获取代码与准备环境

```bash
cd ~/code/LOD2d-CPP
git status --short
git switch main
git pull --ff-only
git rev-parse HEAD
```

基础工具：

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake libeigen3-dev \
  libsuitesparse-dev numactl time gawk
```

`numactl` 不是硬要求。脚本默认 `NUMA_POLICY=auto`：存在 `numactl` 时采用
interleave，不存在时自动退回默认 NUMA 策略，不会像旧脚本那样直接停止。

```bash
chmod +x scripts/run_helmholtz_H_convergence_server.sh
```

## 4. 第一步：小规模正确性检查

```bash
cmake -S . -B build-H-convergence-server \
  -DCMAKE_BUILD_TYPE=Release \
  -DLOD2D_USE_OPENMP=ON \
  -DLOD2D_BUILD_TESTS=ON \
  -DLOD2D_BUILD_BENCHMARKS=ON

cmake --build build-H-convergence-server \
  --target bench_helmholtz_H_convergence -j 32

ctest --test-dir build-H-convergence-server \
  -R helmholtz_H_convergence_smoke --output-on-failure
```

必须看到测试通过后再做 pilot。

## 5. 第二步：线程数和局部求解器 pilot

先比较 16/32/64 线程以及 DirectSaddle/DirectSchur：

```bash
MODE=pilot \
PILOT_THREADS="16 32 64" \
PILOT_SOLVERS="saddle schur" \
K=32 \
PILOT_H_LEVEL=12 \
PILOT_FINE_LEVEL=18 \
PILOT_ELL=5 \
NUMA_POLICY=auto \
bash scripts/run_helmholtz_H_convergence_server.sh
```

查看时间和内存：

```bash
grep -H -E "Elapsed|Percent of CPU|Maximum resident" \
  results/helmholtz_H_convergence_server/time_pilot_*.txt
```

查看数值诊断：

```bash
for f in results/helmholtz_H_convergence_server/summary_pilot_*.csv; do
  echo "===== $f"
  column -s, -t < "$f" | less -S
done
```

选择原则：

1. 所有候选必须通过 `--check`，没有 `.done` 的运行不算结果；
2. 优先最短 wall time；
3. 若速度接近，选择较小 peak RSS；
4. 正式深网格默认从 32 线程开始，64 线程只有在 pilot 明显更快且内存增量可接受时使用；
5. 不要把 DirectSchur 固定当作必然更快；它在大尺度间隔上可能占优，但此前较小间隔实验中 DirectSaddle 更快。

pilot 完成后，把选出的 `THREADS` 和 `SOLVER` 明确写入正式命令。

## 6. 第三步：正式 \(k=32\) 六点实验

假设 pilot 选择 32 线程和 DirectSchur：

```bash
MODE=main \
THREADS=32 \
SOLVER=schur \
K=32 \
FINE_LEVEL=21 \
H_LEVELS="10 11 12 13 14 15" \
ELL=5 \
FINE_REFERENCE_LEVEL=15 \
EXPORT_FIELDS=0 \
NUMA_POLICY=auto \
MIN_AVAILABLE_GIB=200 \
bash scripts/run_helmholtz_H_convergence_server.sh
```

如果 pilot 选择 DirectSaddle，就改成 `SOLVER=saddle`。脚本每个 \(H\) 启动
一个独立进程，并写 `.done` 标记；中断后重跑同一命令会跳过已完成点。
只有需要覆盖完成点时才使用：

```bash
RUN_FORCE=1 MODE=main ... \
  bash scripts/run_helmholtz_H_convergence_server.sh
```

不要在没有保存旧结果的情况下使用 `RUN_FORCE=1`。

正式结果在：

```text
results/helmholtz_H_convergence_server/
  all_results.csv
  summary_main_*.csv
  run_main_*.log
  time_main_*.txt
  main_*.done
  server_metadata.txt
  snapshots/
```

`all_results.csv` 按 \(H\) 排列，并重新计算跨进程的四组相邻收敛率。

## 7. 网格与可视化数据

即使 `EXPORT_FIELDS=0`，每层也保存：

```text
H_L##_coarse_nodes.csv
H_L##_coarse_elements.csv
```

节点表含坐标、同网格 P1 解和 LOD 粗系数；单元表含连接关系、重心、面积和
直径。固定细网格只保存一次：

```text
fine_L21_nodes.csv
fine_L21_elements.csv
```

论文收敛图直接读取 `all_results.csv`：

- 横轴：`H_max`（log-log）或 `coarse_nodes`；
- 纵轴：`p1_energy_abs`、`lod_energy_abs`；
- 图例必须写 “absolute \(k\)-weighted energy error”；
- 不得把这些列标成 relative error；
- 可用 CSV 中的 `p1_energy_rate`、`lod_energy_rate` 标注相邻斜率；
- 论文最终斜率建议另对选定渐近区间做 log-log 最小二乘拟合，并在图注写清拟合点。

若要画解场、粗尺度分量和校正场，可重新执行并启用：

```bash
EXPORT_FIELDS=1 MODE=main ... \
  bash scripts/run_helmholtz_H_convergence_server.sh
```

这会为每个 \(H\) 增加 `H_L##_fields.csv`。level 21 上六份场文件很大，
只画收敛曲线时不要开启。所有大结果留在服务器或归档盘，不提交 Git。

## 8. 结果验收

```bash
column -s, -t < \
  results/helmholtz_H_convergence_server/all_results.csv | less -S

grep -H -E "Elapsed|Maximum resident" \
  results/helmholtz_H_convergence_server/time_main_*.txt

find results/helmholtz_H_convergence_server \
  -name '*.tmp' -o -name '*.done' | sort
```

验收要点：

- 六行正式数据和六个 `.done`；
- `nesting_coordinate_residual < 1e-13`；
- `p1_residual`、`petrov_residual`、`corrector_residual`、
  `constraint_residual < 1e-8`；
- DirectSchur 时 `schur_residual < 1e-8`、`schur_rcond > 1e-14`；
- `fine_energy_abs` 只在 `FINE_REFERENCE_LEVEL=15` 的行有值；
- LOD/P1 绝对能量误差总体下降；
- 结论只覆盖实际完成的 \(H\) 区间，不根据 pilot 作论文结论。

检查磁盘：

```bash
du -sh results/helmholtz_H_convergence_server
du -sh results/helmholtz_H_convergence_server/snapshots/*
```

## 9. 可选 \(k=64\) 升级实验

只有在 \(k=32\) 六点完成、峰值内存和时间均可接受后再执行。建议：

```bash
RESULT_DIR=results/helmholtz_H_convergence_server_k64 \
MODE=main \
THREADS=32 \
SOLVER=schur \
K=64 \
FINE_LEVEL=22 \
H_LEVELS="12 13 14 15 16" \
ELL=6 \
FINE_REFERENCE_LEVEL=16 \
EXPORT_FIELDS=0 \
NUMA_POLICY=auto \
MIN_AVAILABLE_GIB=220 \
bash scripts/run_helmholtz_H_convergence_server.sh
```

这组参数有五个点，`kh=1/(16sqrt(2))`，最细粗网格仍有 `h/H=1/8`。它可能比
\(k=32\) 正式实验昂贵很多；不要直接把它与 \(k=32\) 并行运行。若某点出现
OOM、外部中断或检查失败，该点不是数值结果，不要把残留 `.tmp` 文件并入图。

## 10. 打包回传

```bash
tar -czf helmholtz_H_convergence_server_results.tar.gz \
  results/helmholtz_H_convergence_server

sha256sum helmholtz_H_convergence_server_results.tar.gz \
  > helmholtz_H_convergence_server_results.tar.gz.sha256
```

若包含 `EXPORT_FIELDS=1` 的大场文件，压缩可能很慢；也可以只回传
`all_results.csv`、各 summary/log/time、metadata 及粗网格 CSV，把固定细网格
和场文件单独归档。
