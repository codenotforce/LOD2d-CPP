# Helmholtz 污染效应服务器实验操作手册

> 目标服务器：AMD EPYC 9554，64 个物理核心，377 GiB 内存
> 实验对象：齐次阻抗 Robin 边界、全局兼容 NVB、人工制造解
> 主实验：固定 `kH=1` 和 `kh=1/8`，考察 `k=32,64,128`

## 1. 数值目标与判据

人工解为

```math
u(x,y)=\phi(x)\phi(y)e^{ikx},
\qquad
\phi(t)=16t^2(1-t)^2,
```

并令

```math
f=-\Delta u-k^2u.
```

由于 `phi` 和 `phi'` 在端点同时为零，人工解满足

```math
\partial_nu-iku=0
\qquad\text{on }\partial\Omega.
```

主误差和参考解误差分别为

```math
E_{\mathrm{exact}}(v;k)
=\frac{\|u-v\|_{1,k}}{\|u\|_{1,k}},
\qquad
E_{\mathrm{ref}}(v;k)
=\frac{\|u_h-v\|_{1,k}}{\|u_h\|_{1,k}},
```

其中

```math
\|w\|_{1,k}^2
=\|\nabla w\|_{L^2(\Omega)}^2
+k^2\|w\|_{L^2(\Omega)}^2.
```

同时记录细网格参考误差

```math
E_h(k)=\frac{\|u-u_h\|_{1,k}}{\|u\|_{1,k}}.
```

在固定 `kH` 时，如果标准粗 P1 FEM 的 `E_exact` 和 `E_ref` 随 `k`
增长，而 `E_h` 不增长，则增长不能归因于参考解变差，可解释为数值污染。
LOD 只有在 `E_exact` 和 `E_ref` 都保持有界或下降时，才能称为在已测范围内
没有观察到污染增长。

## 2. 本次性能路径

波数扫描的主耗时是局部 corrector，而不是细网格参考解。优化后的路径包括：

1. trial corrector 只保存和装配一次。当前系数、Robin 参数和拟插值均为实数，
   test 基函数直接由 trial 基函数逐项共轭得到，仍逐 patch 检查伴随残差。
2. benchmark 显式接受 `--threads`，不再只依赖调用者是否正确设置 OpenMP 环境。
3. 保留 `DirectSaddle` 作为 P1 污染实验默认求解器。本地
   `k=16,H=9,h=15,ell=4` 校准中，DirectSaddle 为约 `19.7 s/1.56 GiB`，
   DirectSchur 为约 `33.7 s/1.47 GiB`；误差一致，但 Schur 慢约 70%。
4. 增加有上限的多槽符号缓存和“矩阵数值完全相同才复用分解”的安全选项。
   本地测试中单槽且不复用仍是最佳默认值；服务器脚本会重新试跑，不能把
   本地结论直接外推到 EPYC。
5. 服务器构建默认启用 Release、`-O3 -march=native` 和 LTO；进程固定到
   物理核心，并把大规模不同 `k` 案例串行执行。
6. 默认使用 `numactl --interleave=all`，避免主线程建立的全局网格和稀疏结构
   全部落在单个 NUMA 节点。仍必须用 pilot 与默认 NUMA 策略比较后再决定。

优化没有改变数学离散。`--check` 同时检查 Petrov、primal、adjoint、
插值约束和人工解指标；Schur 校准还检查 Schur 残差、条件估计和 fallback。

## 3. 预期网格规模

主实验使用 `fine-gap=6`，即固定 `kh=1/8`：

| `k` | `L_H` | `L_h` | `ell=ceil(log2 k)` | coarse elements | fine elements |
|---:|---:|---:|---:|---:|---:|
| 32 | 11 | 17 | 5 | 4,096 | 262,144 |
| 64 | 13 | 19 | 6 | 16,384 | 1,048,576 |
| 128 | 15 | 21 | 7 | 65,536 | 4,194,304 |

严格参考扫描使用 `fine-gap=8`，即固定 `kh=1/16`：

| `k` | `L_H` | `L_h` | fine elements |
|---:|---:|---:|---:|
| 32 | 11 | 19 | 1,048,576 |
| 64 | 13 | 21 | 4,194,304 |

第一次服务器实验不要直接运行 `k=256`。它需要 `L_H=17,L_h=23`，
细网格约 1678 万个三角形，并且 `ell=8`；内存和 corrector 总成本不能由
377 GiB 总内存这一项单独保证。

## 4. 获取代码与依赖

在服务器执行：

```bash
cd ~/code/LOD2d-CPP
git pull --ff-only
git status --short
git rev-parse HEAD
```

工作区应保持干净，并记录实际 commit。Ubuntu 22.04 依赖：

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake ninja-build \
  libeigen3-dev libsuitesparse-dev libtbb-dev numactl time
```

记录硬件和 NUMA 拓扑：

```bash
lscpu
free -h
numactl --hardware
```

不要同时启动其他深网格实验。服务器脚本默认要求启动每个案例前至少还有
256 GiB `MemAvailable`；该检查只是启动门槛，不是运行期间的内存上限。

## 5. 第一步：线程和缓存 pilot

```bash
cd ~/code/LOD2d-CPP
chmod +x scripts/run_helmholtz_pollution_server.sh

MODE=pilot \
PILOT_THREADS="8 16 32 64" \
THREADS=32 PILOT_K=32 \
NUMA_POLICY=interleave \
bash scripts/run_helmholtz_pollution_server.sh
```

脚本会测试：

- 单槽、无数值分解复用下的 `8,16,32,64` 线程；
- 32 线程下 `4` 槽符号缓存；
- 32 线程下 `4/8` 槽和完全相同矩阵分解复用。

比较 wall time、CPU 利用率和峰值内存：

```bash
grep -H -E "Elapsed|Percent of CPU|Maximum resident" \
  results/helmholtz_pollution_server/time_pilot_*.txt
```

检查所有数值行：

```bash
column -s, -t \
  < results/helmholtz_pollution_server/all_results.csv | less -S
```

选择 wall time 最短且内存合理的设置。不要只依据 CPU 百分比。缓存槽数增加会
保留更多 SparseLU 符号/数值对象；若提速不足约 10%，使用保守默认：

```bash
THREADS=32
SYMBOLIC_CACHE_SLOTS=1
FACTORIZATION_REUSE=none
```

如果 `16` 与 `32/64` 线程时间接近，优先 `16`，因为大型 patch 的并行
SparseLU 峰值内存近似随并发 worker 数增长。

NUMA 策略也应在 pilot 上比较一次：

```bash
RESULT_DIR=results/helmholtz_pollution_server_numa_default \
MODE=pilot PILOT_THREADS="32" THREADS=32 \
NUMA_POLICY=default \
bash scripts/run_helmholtz_pollution_server.sh
```

只有默认 NUMA 明显更快时才在主实验改为 `NUMA_POLICY=default`。

## 6. 第二步：主扫描，先到 `k=64`

建议在 `tmux` 或批处理系统中运行：

```bash
tmux new -s helmholtz-pollution

MODE=main K_VALUES="32 64" \
THREADS=32 SYMBOLIC_CACHE_SLOTS=1 FACTORIZATION_REUSE=none \
NUMA_POLICY=interleave MIN_AVAILABLE_GIB=256 \
bash scripts/run_helmholtz_pollution_server.sh
```

把 `THREADS`、缓存和 NUMA 参数替换成 pilot 的胜出配置。中断后执行同一命令，
已有 `.done` 的案例会跳过；只有成功写完一行 CSV 并通过 `--check` 才产生
`.done`。

`k=64` 完成后先查看：

```bash
grep -H -E "Elapsed|Maximum resident" \
  results/helmholtz_pollution_server/time_main_k*.txt
free -h
```

经验性升级门槛：

- `k=64` 峰值 RSS 小于 70 GiB：可以继续 `k=128`；
- 70--100 GiB：先把 worker 数降到 `16`，单独运行 `k=128`；
- 大于 100 GiB、发生 swap、或节点同时有其他任务：不要运行 `k=128`。

继续 `k=128`：

```bash
MODE=main K_VALUES="128" \
THREADS=32 SYMBOLIC_CACHE_SLOTS=1 FACTORIZATION_REUSE=none \
NUMA_POLICY=interleave MIN_AVAILABLE_GIB=300 \
bash scripts/run_helmholtz_pollution_server.sh
```

运行期间建议在另一个终端监控：

```bash
watch -n 10 'free -h; ps -C bench_helmholtz_k -o pid,etime,%cpu,%mem,rss,vsz'
```

一旦出现持续 swap 或系统进入内存压力，应终止当前单个 benchmark；不要删除
已经完成案例的 `.done` 和 CSV。

## 7. 第三步：严格参考网格复核

先只运行 `k=32`：

```bash
MODE=strict STRICT_K_VALUES="32" \
THREADS=32 SYMBOLIC_CACHE_SLOTS=1 FACTORIZATION_REUSE=none \
NUMA_POLICY=interleave \
bash scripts/run_helmholtz_pollution_server.sh
```

确认资源可接受后再运行：

```bash
MODE=strict STRICT_K_VALUES="64" \
THREADS=32 SYMBOLIC_CACHE_SLOTS=1 FACTORIZATION_REUSE=none \
NUMA_POLICY=interleave MIN_AVAILABLE_GIB=300 \
bash scripts/run_helmholtz_pollution_server.sh
```

严格扫描的用途是检查主扫描的 LOD 和参考解误差是否对 `V_h` 敏感。标准粗
FEM 解只依赖 `V_H`，其相对人工解误差在两种 `fine-gap` 下应只出现积分级别
差异。

## 8. 验收和停止条件

每行至少检查：

1. `kH=1`；主扫描 `kh=0.125`，严格扫描 `kh=0.0625`。
2. `source=manufactured`、`solver=saddle`、`direct_fallbacks=0`。
3. Petrov、corrector 和 constraint residual 均小于 `1e-8`。
4. `fine_exact_energy_rel` 足够小且不随 `k` 明显增长。
5. FEM 的 `fem_exact_energy_rel` 与 `fem_energy_rel` 是否同时随大 `k` 增长。
6. LOD 的 `lod_exact_energy_rel` 与 `lod_energy_rel` 是否保持有界。

如果精确 FEM 误差增长、参考解相对 FEM 误差也增长，而
`fine_exact_energy_rel` 不增长，则是可信的污染证据。如果只有相对参考解
误差增长而细参考误差也快速增长，应先加深 `V_h`，不能直接下污染结论。

任何残差检查失败、CSV 只有表头、日志出现 `bad_alloc`、进程被 OOM killer
终止，均视为该案例未完成。脚本不会写 `.done`，修正参数后可直接重跑。

## 9. 结果目录和回传

默认结果目录：

```text
results/helmholtz_pollution_server/
```

需要回传：

- `all_results.csv`：所有成功案例合并表；
- `summary_*.csv`：逐案例原始数值行；
- `time_*.txt`：wall time、CPU 和峰值 RSS；
- `run_*.log`：完整命令及错误；
- `server_metadata.txt`：Git、CPU、内存、NUMA 和运行参数；
- `*.done`：断点续跑状态。

打包命令：

```bash
tar -czf helmholtz_pollution_server_results.tar.gz \
  results/helmholtz_pollution_server
```

回传后再把最终大波数数据和结论写入 `DEVELOPMENT.md`。在服务器结果回来前，
不能把 `k<=32` 的本地结论外推成任意大波数下无污染。
