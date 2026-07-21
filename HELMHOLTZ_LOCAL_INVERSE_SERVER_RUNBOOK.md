# Helmholtz 局部逆不等式服务器实验运行手册

> 目标机器：AMD EPYC 9554，64 个物理核心，377 GiB 内存
>
> 实验对象：齐次 Robin 边界 Helmholtz LOD、全局协调非拟一致 NVB 粗网格、固定全局细空间
>
> 主指标：oversampling 匹配分母 `omega_T^ell`，默认 `ell=3`

## 1. 本次性能优化

`bench_helmholtz_local_inverse` 已做以下优化，同时保留默认完整模式：

1. 各粗单元上的局部广义特征值问题使用 OpenMP `dynamic,1` 调度；大 patch 先后差异由动态调度吸收。
2. `--basis=trial` 只分析 trial LOD 空间。当前实系数下 test 与 trial 共轭，完整校准仍会验证共轭残差。
3. `--denominators=matched` 只计算 `omega_T^ell` 分母；`element-matched` 同时计算同单元和匹配分母。
4. 服务器脚本固定 `OMP_PLACES=cores`、`OMP_PROC_BIND=spread`、`OMP_DYNAMIC=FALSE`，避免线程迁移和动态缩减。
5. 所有长任务使用独立 `.done` 文件。只有成功结束才写入该文件，重新运行会跳过已完成案例。

本地 WSL 的同一 `L_h=12` 校准网格上，完整模式用时从单线程 `7.86 s` 降到 8 线程 `2.84 s`；只保留 trial/matched 后为 `2.40 s`。局部逆分析阶段从约 `154 ms` 降到 `27 ms`，匹配分母数值逐位一致。服务器必须重新做线程试跑，不能直接外推本地加速比。

## 2. 获取代码和安装依赖

在服务器执行：

```bash
git clone https://github.com/codenotforce/LOD2d-CPP.git
cd LOD2d-CPP
```

Ubuntu/Debian 依赖：

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake ninja-build \
  libeigen3-dev libsuitesparse-dev libtbb-dev time
```

建议先记录机器拓扑：

```bash
lscpu
free -h
numactl --hardware 2>/dev/null || true
```

不要在同一节点上同时启动多个深细网格案例。脚本内部对粗单元并行，而 `h` 扫描案例之间保持串行，以控制峰值内存。

深 `h` 主实验必须使用 `--denominators=matched`。同单元分母的局部质量矩阵会随 `h/H_T` 变小而出现近零方向；它已在本地实验中表现出明显的阈值敏感性，不应让这一诊断量阻断 oversampling 匹配的 `patch3` 主实验。`MODE=all` 中的低成本完整校准仍使用 `all/all` 检查同单元装配。

## 3. 编译与校准

脚本会自动用 Release、`-O3 -march=native` 和 OpenMP 编译，并运行 `helmholtz_local_inverse_smoke`：

```bash
chmod +x scripts/run_helmholtz_local_inverse_server.sh
MODE=pilot BUILD_JOBS=64 PILOT_THREADS="8 16 32 64" \
  bash scripts/run_helmholtz_local_inverse_server.sh
```

比较：

```bash
grep -H -E "Elapsed|Maximum resident" \
  results/helmholtz_local_inverse_server/time_pilot_t*.txt
```

主任务线程数取 wall time 最短者。默认建议从 `THREADS=32` 开始，因为目标反馈网格通常只有约 16–30 个 coarse patches；64 个物理核心并不意味着64线程一定更快。若32和64线程相近，选32以减少每线程稀疏分解的并发内存。

## 4. 固定粗网格的深 h 扫描

先在同一个等级范围约 `3:8` 的粗网格上只改变全局细层，隔离 `h` 误差：

```bash
THREADS=32 MODE=hscan HSCAN_LEVELS="14 15 16" HSCAN_STEPS=5 \
  bash scripts/run_helmholtz_local_inverse_server.sh
```

对每个粗单元定义

```math
q_T=\max_{t\subset T}\frac{h_t}{H_T},
\qquad q_{\max}=\max_T q_T.
```

主局部逆常数为

```math
C_{\mathrm{inv},\ell}(h)
=\max_{T\in\mathcal T_H}
H_T\sup_{0\ne v\in V_{H,\ell}^{\mathrm{ms}}}
\frac{\|\nabla v\|_{L^2(T)}}
     {\|v\|_{L^2(\omega_T^\ell)}}.
```

`L_h=14,15,16` 预期依次减小 `q_max`；实际比值必须以 CSV 的 `q_max` 为准。最低验收要求：

- 至少最后一个点满足 `q_max <= 0.0625`，即所有粗单元上 `h_T/H_T <= 1/16`；
- 最后两级 `patch3` 最大值相对变化不超过 `2%`；
- `fixed_fine_mesh=1`，`nesting_residual <= 1e-10`；
- Petrov、corrector、constraint 和局部特征值残差通过程序内置 `--check`。

服务器主脚本只输出 `patch3` 行。如果还要研究同单元诊断，应另行在较浅的 `L_h` 上运行 `--denominators=element-matched`，并单独解释质量条件数，不能把它与 `patch3` 的通过/失败混在一起。

如果 `L_h=16` 的峰值 RSS 小于 `120 GiB` 且最后两点尚未进入平台，可以追加：

```bash
THREADS=32 MODE=hscan HSCAN_LEVELS="17" HSCAN_STEPS=5 \
  bash scripts/run_helmholtz_local_inverse_server.sh
```

不要在没有查看 `time_hscan_Lh16.txt` 的情况下直接上 `L_h=17`。

## 5. 最大值反馈加细主实验

使用当前 `patch3` 最大值所在单元驱动下一次 NVB 加细，只标记目标单元，NVB 自动完成协调闭包：

```bash
THREADS=32 MODE=feedback FINE_LEVEL=16 STEPS=6 ELL=3 \
  bash scripts/run_helmholtz_local_inverse_server.sh
```

每轮执行

```math
T_j^*=\mathop{\arg\max}_{T\in\mathcal T_{H,j}}Q_{T,\ell},
\qquad
\mathcal T_{H,j+1}=\operatorname{NVBclose}
\bigl(\operatorname{refine}(\mathcal T_{H,j},\{T_j^*\})\bigr).
```

检查以下轨迹：

1. `L_max-L_min` 和 `grading_ratio=H_max/H_min` 是否总体增大；
2. `neighbor_ratio` 是否保持 NVB 局部形状正则的受控值；
3. `patch3` 最大值是否在 `q_max <= 1/16` 后仍随 grading ratio 持续增长；
4. 最大值位置、等级和边界/过渡区标记是否迁移；
5. 每轮 `fixed_fine_mesh=1` 且 `V_H subset V_h` 的四重嵌套检查通过。

数值实验只能说明已测参数范围内是否观察到增长，不能把平台直接写成网格无关常数的证明。

## 6. 一次执行全部实验

线程试跑后若确定32线程最优，可执行：

```bash
THREADS=32 BUILD_JOBS=64 MODE=all \
  HSCAN_LEVELS="14 15 16" FINE_LEVEL=16 STEPS=6 \
  bash scripts/run_helmholtz_local_inverse_server.sh
```

建议在 `tmux` 或批处理系统中运行。脚本不会自动并行多个大案例。中断后执行相同命令即可，存在 `.done` 的成功案例会被跳过；需要强制重跑时设置 `RUN_FORCE=1`。

## 7. 结果回传

服务器结果位于：

```text
results/helmholtz_local_inverse_server/
```

需要回传的文件包括：

- `summary_*.csv`：主数值表；
- `elements_*.csv`、`mesh_*.csv`：反馈实验逐单元值和网格轨迹；
- `time_*.txt`：wall time 与峰值 RSS；
- `run_*.log`：完整命令和错误；
- `server_metadata.txt`：CPU、内存、Git 版本和运行参数；
- `*.done`：成功完成标记。

回传后再把服务器结果追加到 `DEVELOPMENT.md`。不要覆盖当前本地 WSL 结果目录。

## 8. 深细网格出现 local Hermitian check 失败时

旧版服务器脚本曾在 `hscan` 和 `feedback` 中使用 `element-matched`。若日志命令包含

```text
--denominators=element-matched
```

请先 `git pull` 更新脚本。失败案例不会产生 `.done` 文件，直接重新执行原来的 `MODE=hscan` 或 `MODE=feedback` 命令即可；脚本会覆盖该案例的残缺 CSV。新版报错同时给出 iteration、basis、denominator、Hermitian defect、eigen residual、energy identity error 和最大质量条件数，便于区分主 `patch3` 失败与同单元诊断病态。
