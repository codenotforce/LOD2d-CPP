# 最新 Algorithm 1/2：E1/E2 服务器操作手册（2026-08-24）

本手册对应论文
`helmholtz_lod_certified_amsart_revised.tex` 的 SHA256
`38cc97e4c4997b42f133903febc32eb20e5b475eca9b3d54fec441c0d3736f15`，以及
`schema_version=6`。实验按 factor、pilot、main 三道门串行执行；不要直接启动 main，
也不要并行运行两条 PALOD/HLOD 轨迹。

## 1. 本轮论文审查结论

- E1 的 Algorithm 1 仍是 safeguarded fixed-reference epoch：epoch 内 reference 固定，
  candidate 逐步富化，lazy dual 与结构刷新规则仍在，刷新后继承 `ell`。本轮主要是把
  standard/moving-reference 的边界和术语写清楚，没有改变 E1 主实验的数据生成算法。
  已完成且通过审计的 E1 主实验无需重跑；脚本同时接受该 E1 所用的上一版论文基线。
- E2 的 Algorithm 2 已实质改为 moving-reference singularity-aware LOD。每个非终止
  solve--estimate 步都冻结并富化 candidate，恢复 prospective coarse/candidate 精确匹配，
  检查 candidate Galerkin 适定性，然后立即执行
  `(T_H,T_h,T_c)<-(T_H^+,T_c^+,T_c^+)`。不存在多 H-step epoch、lazy candidate dual、
  level-gap trigger 或 graded reserve。
- 物理区域冻结为 `R_*=0.125`，`ell_S` 是覆盖该物理球的最小 coarse 邻域层数，
  `Omega_F=N^ell(Omega_S)`。coarse 只对 regular 区做完整 Dörfler 标记；candidate 在
  `Omega_F` 和 regular 区分别满足 `theta_c` Dörfler，再做共形及 hierarchy closure。
  candidate 冻结后，matching 循环只能继续加细 proposed coarse，不得再改变 candidate。
- 实验仍保留欠分辨的 pre-asymptotic 点，故 `run.json` 明确标记
  `manuscript_conformance=implementation-study-variant`；不能把这些点说成已满足论文的
  coarse resolution 假设。

## 2. 已完成的实现与性能门

E2 runner 已实现上述 moving-reference 状态机，并自动检查：

- 每个非终止 solve 恰有一次 `moving_reference_closure`；
- `requested_target_gap=0`、`matching_spill=0`，promotion 后 reference=candidate；
- moving-reference 路径 candidate dual 次数为 0；
- coarse full-regular 和 candidate F/R 两侧 Dörfler 均成立；
- `ell` 跨 moving step 继承；
- promoted candidate 解预先计算并复用于下一步 reference stability，避免重复大系统求解。

RT2 candidate flux 已将同一个 patch 约束矩阵原来的两次 rank-revealing 分解
（COD 求 particular solution、FullPivLU 求 kernel）合并为一次 FullPivLU。相同 h12
两步因子实验中：

| 指标 | 优化前 | 优化后 | 变化 |
|---|---:|---:|---:|
| wall time | 35.92 s | 30.84 s | -14.1% |
| candidate RT2 reconstruction | 12.77 s | 5.94 s | -53.5% |
| RT2 patch solve | 7.35 s | 1.76 s | -76.0% |
| peak RSS | 7,310,144 KiB | 7,297,036 KiB | 基本不变 |

所有参考误差、精确误差、candidate 指标和标记数逐项不变，完整 RT2/P2 回归通过。
曾测试整 patch 精确矩阵缓存，命中 `0/12545` 且显著变慢，已删除，不能重新启用。

本地 h12 四步 moving-reference pilot 结果：4 个 solve、3 次 promotion、0 次 dual，
精确相对能量误差
`0.44522 -> 0.36812 -> 0.31518 -> 0.26601`，reference 相对误差
`0.41271 -> 0.33460 -> 0.28368 -> 0.23456`；wall time 62.94 s，peak RSS
7,584,824 KiB，swap 0。该结果是上服务器 pilot 的门禁，不是最终论文收敛阶。

## 3. 获取代码并检查服务器

```bash
ssh shuihan
cd /home/sutai/code/LOD2d-CPP
git fetch --prune origin
git switch codex/palod-streaming-gram
git pull --ff-only origin codex/palod-streaming-gram
git status --short
git rev-parse HEAD
free -h
df -h "$PWD"
```

`git status --short` 必须没有 tracked 修改；不要删除尚未回传的 `results/`。默认使用
16 个 build/patch threads，并关闭局部线性代数的嵌套多线程，避免过度订阅。

## 4. E2 factor

```bash
tmux new -s e2-factor
cd /home/sutai/code/LOD2d-CPP
MODE=e2-revised-factor VALIDATE=1 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$PWD/results/E2-revised-factor-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

factor 依次运行 moving-reference h12 两步和 standard reference-epoch refresh factor。
必须满足：exit status 0、swap 0、无 `Failed`，moving run 为 2 solve/1 promotion/0 dual，
所有 matching spill 为 0，F/R 与 full-regular Dörfler 门通过。

## 5. E2 pilot

factor 通过后再运行：

```bash
tmux new -s e2-pilot
cd /home/sutai/code/LOD2d-CPP
MODE=e2-revised-pilot VALIDATE=1 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$PWD/results/E2-revised-pilot-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

moving pilot 为 8 个 moving steps；standard PALOD 为独立 fixed-reference 对照。E2 moving
的收敛趋势按所有 `SolveAndEstimate` 点检查，因为每一步都会更新 reference；
`relative_reference_energy` 是相对于该步 reference 的量，跨步比较时论文图必须以共同制造解
的 `relative_exact_energy` 为主。不要再要求“同一 epoch 至少三点”，也不要跨 moving step
拟合 stepwise reference error。

pilot 若超时或触发资源门，先汇总：

```bash
grep -H -E 'Elapsed|Percent of CPU|Maximum resident|Swaps|Exit status' "$RESULT_DIR"/logs/*.time
grep -H -E 'state=|stop_reason=|output=' "$RESULT_DIR"/logs/*.stdout
```

然后分析 `iterations.csv` 的 `time_corrector`、`time_theta`、
`time_candidate_flux_*`、`time_reference_riesz` 和 `time_mesh`；不要直接提高上限。

## 6. E2 main

只有 factor 与 pilot 全部通过后才启动：

```bash
tmux new -s e2-main
cd /home/sutai/code/LOD2d-CPP
MODE=e2-revised-main VALIDATE=1 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$PWD/results/E2-revised-main-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

执行顺序为 AFEM、standard PALOD、fixed HLOD、moving-reference PALOD，最慢且最需观察的
moving 方法最后运行。main 配置为 H3/h12、`R_*=0.125`、`theta_H=0.3`、
`theta_c=0.5`、15 个 moving steps。服务器约 366 GiB 时建议：

```bash
systemd-run --user --scope -p MemoryMax=300G -p MemorySwapMax=0 \
  env MODE=e2-revised-main VALIDATE=1 JOBS=16 PATCH_THREADS=16 \
  RESULT_DIR="$PWD/results/E2-revised-main-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  bash scripts/run_helmholtz_adaptive_paper_server.sh
```

## 7. 监控与验收

```bash
ROOT=/path/to/current/result-dir
watch -n 5 'free -h; pgrep -af bench_helmholtz_adaptive_paper'
tail -F "$ROOT"/logs/*.stdout
find "$ROOT/runs" -name '*.done' -print
grep -H -E 'Elapsed|Percent of CPU|Maximum resident|Swaps|Exit status' "$ROOT"/logs/*.time
```

moving E2 stdout 重点查看 `[hybrid-moving-reference]`、`[hybrid-matching]`、
`[helmholtz-model]` 和 `[hybrid-coarse-marking]`。每个 schema-v6 run 必须包含：

```text
iterations.csv
run.json
summary.csv
epoch_history.csv
mesh_manifest.csv
corrector_work.csv
hybrid_reserve.csv
```

`hybrid_reserve.csv` 在 moving 路径中保存的是 `moving_reference_closure`，不是旧版 graded
reserve。`.done` 只有在结构、Dörfler、物理半径、matching、promotion 数和误差净下降门
全部通过后才会创建。

完成后回传完整结果目录而非只回传图片：

```bash
sha256sum -c "$ROOT/SHA256SUMS"
cd "$(dirname "$ROOT")"
tar --zstd -cf "$(basename "$ROOT").tar.zst" "$(basename "$ROOT")"
sha256sum "$(basename "$ROOT").tar.zst" > "$(basename "$ROOT").tar.zst.sha256"
```
