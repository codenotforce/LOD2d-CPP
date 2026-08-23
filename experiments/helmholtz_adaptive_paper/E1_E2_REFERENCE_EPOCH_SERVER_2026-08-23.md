# E1/E2 reference-epoch 主实验服务器操作手册（2026-08-23）

本手册对应当前 schema-v5 reference-epoch 实现和保留的 schema-v4 对照方法。E1 是
localized smooth；E2 是 L-shaped mixed-boundary low-regularity。两组实验不得跨 epoch
拟合 reference error，也不得把 `WorkLimitReached: maximum_H_steps reached` 误判为失败。

## 1. 获取冻结代码

```bash
cd ~/code/LOD2d-CPP
git fetch --prune origin
git switch codex/palod-streaming-gram
git pull --ff-only origin codex/palod-streaming-gram
git status --short
git rev-parse HEAD
```

`git status --short` 必须为空。不要删除服务器上尚未回传的 `results/`；必要时先移动到
仓库外归档。项目同时保留两份论文哈希：schema-v5 使用
`MANUSCRIPT_BASELINE.sha256`，schema-v4 对照使用
`MANUSCRIPT_BASELINE_LEGACY_V4.sha256`，服务器脚本会按配置自动选择。

## 2. 资源与线程

366 GiB 服务器建议：

```bash
export JOBS=32
export PATCH_THREADS=16
export MIN_AVAILABLE_GIB=128
```

所有方法由脚本串行执行；不要同时启动两条 PALOD/HLOD 轨迹。16 个 patch threads
避免小 patch 过度订阅，并为 Eigen/SuiteSparse 和系统保留核。误差--时间比较必须保持
相同 `PATCH_THREADS`、编译器和提交。

## 3. 一次运行 E1 与 E2

### 3.1 先校准 E2 hybrid 的固定物理半径

`R_*=0.25` 在 H3/h12 上产生了单个病态大的 transition corrector patch，服务器运行约
7 小时仍未完成。不要直接重启该配置。先顺序运行 `R_*=0.0625,0.125` 的单 H-step
pilot：

```bash
tmux new -s helmholtz-e2-radius-pilot

MODE=e2-radius-pilot \
JOBS=32 PATCH_THREADS=16 MIN_AVAILABLE_GIB=128 VALIDATE=1 \
RESULT_DIR="$PWD/results/E2-radius-pilot-$(git rev-parse --short HEAD)" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

提交 `2fd2ec3` 的两组 pilot 均已完成。`R_*=0.0625/0.125` 的 wall 为
`13.98/14.84 s`，峰值 RSS 约 `2.095/2.099 GiB`，最大 patch 均为 3360，exact
relative energy 为 `0.100152/0.100157`，skipped work units 为 `2200/6408`。因此正式
E2 hybrid 冻结 `R_*=0.125`。

运行时查看新增的 `[hybrid-matching]`、`[hybrid-preflight]` 和模型阶段日志。每个 hybrid
配置在 corrector 前以 `100000` 个 reference elements 为单 patch 硬上限；超过即明确
失败。若网格、参考层或制造解改变，必须重新运行本节 pilot，不得直接沿用该冻结值。

### 3.2 正式运行

```bash
chmod +x scripts/run_helmholtz_adaptive_paper_server.sh
tmux new -s helmholtz-e1-e2

MODE=e1-e2-main \
JOBS=32 PATCH_THREADS=16 MIN_AVAILABLE_GIB=128 VALIDATE=1 \
RESULT_DIR="$PWD/results/E1-E2-reference-epoch-$(git rev-parse --short HEAD)" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

断开 tmux：`Ctrl-b d`；恢复：

```bash
tmux attach -t helmholtz-e1-e2
```

脚本只为完整通过验收的配置写 `.done`；再次执行同一命令会跳过已完成配置。若只想先拿
E2 headline，可运行：

```bash
MODE=e2-main \
JOBS=32 PATCH_THREADS=16 MIN_AVAILABLE_GIB=128 VALIDATE=1 \
RESULT_DIR="$PWD/results/E2-reference-epoch-$(git rev-parse --short HEAD)" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

E2 顺序为 hybrid PALOD、standard PALOD、AFEM、fixed LOD；即使最后一个 fixed LOD
成本过高，前三条正式主曲线也已保留。只补 fixed LOD 时使用：

```bash
MODE=custom VALIDATE=0 JOBS=32 PATCH_THREADS=16 MIN_AVAILABLE_GIB=128 \
CONFIGS='experiments/helmholtz_adaptive_paper/configs/E2-S-hlod-fixed-k16-H3-h16-ell3-step12-v4.json' \
RESULT_DIR="$PWD/results/E2-fixed-LOD-$(git rev-parse --short HEAD)" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

## 4. 运行中监控

```bash
watch -n 5 'free -h; ps -C bench_helmholtz_adaptive_paper \
  -o pid,etime,%cpu,%mem,rss,vsz,cmd --sort=-rss'
```

查看当前方法和阶段：

```bash
ROOT="results/E1-E2-reference-epoch-$(git rev-parse --short HEAD)"
find "$ROOT/runs" -name '*.done' -print
tail -n 40 "$ROOT"/logs/*.stdout
grep -H -E 'Elapsed|Percent of CPU|Maximum resident|Swaps|Exit status' \
  "$ROOT"/logs/*.time
```

schema-v5 runner 在轨迹结束时原子式写 `iterations.csv`，因此运行中以 stdout 的
`corrector/Gram/RT2/mesh` 分阶段日志和进程 RSS 为准；不要因为 CSV 暂未出现就重启。
swap 持续增长、`MemAvailable<36 GiB` 或 RSS 无界增长时停止当前 benchmark，保留日志，
不要删除已完成 `.done`。

## 5. E1/E2 配置口径

E1 PALOD：`H2/h12`、`ell0=2`、gap guard 4、refresh target gap 6、15 H-steps、最多
4 epochs，ell 跨 epoch 继承，且只在剩余预算至少可形成 4 个新解点时刷新。

E2 制造解固定为

\[
u=\chi(r)r^{2/3}\sin(2\theta/3)+0.05\,\psi(x,y)e^{i\kappa x},
\qquad \kappa=16,
\]

凹角两条边为 homogeneous Dirichlet，外边为 homogeneous Robin。E2 hybrid 与 standard
均从 H-level 3、h-level 12 开始；hybrid `ell<=4`，standard 允许证书自动增长到
`ell<=10`。E2 hybrid 冻结 `hybrid_minimum_physical_radius=0.125`，并在每次
corrector check 选择最小的 (l_s\ge\ell)，使
(B_{R_*}(z)\cap\Omega\subset N^{2l_s}(z))；因此局部加细不会让奇异/过渡区域的
物理范围缩小。hybrid 的 level-gap guard 只检查 regular region 中的正 gap；匹配区
`H=h` 的零 gap 由严格 containment 判据保护。禁止改回全局最小 gap，否则每个 hybrid
epoch 会被构造性的零 gap 立即刷新。

历史 H3/h8、`R_*=0.25` 固定半径探针在 refresh 后得到 `hybrid_l_s=11`、
`hybrid_covered_physical_radius=0.254116...`，wall 约 67 s、峰值 RSS 约 3.16 GiB、
无 swap。h12 已证明 `0.25` 会形成不可接受的大 patch；当前 h12 pilot 已据 preflight
patch count、wall 与 RSS 将生产值校准为 `0.125`。

## 6. 完成验收

```bash
sha256sum -c "$ROOT/SHA256SUMS"
find "$ROOT/runs" -name '*.done' -print
grep -H -E 'state=|stop_reason=|output=' "$ROOT"/logs/*.stdout
grep -H -E 'Elapsed|Percent of CPU|Maximum resident|Swaps|Exit status' \
  "$ROOT"/logs/*.time
```

schema-v5 合法终态为 `Converged` 或带明确资源原因的 `WorkLimitReached`；`Failed` 不可用。
每个 v5 run 必须包含 `iterations.csv/run.json/summary.csv/epoch_history.csv/`
`mesh_manifest.csv/corrector_work.csv`。E2 hybrid 还必须包含
`mesh_E2_final_hybrid_regions.vtu`。

E2 headline gate：

- exact energy error 在同一 epoch 内至少 3--4 个点稳定下降；
- 收敛阶只按 epoch 分段拟合；目标是自适应二维能量误差约 `N^{-1/2}`，不能跨 refresh
  拼接制造该斜率；
- hybrid 的 `skipped_corrector_work_units` 必须为正；
- hybrid 每个 corrector record 必须满足 `hybrid_l_s >= ell`、
  `hybrid_covered_physical_radius >= hybrid_minimum_physical_radius`；
- standard 的大 `ell`、预收敛段和较高成本原样保留；
- mixed boundary 与上述制造解四个参数必须出现在 v5 `run.json.config` 中。

## 7. 绘图

先定位每个配置下唯一的 run 目录：

```bash
find "$ROOT/runs" -name run.json -printf '%h\n'
```

E2 误差、时间、分 epoch 收敛阶和 skipped work：

```bash
python3 tools/visualization/plot_reference_epoch_e1.py \
  --experiment E2 \
  --hybrid <E2-hybrid-run-dir> \
  --standard <E2-standard-run-dir> \
  --afem <E2-afem-run-dir> \
  --fixed-lod <E2-fixed-lod-run-dir> \
  --output figures/paper/E2-S-k16-reference-epoch-main
```

E2 真实网格演化、最终 matching region 和凹角放大：

```bash
python3 tools/visualization/plot_reference_epoch_meshes.py \
  --experiment-label E2 --epoch 0 --checkpoints 4 \
  --run-dir <E2-hybrid-run-dir> \
  --output-dir figures/paper/E2-S-k16-reference-epoch-meshes
```

E1 使用同一脚本，误差图传入 `--experiment E1 --palod ... --fixed-lod ... --ufem ...
--afem ...`；网格图传入 `--experiment-label E1`。

## 8. 回传

```bash
cd results
tar --zstd -cf E1-E2-reference-epoch-$(git -C .. rev-parse --short HEAD).tar.zst \
  E1-E2-reference-epoch-$(git -C .. rev-parse --short HEAD)
sha256sum E1-E2-reference-epoch-*.tar.zst > E1-E2-reference-epoch.SHA256
```

回传 tarball 与 SHA256，或把结果放到单独数据分支；不要只回传 PNG/PDF。必须保留
runtime configs、CSV、VTU、time/stdout、硬件信息和顶层 `SHA256SUMS`。
