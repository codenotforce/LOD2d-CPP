# E1 统一 H4、theta_H=0.3 主实验操作手册（2026-08-24）

本轮只执行 E1。四种方法统一从 `initial_coarse_level=4` 开始；自适应方法
使用 `theta_H=0.3`。`theta_c=0.5` 保持不变，因为它控制 candidate 证书质量，
不应与 coarse Dörfler 参数一起调整。

## 1. 获取代码

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

`git status --short` 不得包含 tracked 修改。历史 `results/` 无需删除。

## 2. 门禁顺序

先运行 factor，再运行 pilot。二者不能同时启动：

```bash
MODE=e1-revised-unified-factor VALIDATE=1 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$PWD/results/E1-unified-H4-theta03-factor-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  scripts/run_helmholtz_adaptive_paper_server.sh

MODE=e1-revised-unified-pilot VALIDATE=0 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$PWD/results/E1-unified-H4-theta03-pilot-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

进入 main 前检查：

- `.time` 中 `Swaps: 0` 且退出码为 0；
- `corrector_work.csv` 中 `active=rebuilt+reused`；
- `corrector_cache_oversized_misses=0`，否则 32768-DoF 缓存准入仍不足；
- `corrector_cache_peak_bytes <= 2147483648`；该字段把 LRU 保留量与
  并行稀疏分解的临时工作集区分开；
- `gram_structure_parallel_threads=16`、`gram_parallel_threads=16`；
- refresh 后 `ell` 继承，且每个新 epoch 至少有 3 个已求解点；
- exact energy error 总体下降，不能通过删掉 pre-asymptotic 点制造收敛。

## 3. tmux 后台主实验

主实验脚本内部按预计耗时从短到长串行执行：

1. AFEM（H4，`theta_H=0.3`，最多 60 个 H-step）；
2. UFEM（H4 到 uniform level 18，共 14 步）；
3. fixed LOD（H4、h18、ell=3、最多 18 个 H-step）；
4. reference-epoch PALOD（H4/h12、最多 20 个 H-step、最多 4 个 epoch）。

```bash
cd /home/sutai/code/LOD2d-CPP
SESSION=e1-h4-theta03
RESULT_DIR="$PWD/results/E1-unified-H4-theta03-main-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)"

tmux new-session -d -s "$SESSION" \
  "cd '$PWD' && systemd-run --user --scope -p MemoryMax=340G -p MemorySwapMax=0 env MODE=e1-revised-unified-main VALIDATE=0 JOBS=16 PATCH_THREADS=16 RESULT_DIR='$RESULT_DIR' bash scripts/run_helmholtz_adaptive_paper_server.sh 2>&1 | tee '$RESULT_DIR.launch.log'"

tmux attach -t "$SESSION"
```

从 tmux 脱离使用 `Ctrl-b d`，不会终止实验。

## 4. 监控

```bash
ROOT=/path/to/current/result-dir
watch -n 5 'free -h; pgrep -af bench_helmholtz_adaptive_paper'
tail -F "$ROOT"/logs/*.stdout
find "$ROOT/runs" -name '*.done' -print
grep -H -E 'Elapsed|Percent of CPU|Maximum resident|Swaps|Exit status' \
  "$ROOT"/logs/*.time
```

PALOD 当前 epoch/H-step：

```bash
RUN=$(find "$ROOT/runs" -path '*PALOD*' -name iterations.csv -print -quit)
tail -n 5 "$RUN"
```

## 5. 绘图与掉阶判据

```bash
python3 tools/visualization/plot_reference_epoch_e1.py \
  --experiment E1 \
  --palod /path/to/PALOD/run-dir \
  --fixed-lod /path/to/fixed-LOD/run-dir \
  --ufem /path/to/UFEM/run-dir \
  --afem /path/to/AFEM/run-dir \
  --output figures/paper/E1-unified-H4-theta03-main
```

图中只含误差--DoF，不绘制 cumulative wall time。JSON audit 按 epoch 分段拟合
PALOD；禁止跨 refresh 拟合。主判据为最后一个具有至少 3 个点的 epoch 的 exact
energy--DoF 指数不低于 `0.4`（理论目标 `0.5`，容差 `0.1`）。

若掉阶，依次检查并处理：

1. 最后 epoch 点数不足：增加 H-step 预算，不改误差口径；
2. reference error 已平台：提高 refresh 后目标 gap 或加深新 reference；
3. localization error 占主导：提高继承的 `ell`/`ell_max`，不得跨 epoch 重置 ell；
4. 仍掉阶：保持 H4 与 `theta_H=0.3`，增加每 epoch 的局部加细点数并重新分段拟合。

第一个 epoch 的三网格图：

```bash
python3 tools/visualization/plot_reference_epoch_meshes.py \
  --run-dir /path/to/PALOD/run-dir \
  --output-dir figures/paper/E1-unified-H4-theta03-main-meshes \
  --epoch 0 --experiment-label E1 --all-checkpoints --page-columns 4
```

同一 epoch 的 reference 每个 H-step 都绘制；JSON 必须报告
`reference_unchanged_within_epoch=true`。
