# E1 revised reference-epoch 主实验操作手册（2026-08-24）

本手册只执行 E1。E2 在论文 Algorithm 2 再次冻结前不得启动。

## 1. 获取代码与构建

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

`git status --short` 不得包含 tracked 修改。历史 `results/` 不要删除。

## 2. 运行顺序

必须串行通过 factor、pilot、main；不能直接启动 main。

```bash
MODE=e1-revised-factor VALIDATE=1 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$PWD/results/E1-revised-factor-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  scripts/run_helmholtz_adaptive_paper_server.sh

MODE=e1-revised-pilot VALIDATE=0 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$PWD/results/E1-revised-pilot-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  scripts/run_helmholtz_adaptive_paper_server.sh

MODE=e1-revised-main VALIDATE=0 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$PWD/results/E1-revised-main-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

建议分别放入独立 `tmux`，但前一阶段生成 `.done` 且科学门通过后才能启动下一阶段。
main 可用：

```bash
systemd-run --user --scope -p MemoryMax=300G -p MemorySwapMax=0 \
  env MODE=e1-revised-main VALIDATE=0 JOBS=16 PATCH_THREADS=16 \
  RESULT_DIR="$PWD/results/E1-revised-main-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  bash scripts/run_helmholtz_adaptive_paper_server.sh
```

## 3. factor/pilot 性能门

进入 main 前必须同时满足：

- `.time` 中 `Swaps: 0`，没有 resource cap 或非结构化退出；
- `corrector_parallel_threads`、`gram_parallel_threads`、
  `candidate_flux_parallel_threads` 在正式工作行达到请求线程数；
- `corrector_work.csv` 的 active 数等于 rebuilt+reused；
- Gram structure/factor/action、corrector、RT2、embedding/mesh commit 均有分阶段计时；
- factor 覆盖首次 reference refresh 后的大网格 corrector/Gram 峰值；
- pilot 至少有一个 refresh 后含三个 `SolveAndEstimate` 点的 epoch，且 exact error 总体下降；
- `ell` 跨 refresh 继承；结构触发的 refresh 不执行无用 candidate dual；
- `mesh_manifest.csv` 覆盖 epoch 0 每个已求解 H-step 的
  coarse/reference/candidate 三元组；所有 reference 行具有同一文件和
  `reference_mesh_version`。

若 factor/pilot 变慢，先比较分阶段累计时间。只有出现新的串行热点、重复 embedding、
Gram 分解回退或内存增长才修改实现；不得通过减少 H-step、删除预收敛点或放宽误差口径
制造主结果。

## 4. 监控

```bash
ROOT=/path/to/current/result-dir
watch -n 5 'free -h; pgrep -af bench_helmholtz_adaptive_paper'
tail -F "$ROOT"/logs/*.stdout
find "$ROOT/runs" -name '*.done' -print
grep -H -E 'Elapsed|Percent of CPU|Maximum resident|Swaps|Exit status' \
  "$ROOT"/logs/*.time
```

脚本冻结 commit、二进制 SHA256 和 patch threads；同一 `RESULT_DIR` 不能混用不同构建。

## 5. 生成论文图与审计

先定位四个 run 目录，然后绘制误差--DoF/时间图：

```bash
python3 tools/visualization/plot_reference_epoch_e1.py \
  --experiment E1 \
  --palod /path/to/PALOD/run-dir \
  --fixed-lod /path/to/fixed-LOD/run-dir \
  --ufem /path/to/UFEM/run-dir \
  --afem /path/to/AFEM/run-dir \
  --output figures/paper/E1-revised-main
```

`E1-revised-main.json` 自动报告：

- PALOD 每个 epoch 的 exact/reference error--DoF 指数；
- tail-4 exact error 指数是否低于预期 `0.5` 超过容差 `0.1`；
- 在 `0.5,0.2,0.1,0.05,0.02,0.01` 等所有方法共同达到的误差目标上，各方法所需最小 DoF；
- PALOD 是否在每个可比目标上使用最少 DoF。

这些是待数据检验的实验命题。若结果为 false/null，必须原样报告，不能删点、跨 epoch
拟合或改目标。

只绘制第一个 epoch 的逐 H-step 三网格图：

```bash
python3 tools/visualization/plot_reference_epoch_meshes.py \
  --run-dir /path/to/PALOD/run-dir \
  --output-dir figures/paper/E1-revised-main-meshes \
  --epoch 0 --experiment-label E1 --all-checkpoints --page-columns 4
```

输出包含多页 PDF、逐页 PNG 和 JSON audit。每列对应同一 H-step，三行依次为
`T_H/T_h/T_c`；reference 每列都画，但同一 epoch 只存一份 VTU。JSON 必须满足
`reference_unchanged_within_epoch=true`，且只有一个 reference version/SHA256。
