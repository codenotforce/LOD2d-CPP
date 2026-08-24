# 新版 Algorithm 1/2：E1/E2 服务器操作手册（2026-08-24）

> 2026-08-24 后续决定：当前只执行 E1，采用
> `E1_REVISED_MAIN_SERVER_2026-08-24.md`；E2 等论文 Algorithm 2 再次冻结后再恢复。

本手册对应 `schema_version=6`、新版制造解和新版 reference-epoch 状态机。所有结果均为
`claim=implementation-study`。必须依次通过 factor、pilot、main 三道门；不得直接启动 main，
也不得并行运行两条 PALOD/HLOD 轨迹。

## 1. 算法口径与已知修正

- E1 使用新版 Algorithm 1：在当前 candidate 上先计算指标、标记和加细，再添加 hierarchy
  closure，整个 candidate 阶段先于分支；结构刷新跳过 dual；每个
  epoch 首解初始化 lazy-dual baseline；收敛前强制 dual；刷新后继承 `ell`；reserve
  trigger/target 为 `4/9`。
- 按既定 implementation-study 口径，实验保留 `kappa H` 很大的预收敛点，不执行论文
  coarse-admissibility 分辨率门；reference/candidate 适定性采用数值稀疏分解检查。因此
  Algorithm 1 的 `manuscript_conformance` 为 `implementation-study-variant`，不是
  `direct`。
- E2 hybrid 使用新版 Algorithm 2：
  `ell_S=min{j:B_R(S) subset N^j(S)}`、`Omega_S=N^ell_S(S)`、
  `Omega_F=N^ell(Omega_S)`；只跳过 `Omega_S` 中的零 corrector；coarse 只在 regular
  区标记；candidate 在 `Omega_F` 与 regular 区分别做 Dörfler 标记，并同样在标记/加细后
  才添加 hierarchy closure。新版定义不再要求 `ell_S>=ell`：固定物理半径只决定
  singular core，oversampling `ell` 独立决定 transition buffer。
- 一般共形 NVB 网格上，“`Omega_F` 精确 H/h matching”与“全部相邻 regular 单元均有
  正的统一 reserve”不总能同时实现。程序因此显式采用可审计的 implementation erratum：
  自动选择最小的 spill-free conformity collar，目标剖面为
  `min(g_tar,max(0,d_graph(T,Omega_F)-b))`。每次 trial/closure 写入
  `hybrid_reserve.csv`，`run.json.algorithm_variant.manuscript_conformance` 明确为
  `implementation-erratum`，不得把它隐去后声称逐字实现论文当前公式。
- E2 coarse marking 优先从 closure-safe regular 子集选点，但所选集合仍必须满足原始完整
  regular 指标质量的 Dörfler 不等式；若安全子集质量不足，回退到原始 full-regular
  标记并按论文触发 structural refresh。h8 校准表明 `theta_H=0.1` 可形成 3 点完整
  epoch；`0.3/0.5` 只有 1--2 点，故主配置冻结为 `0.1`。

## 2. 获取代码并检查服务器

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

`git status --short` 必须没有 tracked 修改；不要删除尚未回传的 `results/`。脚本默认使用
16 个 build/patch threads，并按阶段设置可用内存、磁盘和外部 timeout 门。E1 fixed-LOD
配置因单任务局部直接分解的内存压力保留 `maximum_patch_threads=4`，其余新版 PALOD
配置使用 16 个 patch threads。

## 3. Factor：必须首先运行

E1 factor 覆盖首次 refresh 后的 corrector/Gram 峰值；旧的单 H-step factor 不足以判断
主实验内存。

```bash
tmux new -s e1-factor
MODE=e1-revised-factor VALIDATE=1 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$PWD/results/E1-revised-factor-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

E1 通过后运行 E2 hybrid/standard factor：

```bash
tmux new -s e2-factor
MODE=e2-revised-factor VALIDATE=1 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$PWD/results/E2-revised-factor-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

建议在支持 user systemd 的服务器上直接用下列形式设置运行期内存上限（把 `MODE` 和
`RESULT_DIR` 换成相应阶段）：

```bash
systemd-run --user --scope -p MemoryMax=64G -p MemorySwapMax=0 \
  env MODE=e1-revised-factor VALIDATE=1 JOBS=16 PATCH_THREADS=16 \
  RESULT_DIR="$PWD/results/E1-revised-factor-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  bash scripts/run_helmholtz_adaptive_paper_server.sh
```

如果 user systemd 不可用，保留脚本的 DoF、patch、available-memory 和 timeout 门，不要
用 swap 兜底。

Factor gate：

- exit status 0，`run.json.status != Failed`；
- stop reason 只能是固定 H-step 预算或“剩余预算不足以开启新 epoch”；
- `Swaps: 0`；
- E2 `hybrid_reserve.csv` 中每个 `closure` 行均为 `status=achieved`、
  `matching_spill=0`、`profile_margin_after>=0`；
- E2 每个 corrector 均满足 `hybrid_covered_physical_radius>=0.125` 且 patch 上限不超过
  100000 个 fine elements。

## 4. Pilot：factor 通过后运行

```bash
tmux new -s e1-pilot
MODE=e1-revised-pilot VALIDATE=0 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$PWD/results/E1-revised-pilot-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  scripts/run_helmholtz_adaptive_paper_server.sh

tmux new -s e2-pilot
MODE=e2-revised-pilot VALIDATE=0 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$PWD/results/E2-revised-pilot-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

两条 pilot 仍应串行；第一条结束后再启动第二条。推荐额外限制
`MemoryMax=160G,MemorySwapMax=0`。

Pilot gate：

- E1 必须至少形成一个 refresh 后含 3 个解点的可拟合 epoch，`ell` 不得在 refresh 后
  回退；保留预收敛点，不跨 epoch 拟合。
- E2 hybrid 至少有一个完整 epoch 含 3 个 `SolveAndEstimate` 点；这些点的 exact error
  应总体下降。h8 本地门禁的首个三点段斜率约为 `-0.47`，只作为“可以上服务器”的
  证据，不是最终论文斜率。
- `iterations.csv` 中 candidate F/R 两侧各自达到 `theta_c=0.5`；
  `hybrid_full_regular_doerfler=true`。
- 如果 factor/pilot 触发 resource cap，先分析 `corrector_work.csv`、Gram/RT2 分阶段计时和
  `.time`；不得直接提高上限。

## 5. Main：只在两组 pilot 通过后启动

E1 四方法：

```bash
tmux new -s e1-main
MODE=e1-revised-main VALIDATE=0 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$PWD/results/E1-revised-main-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

E2 四方法（fixed LOD 成本过高时可在前三条完成后单独停止/补跑）：

```bash
tmux new -s e2-main
MODE=e2-revised-main VALIDATE=0 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$PWD/results/E2-revised-main-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

E2 main 中 AFEM/fixed-LOD 配置各保留最多 24 小时内部 wall guard，脚本外部门为
24 小时加 15 分钟，以便收到 `SIGINT` 后结构化写完 CSV/JSON；PALOD 通常应远早于此结束。

推荐 main 使用 `MemoryMax=300G,MemorySwapMax=0`；脚本另要求启动时至少 320 GiB
available memory 和 200 GiB free disk。可执行形式为：

```bash
systemd-run --user --scope -p MemoryMax=300G -p MemorySwapMax=0 \
  env MODE=e1-revised-main VALIDATE=0 JOBS=16 PATCH_THREADS=16 \
  RESULT_DIR="$PWD/results/E1-revised-main-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  bash scripts/run_helmholtz_adaptive_paper_server.sh
```

E1 PALOD 为 H4/h12、12 解点、epoch 安全上限 13；E2 hybrid 为 H3/h12、
`R_*=0.125`、`theta_H=0.1`、15 解点、epoch 安全上限 16；E2 standard 为 H3/h12、
12 解点、target gap 9、epoch 安全上限 13。epoch 上限只防失控，正常终止由 H-step、
资源或 wall-time 门决定。

## 6. 运行中监控

```bash
ROOT=/path/to/current/result-dir
watch -n 5 'free -h; pgrep -af bench_helmholtz_adaptive_paper'
tail -F "$ROOT"/logs/*.stdout
find "$ROOT/runs" -name '*.done' -print
grep -H -E 'Elapsed|Percent of CPU|Maximum resident|Swaps|Exit status' \
  "$ROOT"/logs/*.time
```

stdout 使用 line buffering，可实时看到：

- `[helmholtz-model]`：operators/correctors/basis/factorization；
- `[hybrid-coarse-marking]`：regular/safe/marked mass、collar、是否 closure-safe；
- `[hybrid-refresh-collar-trial]`：spill/deficit/far gap；
- `[hybrid-refresh-closure]`：最终 collar 和 candidate 尺寸；
- `[reference-stability]`：每 epoch 一次的参考求解。

正式 CSV/JSON 在 driver 结构化结束后写出。进程被外部杀死时保留 stdout、`.time` 和整个
未完成目录，不要创建 `.done`，也不要覆盖重跑到同一 `RESULT_DIR`。
脚本会在 `server-build-identity.txt` 冻结 commit、二进制 SHA256 和 patch thread 数；若
同一目录已记录不同 identity，会在读取旧 `.done` 前拒绝运行，防止跨提交混合结果。

新版模式只有在通用结构检查及上述 factor/pilot 科学门（零 swap、reserve/matching、
Dörfler、物理半径、patch cap、`ell` 继承与三点 epoch）自动通过后才创建 `.done`；因此
`.done` 表示该配置通过当前阶段门，而不只是进程退出。

## 7. 完成验收与回传

```bash
sha256sum -c "$ROOT/SHA256SUMS"
find "$ROOT/runs" -name run.json -print
grep -H -E 'state=|stop_reason=|output=' "$ROOT"/logs/*.stdout
```

每个 schema-v6 run 必须包含：

```text
iterations.csv
run.json
summary.csv
epoch_history.csv
mesh_manifest.csv
corrector_work.csv
hybrid_reserve.csv
```

E2 hybrid 还应包含 `mesh_E2_final_hybrid_regions.vtu`。最终收敛率只按 epoch 内的
`SolveAndEstimate` 点拟合；少于 3 点的 epoch 标为 `insufficient_data`，不得跨 refresh
拼接。`time_method_cumulative` 已扣除 exact/reference validation 和网格快照复制，但
wall-time 保护仍使用包含这些开销的原始 wall clock，二者不得混称。

回传完整目录（runtime config、CSV/JSON、VTU、stdout、time、硬件信息、SHA256），不要只
回传图片：

```bash
cd "$(dirname "$ROOT")"
tar --zstd -cf "$(basename "$ROOT").tar.zst" "$(basename "$ROOT")"
sha256sum "$(basename "$ROOT").tar.zst" > "$(basename "$ROOT").tar.zst.sha256"
```
