# Helmholtz 自适应 LOD 论文实验服务器操作手册

本手册用于 practical-paper schema v4 实验。资源 pilot、epoch-2 校准和本机可承受的
R2a 深层 reference gate 已完成；当前只在服务器串行运行 R2a epoch 6/level 16 与
S epoch 3/level 13 校准，再做相邻 reference audit。两项 audit 都通过后才能冻结
R2a/S 的五方法主实验。名称含 `pilot`、`calibration` 或 `reference-adequacy` 的输出
不得写入论文。

## 1. 什么时候能得到论文图表

1. 完成 R2a、`kappa=16` 五方法后，可生成第一组论文候选图：reference energy
   error--DOF、error--累计方法时间、PALOD 的 `eta_H/Theta_loc/ell` 历史以及最终网格。
2. 完成 S、`kappa=16` 五方法后，主体图表才完整。
3. `kappa=8,32` 的代表性三方法是波数补充图；最终重复 3 次后才能冻结论文时间表。

数值轨迹若没有离开预渐近区，`run.json` 会写 `pre_asymptotic`；不只为改变该标签追加
reference refinement。只有独立 adequacy gate 失败才推进显式 reference epoch。共同目标
未达到时在 `summary.csv` 中保留 `not_reached`。

## 2. 已完成的服务器资源基线

AMD EPYC 9554 服务器、同一 Release 二进制：

| worker | R2a 时间/RSS | S 时间/RSS | swap |
|---:|---:|---:|---:|
| 4 | 9.32 s / 247 MiB | 35.94 s / 1.29 GiB | 0 |
| 8 | 9.17 s / 261 MiB | 35.42 s / 1.31 GiB | 0 |
| 16 | 9.24 s / 289 MiB | 35.74 s / 1.35 GiB | 0 |

三组非计时数值轨迹逐字段一致。4 worker 相对 8 worker 只慢约 1.5%，因此冻结
`PATCH_THREADS=4`；16 worker 不再使用。certificate 占方法时间约 85%--89%。所有案例
串行执行，不并行启动多个 benchmark。

脚本现在会拒绝 calibration/custom 模式的非 4 线程值。只有明确标成非生产线程研究时
才可设置 `ALLOW_NONFROZEN_PATCH_THREADS=1`；该输出不得用于冻结资源口径。

## 3. 服务器要求

- Ubuntu 22.04 或更新版本；本地 SSD 工作目录，避免 NFS 构建。
- 建议至少 32 GiB RAM、16 个物理核、10 GiB 可用磁盘；正式 S 扩展建议保留
  64 GiB `MemAvailable`。
- 依赖：

```bash
sudo apt update
sudo apt install -y build-essential cmake g++ python3 \
  libeigen3-dev libsuitesparse-dev libtbb-dev time tmux
```

服务器代理只影响拉取代码；实验本身不联网。

## 4. 获取并核对代码

```bash
git clone https://github.com/codenotforce/LOD2d-CPP.git
cd LOD2d-CPP
git switch codex/reference-epoch-hierarchy   # 合并到 main 后改成 main
git pull --ff-only
git status --short
git rev-parse HEAD
```

工作区必须干净。不要用 `git reset --hard` 清理含有未回传结果的服务器目录。

## 5. correctness smoke

```bash
chmod +x scripts/run_helmholtz_adaptive_paper_server.sh
MODE=smoke PATCH_THREADS=4 JOBS=16 MIN_AVAILABLE_GIB=8 \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

脚本会构建 Release、执行三个核心测试，再串行运行 R1 五方法 smoke。每个成功案例
才创建 `.done`；重复相同命令会跳过已有 `.done`。

## 6. 当前 reference gate 与下一批服务器校准

epoch-2 的独立 12→13 audits 已完成且均失败：R2a 的
`terminal_error_fraction=8.51361`，S 为 `1.15148`。R2a 在本机继续完成 level 13、14、15
的六步校准，但对应相邻 audit 分数仍为 `4.61608/2.93656/2.44423`。level-15 运行 peak
RSS 为 8.40 GiB、wall time 为 24:49、swap 为 0；level 16 转服务器。S level-12 上传结果
为 `success/TrajectoryComplete`，但实际用了 8 patch threads；它只用于决定推进 epoch，
下一层必须恢复冻结的 4 线程。

先运行两个 calibration。脚本串行执行配置；建议至少保留 64 GiB `MemAvailable`：

```bash
CONFIGS='experiments/helmholtz_adaptive_paper/configs/R2a-palod-k16-epoch6-level16-calibration-v4.json experiments/helmholtz_adaptive_paper/configs/S-palod-k16-epoch3-level13-step6-calibration-v4.json' \
MODE=custom PATCH_THREADS=4 JOBS=16 MIN_AVAILABLE_GIB=64 \
RESULT_DIR="$PWD/results/adaptive-paper-reference-epoch-next-server" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

不得只看退出码。两个配置都必须有 `.done`，且各自 `run.json` 必须为
`status=success`、`driver_state=TrajectoryComplete`、
`stop_reason="fixed H-step trajectory complete"`；`iterations.csv` 最终动作必须为
`CompleteTrajectory`。同时检查：

```bash
sha256sum -c results/adaptive-paper-reference-epoch-next-server/SHA256SUMS
grep -H -E 'Elapsed|Percent of CPU|Maximum resident|Swaps|Exit status' \
  results/adaptive-paper-reference-epoch-next-server/logs/*.time
grep -H -E 'state=|convergence_regime=|reference_cache=' \
  results/adaptive-paper-reference-epoch-next-server/logs/*.stdout
cat results/adaptive-paper-reference-epoch-next-server/server-build-identity.txt
```

两个 calibration 均成功后，分别运行相邻 reference audit。下面的 `find` 每个目录必须
恰好找到一份 `iterations.csv`；若为空或多于一份，先停下检查，不能猜 run ID：

```bash
CALIB="$PWD/results/adaptive-paper-reference-epoch-next-server"
R2A_SOURCE=$(find "$CALIB/runs/R2a-palod-k16-epoch6-level16-calibration-v4" \
  -mindepth 2 -maxdepth 2 -type f -name iterations.csv -print)
S_SOURCE=$(find "$CALIB/runs/S-palod-k16-epoch3-level13-step6-calibration-v4" \
  -mindepth 2 -maxdepth 2 -type f -name iterations.csv -print)
[[ $(printf '%s\n' "$R2A_SOURCE" | grep -c .) -eq 1 ]]
[[ $(printf '%s\n' "$S_SOURCE" | grep -c .) -eq 1 ]]

BUILD_DIR="$PWD/build-adaptive-paper-server" \
RESULT_DIR="$PWD/results/R2a-reference-adequacy-epoch6-level16-to17-v4" \
SOURCE_ITERATIONS="$R2A_SOURCE" \
TEMPLATE="$PWD/experiments/helmholtz_adaptive_paper/configs/R2a-palod-k16-epoch6-level16-reference-audit-v4.json" \
JOBS=16 scripts/run_helmholtz_reference_adequacy_server.sh

BUILD_DIR="$PWD/build-adaptive-paper-server" \
RESULT_DIR="$PWD/results/S-reference-adequacy-epoch3-level13-to14-v4" \
SOURCE_ITERATIONS="$S_SOURCE" \
TEMPLATE="$PWD/experiments/helmholtz_adaptive_paper/configs/S-palod-k16-epoch3-level13-reference-audit-v4.json" \
JOBS=16 scripts/run_helmholtz_reference_adequacy_server.sh
```

逐目录运行 `sha256sum -c SHA256SUMS` 并查看 `reference_adequacy.json`。只有两项
`terminal_error_fraction<=0.25` 才冻结 reference epoch；否则继续下一显式 epoch。
audit 只求相邻两个 reference FEM，绝不重跑 PALOD。

`run.json` 的平台判据仍只用于诊断：三点窗几何平均比不小于 0.9 且整体波动不超过 15%
才报告 `plateau_observed=true`；reference error 不反馈 MARK/STOP。真正资源上限才使用
`censored_*` 状态。

## 7. 运行时监控和停止条件

在 `tmux` 中运行脚本，另开窗口监控：

```bash
watch -n 5 'free -h; ps -C bench_helmholtz_adaptive_paper \
  -o pid,etime,%cpu,%mem,rss,vsz,cmd --sort=-rss'
```

立即停止当前 benchmark 的条件：

- swap 持续增长；
- `MemAvailable` 低于 10% 总内存；
- RSS 连续增长但迭代日志/CPU 长时间无进展；
- 文件系统剩余空间低于 5 GiB。

只终止 benchmark 进程；不要删除已完成案例的 `.done`、run 目录或共享
`reference-cache`。重新执行相同命令会继续未完成案例。

## 8. reference cache 与计时口径

脚本把五方法指向同一个 `reference-cache`。缓存键绑定 reference 网格、边界、完整
算子、载荷、求解器格式、Git commit 和 benchmark 二进制 SHA-256。损坏或输入变化会
自动 miss 并重算。`run.json` 记录 cache hit/key。

reference 解和流式 reference-error 评估始终从 method time 排除。候选向量每步评估后
立即释放，不再随迭代数累积保存在 journal 中。

## 9. 正式 E1 运行

当前 E1 被 reference adequacy gate 明确阻塞。只有 R2a 和 S 的最终 audit 都通过后，才生成/
提交冻结的十个 E1 配置。收到配置提交后使用：

```bash
CONFIGS='experiments/helmholtz_adaptive_paper/configs/<frozen-1>.json ...' \
MODE=custom PATCH_THREADS=4 JOBS=16 MIN_AVAILABLE_GIB=64 \
RESULT_DIR="$PWD/results/adaptive-paper-e1" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

不要手工修改冻结配置。脚本只在结果目录生成 runtime copy，将 `git_commit` 和
`build_hash` 替换成服务器实际值；原模板不变。

## 10. 回传结果

```bash
cd results
tar --zstd -cf adaptive-paper-reference-epoch-$(git -C .. rev-parse --short HEAD).tar.zst \
  adaptive-paper-reference-epoch-next-server \
  R2a-reference-adequacy-epoch6-level16-to17-v4 \
  S-reference-adequacy-epoch3-level13-to14-v4
sha256sum adaptive-paper-reference-epoch-*.tar.zst \
  > adaptive-paper-reference-epoch.SHA256
```

回传 `.tar.zst` 和 `.SHA256`。优先作为 GitHub Release asset 上传；不要把解包后的结果
目录提交到代码分支，因为 Git 文本换行转换会使原始 `SHA256SUMS` 失效。新脚本在结果目录
内写相对路径清单，解包后可直接运行 `sha256sum -c SHA256SUMS`。不要只回传 PDF；必须包含 runtime configs、
`run.json`、CSV、VTU、`.time`、`.stdout`、硬件信息和顶层 `SHA256SUMS`。
