# Helmholtz 自适应 LOD 论文实验服务器操作手册

本手册用于 practical-paper schema v3 实验。4/8/16 worker 资源 pilot 已完成；当前
只运行延长校准轨迹，收到结果并冻结 `practical_stop_tolerance`/epoch 判据后，才运行
R2a/S 的五方法主实验。名称含 `pilot` 或 `calibration` 的输出不得写入论文。

## 1. 什么时候能得到论文图表

1. 完成 R2a、`kappa=16` 五方法后，可生成第一组论文候选图：reference energy
   error--DOF、error--累计方法时间、PALOD 的 `eta_H/Theta_loc/ell` 历史以及最终网格。
2. 完成 S、`kappa=16` 五方法后，主体图表才完整。
3. `kappa=8,32` 的代表性三方法是波数补充图；最终重复 3 次后才能冻结论文时间表。

数值轨迹若没有离开预渐近区，`run.json` 会写 `pre_asymptotic`；不为改变该标签追加
reference refinement。共同目标未达到时在 `summary.csv` 中保留 `not_reached`。

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

## 6. R2a/S 延长校准轨迹

schema v3 的 `fixed_work_horizon` 忽略 `U_prac` 停止，仅按冻结的 H-step 上限取得
完整校准轨迹；reference error 只进入报告，不反馈 MARK/STOP。运行 R2a 6 步、S 5 步：

```bash
RESULT_DIR="$PWD/results/adaptive-paper-extended-calibration-v3" \
MODE=calibration PATCH_THREADS=4 JOBS=16 MIN_AVAILABLE_GIB=16 \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

查看时间和内存：

```bash
grep -H -E 'Elapsed|Percent of CPU|Maximum resident|Swaps' \
  results/adaptive-paper-extended-calibration-v3/logs/*.time
grep -H 'reference_cache=' \
  results/adaptive-paper-extended-calibration-v3/logs/*.stdout
```

`run.json` 应显示 `trajectory_policy=fixed_work_horizon`，并报告最后一步误差比、对数改善、
连续平台步数和 `plateau_observed`。单次误差不降会清零平台计数，不能触发平台结论。

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

在审阅两个 pilot 后再生成/提交冻结的十个 E1 配置。收到配置提交后使用：

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
tar --zstd -cf adaptive-paper-server-$(git -C .. rev-parse --short HEAD).tar.zst \
  adaptive-paper-extended-calibration-v3
sha256sum adaptive-paper-server-*.tar.zst > adaptive-paper-server.SHA256
```

回传 `.tar.zst` 和 `.SHA256`。优先作为 GitHub Release asset 上传；不要把解包后的结果
目录提交到代码分支，因为 Git 文本换行转换会使原始 `SHA256SUMS` 失效。不要只回传 PDF；必须包含 runtime configs、
`run.json`、CSV、VTU、`.time`、`.stdout`、硬件信息和顶层 `SHA256SUMS`。
