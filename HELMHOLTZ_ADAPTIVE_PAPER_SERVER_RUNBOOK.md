# Helmholtz 自适应 LOD 论文实验服务器操作手册

本手册用于新版 practical-paper v2 实验。当前只允许先跑资源 pilot；收到 pilot
结果并冻结 E1 配置后，才运行 R2a/S 的五方法主实验。不要把名称含
`resource-pilot` 的输出写入论文。

## 1. 什么时候能得到论文图表

1. 完成 R2a、`kappa=16` 五方法后，可生成第一组论文候选图：reference energy
   error--DOF、error--累计方法时间、PALOD 的 `eta_H/Theta_loc/ell` 历史以及最终网格。
2. 完成 S、`kappa=16` 五方法后，主体图表才完整。
3. `kappa=8,32` 的代表性三方法是波数补充图；最终重复 3 次后才能冻结论文时间表。

数值轨迹若没有离开预渐近区，`run.json` 会写 `pre_asymptotic`；不为改变该标签追加
reference refinement。共同目标未达到时在 `summary.csv` 中保留 `not_reached`。

## 2. 当前本机资源基线

同一 Release 路径、16 个 OpenMP worker：

| pilot | reference | H 步 | 方法时间 | 峰值 RSS | 结果 |
|---|---:|---:|---:|---:|---|
| R2a/PALOD/k16 | level 10，1089 nodes | 3 | 13.5 s | 293 MiB | 达到 0.02；0.01 未达到 |
| S/PALOD/k16 | level 10，3201 nodes | 2 | 72.4 s | 1.37 GiB | 稳定下降，尚未达到 0.1 |

R2a 不属于大内存实验。S 的长轨迹预计明显更贵，因此先在服务器运行相同 pilot，
再决定正式 H-step/reference 上限。所有案例串行执行，不并行启动多个 benchmark。

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
MODE=smoke PATCH_THREADS=8 JOBS=16 MIN_AVAILABLE_GIB=8 \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

脚本会构建 Release、执行三个核心测试，再串行运行 R1 五方法 smoke。每个成功案例
才创建 `.done`；重复相同命令会跳过已有 `.done`。

## 6. R2a/S 资源 pilot

先比较 8 和 16 个 patch worker。使用不同结果目录，避免覆盖：

```bash
RESULT_DIR="$PWD/results/adaptive-paper-pilot-t8" \
MODE=pilot PATCH_THREADS=8 JOBS=16 MIN_AVAILABLE_GIB=16 \
  scripts/run_helmholtz_adaptive_paper_server.sh

RESULT_DIR="$PWD/results/adaptive-paper-pilot-t16" \
MODE=pilot PATCH_THREADS=16 JOBS=16 MIN_AVAILABLE_GIB=16 \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

查看时间和内存：

```bash
grep -H -E 'Elapsed|Percent of CPU|Maximum resident|Swaps' \
  results/adaptive-paper-pilot-t*/logs/*.time
grep -H 'reference_cache=' \
  results/adaptive-paper-pilot-t*/logs/*.stdout
```

选择 wall time 较短且无 swap 的设置。若 16 worker 提速不足 10%，正式运行使用 8，
因为每个 worker 会持有局部 saddle 分解和装配缓冲。

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
MODE=custom PATCH_THREADS=8 JOBS=16 MIN_AVAILABLE_GIB=64 \
RESULT_DIR="$PWD/results/adaptive-paper-e1" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

不要手工修改冻结配置。脚本只在结果目录生成 runtime copy，将 `git_commit` 和
`build_hash` 替换成服务器实际值；原模板不变。

## 10. 回传结果

```bash
cd results
tar --zstd -cf adaptive-paper-server-$(git -C .. rev-parse --short HEAD).tar.zst \
  adaptive-paper-pilot-t8 adaptive-paper-pilot-t16
sha256sum adaptive-paper-server-*.tar.zst > adaptive-paper-server.SHA256
```

回传 `.tar.zst` 和 `.SHA256`。不要只回传绘图后的 PDF；必须包含 runtime configs、
`run.json`、CSV、VTU、`.time`、`.stdout`、硬件信息和顶层 `SHA256SUMS`。
