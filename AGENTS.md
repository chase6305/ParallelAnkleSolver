# 仓库约定

本仓库用于双足机器人脚踝与腿部并联机构求解。核心为 C++17/Eigen，YAML 加载器独立可选。保留几何参数来源和单位；旋转脚踝关节顺序为 `[pitch, roll]`，电机顺序为接线顺序 `[motor1, motor2]`。其他机构必须分别定义任务坐标、转动顺序、执行器坐标及雅可比单位，各模型分别使用自己的坐标约定。新增机型用独立配置，不覆盖旧模型的参数或回归数据。

修改求解算法时优先验证杆长闭环、IK/FK 往返、速度/力矩功率守恒及失败状态。验证时核对参数、方向和旋转中心。

默认验证入口为 `./scripts/check.sh`。

# 多对话开发的资源约定

用户同时运行多个开发对话，优先保证整机稳定、桌面响应和已有任务进度，避免 OOM 和 CPU 争抢。

- 启动编译、批量测试、基准测试、训练或大型数据处理前，检查 `free -h`、`uptime` 和相关进程；按整机当前余量分配资源，不把每个对话都当作独占机器。
- 编译默认最多 2 个 worker；测试默认串行。显式使用 `cmake --build ... --parallel 2`、`make -j2`、`ninja -j2` 等，禁止裸 `-j`、`-j$(nproc)` 和自动占满 CPU。ROS/colcon 同时限制包级和包内并行，默认包级 1、包内 2。
- 对计算任务显式设置 `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`，以及适用的 `CMAKE_BUILD_PARALLEL_LEVEL=2 CTEST_PARALLEL_LEVEL=1 MAX_JOBS=2 RAYON_NUM_THREADS=2`；同时限制程序自身的 worker 数，避免嵌套线程池放大并发。
- 编译、批量测试和其他重任务使用同一个跨对话排他锁 `$HOME/.cache/codex-heavy-work.lock`，参与约定的对话在整机上同时只运行一项重任务。示例：

  ```bash
  flock -n -E 75 "$HOME/.cache/codex-heavy-work.lock" env CMAKE_BUILD_PARALLEL_LEVEL=2 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 nice -n 10 cmake --build build --parallel 2
  ```

  重任务全程持有锁；不得在锁内后台启动任务后立即退出。取锁失败表示资源正忙，应继续轻量工作或稍后重试，不能绕过锁，也不能把取锁失败报告成构建/测试失败。获得锁后仍需检查资源，锁不控制未采用本约定的进程。
- 为系统和其他对话保留至少 `max(12 GiB, 总内存的 20%)` 可用内存。若低于该值，或 1 分钟 load average 达到可用逻辑核数的 75%，先推迟新重任务；同时参考 CPU 实际利用率、swap 活动和 PSI 判断压力。负载升高时优先降低自己任务的并发或暂停自己的重任务。
- 大数据先估算峰值内存，优先分块、流式处理和小样本验证。只运行与本次改动相关的必要检查，避免重复构建、同时跑多套完整测试和无期限压力测试；长任务设置合适的超时并关注执行期间的资源变化。
- 未经用户明确要求，不终止、暂停或修改其他对话及用户已有进程，不调整全局内核、swap 或系统级 CPU/内存限制。结束自己启动的临时服务和闲置后台任务。
- 这些是后续任务读取的协作规则，不是操作系统硬配额；已有会话是否重新读取规则不能假定。汇报时清楚区分规则已保存与所有运行中对话已实际限流。
