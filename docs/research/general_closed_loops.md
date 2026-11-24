# 通用闭环机构：论文、开源实现与接入建议

核查日期：2026-11-12。以下仅核查作者论文、项目仓库、官方文档与源码；未安装依赖、克隆项目或运行其算法。因此“支持”表示公开材料和代码提供相应接口，不表示已在本仓库复现，也不代表任意机构存在唯一、全局可求的运动学解。

针对“适配不同腿部串并联构型”，优先研究 **Gepetto 的 toolbox-parallel-robots + example-parallel-robots**：它们同时提供主/被动关节划分、闭环运动学工具和多种腿部模型。需要图结构编译时研究 RoboCompiler；需要高效嵌入闭环动力学时研究 GRBDA。Pinocchio 的约束动力学论文提供基础算法，但约束动力学与给定电机位置求机构装配的 FK 是不同问题。

## 1. Gepetto 工具箱、模型库与 ClosedLoopMotion

[Optimal Control of Walkers with Parallel Actuation](https://arxiv.org/abs/2504.00642) 将闭环约束及其导数直接纳入腿式机器人最优控制，而不是先用近似串联腿替代原机构。[toolbox-parallel-robots](https://github.com/Gepetto/toolbox-parallel-robots) 提供装配、投影、主/被动关节划分、闭环 FK/IK 和闭环雅可比；[模型库](https://github.com/Gepetto/example-parallel-robots) 包含 Cassie、Digit 单腿/双腿、Kangaroo、WL16、五杆和 Delta，且作者明确说明模型基于公开信息，不替代官方精确模型。

它适合先做跨构型参考后端。[雅可比实现](https://github.com/Gepetto/toolbox-parallel-robots/blob/d7b093e/toolbox_parallel_robots/jacobian.py) 将闭环约束列分为电机与被动部分，通过伪逆构造局部速度映射；这不自动保证奇异点可逆或解唯一。[FK 源码](https://github.com/Gepetto/toolbox-parallel-robots/blob/d7b093e/toolbox_parallel_robots/forward_kinematics.py) 采用 SciPy SLSQP 或 CasADi/Ipopt，依赖初值；其中部分失败路径仍返回调试解，接入时必须独立检查闭环残差和状态。

[ClosedLoopMotion](https://github.com/LudovicDeMatteis/ClosedLoopMotion) 是对应运动生成示例，其 [环境文件](https://github.com/LudovicDeMatteis/ClosedLoopMotion/blob/60516ee/flake.nix) 还依赖特定 Crocoddyl/sobec 分支，负担明显高于单独使用工具箱。其 README 的 ICRA 投稿说明不能单独用来认定最终发表会议。

| 项目 | 核查到的提交前缀 | 授权与依赖 |
| --- | --- | --- |
| toolbox-parallel-robots | [`d7b093e`](https://github.com/Gepetto/toolbox-parallel-robots/commit/d7b093e) | [BSD-3-Clause](https://github.com/Gepetto/toolbox-parallel-robots/blob/d7b093e/LICENSE)；[包元数据](https://github.com/Gepetto/toolbox-parallel-robots/blob/d7b093e/pyproject.toml) 要求 Python ≥3.8、Pinocchio ≥3.7；源码另使用 NumPy/SciPy，部分模块需要 CasADi/Ipopt。 |
| example-parallel-robots | [`147faa5`](https://github.com/Gepetto/example-parallel-robots/commit/147faa5) | [BSD-3-Clause](https://github.com/Gepetto/example-parallel-robots/blob/147faa5/LICENSE)；使用前仍需核对所选模型精度与模型内资产说明。 |
| ClosedLoopMotion | [`60516ee`](https://github.com/LudovicDeMatteis/ClosedLoopMotion/commit/60516ee) | 核查时未发现顶层 LICENSE，不能将相邻项目的 BSD 授权自动套用到该仓库；先作为论文复现参考。 |

## 2. RoboCompiler：从物理图到局部运动/力接口

[RoboCompiler 论文 v2](https://arxiv.org/abs/2609.35717v2) 于 2026-09-29 更新，并直接链接 [作者代码](https://github.com/Mehdi-Heydari-Shahna/RoboCompiler)，[代码许可证为 Apache-2.0](https://github.com/Mehdi-Heydari-Shahna/RoboCompiler/blob/main/LICENSE)。输入包含刚体、关节、连接坐标系、惯量和执行器端口；输出闭环路径、解析残差雅可比、局部独立速度提升映射，以及虚功一致的驱动力和约化动力学接口。示例包括闭链腿 Kangaroo 与六 UPS Stewart 平台。

这是通用模型接口的强参考，但它明确依赖局部装配分支、秩条件及延拓路径，不能解释为任意输入图的全局解析解。[README](https://github.com/Mehdi-Heydari-Shahna/RoboCompiler/blob/main/README.md) 提供分散的研究包入口，参考环境为 Linux/Python 3.12，主要流程固定 Pinocchio 3.8.0、MuJoCo 3.3.7；Isaac Sim 需另行部署，挖掘机完整模型未包含。宜先阅读 Kangaroo/Stewart 的图记录与求解状态设计，再评估适配；本次未取得可固定的提交 SHA，以上为核查日 main 快照。

## 3. Pinocchio：通用约束动力学基础

[Proximal and Sparse Resolution of Constrained Dynamic Equations](https://www.roboticsproceedings.org/rss17/p017.html) 是 RSS **2021** 论文，研究带闭环/接触约束的正向动力学、奇异约束下的近端迭代和稀疏分解。[LCABA 论文项目页](https://simple-robotics.github.io/publications/lcaba/) 给出进一步面向内部环路的递归正向动力学算法。两者都公开关联 [Pinocchio](https://github.com/stack-of-tasks/pinocchio)。

其 [四杆闭环示例](https://github.com/stack-of-tasks/pinocchio/blob/501dec5/examples/simulation-closed-kinematic-chains.py) 展示开树加闭环约束的建模方式。约束模型/雅可比很适合做参考适配器，但这些动力学算法本身不是“电机位置 → 全部被动关节位置”的全局装配器。切断处的约束自由度必须匹配实际关节，不能把所有闭环都刚性焊接。核查提交前缀为 [`501dec5`](https://github.com/stack-of-tasks/pinocchio/commit/501dec5)，[BSD-2-Clause](https://github.com/stack-of-tasks/pinocchio/blob/501dec5/LICENSE)。[构建选项](https://github.com/stack-of-tasks/pinocchio/blob/501dec5/CMakeLists.txt) 可裁剪解析器、碰撞与 Python 接口；直接引入完整 Python 环境仍明显重于当前 Eigen 核心。

## 4. GRBDA 与 URDF+：局部闭环作为广义关节

[A Propagation Perspective on Recursive Forward Dynamics for Systems with Kinematic Loops](https://arxiv.org/abs/2311.13732) 将约束嵌入解释为刚体簇之间的广义关节与运动/力子空间，针对齿轮、差动与四杆等执行器局部机构很贴近腿部应用。[URDF+ 论文](https://arxiv.org/abs/2411.19753) 进一步提供环路关节描述和分组解析方法。[GRBDA 开源代码](https://github.com/ROAM-Lab-ND/generalized_rbda) 提供 C++ 约束嵌入动力学和相应模型接口。

优势是保留真实执行器子机构并利用局部环路结构；其重点是动力学传播，不应当成无条件适用所有拓扑的通用 FK/IK 包。簇的内部运动学、所选独立变量及分支仍必须可用，跨越大范围的非局部环路也可能扩大簇和计算量。核查提交前缀为 [`e68a6ec`](https://github.com/ROAM-Lab-ND/generalized_rbda/commit/e68a6ec)，[MIT](https://github.com/ROAM-Lab-ND/generalized_rbda/blob/e68a6ec/LICENSE)。[README](https://github.com/ROAM-Lab-ND/generalized_rbda/blob/e68a6ec/README.md) 与 [CMake](https://github.com/ROAM-Lab-ND/generalized_rbda/blob/e68a6ec/CMakeLists.txt) 显示 C++17、Eigen、CasADi 和 URDF 解析依赖，因此建议在需要完整动力学时再接入。

## 5. PlaCo：任务约束形式的实用数值运动学

[PlaCo 论文](https://arxiv.org/abs/2511.06141) 对应 [C++/Python 开源框架](https://github.com/Rhoban/placo)。其 [闭环运动学文档](https://placo.readthedocs.io/en/latest/kinematics/loop_closures.html) 明确展示将机构切成树，再用硬相对位置约束闭环，并通过轴掩码表达平面约束；可以分别配置电机关节控制和任务空间控制的求解器。它对新增不同几何机构的原型验证很方便，但属于依赖初值的局部数值求解，不枚举全部装配分支。

[闭环动力学文档](https://placo.readthedocs.io/en/latest/dynamics/loop_closures.html) 还要求为闭环任务添加内部力，并将被动自由度约束为零驱动力矩；仅添加几何任务并不等同于完成正确动力学建模。核查提交前缀为 [`93adb57`](https://github.com/Rhoban/placo/commit/93adb57)，[MIT](https://github.com/Rhoban/placo/blob/93adb57/LICENSE)。其 [Python 包依赖](https://github.com/Rhoban/placo/blob/93adb57/pyproject.toml) 包含 Pinocchio、URDF 解析、Meshcat 等，适合作为可选原型后端，避免成为本库基本编译依赖。

## 适配到本仓库的建议

以下为根据上述接口作出的工程适配路线。当前已实现欧氏等式残差/雅可比回调、局部求解、尺度/边界和秩/相容性诊断，接口及验证范围见[后端说明](../backends.md)。外部框架适配器、自动图编译和非欧氏配置积分尚未接入；下面相关条目仍是后续建议：

1. **先增加可选参考适配器。** 用 Gepetto 模型中的五杆、Digit 或 Kangaroo 建立外部基准，传入开树模型、闭环约束和执行器选择。外部 Python/Pinocchio 环境独立安装，本库默认 C++17/Eigen 路径保持轻量。
2. **核心先抽象约束，而不是泛化固定 2×2 脚踝公式。** 原生模型接口应能返回约束残差、约束雅可比、变量/执行器映射和配置积分操作。对于姿态四元数，配置导数与广义速度不是同一向量，不能直接套用欧氏加法。
3. **显式保存秩与兼容性信息。** 在局部有效的独立/被动划分中，求解 `J_passive * D = -J_independent`，再验证完整约束残差和 `J * E ≈ 0`。伪逆只能给出指定意义下的最小二乘解；欠约束、矛盾约束或秩变化时，应返回诊断并要求重新选变量或给出次级目标。
4. **独立速度不总是电机坐标。** 对直线缸、耦合电机和复杂差动，还应提供执行器端口函数及其微分；以虚功关系验证力矩/推力映射。非方阵情况下需要驱动力分配或可实现性求解，不能直接求逆。
5. **将“通用”限定为明确的约束类型和局部分支。** 用闭环误差、杆长/铰链条件、正逆解往返及功率守恒验证每种构型；增加独立分支标识、初值/延拓和失败状态后，再考虑动力学与自动图编译。

这些局部速度映射的依据可以直接对照 [Gepetto 的雅可比实现](https://github.com/Gepetto/toolbox-parallel-robots/blob/d7b093e/toolbox_parallel_robots/jacobian.py) 以及 [RoboCompiler 的局部秩与提升映射推导](https://arxiv.org/html/2609.35717v2#S3.SS2)。它们都不提供“仅靠任意 URDF，便能保证得到任意机构唯一全局解”的结论。

## 未作为首选的相关方向

- **Drake** 可作为显式约束与优化工具的后续候选，但 [官方 `SetEnableLoopTopology` 文档](https://drake.mit.edu/pydrake/pydrake.multibody.plant.html#pydrake.multibody.plant.MultibodyPlant.SetEnableLoopTopology) 在核查时仍说明自动闭环拓扑处理处于开发中，不能据此承诺自动导入任意闭链图。
- **JaxSim** 的 [作者 README](https://github.com/gbionics/jaxsim) 已明确列出闭链动力学支持，不能笼统说它只支持开树；但本次没有进一步核实各闭环约束类型和装配接口，暂不把它列为通用腿部 FK 方案。其定位是 JAX 可微仿真，依赖与编译模式也不同于当前 C++ 小型求解库。
- 未核实到与本任务匹配、且论文与开源代码能对应的 **CLeaR** 项目；不将同名的数据关联或多机器人通行论文混入闭链机构求解候选。
