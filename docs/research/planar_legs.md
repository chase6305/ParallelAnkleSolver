# 平面并联腿、闭链膝与可复用求解模型

核查日期：2026-11-12。范围限于作者论文、作者/机构仓库、实际源码及其许可证；未下载大型模型、未导入第三方实现、未运行外部代码。

本轮已独立实现并验证通用平面五杆后端 `PlanarFiveBar`，实际接口见[后端说明](../backends.md)。平面五杆的任务坐标是足端二维位置；闭链膝通常是一个电机到一个关节的非线性传动；Cassie 一类机构还包含被动自由度与弹性约束。这些模型共享数值验证原则，几何闭环、坐标和分支分别定义。下列资料为后续扩展提供参考，外部程序未在本仓库运行。

## 1. 通用平面五杆：polypantograph

- 作者技术报告：[Kinematic Modeling, Analysis, and Real-Time Implementation of a Planar Five-Bar Pantograph Robot](https://doi.org/10.5281/zenodo.20616265)，Abolfazl Mohebbi，2026。
- [作者仓库](https://github.com/polar-polymtl/polypantograph)。README 明确声明它是该报告及另一篇动力学报告的配套软件；Zenodo 将运动学文档分类为 `Other`，本记录不把它描述为已核实的同行评审论文。
- 核查版本：[9c7dceb2845575bc38132d32fe1c9e792b28d8a2](https://github.com/polar-polymtl/polypantograph/tree/9c7dceb2845575bc38132d32fe1c9e792b28d8a2)。
- 代码入口：[fkine.py](https://github.com/polar-polymtl/polypantograph/blob/9c7dceb2845575bc38132d32fe1c9e792b28d8a2/fkine.py)、[ikine.py](https://github.com/polar-polymtl/polypantograph/blob/9c7dceb2845575bc38132d32fe1c9e792b28d8a2/ikine.py)、[jacobians.py](https://github.com/polar-polymtl/polypantograph/blob/9c7dceb2845575bc38132d32fe1c9e792b28d8a2/jacobians.py)。
- 许可证：[PolyForm Noncommercial 1.0.0](https://github.com/polar-polymtl/polypantograph/blob/9c7dceb2845575bc38132d32fe1c9e792b28d8a2/LICENSE)。源码头部也标明仅限非商业使用，仓库说明商业用途需要作者另行授权；不能按 MIT/BSD 类项目直接引入商业求解库。

适合研究的内容是：不等长两支链、显式正解装配分支、每侧逆解工作模式、约束雅可比 A/B、奇异性和条件数。它研究的是两个主动转动关节驱动足端二维位置的五杆闭链，不能仅因为名称含 pantograph 就与被动比例缩放连杆等同。

## 2. 真实五杆腿的特例：Stanford Doggo

- 论文：[Stanford Doggo: An Open-Source, Quasi-Direct-Drive Quadruped](https://arxiv.org/abs/1905.04254)，ICRA 2019。论文参考文献 [20] 直接指向下列项目仓库，论文与项目关联已核实。
- [项目仓库](https://github.com/Nate711/StanfordDoggoProject)，核查版本 [ce5716b5941739a006ec604669bca5676f265dc1](https://github.com/Nate711/StanfordDoggoProject/tree/ce5716b5941739a006ec604669bca5676f265dc1)，[顶层 MIT 许可证](https://github.com/Nate711/StanfordDoggoProject/blob/ce5716b5941739a006ec604669bca5676f265dc1/LICENSE)。
- 项目固定的 ODrive 子模块版本为 [ec25d02f7036b0d7eff010fdc516bd0054a86c56](https://github.com/Nate711/ODrive/tree/ec25d02f7036b0d7eff010fdc516bd0054a86c56)。实际运动学/力映射代码在 [Firmware/MotorControl/controller.cpp](https://github.com/Nate711/ODrive/blob/ec25d02f7036b0d7eff010fdc516bd0054a86c56/Firmware/MotorControl/controller.cpp)，该目录有独立的 [MIT 许可证](https://github.com/Nate711/ODrive/blob/ec25d02f7036b0d7eff010fdc516bd0054a86c56/Firmware/LICENSE)。

Doggo 的两个驱动轴**同轴**，代码使用两电机角的平均值和半差定义虚拟腿方向 `theta` 与开合角 `gamma`。控制器中的足端位置、腿长导数及力到力矩的映射利用了该对称几何。它可作为通用五杆后端的**同轴、等长支链回归特例**，不能直接当成支持任意基座间距、任意四根杆长的通用算法。

该固件主要提供特定分支的前向几何与控制映射，并未在这里核实到完整的通用五杆多分支 IK 接口。另一个 `Nate711/Doggo` 控制程序子模块没有查到根目录 LICENSE；顶层项目许可证不应自动推广为所有子模块都使用同一许可证。

## 3. 多种闭链与复杂腿：LinkageKinematicModels

- [作者仓库](https://github.com/LionelBirglen/LinkageKinematicModels)，Lionel Birglen，Polytechnique Montréal。
- 核查版本：[031b9b711c03387082e2eadf9e231bcc9245fcbb](https://github.com/LionelBirglen/LinkageKinematicModels/tree/031b9b711c03387082e2eadf9e231bcc9245fcbb)。
- [GNU AGPL v3 许可证](https://github.com/LionelBirglen/LinkageKinematicModels/blob/031b9b711c03387082e2eadf9e231bcc9245fcbb/LICENSE)。这是强 copyleft 软件，不作为本库的默认代码依赖或复制来源。
- 通用五杆入口：[fivebar_direct_kinematics.m](https://github.com/LionelBirglen/LinkageKinematicModels/blob/031b9b711c03387082e2eadf9e231bcc9245fcbb/Matlab/fivebar_direct_kinematics.m)、[fivebar_inverse_kinematics.m](https://github.com/LionelBirglen/LinkageKinematicModels/blob/031b9b711c03387082e2eadf9e231bcc9245fcbb/Matlab/fivebar_inverse_kinematics.m)。前者返回两个圆交点装配解；后者枚举最多四个主动关节解。该 MATLAB 接口使用角度制，不能直接按本库 rad 接口代入。

仓库还包含四杆、滑块曲柄、Stephenson III 和 Q2 腿。其 README 将 Q2 实现明确关联到：

1. [Kinematic Analysis of a Biocompatible Lower Limb Model](https://doi.org/10.1007/978-3-031-45705-0_29)，Birglen、Hely，2023。
2. [A Compliant Q2 Linkage to Model the Human Lower Limb Motion](https://doi.org/10.1007/978-3-031-95489-4_11)，Birglen，2025。

Q2 腿是一个转动输入与一个伸缩输入共同驱动的复杂多闭环机构。仓库说明正解返回 8 个固定槽位装配模式并分别标记有效性；Q2 逆解仍列为后续工作。这能支持“后端能力声明”的设计：不能因为某个模型有 FK 和几何显示，就统一宣称其已有 IK、速度映射和任意分支追踪。仓库 README 的语言功能表与部分说明存在版本不同步，使用时需逐个检查实际文件。

## 4. 真实机器人参考的适用边界

| 平台与一手材料 | 已核实的可用内容 | 与通用两主动自由度五杆的区别 |
| --- | --- | --- |
| [Ascento 机构论文](https://arxiv.org/abs/2005.11435)、[LQR-Assisted Whole-Body Control of a Wheeled Bipedal Robot with Kinematic Loops](https://arxiv.org/abs/2005.11431) | 后一篇明确研究每条腿的四杆闭环及闭环动力学表达 | 每条腿只有一个腿部执行器，另外一个执行器驱动车轮；不是两个腿部执行器控制平面足端位置。本次未核实到论文官方独立 FK/IK 库和相应代码许可证，不提供未经确认的代码配对。 |
| [Cassie 设计论文](https://leggedrobotlocomotion.wordpress.com/wp-content/uploads/2017/07/abatea.pdf)、[官方运动学说明](https://github.com/agilityrobotics/cassie-doc/wiki/Kinematic-Model) | 真实腿部机构及关节/闭环定义 | 包含被动连接、弹性元件与空间连接杆，不能将被动角误当作独立电机输入。 |
| [Agility 提供的 Cassie MJCF 在 MuJoCo Menagerie 中的模型](https://github.com/google-deepmind/mujoco_menagerie/tree/370c43965323c0944b433b202aac5058b5f907a9/agility_cassie) | [README](https://github.com/google-deepmind/mujoco_menagerie/blob/370c43965323c0944b433b202aac5058b5f907a9/agility_cassie/README.md) 明确来源；[模型](https://github.com/google-deepmind/mujoco_menagerie/blob/370c43965323c0944b433b202aac5058b5f907a9/agility_cassie/cassie.xml) 保留 Achilles/plantar 连杆的四个 connect 闭环；[MIT 许可证](https://github.com/google-deepmind/mujoco_menagerie/blob/370c43965323c0944b433b202aac5058b5f907a9/agility_cassie/LICENSE) 已核实 | 适合作为独立仿真验证模型；不是本次已核实的通用解析运动学库，也不声称是上述 2017 设计论文的配套实现。 |
| [LimX Oli 官方描述](https://github.com/limxdynamics/humanoid-description/tree/4c0a3c71200b810a52831c7c2b2de0ed6ce819a4)、[官方并联关节映射说明](https://github.com/limxdynamics/humanoid-rl-deploy-python/blob/main/doc/parallel_joint_mapping_en.md) | 官方描述以 A/B 电机与等式约束建模踝部和腰部；[Apache-2.0 许可证](https://github.com/limxdynamics/humanoid-description/blob/4c0a3c71200b810a52831c7c2b2de0ed6ce819a4/LICENSE) 已核实 | 是空间 pitch/roll 与 A/B 电机的耦合，不是平面五杆足端模型。映射说明中 PR↔AB 转换由 SDK 内部处理；不能据此声称 SDK 求解器实现已经公开。 |

Digit 在本次定向核查中没有形成比上述条目更清楚的“论文—独立运动学实现—许可证”证据链，因此暂不列作可直接集成的求解代码来源。

## 5. 适合本仓库的独立实现路线

以下说明当前五杆实现及后续机型接入要求。原生实现从杆长与圆交点独立推导，未移植所列外部程序。

### 平面五杆后端

`FiveBarParams` 包含两个二维固定支点、左右各自的近端/远端杆长，支持非同轴和不等长支链。任务坐标为二维位置；全部长度采用调用方选定的同一单位。电机角按全局 +X 的逆时针方向定义，编码器正方向与零位由调用方显式转换。

接口分别提供：

- `ForwardAll(motor_angles)`：枚举足端两个装配解，退化/重合情况有明确状态。
- `InverseAll(foot_position)`：枚举每侧肘形组合，去重相切解。
- `Analyze(foot_position, motor_angles)`：闭合配置的 A/B 约束雅可比、工作/装配分支及两类奇异性诊断。

调用方从 `solutions` 中按返回的 `branch` 选择装配/工作模式；需要被动关节角时，由适配层从连接点几何恢复。连续轨迹结合上一帧、关节行程和奇异性决定分支。

**装配模式与工作模式应分开记录。** 固定两个主动角后，左右肘端的两圆通常给出两个足端交点，这是正解装配模式。固定足端位置后，每侧二连杆通常各有两个肘形，组合成逆解工作模式。给一个统一的 `elbow_up` 布尔值不足以表达所有解。

设闭环为 `phi(x, q_motor) = 0`，定义 `A = d(phi)/d(x)`、`B = -d(phi)/d(q_motor)`，则 `A x_dot = B q_motor_dot`。两者的奇异性有不同几何含义，应保留 A/B 的诊断信息，避免只返回最终雅可比的行列式而丢失原因。具体符号约定需在本库接口中固定，并以功率守恒测试验证力映射。

### 与其他机构共用的层

- 四杆膝传动使用独立的标量输入/输出模型，返回非线性传动比及分支；不要把长度运动输出硬解释为 pitch/roll。
- 现有双电机空间脚踝保留二维角度任务空间，未来可进一步参数化电机轴与转动链。
- 被动/弹性闭环可使用独立坐标与从属坐标的约束接口；弹性自由度的力学方程不应被纯几何求解器消除。
- 共用状态枚举、分支 ID、奇异性阈值、收敛信息和线性求解工具；各模型声明实际支持的 FK/IK/雅可比能力。

### 机型接入的几何验证要求

1. 左右不等长、非零基座间距和旋转后的基座，分别验证每条杆长度与闭环残差。
2. FK/IK 所有有效分支往返；测试相切时去重、两圆重合、不可达和近奇异构型。
3. 同轴等长情形与 Doggo 的虚拟腿解析关系交叉检查，作为特例回归。
4. 雅可比与独立有限差分比较，速度和力映射验证功率守恒。
5. 连续路径保持分支；不能跨越奇异点时返回明确状态，不静默跳到另一个装配解。

当前原生后端已覆盖不对称、基座变换、多分支、相切、重合、不可达、差分导数及 Doggo 几何回归。实际机型的编码器标定、关节行程和连续分支选择仍由适配层提供。未复制非商业或 AGPL 仓库代码。
