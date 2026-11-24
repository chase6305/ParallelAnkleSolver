# 串并联混合机构的建模与拓扑雅可比

核查日期：2026-11-12。以下是作者公开资料的核查结果；尚未在本仓库复现这些第三方程序。

## TriP：模块化串并联混合建模

- 论文：Baumgärtner、Miller，2022，[JOSS 7(71), 3967](https://joss.theoj.org/papers/10.21105/joss.03967)。
- 作者代码：[TriPed-Robot/trip_kinematics](https://github.com/TriPed-Robot/trip_kinematics)。
- 代码许可证：[MIT](https://github.com/TriPed-Robot/trip_kinematics/blob/main/LICENSE)。论文许可证与代码许可证分别核查。
- 作者文档：[TriP documentation](https://trip-kinematics.readthedocs.io/en/main/)。

论文将并联子机构建模为 `KinematicGroup`，串联连接用 `Transformation` 表示，并把两者组合成机器人。作者提供正/逆运动学、雅可比与位姿目标接口，示例包含 TriPed 和挖掘机。论文明确定位于机构建模与设计研究，不以硬实时控制为目标。

对本仓库的启发是按机构族提供独立求解器，再组合成腿部运动链。一个子机构可以有解析法或数值法；模型需要明确主动关节、被动关节、装配分支以及任务坐标。上述组合方式是本仓库拟采用的设计思路，不表示已经接入 TriP。

## ACRoD：从拓扑构建闭环雅可比

- 论文：Jacob、Datta，2024，[JOSS 9(95), 5927](https://joss.theoj.org/papers/10.21105/joss.05927)。
- 作者代码：[suneeshjacob/ACRoD](https://github.com/suneeshjacob/ACRoD)。
- 代码许可证：[GPL-3.0](https://github.com/suneeshjacob/ACRoD/blob/main/LICENSE)。
- 作者说明：[README](https://github.com/suneeshjacob/ACRoD/blob/main/README.md)。

ACRoD 根据机器人拓扑矩阵，生成给定构型、指定末端点的平面或空间雅可比。作者将其用于非冗余机构的尺寸设计和拓扑分析。示例包括平面串并联混合链及空间 RSSR-SSR 机构。它可以帮助交叉检查约束雅可比的构造，但不能据此宣称已获得任意闭链机构的全局 FK/IK。

本仓库首先独立实现几何闭环与约束分块的数值接口；ACRoD 可作为离线验证候选。当前没有复制其源码，也没有把它作为核心库依赖。
