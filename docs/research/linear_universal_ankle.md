# 两条线性驱动支链与中心万向节：空间脚踝组合示例

示例文件：`examples/linear_universal_ankle_demo.cpp`。它组合现有 `SolvePrismaticLength`、`SolveEqualities` 和 `DifferentiateImplicitly`，没有新增核心求解 API，也没有引入外部机器人求解库。构型依据公开的中心 U 约束支链与两条可变长支链；尺寸为本仓库自主给定的 SI 测试夹具，**不是 HyRoDyn、RH5 或其他硬件的导入尺寸**。

## 文献与源代码对应关系

HyRoDyn 官方固定提交 [`d2fc248eb0270b1251a6233208c86df9b5776a8e`](https://github.com/dfki-ric/hyrodyn/tree/d2fc248eb0270b1251a6233208c86df9b5776a8e) 使用模块名 `2SPU+1U`。但其 [SPU2U1.hpp](https://github.com/dfki-ric/hyrodyn/blob/d2fc248eb0270b1251a6233208c86df9b5776a8e/src/SPU2U1.hpp) 按基座到平台的方向实际描述 **U-P-S** 支链，并将中心万向节的两个角作为独立坐标。这种反向命名不能用于推断 URDF 中关节的实际顺序。

其 [SPU2U1.cpp](https://github.com/dfki-ric/hyrodyn/blob/d2fc248eb0270b1251a6233208c86df9b5776a8e/src/SPU2U1.cpp) 的几何逆模型明确采用 `[roll,pitch]`、`Rx*Ry`、旋转平台连接点及两点距离求驱动长度。本文独立推导相同的约化长度关系，不复制其实现。HyRoDyn 全模型还恢复支链被动关节，当前示例不包含这一层。框架的模块化解析/数值方法可对照 [DFKI 的 2022 论文说明](https://robotik.dfki-bremen.de/en/research/publications/12827)；外部代码采用 [BSD-3-Clause](https://github.com/dfki-ric/hyrodyn/blob/d2fc248eb0270b1251a6233208c86df9b5776a8e/LICENSE)。

## 坐标与自主夹具

中心 U 关节的两根轴正交且共点于 `O`。基座坐标系中第一根轴为 X；第二根为经第一关节带动后的 Y 轴。平台零位坐标系与基座同向，采用主动旋转、列向量：

```text
q = [roll, pitch] rad
R(q) = Rx(roll) Ry(pitch)
C_i(q) = O + R(q) p_i
d_i(q) = C_i(q) - A_i
ell_i(q) = ||d_i(q)||
```

旋转脚踝后端使用 `[pitch,roll]`、`Ry*Rx`。跨后端移植时应重新推导有限旋转和雅可比。中心 U 已通过姿态参数化约束平台运动，平台原点不允许自由平移或独立第三轴转动。

| 参数 | 值，单位 m |
| --- | --- |
| `O` | `(0,0,0)` |
| 平台连接点 `p1/p2` | `(.04,.03,0)` / `(.04,-.03,0)` |
| 基座连接点 `A1/A2` | `(.04,.03,.20)` / `(.04,-.03,.20)` |
| 两支链长度范围 | `[.15,.25]` |
| 本示例局部关节区间 | roll、pitch 均为 `[-.3,.3] rad` |

两执行器输出为**长度**。这里只表示理想支链长度；真实编码器的符号、零位、螺距、减速比及伸缩缸参考长度需要另外标定。特别是绝对长度与相对零位伸长量不能混用。

## 几何逆解、闭环正解与导数

给定 q 时，分别向 `SolvePrismaticLength` 提供基座点 `A_i`、空间点 `C_i(q)` 和行程，得到 `ell_i` 与点位置梯度。记 `e_x/e_y` 为单位轴，连接点导数为：

```text
C_i,roll  = e_x cross (R p_i)
C_i,pitch = Rx [e_y cross (Ry p_i)]
d(ell_i)/dq = (d_i/ell_i)^T [C_i,roll  C_i,pitch]
```

给定两条长度时，局部 FK 以 q 为两个未知量，使用原始杆长约束：

```text
F_i(q,ell) = (||C_i(q)-A_i||^2 - ell_i^2) / 2 = 0
F_q[i,:]   = d_i^T [C_i,roll  C_i,pitch]
F_ell      = -diag(ell1,ell2)
```

`SolveEqualities` 接收 F 和 `F_q`、有效区间及 warm start。F 的单位为 m²，残差尺度取 `.20^2 m²`，角度变量尺度取 `.3 rad`。将执行器长度单独交给 `DifferentiateImplicitly`，可在已闭合且相关分块满秩时得到：

```text
F_ell (dell/dq) = -F_q
F_q (dq/dell)   = -F_ell
```

对线性支链力 `f`（N），虚功关系为 `tau_q=J_ell_q^T f`（Nm）、`ell_dot=J_ell_q q_dot`（m/s）。长度为零、支链行程不满足或 `F_q` 秩亏时，不能继续使用唯一逆映射。

## 独立核对与有限范围

示例具有可供 CTest 使用的退出状态和以下有限自检：

- 零位解析值：`ell=[.20,.20]`，`J_ell_q=[[-.03,.04],[.03,.04]] m/rad`。
- 非零双轴姿态使用独立标量 `Rx*Ry` 展开和距离公式核对长度，并以其中心差分核对雅可比。
- 六个有限姿态、零初值、另一个初值及上一姿态 warm start 核对局部 FK 与原始长度闭环。
- `q_dot=[.2,-.1] rad/s`、`f=[100,150] N` 的零位解析值为 `ell_dot=[-.010,.002] m/s`、`tau_q=[1.5,10] Nm`，两侧功率均为 `-.7 W`；其他采样点也验证非对角速度/力映射。
- 超行程输入、区间外初值和区间内不可实现的长度目标均须明确失败。非法初值须在首次求值前返回 `InvalidInput`；行程内的不可行目标须进入几何失败处理。

示例的不可行目标为 `[.15,.25] m`。它的不可行性可独立证明：反三角不等式给出 `|ell1-ell2| <= ||R(p1-p2)-(A1-A2)|| = .12 |sin(roll/2)|`。在本角域内其上界约为 `.01793 m`，小于请求的 `.10 m` 长度差；不依赖局部迭代是否找到解来认定不可行。

这是理想关节、刚性连接的约化空间几何例子。它不枚举所有装配模式，不模拟硬件，不恢复全部 U/S 轴承角，不验证轴承角限位、碰撞或柔顺性。局部 FK 的成功仅表示在本区间和初值附近找到闭合解，不能作为全空间或全类型机构的唯一性证明。

当前可选 `second_order.h` 的公共接口绑定 `PMAnkle`，不能直接用于本示例的线性执行器。若后续需要长度传动 Hessian 与预载刚度，须单独推导 `ell_i(q)` 的二阶导数、以 N/m 与 Nm/rad 区分执行器/任务刚度，并进行独立势能验证后再接入。
