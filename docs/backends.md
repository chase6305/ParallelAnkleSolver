# 原生求解后端与坐标约定

旋转脚踝、五杆、支链与显式约束后端使用 Eigen。各后端分别定义输入、输出、分支和单位。

## 曲柄—连杆脚踝

接口见 `parallel_ankle/parallel_ankle.h`。任务坐标为 `[pitch, roll] rad`，电机输出为接线顺序 `[motor1, motor2] rad`。默认雅可比为 `d(theta)/d(q)`，角度单位相消；模型的所有长度统一缩放不改变该雅可比。

该模型使用固定的 X 轴曲柄平面和 `Ry(pitch) Rx(roll)` 运动链，允许两个旋转中心不重合。该模型求解指定物理侧的固定曲柄分支。主值分支之外、物理侧变化或输出角无法表达所选构型时返回失败。显式三参数 `ComputeJacobian(pitch, roll, eps)` 为迁移兼容的数值差分接口，默认双参数版本为解析导数。

可选 `second_order.h` 自由函数返回同一分支的两组角度 Hessian，并由电机恢复预载映射完整静态刚度。该接口采用本模型的任务坐标和电机符号。计算失败和矩阵谱标志的含义见[二阶说明](research/second_order_transmission.md)。

## 一般平面五杆

接口见 `parallel_ankle/five_bar.h`，示例为 `five_bar_demo`。两个固定支点 `A1/A2` 可以非同轴、非水平、非对称；每侧近端杆和远端杆长度分别配置。主动角为从全局 +X 方向计量的逆时针绝对角，单位 rad。此后端没有隐含的编码器符号和零位补偿，编码器转换由调用方显式完成。

```cpp
parallel_ankle::FiveBarParams parameters;
parameters.pivots = {Eigen::Vector2d(-0.4, 0.1), Eigen::Vector2d(0.5, -0.1)};
parameters.proximal = {0.8, 0.7};
parameters.distal = {1.0, 0.9};
parallel_ankle::PlanarFiveBar leg(parameters);
const auto inverse = leg.InverseAll({0.1, 1.0});
// 逐一检查 inverse.solutions[0..count)，显式选择 working/assembly 分支。
```

`InverseAll` 返回最多四组工作模式。`ForwardAll` 对给定电机角返回最多两个足端位置，分别对应远端两圆的装配交点；相切解去重。两个圆完全重合时存在无限多个足端位置，返回 `Singular/count=0`，不会凭空选取一个点。某个有限逆解可能恰好落在无限正解集合中，因此不能要求每组有效 IK 都存在唯一有限 FK。

`branch.assembly` 为 `cross(E2-E1, P-E1)` 的符号，`branch.working[i]` 为 `cross(Ei-Ai, P-Ei)` 的符号，零表示分支边界。连续运动应保留所选分支并检查奇异性；仅按角度最近选择不保证不会在奇异点改变装配模式。

`Analyze` 暴露约束分块：

```text
A x_dot = B theta_dot
motor_jacobian = B^(-1) A
```

A 的每行是对应远端杆的单位方向；B 的对角元为该方向与曲柄端速度方向的点积。`actuator_singular` 和 `assembly_singular` 分别报告曲柄/远端杆共线与两远端杆共线。仅发生装配奇异时，`motor_jacobian` 可能仍有限，但不可据此对其求逆。

几何长度单位由参数选择。雅可比单位是 `rad/长度单位`；电机力矩单位为 Nm 时，`J^T tau_motor` 的单位为 `Nm/长度单位`，只有长度为 m 时才是 N。以 mm 输入时应保留这一单位差异，不将 `Nm/mm` 直接标作 N。

输出点与电机角转换为 double 后再次验证杆长。极大平移、极端长短杆比或精度不足时可返回 `NumericalFailure`，而不会把不闭合的浮点近似当作成功。

## 任意轴旋转与直线支链

接口见 `parallel_ankle/transmissions.h`，示例为 `transmissions_demo`。这是可组合的单支链原语，不是完整空间平台的正解器。

`SolveRotaryRod` 使用固定电机中心 A、正交单位基 u/v、曲柄长度 L 和连杆长度 R，求解：

```text
B(gamma) = A + L (u cos(gamma) + v sin(gamma))
||C - B(gamma)|| = R
```

连接点 C 可以处于三维空间任意位置。返回最多两个模 2π 的相位族，并分别给出分支符号和局部梯度 `d(gamma)/dC`；电机方向与零位按模 2π 相位校准。它不会枚举无限多个绝对编码器绕数，也没有隐含的转角行程选择。返回的相位为 `[-pi, pi)`，跨区间边界时应做局部解缠后再作差。

相切时存在一个几何相位，但传动导数不成立；`success()` 可为 true，`derivative_valid` 为 false。连接点位于轴线、方程与相位无关且闭合时，返回 `Continuum/count=0`。不能可靠分类或不能用 double 表示所选分支时返回 `NumericalFailure`，不凭空输出零角。

`SolvePrismaticLength` 返回固定点到连接点的距离、方向梯度和行程可达性。若允许零行程，零长度几何上有效但方向不唯一，因此单独报告 `Singular` 且不提供有效梯度。轴承角限位、螺距/编码器标定和多条支链的平台闭环均由上层模型提供。

可将这些支链的原始杆长方程组合成下述等式约束问题，求多支链脚踝、膝关节或串并联腿的局部位姿。每个模型仍须明确转动顺序、连接坐标系、独立变量和分支。

`linear_universal_ankle_demo` 提供一个完整组合：两条 UPS 长度驱动与中心万向节约束平台，坐标为 `[roll,pitch]`、`Rx*Ry`，长度用 m。示例从长度 IK、带区间和初值的局部 FK，连接到双向隐式雅可比、N/Nm 虚功与独立几何/差分核对。它使用自主夹具，不导入外部硬件尺寸；构型来源和范围见[空间脚踝示例说明](research/linear_universal_ankle.md)。

## 显式等式约束

接口见 `parallel_ankle/constraints.h`，位于 `parallel_ankle::constraints`。它为用户提供的局部欧氏方程 `F(x)=0` 求解，不自动推断机器人拓扑或主动/被动关节。示例 `constraints_demo` 使用两条定长杆约束求两个圆交点之一。

`EqualityProblem` 明确未知量数和方程数，分别提供残差及雅可比回调。回调必须写入全部输出；配置不可求值时返回 false，错误尺寸或非有限输出是回调契约错误。接口记录求值次数、迭代次数、秩、尺度化残差、步长和失败消息。

局部迭代采用带阻尼 SVD 步、回溯和可选的变量区间投影。初值为 warm start，必须满足给定边界。欠约束或冗余系统可以得到小残差，但 `success()` 不表示全局或局部唯一；同时读取 `rank/full_column_rank`。该接口不枚举装配模式，不包含四元数积分或碰撞检查。

### 尺度与容差

`residual_scales` 是每个方程的单位尺度，`variable_scales` 是每个未知量的单位尺度：

```text
r_scaled = F / residual_scales
J_scaled = diag(1/residual_scales) F_x diag(variable_scales)
```

归一化雅可比同时应用行和列尺度。若单独的除法会溢出或落入次正规区间，求解器合并两种尺度的二进制指数，再计算最终元素。

归一化雅可比的非零项均为次正规数时，求解器先统一缩放再分解，以保留相对秩判断。残差范数仍按调用方配置的方程尺度报告。

例如杆长平方方程以 m² 为残差单位，其尺度可以采用特征杆长平方；位置未知量的尺度可采用特征杆长。角度未知量使用相应的 rad 尺度。默认全为 1，不能据此承诺混合单位问题自动具有良好条件数。

### 隐式微分

对已闭合的配置 `F(x,y)=0`，将 `A=F_x`、`B=F_y` 交给 `DifferentiateImplicitly`：

```text
A (dx/dy) = -B
```

A 必须满列秩。冗余但相容的方程可以使用；欠约束分块和不相容右端项返回失败，导数为 NaN。相容性逐个独立方向作真正的相对检查，微小的不相容列不会因其他列很大或整体缩放而自动被忽略。精确零右端项单独处理。

返回值转换回用户变量单位后还会重新验算相容性。若有效方向在 double 中下溢为零或反缩放后失真，返回 `NumericalFailure` 并保留 NaN 导数，避免只验证内部尺度化解。某个新迭代点的回调或数值计算失败时，上一点的秩诊断同时失效。

这是在给定变量分块下的局部导数；秩变化时通常需要重新选择独立变量或重新建立分支，不能把伪逆结果作为唯一物理解继续使用。

## 接入另一种腿部机构

为新机型保留一份独立的建模记录，至少包含：

| 内容 | 需要明确的量 |
| --- | --- |
| 任务与执行器坐标 | 变量名称、顺序、单位，以及编码器符号、零位和绕数约定 |
| 连接几何 | 各点所在坐标系、关节轴、转动顺序、固定长度及来源版本 |
| 闭环与自由度 | 主动/被动关节、每种闭合关节实际约束的自由度、弹性变量 |
| 有效域 | 装配/工作分支、行程、局部坐标区间和可用初值 |
| 数值尺度 | 每个残差和未知量的单位尺度，以及接受闭合误差的依据 |

先按几何选择现有后端：平面两转动驱动五杆使用 `PlanarFiveBar`；已知平台姿态下的固定轴曲柄或长度支链使用 `transmissions.h`。完整空间机构由连接点变换和支链约束组成，可从 [线性驱动空间脚踝示例](research/linear_universal_ankle.md) 的 IK、局部 FK 和双向导数流程开始。

对耦合膝踝或多闭环腿，先写出 `F(x,u)=0`，明确给定坐标 `u` 和待解坐标 `x`。`SolveEqualities` 求局部闭合点，再以 `F_x/F_u` 进行隐式微分。全腿耦合保留相应维数和交叉项；中心关节方向约束通过姿态参数化或闭环方程表达，被动弹性另有力学模型。框架适配参考见[一般闭环资料](research/general_closed_loops.md)。

接入验证使用独立的连接点/杆长重建、已知姿态、导数差分及虚功核对，并覆盖目标分支、不可达和奇异配置。连续反馈先以上一帧作初值，再检查当前分支和残差。通过这些模型级检查后，才将该机型列为已适配；论文索引和自主夹具示例分别保留其研究或示范状态。
