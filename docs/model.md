# 机构模型与映射

两个固定电机中心为 A1、A2；曲柄端为 B1、B2；脚部刚体上的连接点为 C1、C2。AB 为曲柄，BC 为定长连杆。电机驱动曲柄在基坐标系的 YZ 平面内旋转。

定义任务角 `q = [pitch, roll]`，电机角 `theta = [theta1, theta2]`，角度均为 rad。`Ry` 绕基坐标系 Y 轴，`Rx` 绕 X 轴；给定零位连接点 `C0`、零位旋转中心 `O_roll` 和 `O_pitch`：

```text
C = Ry(pitch) Rx(roll) C0
  + Ry(pitch) (I - Rx(roll)) O_roll
  + (I - Ry(pitch)) O_pitch
```

因此 pitch 与 roll 中心可以不重合。中心重合于原点时退化为 `C = Ry Rx C0`。

## 单执行器解析分支

对每个执行器，令 `r = C - A`：

```text
a = r.y
b = -r.z
c = (l_rod² - l_bar² - ||r||²) / (2 l_bar)
D = a² + b²
sine = (b c + |a| sqrt(D - c²)) / D
cosine = (-|a| c + b sqrt(D - c²)) / D
gamma = atan2(sine, cosine)
theta = direction * gamma + alpha
```

`direction` 为 `is_left_motor ? +1 : -1`。固定解析分支要求 `cosine >= 0`，并要求连接点仍处于零位确定的曲柄物理侧。仅对 sine 取 asin 会丢失余弦符号，导致伪解；这里同时检查余弦、杆长和工作根方向。普通路径先用最大几何尺度归一化，避免长度平方溢出。

长曲柄旁的短连杆可能在平方差中被消去，即使微小开角可用 double 表示。普通公式未通过物理复核时，使用长精度的 YZ 圆与球面截面、稳定半角恢复同一个根：

```text
rho = hypot(r.y, r.z)
H = sqrt((l_rod - |r.x|) (l_rod + |r.x|))
d = |l_bar - rho|, S = l_bar + rho
opening = 2 atan2(sqrt((H-d)(H+d)), sqrt((S-H)(S+H)))
gamma = canonical(atan2(r.z, |r.y|) + opening)
```

首先检查平面距离与三角形不等式，再复核主值余弦和原始三维杆长；不把负判别式无条件夹到零。投影半径为零时相位不唯一，返回失败。

成功输出须满足杆长相对误差 `512 * epsilon(double)`，容差以连杆自身长度计，不随长曲柄放大。还核验选定根的不变量 `W=|r.y| sin(gamma)-r.z cos(gamma) >= 0`，切点只容许同一连杆尺度的浮点舍入带。对编码器输出解码后的角度再次执行两项检查，并要求与选定角的差不超过 `1e-12 rad`。因此相位被零偏舍入到另一个仍闭合的根，也不会报告成功。核心与可选二阶模块共用同一私有几何求解器。

工作根判定先用保守舍入界判断能否明确接受或拒绝，只有临界情况才调用长精度融合乘加。该快速路径保留相同判定阈值；普通输入可避免本机软件融合乘加的额外成本。

几何约束为 `||B - A|| = l_bar` 和 `||C - B|| = l_rod`。测试直接重建 B、C 的坐标并检查杆长，避免仅用 IK/FK 自洽来判断正确性。

## 正解与雅可比

逆解函数记为 `theta = f(q)`。默认雅可比从杆长闭环隐式求导。令 `s` 为固定曲柄侧（±1），`L` 为曲柄长度，`r = C - A`：

```text
u(gamma) = [0, s cos(gamma), sin(gamma)]
u'(gamma) = [0, -s sin(gamma), cos(gamma)]
v = r - L u
d(gamma)/dq = v dot (dC/dq) / (L r dot u')
J_motor = direction * d(gamma)/dq
```

当分母接近零时，执行器处于传动奇异位；这与两路雅可比组成的 2×2 矩阵秩亏是不同的检查。两个旋转中心不重合时：

```text
rolled = Rx(roll) (C0 - O_roll)
C_pitch = e_y cross (C - O_pitch)
C_roll = Ry(pitch) (e_x cross rolled)
```

独立的前向差分参考默认步长 `1e-6 rad`：

```text
J[:, k] = (f(q + eps e_k) - f(q)) / eps
```

数值版本要求所有差分点可达，并拒绝无法在浮点数中产生有效扰动的步长。正解求 `f(q) - theta_target = 0`，每步通过带列主元的线性求解得到 Newton 增量，并尝试最多 12 次步长回溯；只有残差下降且逆解可达时才接受新点。超出当前分支电机角区间 `alpha ± pi/2` 的目标在迭代前返回不可达。端点按与 IK 相同的可表示加法比较，避免先减零偏造成端点误拒绝；残差使用稳定范数，极小非零误差不会因平方下溢被认作零。

收敛状态包括成功、非法输入、不可达迭代/差分点、奇异雅可比、无下降步长、达到迭代上限。初值默认为零，也可传入上一帧姿态。

## 速度、力矩与增益

```text
theta_dot = J q_dot
q_dot     = solve(J, theta_dot)
tau_task  = J^T tau_motor
tau_motor = solve(J^T, tau_task)

K_motor = J^(-T) K_task J^(-1)
D_motor = J^(-T) D_task J^(-1)
```

力矩关系满足 `tau_task · q_dot = tau_motor · theta_dot`。增益接口输入 pitch、roll 两个非负对角元素，返回完整的 2×2 电机矩阵。

刚度公式是无预载、局部微小位移下的线性映射，未包含姿态相关雅可比在预载下产生的额外刚度项。库保留交叉项，调用方需要按完整矩阵使用。

`second_order.h` 中的独立接口提供电机角 Hessian 和带预载静态映射：

```text
H_i = d²(theta_i)/dq²
G = -sum_i(tau_motor_i H_i)
K_task = J^T K_motor J + G
K_motor = J^(-T) (K_task - G) J^(-1)
```

这里力矩采用恢复力约定 `tau_motor = -gradient(U_motor)`，预载必须与当前姿态、坐标方向和单位一致。对称不定结果保留并报告谱标志；静态映射不包含运动中控制器的速度、延迟等项。完整推导和验证范围见[二阶说明](research/second_order_transmission.md)。
