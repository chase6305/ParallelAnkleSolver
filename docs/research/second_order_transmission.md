# Second-order transmission derivatives and loaded static stiffness

**Optional API, added 2026-11-12.** `parallel_ankle/second_order.h` exposes requested-only Hessians and a loaded static stiffness map for `PMAnkle`. The existing `GetMotorStiffnessMatrixFromTaskStiffness()` interface retains its documented unloaded, local congruence mapping. This independently derived extension does not change that interface or establish hardware validation.

## Convention and static scope

Let task coordinates be `q`, motor coordinates be `theta = f(q)`, and define:

```text
J   = d(theta)/dq
H_i = d²(theta_i)/dq²
```

Each `H_i` is a symmetric task-space matrix on a smooth, fixed assembly branch. For a conservative motor potential `U_m(theta)`, use the **restoring torque** and positive stiffness convention:

```text
tau_motor = -gradient_theta U_m
K_motor  = Hessian_theta U_m = -d(tau_motor)/d(theta)
U_task(q) = U_m(f(q))
tau_task = -gradient_q U_task = J^T tau_motor
K_task   = Hessian_q U_task = -d(tau_task)/dq
```

Differentiating the torque relation, including the variation of `J`, gives:

```text
d(tau_task)/dq = -J^T K_motor J + sum_i(tau_motor_i H_i)

K_task = J^T K_motor J - sum_i(tau_motor_i H_i)
```

For a square, invertible `J`, the corresponding reverse mapping is:

```text
K_motor = J^(-T) [K_task + sum_i(tau_motor_i H_i)] J^(-1)
```

The torque in the geometric correction is the motor restoring preload **at the same operating point**, in the same signed motor-coordinate convention as `J` and `H_i`. A desired task preload can first be mapped through `J^T tau_motor = tau_task`. The matrices can be indefinite under some preloads; this formula does not guarantee that a positive-gain motor controller can realize an arbitrary requested task stiffness.

This is a static conservative tangent-stiffness result. It does not include moving-state feedback terms, friction, backlash, compliance of the links, or delay. The differential-actuation paper derives additional preload and velocity-dependent contributions when transferring a moving PD controller; its full result must not be replaced by this static specialization. [Primary paper, section II-D, equations 25–30](https://arxiv.org/html/2503.22459v2).

## Exact scalar energy check

Choose a nonlinear transmission and a quadratic motor potential:

```text
theta = q²
U_m(theta) = (k/2) (theta - 2)²
```

At `q=1`, the quantities are:

```text
theta=1, J=2, H=2, tau_motor=k, K_motor=k
```

Directly composing the potential gives:

```text
U_task(q) = (k/2) (q² - 2)²
dU_task/dq = 2k q(q² - 2)
d²U_task/dq² = 6k q² - 4k
K_task(q=1) = 2k
```

The transmission formula agrees: `J² K_motor - tau_motor H = 4k - 2k = 2k`. An incorrect plus sign in the task-space relation would produce `6k`. Conversely, the loaded inverse map returns `(2k + 2k)/4 = k`; omitting preload would return `k/2`.

## Implicit Hessian of a rotary crank and rigid rod

Assume a fixed motor center `A`, a fixed crank length `L`, a fixed rod length `R`, and a unit circular crank direction `u(gamma)` perpendicular to its fixed rotation axis. Let the foot attachment position be `C(q)`:

```text
r = C(q) - A
v = r - L u(gamma)
F(q,gamma) = v dot v - R² = 0
u'  = du/dgamma
u'' = -u
D = L (r dot u')
```

For each task coordinate `q_j`, write `C_j = dC/dq_j`, `C_jk = d²C/(dq_j dq_k)`, and `gamma_j = dgamma/dq_j`. The first derivative is:

```text
gamma_j = (v dot C_j) / D
```

The required partial derivatives, holding `gamma` fixed when differentiating with respect to `q`, are:

```text
F_gamma       = -2D
F_jk          = 2(C_j dot C_k + v dot C_jk)
F_j,gamma     = -2L(C_j dot u')
F_gamma,gamma = 2L(r dot u)
```

Differentiating `F(q,gamma(q))=0` again therefore gives:

```text
gamma_jk = [C_j dot C_k + v dot C_jk
            - L((C_j dot u') gamma_k + (C_k dot u') gamma_j)
            + L(r dot u) gamma_j gamma_k] / D
```

For the motor convention `theta_i = direction_i * gamma_i + offset_i`, with a constant sign `direction_i` in `{+1,-1}`:

```text
J_i   = direction_i * gamma_gradient_i
H_i   = direction_i * gamma_Hessian_i
```

When reusing a motor-space Jacobian row, recover `gamma_j = direction_i * J_i[j]` before evaluating the crank-angle Hessian. The symmetry of the expression is explicit in `j,k`.

The transmission denominator `D` must be nonzero. Hessians generally amplify errors more strongly than first derivatives near this singularity. Branch feasibility, angle representation, closure residuals, and finite outputs must all remain checked. Numerical normalization can reuse the current relative-coordinate scaling, treating the scale as a fixed numerical preconditioner within one derivative evaluation.

For a different parameterized actuator, the general scalar implicit rule remains:

```text
gamma_jk = -(F_jk + F_j,gamma gamma_k + F_k,gamma gamma_j
             + F_gamma,gamma gamma_j gamma_k) / F_gamma
```

The simplified circular-crank expression above does not cover a moving motor center, a configuration-dependent motor axis, an axial crank component, variable rod lengths, or arbitrary knee/ankle coupling without adding their derivatives.

## Rotary ankle with separated roll and pitch centers

Let `O_r`, `O_p` denote the zero-pose roll and pitch centers. The current relative-coordinate transform can be written:

```text
rolled = Rx(roll) (C0 - O_r)
pitch_relative = Ry(pitch) [rolled + (O_r - O_p)]
C = pitch_relative + O_p
```

With fixed unit axes `e_x`, `e_y`, its first derivatives are:

```text
C_p = e_y cross pitch_relative
C_r = Ry(pitch) [e_x cross rolled]
```

The second derivatives are:

```text
C_pp = e_y cross C_p
C_pr = C_rp = e_y cross C_r
C_rr = Ry(pitch) [e_x cross (e_x cross rolled)]
```

These depend only on relative coordinates, so a common representable translation of every geometric point leaves the angles, Jacobian, and Hessians unchanged. Forming a large world-space `C` before subtracting the centers would unnecessarily lose precision.

## Optional API and validation

Link `ParallelAnkle::core` and include `parallel_ankle/second_order.h`:

```cpp
const auto H = parallel_ankle::ComputeMotorAngleHessians(ankle, pitch, roll);
const auto loaded = parallel_ankle::MapLoadedStaticTaskStiffness(
    ankle, pitch, roll, K_task, tau_motor_at_pose);
// loaded.motor_stiffness: full symmetric 2x2 matrix
// loaded.geometric_stiffness = -sum_i(tau_motor_at_pose[i] * H[i])
```

`H[i]` differentiates motor `i` with respect to `[pitch,roll]`, including that
motor's sign. Constant encoder offsets do not change its derivatives. The
implementation uses the same ideal geometric angle as the core Jacobian,
rather than recovering it from rounded encoder output; this matters near a
transmission fold. Ordinary IK/FK/Jacobian calls do not compute Hessians.

The task matrix must be finite and symmetric. Symmetric indefinite matrices
are preserved; relative spectral flags report positive semidefiniteness of the
task and motor matrices. They do not certify passivity or closed-loop stability.
Scaling is used only for symmetry and spectral diagnostics. Existing symmetric
entries, including representable tiny diagonal or off-diagonal values, are not
reconstructed from a globally normalized matrix. Symmetrization keeps the
diagonal and uses a stable midpoint for the off-diagonal pair, avoiding both
large-value overflow and premature subnormal halving.
The inverse congruence is accumulated in `long double` until both factors have
been applied, then checked against the finite `double` output range before
narrowing. This prevents a small entry from disappearing in an intermediate
product even when the final entry is representable. The preload correction
retains its checked `double` arithmetic contract.
Nonfinite inputs, unsupported branches, singular inverse maps and arithmetic
overflow fail explicitly. A rank-deficient assembled Jacobian can still have
well-defined individual angle Hessians, but the inverse stiffness map rejects it.

The tests compare Hessians with central differences of the analytic Jacobian
at multiple steps, check mixed symmetry, both motor signs, separated centers,
uniform scaling and large representable common translations. An independent
second difference of a composed scalar motor potential checks the preload
sign and off-diagonal task curvature. Other cases cover zero-preload reduction,
indefinite matrices at normal and tiny scales, invalid inputs, transmission
singularities, assembled rank loss, and offset invariance near a fold.
Identity-Jacobian fixtures additionally check mixed-scale and subnormal matrix
entries, transpose-invariant averaging, and neighboring huge off-diagonal values.
A diagonally scaled Jacobian and its reversed scaling independently check that
the two congruence factors preserve a representable subnormal coupling; an
actually out-of-range result must still fail explicitly.

Run `loaded_stiffness_demo` for a full matrix example. No performance measurement
for this optional calculation is implied by the first-order Jacobian benchmark.
This API is specific to `PMAnkle`; the general formula above can guide another
model's own second derivatives but does not provide them automatically.
