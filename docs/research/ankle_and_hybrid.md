# Parallel ankle and hybrid-leg methods: verified sources

Checked on **2026-11-12**. This review inspected primary paper pages, author/institution pages, repository trees, selected source files, and license files. It did not clone or build these projects, import their implementations, or validate their runtime results. A public paper, a public repository, and an explicitly licensed implementation are different evidence levels.

## Selection overview

| Resource | Topology or abstraction | Evidence level | Proposed role |
| --- | --- | --- | --- |
| [Optimal ankle design, 2509.16469](https://arxiv.org/abs/2509.16469) | SPU / RSU, central universal joint | Open paper; dedicated code not verified | Independent geometric models and feasible parameterization |
| [Nonlinear coupled leg, 2610.16604v2](https://arxiv.org/html/2610.16604v2) | Double-four-bar knee, 2RSU+1U ankle, knee/ankle transmission | Open paper; dedicated code not verified | Independent coupled transmission model |
| [Differential actuation, 2503.22459v2](https://arxiv.org/html/2503.22459v2) | Four-bar knee and two-axis ankle | Open paper; release promised, dedicated repository not verified | Higher derivatives and loaded impedance transfer |
| [HyRoDyn](https://github.com/dfki-ric/hyrodyn) | Configurable serial/parallel submechanisms and implicit constraints | Verified licensed code | Optional broad-topology backend |
| [Toolbox Parallel Robots](https://github.com/Gepetto/toolbox-parallel-robots) | Pinocchio models plus bilateral constraints | Verified licensed code | Independent numerical cross-check or optional adapter |
| [Wearable 3-RRR ankle](https://github.com/YufeiZhang0601/wearable-3rrr-ankle-rehab) | Three-axis rehabilitation ankle, scalar rod constraints | Verified licensed code | Different-topology example and parameter extraction reference |

## Papers: independently implementable methods

### A Framework for Optimal Ankle Design of Humanoid Robots

The paper derives geometric inverse kinematics for SPU and RSU ankles and distinguishes physical choices among mathematical roots. RSU geometry includes motor-axis yaw, attachment locations, crank lengths, and rod lengths; the operational-region parameterization constrains feasible designs. It also evaluates configurations through Jacobian-based metrics. The implementation described in the paper uses Grasshopper, Python, and the Tunny optimizer. Its Tunny GitHub citation is an optimization dependency, not a released ankle solver. No dedicated author code repository was verified in this review. Reuse the published mathematical ideas through an independent implementation, retaining explicit branch choices and feasibility checks; do not label a newly written implementation as the authors' released code. [Primary paper and equations](https://arxiv.org/html/2509.16469v1).

### Design and Analysis of a Nonlinear Parallel Mechanism Coupling for Dynamic Humanoid Legs

The latest verified revision is **v2, 2026-10-23**. It combines a nonlinear double-four-bar knee, a 2RSU+1U ankle, and a belt transmission crossing the knee. The combined independent coordinates have three components: knee, ankle pitch, and ankle roll. A complete transmission therefore needs a coupled 3×3 derivative rather than a standalone ankle's 2×2 Jacobian. The paper derives knee geometry analytically and describes Newton-Raphson resolution of the ankle's position constraints. No dedicated source repository is linked or verified here. It is a reference for an independently implemented coupled model; HyRoDyn's related authorship and citations do not prove that its repository contains this particular hardware implementation. [Primary paper, sections III–IV](https://arxiv.org/html/2610.16604v2).

### Control of Humanoid Robots with Parallel Mechanisms using Differential Actuation Models

Revision **v2, 2025-10-08** develops compact nonlinear transmission maps for four-bar knees and two-axis ankles, with first and second derivatives. It discusses impedance transfer and numerical recovery of serial coordinates from motor encoders using temporal continuity. Separating geometry and transmission from the serial robot model makes the approach applicable to multiple ankle mechanisms. The paper says that GPU implementations will be released as open source; this statement is not a verified release or a software license. No dedicated implementation repository was established in this review. Its analytical derivatives and preload-dependent impedance terms are candidates for independent implementation and physically checked regression tests. [Primary paper, method and implementation](https://arxiv.org/html/2503.22459v2).

## Licensed implementations

### HyRoDyn: broad C++ submechanism library

Verified commit: [`d2fc248eb0270b1251a6233208c86df9b5776a8e`](https://github.com/dfki-ric/hyrodyn/tree/d2fc248eb0270b1251a6233208c86df9b5776a8e), dated 2026-10-08. The root [LICENSE](https://github.com/dfki-ric/hyrodyn/blob/d2fc248eb0270b1251a6233208c86df9b5776a8e/LICENSE) is BSD-3-Clause.

The repository contains actual 2SPRR+1U and 2SPU+1U modules, analytical closure Jacobians and derivatives, numerical forward geometry, robot URDF/YAML examples, and a general numerical-constraint module. Its README explicitly associates the software with the 2019 analytical workbench and the 2022 hybrid numerical/analytical paper. The inspected SPRR forward solver starts at zero, directly inverts its Jacobian, and aborts at twenty iterations; it is not an all-assembly-mode solver. Its ankle convention uses `[roll,pitch]` and `Rx*Ry`, so an adapter must convert coordinates and conventions deliberately. Swapping two scalar angles alone does not equate `Rx*Ry` and `Ry*Rx`; first compare the actual foot rotation and attachment-point positions under both models. This makes it a useful optional backend, while our API should retain explicit statuses, warm starts, branch validation, and failure isolation. [Source and README](https://github.com/dfki-ric/hyrodyn/tree/d2fc248eb0270b1251a6233208c86df9b5776a8e).

Specific evidence:

- [SPRR geometry and Jacobian implementation](https://github.com/dfki-ric/hyrodyn/blob/d2fc248eb0270b1251a6233208c86df9b5776a8e/src/SPRR2U1.cpp), [SPU interface](https://github.com/dfki-ric/hyrodyn/blob/d2fc248eb0270b1251a6233208c86df9b5776a8e/src/SPU2U1.hpp).
- [General implicit-constraint interface](https://github.com/dfki-ric/hyrodyn/blob/d2fc248eb0270b1251a6233208c86df9b5776a8e/src/NumericalLoopConstraints.hpp), [assembly dispatch](https://github.com/dfki-ric/hyrodyn/blob/d2fc248eb0270b1251a6233208c86df9b5776a8e/src/SubmechanismsAssembly.cpp).
- [Leg submechanism configuration](https://github.com/dfki-ric/hyrodyn/blob/d2fc248eb0270b1251a6233208c86df9b5776a8e/robot/leg/urdf/submechanisms.yml) actually selects `rrPr` knee and `2SPRR+1U` ankle modules.
- [Official institutional description](https://robotik.dfki-bremen.de/en/research/softwaretools/hyrodyn) lists additional supported types and numerical fallback. This is capability documentation, not runtime verification.
- [R3US2_ANKLE.hpp](https://github.com/dfki-ric/hyrodyn/blob/d2fc248eb0270b1251a6233208c86df9b5776a8e/src/R3US2_ANKLE.hpp) retains a zlib-license file header while the root license is BSD-3-Clause. Preserve file-level notices and resolve this discrepancy before copying that module; the review imports no files.

### Toolbox Parallel Robots: generic constraint-space reference

Verified commit: [`d7b093e92775061e82c548de198289d1aaf9cb90`](https://github.com/Gepetto/toolbox-parallel-robots/tree/d7b093e92775061e82c548de198289d1aaf9cb90), dated 2026-10-27. The [LICENSE](https://github.com/Gepetto/toolbox-parallel-robots/blob/d7b093e92775061e82c548de198289d1aaf9cb90/LICENSE) is BSD-3-Clause.

This library treats closed chains as an open Pinocchio model augmented with bilateral constraints and an actuation partition. Public modules implement closure projection, forward/inverse kinematics, inverse dynamics, and closed-loop Jacobians. The inspected Jacobian code partitions motor and free coordinates and uses a pseudoinverse. One forward solver uses CasADi/Ipopt, a preceding configuration, constraints, and position limits. These are local numerical methods; source inspection establishes no complete assembly-mode enumeration guarantee. Use it for a numerical backend or cross-validation rather than claiming that it is a drop-in analytical ankle model. This review did not establish it as the released implementation of arXiv:2503.22459, despite overlapping authors. [Jacobian source](https://github.com/Gepetto/toolbox-parallel-robots/blob/d7b093e92775061e82c548de198289d1aaf9cb90/toolbox_parallel_robots/jacobian.py), [forward solver](https://github.com/Gepetto/toolbox-parallel-robots/blob/d7b093e92775061e82c548de198289d1aaf9cb90/toolbox_parallel_robots/forward_kinematics.py).

### Wearable 3-RRR ankle: a different-topology example

Verified commit: [`fe42bf73a295e17db841e6325066b2f79c371d1c`](https://github.com/YufeiZhang0601/wearable-3rrr-ankle-rehab/tree/fe42bf73a295e17db841e6325066b2f79c371d1c), dated 2026-05-12. The [LICENSE](https://github.com/YufeiZhang0601/wearable-3rrr-ankle-rehab/blob/fe42bf73a295e17db841e6325066b2f79c371d1c/LICENSE) is MIT.

The repository contains a rehabilitation mechanism, URDF, parameter files, a report, and Python inverse kinematics. The selected IK implementation uses arbitrary rotary axes and scalar rod-length constraints, samples motor bounds, bisects detected sign changes, and selects a root near the previous solution. This is useful for parameter conventions and continuity-aware root selection. It does not establish analytical FK or Jacobians, and sign-change bracketing can miss tangent or closely spaced roots unless a sample happens to satisfy tolerance. It is a different topology and application from a dynamic biped's two-axis ankle; reuse should preserve that distinction rather than treating its demonstration model as universal leg support. [Inspected IK source](https://github.com/YufeiZhang0601/wearable-3rrr-ankle-rehab/blob/fe42bf73a295e17db841e6325066b2f79c371d1c/rh1/parallel_ik.py).

## Public code without a verified reuse license

The full repository trees and GitHub metadata were checked for these snapshots. Neither had a detected license or a root license file. They are useful discovery references, but this review does not classify them as licensed reusable implementations or import their code.

| Repository snapshot | Visible contents | Evidence limit |
| --- | --- | --- |
| [RRGGZZ/parallel-ankle-sim, `afaa892…`](https://github.com/RRGGZZ/parallel-ankle-sim/tree/afaa892bc4cbfd65fbbe615fce1deefb8d5fabf1) | PNDbotics analytical IK/FK/Jacobian module, MuJoCo/Pinocchio comparison scripts | No verified license; no author-paper association established |
| [Wangxinhui-bot/parallel_ankle_kinematics, `ad6cf6f…`](https://github.com/Wangxinhui-bot/parallel_ankle_kinematics/tree/ad6cf6fc9c2cd0fe3dda779e467d52d5e12a6539) | A single C++ ankle demo | No verified license or linked publication |

## Recommended architecture for this repository

The following is our engineering inference from the reviewed material, not a claim that these capabilities already exist here.

1. Keep each mechanism geometry in a versioned configuration. Introduce separate geometric families for rotary-crank/rod and linear-actuator ankles, with explicit axes, attachment frames, actuator sign/offset, travel bounds, and assembly branch. A label such as `2RSS+U` alone is insufficient to infer frames or physical constraints.
2. Independently implement the elementary scalar closure equations and their implicit derivatives in the Eigen core. Enumerate valid elementary roots before branch/limit selection, then validate the original vector or rod constraints. Root selection should use an explicit assembly mode and optional previous state.
3. Add a generic equality-constraint interface for an unfamiliar topology: residual `Phi(x,y)`, partial Jacobians, and a local solve for dependent coordinates `x` at prescribed independent coordinates `y`. In a regular square partition, the derivative follows from `Phi_x * dx/dy = -Phi_y`; report rank failure and residual inconsistency explicitly. A pseudoinverse must not silently convert an incompatible mechanism into a successful one.
4. Represent knee/ankle coupling as a separate composition layer with a full transmission Jacobian. Preserve cross terms in velocity and force mappings, then extend derivatives and loaded impedance only with dedicated validation.
5. Offer HyRoDyn as an optional adapter after a small, pinned integration experiment; keep RBDL/URDF parsing out of the lightweight core. Treat process-aborting paths, convention conversions, file-level notices, and local-solver limits as concrete adapter work. Toolbox Parallel Robots is an alternative independent numerical reference.

No reviewed method provides a verified guarantee that arbitrary leg topology, every assembly mode, every singular configuration, and all collision constraints are solved globally. Each backend needs its own supported-topology declaration and physical regression fixtures.
