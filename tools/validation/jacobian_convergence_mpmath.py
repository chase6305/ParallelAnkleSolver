#!/usr/bin/env python3
"""Independent high-precision geometry check; run only under the shared work lock.

Consumes <prefix>.meta.json and <prefix>.csv from jacobian_convergence_probe.
Every JSON/CSV double is reconstructed as its exact binary-float rational.
No NumPy, SciPy, YAML or production solver bindings are used.
"""

import argparse
import csv
import json
import math
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prefix", type=Path)
    parser.add_argument("--dps", type=int, default=90, choices=range(75, 101))
    args = parser.parse_args()

    # This optional dependency is imported only when the probe is executed.
    import mpmath as mp

    mp.mp.dps = args.dps
    prefix = str(args.prefix)
    metadata = json.loads(Path(prefix + ".meta.json").read_text())

    def exact(value):
        if isinstance(value, int):
            return mp.mpf(value)
        numerator, denominator = float(value).as_integer_ratio()
        return mp.mpf(numerator) / denominator

    def vector(key):
        return [exact(value) for value in metadata[key]]

    def add(first, second):
        return [a + b for a, b in zip(first, second)]

    def subtract(first, second):
        return [a - b for a, b in zip(first, second)]

    def rotate_x(point, angle):
        x, y, z = point
        return [x, mp.cos(angle) * y - mp.sin(angle) * z, mp.sin(angle) * y + mp.cos(angle) * z]

    def rotate_y(point, angle):
        x, y, z = point
        return [mp.cos(angle) * x + mp.sin(angle) * z, y, -mp.sin(angle) * x + mp.cos(angle) * z]

    pitch = exact(metadata["pitch"])
    roll = exact(metadata["roll"])
    roll_center = vector("r_roll_center")
    pitch_center = vector("r_pitch_center")
    functions = []
    baseline_errors = []
    for motor in (1, 2):
        anchor = vector(f"r_a{motor}")
        foot_zero = vector(f"r_c{motor}")
        bar = exact(metadata[f"l_bar{motor}"])
        rod = exact(metadata[f"l_rod{motor}"])
        offset = exact(metadata[f"alpha_{motor}_rad"])
        direction = 1 if metadata[f"is_left_motor{motor}"] else -1
        side = 1 if foot_zero[1] > anchor[1] else -1
        baseline_gamma = direction * (exact(metadata["theta"][motor - 1]) - offset)

        def geometry_candidates(qp, qr, anchor=anchor, foot_zero=foot_zero, bar=bar, rod=rod, side=side):
            local = rotate_x(subtract(foot_zero, roll_center), qr)
            delta = add(rotate_y(add(local, subtract(roll_center, pitch_center)), qp),
                        subtract(pitch_center, anchor))
            # General crank-circle/sphere intersection: no source asin formula.
            a, b = side * delta[1], delta[2]
            projected = mp.sqrt(a * a + b * b)
            rhs = (sum(value * value for value in delta) + bar * bar - rod * rod) / (2 * bar)
            discriminant = projected * projected - rhs * rhs
            if discriminant < 0:
                raise ValueError("high-precision geometry reports an unreachable target")
            center = mp.atan2(b, a)
            opening = mp.atan2(mp.sqrt(discriminant), rhs)
            return center + opening, center - opening

        candidates = geometry_candidates(pitch, roll)
        distances = [abs(mp.atan2(mp.sin(value - baseline_gamma), mp.cos(value - baseline_gamma)))
                     for value in candidates]
        chosen = 0 if distances[0] <= distances[1] else 1
        winding = mp.nint((baseline_gamma - candidates[chosen]) / (2 * mp.pi))
        baseline_errors.append(distances[chosen])

        def motor_angle(qp, qr, geometry=geometry_candidates, chosen=chosen, winding=winding,
                        direction=direction, offset=offset):
            return direction * (geometry(qp, qr)[chosen] + 2 * mp.pi * winding) + offset

        functions.append(motor_angle)

    reference = [[mp.diff(lambda qp: function(qp, roll), pitch),
                  mp.diff(lambda qr: function(pitch, qr), roll)]
                 for function in functions]
    norm_reference = mp.sqrt(sum(value * value for row in reference for value in row))
    print(f"mpmath dps={args.dps}; exact binary-float geometry and pose")
    print("independent geometric atan2 branch selected only by base motor angle")
    print("reference Jacobian:")
    for row in reference:
        print("  " + "  ".join(mp.nstr(value, 30) for value in row))
    print("base crank-angle differences:", ", ".join(mp.nstr(value, 12) for value in baseline_errors))

    rows = list(csv.DictReader(Path(prefix + ".csv").open(newline="")))
    output_rows = []
    print("method     eps          relative_Frobenius_error_to_mp    stencil_status")
    for row in rows:
        values = [float(row[key]) for key in ("J00", "J01", "J10", "J11")]
        if row["status"] != "ok" or not all(math.isfinite(value) for value in values):
            print(f"{row['method']:10s} {row['eps']:12s} {'n/a':32s} {row['status']}")
            output_rows.append({"method": row["method"], "eps": row["eps"], "dps": args.dps,
                                "status": row["status"], "relative_frobenius_error_to_mp": "nan"})
            continue
        error = mp.sqrt(sum((exact(values[2 * r + c]) - reference[r][c]) ** 2
                            for r in range(2) for c in range(2)))
        relative = error / (norm_reference if norm_reference else 1)
        text = mp.nstr(relative, 16)
        print(f"{row['method']:10s} {row['eps']:12s} {text:32s} {row['status']}")
        output_rows.append({"method": row["method"], "eps": row["eps"], "dps": args.dps,
                            "status": row["status"], "relative_frobenius_error_to_mp": text})
    output_path = Path(prefix + ".mp.csv")
    with output_path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=["method", "eps", "dps", "status", "relative_frobenius_error_to_mp"])
        writer.writeheader()
        writer.writerows(output_rows)
    print("Wrote", output_path)


if __name__ == "__main__":
    main()
