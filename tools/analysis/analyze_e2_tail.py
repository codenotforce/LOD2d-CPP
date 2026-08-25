#!/usr/bin/env python3
"""Audit moving-PALOD exact-error tail rates without crossing method refreshes."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def fit(points: list[tuple[int, float]], count: int) -> dict[str, float | int] | None:
    selected = points[-count:]
    if len(selected) < count:
        return None
    x = [math.log(dof) for dof, _ in selected]
    y = [math.log(error) for _, error in selected]
    x_mean = sum(x) / count
    y_mean = sum(y) / count
    ss_x = sum((value - x_mean) ** 2 for value in x)
    ss_y = sum((value - y_mean) ** 2 for value in y)
    covariance = sum(
        (x_value - x_mean) * (y_value - y_mean)
        for x_value, y_value in zip(x, y)
    )
    if ss_x <= 0.0 or ss_y <= 0.0:
        return None
    return {
        "count": count,
        "exponent": -covariance / ss_x,
        "pearson_r": covariance / math.sqrt(ss_x * ss_y),
        "first_dof": selected[0][0],
        "last_dof": selected[-1][0],
        "first_error": selected[0][1],
        "last_error": selected[-1][1],
    }


def locate_iterations(run_dir: Path) -> Path:
    direct = run_dir / "iterations.csv"
    if direct.is_file():
        return direct
    matches = list(run_dir.glob("**/iterations.csv"))
    if len(matches) != 1:
        raise SystemExit(
            f"expected one iterations.csv below {run_dir}, found {len(matches)}"
        )
    return matches[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--gate", choices=("none", "medium", "main"), default="none")
    arguments = parser.parse_args()

    iterations = locate_iterations(arguments.run_dir)
    by_dof: dict[int, float] = {}
    with iterations.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row.get("action") != "SolveAndEstimate":
                continue
            try:
                dof = int(row["N_H"])
                error = float(row["relative_exact_energy"])
            except (KeyError, TypeError, ValueError):
                continue
            if dof > 0 and error > 0.0 and math.isfinite(error):
                by_dof[dof] = error
    points = sorted(by_dof.items())
    if len(points) < 4:
        raise SystemExit("fewer than four finite exact-error observations")

    fits = {str(count): fit(points, count) for count in (4, 6, 8, 12)}
    late_ratios = [
        points[index][1] / points[index - 1][1]
        for index in range(max(1, len(points) - 7), len(points))
    ]
    report = {
        "iterations": str(iterations),
        "observations": len(points),
        "final_dof": points[-1][0],
        "final_exact_relative_energy": points[-1][1],
        "maximum_last8_error_ratio": max(late_ratios, default=1.0),
        "fits": fits,
        "gate": arguments.gate,
        "accepted": True,
        "reasons": [],
    }

    reasons: list[str] = []
    if arguments.gate == "medium":
        fit6, fit8 = fits["6"], fits["8"]
        if len(points) < 12 or fit6 is None or fit8 is None:
            reasons.append("medium gate needs at least twelve observations")
        else:
            if fit6["exponent"] <= 0.0 or fit8["exponent"] <= 0.0:
                reasons.append("medium tail does not decay")
            if fit6["pearson_r"] > -0.98:
                reasons.append("final-six log-log correlation is weaker than -0.98")
        if report["maximum_last8_error_ratio"] > 1.02:
            reasons.append("a late exact-error increase exceeds two percent")
    elif arguments.gate == "main":
        fit8, fit12 = fits["8"], fits["12"]
        if len(points) < 16 or fit8 is None or fit12 is None:
            reasons.append("main gate needs at least sixteen observations")
        else:
            if fit8["exponent"] < 0.45:
                reasons.append("final-eight exponent is below 0.45")
            if abs(fit8["exponent"] - fit12["exponent"]) > 0.15:
                reasons.append("final-eight and final-twelve exponents differ by more than 0.15")
            if fit8["pearson_r"] > -0.98:
                reasons.append("final-eight log-log correlation is weaker than -0.98")
        if report["maximum_last8_error_ratio"] > 1.02:
            reasons.append("a late exact-error increase exceeds two percent")

    report["reasons"] = reasons
    report["accepted"] = not reasons
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    print(text, end="")
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(text, encoding="utf-8")
    return 0 if not reasons else 2


if __name__ == "__main__":
    raise SystemExit(main())
