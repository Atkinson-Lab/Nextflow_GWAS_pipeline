#!/usr/bin/env python3
"""
Extract credible sets from fine-mapping results.
"""

import argparse
import pandas as pd
import numpy as np


def extract_credible_sets(
    input_path: str, output_path: str, pip_threshold: float = 0.95
) -> None:
    """
    Extract credible sets from fine-mapping output.

    A credible set contains variants whose cumulative PIP exceeds the threshold,
    representing the set of variants that likely contains the causal variant.
    """

    df = pd.read_csv(input_path, sep="\t")

    # Detect PIP column
    pip_col = None
    for col in df.columns:
        if "PIP" in col.upper() or "POSTERIOR" in col.upper():
            pip_col = col
            break

    if pip_col is None:
        raise ValueError("Could not detect PIP column")

    # Detect credible set column if present
    cs_col = None
    for col in df.columns:
        if "CS" in col.upper() or "CREDIBLE" in col.upper():
            cs_col = col
            break

    # If CS column exists, use it
    if cs_col is not None:
        df_cs = df[df[cs_col].notna() & (df[cs_col] != -1)].copy()
    else:
        # Otherwise, create credible sets based on PIP threshold
        df = df.sort_values(pip_col, ascending=False)
        df["cumsum_pip"] = df[pip_col].cumsum()

        # Find variants until cumulative PIP exceeds threshold
        df_cs = df[df["cumsum_pip"] <= pip_threshold + df[pip_col]].copy()
        df_cs["CS"] = 1  # Single credible set

    # Add additional info
    df_cs["in_credible_set"] = True

    # Output columns
    output_cols = ["SNP", "CHR", "POS", pip_col, "in_credible_set"]
    if cs_col:
        output_cols.append(cs_col)

    # Filter to existing columns
    output_cols = [c for c in output_cols if c in df_cs.columns]

    df_cs[output_cols].to_csv(output_path, sep="\t", index=False)

    print(f"Credible set size: {len(df_cs)}")
    print(f"Mean PIP in credible set: {df_cs[pip_col].mean():.3f}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract credible sets from fine-mapping"
    )
    parser.add_argument("--input", required=True, help="Fine-mapping results")
    parser.add_argument("--output", required=True, help="Credible set output")
    parser.add_argument(
        "--pip-threshold",
        type=float,
        default=0.95,
        help="PIP threshold for credible sets",
    )
    args = parser.parse_args()

    extract_credible_sets(args.input, args.output, args.pip_threshold)


if __name__ == "__main__":
    main()
