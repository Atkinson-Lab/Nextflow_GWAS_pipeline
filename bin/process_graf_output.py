#!/usr/bin/env python3
"""
Process GRAF-ANC ancestry inference output to standardized format.
"""

import argparse
import pandas as pd


def process_graf(input_path: str, output_path: str, pcs_path: str) -> None:
    """
    Process GRAF output to standardized ancestry calls.

    GRAF-ANC outputs 8 major ancestry groups plus heterogeneous:
    - EUR: European
    - AFR: African
    - EAS: East Asian
    - SAS: South Asian
    - AMR: Native American
    - MID: Middle Eastern
    - AAC: African American (admixed)
    - AHI: American Hispanic (admixed)
    - HET: Heterogeneous/Other
    """

    # Read GRAF output
    df = pd.read_csv(input_path, sep="\t")

    # Standardize column names
    # GRAF typically outputs: Subject, GD1, GD2, Ancestry, Prob_EUR, Prob_AFR, etc.
    column_rename = {
        "Subject": "SAMPLE_ID",
        "Self-reported ancestry": "SELF_REPORTED",
    }

    # Find ancestry probability columns
    prob_cols = [c for c in df.columns if c.startswith("Prob_") or c.startswith("P_")]

    # Determine ancestry assignment
    if "Ancestry" in df.columns:
        df["ANCESTRY"] = df["Ancestry"]
    else:
        # Assign based on highest probability
        prob_df = df[prob_cols]
        df["ANCESTRY"] = prob_df.idxmax(axis=1).str.replace("Prob_", "").str.replace("P_", "")

    # Create output dataframe
    output_df = pd.DataFrame(
        {
            "SAMPLE_ID": df.iloc[:, 0],  # First column is usually sample ID
            "ANCESTRY": df["ANCESTRY"],
        }
    )

    # Add probability columns if present
    for col in prob_cols:
        ancestry = col.replace("Prob_", "").replace("P_", "")
        output_df[f"PROB_{ancestry}"] = df[col]

    # Write ancestry calls
    output_df.to_csv(output_path, sep="\t", index=False)

    # Write PCs if available
    pc_cols = [c for c in df.columns if c.startswith("PC") or c.startswith("GD")]
    if pc_cols:
        pcs_df = pd.DataFrame({"SAMPLE_ID": df.iloc[:, 0]})
        for i, col in enumerate(pc_cols, 1):
            pcs_df[f"PC{i}"] = df[col]
        pcs_df.to_csv(pcs_path, sep="\t", index=False)

    # Summary statistics
    print("Ancestry distribution:")
    print(output_df["ANCESTRY"].value_counts())


def main():
    parser = argparse.ArgumentParser(
        description="Process GRAF ancestry output"
    )
    parser.add_argument("--input", required=True, help="GRAF output file")
    parser.add_argument("--output", required=True, help="Standardized ancestry calls")
    parser.add_argument("--pcs", required=True, help="Principal components output")
    args = parser.parse_args()

    process_graf(args.input, args.output, args.pcs)


if __name__ == "__main__":
    main()
