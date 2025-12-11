#!/usr/bin/env python3
"""
Validate and standardize the input samplesheet for the ancestry-aware GWAS pipeline.
"""

import argparse
import csv
import sys
from pathlib import Path


def validate_samplesheet(samplesheet_path: str, output_path: str) -> None:
    """
    Validate samplesheet and write standardized output.

    Expected columns:
    - sample_id: Unique sample identifier
    - cohort: Optional cohort name
    - bed/bim/fam: PLINK1 files OR
    - pgen/pvar/psam: PLINK2 files OR
    - bgen/sample/bgi: BGEN files OR
    - vcf/vcf_index: VCF files
    """
    required_columns = ["sample_id"]
    valid_format_groups = [
        {"bed", "bim", "fam"},
        {"pgen", "pvar", "psam"},
        {"bgen", "sample"},
        {"vcf"},
    ]

    errors = []
    validated_rows = []

    with open(samplesheet_path, "r") as f:
        reader = csv.DictReader(f)
        headers = set(reader.fieldnames or [])

        # Check required columns
        for col in required_columns:
            if col not in headers:
                errors.append(f"Missing required column: {col}")

        if errors:
            print("\n".join(errors), file=sys.stderr)
            sys.exit(1)

        for i, row in enumerate(reader, start=2):
            row_errors = []

            # Check sample_id
            if not row.get("sample_id", "").strip():
                row_errors.append(f"Row {i}: sample_id is empty")
                continue

            # Determine format and validate files exist
            format_found = False
            for format_group in valid_format_groups:
                has_all = all(
                    row.get(col, "").strip() for col in format_group
                )
                if has_all:
                    format_found = True
                    # Validate file paths
                    for col in format_group:
                        file_path = row.get(col, "").strip()
                        if file_path and not Path(file_path).exists():
                            row_errors.append(
                                f"Row {i}: File not found: {file_path}"
                            )
                    break

            if not format_found:
                row_errors.append(
                    f"Row {i}: No valid genotype format found. "
                    "Need bed/bim/fam, pgen/pvar/psam, bgen/sample, or vcf"
                )

            if row_errors:
                errors.extend(row_errors)
            else:
                validated_rows.append(row)

    if errors:
        print("Samplesheet validation errors:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        sys.exit(1)

    # Write validated samplesheet
    if validated_rows:
        with open(output_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=validated_rows[0].keys())
            writer.writeheader()
            writer.writerows(validated_rows)

        print(f"Validated {len(validated_rows)} samples")
    else:
        print("No valid samples found", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Validate ancestry-aware GWAS samplesheet"
    )
    parser.add_argument("samplesheet", help="Input samplesheet CSV")
    parser.add_argument("output", help="Output validated samplesheet CSV")
    args = parser.parse_args()

    validate_samplesheet(args.samplesheet, args.output)


if __name__ == "__main__":
    main()
