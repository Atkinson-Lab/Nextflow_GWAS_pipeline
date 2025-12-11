#!/usr/bin/env python3
"""
Generate comprehensive HTML report for the ancestry-aware GWAS pipeline.
"""

import argparse
import json
from pathlib import Path
from datetime import datetime

HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Ancestry-Aware GWAS Pipeline Report</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            line-height: 1.6;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
        }}
        .header h1 {{
            margin: 0;
            font-size: 2.5em;
        }}
        .section {{
            background: white;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .section h2 {{
            color: #333;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #667eea;
            color: white;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
        .stat-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }}
        .stat-card .value {{
            font-size: 2em;
            font-weight: bold;
            color: #667eea;
        }}
        .stat-card .label {{
            color: #666;
            margin-top: 5px;
        }}
        .ancestry-badge {{
            display: inline-block;
            padding: 4px 12px;
            border-radius: 20px;
            font-size: 0.9em;
            font-weight: 500;
        }}
        .EUR {{ background: #e3f2fd; color: #1565c0; }}
        .AFR {{ background: #fff3e0; color: #e65100; }}
        .EAS {{ background: #e8f5e9; color: #2e7d32; }}
        .SAS {{ background: #fce4ec; color: #c2185b; }}
        .AMR {{ background: #f3e5f5; color: #7b1fa2; }}
        .MID {{ background: #e0f7fa; color: #00838f; }}
        .footer {{
            text-align: center;
            color: #666;
            margin-top: 30px;
            padding: 20px;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Ancestry-Aware GWAS Pipeline Report</h1>
        <p>Generated: {timestamp}</p>
    </div>

    <div class="section">
        <h2>Analysis Summary</h2>
        <div class="stat-grid">
            <div class="stat-card">
                <div class="value">{n_traits}</div>
                <div class="label">Traits Analyzed</div>
            </div>
            <div class="stat-card">
                <div class="value">{n_ancestries}</div>
                <div class="label">Ancestry Groups</div>
            </div>
            <div class="stat-card">
                <div class="value">{n_significant}</div>
                <div class="label">Significant Loci</div>
            </div>
            <div class="stat-card">
                <div class="value">{n_credible_sets}</div>
                <div class="label">Credible Sets</div>
            </div>
        </div>
    </div>

    <div class="section">
        <h2>GWAS Results by Ancestry</h2>
        {gwas_table}
    </div>

    <div class="section">
        <h2>Meta-Analysis Results</h2>
        {meta_table}
    </div>

    <div class="section">
        <h2>Heritability Estimates</h2>
        {h2_table}
    </div>

    <div class="section">
        <h2>PRS Performance</h2>
        {prs_table}
    </div>

    <div class="footer">
        <p>Ancestry-Aware GWAS Pipeline v1.0.0</p>
        <p>Report generated using Nextflow</p>
    </div>
</body>
</html>
"""


def generate_report(
    gwas_results: str,
    meta_results: str,
    fm_results: str,
    h2_results: str,
    prs_results: str,
    output_path: str,
) -> None:
    """Generate HTML report."""

    # Placeholder data for demonstration
    report_data = {
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "n_traits": "N/A",
        "n_ancestries": "9",
        "n_significant": "N/A",
        "n_credible_sets": "N/A",
        "gwas_table": "<p>GWAS results will be displayed here</p>",
        "meta_table": "<p>Meta-analysis results will be displayed here</p>",
        "h2_table": "<p>Heritability estimates will be displayed here</p>",
        "prs_table": "<p>PRS performance metrics will be displayed here</p>",
    }

    # Parse results files and populate tables
    # This is a simplified version - in practice, would parse actual result files

    html = HTML_TEMPLATE.format(**report_data)

    with open(output_path, "w") as f:
        f.write(html)

    print(f"Report generated: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate pipeline report"
    )
    parser.add_argument("--gwas-results", required=True)
    parser.add_argument("--meta-results", default="NA")
    parser.add_argument("--fm-results", default="NA")
    parser.add_argument("--h2-results", default="NA")
    parser.add_argument("--prs-results", default="NA")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    generate_report(
        args.gwas_results,
        args.meta_results,
        args.fm_results,
        args.h2_results,
        args.prs_results,
        args.output,
    )


if __name__ == "__main__":
    main()
