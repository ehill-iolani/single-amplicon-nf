#!/usr/bin/env python3
"""
Combine the per-sample outputs into the run-level report:
  consensus_summary.tsv   one row per sample (including samples that got no consensus)
  all_consensus.fasta     every sample's consensus, header carries the QC status
  run_qc_summary.html     the summary table, status-coloured

Stdlib only. Per-sample inputs are the small TSVs the pipeline steps emit; a
sample with no qc row (too few reads) is reported as FAIL / too_few_reads.
"""
import argparse
import csv
import html


def read_rows(paths):
    """{sample: row-dict} over files that each hold a header and one row."""
    rows = {}
    for p in paths:
        with open(p) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                rows[row["sample"]] = row
    return rows


def read_fasta(path):
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                header = line[1:]
            else:
                seq.append(line)
    return header, "".join(seq)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", nargs="+", required=True)
    ap.add_argument("--purity", nargs="*", default=[])
    ap.add_argument("--qc", nargs="*", default=[])
    ap.add_argument("--consensus", nargs="*", default=[])
    ap.add_argument("--min-reads", type=int, default=40)
    ap.add_argument("--min-dominant-frac", type=float, default=0.8)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-html", required=True)
    args = ap.parse_args()

    counts = read_rows(args.counts)
    purity = read_rows(args.purity)
    qc = read_rows(args.qc)
    consensus = {}
    for p in args.consensus:
        header, seq = read_fasta(p)
        if header:
            consensus[header.split()[0]] = seq

    cols = [
        "sample", "status", "reasons", "reads_after_trim", "reads_used",
        "dominant_fraction", "n_clusters", "consensus_length", "n_ambiguous",
        "mean_depth", "mapped_fraction", "median_identity",
    ]

    summary = []
    for sample in sorted(counts):
        c, q, pu = counts[sample], qc.get(sample), purity.get(sample)
        reasons = []
        if q is None:
            status = "FAIL"
            reasons.append(f"too_few_reads({c['reads_used']}<{args.min_reads})")
        else:
            status = q["qc_status"]
            if q["qc_reasons"]:
                reasons.append(q["qc_reasons"])
        if pu is not None and float(pu["dominant_fraction"]) < args.min_dominant_frac:
            if status == "PASS":
                status = "WARN"
            reasons.append(f"mixed_reads(dominant={float(pu['dominant_fraction']):.2f}<{args.min_dominant_frac:g})")
        row = {
            "sample": sample,
            "status": status,
            "reasons": ";".join(reasons),
            "reads_after_trim": c["reads_after_trim"],
            "reads_used": c["reads_used"],
            "dominant_fraction": pu["dominant_fraction"] if pu else "NA",
            "n_clusters": pu["n_clusters"] if pu else "NA",
            "consensus_length": q["consensus_length"] if q else "NA",
            "n_ambiguous": q["n_ambiguous"] if q else "NA",
            "mean_depth": q["mean_depth"] if q else "NA",
            "mapped_fraction": q["mapped_fraction"] if q else "NA",
            "median_identity": q["median_identity"] if q else "NA",
        }
        summary.append(row)

    with open(args.out_summary, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(summary)

    with open(args.out_fasta, "w") as out:
        for row in summary:
            seq = consensus.get(row["sample"])
            if seq:
                out.write(f">{row['sample']} length={len(seq)} mean_depth={row['mean_depth']} status={row['status']}\n{seq}\n")

    colour = {"PASS": "#1a7f37", "WARN": "#9a6700", "FAIL": "#cf222e"}
    with open(args.out_html, "w") as out:
        n = {s: sum(1 for r in summary if r["status"] == s) for s in colour}
        out.write("<!doctype html><meta charset=utf-8><title>Run summary</title>"
                  "<style>body{font:14px system-ui;margin:2em}table{border-collapse:collapse}"
                  "td,th{border:1px solid #d0d7de;padding:4px 8px;text-align:left}th{background:#f6f8fa}</style>")
        out.write(f"<h1>Single-amplicon run summary</h1><p>{len(summary)} samples: "
                  f"{n['PASS']} PASS, {n['WARN']} WARN, {n['FAIL']} FAIL</p><table><tr>")
        out.write("".join(f"<th>{html.escape(c)}</th>" for c in cols) + "</tr>")
        for row in summary:
            out.write("<tr>")
            for c in cols:
                style = f' style="color:{colour[row[c]]};font-weight:600"' if c == "status" else ""
                out.write(f"<td{style}>{html.escape(str(row[c]))}</td>")
            out.write("</tr>")
        out.write("</table>")


if __name__ == "__main__":
    main()
