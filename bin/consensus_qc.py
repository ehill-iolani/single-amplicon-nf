#!/usr/bin/env python3
"""
QC for one sample's consensus, from its reads mapped back to it (minimap2 PAF).

Stdlib only, so it also runs stand-alone. Writes
  <prefix>.qc.tsv         one row: length, ambiguous bases, depth, identity, status
  <prefix>.coverage.tsv   depth binned along the consensus (for plotting)

Status:
  FAIL  mean depth < --min-depth, or fewer than --min-primary-frac of the reads
        map to the consensus, or there is no consensus sequence
  WARN  part of the consensus is covered below --min-depth, or the consensus
        contains ambiguous bases
  PASS  otherwise
"""
import argparse
import statistics

COVERAGE_BINS = 200
# share of positions that must reach --min-depth before the sample counts as covered
MIN_BREADTH = 0.95


def read_fasta_seq(path):
    """Concatenated sequence of the first record (a consensus is one record)."""
    seq, seen = [], False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seen:
                    break
                seen = True
            else:
                seq.append(line.strip())
    return "".join(seq)


def read_paf(path):
    """Primary alignments only, as (tstart, tend, matches, block_len, read_id).
    minimap2 marks the primary alignment of each segment tp:A:P (tp:A:I for an
    inversion); a read split into several P segments is a chimera/supplementary
    signal and is counted separately."""
    segments = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            tags = {t[:2]: t[5:] for t in f[12:] if len(t) > 5}
            if tags.get("tp") not in ("P", "I"):
                continue
            segments.append((int(f[7]), int(f[8]), int(f[9]), int(f[10]), f[0]))
    return segments


def depth_array(length, segments):
    diff = [0] * (length + 1)
    for tstart, tend, _, _, _ in segments:
        diff[max(0, tstart)] += 1
        diff[min(length, tend)] -= 1
    depth, run = [], 0
    for i in range(length):
        run += diff[i]
        depth.append(run)
    return depth


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--consensus", required=True)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--n-reads", type=int, required=True, help="reads that were mapped (denominator for mapped_fraction)")
    ap.add_argument("--min-depth", type=float, default=30)
    ap.add_argument("--min-primary-frac", type=float, default=0.7)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    seq = read_fasta_seq(args.consensus).upper()
    length = len(seq)
    n_ambiguous = sum(1 for b in seq if b not in "ACGT")

    segments = read_paf(args.paf) if length else []
    reads_mapped = len({s[4] for s in segments})
    segs_per_read = {}
    for s in segments:
        segs_per_read[s[4]] = segs_per_read.get(s[4], 0) + 1
    reads_split = sum(1 for n in segs_per_read.values() if n > 1)

    depth = depth_array(length, segments) if length else []
    mean_depth = sum(depth) / length if length else 0.0
    min_depth_seen = min(depth) if depth else 0
    breadth = sum(1 for d in depth if d >= args.min_depth) / length if length else 0.0
    identities = [m / b for _, _, m, b, _ in segments if b > 0]
    median_identity = statistics.median(identities) if identities else 0.0
    mapped_fraction = reads_mapped / args.n_reads if args.n_reads else 0.0

    reasons = []
    if not length:
        status = "FAIL"
        reasons.append("no_consensus")
    else:
        if mean_depth < args.min_depth:
            reasons.append(f"low_depth({mean_depth:.0f}x<{args.min_depth:g}x)")
        if mapped_fraction < args.min_primary_frac:
            reasons.append(f"low_mapped_fraction({mapped_fraction:.2f}<{args.min_primary_frac:g})")
        status = "FAIL" if reasons else "PASS"
        if status == "PASS":
            if breadth < MIN_BREADTH:
                status = "WARN"
                reasons.append(f"uneven_coverage({breadth:.2f} of positions >= {args.min_depth:g}x)")
            if n_ambiguous:
                status = "WARN"
                reasons.append(f"ambiguous_bases({n_ambiguous})")

    cols = [
        "sample", "consensus_length", "n_ambiguous", "reads_mapped", "mapped_fraction",
        "reads_split", "mean_depth", "min_depth", "breadth_at_min_depth",
        "median_identity", "qc_status", "qc_reasons",
    ]
    row = [
        args.sample, length, n_ambiguous, reads_mapped, f"{mapped_fraction:.4f}",
        reads_split, f"{mean_depth:.1f}", min_depth_seen, f"{breadth:.4f}",
        f"{median_identity:.4f}", status, ";".join(reasons),
    ]
    with open(f"{args.out_prefix}.qc.tsv", "w") as out:
        out.write("\t".join(cols) + "\n")
        out.write("\t".join(str(v) for v in row) + "\n")

    with open(f"{args.out_prefix}.coverage.tsv", "w") as out:
        out.write("sample\tstart\tend\tmean_depth\n")
        if length:
            n_bins = min(COVERAGE_BINS, length)
            for i in range(n_bins):
                start, end = i * length // n_bins, (i + 1) * length // n_bins
                out.write(f"{args.sample}\t{start}\t{end}\t{sum(depth[start:end]) / (end - start):.1f}\n")


if __name__ == "__main__":
    main()
