#!/usr/bin/env python3
"""
Trim a polished consensus back to the part its own reads support.

spoa (-r 0) can run the consensus out past the amplicon on a few reads' extra
sequence: chimeric reads that kept ONT adapter / native-barcode flank through
cutadapt (validation on Nishii 2026, PRJEB85706: 45 of 88 consensus carried
> 20 bp of it, up to 263 bp, from ~1% of the reads). racon and medaka polish
that tail but can't remove it. Those bases are covered by a handful of reads
while the rest of the amplicon is covered by nearly all of them, so each end
is cut inward to the first position whose depth reaches --frac x the median
depth. Only the ends are cut, never the middle, and nothing is cut when no
position reaches the threshold (no alignments: CONSENSUS_QC fails it anyway).

--frac is a share of the median, not --min-depth: reads that span the whole
amplicon give a flat depth plateau, so the cut doesn't depend on how deep the
sample was sequenced. It is kept low (0.2) because the last few real bases
dip a little below the plateau (read ends soft-clip), and the tail sits at ~1%.

Depth is counted from the primary and supplementary alignments in the SAM
(minimap2 -a of the consensus's reads against it), like CONSENSUS_QC counts
the PAF's tp:P segments. Stdlib only.

Writes
  <prefix>.consensus.fasta   the trimmed consensus, header kept + ends_trimmed=<5'>,<3'>
  <prefix>.end_trim.tsv      one row: lengths, depths, bases cut and the cut sequence
"""
import argparse
import gzip
import re
import statistics

CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")
REF_OPS = set("MDN=X")


def read_fasta(path):
    """(header without '>', sequence) of the first record."""
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    break
                header = line[1:]
            else:
                seq.append(line.strip())
    return header, "".join(seq)


def alignment_spans(path):
    """(start, end) on the reference, 0-based half-open, of every primary and
    supplementary alignment. Unmapped (4) and secondary (256) records are skipped."""
    opener = gzip.open if path.endswith(".gz") else open
    spans = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.split("\t", 6)
            if len(f) < 6:
                continue
            flag = int(f[1])
            if flag & 4 or flag & 256 or f[5] == "*":
                continue
            start = int(f[3]) - 1
            span = sum(int(n) for n, op in CIGAR_OP.findall(f[5]) if op in REF_OPS)
            spans.append((start, start + span))
    return spans


def depth_array(length, spans):
    diff = [0] * (length + 1)
    for start, end in spans:
        diff[max(0, min(length, start))] += 1
        diff[max(0, min(length, end))] -= 1
    depth, run = [], 0
    for i in range(length):
        run += diff[i]
        depth.append(run)
    return depth


def trim_bounds(depth, frac):
    """(start, end, median, threshold): the kept slice [start:end] runs from the
    first to the last position with depth >= frac x median. Untrimmed when no
    position qualifies."""
    if not depth:
        return 0, 0, 0.0, 0.0
    median = statistics.median(depth)
    threshold = frac * median
    ok = [i for i, d in enumerate(depth) if d >= threshold and d > 0]
    if not ok:
        return 0, len(depth), median, threshold
    return ok[0], ok[-1] + 1, median, threshold


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--consensus", required=True)
    ap.add_argument("--sam", required=True)
    ap.add_argument("--frac", type=float, default=0.2, help="share of the median depth an end base needs to be kept")
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    header, seq = read_fasta(args.consensus)
    depth = depth_array(len(seq), alignment_spans(args.sam)) if seq else []
    start, end, median, threshold = trim_bounds(depth, args.frac)
    if not seq:
        start = end = 0
    cut5, cut3 = start, len(seq) - end
    kept = seq[start:end]

    with open(f"{args.out_prefix}.consensus.fasta", "w") as out:
        if header is not None:
            # one line, like racon writes it
            out.write(f">{header} ends_trimmed={cut5},{cut3}\n{kept}\n")

    cols = ["sample", "untrimmed_length", "consensus_length", "trimmed_start", "trimmed_end",
            "median_depth", "trim_depth", "trimmed_start_seq", "trimmed_end_seq"]
    row = [args.sample, len(seq), len(kept), cut5, cut3, f"{median:g}", f"{threshold:.1f}",
           seq[:start], seq[end:]]
    with open(f"{args.out_prefix}.end_trim.tsv", "w") as out:
        out.write("\t".join(cols) + "\n")
        out.write("\t".join(str(v) for v in row) + "\n")


if __name__ == "__main__":
    main()
