#!/usr/bin/env python3
"""
Builds the synthetic dataset under tests/data/ (deterministic: fixed seed).

Four samples, each chosen to land on a different outcome:
  sample_a  300 reads of amplicon A                    -> PASS
  sample_b  300 reads of amplicon B (unrelated to A)   -> PASS
  sample_c  225 reads of A + 75 of B (75/25 mix)       -> WARN mixed_reads, consensus = A
  sample_d  20 reads of A                              -> FAIL too_few_reads (< --min_reads)

Reads carry ~6% ONT-like error (2% each substitution / insertion / deletion),
random orientation, and are split over two part-files to exercise MERGE_FASTQ.
"""
import gzip
import os
import random

random.seed(7)
HERE = os.path.dirname(os.path.abspath(__file__))
COMP = str.maketrans("ACGT", "TGCA")


def rand_seq(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def revcomp(s):
    return s.translate(COMP)[::-1]


def noisy(seq, sub=0.02, ins=0.02, dele=0.02):
    out = []
    for b in seq:
        r = random.random()
        if r < dele:
            continue
        out.append(random.choice([x for x in "ACGT" if x != b]) if r < dele + sub else b)
        if random.random() < ins:
            out.append(random.choice("ACGT"))
    return "".join(out)


def read(seq, i):
    s = noisy(seq)
    if random.random() < 0.5:
        s = revcomp(s)
    return f"@read{i}\n{s}\n+\n{'0' * len(s)}\n"  # Q15


def write_sample(name, reads):
    d = os.path.join(HERE, "reads", name)
    os.makedirs(d, exist_ok=True)
    half = len(reads) // 2
    for part, chunk in enumerate((reads[:half], reads[half:])):
        with gzip.open(os.path.join(d, f"part_{part}.fastq.gz"), "wt") as fh:
            fh.writelines(chunk)


A, B = rand_seq(650), rand_seq(650)

write_sample("sample_a", [read(A, i) for i in range(300)])
write_sample("sample_b", [read(B, i) for i in range(300)])
mix = [read(A, i) for i in range(225)] + [read(B, 1000 + i) for i in range(75)]
random.shuffle(mix)
write_sample("sample_c", mix)
write_sample("sample_d", [read(A, i) for i in range(20)])

with open(os.path.join(HERE, "samplesheet.csv"), "w") as fh:
    fh.write("sample,fastq\n")
    for s in ("sample_a", "sample_b", "sample_c", "sample_d"):
        fh.write(f"{s},tests/data/reads/{s}/*.fastq.gz\n")
