"""Deterministic synthetic mpileup data for the ebscore tests (no BAMs needed).

Writes into tests/data/synthetic/:
- ref.tsv           chrom, pos, ref base of the synthetic genome
- pon.mpileup       raw `samtools mpileup` text of 6 PON normals
- pair.mpileup      raw mpileup of normal (sample 1) + tumor (sample 2)
- variants.tsv      ANNOVAR-style variant table (Chr Start End Ref Alt + scenario)
- variants.vcf      the same variants as VCF (the multi-allelic sites as one record with two ALTs)

Scenarios cover: background noise, a strand-biased PON artifact, real SNVs, a strong SNV (400x tumor; with only 6 PON
normals at 20-80x it still scores EB < 16, so the sf/legacy Fisher split is tested on the large-PON fixture in
conftest.py), insertions, deletions, a position without PON coverage, a position without tumor coverage, a
multi-allelic site and an MNV. Read tokens include start markers with special MAPQ chars (^$, ^+, ^-, ^5) and end
markers.

Run: python tests/data/make_synthetic.py   (any env with numpy)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

OUT = Path(__file__).parent / "synthetic"
SEED = 20261122
N_PON = 6
BASES = "ACGT"
CHROMS = {"chr7": (1000, 260), "chr12": (5000, 140)}  # start, length


def make_tokens(rng, depth, ref, events):
    """events: list of (kind, payload, prob, fwd_share). kind in snv|ins|del."""
    tokens = []
    for _ in range(depth):
        fwd = rng.random() < 0.5
        tok = ""
        if rng.random() < 0.06:
            tok += "^" + chr(int(rng.integers(33, 75)))  # MAPQ char, may be $ + - , . or a digit
        base = "." if fwd else ","
        suffix = ""
        u = rng.random()
        acc = 0.0
        for kind, payload, prob, fwd_share in events:
            acc += prob
            if u < acc:
                # strand-biased events: re-draw the strand
                fwd = rng.random() < fwd_share
                base = "." if fwd else ","
                if kind == "snv":
                    base = payload if fwd else payload.lower()
                elif kind == "ins":
                    suffix = f"+{len(payload)}{payload if fwd else payload.lower()}"
                elif kind == "del":
                    suffix = f"-{len(payload)}{payload if fwd else payload.lower()}"
                break
        tok += base + suffix
        if rng.random() < 0.04:
            tok += "$"
        tokens.append(tok)
    if rng.random() < 0.1 and depth:
        tokens.append("*" if rng.random() < 0.5 else "#")
    return tokens


def sample_field(rng, depth, ref, events):
    if depth == 0:
        return "0\t*\t*"
    tokens = make_tokens(rng, depth, ref, events)
    quals = "".join(chr(int(rng.integers(35, 75))) for _ in tokens)
    return f"{len(tokens)}\t{''.join(tokens)}\t{quals}"


def noise(rng, ref):
    alts = [b for b in BASES if b != ref]
    return [("snv", alts[int(rng.integers(0, 3))], 0.004, 0.5), ("ins", "A", 0.001, 0.5), ("del", "C", 0.001, 0.5)]


def main():
    rng = np.random.default_rng(SEED)
    OUT.mkdir(parents=True, exist_ok=True)
    genome = {c: "".join(rng.choice(list(BASES), size=n)) for c, (s, n) in CHROMS.items()}

    # scenario per (chrom, offset): kind, details
    plan = {}
    variants = []  # Chr Start End Ref Alt scenario (+ vcf fields)

    def ref_at(chrom, pos):
        return genome[chrom][pos - CHROMS[chrom][0]]

    def add(chrom, off, kind, **kw):
        plan[(chrom, CHROMS[chrom][0] + off)] = (kind, kw)

    for chrom, (start, n) in CHROMS.items():
        offs = rng.permutation(np.arange(5, n - 5))
        it = iter(int(o) for o in offs)
        for _ in range(8):
            add(chrom, next(it), "real_snv")
        for _ in range(4):
            add(chrom, next(it), "artifact_snv")
        add(chrom, next(it), "strong_snv")
        for _ in range(3):
            add(chrom, next(it), "ins")
        for _ in range(3):
            o = next(it)
            add(chrom, o, "del")  # o is the anchor base
            plan[(chrom, start + o + 1)] = ("after_del", {})
        add(chrom, next(it), "pon_gap")
        add(chrom, next(it), "tumor_gap")
        add(chrom, next(it), "multi")
        add(chrom, next(it), "mnv")
        for _ in range(4):
            add(chrom, next(it), "noise_call")  # a "variant" that is pure noise

    pon_lines, pair_lines, ref_lines = [], [], []
    for chrom, (start, n) in CHROMS.items():
        for pos in range(start, start + n):
            ref = ref_at(chrom, pos)
            ref_lines.append(f"{chrom}\t{pos}\t{ref}")
            kind, kw = plan.get((chrom, pos), ("bg", {}))
            alts = [b for b in BASES if b != ref]
            pon_ev = noise(rng, ref)
            tum_ev = noise(rng, ref)
            pon_depths = [int(rng.integers(20, 80)) for _ in range(N_PON)]
            n_depth, t_depth = int(rng.integers(30, 90)), int(rng.integers(40, 150))
            if kind in ("real_snv", "noise_call", "multi", "mnv"):
                alt = alts[int(rng.integers(0, 3))]
                if kind == "real_snv":
                    tum_ev = [("snv", alt, 0.25, 0.5)] + tum_ev
                elif kind == "multi":
                    alt2 = [b for b in alts if b != alt][0]
                    tum_ev = [("snv", alt, 0.2, 0.5), ("snv", alt2, 0.1, 0.5)] + tum_ev
                    variants.append((chrom, pos, pos, ref, alt2, "multi"))
                if kind == "mnv":
                    nxt = ref_at(chrom, pos + 1)
                    alt_n = [b for b in BASES if b != nxt][0]
                    variants.append((chrom, pos, pos + 1, ref + nxt, alt + alt_n, "mnv"))
                else:
                    variants.append((chrom, pos, pos, ref, alt, kind))
            elif kind == "artifact_snv":
                alt = alts[int(rng.integers(0, 3))]
                pon_ev = [("snv", alt, 0.05, 0.9)] + pon_ev
                tum_ev = [("snv", alt, 0.06, 0.9)] + tum_ev
                variants.append((chrom, pos, pos, ref, alt, kind))
            elif kind == "strong_snv":
                alt = alts[int(rng.integers(0, 3))]
                t_depth = 400
                tum_ev = [("snv", alt, 0.45, 0.5)] + tum_ev
                variants.append((chrom, pos, pos, ref, alt, kind))
            elif kind == "ins":
                seq = "".join(rng.choice(list(BASES), size=int(rng.integers(1, 4))))
                tum_ev = [("ins", seq, 0.2, 0.5)] + tum_ev
                variants.append((chrom, pos, pos, "-", seq, kind))
            elif kind == "del":
                dref = ref_at(chrom, pos + 1)
                tum_ev = [("del", dref, 0.2, 0.5)] + tum_ev
                variants.append((chrom, pos + 1, pos + 1, dref, "-", kind))
            elif kind == "pon_gap":
                pon_depths = [0] * N_PON
                alt = alts[0]
                tum_ev = [("snv", alt, 0.3, 0.5)] + tum_ev
                variants.append((chrom, pos, pos, ref, alt, kind))
            elif kind == "tumor_gap":
                t_depth = 0
                variants.append((chrom, pos, pos, ref, alts[1], kind))
            if any(pon_depths):
                fields = [sample_field(rng, d, ref, pon_ev) for d in pon_depths]
                pon_lines.append("\t".join([chrom, str(pos), ref, *fields]))
            fields = [sample_field(rng, n_depth, ref, noise(rng, ref)), sample_field(rng, t_depth, ref, tum_ev)]
            pair_lines.append("\t".join([chrom, str(pos), ref, *fields]))

    (OUT / "ref.tsv").write_text("\n".join(ref_lines) + "\n")
    (OUT / "pon.mpileup").write_text("\n".join(pon_lines) + "\n")
    (OUT / "pair.mpileup").write_text("\n".join(pair_lines) + "\n")
    order = {c: i for i, c in enumerate(CHROMS)}
    variants.sort(key=lambda v: (order[v[0]], v[1], v[4]))
    with open(OUT / "variants.tsv", "w") as fh:
        fh.write("Chr\tStart\tEnd\tRef\tAlt\tscenario\n")
        for v in variants:
            fh.write("\t".join(map(str, v)) + "\n")
    with open(OUT / "variants.vcf", "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for c, (s, n) in CHROMS.items():
            fh.write(f"##contig=<ID={c},length={s + n + 100}>\n")
        fh.write('##INFO=<ID=SC,Number=1,Type=String,Description="test scenario">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        records = []  # [chrom, pos, ref, [alts], kind]; same chrom/pos/ref -> one multi-allelic record
        for chrom, s, e, ref, alt, kind in variants:
            if ref == "-":  # insertion after s
                anchor = ref_at(chrom, s)
                rec = (s, anchor, anchor + alt)
            elif alt == "-":  # deletion of s..e, anchor s-1
                anchor = ref_at(chrom, s - 1)
                rec = (s - 1, anchor + ref, anchor)
            else:
                rec = (s, ref, alt)
            if records and records[-1][:3] == [chrom, rec[0], rec[1]]:
                records[-1][3].append(rec[2])
            else:
                records.append([chrom, rec[0], rec[1], [rec[2]], kind])
        for chrom, pos, ref, alts, kind in records:
            fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{','.join(alts)}\t.\tPASS\tSC={kind}\n")
    print(f"{len(pon_lines)} PON lines, {len(pair_lines)} pair lines, {len(variants)} variants -> {OUT}")


if __name__ == "__main__":
    main()
