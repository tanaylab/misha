"""Writes the source SAM files of the BAM fixtures."""
import random
import sys
from pathlib import Path

out = Path(sys.argv[1])  # writes single.src.sam, paired.src.sam and empty.src.sam
out.mkdir(parents=True, exist_ok=True)
rng = random.Random(20261124)

HEADER = [
    "@HD\tVN:1.6\tSO:unsorted",
    "@SQ\tSN:chr1\tLN:500000",
    "@SQ\tSN:chr2\tLN:300000",
    "@SQ\tSN:X\tLN:200000",          # alias of chrX
    "@SQ\tSN:chrUn_absent\tLN:100000",  # not in the database
]


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def rec(name, flag, rname, pos, mapq, cigar, rnext="*", pnext=0, tlen=0, s=None, qual="*"):
    if s is None:
        # query length of the CIGAR
        n, q = "", 0
        for ch in cigar:
            if ch.isdigit():
                n += ch
            else:
                if ch in "MIS=X":
                    q += int(n)
                n = ""
        s = "A" * q if q else "*"
    return "\t".join(map(str, [name, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, s, qual]))


single = list(HEADER) + [
    rec("f1", 0, "chr1", 1001, 60, "50M", s=seq(50), qual="I" * 50),
    rec("r1", 16, "chr1", 2001, 60, "50M", s=seq(50), qual="I" * 50),
    rec("fdup1", 0, "chr1", 3001, 60, "50M"),
    rec("fdup2", 0, "chr1", 3001, 30, "50M"),
    rec("fdup3", 1024, "chr1", 3001, 60, "50M"),       # marked duplicate: still imported
    rec("rdup1", 16, "chr1", 4001, 60, "50M"),          # 5' end 4050
    rec("rdup2", 16, "chr1", 4011, 60, "40M"),          # 5' end 4050
    rec("rdup3", 16, "chr1", 4001, 60, "40M10S"),       # soft clip at 5' end: 4040
    rec("sc1", 0, "chr1", 5001, 60, "5S45M"),
    rec("sc2", 16, "chr1", 5101, 60, "45M5S"),
    rec("hc1", 16, "chr1", 5201, 60, "10H40M"),
    rec("ins1", 16, "chr1", 6001, 60, "20M5I25M"),
    rec("del1", 16, "chr1", 6101, 60, "20M5D30M"),
    rec("spl1", 16, "chr1", 6201, 60, "20M1000N30M"),
    rec("eqx1", 16, "chr1", 6301, 60, "10=1X39="),
    rec("pad1", 16, "chr1", 6401, 60, "20M2P30M"),
    rec("mix1", 16, "chr1", 6501, 60, "3S10M2I5M3D10M100N10M4S2H"),
    rec("unm1", 4, "chr1", 7001, 0, "*"),                # unmapped, placed
    rec("unm2", 4, "*", 0, 0, "*", s=seq(30)),            # unmapped, unplaced
    rec("unm3", 20, "chr2", 7001, 0, "50M"),              # unmapped with a CIGAR
    rec("sec1", 256, "chr1", 8001, 60, "50M"),
    rec("qcf1", 512, "chr1", 8101, 60, "50M"),
    rec("sup1", 2048, "chr1", 8201, 60, "30M20S"),
    rec("sup2", 2064, "chr1", 8301, 60, "20S30M"),
    rec("pf1", 65, "chr1", 8401, 60, "50M"),              # paired flags in a single-end import
    rec("pf2", 129, "chr1", 8501, 60, "50M"),
    rec("pf3", 99, "chr1", 8601, 60, "50M", "=", 8701, 150),
    rec("pf4", 147, "chr1", 8701, 60, "50M", "=", 8601, -150),
    rec("mq0", 0, "chr1", 9001, 0, "50M"),
    rec("mq5", 16, "chr1", 9101, 5, "50M"),
    rec("mq29", 0, "chr1", 9201, 29, "50M"),
    rec("mq30", 16, "chr1", 9301, 30, "50M"),
    rec("mq255", 0, "chr1", 9401, 255, "50M"),
    rec("nocig", 0, "chr1", 9501, 60, "*", s=seq(50)),    # mapped without a CIGAR: unmapped
    rec("allclip", 0, "chr1", 9601, 60, "50S"),           # no aligned base: unmapped
    rec("noseq", 16, "chr1", 9701, 60, "50M", s="*"),
    rec("long1", 16, "chr1", 9801, 60, "150M", s=seq(150), qual="F" * 150),
    rec("end1", 0, "chr1", 500000, 60, "1M"),             # last base, forward
    rec("end2", 16, "chr1", 500000, 60, "1M"),            # last base, reverse
    rec("end3", 0, "chr1", 499990, 60, "50M"),            # forward, runs past the end
    rec("end4", 16, "chr1", 499990, 60, "50M"),           # reverse, 5' end past the end: unmapped
    rec("end5", 16, "chr1", 499951, 60, "50M"),           # reverse, 5' end on the last base
    rec("past1", 0, "chr1", 500001, 60, "50M"),           # POS past the end: unmapped
    rec("first", 0, "chr1", 1, 60, "50M"),
    rec("firstr", 16, "chr1", 1, 60, "50M"),
    rec("c2f", 0, "chr2", 101, 60, "50M"),
    rec("c2r", 16, "chr2", 201, 60, "50M"),
    rec("xf", 0, "X", 101, 60, "50M"),
    rec("xr", 16, "X", 201, 60, "50M"),
    rec("ab1", 0, "chrUn_absent", 101, 60, "50M"),       # chromosome not in the database
    rec("ab2", 16, "chrUn_absent", 201, 60, "50M"),
    rec("ab3", 256, "chrUn_absent", 301, 60, "50M"),     # filtered flag on an absent chromosome
    rec("same1", 0, "chr2", 1001, 60, "50M"),             # forward and reverse at one point
    rec("same2", 16, "chr2", 952, 60, "50M"),             # 5' end 1001 - 1 = 1000
    # 70000 CIGAR operations: BAM keeps such a CIGAR in a CG tag, with a placeholder in the record
    rec("cg1", 16, "chr1", 100001, 60, "1M1D" * 35000),
]
# bulk: enough records to span several BGZF blocks, with records across block boundaries
for i in range(2500):
    rname, ln = rng.choice([("chr1", 500000), ("chr2", 300000), ("X", 200000), ("chrUn_absent", 100000)])
    flag = rng.choice([0, 16, 0, 16, 1024, 1040, 4, 256, 2048])
    L = rng.choice([36, 50, 75])
    pos = rng.randint(1, ln - 10)
    cig = rng.choice([f"{L}M", f"5S{L-5}M", f"{L-10}M10S", f"10M{rng.randint(1,500)}N{L-10}M", f"10M2I{L-12}M", f"10M3D{L-10}M"])
    single.append(rec(f"b{i}", flag, rname, pos, rng.choice([0, 10, 30, 60]), cig if not (flag & 4) else "*"))

(out / "single.src.sam").write_text("\n".join(single) + "\n")


def pair(name, rname, pos1, pos2, L1=50, L2=50, mapq1=60, mapq2=60, first_rev=False, extra=0, tlen=None, rnext="=",
         cig1=None, cig2=None):
    """A proper pair: first mate at pos1, second at pos2 (1-based)."""
    left, right = min(pos1, pos2), max(pos1 + L1, pos2 + L2)
    t = right - left if tlen is None else tlen
    f1 = 1 | 2 | 64 | extra | (16 if first_rev else 32)
    f2 = 1 | 2 | 128 | extra | (32 if first_rev else 16)
    s1 = t if pos1 <= pos2 else -t
    return [rec(name, f1, rname, pos1, mapq1, cig1 or f"{L1}M", rnext, pos2, s1),
            rec(name, f2, rname, pos2, mapq2, cig2 or f"{L2}M", rnext, pos1, -s1)]


paired = list(HEADER)
paired += pair("p1", "chr1", 101, 251)
paired += pair("p2", "chr1", 101, 251)                     # duplicate fragment
paired += pair("p3", "chr1", 1101, 1001, first_rev=True)   # first mate on the reverse strand
paired += pair("p4", "chr1", 2001, 2101, mapq1=5)          # first mate MAPQ 5
paired += pair("p5", "chr1", 2201, 2301, mapq2=5)          # second mate MAPQ 5
paired += pair("p6", "chr1", 3001, 3001)                   # both mates at one position
paired += pair("p7", "chr1", 4001, 9001)                   # TLEN 5050 > default max.fraglen
paired += pair("p8", "chr1", 10001, 10001 + 1900)          # TLEN 1950
paired += pair("p9", "chr1", 11001, 11101, extra=1024)     # marked duplicate
paired += pair("p10", "chr1", 12001, 12101, extra=512)     # QC fail
paired += pair("p11", "chr1", 13001, 13101, extra=256)     # secondary
paired += pair("p12", "chr1", 14001, 14101, extra=2048)    # supplementary
paired += pair("p13", "chr1", 15001, 15051, cig1="5S45M", cig2="45M5S")
paired += pair("p14", "chr1", 16001, 16101, tlen=0)        # TLEN 0
paired += pair("p15", "chr1", 499901, 499951)              # fragment to the end of chr1
paired += pair("p16", "chr1", 499951, 499981, L2=50)       # fragment past the end of chr1
paired += pair("p17", "chr2", 101, 201)
paired += pair("p18", "X", 101, 201)
paired += pair("p19", "chrUn_absent", 101, 201)
paired += pair("p20", "chr1", 17001, 17101, first_rev=True, extra=1024)
paired += [
    rec("np1", 65, "chr1", 18001, 60, "50M", "=", 18101, 150),     # not proper
    rec("np1", 129, "chr1", 18101, 60, "50M", "=", 18001, -150),
    rec("mu1", 73, "chr1", 19001, 60, "50M", "=", 19001, 0),       # mate unmapped
    rec("mu1", 133, "chr1", 19001, 0, "*", "=", 19001, 0),
    rec("mu2", 67 | 8, "chr1", 19101, 60, "50M", "=", 19101, 0),   # proper flag but mate unmapped
    rec("mu2", 133 | 2, "chr1", 19101, 0, "*", "=", 19101, 0),
    rec("un1", 77, "*", 0, 0, "*", "*", 0, 0, seq(50)),           # both unmapped
    rec("un1", 141, "*", 0, 0, "*", "*", 0, 0, seq(50)),
    rec("xc1", 97, "chr1", 20001, 60, "50M", "chr2", 20001, 0),   # mates on two chromosomes
    rec("xc1", 145, "chr2", 20001, 60, "50M", "chr1", 20001, 0),
    rec("pn0", 99, "chr1", 21001, 60, "50M", "*", 0, 150),        # proper, PNEXT 0
    rec("se1", 0, "chr1", 22001, 60, "50M"),                      # single-end record
    rec("se2", 16, "chr1", 22101, 60, "50M"),
    rec("ncig", 99, "chr1", 23001, 60, "*", "=", 23101, 150, seq(50)),  # first mate without CIGAR
    rec("ncig", 147, "chr1", 23101, 60, "50M", "=", 23001, -150),
]
for i in range(1000):
    rname, ln = rng.choice([("chr1", 500000), ("chr2", 300000), ("X", 200000), ("chrUn_absent", 100000)])
    p1 = rng.randint(1, ln - 3000)
    p2 = p1 + rng.randint(0, 2500)
    paired += pair(f"bp{i}", rname, p1, p2, first_rev=rng.random() < 0.5, mapq1=rng.choice([0, 10, 30, 60]),
                   extra=rng.choice([0, 0, 0, 1024, 256]))
(out / "paired.src.sam").write_text("\n".join(paired) + "\n")

(out / "empty.src.sam").write_text("\n".join(HEADER) + "\n")
