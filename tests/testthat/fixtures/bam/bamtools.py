"""Small BAM utilities for the fixtures: re-blocking a BAM's BGZF stream and writing raw records."""
import gzip
import random
import struct
import sys
import zlib

EOF_BLOCK = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def bgzf_block(data: bytes) -> bytes:
    c = zlib.compressobj(6, zlib.DEFLATED, -15)
    cdata = c.compress(data) + c.flush()
    bsize = 18 + len(cdata) + 8 - 1
    return (b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", bsize) + cdata +
            struct.pack("<II", zlib.crc32(data) & 0xffffffff, len(data)))


def bgzf(data: bytes, rng=None, max_block=65280, empty_every=0) -> bytes:
    out = []
    i = 0
    n = 0
    while i < len(data):
        size = rng.randint(1, max_block) if rng else max_block
        out.append(bgzf_block(data[i:i + size]))
        i += size
        n += 1
        if empty_every and n % empty_every == 0:
            out.append(bgzf_block(b""))
    out.append(EOF_BLOCK)
    return b"".join(out)


def header(refs, text=""):
    t = text.encode()
    h = b"BAM\x01" + struct.pack("<i", len(t)) + t + struct.pack("<i", len(refs))
    for name, ln in refs:
        nb = name.encode() + b"\x00"
        h += struct.pack("<i", len(nb)) + nb + struct.pack("<i", ln)
    return h


OPS = {c: i for i, c in enumerate("MIDNSHP=XB")}


def record(qname, flag, refid, pos, mapq, cigar_ops, next_refid=-1, next_pos=-1, tlen=0, l_seq=None):
    """cigar_ops: list of (len, opcode)"""
    if l_seq is None:
        l_seq = sum(n for n, op in cigar_ops if op in (0, 1, 4, 7, 8))
    qn = qname.encode() + b"\x00"
    body = struct.pack("<iiBBHHHiiii", refid, pos, len(qn), mapq, 4680, len(cigar_ops), flag, l_seq, next_refid, next_pos, tlen)
    body += qn + b"".join(struct.pack("<I", (n << 4) | op) for n, op in cigar_ops)
    body += b"\x11" * ((l_seq + 1) // 2) + b"\xff" * l_seq
    return struct.pack("<i", len(body)) + body


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "reblock":  # reblock in.bam out.bam seed
        data = gzip.decompress(open(sys.argv[2], "rb").read())
        open(sys.argv[3], "wb").write(bgzf(data, random.Random(int(sys.argv[4])), max_block=8000, empty_every=5))
    elif cmd == "edge":  # edge out.bam
        refs = [("chr1", 500000), ("chr2", 300000), ("X", 200000), ("chrUn_absent", 100000)]
        recs = [
            record("nocig_mapped", 0, 0, 1000, 60, [], l_seq=50),            # mapped, no CIGAR
            record("nocig_noseq", 16, 0, 1100, 60, [], l_seq=0),
            record("bop", 16, 0, 2000, 60, [(20, 0), (5, 9), (30, 0)]),       # B operation
            record("badop", 16, 0, 2100, 60, [(20, 0), (5, 12), (30, 0)]),    # operation code 12
            record("badop_f", 0, 0, 2200, 60, [(50, 0), (3, 15)]),
            record("nopos", 0, 0, -1, 60, [(50, 0)]),                         # placed on chr1 without a position
            record("mq255", 16, 1, 3000, 255, [(50, 0)]),
            record("zero_span", 16, 0, 3100, 60, [(10, 4), (5, 1), (10, 4)]), # no reference base
            record("big", 16, 0, 4000, 60, [(100000, 3), (10, 0)], l_seq=10), # N of 100000
            record("hugeN", 16, 0, 5000, 60, [(0xfffffff, 3)] * 4097 + [(10, 0)], l_seq=10),  # span >= 2^40
            record("neg_tlen", 99, 0, 6000, 60, [(50, 0)], 0, 6100, -150),
            record("mate_unplaced", 99, 0, 6200, 60, [(50, 0)], -1, -1, 150),
            record("ok", 0, 2, 100, 60, [(50, 0)]),
        ]
        open(sys.argv[2], "wb").write(bgzf(header(refs) + b"".join(recs)))
    elif cmd == "corrupt":  # corrupt outdir: invalid files that htslib (samtools) refuses too
        out = sys.argv[2]
        refs = [("chr1", 500000)]
        ok = record("ok", 0, 0, 100, 60, [(50, 0)])
        data = header(refs) + ok + ok
        open(f"{out}/bad-truncated-record.bam", "wb").write(bgzf(data[:-7]))
        open(f"{out}/bad-cigar-length.bam", "wb").write(bgzf(header(refs) + ok + record("bad", 0, 0, 200, 60, [(40, 0)], l_seq=50)))
        open(f"{out}/bad-refid.bam", "wb").write(bgzf(header(refs) + ok + record("bad", 0, 0, 200, 60, [(50, 0)], next_refid=3)))
        # a header whose n_ref (2^31 - 1) is more references than the file holds
        hdr = bytearray(header(refs))
        struct.pack_into("<i", hdr, 8, 0x7fffffff)
        open(f"{out}/bad-n-ref.bam", "wb").write(bgzf(bytes(hdr) + ok))
