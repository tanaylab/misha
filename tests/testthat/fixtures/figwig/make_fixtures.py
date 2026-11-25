# Writes the figwig fixtures of test-gtrack.import_figwig.R next to this script, with figwig's
# own writer, and the values figwig's own reader reads from them.
# Run with: python3 -I make_fixtures.py   (figwig == 0.3.2, numpy)
#
# The format files follow figwig's tests/format_fixtures.py (MIT, Jacob Schreiber): tracks built
# so that the chunks take every kind (dense, sparse, runs, matrix) and element type (u1 u2 u4 i1
# i2 i4 f4 f8 and the "one" code), in both orders and both codecs, a stored chunk in a compressed
# file, several runs per group, dead bytes left by an append, blocks of track order, and
# chunks written from entries. The script prints the chunk kinds each file holds.
#
# <name>.values.tsv: every run of equal values other than 0 that figwig's reader reads (track,
# chrom, start, end, value; NaN as "nan"), on chromosomes chr1 (300 bp), chr2 (100 bp) and chr3
# (64 bp), or chrA and chrB for blocks.fw. zlib.values.tsv holds the values of zlib.fw, stored.fw,
# tracks.fw and tracks_zlib.fw, which are the same; appended.fw holds them and track "late", 9 at
# chr1:65. blocks.values.tsv leaves out track "ones" of blocks.fw, 1 on bases 0 to 65,999 of chrA.

import collections
import os
import struct
import zlib

import numpy
from figwig import FigWigReader, FigWigWriter

HERE = os.path.dirname(os.path.abspath(__file__))
CHROMS = {"chr1": 300, "chr2": 100, "chr3": 64}
TILE = 64
N_CELLS = 120


def path(name):
    return os.path.join(HERE, name)


NAMED = ["ones", "random", "random2", "runs", "runs2", "neg8", "neg16", "neg32", "u16", "u32",
         {"name": "wide", "dtype": "float64"}, {"name": "signal", "dtype": "float32", "assay": "ATAC"},
         "dense16", "runs_int", "empty", {"name": "dense_f32", "dtype": "float32"}, "dense_i16",
         {"name": "dense_f64", "dtype": "float64"}]
CELLS = [{"name": f"cell{c}", "cell": c} for c in range(N_CELLS)]
SITES = [{"name": f"site{c}"} for c in range(N_CELLS)]
TRACKS = NAMED + CELLS + SITES
GROUPS = [["ones"], ["random", "random2"], ["runs", "runs2"], ["neg8"], ["neg16"], ["neg32"], ["u16"], ["u32"],
          ["wide", "signal"], ["dense16"], ["runs_int"], ["empty"], ["dense_f32"], ["dense_i16"], ["dense_f64"],
          [c["name"] for c in CELLS],
          [s["name"] for s in SITES]]


def values_of(length, rng):
    v = numpy.zeros((len(TRACKS), length))
    v[0, ::9] = 1                                                       # sparse, "one"
    v[1, :] = rng.integers(0, 256, length)                              # dense uint8, random: stored
    v[2, :] = rng.integers(0, 256, length)
    v[3, 20:200] = 2.5                                                  # runs, float32
    v[4, 30:40] = 1.5
    v[5, ::5] = -rng.integers(1, 100, (length + 4) // 5)                # int8
    v[6, ::6] = -rng.integers(200, 30000, (length + 5) // 6)            # int16
    v[7, ::7] = -rng.integers(40000, 2**31 - 1, (length + 6) // 7)      # int32
    v[8, ::4] = rng.integers(256, 65536, (length + 3) // 4)             # uint16
    v[9, ::8] = rng.integers(65536, 2**32, (length + 7) // 8)           # uint32
    v[10, 10] = 2**24 + 1                                               # float64 track: f8 chunks
    v[10, 11] = 0.5
    v[10, 12] = 0.1
    v[11, ::3] = rng.normal(size=(length + 2) // 3).astype(numpy.float32)
    v[11, 7] = numpy.nan                                                # float32 with NaN
    v[12, :] = rng.integers(256, 65536, length)                         # dense uint16
    v[13, 5:60] = 3                                                     # runs, uint8
    # 14, "empty", stays empty
    v[15, :] = rng.random(length).astype(numpy.float32)                 # dense f4
    v[16, :] = rng.integers(-32768, 32768, length)                      # dense i2
    v[17, :] = rng.normal(size=length) * 1e6                            # dense f8
    # few of many cells have a value in a tile: matrix chunks
    first = len(NAMED)
    for c in range(0, N_CELLS, 30):
        v[first + c, (c * 11) % length] = c // 30 + 1
    for c in range(0, N_CELLS, 30):
        v[first + N_CELLS + c, (c * 7) % length] = 1
    return v


def write(name, level=6, zooms=2, order="positions", interleave=True):
    rng = numpy.random.default_rng(0)
    with FigWigWriter(path(name), CHROMS, TRACKS, genome="toy", attributes={"made": "misha", "n": [1, 2]},
                      tile=TILE, groups=GROUPS, level=level, zooms=zooms, order=order, n_jobs=1) as writer:
        for chrom, length in CHROMS.items():
            values = values_of(length, rng)
            if interleave and order == "positions":
                # each group in turn, two calls each, so that a group's chunks lie in several runs
                half = length // 2 // TILE * TILE or length
                for group in GROUPS:
                    ids = [writer._track_ids[n] for n in group]
                    writer.write(chrom, [0], values[None, ids, :half], tracks=group)
                    if half < length:
                        writer.write(chrom, [half], values[None, ids, half:], tracks=group)
            else:
                writer.write(chrom, [0], values[None])


def write_appended():
    with open(path("zlib.fw"), "rb") as src, open(path("appended.fw"), "wb") as dst:
        dst.write(src.read())
    with FigWigWriter(path("appended.fw"), None, ["late"], append=True) as writer:
        v = numpy.zeros((1, 1, CHROMS["chr1"]))
        v[0, 0, 65] = 9
        writer.write("chr1", [0], v)


def write_blocks():
    chroms = {"chrA": 70_000, "chrB": 10}
    with FigWigWriter(path("blocks.fw"), chroms, ["t", "u", "ones"], order="tracks", level=6, zooms=0) as writer:
        rng = numpy.random.default_rng(1)
        v = numpy.zeros((1, 3, chroms["chrA"]))
        v[0, 0, rng.choice(chroms["chrA"], 2_000, replace=False)] = rng.integers(1, 1000, 2_000)
        v[0, 1, :5] = [1, numpy.nan, 3, 0, 5]
        v[0, 2, :66_000] = 1                              # 66,000 entries: two blocks of the one code
        writer.write("chrA", [0], v)
        w = numpy.zeros((1, 3, chroms["chrB"]))
        w[0, 0, 9] = 7
        writer.write("chrB", [0], w)


def write_entries():
    tracks = ["steps", "random", {"name": "gaps", "dtype": "float32"}, "scattered"]
    with FigWigWriter(path("entries.fw"), CHROMS, tracks, tile=TILE,
                      groups=[[t] for t in ["steps", "random", "gaps", "scattered"]], level=6, zooms=1) as writer:
        rng = numpy.random.default_rng(2)
        for chrom, length in CHROMS.items():
            v = numpy.zeros((4, length))
            v[0] = numpy.repeat(rng.integers(0, 4, -(-length // 8)), 8)[:length]
            v[1] = rng.integers(1, 200, length)
            v[2] = numpy.repeat(rng.integers(0, 3, -(-length // 16)), 16)[:length]
            v[2, 20:45] = numpy.nan
            v[3, ::11] = 3
            which, positions = numpy.nonzero(v)
            writer.write_sparse(chrom, positions, which, v[which, positions])


def write_semantics():
    # attributes of every JSON kind, names to sanitize, and chromosomes: "1" for chr1 (an alias),
    # chr2 shorter than the database's, chrUn_x not in the database, no chr3
    tracks = [
        {"name": "atac.+", "assay": "ATAC", "rep": 1, "description": "from the file", "nested": {"k": [1, 2.5], "s": "é"},
         "flag": True, "none": None, "empty": "", "created.by": "someone", "binsize": 10, "type": "x", "uni": "café"},
        {"name": "atac.-", "assay": "ATAC", "rep": 2, "dtype": "float32"},
        {"name": "H3K27ac-rep1", "assay": "ChIP", "rep": "2"},
        {"name": "1x", "assay": "ChIP"},
    ]
    chroms = {"1": 300, "chr2": 90, "chrUn_x": 50}
    with FigWigWriter(path("semantics.fw"), chroms, tracks, tile=16, groups=2) as writer:
        rng = numpy.random.default_rng(3)
        for chrom, length in chroms.items():
            v = numpy.zeros((1, 4, length))
            v[0, 0, rng.choice(length, length // 4, replace=False)] = rng.integers(1, 5, length // 4)
            v[0, 1] = numpy.round(rng.random(length), 2)
            v[0, 1, 10:20] = numpy.nan
            v[0, 2, 30:70] = 2
            v[0, 3, ::7] = 1
            writer.write(chrom, [0], v)


def write_numattrs():
    # number attributes written as floats, for select matching 2, 2.0 and "2.0" alike
    tracks = [{"name": "a", "rep": 2.0, "frac": 0.1}, {"name": "b", "rep": 3, "frac": 0.25}]
    chroms = {"chr1": 100}
    with FigWigWriter(path("numattrs.fw"), chroms, tracks) as writer:
        writer.write("chr1", [0], numpy.ones((1, 2, 100), dtype=numpy.uint8))


def write_nulgenome():
    # a genome name holding a NUL byte, which an R string cannot
    with FigWigWriter(path("nulgenome.fw"), {"chr1": 100}, ["a"], genome="hg\x0038") as writer:
        writer.write("chr1", [0], numpy.ones((1, 1, 100), dtype=numpy.uint8))


def runs_of(reader):
    rows = []
    for t, rec in enumerate(reader.tracks):
        for chrom, size in reader.chrom_sizes.items():
            x = reader.read([chrom], [0], size, tracks=[t], dtype="float64")[0, 0]
            i = 0
            while i < size:
                j = i + 1
                while j < size and (x[j] == x[i] or (numpy.isnan(x[j]) and numpy.isnan(x[i]))):
                    j += 1
                if x[i] != 0:
                    rows.append(f"{rec['name']}\t{chrom}\t{i}\t{j}\t{'nan' if numpy.isnan(x[i]) else repr(float(x[i]))}")
                i = j
    return rows


def chunk_kinds(name):
    """(kind, element code, stored) of every chunk or block of a file, from its index."""
    with open(path(name), "rb") as f:
        data = f.read()
    flags, = struct.unpack_from("<H", data, 6)
    n_groups, = struct.unpack_from("<I", data, 24)
    index, = struct.unpack_from("<Q", data, 48)
    optional, = struct.unpack_from("<I", data, 100)
    kinds = collections.Counter()
    for g in range(n_groups):
        rec, = struct.unpack_from("<Q", data, index + 8 * g)
        if flags & 1:
            n_blocks, = struct.unpack_from("<I", data, rec)
            for b in range(n_blocks):
                kind, = struct.unpack_from("<I", data, rec + 8 + 40 * b + 32)
                kinds[("block", (kind >> 3) & 15, 0)] += 1
        else:
            n_runs, = struct.unpack_from("<I", data, rec)
            p = rec + 8 + 32 * n_runs
            n_live, n_bytes = struct.unpack_from("<II", data, p)
            live = zlib.decompress(data[p + 8:p + 8 + n_bytes])
            for k in live[4 * n_live:5 * n_live]:
                kinds[(["zero", "dense", "sparse", "matrix", "runs"][k & 7], (k >> 3) & 15, k >> 7)] += 1
    return kinds, optional


def main():
    write("zlib.fw")
    write("stored.fw", level=0, zooms=0, interleave=False)
    write("tracks.fw", order="tracks", level=0, zooms=1)
    write("tracks_zlib.fw", order="tracks", level=6)
    write_appended()
    write_blocks()
    write_entries()
    write_semantics()
    write_numattrs()
    write_nulgenome()

    codes = ["u1", "u2", "u4", "i1", "i2", "i4", "f4", "one", "f8"]
    format_rows = None
    for name in ["zlib.fw", "stored.fw", "tracks.fw", "tracks_zlib.fw", "appended.fw", "blocks.fw", "entries.fw",
                 "semantics.fw"]:
        rows = runs_of(FigWigReader(path(name)))
        if name == "zlib.fw":
            format_rows = rows
        elif name in ("stored.fw", "tracks.fw", "tracks_zlib.fw"):
            assert rows == format_rows, name
        elif name == "appended.fw":
            assert rows == format_rows + ["late\tchr1\t65\t66\t9.0"], name
        if name in ("zlib.fw", "blocks.fw", "entries.fw", "semantics.fw"):
            if name == "blocks.fw":
                ones = [r for r in rows if r.startswith("ones\t")]
                assert ones == ["ones\tchrA\t0\t66000\t1.0"]
                rows = [r for r in rows if not r.startswith("ones\t")]
            with open(path(name[:-3] + ".values.tsv"), "w") as f:
                f.write("track\tchrom\tstart\tend\tvalue\n")
                f.write("\n".join(rows) + "\n")
        kinds, optional = chunk_kinds(name)
        shown = ", ".join(f"{k} {codes[c]}{' stored' if s else ''}: {n}" for (k, c, s), n in sorted(kinds.items()))
        print(f"{name} ({os.path.getsize(path(name))} bytes, optional bits {optional}): {shown}")


main()
