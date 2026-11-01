# Writes the Parquet fixtures of test-gtrack.import_parquet.R next to this script.
# Run with: python3 -I make_fixtures.py   (pyarrow >= 15, pandas)
#
# Common data: chr1 (1000 bp) and chr2 (600 bp) in 32 bp bins, the last bin of each short.
# a (float32) = k + 0.25 and b (float64) = chrom_index * 1000 + k, where k = start / 32.

import os

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

HERE = os.path.dirname(os.path.abspath(__file__))
CHROMS = [("chr1", 1000), ("chr2", 600)]
BIN = 32


def path(name):
    return os.path.join(HERE, name)


def bins(chroms=CHROMS):
    rows = []
    for ci, (chrom, size) in enumerate(chroms):
        for start in range(0, size, BIN):
            k = start // BIN
            rows.append((chrom, start, min(start + BIN, size), k + 0.25, ci * 1000.0 + k))
    return rows


def frame(rows, cols=("chrom", "start", "end", "a", "b")):
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "a", "b"])
    df["start"] = df["start"].astype("int64")
    df["end"] = df["end"].astype("int64")
    df["a"] = df["a"].astype("float32")
    df["b"] = df["b"].astype("float64")
    return df[list(cols)]


def table(df):
    return pa.Table.from_pandas(df, preserve_index=False)


def write(name, df, **kw):
    pq.write_table(table(df), path(name), **kw)


rows = bins()
df = frame(rows)

df.to_parquet(path("pandas_default.parquet"), index=False)
write("plain_uncompressed.parquet", df, use_dictionary=False, compression="none")
write("gzip_v2.parquet", df, compression="gzip", data_page_version="2.0")

# a is null where k % 7 == 0
nulls = table(df)
a = pa.array([None if r[1] // BIN % 7 == 0 else r[3] for r in rows], type=pa.float32())
nulls = nulls.set_column(nulls.schema.get_field_index("a"), "a", a)
pq.write_table(nulls, path("nulls.parquet"))

# row groups of 20 rows (the middle one holds both chr1 and chr2), several pages per chunk
write("mixed_rowgroups.parquet", df, row_group_size=20, data_page_size=64, write_batch_size=4)
# the dictionary outgrows its 64 byte limit: dictionary pages first, then PLAIN pages
write("dict_fallback.parquet", df, dictionary_pagesize_limit=64, write_batch_size=8)
write("byte_stream_split.parquet", df, use_dictionary=False, use_byte_stream_split=["start", "end", "a", "b"])

# per-chromosome shards, the value columns in a different order
write("shard_chr1.parquet", frame([r for r in rows if r[0] == "chr1"]))
write("shard_chr2.parquet", frame([r for r in rows if r[0] == "chr2"], ("chrom", "start", "end", "b", "a")))

write("alias.parquet", frame([(r[0][3:],) + r[1:] for r in rows]))
write("unknown_chrom.parquet", frame([r for r in rows if r[0] == "chr1"] + [("scaffold_9", s, s + 10, 1.0, 2.0) for s in (0, 10, 20)]))

write("zstd.parquet", df, compression="zstd")
write("delta.parquet", df, use_dictionary=False, column_encoding={"start": "DELTA_BINARY_PACKED"})

# chr1 row k = 3 twice
write("overlap.parquet", frame(rows[:4] + rows[3:]))
# the last chr2 row ends at 640, past the end of chr2
write("beyond_end.parquet", frame(rows[:-1] + [rows[-1][:2] + (640,) + rows[-1][3:]]))

# chr1 in 50 bp intervals, value k + 0.5: not aligned to 32 bp bins
unaligned = pd.DataFrame({
    "chrom": "chr1",
    "start": pd.Series(range(0, 1000, 50), dtype="int64"),
    "end": pd.Series(range(50, 1050, 50), dtype="int64"),
    "value": [k + 0.5 for k in range(20)],
})
pq.write_table(table(unaligned), path("unaligned.parquet"))
unaligned.to_csv(path("unaligned.tsv"), sep="\t", index=False)

with_name = df.copy()
with_name["name"] = ["r%d" % i for i in range(len(df))]
write("with_name.parquet", with_name)

write("single.parquet", frame(rows, ("chrom", "start", "end", "a")))


# plain_uncompressed.parquet with fields of every Thrift compact type that parquet.thrift does
# not have (ids 100 and up, as a newer writer might add) at the start of FileMetaData, before the
# fields a reader needs: it must skip them exactly
def varint(v):
    out = b""
    while v >= 0x80:
        out += bytes([v & 0x7F | 0x80])
        v >>= 7
    return out + bytes([v])


def field(fid, ctype):  # the long form of a field header: type, then the zigzag id
    return bytes([ctype]) + varint(fid << 1)


def binary(s):
    return varint(len(s)) + s


nested = (  # a struct holding a list of two structs and a false bool
    field(1, 9) + bytes([0x2C]) + field(1, 5) + varint(4) + b"\x00" + b"\x00"
    + field(2, 2) + b"\x00"
)
unknown = b"".join([
    field(100, 1),                                          # bool true (no value bytes)
    field(101, 3) + b"\x7f",                                # i8
    field(102, 4) + varint(300),                            # i16
    field(103, 7) + b"\x00" * 7 + b"\x40",                  # double
    field(104, 8) + binary(b"a binary value"),              # binary
    field(105, 9) + bytes([0x35]) + varint(2) * 3,          # list<i32>, short form
    field(106, 9) + bytes([0xF6]) + varint(20) + varint(6) * 20,  # list<i64>, long form
    field(107, 10) + bytes([0x21]) + b"\x01\x02",           # set<bool>
    field(108, 11) + varint(2) + bytes([0x8C]) + (binary(b"k") + nested) * 2,  # map<binary, struct>
    field(109, 11) + b"\x00",                               # empty map
    field(110, 12) + nested,                                # struct
])
data = open(path("plain_uncompressed.parquet"), "rb").read()
n = int.from_bytes(data[-8:-4], "little")
footer = data[-8 - n:-8]
# FileMetaData starts with version (field 1, i32) in the short form, whose id is relative to the
# field before it: write that header in the long form after the unknown fields
assert footer[0] == 0x15
footer = unknown + field(1, 5) + footer[1:]
with open(path("unknown_fields.parquet"), "wb") as fh:
    fh.write(data[:-8 - n] + footer + len(footer).to_bytes(4, "little") + b"PAR1")
assert pq.read_table(path("unknown_fields.parquet")).equals(pq.read_table(path("plain_uncompressed.parquet")))
