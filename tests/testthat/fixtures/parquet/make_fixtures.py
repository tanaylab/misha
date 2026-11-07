# Writes the Parquet fixtures of test-gtrack.import_parquet.R next to this script.
# Run with: python3 -I make_fixtures.py   (pyarrow >= 15, pandas)
#
# Common data: chr1 (1000 bp) and chr2 (600 bp) in 32 bp bins, the last bin of each short.
# a (float32) = k + 0.25 and b (float64) = chrom_index * 1000 + k, where k = start / 32.

import datetime
import decimal
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
# chr2 rows with a value column the chr1 shard does not have
extra = frame([r for r in rows if r[0] == "chr2"])
extra["c"] = extra["b"]
write("extra_column.parquet", extra)

# shard_chr2 with the first page of column a overwritten: the footer and the chrom column are
# intact, so the file opens and the error comes from inside the import
data = bytearray(open(path("shard_chr2.parquet"), "rb").read())
col_a = pq.ParquetFile(path("shard_chr2.parquet")).metadata.row_group(0).column(4)
assert col_a.path_in_schema == "a"
first_page = col_a.dictionary_page_offset or col_a.data_page_offset
data[first_page:first_page + 8] = b"\xff" * 8
open(path("corrupt_page.parquet"), "wb").write(bytes(data))

write("alias.parquet", frame([(r[0][3:],) + r[1:] for r in rows]))
write("unknown_chrom.parquet", frame([r for r in rows if r[0] == "chr1"] + [("scaffold_9", s, s + 10, 1.0, 2.0) for s in (0, 10, 20)]))

write("zstd.parquet", df, compression="zstd")
write("delta.parquet", df, use_dictionary=False, column_encoding={"start": "DELTA_BINARY_PACKED"})

# chr1 row k = 3 twice
write("overlap.parquet", frame(rows[:4] + rows[3:]))
# the last chr2 row ends at 640, past the end of chr2
write("beyond_end.parquet", frame(rows[:-1] + [rows[-1][:2] + (640,) + rows[-1][3:]]))

# Rows not aligned to 32 bp bins: chr1 in 50 bp intervals with value k + 0.5, except that
# [250, 300) is missing (a gap) and [500, 550) holds NaN (a value, not a null); chr2 has two
# short rows. So some bins are only partly covered.
unaligned_rows = [("chr1", s, s + 50, s / 50 + 0.5) for s in range(0, 1000, 50) if s != 250]
unaligned_rows = [r if r[1] != 500 else r[:3] + (float("nan"),) for r in unaligned_rows]
unaligned_rows += [("chr2", 10, 20, 7.5), ("chr2", 100, 130, -2.25)]
unaligned = pa.table({
    "chrom": pa.array([r[0] for r in unaligned_rows]),
    "start": pa.array([r[1] for r in unaligned_rows], pa.int64()),
    "end": pa.array([r[2] for r in unaligned_rows], pa.int64()),
    "value": pa.array([r[3] for r in unaligned_rows], pa.float64()),
})
assert unaligned.column("value").null_count == 0
pq.write_table(unaligned, path("unaligned.parquet"))
unaligned.to_pandas().to_csv(path("unaligned.tsv"), sep="\t", index=False, na_rep="NaN")

with_name = df.copy()
with_name["name"] = ["r%d" % i for i in range(len(df))]
write("with_name.parquet", with_name)

write("single.parquet", frame(rows, ("chrom", "start", "end", "a")))

# Physical number types whose converted / logical type changes their meaning. u32 (INT32,
# UINT_32) and u64 (INT64, UINT_64) hold values past the signed range; i8 is INT32 with INT_8;
# dec is DECIMAL(5, 2) stored as INT32 (k + 0.25 is stored as 100 * k + 25); ts (INT64,
# TIMESTAMP with no converted type) and d (INT32, DATE) are not numbers to misha.
annotated = pa.table({
    "chrom": pa.array([r[0] for r in rows]),
    "start": pa.array([r[1] for r in rows], pa.int64()),
    "end": pa.array([r[2] for r in rows], pa.int64()),
    "u32": pa.array([3_000_000_000 + r[1] // BIN * 256 for r in rows], pa.uint32()),
    "u64": pa.array([2**63 + r[1] // BIN * 2**40 for r in rows], pa.uint64()),
    "i8": pa.array([-(r[1] // BIN) for r in rows], pa.int8()),
    "dec": pa.array([decimal.Decimal(r[1] // BIN) + decimal.Decimal("0.25") for r in rows], pa.decimal128(5, 2)),
    "ts": pa.array([datetime.datetime(2020, 1, 1)] * len(rows), pa.timestamp("ns")),
    "d": pa.array([datetime.date(2020, 1, 1)] * len(rows), pa.date32()),
})
pq.write_table(annotated, path("annotated.parquet"), store_decimal_as_integer=True)

# a column whose name is not a valid track name
bad_name = df.copy()
bad_name["x * 2"] = df["b"] * 2
write("bad_name.parquet", bad_name)


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

# a file with the columns of single.parquet and no row groups at all (pyarrow's ParquetFile.read()
# gives an empty table)
with pq.ParquetWriter(path("no_rowgroups.parquet"), table(frame(rows, ("chrom", "start", "end", "a"))).schema):
    pass
assert pq.ParquetFile(path("no_rowgroups.parquet")).metadata.num_row_groups == 0

# a is entirely null in the chr2 row group (one row group per chromosome): the column chunk has no
# values, and with dictionary encoding its dictionary page has none either
allnull = df.copy()
allnull.loc[allnull["chrom"] == "chr2", "a"] = float("nan")
write("allnull_dict.parquet", allnull, row_group_size=32)
write("allnull_plain.parquet", allnull, row_group_size=32, use_dictionary=False)
for name in ("allnull_dict.parquet", "allnull_plain.parquet"):
    assert pq.ParquetFile(path(name)).metadata.num_row_groups == 2

# single.parquet with the num_values of its first column chunk (the chrom column) set to 2**40 while
# the row group has 51 rows: the reader must refuse it before decoding anything. The footer
# holds num_values as the zigzag varint 102 (0x66) after the field header 0x16; the first such
# pair that matters is the second one (the first is in the schema); the asserts below check it
data = open(path("single.parquet"), "rb").read()
n = int.from_bytes(data[-8:-4], "little")
footer = data[-8 - n:-8]
i = footer.index(b"\x16\x66", footer.index(b"\x16\x66") + 1)
footer = footer[:i + 1] + varint(2 * 2**40) + footer[i + 2:]
with open(path("huge_num_values.parquet"), "wb") as fh:
    fh.write(data[:-8 - n] + footer + len(footer).to_bytes(4, "little") + b"PAR1")
md = pq.ParquetFile(path("huge_num_values.parquet")).metadata
assert md.row_group(0).num_rows == 51 and md.row_group(0).column(0).num_values == 2**40
assert md.row_group(0).column(1).num_values == 51

# single.parquet with a null in rows 16 to 31 of a, written uncompressed and PLAIN, so that the
# definition levels of a are the RLE runs (16 x 1, 16 x 0, 19 x 1). Their values are then set to
# 2 and 0: the runs still count 32 values, which the page holds, but a level of 2 does not exist
# (pyarrow refuses the file: "Malformed levels")
single = pq.read_table(path("single.parquet"))
a = pa.array([None if 16 <= i < 32 else v for i, v in enumerate(single.column("a").to_pylist())], type=pa.float32())
pq.write_table(single.set_column(3, "a", a), path("bad_def_level.parquet"), use_dictionary=False, compression="none")
data = open(path("bad_def_level.parquet"), "rb").read()
levels = bytes.fromhex("06 00 00 00 20 01 20 00 26 01")  # their 4-byte length, then the runs
assert data.count(levels) == 1
with open(path("bad_def_level.parquet"), "wb") as fh:
    fh.write(data.replace(levels, bytes.fromhex("06 00 00 00 20 02 20 00 26 00")))

# single.parquet with 1601 rows in its footer, one more than the bases of the test genome
# (chr1 1000 + chr2 600 bp): num_rows of the file and of the row group, and num_values of every
# column chunk, all the field header 0x16 then the zigzag varint 102 (51)
data = open(path("single.parquet"), "rb").read()
n = int.from_bytes(data[-8:-4], "little")
footer = data[-8 - n:-8]
assert footer.count(b"\x16\x66") == 6
footer = footer.replace(b"\x16\x66", b"\x16" + varint(2 * 1601))
with open(path("more_rows_than_bases.parquet"), "wb") as fh:
    fh.write(data[:-8 - n] + footer + len(footer).to_bytes(4, "little") + b"PAR1")
md = pq.ParquetFile(path("more_rows_than_bases.parquet")).metadata
assert md.num_rows == 1601 and md.row_group(0).num_rows == 1601
assert all(md.row_group(0).column(k).num_values == 1601 for k in range(4))

# Many chromosomes per row group: chr1 .. chr41 of 100 + 7 * k bp (k = 1 .. 41), chr41 without rows.
# Odd k: rows that are the 32 bp bins; even k: 50 bp rows, without [100, 150). value = k + start / 1024.
# many_chroms.parquet holds them in order in row groups of 37 rows (each holds several
# chromosomes, and chromosomes span row groups); many_chroms_shuffled.parquet in random order in
# row groups of 40 rows (each holds rows of half the chromosomes or more). many_chroms.tsv holds
# the same rows for gtrack.import.
many = []
for k in range(1, 41):
    size = 100 + 7 * k
    step = 32 if k % 2 else 50
    for start in range(0, size, step):
        if step == 50 and start == 100:
            continue
        many.append(("chr%d" % k, start, min(start + step, size), k + start / 1024))
many = pd.DataFrame(many, columns=["chrom", "start", "end", "value"])
pq.write_table(table(many), path("many_chroms.parquet"), row_group_size=37)
pq.write_table(table(many.sample(frac=1, random_state=1)), path("many_chroms_shuffled.parquet"), row_group_size=40)
many.to_csv(path("many_chroms.tsv"), sep="\t", index=False)
for name, min_chroms in (("many_chroms.parquet", 3), ("many_chroms_shuffled.parquet", 20)):
    pf = pq.ParquetFile(path(name))
    assert pf.metadata.num_row_groups > 5
    assert all(len(set(pf.read_row_group(i).column("chrom").to_pylist())) >= min_chroms for i in range(pf.metadata.num_row_groups - 1))

# single.parquet with start DOUBLE and end FLOAT (arrow::write_parquet writes misha's data frames
# with double start and end); then with a start of 160.5, and with a NaN start (a value, not a null)
float_coords = table(frame(rows, ("chrom", "start", "end", "a")))
float_coords = float_coords.set_column(1, "start", float_coords.column("start").cast(pa.float64()))
float_coords = float_coords.set_column(2, "end", float_coords.column("end").cast(pa.float32()))
pq.write_table(float_coords, path("float_coords.parquet"))
for name, bad in (("float_coords_fraction.parquet", 160.5), ("float_coords_nan.parquet", float("nan"))):
    starts = float_coords.column("start").to_pylist()
    starts[5] = bad
    pq.write_table(float_coords.set_column(1, "start", pa.array(starts, pa.float64())), path(name))
assert pq.read_table(path("float_coords_nan.parquet")).column("start").null_count == 0

# single.parquet with a chromosome name of 1025 bytes in every row: one dictionary entry, which
# the reader would otherwise copy into every row
long_chrom = frame(rows, ("chrom", "start", "end", "a"))
long_chrom["chrom"] = "c" * 1025
write("long_chrom.parquet", long_chrom)
assert pq.ParquetFile(path("long_chrom.parquet")).metadata.row_group(0).column(0).dictionary_page_offset is not None

# Empty files, as shards with no rows: pyarrow writes one row group of no rows (column chunks with
# data_page_offset 0); pandas also types the chrom column of an empty frame as null.
pq.write_table(table(frame(bins())).slice(0, 0), path("empty_pyarrow.parquet"))
frame([]).to_parquet(path("empty_pandas.parquet"), index=False)
