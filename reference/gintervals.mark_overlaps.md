# Mark overlapping intervals with a group ID

Mark overlapping intervals with a group ID

## Usage

``` r
gintervals.mark_overlaps(
  intervals,
  group_col = "overlap_group",
  unify_touching_intervals = TRUE
)
```

## Arguments

- intervals:

  intervals set

- group_col:

  name of the column to store the overlap group IDs (default:
  "overlap_group")

- unify_touching_intervals:

  if 'TRUE', touching one-dimensional intervals are unified

## Value

The intervals set with an additional column containing group IDs from
gintervals.canonic mapping. All overlapping intervals will have the same
group ID.

## Examples

``` r

gdb.init_examples()
# Create sample overlapping intervals
intervs <- data.frame(
    chrom = "chr1",
    start = c(11000, 100, 10000, 10500),
    end = c(12000, 200, 13000, 10600),
    data = c(10, 20, 30, 40)
)

# Mark overlapping intervals
intervs_marked <- gintervals.mark_overlaps(intervs)

# Use custom column name
intervs_marked <- gintervals.mark_overlaps(intervs, group_col = "my_groups")
```
