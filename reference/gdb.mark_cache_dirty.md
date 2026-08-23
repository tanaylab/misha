# Mark cached track list as dirty

When tracks or interval sets are modified outside of misha (e.g. files
copied manually), the cached inventory may become out of date. Calling
this helper marks the cache as dirty so the next
[`gsetroot()`](https://tanaylab.github.io/misha/reference/gdb.init.md)
forces a rescan.

## Usage

``` r
gdb.mark_cache_dirty()
```

## Value

Invisible `TRUE` if the dirty flag was written, `FALSE` otherwise.

## See also

[`gdb.reload`](https://tanaylab.github.io/misha/reference/gdb.reload.md),
[`gsetroot`](https://tanaylab.github.io/misha/reference/gdb.init.md)
