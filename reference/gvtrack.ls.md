# Returns a list of virtual track names

Returns a list of virtual track names.

## Usage

``` r
gvtrack.ls(
  pattern = "",
  ignore.case = FALSE,
  perl = FALSE,
  fixed = FALSE,
  useBytes = FALSE
)
```

## Arguments

- pattern, ignore.case, perl, fixed, useBytes:

  see 'grep'

## Value

An array that contains the names of virtual tracks.

## Details

This function returns a list of virtual tracks that exist in current R
environment that match the pattern (see 'grep'). If called without any
arguments all virtual tracks are returned.

## See also

[`grep`](https://rdrr.io/r/base/grep.html),
[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md),
[`gvtrack.rm`](https://tanaylab.github.io/misha/reference/gvtrack.rm.md)

## Examples

``` r

gdb.init_examples()
gvtrack.create("vtrack1", "dense_track", "max")
gvtrack.create("vtrack2", "dense_track", "quantile", 0.5)
gvtrack.ls()
#>  [1] "vtrack1"            "vtrack2"            "vtrack3"           
#>  [4] "vtrack4"            "cov"                "motif_score"       
#>  [7] "max_motif_score"    "max_motif_pos"      "cg_count"          
#> [10] "cg_frac"            "at_pos"             "at_neg"            
#> [13] "at_both"            "g_frac"             "c_frac"            
#> [16] "masked_count"       "masked_frac"        "gc"                
#> [19] "value_track"        "value_track_max"    "spatial_pwm"       
#> [22] "regular_pwm"        "spatial_extended"   "window_pwm"        
#> [25] "window_spatial_pwm" "pwm_score"          "edist_above"       
#> [28] "edist_below"        "edist_below_all"    "edist_max2"        
gvtrack.ls(pattern = "*2")
#> [1] "vtrack2"    "edist_max2"
```
