# Deletes a virtual track

Deletes a virtual track.

## Usage

``` r
gvtrack.rm(vtrack = NULL)
```

## Arguments

- vtrack:

  virtual track name

## Value

None.

## Details

This function deletes a virtual track from current R environment.

## See also

[`gvtrack.create`](https://tanaylab.github.io/misha/reference/gvtrack.create.md),
[`gvtrack.ls`](https://tanaylab.github.io/misha/reference/gvtrack.ls.md)

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
gvtrack.rm("vtrack1")
gvtrack.ls()
#>  [1] "vtrack2"            "vtrack3"            "vtrack4"           
#>  [4] "cov"                "motif_score"        "max_motif_score"   
#>  [7] "max_motif_pos"      "cg_count"           "cg_frac"           
#> [10] "at_pos"             "at_neg"             "at_both"           
#> [13] "g_frac"             "c_frac"             "masked_count"      
#> [16] "masked_frac"        "gc"                 "value_track"       
#> [19] "value_track_max"    "spatial_pwm"        "regular_pwm"       
#> [22] "spatial_extended"   "window_pwm"         "window_spatial_pwm"
#> [25] "pwm_score"          "edist_above"        "edist_below"       
#> [28] "edist_below_all"    "edist_max2"        
```
