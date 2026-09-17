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
#>  [7] "max_motif_score"    "max_motif_pos"      "potts_max"         
#> [10] "potts_count"        "potts_p"            "potts_m"           
#> [13] "cg_count"           "cg_frac"            "at_pos"            
#> [16] "at_neg"             "at_both"            "g_frac"            
#> [19] "c_frac"             "masked_count"       "masked_frac"       
#> [22] "gc"                 "value_track"        "value_track_max"   
#> [25] "spatial_pwm"        "regular_pwm"        "spatial_extended"  
#> [28] "window_pwm"         "window_spatial_pwm" "pwm_score"         
#> [31] "edist_above"        "edist_below"        "edist_below_all"   
#> [34] "edist_max2"        
gvtrack.rm("vtrack1")
gvtrack.ls()
#>  [1] "vtrack2"            "vtrack3"            "vtrack4"           
#>  [4] "cov"                "motif_score"        "max_motif_score"   
#>  [7] "max_motif_pos"      "potts_max"          "potts_count"       
#> [10] "potts_p"            "potts_m"            "cg_count"          
#> [13] "cg_frac"            "at_pos"             "at_neg"            
#> [16] "at_both"            "g_frac"             "c_frac"            
#> [19] "masked_count"       "masked_frac"        "gc"                
#> [22] "value_track"        "value_track_max"    "spatial_pwm"       
#> [25] "regular_pwm"        "spatial_extended"   "window_pwm"        
#> [28] "window_spatial_pwm" "pwm_score"          "edist_above"       
#> [31] "edist_below"        "edist_below_all"    "edist_max2"        
```
