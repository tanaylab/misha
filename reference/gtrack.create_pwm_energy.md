# Creates a new track from PSSM energy function

Creates a new track from PSSM energy function.

## Usage

``` r
gtrack.create_pwm_energy(
  track = NULL,
  description = NULL,
  pssmset = NULL,
  pssmid = NULL,
  prior = NULL,
  iterator = NULL
)
```

## Arguments

- track:

  track name

- description:

  a character string description

- pssmset:

  name of PSSM set: 'pssmset.key' and 'pssmset.data' must be presented
  in 'GROOT/pssms' directory

- pssmid:

  PSSM id

- prior:

  prior

- iterator:

  track expression iterator for the newly created track

## Value

None.

## Details

This function creates a new track with values of a PSSM energy function.
PSSM parameters (nucleotide probability per position and pluralization)
are determined by 'pssmset' key and data files ('pssmset.key' and
'pssmset.data'). These two files must be located in 'GROOT/pssms'
directory. The type of the created track is determined by the type of
the iterator. 'description' is added as a track attribute.

## See also

[`gtrack.create`](https://tanaylab.github.io/misha/reference/gtrack.create.md),
[`gtrack.2d.create`](https://tanaylab.github.io/misha/reference/gtrack.2d.create.md),
[`gtrack.create_sparse`](https://tanaylab.github.io/misha/reference/gtrack.create_sparse.md),
[`gtrack.smooth`](https://tanaylab.github.io/misha/reference/gtrack.smooth.md),
[`gtrack.modify`](https://tanaylab.github.io/misha/reference/gtrack.modify.md),
[`gtrack.rm`](https://tanaylab.github.io/misha/reference/gtrack.rm.md),
[`gtrack.info`](https://tanaylab.github.io/misha/reference/gtrack.info.md),
[`gdir.create`](https://tanaylab.github.io/misha/reference/gdir.create.md)

## Examples

``` r
# \donttest{
gdb.init_examples()
gtrack.create_pwm_energy("pwm_energy_track", "Test track", "pssm",
    3, 0.01,
    iterator = 100
)
gextract("pwm_energy_track", gintervals(1, 0, 1000))
#>    chrom start  end pwm_energy_track intervalID
#> 1   chr1     0  100        -53.42195          1
#> 2   chr1   100  200        -52.64308          1
#> 3   chr1   200  300        -52.11026          1
#> 4   chr1   300  400        -51.71758          1
#> 5   chr1   400  500        -49.75199          1
#> 6   chr1   500  600        -44.75133          1
#> 7   chr1   600  700        -49.11913          1
#> 8   chr1   700  800        -48.81431          1
#> 9   chr1   800  900        -48.35187          1
#> 10  chr1   900 1000        -49.70630          1
# }
```
