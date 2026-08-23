# Create directories needed for track creation

This function creates the directories needed for track creation. For
example, if the track name is 'proj.sample.my_track', this function
creates the directories 'proj' and 'sample'. Use this function with
caution - a long track name may create a deep directory structure.

## Usage

``` r
gtrack.create_dirs(track, mode = "0777")
```

## Arguments

- track:

  name of the track

- mode:

  see 'dir.create'

## Value

None.

## Examples

``` r

gdb.init_examples()

# This creates the directories 'proj' and 'sample'
gtrack.create_dirs("proj.sample.my_track")
```
