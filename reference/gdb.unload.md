# Unloads the genome database

Unloads the currently active genome database and clears the session
state.

## Usage

``` r
gdb.unload()
```

## Value

None.

## Details

Resets the misha session to an uninitialized state: the database root,
the working directory, the genome intervals, the virtual tracks, the
loaded datasets, the chromosome aliases and the track / intervals-set
auto-completion symbols are all cleared. After calling this function
[`gdb.init`](https://tanaylab.github.io/misha/reference/gdb.init.md) (or
[`gsetroot`](https://tanaylab.github.io/misha/reference/gdb.init.md))
must be called again before any other misha function is used.

Package internals are preserved, so the package itself stays loaded;
this is not the same as detaching the package.

## See also

[`gdb.init`](https://tanaylab.github.io/misha/reference/gdb.init.md),
[`gsetroot`](https://tanaylab.github.io/misha/reference/gdb.init.md),
[`gdb.reload`](https://tanaylab.github.io/misha/reference/gdb.reload.md)

## Examples

``` r

gdb.init_examples()
gdb.unload()
gdb.init_examples()
```
