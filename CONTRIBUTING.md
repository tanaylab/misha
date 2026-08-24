# Contributing to misha

## First-time setup

After cloning, install the repo’s git hooks (one-time, per clone):

``` bash
tools/install-hooks.sh
```

This wires up:

- **pre-commit** — runs `styler` on staged R files (matches the CI style
  pass) and refuses to commit anything under `dev/`.
- **pre-push** — runs
  [`pkgbuild::compile_dll()`](https://pkgbuild.r-lib.org/reference/compile_dll.html)
  when pushed commits touch `src/`, and `devtools::document()` when they
  touch `R/` (aborts the push if the regenerated `man/` or `NAMESPACE`
  differ from what’s committed).

## Branch naming

- `fix/...` — bug fixes
- `feat/...` — new features
- `refactor/...` — code restructuring without behavior change
- `docs/...` — documentation only
- `test/...` — tests only

## Tests

The full test suite is run on the Tanay Lab cluster (some tests require
cluster-only data and resources). For external contributions, run the
subset relevant to your change locally and rely on CI for the rest.

``` r

# run a single test file
testthat::test_file("tests/testthat/test-foo.R")

# run everything available locally, in parallel (some tests will skip
# without cluster-only data)
Sys.setenv(TESTTHAT_PARALLEL = "TRUE")
devtools::test()
```

### The shared test database

Most tests need the lab’s 5.9 TB test database on NFS, which is far too
large to copy. They do not open it directly. Every test process builds a
private **overlay** of it under
[`tempdir()`](https://rdrr.io/r/base/tempfile.html) - directories misha
writes into (`tracks`, `tracks/test`, `tracks/temp`, …) are real
directories, and the tracks and interval sets inside them are symlinks.
The overlay reads exactly like the shared database and costs about 0.1 s
and 700 KB, and every track, interval set or attribute a test writes
lands in [`tempdir()`](https://rdrr.io/r/base/tempfile.html) instead of
in lab data.

Concurrent runs are therefore fine: two worktrees, two people, or a
local run alongside CI. Nothing in the suite may root a session at
`shared_test_db_path()`; use one of

- `test_db_root()` - this process’s overlay, shared by every file in it.
- `create_isolated_test_db()` - a fresh overlay for the calling file,
  removed when that file finishes.
- `load_test_db()` - `test_db_root()` plus an emptied `tracks/temp`.

Nothing needs cleaning up in the shared database: the overlays never
write to it, and a killed run leaves only symlinks under `TMPDIR`.

Overlays and the databases some tests build from scratch live under
[`tempdir()`](https://rdrr.io/r/base/tempfile.html), so point `TMPDIR`
at a filesystem with several GB free before starting R:

``` sh
TMPDIR=/path/with/space R -e 'Sys.setenv(TESTTHAT_PARALLEL="TRUE"); devtools::test()'
```

On a shared machine `/tmp` often sits on a full root filesystem. Running
out of space mid-run leaves partially copied test databases, and the
resulting failures point at whichever test file happened to use them
rather than at the disk — expect errors like
`Interval test.fixedbin does not exist`. Failures that move between
files across runs are a symptom of this, not of flaky tests.

GitHub Actions runs `R CMD check` on every push and PR — treat that as
the authoritative signal.

## Other dev commands

Compile after C++ changes:

``` r

pkgbuild::compile_dll(debug = FALSE)
```

Regenerate documentation:

``` r

devtools::document()
```
