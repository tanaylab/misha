# Contributing to misha

## First-time setup

After cloning, install the repo's git hooks (one-time, per clone):

```bash
tools/install-hooks.sh
```

This wires up:

- **pre-commit** — runs `styler` on staged R files (matches the CI style
  pass) and refuses to commit anything under `dev/`.
- **pre-push** — runs `pkgbuild::compile_dll()` when pushed commits touch
  `src/`, and `devtools::document()` when they touch `R/` (aborts the push
  if the regenerated `man/` or `NAMESPACE` differ from what's committed).

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

```r
# run a single test file
testthat::test_file("tests/testthat/test-foo.R")

# run everything available locally, in parallel (some tests will skip
# without cluster-only data)
Sys.setenv(TESTTHAT_PARALLEL = "TRUE")
devtools::test()
```

A parallel run gives each test file its own copy of the test database under
`tempdir()`, so point `TMPDIR` at a filesystem with several GB free before
starting R:

```sh
TMPDIR=/path/with/space R -e 'Sys.setenv(TESTTHAT_PARALLEL="TRUE"); devtools::test()'
```

On a shared machine `/tmp` often sits on a full root filesystem. Running out
of space mid-run leaves partially copied test databases, and the resulting
failures point at whichever test file happened to use them rather than at the
disk — expect errors like `Interval test.fixedbin does not exist`. Failures
that move between files across runs are a symptom of this, not of flaky tests.

GitHub Actions runs `R CMD check` on every push and PR — treat that as the
authoritative signal.

## Other dev commands

Compile after C++ changes:

```r
pkgbuild::compile_dll(debug = FALSE)
```

Regenerate documentation:

```r
devtools::document()
```
