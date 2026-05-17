> **For AI coding agents working on the misha source code.** This is the
> *development* guide - workflow, hooks, build, test, release, contribution
> conventions. It is **not** user documentation. If you (or whoever invoked
> you) want to know how to *use* misha, read the
> [README](README.md) and the
> [package docs](https://tanaylab.github.io/misha/) instead.

## Initial Setup

After cloning, install the repo's git hooks (one-time, per clone):

tools/install-hooks.sh

This wires up a pre-commit hook (styler on staged R files, blocks `dev/`)
and a pre-push hook (`pkgbuild::compile_dll` on `src/` changes,
`devtools::document()` on `R/` changes).

## Running current R code

To run the current R code, use the following command:

R -e "devtools::load_all(export_all = FALSE)"

Do **not** use `library(name_of_package)` to load the package, this will load the package from the installed version, not the development version.

## Compilation Instructions

To compile the package after making C++ changes:

R -e 'pkgbuild::compile_dll(debug = FALSE)'

Or (for safety):

R -e 'pkgbuild::clean_dll(); pkgbuild::compile_dll(debug = FALSE)'

## Testing Instructions

To run specific tests:

R -e 'devtools::test(filter = "name_of_test_file_without_extension")'

To run a specific test in a file do (e.g.):

R -e "testthat::test_file('tests/testthat/test-pwm-sliding-window.R', desc = 'PWM sliding window MAX_POS mode works with iterator=20 and shifts')"

To run all tests in parallel (preferred - serial mode is dramatically slower):

R -e 'Sys.setenv(TESTTHAT_PARALLEL = "TRUE"); devtools::test()'

To check the package:

R -e "devtools::check()"

**Always prefer running all tests in parallel and not doing it one by one. It would be significantly slower without parallelism**

## Adding debug prints in C++ code

To add debug prints in C++ code, use REprintf.

Remove it afterwards, and compile the package with `R -e "pkgbuild::clean_dll(); pkgbuild::compile_dll(debug = F)"`.


## Documentation

R -e 'styler::style_pkg(indent_by = 4); devtools::document()'

## Adding new functions

Add @export to roxygen documentation and run documentation. Add to _pkgdown.yml. Test with `R -e 'devtools::test(filter = "name_of_test_file_without_extension")'` and run `R -e "devtools::check()"`.

## `dev/` Directory Structure

**Note:** `dev/` is an internal scratch directory for development artifacts (notes, benchmarks, scratch tests, debugging scripts, archived material, step-by-step skill guides). It is gitignored from the GitHub repo and only exists in local working checkouts that pull it explicitly. If you are reading this file on GitHub, the sections below reference `dev/` files that are not visible there - the inline summaries here are sufficient for contributing without `dev/` access.

Subdirectory layout (when present):

- `dev/notes/` - markdown notes, plans, design docs.
- `dev/benchmarks/` - performance benchmarks.
- `dev/tests/` - regression / scratch tests.
- `dev/debugging/` - bug repro scripts.
- `dev/archive/` - historical material.
- `dev/skills/` - step-by-step guides for recurring dev tasks (e.g. adding a new vtrack, publishing a release). Internal-only; if you have a checkout with `dev/` available, consult these before the corresponding task.

## Versioning Policy

Version bumps:
- **MINOR** - new exported functions, default changes, return-shape changes, or new on-disk formats.
- **PATCH** - bug fixes (even when output changes), perf work, new optional arguments with backward-compatible defaults, internal work.

NEWS bullets for output-changing fixes start with `**Behavior fix:**`; deliberate breaks start with `**Breaking:**`. Keep NEWS bullets short and user-facing - one or two lines for someone deciding whether to upgrade. Detailed rationale, benchmarks, and root-cause writeups belong in the commit message, not in NEWS.

## Publishing a Release

To publish a new misha release (conda package on anaconda.org):

1. **Update `NEWS.md`** with a new section at the top for the new version, per the policy above.
2. **Pick the bump type** (MINOR vs PATCH) and update the `Version:` field in `DESCRIPTION`.
3. **Commit, push to master.**
4. **Create a tag and GitHub release:**
   ```bash
   git tag v<VERSION>
   git push origin v<VERSION>
   gh release create v<VERSION> --title "<VERSION>" --notes "<release notes>"
   ```
   Pushing a `v*` tag triggers `.github/workflows/conda-publish.yml`, which builds the conda package on Ubuntu and macOS and uploads it to `anaconda.org/aviezerl/misha`.

Alternatively, trigger a manual conda build without tagging:
```bash
gh workflow run conda-publish.yml -f version=<VERSION>
```

## Branch Naming Convention

Use the following prefixes for branch names:

- `fix/` - bug fixes (e.g. `fix/numeric-stability-v5.5.2`)
- `feat/` - new features (e.g. `feat/pwm-sliding-window`)
- `refactor/` - code refactoring without behavior changes
- `docs/` - documentation-only changes
- `test/` - test-only changes

Never use personal name prefixes (e.g. `aviezerl/...`).

## GitHub PR Comments

When posting GitHub PR comments, do not include the "Generated with Claude Code" footer or reaction request.

## General Guidelines

- **Clarity and Conciseness**: Strive for clarity and conciseness in all development activities, from code to documentation.
- **Modularity**: Keep components modular and well-defined.
- **Testing**: All new features should be accompanied by tests.
- **Documentation**: Document new features and changes in the relevant files.

## Notes Policy

When adding new notes to the `dev/notes` directory, please adhere to the following guidelines:

- **Dated Notes**: If a note pertains to a specific date (e.g., a code review, a plan, a summary of work), include the date in the filename in `YYYY-MM-DD` format. For example, `2025-10-29_my-new-feature-plan.md`.
- **Categorization**: Place the note in the appropriate subdirectory (`features`, `optimizations`, `bug-fixes`, `specs`, `archive`).
