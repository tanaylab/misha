
<!-- README.md is generated from README.Rmd. Please edit that file -->

# misha

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/misha)](https://CRAN.R-project.org/package=misha)
[![R-CMD-check](https://github.com/tanaylab/misha/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tanaylab/misha/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `misha` package is a toolkit for analysis of genomic data. it
implements an efficient data structure for storing genomic data, and
provides a set of functions for data extraction, manipulation and
analysis.

## Installation

You can install the released version of misha from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("misha")
```

Or from conda:

``` bash
conda install -c aviezerl r-misha
```

And the development version from GitHub with:

``` r
remotes::install_github("tanaylab/misha")
```

## Usage

See the
[Genomes](https://tanaylab.github.io/misha/articles/Genomes.html)
vignette for instructions on how to create a misha database for common
genomes.

See the [user
manual](https://tanaylab.github.io/misha/articles/Manual.html) for more
usage details.

### Using misha with an LLM agent

**Drop-in prompt (no clone needed).** Paste the block below into your
agent at the start of a misha task. It points the agent at the raw
files on GitHub, so it works without a local checkout:

````
Before writing any misha code, fetch and read:

- https://raw.githubusercontent.com/tanaylab/misha/master/agent-guides/misha-core.md  (mandatory: concepts + everyday recipes)
- https://raw.githubusercontent.com/tanaylab/misha/master/agent-guides/misha-anti-patterns.md  (silent footguns; cross-referenced from core)
- https://raw.githubusercontent.com/tanaylab/misha/master/agent-guides/misha-advanced.md  (consult on demand: 2D / Hi-C, PWM, import/export, new genomes)

Follow the conventions in those files. When you hit a recipe with an
"Avoid:" block, treat it as a hard rule.
````

For agents (Claude Code, Copilot, Cursor, etc.) writing misha analysis
code in a downstream project, point them at the maintained agent guides
in this repo:

- [`agent-guides/misha-core.md`](https://github.com/tanaylab/misha/blob/master/agent-guides/misha-core.md)
  - concepts, bootstrap, and the everyday recipes (intervals,
  annotation, distance, extract, vtracks, gscreen, gdist,
  gtrack.create). Start here.
- [`agent-guides/misha-advanced.md`](https://github.com/tanaylab/misha/blob/master/agent-guides/misha-advanced.md)
  - 2D / Hi-C pile-ups, insulation, sequence and PWM tracks, bulk
  import/export, new genomes and cross-species.
- [`agent-guides/misha-anti-patterns.md`](https://github.com/tanaylab/misha/blob/master/agent-guides/misha-anti-patterns.md)
  - silent footguns referenced inline from the above.
- [`agent-guides/skills/`](https://github.com/tanaylab/misha/tree/master/agent-guides/skills)
  - deep playbooks for specific tasks. Currently: [`importing-tracks`](https://github.com/tanaylab/misha/blob/master/agent-guides/skills/importing-tracks/SKILL.md) (format chooser across all `gtrack.*import*` variants + pre/post-import validation). Load when the task specifically calls for one of these.

The core guide is ~4k words and targets a system-prompt-sized context.
For Claude Code-style setups, dropping `misha-core.md` (or all three)
into the project’s `CLAUDE.md` / `AGENTS.md` is the intended use.

#### Running scripts from old versions of misha (\< 4.2.0)

Starting in `misha` 4.2.0, the package no longer stores global variables
such as `ALLGENOME` or `GROOT`. Instead, these variables are stored in a
special environment called `.misha`. This means that scripts written for
older versions of `misha` will no longer work. To run such scripts,
either add a prefix of `.misha$` to all those variables
(`.misha$ALLGENOME` instead of `ALLGENOME`), or run the following
command before running the script:

``` r
ALLGENOME <<- .misha$ALLGENOME
GROOT <<- .misha$GROOT
ALLGENOME <<- .misha$ALLGENOME
GINTERVID <<- .misha$GINTERVID
GITERATOR.INTERVALS <<- .misha$GITERATOR.INTERVALS
GROOT <<- .misha$GROOT
GWD <<- .misha$GWD
GTRACKS <<- .misha$GTRACKS
```
