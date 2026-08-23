# Runs R commands on a cluster

Runs R commands on a cluster that supports SGE.

## Usage

``` r
gcluster.run(
  ...,
  opt.flags = "",
  max.jobs = 400,
  debug = FALSE,
  R = "R",
  control_dir = NULL
)
```

## Arguments

- ...:

  R commands

- opt.flags:

  optional flags for qsub command

- max.jobs:

  maximal number of simultaneously submitted jobs

- debug:

  if 'TRUE', additional reports are printed

- R:

  command that launches R

- control_dir:

  directory where the control files are stored. Note that this directory
  should be accessible from all nodes. If 'NULL', a temporary directory
  would be created under the current misha database.

## Value

Return value ('retv') is a list, such that 'retv\[\[i\]\]' represents
the result of the run of command number 'i'. Each result consists of 4
fields that can be accessed by 'retv\[\[i\]\]\$FIELDNAME':

|  |  |
|----|----|
| *FIELDNAME* | *DESCRIPTION* |
| exit.status | Exit status of the command. Possible values: 'success', 'failure' or 'interrupted'. |
| retv | Return value of the command. |
| stdout | Standard output of the command. |
| stderr | Standard error of the command. |

## Details

This function runs R commands on a cluster by distributing them among
cluster nodes. It must run on a machine that supports Sun Grid Engine
(SGE). The order in which the commands are executed can not be
guaranteed, therefore the commands must be inter-independent.

Optional flags to 'qsub' command can be passed through 'opt.flags'
parameter. Users are strongly recommended to use only '-l' flag as other
flags might interfere with those that are already used (-terse, -S, -o,
-e, -V). For additional information please refer to the manual of
'qsub'.

The maximal number of simultaneously submitted jobs is controlled by
'max.jobs'.

Set 'debug' argument to 'TRUE to allow additional report prints.

'gcluster.run' launches R on the cluster nodes to execute the commands.
'R' argument specifies how R executable should be invoked.

## Examples

``` r
# \donttest{
gdb.init_examples()
# Run only on systems with Sun Grid Engine (SGE)
if (FALSE) {
    v <- 17
    gcluster.run(
        gsummary("dense_track + v"),
        {
            intervs <- gscreen("dense_track > 0.1", gintervals(1, 2))
            gsummary("sparse_track", intervs)
        },
        gsummary("rects_track")
    )
}
# }
```
