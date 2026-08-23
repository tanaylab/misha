# Toolkit for analysis of genomic data

'misha' package is intended to help users to efficiently analyze genomic
data achieved from various experiments.

## Details

For a complete list of help resources, use
[`library(help = "misha")`](https://rdrr.io/r/base/library.html).

The following options are available for the package. Use 'options'
function to alter the value of the options.

|  |  |  |
|----|----|----|
| NAME | DEFAULT | DESCRIPTION |
| gmax.data.size | AUTO | Auto-configured based on system RAM and processes. |
|  |  | Formula: min((RAM \* 0.7) / gmax.processes, 10GB). |
|  |  | Controls max in-memory result buffer size for streaming/sampling. |
|  |  | Operations like gextract, gscreen use this as upper bound. |
| gbig.intervals.size | 1000000 | Threshold for converting interval sets to disk-based |
|  |  | "big" format. When interval count exceeds this, intervals |
|  |  | are stored per-chromosome instead of all in memory. |
|  |  | Note: Independent of gmax.data.size (different purposes). |
| gmax.mem.usage | 10000000 | Maximal memory consumption of all child |
|  |  | processes in KB before the limiting algorithm is invoked. |
| gmax.processes | AUTO | Auto-configured to min(32, 70% of CPU cores). |
|  |  | Maximal number of processes for multitasking. |
| gmax.processes2core | 2 | Maximal number of processes per CPU core for multitasking |
| gseq.extract.probe.usec | 100 | gseq.extract times its first reads |
|  |  | and distributes across processes when the measured cost per |
|  |  | interval reaches this many microseconds. 0 always distributes, |
|  |  | a large value never does. |
| gmin.scope4process | 10000 | Minimal scope range (for 2D: surface) assigned to a process |
|  |  | in multitasking mode. |
| gbuf.size | 1000 | Size of track expression values buffer. |
| gmultitask.max.records.factor | 64 | Multiplier to inflate multitask |
|  |  | max_records estimates to avoid under-allocation. |
| gtrack.chunk.size | 100000 | Chunk size in bytes of a 2D track. If '0' chunk size |
|  |  | is unlimited. |
| gtrack.num.chunks | 0 | Maximal number of 2D track chunks simultaneously stored |
|  |  | in memory. |
| gmultitasking | TRUE | Enable/disable automatic parallelization. Some functions |
|  |  | may choose single-threaded mode for very small workloads. |
| gpermissions.umask | "0007" | umask applied around misha's own writes |
|  |  | into a database: 660 files, 770 directories - group-writable, |
|  |  | no world access - matching what the C++ layer enforces. Set to |
|  |  | "0002" for the world-readable 664/775 that misha produced before |
|  |  | 5.11.21, or to NULL to leave the process umask untouched. The |
|  |  | umask is restored as soon as the write returns. |

More information about the options can be found in 'User manual' of the
package.

## See also

Useful links:

- <https://tanaylab.github.io/misha/>

- <https://github.com/tanaylab/misha>

- Report bugs at <https://github.com/tanaylab/misha/issues>

## Author

**Maintainer**: Aviezer Lifshitz <aviezer.lifshitz@weizmann.ac.il>

Authors:

- Aviezer Lifshitz <aviezer.lifshitz@weizmann.ac.il>

- Misha Hoichman <misha@hoichman.com>

- Eitan Yaffe <eitan.yaffe@weizmann.ac.il>

- Amos Tanay <amos.tanay@weizmann.ac.il>

Other contributors:

- Weizmann Institute of Science \[copyright holder\]
