#!/usr/bin/env bash
# Split a task's CPUs between the tools in an align/dedup pipe; the only place these ratios live.
# usage: eval "$(threads.sh CPUS [--no-aligner])"
# sets:
#   BWA_THREADS    bwa-mem3 -t: 2/3 of the CPUs, as the bwa-mem3 docs suggest (0 with --no-aligner)
#   READ_THREADS   input decoding alongside the aligner: bazam -n, or samtools merge -@
#   SORT_THREADS   samtools sort -@ while upstream is still running (spill compression);
#                  one CPU is left for dupblaster, which has no threads option and uses ~1-1.5 cores
#   WRITE_THREADS  samtools view -@ for the final CRAM; sort holds its output until its input ends,
#                  so this mostly runs after alignment has finished and can use nearly every CPU
set -euo pipefail

cpus=$1
aligner=${2:-}
max() { echo $(($1 > $2 ? $1 : $2)); }

if [ "$aligner" = --no-aligner ]; then
  bwa=0
else
  bwa=$(max 1 $((cpus * 2 / 3)))
fi
io=$(max 1 $((cpus - bwa)))
read=$(max 1 $((io / 2)))

echo "BWA_THREADS=$bwa"
echo "READ_THREADS=$read"
echo "SORT_THREADS=$(max 1 $((io - read - 1)))"
echo "WRITE_THREADS=$(max 1 $((cpus - 1)))"
