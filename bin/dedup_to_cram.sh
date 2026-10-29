#!/usr/bin/env bash
# Read a query-grouped SAM/BAM stream on stdin, mark duplicates with dupblaster, coordinate-sort,
# and write SG_ID.cram (CRAM 3.0, + .crai) and SG_ID.markduplicates-metrics.
# usage: dedup_to_cram.sh SG_ID REF_FASTA SORT_THREADS WRITE_THREADS [DUPBLASTER_ARGS...]
set -euo pipefail

id=$1 fasta=$2 sort_threads=$3 write_threads=$4
shift 4

mkdir -p tmp
dupblaster --metrics-prefix metrics --tmp-dir tmp "$@" \
  | samtools sort -u -@"$sort_threads" -T tmp/sort \
  | samtools view --write-index -@"$write_threads" -T "$fasta" -O cram,version=3.0 -o "$id.cram" -

check_alignments.sh "$id.cram" "$fasta"
mv metrics.duplicate-metrics.tsv "$id.markduplicates-metrics"
