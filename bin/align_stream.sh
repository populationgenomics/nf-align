#!/usr/bin/env bash
# Align one unit of input with bwa-mem3, writing uncompressed BAM to stdout (input read order, query-grouped).
# usage: align_stream.sh KIND SHARD INDEX_PREFIX SG_ID BWA_THREADS READ_THREADS FILES...
#   fastq:     FILES = R1 R2
#   fastq_ora: FILES = R1.ora R2.ora ORA_REFERENCE_TARBALL
#   cram|bam:  FILES = READS READS_INDEX REF_FASTA REF_FAI    (SHARD = "i,N" for bazam, or empty)
set -euo pipefail

kind=$1 shard=$2 index=$3 id=$4 bwa_threads=$5 read_threads=$6
shift 6
# -K fixes the batch size, so output doesn't change with the thread count (insert size inference is per batch)
bwa=(bwa-mem3 mem -t "$bwa_threads" -K 100000000 --bam -R "@RG\tID:${id}\tSM:${id}\tPL:ILLUMINA")

case $kind in
  fastq)
    "${bwa[@]}" "$index" "$1" "$2"
    ;;
  fastq_ora)
    # stream through FIFOs rather than decompress to disk (~200 GB for a 30x genome). bwa-mem3 can't see a
    # writer fail, so it would align a truncated stream and exit 0: wait on each orad to collect its status.
    tar -xf "$3"
    mkfifo R1.fq R2.fq
    orad -c --ora-reference oradata_homo_sapiens "$1" > R1.fq &
    r1_pid=$!
    orad -c --ora-reference oradata_homo_sapiens "$2" > R2.fq &
    r2_pid=$!
    # if bwa-mem3 exits early, an orad blocked opening its FIFO still holds our stdout, and the pipe never closes
    trap 'kill "$r1_pid" "$r2_pid" 2> /dev/null || true' EXIT
    "${bwa[@]}" "$index" R1.fq R2.fq
    wait "$r1_pid"
    wait "$r2_pid"
    ;;
  cram|bam)
    # bazam wants foo.crai / foo.bai, not foo.cram.crai / foo.bam.bai
    [ -e "${1%m}i" ] || ln -s "$(realpath "$2")" "${1%m}i"
    # interleaved pairs on stdin (-p -); pipefail catches a failure on either side
    bazam -Xmx16g -Dsamjdk.reference_fasta="$3" -n"$read_threads" -bam "$1" ${shard:+-s "$shard"} \
      | "${bwa[@]}" -p "$index" -
    ;;
  *)
    echo "align_stream.sh: unknown input kind '$kind'" >&2
    exit 1
    ;;
esac
