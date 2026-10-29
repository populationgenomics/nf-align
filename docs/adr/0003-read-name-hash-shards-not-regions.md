# CRAM/BAM inputs are sharded by read-name hash, not by genomic region

CRAM/BAM inputs are split into `realignment_shards` pieces with `bazam -s i,N`, which picks read pairs by
a hash of the read name. The original request said "divided roughly by regions". Hash shards were chosen
because both mates always land in the same shard, unmapped and unplaced reads are covered with no special
case, shards come out close to equal in size, and no region list per reference build is needed.

## Considered Options

- **Region shards (`bazam -L`, or contig lists).** Pairs that span a boundary, or whose mates sit on
  different chromosomes, can be emitted twice or dropped. Unmapped reads would need a shard of their own.
  The shards would also be uneven.
