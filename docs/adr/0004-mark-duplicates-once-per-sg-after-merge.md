# Duplicates are marked once per SG, after the merge

When an SG has several alignment units, each unit is aligned and name-sorted on its own. The units are
then merged (`samtools merge -n`) into one query-grouped stream, and dupblaster marks duplicates on that
stream. Marking per unit would miss duplicates whose copies landed in different lanes or shards. The
price is a single merge task per SG that handles the whole SG's reads.
