This is the folder where reference genomes would go (not add to GitHub due to size).

These are the same reference genomes and gff3 files that were deposited into MycoCosm.

For simplicity, primary contigs and haplotigs were separated into two different files.

Annotations for primary contigs and haplotigs were also separated into two separate files.

Reference genomes were indexed like this:

module load bwa/0.7.15 samtools/1.5
bwa index ref.fa
samtools faidx ref.fa
