# Information about reference genomes
This is the folder where reference genomes would go (not added to GitHub due to size).

These are the same reference genomes and gff3 files that were deposited into MycoCosm and NCBI (https://www.ncbi.nlm.nih.gov/bioproject/398546).

For simplicity, primary contigs and haplotigs were separated into two different files.

Annotations for primary contigs and haplotigs were also separated into two separate files.

Reference genomes were indexed like this:
module load bwa/0.7.15 samtools/1.5
bwa index ref.fa
samtools faidx ref.fa

All files can be found at this link: https://s3.msi.umn.edu/1990_2015_pca_project/ref_genome_files.tar.gz
