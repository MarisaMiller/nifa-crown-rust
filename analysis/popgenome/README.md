# Details for preparing files for PopGenome analysis
In this directory are scripts and files needed to perform analysis with the PopGenome R package.

The 1990, 2015, and combined directories are linked here from the ../snp_calling/freebayes location as the vcf files separated by chromosome were generated during that analysis. See the freebayes_pipeline.sh script for details.

Here are the commands used to generate the gff and reference files split by contig (I have tabix and bioawk installed on MSI and in my $PATH):
```
#For preparing GFF files
ln -s ../../reference_genomes/Puccinia_coronata_avenae_12SD80.p.gff3
(grep ^"#" Puccinia_coronata_avenae_12SD80.p.gff3; grep -v ^"#" Puccinia_coronata_avenae_12SD80.p.gff3 | sort -k1,1 -k4,4n) | bgzip > Puccinia_coronata_avenae_12SD80.p.gff3.gz
tabix -p gff Puccinia_coronata_avenae_12SD80.p.gff3.gz

mkdir GFF
for contig in $(bioawk -c gff '{print $seqname}' Puccinia_coronata_avenae_12SD80.p.gff3.gz | sort | uniq)
do
#Make sure the gff files end with gff, not gff3, or else they can't be loaded by PopGenome
tabix -h Puccinia_coronata_avenae_12SD80.p.gff3.gz $contig > ./GFF/${contig}.gff
done

mkdir FASTA
ln -s ../../reference_genomes/Puccinia_coronata_avenae_12SD80.primary.fa

while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.fa
        echo $line > ./FASTA/$outfile
    else
        echo $line >> ./FASTA/$outfile
    fi
done < Puccinia_coronata_avenae_12SD80.primary.fa
```

An R bug occurs if GFF files have places in the last column that contain single quotes (see this blog for more information on this issue http://kbroman.org/blog/2017/08/08/eof-within-quoted-string/), and then they can't be successfully read in with PopGenome. I notified the author of PopGenome about the bug, but to get around this issue for now I will remove the offending quotes in the 53 GFF files that have them. This might make readability of certain information in the INFO column of the gff files worse, for example 5'-3' will now be 5-3, but it is the only workaround for now.
```
for F in $(grep -l "'" ./GFF/*.gff)
do
sed -i -e "s/'//g" "$F"
done
```

Here is a version that only removes files with odd numbers of quotes, which I later decided not to use.
```
for F in $(grep -l "'" ./GFF/*.gff)
do
if [[ $(( $(grep -o "'" "$F" | wc -l) % 2 )) -eq 1 ]]
then
sed -i -e "s/'//g" "$F"
fi
done
```

Make directories with bgzip compressed and tabix indexed split vcf files
```
mkdir 1990_gz
cp 1990/* 1990_gz/
cd 1990_gz
for vcf in $(ls *.vcf)
do
bgzip $vcf
tabix -p vcf $vcf.gz
done

mkdir 2015_gz
cp 2015/* 2015_gz/
cd 2015_gz
for vcf in $(ls *.vcf)
do
bgzip $vcf
tabix -p vcf $vcf.gz
done
```

## To get effector, BUSCO, and haustorially expressed gene lists 
For BUSCO list, BUSCO v2.0 was run in protein mode with all predicted proteins on primary contigs in 12SD80. The basidiomycota BUSCO set was used for comparison.

For haustorially expressed genes and effectors, the 12SD80 primary contig genes from Cluster 4 and 5 from Miller et al, mBio, 2018 were used.

For effectors, all effectors predicted with EffectorP on 12SD80 primary contigs in the same paper mentioned above were used.

Then to get gff for just effector gene bodies, I used these commands:
```
grep -f ./gene_lists/effectors_list_sd80_p.txt ../../reference_genomes/Puccinia_coronata_avenae_12SD80.p.gff3 | grep "gene" > ./gene_lists/effector_genes_sd80_p.gff3

grep -f ./gene_lists/haus_expr_secreted_effectors_list_sd80_p.txt ../../reference_genomes/Puccinia_coronata_avenae_12SD80.p.gff3 | grep "gene" > ./gene_lists/haus_expr_secreted_effectors_sd80_p.gff3 
```
