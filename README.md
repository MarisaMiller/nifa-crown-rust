# NIFA Oat Crown Rust Project
This repository contains all scripts and files to re-create the analysis for the NIFA Oat Crown Rust project. Please forgive any clumsy git use or scripting, this is my first attempt at organizing a bioinformatics project this way.

Some files were not uploaded here due to size limitations. Notes have been added to directories to indicate this when possible. Links to large files have been provided whenever possible (such as annotation files), and links to raw data at NCBI have been included so analysis can be replicated from the beginning.

Note, see README files and/or job scripts in all directories as well for detailed commands, parameters, and versions.

In this project, I re-sequenced 60 isolates of *Puccinia coronata* f. sp. *avenae* (*Pca*), using paired-end sequencing (100bp PE). This sequencing was done by the UMGC sequencing core (100 ng of DNA used for TruSeq Nano DNA procedure and a 350 bp insert size). 2 batches of 30 libraries were multiplexed and sequenced in 3.5 lanes (HiSeq 2500, High Output Mode, 125 bp paired-end reads) at the University of Minnesota Genomics Center (UMGC) (MN, USA).

Data was deposited by UMGC into the relevant directories for our group, and subsequently I copied all raw FASTQ sequences to the raw_data directory. 12SD80 and 12NC29 sequencing data generated in the *Pca* genome sequencing project were also placed in the raw_data directory. Then, the rename.py script in the scripts directory was used to rename the files.

## General workflow to reproduce my analysis.
1. Download raw data from NCBI (see the BioProject/SRA information in the raw_data directory). 
1. Run the QC_mapping_pipeline.sh script in the scripts directory (script has detailed comments and instructions). 
2. Examine the quality stats in the analysis/trimmed_reads directory to determine if all samples can be used for further analysis. Run the count_reads.sh script in analysis/trimmed_reads to get counts of surviving reads after trimming.
3. Examine the mapping statistics with the generate_stats.sh script in analysis/mapped_reads. Then, use bedtools to generate per-sample average coverage (in the coverage_stats.sh script) and then use the mapping_genome_coverage.R script in the analysis/mapped_reads directory to generate plots.
4. Calculate coverage across the binned genome with the genome_coverage.pbs script in analysis/genome_wide_coverage/. Then, use the window_genome_coverage.R in the same directory to refine the analysis and plot the results.
5. Call SNPs with FreeBayes, perform variant filtering, and generate statistics. This is outlined in analysis/snp_calling/freebayes/freebayes_pipeline.sh, and in this script are also various steps to get vcf files ready for comparison with other variant callers and to combine vcfs.
6. To compare FreeBayes to other callers, follow the scripts in analysis/snp_calling/gatk (Gatk_pipeline.sh) and analysis/snp_calling/samtools (samtools_pipeline.sh). Then, use the various scripts in the analysis/snp_calling/compare_gatk_freebayes_samtools/ directory to perform analysis of shared and unique variants (compare_variant_callers.sh and get_unique_variants.sh). To overlap variants to genes and repeats, and get surrounding GC content, use the variant_featureOverlap_andGC.sh script. Finally, plot results with UpSet_variant_caller_and_year_comparison.R, genomic_neighborhood_shared_unique_variants.R, and gc_shared_unique_variants.R.
7. Run kWIP using the kwip.sh script in analysis/kWIP. I have an installation of kWIP in my scripts directory, but you will need to install yourself. Plot results with the img_mod_from_original.R script.
8. Assess isolate contamination with vcfR. See the check_contamination_script.R script in analysis/dapc_poppr/poppr_analysis/.
9. Annotate variants using Annovar with the annotate_variants.sh script in analysis/snp_calling/variant_annotation/
10. Perform DAPC analysis and generate phylogenetic tree with the dapc_script.R script in analysis/dapc_poppr/adegenet_dapc_analysis/
11. Perform LD analysis with two different measures: one method using PoppR in analysis/dapc_poppr/dapc_poppr/poppr_analysis/poppr_LD_script.R and the other method using vcftools in analysis/r_squared/. First run the shell script calc_R2.sh and then plot the results with r2_analysis.R.
12. Use PopGenome to perform population genomic analyses. See the README and the popGenome_analysis.R in analysis/popgenome/ for details of all analyses.
13. Lastly, thoroughly characterize effectors with the scripts in analysis/effector_analysis/. See the main script general_effector_analysis.R, but also see the scripts in the subdirectories to generate files needed for main analysis. This include the blast.pbs scipt in effector_analysis/blast/, the effector_coverage.pbs and effector_coverage.R scripts in effector_analysis/pav_analysis/, the subset_alignments.pbs in effector_analysis/subset_alignments/, and the effector_variants.sh script in effector_analysis/variation_analysis/. First run all supporting scripts in these directories and the run the summary script in analysis/effector_analysis/.

## Other analyses and information
Analysis of phenotypic data is in the isolate_phenotyping directory, and information for upload of raw data to NCBI is in the ncbi_upload directory.
Script and data for generating the plots of oat and buckthorn distribution is in the buckthorn_oat_analysis folder in the analysis directory.
