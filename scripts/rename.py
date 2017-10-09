#!/usr/bin/python

#This script is for renaming Illumina files received from UMGC at UMN based on a sample name file
#The data files received from UMGC may be in a format like 2015_1_S1_R1_001.fastq, but you want a format like descriptive_name_R1(or R2).fastq
#So, by providing a sample name file in the format below, you can rename the R1 and R2 ffiles accordingly
#1,15ND19-2
#2,15ND19-5
#3,15ND20-3
#etc...

#The format changes occasionally for UMGC files so this script might not work for you, be sure to check the latest format

#Make sure to module load python3 on MSI
#As of Octoboer 2017, this was Python 3.4.5 |Anaconda 2.1.0
#Run this in the same directory as fastq files, or use os.chdir to direct to the proper location

import glob, re, os

sample_file = open("sample_names.txt")

replace_key = {}

for line in sample_file:
	fields = line.rstrip("\n").split(",")
	replace_key[fields[0]] = fields[1]

for filename in glob.glob('*.fastq'):
	sample_code = re.search("\d+_(\d+)_S\d+_(R\d)_001\.fastq", filename)
	if sample_code.group(1) in replace_key.keys():
		new_name = re.sub(r"\d+_\d+_S\d+_R\d_001(\.fastq)", replace_key.get(sample_code.group(1)), filename)
		os.rename(filename, new_name + "_" + sample_code.group(2) + ".fastq")

sample_file.close()
