In this directory would be gzip compressed fastq files that contain adapter-trimmed reads after running the QC_mapping_pipeline.sh script. These were not uploaded due to size limitations.

There are also plots of the quality of each library and also statistics. Also included is the script used to get counts of surviving PE reads after trimming.

Subsequently, all of the trimmed read files were uploaded to s3 (see the s3.sh script) and removed from MSI primary storage, and can be accessed at the following link: https://s3.msi.umn.edu/nifa_crown_rust/trimmed_data.tar.gz
