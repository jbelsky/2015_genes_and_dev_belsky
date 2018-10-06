#!/bin/bash

# Enter the working directory
dir=/data/data2/jab112/2014_mnase_manuscript/scripts/2014_belsky_genes_and_dev/
cd $dir

# Clear the bin
rm -r bin/*/

# Compile the script if needed
javac -d bin/ -sourcepath src/ src/acsmotif/GetFeatureMotifSeq.java

# Set the parameters
input_feature="/data/data2/jab112/2014_mnase_manuscript/datasets/jab112_yeast_feature_files/replication_origins/oridb_acs_feature_file_curated_798_sites_timing_whitehouse_raw_oem.csv"
output_feature="/data/data2/jab112/2014_mnase_manuscript/datasets/jab112_yeast_feature_files/replication_origins/oridb_acs_feature_file_curated_798_sites_timing_whitehouse_raw_oem_acs_seq.csv"
fasta_dir="/data/data2/jab112/2014_mnase_manuscript/datasets/jab112_yeast_feature_files/genome/sacCer2_sgdR61_chr"
win=33;

# Run the program
java -cp bin acsmotif.GetFeatureMotifSeq $input_feature $output_feature $fasta_dir $win
