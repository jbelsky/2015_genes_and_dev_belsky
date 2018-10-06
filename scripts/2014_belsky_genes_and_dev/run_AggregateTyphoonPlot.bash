#!/bin/bash

# Enter the working directory
dir=/data/data2/jab112/2014_mnase_manuscript/scripts/2014_belsky_genes_and_dev/
cd $dir

# Clear the bin/
rm -r bin/*/

# Compile the script if needed
javac -d bin/ -sourcepath src/ -cp lib/sam-1.67.jar src/twodimplot/AggregateTyphoonPlot.java

# Set the dir
dir="/data/data2/jab112/2014_mnase_manuscript/figures/supplementary/sfigure_6/sfigure_6_datasets/"

# Set the file names
feature_file=${dir}right_mcm_cluster_oridb_feature_file_228_sites.csv
bam_file="/data/illumina_pipeline/aligned_experiments/DM243/dm243.bam"
output_file=${dir}right_mcm_cluster_oridb_feature_file_228_sites_agg_typhoon_plot.csv

# Run the program
java -cp bin/:lib/sam-1.67.jar twodimplot.AggregateTyphoonPlot $bam_file $feature_file $output_file
