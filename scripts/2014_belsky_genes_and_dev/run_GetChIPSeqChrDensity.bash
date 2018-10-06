#!/bin/bash

# Enter the working directory
dir=/data/data2/jab112/2014_mnase_manuscript/scripts/2014_belsky_genes_and_dev/
cd $dir

# Compile the script if needed
javac -d bin/ -sourcepath src/ -cp lib/sam-1.67.jar src/mnasedens/ChIPSeqChrDensitySignal.java

# Set the parameters
bam_file_name=/data/illumina_pipeline/aligned_experiments/DM272/dm272.bam
output_header=/data/data2/jab112/2014_mnase_manuscript/datasets/chip_data/input/input_chip_seq_dm272_density_signal_chr_
output_distr=/data/data2/jab112/2014_mnase_manuscript/datasets/chip_data/input/signal_distribution/input_chip_seq_dm272_dist_output.csv
shift_var=75
bw=30

# Run the program
java -cp bin:lib/sam-1.67.jar mnasedens.ChIPSeqChrDensitySignal $bam_file_name $output_header $output_distr $shift_var $bw
