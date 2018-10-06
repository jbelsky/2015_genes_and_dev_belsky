#!/bin/bash

# Enter the working directory
dir=/data/data2/jab112/2014_mnase_manuscript/scripts/2014_belsky_genes_and_dev/
cd $dir

# Compile the script if needed
javac -d bin/ -sourcepath src/ -cp lib/sam-1.67.jar src/mnasedens/GetMNaseChrDensity.java

# Set the dm_id
dm_id=(466 467 468 469)

# Set the header
header=("g1_nuc_cdc6" "g1_nuc_cdc6" "g1_nuc_wt_37" "g1_nuc_wt_37")

# Set the non-changing parameters
bam_file_dir=/data/illumina_pipeline/aligned_experiments/
output_file_dir=/data/data2/jab112/2014_mnase_manuscript/datasets/nuc_data/raw_data2/
output_file_distr_dir=/data/data2/jab112/2014_mnase_manuscript/datasets/nuc_data/signal_distribution2/
frag_low=150
frag_high=175
bw=20

# Iterate through each dm_id
for(( i=0; i<${#dm_id[@]}; i++ ))
do

	# Get the dm_id
	id=${dm_id[$i]}

	# Set the parameters
	bam_file_name=${bam_file_dir}DM${id}/dm${id}.bam
	output_header=${output_file_dir}${header[$i]}_dm${id}_nuc_${frag_low}_${frag_high}_density_signal_chr_
	output_distr=${output_file_distr_dir}${header[$i]}_dm${id}_nuc_${frag_low}_${frag_high}_dist_output.csv

	# Run the program
	java -cp bin:lib/sam-1.67.jar mnasedens.GetMNaseChrDensity $bam_file_name $output_header $output_distr $frag_low $frag_high $bw

done
