#!/bin/bash

# Enter the working directory
dir=/data/data2/jab112/2014_mnase_manuscript/scripts/2014_belsky_genes_and_dev/
cd $dir

# Remove the bin
rm -r bin/*/

# Compile the script if needed
javac -d bin/ -sourcepath src/ -cp lib/sam-1.67.jar src/mnasedens/GetFragmentDistributionOverFeature.java

# Set the parameters
bam_file_dir="/data/illumina_pipeline/aligned_experiments/"
bam_file_id=(243)
work_dir="/data/data2/jab112/2014_mnase_manuscript/figures/supplementary/sfigure_6/sfigure_6_datasets/"
feature_file_name=("left_mcm_cluster_oridb_feature_file_168_sites"
		   "right_mcm_cluster_oridb_feature_file_228_sites"
		  )

# Iterate through each
for d in ${bam_file_id[@]}
do
	
	# Get the bam_file
	bam_file=${bam_file_dir}DM${d}/dm${d}.bam

	# Iterate through the feature file names
	for f in ${feature_file_name[@]}
	do

		echo -e "Finding the frag distribution for ${f} over dm${d}.bam...\r"

		# Get the feature file
		feature_file=${work_dir}${f}.csv

		# Set the output file
		output_file=${work_dir}${f}_nucleosome_dm${d}_fragment_length_distribution.csv

		# Run the script
		java -cp bin:lib/sam-1.67.jar mnasedens.GetFragmentDistributionOverFeature $feature_file $bam_file $output_file -91 164

	done

done

echo -e "\n\tComplete!"
