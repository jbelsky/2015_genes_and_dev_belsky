# Clear the workspace
rm(list = ls())
graphics.off()

# Get the chr lengths
bai_file.v = system(command = "samtools idxstats /data/illumina_pipeline/aligned_experiments/DM242/dm242.bam",
		    intern = T
		   )

# Convert to a list
bai_file.l = strsplit(bai_file.v, split = "\t")

# Convert to a dataframe
bai_file.df = as.data.frame(matrix(unlist(bai_file.l), ncol = 4, byrow = T))
colnames(bai_file.df) = c("chr", "length", "read_num", "unaligned")

# Convert the length into a numeric
bai_file.df$length = as.numeric(bai_file.df$length)

# Get the chromosome length
chr_length.v = bai_file.df$length[1:16]

# Get the percentage of each chr length
chr_percent.v = chr_length.v / sum(chr_length.v)

# Select the total number of positions for the signal distribution
total_position = 2000

# Create the output file
output_file = paste("/data/data2/jab112/2014_mnase_manuscript/datasets/",
		    "genome_background_feature_file_2000_sites.csv", sep = ""
		   )

# Get the number of positions sampled for each chromosome
chr_samp.v = round(total_position * chr_percent.v)

if(sum(chr_samp.v) != total_position){

	# Get the diff
	diff = sum(chr_samp.v) - total_position
	
	# Randomly modify the count from one chromosome
	chr_change = sample(1:16, 1)

	chr_samp.v[chr_change] = chr_samp.v[chr_change] - diff

}

# Make the storage dataframe
feature.df = data.frame(name = paste("sample_genome_", 1:total_position, sep = ""),
			chr = rep(1:16, chr_samp.v),
			pos = 0,
			strand = "+"
		       )

# Iterate through each chr
for(i in 1:16){

	cat("Sampling chr ", i, "...\r", sep = "")

	# Get indices based on the total_position and chr_percent.v
	sample_pos.v = sample(500:(chr_length.v[i] - 500), chr_samp.v[i], replace = FALSE)

	# Update the position
	feature.df[which(feature.df$chr == i),"pos"] = sort(sample_pos.v)


}

cat("\n\tComplete!\n")

# Write the output table
write.table(feature.df, file = output_file, sep = ",", col.names = T, row.names = F, quote = F)
