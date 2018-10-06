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
bai_file.df$length = as.numeric(bai_file.df$read_num)

# Get the percentage of each chr length
chr_percent.v = bai_file.df$length[1:16] / sum(bai_file.df$length[1:16])

# Select the total number of positions for the signal distribution
total_position = 1.2E6

# Create the output file
output_file = paste("/data/data2/jab112/2014_mnase_manuscript/datasets/subnuc_data/signal_distribution/",
		    "g2_dm242_subnuc_0_120_dist_output.csv", sep = ""
		   )

# Create the header file
input_header = paste("/data/data2/jab112/2014_mnase_manuscript/datasets/subnuc_data/quant_normalize/",
		     "g2_dm242_subnuc_0_120_density_signal_chr_", sep = ""
		    )

# Iterate through each chr
for(i in 1:16){

	cat("Sampling chr ", i, "...\r", sep = "")

	# Load the signal distribution file
	sig.df = read.csv(paste(input_header, i, "_quant_normalize.csv", sep = ""))

	# Get indices based on the total_position and chr_percent.v
	sample_idx = sample(nrow(sig.df), round(total_position * chr_percent.v[i]), replace = FALSE)

	# Output the positions that are greater than 0
	sig.v = sig.df$count[sample_idx]

	# Write the output
	write(sig.v[which(sig.v > 0)], file = output_file, ncolumns = 1, append = T)

}

cat("\n\tComplete!\n")
