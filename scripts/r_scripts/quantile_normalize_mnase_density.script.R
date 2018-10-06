# Clear the workspace
rm(list = ls())
graphics.off()

# Set the working directory
setwd("/data/data2/jab112/2014_mnase_manuscript/datasets/chip_data/orc/")

# Set the header
header.v = c("input_chip_seq_dm272"
	    )

# Set the reference file name
ref_dist = "signal_distribution/dm265_dist_output.csv"

# Load the reference distribution
ref_sample_dist.v = scan(ref_dist)

# Iterate through each sample
for(h in header.v){

	cat("Quantile-normalizing ", h, "...\n", sep = "")

	# Set the sample distribution
	sample_dist = paste("/data/data2/jab112/2014_mnase_manuscript/datasets/chip_data/input/signal_distribution/",
			    "input_chip_seq_dm272_dist_output.csv", sep = "")

	# Set the sample header
	sample_input_header = paste("/data/data2/jab112/2014_mnase_manuscript/datasets/chip_data/input/raw_data/",
				    h, "_density_signal_chr_", sep = "")

	# Set the sample output header
	sample_output_header = paste("quant_normalize/", h, "_orc_density_signal_chr_", sep = "")

	##########################################################################

	# Load the distributions
	sample_to_quant_dist.v = scan(sample_dist)

	# Get the cdf of the reference and sample distribution
	sample.cdf = ecdf(sample_to_quant_dist.v)

	# Iterate through each chr
	for(c in 1:16){

		cat("\tChr ", c, "...\r", sep = "")

		# Load the signal of the sample to quantile normalize
		raw_sig.df = read.csv(paste(sample_input_header, c, ".csv", sep = ""))

		# Get the quantile distribution of the raw_signal
		quant.v = sample.cdf(raw_sig.df$signal)

		# Quantile-Normalize
		chr_norm_count.v = quantile(ref_sample_dist.v, probs = quant.v, names = F)

		# Update the data frame
		raw_sig.df$signal = chr_norm_count.v

		# Write the output table
		write.table(raw_sig.df,
			    file = paste(sample_output_header, c, "_quant_normalize.csv", sep = ""),
			    sep = ",", col.names = T, row.names = F, quote = F
			   )

	}

	cat("\n\tComplete!\n")

}
