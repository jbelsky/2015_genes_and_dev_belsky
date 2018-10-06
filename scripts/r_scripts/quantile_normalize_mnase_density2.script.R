# Clear the workspace
rm(list = ls())
graphics.off()

# Set the working directory
setwd("/data/data2/jab112/2014_mnase_manuscript/datasets/nuc_data/")

# Set the header
header.v = c("g1_nuc_cdc6_dm466",	
	     "g1_nuc_cdc6_dm467",	
	     "g1_nuc_wt_37_dm468",	
	     "g1_nuc_wt_37_dm469"
	    )

# Set the reference file name
ref_dist = "g2_dm242"

# Iterate through each sample
for(h in header.v){

	cat("Quantile-normalizing ", h, "...\n", sep = "")

	# Set the sample distribution
	reference_input_header = 
		paste("raw_data2/", ref_dist, 
		      "_nuc_150_175_density_signal_chr_", sep = "")

	# Set the sample header
	sample_input_header = paste("raw_data2/", h, "_nuc_150_175_density_signal_chr_", sep = "")

	# Set the sample output header
	sample_output_header = paste("quant_normalize2/", h, "_nuc_150_175_density_signal_chr_", sep = "")

	##########################################################################

	# Iterate through each chr
	for(c in 1:16){

		cat("\tChr ", c, "...\r", sep = "")

		# Load the reference distribution
		ref.df = read.csv(paste(reference_input_header, c, ".csv", sep = ""))

		# Load the signal of the sample to quantile normalize
		sample.df = read.csv(paste(sample_input_header, c, ".csv", sep = ""))
		
		# Get the cdf of the sample distribution
		sample.cdf = ecdf(sample.df$signal)

		# Get the quantile distribution of the raw_signal
		quant.v = sample.cdf(sample.df$signal)

		# Quantile-Normalize
		chr_norm_count.v = quantile(ref.df$signal, probs = quant.v, names = F)

		# Update the data frame
		sample.df$signal = chr_norm_count.v

		# Write the output table
		write.table(sample.df,
			    file = paste(sample_output_header, c, "_quant_normalize.csv", sep = ""),
			    sep = ",", col.names = T, row.names = F, quote = F
			   )

	}

	cat("\n\tComplete!\n")

}
