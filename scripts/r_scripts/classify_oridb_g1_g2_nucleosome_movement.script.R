# Clear the workspace
rm(list = ls())
graphics.off()

# Load the functions
source("/data/data2/jab112/2014_mnase_manuscript/datasets/r_scripts/make_dataset_heatmap_functions.R")

# Create the function
get_sample_from_density = function(nuc.v, rel_pos, win = 100, set_num = 500, sample_num_per_set = 100){

	# Set the upper and lower pos
	pos_low = rel_pos - win
	pos_high = rel_pos + win

	# Get the density range
	nuc_sig_range.v = nuc.v[as.character(pos_low:pos_high)]

	# Normalize the density to make a probability
	nuc_sig_range.v = nuc_sig_range.v / sum(nuc_sig_range.v)

	# Get the total number of samples from the distribution
	nuc_samp.v = sample(x = -win:win, size = set_num * sample_num_per_set, replace = T, prob = nuc_sig_range.v)	

	# Convert to a sample matrix
	nuc_samp.m = matrix(nuc_samp.v, nrow = set_num, ncol = sample_num_per_set, byrow = T)

	# Calculate the median position for each sample
	nuc_median.v = apply(nuc_samp.m, 1, median)

	# Return the nuc_median
	return(nuc_median.v)

}


# Set the work dir
work_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_4/figure_4_datasets/"

# Set the file name types
file_name_type.v = c("g2_dm242", "g2_dm356",
		     "g1_dm243", "g1_dm354"
		    )

# Set the file name footer
footer = "nuc_150_175_density_signal_around_oridb_acs_feature_file_curated_798_sites_win_500bp.csv"

# Set the output file name
output_file_name = "g1_g2_nuc_movement_around_oridb_acs_feature_file_curated_798_sites.csv"

# Set up the storage matrix
mat.l = vector("list", 2)
names(mat.l) = c("g2", "g1")

# Enter into the matrix list
for(i in 1:2){

	# Get the output density matrix
	mat.l[[i]] = average_matrices(convert_csv_to_matrix(paste(work_dir, file_name_type.v[2*i - 1], "_", footer, sep = "")),
				      convert_csv_to_matrix(paste(work_dir, file_name_type.v[2*i], "_", footer, sep = ""))
				     )

}

# Load the g2 nucleosome positions
g2_nuc_peaks.df = read.csv(paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_4/figure_4_datasets/",
				 "g2_nuc_left_right_rel_pos_around_oridb_acs_feature_file_curated_798_sites.csv", sep = ""
			        )
			  )

# Set the extra columns
g2_nuc_peaks.df = g2_nuc_peaks.df[,1:6]

# Set the quantiles
quantile.v = c(0.025, 0.05, 0.1, 0.25, 0.50, 0.75, 0.9, 0.95, 0.975)

# Set the column names
col_names = paste(rep(c("left", "right"), each = length(quantile.v)), quantile.v, sep = "")

# Initialize the columns
g2_nuc_peaks.df[,col_names] = 0

# Iterate through each nucleosome
for(i in 1:nrow(g2_nuc_peaks.df)){

	cat("Finding the nucleosome positions for ", g2_nuc_peaks.df$name[i], "...\r", sep = "")

	for(n in 1:2){

		# Set the type of nucleosome
		nuc_pos_type = c("left", "right")[n]
		nuc_output_type.l = list(left = paste("left", quantile.v, sep = ""),
					 right = paste("right", quantile.v, sep = "")			
					)

		# Get the G2 position
		g2_pos = g2_nuc_peaks.df[i,nuc_pos_type]

		if(!is.na(g2_pos)){

			# Get the median samples
			g2_sample_nuc_position.v = 
				get_sample_from_density(mat.l$g2[i,], g2_pos, 
						    	win = 100, set_num = 500, sample_num_per_set = 100
						       )

			g1_sample_nuc_position.v = 
				get_sample_from_density(mat.l$g1[i,], g2_pos, 
							win = 100, set_num = 500, sample_num_per_set = 100
						       )

			# Get the distribution of sample nucleosome differences
			g1_g2_diff_dist.v = g1_sample_nuc_position.v - g2_sample_nuc_position.v

			# Find the number of sites that are greater in the g1 sample
			g2_nuc_peaks.df[i,nuc_output_type.l[[n]]] = quantile(g1_g2_diff_dist.v, probs = quantile.v)

		}else{

			g2_nuc_peaks.df[i,nuc_output_type.l[[n]]] = NA

		}

	}

}

cat("\n")

g2_nuc_peaks.df = read.csv("/data/data2/jab112/2014_mnase_manuscript/figures/figure_4/figure_4_datasets/g1_g2_nuc_movement_around_oridb_acs_feature_file_curated_798_sites.csv")

# Set the movement threshold
rel_pos_thresh = 0

# Get the indices corresponding to each type
up_idx = which(g2_nuc_peaks.df$left0.9 < 0)
down_idx = which(g2_nuc_peaks.df$right0.1 > 0)

# Check if there are any idx in both groups
both_idx = up_idx[which(up_idx %in% down_idx)]

if(any(both_idx)){

	# Get the absolute value median for each group
	both_left_med.v = abs(g2_nuc_peaks.df$left0.5[both_idx])
	both_right_med.v = abs(g2_nuc_peaks.df$right0.5[both_idx])

	# Compare the medians
	remove_right = both_idx[which(both_left_med.v >= both_right_med.v)]
	remove_left = both_idx[which(both_left_med.v < both_right_med.v)]

	if(any(remove_left)){

		up_idx = up_idx[!up_idx %in% remove_left]

	}

	if(any(remove_right)){
		
		down_idx = down_idx[!down_idx %in% remove_right]
	
	}

}

static_idx = (1:nrow(g2_nuc_peaks.df))[-c(up_idx, down_idx)]

# Classify each type
g2_nuc_peaks.df[up_idx,"type"] = "left_movement"
g2_nuc_peaks.df[down_idx,"type"] = "right_movement"
g2_nuc_peaks.df[static_idx,"type"] = "static"

# Write the output
write.table(g2_nuc_peaks.df, file = paste(work_dir, output_file_name, sep = ""),
	    sep = ",", col.names = T, row.names = F, quote = F
	   )
# Load the functions
function_files.v = list.files(path = "/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions/", full.names = T)
for(i in 1:length(function_files.v)){
	source(function_files.v[i])
}
