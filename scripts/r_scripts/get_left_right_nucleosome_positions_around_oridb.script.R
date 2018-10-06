# Clear the workspace
rm(list = ls())
graphics.off()



# Open the libraries
library(GenomicRanges)
library(nucleR)

#############################################
# Load the functions
function_files.v = list.files(path = "/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions/", full.names = T)
for(i in 1:length(function_files.v)){
	source(function_files.v[i])
}










#############################################
# Data Processing

# Set the parameters
min_thresh = 0.25
peak_width = 150

# Set the work dir
work_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_4/figure_4_datasets/"

# Set the file name types
file_name_type.v = c("g1_dm243", "g1_dm354")

# Set the file name footer
footer = "nuc_150_175_density_signal_around_oridb_acs_feature_file_curated_798_sites_win_500bp.csv"

# Load the feature file
feature_file_name = paste("/data/data2/jab112/2014_mnase_manuscript/datasets/jab112_yeast_feature_files/",
			  "replication_origins/oridb_acs_feature_file_curated_798_sites_timing_oem.csv", sep = ""
			 )

# Write the output name
output_file_name = "g1_nuc_left_right_rel_pos_around_oridb_acs_feature_file_curated_798_sites.csv"

# Get the output density matrix
nuc.m = average_matrices(convert_csv_to_matrix(paste(work_dir, file_name_type.v[1], "_", footer, sep = "")),
		         convert_csv_to_matrix(paste(work_dir, file_name_type.v[2], "_", footer, sep = ""))
			)

# Load the feature file
oridb.df = read.csv(feature_file_name)

# Set the storage
nuc_peaks.df = data.frame(oridb.df[,1:4],
			  left = numeric(nrow(nuc.m)),
			  right = numeric(nrow(nuc.m))
			 )

# Set the peak window
peak_win = (peak_width - 1) / 2

for(j in 1:nrow(nuc.m)){

	# Set the cov.v
	cov.v = nuc.m[j,]

	# Get the peaks
	nuc.df = get_mod_peaks(cov.v, 0, min_thresh = min_thresh, peak_width = peak_width)

	# Enter the nucleosome position into the output
	nuc_peaks.df[j,c("left", "right")] = c(max(get_nuc_peak(nuc.df, -91)), 
					       min(get_nuc_peak(nuc.df, 148))
					      )

}

# Write the output
write.table(nuc_peaks.df, file = paste(work_dir, output_file_name, sep = ""),
	    sep = ",", col.names = T, row.names = F, quote = F
	   )
