#############################################
# Set the Filenames
feature.fn = paste("/data/data2/jab112/2014_mnase_manuscript/datasets/jab112_yeast_feature_files/replication_origins/",
		   "oridb_acs_feature_file_curated_798_sites_timing_whitehouse_raw_oem_acs_seq.csv", sep = ""
		  )

oridb_subnuc_density.fn = paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/",
				"oridb_acs_feature_file_curated_798_sites_left_win_50bp_right_win_150bp.csv", 
				sep = ""
			       )

g1_agg_typhoon_plot.fn = paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_5/figure_5_datasets/",
			       "g1_and_g2_subnuc_peak_acs_398_subset_aggregate_typhoon_plot_g1_dm243_",
			       "win_400_high_frag_250.csv",
			       sep = ""
			      )

# Set the work dir
work.dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/"

# Set the file name types
datasets.fn.v = c("orc_chip_seq_dm265", "orc_chip_seq_dm287",
		  "mcm_chip_seq_dm82", "mcm_chip_seq_dm282"
		 )

# Set the file name footer
footer = "signal_around_oridb_acs_feature_file_curated_798_sites_win_500bp.csv"

# Load the peak_mat.l
load("/data/data2/jab112/2014_mnase_manuscript/figures/figure_5/figure_5_datasets/orc_and_mcm_peak_mat.l")


#############################################
# Set the Parameters
orc_thresh = 0.8
mcm_thresh = 0.5




################################################################################################################
# R Functions

# Load the functions
function_files.v = list.files(path = "/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions/", full.names = T)
for(i in 1:length(function_files.v)){
	source(function_files.v[i])
}






#############################################
# Data Processing

# Load the feature file name
feature.df = read.csv(feature.fn)

if(0){
# Create the output peak list
peak_mat.l = vector("list", 2)
names(peak_mat.l) = c("orc", "mcm")

# Create the dataset matrix
for(i in 1:length(peak_mat.l)){

	# Get the output density matrix
	mat.m = average_matrices(convert_csv_to_matrix(paste(work.dir, datasets.fn.v[2*i - 1], "_", footer, sep = "")),
				 convert_csv_to_matrix(paste(work.dir, datasets.fn.v[2*i], "_", footer, sep = ""))
				)

	# Enter in the row names
	rownames(mat.m) = feature.df$name

	# Get the output peak file name
	peak_mat.l[[i]] = get_chip_peak_over_mat(mat.m, threshold_peak = c(orc_thresh, mcm_thresh)[i])

}
}

# Get the plot idx
oridb_idx.l = get_oridb_idx(oridb_subnuc_density.fn)

# Get the idx 
idx = sort(unlist(oridb_idx.l[c("g1", "g2")]))





# Load the typhoon plot
typhoon_plot_mat.m = convert_typhoon_plot_to_matrix(g1_agg_typhoon_plot.fn)














###################################
# Plotting

# Set up the par
par(mar = c(3.1, 3.1, 2, 1), mgp = c(2, 0.75, 0))

# Make the dens_dot_plot
dens_dot_plot(typhoon_plot_mat.m,
	      z_max = 700, lowCol = "white", highCol = "#6B6B6B",
	      x_axt = "n", y_axt = "n",
	      x_label = "Relative distance from ACS (bp)",
	      y_label = "Fragment length (bp)"
	     )
axis(1, at = c(-400.5, -200, 0, 200, 400.5), labels = seq(-400, 400, 200))
axis(2, at = c(0.5, 50, 100, 150, 200, 250.5), labels = seq(0, 250, 50))

# Set the y_axis_midpoints
y_axis_mid.v = c(150, 210)

# Enter the points
for(i in 1:2){

	# Get the x_axis points
	chip_peaks.v = peak_mat.l[[i]][idx,"rel_pos"]
	chip_peaks.v = chip_peaks.v[which(!is.na(chip_peaks.v))]

	# Create the y_axis_scatter
	y = rnorm(length(chip_peaks.v), mean = y_axis_mid.v[i], sd = 10)
	y[which(y > y_axis_mid.v[i] + 25)] = y_axis_mid.v[i] + 25
	y[which(y < y_axis_mid.v[i] - 25)] = y_axis_mid.v[i] - 25

	# Enter in the points
	points(x = chip_peaks.v, y = y, col = c("#00640095", "#6A287E95")[i], pch = 19, cex = 0.5)

}

# Add in the legend
legend("bottomright", legend = c("Mcm2-7 peak", "ORC peak"), 
		      col = c("#6A287E", "#006400"), pch = 19, bty = "n", horiz = F, cex = 1)
