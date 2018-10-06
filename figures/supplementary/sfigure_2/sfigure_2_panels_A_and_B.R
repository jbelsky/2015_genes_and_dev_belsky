##################################################
# Set the file parameters

# Set the file names
file_name_type.v = c("g1_dm243", "g1_dm354",
		     "g2_dm242", "g2_dm356",
		     "g1_cdc6_dm479", "g1_cdc6_dm481"
		    )

# Set teh total density signal file
total_density_feature.fn = 
	paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/",
	      "oridb_acs_feature_file_curated_798_sites_left_win_50bp_right_win_150bp.csv",
	      sep = ""
	     )

# Set the directory
oridb_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/"
abf1_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/supplementary/sfigure_1/abf1_subnuc/"

# Set the footer
oridb_footer = "signal_around_oridb_acs_feature_file_curated_798_sites_win_500bp.csv"
abf1_footer = "signal_around_abf1_feature_file_win_500bp.csv"


#######################################
# Data Processing

# Get the plot idx
oridb_idx.l = get_oridb_idx(total_density_feature.fn)

# Create the storage matrix list
oridb_mat.l = vector("list", 3)
names(oridb_mat.l) = c("g1", "g2", "g1_cdc6")

# Get the average signal over each subset
for(i in 1:3){

	# Create the average signal
	mat.m = average_matrices(
		convert_csv_to_matrix(paste(oridb_dir, file_name_type.v[2*i-1], "_", oridb_footer, sep = "")),
		convert_csv_to_matrix(paste(oridb_dir, file_name_type.v[2*i], "_", oridb_footer, sep = ""))
		)

	# Get the average signal over each footprint
	oridb_mat.l[[i]] = apply(mat.m[oridb_idx.l$g2,], 2, mean)

}

# Create the storage matrix list
abf1_mat.l = vector("list", 3)
names(abf1_mat.l) = c("g1", "g2", "g1_cdc6")

# Get the average signal over each subset
for(i in 1:3){

	# Create the average signal
	mat.m = average_matrices(
		convert_csv_to_matrix(paste(abf1_dir, file_name_type.v[2*i-1], "_", abf1_footer, sep = "")),
		convert_csv_to_matrix(paste(abf1_dir, file_name_type.v[2*i], "_", abf1_footer, sep = ""))
		)

	# Get the average signal over each footprint
	abf1_mat.l[[i]] = apply(mat.m, 2, mean)

}






#######################################
# Set the specific functions
make_base_plot = function(x_lab, y_max){

	# Make the aggregate plot over the ACS
	plot(0, 0, type = "n",
	     xlim = c(-250, 250), xaxs = "i", xaxt = "n",
	     ylim = c(0, y_max),
	     xlab = x_lab,
	     ylab = "Average footprint density"
	    )
	axis(1, at = seq(-250, 250, 250))
	axis(1, at = c(-125, 125), labels = F)

	legend("topleft", legend = c("G1 WT", expression(paste("G1 ", italic("cdc6-1"), sep = ""), "G2")), lwd = 3, 
	       col = c("red", "purple", "darkgreen"), bty = "n"
	      )
	# legend("topright", legend = "G2", lwd = 3, col = "darkgreen", bty = "n", cex = 1.25)

}


#############################################
# Plotting

# Set up the plot
par(cex.axis = 1, cex.lab = 1, mar = c(4.1, 4.1, 1, 2.1), mgp = c(2.75, 1, 0))

# Split the screen in half
fp_scr.m = matrix(c(0, 0.5, 0, 1,
		    0.5, 1, 0, 1
		   ),
		   ncol = 4, byrow = T
		  )

# Split the screen
fp_scr.s = split.screen(fp_scr.m)

############################################
# Make the oridb plot

# Move to screen 1
screen(fp_scr.s[1])

# Make the aggregate plot over the ACS
make_base_plot("Relative distance from ACS (bp)", 4)

# Add in the lines
lines(-500:500, oridb_mat.l$g1, col = "red", lwd = 2)
lines(-500:500, oridb_mat.l$g2, col = "darkgreen", lwd = 2)
lines(-500:500, oridb_mat.l$g1_cdc6, col = "purple", lwd = 2)




###################################
# Make the abf1 plot

# Move to screen 2
screen(fp_scr.s[2])

# Make the aggregate plot over the ACS
make_base_plot("Relative distance from Abf1 (bp)", 8.5)

# Add in the lines
lines(-500:500, abf1_mat.l$g1, col = "red", lwd = 2)
lines(-500:500, abf1_mat.l$g2, col = "darkgreen", lwd = 2)
lines(-500:500, abf1_mat.l$g1_cdc6, col = "purple", lwd = 2)
