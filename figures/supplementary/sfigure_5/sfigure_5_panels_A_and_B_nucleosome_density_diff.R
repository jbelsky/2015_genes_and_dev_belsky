#############################################
# Set the Filenames
feature.fn = paste("/data/data2/jab112/2014_mnase_manuscript/datasets/jab112_yeast_feature_files/replication_origins/",
		   "oridb_acs_feature_file_curated_798_sites_timing_whitehouse_raw_oem_acs_seq.csv", sep = ""
		  )

oridb_subnuc_density.fn = paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/",
				"oridb_acs_feature_file_curated_798_sites_left_win_50bp_right_win_150bp.csv", 
				sep = ""
			       )


# Get the mcm single density
mcm_datasets.fn.v = paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/",
			  "mcm_chip_seq_dm", c(82, 282), "_signal_around_oridb_acs_feature_file_curated_798_sites_win_500bp.csv",
			  sep = ""
			 )

# Set the work dir
work_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_4/figure_4_datasets/"

# Set the file name types
file_name_type.v = c("g1_dm243", "g1_dm354",
		     "g2_dm242", "g2_dm356"
		    )

# Set the file name footer
footer = "nuc_150_175_density_signal_around_oridb_acs_feature_file_curated_798_sites_win_500bp.csv"

# Load the peak file
load("/data/data2/jab112/2014_mnase_manuscript/figures/figure_5/figure_5_datasets/orc_and_mcm_peak_mat.l")

# Set the x_win
x_win = 400

#############################################
# Data Processing

# Load the feature file name
feature.df = read.csv(feature.fn)

# Get the plot idx
oridb_idx.l = get_oridb_idx(oridb_subnuc_density.fn)

# Get the idx 
oridb_subset_idx = sort(unlist(oridb_idx.l[c("g1", "g2")]))

# Subset on the oridb_subset_idx
feature.df = feature.df[oridb_subset_idx,]

# Get the mcm density
mcm.m = average_matrices(convert_csv_to_matrix(mcm_datasets.fn.v[1]),
			 convert_csv_to_matrix(mcm_datasets.fn.v[2])
			)[oridb_subset_idx,]

# Get the Mcm sum in a 100 bp window around the left (-90) and right (160) nucleosomes
left_mcm_sig.v = apply(mcm.m[,as.character(-190:10)], 1, sum)
right_mcm_sig.v = apply(mcm.m[,as.character(60:260)], 1, sum)

# Get the left and right mcm indices
left_idx = which(left_mcm_sig.v > right_mcm_sig.v)
right_idx = which(left_mcm_sig.v < right_mcm_sig.v)










#############################################
# Nucleosome Data Processing

# Create the storage matrix list
mat.l = vector("list", 2)
names(mat.l) = c("g1", "g2")

# Create the dataset matrix
for(i in 1:2){

	# Get the output density matrix
	mat.l[[i]] = average_matrices(convert_csv_to_matrix(paste(work_dir, file_name_type.v[2*i - 1], "_", footer, sep = "")),
				      convert_csv_to_matrix(paste(work_dir, file_name_type.v[2*i], "_", footer, sep = ""))
				     )

	# Subset on the G1 and G2 idx
	mat.l[[i]] = mat.l[[i]][oridb_subset_idx,]

}

# Plotting
nuc_scr.m = matrix(c(0, 0.5, 0.7, 0.9,
		     0, 0.5, 0, 0.7,
		     0.5, 1, 0.7, 0.9,
		     0.5, 1, 0, 0.7
		    ), ncol = 4, byrow = T
		  )

# Split the screen
nuc_scr.s = split.screen(nuc_scr.m)

# Open the screen
screen(nuc_scr.s[1])

# Set the mar
par(mar = c(0, 4.1, 0, 2.1))

# Set the plot area
set_chromatin_schematic(x_start = -400, x_end = 400, y_start = 0, y_end = 2)

# Get the nucleosome peaks
nuc_peaks.df = get_mod_peaks(apply(mat.l[[1]][left_idx,], 2, mean), x_mid = 0, peak_width = 150,
			     min_thresh = 0.5
			    )

# Make the nucleosome
plot_nucleosome(nuc_peaks.df[2:3,], y_max = 1.5, yh = 0.25, y0 = 1.5)
plot_nucleosome(nuc_peaks.df[2:3,], y_max = 1.5, yh = 0.25, y0 = 0.5)

# Add the Mcm to the G1
mcm_pos.v = nuc_peaks.df$pos[2] + 65 + c(17, 51)
plot_mcm(mcm_pos.v[1], x_w = 17, y0 = 1.5, yh = 0.27, obj_col = "purple")
plot_mcm(mcm_pos.v[2], x_w = 17, y0 = 1.5, yh = 0.27, obj_col = "purple")

# Add the text
text(x = -250, y = c(1.5, 0.5), labels = c("G1", "G2"), cex = 1.25)

# Open the screen
screen(nuc_scr.s[2])
par(mar = c(5.1, 4.1, 0.85, 2.1))

# Make the nucleosome density plot
plot(0, 0, type = "n",
     xlim = c(-400, 400), xaxs = "i", xaxt = "n",
     ylim = c(0, 2), yaxt = "n",
     xlab = "Relative distance from ACS (bp)",
     ylab = "Average nucleosome density"
    )
axis(1, at = seq(-400, 400, 200))
axis(2, at = 0:2)
axis(2, at = c(0.5, 1.5), labels = F)

# Set the legend
legend("topleft", legend = "G1", lwd = 2, col = "red", bty = "n")
legend("topright", legend = "G2", lwd = 2, col = "darkgreen", bty = "n")

lines(-500:500, apply(mat.l[[1]][left_idx,], 2, mean), col = "red")
lines(-500:500, apply(mat.l[[2]][left_idx,], 2, mean), col = "darkgreen")





# Open the screen
screen(nuc_scr.s[3])

# Set the mar
par(mar = c(0, 4.1, 0, 2.1))

# Set the plot area
set_chromatin_schematic(x_start = -400, x_end = 400, y_start = 0, y_end = 2)

# Get the nucleosome peaks
nuc_peaks.df = get_mod_peaks(apply(mat.l[[1]][right_idx,], 2, mean), x_mid = 0, peak_width = 150,
			     min_thresh = 0.5
			    )

# Make the nucleosome
plot_nucleosome(nuc_peaks.df[2:3,], y_max = 1.5, yh = 0.25, y0 = 1.5)
plot_nucleosome(nuc_peaks.df[2:3,], y_max = 1.5, yh = 0.25, y0 = 0.5)

# Add the Mcm to the G1
mcm_pos.v = nuc_peaks.df$pos[3] - 65 - c(17, 51)
plot_mcm(mcm_pos.v[1], x_w = 17, y0 = 1.5, yh = 0.27, obj_col = "purple")
plot_mcm(mcm_pos.v[2], x_w = 17, y0 = 1.5, yh = 0.27, obj_col = "purple")

# Add the text
text(x = -250, y = c(1.5, 0.5), labels = c("G1", "G2"), cex = 1.25)

# Open the screen
screen(nuc_scr.s[4])
par(mar = c(5.1, 4.1, 0.85, 2.1))

# Make the nucleosome density plot
plot(0, 0, type = "n",
     xlim = c(-400, 400), xaxs = "i", xaxt = "n",
     ylim = c(0, 2), yaxt = "n",
     xlab = "Relative distance from ACS (bp)",
     ylab = "Average nucleosome density"
    )
axis(1, at = seq(-400, 400, 200))
axis(2, at = 0:2)
axis(2, at = c(0.5, 1.5), labels = F)

# Set the legend
legend("topleft", legend = "G1", lwd = 2, col = "red", bty = "n")
legend("topright", legend = "G2", lwd = 2, col = "darkgreen", bty = "n")

lines(-500:500, apply(mat.l[[1]][right_idx,], 2, mean), col = "red")
lines(-500:500, apply(mat.l[[2]][right_idx,], 2, mean), col = "darkgreen")
