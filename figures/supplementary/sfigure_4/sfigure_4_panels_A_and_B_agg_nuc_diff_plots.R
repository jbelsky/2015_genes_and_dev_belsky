#############################################
# Input file names

# Set the work dir
work_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_4/figure_4_datasets/"

# Get the total density feature file name
total_density_feature_file_name = paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/",
					"oridb_acs_feature_file_curated_798_sites_left_win_50bp_right_win_150bp.csv", sep = "")

# Set the feature file names
nuc_feature.fn = paste(work_dir, "g1_g2_nuc_movement_around_oridb_acs_feature_file_curated_798_sites.csv", sep = "")

# Oridb feature file name
oridb.fn = paste("/data/data2/jab112/2014_mnase_manuscript/datasets/jab112_yeast_feature_files/replication_origins/",
		 "oridb_acs_feature_file_curated_798_sites_timing_whitehouse_raw_oem_acs_seq.csv", sep = ""
		)


# Set the file name types
file_name_type.v = c("g1_dm243", "g1_dm354",
		     "g2_dm242", "g2_dm356",
		     "g1_cdc6_dm479", "g1_cdc6_dm481"
		    )

# Set the file name footer
footer = "nuc_150_175_density_signal_around_oridb_acs_feature_file_curated_798_sites_win_500bp.csv"









#############################################
# Data Processing


# Get the G1 and G2 subsets
oridb_idx.l = get_oridb_idx(total_density_feature_file_name)

# Get the oridb_subset
oridb_subset_idx = sort(unlist(oridb_idx.l[c("g1", "g2")]))

# Create the storage matrix list
mat.l = vector("list", 3)
names(mat.l) = c("g1", "g2", "cdc6")

# Set the feature file names
nuc_type.df = read.csv(nuc_feature.fn)
oridb.df = read.csv(oridb.fn)

# Subset on the oridb_subset_idx
nuc_type.df = nuc_type.df[oridb_subset_idx,]
oridb.df = oridb.df[oridb_subset_idx,]

# Create the dataset matrix
for(i in 1:3){

	# Get the output density matrix
	mat.l[[i]] = average_matrices(convert_csv_to_matrix(paste(work_dir, file_name_type.v[2*i - 1], "_", footer, sep = "")),
				      convert_csv_to_matrix(paste(work_dir, file_name_type.v[2*i], "_", footer, sep = ""))
				     )

	# Subset on the G1 and G2 idx
	mat.l[[i]] = mat.l[[i]][oridb_subset_idx,]

}

# Set the reference distribution
ref.v = as.vector(mat.l[[2]])

# Quantile normalize the first and third matrices
for(i in c(1,3)){

	# Get the matrix
	mat.m = mat.l[[i]]

	# Get the cdf function of the mat
	mat.cdf = ecdf(as.vector(mat.m))

	# Get the cdf of the mat
	mat.v = mat.cdf(as.vector(mat.m))

	# Quantile normalize
	mat_quant.v = quantile(ref.v, mat.v)

	# Remake the matrix
	mat.l[[i]] = matrix(mat_quant.v, ncol = 1001)

}

# Get the idx
left = which(nuc_type.df$type == "left_movement")
right = which(nuc_type.df$type == "right_movement")
static = which(nuc_type.df$type == "static")

# Set the idx.l
idx.l = list(left, right, static)

##########################################
# Plotting

# Split the overall screen
overall_plot_agg.m = matrix(c(0, 1, 0.5, 1,
			      0, 1, 0, 0.5),
			    ncol = 4, byrow = T
			   )

# Open the screen
overall_plot_agg.s = split.screen(overall_plot_agg.m)

# Split the bottom and top
dens_split.m = matrix(c(0, 1, 0.633, 0.9,
			0, 1, 0.367, 0.633,
			0, 1, 0.1, 0.367,
			0, 1, 0.9, 1,
			0, 1, 0, 0.1,
			0.05, 0.1, 0.1, 0.9
		       ),
		      ncol = 4, byrow = T
		     )

# Set the plot index
plot_idx.l = list(c(1, 2), c(3, 2))

for(a in 1:2){

	# Open the screen
	screen(overall_plot_agg.s[a])

	# Split the current screen
	dens_split.s = split.screen(dens_split.m)

	for(i in 1:3){

		# Open the screen
		screen(dens_split.s[i])	

		# Make the plot
		par(mar = c(0.5, 4.1, 0.5, 2.1), mgp = c(3, 0.5, 0))

		plot(0, 0, type = "n",
		     xlim = c(-500, 500), xaxs = "i", xaxt = "n",
		     ylim = c(0, 2.5), yaxs = "i", yaxt = "n",
		     ann = F
		    )
		axis(1, at = seq(-500, 500, 250), labels = F)
		axis(2, at = 0:2)

		text(x = -475, y = 2, labels = c("Upstream", "Downstream", "Static")[i], adj = 0)

		for(j in plot_idx.l[[a]]){
			lines(-500:500, apply(mat.l[[j]][idx.l[[i]],], 2, mean), 
			      col = c("red", "darkgreen", "red")[j]
			     )
		}


	}

	# Y-axis Label
	screen(dens_split.s[4])
	par(mar = c(0, 4.1, 0, 2.1))
	set_chromatin_schematic()
	legend("left", legend = c("G1 WT", expression(paste("G1 ", italic("cdc6-1"))))[a], col = "red", lwd = 2, lty = 1, bty = "n", xjust = 0)
	legend("right", legend = "G2 WT", col = "darkgreen", lwd = 2, lty = 1, bty = "n", xjust = 1)

	# X-axis Label
	screen(dens_split.s[5])
	par(mar = c(0, 4.1, 0, 2.1))
	set_chromatin_schematic(x_start = -500, x_end = 500)
	text(x = 0.5, y = 0.4, labels = "Relative distance from ACS (bp)")
	axis(1, at = seq(-500, 500, 500), labels = T, tick = F, pos = 1.5)

	# Y-axis Label
	screen(dens_split.s[6])
	par(mar = rep(0, 4))
	set_chromatin_schematic()
	text(x = 0.5, y = 0.5, labels = "Average nucleosome density", srt = 90)

}
