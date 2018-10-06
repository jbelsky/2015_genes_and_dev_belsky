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

# Set the plot idx
plot_idx = c(static, right, left)

# Get the lengths of each data set
idx_length.v = c(length(static), length(right), length(left))

# Set the y-axis
y_axis = c(idx_length.v[1]/2, idx_length.v[1] + (idx_length.v[2]/2),
	   sum(idx_length.v[1:2]) + (idx_length.v[3]/2)
	  )

# Create the difference lists
diff.l = vector("list", 2)

# Create the difference matrix
diff.l[[1]] = log2(mat.l$g1 + 1) - log2(mat.l$g2 + 1)
diff.l[[2]] = log2(mat.l$cdc6 + 1) - log2(mat.l$g2 + 1)




##########################################
# Plotting

# Split the bottom and top
dens_split.m = matrix(c(0, 1, 0.5, 1,
			0, 1, 0, 0.5),
		      ncol = 4, byrow = T
		     )

dens_split.s = split.screen(dens_split.m)

for(i in 1:2){

	# Open the screen
	screen(dens_split.s[i])	

	# Make the plot
	par(mar = c(4.1, 4.1, 2, 2.1))
	dens_dot_plot(diff.l[[i]][plot_idx,], z_min = -2, z_max = 2,
		      lowCol = "darkgreen", medCol = "white", highCol = "red",
		      plot_title = "",
		      x_label = "", y_label = "",
		      x_axt = "n", y_axt = "n")
	if(i == 1){
		title(main = "G1 WT", adj = 0, col.main = "red", font.main = 1, cex.main = 1)
	}else{
		title(main = "G1", adj = 0, col.main = "red", line = 1.25, cex.main = 1, font.main = 1)
		title(main = "cdc6-1", font.main = 3, adj = 0, col.main = "red", cex.main = 1, line = 0.25)
	}
# expression(paste("G1\n", italic("cdc6-1"), sep = "")))[i], adj = 0, col.main = "red", font.main = 1, cex.main = 1)
	title(main = "G2 WT", adj = 1, col.main = "darkgreen", font.main = 1, cex.main = 1)
	title(main = expression(paste(log[2], " nucleosome ratio")), font.main = 1, cex.main = 0.95, line = 0.9)
	title(xlab = "Relative distance from ACS (bp)", line = 2.5)
	title(ylab = "OriDB origin by nucleosome shift", line = 2.5)
	axis(1, at = c(-500.5, -250, 0, 250, 500.5), labels = c(-500, "", 0, "", 500))

	axis(4, at = y_axis, tick = F, pos = 440,
		labels = idx_length.v
	    )
	axis(2, at = y_axis, labels = c("Static", "Downstream", "Upstream"), tick = T, cex.axis = 0.9)

	# Set the dividing line
	abline(h = c(idx_length.v[1], sum(idx_length.v[1:2])), col = "black", lty = 2, lwd = 2)

}
