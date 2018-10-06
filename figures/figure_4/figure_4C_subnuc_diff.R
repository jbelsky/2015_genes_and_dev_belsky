#############################################
# Input file names

# Set the work dir
work_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_4/figure_4_datasets/"
dataset_work_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/"

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
footer = "signal_around_oridb_acs_feature_file_curated_798_sites_win_500bp.csv"









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
	mat.l[[i]] = average_matrices(convert_csv_to_matrix(paste(dataset_work_dir, file_name_type.v[2*i - 1], "_", footer, sep = "")),
				      convert_csv_to_matrix(paste(dataset_work_dir, file_name_type.v[2*i], "_", footer, sep = ""))
				     )

	# Subset on the G1 and G2 idx
	mat.l[[i]] = mat.l[[i]][oridb_subset_idx,]

}

# Get the idx
left = which(nuc_type.df$type == "left_movement")
right = which(nuc_type.df$type == "right_movement")
static = which(nuc_type.df$type == "static")

# Set the idx.l
idx.l = list(left, right, static)

##########################################
# Plotting

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

dens_split.s = split.screen(dens_split.m)

for(i in 1:3){

	# Open the screen
	screen(dens_split.s[i])	

	# Make the plot
	par(mar = c(0.5, 4.1, 0.5, 2.1), mgp = c(3, 0.5, 0))

	plot(0, 0, type = "n",
	     xlim = c(-500, 500), xaxs = "i", xaxt = "n",
	     ylim = c(0, 4.5), yaxs = "i",
	     ann = F
	    )
	axis(1, at = seq(-500, 500, 250), labels = F)
	text(x = -475, y = 3.8, labels = c("Upstream", "Downstream", "Static")[i], adj = 0)

	for(j in 1:3){
		lines(-500:500, apply(mat.l[[j]][idx.l[[i]],], 2, mean), 
		      col = c("red", "darkgreen", "purple")[j]
		     )
	}


}

# Y-axis Label
screen(dens_split.s[4])
par(mar = c(0, 4.1, 0, 2.1))
set_chromatin_schematic()
text(x = 0.5, y = 0.2, labels = "Footprint density by nucleosome shift")
par(xpd = T)
legend(x = c(-0.1, 1.1), y = c(1.25, 0.25),
       legend = c("G1 WT", expression(paste("G1 ", italic("cdc6-1"), sep = "")), "G2 WT"), lwd = 2, lty = 1,
       col = c("red", "purple", "darkgreen"), bty = "n", horiz = T, cex = 0.9
      )
par(xpd = F)

# X-axis Label
screen(dens_split.s[5])
par(mar = c(0, 4.1, 0, 2.1))
set_chromatin_schematic(x_start = -500, x_end = 500)
text(x = 0.5, y = 0.2, labels = "Relative distance from ACS (bp)")
axis(1, at = seq(-500, 500, 500), labels = T, tick = F, pos = 1.5)

# SubFigure Title Label
screen(dens_split.s[6])
par(mar = rep(0, 4))
set_chromatin_schematic()
text(x = 0.5, y = 0.5, labels = "Average footprint density", srt = 90)
