# Clear the workspace
graphics.off()
rm(list = ls())

#############################################
# Load the functions
function_files.v = list.files(path = "/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions/", full.names = T)
for(i in 1:length(function_files.v)){
	source(function_files.v[i])
}



#############################################
# Set the file names

# Set the total density feature file name
total_density_feature_file_name = 
	paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/",
	      "oridb_acs_feature_file_curated_798_sites_left_win_50bp_right_win_150bp.csv", sep = ""
	     )

# Set the file types
dm_file.v = c("g2_dm242", "g2_dm356", "g2_orc1_161_dm261", "g2_orc1_161_dm334")

# Set the dir
dir = "/data/data2/jab112/2014_mnase_manuscript/figures/supplementary/sfigure_1/"
setwd(dir)

#############################################
# Plotting


# Set up the plot
png(file = paste("/data/data2/jab112/2014_mnase_manuscript/figures/genes_and_dev_figures/",
		 "sfig_orcmtCtrl.png", sep = ""),
    width = 8.5, height = 8.5, units = "in", res = 300
   )


# Split the screen into fourths
scr.m = matrix(c(0, 0.5, 0.5, 1,
		 0.5, 1, 0.5, 1,
		 0.5, 1, 0, 0.5,
		 0, 0.05, 0.925, 0.975,
		 0.5, 0.55, 0.925, 0.975,
		 0.5, 0.55, 0.425, 0.475
		),
	       ncol = 4, byrow = T
	      )

# Split the screen
close.screen(all.screens = T)
scr.s = split.screen(scr.m)



####################################
# TSS Nucleosome Plot

# Move to screen 1
screen(scr.s[1])

# Make the nucleosome plot
nuc.l = vector("list", 4)

# Get the nucleosome average plot
for(i in 1:4){

	# Load the file
	nuc.df = read.csv(paste("xu_tss_nuc/", dm_file.v[i], "_xu_tss_nuc.csv", sep = ""))

	# Convert to matrix
	nuc.m = as.matrix(nuc.df[,-(1:4)])	

	# Make the average plot
	nuc.l[[i]] = apply(nuc.m, 2, mean)

}

# Remove the nuc.df and nuc.m
rm(nuc.df, nuc.m)

# Get the G2 and G2 orc1-161 averages
g2.v = (nuc.l[[1]] + nuc.l[[2]]) / 2
g2_orc1.v = (nuc.l[[3]] + nuc.l[[4]]) / 2

par(cex.axis = 1.25, cex.lab = 1.25)

# Make the plot
plot(-500:500, g2.v, type = "l", col = "red", lwd = 2,
     xlim = c(-500, 500), xaxs = "i", xaxt = "n",
     ylim = c(0, 3),
     xlab = "Relative distance from TSS (bp)",
     ylab = "Average nucleosome density"
    )
lines(-500:500, g2_orc1.v, col = "darkgreen", lwd = 2)
axis(1, at = seq(-500, 500, 250))

# Add in the legend
legend("topleft", legend = "WT", lwd = 3, col = "red", bty = "n", cex = 1.25)
legend("topright", legend = expression(italic("orc1-161")), lwd = 3, col = "darkgreen", bty = "n", cex = 1.25)



####################################
# Abf1 Footprint Plot

# Move to screen 2
screen(scr.s[2])

# Add in the Abf1 subnuc
subnuc.l = vector("list", 4)

# Get the nucleosome average plot
for(i in 1:4){

	# Load the file
	subnuc.df = read.csv(paste("abf1_subnuc/", dm_file.v[i], "_signal_around_abf1_feature_file_win_500bp.csv", sep = ""))

	# Convert to matrix
	subnuc.m = as.matrix(subnuc.df[,-(1:4)])	

	# Make the average plot
	subnuc.l[[i]] = apply(subnuc.m, 2, mean)

}

# Get the G2 and G2 orc1-161 averages
g2.v = (subnuc.l[[1]] + subnuc.l[[2]]) / 2
g2_orc1.v = (subnuc.l[[3]] + subnuc.l[[4]]) / 2

# Make the plot
plot(-500:500, g2.v, type = "l", col = "red", lwd = 2,
     xlim = c(-500, 500), xaxs = "i", xaxt = "n",
     ylim = c(0, 9),
     xlab = "Relative distance from Abf1 (bp)",
     ylab = "Average footprint density"
    )
lines(-500:500, g2_orc1.v, col = "darkgreen", lwd = 2)
axis(1, at = seq(-500, 500, 250))

# Add in the legend
legend("topleft", legend = "WT", lwd = 3, col = "red", bty = "n", cex = 1.25)
legend("topright", legend = expression(italic("orc1-161")), lwd = 3, col = "darkgreen", bty = "n", cex = 1.25)





####################################
# Origin Footprint Plot

# Move to screen 3
screen(scr.s[3])

# Get the oridb subset
oridb_idx.l = get_oridb_idx(total_density_feature_file_name)

# Get the indices
g2_idx = oridb_idx.l$g2

# Create the storage matrix list
mat.l = vector("list", 2)
names(mat.l) = c("g2", "g2_orc1_161")

# Set the oridb footprint dir
oridb_footprint_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/"

# Set the footer
footer = "signal_around_oridb_acs_feature_file_curated_798_sites_win_500bp.csv"

# Create the dataset matrix
for(i in 1:2){

	# Get the output density matrix
	mat.m = average_matrices(
		convert_csv_to_matrix(paste(oridb_footprint_dir, dm_file.v[2*i-1], "_", footer, sep = "")),
		convert_csv_to_matrix(paste(oridb_footprint_dir, dm_file.v[2*i], "_", footer, sep = ""))
		)

	# Get the average signal over the G2 footprint
	mat.l[[i]] = apply(mat.m[g2_idx,], 2, mean)

}

# Make the aggregate plot over the ACS
plot(-500:500, mat.l$g2, type = "l", col = "red", lwd = 2,
     xlim = c(-500, 500), xaxs = "i", xaxt = "n", ylim = c(0, 3.5),
     xlab = "Relative distance from ACS (bp)",
     ylab = "Average footprint density"

    )
lines(-500:500, mat.l$g2_orc1_161, col = "darkgreen", lwd = 2)
axis(1, at = seq(-500, 500, 250))

# Add in the legend
legend("topleft", legend = "WT", lwd = 3, col = "red", bty = "n", cex = 1.25)
legend("topright", legend = expression(italic("orc1-161")), lwd = 3, col = "darkgreen", bty = "n", cex = 1.25)







#####################################
# Add the figure labels

# Add in the figure labels
screen(scr.s[4])
add_figure_label("A")

screen(scr.s[5])
add_figure_label("B")

screen(scr.s[6])
add_figure_label("C")

# Close the screens
close.screen(all.screens = T)

# Close the device
dev.off()
