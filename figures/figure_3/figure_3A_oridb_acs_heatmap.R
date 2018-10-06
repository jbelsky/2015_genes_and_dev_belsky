# Set the dir
dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/"

#############################################
# Data Processing

# Set the work dir
work_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/"

# Set the file name types
file_name_type.v = c("g1_dm243", "g1_dm354",
		     "g2_dm242", "g2_dm356"
		    )

# Set the file name footer
footer = "signal_around_oridb_acs_feature_file_curated_798_sites_win_500bp.csv"

# Set the total density signal file
total_density_feature_file_name = 
	paste(work_dir, "oridb_acs_feature_file_curated_798_sites_left_win_50bp_right_win_150bp.csv", sep = "")

# Create the storage matrix list
mat.l = vector("list", 3)
names(mat.l) = c("g1", "g2")

# Create the dataset matrix
for(i in 1:2){

	# Get the output density matrix
	mat.l[[i]] = average_matrices(convert_csv_to_matrix(paste(work_dir, file_name_type.v[2*i - 1], "_", footer, sep = "")),
				      convert_csv_to_matrix(paste(work_dir, file_name_type.v[2*i], "_", footer, sep = ""))
				     )

}

# Get the plot idx
oridb_idx.l = get_oridb_idx(total_density_feature_file_name)

# Set the plot idx
plot_idx = c(sample(oridb_idx.l$g0), sample(oridb_idx.l$g1), sample(oridb_idx.l$g2))

# Get the lengths of each dataset
idx_length.v = c(length(oridb_idx.l$g0), length(oridb_idx.l$g1), length(oridb_idx.l$g2))

# Set the y-axis positions
y_axis = numeric(3)
for(i in 1:length(y_axis)){
	if(i == 1){
		y_axis[i] = idx_length.v[i] / 2
	}else{
		y_axis[i] = sum(idx_length.v[1:(i-1)]) + idx_length.v[i] / 2
	}
}

# Set the dividing lines
div_line = c(idx_length.v[1], sum(idx_length.v[1:2]))

###################################
# Plotting

# Set up the overall screen
ovr_scr.m = matrix(c(0.075, 1, 0.05, 1,
		     0.075, 1, 0, 0.05,
		     0, 0.075, 0.05, 1),
		   ncol = 4, byrow = T
		  )

# Split the screen
ovr_scr.s = split.screen(ovr_scr.m)

# Set the plot title
plot_title.v = c("G1", "G2")

# Set the heatmap_mar
heatmap_mar = c(2, 2, 2, 0.75)

# Open the heatmap screen
screen(ovr_scr.s[1])

# Set up the heatmap scr
heatmap_scr.m = matrix(c(0, 0.5, 0, 1,
			 0.5, 1, 0, 1
			),
		       ncol = 4, byrow = T
		      )
heatmap_scr.s = split.screen(heatmap_scr.m)

# Set the screen
screen(heatmap_scr.s[1])

# Set up the plot
par(mar = heatmap_mar, mgp = c(3, 0.75, 0))

# Make the plot
dens_dot_plot(mat.l[[1]][plot_idx,], z_max = 4, 
	      plot_title = "",
	      x_axt = "n", y_axt = "n")
title(main = plot_title.v[1])
axis(1, at = c(-500.5, -250, 0, 250, 500.5), labels = c(-500, "", 0, "", 500))

axis(2, at = y_axis, 
	labels = c("No\nFootprint", "G1-Only\nFootprint", "G1 & G2\nFootprint"),
	cex.axis = 0.9, pos = -440, tick = F
    )

# Set the dividing line
abline(h = div_line, col = "black", lty = 2, lwd = 2)

# Set the screen
screen(heatmap_scr.s[2])

# Set up the plot
par(mar = heatmap_mar, mgp = c(3, 0.75, 0))

# Make the plot
dens_dot_plot(mat.l[[2]][plot_idx,], z_max = 4, 
	      plot_title = "",
	      x_axt = "n", y_axt = "n")
title(main = plot_title.v[2]) 
axis(1, at = c(-500.5, -250, 0, 250, 500.5), labels = c(-500, "", 0, "", 500))

axis(2, at = y_axis, labels = idx_length.v)

# Set the dividing line
abline(h = div_line, col = "black", lty = 2, lwd = 2)


# Add in the x-label
screen(ovr_scr.s[2])
par(mar = c(0, 0, 0, 0))
plot(0, 0, type = "n", bty = "n",
     xlim = c(0, 1), xaxs = "i", xaxt = "n",
     ylim = c(0, 1), yaxs = "i", yaxt = "n",
     xlab = "", ylab = "", main = ""
    )

# Add in the label
text(x = 0.5, y = 0.5, labels = "Relative distance from ACS (bp)", cex = 1.125)

# Add in the y-label
screen(ovr_scr.s[3])
par(mar = c(2, 0, 2, 0))
plot(0, 0, type = "n", bty = "n",
     xlim = c(0, 1), xaxs = "i", xaxt = "n",
     ylim = c(0, 798), yaxs = "i", yaxt = "n",
     xlab = "", ylab = "", main = ""
    )

# Add in the label
text(x = 0.4, y = 798/2, labels = "Putative OriDB origins", srt = 90, cex = 1.125)

# Add in the y-axis
# par(xpd = T)
# axis(2, at = y_axis_pos, labels = c(length(oridb_idx.l$g0), length(oridb_idx.l$g2)), pos = 1.4, tick = F)
# par(xpd = F)
