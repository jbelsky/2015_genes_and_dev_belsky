# Clear the workspace
graphics.off()
rm(list = ls())

# Load the functions
function_files.v = list.files(path = "/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions/", full.names = T)
for(i in 1:length(function_files.v)){
	source(function_files.v[i])
}

# Set the dir
dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/"


#############################################
# Data Processing

# Set the work dir
work_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/"

# Set the file name types
file_name_type.v = c("g2_dm242", "g2_dm356",
		     "g2_orc1_161_dm261", "g2_orc1_161_dm334",
		     "orc_chip_seq_dm265", "orc_chip_seq_dm287"
		    )

# Set the file name footer
footer = "signal_around_oridb_acs_feature_file_curated_798_sites_win_500bp.csv"

# Set the total density signal file
total_density_feature_file_name = 
	paste(work_dir, "oridb_acs_feature_file_curated_798_sites_left_win_50bp_right_win_150bp.csv", sep = "")

# Create the storage matrix list
mat.l = vector("list", 3)
names(mat.l) = c("g2", "g2_orc1_161", "orc_chip_seq")

# Create the dataset matrix
for(i in 1:3){

	# Get the output density matrix
	mat.l[[i]] = average_matrices(convert_csv_to_matrix(paste(work_dir, file_name_type.v[2*i - 1], "_", footer, sep = "")),
				      convert_csv_to_matrix(paste(work_dir, file_name_type.v[2*i], "_", footer, sep = ""))
				     )

}

# Get the plot idx
oridb_idx.l = get_oridb_idx(total_density_feature_file_name)

# Get the non_idx for figure 2
non_idx = sort(unlist(oridb_idx.l[c("g0", "g1")]))

# Set the plot idx
plot_idx = c(sample(non_idx), sample(oridb_idx.l$g2))

###################################
# Plotting

# Set up the plot
#png(file = paste("/data/data2/jab112/2014_mnase_manuscript/figures/genes_and_dev_figures/",
#		 "fig_oriHeatmap.png", sep = ""),
#    width = 8.5, height = 8.5, units = "in", res = 512
#   )

tiff(file = paste("/data/data2/jab112/2014_mnase_manuscript/figures/genes_and_dev_figures/",
		  "fig_oriHeatmap_lzw.tiff", sep = ""),
     width = 8.5, height = 8.5, units = "in", res = 900, compression = "lzw"
    )

# Close all the screens
close.screen(all.screens = T)

# Set up the overall screen
ovr_scr.m = matrix(c(0.075, 1, 0.05, 1,
		     0.075, 1, 0, 0.05,
		     0, 0.075, 0.05, 1),
		   ncol = 4, byrow = T
		  )

# Split the screen
ovr_scr.s = split.screen(ovr_scr.m)

# Get the y-axis positions
y_axis_pos = c(length(non_idx)/2, length(non_idx) + (length(oridb_idx.l$g2)/2))

# Set the plot title
plot_title.v = c("WT", expression(paste(italic("orc1-161"))), "", "ORC ChIP-seq")

# Set the heatmap_mar
heatmap_mar = c(2, 0.975, 3, 0.975)

# Open the heatmap screen
screen(ovr_scr.s[1])

# Set up the heatmap scr
heatmap_scr.m = matrix(c(0, 0.25, 0, 1,
			 0.25, 0.5, 0, 1,
			 0.5, 0.75, 0, 1,
			 0.75, 1, 0, 1,
			 0, 0.03, 0.95, 1,
			 0.25, 0.28, 0.95, 1,
			 0.5, 0.53, 0.95, 1,
			 0.75, 0.78, 0.95, 1
			),
		       ncol = 4, byrow = T
		      )
heatmap_scr.s = split.screen(heatmap_scr.m)

# Set the screen
screen(heatmap_scr.s[1])

# Set up the plot
par(mar = heatmap_mar)

# Make the plot
dens_dot_plot(mat.l[[1]][plot_idx,], z_max = 4, 
	      plot_title = "",
	      x_axt = "n", y_axt = "n")
title(main = plot_title.v[1], line = 1) 
axis(1, at = c(-500.5, -250, 0, 250, 500.5), labels = c(-500, "", 0, "", 500))
axis(2, at = y_axis_pos, labels = F)

# Set the dividing line
abline(h = length(non_idx), col = "black", lty = 2, lwd = 2)

# Set the screen
screen(heatmap_scr.s[2])

# Set up the plot
par(mar = heatmap_mar)

# Make the plot
dens_dot_plot(mat.l[[2]][plot_idx,], z_max = 4, 
	      plot_title = "",
	      x_axt = "n", y_axt = "n")
title(main = plot_title.v[2], line = 1) 
axis(1, at = c(-500.5, -250, 0, 250, 500.5), labels = c(-500, "", 0, "", 500))
axis(2, at = y_axis_pos, labels = F)

# Set the dividing line
abline(h = length(non_idx), col = "black", lty = 2, lwd = 2)

# Make the difference plot
screen(heatmap_scr.s[3])
par(mar = heatmap_mar)

dens_dot_plot((mat.l[[1]] - mat.l[[2]])[plot_idx,], z_min = -4, z_max = 4, plot_title = "",
	      lowCol = "darkgreen", medCol = "white", highCol = "red", x_axt = "n", y_axt = "n"
	     )
axis(1, at = c(-500.5, -250, 0, 250, 500.5), labels = c(-500, "", 0, "", 500))
axis(2, at = y_axis_pos, labels = F)
abline(h = length(non_idx), col = "black", lty = 2, lwd = 2)

title(main = expression("WT"), col.main = "red", adj = 0.05, line = 1) 
title(main = expression(italic("orc1-161")), col.main = "darkgreen", adj = 0.95, line = 1)

# Make the ORC ChIP-seq plot
screen(heatmap_scr.s[4])
par(mar = heatmap_mar)

dens_dot_plot(mat.l[[3]][plot_idx,], z_max = 3, plot_title = "", x_axt = "n", y_axt = "n", 
	      highCol = "#7F7F7F")
axis(1, at = c(-500.5, -250, 0, 250, 500.5), labels = c(-500, "", 0, "", 500))
axis(2, at = y_axis_pos, labels = F)
abline(h = length(non_idx), col = "black", lty = 2, lwd = 2)

title(main = plot_title.v[4], line = 1) 

# Add in the figure labels
screen(heatmap_scr.s[5])
add_figure_label("A")

screen(heatmap_scr.s[6])
add_figure_label("B")

screen(heatmap_scr.s[7])
add_figure_label("C")

screen(heatmap_scr.s[8])
add_figure_label("D")

# Add in the x-label
screen(ovr_scr.s[2])
par(mar = c(0, 0, 0, 0))
plot(0, 0, type = "n", bty = "n",
     xlim = c(0, 1), xaxs = "i", xaxt = "n",
     ylim = c(0, 1), yaxs = "i", yaxt = "n",
     xlab = "", ylab = "", main = ""
    )

# Add in the label
text(x = 0.5, y = 0.5, labels = "Relative distance from ACS (bp)", cex = 1.25)

# Add in the y-label
screen(ovr_scr.s[3])
par(mar = c(2, 0, 3, 0))
plot(0, 0, type = "n", bty = "n",
     xlim = c(0, 1), xaxs = "i", xaxt = "n",
     ylim = c(0, 798), yaxs = "i", yaxt = "n",
     xlab = "", ylab = "", main = ""
    )

# Add in the label
text(x = 0.4, y = 798/2, labels = "Putative OriDB origins", srt = 90, cex = 1.25)

# Add in the y-axis
par(xpd = T)
axis(2, at = y_axis_pos, labels = c(length(non_idx), length(oridb_idx.l$g2)), pos = 1.4, tick = F)
par(xpd = F)

# Close the screens
close.screen(all.screens = T)

# Close the device
dev.off()
