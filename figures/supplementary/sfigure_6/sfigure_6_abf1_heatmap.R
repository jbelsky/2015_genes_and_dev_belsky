# Clear the workspace
rm(list = ls())
graphics.off()

# Load the functions
function_files.v = list.files(path = "/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions", full.names = T)
for(i in 1:length(function_files.v)){
	source(function_files.v[i])
}

# Create the plot
png(file = paste("/data/data2/jab112/2014_mnase_manuscript/figures/genes_and_dev_figures/",
		 "sfig_abf1Heatmap.png", sep = ""),
    width = 8.5, height = 8.5, units = "in", res = 400
   )

# Set the filenames
abf1_chip.fn = paste("/data/data2/jab112/2014_mnase_manuscript/figures/supplementary/sfigure_6/sfigure_6_datasets/",
		     "abf1_dm290_by_strand_win_500_bw_20.csv", sep = ""
		    )

abf1_nuc.fn = paste("/data/data2/jab112/2014_mnase_manuscript/figures/supplementary/sfigure_6/sfigure_6_datasets/",
		    c("g1_dm243", "g1_dm354"), "_nuc_signal_around_abf1_feature_file_win_500bp.csv", sep = ""
		   )



#########################################
# Data Processing

# Load the ChIP datasets
abf1.l = convert_strand_csv_to_matrix(abf1_chip.fn)

# Get the nucleosome density
abf1.m = average_matrices(convert_csv_to_matrix(abf1_nuc.fn[1]),
			  convert_csv_to_matrix(abf1_nuc.fn[2])
			 )

# Get the average nucleosome signal
abf1_nuc.v = apply(abf1.m, 2, mean)

# Get the average ChIP density
abf1_chip_fwd.v = apply(abf1.l$pos, 2, mean)
abf1_chip_rev.v = apply(abf1.l$neg, 2, mean)

# Get the nucleosome peaks
nuc_peaks.df = get_mod_peaks(abf1_nuc.v, x_mid = 0, peak_width = 150, min_thresh = 0.5)


#########################################
# Plotting

# Set up the plot
scr.m = matrix(c(0, 1, 0.9, 1,
		 0, 1, 0.6, 0.9,
		 0, 1, 0, 0.6
		), ncol = 4, byrow = T
	      )

# Get the screen
scr.s = split.screen(scr.m)

# Make the model plot
screen(scr.s[1])
par(mar = c(0, 8.1, 0, 8.1), cex.lab = 1.25, cex.axis = 1.25)

set_chromatin_schematic(x_start = -400, x_end = 400, y_start = 0, y_end = 1)

# Set the legend
legend(x = c(-400, 400), y = c(0, 1), bty = "n",
       legend = c("Fwd strand\ndensity", "Rev strand\ndensity", "Nucleosome\ndensity"),
       col = c("red", "darkgreen", "blue"),
       lwd = 2,
       horiz = T, cex = 1
      )

# Add in the text
# text(x = 0, y = 1.5, labels = "Abf1", cex = 1.25)

screen(scr.s[2])

# Set the margin
par(mar = c(1, 8.1, 0.5, 8.1))

# Make the plot
plot(0, 0, type = "n", 
     xlim = c(-400, 400), xaxs = "i", xaxt = "n",
     ylim = c(0, 8), yaxt = "n",
     xlab = "", ylab = "Average ChIP density"
    )
axis(1, at = seq(-400, 400, 200), labels = F)
axis(2, at = seq(0, 8, 2))
# axis(2, at = c(1.5, 4.5), labels = F)

# Add the ChIP density
lines(-500:500, abf1_chip_fwd.v, col = "red")
lines(-500:500, abf1_chip_rev.v, col = "darkgreen")

# Get the maximum position
x_pos = -500:500
abf1_max_fwd_pos = x_pos[which(abf1_chip_fwd.v == max(abf1_chip_fwd.v))]
abf1_max_rev_pos = x_pos[which(abf1_chip_rev.v == max(abf1_chip_rev.v))]
abline(v = c(abf1_max_fwd_pos, abf1_max_rev_pos), col = c("red", "darkgreen"), lty = 2)

# Plot the nucleosome
plot_nucleosome(nuc_peaks.df, y_max = 2.5, yh = 0.75, y0 = 7)
rect(xleft = -30, ytop = 7.75, xright = 30, ybottom = 6.25, col = "darkgreen")


par(new = T)
set_chromatin_schematic(x_start = -400, x_end = 400, y_start = 0, y_end = 3)
lines(x_pos, abf1_nuc.v, col = "blue")




# Move to the next screen
screen(scr.s[3])

# Set the margin
par(mar = c(5.1, 8.1, 0, 8.1), cex.lab = 1.25, cex.axis = 1.25)

# Make the plot
merge_coverage_plot(abf1.l$pos[,as.character(-400:400)], 
		    abf1.l$neg[,as.character(-400:400)], 
		    z_max = 10,
		   )
title(xlab = "Relative distance from Abf1 motif (bp)",
      ylab = "Abf1 binding site"
     )
axis(2, at = seq(0, 200, 100))
axis(2, at = seq(50, 250, 100), labels = F)
axis(1, at = seq(-400, 400, 200))




# Close the screens
close.screen(all.screens = T)

# Close the device
dev.off()
