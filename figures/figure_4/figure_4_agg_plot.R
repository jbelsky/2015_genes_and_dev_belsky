# Clear the workspace
graphics.off()
rm(list = ls())

# Load the functions
function_files.v = list.files(path = "/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions/", full.names = T)
for(i in 1:length(function_files.v)){
	source(function_files.v[i])
}

# Set the dir
subpanel_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_4/"

# Set up the plot
# png(file = paste("/data/data2/jab112/2014_mnase_manuscript/figures/genes_and_dev_figures/",
#		 "fig_oriNuc.png", sep = ""),
#    width = 8.5, height = 8.5, units = "in", res = 512
#   ) 
tiff(file = paste("/data/data2/jab112/2014_mnase_manuscript/figures/genes_and_dev_figures/",
		  "fig_oriNuc_lzw.tiff", sep = ""),
     width = 8.5, height = 8.5, units = "in", res = 900, compression = "lzw"
    )
      

# Setup the screen
overall_scr.m = matrix(c(0, 0.5, 0, 1,
			 0.5, 1, 0.5, 1,
			 0.5, 1, 0, 0.5,
			 0, 0.03, 0.95, 1,
			 0, 0.03, 0.45, 0.5,
			 0.5, 0.53, 0.95, 1,
			 0.5, 0.53, 0.45, 0.5
			),
		       ncol = 4, byrow = T
		      )

# Split the overall screen
close.screen(all.screens = T)
overall_scr = split.screen(overall_scr.m)

# Panels A and B: Nuc Difference
screen(overall_scr[1])
source(paste(subpanel_dir, "figure_4A_and_4B_nuc_diff_heatmaps.R", sep = ""))

# Panel C: Subnuc Difference
screen(overall_scr[2])
source(paste(subpanel_dir, "figure_4C_subnuc_diff.R", sep = ""))

# Panel D: Origin Efficiency and Replication Timing
screen(overall_scr[3])
source(paste(subpanel_dir, "figure_4D_compare_nuc_origin_efficiency_and_timing.R", sep = ""))

# Add in the Figure Labels
screen(overall_scr[4])
add_figure_label("A")

screen(overall_scr[5])
add_figure_label("B")

screen(overall_scr[6])
add_figure_label("C")

screen(overall_scr[7])
add_figure_label("D")

# Close all the screens
close.screen(all.screens = T)

# Close the device
dev.off()
