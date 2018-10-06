# Clear the workspace
graphics.off()
rm(list = ls())

# Set the dir
dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_1/"

# Load the functions
function_files.v = list.files(path = "/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions/", full.names = T)
for(i in 1:length(function_files.v)){
	source(function_files.v[i])
}

# Set up the plot
# Original dimension was width = 12.5 x height = 8
tiff(file = paste("/data/data2/jab112/2014_mnase_manuscript/figures/genes_and_dev_figures/",
		  "fig_ARS107_2_lzw.tiff", sep = ""),
     width = 12.5, height = 8, units = "in", res = 900, compression = "lzw"
    )
#png(file = paste("/data/data2/jab112/2014_mnase_manuscript/figures/genes_and_dev_figures/",
#		 "fig_ARS107_2.png", sep = ""),
#    width = 12.5, height = 8, units = "in", res = 300
#   )

# Set up the overall screen
overall_scr.m = matrix(c(0, 0.72, 0, 1,
			 0.72, 1, 0, 1,
			 0, 0.03, 0.95, 1,
			 0.72, 0.75, 0.95, 1
			),
		       ncol = 4, byrow = T 
		      )

# Close all the screens
close.screen(all.screens = T)

# Setup the screen
split.screen(overall_scr.m)

# Enter in each component

# Typhoon Plot
screen(1)
source(paste(dir, "panel_A/fig_1A_ars107_typhoon.R", sep = ""))

# V-Plots
screen(2)
source(paste(dir, "panel_B/fig_1B_acs_aggregate_typhoon_plot.R", sep = ""))

# Figure labels
screen(3)
add_figure_label("A", label_cex_ratio = 12.5/8.5)

screen(4)
add_figure_label("B", label_cex_ratio = 12.5/8.5)

# Close the screens
close.screen(all.screens = T)

# Close the device
dev.off()
