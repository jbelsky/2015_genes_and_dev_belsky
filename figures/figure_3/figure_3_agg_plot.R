# Clear the workspace
graphics.off()
rm(list = ls())

# Load the functions
function_files.v = list.files(path = "/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions/", full.names = T)
for(i in 1:length(function_files.v)){
	source(function_files.v[i])
}

# Set the dir
subpanel_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_3/"

# Set up the plot
# png(file = paste("/data/data2/jab112/2014_mnase_manuscript/figures/genes_and_dev_figures/",
#		 "fig_oriG1Only.png", sep = ""),
#    width = 8.5, height = 8.5, units = "in", res = 524
#   ) 
tiff(file = paste("/data/data2/jab112/2014_mnase_manuscript/figures/genes_and_dev_figures/",
		  "fig_oriG1Only_lzw.tiff", sep = ""),
     width = 8.5, height = 8.5, units = "in", res = 900, compression = "lzw"
    ) 

# Setup the screen
overall_scr.m = matrix(c(0, 1/2, 0, 1,
			 1/2, 1, 0, 1),
		       ncol = 4, byrow = T
		      )

# Split the overall screen
close.screen(all.screens = T)
overall_scr = split.screen(overall_scr.m)

# Enter the first screen and split 2/3 and 1/3
screen(overall_scr[1])
horiz_split_1.m = matrix(c(0, 1, 0.5, 1,
			   0, 1, 0, 0.5,
			   0, 0.06, 0.95, 1,
			   0, 0.06, 0.45, 0.5
			  ),
		       	 ncol = 4, byrow = T
		        )

horiz_split_1.s = split.screen(horiz_split_1.m)

# Panel A: Subnuc and MCM Heatmaps
screen(horiz_split_1.s[1])
source(paste(subpanel_dir, "figure_3A_oridb_acs_heatmap.R", sep = ""))

# Panel B: Origin Efficiency and Replication Timing
screen(horiz_split_1.s[2])
source(paste(subpanel_dir, "figure_3C_compare_origin_efficiency_and_timing.R", sep = ""))

# Add the Figure Label
screen(horiz_split_1.s[3])
add_figure_label("A")

screen(horiz_split_1.s[4])
add_figure_label("C")

# Split the 2nd screen into 2
screen(overall_scr[2])

horiz_split_2.m = matrix(c(0, 1, 0.5, 1,
			   0, 1, 0, 0.5,
			   0, 0.06, 0.96, 1,
			   0, 0.06, 0.46, 0.5
			  ),
		         ncol = 4, byrow = T
		        )

horiz_split_2.s = split.screen(horiz_split_2.m)

# Panel B: ORC Signal Comparison
screen(horiz_split_2.s[1])
source(paste(subpanel_dir, "figure_3B_compare_orc_density.R", sep = ""))

# Panel C: MCM Signal Comparison
screen(horiz_split_2.s[2])
source(paste(subpanel_dir, "figure_3D_compare_mcm_density.R", sep = ""))

# Add the Figure Label
screen(horiz_split_2.s[3])
add_figure_label("B")

screen(horiz_split_2.s[4])
add_figure_label("D")

# Close all the screens
close.screen(all.screens = T)

# Close the device
dev.off()
