# Clear the workspace
graphics.off()
rm(list = ls())

# Load the functions
function_files.v = list.files(path = "/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions/", full.names = T)
for(i in 1:length(function_files.v)){
	source(function_files.v[i])
}

# Set the dir
subpanel_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_5/"

# Set up the plot
# png(file = paste("/data/data2/jab112/2014_mnase_manuscript/figures/genes_and_dev_figures/",
#		 "fig_oriMcm.png", sep = ""),
#    height = 8.5, width = 8.5, units = "in", res = 516
#   ) 
tiff(file = paste("/data/data2/jab112/2014_mnase_manuscript/figures/genes_and_dev_figures/",
		  "fig_oriMcm_lzw.tiff", sep = ""),
     height = 8.5, width = 8.5, units = "in", res = 890, compression = "lzw"
    ) 

# Setup the screen
overall_scr.m = matrix(c(0, 0.45, 0.5, 1,
			 0, 0.45, 0, 0.5,
			 0.45, 1, 0, 1,
			 0, 0.03, 0.95, 1,
			 0, 0.03, 0.45, 0.5,
			 0.45, 0.48, 0.95, 1
			),
		       ncol = 4, byrow = T
		      )

# Split the overall screen
close.screen(all.screens = T)
overall_scr = split.screen(overall_scr.m)

# Panel A: Typhoon Plot with ORC and MCM Peak
screen(overall_scr[1])
source(paste(subpanel_dir, "figure_5_panel_A_agg_typhoon_plot.R", sep = ""))

# Panel B: Read Schematic
screen(overall_scr[2])
source(paste(subpanel_dir, "figure_5_panel_B_create_read_schematic.R", sep = ""))

# Panel C: ORC and Mcm2-7 by Strand
screen(overall_scr[3])
source(paste(subpanel_dir, "figure_5_panel_C_orc_mcm_heatmap.R", sep = ""))

# Add in the Figure Labels
screen(overall_scr[4])
add_figure_label("A")

screen(overall_scr[5])
add_figure_label("B")

screen(overall_scr[6])
add_figure_label("C")

# Close all the screens
close.screen(all.screens = T)

# Close the device
dev.off()
