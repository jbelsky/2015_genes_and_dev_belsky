# Clear the workspace
graphics.off()
rm(list = ls())

# Load the functions
function_files.v = list.files(path = "/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions/", full.names = T)
for(i in 1:length(function_files.v)){
	source(function_files.v[i])
}

# Set the dir
subpanel_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/supplementary/sfigure_4/"

# Set up the plot
png(file = paste("/data/data2/jab112/2014_mnase_manuscript/figures/genes_and_dev_figures/",
		 "sfig_oriNucAgg.png", sep = ""),
    width = 8.5, height = 8.5, units = "in", res = 300
   )

# Create the overall screen
scr.m = matrix(c(0, 0.5, 0, 1,
		 0.5, 1, 0.5, 1,
		 0, 0.05, 0.96, 0.99,
		 0, 0.05, 0.46, 0.49,
		 0.5, 0.55, 0.96, 0.99
		),
	       ncol = 4, byrow = T
	      )

# Close the screens
close.screen(all.screens = T)
scr.s = split.screen(scr.m)

# Make the nucleosome difference plots
screen(scr.s[1])
source(paste(subpanel_dir, "sfigure_4_panels_A_and_B_agg_nuc_diff_plots.R", sep = ""))

# Make the Mcm2-7 signal difference plots
screen(scr.s[2])
source(paste(subpanel_dir, "sfigure_4_panel_C_compare_mcm_density.R", sep = ""))

# Add the figure labels
screen(scr.s[3])
add_figure_label("A")

screen(scr.s[4])
add_figure_label("B")

screen(scr.s[5])
add_figure_label("C")

# Close the screens
close.screen(all.screens = T)

# Close the device
dev.off()
