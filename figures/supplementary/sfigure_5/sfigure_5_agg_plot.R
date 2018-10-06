# Clear the workspace
# rm(list = ls())
graphics.off()

# Load the functions
function_files.v = list.files(path = "/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions", full.names = T)
for(i in 1:length(function_files.v)){
	source(function_files.v[i])
}

# Create the plot
png(file = paste("/data/data2/jab112/2014_mnase_manuscript/figures/genes_and_dev_figures/",
		 "sfig_oriMcmNuc.png", sep = ""),
    width = 8.5, height = 8.5, units = "in", res = 400
   )

# Create the overall screen
scr.m = matrix(c(0, 1, 0.5, 1,
		 0, 0.5, 0, 0.5,
		 0.5, 1, 0, 0.5,
		 0, 0.05, 0.96, 0.99,
		 0.5, 0.55, 0.96, 0.99,
		 0, 0.05, 0.46, 0.49,
		 0.5, 0.55, 0.46, 0.49
		),
	       ncol = 4, byrow = T
	      )

# Close the screens
close.screen(all.screens = T)
scr.s = split.screen(scr.m)

# Make the nucleosome difference plots
screen(scr.s[1])
source(paste("/data/data2/jab112/2014_mnase_manuscript/figures/supplementary/sfigure_5/", 
	     "sfigure_5_panels_A_and_B_nucleosome_density_diff.R", sep = ""
	    )
      )

# Make the fragment distribution plot
screen(scr.s[2])
source(paste("/data/data2/jab112/2014_mnase_manuscript/figures/supplementary/sfigure_5/", 
	     "sfigure_5_panel_C_alt_frag_plot.R", sep = ""))

# Make the footprint difference plots
screen(scr.s[3])
source(paste("/data/data2/jab112/2014_mnase_manuscript/figures/supplementary/sfigure_5/", 
	     "sfigure_5_panel_D_compare_ori_eff_and_timing.R", sep = ""))

# Add the figure labels
screen(scr.s[4])
add_figure_label("A")

screen(scr.s[5])
add_figure_label("B")

screen(scr.s[6])
add_figure_label("C")

screen(scr.s[7])
add_figure_label("D")

# Close the screens
close.screen(all.screens = T)

# Close the device
dev.off()
