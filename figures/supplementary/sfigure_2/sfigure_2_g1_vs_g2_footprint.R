# Clear the workspace
graphics.off()
rm(list = ls())

# Load the functions
function_files.v = list.files(path = "/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions/", full.names = T)
for(i in 1:length(function_files.v)){
	source(function_files.v[i])
}

#############################################
# Plotting

# Set up the plot
png(file = paste("/data/data2/jab112/2014_mnase_manuscript/figures/genes_and_dev_figures/",
		 "sfig_oriSubG1.png", sep = ""),
    width = 8.5, height = 8.5, units = "in", res = 300
   ) 
par(cex.axis = 1.25, cex.lab = 1.25, mar = c(4.1, 4.1, 1, 2.1), mgp = c(2.75, 1, 0))

# Split the screen in half
scr.m = matrix(c(0, 1, 0.5, 1,
		 0.5, 1, 0, 0.5,
		 0, 0.035, 0.95, 1,
		 0.5, 0.53, 0.95, 1,
		 0.5, 0.53, 0.45, 0.5
		),
	       ncol = 4, byrow = T
	      )

# Split the screen
close.screen(all.screens = T)
scr.s = split.screen(scr.m)

# Make the cdc6 panels
screen(scr.s[1])
source("/data/data2/jab112/2014_mnase_manuscript/figures/supplementary/sfigure_2/sfigure_2_panels_A_and_B.R")

# Make the qpcr panel
screen(scr.s[2])
source("/data/data2/jab112/2014_mnase_manuscript/figures/supplementary/sfigure_2/sfigure_2_panel_C_qpcr.R")

##################################
# Add in the figure labels

# Add in the figure labels
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
