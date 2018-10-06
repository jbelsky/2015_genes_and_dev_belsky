# Set the working directory
dataset_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_1/panel_B/"
setwd(dataset_dir)

# Set the type
cell_type.v = c("g2", "g2", "g2_orc1_161")
data_type.v = c("nuc", "subnuc", "nuc")
dm_name.v = c("dm242", "dm242", "dm261")


# Set the fragment width thresholds
low_frag = c(150, 0, 150)
high_frag = c(175, 120, 175)

# Set the minimum threshold
min_thresh.v = c(0.25, 2, 0.25)

# Set the peak storage list
peak_storage.l = vector("list", 3)

for(i in 1:3){

	# Get the G2 nucleosome vector
	mat.df = read.csv(paste("eaton_orc_chipseq_acs_start_220_sites_",
				cell_type.v[i], "_", dm_name.v[i], "_", data_type.v[i], 
				"_", low_frag[i], "_", high_frag[i], "_density.csv", sep = ""))

	# Convert to matrix
	mat.m = as.matrix(mat.df[,-(1:4)])

	# Get the mean coverage
	nuc.v = apply(mat.m, 2, mean)

	# Find the peaks
	peak_storage.l[[i]] = get_mod_peaks(nuc.v, x_mid = 0, min_thresh = min_thresh.v[i], peak_width = 150)

}








# Setup the screen
scr.m = matrix(c(0.05, 1, 0.633, 0.925,
		 0.05, 1, 0.341, 0.633,
		 0.05, 1, 0.05, 0.341,
		 0.05, 1, 0, 0.05,
		 0, 0.075, 0.05, 1,
		 0.05, 1, 0.925, 1
		),
	       ncol = 4, byrow = T
	      )
s_num.v = split.screen(scr.m)

# Set the x_start and x_end
x_start = -400
x_end = 400

# Set the dm_id
dm_id = c("g2", "g2_orc1_161")

# Set the dm_name
dm_name = c("dm242", "dm261")

# Plot labels
plot_label.v = c("WT", expression(paste(italic("orc1-161"))))

# Set the storage matrix
mat.l = vector("list", 2)

# Make the plot
for(i in 1:2){

	# Set up the plot
	screen(s_num.v[i])
	par(mar = c(1, 3, 0.5, 2))

	# Load the aggregate typhoon plots
	mat.df = read.csv(paste(dataset_dir, "eaton_orc_chipseq_acs_start_220_sites_aggregate_typhoon_plot_",
					     dm_id[i], "_", dm_name[i], "_win_400_high_frag_250.csv", sep = ""))

	# Convert to a matrix
	mat.m = as.matrix(mat.df)
	colnames(mat.m) = -400:400

	# Correct the mat to be 30E6
	mat.m = mat.m * (30E6) / get_num_frag(dm_name[i], 0, 250)

	# Make the plot (20 for V-plot)
	dens_dot_plot(mat.m, z_max = 1000, x_axt = "n", y_axt = "n")
	axis(1, at = c(-400.5, -200, 0, 200, 400.5), labels = F)
	axis(2, at = c(0.5, 50, 100, 150, 200, 250.5), labels = c(0, "", 100, "", 200, ""), cex.axis = 1.125)

	# Add the text
	text(x = -390, y = 10, adj = c(0, 0), labels = plot_label.v[i], cex = 1.5)

	# Store in the list
	mat.l[[i]] = mat.m

	# Remove the mat.m and mat.df
	rm(mat.m, mat.df)

	if(i == 1){
		plot_nucleosome(peak_storage.l[[1]], y_max = 1.6, y0 = 225, yh = 15)
		plot_subnucleosome(peak_storage.l[[2]], y_max = 0.7, y_bot = 210, y_top = 240)

	}else{
		plot_nucleosome(peak_storage.l[[3]], y_max = 1.6, y0 = 225, yh = 15)
	}	

}

# Move to the next screen
screen(s_num.v[3])
par(mar = c(1, 3, 0.5, 2))

# Make the merge coverage plot (20 for v-plot)
merge_coverage_plot(mat.l[[1]], mat.l[[2]], z_max = 1000) 
axis(1, at = seq(-400, 400, 200), labels = F)
axis(2, at = c(1, seq(50, 250, 50)), labels = c(0, "", 100, "", 200, ""), cex.axis = 1.125)

# Add in the label
text(x = -390, y = 10, adj = c(0, 0), labels = "WT", cex = 1.5, col = "red")
text(x = 390, y = 10, adj = c(1, 0), labels = "orc1-161", cex = 1.5, col = "green", font = 3)

# Add in the x-label
screen(s_num.v[4])
par(mar = c(0, 3, 0, 2))
plot(0, 0, type = "n", bty = "n",
     xlim = c(x_start, x_end), xaxs = "i", xaxt = "n",
     ylim = c(0, 1), yaxs = "i", yaxt = "n",
     xlab = "", ylab = "", main = ""
    )

# Add in the label
par(xpd = T)
text(x = 0, y = 0.5, labels = "Relative distance from ACS (bp)", cex = 1.25)

# Add in the x-axis
axis(1, at = seq(-400, 400, 200), pos = 1.75, tick = F, cex.axis = 1.125)
par(xpd = F)

# Add in the y-label
screen(s_num.v[5])
par(mar = c(0, 0, 0, 0))
plot(0, 0, type = "n", bty = "n",
     xlim = c(0, 1), xaxs = "i", xaxt = "n",
     ylim = c(0, 1), yaxs = "i", yaxt = "n",
     xlab = "", ylab = "", main = ""
    )

# Add in the label
text(x = 0.5, y = 0.5, labels = "Fragment length (bp)", srt = 90, cex = 1.5)

# Add in the top
screen(s_num.v[6])
par(mar = c(0, 3, 0, 2))
set_chromatin_schematic(x_start = -400, x_end = 400, y_start = 0, y_end = 3)

x_len = seq(-175, 175, length.out = 11)
x_len[7:11] = x_len[7:11] - 5
x_len[8:11] = x_len[8:11] - 5
x_len[11] = x_len[11] - 10

x_len = x_len + 10

text(x = x_len, y = 1,
     labels = unlist(strsplit("TTTTATGTTTA", split = "")), adj = c(0.5, 0),
     cex = c(0.6245, 0.9526, 0.9486, 0.7747, 0.9130, 0.6917, 0.5692, 0.8893, 0.9842, 0.9881, 0.5415)
    )

segments(x0 = -400, y0 = 1.25, x1 = -200, y1 = 1.25, lwd = 2)
segments(x0 = 200, y0 = 1.25, x1 = 400, y1 = 1.25, lwd = 2)
segments(x0 = -400, y0 = 0.5, x1 = 400, y1 = 0.5, col = "#909090", lwd = 2)

text(x = -400, y = 2.5, labels = "Upstream", adj = 0)
text(x = 400, y = 2.5, labels = "Downstream", adj = 1)

par(xpd = T)
arrows(x0 = -410, y0 = 2.5, x1 = -485, y1 = 2.5, length = 0.1, lwd = 2)
arrows(x0 = 410, y0 = 2.5, x1 = 485, y1 = 2.5, length = 0.1, lwd = 2)
text(x = -425, y = 1.25, labels = "5'")
text(x = -425, y = 0.5, labels = "3'", col = "#909090")
text(x = 425, y = 1.25, labels = "3'")
text(x = 425, y = 0.5, labels = "5'", col = "#909090")
par(xpd = F)
	
segments(x0 = 0, y0 = 0, x1 = -175, y1 = 0.4, lty = 3)
segments(x0 = 0, y0 = 0, x1 = 175, y1 = 0.4, lty = 3)
