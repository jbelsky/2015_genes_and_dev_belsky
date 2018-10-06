# Main dir
main_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_1/panel_A/"
setwd(main_dir)

# Load the feature file
feature.df = read.csv("ars107_feature_table.csv")
win = 1000

# Set the axes
feature_chr = feature.df$chr
x_mid = feature.df$mid
x_start = x_mid - win
x_end = x_mid + win

# Setup the screen
scr.m = matrix(c(0, 1, 0.9, 1,
		 0, 1, 0.8, 0.9,
		 0, 1, 0.425, 0.8,
		 0, 1, 0.05, 0.425,
		 0, 1, 0, 0.05,
		 0, 0.05, 0.05, 0.8),
	       ncol = 4, byrow = T
	      )

# Split the screens
s_num.v = split.screen(scr.m)

# Set the x_left parameter
x_left = 4.5

# Set up the plot
screen(s_num.v[1])
par(mar = c(0, x_left, 0, 2))

# Make the gene schematic
make_gene_schematic(feature_chr, x_start, x_end)

# Enter in the features (note these are based off of the subnuc peak position)
text(x = 124493, y = 0.5, labels = "ACS", col = "darkgreen", srt = 90)
text(x = 124720, y = 0.5, labels = "Abf1", col = "darkgreen", srt = 90)

# Set up the chromatin schematic
screen(s_num.v[2])
par(mar = c(0, x_left, 0, 2))

# Make the gene schematic
set_chromatin_schematic(x_start, x_end)

# Set the y-axis pos
y_axis_pos = c(0.75, 0.25)

# Set the plot labels
par(xpd = T)
text(x = x_start - 100, y = y_axis_pos, adj = c(0.5, 0.5), labels = c("WT", expression(italic("orc1-161"))), cex = 1.25)
par(xpd = F)

# Set the dm_type
dm_type = c("g2", "g2_orc1_161")

# Set the dm number
dm_num = c("dm242", "dm261")

# Iterate through each type
for(i in 1:2){

	if(i == 1){

		# Load the density plot
		subnuc.df = read.csv(paste("ars107_", dm_type[i], "_", dm_num[i], "_subnuc_0_120_density.csv", sep = ""))
		
		# Convert to a vector
		subnuc.v = as.numeric(subnuc.df[1,-(1:4)])

		# Get the subnuc peaks
		subnuc_peaks.df = get_mod_peaks(subnuc.v, x_mid, min_thresh = 3)			       
		
	}else{

		subnuc_peaks.df = subnuc_peaks.df[2,]

	}

	# Plot the subnucleosome
	plot_subnucleosome(subnuc_peaks.df, 2, y_axis_pos[i] - 0.225, y_axis_pos[i] + 0.225)

	# Load the density plot
	nuc.df = read.csv(paste("ars107_", dm_type[i], "_", dm_num[i], "_nuc_150_175_density.csv", sep = ""))
	
	# Convert to a vector
	nuc.v = as.numeric(nuc.df[1,-(1:4)])

	# Get the subnuc peaks
	nuc_peaks.df = get_mod_peaks(nuc.v, x_mid, min_thresh = 0.2, peak_width = 150)			       

	# Plot the subnucleosome
	plot_nucleosome(nuc_peaks.df, 5, y_axis_pos[i], 0.175)

}

# Set the screen_idx
screen_idx = s_num.v[c(3, 4)]

# Plot label
plot_label = c("WT", expression(paste(italic("orc1-161"))))

# Iterate through each type
for(i in 1:2){

	# Move to the second screen
	screen(screen_idx[i])
	par(mar = c(2, x_left, 1, 2), cex = 1.25)

	# Load the mat
	mat.df = read.csv(paste("ars107_typhoon_plot_", dm_type[i], "_", dm_num[i], "_win_1000_high_frag_250.csv", sep = ""))

	# Convert to matrix
	mat.m = as.matrix(mat.df)
	colnames(mat.m) = x_start:x_end

	# Get the total fragments
	total_frag = get_num_frag(dm_num[i], 0, 250)

	# Scale to 30E6 reads
	z_max = 10 * total_frag / (30E6)

	# Make the plot
	dens_dot_plot(mat.m, z_max = z_max, x_axt = "n", y_axt = "n", plot_title_line = 0)
	axis(1, at = c(x_start - 0.5, x_mid - 500, x_mid, x_mid + 500, x_end + 0.5), labels = F)
	axis(2, at = c(0.5, 100, 200), labels = c(0, 100, 200))
	axis(2, at = c(50, 150, 250.5), labels = F)

	# Enter the label
	text(x = x_start + 20, y = 10, adj = c(0, 0), labels = plot_label[i], cex = 1.5)

	# if(i == 1){
	#	par(xpd = T)
	#	arrows(x0 = 124850, y0 = 300, x1 = 124850, y1 = 200, col = "red", lwd = 2)
	#	arrows(x0 = 124717, y0 = 300, x1 = 124717, y1 = 120, col = "darkgreen", lwd = 2)
	#	par(xpd = F)
	# }

}

# Add in the x-label
screen(s_num.v[5])
par(mar = c(0, x_left, 0, 2))

# Set up the plot
set_chromatin_schematic(x_start, x_end)

# Add in the label
text(x = x_mid, y = 0.5, labels = "Chr I position (kb)", cex = 1.25)

# Add in the x-axis
axis(1, at = c(x_start, x_mid, x_end), 
	labels = formatC(c(x_start, x_mid, x_end) / 1000, format = "f", digits = 1), 
	pos = 2, tick = F, cex.axis = 1.25
    )

# Add in the y-label
screen(s_num.v[6])
par(mar = c(0, 0, 0, 0))
plot(0, 0, type = "n", bty = "n",
     xlim = c(0, 1), xaxs = "i", xaxt = "n",
     ylim = c(0, 1), yaxs = "i", yaxt = "n",
     xlab = "", ylab = "", main = ""
    )

# Add in the label
text(x = 0.5, y = 0.5, labels = "Fragment length (bp)", srt = 90, cex = 1.5)
