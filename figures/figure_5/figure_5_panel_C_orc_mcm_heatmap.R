#############################################
# Set the Filenames
feature.fn = paste("/data/data2/jab112/2014_mnase_manuscript/datasets/jab112_yeast_feature_files/replication_origins/",
		   "oridb_acs_feature_file_curated_798_sites_timing_whitehouse_raw_oem_acs_seq.csv", sep = ""
		  )

oridb_subnuc_density.fn = paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/",
				"oridb_acs_feature_file_curated_798_sites_left_win_50bp_right_win_150bp.csv", 
				sep = ""
			       )

# Load the peak_mat
load("/data/data2/jab112/2014_mnase_manuscript/figures/figure_5/figure_5_datasets/orc_and_mcm_peak_mat.l")

# Set the work dir
work.dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_5/figure_5_datasets/"

# Set the file name types
datasets.fn.v = c("orc_chip_seq_dm265", "orc_chip_seq_dm287",
		  "mcm_chip_seq_dm82", "mcm_chip_seq_dm282"
		 )

# Set the file name footer
footer = "signal_by_strand_around_oridb_acs_feature_file_curated_798_sites_win_500bp.csv"

# Get the mcm single density
mcm_datasets.fn.v = paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/",
			  "mcm_chip_seq_dm", c(82, 282), "_signal_around_oridb_acs_feature_file_curated_798_sites_win_500bp.csv",
			  sep = ""
			 )

# Get the g1 nucleosome density
g1_nuc_datasets.fn.v = paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_4/figure_4_datasets/",
			     "g1_dm", c(243, 354), "_nuc_150_175_density_signal_around_",
			     "oridb_acs_feature_file_curated_798_sites_win_500bp.csv",
			     sep = ""
			    )


# Set the x_win
x_win = 400

#############################################
# Data Processing

# Load the feature file name
feature.df = read.csv(feature.fn)

# Get the plot idx
oridb_idx.l = get_oridb_idx(oridb_subnuc_density.fn)

# Get the idx 
idx = sort(unlist(oridb_idx.l[c("g1", "g2")]))

# Get the mcm peak over the G1 subset
mcm_peak.df = peak_mat.l[["mcm"]][idx,]

# Get the peak names
# left_names = mcm_peak.df[which(mcm_peak.df$rel_pos <= 35 & mcm_peak.df$signal > 1),"name"]
# right_names = mcm_peak.df[which(mcm_peak.df$rel_pos > 35 & mcm_peak.df$signal > 1),"name"]

# Get the indices from the feature.df file
# left_idx = which(feature.df$name %in% left_names)
# right_idx = which(feature.df$name %in% right_names)

# Get the mcm density
mcm.m = average_matrices(convert_csv_to_matrix(mcm_datasets.fn.v[1]),
			 convert_csv_to_matrix(mcm_datasets.fn.v[2])
			)[idx,]

# Get the Mcm sum in a 100 bp window around the left (-90) and right (160) nucleosomes
left_mcm_sig.v = apply(mcm.m[,as.character(-190:10)], 1, sum)
right_mcm_sig.v = apply(mcm.m[,as.character(60:260)], 1, sum)

# Get the left and right mcm indices
left_idx = which(left_mcm_sig.v > right_mcm_sig.v)
right_idx = which(left_mcm_sig.v < right_mcm_sig.v)

# Create the output peak list
mat.l = vector("list", 2)
names(mat.l) = c("orc", "mcm")

# Create the dataset matrix
for(i in 1:length(mat.l)){

	# Get the output lists for each replicate
	mat1.l = convert_strand_csv_to_matrix(paste(work.dir, datasets.fn.v[2*i - 1], "_", footer, sep = ""))
	mat2.l = convert_strand_csv_to_matrix(paste(work.dir, datasets.fn.v[2*i], "_", footer, sep = ""))

	# Get the output density matrix for each strand
	pos.m = average_matrices(mat1.l$pos, mat2.l$pos)
	neg.m = average_matrices(mat1.l$neg, mat2.l$neg)


	# Store the idx subset as a list in the mat.l
	mat.l[[i]] = list(pos = pos.m[idx,as.character(-400:400)], neg = neg.m[idx,as.character(-400:400)])

}

# Get the nucleosome density
g1_nuc.m = average_matrices(convert_csv_to_matrix(g1_nuc_datasets.fn.v[1]),
			    convert_csv_to_matrix(g1_nuc_datasets.fn.v[2])
			   )[idx,as.character(-400:400)]

# Get the nucleosome density for each subset
g1_left.v = apply(g1_nuc.m[left_idx,], 2, mean)
g1_right.v = apply(g1_nuc.m[right_idx,], 2, mean)

# Get the nucleosome peaks
left_nuc_peaks.df = get_mod_peaks(g1_left.v, x_mid = 0, peak_width = 150, min_thresh = 0.5)

right_nuc_peaks.df = get_mod_peaks(g1_right.v, x_mid = 0, peak_width = 150, min_thresh = 0.5)
right_nuc_peaks.df[4,] = c(347, 0.8881285)

# Set the plot_idx
idx.l = list(left_idx, right_idx)




#############################################
# Plot

# Create the plot idx
plot_idx = c(right_idx, left_idx)

# Create the plot axes
y_axis_pos = c(length(right_idx)/2, length(right_idx) + length(left_idx)/2)


# Outer scr
outer_scr.m = matrix(c(0.05, 1, 0.025, 0.925,
		       0, 0.05, 0.025, 0.925,
		       0.05, 1, 0.925, 1,
		       0.05, 1, 0, 0.045
		      ),
		     ncol = 4, byrow = T
		    )
os.s = split.screen(outer_scr.m)
screen(os.s[1])


# Set up the scr.m
# Model scr
model_scr.m = matrix(c(0, 1, 1/3, 2/3,
		       0, 1, 2/3, 1,
		       0, 1, 0, 1/3
		      ),
		     ncol = 4, byrow = T
		    )

lr_scr.m = matrix(c(0, 0.5, 0, 1,
		    0.5, 1, 0, 1
		   ),
		  ncol = 4, byrow = T
		 )


make_merge_plot = function(chip_idx){

	# Set the par
	par(mar = c(0.5, 1, 0.5, 1), tcl = -0.3, mgp = c(3, 0.25, 0), cex.axis = 0.85)
	merge_coverage_plot(mat.l[[chip_idx]]$pos[plot_idx,], mat.l[[chip_idx]]$neg[plot_idx,], z_max = 5)

	# Add in the axis
	axis(1, at = seq(-400, 400, 200), labels = F)
	if(chip_idx == 1){
		axis(2, at = y_axis_pos, labels = c(length(right_idx), length(left_idx)))
	}else{
		axis(2, at = y_axis_pos, labels = c("Downstream", "Upstream"), tick = F, line = 0)
	}

	# Divide the clusters
	abline(h = length(right_idx) + 0.5, col = "white", lwd = 2)

}

make_chip_plot = function(y_high, chip_idx, oridb_sub, par_b, par_t){

	par(mar = c(par_b, 1, par_t, 1), tcl = -0.3, mgp = c(3, 0.25, 0), cex.axis = 0.85)

	# Make the aggregate plot
	set_chromatin_schematic(x_start = -x_win, x_end = x_win, y_start = 0, y_end = y_high)
	box()

	axis(2, at = seq(0, y_high, 2), labels = T)
	axis(2, at = seq(1, y_high + 1, 2), labels = F)

	axis(1, at = seq(-x_win, x_win, x_win/2), labels = F)

	# Get the forward and reverse strand averages
	fwd.v = apply(mat.l[[chip_idx]]$pos[oridb_sub,], 2, mean)
	rev.v = apply(mat.l[[chip_idx]]$neg[oridb_sub,], 2, mean)

	# Enter in the lines
	lines(-x_win:x_win, fwd.v, col = "red")
	lines(-x_win:x_win, rev.v, col = "darkgreen")

	# Find the maximum amount
	fwd_idx = which(fwd.v == max(fwd.v))
	rev_idx = which(rev.v == max(rev.v))

	# Add in the boundary lines
	x = -x_win:x_win
	fwd_peak = x[fwd_idx]
	rev_peak = x[rev_idx]

	abline(v = c(fwd_peak, rev_peak), lty = 2, col = c("red", "darkgreen"))
	
	return(c(fwd_peak, rev_peak))

}

make_nucleosome_plot = function(oridb_sub, nuc_sub.v, nuc_p.df){

	par(new = T, tcl = -0.3, mgp = c(3, 0.25, 0), cex.axis = 0.85)

	# Set the schematic
	set_chromatin_schematic(x_start = -x_win, x_end = x_win, y_start = 0, y_end = 2)
	lines(-x_win:x_win, nuc_sub.v, type = "l", col = "blue")

	# Enter in the nucleosome peaks
	plot_nucleosome(nuc_p.df, y_max = 1.5, yh = 0.14, y0 = 1.75)

}

plot_orc = function(peak.v){

	# Get the center position
	orc_pos = mean(peak.v)

	# Plot ORC
	rect(xleft = orc_pos - 30, xright = orc_pos + 30, ybottom = 1.75 - 0.125, ytop = 1.75 + 0.125,
	     col = "darkgreen"
	    )

}

add_mcm = function(nuc_pos, mcm_edge_pos){

	# Get the mcm positions
	if(nuc_pos < 0){
		mcm_pos.v = (mcm_edge_pos - 68) + c(17, 51)
	}else{
		mcm_pos.v = (mcm_edge_pos + 68) - c(17, 51)
	}

	# Make the Mcm
	for(m in 1:2){
		plot_mcm(mcm_pos.v[m], x_w = 15, y0 = 1.75, yh = 0.155, obj_col = "purple")
	}

}


# Open the left-right screen
lr_scr.s = split.screen(lr_scr.m)
screen(lr_scr.s[1])

# Split the model screen
ms_scr.s = split.screen(model_scr.m)

screen(ms_scr.s[1])
make_merge_plot(1)

screen(ms_scr.s[2])
orc_peak.v = make_chip_plot(4, 1, idx.l[[1]], par_b = 1, par_t = 2)
title(main = "ORC")
make_nucleosome_plot(idx.l[[1]], g1_left.v, left_nuc_peaks.df)
plot_orc(orc_peak.v)

screen(ms_scr.s[3])
orc_peak.v = make_chip_plot(4, 1, idx.l[[2]], par_b = 2, par_t = 1)
axis(1, at = c(-x_win, 0, x_win), mgp = c(3, 0.25, 0), tcl = -0.3, cex.axis = 0.85)
make_nucleosome_plot(idx.l[[1]], g1_right.v, right_nuc_peaks.df)
plot_orc(orc_peak.v)

screen(lr_scr.s[2])

# Split the model screen
ms_scr.s = split.screen(model_scr.m)

screen(ms_scr.s[1])
make_merge_plot(2)

screen(ms_scr.s[2])
mcm_peak.v = make_chip_plot(6, 2, idx.l[[1]], par_b = 1, par_t = 2)
title(main = "Mcm2-7")
make_nucleosome_plot(idx.l[[1]], g1_left.v, left_nuc_peaks.df)
add_mcm(left_nuc_peaks.df$pos[2], mcm_peak.v[2])

screen(ms_scr.s[3])
mcm_peak.v = make_chip_plot(6, 2, idx.l[[2]], par_b = 2, par_t = 1)
axis(1, at = c(-x_win, 0, x_win), mgp = c(3, 0.25, 0), tcl = -0.3, cex.axis = 0.85)
make_nucleosome_plot(idx.l[[2]], g1_right.v, right_nuc_peaks.df)
add_mcm(right_nuc_peaks.df$pos[3], mcm_peak.v[1])



# Add in the labels
screen(os.s[2])
par(mar = c(2, 0, 2, 0))
set_chromatin_schematic()
text(x = 0.5, y = 0.5, labels = "ChIP density across OriDB sites", srt = 90)
text(x = 0.5, y = c(0.14, 1 - 0.14), labels = "Average ChIP density", srt = 90)

screen(os.s[3])
par(mar = c(0, 0, 0, 0))
set_chromatin_schematic()
legend(x = c(0, 0.25), y = c(0, 1), legend = c("Fwd strand\ndensity"), lwd = 2, col = "red", adj = 0, bg = NA, bty = "n")
legend(x = c(0.31, 0.5), y = c(0, 1), legend = c("Rev strand\ndensity"), lwd = 2, col = "darkgreen", bg = NA, bty = "n")
legend(x = c(0.61, 1), y = c(0, 1), legend = c("Nucleosome\ndensity"), lwd = 2, col = "blue", bg = NA, bty = "n")

screen(os.s[4])
par(mar = rep(0, 4))
set_chromatin_schematic()
text(x = 0.5, y = 0.5, labels = "Relative distance from ACS (bp)")
