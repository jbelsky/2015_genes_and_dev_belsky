#############################################
# Set the file names

# Get the total density feature file name
total_density_feature_file_name = paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/",
					"oridb_acs_feature_file_curated_798_sites_left_win_50bp_right_win_150bp.csv", sep = "")

# Set the feature file names
nuc_type.df = read.csv(paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_4/figure_4_datasets/", 
			     "g1_g2_nuc_movement_around_oridb_acs_feature_file_curated_798_sites.csv", 
			     sep = "")
		      )
oridb.df = read.csv(paste("/data/data2/jab112/2014_mnase_manuscript/datasets/jab112_yeast_feature_files/replication_origins/",
			  "oridb_acs_feature_file_curated_798_sites_timing_whitehouse_raw_oem_acs_seq.csv", sep = "")
		   )

#############################################
# Data Processing

# Get the G1 and G2 subsets
oridb_idx.l = get_oridb_idx(total_density_feature_file_name)

# Get the oridb_subset
oridb_subset_idx = sort(unlist(oridb_idx.l[c("g1", "g2")]))

# Create the storage matrix list
mat.l = vector("list", 2)
names(mat.l) = c("g1", "g2")

# Subset on the oridb_subset_idx
nuc_type.df = nuc_type.df[oridb_subset_idx,]
oridb.df = oridb.df[oridb_subset_idx,]

# Get the idx
left = which(nuc_type.df$type == "left_movement")
right = which(nuc_type.df$type == "right_movement")
static = which(nuc_type.df$type == "static")

# Create the idx.l
idx.l = list(left = left, right = right, static = static)










#############################################
# Plotting

# Set up the plot
scr_eff.m = matrix(c(0, 1, 0, 0.75,
		     0, 1, 0.75, 1),
		   byrow = T, ncol = 4
		  )

# Split the screen
scr_eff.s = split.screen(scr_eff.m)

# Create the efficiency plot
screen(scr_eff.s[1])

# Set the plot parameters
par(mar = c(3.1, 4.1, 0.25, 2.1), mgp = c(2.5, 0.5, 0))

# Set up the plot
plot(0, 0, type = "n", 
     ylab = "Origin efficiency", xlab = "",
     xlim = c(0.5, 3.5), xaxs = "i", xaxt = "n",
     ylim = c(-0.02, 1), yaxt = "n",
    )
axis(2, at = seq(0.0, 1.0, 0.2), cex.axis = 1)
axis(1, at = c(1, 2, 3), labels = c("Upstream", "Downstream", "Static"))

# Set the storage ori.v
ori.l = vector("list", 3)

# Iterate through each dataset
for(i in 1:3){

	# Get the efficiencies
	ori.v = oridb.df$raw_oem[idx.l[[i]]]

	# Subset only on those ori with an efficiency
	ori.v = ori.v[which(ori.v > -1)]

	# Set anything less than 0.05 to 0
	ori.v[which(ori.v < 0.05)] = 0

	# Save the ori.v for statistical comparison later
	ori.l[[i]] = ori.v

	# Add in the Efficiency Points
	points(x = rnorm(length(ori.v), i, 0.075),
	       y = ori.v,
	       col = "#0000FF75", pch = 19, cex = 1
	      )

	# Add in the median
	segments(x0 = i - 0.25, y0 = median(ori.v),
		 x1 = i + 0.25, y1 = median(ori.v),
		 col = "red", lwd = 2
		)

}



# Add in the significance hatches
segments(x0 = 1, y0 = 0.85, x1 = 1, y1 = 0.985, lwd = 1)
segments(x0 = 1, y0 = 0.985, x1 = 3, y1 = 0.985, lwd = 1)
segments(x0 = 3, y0 = 0.985, x1 = 3, y1 = 0.85, lwd = 1)
text(x = 2, y = 0.96, adj = c(0.5, 0), labels = "*", cex = 1.25)

segments(x0 = 2, y0 = 0.85, x1 = 2, y1 = 0.925, lwd = 1)
segments(x0 = 2, y0 = 0.925, x1 = 3, y1 = 0.925, lwd = 1)
text(x = 2.5, y = 0.895, adj = c(0.5, 0), labels = "*", cex = 1.25)





# Create the pie chart
screen(scr_eff.s[2])
par(mar = c(0, 4.1, 0, 2.1))

# Set up the plot
plot(0, 0, type = "n", bty = "n",
     ylab = "", xlab = "",
     xlim = c(0.5, 3.5), xaxs = "i", xaxt = "n",
     ylim = c(0, 1), yaxs = "i", yaxt = "n"
    )

# Add in the legend
par(xpd = T)
text(x = 0.5, y = 0.75, adj = 0.5, labels = "Early", col = "#4DAF4A", font = 1)
text(x = 0.5, y = 0.5, adj = 0.5, labels = "Late", col = "#E41A1C", font = 1)
text(x = 0.5, y = 0.25, adj = 0.5, labels = "Dormant", col = "#999999", font = 1)
par(xpd = F)

# Make the chi_sq output
chi_sq.m = matrix(nrow = 3, ncol = 3)

# Get the proportion of early for each dataset (only if there is an efficiency)
for(i in 1:3){

	# Get the timing vector for each dataset
	ori_sub.df = oridb.df[idx.l[[i]],]

	# Get the proportion of each category
	dormant = length(which(ori_sub.df$raw_oem < 0.05))
	early = length(which(ori_sub.df$raw_oem >= 0.05 & ori_sub.df$timing == "early"))
	late = length(which(ori_sub.df$raw_oem >= 0.05 & ori_sub.df$timing == "late"))

	# Update the chi_sq matrix
	chi_sq.m[i,] = c(early, late, dormant)

	# Add the pie chart
	add_pie_chart_3(h1 = 1/4, h2 = 0, x0 = i, y0 = 0.5, perc.v = c(early, late, dormant)/nrow(ori_sub.df))

}

# Calculate the p-values for each
ori_eff.v = numeric(3)
chi_sq.v = numeric(3)
ori_eff_idx = 1
for(a in 1:2){
	for(b in (a+1):3){

		# Get the origin efficiency difference
		ori_eff.v[ori_eff_idx] = wilcox.test(ori.l[[a]], ori.l[[b]],
						     alternative = "two.sided", mu = 0, paired = FALSE)$p.value
		
		# Get the chi-square difference
		chi_sq.v[ori_eff_idx] = chisq.test(chi_sq.m[c(a,b),])$p.value

		# Update the index
		ori_eff_idx = ori_eff_idx + 1

	}
}
