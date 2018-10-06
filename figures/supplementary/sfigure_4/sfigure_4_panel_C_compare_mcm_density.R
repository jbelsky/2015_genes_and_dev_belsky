#############################################
# Set the file names

# Get the total density feature file name
total_density_feature_file_name = paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/",
					"oridb_acs_feature_file_curated_798_sites_left_win_50bp_right_win_150bp.csv", sep = "")

# Load the file containing the total mcm dataset
mcm_total.df = read.csv(paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_3/figure_3_datasets/",
			      "oridb_acs_feature_file_curated_798_sites_mcm_chip_seq_left_win_400bp_right_win_400bp.csv", 
			      sep = ""
			     )
		       )

# Set the feature file names
nuc_type.df = read.csv(paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_4/figure_4_datasets/", 
			     "g1_g2_nuc_movement_around_oridb_acs_feature_file_curated_798_sites.csv", 
			     sep = "")
		      )

#############################################
# Data Processing

# Get the G1 and G2 subsets
oridb_idx.l = get_oridb_idx(total_density_feature_file_name)

# Get the oridb_subset
oridb_subset_idx = sort(unlist(oridb_idx.l[c("g1", "g2")]))

# Subset on the oridb_subset_idx
nuc_type.df = nuc_type.df[oridb_subset_idx,]
mcm_total.df = mcm_total.df[oridb_subset_idx,]

# Get the idx
left = which(nuc_type.df$type == "left_movement")
right = which(nuc_type.df$type == "right_movement")
static = which(nuc_type.df$type == "static")

# Create the idx.l
idx.l = list(left = left, right = right, static = static)

# Create the average signal
mcm_total.df[,"mcm"] = (mcm_total.df$mcm_chip_seq_dm82 + mcm_total.df$mcm_chip_seq_dm282) / 2
mcm_total.df[,"input"] = (mcm_total.df$input_chip_seq_dm271_mcm + mcm_total.df$input_chip_seq_dm272_mcm) / 2
mcm_total.df[,"log2_diff"] = log2(mcm_total.df$mcm + 1) - log2(mcm_total.df$input + 1)


############################################
# Plotting

# Set the mar
par(mar = c(4.1, 4.1, 3, 1), mgp = c(2.5, 1, 0))

# Set up the plot
plot(0, 0, type = "n", 
     ylab = expression(paste(log[2], " ratio Mcm2-7 ChIP-seq to input")), xlab = "Nucleosome shift class",
     xlim = c(0.5, 3.5), xaxs = "i", xaxt = "n",
     ylim = c(-1, 4.5)
    )
axis(1, at = 1:3, labels = c("Upstream", "Downstream", "Static"))
par(mgp = c(3, 1, 0))

# Set the x_mid
x_mid.v = 1:3

# Create the mcm.l for statistical comparison later
mcm.l = vector("list", 3)

# Iterate through each datapoint
for(i in 1:3){

	# Get the mcm_subset.df
	mcm_subset.df = mcm_total.df[idx.l[[i]],]

	# Get the log2 subset of points with sufficient sampling
	mcm_subset.v = mcm_subset.df$log2_diff[which(mcm_subset.df$mcm >= 40 & mcm_subset.df$input >= 40)]

	# Store the subset for statistical comparison later
	mcm.l[[i]] = mcm_subset.v

	# Add in the G2 efficiency points
	points(x = rnorm(length(mcm_subset.v), x_mid.v[i], 0.075),
	       y = mcm_subset.v,
	       col = "#7F7F7F75", pch = 19, cex = 1
	      )

	segments(x0 = x_mid.v[i] - 0.25, y0 = median(mcm_subset.v),
		 x1 = x_mid.v[i] + 0.25, y1 = median(mcm_subset.v),
		 col = "red", lwd = 2)

}

y_low = 4
y_outer = 4.4
y_right = 4.2

# Add in the significance hatches
segments(x0 = 1, y0 = y_low, x1 = 1, y1 = y_outer, lwd = 1)
segments(x0 = 2, y0 = y_low, x1 = 2, y1 = y_right, lwd = 1)
segments(x0 = 3, y0 = y_low, x1 = 3, y1 = y_outer, lwd = 1)

# Connect each of the segments
segments(x0 = 1, y0 = y_outer, x1 = 3, y1 = y_outer, lwd = 1)
segments(x0 = 2, y0 = y_right, x1 = 3, y1 = y_right, lwd = 1)

text(x = c(2, 2.5), y = c(y_outer, y_right) + 0.1, labels = "*", cex = 1.25)

# Create the p-value vector
mcm_compare.v = numeric(3)

# Set the idx
idx = 1

for(a in 1:2){
	for(b in (a+1):3){
		
		# Get the Mcm2-7 difference
		mcm_compare.v[idx] = wilcox.test(mcm.l[[a]], mcm.l[[b]],
						 alternative = "two.sided", mu = 0, paired = FALSE)$p.value

		# Update the idx
		idx = idx + 1
		
	}

}
