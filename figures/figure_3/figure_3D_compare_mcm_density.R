#############################################
# Set the file names

# Get the total density feature file name
total_density_feature_file_name = paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/",
					"oridb_acs_feature_file_curated_798_sites_left_win_50bp_right_win_150bp.csv", sep = "")

# Load the file containing the total mcm dataset
chip_total.df = read.csv(paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_3/figure_3_datasets/",
			       "oridb_acs_feature_file_curated_798_sites_mcm_chip_seq_left_win_400bp_right_win_400bp.csv", 
			      sep = ""
			     )
		       )

#############################################
# Data Processing

# Get the G1 and G2 subsets
oridb_idx.l = get_oridb_idx(total_density_feature_file_name)

# Create the average signal
chip_total.df[,"chip"] = (chip_total.df$mcm_chip_seq_dm82 + chip_total.df$mcm_chip_seq_dm282) / 2
chip_total.df[,"input"] = (chip_total.df$input_chip_seq_dm271_mcm + chip_total.df$input_chip_seq_dm272_mcm) / 2
chip_total.df[,"log2_diff"] = log2(chip_total.df$chip + 1) - log2(chip_total.df$input + 1)


############################################
# Plotting





# Set the mar
par(mar = c(3.1, 4.1, 1, 1))

# Set up the plot
plot(0, 0, type = "n", 
     ylab = expression(paste(log[2], " ratio Mcm2-7 ChIP-seq to input")), xlab = "",
     xlim = c(0.5, 3.5), xaxs = "i", xaxt = "n",
     ylim = c(-2.1, 5.1), yaxs = "i"
    )
par(mgp = c(3, 1.75, 0))
axis(1, at = 1:3, labels = c("G1 & G2\nFootprint", "G1-Only\nFootprint", "No\nFootprint"))
par(mgp = c(3, 1, 0))

# Set the x_mid
x_mid.v = 1:3

# Create the processed orc list
chip.l = vector("list", length = 3)
names(chip.l) = c("g2", "g1", "g0")

# Iterate through each datapoint
for(i in 1:3){

	# Get the ChIP Subset
	chip_subset.df = chip_total.df[oridb_idx.l[[i]],]

	# Get the log2 subset of points with sufficient sampling
	chip_subset.v = chip_subset.df$log2_diff[which(chip_subset.df$chip >= 40 & chip_subset.df$input >= 40)]

	# Enter the subset into the list
	chip.l[[i]] = chip_subset.v

	# Get the interquartile range
	q1 = quantile(chip_subset.v, 0.25)
	q3 = quantile(chip_subset.v, 0.75)
	iqr = q3 - q1

	# Eliminate outliers from plot
	chip_subset.v = chip_subset.v[which(chip_subset.v >= (q1 - 3 * iqr) &
					    chip_subset.v <= (q3 + 3 * iqr)
					   )
				     ]

	# Add in the G2 efficiency points
	points(x = rnorm(length(chip_subset.v), x_mid.v[i], 0.075),
	       y = chip_subset.v,
	       col = "#7F7F7F75", pch = 19, cex = 1
	      )
	segments(x0 = x_mid.v[i] - 0.25, y0 = median(chip_subset.v),
		 x1 = x_mid.v[i] + 0.25, y1 = median(chip_subset.v),
		 col = "red", lwd = 2)

}

y_low = 4.4
y_outer = y_low + 0.4
y_left = y_low + 0.1
y_right = y_low + 0.2

# Add in the significance hatches
segments(x0 = 1, y0 = y_low, x1 = 1, y1 = y_outer, lwd = 1)
segments(x0 = 2, y0 = y_low, x1 = 2, y1 = y_right, lwd = 1)
segments(x0 = 3, y0 = y_low, x1 = 3, y1 = y_outer, lwd = 1)

# Connect each of the segments
segments(x0 = 1, y0 = y_outer, x1 = 3, y1 = y_outer, lwd = 1)
# segments(x0 = 1, y0 = y_left, x1 = 2, y1 = y_left, lwd = 1)
segments(x0 = 2, y0 = y_right, x1 = 3, y1 = y_right, lwd = 1)

text(x = c(2, 2.5), y = c(y_outer, y_right) + 0.1, labels = "*", cex = 1.25)


# Create the p-value vector
chip_compare.v = numeric(3)

# Set the idx
idx = 1

for(a in 1:2){
	for(b in (a+1):3){

		# Get the p-values
		chip_compare.v[idx] = wilcox.test(chip.l[[a]], chip.l[[b]],
						  alternative = "two.sided", mu = 0, paired = FALSE)$p.value

		# Update the idx
		idx = idx + 1		

	}

}
