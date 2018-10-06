# Set the function
compare_ori_eff_and_timing = function(oridb.l, xlab.v, ori_xlab_line = 2){

	###################################################
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
	par(mar = c(3.1, 4.1, 0.25, 2.1))

	# Set up the plot
	plot(0, 0, type = "n", 
	     ylab = "", xlab = "",
	     xlim = c(0.5, 2.5), xaxs = "i", xaxt = "n",
	     ylim = c(-0.02, 1), yaxt = "n",
	    )
	axis(2, at = seq(0.0, 1.0, 0.2), cex.axis = 1)
	par(mgp = c(3, ori_xlab_line, 0))
	axis(1, at = c(1, 2), labels = xlab.v)
	par(mgp = c(3, 2, 0))

	title(ylab = "Origin efficiency", line = 2.25)

	# Set the x_mid
	x_mid.v = 1:2

	# Set the storage ori.v
	ori.l = vector("list", 2)

	# Iterate through each dataset
	for(i in 1:2){

		# Get the efficiencies
		ori.v = oridb.l[[i]]$raw_oem

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
	if(0){
	segments(x0 = 1, y0 = 0.9, x1 = 1, y1 = 0.95, lwd = 1)
	segments(x0 = 1, y0 = 0.95, x1 = 2, y1 = 0.95, lwd = 1)
	segments(x0 = 2, y0 = 0.9, x1 = 2, y1 = 0.95, lwd = 1)
	text(x = 1.5, y = 0.96, adj = c(0.5, 0.5), labels = "*", cex = 1.25)
	}

	# Create the pie chart
	screen(scr_eff.s[2])
	par(mar = c(0, 4.1, 0, 2.1))

	# Set up the plot
	plot(0, 0, type = "n", bty = "n",
	     ylab = "", xlab = "",
	     xlim = c(0.5, 2.5), xaxs = "i", xaxt = "n",
	     ylim = c(0, 1), yaxs = "i", yaxt = "n"
	    )

	# Add in the legend
	par(xpd = T)
	text(x = 0.5, y = 0.75, adj = 0.5, labels = "Early", col = "#4DAF4A", font = 1)
	text(x = 0.5, y = 0.5, adj = 0.5, labels = "Late", col = "#E41A1C", font = 1)
	text(x = 0.5, y = 0.25, adj = 0.5, labels = "Dormant", col = "#999999", font = 1)
	par(xpd = F)

	# Make the chi_sq output
	chi_sq.m = matrix(nrow = 2, ncol = 3)

	# Get the proportion of early for each dataset (only if there is an efficiency)
	for(i in 1:2){

		# Get the timing vector for each dataset
		ori_sub.df = oridb.l[[i]]

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
	ori_eff = wilcox.test(ori.l[[1]], ori.l[[2]],
			      alternative = "two.sided", mu = 0, paired = FALSE)$p.value

	# Get the chi-square difference
	chi_sq = chisq.test(chi_sq.m)$p.value

	return(list(ori_eff = ori_eff, chi_sq = chi_sq))	

}
