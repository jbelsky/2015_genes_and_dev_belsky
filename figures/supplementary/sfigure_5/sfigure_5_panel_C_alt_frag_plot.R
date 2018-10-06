# Load the peak_mat.l file
load("/data/data2/jab112/2014_mnase_manuscript/figures/figure_5/figure_5_datasets/orc_and_mcm_peak_mat.l")

# Get the Mcm peaks for the idx
mcm_peak.df = peak_mat.l[["mcm"]][oridb_subset_idx,]

# Split the plot into 4 quadrants
typ_mat.m = matrix(c(0, 0.5, 0.4, 1,
		     0, 0.5, 0, 0.4,
		     0.5, 1, 0.4, 1,
		     0.5, 1, 0, 0.4
		    ), ncol = 4, byrow = T
		  )

# Get the screen
typ_mat.s = split.screen(typ_mat.m)

# Make the typhoon plot
screen(typ_mat.s[1])
par(mar = c(3.1, 3.1, 2.5, 0.5), mgp = c(1.5, 0.5, 0), cex.axis = 0.75, cex.lab = 0.65)

# Get the typhoon plot
typhoon.m = as.matrix(read.csv(paste("/data/data2/jab112/2014_mnase_manuscript/figures/supplementary/sfigure_5/",
				     "sfigure_5_datasets/left_mcm_cluster_oridb_feature_file_168_sites_agg_typhoon_plot.csv", sep = ""
				    )
			      )
		     )
colnames(typhoon.m) = -400:400

dens_dot_plot(typhoon.m[,as.character(-250:250)], z_max = 250, x_axt = "n", y_axt = "n",
	      lowCol = "white", highCol = "#686B6B"
	     )
axis(2, at = c(0.5, 100, 200), labels = F)
axis(2, at = c(50, 150, 250.5), labels = c(50, 150, 250))
axis(1, at = seq(-200, 200, 200))
axis(1, at = c(-100, 100), labels = F)
title(xlab = "Relative distance from ACS (bp)")
title(ylab = "Fragment length (bp)")
abline(v = c(-91, 164), col = "blue", lty = 2)

par(xpd = T)
mcm_pos.v = -90 + c(-17, 17)
plot_mcm(mcm_pos.v[1], x_w = 17, y0 = 275, yh = 20, obj_col = "purple")
plot_mcm(mcm_pos.v[2], x_w = 17, y0 = 275, yh = 20, obj_col = "purple")

screen(typ_mat.s[2])
par(mar = c(1, 3.1, 0.5, 0.5), mgp = c(1.5, 0.5, 0), cex.axis = 0.75, cex.lab = 0.65)

# Load the fragment lengths
frag.df = read.csv(paste("/data/data2/jab112/2014_mnase_manuscript/figures/supplementary/sfigure_5/sfigure_5_datasets/",
			 "left_mcm_cluster_oridb_feature_file_168_sites_nucleosome_dm243_fragment_length_distribution.csv", sep = ""
			)
		  )

# Get the distributions
frag_idx = which(frag.df$frag_length >= 50 & frag.df$frag_length <= 250)
left_dist.v = rep(50:250, frag.df[frag_idx,2])
right_dist.v = rep(50:250, frag.df[frag_idx,3])

# Make the boxplot
boxplot(list(left_dist.v, right_dist.v), outline = F,
	xlim = c(0, 1), xaxs = "i", at = c(0.3, 0.85), xaxt = "n", boxwex = 0.25,
	names = NA,
	ylim = c(50, 250)
       )
title(ylab = "Fragment length distribution")
text(x = c(0.1, 0.6), y = 165, labels = c("Upstream", "Downstream"), cex = 0.8, srt = 90)



# Make the typhoon plot
screen(typ_mat.s[3])
par(mar = c(3.1, 3.1, 2.5, 0.5), mgp = c(1.5, 0.5, 0), cex.axis = 0.75, cex.lab = 0.65)

# Get the typhoon plot
typhoon.m = as.matrix(read.csv(paste("/data/data2/jab112/2014_mnase_manuscript/figures/supplementary/sfigure_5/",
				     "sfigure_5_datasets/right_mcm_cluster_oridb_feature_file_228_sites_agg_typhoon_plot.csv", sep = ""
				    )
			      )
		     )
colnames(typhoon.m) = -400:400

dens_dot_plot(typhoon.m[,as.character(-250:250)], z_max = 340, x_axt = "n", y_axt = "n",
	      lowCol = "white", highCol = "#686B6B"
	     )
axis(2, at = c(0.5, 100, 200), labels = F)
axis(2, at = c(50, 150, 250.5), labels = c(50, 150, 250))
axis(1, at = seq(-200, 200, 200))
axis(1, at = c(-100, 100), labels = F)
title(xlab = "Relative distance from ACS (bp)")
title(ylab = "Fragment length (bp)")
abline(v = c(-91, 164), col = "blue", lty = 2)

# Get the left mcm peak idx
if(0){
mcm_pos.v = mcm_peak.df$rel_pos[right_idx]
mcm_pos.v = mcm_pos.v[!is.na(mcm_pos.v)]

y = rnorm(length(mcm_pos.v), mean = 210, sd = 10)
y[which(y > 235)] = 235
y[which(y < 185)] = 185
points(x = mcm_pos.v, y = y, col = "#6A287E95", pch = 19, cex = 0.3)
}

# Add in the Mcm
par(xpd = T)
mcm_pos.v = 160 + c(-17, 17)
plot_mcm(mcm_pos.v[1], x_w = 17, y0 = 275, yh = 20, obj_col = "purple")
plot_mcm(mcm_pos.v[2], x_w = 17, y0 = 275, yh = 20, obj_col = "purple")

screen(typ_mat.s[4])
par(mar = c(1, 3.1, 0.5, 0.5), mgp = c(1.5, 0.5, 0), cex.axis = 0.75, cex.lab = 0.65)

# Load the fragment lengths
frag.df = read.csv(paste("/data/data2/jab112/2014_mnase_manuscript/figures/supplementary/sfigure_5/sfigure_5_datasets/", 
			 "right_mcm_cluster_oridb_feature_file_228_sites_nucleosome_dm243_fragment_length_distribution.csv", sep = ""
			)
		  )

# Get the distributions
frag_idx = which(frag.df$frag_length >= 50 & frag.df$frag_length <= 250)
left_dist.v = rep(50:250, frag.df[frag_idx,2])
right_dist.v = rep(50:250, frag.df[frag_idx,3])

# Make the boxplot
boxplot(list(left_dist.v, right_dist.v), outline = F,
	xlim = c(0, 1), xaxs = "i", at = c(0.3, 0.85), xaxt = "n", boxwex = 0.25,
	names = NA,
	ylim = c(50, 250)
       )
title(ylab = "Fragment length distribution")
text(x = c(0.1, 0.6), y = 165, labels = c("Upstream", "Downstream"), cex = 0.8, srt = 90)
