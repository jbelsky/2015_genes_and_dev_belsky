#############################################
# Set the Filenames
feature.fn = paste("/data/data2/jab112/2014_mnase_manuscript/datasets/jab112_yeast_feature_files/replication_origins/",
		   "oridb_acs_feature_file_curated_798_sites_timing_whitehouse_raw_oem_acs_seq.csv", sep = ""
		  )

oridb_subnuc_density.fn = paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/",
				"oridb_acs_feature_file_curated_798_sites_left_win_50bp_right_win_150bp.csv", 
				sep = ""
			       )


# Get the mcm single density
mcm_datasets.fn.v = paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/",
			  "mcm_chip_seq_dm", c(82, 282), "_signal_around_oridb_acs_feature_file_curated_798_sites_win_500bp.csv",
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

# Subset on the oridb_subset_idx
feature.df = feature.df[idx,]

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

######################################################################
par(cex.lab = 1, cex.axis = 1)

# Make the origin efficiency and timing plot
p_val.l = compare_ori_eff_and_timing(list(feature.df[left_idx,], feature.df[right_idx,]),
				     c("Upstream", "Downstream"),
				     ori_xlab_line = 0.5
				    )
title(xlab = "Mcm2-7 Loading Class", line = 2.75)
