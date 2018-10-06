classify_nuc_movement = function(nuc_file_name, rel_pos_thresh = 0, left_col = "left0.975", right_col = "right0.025"){

	# Open the nucleosome data frame
	nuc.df = read.csv(nuc_file_name)

	# Get the indices corresponding to each type
	up_idx = which(nuc.df[,left_col] < -rel_pos_thresh)
	down_idx = which(nuc.df[,right_col] > rel_pos_thresh)

	# Check if there are any idx in both groups
	both_idx = up_idx[which(up_idx %in% down_idx)]

	if(any(both_idx)){

		# Get the absolute value median for each group
		both_left_med.v = abs(nuc.df$left0.5[both_idx])
		both_right_med.v = abs(nuc.df$right0.5[both_idx])

		# Compare the medians
		remove_right = both_idx[which(both_left_med.v >= both_right_med.v)]
		remove_left = both_idx[which(both_left_med.v < both_right_med.v)]

		if(any(remove_left)){

			up_idx = up_idx[!up_idx %in% remove_left]

		}

		if(any(remove_right)){
			
			down_idx = down_idx[!down_idx %in% remove_right]
		
		}

	}

	static_idx = (1:nrow(nuc.df))[-c(up_idx, down_idx)]

	# Classify each type
	nuc.df[up_idx,"type"] = "left_movement"
	nuc.df[down_idx,"type"] = "right_movement"
	nuc.df[static_idx,"type"] = "static"

	# Return the nuc.df
	return(nuc.df)

}
