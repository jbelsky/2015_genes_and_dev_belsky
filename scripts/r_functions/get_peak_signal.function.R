get_mod_peaks = function(cov.v, x_mid, min_thresh = 1, peak_width = 75){

	# Set the peak window
	peak_win = (peak_width - 1) / 2

	# Perform an FFT on the coverage
	cov_fft.v = filterFFT(cov.v, pcKeepComp = 0.02)

	# Find the peaks
	cov_peaks.v = peakDetection(cov_fft.v, threshold = min_thresh, score = F, width = 1)

	if(!is.null(cov_peaks.v)){

		# Convert to GenomicRanges
		peaks.gr = GRanges(seqnames = 1, 
				   ranges = IRanges(start = cov_peaks.v - 1, width = 1)
				  )

		# Enter in the peak position and the signal
		values(peaks.gr)$peak = mid(ranges(peaks.gr))
		values(peaks.gr)$sig = cov_fft.v[values(peaks.gr)$peak]

		neg_idx = numeric()

		# Ensure that the peak is the highest signal in the range
		for(i in 1:length(peaks.gr)){	

			# Get the maximum signal in the peak_width range
			start = values(peaks.gr)$peak[i] - peak_win
			end = values(peaks.gr)$peak[i] + peak_win

			if(start < 1){
				start = 1
			}else if(end > length(cov_fft.v)){
				end = length(cov_fft.v)
			}

			peak_range_sig = max(cov_fft.v[start:end])

			if(peak_range_sig > values(peaks.gr)$sig[i]){

				neg_idx = c(neg_idx, i)

			}

		}

		if(length(neg_idx) > 0){

			# Get the peak listing
			peaks.gr = peaks.gr[-neg_idx]

		}

		# Get the win
		win = (length(cov_fft.v) - 1)/2

		# Get the peaks
		peaks.df = data.frame(pos = values(peaks.gr)$peak - win,
				      sig = values(peaks.gr)$sig
				     )

		# Remove peaks within 40 bp of either end
		peaks.df = peaks.df[which(peaks.df$pos >= (-win + 40) & peaks.df$pos <= (win - 40)),]

	}else{

		peaks.df = data.frame(pos = numeric(), sig = numeric())

	}

	# Adjust the position
	peaks.df$pos = x_mid + peaks.df$pos

	return(peaks.df)

}

# Set the ChIP peak function
get_chip_peak_over_mat = function(chip.m, threshold_peak){

	# Set the output data frame
	output.df = data.frame(name = rownames(chip.m),
			       rel_pos = numeric(nrow(chip.m)),
			       signal = numeric(nrow(chip.m))
			      )

	# Iterate through each row
	for(i in 1:nrow(chip.m)){

		cat(output.df$name[i], ",", i, "\r", sep = "")

		# Get the peaks
		peaks.df = get_mod_peaks(chip.m[i,], 0, min_thresh = threshold_peak, peak_width = 150)	

		if(nrow(peaks.df) > 0){

			# Sort the peaks on signal intensity
			peaks.df = peaks.df[order(peaks.df$sig, decreasing = T),]
			
			# Fill out the output table
			output.df$rel_pos[i] = peaks.df$pos[1]
			output.df$signal[i] = peaks.df$sig[1]

		}else{

			# Set this row as NA
			output.df[i,c("rel_pos", "signal")] = c(NA, NA)

		}


	}

	cat("\n\tComplete!\n")

	# Return the output.df
	return(output.df)

}
# Set the function for finding nucleosome peaks
get_nuc_peak = function(nuc_output.df, rel_pos, nuc_win = 100){

	# Find the idx of the nuc_output within the range
	nuc_idx = which(abs(nuc_output.df$pos - rel_pos) <= nuc_win)

	if(any(nuc_idx)){
		return(nuc_output.df$pos[nuc_idx])
	}else{
		return(NA)
	}

}
