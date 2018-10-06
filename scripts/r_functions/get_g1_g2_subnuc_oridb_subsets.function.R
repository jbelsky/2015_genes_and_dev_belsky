# Functions to load in data

get_oridb_idx = function(input_file_name){

	# Load the oridb data files
	oridb.df = read.csv(input_file_name)

	# Get the average signal for the G2 oridb
	g1_oridb.v = (oridb.df$g1_dm243 + oridb.df$g1_dm354) / 2
	g2_oridb.v = (oridb.df$g2_dm242 + oridb.df$g2_dm356) / 2
	g2_orc1_oridb.v = (oridb.df$g2_orc1_161_dm261 + oridb.df$g2_orc1_161_dm334) / 2

	# Get the threshold for the g2_oridb and g2_orc1_oridb
	g2_oridb_thresh.v = g2_oridb.v * 1.5
	g2_orc1_oridb_thresh.v = g2_orc1_oridb.v * 1.5

	# Find the oridb sites that have at least 150 and are greater than the g2_orc1_qpois.v threshold
	# Threshold of 100, 125, 150
	# Based on bg distribution, corresponds to cdf of 0.7430, 0.7973, 0.8251
	g2_oridb_idx = which(g2_oridb.v >= 100 & g2_oridb.v >= g2_orc1_oridb_thresh.v)
	g2_oridb_non_idx = (1:nrow(oridb.df))[-g2_oridb_idx]

	# Find the g1_oridb_idx
	g1_oridb_idx = which(g1_oridb.v >= 100 & g1_oridb.v >= g2_oridb_thresh.v)

	# Subset only on those g1_oridb_idx that are not already in the g2_oridb_idx
	g1_oridb_idx = g2_oridb_non_idx[which(g2_oridb_non_idx %in% g1_oridb_idx)]
	no_bind_idx = g2_oridb_non_idx[which(!g2_oridb_non_idx %in% g1_oridb_idx)]

	# Return a list of the indices
	idx.l = list(g2 = g2_oridb_idx, g1 = g1_oridb_idx, g0 = no_bind_idx)
	return(idx.l)

}
