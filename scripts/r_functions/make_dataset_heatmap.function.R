# Functions List

convert_csv_to_matrix = function(input_file_name){

	# Read in the input file name
	mat.df = read.csv(input_file_name)

	# Convert to matrix
	mat.m = as.matrix(mat.df[,-(1:4)])

	# Find the window
	win = (ncol(mat.m) - 1)/2

	# Enter the column names
	colnames(mat.m) = -win:win

	# Return the matrix
	return(mat.m)

}

convert_strand_csv_to_matrix = function(input_file_name){

	# Read in the input file name
	mat.df = read.csv(input_file_name)

	# Convert to matrix
	mat.m = as.matrix(mat.df[,-(1:4)])

	# Find the total window
	total_win = ncol(mat.m)/2

	# Create the submatrices
	pos.m = mat.m[,1:total_win]
	neg.m = mat.m[,(total_win+1):ncol(mat.m)]

	# Find the window
	win = (ncol(pos.m) - 1)/2

	# Enter the column names
	colnames(pos.m) = -win:win
	colnames(neg.m) = -win:win

	# Create the output list
	mat.l = list(pos = pos.m, neg = neg.m)

	# Return the matrix list
	return(mat.l)

}

convert_typhoon_plot_to_matrix = function(input_file_name){

	# Read in the input file name
	mat.df = read.csv(input_file_name)

	# Convert to matrix
	mat.m = as.matrix(mat.df)

	# Find the window
	win = (ncol(mat.m) - 1)/2

	# Enter the column names
	colnames(mat.m) = -win:win

	# Return the matrix
	return(mat.m)

}

# Average the input matrices
average_matrices = function(a.m, b.m){

	# Set the output matrix
	output.m = (a.m + b.m) / 2

	# Ensure that no position is less than 0?
	output.m[which(output.m < 0)] = 0

	# Return output.m
	return(output.m)

}
