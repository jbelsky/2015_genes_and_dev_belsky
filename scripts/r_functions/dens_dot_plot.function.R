# Load gplots
library(gplots)

######################################################################################################################################################
# dens_dot_plot: Makes a heatmap of data within a matrix
#
# PARAMETERS:
#	dot.m:		The length-midpoint matrix outputed by make_dens_matrix_GR
#	z_min:		The minimum intensity that should be plotted
#	z_max:		The maximum intensity that should be plotted
#	lowCol:		The color specifying the lowest intensity points
#	highCol:	The color specifying the highest intensity points
#	numColors:	The total number of colors for the plot
#	plot_title:	Title for the plot
#	x_label:	An x-label for the plot
#	y_label:	A y-label for the plot
#	use_row_names:	Whether to number rows as 1:nrow(dot.m) or use the designated rownames(dot.m)
#	plot_title_size: The cex number for the plot title
#	plot_title_line: The line number location for the plot title
#	plot_box:	A Boolean on whether to add a box around the entire plot
#	x_axt, y_axt:	Whether to plot ("s") or suppress ("n") the x and y axes
# RETURNS:
#	A heatmap plot utilizing the image function
dens_dot_plot = function(dot.m, z_min = 0, z_max = 100, 
			 lowCol = "white", medCol = "", highCol = "blue", numColors = 100,
		         plot_title = "", x_label = "", y_label = "", 
			 use_row_names = FALSE,
			 plot_title_line = NA, plot_box = TRUE, x_axt = "s", y_axt = "s"){

	# For points that are either above or below z_max or z_min respectively, set them to
	# the z_min and z_max (otherwise, plot shows arbitrary colors
        dot.m[which(dot.m >= z_max)] = z_max
	dot.m[which(dot.m <= z_min)] = z_min

	# Set the xValues conversion boolean
	xValues_convert = TRUE

	# Get the current column names
	if(!is.null(colnames(dot.m))){

		# Check if the values can be converted to numeric
		if(is.numeric(type.convert(colnames(dot.m)))){

			xValues = as.numeric(colnames(dot.m))

			# Flip the xValues conversion boolean
			xValues_convert = FALSE
		
		}

	}
		
	# If the xValues have to be converted, assume a span of 1
	if(xValues_convert){

		# Get the number of columns
		col_num = ncol(dot.m)
		
		# Get the window
		x_win = (col_num - 1) / 2

		# Set the column names
		xValues = -x_win:x_win
	}

	if(!use_row_names | is.null(rownames(dot.m))){
		yValues = 0.5:(nrow(dot.m) + 0.5)
	}else{
		yValues = as.numeric(rownames(dot.m))
	}

	# Make the colorpanel
	if(nchar(medCol) > 0){
		make_colorpanel = colorpanel(numColors, lowCol, medCol, highCol)
	}else{
		make_colorpanel = colorpanel(numColors, lowCol, highCol)
	}	

	# Make the heatmap utilizing the parameters specified above
        image(xValues, yValues, t(dot.m), col = make_colorpanel, zlim = c(z_min, z_max),
              xlab = x_label, ylab = y_label, xaxt = x_axt, yaxt = y_axt, bty = "n"
        ) 

	# Set the title
	title(main = plot_title, line = plot_title_line)

	# Add a box around the plot
	if(plot_box){ 
        	box(which = "plot", lty = "solid")
	}

}

merge_coverage_plot = function(g1.m, g2.m, z_max = 10, plot_bty = "o"){

	# Get the flank
	flank = (ncol(g1.m) - 1)/2
	
	# Get the frag_high
	frag_high = nrow(g1.m)

	g1.v = as.vector(t(g1.m))
	g2.v = as.vector(t(g2.m))

	g1.v[which(g1.v > z_max)] = z_max
	g2.v[which(g2.v > z_max)] = z_max
	g1.v[which(g1.v < 0)] = 0
	g2.v[which(g2.v < 0)] = 0

	# Create the coordinates for the graph
	x_left = rep((-flank:flank) - 0.5, frag_high)
	y_low = rep((1:frag_high) - 0.5, each = (2 * flank + 1))

	# Create the plot
	plot(0, 0, type = "n", bty = plot_bty,
	     xlim = c(-flank, flank), xaxs = "i", xaxt = "n",
	     ylim = c(0.5, frag_high + 0.5), yaxs = "i", yaxt = "n",
	     ann = F
	    )

	# Add in the colors
	rect(x_left, y_low, x_left + 1, y_low + 1, col = rgb(g1.v, g2.v, 0, z_max, maxColorValue = z_max), border = NA)

}

merge_coverage_plot_2 = function(g1.m, g2.m, z_max = 10, plot_bty = "o"){

	# Get the flank
	flank = (ncol(g1.m) - 1)/2
	
	# Get the frag_high
	frag_high = nrow(g1.m)

	g1.v = as.vector(t(g1.m))
	g2.v = as.vector(t(g2.m))

	g1.v[which(g1.v > z_max)] = z_max
	g2.v[which(g2.v > z_max)] = z_max
	g1.v[which(g1.v < 0)] = 0
	g2.v[which(g2.v < 0)] = 0

	# Create the coordinates for the graph
	x_left = rep((-flank:flank) - 0.5, frag_high)
	y_low = rep((1:frag_high) - 0.5, each = (2 * flank + 1))

	# Create the plot
	plot(0, 0, type = "n", bty = plot_bty,
	     xlim = c(-flank, flank), xaxs = "i", xaxt = "n",
	     ylim = c(0.5, frag_high + 0.5), yaxs = "i", yaxt = "n",
	     ann = F
	    )

	# Add in the colors

	r_channel = g1.v
	g_channel = g2.v
	b_channel = z_max

	r_channel[which(g1.v >= g2.v)] = z_max
	g_channel[which(g2.v >= g1.v)] = z_max
	b_channel[which(g1.v >= g2.v)] = b_channel - g1.v
	b_channel[which(g2.v > g1.v)] = b_channel - g2.v

	# rect(x_left, y_low, x_left + 1, y_low + 1, col = rgb(g1.v, g2.v, 0, z_max, maxColorValue = z_max), border = NA)
	rect(x_left, y_low, x_left + 1, y_low + 1, col = rgb(r_channel, g_channel, b_channel, z_max, maxColorValue = z_max), border = NA)

}
