# Load the nucleR library
library(GenomicRanges)
library(nucleR)

# Create pie chart
add_pie_chart = function(h1 = 1/3, h2 = 0, x0, y0, early_prop){

	# Get the height in inches
	height_plot_1 = grconvertY(h1, from = "npc", to = "inches")
	height_plot_2 = grconvertY(h2, from = "npc", to = "inches")

	# Get the x and y distances
	x_dist = grconvertX(height_plot_1, from = "inches", to = "user") - 
		 grconvertX(height_plot_2, from = "inches", to = "user")

	y_dist = grconvertY(height_plot_1, from = "inches", to = "user") - 
		 grconvertY(height_plot_2, from = "inches", to = "user")

	# Get the early_arc_length
	e_arc = 2 * pi * early_prop	

	# Enter into the pie chart
	e_theta = seq(pi/2, pi/2 - e_arc, length = 1000)

	# Set the parameters
	e_x = c(x0 + x_dist * cos(e_theta), x0)
	e_y = c(y0 + y_dist * sin(e_theta), y0)

	# Make the polygon
	polygon(e_x, e_y, col = "#4DAF4A")

	# Get the late_arc_length
	l_arc = 2 * pi - e_arc

	# Get the l_theta
	l_theta = seq(pi/2 - e_arc, -3 * pi / 2, length = 1000)
	
	# Set the parameters
	l_x = c(x0 + x_dist * cos(l_theta), x0)
	l_y = c(y0 + y_dist * sin(l_theta), y0)

	# Make the polygon
	polygon(l_x, l_y, col = "#E41A1C")
	
}

add_pie_chart_3 = function(h1 = 1/3, h2 = 0, x0, y0, perc.v){

	# Get the height in inches
	height_plot_1 = grconvertY(h1, from = "npc", to = "inches")
	height_plot_2 = grconvertY(h2, from = "npc", to = "inches")

	# Get the x and y distances
	x_dist = grconvertX(height_plot_1, from = "inches", to = "user") - 
		 grconvertX(height_plot_2, from = "inches", to = "user")

	y_dist = grconvertY(height_plot_1, from = "inches", to = "user") - 
		 grconvertY(height_plot_2, from = "inches", to = "user")

	# Set the colors
	colors.v = c("#4DAF4A", "#E41A1C", "#999999")

	# Set the current pi position
	current_pos = 0

	# Iterate through each percentage	
	for(i in 1:length(perc.v)){

		# Get the arc length
		arc = 2 * pi * perc.v[i]	
		
		# Get the end positions
		end_pos = current_pos + arc

		# Get the theta
		theta = seq((2*pi - current_pos) + (pi/2), (2*pi - end_pos) + (pi/2), length = 1000)

		# Set the parameters
		x = c(x0 + x_dist * cos(theta), x0)
		y = c(y0 + y_dist * sin(theta), y0)

		# Make the polygon
		polygon(x, y, col = colors.v[i])

		# Update the position
		current_pos = end_pos

	}

}

# Add Figure Label
add_figure_label = function(label, label_cex_ratio = 1){

	# Set up the plot
	par(mar = rep(0, 4))
	plot(0, 0, type = "n", bty = "n", axes = F)
	text(x = 0, y = 0, labels = label, cex = 2 * label_cex_ratio, font = 2)

}


# Plot the distribution of y-values
plot_sig_dist = function(y.v, x_pos, x_sd = 0.075, point_col = "#0000FF75"){

	# Add in the points
	points(x = rnorm(length(y.v), x_pos, x_sd), y = y.v,
	       col = point_col, pch = 19, cex = 1
	      )

	# Add in the median line
	segments(x0 = x_pos - 0.25, y0 = median(y.v),
		 x1 = x_pos + 0.25, y1 = median(y.v),
		 col = "red", lwd = 2
		)

}
