# Load the file
qpcr.df = read.csv("/data/data2/jab112/2014_mnase_manuscript/figures/supplementary/sfigure_2/qpcr_files/wt_cdc6_orc_mcm_qpcr.csv", 
		   fill = T, row.names = 1
		  )

# Get the letters
letter.v = toupper(letters[1:16])

# Create the storage data frame
stor.df = data.frame(y1 = numeric(16), y2= numeric(16))
rownames(stor.df) = letter.v

# Get the digits
digits.v = formatC(1:10, format = "d", width = 2, flag = "0")

for(i in 1:16){

	# Get the char string
	char_string.v = paste(letter.v[i], digits.v, sep = "")

	# Get the standard curves
	y = qpcr.df[char_string.v[1:8],"cq"]
	x = log10(qpcr.df[char_string.v[1:8],"sq"])

	# Get the linear model coeff
	coeff.v = lm(y ~ x)$coefficients
	b = as.numeric(coeff.v[1])
	m = as.numeric(coeff.v[2])

	# Get the starting quantities for each
	r1 = 10^((qpcr.df[char_string.v[9],"cq"] - b)/m)
	r2 = 10^((qpcr.df[char_string.v[10],"cq"] - b)/m)

	# Store in the dataframe
	stor.df[i,"y1"] = r1
	stor.df[i,"y2"] = r2

}

# Get the row means
stor.df[,"avg"] = rowMeans(stor.df)

# Create a new dataframe
orc.df = data.frame(name = rep(c("ARS1", "ARS305", "Ctrl1", "Ctrl2"), 2),
		    type = rep(c("wt", "cdc6"), each = 4),
		    qpcr_avg = stor.df$avg[c(1:4,9:12)]
		   )

mcm.df = data.frame(name = rep(c("ARS1", "ARS305", "Ctrl1", "Ctrl2"), 2),
		    type = rep(c("wt", "cdc6"), each = 4),
		    qpcr_avg = stor.df$avg[c(5:8,13:16)]
		   )

# Enter into a list
chip.l = list(orc = orc.df, mcm = mcm.df)


##############################################
# Plotting

# Set up the screen
qpcr_scr.m = matrix(c(0.1, 1, 0.85, 1,
		      0.1, 1, 0.425, 0.85,
		      0.1, 1, 0, 0.425,
		      0, 0.1, 0, 0.9
		     ), ncol = 4, byrow = T
		   )

# Set up the screen
qpcr_scr.s = split.screen(qpcr_scr.m)

# Set the color
graph_col.v = c("#000000", "#444444", "#888888", "#CCCCCC")

# Add in the legend
screen(qpcr_scr.s[1])
par(mar = c(0, 2.1, 0, 2.1))
set_chromatin_schematic(x_start = 0, x_end = 10)

# Enter in the legend
legend(x = c(0, 10), y = c(0, 1), legend = c(expression(italic("ARS1")), expression(italic("ARS305")), "Ctrl1", "Ctrl2"),
       fill = graph_col.v, bty = "n", ncol = 2
      )

# Set the screen_idx
s_idx = c(2, 3)

for(i in 1:2){

	# Open the screen	
	screen(qpcr_scr.s[s_idx[i]])

	# Set the par
	par(mar = c(2.1, 2.1, 2.1, 2.1), mgp = c(3, 0.5, 0))

	# Get the data frame of interest
	mat.df = chip.l[[i]]

	# Get the relative position to the first plot
	rel_sig = mat.df$qpcr_avg[1]

	# Get the relative enrichment
	mat.df[,"rel_enrich"] = mat.df$qpcr_avg / rel_sig

	# Set up the plot
	plot(0, 0, type = "n",
	     xlim = c(0, 10), xaxs = "i", xaxt = "n",
	     ylim = c(0, 4), yaxt = "n",
	     ylab = "",
	     xlab = "",
	     main = c("ORC", "Mcm2-7")[i]
	    )
	axis(2, at = 0:4)
	axis(1, at = 2.5, labels = "WT", tick = F)
	axis(1, at = 7.5, labels = "cdc6-1", font.axis = 3, tick = F)

	# Set the positioning
	x_start = c(0, 5)

	# Set the plot_idx
	plot_idx = 1

	# Enter the bar graphs
	for(a in 1:2){

		# Iterate through each set
		for(b in 1:4){

			# Get the values
			x_val = x_start[a] + b
			y_val = mat.df$rel_enrich[plot_idx]

			# Make the rectangle
			rect(xleft = x_val - 0.5, xright = x_val + 0.5,
			     ybottom = 0, ytop = y_val,
			     col = graph_col.v[b]
			    )

			# Increment the plot_idx
			plot_idx = plot_idx + 1
		
		}

	}

}

# Add in the y-axis
screen(qpcr_scr.s[4])
par(mar = rep(0, 4))
set_chromatin_schematic()
text(x = 0.25, y = 0.5, label = "Fold enrichment (IP/input)", srt = 90)
text(x = 0.75, y = 0.5, label = expression(paste("relative to ", italic("ARS1"), " WT", sep = "")), srt = 90)
# text(x = 0.5, y = 0.5, label = expression(paste(plain("Fold enrichment (IP/input)\nrelative to "), italic("ARS1"), plain("WT"), sep = "")), srt = 90)
