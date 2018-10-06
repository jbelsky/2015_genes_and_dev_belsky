# Clear the workspace
graphics.off()
rm(list = ls())

# Load the grid library
library(grid)

# Set the dir
dir = "/data/data2/jab112/2014_mnase_manuscript/figures/supplementary/sfigure_3/"
setwd(dir)

# Load the letter functions
source("pwm_letter_functions.R")

# Load the functions
function_files.v = list.files(path = "/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions", full.names = T)
for(i in 1:length(function_files.v)){
	source(function_files.v[i])
}

# Get the information content
get_position_ic_alt = function(site_prob.v, bg_dist.v){

	# Get the site prob <0
	idx = which(site_prob.v == 0)
	site_prob.v[idx] = 1E-4

	# Calculate the information content at the position
	position_ic = 0
	for(i in 1:4){
		position_ic = position_ic + site_prob.v[i] * log2(site_prob.v[i] / bg_dist.v[i])
	}

	return(position_ic)

}

# Get the information content
get_position_ic = function(site_prob.v, bg_dist){

	# Calculate the information content at the position
	position_ic = 0
	for(i in 1:4){
		position_ic = position_ic + -site_prob.v[i] * log2(site_prob.v[i])
	}

	# Calculate the information obtained at position
	position_ic_ob = bg_dist - position_ic

	return(position_ic_ob)

}


# Get the letter heights
get_nucleotide_heights = function(site_prob.v, bg_dist, x_pos, width){

	# Get the total information content at the position
	position_ic_ob = get_position_ic_alt(site_prob.v, bg_dist)
	
	# Calculate the height of each letter
	letter_height.v = site_prob.v * position_ic_ob

	# Calculate the letter height positions
	letter_height_pos.v = numeric(4)
	letter_height_idx = order(letter_height.v)

	letter_height_pos.v[0] = 0
	for(i in 2:4){
		letter_height_pos.v[i] = letter_height_pos.v[i-1] + letter_height.v[letter_height_idx[i-1]]	
	}
	
	# Get the list containing the data for each of the letters
	letter_info.l = vector("list", 4)

	# Get the x position
	letter_xpos_begin = x_pos - width/2

	# Get the A coordinates
	letter_info.l[[1]] = letterA(letter_xpos_begin, 
				     letter_height_pos.v[which(letter_height_idx == 1)], 
				     letter_height.v[1], 
				     width
				    )

	# Get the C coordinates
	letter_info.l[[2]] = letterC(letter_xpos_begin, 
				     letter_height_pos.v[which(letter_height_idx == 2)], 
				     letter_height.v[2], 
				     width
				    )

	# Get the G coordinates
	letter_info.l[[3]] = letterG(letter_xpos_begin, 
				     letter_height_pos.v[which(letter_height_idx == 3)], 
				     letter_height.v[3], 
				     width
				    )

	# Get the T coordinates
	letter_info.l[[4]] = letterT(letter_xpos_begin, 
				     letter_height_pos.v[which(letter_height_idx == 4)], 
				     letter_height.v[4], 
				     width
				    )

	return(letter_info.l)

}

get_pwm_matrix = function(sequence.v){

	# Set the output matrix
	mat.m = matrix(nrow = nchar(sequence.v[1]), ncol = 4)
	colnames(mat.m) = c("A", "C", "G", "T")
	rownames(mat.m) = 0:(nrow(mat.m) - 1)

	# Convert to a seq matrix
	seq.m = matrix(unlist(strsplit(sequence.v, split = "")), nrow = length(sequence.v), byrow = T)

	# Create the factor function
	get_nuc_percentage = function(x){

		# Create a factor
		seq.v = factor(x, levels = c("A", "C", "G", "T"))

		# Return the precentages
		return(summary(seq.v) / length(seq.v))

	}

	# Enter into the matrix
	mat.m = t(apply(seq.m, 2, get_nuc_percentage))
	rownames(mat.m) = 0:(nrow(mat.m) - 1)

	# Return the mat.m
	return(mat.m)

}



# Set up the plot
png(file = paste("/data/data2/jab112/2014_mnase_manuscript/figures/genes_and_dev_figures/",
		 "sfig_acsPWM.png", sep = ""),
    width = 8.5, height = 8.5, units = "in", res = 300
   ) 

# Set the oridb filename
oridb_file_name = paste("/data/data2/jab112/2014_mnase_manuscript/datasets/jab112_yeast_feature_files/replication_origins/",
			"oridb_acs_feature_file_curated_798_sites_timing_whitehouse_raw_oem_acs_seq.csv", sep = ""
		       )

# Set the total density file name
total_density_feature.fn = 
	paste("/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/",
	      "oridb_acs_feature_file_curated_798_sites_left_win_50bp_right_win_150bp.csv",
	      sep = ""
	     )

# Get the oridb idx
oridb_idx.l = get_oridb_idx(total_density_feature.fn)

# Get the oridb_feature
oridb.df = read.csv(oridb_file_name)

# Calculate the background information
a_bg_prob = 0.30901
c_bg_prob = 0.19167
g_bg_prob = 0.19130
t_bg_prob = 0.30801

# Get the background distribution
bg_dist.v = c(a_bg_prob, c_bg_prob, g_bg_prob, t_bg_prob)

# Set up the new plot
grid.newpage()

# Create the grid space
vp = viewport(x = unit(0.5, "npc"), y = unit(5/6, "npc"), width = unit(1, "npc"), height = unit(1/3, "npc"))
vp2 = viewport(x = unit(0.5, "npc"), y = unit(3/6, "npc"), width = unit(1, "npc"), height = unit(1/3, "npc"))
vp3 = viewport(x = unit(0.5, "npc"), y = unit(1/6, "npc"), width = unit(1, "npc"), height = unit(1/3, "npc"))
v.l = list(vp, vp2, vp3)
plot_vp = plotViewport(margins = c(4, 4, 2, 2))
data_vp = dataViewport(xscale = c(0.5, 33.5), yscale = c(0, 2))

for(i in 1:3){

	# Load the pwm
	pwm.m = get_pwm_matrix(oridb.df$sequence[oridb_idx.l[[i]]])

	pushViewport(v.l[[i]])

	grid.text(label = c("A", "B", "C")[i], x = unit(0.01, "npc"), y = unit(0.975, "npc"), hjust = 0, vjust = 1,
		  gp = gpar(fontsize = 24, fontface = "bold")	
		 )

	# Create the Data VP
	pushViewport(plot_vp)
	pushViewport(data_vp)


	# Add a rectangle around the grid
	grid.rect()

	# Set the axes
	grid.xaxis(at = c(5, 10, 15, 20, 25, 30), label = T)
	grid.yaxis(at = seq(0, 2, 0.5), label = T)

	# Set the axis label
	grid.text(label = "Position (bp)", y = unit(-3, "lines"))
	grid.text(label = "Information Content", x = unit(-3, "lines"), rot = 90)

	# Add the title
	grid.text(label = c("G1 & G2 Footprint", "G1-Only Footprint", "No Footprint")[i], 
		  y = unit(2, "native") + unit(0.5, "lines"),
		  gp = gpar(fontsize = 15)
		 )

	# Iterate through each position
	for(nuc_pos in 0:32){

		# Get the probability distribution at the site
		site_prob_dist.v = as.numeric(pwm.m[as.character(nuc_pos),])

		# Get the information content at the position
		position_info_content = get_position_ic_alt(site_prob_dist.v, bg_dist.v)

		if(position_info_content >= 0){

			# Get the list of letter information at the position
			cur_letter_info.l = get_nucleotide_heights(site_prob.v = site_prob_dist.v,
								   bg_dist = bg_dist.v,
								   x_pos = nuc_pos + 1, 
								   width = 0.75
								  )

			# Plot the letters
			for(i in 1:4){
				grid.polygon(x = unit(cur_letter_info.l[[i]]$x, "native"),
					     y = unit(cur_letter_info.l[[i]]$y, "native"),
					     id = cur_letter_info.l[[i]]$id,
					     gp = gpar(fill = cur_letter_info.l[[i]]$fill)
					    )
			}

		}

	}

	popViewport()
	popViewport()
	popViewport()

}

# Remove the viewports
popViewport(0)

# Close the device
dev.off()
