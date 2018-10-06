# Load the nucleR library
library(GenomicRanges)
library(nucleR)

# Set up the gene
make_gene_schematic = function(feature_chr, feature_start, feature_end, 
			       y_low = 0, y_high = 1, cex_title = 1, bg_type = "white",
			       proteinCoding = T, geneName = T, omit_genes = NA, x_pos_title = 50
			      ){

	# Set up the plot
	plot(0, 0, type = "n", bty = "n", bg = bg_type,
	     xlim = c(feature_start, feature_end), xaxs = "i", xaxt = "n",
	     ylim = c(0, 1), yaxs = "i", yaxt = "n",
	     ann = F
	    )

	# Load the gene dataframe
	gene.df = read.csv(paste("/data/data2/jab112/mnase_data/macalpine_mnase/yeast_metal_analysis/feature_files/",
				 "sacCer2_ucsc_sgdGeneTable.csv", sep = ""))

	# Subset only on protein coding if selected
	if(proteinCoding){

		idx = which(gene.df$name != gene.df$sgd_name)

		gene.df = gene.df[idx,]

	}

	# Omit any gene if necessary
	if(any(!is.na(omit_genes))){
	
		gene.df = gene.df[-which(gene.df$sgd_name %in% omit_genes),]

	}

	# Convert to a GenomicRanges object
	gene.gr = GRanges(seqnames = gene.df$chr,
			  ranges = IRanges(start = gene.df$start, end = gene.df$end),
			  strand = gene.df$strand
			 )
	names(gene.gr) = gene.df$name

	# Create the feature gr
	feature.gr = GRanges(seqnames = feature_chr,	
			     ranges = IRanges(start = feature_start, end = feature_end)
			    )

	# Find the overlaps
	overlaps.df = as.data.frame(as.matrix(findOverlaps(feature.gr, gene.gr)))

	if(any(nrow(overlaps.df))){

		# Enter in the genes
		for(i in 1:nrow(overlaps.df)){
			plot_gene(gene.df[overlaps.df$subjectHits[i],], y_low, y_high, 
				  feature_start, feature_end, cex_title, geneName, x_pos_title)
		}

	}

}

# Plot gene
plot_gene = function(gene.v, ylow, yhigh, xstart, xend, cex_title, geneName, x_pos_title = 50){
	
	# Get y_mid
	ymid = (yhigh + ylow) / 2
	
	# Add in the text
	if(gene.v$strand == "+"){

		# Make the rectangle
		rect(gene.v$start, ymid + 0.1, gene.v$end, yhigh - 0.1, col = "gray")

		if(geneName){
			if(gene.v$start >= xstart){
				text(x = gene.v$start + x_pos_title, y = yhigh - 0.15, adj = c(0, 1),
				     labels = gene.v$sgd_name, font = 3, cex = cex_title)
			}else{
				text(x = gene.v$end - x_pos_title, y = yhigh - 0.15, adj = c(1, 1),
				     labels = gene.v$sgd_name, font = 3, cex = cex_title)
			}
		}
	}else{

		# Make the rectangle
		rect(gene.v$start, ylow + 0.1, gene.v$end, ymid - 0.1, col = "gray")

		if(geneName){
			if(gene.v$end <= xend){
				text(x = gene.v$end - x_pos_title, y = ylow + 0.15, adj = c(0, 1),
				     labels = gene.v$sgd_name, srt = 180, font = 3, cex = cex_title)
			}else{
				text(x = gene.v$start + x_pos_title, y = ylow + 0.15, adj = c(1, 1),
				     labels = gene.v$sgd_name, srt = 180, font = 3, cex = cex_title)
			}
		}
	}

}

# Set up the schematic section
set_chromatin_schematic = function(x_start = 0, x_end = 1, y_start = 0, y_end = 1){

	plot(0, 0, type = "n", bty = "n",
	     xlim = c(x_start, x_end), xaxs = "i", xaxt = "n",
	     ylim = c(y_start, y_end), yaxs = "i", yaxt = "n",
	     ann = F
	    )

}

# Make the nucleosome
plot_nucleosome = function(nuc.df, y_max, y0 = 0.5, yh = 0.2, nuc_col = "#FF0000"){
	
	# Set up the angle vector
	theta = seq(0, 2 * pi, length = 1000)

	# Set the y
	y = y0 + yh * sin(theta)

	for(i in 1:nrow(nuc.df)){

		# Get the position
		x0 = nuc.df$pos[i]

		# Find the coordinates for the nucleosome at each theta position
		x = x0 + 75 * cos(theta)

		# Find the signal color shading
		sig_shade = round(100 * nuc.df$sig[i] / y_max)
		
		if(sig_shade > 99){
			sig_shade = 99
		}

		sig_shade_str = formatC(sig_shade, flag = "0#", format = "d", width = 2)

		# Plot the nucleosome
		polygon(x, y, col = paste(nuc_col, sig_shade_str, sep = ""))

	}


}

# Make the subnucleosome
plot_subnucleosome = function(subnuc.df, y_max, y_bot, y_top){

	for(i in 1:nrow(subnuc.df)){

		# Get the position
		x_subnuc_pos = subnuc.df$pos[i]

		# Find the signal color shading
		sig_shade = round(100 * subnuc.df$sig[i] / y_max)
		
		if(sig_shade > 99){
			sig_shade = 99
		}

		# Enter the schematic for the subnucleosome
		rect(xleft = x_subnuc_pos - 25,  ybottom = y_bot,
		     xright = x_subnuc_pos + 25, ytop = y_top, 
		     col = paste("#006400", sig_shade, sep = "")
		    )

	}

}

# Make the nucleosome
plot_mcm = function(x0, x_w, y0 = 0.5, yh = 0.2, obj_col = "purple"){
	
	# Set up the angle vector
	theta = seq(0, 2 * pi, length = 1000)

	# Set the y
	y = y0 + yh * sin(theta)

	# Find the coordinates for the nucleosome at each theta position
	x = x0 + x_w * cos(theta)

	# Plot the nucleosome
	polygon(x, y, col = obj_col)

}

# Get the fragment length coverage
get_num_frag = function(dm_id, low_frag, high_frag){

	# Load the fragment length distribution
	frag.df = read.csv(paste("/data/data2/jab112/2014_mnase_manuscript/datasets/mnase_paired_end_fragment_length_distribution/", 
			   	 dm_id, "_fragment_length_distribution.csv", sep = ""))

	# Subset on the fragments in range
	idx = which(frag.df$frag_length >= low_frag & frag.df$frag_length <= high_frag)

	# Get the total number of fragments in the range
	total_frag = sum(frag.df$number_fragments[idx])

	return(total_frag)

}
