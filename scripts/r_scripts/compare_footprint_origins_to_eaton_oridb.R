# Clear the workspace
graphics.off()
rm(list = ls())

# Load the library
library(GenomicRanges)

# Load the functions
function_files.v = list.files(path = "/data/data2/jab112/2014_mnase_manuscript/scripts/r_functions/", full.names = T)
for(i in 1:length(function_files.v)){
	source(function_files.v[i])
}

#############################################
# Filenames
work_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/"

total_density.fn = "oridb_acs_feature_file_curated_798_sites_left_win_50bp_right_win_150bp.csv"

dataset_dir = "/data/data2/jab112/2014_mnase_manuscript/datasets/jab112_yeast_feature_files/replication_origins/"

oridb.fn = "oridb_acs_feature_file_curated_798_sites_timing_whitehouse_raw_oem_acs_seq.csv"
eaton_orc_chip.fn = "eaton_GSM424494_wt_G2_orc_chip_combined_peak_regions.csv"

# Get the eaton ACS
eaton_acs.fn = "eaton_acs_filter_rDNA_match_oridb_acs_238_sites.csv"

#############################################
# Data Processing

# Set the work dir
work_dir = "/data/data2/jab112/2014_mnase_manuscript/figures/figure_2/figure_2_datasets/"

# Set the total density signal file
total_density_feature_file_name = paste(work_dir, total_density.fn, sep = "")

# Get the plot idx
oridb_idx.l = get_oridb_idx(total_density_feature_file_name)

# Get the oridb feature file
oridb.df = read.csv(paste(dataset_dir, oridb.fn, sep = ""))

# Get the eaton_orc_chip.fn
eaton_orc.df = read.csv(paste(dataset_dir, eaton_orc_chip.fn, sep = ""))

# Get the eaton_acs.fn
eaton_acs.df = read.csv(paste(dataset_dir, eaton_acs.fn, sep = ""))

# Convert to GenomicRanges
oridb.gr = GRanges(seqnames = oridb.df$chr, ranges = IRanges(start = oridb.df$pos, width = 1))
eaton_orc.gr = GRanges(seqnames = eaton_orc.df$chr, ranges = IRanges(start = eaton_orc.df$start, end = eaton_orc.df$end))
eaton_acs.gr = GRanges(seqnames = eaton_acs.df$chr, ranges = IRanges(start = eaton_acs.df$pos, width = 1))


# Get the overlaps
oridb_eaton_overlap.v = countOverlaps(oridb.gr, eaton_orc.gr)

eaton_overlap.v = countOverlaps(eaton_acs.gr, eaton_orc.gr)
