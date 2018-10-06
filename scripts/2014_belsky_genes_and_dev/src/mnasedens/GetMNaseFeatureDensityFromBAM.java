package mnasedens;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import functions.*;
import net.sf.samtools.SAMFileReader;


public class GetMNaseFeatureDensityFromBAM {

	public static String get_output_header(int w){
		
		// Set the initial output String
		String output_string = "name,chr,pos,strand,";
		
		// Iterate through each win for the output string
		String sep = ",";
		for(int i = -w; i <= w; i++){
			if(i == w){sep = "\n";}
			output_string += (i + sep);
		}
		
		// Return the output string
		return(output_string);
		
	}

	public static void main(String[] args) throws IOException {

		// Set the bam file names
		String[] bam_file_name_arr =
			{"/data/data2/jab112/SRR001091_30bp.bam"
			};
		
		// Set the output file name header
		String[] output_file_header_arr = 
			{"/data/data2/jab112/SRR001091_30bp_around_origins"
			};
		
		// Set the input feature name
		String input_feature_file = 
			"/data/data2/jab112/2014_mnase_manuscript/datasets/jab112_yeast_feature_files/" +
			"replication_origins/oridb_acs_feature_file_curated_798_sites_timing_whitehouse_raw_oem.csv";
									
		// Set the window
		int win = 1000;
		
		// Set the low frag
		int low_frag = 150;
		
		// Set the high frag
		int high_frag = 175;
		
		// Set the bw
		double bw = 20;	
		
		// Set the output file name footer
		String output_file_footer = "_nuc_density_signal_" + low_frag + "_" + high_frag + 
									"_bw_" + bw + "_win_" + win + ".csv";
				
		//////////////////////////////////////////////////////////////////////////////////
		// Begin script here
		
		// Read in the feature file list
		ArrayList<TF> tf_list = TF.read_in_tf_list(input_feature_file);
	
		// Iterate through each file
		for(int id = 0; id < bam_file_name_arr.length; id++){
		
			// Get the read scale factor
			double read_scale_factor = 10.0E6 / 
										(double) BAMInput.get_number_aligned_reads(bam_file_name_arr[id]);
				
			// Open the output buffer
			BufferedWriter output = new BufferedWriter(
					new FileWriter(output_file_header_arr[id] + output_file_footer));
		
			// Write the header
			output.write(get_output_header(win));
		
			// Open the BAM File
			SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name_arr[id]), 
													   new File(bam_file_name_arr[id] + ".bai")
													  );
		
			// Iterate through each feature
			for(int i = 0; i < tf_list.size(); i++){
		
				if((i + 1) % 500 == 0){
					System.out.println("We are on TF " + (i+1) + " of " + tf_list.size());
				}
				
				// Get the subnucleosome coverage density over the region
				double[] feat_subnuc_cov = 
						DensityEst.get_mnase_coverage_density(bam_file, tf_list.get(i), win,
												   			  low_frag, high_frag,
						        							  bw, read_scale_factor
														   	 );
	
				// Write the subnuc output
				tf_list.get(i).write_reads_output(output, feat_subnuc_cov, 1);
		
			}
		
			// Close the output buffer
			output.close();

		}
		
	}
		
}
