package acsmotif;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;


public class MotifScoreAtFeaturePosition {


	public static void main(String[] args) throws IOException {

		// Set the input file name
		String input_feature_file_name = 
			"c:/Users/Jason/Desktop/acs_test_updated.csv";
		
		// Set the output file name
		String output_file_name = 
			"c:/Users/Jason/Desktop/acs_test_updated_new_scores.csv";
		
		// Set the pwm file
		String pwm_file_name = "z:/datasets/jab112_yeast_feature_files/yeast/replication_origins/pwm/ACS.pwm";
		
		// Set the bg file name
		String bg_model_file_name = 
			"z:/datasets/jab112_yeast_feature_files/yeast/genome/yeast_1mer_to_6mer_background_nucleotide_freq.tsv";
				
		// Set the fasta file directory
		String fasta_file_header = 
			"z:/datasets/jab112_yeast_feature_files/yeast/genome/sacCer2_sgdR61_chr";
				
		////////////////////////////////////////////////////////////////////////////////////////////////
				
		// Open the input buffer
		BufferedReader input = new BufferedReader(new FileReader(input_feature_file_name));
		
		// Open the output buffer
		BufferedWriter output = new BufferedWriter(new FileWriter(output_file_name));
		
		// Read in the header
		String line = input.readLine();
		
		// Write the header
		output.write(line + "\n");
		
		// Read in the acs pwm
		double[][] tf_pwm = PWMFunctions.read_in_pwm(pwm_file_name);
		
		// Get the background model
		HashMap<String, Double> bg_model = PWMFunctions.read_in_bg_model(bg_model_file_name);
			
		// Set the chr
		String chr = "0";
		
		// Initialize the chr_seq
		String chr_seq = "";
		
		// Set the DecimalFormat
		DecimalFormat df = new DecimalFormat("#.####");
				
		// Iterate through each line
		while((line = input.readLine()) != null){
			
			// Split on comma
			String[] line_arr = line.split(",");
			
			// Get the chromosome
			String feat_chr = line_arr[1];
			
			// Get the chromosome sequence if a new feat_chr
			if(!feat_chr.equals(chr)){
				chr = feat_chr;
				chr_seq =  GenomeFunctions.read_in_chr_fasta_file(chr, fasta_file_header);
			}

			// Get the position and strand
			int pos = Integer.parseInt(line_arr[2]);
			char str = line_arr[3].charAt(0);
			
			// Find the log_diff_score
			double log_diff_score = PWMFunctions.get_log_score(chr_seq, tf_pwm, pos, str, bg_model);
				
			// Update the score
			line_arr[4] = df.format(log_diff_score);
			
			// Set the output_sep
			String output_sep = ",";
			
			// Write the output
			for(int s = 0; s < line_arr.length; s++){
				if(s == line_arr.length - 1){
					output_sep = "\n";
				}

				output.write(line_arr[s] + output_sep);

			}
								
		}
			
		// Close the input buffer
		input.close();
		
		// Close the output buffer
		output.close();
		
	}

}
