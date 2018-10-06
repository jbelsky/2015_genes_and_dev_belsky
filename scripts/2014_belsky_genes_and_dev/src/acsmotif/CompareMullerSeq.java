package acsmotif;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;


public class CompareMullerSeq {

	public static HashMap<String, String> get_muller_seq(String input_file_name) throws IOException{
		
		// Open the input buffer
		BufferedReader input = new BufferedReader(new FileReader(input_file_name));
		
		// Read in the header
		String line = input.readLine();
		
		// Create the output hash
		HashMap<String, String> seq_hash = new HashMap<String, String>();
		
		// Add each of the sequences to the hash
		while((line = input.readLine()) != null){

			// Split the line on tab
			String[] line_arr = line.split("\t");
		
			// Enter into the hash
			seq_hash.put(line_arr[0], line_arr[1]);
			
		}
		
		// Close the input buffer
		input.close();
		
		// Return the hash
		return(seq_hash);
	
	}
	
	
	public static void main(String[] args) throws IOException {

		// Set the input file name
		String input_file_name = "c:/Users/Jason/Desktop/acs_test_updated_new_scores.csv";
		
		// Set the output file name
		String muller_file_name = 
			"c:/Users/Jason/Desktop/2014_08_15_workspace/muller_confirmed_acs.tsv";
		
		// Set the pwm file
		String pwm_file_name = "z:/datasets/jab112_yeast_feature_files/yeast/replication_origins/pwm/ACS.pwm";
		
		// Set the fasta file directory
		String fasta_file_header = 
			"z:/datasets/jab112_yeast_feature_files/yeast/genome/sacCer2_sgdR61_chr";
				
		////////////////////////////////////////////////////////////////////////////////////////////////
			
		// Read in the acs pwm
		double[][] tf_pwm = PWMFunctions.read_in_pwm(pwm_file_name);
	
		// Get the PWM width
		int pwm_width = tf_pwm[0].length;
		
		// Set the chr
		String chr = "0";
		
		// Initialize the chr_seq
		String chr_seq = "";
		
		// Get the Muller Hash
		HashMap<String, String> muller_hash = get_muller_seq(muller_file_name);
		
		// Open the input buffer
		BufferedReader input = new BufferedReader(new FileReader(input_file_name));
	
		// Read in the header
		String line = input.readLine();
		
		// Split the line on comma
		String[] line_arr = line.split(",");
		
		// Get the relevant indices
		int chr_idx = -1;
		int pos_idx = -1;
		int str_idx = -1;
		int ars_name_idx = -1;
		int ars_start_idx = -1;
		int ars_end_idx = -1;
		
		for(int h = 0; h < line_arr.length; h++){
			if(line_arr[h].equals("chr")){
				chr_idx = h;
			}else if(line_arr[h].equals("pos")){
				pos_idx = h;
			}else if(line_arr[h].equals("strand")){
				str_idx = h;
			}else if(line_arr[h].equals("ars_name")){
				ars_name_idx = h;
			}else if(line_arr[h].equals("ars_start")){
				ars_start_idx = h;
			}else if(line_arr[h].equals("ars_end")){
				ars_end_idx = h;
			}
				
		}
		
		// Iterate through each line
		while((line = input.readLine()) != null){
			
			// Split on comma
			line_arr = line.split(",");
			
			// Get the acs parameters
			String feat_chr = line_arr[chr_idx];
			String ars_name = line_arr[ars_name_idx];
			int pos = Integer.parseInt(line_arr[pos_idx]);
			char str = line_arr[str_idx].charAt(0);
			int ars_start = Integer.parseInt(line_arr[ars_start_idx]);
			int ars_end = Integer.parseInt(line_arr[ars_end_idx]);
			
			// Get the chromosome sequence if a new feat_chr
			if(!feat_chr.equals(chr)){
				chr = feat_chr;
				chr_seq =  GenomeFunctions.read_in_chr_fasta_file(chr, fasta_file_header);
			}


			// Get the sequence
			String origin_acs = "";
			if(str == '+'){
				origin_acs = chr_seq.substring(pos, pos + pwm_width);
			}else{
				origin_acs = GenomeFunctions.reverse_complement_sequence(
								chr_seq.substring(pos - pwm_width + 1, pos + 1)
								);
			}
		
			if(muller_hash.containsKey(ars_name)){
				
				if(!origin_acs.contains(muller_hash.get(ars_name))){
					System.out.println(line_arr[0] + "\t" + ars_name + "\t" + origin_acs + "\t" + muller_hash.get(ars_name));
				
					// Get the muller substring
					String muller_substring = muller_hash.get(ars_name);

					// Search across the ARS region for a sequence match
					for(int p = ars_start; p <= ars_end; p++){
						
						// Check a match
						if(chr_seq.substring(p, p + muller_substring.length()).equals(muller_substring)){
							System.out.println("\tMatch at position " + p + " on the '+' strand");
						}
						
						if(GenomeFunctions.reverse_complement_sequence(
								chr_seq.substring(p - muller_substring.length() + 1, p + 1)
								).equals(muller_substring)
						   ){
							System.out.println("\tMatch at position " + p + " on the '-' strand");
						}
						
					}
								
				}
					
			}			
		
		}
		
		// Close the input buffer
		input.close();
			
	}

}
