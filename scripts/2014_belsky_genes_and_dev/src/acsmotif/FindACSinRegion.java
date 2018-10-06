package acsmotif;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;


public class FindACSinRegion {

	public static ArrayList<FoundMotif> find_motif_instances(
				double[][] pwm, String chr_seq,	ARS feature, double min_score, HashMap<String, Double> bg_model){
		
		// Set up the FoundMotif ArrayList
		ArrayList<FoundMotif> motif_list = new ArrayList<FoundMotif>();
		
		// Iterate over every position on both strands to evaluate the motif
		for(int p = feature.getStart(); p <= feature.getEnd(); p++){
			
			// Get the scores on each strand
			
			if(p + pwm[0].length < chr_seq.length()){
				
				double log_diff_score_fwd = PWMFunctions.get_log_score(chr_seq, pwm, p, '+', bg_model);
			
				if(log_diff_score_fwd >= min_score){
					
					// Create a new FoundMotif
					FoundMotif fm = new FoundMotif(feature.getName(), feature.getChr(), 
												   p, p + pwm[0].length - 1, '+', log_diff_score_fwd);
					
					// Add to the motif_list
					motif_list.add(fm);
					
				}
			
			}
			
			if(p - pwm[0].length > 0){

				// Get the score on the reverse strand
				double log_diff_score_rev = PWMFunctions.get_log_score(chr_seq, pwm, p, '-', bg_model);

				if(log_diff_score_rev >= min_score){

					// Create a new FoundMotif
					FoundMotif fm = new FoundMotif(feature.getName(), feature.getChr(), 
							p - pwm[0].length + 1, p, '-', log_diff_score_rev);

					// Add to the motif_list
					motif_list.add(fm);

				}

			}
			
		}
		
		// Sort by score
		Collections.sort(motif_list);
		
		// Return the motif list
		return(motif_list);
		
	}
	
	public static String get_motif_output(ArrayList<FoundMotif> motif_list){
		
		// Initialize the output string
		String output_str = "";
		
		for(int f = 0; f < motif_list.size(); f++){
			
			// Get the current motif
			FoundMotif cm = motif_list.get(f);
			
			output_str += cm.getOutput() + "\n";
			
		}
		
		// Return the output_str
		return(output_str);
		
	}
	
	public static String get_output_string(FoundMotif fm, ARS ars){
		
		// Create the DecimalFormat
		DecimalFormat df = new DecimalFormat("#.####");
		
		// Create the output string
		String output_str = ars.getName() + "," +
							ars.getChr() + "," +
							fm.getMotifStart() + "," +
							fm.getStr() + "," +
							df.format(fm.getScore()) + "," +
							ars.getStart() + "," +
							ars.getEnd() + "," +
							ars.getARSName() + "," +
							ars.getConfidence() + "\n";
		
		// Return the output string
		return(output_str);
							
	}
	
	public static void main(String[] args) throws IOException {

		// Set the input file name
		String input_feature_file = 
			"c:/Users/Jason/Desktop/2014_08_15_workspace/muller_diff.csv";
					
		// Set the output file name
		String output_file_name = 
			"c:/Users/Jason/Desktop/my_elon_out.csv"; 
					
		// Set the pwm file
		String pwm_file_name = "z:/datasets/jab112_yeast_feature_files/yeast/replication_origins/pwm/ACS.pwm";
		
		// Set the bg file name
		String bg_model_file_name = 
			"z:/datasets/jab112_yeast_feature_files/yeast/genome/yeast_1mer_to_6mer_background_nucleotide_freq.tsv";
				
		// Set the fasta file directory
		String fasta_file_header = 
			"z:/datasets/jab112_yeast_feature_files/yeast/genome/sacCer2_sgdR61_chr";
					
		// Set the threshold score
		double thresh_score = 4;

		///////////////////////////////////////////////////////////////////////////////////////
		
		// Get the TF HashMap
		HashMap<String, ArrayList<ARS>> ars_map = ARS.read_in_tf_map(input_feature_file);
			
		// Open the output buffer
		BufferedWriter output = new BufferedWriter(new FileWriter(output_file_name));
		
		// Write the header
		output.write("name,chr,pos,strand,motif_score,ars_start,ars_end,ars_name,oridb_status\n");
		
		// Read in the pwm
		double[][] pwm = PWMFunctions.read_in_pwm(pwm_file_name);
				
		// Get the background model
		HashMap<String, Double> bg_model = PWMFunctions.read_in_bg_model(bg_model_file_name);
		
		// Set the counter for number of ARSs without a high-scoring ACS
		int NA_ARS = 0;
		
		// Iterate through each TF
		for(int c = 1; c <= 16; c++){
			
			// Get the tf_list on the chr
			ArrayList<ARS> ars_list = ars_map.get(Integer.toString(c));
			
			// Iterate through each feature
			if(ars_list != null){
				
				// Get the chromosome sequence
				String chr_seq =  GenomeFunctions.read_in_chr_fasta_file(Integer.toString(c), fasta_file_header);	
		
				// Iterate through each feature
				for(int f = 0; f < ars_list.size(); f++){
					
					// Find the motif instances
					ArrayList<FoundMotif> motif_instance_list = 
						find_motif_instances(pwm, chr_seq, ars_list.get(f), thresh_score, bg_model);
						
					// Create the output
					if(motif_instance_list.size() > 0){
						output.write(get_output_string(motif_instance_list.get(0), ars_list.get(f)));
						
						for(int a = 0; a < motif_instance_list.size(); a++){
							FoundMotif fm = motif_instance_list.get(a);
							System.out.println(ars_list.get(f).getName() + "\t" + fm.getMotifStart() + "\t" + fm.getStr() + "\t" + fm.getScore());
						}
						
					}else{
						System.out.print(ars_list.get(f).getName() + " has no high scoring ACS\n");
						NA_ARS++;
					}
					
				}
				
			}
			
		}
		
		// Close the output buffer
		output.close();
		
		System.out.println("There are " + NA_ARS + " ARS regions without a high-scoring ACS");
		
	}

}
