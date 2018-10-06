package mnasedens;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import functions.TF;

public class GetDensitySignalAroundFeature {
	
	public static HashMap<String, BufferedWriter> get_output_buffer_hash(
		String[] input_files, String output_file_dir, String footer) throws IOException{
		
		// Set up the storage hash
		HashMap<String, BufferedWriter> output_buffer_hash = new HashMap<String, BufferedWriter>();
		
		// Iterate through each dataset to open the output buffer
		for(String f : input_files){
			
			// Get the output file
			String output_file = output_file_dir + f + footer;
			
			// Set up the new output buffer
			output_buffer_hash.put(f, new BufferedWriter(new FileWriter(output_file)));
			
		}
	
		// Return the output buffer hash
		return(output_buffer_hash);
		
	}
	
	public static void write_output_header(HashMap<String, BufferedWriter> output_buffer_hash, int win) throws IOException{
		
		// Get the iterator over the keyset
		Iterator<String> itr = output_buffer_hash.keySet().iterator();
		
		// Write the output for each dataset
		while(itr.hasNext()){
			
			// Get the BufferedWriter
			BufferedWriter output = output_buffer_hash.get(itr.next());
			
			// Write the output subnucleosome header
			output.write("name,chr,pos,strand,");
			String sep = ",";
			for(int i = -win; i <= win; i++){
				if(i == win){sep = "\n";}
				output.write(i + sep);
			}
			
		}

	}
	
	public static String get_feature_output(TF feature, double[] vector, int win){
		
		// Get the double array subset
		double[] sig = new double[2*win + 1];
		
		// Set the sig_idx
		int sig_idx = 0;
		
		// Get the signal around the feature
		for(int p = feature.getPos() - win; p <= feature.getPos() + win; p++){
			
			// If p is negative or above the matrix, enter a -1
			if(p < 0 || p >= vector.length){
				sig[sig_idx++] = -1;
			}else{
				sig[sig_idx++] = vector[p];
			}
			
		}
		
		// If the feature is on the negative strand, reverse complement the sig
		if(feature.getStrand() == '-'){
			sig = reverse_complement_signal_vector(sig);
		}
		
		// Initialize the output string
		String output_str = feature.toString() + ",";
		
		// Set up the DecimalFormat
		DecimalFormat df = new DecimalFormat("#.####");
		
		// Initialize the sep
		String sep = ",";
		
		// Add to the output_str
		for(int i = 0; i < sig.length; i++){
			if(i == sig.length -1){
				sep = "\n";
			}
			output_str += (df.format(sig[i]) + sep);
		}

		// Return the output string
		return(output_str);
		
	}
	
	public static double[] reverse_complement_signal_vector(double[] vec){
	
		// Create the new vector
		double[] rc_vec = new double[vec.length];
		
		// Make the reverse complement
		for(int i = 0; i < vec.length; i++){
			rc_vec[(vec.length - 1) - i] = vec[i];
		}
		
		// Return the reverse complemented vector
		return(rc_vec);
		
	}

	public static void close_output_buffers(HashMap<String, BufferedWriter> output_buffer_hash) throws IOException{
		
		// Get the iterator over the keyset
		Iterator<String> itr = output_buffer_hash.keySet().iterator();
		
		// Write the output for each dataset
		while(itr.hasNext()){
			
			// Get the BufferedWriter
			output_buffer_hash.get(itr.next()).close();
			
		}
		
	}
	
	public static void main(String[] args) throws IOException {

		// Set the string headers
		String[] dataset_header = {"g2_dm242",
					   "g2_dm356"
					  };
				
		// Set the directory
		String dataset_dir = "/data/data2/jab112/2014_mnase_manuscript/datasets/nuc_data/quant_normalize/";
		
		// Set the footers
		String[] input_footer1 = {"nuc_150_175_density_signal_chr",
					  "nuc_150_175_density_signal_chr"
					 };
		
		String input_footer2 = "quant_normalize.csv";
		
		// Set the input feature file name
		String input_feature_file = 
			// "/data/data2/jab112/2014_mnase_manuscript/datasets/jab112_yeast_feature_files/" + 
			// "replication_origins/oridb_acs_feature_file_curated_798_sites_timing_whitehouse_raw_oem.csv";
			// "tf/macisaac_sacCer2/abf1.csv";
			"/home/jab112/feature_files/yeast/replication_origins/eaton_origins/" +
			"eaton_acs_filter_rDNA_match_oridb_acs_238_sites.csv";	
					
		// Set the window
		int win = 500;
		
		// Set the output file dir
		String output_dir = "/home/jab112/jb_temp/";
				
		// Set the output file footer
		String output_footer = 
			"_nuc_signal_around_eaton_238_sites_feature_file_win_" + 
			+ win + "bp.csv";				  
		
		// Set a reference bam file name
		String bam_file_name = "/data/illumina_pipeline/aligned_experiments/DM242/dm242.bam";
		
		////////////////////////////////////////////////////////////////////////////////////
		
		// Set the file names
		String[] input_file_name_arr = new String[dataset_header.length];
		for(int h = 0; h < input_file_name_arr.length; h++){
			input_file_name_arr[h] = dataset_dir +  
									 dataset_header[h] + "_" + input_footer1[h];
		}
		
		// Open the output buffers
		HashMap<String, BufferedWriter> output_hash = get_output_buffer_hash(dataset_header, output_dir, output_footer);
		
		// Write the headers for the output
		write_output_header(output_hash, win);
		
		// Get the TF HashMap
		HashMap<String, ArrayList<TF>> tf_map = TF.read_in_tf_map(input_feature_file);
		
		// Iterate through each TF
		for(int c = 1; c <= 16; c++){
			
			// Get the tf_list on the chr
			ArrayList<TF> tf_list = tf_map.get(Integer.toString(c));
			
			// Iterate through each feature
			if(tf_list != null){
						
				// Get the double array
				System.out.println("Reading in the data for chr " + c + "...");
				double[][] input_data = 
					GetTotalQuantNormSignalAroundFeature.get_chr_density(
						input_file_name_arr, input_footer2,	bam_file_name, Integer.toString(c)
						);
					
				// Iterate through each feature
				for(int f = 0; f < tf_list.size(); f++){
				
					// Get the feature
					TF feature = tf_list.get(f);
					
					// Iterate through each dataset
					for(int d = 0; d < dataset_header.length; d++){
					
						// Get the output buffer
						BufferedWriter output = output_hash.get(dataset_header[d]);
						
						// Write the output
						output.write(get_feature_output(feature, input_data[d], win));
						
					}
					
				}
					
			}
			
		}

		// Close the output buffer
		close_output_buffers(output_hash);
		
	}
		
}
