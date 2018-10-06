package acsmotif;

import java.io.BufferedWriter;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class GetFeatureMotifSeq {

	public static void main(String[] args) throws IOException {

		// Set the input file name
		String input_file_name = args[0];
		
		// Set the output file name
		String output_file_name = args[1];
		
		// Set the fasta genome file directory
		String fasta_dir =  args[2];
		
		// Set the window
		int win = Integer.parseInt(args[3]);

		//////////////////////////////////////////////
		
		// Read in the input filename
		BufferedReader input = new BufferedReader(new FileReader(input_file_name));
		
		// Write the output
		BufferedWriter output = new BufferedWriter(new FileWriter(output_file_name));
				
		// Read in the header
		String line = input.readLine();
		
		// Output the header
		output.write(line + ",sequence\n");
		
		// Set the chr
		String chr = "0";
		
		// Initialize the sequence
		String seq = "";
		String chr_seq = "";
						
		// Iterate through each feature
		while((line = input.readLine()) != null){
			
			// Split on comma
			String[] line_arr = line.split(",");
			
			// Get the chromosome, position, and strand
			String oridb_chr = line_arr[1];
			int pos = Integer.parseInt(line_arr[2]);
			char strand = line_arr[3].charAt(0);
			
			// Update the chromsome if necessary
			if(!chr.equals(oridb_chr)){
				
				// Update the chr
				chr = oridb_chr;
				
				// Get the chromosome sequence
				chr_seq =  GenomeFunctions.read_in_chr_fasta_file(chr, fasta_dir);
						
			}
		
			// Get the sequence from -100 to 101
			if(strand == '+'){
						
				seq = chr_seq.substring(pos, pos + win);
						
			}else{
						
				// Adjust the position to the midpoint of the ACS
				seq = GenomeFunctions.reverse_complement_sequence(chr_seq.substring(pos - win + 1, pos + 1));
					
			}
	
			// Output the results
			output.write(line + "," + seq + "\n");
								
		}
		
		// Close the buffers
		input.close();
		output.close();

	}

}
