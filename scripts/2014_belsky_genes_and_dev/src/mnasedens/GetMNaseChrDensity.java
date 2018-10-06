package mnasedens;
import functions.BAMInput;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Random;

import net.sf.samtools.SAMFileReader;

public class GetMNaseChrDensity {

	public static void write_chr_output(String output_name, BufferedWriter dist_output, 
										double[] sig, String chr_name) throws IOException{
		
		// Set the output format
		DecimalFormat df = new DecimalFormat("#.####");
				
		// Create the random
		Random r = new Random();
		
		// Write the output
		BufferedWriter output = new BufferedWriter(new FileWriter(output_name));
		
		// Write the header
		output.write("chr,pos,signal\n");
		
		// Iterate through the positions to write the output
		for(int p = 1; p < sig.length; p++){

			// Write the output
			output.write(chr_name + "," + 
						 p + "," + 
						 df.format(sig[p]) + "\n"
						);
			
			if(sig[p] > 0){
				
				if(r.nextDouble() < 0.1){
					
					dist_output.write(df.format(sig[p]) + "\n");
					
				}
					
			}
			
		}
		
		// Close the output buffer
		output.close();
		
	}
	
	
	
	public static void main(String[] args) throws IOException {

		// Set the bam file name
		String bam_file_name = args[0];
	
		// Set up the output file name
		String output_file_name_header = args[1];
		
		String output_distribution_file_name = args[2];
		
		// Set the fragment length thresholds
		int frag_length_low = Integer.parseInt(args[3]);
		int frag_length_high = Integer.parseInt(args[4]);
		
		// Set the bw
		double bw = Double.parseDouble(args[5]); 
				
		//////////////////////////////////////////////////////////////////////////////////////////////
		
		// Get the output distribution
		BufferedWriter dist_output = new BufferedWriter(new FileWriter(output_distribution_file_name));
					
		// Get the bam file
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));
			
		// Iterate through each chr
		for(int chr = 1; chr <= 16; chr++){
			
			// Get the coverage
			double[] cov = 
				BAMInput.read_in_paired_sam_file_density(bam_file_name, Integer.toString(chr), 
														 bw, frag_length_low, frag_length_high
														);
			
			// Write the output
			write_chr_output(output_file_name_header + chr + ".csv",
							 dist_output, cov, Integer.toString(chr)
							);
			
		}
			
		// Close the bam file
		bam_file.close();
		
		// Close the dist output
		dist_output.close();
		
		System.out.println("\tComplete!");
				
	}

}
