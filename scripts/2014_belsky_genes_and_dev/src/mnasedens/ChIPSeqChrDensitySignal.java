package mnasedens;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import functions.BAMInput;
import net.sf.samtools.SAMFileReader;

public class ChIPSeqChrDensitySignal{

	public static void main(String[] args) throws IOException {

		// Set the bam file name
		String bam_file_name = args[0];
	
		// Set up the output file name
		String output_file_name_header = args[1];
		
		// Set the output file name distribution
		String output_distribution_file_name = args[2];
		
		// Set the fragment length thresholds
		int shift = Integer.parseInt(args[3]);
				
		// Set the bw
		double bw = Double.parseDouble(args[4]); 
				
		//////////////////////////////////////////////////////////////////////////////////////////////
		
		// Get the output distribution
		BufferedWriter dist_output = new BufferedWriter(new FileWriter(output_distribution_file_name));
					
		// Get the bam file
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));
			
		// Iterate through each chr
		for(int chr = 1; chr <= 16; chr++){
			
			// Set the chromosome, start, and end positions
			int start_pos = 0;
			int end_pos = BAMInput.get_chr_length(bam_file_name, Integer.toString(chr));
				
			System.out.println(end_pos);
			
			// Get the coverage
			double[] cov = BAMInput.read_in_sam_file_density(bam_file_name, Integer.toString(chr), 
															 start_pos, end_pos, shift, bw
															);
	
			// Write the output
			GetMNaseChrDensity.write_chr_output(output_file_name_header + chr + ".csv",
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
