package mnasedens;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import functions.TF;

public class GetFragmentDistributionOverFeature {

	public static void main(String[] args) throws IOException {

		////////////////////////////////////////////////////////////////////////////////////////
		// Enter the input parameters
		
		// Set the input feature file
		String input_feature_file = args[0];
		
		// Set the bam file name
		String input_bam_file = args[1];
		
		// Set the output file name
		String output_file_name = args[2];
		
		// Enter the positions into the array
		int[] pos_array = new int[2];
		pos_array[0] = Integer.parseInt(args[3]);
		pos_array[1] = Integer.parseInt(args[4]);
		
		////////////////////////////////////////////////////////////////////////////////////////
		// Data Processing
		
		// Read in the feature file
		ArrayList<TF> tf_list = TF.read_in_tf_list(input_feature_file);
	
		// Set up the fragment length distribution
		int[][] frag_dist = new int[1001][2];
		
		// Open the BAM File
		SAMFileReader bam_file = new SAMFileReader(new File(input_bam_file), new File(input_bam_file + ".bai"));
		
		// Iterate through each TF
		for(TF t : tf_list){
			
			// Get the feature position
			int feature_pos = t.getPos();
			
			// Set the multiplier
			int mult = 1;
			if(t.getStrand() == '-'){
				mult = -1;
			}
			
			// Iterate through each pos array
			for(int i = 0; i < 2; i++){
				
				// Get the position
				int pos = feature_pos + mult * pos_array[i];
				
				// Get the bam_itr
				SAMRecordIterator bam_itr = bam_file.queryOverlapping(t.getChr(), pos - 250, pos);
			
				// Iterate through each record
				while(bam_itr.hasNext()){
				
					// Get the SAMRecord
					SAMRecord read = bam_itr.next();
				
					// Get the read start and end coordinates
					int start = read.getAlignmentStart();
					int width = read.getInferredInsertSize();
					int end = start + width - 1;
				
					// If the read overlaps the position, enter into the frag_dist
					if((start <= pos) && (end >= pos)){
						frag_dist[width][i]++;
					}
					
				}
			
				// Close the itr
				bam_itr.close();
								
			}
			
		}
		
		// Close the bam file
		bam_file.close();
		
		// Write the output
		BufferedWriter output = new BufferedWriter(new FileWriter(output_file_name));
		
		// Create the header
		output.write("frag_length," + pos_array[0] + "," + pos_array[1] + "\n");
		
		// Iterate through the frag_length
		for(int i = 20; i < frag_dist.length; i++){
			output.write(i + "," + frag_dist[i][0] + "," + frag_dist[i][1] + "\n");
		}
		
		// Close the output buffer
		output.close();
		
	}

}
