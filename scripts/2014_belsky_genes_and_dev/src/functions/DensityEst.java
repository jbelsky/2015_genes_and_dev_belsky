package functions;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;


public class DensityEst {
	
	public static void add_to_chr_array(double[] chr_array, double val, double pos, double bw){
		
		double normal_a = 1 / Math.sqrt(2 * Math.PI * Math.pow(bw, 2));
		double normal_b = -1 / (2 * Math.pow(bw, 2));
		
		for(int p = (int) Math.max(0, pos - 2 * bw);
			p <= (int) Math.min(chr_array.length - 1, pos + 2 * bw);
			p++){
			chr_array[p] += val * normal_a * Math.exp(normal_b * Math.pow(pos - p, 2));
		}	
		
	}
	
	public static double[] get_chip_seq_coverage_density(SAMFileReader bam_file, TF t, int win,
	        											 int shift, double bw, double read_val
														){

		// Set up the array around the ACS position
		double[] storage = new double[2 * win + 1];
		
		// Get the SAM Record Iterator
		SAMRecordIterator bam_itr = bam_file.queryOverlapping(t.getChr(),
								  							  t.getPos() - win - shift,
								  							  t.getPos() + win + shift
								 							 );
		
		// Get the array starting position
		int arr_start_pos = t.getPos() - win;
		
		// Set the read pos
		int read_mid = 0;
		
		while(bam_itr.hasNext()){
		
			// Get each read (the SAMRecord object has functions to get at the SAM attributes)
			SAMRecord r = bam_itr.next();
		
			// Get the read midpoint
			if(r.getReadNegativeStrandFlag()){
				read_mid = r.getAlignmentEnd() - shift;
			}else{
				read_mid = r.getAlignmentStart() + shift;
			}
			
			// Add the read to the array
			DensityEst.add_to_chr_array(storage, read_val, read_mid - arr_start_pos, bw);
		
		}
		
		// Close the bam iterator
		bam_itr.close();
		
		// Return the storage vector
		return(storage);
		
	}
	
	public static double[] get_mnase_coverage_density(SAMFileReader bam_file, TF t, int win,
													  int frag_low, int frag_high,
													  double bw, double read_val
													 ){

		// Set up the array around the ACS position
		double[] storage = new double[2 * win + 1];

		// Get the SAM Record Iterator
		SAMRecordIterator bam_itr = bam_file.queryOverlapping(t.getChr(),
				t.getPos() - 2 * win,
				t.getPos() + win
				);

		// Get the array starting position
		int arr_start_pos = t.getPos() - win;

		while(bam_itr.hasNext()){

			// Get each read (the SAMRecord object has functions to get at the SAM attributes)
			SAMRecord r = bam_itr.next();

			// Get the read midpoint
			int read_start = r.getAlignmentStart();
			int read_width = r.getInferredInsertSize();
			int read_mid = read_start + (read_width/2);
			
			// Check if read falls within the window coordinates
			if((read_width >= frag_low) && (read_width <= frag_high)){
			
				// Add the read to the array
				add_to_chr_array(storage, read_val, read_mid - arr_start_pos, bw);
				
			}

		}

		// Close the bam iterator
		bam_itr.close();

		// Return the storage vector
		return(storage);

	}
	
}
