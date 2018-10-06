package acsmotif;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;

public class ARS implements Comparable<ARS> {

	protected String oridb_name;
	protected String chr;
	protected int start;
	protected int end;
	protected String ars_name;
	protected String confidence;
			
	public ARS(String n, String c, int st, int en, String ars_n, String type){
		oridb_name = n;
		ars_name = ars_n;
		chr = c;
		start = st;
		end = en;
		confidence = type;
	}
	
	// Return functions
	public String getName(){
		return(oridb_name);
	}
	
	public String getChr(){
		return(chr);
	}
		
	public int getStart(){
		return(start);
	}
	
	public int getEnd(){
		return(end);
	}
	
	public String getARSName(){
		return(ars_name);
	}
	
	public String getConfidence(){
		return(confidence);
	}
	
	// Static Functions
	public static HashMap<String, ArrayList<ARS>> read_in_tf_map(String filename) 
			throws IOException{
		
		// Get the Buffer for the filename and read in the header
		BufferedReader input = new BufferedReader(new FileReader(filename));
		
		// Read in the header
		String line = input.readLine();
				
		// Split the header
		String[] line_arr = line.split(",");
		
		// Find the indices corresponding to ars_start, ars_end, and oridb_status
		int ars_start_idx = -1;
		int ars_end_idx = -1;
		int ars_name_idx = -1;
		int oridb_status_idx = -1;

		for(int i = 0; i < line_arr.length; i++){
			if(line_arr[i].equals("ars_start")){
				ars_start_idx = i;
			}else if(line_arr[i].equals("ars_end")){
				ars_end_idx = i;
			}else if(line_arr[i].equals("ars_name")){
				ars_name_idx = i;
			}else if(line_arr[i].equals("oridb_status")){
				oridb_status_idx = i;
			}
		}
		
		// Set up the HashMap
		HashMap<String, ArrayList<ARS>> ars_map = new HashMap<String, ArrayList<ARS>>();
				
		while((line = input.readLine()) != null){
			
			// Split the line
			line_arr = line.split(",");
			String chr = line_arr[1];
			
			if(!ars_map.containsKey(chr)){
				ArrayList<ARS> acs_list = new ArrayList<ARS>();
				ars_map.put(chr, acs_list);
			}
			
			// Get the ARS characteristics
			String ars_name = line_arr[0];
			int start = Integer.parseInt(line_arr[ars_start_idx]);
			int end = Integer.parseInt(line_arr[ars_end_idx]);
			String arsName = line_arr[ars_name_idx];
			String ars_confidence = line_arr[oridb_status_idx];
			
			// Create a new ARS object
			ARS ars_obj = new ARS(ars_name, chr, start, end, arsName, ars_confidence);
						
			// Enter the ARS into the map
			ars_map.get(chr).add(ars_obj);						
		
		}
		
		// Close the Buffer
		input.close();
		
		return(ars_map);
		
	}
	

	
	static Comparator<ARS> motifCompare = new Comparator<ARS>(){

		public int compare(ARS t1, ARS t2) {
			
			if(t1.getStart() < t2.getStart())
				return -1;
			else if(t1.getStart() > t2.getStart())
				return 1;
			else
				return 0;
		}
		
	};
	
	static Comparator<Object> posCompare = new Comparator<Object>(){

		public int compare(Object o1, Object o2) {
			ARS tf = (ARS) o1;
			int pos = (Integer) o2;
			
			if(tf.getStart() < pos){
				return(-1);
			}else if(tf.getStart() > pos){
				return(1);
			}else{
				return(0);
			}
			
		}
		
	};

	@Override
	public int compareTo(ARS o) {
		// TODO Auto-generated method stub
		return(this.start - o.getStart());
	}
	
	
}