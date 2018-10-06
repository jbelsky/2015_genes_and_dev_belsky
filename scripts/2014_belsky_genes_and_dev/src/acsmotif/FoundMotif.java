package acsmotif;
import java.text.DecimalFormat;


public class FoundMotif implements Comparable<FoundMotif> {

	private String tf_name;
	private String chr;
	private int start;
	private int end;
	private char str;
	private double score;
	
	public FoundMotif(String t, String chromosome, int s, int e, char c, double sc){
		tf_name = t;
		chr = chromosome;
		start = s;
		end = e;
		str = c;
		score = sc;
	}

	public String getOutput(){
		
		DecimalFormat df = new DecimalFormat("#.####");
		
		return(tf_name + "," + 
			   chr + "," + 
			   getMotifStart() + "," +
			   str + "," +
			   start + "," + 
			   end + "," + 
			   df.format(score));
	}
	
	public String getTFName(){
		return(tf_name);
	}
	
	public int getStart(){
		return(start);
	}
	
	public int getEnd(){
		return(end);
	}
	
	public char getStr(){
		return(str);
	}
	
	public double getScore(){
		return(score);
	}

	public int getMotifStart(){
		if(str == '+')
			return(start);
		else
			return(end);
	}
	
	public int compareTo(FoundMotif o) {
		
		if(score > o.score){
			return(-1);
		}else if(score < o.score){
			return(1);
		}else{
			return(0);
		}
	}
	
}
