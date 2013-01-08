package alignment;

import java.io.Serializable;
import java.util.ArrayList;

public class GeneInfo implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	ArrayList <String> GOterms;
	String location;
	
	public void addGoTerm(String GOterm){
		if(GOterms == null) GOterms = new ArrayList<String>();
		GOterms.add(GOterm);
	}
	
	public void containsGOTerm(String GOterm){
		GOterms.contains(GOterm);
	}
	
	public ArrayList<String> getGOterms(){
		return this.GOterms;
	}
	
	public void setLocaiton(String location){
		this.location = location;
	}
	
	public String getLocation(){
		return this.location;
	}
	
	public String getGeneInfo(){
		String info = "";
		if(location != null)
			info = location;
		if(GOterms != null && GOterms.size()> 0){
			info +="   GOTERMS: ";
		for(int i = 0; i < GOterms.size(); i++){
			info = info +"  "+GOterms.get(i);
			
		}
		}
		return info;
		
	}
	
	
	
}
