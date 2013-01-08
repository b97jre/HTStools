
package alignment;

import java.io.FileReader;
import java.io.Serializable;
import java.util.Hashtable;

import extra.ExtendedReader;

import general.Functions;

public class GenesHT extends Hashtable<String,Gene>  implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	
	public GenesHT(){}
	
	
	public void addGene(Gene G){
		if(!containsGene(G.getName()))
			this.put(G.getName(),G);
		else{
			System.out.println("Gene "+G.getName() +" already exist in hashtable");
		}
	}

	
	public Gene getGene(String GeneName){
		if(this.containsKey(GeneName)) return this.get(GeneName) ;
		return null;
	}

	public boolean containsGene(String GeneName){
		if(this.containsKey(GeneName)) return true;
		return false;
		
	}
	
	
	public String[] getGeneNames(){
		String[] Names = new String[this.size()];
		for(int i = 0; i < this.size(); i++){
			Names[i] = this.get(i).getName();
		}
		return Names;
	}
	
	
	public void parseFpkm_tracking_file(String dir,String file){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(dir+"/"+file));
			ER.readLine();// remove first row;
			
		}catch(Exception E){
			E.printStackTrace();
		}
		
	}

	
}
	
	