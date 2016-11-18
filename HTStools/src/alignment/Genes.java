package alignment;

import java.io.Serializable;
import java.util.ArrayList;

import general.Functions;

public class Genes extends ArrayList<Gene>  implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	
	public Genes(){}
	
	
	public void setGene(){
		Gene A = this.get(1);
		add(A);
		
	}
	
	
	
	public int[] findGene(String GeneName){
		int[] positions = null;
		if(this.size() == 0) return null;
		for(int i = 0; i < this.size();i++){
			if(GeneName.compareTo(this.get(i).getName()) == 0)
				positions = Functions.addInt(positions,i);
		}
		return positions;
	}
	
	public Gene[] getGenes(String GeneName){
		int[] geneNr = findGene(GeneName);
		if(geneNr == null) return null;
		Gene[] genes = new Gene[geneNr.length];
		for(int i = 0; i < geneNr.length;i++)
			genes[i] = this.get(geneNr[i]);
		return genes;
	}
	
	
	public String[] getGeneNames(){
		String[] Names = new String[this.size()];
		for(int i = 0; i < this.size(); i++){
			Names[i] = this.get(i).getName();
		}
		return Names;
	}
}
	
	
	
	