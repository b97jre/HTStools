package alignment;

import java.io.Serializable;
import java.util.ArrayList;


public class Intergenic extends Gene implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	
	protected ArrayList <Hit> hits;
	
	Intergenic(){}
	
	Intergenic(int left, int right,String ID,String Name,String description){
		this.ID = ID;
		this.Name = Name;
		this.left = left;
		this.right = right;
		this.description = description;
		kind = this.intergenic;
	}
	
	public void setInfo(Gene otherGene){
		this.ID = otherGene.ID;
		this.Name = otherGene.Name;
		this.plusStrand = otherGene.plusStrand;
		this.left = otherGene.left;
		this.right = otherGene.right;
		this.description = otherGene.description;
		kind = this.intergenic;

	}
	
}