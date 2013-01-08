package alignment;

import java.io.Serializable;
import java.util.ArrayList;


public class CodingGene extends Gene implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	protected ArrayList <mRNA> mRNAs;
	
	CodingGene(){
		kind = this.mRNA;

	}
	
	CodingGene(int left, int right,boolean plusStrand,String ID,String Name,String description){
		this.ID = ID;
		this.Name = Name;
		this.plusStrand = plusStrand;
		this.left = left;
		this.right = right;
		this.description = description;
		kind = this.mRNA;
		}
	
	CodingGene(int left, int right,String ID,String Name,String description){
		this.ID = ID;
		this.Name = Name;
		this.left = left;
		this.right = right;
		this.description = description;
		kind = this.mRNA;
		}
	
	
	public int getKind(boolean plusStrand){
		return kind;
	}
		
			
	public void print(){
		if(plusStrand){
			System.out.println(ID+"\t"+left+"\t"+right+"\t+\t"+description);
		}
		else
			System.out.println(ID+"\t"+left+"\t"+right+"\t-\t"+description);
	}
	

	
	
	public boolean addmRNA(mRNA newmRNA){
		if(newmRNA.parent.compareTo(ID) == 0){
			if(this.mRNAs == null)
				this.mRNAs = new ArrayList<mRNA>();
			this.mRNAs.add(newmRNA);
			return true;
		}
		return false;
	}
	
	public boolean addExon(Exon newExon){
		if(mRNAs != null){
		for(int i = mRNAs.size()-1;i > -1; i--){
			if(mRNAs.get(i).addExon(newExon)) return true;
		}
		}
		return false;
	}

	public ArrayList<mRNA> getMRNAs() {
		return mRNAs;
	}

	public void setMRNAs(ArrayList<mRNA> as) {
		mRNAs = as;
	}


}