package alignment;

import general.ExtendedWriter;
import general.Functions;
import general.RNAfunctions;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;


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

	public void printPremRNA( ExtendedWriter EW, int[] sequence){
		if(this.mRNAs != null){
			for(int i = 0; i < this.mRNAs.size();i++){
				mRNAs.get(i).printPremRNA(EW,sequence);
			}
		}
	}

	public void printCodingRNA( ExtendedWriter EW, int[] sequence,String contigName){
		if(this.mRNAs != null){
			for(int i = 0; i < this.mRNAs.size();i++){
				mRNAs.get(i).printCodingRNA(EW,sequence,contigName);
			}
		}
	}

	public void printmRNA( ExtendedWriter EW, int[] sequence,String contigName){
		if(this.mRNAs != null){
			for(int i = 0; i < this.mRNAs.size();i++){
				mRNAs.get(i).printmRNA(EW,sequence,contigName);
			}
		}
	}

	public void  printUpstream_5UTR_FirstIntron_Sequence(ExtendedWriter EW, int[] sequence, int  upstreamLength){
		if(this.mRNAs != null){
			for(int i = 0; i < this.mRNAs.size();i++){
				mRNAs.get(i).printUpstream_5UTR_FirstIntron_Sequence(EW,sequence,upstreamLength);
			}
		}
	}

	
	

	public void printmRNA(ExtendedWriter EW, int[] sequence,String contigName,Hashtable <Integer,StructuralVariation> SVs,String Sample){
		if(this.mRNAs != null){
			for(int i = 0; i < this.mRNAs.size();i++){
				mRNAs.get(i).printmRNA(EW,sequence,contigName,SVs,Sample);
			}
		}
	}

	public void printPersonalmRNAInfo( String ContigName,Hashtable <Integer,
			StructuralVariation> SVs, ArrayList<Integer> SVorder, ArrayList<String> samples,ExtendedWriter info){
		if(this.mRNAs != null){
			for(int i = 0; i < this.mRNAs.size();i++){
				mRNAs.get(i).printPersonalmRNAInfo(ContigName,SVs,SVorder,samples,info);
			}
		}

	}

	public void printPersonalmRNAInfo( String ContigName,Hashtable <Integer,
			StructuralVariation> SVs, ArrayList<Integer> SVorder, String sample,ExtendedWriter info){
		if(this.mRNAs != null){
			for(int i = 0; i < this.mRNAs.size();i++){
				mRNAs.get(i).printPersonalmRNAInfo(ContigName,SVs,SVorder,sample,info);
			}
		}

	}


	public void printPersonalmRNA(ExtendedWriter[] EWs, int[] sequence,String contigName,Hashtable <Integer,StructuralVariation> SVs,ArrayList<Integer> SVorder, String[] samples){
		if(this.mRNAs != null){
			for(int i = 0; i < this.mRNAs.size();i++){
				mRNAs.get(i).printPersonalmRNA(EWs,sequence,contigName,SVs,SVorder,samples);
			}
		}

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


	public boolean add5UTR(FUTR newExon, String parent){
		if(mRNAs != null){
			for(int i = mRNAs.size()-1;i > -1; i--){
				if(mRNAs.get(i).isParent(parent)){
					mRNAs.get(i).add5UTR(newExon);
					return true;
				}
			}
		}
		return false;
	}


	public boolean addCDS(CDS newExon,String parent){
		if(mRNAs != null){
			for(int i = mRNAs.size()-1;i > -1; i--){
				if(mRNAs.get(i).isParent(parent)){
					mRNAs.get(i).addCDS(newExon);
					return true;
				}
			}
		}
		return false;
	}

	public boolean add3UTR(TUTR newExon, String parent){
		if(mRNAs != null){
			for(int i = mRNAs.size()-1;i > -1; i--){
				if(mRNAs.get(i).isParent(parent)){
					mRNAs.get(i).add3UTR(newExon);
					return true;
				}
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