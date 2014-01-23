package Sequence;

import general.ExtendedWriter;

public class Pair {
	
	String refName1;
	String end1;
	String end2;
	String refName2 ;
	int code;
	int nrOfHits;
	
	
	
	
	Pair(FastaSequences ref1, String end1, FastaSequences ref2,String end2, int code){
		this.refName1 = ref1.Name;
		this.end1 = end1;
		this.refName2 = ref2.Name;
		this.end2 = end2; 
		this.code = code;
		this.nrOfHits = 1;
		
	}
	public boolean isTheSame(FastaSequence ref1, FastaSequence ref2){
		if(this.refName1.compareTo(ref1.Name) == 0 && this.refName2.compareTo(ref2.Name) == 0) return true;
		if(this.refName1.compareTo(ref2.Name) == 0 && this.refName2.compareTo(ref1.Name) == 0) return true;
		return false;
	}
	
	public void addHit(){
		nrOfHits++;
	}
	
	public void printPair(){
		System.out.println(this.refName1+"\t"+this.refName2+"\t"+this.code+"\t"+this.nrOfHits);
	}
	
	public void printPair(ExtendedWriter EW){
		EW.println(this.refName1+"\t"+this.refName2+"\t"+this.code+"\t"+this.nrOfHits);
	}
	
	
	
	
}
