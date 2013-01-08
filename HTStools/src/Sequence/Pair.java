package Sequence;

import general.ExtendedWriter;

public class Pair {
	
	FastaSequences ref1;
	String end1;
	String end2;
	FastaSequences ref2;
	int code;
	int nrOfHits;
	
	
	
	
	Pair(FastaSequences ref1, String end1, FastaSequences ref2,String end2, int code){
		this.ref1 = ref1;
		this.end1 = end1;
		this.ref2 = ref2;
		this.end2 = end2; 
		this.code = code;
		this.nrOfHits = 1;
		
	}
	public boolean isTheSame(FastaSequence ref1, FastaSequence ref2){
		if(this.ref1.Name.compareTo(ref1.Name) == 0 && this.ref2.Name.compareTo(ref2.Name) == 0) return true;
		if(this.ref1.Name.compareTo(ref2.Name) == 0 && this.ref2.Name.compareTo(ref1.Name) == 0) return true;
		return false;
	}
	
	public void addHit(){
		nrOfHits++;
	}
	
	public void printPair(){
		System.out.println(this.ref1.Name+"\t"+this.ref2.Name+"\t"+this.code+"\t"+this.nrOfHits);
	}
	
	public void printPair(ExtendedWriter EW){
		EW.println(this.ref1.Name+"\t"+this.ref2.Name+"\t"+this.code+"\t"+this.nrOfHits);
	}
	
	
	
	
}
