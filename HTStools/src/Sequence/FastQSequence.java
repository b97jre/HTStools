package Sequence;

import java.io.Serializable;
import java.util.ArrayList;

import general.ExtendedReader;
import general.ExtendedWriter;
import general.RNAfunctions;


public class FastQSequence extends Object implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	String name;
	String sequence;
	String quality;
	String otherName;
	
	
//	@HWI-ST167:1:1101:20376:1996#0/1
//	TAGTGATAACCATCCTATTCATTTTGCAAAGTTGATATATTACTTTCTGTTCATATATTTGTATTTGGTTGGACTTCTGTGAACATGAAAACACCGAATTC
//	+
//	bbbeeeeegggggiiiiiiiiiiihiiiiiiiihhiiiiiiiiiiiiiiiiiiihhiihiiiihihhihiihhiighiihihhhhfcggbddeeebccccb
	
	
	public double[] nrOfHits;
	
	

	public FastQSequence(){
		this.name = "not present";
		this.sequence = null;
		this.otherName = null;
		this.quality = null;
		this.nrOfHits = null;
		
	}
	
	
	FastQSequence(String Name, String Sequence,String otherName, String Quality, int exp){
		this.name = Name;
		this.sequence = Sequence;
		this.otherName = otherName;
		this.quality = Quality;
		this.nrOfHits = new double[exp];
	}
	

	public boolean addInfo(ExtendedReader ER){
		this.name = ER.readLine();
		this.sequence = ER.readLine();
		this.otherName = ER.readLine();
		this.quality = ER.readLine();
		return this.QC();
	}

	public boolean addInfo2(ExtendedReader ER){
		this.name = ER.readLine();
		this.sequence = ER.readLine();
		this.otherName = ER.readLine();
		this.quality = ER.readLine().substring(0,this.sequence.length());
		return this.QC();
	}
	
	
	
	public boolean addInfo(ExtendedReader ER, int length){
		this.name = ER.readLine();
		this.sequence = ER.readLine();
		this.otherName = ER.readLine();
		this.quality = ER.readLine();
		return this.QC(length);
	}
	
	
	public boolean QC(){
		if(this.sequence == null) return false;
		if(this.quality == null) return false;
		if(this.sequence.length() != this.quality.length()){
			System.out.println(this.name+ " does not have the same sequence and quality length or sequenceLength is less than 40");
			return false;
		}
		
		return true;
	}
	
	public boolean QC(int length){
		if(this.sequence == null) return false;
		if(this.quality == null) return false;
		if(this.sequence.length() != this.quality.length()){
			System.out.println(this.name+ " does not have the same sequence and quality length");
			return false;
		}
		if(this.sequence.length() < length) return false;
		return true;
	}


	
	public boolean isPair(FastQSequence reverse){
		String[] Name1 = this.name.split("/");
		String[] Name2 = reverse.name.split("/");
		if(Name1[0].compareTo(Name2[0]) == 0) return true;
		return false;
	}
	
	
	
	
	
	
	public void printFasta(ExtendedWriter EW){
		EW.println(">"+this.name.substring(1));
		EW.println(this.sequence);
	}
	
	public void printFasta(ExtendedWriter EW, String ExtraName){
		EW.println(">"+ExtraName+" "+this.name.substring(1));
		EW.println(this.sequence);
	}

	public void printFastQ(ExtendedWriter EW){
		EW.println(this.name);
		EW.println(this.sequence);
		EW.println(this.otherName);
		EW.println(this.quality);
	}
	
	public void printFastQreverse(ExtendedWriter EW){

		String[] Name1 = this.name.split("/");
		if(Name1[1].contains("1")) 
			EW.println(Name1[1]+"/2");
		else
			EW.println(Name1[1]+"/1");
		EW.println(this.sequence);
		EW.println(this.otherName);
		EW.println(this.quality);
	}

	
	public void printFastQ(ExtendedWriter EW, int min){
		if(this.sequence.length() >= min){
			EW.println(this.name);
			EW.println(this.sequence);
			EW.println(this.otherName);
			EW.println(this.quality);
		}
	}
	
	public void printFastQ(ExtendedWriter EW, int min, int max){
		if(this.sequence.length() >= min && this.sequence.length() <= max){
			EW.println(this.name);
			EW.println(this.sequence);
			EW.println(this.otherName);
			EW.println(this.quality);
		}
	}
	
	
	public void printNrOfHits(ExtendedWriter EW){
		EW.print(name +"\t"+sequence.length());
		for(int i =0 ; i < nrOfHits.length; i++)	
			EW.print("\t"+nrOfHits[i]);
		EW.println();
	}
	
	public void addHit(double nrOfHits){
		this.nrOfHits[0] += nrOfHits;
	}
	
	public void addHit(double nrOfHits,int exp){
		this.nrOfHits[exp] += nrOfHits;
	}
	
	public String getSequence() {
		return sequence;
	}

	public int[] getIntSequence() {
		return RNAfunctions.RNAString2Int(sequence);
	}
	public void setSequence(int[] sequence) {
		this.sequence = RNAfunctions.RNAInt2String(sequence);
	}
	
	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public int[][] getORFinfo(){
		return FastaSequence.findLongestORFinfo(this.sequence);
	}
	
	public int removePrimer(FastaSequences primers, double cutoff){
		float bestMatch = -1;
		int bestPrimer = -1;
		int nrOfPrimers = primers.size();
		for(int i = 0; i < nrOfPrimers;i++){
			
			int location = findPrimer(cutoff,0,10,2,primers.get(i).Sequence);
			if(location > -1){
				bestMatch = location;
				i = nrOfPrimers;
			}
		}
		if(bestMatch > -1)
			removePrimer((int)bestMatch);
		return this.sequence.length();
	}
	
	
	
	
	public int findPrimer(double cutoff, int start, int maxLength, int minMatches ,int[] primer){
		
		int searchLength = Math.min(maxLength, primer.length);
		
		int missmatches = (int)((1-cutoff)*(double)searchLength);
		
		int location  = start;
		int[] Sequence = RNAfunctions.RNAString2Int(this.sequence);
		while(location + minMatches < Sequence.length){
			int mm = 0;
			int m = 0;
			int pointer = 0;
			while(location+pointer < Sequence.length && pointer < searchLength && mm <= missmatches){
				if(Sequence[location+pointer] != primer[pointer]) mm++;
				else m++;
				pointer++;
			}
			if(((double)m/(double)(m+mm)) >= cutoff ) 
				return location-1;
			location++;
		}
		return -1;
	}
	
	
	
	public void removePrimer(int position){
		this.sequence = this.sequence.substring(0,position);
		this.quality = this.quality.substring(0, position);
	}
	
	
// Not finished yet Using SeqPrep instead.
	public FastQSequence joinPairs(FastQSequence PairedEnd, int length,int start){
		FastQSequence joined;
		int[] seq = RNAfunctions.RNAString2Int(this.sequence);
		int[] seq2 = RNAfunctions.getReverseComplement(RNAfunctions.RNAString2Int(PairedEnd.sequence));
		
		
//		int joinedSeqLocation = joinSequences(seq, seq2, length,start);
		
		
		return null;
	}
	
	private int findJoinedSequences(int[] seq, int[] seq2, int length, int start){
		boolean found = false;
		int position = start;
		while(position < seq.length-length && found ==  false){
			int count = 0;
			boolean same = true;
			while(same && position+count < seq.length){
				if(seq[position+count] == seq2[count])
					count++;
				else
					same = false;
			}
			if(same){
				return position;
			}
			position++;
		}
		return -1;
	}
	
	
	
	// Both arugemnts below are not tested and are not ready without testing.
	private String joinSequences(int position, int[] seq, int[] seq2){ 				
		int[] joinedSequences = new int[position+seq2.length];
		for(int i = 0; i < position;i++){
			joinedSequences[i] = seq[i];
		}	
		for(int i = 0; i < seq2.length;i++){
			joinedSequences[i+position] = seq2[i];
		}	
		return RNAfunctions.DNAInt2String(joinedSequences);
	}
	
	private String joinQuality(int position, char[] Quality, char[] otherQuality){
		char[] joinedQuality = new char[position+otherQuality.length];
		for(int i = 0; i < position; i++){
			joinedQuality[i] = Quality[i];
		}
		for(int i= position; i < Quality.length; i++){
			joinedQuality[i] = 'h';
		}
		for(int i = 0; i < otherQuality.length- (Quality.length-position);i++){
			joinedQuality[Quality.length+i] = otherQuality[otherQuality.length-Quality.length-position-i];
			
		}
		return new String(joinedQuality);
		
	}



	

}
