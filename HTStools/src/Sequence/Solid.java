package Sequence;

import java.io.Serializable;
import java.util.ArrayList;

import alignment.Hit;

import general.ExtendedWriter;
import general.RNAfunctions;


public class Solid extends Object implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public String Name;
	public String Sequence;
	public ArrayList <Hit> hits; 
	String infoLine;
	String SequenceLine;
	int[] tag;
	String lastTag;
	int[] colorCode;
	


	
	Solid(){
		this.Name = "temp";
	}

	Solid(String Name, String Sequence){
		this.Name = Name;
		String temp = Name.substring(1);
		String[] tempArray = temp.split("_");
		tag = new int[tempArray.length-1];
		for(int i = 0; i < tag.length; i++){
			tag[i]= Integer.parseInt(tempArray[i]);
		}
		this.lastTag = tempArray[tempArray.length-1];
		this.Sequence = Sequence;
		hits = new ArrayList<Hit>();
		
	}

	Solid(String Name){
		this.Name = Name;
		String temp = Name.substring(1);
		String[] tempArray = temp.split("_");
		tag = new int[tempArray.length-1];
		for(int i = 0; i < tag.length; i++){
			tag[i]= Integer.parseInt(tempArray[i]);
		}
		this.lastTag = tempArray[tempArray.length-1];
		
		hits = new ArrayList<Hit>();
	}	

	
	
	public void addSequence(FastaSequence Fseq){
		String seq = Fseq.getStringSequence();
		this.colorCode = RNAfunctions.Fasta2CFasta(seq);
		String colorSeq = seq.substring(0,1);
		for(int i = 0; i < colorCode.length; i++){
			colorSeq = colorSeq + colorCode[i];
		}
		this.Sequence = colorSeq;
	}
	
	
	public void mapHits(FastaSequences FS){
		if(hits.size() > 0){
			for(int i = 0; i< hits.size();i++){
				hits.get(i).mapHit(FS,1.0/(double)hits.size());
			}
		}
	}

	
	public int countHits(){
		if(hits.size() > 0){
			for(int i = 0; i< hits.size();i++){
				hits.get(i).nrOfHits = this.hits.size();
			}
		}
		return this.hits.size();
	}
	
	
	public void mapHits(FastaSequences FS, int exp){
		if(hits.size() > 0){
			for(int i = 0; i< hits.size();i++){
				hits.get(i).mapHit(FS, 1.0/(double)hits.size(), exp);
			}
		}
	}

	
	public boolean sameLocation(String chromosome, int location){
		for(int i = 0; i < hits.size(); i++){
			if(hits.get(i).sameLocation(chromosome,location))
				return true;
		}
		return false;
	}

	public void setColorCode(){
		if(this.colorCode == null && this.Sequence != null){
			this.colorCode = setColorCode(Sequence.substring(1));
		}
		if(this.Sequence == null)System.out.println(this.Name);
	}

	public static int[][] setColorCodes(String[] colorCodeStrings){
		int[][] colorCodes = new int[colorCodeStrings.length][];
		for(int i = 0; i < colorCodeStrings.length; i++){
			colorCodes[i] = setColorCode(colorCodeStrings[i]);
		}
		return colorCodes;
	}

	public static int[] setColorCode(String colorCodeString){
		char[] ccode = colorCodeString.toCharArray();
		int[] colorCode = new int[ccode.length];
		for(int i = 0; i < ccode.length; i++){
			colorCode[i] = (int)ccode[i]-48;
		}
		return colorCode;
	}
	

	public int removePrimer(CfastaSequences primers, double cutoff){
		float bestMatch = -1;
		int bestPrimer = -1;
		int nrOfPrimers = primers.size();
		for(int i = 0; i < nrOfPrimers;i++){
			
			int location = findPrimer(cutoff,0,10,2,primers.get(i).colorCode);
			if(location > -1){
				bestMatch = location;
				i = nrOfPrimers;
			}
		}
		if(bestMatch > -1)
			removePrimer((int)bestMatch);
		return this.colorCode.length;
		
	}
	
	
	public int findPrimer(double cutoff, int start, int maxLength, int minMatches ,int[] primer){
		
		int searchLength = Math.min(maxLength, primer.length);
		
		int missmatches = (int)((1-cutoff)*(double)searchLength);
		
		if(this.colorCode == null){
			setColorCode();
		}	
		int location  = start;
		while(location + minMatches < colorCode.length){
			int mm = 0;
			int m = 0;
			int pointer = 0;
			while(location+pointer < colorCode.length && pointer < searchLength && mm <= missmatches){
				if(colorCode[location+pointer] != primer[pointer]) mm++;
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
		int[] newColorCode = new int[position]; 
		if(this.colorCode == null){
			setColorCode();
		}	
		for(int i = 0; i < position; i++){
			newColorCode[i] = this.colorCode[i];
		}
		
		this.colorCode = newColorCode;
		
		this.Sequence = this.Sequence.substring(0,position);
		
	}

	
	
	
	
	public double compareTo(Solid otherSequence){
		double location = 0;
		int pointer = 0;
		while(location == 0 && pointer < tag.length){
			if(this.tag[pointer] == otherSequence.tag[pointer]) location = 0;
			else if(this.tag[pointer] < otherSequence.tag[pointer]) location = -1;
			else location = 1; 
			pointer++;
		}
		if(location != 0)
			return location;
		return this.lastTag.compareTo(otherSequence.lastTag);
	}
	
	
	public void addInfo(String infoLine, String SequenceLine){
		this.infoLine = infoLine;
		this.SequenceLine = SequenceLine;
	}
	
	public void printInfoEasy(ExtendedWriter EW){
		EW.println(this.infoLine);
		EW.println(this.SequenceLine);
	}
	
	public void printcFasta(ExtendedWriter EW){
		EW.print(Name);
		for(int i = 0; i < this.hits.size(); i++){
			EW.print(",");
			this.hits.get(i).print(EW);
		}
		EW.println();
		EW.println(this.Sequence);
	}

	public boolean printcFastaSize(ExtendedWriter EW, int length){
		
		if(this.colorCode.length >=length){
			printcFasta(EW);
			return true;
		}
		return false;
	}

	public boolean printcFastaSize(ExtendedWriter EW, int min, int max){
		
		if(this.colorCode.length >=min && this.colorCode.length < max){
			printcFasta(EW);
			return true;
		}
		return false;
	}
	
	
	public void printrMapperHits(ExtendedWriter EW){
		for(int i = 0; i < this.hits.size(); i++){
		EW.print(Name);
		this.hits.get(i).printRmapper(EW);
		//System.out.println(this.hits.get(i).nrOfReads);
		EW.println("\t"+this.Sequence);
		}
	}
	
	public void printRedundantrMapperHits(ExtendedWriter EW, int cutoff){
		for(int i = 0; i < this.hits.size(); i++){
			if(this.hits.get(i).nrOfReads > cutoff){
				EW.print(Name);
				this.hits.get(i).printRmapper(EW);
				System.out.println(this.hits.get(i).nrOfReads);
				EW.println("\t"+this.Sequence);
			}
		}
	}
	
	public String getFastaHitInfo(int hitNumber){
		return this.Name+","+this.hits.get(hitNumber).getFastaInfo();
	}
	
	
	
	public void printcFasta(ExtendedWriter EW, int kind){
		if(contains(kind)){
			EW.print(Name);
			for(int i = 0; i < this.hits.size(); i++){
				if(this.hits.get(i).getKind()== kind){
				EW.print(",");
				this.hits.get(i).print(EW);
				}
		}
		EW.println();
		EW.println(this.Sequence);
		}
	}


	public int[] countKinds(int[] count){

		/*	1=mRNA
			2=ncRNA
			3=intergenic
			4=antisense
			5=repeats
			6=DIRS
		 */
		
		for(int i = 0; i < this.hits.size(); i++){
			count[this.hits.get(i).getKind()]++;
		}
		return count;	
	}
	
	public double[] countNormalizedKinds(double[] count){

		/*	1=mRNA
			2=ncRNA
			3=intergenic
			4=antisense
			5=repeats
			6=DIRS
		 */
		
		for(int i = 0; i < this.hits.size(); i++){
			count[this.hits.get(i).getKind()] =count[this.hits.get(i).getKind()] +((double)1/this.hits.size());
		}
		return count;	
	}


	
	
	private boolean contains(int kind){
		for(int i = 0; i < this.hits.size(); i++){
			if(this.hits.get(i).getKind() == kind)
				return true;
		}
		return false;
	
		
	}
	
		
		
	
	
	public void addHit(int start, int length, int missmatches, boolean plusStrand){
		this.hits.add(new Hit(start,length,missmatches,plusStrand));
	}

	public void addHit(int start, int length, int missmatches, boolean plusStrand,String chromosome){
		this.hits.add(new Hit(start,length,missmatches,plusStrand, chromosome));
	}
	
	public void addHit(int start, int length, boolean plusStrand,String chromosome, String editString, int score){
		if(hits.size() > 0)
			if(this.hits.get(0).score < score){
				this.hits = new ArrayList<Hit>();
			}
			else if(this.hits.get(0).score > score)
				return;
		this.hits.add(new Hit(start,length,plusStrand, chromosome,editString,score,this));
	}

	
	
	public void addHit(Hit newHit){
		if(hits.size() > 0)
			if(this.hits.get(0).score < newHit.score){
				this.hits = new ArrayList<Hit>();
			}
			else if(this.hits.get(0).score > newHit.score)
				return;
		newHit.solidSequence = this;
		this.hits.add(newHit);
	}
	
	
	public void addHits(Solid cfasta){
		ArrayList<Hit> moreHits= cfasta.hits;
		for(int i = 0;i < moreHits.size();i++)
			this.hits.add(moreHits.get(i));
	}
	
	
	

	
}
