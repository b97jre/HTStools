package Blast;

import alignment.Chromosome;
import alignment.Gene;
import Sequence.FastaSequences;
import Sequence.Solid;
import general.ExtendedWriter;

public class BlastHit {
	//trinity_4_1_1	trinity_4_1_1	100.00	79	0	0	1	79	1	79	1e-27	 118
	protected String queryName;
	protected String hitName;
	protected double similarity;
	protected int length;
	protected int missmatches;
	protected int gaps;
	protected int queryStart;
	protected int queryStop;
	protected int hitStart;
	protected int hitStop;
	protected double Evalue;
	protected double score;
	protected Gene specificGene;
	public int nrOfHits;



	public BlastHit(int score, Gene query, String HitName){
		this.score = score;	
		this.specificGene = query;
		this.hitName = HitName;
	}


	public BlastHit(
			String	queryName,
			String hitName,
			double similarity,
			int length,
			int missmatches,
			int gaps,
			int queryStart,
			int queryStop,
			int hitStart,
			int hitStop,
			double Evalue,
			double score
	){
		this.hitName = hitName;
		this.queryName =queryName;
		this.similarity = similarity;
		this.length = length;
		this.missmatches = missmatches;
		this.gaps = gaps;
		this.queryStart = queryStart;
		this.queryStop = queryStop;
		this.hitStart =hitStart;
		this.hitStop = hitStop;
		this.Evalue = Evalue;
		this.score = score;
	}


	public void printHit(ExtendedWriter EW){

		EW.print(queryName);
		EW.print("\t"+hitName);
		EW.print("\t"+similarity);
		EW.print("\t"+length);
		EW.print("\t"+missmatches);
		EW.print("\t"+gaps);
		EW.print("\t"+queryStart);
		EW.print("\t"+queryStop);
		EW.print("\t"+hitStart);
		EW.print("\t"+hitStop);
		EW.print("\t"+Evalue);
		EW.println("\t"+score);

	}

	public void printBestHit(ExtendedWriter EW,double totalScore){

		EW.print(queryName);
		EW.print("\t"+hitName);
		EW.print("\t"+score);
		EW.println("\t"+(score/totalScore));

	}

	public void printBestLength(ExtendedWriter EW,double bestLength){

		EW.print(queryName);
		EW.print("\t"+hitName);
		EW.print("\t"+length);
		EW.println("\t"+((double)length/bestLength));

	}
	
	
	public void printNoHit(ExtendedWriter EW,double totalScore){

		EW.print(queryName);
		EW.print("\tNA");
		EW.print("\t0");
		EW.println("\t0");

	}
	
	

	public void printHit(){

		System.out.print(queryName);
		System.out.print("\t"+hitName);
		System.out.print("\t"+similarity);
		System.out.print("\t"+length);
		System.out.print("\t"+missmatches);
		System.out.print("\t"+gaps);
		System.out.print("\t"+queryStart);
		System.out.print("\t"+queryStop);
		System.out.print("\t"+hitStart);
		System.out.print("\t"+hitStop);
		System.out.print("\t"+Evalue);
		System.out.println("\t"+score);

	}

	
	public static BlastHit reversed(BlastHit BH){
		BlastHit rev = new BlastHit(BH.hitName,BH.queryName,BH.similarity,BH.length,BH.missmatches,BH.gaps,BH.queryStart,BH.queryStop,BH.hitStart, BH.hitStop, BH.Evalue, BH.score);
		return rev;

	}

	public void merge(BlastHit BH, int penalty, int length){
//		System.out.println("Before");
//		this.printHit();
//		BH.printHit();		
//		System.out.println();

		this.length = this.length+ BH.length+length;
		this.missmatches = missmatches + BH.missmatches;
		this.gaps = gaps + BH.gaps+length;
		this.similarity = (double)(this.length-this.missmatches-this.gaps )/(double)this.length*100;
		this.Evalue = this.Evalue*BH.Evalue;
		this.score = score+BH.score-length*penalty;

		if(this.hitStart < this.hitStop){
			if (this.hitStart < BH.hitStart){
				this.hitStop = BH.hitStop;
				this.queryStop = BH.queryStop;
			}
			else{
				this.hitStart = BH.hitStart;
				this.queryStart = BH.queryStart;
			}
		}else{
			if (this.hitStart < BH.hitStart){
				this.hitStart = BH.hitStart;
				this.queryStart = BH.queryStart;
			}
			else{
				this.hitStop = BH.hitStop;
				this.queryStop = BH.queryStop;
			}

		}
//		System.out.println("After");
//		this.printHit();
//		System.out.println();
	}

	public boolean isForward(){
		if(this.hitStop>this.hitStart){
			if(this.queryStop>this.queryStart)return true;
			return false;
		}
		if(this.hitStop<this.hitStart){
			if(this.queryStop<this.queryStart)return true;
			return false;
		}
		System.out.println("this should never happen in isForward()");
		return false;

	}


	public String getQueryName() {
		return queryName;
	}


	public void setQueryName(String queryName) {
		this.queryName = queryName;
	}


	public String getHitName() {
		return hitName;
	}


	public void setHitName(String hitName) {
		this.hitName = hitName;
	}


	public double getSimilarity() {
		return similarity;
	}


	public void setSimilarity(double similarity) {
		this.similarity = similarity;
	}


	public int getLength() {
		return length;
	}


	public void setLength(int length) {
		this.length = length;
	}


	public int getMissmatches() {
		return missmatches;
	}


	public void setMissmatches(int missmatches) {
		this.missmatches = missmatches;
	}


	public int getGaps() {
		return gaps;
	}


	public void setGaps(int gaps) {
		this.gaps = gaps;
	}


	public int getQueryStart() {
		return queryStart;
	}


	public void setQueryStart(int queryStart) {
		this.queryStart = queryStart;
	}


	public int getQueryStop() {
		return queryStop;
	}


	public void setQueryStop(int queryStop) {
		this.queryStop = queryStop;
	}


	public int getHitStart() {
		return hitStart;
	}


	public void setHitStart(int hitStart) {
		this.hitStart = hitStart;
	}


	public int getHitStop() {
		return hitStop;
	}


	public void setHitStop(int hitStop) {
		this.hitStop = hitStop;
	}


	public double getEvalue() {
		return Evalue;
	}


	public void setEvalue(double evalue) {
		Evalue = evalue;
	}


	public double getScore() {
		return score;
	}


	public void setScore(double score) {
		this.score = score;
	}


	public Gene getSpecificGene() {
		return specificGene;
	}


	public void setSpecificGene(Gene specificGene) {
		this.specificGene = specificGene;
	}


	public int getNrOfHits() {
		return nrOfHits;
	}


	public void setNrOfHits(int nrOfHits) {
		this.nrOfHits = nrOfHits;
	}




	//	
	//	public BlastHit(int start, int length, int missmatches, boolean plusStrand){
	//		
	//		this.start = start;
	//		this.length = length;
	//		this.plusStrand = plusStrand;
	//		this.missmatches = missmatches;
	//	}
	//	
	//	public BlastHit(int start, int length, int missmatches, boolean plusStrand, String chromosome){
	//		
	//		this.start = start;
	//		this.length = length;
	//		this.plusStrand = plusStrand;
	//		this.missmatches = missmatches;
	//		this.chromosome = chromosome;
	//	}
	//	
	//	public BlastHit(int start, int length, boolean plusStrand, String chromosome,String readSeq, int score){
	//		this.start = start;
	//		this.length = length;
	//		this.plusStrand = plusStrand;
	//		this.chromosome = chromosome;
	//		this.readSeaq = readSeq;
	//		this.score = score;
	//	
	//	}
	//	
	//	public BlastHit(int start, int length, boolean plusStrand, String chromosome,String readSeq, int score, Solid SolidSequence){
	//		this.start = start;
	//		this.length = length;
	//		this.plusStrand = plusStrand;
	//		this.chromosome = chromosome;
	//		this.readSeaq = readSeq;
	//		this.score = score;
	//		this.solidSequence = SolidSequence;
	//	}
	//
	//	
	//	public BlastHit(String chromosome, boolean plusStrand,int start, int contigend,int readStart, int readEnd,
	//			int length, int score,String readSeq){
	//
	//		this.chromosome = chromosome;
	//		this.plusStrand = plusStrand;
	//		this.start = start;
	//		this.contigend = contigend;
	//		this.readstart = readStart;
	//		this.readend = readEnd;
	//		this.length = length;
	//		this.score = score;
	//		this.readSeaq = readSeq;
	//	}
	//
	//
	//	public void printRmapper(ExtendedWriter EW){
	//		//#FORMAT: readname contigname strand contigstart contigend readstart readend readlength score editstring readsequence
	//		//>856_759_1307_F3        adaptor + bc1 + P2      +       1       25      15      39      50      208     x18x2x5 A22301003111312033020103031311231101032002322003120
	//		String strand = "+";
	//		if(!plusStrand) strand = "-";
	//			
	//			EW.print("\t"+this.chromosome+
	//					"\t"+strand+
	//					"\t"+this.start+
	//					"\t"+this.contigend+
	//					"\t"+this.readstart+
	//					"\t"+this.readend+
	//					"\t"+this.length+
	//					"\t"+this.score+
	//					"\t"+this.readSeaq
	//					);
	//	}
	//	
	//	public void printHit(ExtendedWriter EW){
	//		//#FORMAT: readname contigname strand contigstart contigend readstart readend readlength score editstring readsequence
	//		//>856_759_1307_F3        adaptor + bc1 + P2      +       1       25      15      39      50      208     x18x2x5 A22301003111312033020103031311231101032002322003120
	//		String strand = "+";
	//		if(!plusStrand) strand = "-";
	//			EW.print(this.solidSequence.Name+
	//					"\t"+this.chromosome+
	//					"\t"+strand+
	//					"\t"+this.start+
	//					"\t"+this.contigend+
	//					"\t"+this.readstart+
	//					"\t"+this.readend+
	//					"\t"+this.length+
	//					"\t"+this.score+
	//					"\t"+this.readSeaq
	//			);
	//			if(this.solidSequence.Sequence != null)
	//				EW.println("\t"+this.solidSequence.Sequence);
	//			else
	//				EW.println();
	//	}
	//	
	//	
	//
	//	
	//	
	//	
	//	public void printHitOLD(ExtendedWriter EW){
	//		if(this.plusStrand)
	//			EW.println(
	//				this.solidSequence.Name+"\t"+
	//				"+\t"+
	//				this.start+"\t"+
	//				this.length+"\t"+
	//				this.score+"\t"+
	//				this.readSeaq);		
	//		else
	//			EW.println(
	//					this.solidSequence.Name+"\t"+
	//					"-\t"+
	//					this.start+"\t"+
	//					this.length+"\t"+
	//					this.score+"\t"+
	//					this.readSeaq);		
	//	}
	//	
	//	
	//	public void printSurrounding(Chromosome C, int cutoff, int width, String GeneName, ExtendedWriter EW){
	//		if(this.nrOfReads > cutoff){
	//			C.printSurroundingSequence(start, length, plusStrand, width, GeneName ,nrOfReads ,  EW);
	//		}
	//	}
	//	
	//	
	//	public double getWeightedHit(){
	//		return 1/(double)nrOfHits;
	//		
	//	}
	//	
	//
	//	
	//
	//	public boolean sameLocation(String chromosome, int start){
	//		if(this.chromosome.compareTo(chromosome ) == 0 && this.start == start) return true;
	//		return false;
	//	}
	//	
	//	public boolean sameLocation(boolean plusStrand, int start){
	//		if(this.plusStrand == plusStrand && this.start == start) return true;
	//		return false;
	//	}
	//	
	//	
	//	public boolean isAntisense(int surrounding, int otherLocation, boolean OtherPlusStrand){
	//		if(this.plusStrand != OtherPlusStrand){
	//			int location = this.getLocation();
	//			if(otherLocation > location-surrounding && otherLocation < location+surrounding)
	//				return true;
	//		}
	//		return false;
	//	}
	//	
	//	public int getLocation(){
	//		if(this.plusStrand)
	//			return this.start+this.length/2;
	//		return this.start-this.length/2;
	//	}
	//
	//	
	//	
	//	public void mapHit(FastaSequences FS, double  weight){
	//		
	//		FS.get(Integer.parseInt(chromosome)-1).nrOfHits[0] += weight;
	//	}
	//	
	//	public void mapHit(FastaSequences FS, double  weight, int exp){
	//		
	//		FS.get(Integer.parseInt(chromosome)-1).nrOfHits[exp] += weight;
	//	}
	//	
	//	
	//	public int getKind(){
	//		if(specificGene != null){
	//			int kind = specificGene.kind;
	//			if(kind ==  intergenic)
	//				return kind;
	//			if(kind == repeat){
	//				if(specificGene.Name.indexOf("DIRS-1") > -1){
	//					return 6;
	//				}
	//				else{
	//					return kind;
	//				}
	//			}	
	//			if(kind == ncRNA)
	//				return ncRNA;
	//			if(kind == mRNA){
	//				if(specificGene.getName().indexOf("_RTE") >-1 || specificGene.getName().indexOf("_TE") > -1)
	//						return repeat;
	//			}
	//		
	//			if(this.plusStrand != specificGene.plusStrand)
	//				return antisense;
	//			else return specificGene.kind;
	//		}
	//		return -1;
	//	}
	//	
	//	public void print(ExtendedWriter EW){
	//		if(chromosome != null){
	//			EW.print(this.chromosome+"_");
	//		}
	//		EW.print(this.start+"."+this.missmatches);
	//		if(this.length > 0){
	//			EW.print("."+this.length);
	//		}
	//			
	//	}
	//
	//	
	//	public String getFastaInfo(){
	//		String info = "";
	//		if(chromosome != null){
	//			info = this.chromosome+"_";
	//		}
	//		info += this.start+"."+this.missmatches;
	//		if(this.length > 0){
	//			info += "."+this.length;
	//		}
	//		return info;
	//	}

}
