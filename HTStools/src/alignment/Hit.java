package alignment;

import Sequence.FastaSequences;
import general.ExtendedWriter;

public class Hit {
	protected int start;
	protected int length;
	protected int missmatches;
	protected boolean plusStrand;
	protected String chromosome;
	protected int kind;
	protected String readSeaq;
	protected Gene specificGene;
	public int score;
	protected int contigend;// = Integer.parseInt(info[4]);      //		25
	protected int readstart;//  = Integer.parseInt(info[5]);     // 		24
	protected int readend;//  = Integer.parseInt(info[6]);     // 		24

	
	public int nrOfHits;
	public int nrOfReads = 1;

	private final int mRNA = 1;
	private final int ncRNA = 2;
	private final int intergenic = 3;
	private final int antisense = 4;
	private final int repeat = 5;
	

	public boolean hasSameLocation(Hit otherHit){
		if(this.start == otherHit.start &&
				this.length == otherHit.length&&
				this.plusStrand == otherHit.plusStrand 
		)
			return true;
		return false;
	}
	
	public Hit(int score, Gene query, String HitName){
			this.score = score;
			this.specificGene = query;
			this.chromosome = HitName;
	}
	

	
	public Hit(int start, int length, int missmatches, boolean plusStrand){
		
		this.start = start;
		this.length = length;
		this.plusStrand = plusStrand;
		this.missmatches = missmatches;
	}
	
	public Hit(int start, int length, int missmatches, boolean plusStrand, String chromosome){
		
		this.start = start;
		this.length = length;
		this.plusStrand = plusStrand;
		this.missmatches = missmatches;
		this.chromosome = chromosome;
	}
	
	public Hit(int start, int length, boolean plusStrand, String chromosome,String readSeq, int score){
		this.start = start;
		this.length = length;
		this.plusStrand = plusStrand;
		this.chromosome = chromosome;
		this.readSeaq = readSeq;
		this.score = score;
	
	}
	

	
	public Hit(String chromosome, boolean plusStrand,int start, int contigend,int readStart, int readEnd,
			int length, int score,String readSeq){

		this.chromosome = chromosome;
		this.plusStrand = plusStrand;
		this.start = start;
		this.contigend = contigend;
		this.readstart = readStart;
		this.readend = readEnd;
		this.length = length;
		this.score = score;
		this.readSeaq = readSeq;
	}


	public void printRmapper(ExtendedWriter EW){
		//#FORMAT: readname contigname strand contigstart contigend readstart readend readlength score editstring readsequence
		//>856_759_1307_F3        adaptor + bc1 + P2      +       1       25      15      39      50      208     x18x2x5 A22301003111312033020103031311231101032002322003120
		String strand = "+";
		if(!plusStrand) strand = "-";
			
			EW.print("\t"+this.chromosome+
					"\t"+strand+
					"\t"+this.start+
					"\t"+this.contigend+
					"\t"+this.readstart+
					"\t"+this.readend+
					"\t"+this.length+
					"\t"+this.score+
					"\t"+this.readSeaq
					);
	}
	
	
	

	
	
	public void printSurrounding(Chromosome C, int cutoff, int width, String GeneName, ExtendedWriter EW){
		if(this.nrOfReads > cutoff){
			C.printSurroundingSequence(start, length, plusStrand, width, GeneName ,nrOfReads ,  EW);
		}
	}
	
	
	public double getWeightedHit(){
		return 1/(double)nrOfHits;
		
	}
	

	

	public boolean sameLocation(String chromosome, int start){
		if(this.chromosome.compareTo(chromosome ) == 0 && this.start == start) return true;
		return false;
	}
	
	public boolean sameLocation(boolean plusStrand, int start){
		if(this.plusStrand == plusStrand && this.start == start) return true;
		return false;
	}
	
	
	public boolean isAntisense(int surrounding, int otherLocation, boolean OtherPlusStrand){
		if(this.plusStrand != OtherPlusStrand){
			int location = this.getLocation();
			if(otherLocation > location-surrounding && otherLocation < location+surrounding)
				return true;
		}
		return false;
	}
	
	public int getLocation(){
		if(this.plusStrand)
			return this.start+this.length/2;
		return this.start-this.length/2;
	}

	
	
	public void mapHit(FastaSequences FS, double  weight){
		
		FS.get(Integer.parseInt(chromosome)-1).nrOfHits[0] += weight;
	}
	
	public void mapHit(FastaSequences FS, double  weight, int exp){
		
		FS.get(Integer.parseInt(chromosome)-1).nrOfHits[exp] += weight;
	}
	
	
	public int getKind(){
		if(specificGene != null){
			int kind = specificGene.kind;
			if(kind ==  intergenic)
				return kind;
			if(kind == repeat){
				if(specificGene.Name.indexOf("DIRS-1") > -1){
					return 6;
				}
				else{
					return kind;
				}
			}	
			if(kind == ncRNA)
				return ncRNA;
			if(kind == mRNA){
				if(specificGene.getName().indexOf("_RTE") >-1 || specificGene.getName().indexOf("_TE") > -1)
						return repeat;
			}
		
			if(this.plusStrand != specificGene.plusStrand)
				return antisense;
			else return specificGene.kind;
		}
		return -1;
	}
	
	public void print(ExtendedWriter EW){
		if(chromosome != null){
			EW.print(this.chromosome+"_");
		}
		EW.print(this.start+"."+this.missmatches);
		if(this.length > 0){
			EW.print("."+this.length);
		}
			
	}

	
	public String getFastaInfo(){
		String info = "";
		if(chromosome != null){
			info = this.chromosome+"_";
		}
		info += this.start+"."+this.missmatches;
		if(this.length > 0){
			info += "."+this.length;
		}
		return info;
	}
	
}
