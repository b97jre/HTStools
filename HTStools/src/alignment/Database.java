package alignment;

import general.Functions;
import general.RNAfunctions;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;

import Sequence.Solid;

import general.ExtendedReader;
import general.ExtendedWriter;




public class Database implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//private String GenomeName;
	private String Name;
	private String ID;

	private int length;
	private int[] sequence;
	Hashtable <Integer,StructuralVariation> SVs;
	ArrayList <Exon> filters;
	protected ArrayList <Gene> codingGenes;
	protected ArrayList <Gene> ncRNAs;
	protected ArrayList <Gene> repeats;
	protected ArrayList <Gene> intergenicRegions;

	private double coverage;
	private double coverage2;




	public void addGFF3LineInfo(String[] columns){

		if(columns[2].indexOf("gene") == 0)addGene(columns);
		else if(columns[2].indexOf("mRNA") == 0)addmRNA(columns);
		else if(columns[2].indexOf("exon") == 0)addExon(columns);
		else if(columns[2].indexOf("chromosome") == 0)getChromosomeInfo(columns);
		else if(columns[2].indexOf("dr") ==0)addNCRNA(columns);
		else if(columns[2].indexOf("Dictyostelium discoideum complex repeat") > 0)addRepeat(columns);

	}



	protected ArrayList <Hit> hits;



	public static void main(String[] args) {
		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}
		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);
		String Chromosome = Functions.getValue(T, "-c", "1");
		String gffDir = Functions.getValue(T, "-gffDir", ".");
		String gffFile = Functions.getValue(T, "-gffFile", "chromosome_"+Chromosome+".gff");
		Chromosome test = new Chromosome(gffDir,gffFile,Chromosome);
		int start = Integer.parseInt(Functions.getValue(T, "-start", "300"));
		int length2 = Integer.parseInt(Functions.getValue(T, "-length", "21"));
		int width = Integer.parseInt(Functions.getValue(T, "-width", "200"));
		boolean plusStrand = Boolean.parseBoolean(Functions.getValue(T, "-plusStrand", "true"));
		test.printHitSequence(start, length2, plusStrand, width,"genome");

	}	





	Database(String Name, ExtendedReader ER){
		this.setName(Name);
		readChromosomeSequence(ER);
	}

	Database(String Name){
		this.setName(Name);
		this.codingGenes = new ArrayList<Gene>();
		this.ncRNAs = new ArrayList<Gene>();
		this.repeats = new ArrayList<Gene>();

	}

	public void initiateAnnotation(){
		this.codingGenes = new ArrayList<Gene>();
		this.ncRNAs = new ArrayList<Gene>();
		this.repeats = new ArrayList<Gene>();

	}



	public boolean isDatabase(String name){
		if(this.Name.indexOf(name) ==0) return true;
		return false;
	}




	public void addCoverage(int start, int stop, String info){
		this.coverage = (double)(stop-start) / (double)(this.length);

		int ind1 = info.indexOf("cov \"") +5;
		int ind2 = info.lastIndexOf("\"");
		String temp = info.substring(ind1,ind2);
		this.coverage2 = Double.parseDouble(temp);
	}

	public void printCoverage(){
		System.out.println(this.Name +"\t"+ this.coverage+"\t"+ this.coverage2);
	}


	public void printDistribution(ExtendedWriter EW, int cutoff){
		int[] plusStrand = new int[this.length];
		int[] minusStrand = new int[this.length];
		//		EW.println(this.Name);
		//		EW.println("cutoff:"+ cutoff);
		try{
			if(this.hits != null){

				for(int i = 0; i< this.hits.size();i++ ){
					int start = this.hits.get(i).start;
					int stop = this.hits.get(i).contigend;
					boolean strand = this.hits.get(i).plusStrand;
					if(strand){
						for(int j = start; j <= stop; j++){
							plusStrand[j-1]++;	
						}
					}
					else{
						for(int j = start; j <= stop; j++){
							minusStrand[j-1]++;
						}
					}
				}
			}
		}
		catch(Exception E){E.printStackTrace();}
		//		EW.println(this.Name);
		//		EW.println("cutoff:"+ cutoff);
		//EW.println(" \t"+experiment+"\t"+experiment);
		int total = 0;
		EW.println("Loc\tPlusstrand\tMinusStrand\tDifference");
		for(int i =0 ; i < plusStrand.length;i++){
			total +=(plusStrand[i] + minusStrand[i]);
			if(plusStrand[i] > cutoff || minusStrand[i] >cutoff )
				EW.println((i+1)+"\t"+plusStrand[i]+"\t-"+minusStrand[i]+"\t"+(plusStrand[i]-minusStrand[i]));
			else if(i < plusStrand.length-1 && (plusStrand[i+1] > cutoff || minusStrand[i+1] >cutoff ) )
				EW.println((i+1)+"\t0\t0\t0");
			else if(i > 0 && (plusStrand[i-1] > cutoff || minusStrand[i-1] >cutoff ) )
				EW.println((i+1)+"\t0\t0\t0");
		}
		EW.println();
		//		System.out.print(total/plusStrand.length+"\t");
	}

	public void printOverallDistribution(ExtendedWriter EW, int cutoff){
		int[] plusStrand = new int[this.length];
		int[] minusStrand = new int[this.length];
		//		EW.println(this.Name);
		//		EW.println("cutoff:"+ cutoff);
		try{
			if(this.hits != null){
				for(int i = 0; i< this.hits.size();i++ ){
					int start = this.hits.get(i).start;
					int stop = this.hits.get(i).contigend;
					boolean strand = this.hits.get(i).plusStrand;
					if(strand){
						for(int j = start; j <= stop; j++){
							plusStrand[j-1]++;	
						}
					}
					else{
						for(int j = start; j <= stop; j++){
							minusStrand[j-1]++;
						}
					}
				}
			}
		}
		catch(Exception E){E.printStackTrace();}
		int plusStrandMax = 0;
		int negStrandMax = 0;
		int diffMax;
		double plusStrandMean;
		double plusStrandStd;
		double negStrandMean;
		double negStrandStd;


		int plusTotal = 0;
		int negTotal = 0;
		for(int i =0 ; i < plusStrand.length;i++){
			plusTotal +=plusStrand[i];
			negTotal +=minusStrand[i];
			if(plusStrand[i] > plusStrandMax) plusStrandMax = plusStrand[i];
			if(minusStrand[i] > negStrandMax) negStrandMax  = minusStrand[i];
		}
		plusStrandMean = Functions.getMean(plusStrand);
		negStrandMean = Functions.getMean(minusStrand);
		plusStrandStd = Functions.getSD(plusStrand, plusStrandMean);
		negStrandStd = Functions.getSD(minusStrand, negStrandMean);

		EW.println(this.Name+"\t"+plusStrandMax+"\t"+negStrandMax+"\t"+plusStrandMean+"\t"+plusStrandStd+"\t"+negStrandMean+"\t"+negStrandStd);
	}

	public void printMaxSequence(ExtendedWriter EW, int length, int cutoff){
		int[] plusStrand = new int[this.length];
		int[] minusStrand = new int[this.length];
		//		EW.println(this.Name);
		//		EW.println("cutoff:"+ cutoff);
		try{
			if(this.hits != null){
				for(int i = 0; i< this.hits.size();i++ ){
					int start = this.hits.get(i).start;
					int stop = this.hits.get(i).contigend;
					boolean strand = this.hits.get(i).plusStrand;
					System.out.println(Math.abs(start-stop)+1);
					if( Math.abs(start-stop)== length-1){
						if(strand){
							plusStrand[start-1]++;	
						}
						else{
							minusStrand[start-1]++;
						}
					}
				}
			}
		}
		catch(Exception E){E.printStackTrace();}
		int plusStrandMax = 0;
		int negStrandMax = 0;
		int plusStrandMaxLocation = 0;
		int negStrandMaxLocation = 0;


		for(int i =0 ; i < plusStrand.length;i++){
			if(plusStrand[i] > plusStrandMax){
				plusStrandMax = plusStrand[i];
				plusStrandMaxLocation = i;
			}
			if(minusStrand[i] > negStrandMax){
				negStrandMax  = minusStrand[i];
				negStrandMaxLocation = i;
			}
		}
		int[] plus = new int[length];
		int[] neg = new int[length];

		if(plusStrandMax > cutoff){
			for(int i = 0 ; i < length; i++ ){
				plus[i] = this.sequence[i+plusStrandMaxLocation];
			}
		}
		if(negStrandMax > cutoff){
			for(int i = 0 ; i < length; i++ ){
				neg[i] = this.sequence[negStrandMaxLocation+i];
			}
			neg = RNAfunctions.getReverseComplement(neg);
		}




		if(plusStrandMax > cutoff){
			EW.println(">"+this.Name+"_plustStrand_[nrOfreads:"+plusStrandMax+"]_(length:"+length+")_{start:"+(plusStrandMaxLocation+1)+"}");
			EW.println(RNAfunctions.DNAInt2String(plus));
		}
		if(negStrandMax > cutoff){
			EW.println(">"+this.Name+"_negStrand_[nrOfreads:"+negStrandMax+"]_(length:"+length+")_{start:"+(negStrandMaxLocation+length)+"}");
			EW.println(RNAfunctions.DNAInt2String(neg));
		}

	}



	private void readChromosomeSequence(ExtendedReader ER){
		int pointer = 0;
		int[] sequence = new int[20000000];
		while(ER.more() && ER.lookAhead() != '>'){
			int[] subSequence = RNAfunctions.RNAString2Int(ER.readLine());
			for(int i = 0; i < subSequence.length;i++){
				sequence[pointer] = subSequence[i];
				pointer++;
			}
		}


		int[] sequence2 = new int[pointer];
		this.length = pointer;
		for(int i = 0; i < pointer; i++){
			sequence2[i] = sequence[i];
		}

		this.sequence = sequence2;
	}



	public void getChromosomeSequenceSize(ExtendedReader ER){
		this.length = 0;
		while(ER.more() && ER.lookAhead() != '>'){
			this.length += ER.readLine().length();
		}
	}


	public void mapSolidSequence2Database(Solid hit){
		for(int i = 0; i < hit.hits.size(); i++){
			if(hit.hits.get(i).chromosome.compareTo(this.Name) == 0){
				this.addHit(hit.hits.get(i));
			}
		}
	}


	public void addHit(Hit newHit){
		if(hits == null)
			hits = new ArrayList<Hit>();
		hits.add(newHit);
	}

	public int getNrOfHits(){
		return this.hits.size();
	}

	public void addVCFinfo(ArrayList<String> Samples, String[] VCFinfo){
		//		0			1		2     3		  4			5		6		7		8		9				10			11			12				13				14				15				16		 
		// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
		// scaffold_1      5       .       A       	      56.41   .       AC=3;AF=0.188;AN=16;BaseQRankSum=-1.804;DP=72;Dels=0.00;FS=0.000;HaplotypeScore=2.4206;MLEAC=2;MLEAF=0.125;MQ=26.05;MQ0=11;MQRankSum=0.618;QD=2.69;ReadPosRankSum=1.053 GT:AD:DP:GQ:PL  0|0:14,0:14:30:0,30,268 0|0:1,0:1:3:0,3,29      0/1:13,5:18:72:72,0,247 0|0:14,0:14:39:0,39,316 0|0:13,0:13:39:0,39,328 0|0:4,0:4:3:0,3,25      1|1:2,1:3:3:23,3,0      0|0:5,0:5:3:0,3,25
		if(this.SVs == null){
			this.SVs = new Hashtable<Integer,StructuralVariation>();
		}
		Integer Location = Integer.decode(VCFinfo[1]);
		if(VCFinfo[5].compareTo(".")!=0)
			this.SVs.put(Location, new StructuralVariation(Samples,VCFinfo));
	}

	public void addVCFphase(String sample, String[] VCFinfo, ExtendedWriter EW,int[] info){
		//		0			1		2     3		  4			5		6		7		8		9				10			11			12				13				14				15				16		 
		// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
		// scaffold_1      5       .       A       	      56.41   .       AC=3;AF=0.188;AN=16;BaseQRankSum=-1.804;DP=72;Dels=0.00;FS=0.000;HaplotypeScore=2.4206;MLEAC=2;MLEAF=0.125;MQ=26.05;MQ0=11;MQRankSum=0.618;QD=2.69;ReadPosRankSum=1.053 GT:AD:DP:GQ:PL  0|0:14,0:14:30:0,30,268 0|0:1,0:1:3:0,3,29      0/1:13,5:18:72:72,0,247 0|0:14,0:14:39:0,39,316 0|0:13,0:13:39:0,39,328 0|0:4,0:4:3:0,3,25      1|1:2,1:3:3:23,3,0      0|0:5,0:5:3:0,3,25

		Integer Location = Integer.decode(VCFinfo[1]);
		if(this.SVs.containsKey(Location)){
			this.SVs.get(Location).addPhase(EW,VCFinfo,sample,info);
			EW.print(this.Name+"\t"+Location+"\t");
			SVs.get(Location).printSample(EW, sample);
			EW.println();
		}

	}

	public void markSubset(String[] VCFinfo){
		//		0			1		2     3		  4			5		6		7		8		9				10			11			12				13				14				15				16		 
		// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
		// scaffold_1      5       .       A       	      56.41   .       AC=3;AF=0.188;AN=16;BaseQRankSum=-1.804;DP=72;Dels=0.00;FS=0.000;HaplotypeScore=2.4206;MLEAC=2;MLEAF=0.125;MQ=26.05;MQ0=11;MQRankSum=0.618;QD=2.69;ReadPosRankSum=1.053 GT:AD:DP:GQ:PL  0|0:14,0:14:30:0,30,268 0|0:1,0:1:3:0,3,29      0/1:13,5:18:72:72,0,247 0|0:14,0:14:39:0,39,316 0|0:13,0:13:39:0,39,328 0|0:4,0:4:3:0,3,25      1|1:2,1:3:3:23,3,0      0|0:5,0:5:3:0,3,25

		Integer Location = Integer.decode(VCFinfo[1]);
		if(this.SVs.containsKey(Location)){
			this.SVs.get(Location).SkellyPrint = true;
		}

	}


	public void comparePhasedVCFinfo(Database otherDatabase, String sample, ExtendedWriter EW){
		//		0			1		2     3		  4			5		6		7		8		9				10			11			12				13				14				15				16		 
		// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
		// scaffold_1      5       .       A       	      56.41   .       AC=3;AF=0.188;AN=16;BaseQRankSum=-1.804;DP=72;Dels=0.00;FS=0.000;HaplotypeScore=2.4206;MLEAC=2;MLEAF=0.125;MQ=26.05;MQ0=11;MQRankSum=0.618;QD=2.69;ReadPosRankSum=1.053 GT:AD:DP:GQ:PL  0|0:14,0:14:30:0,30,268 0|0:1,0:1:3:0,3,29      0/1:13,5:18:72:72,0,247 0|0:14,0:14:39:0,39,316 0|0:13,0:13:39:0,39,328 0|0:4,0:4:3:0,3,25      1|1:2,1:3:3:23,3,0      0|0:5,0:5:3:0,3,25


		if(this.SVs!=null && otherDatabase.SVs != null){
			List<Integer> sortedKeys=new ArrayList<Integer>(this.SVs.keySet());
			Collections.sort(sortedKeys);
			List<Integer> OtherSortedKeys=new ArrayList<Integer>(otherDatabase.SVs.keySet());
			Collections.sort(OtherSortedKeys);

			int currentPhase = 0;

			int nrOfSNPs = 0;
			int nrOfIdentical = 0;
			int start = sortedKeys.get(0);
			for(int i = 0; i < sortedKeys.size();i++){
				if(otherDatabase.SVs.containsKey(sortedKeys.get(i))){
					StructuralVariationSample SV1 = SVs.get(sortedKeys.get(i)).getSample(sample);
					StructuralVariationSample SV2 = otherDatabase.SVs.get(sortedKeys.get(i)).getSample(sample);

					if(!SV1.isHomozygous()){
						if(SV1.phaseNr != currentPhase){
							int stop = sortedKeys.get(i-1);

							if(nrOfSNPs != 1){
								if(nrOfIdentical*2 < nrOfSNPs){
									nrOfIdentical = nrOfSNPs - nrOfIdentical;
								}
								EW.println(this.Name+"\t"+currentPhase+"\t"+start+"\t"+stop+"\t"+nrOfSNPs+"\t"+nrOfIdentical);
							}
							start = sortedKeys.get(i);
							currentPhase=SV1.phaseNr;
							nrOfSNPs = 1; 
							nrOfIdentical = 0;
							if(SV1.isSame(SV2)){
								nrOfIdentical++;
							}
						}
					}else{
						nrOfSNPs++;
						if(SV1.isSame(SV2)){
							nrOfIdentical++;
						}
					}
				}
			}
		}
	}


	public void trimPhasedVCFinfo(Database otherDatabase, String sample){
		//		0			1		2     3		  4			5		6		7		8		9				10			11			12				13				14				15				16		 
		// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
		// scaffold_1      5       .       A       	      56.41   .       AC=3;AF=0.188;AN=16;BaseQRankSum=-1.804;DP=72;Dels=0.00;FS=0.000;HaplotypeScore=2.4206;MLEAC=2;MLEAF=0.125;MQ=26.05;MQ0=11;MQRankSum=0.618;QD=2.69;ReadPosRankSum=1.053 GT:AD:DP:GQ:PL  0|0:14,0:14:30:0,30,268 0|0:1,0:1:3:0,3,29      0/1:13,5:18:72:72,0,247 0|0:14,0:14:39:0,39,316 0|0:13,0:13:39:0,39,328 0|0:4,0:4:3:0,3,25      1|1:2,1:3:3:23,3,0      0|0:5,0:5:3:0,3,25

		if(this.SVs!=null && otherDatabase.SVs != null){
			List<Integer> sortedKeys=new ArrayList<Integer>(this.SVs.keySet());
			Collections.sort(sortedKeys);
			int currentPhase = SVs.get(sortedKeys.get(0)).getPhase(sample);
			for(int i = 0; i < sortedKeys.size();i++){
				if(otherDatabase.SVs.containsKey(sortedKeys.get(i))){
					if(SVs.get(sortedKeys.get(i)).getPhase(sample) != currentPhase){
						if(SVs.get(sortedKeys.get(i)).isPhased(sample)){
							System.out.println("Unphasing "+this.Name +"\t"+sortedKeys.get(i));
						}
						SVs.get(sortedKeys.get(i)).unPhase(sample);
						currentPhase = SVs.get(sortedKeys.get(i)).getPhase(sample);
					}
				}else{
					System.out.println("removing "+this.Name +"\t"+sortedKeys.get(i));
					SVs.remove(sortedKeys.get(i));
				}
			}
		}else if(this.SVs!=null){
			this.SVs =null;
		}

	}


	public void compareVCFinfo(Database otherDatabase,ArrayList<String> samples, ExtendedWriter EW1,ExtendedWriter EW2){
		//		0			1		2     3		  4			5		6		7		8		9				10			11			12				13				14				15				16		 
		// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
		// scaffold_1      5       .       A       	      56.41   .       AC=3;AF=0.188;AN=16;BaseQRankSum=-1.804;DP=72;Dels=0.00;FS=0.000;HaplotypeScore=2.4206;MLEAC=2;MLEAF=0.125;MQ=26.05;MQ0=11;MQRankSum=0.618;QD=2.69;ReadPosRankSum=1.053 GT:AD:DP:GQ:PL  0|0:14,0:14:30:0,30,268 0|0:1,0:1:3:0,3,29      0/1:13,5:18:72:72,0,247 0|0:14,0:14:39:0,39,316 0|0:13,0:13:39:0,39,328 0|0:4,0:4:3:0,3,25      1|1:2,1:3:3:23,3,0      0|0:5,0:5:3:0,3,25


		if(this.SVs!=null && otherDatabase.SVs != null){
			List<Integer> sortedKeys=new ArrayList<Integer>(this.SVs.keySet());
			Collections.sort(sortedKeys);
			List<Integer> OtherSortedKeys=new ArrayList<Integer>(otherDatabase.SVs.keySet());
			Collections.sort(OtherSortedKeys);

			for(int i = 0; i < sortedKeys.size();i++){
				if(!otherDatabase.SVs.containsKey(sortedKeys.get(i))){
					EW1.print(this.Name+"\t"+sortedKeys.get(i)+"\t");
					SVs.get(sortedKeys.get(i)).printSamples(EW1,samples);
					EW1.println();
				}
			}
			for(int i = 0; i < OtherSortedKeys.size();i++){
				if(!this.SVs.containsKey(OtherSortedKeys.get(i))){
					EW2.print(otherDatabase.Name+"\t"+OtherSortedKeys.get(i)+"\t");
					otherDatabase.SVs.get(OtherSortedKeys.get(i)).printSamples(EW2,samples);
					EW2.println();
				}
			}

		}
	}



	public String[] addVCFSamples(ExtendedReader ER, ArrayList<String> newSamples,String [] VCFinfo){

		if(SVs == null){
			while(VCFinfo[0].compareTo(this.Name) == 0)
				VCFinfo = ER.readLine().split("\t");
			System.out.println("Finished");
			System.out.println(this.Name +"\t"+VCFinfo[0]);
			return VCFinfo;
		}
		//		0			1		2     3		  4			5		6		7		8		9				10			11			12				13				14				15				16		 
		// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
		// scaffold_1      5       .       A       T       56.41   .       AC=3;AF=0.188;AN=16;BaseQRankSum=-1.804;DP=72;Dels=0.00;FS=0.000;HaplotypeScore=2.4206;MLEAC=2;MLEAF=0.125;MQ=26.05;MQ0=11;MQRankSum=0.618;QD=2.69;ReadPosRankSum=1.053 GT:AD:DP:GQ:PL  0|0:14,0:14:30:0,30,268 0|0:1,0:1:3:0,3,29      0/1:13,5:18:72:72,0,247 0|0:14,0:14:39:0,39,316 0|0:13,0:13:39:0,39,328 0|0:4,0:4:3:0,3,25      1|1:2,1:3:3:23,3,0      0|0:5,0:5:3:0,3,25
		List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
		Collections.sort(sortedKeys);
		int pointer = 0;
		Integer Location = null;

		//First run
		if(VCFinfo[0].compareTo(this.Name) != 0){
			System.out.println("Finished");
			//System.out.println(this.Name +"\t"+VCFinfo[0]);
			return VCFinfo;
		}
		Location = Integer.decode(VCFinfo[1]);
		while(pointer < sortedKeys.size() && Location > sortedKeys.get(pointer)){
			pointer++;
		}
		if(pointer < sortedKeys.size() && Location.compareTo(sortedKeys.get(pointer)) == 0){
			//System.out.println("adding sample: "+Location+"\t"+sortedKeys.get(pointer));
			this.SVs.get(Location).addSamples(newSamples,VCFinfo);
			VCFinfo = ER.readLine().split("\t");
		}
		try{
			while(ER.more()){

				//				System.out.println(this.Name +"\t"+VCFinfo[0]);

				if(VCFinfo[0].compareTo(this.Name) != 0){
					System.out.println("Finished");
					//System.out.println(this.Name +"\t"+VCFinfo[0]);
					return VCFinfo;
				}
				Location = Integer.decode(VCFinfo[1]);
				while(pointer < sortedKeys.size() && Location > sortedKeys.get(pointer)){
					pointer++;
				}
				if(pointer < sortedKeys.size() && Location.compareTo(sortedKeys.get(pointer)) == 0 ){
					//System.out.println("adding sample: "+Location+"\t"+sortedKeys.get(pointer));
					this.SVs.get(Location).addSamples(newSamples,VCFinfo);
				}
				//System.out.println("testing "+ Location);
				VCFinfo = ER.readLine().split("\t");
			}

		}
		catch(Exception E){
			//		System.out.println(Location+"\t"+sortedKeys.get(pointer));
			E.printStackTrace();

		}
		return null;
	}


	public void addPhase(boolean readPhased){
		if(SVs!=null){
			List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
			SVs.get(sortedKeys.get(0)).addPhase(null,readPhased);
			for(int i = 1; i < sortedKeys.size();i++){
				SVs.get(sortedKeys.get(i)).addPhase(SVs.get(sortedKeys.get(i-1)),readPhased);
			}
		}



	}


	//	public String[] checkPhased(){
	//		
	//		List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
	//		Collections.sort(sortedKeys);
	//		int pointer = 0;
	//		Integer Location = null;
	//		int phase = 0;
	//		for(int i =0 ; i < sortedKeys.size();i++){
	//			phase = SVs.get(i).getPhase(phase);	
	//		}
	//		if(VCFinfo[0].compareTo(this.Name) != 0){
	//			System.out.println("Finished");
	//			//System.out.println(this.Name +"\t"+VCFinfo[0]);
	//			return VCFinfo;
	//		}
	//		Location = Integer.decode(VCFinfo[1]);
	//		while(pointer < sortedKeys.size() && Location > sortedKeys.get(pointer)){
	//			pointer++;
	//		}
	//		if(pointer < sortedKeys.size() && Location.compareTo(sortedKeys.get(pointer)) == 0){
	//			//System.out.println("adding sample: "+Location+"\t"+sortedKeys.get(pointer));
	//			this.SVs.get(Location).addSamples(newSamples,VCFinfo);
	//			VCFinfo = ER.readLine().split("\t");
	//		}
	//		try{
	//			while(ER.more()){
	//
	////				System.out.println(this.Name +"\t"+VCFinfo[0]);
	//				
	//				if(VCFinfo[0].compareTo(this.Name) != 0){
	//					System.out.println("Finished");
	//					//System.out.println(this.Name +"\t"+VCFinfo[0]);
	//					return VCFinfo;
	//				}
	//				Location = Integer.decode(VCFinfo[1]);
	//				while(pointer < sortedKeys.size() && Location > sortedKeys.get(pointer)){
	//					pointer++;
	//				}
	//				if(pointer < sortedKeys.size() && Location.compareTo(sortedKeys.get(pointer)) == 0 ){
	//					//System.out.println("adding sample: "+Location+"\t"+sortedKeys.get(pointer));
	//					this.SVs.get(Location).addSamples(newSamples,VCFinfo);
	//				}
	//				//System.out.println("testing "+ Location);
	//				VCFinfo = ER.readLine().split("\t");
	//			}
	//			
	//		}
	//		catch(Exception E){
	//	//		System.out.println(Location+"\t"+sortedKeys.get(pointer));
	//			E.printStackTrace();
	//			
	//		}
	//		return null;
	//	}




	public void addBEDfilterInfo( String[] BEDinfo){
		//	0			1		2     3		  4			5		6		7		8		9			 
		// #CHROM(0)  Start     Stop     Name     .     .    INFO  type    Something  XTR 	AInfo   
		// scaffold_1      767     2124    PAC:20891551.exon.3     .       -       phytozome8_0    exon    .       ID=PAC:20891551.exon.3;Parent=PAC:20891551;pacid=20891551


		if(this.filters == null){
			this.filters = new ArrayList<Exon>();
		}
		Exon newFilter = new Exon(BEDinfo);
		this.filters.add(newFilter);
	}

	public void sortBEDfilters(){
		if(this.filters != null)
			Collections.sort(this.filters);
	}

	public void mergeBEDfilters(){
		if(this.filters == null) return;
		Collections.sort(this.filters);
		for(int i = 0; i< this.filters.size()-1;i++){
			if(filters.get(i).join(filters.get(i+1))){
				filters.remove(i+1);
				i--;
			}

		}

	}

	public void splitBEDfilters(){

		if(this.filters == null) return;
		ArrayList <Exon> NewFilters = new ArrayList <Exon>();
		boolean noMoreSplits = false;
		while(!noMoreSplits){
			noMoreSplits = true;
			Collections.sort(this.filters);
			boolean noSplits = true;
			int i = 0;
			while(i < this.filters.size()-1 && noSplits){
				if(filters.get(i).overlaps(filters.get(i+1))){
					ArrayList <Exon> tempFilters = filters.get(i).split(filters.get(i+1));
					filters.remove(i);
					filters.remove(i);
					for(int j = 0; j< tempFilters.size();j++ ){
						filters.add(tempFilters.get(j));
					}
					noSplits = false;
				}
				else{
					i++;	
				}
			}
			noMoreSplits = noSplits;
		}

	}

	public void countBEDfilters(int size, double fractionCutoff, ExtendedWriter EW , ExtendedWriter BedEW){
		if(this.filters == null)return;
		int bins = this.length/size+1; 
		int[]counts = new int[bins];  
		double[] fraction = new double[bins];  
		int pointer = 0;
		Collections.sort(this.filters);
		for(int j = 0; j < this.length; j++ ){
			if(this.filters.size() > pointer &&this.filters.get(pointer).left <=j){
				while(this.filters.get(pointer).right > j){
					counts[j/size]++;
					j++;
				}
				pointer++;
			}
		}

		for(int i = 0; i < counts.length;i++){
			fraction[i] = (double)counts[i]/(double)size;
			EW.println(this.Name+"\t"+i*size+"\t"+counts[i]+"\t"+fraction[i] );
		}


		int  count = 0;
		for(int i = 0; i < counts.length;i++){
			if(fraction[i] >fractionCutoff){ 
				BedEW.print(this.Name+"\t"+i*size+"\t");
				while(i < counts.length && fraction[i] >fractionCutoff){i++;}
				BedEW.println((i*size+size)+"\tFractionAbove"+count);
				count++;
			}
		}



	}

	public void printpreMRNAs(ExtendedWriter EW){
		if(this.codingGenes != null)
			for(int i = 0; i< this.codingGenes.size();i++){
				codingGenes.get(i).printPremRNA(EW,this.sequence);
			}
	}

	//	public void printCodingRNAs(ExtendedWriter EW){
	//		if(this.codingGenes != null)
	//			for(int i = 0; i< this.codingGenes.size();i++){
	//				codingGenes.get(i).printCodingRNA(EW,this.sequence,this.Name);
	//			}
	//	}

	public void printHaploContig(ExtendedWriter EW,String sample,boolean father,ExtendedWriter SNPinformation ){
		String newName = this.Name+"_"+sample+"_phased_";
		if(father) newName+="father.fa";
		else newName+= "mother.fa";
		List<Integer> sortedKeys = null;
		if(SVs != null){
			sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
		}
		EW.println(">"+newName);
		System.out.println(">"+newName);
		int newLine = 100;
		int pointer = 0;
		for(int i= 0; i < this.length;i++){
			if(i == newLine){
				newLine = newLine+100;
				EW.println();
			}
			if(sortedKeys != null && pointer < sortedKeys.size() && sortedKeys.get(pointer) == i+1){
				char[] SNP = this.SVs.get(sortedKeys.get(pointer)).getPhasedSNP(EW, sample,father);
				for(int j = 0; j < SNP.length; j++){
					EW.print(SNP[j]);
				}
				SNPinformation.print((i+1)+"\t"+RNAfunctions.DNAInt2char(this.sequence[i])+"\t");
				for(int j = 0; j < SNP.length; j++){
					SNPinformation.print(SNP[j]);
				}
				SNPinformation.println();
				pointer++;
			}
			else
				EW.print(RNAfunctions.DNAInt2char(this.sequence[i]));
		}
	}
	public void printVCFinfo(ExtendedWriter EW,String sample){
		int[] vcfInfo = new int[5];
		List<Integer> sortedKeys = null;
		if(SVs != null){
			sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
			for(int i = 0; i < sortedKeys.size();i++){
				vcfInfo = SVs.get(sortedKeys.get(i)).getSampleVCFinfo(vcfInfo,sample);
			}
		}
		
		EW.println(this.Name+"\t"+vcfInfo[0]+"\t"+vcfInfo[1]+"\t"+vcfInfo[2]+
				"\t"+vcfInfo[3]+"\t"+vcfInfo[4]+"\t"+sample);

	}



	public void printNNContig(ExtendedWriter EW,String sample){
		List<Integer> sortedKeys = null;
		if(SVs != null){
			sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
		}
		EW.println(">"+this.Name);
		int newLine = 100;
		int pointer = 0;
		for(int i= 0; i < this.length;i++){
			if(i == newLine){
				newLine = newLine+100;
				EW.println();
			}
			if(sortedKeys != null &&pointer < sortedKeys.size() && sortedKeys.get(pointer) == i+1){
				char[] SNP = this.SVs.get(sortedKeys.get(pointer)).getNNSNP(EW, sample);
				for(int j = 0; j < SNP.length; j++){
					EW.print(SNP[j]);
				}
			}else
				EW.print(RNAfunctions.DNAInt2char(this.sequence[i]));
		}

	}


	public void printmRNAs(ExtendedWriter EW){
		List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
		Collections.sort(sortedKeys);
		if(this.codingGenes != null)
			for(int i = 0; i< this.codingGenes.size();i++){
				codingGenes.get(i).printmRNA(EW,this.sequence,this.Name);
			}
	}

	public void printPeronsalmRNAs(ExtendedWriter[] EWs, String[] Samples){
		ArrayList<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
		Collections.sort(sortedKeys);
		if(this.codingGenes != null)
			for(int i = 0; i< this.codingGenes.size();i++){
				codingGenes.get(i).printPersonalmRNA(EWs,this.sequence,this.Name,this.SVs,sortedKeys,Samples);
			}
	}




	public void printPersonalmRNAsInfo(ArrayList<String> Samples, ExtendedWriter Info){
		if(this.SVs != null){
			ArrayList<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
			if(this.codingGenes != null)
				for(int i = 0; i< this.codingGenes.size();i++){
					codingGenes.get(i).printPersonalmRNAInfo(this.Name,this.SVs,sortedKeys,Samples, Info);

				}
		}
	}

	public void printPersonalmRNAsInfo(String Sample, ExtendedWriter Info){
		if(this.SVs != null){
			ArrayList<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
			if(this.codingGenes != null)
				for(int i = 0; i< this.codingGenes.size();i++){
					codingGenes.get(i).printPersonalmRNAInfo(this.Name,this.SVs,sortedKeys,Sample, Info);

				}
		}
	}



	public void printFilters(){
		if(this.filters != null)
			for(int i = 0; i< this.filters.size();i++){
				filters.get(i).printBED(this.Name);
			}
	}

	public void printFilters(ExtendedWriter EW){
		if(this.filters != null)
			for(int i = 0; i< this.filters.size();i++){
				filters.get(i).printBED(this.Name,EW);
			}
	}


	public void filterVCFinfoOutside(){
		if(SVs!=null){
			if(this.filters==null){
				this.SVs=null;
				return;
			}
			List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
			int pointer = 0;
			for(int i = 0; i < sortedKeys.size();i++)
			{
				if(this.filters.get(this.filters.size()-1).right<sortedKeys.get(i)){
					SVs.remove(sortedKeys.get(i));
				} 
				else if(sortedKeys.get(i) < this.filters.get(pointer).left){
					SVs.remove(sortedKeys.get(i));
					//System.out.println("removing "+sortedKeys.get(i));
				}
				else if(sortedKeys.get(i) > this.filters.get(pointer).right){
					i--;
					pointer++;
				}
			}
		}
	}

	public void filterVCFinfoInside(){
		if(SVs!=null){
			if(this.filters==null){
				this.SVs=null;
				return;
			}
			List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
			int pointer = 0;
			for(int i = 0; i < sortedKeys.size();i++)
			{
				if(this.filters.get(this.filters.size()-1).right<sortedKeys.get(i)){
					i=sortedKeys.size();
				} 
				else if(sortedKeys.get(i) <= this.filters.get(pointer).right && sortedKeys.get(i) >= this.filters.get(pointer).left){
					SVs.remove(sortedKeys.get(i));
					//System.out.println("removing "+sortedKeys.get(i));
				}
				else if(sortedKeys.get(i) > this.filters.get(pointer).right){
					i--;
					pointer++;
				}
			}
		}
	}

	public void annotateVCFinfo(ExtendedWriter EW){
		if(SVs!=null){
			if(this.filters==null){
				this.SVs=null;
				return;
			}
			List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
			int pointer = 0;
			for(int i = 0; i < sortedKeys.size();i++)
			{
				if(sortedKeys.get(i) > this.filters.get(this.filters.size()-1).right) 
					return;
				else if(sortedKeys.get(i) >= this.filters.get(pointer).left && sortedKeys.get(i) <= this.filters.get(pointer).right){
					EW.println(this.Name+"\t"+sortedKeys.get(i)+"\t"+this.filters.get(pointer).name);

				}
				else if(sortedKeys.get(i) > this.filters.get(pointer).right){
					i--;
					pointer++;
				}
			}
		}
	}


	public void SkellyFormat(ExtendedWriter EW,String sample){
		if(SVs!=null){
			if(this.filters==null){
				this.SVs=null;
				return;
			}
			List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
			int pointer = 0;
			for(int i = 0; i < sortedKeys.size();i++)
			{
				if(sortedKeys.get(i) > this.filters.get(this.filters.size()-1).right) 
					return;
				else if(sortedKeys.get(i) >= this.filters.get(pointer).left && sortedKeys.get(i) <= this.filters.get(pointer).right){
					if(SVs.get(sortedKeys.get(i)).isHeterozygous(sample) && SVs.get(sortedKeys.get(i)).isPhased(sample))
						EW.println(this.filters.get(pointer).name+"_"+SVs.get(sortedKeys.get(i)).getPhase(sample)+"\t"+this.Name+"_"+sortedKeys.get(i)
								+"\t"+SVs.get(sortedKeys.get(i)).getMotherCount(sample)
								+"\t"+SVs.get(sortedKeys.get(i)).getFatherCount(sample));
				}
				else if(sortedKeys.get(i) > this.filters.get(pointer).right){
					i--;
					pointer++;
				}
			}
		}
	}



	public void SplitFilterReadPhased(String sample){
		if(SVs!=null){
			if(this.filters==null){
				this.SVs=null;
				return;
			}
			List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
			int pointer = 0;
			int count = 0; 
			String filterName = null;
			ArrayList <Exon> AL = new ArrayList<Exon>();
			for(int i = 0; i < sortedKeys.size();i++)
			{
				if(sortedKeys.get(i) > this.filters.get(this.filters.size()-1).right) 
					return;
				else if(sortedKeys.get(i) > this.filters.get(pointer).right){

					i--;
					pointer++;
					count = 0;
				}
				else if(sortedKeys.get(i) >= this.filters.get(pointer).left && sortedKeys.get(i) <= this.filters.get(pointer).right){
					if(SVs.get(sortedKeys.get(i)).isHeterozygous(sample)){
						if(filterName == null){
							filterName = this.filters.get(pointer).name;
						}else{
							if(!SVs.get(sortedKeys.get(i)).Samples.get(sample).phased){
								count++;
							}

						}
					}
					//EW.println(this.filters.get(pointer).name+"_"+count+"\t"+this.Name+"_"+sortedKeys.get(i)+"\t"+SVs.get(sortedKeys.get(i)).getMotherCount(sample)+"\t"+SVs.get(sortedKeys.get(i)).getFatherCount(sample));
				}
			}
		}
	}





	public void removeHomozygous(ArrayList<String> samples){
		if(SVs!=null){
			List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
			int pointer = 0;
			for(int i = 0; i < sortedKeys.size();i++)
			{
				if(SVs.get(sortedKeys.get(i)).isHomozygous(samples))
					SVs.remove(sortedKeys.get(i));
			}
		}
	}

	public void removeHeterozygous(ArrayList<String> samples){
		if(SVs!=null){
			List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
			int pointer = 0;
			for(int i = 0; i < sortedKeys.size();i++)
			{
				if(!SVs.get(sortedKeys.get(i)).isHomozygous(samples))
					SVs.remove(sortedKeys.get(i));
			}
		}
	}

	public void printVCFinfoSamples(ArrayList<String> samples,ExtendedWriter EW){
		if(SVs!=null){
			List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
			for(int i = 0; i < sortedKeys.size();i++)
			{
				EW.print(this.Name+"\t"+sortedKeys.get(i)+"\t");
				SVs.get(sortedKeys.get(i)).printSamples(EW, samples);
				EW.println();
			}
		}
	}

	
	
	public void removeCounts(String sample){
		if(SVs!=null){
			List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
			for(int i = 0; i < sortedKeys.size();i++)
			{
				SVs.get(sortedKeys.get(i)).removeCounts(sample);
			}
		}
		
	}

	public void addCounts(Integer location, String sample, int count, boolean mother){
		if(SVs.containsKey(location))
				SVs.get(location).addCounts(sample,count,mother);
	}
	
	
	
	public void parseMpileUpFile(ExtendedReader ER,String sample, String sep, String parent){
		if(SVs!=null){
			List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
			for(int i = 0; i < sortedKeys.size();i++){
				
				
			}
		}
	}
	
	
	public void printVCFinfoSamples(ArrayList<String> samples,ExtendedWriter EW,String extraInfo){
		if(SVs!=null){
			List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
			for(int i = 0; i < sortedKeys.size();i++)
			{
				EW.print(this.Name+extraInfo+"\t"+sortedKeys.get(i)+"\t");
				SVs.get(sortedKeys.get(i)).printSamples(EW, samples);
				EW.println();
			}
		}
	}

	public void printVCFHetinfoSamples(String sample,ExtendedWriter EW){
		if(SVs!=null){
			List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
			for(int i = 0; i < sortedKeys.size();i++)
			{
				if(!SVs.get(sortedKeys.get(i)).getSample(sample).isHomozygous()){
					EW.print(this.Name+"\t"+sortedKeys.get(i)+"\t");
					SVs.get(sortedKeys.get(i)).printSample(EW, sample);
					EW.println();
				}
			}
		}
	}

	public void printVCFHetinfoSamplesSpecial(String sample,ExtendedWriter EW,String extraInfo){
		if(SVs!=null){
			List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
			for(int i = 0; i < sortedKeys.size();i++)
			{
				if(!SVs.get(sortedKeys.get(i)).getSample(sample).isHomozygous()){
					EW.print(this.Name+extraInfo+"\t"+sortedKeys.get(i)+"\t");
					SVs.get(sortedKeys.get(i)).printSample(EW, sample);
					EW.println();
				}
			}
		}
	}



	public void printVCFinfoSamplesRfriendly(ArrayList<String> samples,ExtendedWriter EW){
		if(SVs!=null){
			List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
			Collections.sort(sortedKeys);
			for(int i = 0; i < sortedKeys.size();i++)
			{
				EW.print(this.Name+"\t"+sortedKeys.get(i)+"\t");
				SVs.get(sortedKeys.get(i)).printSamplesRfriendly(EW, samples);
				EW.println();
			}
		}
	}




	//	public void printVCFinfoRfriendly(String sample,ExtendedWriter EW){
	//		if(SVs!=null){
	//			List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
	//			Collections.sort(sortedKeys);
	//			for(int i = 0; i < sortedKeys.size();i++)
	//			{
	//				EW.print(this.Name+"\t"+sortedKeys.get(i)+"\t");
	//				SVs.get(sortedKeys.get(i)).printSamples(EW, sample);
	//				EW.println();
	//			}
	//		}
	//	}
	//
	//
	//		public void printVCFinfoSampleDistribution(String sample,ExtendedWriter EW, int stepSize){
	//			if(SVs!=null){
	//				int start = 0; 
	//				int stop = stepSize;
	//				List<Integer> sortedKeys=new ArrayList<Integer>(SVs.keySet());
	//				Collections.sort(sortedKeys);
	//				for(int i = 0; i < sortedKeys.size();i++)
	//				{	
	//					while
	//					EW.print(this.Name+"\t"+sortedKeys.get(i)+"\t");
	//					SVs.get(sortedKeys.get(i)).printSample(EW, sample);
	//					EW.println();
	//				}
	//			}
	//		}
	//
	//	

	public void compareDistribution(Database otherRun,ExtendedWriter ER){

		double nrOfHits = this.getNrOfHits();
		double otherNrOfHits = otherRun.getNrOfHits();

		double difference = 0;
		if(nrOfHits > 10 || otherNrOfHits > 10){
			if(nrOfHits < 10)
				nrOfHits = 10;
			else if(otherNrOfHits < 10)
				otherNrOfHits = 10;
			else
				difference = Math.log(otherNrOfHits/nrOfHits);
		}
		if(difference != 0)
			ER.println(this.Name+","+difference+","+nrOfHits+","+otherNrOfHits);
	}

	private void sortHits(){
		for(int i = 1; i < this.hits.size();i++){
			int location = findHit(i, this.hits.get(i));
			if(location != -1){
				this.hits.get(location).nrOfReads++;
			}
		}
		this.hits.trimToSize();
	}

	public void removeNonRedundantHits(int cutoff){
		findRedundancy();
		removeBelowCutoff(cutoff);
	}

	public void removeAntisenseHits(int cutoff, int surrounding){
		if(this.hits != null){
			ArrayList <Hit> newHits = new ArrayList<Hit>();
			for(int i = 0; i < this.hits.size();i++){
				newHits.add(this.hits.get(i));
			}
			findRedundancy();
			removeDuplicates();
			findNonAntisenseHits(cutoff, surrounding);
			for(int j = 0; j < newHits.size();j++){
				boolean found = false;
				int pointer = 0;
				int start = newHits.get(j).start;
				boolean plusStrand = newHits.get(j).plusStrand;
				while(pointer < this.hits.size() && !found){
					if(this.hits.get(pointer).sameLocation(plusStrand, start))
						found = true;
					pointer++;
				}
				if(!found){
					newHits.remove(j);
					j--;
				}
			}
			this.hits = newHits;
		}
	}

	public void printHits(ExtendedWriter EW){
		if(this.hits != null){
			for(int i = 0; i < this.hits.size();i++){
				this.hits.get(i).printHit(EW);
			}
		}
	}


	public boolean printHits(String Dir, String file){
		try{
			ExtendedWriter EW = null;
			EW = new ExtendedWriter(new FileWriter(Dir+"/"+file+"."+this.Name+".rmapper"));
			for(int i = 0; i < this.hits.size();i++){
				this.hits.get(i).printHit(EW);
			}
			EW.flush();
			EW.close();
		}
		catch(Exception E){
			E.printStackTrace(); 
		}
		return false;
	}

	private void findRedundancy(){
		if(this.hits != null){
			for(int i = 1; i < this.hits.size();i++){
				findSameHits(i, this.hits.get(i));
			}
		}
	}

	private void removeBelowCutoff(int cutoff){
		if(this.hits != null){
			for(int i = 1; i < this.hits.size();i++){
				if(this.hits.get(i).nrOfReads < cutoff){
					this.hits.remove(i);
					i--;
				}
			}
			this.hits.trimToSize();
		}
	}


	private void findNonAntisenseHits(int cutoff,int surrounding){
		if(this.hits != null){
			ArrayList <Hit> newHits = new ArrayList<Hit>();
			for(int i = 0; i < this.hits.size();i++){
				int nrOfReads = this.hits.get(i).nrOfReads;
				boolean plusStrand = this.hits.get(i).plusStrand;
				int location = this.hits.get(i).getLocation();
				int nrOfAntisense = 0;
				for(int j = 0; j < this.hits.size();j++){
					if(this.hits.get(j).isAntisense( surrounding, location, plusStrand))
						nrOfAntisense += this.hits.get(j).nrOfReads;
				}
				if(nrOfAntisense == 0 || nrOfReads/nrOfAntisense > cutoff)newHits.add(this.hits.get(i));
			}
			this.hits = newHits;
		}
	}



	private void removeDuplicates(){
		if(this.hits != null){
			int count = 1;
			for(int i = 1; i < this.hits.size();i++){
				if(i*100/this.hits.size()>count){
					count++;
				}
				int location = findHit(i, this.hits.get(i));
				if(location != -1){
					this.hits.remove(i);
					i--;
				}
			}		
			this.hits.trimToSize();
		}
	}

	public int countLocations(){
		removeDuplicates();
		if(this.hits != null)
			return this.hits.size();
		return 0;
	} 




	public int findHit(int stop, Hit newHit){
		for(int i = 0; i < stop;i++){
			if(this.hits.get(i).hasSameLocation(newHit)){
				newHit.nrOfReads++;
				return i;
			}
		}
		return -1;
	}

	public void findSameHits(int stop, Hit newHit){
		for(int i = 0; i < stop;i++){
			if(this.hits.get(i).hasSameLocation(newHit)){
				newHit.nrOfReads++;
				this.hits.get(i).nrOfReads++;
			}
		}
	}






	public void printHitSequence(int start, int length, boolean plusStrand, int width, String geneName){
		int[] surrSequence = null;
		if(plusStrand){
			int seqStart = start-width;
			int seqEnd = start+length+width;
			if(seqStart < 0)seqStart = 0;
			if(seqEnd >= this.sequence.length) seqEnd = this.sequence.length;
			surrSequence = Functions.getSubarray(this.sequence, seqStart, seqEnd);
			System.out.println(">"+this.Name+"_"+start+"_"+length+"_+_"+width+"_nt_upstreamAndDownstream");
		}
		else{
			int seqStart = start-width;
			int seqEnd = start+length+width;
			if(seqStart < 0)seqStart = 0;
			if(seqEnd >= this.sequence.length) seqEnd = this.sequence.length;
			surrSequence = RNAfunctions.getComplementary(
					Functions.getSubarray(this.sequence, seqStart, seqEnd));
			System.out.println(">"+this.ID+"_"+start+"_"+length+"_-_"+width+"_nt_upstreamAndDownstream");
		}
		System.out.println(RNAfunctions.RNAInt2String(surrSequence));


		if(plusStrand){
			int seqStart = start;
			int seqEnd = start+length;
			if(seqStart < 0)seqStart = 0;
			if(seqEnd >= this.sequence.length) seqEnd = this.sequence.length;
			surrSequence = Functions.getSubarray(this.sequence, seqStart, seqEnd);
			System.out.println(">"+this.ID+"_"+start+"_"+length+"_+");
		}
		else{
			int seqStart = start;
			int seqEnd = start+length;
			if(seqStart < 0)seqStart = 0;
			if(seqEnd >= this.sequence.length) seqEnd = this.sequence.length;
			surrSequence = RNAfunctions.getComplementary(
					Functions.getSubarray(this.sequence, seqStart, seqEnd));
			System.out.println(">"+this.ID+"_"+(start+length)+"_"+length+"_-");
		}
		System.out.println(RNAfunctions.RNAInt2String(surrSequence));
	}

	public void printSurroundingSequence(int start, int length, boolean plusStrand, int width, String geneName,int nrOfReads, ExtendedWriter EW){
		int[] surrSequence = null;
		if(plusStrand){
			int seqStart = start-width-1;
			int seqEnd = start+length+width-1;
			if(seqStart < 0)seqStart = 0;
			if(seqEnd >= this.sequence.length) seqEnd = this.sequence.length;
			surrSequence = Functions.getSubarray(this.sequence, seqStart, seqEnd);
			EW.println(">"+this.Name+","+geneName+","+nrOfReads+","+start+","+length+",+,"+width+"_nt_upstreamAndDownstream");
			System.out.println(">"+this.Name+","+geneName+","+nrOfReads+","+start+","+length+",+,"+width+"_nt_upstreamAndDownstream");

		}
		else{
			int seqStart = start-width-1;
			int seqEnd = start+length+width-1;
			if(seqStart < 0)seqStart = 0;
			if(seqEnd >= this.sequence.length) seqEnd = this.sequence.length;
			surrSequence = RNAfunctions.getComplementary(
					Functions.getSubarray(this.sequence, seqStart, seqEnd));
			EW.println(">"+this.Name+","+geneName+","+nrOfReads+","+(start+21-1)+","+length+",-,"+width+"_nt_upstreamAndDownstream");
			System.out.println(">"+this.Name+","+geneName+","+nrOfReads+","+(start+21-1)+","+length+",-,"+width+"_nt_upstreamAndDownstream");
		}
		EW.println(RNAfunctions.RNAInt2String(surrSequence));
		System.out.println(RNAfunctions.RNAInt2String(surrSequence));

	}

	public void printSolidSequence(Solid hit, ExtendedWriter EW){
		int USlength = 300;
		for(int i = 0; i < hit.hits.size(); i++){
			if(hit.hits.get(i).chromosome.compareTo(this.Name) == 0){
				int start = hit.hits.get(i).start;
				int length = hit.hits.get(i).length;
				int[] surrSequence = null;
				if(start > 0 ){
					int seqStart = start-USlength;
					int seqEnd = start+length+USlength;
					if(seqStart < 0)seqStart = 0;
					if(seqEnd >= this.sequence.length) seqEnd = this.sequence.length;
					surrSequence = Functions.getSubarray(this.sequence, seqStart, seqEnd);
				}
				else{
					int seqStart = this.length+start-length-USlength;
					int seqEnd = this.length+start+USlength;
					if(seqStart < 0)seqStart = 0;
					if(seqEnd >= this.sequence.length) seqEnd = this.sequence.length;
					surrSequence = RNAfunctions.getComplementary(
							Functions.getSubarray(this.sequence, seqStart, seqEnd));
				}
				EW.println(hit.getFastaHitInfo(i)+"_"+USlength+"_nt_upstreamAndDownstream");
				EW.println(RNAfunctions.RNAInt2String(surrSequence));
			}
		}
	}



	private void getChromosomeInfo(String[] columns){
		int left = Integer.parseInt(columns[3]);
		int right = Integer.parseInt(columns[4]);
		this.setLength(right - left+1);

		String [] extra = columns[8].split(";");
		if(extra != null){
			for(int i = 0; i < extra.length; i++){
				if(extra[i].indexOf("ID=") == 0){
					String[] IDs = extra[i].split("=");

					this.ID = IDs[1];

				}
				if(extra[i].indexOf("Name=") == 0){
					String[] IDs = extra[i].split("=");
					this.Name = IDs[1];
				}
			}
		}


	}

	public void setName(String name) {
		Name = name;
	}

	public String getName() {
		return Name;
	}

	public void setLength(int length) {
		this.length = length;
	}

	public int getLength() {
		return length;
	}

	public void setSequence(int[] sequence) {
		this.sequence = sequence;
	}

	public int[] getSequence() {
		return sequence;
	}


	void addGene(String[] columns){
		int left = Integer.parseInt(columns[3]);
		int right = Integer.parseInt(columns[4]);
		boolean plusStrand = true;
		if(columns[6].indexOf("-") > -1)
			plusStrand = false;
		String [] extra = columns[8].split(";");
		String ID = "wrong";
		String Name = "";
		String description = "";
		if(extra != null){
			for(int i = 0; i < extra.length; i++){
				if(extra[i].indexOf("ID=") == 0){
					String[] IDs = extra[i].split("=");
					ID = IDs[1];
				}
				if(extra[i].indexOf("Name=") == 0){
					String[] IDs = extra[i].split("=");
					Name = IDs[1];
				}
				if(extra[i].indexOf("description=") == 0){
					String[] IDs = extra[i].split("=");
					description = IDs[1];
				}
			}
		}
		if(ID.indexOf("wrong")==-1){
			this.ncRNAs.add(new ncRNA(left, right, plusStrand,ID,Name,description));
		}
		else{
			System.out.println("Something wrong when adding gene!");
			for(int i = 0; i < columns.length; i++){
				System.out.print(columns[i]+"\t");
			}
			System.out.println();
		}
	}

	private void addRepeat(String[] columns){
		if(columns [0].compareTo(this.ID) == 0){
			int left = Integer.parseInt(columns[3]);
			int right = Integer.parseInt(columns[4]);
			boolean plusStrand = true;
			if(columns[6].indexOf("-") > -1)
				plusStrand = false;
			String [] extra = columns[8].split("\"");
			String ID = "wrong";
			if(extra != null){
				ID = extra[1];
			}
			if(ID.indexOf("wrong")==-1)this.repeats.add(new Repeat(left, right, plusStrand,ID));
			else{
				System.out.println("Something wrong when adding repeat!");
				for(int i = 0; i < columns.length; i++){
					System.out.print(columns[i]+"\t");
				}
				System.out.println();
			}
		}

	}


	private void addNCRNA(String[] columns){
		if(columns[0].compareTo("Dictyostelium_discoideum_chromosome_"+this.Name) == 0){
			int left = Integer.parseInt(columns[3]);
			int right = Integer.parseInt(columns[4]);
			boolean plusStrand = true;
			if(columns[6].indexOf("-") > -1)
				plusStrand = false;
			String ID = columns[8];
			String Name = columns[8];
			String description = "GC rich ncRNA";
			this.ncRNAs.add(new ncRNA(left, right, plusStrand,ID,Name,description));
		}
	}



	private boolean addmRNA(String[] columns){
		int left = Integer.parseInt(columns[3]);
		int right = Integer.parseInt(columns[4]);
		boolean plusStrand = true;
		if(columns[6].indexOf("-") > -1)
			plusStrand = false;
		String [] extra = columns[8].split(";");
		String ID = "wrong";
		String Name = "";
		String description = "";
		String parent = "";
		if(extra != null){
			for(int i = 0; i < extra.length; i++){
				if(extra[i].indexOf("ID=") == 0){
					String[] IDs = extra[i].split("=");
					ID = IDs[1];
				}
				if(extra[i].indexOf("Name=") == 0){
					String[] IDs = extra[i].split("=");
					Name = IDs[1];
				}
				if(extra[i].indexOf("Parent=") == 0){
					String[] IDs = extra[i].split("=");
					parent = IDs[1];
				}
			}
		}
		if(ID.indexOf("wrong")==-1){
			mRNA newmRNA = new mRNA(left, right, plusStrand,ID,Name,parent,description);
			for(int i = this.ncRNAs.size()-1;i > -1; i--){
				if(this.ncRNAs.get(i).isParent(newmRNA)){
					CodingGene newCG = new CodingGene();
					newCG.setInfo(ncRNAs.get(i));
					newCG.addmRNA(newmRNA);
					this.codingGenes.add(newCG);
					this.ncRNAs.remove(i);
					this.ncRNAs.trimToSize();
					i--;
					return true;
				}
			}
			for(int i = this.codingGenes.size()-1;i > -1; i--){
				CodingGene temp = (CodingGene)this.codingGenes.get(i);
				if(temp.addmRNA(newmRNA)){
					this.codingGenes.remove(i);
					this.codingGenes.add(i,temp);
					return true;
				}
			}

		}
		else{
			System.out.println("Something wrong when adding mRNA!");
			for(int i = 0; i < columns.length; i++){
				System.out.print(columns[i]+"\t");
			}
			System.out.println();
		}
		return false;

	}

	private boolean addExon(String[] columns){
		int left = Integer.parseInt(columns[3]);
		int right = Integer.parseInt(columns[4]);
		String parent = "wrong";
		String [] extra = columns[8].split(";");
		String ID = "wrong";
		String Name = "";
		String description = "";

		if(extra != null){
			for(int i = 0; i < extra.length; i++){
				if(extra[i].indexOf("ID=") == 0){
					String[] IDs = extra[i].split("=");
					ID = IDs[1];
				}
				if(extra[i].indexOf("Name=") == 0){
					String[] IDs = extra[i].split("=");
					Name = IDs[1];
				}
				if(extra[i].indexOf("Parent=") == 0){
					String[] IDs = extra[i].split("=");
					parent = IDs[1];
				}
			}
		}
		if(parent.indexOf("wrong")==-1){
			Exon newExon = new Exon(left, right, parent);
			codingGenes.trimToSize();
			for(int i = this.codingGenes.size()-1;i > -1; i--){
				CodingGene temp = (CodingGene)this.codingGenes.get(i);
				if(temp.addExon(newExon)){
					this.codingGenes.remove(i);
					this.codingGenes.add(i,temp);
					return true;
				}
			}
		}
		else{
			System.out.println("Something wrong when adding exon!");
			for(int i = 0; i < columns.length; i++){
				System.out.print(columns[i]+"\t");
			}
			System.out.println();
		}
		/*			System.out.println("Something wrong when adding exon!");
			for(int i = 0; i < columns.length; i++){
				System.out.print(columns[i]+"\t");
			}
			System.out.println();

		 */
		return false;

	}



}




