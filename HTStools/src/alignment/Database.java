package alignment;

import general.Functions;
import general.RNAfunctions;
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
		int[] sequence = new int[10000000];
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
		// scaffold_1      5       .       A       T       56.41   .       AC=3;AF=0.188;AN=16;BaseQRankSum=-1.804;DP=72;Dels=0.00;FS=0.000;HaplotypeScore=2.4206;MLEAC=2;MLEAF=0.125;MQ=26.05;MQ0=11;MQRankSum=0.618;QD=2.69;ReadPosRankSum=1.053 GT:AD:DP:GQ:PL  0|0:14,0:14:30:0,30,268 0|0:1,0:1:3:0,3,29      0/1:13,5:18:72:72,0,247 0|0:14,0:14:39:0,39,316 0|0:13,0:13:39:0,39,328 0|0:4,0:4:3:0,3,25      1|1:2,1:3:3:23,3,0      0|0:5,0:5:3:0,3,25

		
		if(this.SVs == null){
			this.SVs = new Hashtable<Integer,StructuralVariation>();
		}
		Integer Location = Integer.decode(VCFinfo[1]);
		if(VCFinfo[5].compareTo(".")!=0)
			this.SVs.put(Location, new StructuralVariation(Samples,VCFinfo));
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



	public void addBEDfilterInfo( String[] BEDinfo){
		//	0			1		2     3		  4			5		6		7		8		9			 
		// #CHROM(0)  Start     Stop     Name     .     .    INFO  type    Something  XTR 	AInfo   
		// scaffold_1      767     2124    PAC:20891551.exon.3     .       -       phytozome8_0    exon    .       ID=PAC:20891551.exon.3;Parent=PAC:20891551;pacid=20891551


		if(this.filters == null){
			this.filters = new ArrayList <Exon>();
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




	public void printFilters(){
		if(this.filters != null)
			for(int i = 0; i< this.filters.size();i++){
				filters.get(i).printBED(this.Name);
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
					System.out.println("removing "+sortedKeys.get(i));
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
				if(sortedKeys.get(i) <= this.filters.get(pointer).right && sortedKeys.get(i) >= this.filters.get(pointer).left){
					SVs.remove(sortedKeys.get(i));
					System.out.println("removing "+sortedKeys.get(i));
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



}




