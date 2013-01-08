package alignment;

import general.Functions;
import general.IOTools;
import general.RNAfunctions;
import general.energy;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;

import Sequence.FastaSequence;
import Sequence.Solid;

import sun.tools.tree.ThisExpression;
import general.ExtendedReader;
import general.ExtendedWriter;




public class Contig implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//private String GenomeName;
	private String Name;
	private String ID;
	
	protected ArrayList <Gene> codingGenes;
	protected ArrayList <Gene> ncRNAs;
	protected ArrayList <Gene> repeats;
	protected ArrayList <Gene> intergenicRegions;
	
	private int length;
	private int[] sequence;
	public FastaSequence seq;
	
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
			Contig test = new Contig(gffDir,gffFile,Chromosome);
			int start = Integer.parseInt(Functions.getValue(T, "-start", "300"));
			int length2 = Integer.parseInt(Functions.getValue(T, "-length", "21"));
			int width = Integer.parseInt(Functions.getValue(T, "-width", "200"));
			boolean plusStrand = Boolean.parseBoolean(Functions.getValue(T, "-plusStrand", "true"));
			test.printHitSequence(start, length2, plusStrand, width,"genome");
			
	}	
	
	Contig(String Genomedir, String fileName,String Name){
		this.setName(Name);
		this.codingGenes = new ArrayList<Gene>();
		//printGenes();
		
	}

	Contig(String Name){
		this.setName(Name);
		codingGenes = new ArrayList <Gene>(); 
		
	}


	public void sortGenes(){
		this.codingGenes = sortGenes(this.codingGenes);
		this.ncRNAs = sortGenes(this.ncRNAs);
		this.repeats = sortGenes(this.repeats);
		this.intergenicRegions = addIntergenic();
	}
	
	
	public void addgetORFinfo(FastaSequence FS){
		Gene A = new Gene();
		String[] info = FS.getName().split("\\ ");
		int start =  Integer.parseInt(info[1].substring(1));
		int stop =  Integer.parseInt(info[3].substring(0,info[3].length()-1));
		if(start < stop){
			A.Name = info[0];
			A.left = start;
			A.right = stop;
			A.plusStrand = true;
		}
		else{
			A.left = stop;
			A.right = start;
			A.plusStrand = false;
		}
		A.setFastaSeq(FS);
		this.codingGenes.add(A);
	}
	
	public void printORFs(ExtendedWriter EW){
		for(int i = 0; i < this.codingGenes.size();i++){
			this.codingGenes.get(i).getFastaSeq().printFasta(EW);
		}
	}
	
	public void printORFNames(ExtendedWriter EW){
		for(int i = 0; i < this.codingGenes.size();i++){
			this.codingGenes.get(i).getFastaSeq().printName(EW);
		}
	}
	
	
	
	public int sortORFs(){

		ArrayList <Gene> SortedCodingGenes = new ArrayList<Gene>();
		while(this.codingGenes.size()>0){
			int pointer = 0;
			int maxSize = 0;
			for(int i = 0; i < this.codingGenes.size();i++){
				if(this.codingGenes.get(i).size() > maxSize){
					maxSize = this.codingGenes.get(i).size();
					pointer = i;
				}
			}
			Gene NewGene = this.codingGenes.get(pointer);
			this.codingGenes.remove(pointer);	
			boolean add =true;
			for(int i = 0; i < SortedCodingGenes.size();i++){
				if(SortedCodingGenes.get(i).isOverlapping(NewGene)) add = false;
			}
			if(add) {
				SortedCodingGenes.add(NewGene);
			}
		}
		this.codingGenes = SortedCodingGenes;
		return this.codingGenes.size();
	}
	
	
	public int longestORFs(){

		ArrayList <Gene> SortedCodingGenes = new ArrayList<Gene>();
		int pointer = 0;
		int maxSize = 0;
		for(int i = 0; i < this.codingGenes.size();i++){
			if(this.codingGenes.get(i).size() > maxSize){
				maxSize = this.codingGenes.get(i).size();
				pointer = i;
			}
		}
		Gene NewGene = this.codingGenes.get(pointer);
		SortedCodingGenes.add(NewGene);
		this.codingGenes = SortedCodingGenes;
		return this.codingGenes.size();
	}
	
	
	
	public int getNrOfHits(){
		int total =0;
		for(int i = 0; i < codingGenes.size(); i++){
			total += codingGenes.get(i).getNrOfHits();
		}
		for(int i = 0; i < ncRNAs.size(); i++){
			total += ncRNAs.get(i).getNrOfHits();
		}
		for(int i = 0; i < repeats.size(); i++){
			total += repeats.get(i).getNrOfHits();
		}
		for(int i = 0; i < intergenicRegions.size(); i++){
			total += intergenicRegions.get(i).getNrOfHits();
		}
		return total;
	}


	public void compareDistribution(Contig otherRun,ExtendedWriter ER){
		
		for(int i = 0; i < codingGenes.size(); i++){
			if(codingGenes.get(i).Name.compareTo(otherRun.codingGenes.get(i).Name ) != 0)
				System.out.println(codingGenes.get(i).Name+ "  "+ otherRun.codingGenes.get(i).Name);
			
			codingGenes.get(i).compareNROfHits(this.Name,otherRun.codingGenes.get(i).getNrOfHits(), ER);
		}
		for(int i = 0; i < ncRNAs.size(); i++){
			if(ncRNAs.get(i).Name.compareTo(otherRun.ncRNAs.get(i).Name ) != 0)
				System.out.println(ncRNAs.get(i).Name+ "  "+ otherRun.ncRNAs.get(i).Name);
			ncRNAs.get(i).compareNROfHits(this.Name,otherRun.ncRNAs.get(i).getNrOfHits(), ER);
		}
		for(int i = 0; i < repeats.size(); i++){
			repeats.get(i).compareNROfHits(this.Name,otherRun.repeats.get(i).getNrOfHits(), ER);
			if(repeats.get(i).Name.compareTo(otherRun.repeats.get(i).Name ) != 0)
				System.out.println(repeats.get(i).Name+ "  "+ otherRun.repeats.get(i).Name);
		}
		for(int i = 0; i < intergenicRegions.size(); i++){
			intergenicRegions.get(i).compareNROfHits(this.Name,otherRun.intergenicRegions.get(i).getNrOfHits(), ER);
			if(intergenicRegions.get(i).Name.compareTo(otherRun.intergenicRegions.get(i).Name ) != 0)
				System.out.println(intergenicRegions.get(i).Name+ "  "+ otherRun.intergenicRegions.get(i).Name);
		}
	}
	
	
	public void addGFF3info(String dir, String fileName){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));

			while(ER.more()){
				if((char)ER.lookAhead() == '#'){
					String GFF3line = ER.readLine();
					if(GFF3line.indexOf("FASTA") != -1) readChromosomeSequence(ER);
				}
				else{
					readInfo(ER);
				}
			}
		}catch(Exception E){E.printStackTrace();}
	}
	
	private void readChromosomeSequence(ExtendedReader ER){
		int pointer = 0;
		int[] sequence = new int[this.getLength()];
		ER.skipLine();
//		System.out.println("LŠser in kromosom sekvensen");
		while(ER.more()){
			int[] subSequence = RNAfunctions.RNAString2Int(ER.readLine());
			for(int i = 0; i < subSequence.length;i++){
				sequence[pointer] = subSequence[i];
				pointer++;
			}
		}
		this.sequence = sequence;
//		System.out.println(pointer);
//		System.out.println("Name: "+this.Name);
//		System.out.println("length: "+this.getLength());
//		System.out.println("Klar");
	}
	
	
	public void mapSolidSequence(Solid hit){
		for(int i = 0; i < hit.hits.size(); i++){
			if(hit.hits.get(i).chromosome.compareTo(this.Name) == 0){
				hit.hits.get(i).specificGene = this.getKind(hit.hits.get(i));
			}
		}
	}

	
	public void mapSolidSequenceChromosome(Solid hit){
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
		for(int i = 1; i < this.hits.size();i++){
			int location = findHit(i, this.hits.get(i));
			if(location != -1){
				this.hits.remove(i);
				i--;
			}
		}
		this.hits.trimToSize();
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
	

	public void print3UTR(int start, int length, boolean plusStrand, String geneName, ExtendedWriter EW){
		int[] surrSequence = null;
		if(plusStrand){
			int seqStart = start-1;
			int seqEnd = start+length+1;
			if(seqEnd >= this.sequence.length) seqEnd = this.sequence.length;
			surrSequence = Functions.getSubarray(this.sequence, seqStart, seqEnd);
			EW.println(">"+geneName+"UTR,"+this.Name+","+start+","+length+"nt ,+");
			//System.out.println(">"+this.Name+","+geneName+","+nrOfReads+","+start+","+length+",+,"+width+"_nt_upstreamAndDownstream");

		}
		else{
			int seqStart = start-length;
			int seqEnd = start;
			if(seqStart < 0)seqStart = 0;
			if(seqEnd >= this.sequence.length) seqEnd = this.sequence.length;
			surrSequence = RNAfunctions.getComplementary(
					Functions.getSubarray(this.sequence, seqStart, seqEnd));
			EW.println(">"+geneName+"UTR,"+this.Name+","+start+","+length+"nt ,-");
			//System.out.println(">"+this.Name+","+geneName+","+nrOfReads+","+(start+21-1)+","+length+",-,"+width+"_nt_upstreamAndDownstream");
		}
		EW.println(RNAfunctions.RNAInt2String(surrSequence));
		System.out.println(RNAfunctions.RNAInt2String(surrSequence));
		
	}
	
	
	public void print3UTRtab(int start, int length, boolean plusStrand, String geneName, ExtendedWriter EW){
		int[] surrSequence = null;
		if(plusStrand){
			int seqStart = start-1;
			int seqEnd = start+length+1;
			if(seqEnd >= this.sequence.length) seqEnd = this.sequence.length;
			surrSequence = Functions.getSubarray(this.sequence, seqStart, seqEnd);
			EW.print(geneName+"\t1000\t");
			//System.out.println(">"+this.Name+","+geneName+","+nrOfReads+","+start+","+length+",+,"+width+"_nt_upstreamAndDownstream");

		}
		else{
			int seqStart = start-length;
			int seqEnd = start;
			if(seqStart < 0)seqStart = 0;
			if(seqEnd >= this.sequence.length) seqEnd = this.sequence.length;
			surrSequence = RNAfunctions.getComplementary(
					Functions.getSubarray(this.sequence, seqStart, seqEnd));
			EW.print(geneName+"\t1000\t");
			//System.out.println(">"+this.Name+","+geneName+","+nrOfReads+","+(start+21-1)+","+length+",-,"+width+"_nt_upstreamAndDownstream");
		}
		EW.println(RNAfunctions.RNAInt2String(surrSequence));
//		System.out.println(RNAfunctions.RNAInt2String(surrSequence));
		
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
	
	
	
	public void printDistributionSolidSequence(ExtendedWriter EW){
		for(int i = 0; i < this.ncRNAs.size(); i++){
			if(this.ncRNAs.get(i).getNrOfHits() > 0 ){
				this.ncRNAs.get(i).printNrOfDirectedHits(EW,this.Name, this.length);
			}
		}
		for(int i = 0; i < this.repeats.size(); i++){
			if(this.repeats.get(i).getNrOfHits() > 0 ){
				this.repeats.get(i).printNrOfDirectedHits(EW,this.Name, this.length);
			}
		}
		
		for(int i = 0; i < this.intergenicRegions.size(); i++){
			if(this.intergenicRegions.get(i).getNrOfHits() > 0 ){
				this.intergenicRegions.get(i).printNrOfDirectedHits(EW,this.Name, this.length);
			}
		}
		
		for(int i = 0; i < this.codingGenes.size(); i++){
			if(this.codingGenes.get(i).getNrOfHits() > 0 ){
				this.codingGenes.get(i).printNrOfDirectedHits(EW,this.Name, this.length);
			}
		}
	}

	public void printHits(String dir, String file){
		if(!IOTools.isDir(dir+"/ncRNAs"))
			IOTools.mkDir(dir+"/ncRNAs");
		for(int i = 0; i < this.ncRNAs.size(); i++){
			if(this.ncRNAs.get(i).getNrOfHits() > 0 ){
				if(!this.ncRNAs.get(i).printpremiRNAstructures(dir+"/ncRNAs", file)){
					this.ncRNAs.remove(i);
					i--;
					
				}
			}
		}
		if(!IOTools.isDir(dir+"/repeats"))
			IOTools.mkDir(dir+"/repeats");
		for(int i = 0; i < this.repeats.size(); i++){
			if(this.repeats.get(i).getNrOfHits() > 0 ){
				if(!this.repeats.get(i).printpremiRNAstructures(dir+"/repeats", file)){
					this.repeats.remove(i);
					i--;
				}
			}
		}
		
		if(!IOTools.isDir(dir+"/intergenic"))
			IOTools.mkDir(dir+"/intergenic");
		for(int i = 0; i < this.intergenicRegions.size(); i++){
			if(this.intergenicRegions.get(i).getNrOfHits() > 0 ){
				if(!this.intergenicRegions.get(i).printpremiRNAstructures(dir+"/intergenic", file)){
					this.intergenicRegions.remove(i);
					i--;
				}
			}
		}
		
		if(!IOTools.isDir(dir+"/coding"))
			IOTools.mkDir(dir+"/coding");
		for(int i = 0; i < this.codingGenes.size(); i++){
			if(this.codingGenes.get(i).getNrOfHits() > 0 ){
				if(!this.codingGenes.get(i).printpremiRNAstructures(dir+"/coding", file)){
					this.codingGenes.remove(i);
					i--;
				}
			}
		}
	}
	

	
	
	
	
	
	
	private Gene getKind(Hit newHit){
		int location = newHit.start;
		int length = newHit.length;
		String Name = newHit.getFastaInfo();
		if(location < 0)
			location = this.length +location - length/2;
		else
			location = location + length/2;
		double loc = findLocation(location,this.repeats);
		if(loc > -0.4 &&(int)loc - (int)(loc + 0.5) == 0){
			this.repeats.get((int)loc).addHit(newHit);
			return this.repeats.get((int)loc);
		}
		loc = findLocation(location,this.ncRNAs);
		if(loc > -0.4 &&(int)loc - (int)(loc + 0.5) == 0){
			this.ncRNAs.get((int)loc).addHit(newHit);
			return this.ncRNAs.get((int)loc);
		}
		loc = findLocation(location,this.codingGenes);
		if(loc > -0.4 &&(int)loc - (int)(loc + 0.5) == 0){
			this.codingGenes.get((int)loc).addHit(newHit);
			return this.codingGenes.get((int)loc);
		}
		loc = findLocation(location,this.intergenicRegions);
		if(loc > -0.4 &&(int)loc - (int)(loc + 0.5) == 0){
			this.intergenicRegions.get((int)loc).addHit(newHit);
			return this.intergenicRegions.get((int)loc);
		}
		else{
			System.out.println("Something is wrong with find location in Chromsome");
			System.out.println(location + " chromosome "+ this.Name);
			System.out.println(Name);
		}
		return null;
	}
	
	
	private static double findLocation(int location,ArrayList<Gene> sortedGenes){
		if(sortedGenes == null || sortedGenes.size() == 0)
			return -0.5;
		int leftPointer = 0;
		int rightPointer = sortedGenes.size()-1;
		int smallestLeft = sortedGenes.get(leftPointer).left;
		int largestRight = sortedGenes.get(rightPointer).right;
		if(location > largestRight){
			return (double)rightPointer+0.5;
		} 
		if(location < smallestLeft){
			return (double)leftPointer-0.5;
		} 
		while(rightPointer - leftPointer > 1){
			int middle = (rightPointer+leftPointer)/2;
			if(sortedGenes.get(middle).left> location) rightPointer = middle;
			else if(sortedGenes.get(middle).right < location) leftPointer = middle;
			else return middle;
		}
		if(sortedGenes.get(leftPointer).right> location)return leftPointer;
		else if(sortedGenes.get(rightPointer).left < location)return rightPointer;
		else return (double)leftPointer+0.5;
	}
	
	private static ArrayList <Gene> sortGenes(ArrayList <Gene> unsortedGenes ){
		ArrayList <Gene> sortedGenes = new ArrayList<Gene>();
		for(int i = 0; i < unsortedGenes.size(); i++){
			double pointer = findLocation((unsortedGenes.get(i).left+unsortedGenes.get(i).right)/2,sortedGenes);
			int location = (int)(pointer+0.5);
			sortedGenes.add(location, unsortedGenes.get(i));
			
		}
		return sortedGenes;
	}


	
	
	private  ArrayList <Gene> addIntergenic(){

		ArrayList <Gene> sortedGenes = new ArrayList<Gene>();
		int codingPointer = 0;
		int nextLeft = 0;
		int count = 0;
		while( codingPointer < codingGenes.size() ){
			boolean overlaping = true;
			while(overlaping){
				overlaping = false;
				while(codingPointer < codingGenes.size() && codingGenes.get(codingPointer).left < nextLeft){
					while (codingPointer < codingGenes.size() &&codingGenes.get(codingPointer).right <= nextLeft)codingPointer++;
					if(codingPointer < codingGenes.size() && codingGenes.get(codingPointer).left < nextLeft) nextLeft = codingGenes.get(codingPointer).right;
					overlaping = true;
				}
			}
			int right = this.length;
			int left = nextLeft;
			if(codingPointer < codingGenes.size() &&right > codingGenes.get(codingPointer).left){
				right = codingGenes.get(codingPointer).left;
				nextLeft = codingGenes.get(codingPointer).right;
			}
			sortedGenes.add(new Intergenic(left, right, "intergenic_"+count,"intergenic_"+count, "test"));
			count++;
		}
		if(sortedGenes.size() < 1)
			sortedGenes.add(new Intergenic(0, this.length, "intergenic_"+count,"intergenic_"+count, "test"));
		
		
		
		return sortedGenes;
	}

	private void readInfo(ExtendedReader ER){
		String GFF3line = ER.readLine();
		String[] columns = GFF3line.split("\t");
		if(columns.length > 2){
		if(columns[2].indexOf("gene") == 0)addGene(columns);
		else if(columns[2].indexOf("mRNA") == 0)addmRNA(columns);
		else if(columns[2].indexOf("exon") == 0)addExon(columns);
		else if(columns[2].indexOf("chromosome") == 0)getChromosomeInfo(columns);
		else if(columns[2].indexOf("dr") ==0)addNCRNA(columns);
		else if(columns[2].indexOf("Dictyostelium discoideum complex repeat") > 0)addRepeat(columns);
		
		}
		else
			System.out.println(GFF3line);
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
	
	
	private void addGene(String[] columns){
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
		if(ID.indexOf("wrong")==-1)this.ncRNAs.add(new ncRNA(left, right, plusStrand,ID,Name,description));
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
			System.out.println("Something wrong when adding gene!");
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
		if(columns[8].indexOf("Parent=") == 0){
				String[] IDs = columns[8].split("=");
				parent = IDs[1];
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
			System.out.println("Something wrong when adding gene!");
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




