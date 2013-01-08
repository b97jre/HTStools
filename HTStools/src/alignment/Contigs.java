package alignment;

import general.Functions;
import general.IOTools;

import java.io.FileWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;

import Sequence.CfastaSequences;
import Sequence.FastaSequence;
import Sequence.FastaSequences;
import Sequence.Solid;

import general.ExtendedWriter;




public class Contigs implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//private String GenomeName;
	private String Name;
	protected Hashtable <String,Contig> Chromosomes;

	public static void main(String[] args) {
		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}
		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);

		Contigs.run(T);

	}

	public static void run(Hashtable<String,String> T){
		String dir = Functions.getValue(T, "-d",".");
		String ContigFile = Functions.getValue(T, "-f",".");
		String OrfFile = Functions.getValue(T, "-i",".");
		String outFile = Functions.getValue(T, "-o",ContigFile+"_"+OrfFile+"_pruned.fa");

		Contigs A = new Contigs();
//		System.out.println("Adding contigs");
//		A.addContigsName(dir,ContigFile);
		System.out.println("Adding ORFs");
		A.addGetORFinfo(dir,OrfFile);
		System.out.println("removing overlapping ORFs");
		int NrOfORFS = 0;
		if(!T.containsKey("-longest"))
			NrOfORFS = A.sortContigs();
		else
			NrOfORFS = A.longestContigs();
		Hashtable<String, Gene> ORFs = new Hashtable<String, Gene>(NrOfORFS);
		System.out.println("printing ORFs");
		A.printORFNames(dir, outFile);
	}





	Contigs(){
		this.Chromosomes = new Hashtable <String,Contig>();
	}

	public Contigs(String Name){
		this.Name = Name;
		this.Chromosomes = new Hashtable <String,Contig>();
	}




	public void addContigs(String dir, String file){
		FastaSequences FS = new FastaSequences(dir,file);
		for(int i = 0; i < FS.size(); i++){
			Contig B = new Contig(FS.get(i).getName());
			B.seq = FS.get(i);
			String orfName =  FS.get(i).getName();
			if(orfName.indexOf("\\ ")>-1){
				orfName = orfName.split("\\ ")[0];
				if(orfName.indexOf("_") != orfName.lastIndexOf("_")){
					int index =  orfName.indexOf("_", orfName.indexOf("_")+1);
					orfName = orfName.substring(0,index);
				}
			}
			Chromosomes.put(orfName, B);
		}
	}

	public void addContigsName(String dir, String file){
		FastaSequences FS = new FastaSequences();
		FS.getFastaFileNames(dir,file);
		for(int i = 0; i < FS.size(); i++){
			Contig B = new Contig(FS.get(i).getName());
			B.seq = FS.get(i);
			String orfName =  FS.get(i).getName();
			if(orfName.indexOf("\\ ")>-1){
				orfName = orfName.split("\\ ")[0];
				if(orfName.indexOf("_") != orfName.lastIndexOf("_")){
					int index =  orfName.indexOf("_", orfName.indexOf("_")+1);
					orfName = orfName.substring(0,index);
				}
			}
			Chromosomes.put(orfName, B);
		}
	}

	public int sortContigs(){
		int numberOfORFs = 0;
		
		for(Contig val : this.Chromosomes.values()){
			numberOfORFs += val.sortORFs();
		}
		return numberOfORFs;
	}


	public int  longestContigs(){
		int numberOfORFs = 0;
		for(Contig val : this.Chromosomes.values()){
			numberOfORFs += val.longestORFs();
		}
		return numberOfORFs;
	}

	public void addGetORFinfo(String dir, String file){
		FastaSequences FS = new FastaSequences();
		FS.getFastaFileNames(dir,file);
		System.out.println("Find fasta sequences with orfs");
		addORF(FS);
	}



	public void addORF(FastaSequences FS){
		int pointer = 0;
		int step  = 10000;
		int itemp = 0;
			
		for(int j = 0; j < FS.size(); j++){
			String orfName =  FS.get(j).getName();
			
			
			if(orfName.indexOf(" ") > -1){
				orfName = orfName.split(" ")[0];
				if(orfName.indexOf("_") != orfName.lastIndexOf("_")){
					int index =  orfName.indexOf("_", orfName.indexOf("_")+1);
					orfName = orfName.substring(0,index);
				}
			}

			boolean found = false;
//			System.out.println(orfName);
//			System.out.println(orfName.indexOf(" "));

			try{
				if(this.Chromosomes.containsKey(orfName)){
					this.Chromosomes.get(orfName).addgetORFinfo(FS.get(j));
				}
				else{
					this.Chromosomes.put(orfName,new Contig(orfName));
					this.Chromosomes.get(orfName).addgetORFinfo(FS.get(j));
				}
			}catch(Exception E){
				System.out.println("Something wrong with these two sequences?");
				System.out.println(this.Chromosomes.get(itemp).getName());
				System.out.println(FS.get(j).getName());
				System.out.println();
			}

			if(!found){
//				System.out.println();
//				System.out.println(FS.get(j).getName()+" not found adding it to the system");


			}
			if(j == step){
				System.out.println(step+" orfs processed");
				step = step + 10000;
			}

		}

	}











	private void printORFNames(String dir, String file) {
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(dir+"/"+file));
			for(Contig val : this.Chromosomes.values()){
				
				val.printORFNames(EW);
			}
			EW.flush();
			EW.close();

		}
		catch(Exception E){E.printStackTrace();}

	}

	private void removeNonRedundantHits(int cutoff){
		for(int j = 0; j < this.Chromosomes.size(); j++){
			this.Chromosomes.get(j).removeNonRedundantHits(cutoff);
		}
	}

	private void removeAntisenseHits(int cutoff, int surrounding){
		for(int j = 0; j < this.Chromosomes.size(); j++){
			this.Chromosomes.get(j).removeAntisenseHits(cutoff,surrounding);
		}
	}

	private void printChromosomalHits2(String Dir, String file){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(Dir+"/"+file));
			for(int j = 0; j < this.Chromosomes.size(); j++){
				this.Chromosomes.get(j).printHits(EW);
			}
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}
	}


	private int getNrOfHits(){
		int total = 0;
		for(int  i = 0; i < this.Chromosomes.size(); i++){
			total += this.Chromosomes.get(i).getNrOfHits();
		}
		return total;


	}

	public void compareDifferentRuns(Contigs otherRun, ExtendedWriter ER,String exp1, String exp2){
		int nrOfHits = this.getNrOfHits();
		System.out.println("Number of hits in first genome : "+nrOfHits);
		int otherNrOfHits = otherRun.getNrOfHits();
		System.out.println("Number of hits in second genome: "+otherNrOfHits);

		double ratio = (double)nrOfHits/(double)otherNrOfHits;
		System.out.println("The ratio is: "+ratio);
		ER.println("Name,chromosome,ratio,"+exp1+","+exp2+",location,nrOfHits in "+ exp1+" = "+nrOfHits+",nrOfHits in "+ exp2+" = "+otherNrOfHits);

		for(int  i = 0; i < this.Chromosomes.size(); i++){
			this.Chromosomes.get(i).compareDistribution(otherRun.Chromosomes.get(i),ER);
		}


	}



	private void printDistributionOfHits(String dir,String file){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(dir+"/"+file+".distribution"));
			for(int j = 0; j < this.Chromosomes.size(); j++){
				this.Chromosomes.get(j).printDistributionSolidSequence(EW);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}

	private void printCodingGenes(String dir,String file, int upstream, int downStream){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(dir+"/"+file+".fa"));
			for(int j = 0; j < this.Chromosomes.size(); j++){
				this.Chromosomes.get(j).printDistributionSolidSequence(EW);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}


	private void printChromosomalHits(String dir,String file){
		try{

			for(int j = 0; j < this.Chromosomes.size(); j++){
				if(!IOTools.isDir(dir+"/"+this.Chromosomes.get(j).getName()))
					IOTools.mkDir(dir+"/"+this.Chromosomes.get(j).getName());
				this.Chromosomes.get(j).printHits(dir+"/"+this.Chromosomes.get(j).getName(),file);
			}
		}catch(Exception E){E.printStackTrace();}
	}

	private static void printDistributionOfHits(String dir,String file,String[] experiments, int[][] distribution){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(dir+"/"+file+".distribution.csv"));
			EW.println(" ,ncRNAs,repeats,intergenic regions,coding genes");
			for(int j = 0; j < experiments.length; j++){
				EW.print(experiments[j]);
				for(int i = 0 ; i < distribution[j].length; i++)
					EW.print(","+distribution[j][i]);
				EW.println();
			}

			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}






}





