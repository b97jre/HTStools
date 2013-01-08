package alignment;

import general.Functions;
import general.IOTools;

import java.io.FileWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;

import Sequence.CfastaSequences;
import Sequence.Solid;

import general.ExtendedWriter;




public class Genome implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//private String GenomeName;
	private String Name;
	protected ArrayList <Chromosome> Chromosomes;
	protected CfastaSequences solidSequences;

	public static void main(String[] args) {
		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}
		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);

		Genome.run(T);

	}


	public Genome(String Name, String gffDir, String[] chromosomes, String ncRNAfile, String repeatFile){

		this.Name = Name;
		// add chromosomes
		this.Chromosomes = new ArrayList<Chromosome>();
		for(int i = 0; i< chromosomes.length; i++){
			String gffFile = "chromosome_"+chromosomes[i]+".gff";
			Chromosome A = new Chromosome(gffDir,gffFile,chromosomes[i]);
			this.Chromosomes.add(A);
		}
		// 
		for(int i = 0; i< this.Chromosomes.size(); i++){
			//System.out.println(this.Chromosomes.get(i).getName());
			//System.out.println(B.Chromosomes.get(i).Genes.size());
			if(ncRNAfile != null){
				this.Chromosomes.get(i).addGFF3info(gffDir,ncRNAfile);
			}
			if(repeatFile != null){
				this.Chromosomes.get(i).addGFF3info(gffDir, repeatFile);
			}

			//System.out.println(B.Chromosomes.get(i).Genes.size());
			this.Chromosomes.get(i).sortGenes();
			//System.out.println(B.Chromosomes.get(i).Genes.size());
			//System.out.println();
		}


	}


	public void addSolidSequences(String solidDir, String solidFile){
		this.loadSolidSequences(solidDir, solidFile);
		this.mapSolidSequences();
	}

	public void addRmapperSequences(String rmapperDir, String rmapperFile){
		this.loadSolidSequences(rmapperDir, rmapperFile);
		this.mapSolidSequences();
	}





	public static void run(Hashtable<String,String> T){
		if(T.containsKey("-countGroups")){
			countGroups(T);
			return;
		}

		String solidDir = Functions.getValue(T, "-solidDir", ".");
		String solidFile = Functions.getValue(T, "-solidFile", "sequences.cfasta");
		String rmapperDir = Functions.getValue(T, "-rmapperDir", ".");
		String rmapperFile = Functions.getValue(T, "-rmapperFile", "sequences.rmapper");
		String gffDir = Functions.getValue(T, "-gffDir", solidDir);
		String[] chromosomes = null; 
		if(T.containsKey("-chromosomes"))
			chromosomes = T.get("-chromosomes").split(" ");
		String[] experiments = null;

		if(T.containsKey("-experiments"))
			experiments = T.get("-experiments").split(" ");

		String ncRNAfile,repeatFile;
		ncRNAfile = repeatFile = null;
//		ncRNAfile = Functions.getValue(T, "-ncRNAFile", "090915_GClibrary_annotation.gff");
		repeatFile = Functions.getValue(T, "-repeatsFile", "dicty_repeats.gff");

		if(experiments != null){
			int[][] distribution = new int[experiments.length][4];
			for(int i = 0; i < experiments.length; i++){

				String unique20merDir = solidDir+"/"+experiments[i]+"/unique20mersDir";

				Genome B = new Genome("temp", gffDir, chromosomes, ncRNAfile, repeatFile);
				if(T.containsKey("-solidDir")){
					B.loadSolidSequences(solidDir, solidFile);
				}
				if(T.containsKey("-rmapperDir")){
					B.loadRmapperSequences(rmapperDir, rmapperFile);
				}
				if(T.containsKey("-map")){
					B.mapSolidSequences();
				}
				if(T.containsKey("-removeNonRedundant")){
					//glob/solid/results/d0h/21/d0h.
					int startLength = Integer.parseInt(Functions.getValue(T, "-NucleotideLengthStart", "15"));
					int stopLength = Integer.parseInt(Functions.getValue(T, "-NucleotideLengthStop", "35"));
					int nrOfHitsStart = Integer.parseInt(Functions.getValue(T, "-nrOfHitsStart", "1"));
					int nrOfHitsStop = Integer.parseInt(Functions.getValue(T, "-nrOfHitsStop", "10"));
					int cutoff = Integer.parseInt(Functions.getValue(T, "-cutoff", "10"));
					String filterFile = Functions.getValue(T, "-queryFile", "dicty_genome.fa");
					String resultDir = Functions.getValue(T, "-resultDir", ".");
					for(int j = nrOfHitsStart; j< nrOfHitsStop; j++){
						for(int k = startLength; k< stopLength;k++){
							Genome C = new Genome("temp", gffDir, chromosomes, ncRNAfile, repeatFile);
							String expDir = resultDir+"/"+experiments[i]+"/"+k;
							String expFile = experiments[i]+"."+filterFile+".rmapper."+j;
							String extra = Functions.getValue(T, "-extra", "");
							if(extra.length() > 1)
								expFile+="."+extra;
							C.loadRmapperSequences(expDir, expFile);
							C.mapSolidSequencesonChromosome();
							C.removeNonRedundantHits(cutoff);
							C.printRedundantHits(expDir, expFile+".redundant",cutoff);
						}
					}
				}

				if(T.containsKey("-removeAntisense")){
					//glob/solid/results/d0h/21/d0h.
					int startLength = Integer.parseInt(Functions.getValue(T, "-NucleotideLengthStart", "15"));
					int stopLength = Integer.parseInt(Functions.getValue(T, "-NucleotideLengthStop", "35"));
					int nrOfHitsStart = Integer.parseInt(Functions.getValue(T, "-nrOfHitsStart", "1"));
					int nrOfHitsStop = Integer.parseInt(Functions.getValue(T, "-nrOfHitsStop", "10"));
					int cutoff = Integer.parseInt(Functions.getValue(T, "-cutoff", "10"));
					String filterFile = Functions.getValue(T, "-queryFile", "dicty_genome.fa");
					String resultDir = Functions.getValue(T, "-resultDir", ".");
					int surrounding = Integer.parseInt(Functions.getValue(T, "-surrounding", "200"));
					for(int j = nrOfHitsStart; j< nrOfHitsStop; j++){
						for(int k = startLength; k< stopLength;k++){
							Genome C = new Genome("temp", gffDir, chromosomes, ncRNAfile, repeatFile);
							String expDir = resultDir+"/"+experiments[i]+"/"+k;
							String expFile = experiments[i]+"."+filterFile+".rmapper."+j;
							String extra = Functions.getValue(T, "-extra", "");
							if(extra.length() > 1)
								expFile+="."+extra;
							C.loadRmapperSequences(expDir, expFile);
							C.mapSolidSequencesonChromosome();
							C.removeAntisenseHits(cutoff, surrounding);
							C.printChromosomalHits2(expDir, expFile+".noAntisense");
						}
					}
				}

				if(T.containsKey("-printPremiRNAs")){
					//glob/solid/results/d0h/21/d0h.
					int startLength = Integer.parseInt(Functions.getValue(T, "-NucleotideLengthStart", "15"));
					int stopLength = Integer.parseInt(Functions.getValue(T, "-NucleotideLengthStop", "35"));
					int nrOfHitsStart = Integer.parseInt(Functions.getValue(T, "-nrOfHitsStart", "1"));
					int nrOfHitsStop = Integer.parseInt(Functions.getValue(T, "-nrOfHitsStop", "10"));
					int cutoff = Integer.parseInt(Functions.getValue(T, "-cutoff", "10"));
					String filterFile = Functions.getValue(T, "-queryFile", "dicty_genome.fa");
					String resultDir = Functions.getValue(T, "-resultDir", ".");
					int surrounding = Integer.parseInt(Functions.getValue(T, "-surrounding", "200"));
					for(int j = nrOfHitsStart; j< nrOfHitsStop; j++){
						for(int k = startLength; k< stopLength;k++){
							Genome C = new Genome("temp", gffDir, chromosomes, ncRNAfile, repeatFile);
							String expDir = resultDir+"/"+experiments[i]+"/"+k;
							String expFile = experiments[i]+"."+filterFile+".rmapper."+j;
							String extra = Functions.getValue(T, "-extra", "");
							if(extra.length() > 1)
								expFile+="."+extra;
							C.loadRmapperSequences(expDir, expFile);
							C.mapSolidSequences();
							C.printDistributionOfHits(expDir,expFile);
							//						C.printHits(expDir,expFile);
							C.printSurrounding(expDir,expFile);
						}
					}
				}


				if(T.containsKey("-surrounding")){
					B.printSolidSequence(unique20merDir,solidFile+".surrounding.fa");
				}

				if(T.containsKey("-printGroups")){
					B.solidSequences.printGroups(unique20merDir,solidFile);
				}

				if(T.containsKey("-countNormGroups")){
					B.solidSequences.countNormalizedGroups();
				}




				if(T.containsKey("-distribution") && T.containsKey("-map")){
					distribution[i] = B.DistributionOfHits(distribution[i]);
				}

				if(T.containsKey("-maxNrOfHits")){
					for(int j = 1; j < Integer.parseInt(Functions.getValue(T, "-maxNrOfHits")); j++){
						CfastaSequences temp = new CfastaSequences();
						temp = CfastaSequences.getNrOfHits(B.solidSequences,j);
						temp.printNrOfHits(unique20merDir, solidFile,Integer.parseInt(Functions.getValue(T, "-maxNrOfHits")));
					}
				}

				if(T.containsKey("-UpToNrOfHits")){
					int nrOfHits = Integer.parseInt(Functions.getValue(T, "-UpToNrOfHits"));
					CfastaSequences temp = new CfastaSequences();
					temp = CfastaSequences.getFewerHits(B.solidSequences,nrOfHits);
					temp.printNrOfHits(unique20merDir, solidFile,nrOfHits);
				}



			}
			printDistributionOfHits(solidDir,"distributionOfExperiments.txt", experiments,distribution);
		}
		else{
			Genome B = new Genome("temp", gffDir, chromosomes, ncRNAfile, repeatFile);
			B.print3UTRs(gffDir,"3UTRSequences", 200);
			B.print3UTRsTab(gffDir,"3UTRSequences", 200);


		}
	}





	Genome(){
		this.Chromosomes = new ArrayList<Chromosome>();
	}

	public static void countGroups(Hashtable<String,String> T ){	
		String[] experiments = null; 
		experiments = T.get("-experiments").split(" ");

		int startLength = Integer.parseInt(Functions.getValue(T, "-NucleotideLengthStart", "15"));
		int stopLength = Integer.parseInt(Functions.getValue(T, "-NucleotideLengthStop", "35"));
		int nrOfHitsStart = Integer.parseInt(Functions.getValue(T, "-nrOfHitsStart", "1"));
		int nrOfHitsStop = Integer.parseInt(Functions.getValue(T, "-nrOfHitsStop", "1"));
		if(T.containsKey("-experiments"))
			experiments = T.get("-experiments").split(" ");
		int[][][][] counts = new int[experiments.length][][][];
		for(int i = 0; i < experiments.length; i++){
			System.out.println(experiments[i]);
			int[][][] expCount = new int[stopLength-startLength][][];
			for(int j = startLength; j < stopLength; j++){
				System.out.println("length:\t"+(j));
				System.out.print(" \t");
				for(int k = nrOfHitsStart; k < nrOfHitsStop; k++){
					System.out.print(k+"\t");
				}
				System.out.println("Total");
				int[][] readsCount = new int[nrOfHitsStop-nrOfHitsStart][]; 
				for(int k = nrOfHitsStart; k < nrOfHitsStop; k++){
					readsCount[k-nrOfHitsStart] = countGroups(experiments[i],j,k,T);
				}

				System.out.print("DIRS-1\t");
				for(int k = 0; k < nrOfHitsStop-nrOfHitsStart; k++){
					System.out.print(readsCount[k][6]+"\t");
				}
				System.out.println();

				System.out.print("Complex repeats\t");
				for(int k = 0; k < nrOfHitsStop-nrOfHitsStart; k++){
					System.out.print(readsCount[k][5]+"\t");
				}
				System.out.println();

				System.out.print("ncRNA\t");
				for(int k = 0; k < nrOfHitsStop-nrOfHitsStart; k++){
					System.out.print(readsCount[k][2]+"\t");
				}
				System.out.println();

				System.out.print("intergenic\t");
				for(int k = 0; k < nrOfHitsStop-nrOfHitsStart; k++){
					System.out.print(readsCount[k][3]+"\t");
				}
				System.out.println();

				System.out.print("mRNA\t");
				for(int k = 0; k < nrOfHitsStop-nrOfHitsStart; k++){
					System.out.print(readsCount[k][1]+"\t");
				}
				System.out.println();

				System.out.print("antisenseRNA\t");
				for(int k = 0; k < nrOfHitsStop-nrOfHitsStart; k++){
					System.out.print(readsCount[k][4]+"\t");
				}
				System.out.println();


				expCount[j-startLength] = readsCount;
			}
			counts[i]=expCount;
		}

		for(int i = 0; i < experiments.length; i++){
			System.out.println(experiments[i]);
			for(int j = 0; j < stopLength-startLength; j++){
				System.out.print(" \t");
				for(int k = 1; k < 31; k++){
					System.out.print(k+"\t");
				}
				System.out.println("Total");

				System.out.print("DIRS-1\t");
				for(int k = 1; k < 31; k++){
					System.out.print(counts[i][j][k][6]+"\t");
				}
				System.out.println();

				System.out.print("Complex repeats\t");
				for(int k = 1; k < 31; k++){
					System.out.print(counts[i][j][k][5]+"\t");
				}
				System.out.println();

				System.out.print("ncRNA\t");
				for(int k = 1; k < 31; k++){
					System.out.print(counts[i][j][k][2]+"\t");
				}
				System.out.println();

				System.out.print("intergenic\t");
				for(int k = 1; k < 31; k++){
					System.out.print(counts[i][j][k][3]+"\t");
				}
				System.out.println();

				System.out.print("mRNA\t");
				for(int k = 1; k < 31; k++){
					System.out.print(counts[i][j][k][1]+"\t");
				}
				System.out.println();

				System.out.print("antisenseRNA\t");
				for(int k = 1; k < 31; k++){
					System.out.print(counts[i][j][k][4]+"\t");
				}
				System.out.println();

			}
		}
		/*	1=mRNA
			2=ncRNA
			3=intergenic
			4=antisense
			5=repeats
			6=DIRS
		 */



	}


	private void loadSolidSequences(String dir, String file){
		this.solidSequences = new CfastaSequences();
		this.solidSequences.addSolidSequences(dir,file);
	}

	private void loadRmapperSequences(String dir, String file){
		this.solidSequences = new CfastaSequences();
		this.solidSequences.addRmapperSequences(dir,file);
		this.solidSequences.countHits();
	}

	private void removeNonredundantRmapperSequences(String dir, String file){
		this.solidSequences = new CfastaSequences();
		this.solidSequences.addRmapperSequences(dir,file);




	}


	public static int[] countGroups(String experiment, int SequenceLength, int nrOfHits, Hashtable<String,String> T){

		String resultDir = Functions.getValue(T, "-resultDir", ".");
		String gffDir = Functions.getValue(T, "-gffDir", resultDir);
		String filterFile = Functions.getValue(T, "-queryFile", "dicty_genome.fa");

		String[] chromosomes = null; 
		if(T.containsKey("-chromosomes"))
			chromosomes = T.get("-chromosomes").split(" ");
		String ncRNAfile,repeatFile;
		ncRNAfile = repeatFile = null;
		//ncRNAfile = Functions.getValue(T, "-ncRNAFile", "090915_GClibrary_annotation.gff");
		repeatFile = Functions.getValue(T, "-repeatsFile", "dicty_repeats.gff");
		String finalDir = resultDir+"/"+experiment+"/"+SequenceLength+"/";

		String finalFile = experiment+"."+filterFile+".rmapper."+nrOfHits; 
		String extra = Functions.getValue(T, "-extra", "");
		if(extra.length() > 1)
			finalFile+="."+extra;

		Genome B = new Genome("temp", gffDir, chromosomes, ncRNAfile, repeatFile);
		B.loadRmapperSequences(finalDir, finalFile);
		B.mapSolidSequences();
		return B.solidSequences.countGroups();
	}




	private void mapSolidSequences(){
		for(int  i = 0 ; i < this.solidSequences.size(); i++){
			Solid hit = this.solidSequences.get(i);
			for(int j = 0; j < this.Chromosomes.size(); j++){
				this.Chromosomes.get(j).mapSolidSequence(hit);
			}
		}
	}

	private void mapSolidSequencesonChromosome(){
		for(int  i = 0 ; i < this.solidSequences.size(); i++){
			Solid hit = this.solidSequences.get(i);
			for(int j = 0; j < this.Chromosomes.size(); j++){
				this.Chromosomes.get(j).mapSolidSequenceChromosome(hit);
			}
		}
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

	public void compareDifferentRuns(Genome otherRun, ExtendedWriter ER,String exp1, String exp2){
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



	private void printSolidSequence(String dir, String fileName){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(dir+"/"+fileName));
			for(int  i = 0 ; i < this.solidSequences.size(); i++){
				Solid hit = this.solidSequences.get(i);
				for(int j = 0; j < this.Chromosomes.size(); j++){
					this.Chromosomes.get(j).printSolidSequence(hit,EW);
				}
			}
			EW.flush(); 
			EW.close();
		}catch(Exception E){E.printStackTrace();}
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


	private void print3UTRs(String dir,String file, int downStream){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(dir+"/"+file+".fa"));
			for(int j = 0; j < this.Chromosomes.size(); j++){
				this.Chromosomes.get(j).print3UTRs(EW, downStream);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}

	private void print3UTRsTab(String dir,String file, int downStream){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(dir+"/"+file+".tab"));
			for(int j = 0; j < this.Chromosomes.size(); j++){
				this.Chromosomes.get(j).print3UTRsTab(EW, downStream);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}

	private void printHits(String dir,String file){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(dir+"/"+file));

			this.solidSequences.printRmapperHits(EW);
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}

	private void printRedundantHits(String dir,String file, int cutoff){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(dir+"/"+file));
			this.solidSequences.printRedundantRmapperHits(EW,cutoff);
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


	private void printSurrounding(String dir,String file){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(dir+"/"+file+".surrounding"));
			for(int j = 0; j < this.Chromosomes.size(); j++){
				this.Chromosomes.get(j).printSurrounding(EW,200);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(dir+"/"+file+".fa"));
			for(int j = 0; j < this.Chromosomes.size(); j++){
				this.Chromosomes.get(j).printSurrounding(EW,0);
			}
			EW.flush();
			EW.close();
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


	private int[] DistributionOfHits(int[] distribution){
		try{
			for(int j = 0; j < this.Chromosomes.size(); j++){
				distribution = this.Chromosomes.get(j).DistributionSolidSequence(distribution);
			}
		}catch(Exception E){E.printStackTrace();}
		return distribution;
	}



	public void readGenomeFile(String dir, String[] chromosomeNames){

		for(int i = 0; i< chromosomeNames.length; i++){
			String filename = "chromosome_"+chromosomeNames[i]+".gff";
			this.Chromosomes.add(new Chromosome(dir, filename,chromosomeNames[i]));
		}






	}





}




