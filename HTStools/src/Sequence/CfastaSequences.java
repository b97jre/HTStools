package Sequence;

import general.Functions;
import general.IOTools;
import general.RNAfunctions;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;

import alignment.Hit;

import sun.tools.tree.ThisExpression;
import general.ExtendedReader;
import general.ExtendedWriter;




public class CfastaSequences extends ArrayList <Solid> implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//private String GenomeName;
	private String Name;
	
	

	public CfastaSequences(){
		this.Name = "temp";
	}
	
	public CfastaSequences(String Name){
		this.Name = Name;
	}
	
	
	public CfastaSequences(String dir, String file){
		addSolidSequences(dir,file);
	}
	
	
	public static void main(String[] args) {
			int length = args.length;
			for (int i = 0; i < length; i++){
				args[i] = args[i].trim();
				System.out.print(args[i]+" ");
			}
			System.out.println();
			Hashtable<String,String> T = Functions.parseCommandLine(args);
			if(T.containsKey("-initial"))
				runEntireFiltering(T);
			else
				run(T);
	}
	

	public static void run(Hashtable<String,String> T){
		if(T.containsKey("-initial"))runEntireFiltering(T);
		else{
			String[] experiments = null;
			if(T.containsKey("-experiments"))
				experiments = T.get("-experiments").split(" ");

			if(experiments == null){
				runExperiments(T,"");
			}
			else{ 
				for(int i = 0; i < experiments.length; i++){
					runExperiments(T,experiments[i]);
				}
			}
			
			if(T.containsKey("-primer")){
				int[][] distributions = null;
				String solidDir = Functions.getValue(T, "-solidDir", ".");
				String solidFile = Functions.getValue(T, "-solidFile", "test.cfa");
				String specificDir = Functions.getValue(T, "-specificDir", "");

				if(experiments == null){
					distributions = new int[1][];
					CfastaSequences C1 = new CfastaSequences();
					distributions[0] = C1.getSizeDistribution(solidDir,solidFile,
							Integer.parseInt(Functions.getValue(T, "-max", "500000")),
							Functions.getValue(T, "-primer", "GCGGAACCGGCATGTCGTC").split(" "));
					
				}
				else{
					distributions = new int[experiments.length][];
					
					for(int i = 0; i < experiments.length; i++){
						CfastaSequences C1 = new CfastaSequences();
						distributions[i] = C1.getSizeDistribution(solidDir+"/"+experiments[i]+"/"+specificDir,solidFile,
								Integer.parseInt(Functions.getValue(T, "-max", "500000")),
								Functions.getValue(T, "-primer", "GCGGAACCGGCATGTCGTC").split(" "));
						//C1.printSequences(solidDir,solidFile+".3adapter");
					}
				}
				try{
					if(!IOTools.isDir(solidDir+"/"+specificDir+"/"))
						IOTools.mkDir(solidDir+"/"+specificDir+"/");
					ExtendedWriter EW = new ExtendedWriter(new FileWriter(solidDir+"/"+specificDir+"/"+solidFile+".sizeDistribution.csv"));
					EW.print(" ");
					for(int i = 0; i < experiments.length; i++){
						EW.print(","+experiments[i]);
					}
					EW.println();
					for(int j = 0; j < distributions[0].length;j++){
						EW.print(j);
						for(int i = 0; i < experiments.length; i++){
						EW.print(","+distributions[i][j]);
						}
						EW.println();
						
					}
					EW.flush();
					EW.close();
				}
				catch(Exception E){E.printStackTrace();}
			}
			
		}

		
		if(T.containsKey("-RmapperFile")){
			CfastaSequences A = new CfastaSequences();
			String solidDir = Functions.getValue(T, "-Rdir", ".");
			String solidFile = Functions.getValue(T, "-RmapperFile", "test.cfa");
			if(T.containsKey("-split"))
				A.splitRmapperAdapterSequences(solidDir,solidFile,35);
			if(T.containsKey("-count"))
				A.countRmapperAdapterSequences(solidDir,solidFile,35);
			if(T.containsKey("-count"))
				A.countRmapperAdapterSequences(solidDir,solidFile,35);
			
			else if(T.containsKey("-convert"))
				A.convertRmapperSequences(solidDir,solidFile);
			else if(T.containsKey("-convert2Fasta"))
				A.convertRmapper2Fasta(solidDir,solidFile);
			else if(T.containsKey("-best")){
				String filterFile = Functions.getValue(T, "-queryFile", "dicty_genome.fa");
				printBestRmapperHits(solidDir,solidFile,10000,filterFile);
			}
			else if(T.containsKey("-ordered"))
				if(T.containsKey("-experiments")){
					String[] experiments = T.get("-experiments").split(" ");
					String filterFile = Functions.getValue(T, "-queryFile", "dicty_genome.fa");
				
					int start = 15; 
					int stop = 35;
					int[][][] distribution = new int[experiments.length][stop-start][];
					for(int i = 0; i < experiments.length; i++){
						String tempDir = solidDir +"/"+experiments[i];
						for(int j = start; j < stop; j++){
							String specificDir = tempDir+"/"+j+"/";
							distribution[i][j-start] = printOrderedRmapperHits(specificDir,experiments[i]+"."+filterFile+".rmapper");
						}
					}
					
					for(int i = 0; i < experiments.length; i++){
						System.out.println("Experiment: "+ experiments[i]+"(Length of sequence on rows and nr of Hits per sequence columns)");
						System.out.print("  \t");
						for(int k = 0; k < distribution[i][0].length-1;k++){
							System.out.print((k+1)+"\t");
						}
						System.out.println(">"+distribution[i][0].length);
						for(int j = start; j < stop; j++){
							System.out.print(j);
							for(int k = 0; k < distribution[i][j-start].length;k++){
								System.out.print("\t"+ distribution[i][j-start][k]);
							}
							System.out.println();
						}
						System.out.println();
					}
				}
				else
					printOrderedRmapperHits(solidDir,solidFile);
		}
		
		if(T.containsKey("-cFastaFile")){
			CfastaSequences A = new CfastaSequences();
			String cFastaFileName = Functions.getValue(T, "-cFastaFile", "dicty_genome.fa");
			String cFastaDir = Functions.getValue(T, "-cFastaDir", ".");
			String dir = Functions.getValue(T, "-dir", ".");
			String cFastaOut = Functions.getValue(T, "-cFastaOut", ".");
			String filterFile = Functions.getValue(T, "-queryFile", "dicty_genome.fa");
			if(T.containsKey("-split"))
				A.splitBySizeCfastaSequences(dir,cFastaDir,cFastaFileName,cFastaOut,10,35,filterFile);
		}
		
		
		if(T.containsKey("-filter")){
			String solidDir = Functions.getValue(T, "-dir", ".");
			String solidFile = Functions.getValue(T, "-SF", "test.cfa");
			String filterFile = Functions.getValue(T, "-FF", "test.cfa");
			String outFile = Functions.getValue(T, "-OF", "test.cfa");
			
			
			FilterDictySequences(solidDir,solidFile,
					solidDir,filterFile,
					solidDir,outFile,100000);
				
		}
		
		
		if(T.containsKey("-split2")){
			int[][] distributions = null;
			String solidDir = Functions.getValue(T, "-solidDir", ".");
			String solidFile = Functions.getValue(T, "-solidFile", "test.cfa");
			String specificDir = Functions.getValue(T, "-specificDir", "");
			String[] experiments = null;
			if(T.containsKey("-experiments"))
				experiments = T.get("-experiments").split(" ");

			if(experiments == null){
				distributions = new int[1][];
				CfastaSequences C1 = new CfastaSequences();
				distributions[0] = C1.getSizeDistribution(solidDir,solidFile,
						Integer.parseInt(Functions.getValue(T, "-max", "500000")),
						Functions.getValue(T, "-primer", "GCGGAACCGGCATGTCGTC").split(" "));
				
			}
			else{
				distributions = new int[experiments.length][];
				
				for(int i = 0; i < experiments.length; i++){
					CfastaSequences C1 = new CfastaSequences();
					distributions[i] = C1.getSizeDistribution(solidDir+"/"+experiments[i]+"/"+specificDir,solidFile,
							Integer.parseInt(Functions.getValue(T, "-max", "500000")),
							Functions.getValue(T, "-primer", "GCGGAACCGGCATGTCGTC").split(" "));
					//C1.printSequences(solidDir,solidFile+".3adapter");
				}
			}
			try{
				if(!IOTools.isDir(solidDir+"/"+specificDir+"/"))
					IOTools.mkDir(solidDir+"/"+specificDir+"/");
				ExtendedWriter EW = new ExtendedWriter(new FileWriter(solidDir+"/"+specificDir+"/"+solidFile+".sizeDistribution.csv"));
				EW.print(" ");
				for(int i = 0; i < experiments.length; i++){
					EW.print(","+experiments[i]);
				}
				EW.println();
				for(int j = 0; j < distributions[0].length;j++){
					EW.print(j);
					for(int i = 0; i < experiments.length; i++){
					EW.print(","+distributions[i][j]);
					}
					EW.println();
					
				}
				EW.flush();
				EW.close();
			}
			catch(Exception E){E.printStackTrace();}
		}
		
		else if(T.containsKey("-filter2")){
			String dir = Functions.getValue(T, "-d", ".");
			String suffix = Functions.getValue(T, "-suffix", "gmapper");
			parseGmapperSequences(dir,suffix);
		}

	}
		
	
	
	
	public void mapSequences(FastaSequences FS){
		for(int i = 0; i < this.size();i++)
			this.get(i).mapHits(FS);
	}
	
	public void mapSequences(FastaSequences FS, int exp){
		for(int i = 0; i < this.size();i++)
			this.get(i).mapHits(FS, exp);
	}
	
	public int getNrOfHits(){
		int total = 0;
		for(int i = 0; i < this.size();i++)
			total += this.get(i).countHits();
		return total;
	}
	
		
	
	public static void runEntireFiltering(Hashtable<String,String> T){
		String solidDir = Functions.getValue(T, "-solidDir", ".");
		String genomeDir = Functions.getValue(T, "-genomeDir", solidDir+"/genome");
		String mirBaseDir = Functions.getValue(T, "-mirBaseDir", solidDir+"/miRBase");
		String longChromosomeFile = Functions.getValue(T, "-longChromosomeFile", "ugc_SREK_Sel_CW_solid0105_20090427_ugc_65__71_F3.csfasta_extend.ma.50.2");
		String shortChromosomeFile = Functions.getValue(T, "-shortChromosomeFile", "ugc_SREK_Sel_CW_solid0105_20090427_ugc_65__71_F3.csfasta.ma.50.2");
		String longmiRBaseFile = Functions.getValue(T, "-longChromosomeFile", "ugc_SREK_Sel_CW_solid0105_20090427_ugc_65__71_F3.csfasta_extend.ma.50.4");
		String shortmiRBaseFile = Functions.getValue(T, "-shortChromosomeFile", "ugc_SREK_Sel_CW_solid0105_20090427_ugc_65__71_F3.csfasta.ma.50.2");
		
		String unique20mersFile = Functions.getValue(T, "-unique20mersFile", "unique20merSequences");
		String unique20mersDir = Functions.getValue(T, "-unique20mersDir", solidDir+"/unique20mersDir");
		
		System.out.println("running cfasta");
		String[] chromosomes = T.get("-chromosomes").split(" ");
		
		// aligning long sequences on the chromosomes
		joinLong(genomeDir,chromosomes,longChromosomeFile);
		
		// removing sequenses in short miRBase step that are aligned in longer locations in the chromosome or in the miRBASE library  
		longHitsFilter(mirBaseDir,shortmiRBaseFile, longmiRBaseFile, unique20mersDir, unique20mersFile, genomeDir, longChromosomeFile);
		
		//removing sequences from the short miRBase sequences that exist as longer mirDatabase sequences)
		ncRNARepeatsFilter(unique20mersDir, unique20mersFile, chromosomes);
		
		// join all sequences from the joined long library and from the short filtered sequences of the different libraries
		joinSequences(unique20mersDir, unique20mersFile, genomeDir, longChromosomeFile, chromosomes);
			
		// printing single hits
		printSingleHits(unique20mersDir,unique20mersFile+".extended_filter.ncRNA_repeats_filter.joined");
		
		//printing distribution
		printDistribution(unique20mersDir,unique20mersFile+".extended_filter.ncRNA_repeats_filter.joined");
		

	}

	private static void printSingleHits(String dir,String file){
		CfastaSequences C1 = new CfastaSequences();
		C1.addSolidSequences(dir,file);
		C1 = singleHits(C1);
		C1.printSingle(dir, file);
	}

	private static void printDistribution(String dir,String file){
		CfastaSequences C1 = new CfastaSequences();
		C1.addSolidSequences(dir,file);
		int[] dist = getHitDistribution(C1);
		printDistribution(dist,0,1);
	}

	
	
	
	
	private static void longHitsFilter(String mirBaseDir,String shortmiRBaseFile,String longmiRBaseFile,
			String unique20mersDir, String unique20mersFile,String genomeDir,String longChromosomeFile){
		//removing sequences from the short miRBase sequences that exist as longer mirDatabase sequences
		FilterDictySequences(mirBaseDir,shortmiRBaseFile,
				mirBaseDir,longmiRBaseFile,
				unique20mersDir,"ncRNA_repeats."+unique20mersFile+".extended_filter");

		//removing sequences from the short miRBase sequences that exist as longer longer genome sequences
		FilterDictySequences(unique20mersDir,"ncRNA_repeats."+unique20mersFile+".extended_filter",
				genomeDir, longChromosomeFile+".joined",
				unique20mersDir,"ncRNA_repeats."+unique20mersFile+".extended_filter");
		// finished
	}

	private static void ncRNARepeatsFilter(String unique20mersDir, String unique20mersFile,String [] chromosomes){
		CfastaSequences filter = new CfastaSequences();
		filter.addSolidSequencesSimplest(unique20mersDir, "ncRNA_repeats."+unique20mersFile+".extended_filter");
		for(int i = 0; i< chromosomes.length; i++){
			FilterRepeatDictySequences(unique20mersDir,chromosomes[i]+"."+unique20mersFile+".extended_filter",
					"ncRNA_repeats."+unique20mersFile+".extended_filter", filter,"ncRNA_repeats_filter");
		}
	}

	private static void joinSequences(String unique20mersDir, String unique20mersFile,String genomeDir,String longChromosomeFile,String[] chromosomes){
		ArrayList <Solid> sequences = new ArrayList<Solid>();
		for(int i = 0; i< chromosomes.length; i++){
			CfastaSequences C1 = new CfastaSequences();
			C1.addSolidSequences(unique20mersDir, chromosomes[i]+"."+unique20mersFile+".extended_filter.ncRNA_repeats_filter",chromosomes[i]);
			sequences = joinSequences(sequences,C1);
		}
		CfastaSequences C1 = new CfastaSequences();
		C1.addSolidSequences(genomeDir,longChromosomeFile+".joined");
		sequences = joinSequences(sequences,C1);

		CfastaSequences joined = new CfastaSequences();
		joined.addAll(sequences);
		joined.printJoinedInfo(unique20mersDir, unique20mersFile+".extended_filter.ncRNA_repeats_filter", chromosomes);
	}


	private static void joinLong(String genomeDir,String[] chromosomes, String longChromosomeFile){
		ArrayList <Solid> sequences = new ArrayList<Solid>();
		for(int i = 0; i< chromosomes.length; i++){
			CfastaSequences C1 = new CfastaSequences();
			C1.addSolidSequences(genomeDir+"/"+chromosomes[i],longChromosomeFile,chromosomes[i]);
			sequences = joinSequences(sequences,C1);
		}
		CfastaSequences joined = new CfastaSequences();
		joined.addAll(sequences);
		joined.printJoinedInfo(genomeDir, longChromosomeFile, chromosomes);
	}
	


	public static void runExperiments(Hashtable<String,String> T, String experiment){	
		String solidDir = Functions.getValue(T, "-solidDir", ".");
		String solidFile = Functions.getValue(T, "-solidFile", "test.cfa");
		String longChromosomeFile = Functions.getValue(T, "-longChromosomeFile", "ugc_SREK_Sel_CW_solid0105_20090427_ugc_65__71_F3.csfasta_extend.ma.50.2");
		String shortChromosomeFile = Functions.getValue(T, "-shortChromosomeFile", "ugc_SREK_Sel_CW_solid0105_20090427_ugc_65__71_F3.csfasta.ma.50.2");
		String longmiRBaseFile = Functions.getValue(T, "-longmiRBaseFile", "ugc_SREK_Sel_CW_solid0105_20090427_ugc_65__71_F3.csfasta_extend.ma.50.4");
		String shortmiRBaseFile = Functions.getValue(T, "-shortmiRBaseFile", "ugc_SREK_Sel_CW_solid0105_20090427_ugc_65__71_F3.csfasta.ma.50.2");
		String unique20mersFile = Functions.getValue(T, "-unique20mersFile", "unique20merSequences");



		String genomeDir = Functions.getValue(T, "-genomeDir", solidDir+"/"+experiment+"/genome");
		String mirBaseDir = Functions.getValue(T, "-mirBaseDir", solidDir+"/"+experiment+"/miRBase");
		String unique20mersDir = Functions.getValue(T, "-unique20mersDir", solidDir+"/"+experiment+"/unique20mersDir");

		System.out.println("running cfasta");
		String[] chromosomes = null;
		if(T.containsKey("-chromosomes"))
			chromosomes = T.get("-chromosomes").split(" ");

		System.out.println("running cfasta");


		// aligning long sequences on the chromosomes
		if(T.containsKey("-joinLong"))
			joinLong(genomeDir,chromosomes,longChromosomeFile);

		// removing sequenses in short miRBase step that are aligned in longer locations in the chromosome or in the miRBASE library  
		if(T.containsKey("-longHitsFilter"))
			longHitsFilter(mirBaseDir,shortmiRBaseFile, longmiRBaseFile, unique20mersDir, unique20mersFile, genomeDir, longChromosomeFile);

		//removing sequences from the short miRBase sequences that exist as longer mirDatabase sequences)
		if(T.containsKey("-ncRNARepeatsFilter"))
			ncRNARepeatsFilter(unique20mersDir, unique20mersFile, chromosomes);

		// join all sequences from the joined long library and from the short filtered sequences of the different libraries
		if(T.containsKey("-joinSequences"))
			joinSequences(unique20mersDir, unique20mersFile, genomeDir, longChromosomeFile, chromosomes);

		// printing single hits
		if(T.containsKey("-printSingleHits"))
			printSingleHits(unique20mersDir,unique20mersFile+".extended_filter.ncRNA_repeats_filter.joined");

		//printing distribution
		if(T.containsKey("-printDistribution"))
			printDistribution(unique20mersDir,unique20mersFile+".extended_filter.ncRNA_repeats_filter.joined");






		System.out.println("finished");

	}



	
	public void findPrimers(String primerSequence){
		int [] primer = Solid.setColorCode(primerSequence);
		int [] counter = new int[50];
		int start = 0;
		int count = 0;
		int maxLength = 15;
		int minMatches = 3;
		for(int i = 0; i < this.size(); i++){
			int location = this.get(i).findPrimer(2, start, maxLength, minMatches, primer);
			if(location>0){
				count++;
				counter[location]++;
			}
		}
		
		for(int i = start; i < counter.length; i++){
			System.out.println(i + "\t"+ counter[i]) ;
		}
		
		System.out.println(count);
		System.out.println(this.size());
	}

	
	public int[] findPrimers(String[] primerSequence, int[] counter,CfastaSequences notFound, int minLength){
		int [][] primers = RNAfunctions.fasta2CFasta(primerSequence);
		int start = 0;
		int count = 0;
		int maxLength = 15;
		int minMatches = 3;
		System.out.println("Sequences before :"+this.size());
		for(int i = 0; i < this.size(); i++){
			boolean found = false; 
			int j = 0;
			while(!found && j < primers.length){
				int location = this.get(i).findPrimer(2, start, maxLength, minMatches, primers[j]);
				if(location < minLength){
					counter[location]++;
					found = true;
					this.get(i).removePrimer(location);
				}
				j++;
			}
			if(!found){
				notFound.add(this.get(i));
				this.remove(i);
				i--;
				count++;
			}	
		}
		System.out.println("Sequences after +:"+this.size());
		System.out.println("Sequences removed +:"+count);
		return counter;
		
	}
	
	
	
	private static void FilterDictySequences(String solidDir, String solidFile, String filterDir, String filterFile, String outDir, String outFile){
		CfastaSequences C1 = new CfastaSequences();
		C1.addSolidSequencesSimple(solidDir, solidFile);
		CfastaSequences C2 = new CfastaSequences();
		C2.addSolidSequencesSimple(filterDir, filterFile);

		ArrayList <Solid> uniqueC1 = C1.getUniqueSequences(C2);

		CfastaSequences filtered = new CfastaSequences();
		filtered.setSolidSequences(uniqueC1);

		filtered.printInfoEasy(solidDir,solidFile,filterDir,filterFile,
				C1.size(),C2.size(),filtered.size()
				,outDir,outFile);
	}

	

	private static void FilterDictySequences(String solidDir, String solidFile, String filterDir, String filterFile, String outDir, String outFile, int max){
		try{
		if(!IOTools.isDir(outDir)) 
			IOTools.mkDir(outDir);
		ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile));
		boolean more = true;
		int start = 0;
		while(more){
			CfastaSequences C1 = new CfastaSequences();
			more = C1.addSolidSequencesSimple(solidDir, solidFile,start, start+max);
			int filterStart = 0;
			boolean moreFilter = true;
			
			while(moreFilter){
				CfastaSequences C2 = new CfastaSequences();
				moreFilter = C2.addSolidSequencesSimple(filterDir, filterFile,filterStart, filterStart+max);
				C1.getUniqueSequences2(C2);
				filterStart = filterStart+max;
			}
			for(int i = 0; i < C1.size(); i++){
				C1.get(i).printInfoEasy(EW);
			}
			start =  start+max;
			System.out.println();
			System.out.println("number of sequences read :"+start);
			System.out.println();
			EW.flush();
		}
		EW.close();
		}
		catch(Exception E){E.printStackTrace();}

	}

	
	
	
	
	private static void FilterRepeatDictySequences(String Dir, String solidFile,String filterFile ,CfastaSequences C2, String extension){
		
		CfastaSequences C1 = new CfastaSequences();
		C1.addSolidSequencesSimple(Dir, solidFile);

		ArrayList <Solid> uniqueC1 = C1.getUniqueSequences(C2);

		CfastaSequences filtered = new CfastaSequences();
		filtered.setSolidSequences(uniqueC1);

		filtered.printInfoEasy(Dir,solidFile,Dir,filterFile,
				C1.size(),C2.size(),filtered.size()
				,Dir,solidFile+"."+extension);
	}
	
	
	
		
	
	private void setSolidSequences(ArrayList<Solid> sequences){
		while(this.size() > 0)
			this.remove(1);
		this.addAll(sequences);
	}
	
	private boolean addSolidSequencesSimple(String dir, String fileName, int start, int stop){
		
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));
			int count = 0;
			while(ER.more() && count < stop){
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					if(count >= start && count < stop)
						addcFastaEasy(ER);
					else{
						ER.skipLine();
						ER.skipLine();
					
					}
					count++;
				}
			}
			ER.close();
			if(count < stop)
				return false;
			return true;
		}
		catch(Exception E){E.printStackTrace();}
		return false;
		
	}
	

	
	public static boolean removePrimersDir(String inDir, String outDir, CfastaSequences primers, String suffix,int lowCutoff, int highCutoff){
		
		if(!IOTools.isDir(outDir))
			IOTools.mkDir(outDir);

		ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir,suffix);

		if(fileNames.isEmpty()){
			System.out.println("No "+ suffix+"  files in folder :"+inDir);
			return false;
		}
		else{
			if(!IOTools.isDir(outDir+"/reports"))
				IOTools.mkDir(outDir+"/reports");
			if(!IOTools.isDir(outDir+"/scripts"))
				IOTools.mkDir(outDir+"/scripts");
			try{
			ExtendedWriter EW3 = new ExtendedWriter(new FileWriter(outDir+"/info.txt"));

			EW3.println("Sequences above length "+lowCutoff+" and below length "+highCutoff+" is put i folder short");
			EW3.println("Sequences above length "+highCutoff+" is put i folder long");
			EW3.println();
			EW3.println("Distribution of sequences are:");

			EW3.print("cutoff SequenceFile ");
			for(int i = 0; i< 35; i++){
				EW3.print(i+" ");
			}
			EW3.println();
			
			
			for(int i = 0; i < fileNames.size(); i++){	
				System.out.println("removing primers in file "+fileNames.get(i));
				removePrimers(EW3,inDir, outDir, fileNames.get(i),primers,lowCutoff,highCutoff);
			}
			EW3.flush();
			EW3.close();
			}
			catch(Exception E){E.printStackTrace();}

			
		}
		return true;
	}
	
	public static boolean removePrimers(ExtendedWriter EW3, String inDir, String outDir, String SequenceFile,CfastaSequences primers,int lowCutoff, int highCutoff){
		
		try{
			//System.out.println();
			double cutoff = 0.85;
//			for(double cutoff = 0.0; cutoff < 1.0; cutoff = cutoff+0.1){
				ExtendedReader ER = new ExtendedReader(new FileReader (inDir+"/"+SequenceFile));
				if(!IOTools.isDir(outDir)) 	IOTools.mkDir(outDir);
				if(!IOTools.isDir(outDir+"/short")) IOTools.mkDir(outDir+"/short");
				if(!IOTools.isDir(outDir+"/long")) 	IOTools.mkDir(outDir+"/long");
				if(!IOTools.isDir(outDir+"/other")) 	IOTools.mkDir(outDir+"/other");
				
				ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/short/"+SequenceFile));
				ExtendedWriter EW2 = new ExtendedWriter(new FileWriter(outDir+"/long/"+SequenceFile));
				ExtendedWriter EW4 = new ExtendedWriter(new FileWriter(outDir+"/other/"+SequenceFile));
				
				int count = 0;
				int[] dist = new int[36];
				while(ER.more()){
					if((char)ER.lookAhead() == '#'){
						ER.skipLine();
					}
					else{
						String infoLine = ER.readLine();
						String SequenceLine = ER.readLine();
						Solid Sequence = new Solid(infoLine, SequenceLine);
						dist[Sequence.removePrimer(primers, cutoff)]++;
						if(Sequence.printcFastaSize(EW,17,32));
						else if(Sequence.printcFastaSize(EW2,32));
						else Sequence.printcFasta(EW4);
					}
				}
				int total = 0;
				EW3.print(cutoff+" "+SequenceFile+" ");
				for(int i = 0; i< dist.length; i++){
					EW3.print(dist[i]+" ");
					total = total + dist[i];
				}
				
				EW3.println(total);
			ER.close();
			EW.flush();
			EW.close();
			EW2.flush();
			EW2.close();
//			}
		
		}
		catch(Exception E){E.printStackTrace();}
		return false;
		
	}
	
	
	
	
	private void addSolidSequencesSimple(String dir, String fileName){
		
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));
			int count = 100000;
			while(ER.more()){
				
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					addcFastaEasy(ER);
				}
				if(this.size() > count){
					System.out.println("nr Of Sequences Read" + count);
					count= count + 100000;
				}
			}
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}
	
	public void countHits(){
		int[] NrOfHits = new int[101];
		for(int i = 0; i < this.size(); i++ ){
			NrOfHits[this.get(i).countHits()]++;
		}
		System.out.println("Nr of hits per read distribution:");
		for(int i = 0; i < NrOfHits.length;i++){
			if(NrOfHits[i] > 0)
				System.out.println(i+"\t"+NrOfHits[i]);
		}
	}
	
	public void addSolidSequences(String dir, String fileName){
		
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));
			int count = 100000;
			while(ER.more()){
				
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					addcFasta(ER);
				}
				if(this.size() > count){
					System.out.println("nr Of Sequences Read" + count);
					count= count + 100000;
				}
			}
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}

	public static void parseGmapperSequences(String dir, String suffix){
		ArrayList <String> fileNames = IOTools.getSequenceFiles(dir,suffix);
		for(int i = 0; i < fileNames.size(); i++){
			try{
				ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileNames.get(i)));
			String[] info = fileNames.get(i).split("\\.");
			int count = 100000;
			while(ER.more()){
				
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					extractGmapperInfo(ER, info[0], info[1]);
				}
			}
			ER.close();
		}catch(Exception E){E.printStackTrace();}
		}
	}
	
	
	public void addRmapperSequences(String dir, String fileName){
		
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));
			int count = 100000;
			while(ER.more()){
				
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					addRmapperHit(ER);
				}
				if(this.size() > count){
					System.out.println("nr Of Sequences Read" + count);
					count= count + 100000;
				}
			}
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}
	


	void splitRmapperAdapterSequences(String dir, String fileName, int maxLength){	
		try{
		
		ExtendedWriter[] EW = new ExtendedWriter[maxLength];
		ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));
		
		for(int i = 0; i < maxLength; i++){
			EW[i] = new ExtendedWriter(new FileWriter (dir+"/"+fileName+"."+i+".csfasta"));
		}
		while(ER.more()){
			if((char)ER.lookAhead() == '#'){
				ER.skipLine();
			}
			else
			RmapperRemoveAdapter(ER,EW,maxLength);
		}
		
		ER.close();
		for(int i = 0; i < maxLength; i++){
			EW[i].flush();
			EW[i].close();
		}
		}
			catch(Exception E){E.printStackTrace();}
	
	}

	
	void splitBySizeCfastaSequences(String dir,String fileDir, String fileName,String fileOut, int minLength,  int maxLength, String databaseFile){	
		try{

			ExtendedWriter[] EW = new ExtendedWriter[maxLength];
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/reads/"+fileDir+"/"+fileName));

			for(int i = minLength; i < maxLength; i++){
				if(!IOTools.isDir(dir+"/reads/"+fileDir+"/"+i))IOTools.mkDir(dir+"/reads/"+fileDir+"/"+i);
				EW[i] = new ExtendedWriter(new FileWriter (dir+"/reads/"+fileDir+"/"+i+"/"+fileOut+".csfasta"));
			}
			while(ER.more()){
				if((char)ER.lookAhead() == '>'){
					String name = ER.readLine();
					String Sequence = ER.readLine();
					if(Sequence.length()-1 >= minLength && Sequence.length()-1 < maxLength){
						EW[Sequence.length()-1].println(name);
						EW[Sequence.length()-1].println(Sequence);
					}
				}
			}
			ER.close();
			String cat = new String("cat");
			for(int i = minLength; i < maxLength; i++){
				EW[i].flush();
				EW[i].close();
				splitNrOfReadsCfastaSequences(dir,fileDir,fileOut,"csfasta",10000,i,databaseFile);
				cat += " results/"+fileDir+"/"+i+"/"+fileOut+"."+databaseFile+".rmapper";
			}
			cat += " >results/"+fileDir+"/"+fileOut+"_"+minLength+"_"+(maxLength-1)+"."+databaseFile+".rmapper";
			
			System.out.println(cat);
		}



		catch(Exception E){E.printStackTrace();}

	}


	void splitNrOfReadsCfastaSequences(String dir,String fileDir ,String fileName,String extension, int nrOfReads, int sequenceLength, String DatabaseFile){	
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/reads/"+fileDir+"/"+sequenceLength+"/"+fileName+"."+extension));
			String cat = new  String("cat");
			String rm  = new String("rm");
			int count = 0;
			int start = 0;
			while(ER.more()){
				if(count == 0){
					ExtendedWriter EW = new ExtendedWriter(new FileWriter (dir+"/reads/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+extension));
					count = nrOfReads;
					while (ER.more() && count > 0){
						if((char)ER.lookAhead() == '>'){
							String name = ER.readLine();
							String Sequence = ER.readLine();
							count--;
							EW.println(name);
							EW.println(Sequence);
						}
						else
							ER.skipLine();
					}
					EW.flush();
					EW.close();
/*
 //uppmax
  					if(sequenceLength > 24)
						System.out.println("bin/SHRiMP_1_3_1/bin/rmapper-cs  -s 1111001111 -h 80% -v 85% -n 2 " +
								"glob/solid/reads/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+extension+
								" glob/solid/genome/"+  DatabaseFile+ 
								" >glob/solid/results/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+DatabaseFile+".rmapper");
					else if(sequenceLength > 20)
						System.out.println("bin/SHRiMP_1_3_1/bin/rmapper-cs  -s 1111001111,0011111111,1100111111,1111110011,,1111111100   -h 80% -v 85% -n 2 " +
								"glob/solid/reads/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+extension+
								" glob/solid/genome/"+  DatabaseFile+ 
								" >glob/solid/results/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+DatabaseFile+".rmapper");
					else if(sequenceLength > 18)
						System.out.println("bin/SHRiMP_1_3_1/bin/rmapper-cs  -s 11100111 -h 75% -v 85% -n 2 -R  " +
								"glob/solid/reads/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+extension+
								" glob/solid/genome/"+  DatabaseFile+ 
								" >glob/solid/results/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+DatabaseFile+".rmapper");
					else
						System.out.println("bin/SHRiMP_1_3_1/bin/rmapper-cs  -s 1111111111111 -h 100% -v 100% -n 1 -R  " +
								"glob/solid/reads/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+extension+
								" glob/solid/genome/"+  DatabaseFile+ 
								" >glob/solid/results/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+DatabaseFile+".rmapper");

					cat += " glob/solid/results/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+DatabaseFile+".rmapper.best ";
					rm += " glob/solid/results/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+DatabaseFile+".rmapper " + 
					      " glob/solid/results/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+DatabaseFile+".rmapper.best";
*/

//local
				if(sequenceLength > 18)
					System.out.println("bin/rmapper-cs  -s 11111111111111111111 -h 75% -v 85% -n 2 -R  " +
							"reads/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+extension+
							" genome/"+  DatabaseFile+ 
							" >results/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+DatabaseFile+".rmapper");
				else if(sequenceLength > 13)
					System.out.println("bin/rmapper-cs  -s 1111111111111 -h 100% -v 100% -n 1 -R  " +
							"reads/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+extension+
							" genome/"+  DatabaseFile+ 
							" >results/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+DatabaseFile+".rmapper");

				cat += " results/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+DatabaseFile+".rmapper.best ";
				rm += " results/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+DatabaseFile+".rmapper " + 
				      " results/"+fileDir+"/"+sequenceLength+"/"+fileName+"_"+start+"_"+(start+nrOfReads-1)+"."+DatabaseFile+".rmapper.best";

					
					
					start = start+nrOfReads;

				}
				
				
				
			}
			ER.close();
			if(!IOTools.isDir(dir+"/results/"+fileDir))IOTools.mkDir(dir+"/results/"+fileDir);
			if(!IOTools.isDir(dir+"/results/"+fileDir+"/"+sequenceLength))IOTools.mkDir(dir+"/results/"+fileDir+"/"+sequenceLength);
			System.out.println("java -Xmx500M -jar bin/solid.jar -p CFASTA -Rdir results/"+fileDir+"/"+sequenceLength+
					" -RmapperFile "+fileName + " -best -queryFile "+DatabaseFile);
			cat = cat+" >results/"+fileDir+"/"+sequenceLength+"/"+fileName+"."+DatabaseFile+".rmapper";
			System.out.println(cat);
			System.out.println(rm);
		}
		catch(Exception E){E.printStackTrace();}

	}
	
	
	void countRmapperAdapterSequences(String dir, String fileName, int maxLength){	
		try{

			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(dir+"/"+fileName+".redundantHits"));
			int[] lengths = new int[maxLength];

			while(ER.more()){
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else
					RmapperGetAdapterSizes(ER,lengths,maxLength);
			}

			ER.close();
			for(int i = 0; i < maxLength; i++){
				EW.println(i+"\t"+lengths[i]);
			}
			EW.flush();
			EW.close();		
		}
		catch(Exception E){E.printStackTrace();}

	}

	
	
	void convertRmapperSequences(String dir, String fileName){	
		try{
			String cfastaFileName = null;
			if(fileName.endsWith(".rmapper")){
				cfastaFileName = fileName.substring(0,fileName.length()-".rmapper".length())+".csfasta";
			}
			else
				cfastaFileName = fileName+".csfasta";
				
			ExtendedWriter EW = new ExtendedWriter(new FileWriter (dir+"/"+cfastaFileName));
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));
			while(ER.more()){
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else
					ConvertRmapper2csfasta(ER,EW);
			}
			ER.close();
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}

	}

	void convertRmapper2Fasta(String dir, String fileName){	
		try{
			String cfastaFileName = null;
			if(fileName.endsWith(".rmapper")){
				cfastaFileName = fileName.substring(0,fileName.length()-".rmapper".length())+".fa";
			}
			else
				cfastaFileName = fileName+".fa";
				
			ExtendedWriter EW = new ExtendedWriter(new FileWriter (dir+"/"+cfastaFileName));
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));
			while(ER.more()){
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else
					ConvertRmapper2fasta(ER,EW);
			}
			ER.close();
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}

	}
	
	
	static void printBestRmapperHits(String Dir, String fileName, int size,String filter){
		try{
			int start = 0; 
			int stop = size -1;
			String extension = "_"+start+"_"+stop+"."+filter+".rmapper";
			System.out.println(Dir+"/"+fileName+extension);
			while(IOTools.fileExists(Dir,fileName+extension)){
				CfastaSequences CFsequences  = new CfastaSequences();

				ExtendedReader ER = new ExtendedReader(new FileReader (Dir+"/"+fileName+extension));
				while(ER.more()){
					if((char)ER.lookAhead() == '#'){
						ER.skipLine();
					}
					else
						CFsequences.addRmapperHit(ER);
				}
				ER.close();

				ExtendedWriter EW = new ExtendedWriter(new FileWriter (Dir+"/"+fileName+extension+".best"));
				CFsequences.printRmapperHits(EW);
				EW.flush();
				EW.close();
				
				start = start + size;
				stop = stop + size;
				extension = "_"+start+"_"+stop+"."+filter+".rmapper";
			}
		}
		catch(Exception E){E.printStackTrace();}
	}
	
	static int[] printOrderedRmapperHits(String Dir, String fileName){
		try{
			int count = 0;
			ExtendedWriter [] EWs = new ExtendedWriter[31];
			int[] counts = new int[31];
			for(int i = 1; i < 31; i++)
				EWs[i] = new ExtendedWriter(new FileWriter (Dir+"/"+fileName+"."+i));

			EWs[0] = new ExtendedWriter(new FileWriter (Dir+"/"+fileName+".rest"));

			CfastaSequences CFsequences  = new CfastaSequences();

			ExtendedReader ER = new ExtendedReader(new FileReader (Dir+"/"+fileName));
			Solid oldHit = null;
			while(ER.more()){
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					Solid newHit = CFsequences.addRmapperHit(ER,oldHit);
					if(oldHit != null && newHit.compareTo(oldHit) != 0){
						count++;
						//System.out.println(count + "\t"+oldHit.hits.size() );
						if(oldHit.hits.size() >30 ){
							oldHit.printrMapperHits(EWs[0]);
							counts[30]++;
						}
						else{
							oldHit.printrMapperHits(EWs[oldHit.hits.size()]);
							counts[oldHit.hits.size()-1]++;
						}
					}
					oldHit = newHit;
				}
			}

			ER.close();
			for(int i = 0; i <31; i++){
				EWs[i].flush();
				EWs[i].close();
				System.out.println(i+"\t"+counts[i]);
			}
			
			
			return counts;
		}
		catch(Exception E){E.printStackTrace();}
		return null;
	}

	
	
	
	void addSolidSequences(String dir, String fileName, int max){
		
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));
			int count = 100000;
			while(ER.more() && this.size() < max){
				
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					addcFasta(ER,"wrong");
				}
				if(this.size() > count){
					System.out.println("nr Of Sequences Read" + count);
					count= count + 100000;
				}
			}
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}

	public static void printSizeDistribution(String dir, String fileName, int max){
		CfastaSequences A = new CfastaSequences();
		int[] counter = A.getSizeDistribution(dir, fileName, max);
		for(int i = 0; i < counter.length;i++){
			System.out.println(i-1+" "+counter[i]);
		}
		System.out.println();
	} 
	
	private int[] getSizeDistribution(String dir, String fileName, int max){
		int count = 0;
		try{
			
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));
			
			int[] counter = new int[max];
			while(ER.more()){
				while(ER.more()){
					if((char)ER.lookAhead() == '#'){
						ER.skipLine();
					}
					else{
						int length = getLength(ER);
						if(length < max)
							counter[length]++;
					}
				}
			}
			ER.close();
			return counter;
			}catch(Exception E){E.printStackTrace();}
			return null;
	}

	
	private void getSizeDistribution(String dir, String fileName, String outFileName, int size){
		int count = 0;
		try{
			
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));
			ExtendedWriter EW = new ExtendedWriter(new FileWriter (dir+"/"+outFileName));
			
			int[] counter = new int[size];
			while(ER.more()){
				while(ER.more()){
					if((char)ER.lookAhead() == '#'){
						ER.skipLine();
					}
					else{
						int length = getLength(ER);
						if(length < size)
							counter[length]++;
					}
				}
			}
			ER.close();
			}catch(Exception E){E.printStackTrace();}
	}
	
	
	private int[] getSizeDistribution(String dir, String fileName, int max, String[] primerSequences){
		int count = 0;
		try{
			
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));
			
			int[] counter = new int[50];
			while(ER.more()){
				while(ER.more() && this.size() < max){
					if((char)ER.lookAhead() == '#'){
						ER.skipLine();
					}
					else{
						addcFasta(ER,"wrong");
					}
				}
				
				System.out.println("nr Of Sequences Read: " + count);
				System.out.println("finding primers");
				CfastaSequences notFound = new CfastaSequences();
				counter = findPrimers(primerSequences,counter,notFound,10);
				
				printSequences(dir,fileName+"_"+count+"_"+(count+max)+".3adapter");
				notFound.printSequences(dir,fileName+"_"+count+"_"+(count+max)+".notFound");
				
				count = count + max;
			}
	/*		
			for(int i = 0; i < counter.length; i++){
				System.out.println(i + "\t"+ counter[i]) ;
			}
			
			System.out.println(this.size());
*/			ER.close();
			return counter;
			}catch(Exception E){E.printStackTrace();}
			return null;
			
	}
	
	
	public static CfastaSequences convertFasta2CSfasta(FastaSequences Fseq){
		CfastaSequences CSseq = new  CfastaSequences();
		for(int i = 0; i < Fseq.size();i++){
			Solid temp = new Solid();
			temp.addSequence(Fseq.get(i));
			CSseq.add(temp);
		}
		int test = 0;
		int size = CSseq.size();
		return CSseq;
	}

	
	private int[] splitFile(String dir, String fileName, int maxSize){
		int count = 0;
		try{
			
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));
			
			while(ER.more()){
				while(ER.more() && this.size() < maxSize){
					if((char)ER.lookAhead() == '#'){
						ER.skipLine();
					}
					else{
						addcFastaEasy(ER);
					}
				}
				
				System.out.println("nr Of Sequences Read: " + count);
				System.out.println("finding primers");
				
				printSequences(dir,fileName+"_"+count+"_"+(count+this.size())+".3adapter");
				count = count + maxSize;
			}
			ER.close();
			}catch(Exception E){E.printStackTrace();}
			return null;
			
	}
	
	
	
	
	private void addSolidSequences(String dir, String fileName, String chromsome){
		
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));
			int count = 100000;
			while(ER.more()){
				
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					addcFasta(ER, chromsome);
				}
				if(this.size() > count){
					System.out.println("nr Of Sequences Read" + count);
					count= count + 100000;
				}
			}
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}


	private void addSolidSequences(ArrayList<Solid> sequences, String dir, String fileName){
		
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));
			int count = 100000;
			while(ER.more()){
				
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					addcFasta(ER);
				}
				if(this.size() > count){
					System.out.println("nr Of Sequences Read" + count);
					count= count + 100000;
				}
			}
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}
	
	
	private void addSolidSequencesSimplest(String dir, String fileName){
		
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));
			int count = 100000;
			while(ER.more()){
				
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					addcFastaSimplest(ER);
				}
				if(this.size() > count){
					System.out.println("nr Of Sequences Read" + count);
					count= count + 100000;
				}
			}
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}
	
	private void printInfoEasy(String cfastaDir , String cfastaFile, String compareFastaFile, String compareDir, int nrOfSequences, int nrOfcomparingSequences, int nrOfuniqueSequences, String outDir,String outFile){
		System.out.print("Printing sequences to file "+outDir+"/"+outFile+  "  .......");
		if(!IOTools.isDir(outDir)) 
		IOTools.mkDir(outDir);
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile));
			EW.println("# original file =  "+ cfastaDir +"/"+cfastaFile );
			EW.println("# compare file = "+ compareDir +"/"+compareFastaFile );
			EW.println("# Nr of sequences in original file :"+nrOfSequences );
			EW.println("# Nr of sequences in compare file :"+nrOfcomparingSequences );
			EW.println("# Nr of unique sequences in original file :"+nrOfuniqueSequences );
			EW.println("# Nr of unique sequences in compare file :"+(nrOfcomparingSequences - (nrOfSequences - nrOfuniqueSequences)));
			EW.println("# Nr of same seqences in both files :"+(nrOfSequences - nrOfuniqueSequences) );
			
			for(int i = 0; i < this.size(); i++){
				this.get(i).printInfoEasy(EW);
			}
			EW.flush();
			EW.close();
			
		}catch(Exception E){E.printStackTrace();}
		System.out.println("finished");
	}


	private void printInfo(String cfastaDir , String cfastaFile, String compareFastaFile, String compareDir, int nrOfSequences, int nrOfcomparingSequences, int nrOfuniqueSequences, String outDir,String outFile){
		System.out.print("Printing sequences to file "+outDir+"/"+outFile+  "  .......");
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile));
			EW.println("# original file =  "+ cfastaDir +"/"+cfastaFile );
			EW.println("# compare file = "+ compareDir +"/"+compareFastaFile );
			EW.println("# Nr of sequences in original file :"+nrOfSequences );
			EW.println("# Nr of sequences in compare file :"+nrOfcomparingSequences );
			EW.println("# Nr of unique sequences in original file :"+nrOfuniqueSequences );
			EW.println("# Nr of unique sequences in compare file :"+(nrOfcomparingSequences - (nrOfSequences - nrOfuniqueSequences)));
			EW.println("# Nr of same seqences in both files :"+(nrOfSequences - nrOfuniqueSequences) );
			
			for(int i = 0; i < this.size(); i++){
				this.get(i).printcFasta(EW);
			}
			EW.flush();
			EW.close();
			
		}catch(Exception E){E.printStackTrace();}
		System.out.println("finished");
	}

	private void printJoinedInfo(String outDir,String outFile,String[] chromosomes ){
		System.out.print("Printing sequences to file "+outDir+"/"+outFile+  ".joined  .......");
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile+".joined"));
			EW.println("# joined sequences from multiple files:");
			for(int i = 0; i < chromosomes.length; i++){
				EW.println("# joined file =  "+ outDir +"/"+ chromosomes[i]+"/"+outFile );
			}
			EW.println("#Total number of unique sequences"+  this.size());
			for(int i = 0; i < this.size(); i++){
				this.get(i).printcFasta(EW);
			}
			EW.flush();
			EW.close();
			
		}catch(Exception E){E.printStackTrace();}
		System.out.println("finished");
	}
	
	private void printInfo(String outDir,String outFile){
		System.out.print("Printing sequences to file "+outDir+"/"+outFile+  ".joined  .......");
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile+".joined"));
			for(int i = 0; i < this.size(); i++){
				this.get(i).printcFasta(EW);
			}
			EW.flush();
			EW.close();
			
		}catch(Exception E){E.printStackTrace();}
		System.out.println("finished");
	}
	
	private void printSequences(ExtendedWriter EW){
		for(int i = 0; i < this.size(); i++){
			this.get(i).printcFasta(EW);
		}

		
	}
	
	public void printRmapperHits(ExtendedWriter EW){
		for(int i = 0; i < this.size(); i++){
			this.get(i).printrMapperHits(EW);
		}
	}
	
	public void printRedundantRmapperHits(ExtendedWriter EW, int cutoff){
		for(int i = 0; i < this.size(); i++){
			this.get(i).printRedundantrMapperHits(EW, cutoff);
		}
	}

	
	
	private void printSequences(String outDir,String outFile){
		System.out.print("Printing sequences to file "+outDir+"/"+outFile+"  .......");
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile));
			printSequences(EW);
			EW.flush();
			EW.close();
			
		}catch(Exception E){E.printStackTrace();}
		System.out.println("finished");
	}

	private void printSingle(String outDir,String outFile){
		System.out.print("Printing sequences to file "+outDir+"/"+outFile+  ".single  .......");
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile+".single"));
			EW.println("# fasta sequence "+outDir +"/"+ outFile +" with one unique hit :");
			EW.println("#Total number of sequences with one hit: "+  this.size());
			printSequences(EW);
			EW.flush();
			EW.close();
			
		}catch(Exception E){E.printStackTrace();}
		System.out.println("finished");
	}

	public void printNrOfHits(String outDir,String outFile, int nrOfHits){
		System.out.print("Printing sequences to file "+outDir+"/"+outFile+  "."+nrOfHits+"Hits  .......");
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile+"."+nrOfHits+"Hits"));
			EW.println("# cfasta sequences "+outDir +"/"+ outFile +" with "+ nrOfHits+" hits :");
			EW.println("#Total number of sequences with "+ nrOfHits+" hits: "+  this.size());
			printSequences(EW);
			EW.flush();
			EW.close();
			
		}catch(Exception E){E.printStackTrace();}
		System.out.println("finished");
	}
	
	
	
/*	private final int mRNA = 1;
	private final int ncRNA = 2;
	private final int intergenic = 3;
	private final int antisense = 4;
*/
	
/*	1=mRNA
	2=ncRNA
	3=intergenic
	4=antisense
	5=repeats
*/
	
	public void printGroups(String outDir,String outFile){
		for(int i = 1; i < 6; i++){
			System.out.println("Printing different groupsequences to file "+outDir+"/"+outFile+  ".###  .......");
			try{
				ExtendedWriter EW = null;
				if(i == 1){
					System.out.print("Printing mRNAsequences to file "+outDir+"/"+outFile+  ".mRNA  .......");
					EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile+".mRNA"));
				}
				if(i == 2){
					System.out.print("Printing ncRNA sequences to file "+outDir+"/"+outFile+  ".ncRNA  .......");
					EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile+".ncRNA"));
				}
				if(i == 3){
					System.out.print("Printing intergenic sequences to file "+outDir+"/"+outFile+  ".intergenic  .......");
					EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile+".intergenic"));
				}
				if(i == 4){
					System.out.print("Printing antisense sequences to file "+outDir+"/"+outFile+  ".antisense  .......");
					EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile+".antisense"));
				}
				if(i == 5){
					System.out.print("Printing repeats sequences to file "+outDir+"/"+outFile+  ".repeats  .......");
					EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile+".repeats"));
				}

				//EW.println("# fasta sequence "+outDir +"/"+ outFile +" with one unique hit :");
				//EW.println("#Total number of sequences with one hit: "+  this.size());
				for(int j = 0; j < this.size(); j++){
					this.get(j).printcFasta(EW,i);
				}
				EW.flush();
				EW.close();

			}catch(Exception E){E.printStackTrace();}
			System.out.println("finished");
		}
	}


	public int[]  countGroups(){
		int counts[] = new int[7];
		
		/*	1=mRNA
			2=ncRNA
			3=intergenic
			4=antisense
			5=repeats
			6=DIRS-1
		*/
		
			for(int j = 0; j < this.size(); j++){
				this.get(j).countKinds(counts);
			}
		return counts;
	}
	
	public double[]  countNormalizedGroups(){
		double counts[] = new double[7];
		
		/*	1=mRNA
			2=ncRNA
			3=intergenic
			4=antisense
			5=repeats
			6=DIRS-1
		*/
		
			for(int j = 0; j < this.size(); j++){
				counts = this.get(j).countNormalizedKinds(counts);
			}
			
			
			for(int i = 0; i < counts.length ; i++){
				System.out.println(i+"\t"+counts[i]);
			}
			
			
		return counts;
	}

	
	
	

	
	private void addcFasta(ExtendedReader ER){
		String InfoLine = ER.readLine();
		//>853_36_1821_F3,187_50.1,155_50.1,149_50.1,137_50.1,130_50.1,24_50.1,18_50.1
		
		String Sequence = ER.readLine();
		//T1002123.321003210102333300010033020103031311211220
		
		String[] hits = InfoLine.split(",");
		Solid newHit = new Solid(hits[0],Sequence);
		
		
		if(hits.length > 1){
			for(int i = 1; i < hits.length; i++){
				String[] info = hits[i].split("\\.");
				int index = info[0].indexOf("_");
//				String fastaFile = ""; 
				int location = 0;
				String chromosome = "wrong";
				if(index > -1){
					chromosome = info[0].substring(0,index);
					location = Integer.parseInt(info[0].substring(index+1));
				}
				else{
					location = Integer.parseInt(info[0]);
				}
				int mm = Integer.parseInt(info[1]);
				int length = 20;
				if(info.length > 2)
					length = Integer.parseInt(info[2]);

				boolean plusStrand = true;
				if(location < 0){
//					location = this.length+location;
					plusStrand = false;
				}
				if(chromosome.compareTo("wrong") == 0)
					newHit.addHit(location, length, mm, plusStrand);
				else
					newHit.addHit(location, length, mm, plusStrand,chromosome);
			}
		}
		double location = findLocation(newHit,this);
		int pointer = (int)(location+0.5);
		this.add(pointer,newHit);
	}

	private int addRmapperHit(ExtendedReader ER){
		//#FORMAT: readname contigname strand contigstart contigend readstart readend readlength score editstring readsequence
		//>856_759_1307_F3        adaptor + bc1 + P2      +       1       25      15      39      50      208     x18x2x5 A22301003111312033020103031311231101032002322003120
		String[] info = ER.readLine().split("\t");

		String readname = info[0];    //>855_1396_1408_F3
		String contigname = info[1];  // adaptor + bc1 + P2 
		String strand = info[2];      //		+ 
		int contigstart = Integer.parseInt(info[3]);    // 		2
		int contigend = Integer.parseInt(info[4]);      //		25
		int readstart  = Integer.parseInt(info[5]);     // 		24
		int readend  = Integer.parseInt(info[6]);     // 		24
		int readlength = Integer.parseInt(info[7]);     //      	50
		int score = Integer.parseInt(info[8]);          //      	184
		String editstring = info[9];  //     	x10x2x8x4
		String readsequence = null;
		if(info.length > 10)
			readsequence = info[10];//      A12331211230122022221321230201030323212312003320120

		Solid newHit = new Solid(readname,readsequence);

		double location = findLocation(newHit,this);
		boolean plusStrand = true;
		if(strand.compareTo("+") != 0){
			plusStrand = false;
		}


		Hit hit = new Hit(contigname, plusStrand, contigstart, contigend, readstart,
				readend, readlength, score,editstring);
		int pointer  = (int)location;
		if(location - (int)(location+0.5) == 0)
			this.get(pointer).addHit(hit);
		else{
			newHit.addHit(hit);
			pointer = (int)(location+0.5);
			this.add(pointer,newHit);

		}
		return pointer;
	}

	private int addRmapperHit(ExtendedReader ER, int pointer){
		//#FORMAT: readname contigname strand contigstart contigend readstart readend readlength score editstring readsequence
		//>856_759_1307_F3        adaptor + bc1 + P2      +       1       25      15      39      50      208     x18x2x5 A22301003111312033020103031311231101032002322003120
		String[] info = ER.readLine().split("\t");

		String readname = info[0];    //>855_1396_1408_F3
		String contigname = info[1];  // adaptor + bc1 + P2 
		String strand = info[2];      //		+ 
		int contigstart = Integer.parseInt(info[3]);    // 		2
		int contigend = Integer.parseInt(info[4]);      //		25
		int readstart  = Integer.parseInt(info[5]);     // 		24
		int readend  = Integer.parseInt(info[6]);     // 		24
		int readlength = Integer.parseInt(info[7]);     //      	50
		int score = Integer.parseInt(info[8]);          //      	184
		String editstring = info[9];  //     	x10x2x8x4
		String readsequence = null;
		if(info.length > 10)
			readsequence = info[10];//      A12331211230122022221321230201030323212312003320120

		boolean plusStrand = true;
		if(strand.compareTo("+") != 0){
			plusStrand = false;
		}
		Solid newHit = new Solid(readname,readsequence);
		Hit hit = new Hit(contigname, plusStrand, contigstart, contigend, readstart,
				readend, readlength, score,editstring);

		if(pointer != -1 && this.get(pointer).compareTo(newHit) == 0){
			this.get(pointer).addHit(hit);
		}
		else{



			double location = findLocation(newHit,this);

			pointer  = (int)location;
			if(location - (int)(location+0.5) == 0)
				this.get(pointer).addHit(hit);
			else{
				newHit.addHit(hit);
				pointer = (int)(location+0.5);
				this.add(pointer,newHit);

			}
		}
		return pointer;
	}

	private Solid addRmapperHit(ExtendedReader ER, Solid oldHit){
		//#FORMAT: readname contigname strand contigstart contigend readstart readend readlength score editstring readsequence
		//>856_759_1307_F3        adaptor + bc1 + P2      +       1       25      15      39      50      208     x18x2x5 A22301003111312033020103031311231101032002322003120
		String[] info = ER.readLine().split("\t");

		String readname = info[0];    //>855_1396_1408_F3
		String contigname = info[1];  // adaptor + bc1 + P2 
		String strand = info[2];      //		+ 
		int contigstart = Integer.parseInt(info[3]);    // 		2
		int contigend = Integer.parseInt(info[4]);      //		25
		int readstart  = Integer.parseInt(info[5]);     // 		24
		int readend  = Integer.parseInt(info[6]);     // 		24
		int readlength = Integer.parseInt(info[7]);     //      	50
		int score = Integer.parseInt(info[8]);          //      	184
		String editstring = info[9];  //     	x10x2x8x4
		String readsequence = null;
		if(info.length > 10)
			readsequence = info[10];//      A12331211230122022221321230201030323212312003320120

		boolean plusStrand = true;
		if(strand.compareTo("+") != 0){
			plusStrand = false;
		}
		Solid newHit = new Solid(readname,readsequence);
		Hit hit = new Hit(contigname, plusStrand, contigstart, contigend, readstart,
				readend, readlength, score,editstring);

		if(oldHit != null && oldHit.compareTo(newHit) == 0){
			oldHit.addHit(hit);
		}
		else{
				newHit.addHit(hit);
				return newHit;
		}
		return oldHit;
	}

	private int addRmapperMoreThanHit(ExtendedReader ER, int pointer,int hits){
		//#FORMAT: readname contigname strand contigstart contigend readstart readend readlength score editstring readsequence
		//>856_759_1307_F3        adaptor + bc1 + P2      +       1       25      15      39      50      208     x18x2x5 A22301003111312033020103031311231101032002322003120
		String[] info = ER.readLine().split("\t");

		String readname = info[0];    //>855_1396_1408_F3
		String contigname = info[1];  // adaptor + bc1 + P2 
		String strand = info[2];      //		+ 
		int contigstart = Integer.parseInt(info[3]);    // 		2
		int contigend = Integer.parseInt(info[4]);      //		25
		int readstart  = Integer.parseInt(info[5]);     // 		24
		int readend  = Integer.parseInt(info[6]);     // 		24
		int readlength = Integer.parseInt(info[7]);     //      	50
		int score = Integer.parseInt(info[8]);          //      	184
		String editstring = info[9];  //     	x10x2x8x4
		String readsequence = null;
		if(info.length > 10)
			readsequence = info[10];//      A12331211230122022221321230201030323212312003320120

		boolean plusStrand = true;
		if(strand.compareTo("+") != 0){
			plusStrand = false;
		}
		Solid newHit = new Solid(readname,readsequence);
		Hit hit = new Hit(contigname, plusStrand, contigstart, contigend, readstart,
				readend, readlength, score,editstring);

		if(pointer != -1 && this.get(pointer).compareTo(newHit) == 0){
			this.get(pointer).addHit(hit);
		}
		else{
			if(pointer != -1 &&this.get(pointer).hits.size() < hits)
				this.remove(pointer);

			double location = findLocation(newHit,this);

			pointer  = (int)location;
			if(location - (int)(location+0.5) == 0)
				this.get(pointer).addHit(hit);
			else{
				newHit.addHit(hit);
				pointer = (int)(location+0.5);
				this.add(pointer,newHit);

			}
		}
		return pointer;
	}

	
	
	
	
	

	
	private static void extractGmapperInfo(ExtendedReader ER, String miRNAName, String sampleName ){
		String[] info = ER.readLine().split("\t");
		
		
		String readname = info[0];    //>855_1396_1408_F3
		String contigname = info[1];  // adaptor + bc1 + P2 
		String strand = info[2];      //		+ 
		int contigstart = Integer.parseInt(info[3]);    // 		2
		int contigend = Integer.parseInt(info[4]);      //		25
		int readstart  = Integer.parseInt(info[5]);     // 		24
		int readend  = Integer.parseInt(info[6]);     // 		24
		int readlength = Integer.parseInt(info[7]);     //      	50
		int score = Integer.parseInt(info[8]);          //      	184
		String editstring = info[9];  //     	x10x2x8x4
		String readsequence = null;
		if(info.length > 10)
			readsequence = info[10];//      A12331211230122022221321230201030323212312003320120

		System.out.println(miRNAName+"\t"+sampleName+"\t"+contigstart+"\t"+readlength);
	
	
	}	
	private void RmapperRemoveAdapter(ExtendedReader ER,ExtendedWriter[] EW, int maxLength){
//#FORMAT: readname contigname strand contigstart contigend readstart readend readlength score editstring readsequence
//>856_759_1307_F3        adaptor + bc1 + P2      +       1       25      15      39      50      208     x18x2x5 A22301003111312033020103031311231101032002322003120
		String[] info = ER.readLine().split("\t");
		
		
		String readname = info[0];    //>855_1396_1408_F3
		String contigname = info[1];  // adaptor + bc1 + P2 
		String strand = info[2];      //		+ 
		int contigstart = Integer.parseInt(info[3]);    // 		2
		int contigend = Integer.parseInt(info[4]);      //		25
		int readstart  = Integer.parseInt(info[5]);     // 		24
		int readend  = Integer.parseInt(info[6]);     // 		24
		int readlength = Integer.parseInt(info[7]);     //      	50
		int score = Integer.parseInt(info[8]);          //      	184
		String editstring = info[9];  //     	x10x2x8x4
		String readsequence = null;
		if(info.length > 10)
			readsequence = info[10];//      A12331211230122022221321230201030323212312003320120

		if(readstart-contigstart > -1 && readstart-contigstart < maxLength){
			String sequence = readsequence.substring(0,readstart+1-contigstart);
			EW[readstart-contigstart].println(readname);
			EW[readstart-contigstart].println(sequence);
		}
	}
	
	private void RmapperGetAdapterSizes(ExtendedReader ER,int[] lengths, int maxLength){
		//#FORMAT: readname contigname strand contigstart contigend readstart readend readlength score editstring readsequence
		//>856_759_1307_F3        adaptor + bc1 + P2      +       1       25      15      39      50      208     x18x2x5 A22301003111312033020103031311231101032002322003120

		String[] info = ER.readLine().split("\t");


		int contigstart = Integer.parseInt(info[3]);    // 		2

		int readstart  = Integer.parseInt(info[5]);     // 		24


		if(readstart-contigstart > -1 && readstart-contigstart < maxLength){
			lengths[readstart-contigstart]++;
		}
	}

	
	
	
	private void ConvertRmapper2csfasta(ExtendedReader ER,ExtendedWriter EW){
		//#FORMAT: readname contigname strand contigstart contigend readstart readend readlength score editstring readsequence
		//>856_759_1307_F3        adaptor + bc1 + P2      +       1       25      15      39      50      208     x18x2x5 A22301003111312033020103031311231101032002322003120
		
		String[] info = ER.readLine().split("\t");


		String readname = info[0];    //>855_1396_1408_F3
		String contigname = info[1];  // adaptor + bc1 + P2 
		String strand = info[2];      //		+ 
		int contigstart = Integer.parseInt(info[3]);    // 		2
		int contigend = Integer.parseInt(info[4]);      //		25
		int readstart  = Integer.parseInt(info[5]);     // 		24
		int readend  = Integer.parseInt(info[6]);     // 		24
		int readlength = Integer.parseInt(info[7]);     //      	50
		int score = Integer.parseInt(info[8]);          //      	184
		String editstring = info[9];  //     	x10x2x8x4
		String readsequence = null;
		if(info.length > 10){
			readsequence = info[10];//      A12331211230122022221321230201030323212312003320120
			EW.println(readname);
			EW.println(readsequence);
		}
	}
	
	
	private void ConvertRmapper2fasta(ExtendedReader ER,ExtendedWriter EW){
		//#FORMAT: readname contigname strand contigstart contigend readstart readend readlength score editstring readsequence
		//>856_759_1307_F3        adaptor + bc1 + P2      +       1       25      15      39      50      208     x18x2x5 A22301003111312033020103031311231101032002322003120
		
		String[] info = ER.readLine().split("\t");


		String readname = info[0];    //>855_1396_1408_F3
		String contigname = info[1];  // adaptor + bc1 + P2 
		String strand = info[2];      //		+ 
		int contigstart = Integer.parseInt(info[3]);    // 		2
		int contigend = Integer.parseInt(info[4]);      //		25
		int readstart  = Integer.parseInt(info[5]);     // 		24
		int readend  = Integer.parseInt(info[6]);     // 		24
		int readlength = Integer.parseInt(info[7]);     //      	50
		int score = Integer.parseInt(info[8]);          //      	184
		String editstring = info[9];  //     	x10x2x8x4
		String readsequence = null;
		if(info.length > 10){
			readsequence = info[10];//      A12331211230122022221321230201030323212312003320120
			EW.println(readname);
			String fastaSequence = RNAfunctions.cfasta2fasta(readsequence);
			EW.println(fastaSequence);
		}
	}
	
	
	private int getLength(ExtendedReader ER){
		String InfoLine = ER.readLine();
		String Sequence = ER.readLine();
		return Sequence.length();
		
	}
	
//	private int getCFastaSequence(ExtendedReader ER){
//		String InfoLine = ER.readLine();
//		String Sequence = ER.readLine();
//		CfastaSequence CSF = new 
//	}
	
	private void addcFasta(ExtendedReader ER,String chromosome){
		String InfoLine = ER.readLine();
		String Sequence = ER.readLine();

		String[] hits = InfoLine.split(",");
		Solid newHit = new Solid(hits[0],Sequence);
		if(hits.length > 1){
			for(int i = 1; i < hits.length; i++){
				String[] info = hits[i].split("\\.");
				int index = info[0].indexOf("_");
				//							String fastaFile = ""; 
				int location = 0;
				if(index > -1){
					//chromosome = info[0].substring(0,index);
					location = Integer.parseInt(info[0].substring(index+1));
				}
				else{
					location = Integer.parseInt(info[0]);
				}
				int mm = Integer.parseInt(info[1]);
				int length = 20;
				if(info.length > 2)
					length = Integer.parseInt(info[2]);

				boolean plusStrand = true;
				if(location < 0){
					//								location = this.length+location;
					plusStrand = false;
				}
				newHit.addHit(location, length, mm, plusStrand,chromosome);
			}
		}
		double location = findLocation(newHit,this);
		int pointer = (int)(location+0.5);
		this.add(pointer,newHit);
		

	}

	
	
	private void addcFastaEasy(ExtendedReader ER){
		String InfoLine = ER.readLine();
		//>853_36_1821_F3,187_50.1,155_50.1,149_50.1,137_50.1,130_50.1,24_50.1,18_50.1
		
		String SequenceLine = ER.readLine();
		//T1002123.321003210102333300010033020103031311211220
		
		String[] hits = InfoLine.split(",");
		Solid newHit = new Solid(hits[0]);
		newHit.addInfo(InfoLine, SequenceLine);
		
		double location = findLocation(newHit,this);
		int pointer = (int)(location+0.5);
		this.add(pointer,newHit);
		
	}
	
	
	private void addcFastaSimplest(ExtendedReader ER){
		String InfoLine = ER.readLine();
		//>853_36_1821_F3,187_50.1,155_50.1,149_50.1,137_50.1,130_50.1,24_50.1,18_50.1
		
		ER.skipLine();
		//T1002123.321003210102333300010033020103031311211220
		
		String[] hits = InfoLine.split(",");
		Solid newHit = new Solid(hits[0]);
		
		double location = findLocation(newHit,this);
		int pointer = (int)(location+0.5);
		this.add(pointer,newHit);
		
	}
	
	private static double findLocation(Solid newSequence,ArrayList<Solid> sortedSequences){
		if(sortedSequences == null || sortedSequences.size() == 0)
			return -0.5;
		int leftPointer = 0;
		int rightPointer = sortedSequences.size()-1;
		
		if(newSequence.compareTo(sortedSequences.get(rightPointer)) > 0){
			return (double)rightPointer+0.5;
		} 
		if(newSequence.compareTo(sortedSequences.get(leftPointer)) < 0){
			return (double)leftPointer-0.5;
		} 
		
		while(rightPointer - leftPointer > 1){
			int middle = (rightPointer+leftPointer)/2;
			
			double location = newSequence.compareTo(sortedSequences.get(middle));
			if(location < 0) rightPointer = middle;
			else if(location > 0) leftPointer = middle;
			else return middle;
		}
		
		double location = newSequence.compareTo(sortedSequences.get(leftPointer));
		if(location <0) System.out.println("Something is very wrong with the finder");
		else if(location == 0) return leftPointer;

		location = newSequence.compareTo(sortedSequences.get(rightPointer));
		if(location > 0) System.out.println("Something is very wrong with the finder");
		else if(location == 0) return rightPointer;
		
		return (double)leftPointer+0.5;
	}
	
	
	
	public static CfastaSequences singleHits(CfastaSequences unsortedSequences ){
		CfastaSequences singleHits = new CfastaSequences(unsortedSequences.getName());
		for(int i = 0; i < unsortedSequences.size(); i++){
			if(unsortedSequences.get(i).hits.size() ==1)
				singleHits.add(unsortedSequences.get(i));
		}
		return singleHits;
	}

	
	public static CfastaSequences getFewerHits(CfastaSequences unsortedSequences, int nrOfHits ){
		CfastaSequences singleHits = new CfastaSequences(unsortedSequences.getName());
		for(int i = 0; i < unsortedSequences.size(); i++){
			if(unsortedSequences.get(i).hits.size() <= nrOfHits)
				singleHits.add(unsortedSequences.get(i));
		}
		return singleHits;
	}

	public static CfastaSequences getNrOfHits(CfastaSequences unsortedSequences, int nrOfHits ){
		CfastaSequences singleHits = new CfastaSequences(unsortedSequences.getName());
		for(int i = 0; i < unsortedSequences.size(); i++){
			if(unsortedSequences.get(i).hits.size() == nrOfHits)
				singleHits.add(unsortedSequences.get(i));
		}
		return singleHits;
	}
	
	
	
	public static int[] getHitDistribution(ArrayList <Solid> unsortedSequences){
		int[] nrOfHits = new int[100];
		for(int i = 0; i < unsortedSequences.size(); i++){
			nrOfHits[unsortedSequences.get(i).hits.size()]++;
		}
		
		return nrOfHits;
	}
	
	public static void printDistribution(int[] Distribution, int start, int step){
		for(int i = 0; i< Distribution.length; i++)
			System.out.println((start+i*step)+"\t"+Distribution[i]);
	}

	
	
	private static ArrayList <Solid> sortSequences(ArrayList <Solid> unsortedSequences ){
		ArrayList <Solid> sortedSequences = new ArrayList<Solid>();
		for(int i = 0; i < unsortedSequences.size(); i++){
			double pointer = findLocation(unsortedSequences.get(i),sortedSequences);
			int location = (int)(pointer+0.5);
			sortedSequences.add(location, unsortedSequences.get(i));
		}
		return sortedSequences;
	}
	
	
	private static ArrayList <Solid> joinSequences(ArrayList <Solid> ss1 ,ArrayList <Solid> ss2){
		
		 System.out.println("Sort sequences");
		 int ss1Pointer = 0;
		 int ss2Pointer = 0;
		 int ss1Size = ss1.size();
		 int ss2Size = ss2.size();
		 ArrayList <Solid> JoinSequences = new ArrayList <Solid>(); 
		 int ss1single = 0;
		 int ss2single = 0;
		  int same = 0;
		 
		 while(ss1Pointer < ss1Size && ss2Pointer < ss2Size){
			while(ss1Pointer < ss1Size && ss1.get(ss1Pointer).compareTo(ss2.get(ss2Pointer)) < 0){
				JoinSequences.add(ss1.get(ss1Pointer));
					 ss1Pointer++;
					 ss1single++;
			}
			while(ss1Pointer < ss1Size && ss2Pointer < ss2Size && ss1.get(ss1Pointer).compareTo(ss2.get(ss2Pointer)) > 0){
				JoinSequences.add(ss2.get(ss2Pointer));
				ss2Pointer++;
				ss2single++;
			}
			if(ss1Pointer < ss1Size && ss2Pointer < ss2Size && ss1.get(ss1Pointer).compareTo(ss2.get(ss2Pointer)) == 0){
				Solid joinedHit =  ss1.get(ss1Pointer);
				joinedHit.addHits(ss2.get(ss2Pointer));
				JoinSequences.add(joinedHit);
				 ss2Pointer++;
				 ss1Pointer++;
				 same++;
			 }
		 }
		 while(ss1Pointer < ss1Size){
			 JoinSequences.add(ss1.get(ss1Pointer));
			ss1Pointer++;
			ss1single++;
		}
		 
		while(ss2Pointer < ss2Size){
			JoinSequences.add(ss2.get(ss2Pointer));
			ss2Pointer++;
			ss2single++;
		}
		 
		 
		 System.out.println("Number of sequences in cfastaFile1 =        "+ ss1Size);
		 System.out.println("Number of sequences in cfastaFile2 =        "+ ss2Size);
		 System.out.println("Number of unique sequences in cfastaFile1 = "+ ss1single);
		 System.out.println("Number of unique sequences in cfastaFile2 = "+ ss2single);
		 System.out.println("same = "+ same);
		 return JoinSequences;
		 
		}
	private static ArrayList <Solid> replaceSequences(ArrayList <Solid> ss1 ,ArrayList <Solid> ss2){
		
		 System.out.println("Replace sequences");
		 int ss1Pointer = 0;
		 int ss2Pointer = 0;
		 int ss1Size = ss1.size();
		 int ss2Size = ss2.size();
		 ArrayList <Solid> JoinSequences = new ArrayList <Solid>(); 
		 int ss1single = 0;
		 int ss2single = 0;
		  int same = 0;
		 
		 while(ss1Pointer < ss1Size && ss2Pointer < ss2Size){
			while(ss1Pointer < ss1Size && ss1.get(ss1Pointer).compareTo(ss2.get(ss2Pointer)) < 0){
				JoinSequences.add(ss1.get(ss1Pointer));
					 ss1Pointer++;
					 ss1single++;
			}
			while(ss1Pointer < ss1Size && ss2Pointer < ss2Size && ss1.get(ss1Pointer).compareTo(ss2.get(ss2Pointer)) > 0){
				JoinSequences.add(ss2.get(ss2Pointer));
				ss2Pointer++;
				ss2single++;
			}
			if(ss1Pointer < ss1Size && ss2Pointer < ss2Size && ss1.get(ss1Pointer).compareTo(ss2.get(ss2Pointer)) == 0){
				Solid joinedHit =  ss2.get(ss2Pointer);
				JoinSequences.add(joinedHit);
				 ss2Pointer++;
				 ss1Pointer++;
				 same++;
			 }
		 }
		 while(ss1Pointer < ss1Size){
			 JoinSequences.add(ss1.get(ss1Pointer));
			ss1Pointer++;
			ss1single++;
		}
		 
		while(ss2Pointer < ss2Size){
			JoinSequences.add(ss2.get(ss2Pointer));
			ss2Pointer++;
			ss2single++;
		}
		 
		 
		 System.out.println("Number of sequences in cfastaFile1 =        "+ ss1Size);
		 System.out.println("Number of sequences in cfastaFile2 =        "+ ss2Size);
		 System.out.println("Number of unique sequences in cfastaFile1 = "+ ss1single);
		 System.out.println("Number of unique sequences in cfastaFile2 = "+ ss2single);
		 System.out.println("same = "+ same);
		 return JoinSequences;
		 
		}

	
	private static ArrayList <Solid> NANDSequences(ArrayList <Solid> ss1 ,ArrayList <Solid> ss2){
	
	 System.out.println("Sort sequences");
	 int ss1Pointer = 0;
	 int ss2Pointer = 0;
	 int ss1Size = ss1.size();
	 int ss2Size = ss2.size();
	 ArrayList <Solid> NandSequences = new ArrayList <Solid>(); 
	 int ss1single = 0;
	 int ss2single = 0;
	  int same = 0;
	 
	 while(ss1Pointer < ss1Size && ss2Pointer < ss2Size){
		while(ss1Pointer < ss1Size && ss1.get(ss1Pointer).compareTo(ss2.get(ss2Pointer)) < 0){
				 NandSequences.add(ss1.get(ss1Pointer));
				 ss1Pointer++;
				 ss1single++;
		}
		while(ss2Pointer < ss2Size && ss1.get(ss1Pointer).compareTo(ss2.get(ss2Pointer)) > 0){
			NandSequences.add(ss2.get(ss2Pointer));
			ss2Pointer++;
			ss2single++;
		}
		if(ss1Pointer < ss1Size && ss2Pointer < ss2Size && ss1.get(ss1Pointer).compareTo(ss2.get(ss2Pointer)) == 0){
			 ss2Pointer++;
			 ss1Pointer++;
			 same++;
		 }
	 }
	 while(ss1Pointer < ss1Size){
		NandSequences.add(ss1.get(ss1Pointer));
		ss1Pointer++;
		ss1single++;
	}
	while(ss2Pointer < ss2Size){
		NandSequences.add(ss2.get(ss2Pointer));
		ss2Pointer++;
		ss2single++;
	}
	 
	 
	 System.out.println("Number of sequences in cfastaFile1 =        "+ ss1Size);
	 System.out.println("Number of sequences in cfastaFile2 =        "+ ss2Size);
	 System.out.println("Number of unique sequences in cfastaFile1 = "+ ss1single);
	 System.out.println("Number of unique sequences in cfastaFile2 = "+ ss2single);
	 System.out.println("same = "+ same);
	 return NandSequences;
	 
	}

	private ArrayList <Solid> getUniqueSequences(ArrayList <Solid> otherSequences){
		
		 System.out.println("Sort sequences");
		 int ss1Pointer = 0;
		 int ss2Pointer = 0;
		 int ss1Size = this.size();
		 int ss2Size = otherSequences.size();
		 ArrayList <Solid> UniqueSequences = new ArrayList <Solid>(); 
		 int ss1single = 0;
		 int ss2single = 0;
		 
		 while(ss1Pointer < ss1Size && ss2Pointer < ss2Size){
			 
			while(ss1Pointer < ss1Size && this.get(ss1Pointer).compareTo(otherSequences.get(ss2Pointer)) < 0){
				UniqueSequences.add(this.get(ss1Pointer));
					 ss1Pointer++;
					 ss1single++;
			}
			while(ss2Pointer < ss2Size &&ss1Pointer < ss1Size  &&this.get(ss1Pointer).compareTo(otherSequences.get(ss2Pointer)) > 0){
				ss2Pointer++;
				 ss2single++;
			}
			if(ss1Pointer < ss1Size && ss2Pointer < ss2Size && this.get(ss1Pointer).compareTo(otherSequences.get(ss2Pointer)) == 0){
				 ss2Pointer++;
				 ss1Pointer++;
			 }
		 }
		 while(ss1Pointer < ss1Size ){
			UniqueSequences.add(this.get(ss1Pointer));
			ss1Pointer++;
			ss1single++;
		}

		 while(ss2Pointer < ss2Size ){
				ss2Pointer++;
				ss2single++;
			}
		 
		 System.out.println("Number of sequences in cfastaFile1 =        "+ ss1Size);
		 System.out.println("Number of sequences in cfastaFile2 =        "+ ss2Size);
		 System.out.println("Number of unique sequences in cfastaFile1 = "+ ss1single);
		 System.out.println("Number of unique sequences in cfastaFile2 = "+ ss2single);
		 return UniqueSequences;
		 
		}

	private void getUniqueSequences2(ArrayList <Solid> otherSequences){
		
		 System.out.println("Sort sequences");
		 int ss1Pointer = 0;
		 int ss2Pointer = 0;
		 int ss1Size = this.size();
		 int ss2Size = otherSequences.size();
		 CfastaSequences UniqueSequences = new CfastaSequences(this.getName()); 
		 int ss1single = 0;
		 int ss2single = 0;
		 
		 while(ss1Pointer < ss1Size && ss2Pointer < ss2Size){
			 
			while(ss1Pointer < ss1Size && this.get(ss1Pointer).compareTo(otherSequences.get(ss2Pointer)) < 0){
				UniqueSequences.add(this.get(ss1Pointer));
					 ss1Pointer++;
					 ss1single++;
			}
			while(ss2Pointer < ss2Size &&ss1Pointer < ss1Size  &&this.get(ss1Pointer).compareTo(otherSequences.get(ss2Pointer)) > 0){
				ss2Pointer++;
				 ss2single++;
			}
			if(ss1Pointer < ss1Size && ss2Pointer < ss2Size && this.get(ss1Pointer).compareTo(otherSequences.get(ss2Pointer)) == 0){
				 ss2Pointer++;
				 ss1Pointer++;
			 }
		 }
		 while(ss1Pointer < ss1Size ){
			UniqueSequences.add(this.get(ss1Pointer));
			ss1Pointer++;
			ss1single++;
		}

		 while(ss2Pointer < ss2Size ){
				ss2Pointer++;
				ss2single++;
			}
		 
		 System.out.println("Number of sequences in cfastaFile1 =        "+ ss1Size);
		 System.out.println("Number of sequences in cfastaFile2 =        "+ ss2Size);
		 System.out.println("Number of unique sequences in cfastaFile1 = "+ ss1single);
		 System.out.println("Number of unique sequences in cfastaFile2 = "+ ss2single);
		CfastaSequences UniqueSequenceArray = new CfastaSequences();
		UniqueSequenceArray.addAll(UniqueSequences);
		 
		}
	
	
	
	private ArrayList <Solid> removeDublettesSequences(ArrayList <Solid> otherSequences){
		
		 System.out.println("Sort sequences");
		 int ss1Pointer = 0;
		 int ss2Pointer = 0;
		 int ss1Size = this.size();
		 int ss2Size = otherSequences.size();
		 ArrayList <Solid> UniqueSequences = new ArrayList <Solid>(); 
		 int ss1single = 0;
		 int ss2single = 0;
		 
		 while(ss1Pointer < ss1Size && ss2Pointer < ss2Size){
			 
			while(ss1Pointer < ss1Size && this.get(ss1Pointer).compareTo(otherSequences.get(ss2Pointer)) < 0){
				UniqueSequences.add(this.get(ss1Pointer));
					 ss1Pointer++;
					 ss1single++;
			}
			while(ss2Pointer < ss2Size &&ss1Pointer < ss1Size  &&this.get(ss1Pointer).compareTo(otherSequences.get(ss2Pointer)) > 0){
				UniqueSequences.add(otherSequences.get(ss2Pointer));
				ss2Pointer++;
				ss2single++;
			}
			if(ss1Pointer < ss1Size && ss2Pointer < ss2Size && this.get(ss1Pointer).compareTo(otherSequences.get(ss2Pointer)) == 0){
				UniqueSequences.add(otherSequences.get(ss2Pointer));
				 ss2Pointer++;
				 ss1Pointer++;
			 }
		 }
		 while(ss1Pointer < ss1Size ){
			UniqueSequences.add(this.get(ss1Pointer));
			ss1Pointer++;
			ss1single++;
		}

		 while(ss2Pointer < ss2Size ){
				UniqueSequences.add(otherSequences.get(ss2Pointer));
				ss2Pointer++;
				ss2single++;
			}
		 
		 
		 System.out.println("Number of sequences in cfastaFile1 =        "+ ss1Size);
		 System.out.println("Number of sequences in cfastaFile2 =        "+ ss2Size);
		 System.out.println("Number of unique sequences in cfastaFile1 = "+ ss1single);
		 System.out.println("Number of unique sequences in cfastaFile2 = "+ ss2single);
		 return UniqueSequences;
		 
		}
	
	
	private static ArrayList <Solid> ANDSequences(ArrayList <Solid> ss1 ,ArrayList <Solid> ss2){
		
		 int ss1Pointer = 0;
		 int ss2Pointer = 0;
		 int ss1Size = ss1.size();
		 int ss2Size = ss2.size();
		 ArrayList <Solid> AndSequences = new ArrayList <Solid>(); 
		 int same = 0;
		 
		 while(ss1Pointer < ss1Size && ss2Pointer < ss2Size){
			 while(ss1Pointer < ss1Size && ss1.get(ss1Pointer).compareTo(ss2.get(ss2Pointer)) < 0){
				 ss1Pointer++;
			 }
			 while(ss2Pointer < ss2Size && ss1.get(ss1Pointer).compareTo(ss2.get(ss2Pointer)) > 0){
				 ss2Pointer++;
			 }
			 if(ss1Pointer < ss1Size && ss2Pointer < ss2Size && ss1.get(ss1Pointer).compareTo(ss2.get(ss2Pointer)) == 0){
				 AndSequences.add(ss1.get(ss1Pointer));
				 ss2Pointer++;
				 ss1Pointer++;
				 same++;
			 }
			 else
				 System.out.println("I should not end up here ever");
		 }
		 System.out.println("Number of sequences in cfastaFile1 =      "+ ss1Size);
		 System.out.println("Number of sequences in cfastaFile2 =      "+ ss2Size);
		 System.out.println("Number of sequences that where the same = "+ same);


		 return AndSequences;
	}

	
	
	public void setName(String name) {
		Name = name;
	}

	public String getName() {
		return Name;
	}

	


}




