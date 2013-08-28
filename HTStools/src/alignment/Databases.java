package alignment;

import general.Functions;
import general.IOTools;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;

import sun.org.mozilla.javascript.internal.EcmaError;

import com.sun.org.apache.xpath.internal.functions.Function;

import Sequence.CfastaSequences;
import Sequence.Solid;

import general.ExtendedReader;
import general.ExtendedWriter;




public class Databases implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//private String GenomeName;
	private String Name;
	protected ArrayList <Database> Databases2;
	protected ArrayList <String> samples;
	protected CfastaSequences solidSequences;

	public static void main(String[] args) {
		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}

		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);
		Databases.run(T);

	}




	Databases(String Name,String DatabasesDir,String DatabasesFile ){
		this.Name = Name;
		this.Databases2 = new ArrayList<Database>();
		//System.out.print("Reads database files .......");
		this.loadDatabasesFile(DatabasesDir, DatabasesFile);
		//System.out.println("finished");
	}	

	Databases(String DatabasesFile ){
		this.Name = "temp";
		this.Databases2 = new ArrayList<Database>();
		//System.out.print("Reads database files .......");
		this.loadDatabasesFile(DatabasesFile);
		//System.out.println("finished");
	}	

	public Databases(){
		this.Databases2 = new ArrayList<Database>();
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

		String dir = Functions.getValue(T, "-d", IOTools.getCurrentPath());
		String solidDir = Functions.getValue(T, "-solidDir", ".");
		String solidFile = Functions.getValue(T, "-solidFile", "sequences.cfasta");
		String resultDir = Functions.getValue(T, "-resultDir", ".");
		String resultFile = Functions.getValue(T, "-resultFile", "test.out.vcf");
		String DatabasesDir = Functions.getValue(T, "-DatabasesDir", resultDir);
		String DatabasesFile = Functions.getValue(T, "-DatabasesFile", resultDir);
		String outDir = Functions.getValue(T, "-outDir", resultDir);


		if(T.containsKey("-loadVCFinfo")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);

			String VCF = Functions.getValue(T, "-VCF");
			A.loadVCFFile(VCF);



			A.printVCFSample(outDir, VCF+".Sample.vcf",Functions.getValue(T, "-samples"),VCF);

		}

		if(T.containsKey("-printVCFinfoRfriendly")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);

			String VCF = Functions.getValue(T, "-VCF");
			A.loadVCFFile(VCF);



			A.printVCFSampleRfriendly(outDir, VCF+".Rfriendly",Functions.getValue(T, "-samples"));

		}



		if(T.containsKey("-addVCFinfo")){
			
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			
			file = dir+"/"+file;
			A.getDatabasesSizes(file);

			String VCF = Functions.getValue(T, "-VCF");
			A.loadVCFFile(VCF);

			String addVCFFiles= Functions.getValue(T, "-addVCFfiles");

			String[] files = addVCFFiles.split(",");

			for(int i = 0; i < files.length;i++){
				String[] sampleInfo = files[i].substring(files[i].lastIndexOf("/")+1).split("\\.");
				String sampleName = sampleInfo[0]+"_"+sampleInfo[1]+"_"+sampleInfo[2];
				A.addVCFInfo(files[i],sampleName);
			}
			
			
			A.printVCFSample(outDir, VCF+".Sample.vcf",Functions.getValue(T, "-samples"),VCF);

		}



		if(T.containsKey("-filterVCFinfo")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);

			String VCF = Functions.getValue(T, "-VCF");
			A.loadVCFFile(VCF);

			String BED = Functions.getValue(T, "-BED");
			A.loadFilterBedFile(BED);
			A.mergeFilters();
			
			
			if(!T.containsKey("-complement"))
				A.filterVCFfiles();
			else
				A.filterVCFfilesComplement();
				
			A.printVCFSample(outDir, VCF.substring(0,VCF.lastIndexOf('.'))+"."+BED.substring((BED.lastIndexOf("/")+1),BED.lastIndexOf('.'))+".vcf",Functions.getValue(T, "-samples"),VCF);

		}


	
		
		
		
		
		if(T.containsKey("-annotateSNPs")){
			String file = Functions.getValue(T, "-R", ".");
			String FullPathfile = dir+"/"+file;

			Databases  A= new Databases();
			A.getDatabasesSizes(FullPathfile);

			String VCF = Functions.getValue(T, "-VCF");
			String FullPathVCF = dir+"/"+VCF;
			A.loadVCFFile(FullPathVCF);

			String BED = Functions.getValue(T, "-BED");
			String FullPathBED = dir+"/"+BED;
			
			A.loadFilterBedFile(FullPathBED);
			A.mergeFilters();

			A.annotateVCFfiles(VCF.substring(0,VCF.lastIndexOf('.'))+"."+BED.substring(0,BED.lastIndexOf('.'))+".VCFannotaion");

			//A.printVCFSample(outDir, VCF.substring(0,VCF.lastIndexOf('.'))+"."+BED.substring(0,BED.lastIndexOf('.'))+".vcf",Functions.getValue(T, "-samples"));

		}




		if(T.containsKey("-onlyHeterozygous")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);

			String VCF = Functions.getValue(T, "-VCF");
			A.loadVCFFile(VCF);

			A.removeHomozygous(Functions.getValue(T, "-sample"));
			A.printVCFSample(outDir, VCF.substring(0,VCF.lastIndexOf('.'))+".heterozygous.vcf",null,VCF);

		}

		if(T.containsKey("-onlyHomozygous")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);

			String VCF = Functions.getValue(T, "-VCF");
			A.loadVCFFile(VCF);

			A.removeHeterozygous(Functions.getValue(T, "-sample"));
			A.printVCFSample(outDir, VCF.substring(0,VCF.lastIndexOf('.'))+".homozygous.vcf",null,VCF);

		}



		if(T.containsKey("-sortBed")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);
			String BED = Functions.getValue(T, "-BED");
			A.loadFilterBedFile(BED);
			A.printFilters();
			A.sortFilters();
			A.printFilters();

			A= new Databases();
			A.getDatabasesSizes(file);
			BED = Functions.getValue(T, "-BED");
			A.loadFilterBedFile(BED);
			A.printFilters();
			A.mergeFilters();
			A.printFilters();


		}



		if(T.containsKey("-loadGFF")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);

			String GTF = Functions.getValue(T, "-GFF", ".");
			A.loadGTFFile(GTF);
			A.printCoverage();
		}

		if(T.containsKey("-coverage")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);

			String GTF = Functions.getValue(T, "-GFF", ".");
			A.loadGTFFile(GTF);
			A.printCoverage();

			String VCF = Functions.getValue(T, "-VCF");
			A.loadVCFFile(VCF);
		}


		if(T.containsKey("-printDistribution")){
			Databases B = new Databases("temp", DatabasesDir, DatabasesFile);

			int cutoff = Integer.parseInt(Functions.getValue(T, "-x", "1"));
			B.loadRmapperSequences(resultDir, resultFile);
			B.mapSolidSequences();
			B.printDistribution(cutoff, outDir,resultFile);

		}



		//runPrograms(DatabasesDir, DatabasesFile,solidDir, solidFile,resultDir, resultFile,outDir,
		//			"temp",0,0,T);
	}



	public static void runPrograms(String DatabasesDir, String DatabasesFile, 
			String solidDir, String solidFile, String resultDir, String resultFile, String outDir,
			String experiment, int nucleotideLength,int nrOfHits,
			Hashtable<String,String> T){

		int cutoff = Integer.parseInt(Functions.getValue(T, "-cutoff", "100"));

		Databases B = new Databases(experiment, DatabasesDir, DatabasesFile);
		if(T.containsKey("-solidDir")){
			B.loadSolidSequences(solidDir, solidFile);
		}
		if(T.containsKey("-resultDir")){
			if(!experiment.contains("temp")){
				resultDir = resultDir+"/"+experiment;
				if(nucleotideLength > 0)
					resultDir = resultDir+"/"+nucleotideLength;
				if(nrOfHits > 0){
					resultFile = experiment+"."+DatabasesFile+".rmapper."+nrOfHits;
					if(T.containsKey("-extra"))
						resultFile += "."+ Functions.getValue(T, "-extra", "rmapper");
				}
			}
			B.loadRmapperSequences(resultDir, resultFile);
		}

		B.mapSolidSequences();

		if(T.containsKey("-countGroups")){
			B.solidSequences.countGroups();
		}

		if(T.containsKey("-countlocations")){
			int total = B.countLocations();
			System.out.println(experiment+"\t"+nucleotideLength+"\t"+nrOfHits+"\t"+B.getNrOfSequences()+"\t"+B.getNrOfHits2()+"\t"+total);
		}

		if(T.containsKey("-printDistribution")){
			B.printDistribution(cutoff, outDir,resultFile);
		}

		if(T.containsKey("-printMaxSequences")){
			int length = Integer.parseInt(Functions.getValue(T, "-length", "21"));
			B.printMaxSequences(outDir,resultFile,length, cutoff);
		}

	}


	private void loadSolidSequences(String dir, String file){
		this.solidSequences = new CfastaSequences();
		this.solidSequences.addSolidSequences(dir,file);
	}

	private void loadRmapperSequences(String dir, String file){
		//System.out.print("reading rmapper files.......");
		this.solidSequences = new CfastaSequences();
		this.solidSequences.addRmapperSequences(dir,file);
		//System.out.println("finished");
		//this.solidSequences.countHits();
	}

	private void removeNonredundantRmapperSequences(String dir, String file){
		this.solidSequences = new CfastaSequences();
		this.solidSequences.addRmapperSequences(dir,file);




	}


	public static int[] countGroups(String experiment, int SequenceLength, int nrOfHits, Hashtable<String,String> T){

		String resultDir = Functions.getValue(T, "-resultDir", ".");
		String resultFile = Functions.getValue(T, "-resultFile", ".");
		String DatabasesDir = Functions.getValue(T, "-DatabasesDir", resultDir);
		String DatabasesFile = Functions.getValue(T, "-DatabasesFile", resultDir);

		String[] chromosomes = null; 
		String finalDir = resultDir+"/"+experiment+"/"+SequenceLength+"/";

		String finalFile = experiment+"."+DatabasesFile+".rmapper."+nrOfHits; 
		String extra = Functions.getValue(T, "-extra", "");

		if(extra.length() > 1)
			finalFile+="."+extra;
		Databases C = new Databases("temp",DatabasesDir, DatabasesFile);
		C.loadRmapperSequences(finalDir, finalFile);
		C.mapSolidSequences();
		return C.solidSequences.countGroups();
	}




	private void mapSolidSequences(){
		for(int  i = 0 ; i < this.solidSequences.size(); i++){
			Solid hit = this.solidSequences.get(i);
			for(int j = 0; j < this.Databases2.size(); j++){
				this.Databases2.get(j).mapSolidSequence2Database(hit);
			}
		}
	}



	private int getNrOfHits(){
		int total = 0;
		for(int  i = 0; i < this.Databases2.size(); i++){
			total += this.Databases2.get(i).getNrOfHits();
		}
		return total;
	}


	private int getNrOfHits2(){
		int total = 0;
		return this.solidSequences.getNrOfHits();
	}

	private int getNrOfSequences(){
		return this.solidSequences.size();
	}


	private int countLocations(){
		int total = 0;
		for(int  i = 0; i < this.Databases2.size(); i++){
			total += this.Databases2.get(i).countLocations();
		}
		return total;
	}



	private void printDistribution(int cutoff,String outDir,String resultFile){

		int total = 0;
		try{
			if(!IOTools.isDir(outDir))
				IOTools.mkDir(outDir);
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+resultFile+".distribution.overview"));
			for(int  i = 0; i < this.Databases2.size(); i++){
				this.Databases2.get(i).printOverallDistribution(EW,cutoff);
			}
			EW.flush();
			EW.close();

			for(int  i = 0; i < this.Databases2.size(); i++){
				String Name = IOTools.fixFileName(this.Databases2.get(i).getName());
				EW = new ExtendedWriter(new FileWriter(outDir+"/"+Name+"."+resultFile+".distribution"));
				this.Databases2.get(i).printDistribution(EW,cutoff);
				int length = this.Databases2.get(i).getLength();
				if(Name.indexOf("miRNAstart=") > -1){
					//int A = Integer.parseInt(Name.substring(Name.indexOf("miRNAstart=")+11));
					if(Name.indexOf(",") != -1)
						System.out.println("plot.miRNAdist(dir,experiments,\""+ Name+".distribution\",21,0,"+length+",dir,\""+Name.substring(0,Name.indexOf(","))+"\")");
					else{
						System.out.println("plot.miRNAdist(dir,experiments,\""+ Name+".distribution\",21,0,"+length+",dir,\""+Name.substring(0,Name.indexOf("_"))+"\")");
					}

				}
				else
					System.out.println("plot.Distribution.infestans(dir,experiments,\""+ Name+"\","+length+",dir,legendInfo)");


				EW.flush();
				EW.close();
			}	
		}
		catch(Exception E){E.printStackTrace();}

	}



	private void printVCFSampleRfriendly(String outDir,String resultFile,String sampleString){

		int total = 0;
		try{
			if(!IOTools.isDir(outDir))
				IOTools.mkDir(outDir);
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+resultFile));
			// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
			EW.print("CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
			if(sampleString != null){
				this.samples = new ArrayList<String>();
				String[] sampleArray = sampleString.split(","); 
				for(int i = 0; i < sampleArray.length; i++){
					samples.add(sampleArray[i]);
				}
			}
			for(int i  = 0 ; i < this.samples.size(); i++){
				EW.print("\t"+this.samples.get(i)+"_Call\t"+this.samples.get(i)+"_Count1\t"+this.samples.get(i)+"_Count2\t"+this.samples.get(i)+"_Total\t"+this.samples.get(i)+"_fraction");
			}	
			EW.println();
			for(int  i = 0; i < this.Databases2.size(); i++){
				this.Databases2.get(i).printVCFinfoSamplesRfriendly(this.samples, EW);
			}
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}

	}


	private void printHeader(ExtendedWriter EW, String inFile){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(inFile));
			while(ER.more() && ER.lookAhead() == '#'){
				EW.println(ER.readLine());
			}
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}
	
	
	private void printVCFSample(String outDir,String resultFile,String sampleString, String inFile){

		try{
			if(!IOTools.isDir(outDir))
				IOTools.mkDir(outDir);
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+resultFile));
			printHeader(EW,inFile);
			
			// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
			EW.print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
			if(sampleString != null){
				this.samples = new ArrayList<String>();
				String[] sampleArray = sampleString.split(","); 
				for(int i = 0; i < sampleArray.length; i++){
					samples.add(sampleArray[i]);
				}
			}
			for(int i  = 0 ; i < this.samples.size(); i++){
				EW.print("\t"+this.samples.get(i));
			}
			EW.println();
			for(int  i = 0; i < this.Databases2.size(); i++){
				this.Databases2.get(i).printVCFinfoSamples(this.samples, EW);
			}
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}

	}

	private void printFilters(){
		for(int  i = 0; i < this.Databases2.size(); i++){
			this.Databases2.get(i).printFilters();
		}

	}

	private void sortFilters(){
		int total = 0;
		for(int  i = 0; i < this.Databases2.size(); i++){
			this.Databases2.get(i).sortBEDfilters();
		}

	}

	private void mergeFilters(){
		int total = 0;
		for(int  i = 0; i < this.Databases2.size(); i++){
			this.Databases2.get(i).mergeBEDfilters();
		}

	}

	private void filterVCFfiles(){
		for(int  i = 0; i < this.Databases2.size(); i++){
			this.Databases2.get(i).filterVCFinfoOutside();
		}
	}

	
	private void filterVCFfilesComplement(){
		for(int  i = 0; i < this.Databases2.size(); i++){
			this.Databases2.get(i).filterVCFinfoInside();
		}
	}
	
	
	private void annotateVCFfiles(String VCFfile){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(VCFfile));
			for(int  i = 0; i < this.Databases2.size(); i++){
				this.Databases2.get(i).annotateVCFinfo(EW);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){
			E.printStackTrace();
		}
	}




	private void removeHomozygous(String sampleString){
		ArrayList<String> samples2 =null;
		if(sampleString != null){
			samples2 = new ArrayList<String>();
			String[] sampleArray = sampleString.split(","); 
			for(int i = 0; i < sampleArray.length; i++){
				samples2.add(sampleArray[i]);
			}
		}else{
			samples2 = this.samples;
		}
		for(int  i = 0; i < this.Databases2.size(); i++){
			this.Databases2.get(i).removeHomozygous(samples2);
		}

	}

	private void removeHeterozygous(String sampleString){
		ArrayList<String> samples2 =null;
		if(sampleString != null){
			samples2 = new ArrayList<String>();
			String[] sampleArray = sampleString.split(","); 
			for(int i = 0; i < sampleArray.length; i++){
				samples2.add(sampleArray[i]);
			}
		}else{
			samples2 = this.samples;
		}
		for(int  i = 0; i < this.Databases2.size(); i++){
			this.Databases2.get(i).removeHeterozygous(samples2);
		}

	}



	//
	//	private void printVCFSampleExons(String outDir,String resultFile,String sample,){
	//
	//		int total = 0;
	//		try{
	//			if(!IOTools.isDir(outDir))
	//				IOTools.mkDir(outDir);
	//			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+resultFile));
	//			for(int  i = 0; i < this.Databases2.size(); i++){
	//				this.Databases2.get(i).printVCFinfoSample(sample, EW);
	//			}
	//			EW.flush();
	//			EW.close();
	//		}
	//		catch(Exception E){E.printStackTrace();}
	//
	//	}
	//


	private void printVCFDistribuionSample(String outDir,String resultFile,String sample){

		int total = 0;
		try{
			if(!IOTools.isDir(outDir))
				IOTools.mkDir(outDir);
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+resultFile));
			ArrayList<String> samples = new ArrayList<String>();
			for(int  i = 0; i < this.Databases2.size(); i++){
				this.Databases2.get(i).printVCFinfoSamples(samples, EW);
			}
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}

	}



	private void printMaxSequences(String outDir,String resultFile, int length, int cutoff){

		try{
			if(!IOTools.isDir(outDir))
				IOTools.mkDir(outDir);
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+resultFile+".maxSequences.overview"));
			for(int  i = 0; i < this.Databases2.size(); i++){
				this.Databases2.get(i).printMaxSequence(EW, length, cutoff);
			}
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}

	}


	public void compareDifferentRuns(Databases otherRun, ExtendedWriter ER,String exp1, String exp2){
		int nrOfHits = this.getNrOfHits();
		System.out.println("Number of hits in first genome : "+nrOfHits);
		int otherNrOfHits = otherRun.getNrOfHits();
		System.out.println("Number of hits in second genome: "+otherNrOfHits);

		double ratio = (double)nrOfHits/(double)otherNrOfHits;
		System.out.println("The ratio is: "+ratio);
		ER.println("Name,chromosome,ratio,"+exp1+","+exp2+",location,nrOfHits in "+ exp1+" = "+nrOfHits+",nrOfHits in "+ exp2+" = "+otherNrOfHits);

		for(int  i = 0; i < this.Databases2.size(); i++){
			this.Databases2.get(i).compareDistribution(otherRun.Databases2.get(i),ER);
		}
	}

	public void loadDatabasesFile(String dir, String fileName){
		loadDatabasesFile(dir+"/"+fileName);
	}


	public void loadDatabasesFile( String fileName){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(fileName));
			int count = 0; 
			int progress = 1000;
			while(ER.more()){
				while(ER.lookAhead() != '>')ER.skipLine();
				String DatabaseName = ER.readLine().substring(1);
				this.Databases2.add(new Database(DatabaseName,ER));
				count++;
				if(count > progress){
					System.out.println("sequences read :"+ progress);
					progress += 1000;
				}
			}
		}catch(Exception E){E.printStackTrace();}
	}

	public void getDatabasesSizes(String fileName){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(fileName));
			int count = 0; 
			int progress = 1000;
			while(ER.more()){
				while(ER.lookAhead() != '>')ER.skipLine();
				String DatabaseName = ER.readLine().substring(1);
				Database A = new Database(DatabaseName);
				A.getChromosomeSequenceSize(ER);
				this.Databases2.add(A);
				count++;
				if(count > progress){
					System.out.println("sequences read :"+ progress);
					progress += 1000;
				}
			}
		}catch(Exception E){E.printStackTrace();}
	}


	public void loadGTFFile(String dir, String fileName){
		loadGTFFile(dir+"/"+fileName);
	}

	public void loadGTFFile(String fileName){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(fileName));
			int count = 0;
			int progress = 1000;
			while(ER.more()){
				count++;
				readGTFInfo(ER);
				if(count > progress ){
					progress += 1000;
					System.out.println("gtf transcripts read: "+progress );
				}
			}
		}catch(Exception E){E.printStackTrace();}
	}


	private void readGTFInfo(ExtendedReader ER){
		String GFF3line = ER.readLine();
		String[] columns = GFF3line.split("\t");
		if(columns.length > 2){
			if(columns[2].indexOf("transcript") == 0){
				addCoverage(columns);
			}
		}
		else
			System.out.println(GFF3line);
	}


	public void loadVCFFile(String fileName){
		//	0			1		2     3		  4			5		6		7		8		9				10			11			12				13				14				15				16		 
		// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
		// scaffold_1      5       .       A       T       56.41   .       AC=3;AF=0.188;AN=16;BaseQRankSum=-1.804;DP=72;Dels=0.00;FS=0.000;HaplotypeScore=2.4206;MLEAC=2;MLEAF=0.125;MQ=26.05;MQ0=11;MQRankSum=0.618;QD=2.69;ReadPosRankSum=1.053 GT:AD:DP:GQ:PL  0|0:14,0:14:30:0,30,268 0|0:1,0:1:3:0,3,29      0/1:13,5:18:72:72,0,247 0|0:14,0:14:39:0,39,316 0|0:13,0:13:39:0,39,328 0|0:4,0:4:3:0,3,25      1|1:2,1:3:3:23,3,0      0|0:5,0:5:3:0,3,25


		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(fileName));
			int count = 0;
			int DatabasePointer = 0;
			int progress = 1000;
			String info = ER.readLine();
			while(ER.more() && ER.lookAhead() == '#'){
				info = ER.readLine();
			}
			String [] infoArray = info.split("\t");
			this.samples = new ArrayList<String>();
			for(int i = 9 ; i < infoArray.length;i++){
				samples.add(infoArray[i]);
				System.out.println(infoArray[i]);
			}
			System.out.print("Now gathering vcf data from :"+this.Databases2.get(DatabasePointer).getName());
			while(ER.more()){
				String[] LineInfo = ER.readLine().split("\t");
				while(DatabasePointer < this.Databases2.size() && this.Databases2.get(DatabasePointer).getName().compareTo(LineInfo[0])!=0){
					System.out.println("Finished");
					DatabasePointer++;
					System.out.print("Now gathering vcf data for :"+this.Databases2.get(DatabasePointer).getName());
				}
				//				if(DatabasePointer > this.Databases2.size())System.out.println(LineInfo[0]+"\t"+LineInfo[1]);
				//				else{
				this.Databases2.get(DatabasePointer).addVCFinfo(samples, LineInfo);
				//				}
				count++;
				if(count > progress ){
					progress += 10000;
					System.out.print(".");
				}
			}
			System.out.println("Finished");
		}catch(Exception E){E.printStackTrace();}
	}



	public void addVCFInfo(String fileName,String sampleName){
		//	0			1		2     3		  4			5		6		7		8		9				10			11			12				13				14				15				16		 
		// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
		// scaffold_1      5       .       A       T       56.41   .       AC=3;AF=0.188;AN=16;BaseQRankSum=-1.804;DP=72;Dels=0.00;FS=0.000;HaplotypeScore=2.4206;MLEAC=2;MLEAF=0.125;MQ=26.05;MQ0=11;MQRankSum=0.618;QD=2.69;ReadPosRankSum=1.053 GT:AD:DP:GQ:PL  0|0:14,0:14:30:0,30,268 0|0:1,0:1:3:0,3,29      0/1:13,5:18:72:72,0,247 0|0:14,0:14:39:0,39,316 0|0:13,0:13:39:0,39,328 0|0:4,0:4:3:0,3,25      1|1:2,1:3:3:23,3,0      0|0:5,0:5:3:0,3,25


		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(fileName));
			int count = 0;
			int DatabasePointer = 0;
			int progress = 1000;
			String info = ER.readLine();
			while(ER.more() && ER.lookAhead() == '#'){
				info = ER.readLine();
			}
			String [] infoArray = info.split("\t");

			ArrayList<String> newSamples = new ArrayList<String>();
			for(int i = 9 ; i < infoArray.length;i++){
				while(Functions.contains(samples, sampleName)) sampleName = sampleName+"_1";
				samples.add(sampleName);
				newSamples.add(sampleName);
				System.out.println(sampleName);
			}
			infoArray = ER.readLine().split("\t");
			while(ER.more()){

				System.out.println("Now gathering vcf data from :"+this.Databases2.get(DatabasePointer).getName());
				infoArray = this.Databases2.get(DatabasePointer).addVCFSamples(ER, newSamples,infoArray);
				DatabasePointer++;
				//				}
				count++;
				if(count > progress ){
					progress += 10000;
					System.out.print(".");
				}
			}
			System.out.println("Finished");
		}catch(Exception E){E.printStackTrace();}
	}


	public int findPointer(String name){
		int DatabasePointer = 0;
		while(DatabasePointer < this.Databases2.size() && this.Databases2.get(DatabasePointer).getName().compareTo(name) != 0){
			DatabasePointer++;
		}
		if(DatabasePointer == this.Databases2.size()){
			System.out.println("Could not find "+name);
			return -1;
		}
		return DatabasePointer;



	} 

	public void loadFilterBedFile(String fileName){
		//	0			1		2     3		  4			5		6		7		8		9			 
		// #CHROM(0)  Start     Stop     Name     .     .    INFO  type    Something  XTR 	AInfo   
		// scaffold_1      767     2124    PAC:20891551.exon.3     .       -       phytozome8_0    exon    .       ID=PAC:20891551.exon.3;Parent=PAC:20891551;pacid=20891551
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(fileName));
			String currentDatabase = "";
			int count = 0;
			int DatabasePointer = 0;
			int progress = 1000;
			while(ER.more() && ER.lookAhead() == '#'){
				ER.readLine();
			}
			String[] LineInfo = ER.readLine().split("\t");
			currentDatabase = LineInfo[0];
			DatabasePointer = findPointer(currentDatabase);
			System.out.print("Now gathering filter data for :"+this.Databases2.get(DatabasePointer).getName());

			while(ER.more()){
				if(currentDatabase.compareTo(LineInfo[0])!=0){
					System.out.println("Finished");
					currentDatabase = LineInfo[0];
					
					DatabasePointer = findPointer(currentDatabase);
					System.out.println(currentDatabase+"\t=  "+ this.Databases2.get(DatabasePointer).getName());
					System.out.print("Now gathering filter data for :"+this.Databases2.get(DatabasePointer).getName());						
				}
				//				if(DatabasePointer > this.Databases2.size())System.out.println(LineInfo[0]+"\t"+LineInfo[1]);
				//				else{
				if(DatabasePointer!= -1)
					this.Databases2.get(DatabasePointer).addBEDfilterInfo(LineInfo);
				else{
					System.out.print("No data for  :"+currentDatabase);
					
				}
				LineInfo = ER.readLine().split("\t");
				//				}
				count++;
				if(count > progress ){
					progress += 10000;
					System.out.print(".");
				}
			}
			
			if(currentDatabase.compareTo(LineInfo[0])!=0){
				System.out.println("Finished");
				currentDatabase = LineInfo[0];
				DatabasePointer = findPointer(currentDatabase);
				System.out.print("Now gathering filter data for :"+this.Databases2.get(DatabasePointer).getName());						
			}

			this.Databases2.get(DatabasePointer).addBEDfilterInfo(LineInfo);
			System.out.println("Finished");
		}catch(Exception E){E.printStackTrace();}
	}




	private void addCoverage(String[] columns){
		int databaseNr = findDatabase(columns[0]);
		if(databaseNr > -1)this.Databases2.get(databaseNr).addCoverage(Integer.parseInt(columns[3]), Integer.parseInt(columns[4]), columns[8]);
		else System.out.println(columns[0] +" not found in Database");

	}
	private int findDatabase(String name){
		int pointer = -1;
		for(int i = 0; i < this.Databases2.size(); i++){
			if(this.Databases2.get(i).isDatabase(name)){
				pointer = i;
				i = this.Databases2.size();
			}
		}
		return pointer;
	}
	private void printCoverage(){
		int pointer = -1;
		for(int i = 0; i < this.Databases2.size(); i++){
			this.Databases2.get(i).printCoverage();
		}
	}



}




