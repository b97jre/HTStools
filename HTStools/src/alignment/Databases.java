package alignment;

import general.Functions;
import general.IOTools;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;

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

		String solidDir = Functions.getValue(T, "-solidDir", ".");
		String solidFile = Functions.getValue(T, "-solidFile", "sequences.cfasta");
		String resultDir = Functions.getValue(T, "-resultDir", ".");
		String resultFile = Functions.getValue(T, "-resultFile", ".");
		String DatabasesDir = Functions.getValue(T, "-DatabasesDir", resultDir);
		String DatabasesFile = Functions.getValue(T, "-DatabasesFile", resultDir);
		String outDir = Functions.getValue(T, "-outDir", resultDir);
		
		
		if(T.containsKey("-coverage")){
			String file = Functions.getValue(T, "-i", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);
			
			String GTF = Functions.getValue(T, "-g", ".");
			A.loadGTFFile(GTF);
			A.printCoverage();
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




