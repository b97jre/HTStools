package variousTools;

import general.ExtendedWriter;
import general.Functions;
import general.IOTools;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Hashtable;

public class BWA {

	String refFile;
	String time;
	String projectDir;
	String suffix;
	String split;
	String[] sep; 
	boolean files;
	String forward;
	String reverse;


	int missmatch;
	int seedLength; 
	int nrOfHits;
	int percentage;


	boolean strandSpecifik;
	boolean compress;


	public BWA(){
		this.refFile = "Database";
		this.split = ".";

		projectDir = time =  null;
	}

	public static void main(String []args){

		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}
		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);
		BWA BWA = new BWA();
		BWA.run(T);
	}

	public void run(Hashtable<String,String> T){

		String inDir, outDir, logDir;
		inDir = outDir = logDir = null;
		boolean allPresent = true;


		String timeStamp = Functions.getValue(T, "-TS", Functions.getDateTime());
		SBATCHinfo sbatch = new SBATCHinfo();
		if(!sbatch.addSBATCHinfo(T)){
			allPresent = false;
			return;
		}


		this.files = false;
		if(T.containsKey("-i"))
			inDir= Functions.getValue(T, "-i", ".");
		else if( T.containsKey("-f1") && T.containsKey("-f2")){
			forward= Functions.getValue(T, "-f1", ".");
			reverse= Functions.getValue(T, "-f2", ".");
			this.files=true;
		}
		else{
			System.out.println("must contain inDirectory -i or f1 and f2");
			allPresent = false;
		}


		if(T.containsKey("-o"))
			outDir= Functions.getValue(T, "-o", ".");
		else{
			System.out.println("must contain outDirectory -o");
			allPresent = false;
		}

		projectDir= Functions.getValue(T, "-pDir", IOTools.getCurrentPath());

		if(T.containsKey("-refFile")){
			refFile= Functions.getValue(T, "-refFile", "."); 
		}
		else{
			System.out.println("must contain referenceFile -refFile");
			allPresent = false;
		}


		if(T.containsKey("-time"))
			time = Functions.getValue(T, "-time", ".");
		else{
			System.out.println("must contain likely time -time");
			allPresent = false;
		}

		suffix = Functions.getValue(T,"-suffix","fastq");
		String seperator = Functions.getValue(T,"-sep","1."+suffix+" 2."+suffix);
		this.sep = seperator.split(" ");
		

		if(T.containsKey("-strandSpecific"))
			this.strandSpecifik=true;
		else
			this.strandSpecifik=false;
		
		compress = true;



		//		if(T.containsKey("-searchspace"))
		//			findOptimalValues(T,sbatch, timeStamp, outDir);
		//		else 
		if(T.containsKey("-build")){

		}
		else if(allPresent){
			if(!IOTools.isDir(inDir)){
				if(!IOTools.isDir(projectDir+"/"+inDir))
					inDir = projectDir+"/"+inDir;
				else {
					System.out.println("Neither "+ inDir +" nor "+projectDir+"/"+inDir+ "was found");
					return;
				}
			}
			if(!IOTools.fileExists(refFile)){
				if(!IOTools.fileExists(projectDir+"/"+refFile))
					refFile = projectDir+"/"+refFile;
				else {
					System.out.println("Neither "+ refFile +" nor "+projectDir+"/"+refFile+ "was found");
					return;
				}
			
			}
			else{
				if(refFile.indexOf("/")!= 0)
					refFile = IOTools.getCurrentPath()+"/"+refFile;
			}
			BWAinitiate(sbatch, timeStamp, inDir, projectDir+"/"+outDir);
		}
		else
			System.out.println("\n\nAborting run because of missing arguments for BWA.");
	}


	//
	//	public void findOptimalValues(Hashtable<String,String> T,SBATCHinfo sbatch, String timeStamp, String outDir){
	//		String forward,reverse;
	//		forward = reverse = null;
	//		boolean allPresent = true;
	//		if(T.containsKey("-1") && T.containsKey("-2")){
	//			forward = Functions.getValue(T, "-1", ".");
	//			reverse = Functions.getValue(T, "-2", ".");
	//		}
	//		else if(T.containsKey("-U") ){
	//			forward = Functions.getValue(T, "-U", ".");
	//		}
	//		else{
	//			System.out.println("must contain sequence file");
	//			allPresent = false;
	//		}
	//		if(allPresent){
	//			try{
	//				if(!IOTools.isDir(projectDir+"/scripts"))
	//					IOTools.mkDir(projectDir+"/scripts");
	//				ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_STAR.sh"));
	//				searchSpace(EW,sbatch, forward, reverse, outDir, timeStamp,8);
	//				EW.flush();
	//				EW.close();
	//			}catch(Exception E){E.printStackTrace();}
	//		}
	//	}
	//

	public void BWAinitiate(SBATCHinfo sbatch, String timeStamp, String inDir, String outDir){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+".BWA.sh"));
			if(files){
				BWAFile(EW,sbatch, timeStamp,  outDir,forward, reverse,forward);
			}
			else {
				BWADir(EW,sbatch, timeStamp, inDir, outDir);
			}

			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	} 




	public void BWADir( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String inDir, String outDir ){

		ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir,suffix);
		if(!fileNames.isEmpty()){
			if(!IOTools.isDir(outDir))
				IOTools.mkDirs(outDir);
			try{
				ArrayList <String[]> pairs = IOTools.findPairs(fileNames,sep);
				if(pairs.size()!= 0){
					for(int i = 0; i < pairs.size(); i++){
						String [] temp =  refFile.split("/");
						String refName = refFile;
						if(temp.length > 1)
							refName = temp[temp.length-1];

						String readsName = pairs.get(i)[0].substring(0,pairs.get(i)[0].indexOf(sep[0]));


						BWAFile(generalSbatchScript,sbatch,timestamp,outDir,inDir+"/"+pairs.get(i)[0], inDir+"/"+pairs.get(i)[1],readsName);
					}
				}
				else{
					System.out.println("Something wrong with the seperators? Asuming that these are single end reads");
					System.out.println(sep[0]);
					System.out.println(sep[1]);
					for(int i = 0; i < fileNames.size(); i++){
						System.out.println(fileNames.get(i));

						String [] temp =  refFile.split("/");
						String refName = refFile;
						if(temp.length > 1)
							refName = temp[temp.length-1];

						String fileName = fileNames.get(i).substring(0, fileNames.get(i).indexOf(this.suffix));

						String newOutDir = outDir+"/"+fileName+"_"+refName;
						BWAFile(generalSbatchScript,sbatch,timestamp,newOutDir,inDir+"/"+fileNames.get(i), null,refName);
					}
				}
			}catch(Exception E){E.printStackTrace();}
		}
		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
		for(int i = 0; i < subDirs.size(); i++){
			BWADir(generalSbatchScript,sbatch,timestamp,inDir+"/"+subDirs.get(i),outDir+"/"+subDirs.get(i));
		}



	}



	public void BWAFile( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp, String outDir,String forward,String reverse,String readsName ){

		if(!IOTools.isDir(outDir))
			IOTools.mkDirs(outDir);
		if(!IOTools.isDir(outDir+"/reports"))
			IOTools.mkDir(outDir+"/reports");
		if(!IOTools.isDir(outDir+"/scripts"))
			IOTools.mkDir(outDir+"/scripts");
		try{

			String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_"+readsName+"_BWA.sbatch";
			generalSbatchScript.println("sbatch "+ sbatchFileName);

			ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
			sbatch.printSBATCHinfo(EW,outDir,timestamp,0,"BWA", time);

			
			EW.println();
			EW.println("cd "+ outDir);

			BWACommand(EW, refFile, projectDir+"/"+forward, projectDir+"/"+reverse,8,this.strandSpecifik,this.suffix,readsName,compress);

			EW.println();
			EW.println();

			EW.flush();
			EW.close();

		}catch(Exception E){E.printStackTrace();}
	}
	
	
	public static void BWACommand( ExtendedWriter EW ,String refFile, String inFile1, String inFile2, int nrOfThreads, boolean strandSpecifik, String suffix,String nameSpace, boolean compress){
		
		EW.println("module load bioinfo-tools");
		EW.println("module load bwa/0.6.2");
		
		String refFileBase = refFile;
		if(refFileBase.indexOf("/") != -1)
			refFileBase = refFile.substring(refFile.lastIndexOf("/")+1);
		
		if(inFile2 != null){
			EW.println("bwa aln -t "+nrOfThreads+" "+ refFile+" "+ inFile1+"> "+refFileBase+"."+nameSpace+".1.sai");
			EW.println("bwa aln -t "+nrOfThreads+" "+ refFile+" "+ inFile2+"> "+refFileBase+"."+nameSpace+".2.sai");
			EW.println("bwa sampe "+ refFile+" "+refFileBase+"."+nameSpace+".1.sai "+refFileBase+"."+nameSpace+".2.sai "+ inFile1+" "+ inFile2+" > "+refFileBase+"."+nameSpace+".pe.sam");
			if(compress){
				SamtoolsSBATCH.sam2bam(EW, refFile+"."+nameSpace+".pe.sam",0,0,true,true,true,true,true);
			}
		}
		else{
			EW.println("bwa aln -t "+nrOfThreads+" "+ refFile+" "+ inFile1+"> "+refFileBase+"."+nameSpace+".sai");
			EW.println("bwa samse "+ refFile+" "+refFileBase+"."+nameSpace+".sai "+ inFile1+" > "+refFileBase+"."+nameSpace+".se.sam");
			if(compress){
				SamtoolsSBATCH.sam2bam(EW, refFileBase+"."+nameSpace+".se.sam",0,0,true,true,true,true,true);
			}
		}
		

		EW.println();
		EW.println("echo DONE");

	}

	

}

