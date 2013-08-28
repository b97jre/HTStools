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

public class STAR {

	String referenceDir;
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
	boolean interactive;
	boolean sam2bam;



	public STAR(){
		this.referenceDir = "Database";
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
		STAR STAR = new STAR();
		STAR.run(T);
	}

	public void run(Hashtable<String,String> T){

		String inDir, outDir, logDir;
		inDir = outDir = logDir = null;
		boolean allPresent = true;
		this.interactive = false;


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


		outDir= Functions.getValue(T, "-o", inDir+"_STAR");

		projectDir= Functions.getValue(T, "-pDir", IOTools.getCurrentPath());

		if(T.containsKey("-refDir")){
			referenceDir= Functions.getValue(T, "-refDir", "."); 
			if(referenceDir.lastIndexOf('/')==referenceDir.length()-1){

				referenceDir = referenceDir.substring(0,referenceDir.length()-1);
			}
		}
		else{
			System.out.println("must contain referenceDirectory -refDir");
			allPresent = false;
		}


		if(T.containsKey("-time"))
			time = Functions.getValue(T, "-time", ".");
		else if(T.containsKey("-interactive")){
			this.interactive=true;
		}else{
			System.out.println("must contain likely time -time or run as interactive -interactive");
			allPresent = false;
		}

		suffix = Functions.getValue(T,"-suffix","fastq");
		String seperator = Functions.getValue(T,"-sep","1."+suffix+" 2."+suffix);
		this.sep = seperator.split(" ");
		this.sam2bam = true;
		

		if(T.containsKey("-strandSpecific"))
			this.strandSpecifik=true;
		else
			this.strandSpecifik=false;



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
			if(referenceDir.indexOf('/') != 0){
				if(IOTools.isDir(projectDir+"/"+referenceDir)) referenceDir =  projectDir+"/"+referenceDir;
				else if(IOTools.isDir(IOTools.getCurrentPath()+"/"+referenceDir))  referenceDir =  IOTools.getCurrentPath()+"/"+referenceDir;
				else {
					System.out.println("Neither "+ referenceDir +" nor "+projectDir+"/"+referenceDir+ "nor "+ IOTools.getCurrentPath()+"/"+referenceDir+ " was found");
					return;
				}
			}
			STAR(sbatch, timeStamp, inDir, projectDir+"/"+outDir);
		}
		else
			System.out.println("\n\nAborting run because of missing arguments for STAR.");
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

	public void STAR(SBATCHinfo sbatch, String timeStamp, String inDir, String outDir){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+".STAR.sh"));
			if(interactive)
				STAR.STARCommandLoadGenome(EW, referenceDir);
			if(files){
				STARFile(EW,sbatch, timeStamp,  outDir,forward, reverse);
			}
			else {
				STARDir(EW,sbatch, timeStamp, inDir, outDir);
			}
			if(interactive)
				STAR.STARCommandRemoveGenome(EW, referenceDir);

			EW.flush();
			EW.close();
			
			System.out.println("Execute the following command to start all the runs:");
			System.out.println("sh "+projectDir+"/scripts/"+timeStamp+".STAR.sh >&"+projectDir+"/scripts/"+timeStamp+".STAR.sh.out ");
		}catch(Exception E){E.printStackTrace();}
	} 




	public void STARDir( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String inDir, String outDir ){

		ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir,suffix);
		if(!fileNames.isEmpty()){
			if(!IOTools.isDir(outDir))
				IOTools.mkDirs(outDir);
			try{
				ArrayList <String[]> pairs = IOTools.findPairs(fileNames,sep);
				if(pairs.size()!= 0){
					for(int i = 0; i < pairs.size(); i++){
						String [] temp =  referenceDir.split("/");
						String refName = referenceDir;
						if(temp.length > 1)
							refName = temp[temp.length-1];

						String readsName = pairs.get(i)[0].substring(0,pairs.get(i)[0].indexOf(sep[0]));

						String newOutDir = outDir+"/"+readsName+"_"+refName;

						STARFile(generalSbatchScript,sbatch,timestamp,newOutDir,inDir+"/"+pairs.get(i)[0], inDir+"/"+pairs.get(i)[1]);
					}
				}
				else{
					System.out.println("Something wrong with the seperators? Asuming that these are single end reads");
					//System.out.println(sep[0]);
					//System.out.println(sep[1]);
					for(int i = 0; i < fileNames.size(); i++){
						System.out.println(fileNames.get(i));

						String [] temp =  referenceDir.split("/");
						String refName = referenceDir;
						if(temp.length > 1)
							refName = temp[temp.length-1];

						String fileName = fileNames.get(i).substring(0, fileNames.get(i).indexOf(this.suffix));
						if(fileName.lastIndexOf(".") == fileName.length()) fileName.subSequence(0, fileName.length()-1);

						String newOutDir = outDir+"/"+fileName+"_"+refName;
						STARFile(generalSbatchScript,sbatch,timestamp,newOutDir,inDir+"/"+fileNames.get(i), null);
					}
				}
			}catch(Exception E){E.printStackTrace();}
		}
		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
		for(int i = 0; i < subDirs.size(); i++){
			STARDir(generalSbatchScript,sbatch,timestamp,inDir+"/"+subDirs.get(i),outDir+"/"+subDirs.get(i));
		}



	}



	public void STARFile( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp, String outDir,String forward,String reverse ){

		if(!IOTools.isDir(outDir))
			IOTools.mkDirs(outDir);
		if(!IOTools.isDir(outDir+"/reports"))
			IOTools.mkDir(outDir+"/reports");
		if(!IOTools.isDir(outDir+"/scripts"))
			IOTools.mkDir(outDir+"/scripts");
		try{

			String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_STAR.sbatch";
			if(!interactive)
				generalSbatchScript.println("sbatch "+ sbatchFileName);
			else
				generalSbatchScript.println("sh "+ sbatchFileName);
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
			sbatch.printSBATCHinfo(EW,outDir,timestamp,0,"STAR", time);

			EW.println();
			EW.println("cd "+ outDir);
			if(reverse != null)
				STARCommand(EW, referenceDir, projectDir+"/"+forward, projectDir+"/"+reverse,8,this.strandSpecifik,this.suffix);
			else
				STARCommand(EW, referenceDir, projectDir+"/"+forward, null,8,this.strandSpecifik,this.suffix);
			if(sam2bam){
				SamtoolsSBATCH.sam2bam(EW, outDir+"/Aligned.out.sam", -1, -1, true, true, true, true, true);
			}
			
			EW.println();
			EW.println();
			EW.println("wait");

			EW.flush();
			EW.close();

		}catch(Exception E){E.printStackTrace();}
	}
	
	
	public static void STARCommand( ExtendedWriter EW ,String refDir, String inFile1, String inFile2, int nrOfThreads, boolean strandSpecifik, String suffix){

		
		String bowtiecommand = "STAR ";
		bowtiecommand += " --genomeDir "+ refDir;
		bowtiecommand += " --readFilesIn "+ inFile1;
		if(inFile2 != null){
			bowtiecommand += " "+ inFile2;
		}
		bowtiecommand += " --runThreadN "+ nrOfThreads;
		bowtiecommand += " --genomeLoad LoadAndKeep";

		if(!strandSpecifik) bowtiecommand += " --outSAMstrandField intronMotif ";
		if(suffix.indexOf("gz")>0) bowtiecommand += " --readFilesCommand zcat ";
		if(suffix.indexOf("bz2")>0) bowtiecommand += " --readFilesCommand bzcat ";
		


		EW.println("echo START");
		EW.println();
		EW.println("echo \""+bowtiecommand+"\" 1>&2");
		EW.println(bowtiecommand);
		EW.println();
		EW.println("echo DONE");

	}
	public static void STARCommandLoadGenome( ExtendedWriter EW ,String refDir){

		
		String bowtiecommand = "STAR ";
		bowtiecommand += " --genomeDir "+ refDir;
		bowtiecommand += " --genomeLoad LoadAndExit";



		EW.println("echo START");
		EW.println();
		EW.println("echo \""+bowtiecommand+"\" 1>&2");
		EW.println(bowtiecommand);
		EW.println();
		EW.println("echo DONE");

	}
	
	public static void STARCommandRemoveGenome( ExtendedWriter EW ,String refDir){

		
		String bowtiecommand = "STAR ";
		bowtiecommand += " --genomeDir "+ refDir;
		bowtiecommand += " --genomeLoad remove";



		EW.println("echo START");
		EW.println();
		EW.println("echo \""+bowtiecommand+"\" 1>&2");
		EW.println(bowtiecommand);
		EW.println();
		EW.println("echo DONE");

	}

	

}

