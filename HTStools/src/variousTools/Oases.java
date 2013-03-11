package variousTools;

import general.ExtendedWriter;
import general.Functions;
import general.IOTools;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import Sequence.FastQSequences;

public class Oases {

	String time;
	String projectDir;

	String suffix;
	String split;
	String[] sep; 

	
	boolean pairedEnd;
	int insertLength;

	int maxKmer;
	int minKmer;
	int step;
	int mergeKmer;
	boolean strandSpecific;
	
	public Oases(){
		this.split = ".";
		projectDir = time =  null;
		
	}
	
	
	public void setValues(int maxKmer, int minKmer, int step,int mergeKmer, String suffix, int insertLength, String projectDir,boolean strandSpecific){
		this.maxKmer = maxKmer;
		this.minKmer = minKmer;
		this.step=step;
		this.mergeKmer=mergeKmer;
		this.suffix=suffix;
		this.insertLength =insertLength;
		this.projectDir = projectDir;
		this.strandSpecific = strandSpecific;
	}

	public static void main(String []args){

		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}
		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);
		Oases oases = new Oases();
		oases.run(T);
	}

	public void run(Hashtable<String,String> T){
		String inDir, outDir;
		inDir = outDir =  null;
		boolean allPresent = true;

		String timeStamp = Functions.getValue(T, "-TS",Functions.getDateTime());
		SBATCHinfo sbatch = new SBATCHinfo();
		if(!sbatch.addSBATCHinfo(T)) allPresent = false;

		if(T.containsKey("-i"))
			inDir= Functions.getValue(T, "-i", ".");
		else{
			System.out.println("must contain inDirectory -i");
			allPresent = false;
		}
		if(T.containsKey("-o"))
			outDir= Functions.getValue(T, "-o", ".");
		else{
			System.out.println("must contain outDirectory -o");
			allPresent = false;
		}

		
		this.pairedEnd = true;
		if(T.containsKey("-singleEnd")){
			this.pairedEnd = false;
		}


		projectDir= Functions.getValue(T, "-pDir", IOTools.getCurrentPath());


		time = Functions.getValue(T, "-time", "00.00.00");
		suffix = Functions.getValue(T,"-suffix","fastq");

		
		this.sep = new String[2];
		sep[0] = "1."+suffix;
		sep[1] = "2."+suffix;
		int insertSize = Integer.parseInt(Functions.getValue(T, "-insertSize", "400"));

		if(projectDir.compareTo("notPresent") != 0){
			inDir = projectDir+"/"+inDir;
			outDir = projectDir+"/"+outDir;
		}
		else{
			projectDir = inDir;
		}

	}



	public void oasesFilePE( ExtendedWriter EW , ArrayList <String[]> pairs, String outDir, String memory,int CPUs){
		// oases.pl --seqType fq --left 13_2/QC/13_2_1.fastq  --right 13_2/QC/13_2_2.fastq --output /bubo/proj/b2011098/private/oases/13_2/output  --min_contig_length 200 --CPU 10 --bflyHeapSpace 10000M
	
		
		
		String catFiles = "cat ";
		for(int i = 0; i < pairs.size()-1;i++)
			catFiles +=pairs.get(i)[0]+" ";
		catFiles +=pairs.get(pairs.size()-1)[0];
		catFiles +=" >"+outDir+"/left."+suffix;

		EW.println(catFiles);
		catFiles = "cat ";
		for(int i = 0; i < pairs.size()-1;i++)
			catFiles +=pairs.get(i)[1]+" ";
		catFiles +=pairs.get(pairs.size()-1)[1];
		catFiles +=" >"+outDir+"/right."+this.suffix;
		
		EW.println(catFiles);
		
		EW.println("cd "+outDir);
		
		String oasescommand = "oases_pipeline.py ";
				

		oasescommand += " -p \"-ins_length "+this.insertLength+"\" ";

		oasescommand +="-m "+this.minKmer+" ";
		oasescommand +="-M "+this.maxKmer+" ";
		oasescommand +="-s "+this.step+" ";
		oasescommand +="-g "+this.mergeKmer+" ";
		if(this.strandSpecific)oasescommand +=oasescommand += " -d \"-"+this.suffix+" -strand_specific -separate left."+this.suffix+" right."+this.suffix+" \"";
		else oasescommand += " -d \"-"+this.suffix+" -separate left."+this.suffix+" right."+this.suffix+" \"";
		
		EW.println(oasescommand);
		
		
	}

	public void oasesFileSR( ExtendedWriter EW , ArrayList <String> files, String outDir, String memory, int CPUs){
		// oases.pl --seqType fq --left 13_2/QC/13_2_1.fastq  --right 13_2/QC/13_2_2.fastq --output /bubo/proj/b2011098/private/oases/13_2/output  --min_contig_length 200 --CPU 10 --bflyHeapSpace 10000M

		String catFiles = "cat ";
		for(int i = 0; i < files.size()-1;i++)
			catFiles +=files.get(i)+" ";
		catFiles +=files.get(files.size()-1);

		catFiles +=" >"+outDir+"/single."+this.suffix;
		
		EW.println(catFiles);
		
		EW.println("cd "+outDir);
		
		String oasescommand = "python oases_pipeline.py ";
				

		oasescommand += " -p \"-ins_length "+this.insertLength+" \"";

		oasescommand +="-m "+this.minKmer+" ";
		oasescommand +="-M "+this.maxKmer+" ";
		oasescommand +="-s "+this.step+" ";
		oasescommand +="-g "+this.mergeKmer+" ";
		if(this.strandSpecific)oasescommand +=oasescommand += " -d \"-"+this.suffix+" -strand_specific single."+this.suffix+" \"";
		else oasescommand +=oasescommand += " -d \"-"+this.suffix+" single."+this.suffix+" \"";
		
		EW.println(oasescommand);
		
	}



	public static String getTime(int nrOfSequences){

		int nrofMB = nrOfSequences/1000000+1;
		int hours = nrofMB*3;
		if(hours > 168)hours = 168;	
		System.out.println("Number of sequences are:"+ nrOfSequences);
		String newTime = hours+":00:00";
		return newTime;

	}

	public static String FileStart(int nrOfSequences,String inDir, String outDir, String timestamp, int count, SBATCHinfo sbatch,ExtendedWriter EW , String time){

		int nrofMB = nrOfSequences/1000000+1;
		int hours = nrofMB*3;
		if(hours > 168)hours = 168;	
		System.out.println("Number of sequences are:"+ nrOfSequences);
		String newTime = "";
		if(time.compareTo("00.00.00") == 0){
			System.out.println("time not set will be set to "+hours +" hours");
			newTime =hours+":00:00";
		}
		else{
			newTime = time;
		}


		String[] split1 = outDir.split("/");
		String memory = "24G";
		if(nrofMB < 36){
			System.out.println("Memory allocated will be 36 GB");
			sbatch.printSBATCHinfoFat(EW,outDir,timestamp,count,"oases_"+split1[split1.length-1], newTime);
			memory = "34G";
		}
		else if(nrofMB < 72){
			System.out.println("Memory allocated will be 72 GB");
			sbatch.printSBATCHinfo72GB(EW,outDir,timestamp,count,"oases_"+split1[split1.length-1], newTime);
			memory = "70G";
		}				
		else{
			sbatch.printSBATCHinfohalvan(EW,outDir,timestamp,count,"oases_"+split1[split1.length-1], newTime,nrofMB);
			memory = "72G";
		}


		EW.println("ulimit -s unlimited");

		EW.println("cd "+inDir);


		EW.println();
		EW.println();

		return memory;


	}


//	public void dirPE( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String inDir, String outDir){
//
//		ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir,suffix);
//		if(fileNames.isEmpty()){
//			System.out.println("No "+suffix+" files in folder :"+inDir);
//			if(!IOTools.isDir(outDir))
//				IOTools.mkDir(outDir);
//		}
//		else{
//			if(!IOTools.isDir(outDir))
//				IOTools.mkDir(outDir);
//			if(!IOTools.isDir(outDir+"/reports"))
//				IOTools.mkDir(outDir+"/reports");
//			if(!IOTools.isDir(outDir+"/scripts"))
//				IOTools.mkDir(outDir+"/scripts");
//			int count = 0;
//			try{
//				ArrayList <String[]> pairs = IOTools.findPairs(fileNames,sep);
//				int nrOfSequences = FastQSequences.countSequencesPE(inDir,pairs);
//				String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_"+count+"_oases.sbatch";
//				generalSbatchScript.println("sbatch "+ sbatchFileName);
//
//				ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
//
//				String memory = FileStart(nrOfSequences,  inDir,  outDir,  timestamp,  count,  sbatch, EW,time);
//				int CPUs = 8;
//				oasesFilePE(EW, pairs,outDir, memory,CPUs);
//
//				EW.flush();
//				EW.close();
//
//			}catch(Exception E){E.printStackTrace();}
//
//			count++;
//
//		}
//		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
//		for(int i = 0; i < subDirs.size(); i++){
//			oasesDirPE(generalSbatchScript,sbatch,timestamp,inDir+"/"+subDirs.get(i),outDir+"/"+subDirs.get(i), insertSize);
//		}
//	}
//
//
//	public void oasesDirSR(ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String inDir, String outDir ){
//
//		ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir,suffix);
//		if(fileNames.isEmpty()){
//			System.out.println("No "+suffix+" files in folder :"+inDir);
//			if(!IOTools.isDir(outDir))
//				IOTools.mkDir(outDir);
//		}
//		else{
//			if(!IOTools.isDir(outDir))
//				IOTools.mkDir(outDir);
//			if(!IOTools.isDir(outDir+"/reports"))
//				IOTools.mkDir(outDir+"/reports");
//			if(!IOTools.isDir(outDir+"/scripts"))
//				IOTools.mkDir(outDir+"/scripts");
//			int count = 0;
//			try{
//				String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_"+count+"_oases.sbatch";
//				generalSbatchScript.println("sbatch "+ sbatchFileName);
//				ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
//
//				int nrOfSequences = FastQSequences.countSequencesSR(inDir,fileNames);
//
//				String memory = oasesFileStart(nrOfSequences,  inDir,  outDir,  timestamp,  count,  sbatch, EW,time );
//
//				EW.println();
//				EW.println();
//				oasesFileSR(EW, fileNames,outDir, memory,8);
//
//				EW.flush();
//				EW.close();
//
//			}catch(Exception E){E.printStackTrace();}
//
//			count++;
//
//		}
//		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
//		for(int i = 0; i < subDirs.size(); i++){
//			oasesDirSR(generalSbatchScript,sbatch,timestamp,inDir+"/"+subDirs.get(i),outDir+"/"+subDirs.get(i));
//		}
//	}
//
//	public void oasesDirSRSingle( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String inDir, String outDir ){
//
//		ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir,suffix);
//		if(fileNames.isEmpty()){
//			System.out.println("No "+suffix+" files in folder :"+inDir);
//			if(!IOTools.isDir(outDir))
//				IOTools.mkDir(outDir);
//		}
//		else{
//			if(!IOTools.isDir(outDir))
//				IOTools.mkDir(outDir);
//			if(!IOTools.isDir(outDir+"/reports"))
//				IOTools.mkDir(outDir+"/reports");
//			if(!IOTools.isDir(outDir+"/scripts"))
//				IOTools.mkDir(outDir+"/scripts");
//			int count = 0;
//			int number = 0;
//			int perRun = 100;
//			while(count < fileNames.size()/perRun+1){ 
//				try{
//					String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_"+count+"_oases.sbatch";
//					generalSbatchScript.println("sbatch "+ sbatchFileName);
//
//					ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
//
//					oasesFileStartSingle(inDir,  outDir,  timestamp,  count,  sbatch, EW);
//
//					for(int i = 0; i < perRun;i++){
//						if(number < fileNames.size()){
//							EW.println();
//							oasesFileSingle(EW, fileNames.get(number), inDir,  outDir);
//							number++;
//						}
//					}
//					EW.flush();
//					EW.close();
//
//				}catch(Exception E){E.printStackTrace();}
//
//				count++;
//			}
//		}
//	}
}
