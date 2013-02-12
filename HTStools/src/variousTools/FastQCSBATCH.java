package variousTools;

import general.ExtendedReader;
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

public class FastQCSBATCH {

	String inDir;
	String outDir; 
	String time;
	String projectDir;
	boolean dependencies;

	public FastQCSBATCH(){
		projectDir = time = inDir = outDir = null;
		dependencies = false;
	}

	public static void main(String []args){

		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}
		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);
		FilterFastqSBATCH filter = new FilterFastqSBATCH();
		filter.run(T);
	}

	public void run(Hashtable<String,String> T){

		boolean allPresent = true;

		String timeStamp = Functions.getDateTime();
		SBATCHinfo sbatch = new SBATCHinfo();
		if(!sbatch.addSBATCHinfo(T)) allPresent = false;

		if(T.containsKey("-i"))
			inDir= Functions.getValue(T, "-i", ".");
		else{
			System.out.println("must contain inDirectory -i");
			allPresent = false;
		}
			projectDir= Functions.getValue(T, "-pDir", IOTools.getCurrentPath());

		if(T.containsKey("-t"))
			time = Functions.getValue(T, "-t", "15:00");
		else{
			System.out.println("must contain likely time -t now set to default 15 minutes");
		}

		if(allPresent)
			FastQC(sbatch, timeStamp);
		else
			System.out.println("\n\nAborting run because of missing arguments for fastQC.");
	}


	public void FastQC(SBATCHinfo sbatch, String timeStamp){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_fastQC.sh"));
			ArrayList <String> samples = IOTools.getDirectories(projectDir+"/"+inDir);
			for(int i = 0; i < samples.size(); i++){
				FastQCSample(EW,sbatch,samples.get(i), timeStamp,i);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	} 




	public void FastQCSample(ExtendedWriter generalSbatchScript, SBATCHinfo sbatch , String sample,String timestamp, int count){

		String finalInDir = projectDir+"/"+inDir+"/"+sample;

		String finalOutDir = finalInDir+"/FastQC";
		if(!IOTools.isDir(finalOutDir))
			IOTools.mkDir(finalOutDir);

		ArrayList <String> fileNames = IOTools.getSequenceFiles(finalInDir, "fastq");

		if(fileNames.isEmpty()){
			System.out.println("No .fastq files in folder :"+finalInDir);
			return;
		}
		else{
			if(!IOTools.isDir(finalOutDir+"/reports"))
				IOTools.mkDir(finalOutDir+"/reports");
			if(!IOTools.isDir(finalOutDir+"/scripts"))
				IOTools.mkDir(finalOutDir+"/scripts");

				for(int i = 0; i < fileNames.size(); i++){	
				String sbatchFileName = finalOutDir+"/scripts/"+timestamp+"_"+i+"_FastQC.sbatch";
				generalSbatchScript.println("sbatch "+ sbatchFileName);
				try{
					ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
					sbatch.printSBATCHinfoCore(EW,finalOutDir,timestamp,count, "FastQC", time);

					EW.println();
					EW.println("module load bioinfo-tools");
					EW.println("module load FastQC");
					EW.println();
					EW.println("cd "+finalInDir);
					EW.println("fastqc -o "+finalOutDir+" "+fileNames.get(i));
				
					EW.flush();
					EW.close();
				}
				catch(Exception E){E.printStackTrace();}
				
			}
		}

	}

	public void FastQCSample(ExtendedWriter EW, String inDir, String fileName){
		EW.println();
		EW.println("#####################################################################");
		EW.println("#FastQC of sample "+inDir+"/"+fileName);
		EW.println();
		if(!dependencies){
			EW.println("module load bioinfo-tools");
			EW.println("module load FastQC");
			EW.println();
			dependencies =true;
		}
		EW.println("cd "+inDir);
		EW.println("fastqc -o fastQC "+fileName);
		EW.println();
		EW.println("#FastQC stop");
		EW.println("#####################################################################");
		EW.println();
	}
	

}
