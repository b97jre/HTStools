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

	public FastQCSBATCH(){
		projectDir = time = inDir = outDir = null;
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

		if(T.containsKey("dataDir"))
			inDir= Functions.getValue(T, "dataDir", ".");
		else{
			System.out.println("must contain inDirectory -dataDir");
			allPresent = false;
		}
		if(T.containsKey("pDir"))
			projectDir= Functions.getValue(T, "pDir", ".");
		else{
			System.out.println("must contain projectDirectory -pDir");
			allPresent = false;
		}

		if(T.containsKey("time"))
			time = Functions.getValue(T, "time", ".");
		else{
			System.out.println("must contain likely time -time");
			allPresent = false;
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

			try{
				String sbatchFileName = finalOutDir+"/scripts/"+timestamp+".sbatch";
				generalSbatchScript.println("sbatch "+ sbatchFileName);

				ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
				sbatch.printSBATCHinfo(EW,finalOutDir,timestamp,count, "FastQC", time);

				EW.println();
				EW.println();
				EW.println("module load bioinfo-tools");
				EW.println("module load FastQC");
				EW.println("cd "+finalInDir);
				EW.println();
				EW.println();

				for(int i = 0; i < fileNames.size(); i++){	
					EW.println("fastqc -o "+finalOutDir+" "+fileNames.get(i)+" &");
					if((i+1)%8 == 0 ){
						EW.println();
						EW.println();
						EW.println("wait");
						EW.println();
						EW.println();
					}
				}
				EW.println();
				EW.println();
				EW.println("wait");
				EW.flush();
				EW.close();
			}catch(Exception E){E.printStackTrace();}
		}

	}


}
