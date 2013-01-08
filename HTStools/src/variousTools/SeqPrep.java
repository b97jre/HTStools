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

public class SeqPrep {

	String codeFile;
	String time;
	String[] sep; 
	String suffix;
	boolean hiseq;

	public SeqPrep(){
		 codeFile = time = null;
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

		String inDir = null;
		String outDir = null; 
		String projectDir = null;

		this.sep = new String[2];
		sep[0] = "1.fastq";
		sep[1] = "2.fastq";
		

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
		if(T.containsKey("-o"))
			outDir= Functions.getValue(T, "-o", ".");
		else{
			System.out.println("must contain outDirectory -o");
			allPresent = false;
		}
		if(T.containsKey("-pDir"))
			projectDir= Functions.getValue(T, "-pDir", ".");
		else{
			System.out.println("must contain projectDirectory -pDir");
			allPresent = false;
		}

		if(T.containsKey("-6"))
			hiseq = true;
		else{
			hiseq = false;
		}
	
		if(T.containsKey("-time"))
			time = Functions.getValue(T, "-time", ".");
		else{
			System.out.println("must contain likely time -time");
			allPresent = false;
		}

		suffix = Functions.getValue(T, "-suffix", "fastq");
		codeFile = Functions.getValue(T, "-codeFile", "/bubo/home/h17/johanr/bin/SeqPrep");
	
		if(allPresent)
			filterFasta(sbatch, timeStamp,projectDir,inDir,outDir);
		else
			System.out.println("\n\nAborting run because of missing arguments for SeqPrep.");
	}


	public void filterFasta(SBATCHinfo sbatch, String timeStamp,String projectDir, String inDir, String outDir){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			if(!IOTools.isDir(projectDir+"/"+outDir))
				IOTools.mkDir(projectDir+"/"+outDir);
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_filter_fastq.sh"));
			filter_fastqSample(EW,sbatch,projectDir+"/"+inDir ,projectDir+"/"+outDir, timeStamp,0);
				EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	} 




	public void filter_fastqSample(ExtendedWriter generalSbatchScript, SBATCHinfo sbatch , String inDir, String outDir,String timestamp, int count2){


		if(!IOTools.isDir(outDir))
			IOTools.mkDir(outDir);
		

		ArrayList <String> samples = IOTools.getDirectories(inDir);
		for(int i = 0; i < samples.size(); i++){
			filter_fastqSample(generalSbatchScript,sbatch,inDir+"/"+samples.get(i) ,outDir+"/"+samples.get(i), timestamp,i);
		}

		ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir,suffix);

		if(fileNames.isEmpty()){
			System.out.println("No "+suffix+" files in folder :"+inDir);
			return;
		}
		else{
			if(!IOTools.isDir(outDir+"/reports"))
				IOTools.mkDir(outDir+"/reports");
			if(!IOTools.isDir(outDir+"/scripts"))
				IOTools.mkDir(outDir+"/scripts");

			try{
				String sbatchFileName = outDir+"/scripts/"+timestamp+".sbatch";
				
				generalSbatchScript.println("sbatch "+ sbatchFileName);
				
				ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
				sbatch.printSBATCHinfo(EW,outDir,timestamp,count2,"filter_fastq_"+count2, time);
				EW.println("cd "+inDir);
				
				ArrayList <String[]> pairs = IOTools.findPairs(fileNames,sep);
				
				for(int i = 0; i < pairs.size();i++){
					System.out.println(pairs.get(i)[0] + " " + pairs.get(i)[1]);
					EW.print(codeFile +" -f "+pairs.get(i)[0]+" -r "+pairs.get(i)[1]);
					EW.print(" -1 "+outDir+"/"+pairs.get(i)[0]+".gz -2 "+outDir+"/"+pairs.get(i)[1]+".gz");
					if(hiseq)EW.print(" -6");
					EW.println(" -s "+ outDir+"/"+pairs.get(i)[2]+"merged.fastq.gz &");
			
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
