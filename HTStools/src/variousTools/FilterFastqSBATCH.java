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

public class FilterFastqSBATCH {

	String codeDir;
	String projectDir;
	String inDir;
	String outDir; 
	String appendix;
	String time;
	String[] sep; 

	public FilterFastqSBATCH(){
		projectDir =appendix =  codeDir = inDir = outDir  = time = null;
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

		
		this.sep = new String[2];
		sep[0] = "_1.fastq";
		sep[1] = "_2.fastq";
		

		boolean allPresent = true;

		String timeStamp = getDateTime();
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

		if(T.containsKey("-time"))
			time = Functions.getValue(T, "-time", ".");
		else{
			System.out.println("must contain likely time -time");
			allPresent = false;
		}

		appendix = Functions.getValue(T, "-a", "filtered");

		if(T.containsKey("-codeDir"))
			codeDir = Functions.getValue(T, "-codeDir", ".");
		else{
			System.out.println("must contain the directory of filter_fastq.pl -codeDir");
			allPresent = false;
		}

		if(allPresent)
			filterFasta(sbatch, timeStamp);
		else
			System.out.println("\n\nAborting run because of missing arguments for cutadapt.");
	}


	public void filterFasta(SBATCHinfo sbatch, String timeStamp){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_filter_fastq.sh"));
			ArrayList <String> samples = getDirectories(projectDir+"/"+inDir);
			for(int i = 0; i < samples.size(); i++){
				filter_fastqSample(EW,sbatch,samples.get(i), timeStamp,i);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	} 




	public void filter_fastqSample(ExtendedWriter generalSbatchScript, SBATCHinfo sbatch , String sample,String timestamp, int count2){

		String finalInDir = projectDir+"/"+inDir+"/"+sample;

		if(!IOTools.isDir(projectDir+"/"+outDir))
			IOTools.mkDir(projectDir+"/"+outDir);
		String finalOutDir = projectDir+"/"+outDir+"/"+sample;
		if(!IOTools.isDir(finalOutDir))
			IOTools.mkDir(finalOutDir);

		ArrayList <String> fileNames = getSequenceFiles(finalInDir);

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
				sbatch.printSBATCHinfo(EW,finalOutDir,timestamp,count2,"filter_fastq_"+sample, time);
				EW.println("cd "+finalInDir);
				
				ArrayList <String[]> pairs = IOTools.findPairs(fileNames,sep);
				
				for(int i = 0; i < pairs.size();i++){
					System.out.println(pairs.get(i)[0] + " " + pairs.get(i)[1]);
					EW.print(codeDir+"/filter_fastq.pl "+ appendix +" "+finalOutDir+" ");
					EW.println(pairs.get(i)[0] + " " + pairs.get(i)[1] +" &");
			
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



	private static String getDateTime() {
		DateFormat dateFormat = new SimpleDateFormat("yyyyMMdd_HH_mm_ss");
		Date date = new Date();
		return dateFormat.format(date);
	}



	public static ArrayList<String> getDirectories(String Dir){
		ArrayList<String> dirs = new ArrayList<String>();
		File dir = new File(Dir);

		String[] children = dir.list();
		File[] files = dir.listFiles();

		if (children == null) {
			return null;
			// Either dir does not exist or is not a directory
		} else {
			for (int i=0; i<children.length; i++) {
				// Get filename of file or directory
				if(files[i].isDirectory()){
					dirs.add(children[i]);
				}
			}
		}
		return dirs;
	}



	public static ArrayList<String> getSequenceFiles(String Dir){
		ArrayList<String> SequenceFiles = new ArrayList<String>();
		File dir = new File(Dir);

		String[] children = dir.list();
		if (children == null) {
			return null;
			// Either dir does not exist or is not a directory
		} else {
			for (int i=0; i<children.length; i++) {
				// Get filename of file or directory
				String filename = children[i];
				if(filename.indexOf(".fastq") > -1){
					SequenceFiles.add(filename);
				}
			}
		}
		return SequenceFiles;
	}

	public  static ArrayList<String> getAdapters(String file){

		ArrayList<String> SequenceFiles = new ArrayList<String>();
		if(file != null){
			try{
				ExtendedReader ER = new ExtendedReader(new FileReader(file));
				while(ER.more()){
					SequenceFiles.add(ER.readLine());
				}
			}catch (Exception E){E.printStackTrace();}
		}
		return SequenceFiles;
	}
}
