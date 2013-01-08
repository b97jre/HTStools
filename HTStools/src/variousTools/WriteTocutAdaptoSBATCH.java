package variousTools;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Hashtable;

import general.ExtendedReader;
import general.ExtendedWriter;
import general.Functions;
import general.IOTools;

public class WriteTocutAdaptoSBATCH {



	String projectDir;
	String inDir;
	String outDir; 
	String projectNumber;
	String time;
	String threePrimeAdaptersFile;
	String allAdaptersFile;
	String seqType;
	int cutoff;
	boolean cutadapt;
	boolean QC;
	String suffix;
	String[] sep;
	int length;


	public static void main(String []args){



		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}
		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);

		WriteTocutAdaptoSBATCH sbatchscript = new WriteTocutAdaptoSBATCH();
		sbatchscript.run(T);


	}





	public WriteTocutAdaptoSBATCH(){
		projectDir = inDir = outDir = projectNumber = time = threePrimeAdaptersFile = allAdaptersFile = null;
	}

	public void run(Hashtable<String,String> T){

		boolean allPresent = true;
		cutadapt = true;
		QC = true;

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
		outDir= Functions.getValue(T, "-o", "cutadapt");

		if(T.containsKey("-time"))
			time = Functions.getValue(T, "-time", ".");
		else{
			System.out.println("must contain likely time -time");
			allPresent = false;
		}

		if(T.containsKey("-3Adapter") || T.containsKey("-Adapter")){
			if(T.containsKey("-3Adapter"))
				threePrimeAdaptersFile= Functions.getValue(T, "-3Adapter", ".");
			if( T.containsKey("-Adapter"))
				allAdaptersFile= Functions.getValue(T, "-Adapter", ".");
		}
		else{
			System.out.println("must contain projectNumber -3Adapter or -Adapter now only QC of sequences will be carried out.");
			cutadapt = false;
		}

		seqType = Functions.getValue(T, "-f", ".fastq");
		cutoff = Integer.parseInt(Functions.getValue(T, "-q", "-1"));
		length = Integer.parseInt(Functions.getValue(T, "-l", "40"));
		suffix = Functions.getValue(T,"-suffix","fastq");
		this.sep = new String[2];
		sep[0] = "1."+suffix;
		sep[1] = "2."+suffix;

		if(allPresent)
			cutAdapt(sbatch, projectDir,inDir, outDir, projectNumber, time, threePrimeAdaptersFile,allAdaptersFile,timeStamp,length);
		else
			System.out.println("\n\nAborting run because of missing arguments for cutadapt.");


	}


	public void cutAdapt(SBATCHinfo sbatch,String projectDir,String inDir, String outDir,String projectNumber, String hours, String threeAdaptersFile, String otherAdaptersFile,String timestamp,int length){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			if(!IOTools.isDir(projectDir+"/"+outDir))
				IOTools.mkDir(projectDir+"/"+outDir);
			if(!IOTools.isDir(projectDir+"/"+outDir+"/QC"))
				IOTools.mkDir(projectDir+"/"+outDir+"/QC");
			ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir, seqType);
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timestamp+"cutadapt.sh"));
			if(!fileNames.isEmpty()){
				if(cutadapt)cutAdaptSample(sbatch, EW,projectDir+"/"+inDir,projectDir+"/"+outDir, hours, threeAdaptersFile, otherAdaptersFile, timestamp);
			}
			ArrayList <String> samples = IOTools.getDirectories(projectDir+"/"+inDir);
			for(int i = 0; i < samples.size(); i++){
				if(cutadapt)cutAdaptSample(sbatch, EW,projectDir+"/"+inDir+"/"+samples.get(i),projectDir+"/"+outDir+"/"+samples.get(i), hours, threeAdaptersFile, otherAdaptersFile, timestamp);
				EW.println();
				EW.println("wait");
			}
			EW.println();
			EW.flush();
			EW.close();


		}catch(Exception E){E.printStackTrace();}
	} 




	public  void cutAdaptSample(SBATCHinfo sbatch,ExtendedWriter generalSbatchScript, String inDir, String outDir, String time, String threeAdaptersFile, String otherAdaptersFile, String timestamp){


		if(!IOTools.isDir(outDir))
			IOTools.mkDir(outDir);
		if(!IOTools.isDir(outDir+"/reports"))
			IOTools.mkDir(outDir+"/reports");
		if(!IOTools.isDir(outDir+"/scripts"))
			IOTools.mkDir(outDir+"/scripts");
		if(QC){
			if(!IOTools.isDir(outDir+"/QC"))
				IOTools.mkDir(outDir+"/QC");
		}

		
		
		ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir, seqType);
		if(fileNames.isEmpty()){
			System.out.println("No .fastq files in folder :"+inDir);
			return;
		}

		ArrayList <String> threeAdapters = getAdapters(threeAdaptersFile);
		ArrayList <String> otherAdapters = getAdapters(otherAdaptersFile);


		try{

			ArrayList <String[]> pairs = IOTools.findPairs(fileNames,sep);
	
			
			for(int i = 0; i < pairs.size(); i++){
				String sbatchFile = outDir+"/scripts/"+i+"_"+timestamp+"_cutadapt.sbatch";
				generalSbatchScript.println("sbatch "+ sbatchFile);
				ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFile));
				String[] split = outDir.split("/");

				sbatch.printSBATCHinfoCore(EW,outDir,timestamp,0,i+"_cutadapt_"+split[split.length-1], time);


				EW.println();
				EW.println();
				for(int k = 0; k < 2; k++){
					EW.print("cutadapt --overlap=30");
					if(cutoff > -1)EW.print(" -q "+cutoff+" ");
					for (int j = 0; j < threeAdapters.size(); j++){
						EW.print(" -a "+threeAdapters.get(j));
					}
					for (int j = 0; j < otherAdapters.size(); j++){
						EW.print(" -b "+otherAdapters.get(j));
					}
					EW.print(" -o "+ outDir+"/"+pairs.get(i)[k]);
					EW.println(" "+inDir+"/"+pairs.get(i)[k]+"&");
					
				}
				EW.println("wait");
				if(QC)
					EW.println("java -jar /bubo/home/h17/johanr/bin/HTStools.jar -p sequenceHandling QC -dir "+outDir+ " -outDir "+ outDir+"/QC -l "+length +" -f1 "+pairs.get(i)[0]+" -f2 "+pairs.get(i)[1]);
				EW.println();
				EW.flush();
				EW.close();
			}		
		}catch(Exception E){E.printStackTrace();}
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
