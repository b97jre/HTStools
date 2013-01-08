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

public class DigiNorm {

	String codeFile;
	String time;
	String[] sep; 
	String suffix;
	boolean hiseq;
	int  C;
	int k;
	int N;
	String x ;

	public DigiNorm(){
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
		DigiNorm filter = new DigiNorm();
		filter.run(T);
	}

	public void run(Hashtable<String,String> T){

		String inDir = null;
		String outDir = null; 
		String projectDir = null;

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
			projectDir= Functions.getValue(T, "-pDir", IOTools.getCurrentPath());

	
			time = Functions.getValue(T, "-time", "2:00:00");
		if(!T.containsKey("-time")){
			System.out.println("Default time is set to 2:00:00");
		}

		suffix = Functions.getValue(T, "-suffix", "fastq");
		codeFile = Functions.getValue(T, "-codeFile", "~/bin/khmer/scripts/normalize-by-median.py");
		C = Integer.parseInt(Functions.getValue(T, "-C", "20"));
		k = Integer.parseInt(Functions.getValue(T, "-k", "20"));
		N = Integer.parseInt(Functions.getValue(T, "-N", "4"));
		x = Functions.getValue(T, "-x", "2e9");
	
		if(allPresent)
			DigiNormTop(sbatch, timeStamp,projectDir,inDir,outDir);
		else
			System.out.println("\n\nAborting run because of missing arguments for DigiNorm.");
	}


	public void DigiNormTop(SBATCHinfo sbatch, String timeStamp,String projectDir, String inDir, String outDir){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			if(!IOTools.isDir(projectDir+"/"+outDir))
				IOTools.mkDir(projectDir+"/"+outDir);
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_filter_fastq.sh"));
			DigiNormDir(EW,sbatch,projectDir+"/"+inDir ,projectDir+"/"+outDir, timeStamp,0);
				EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	} 




	public void DigiNormDir(ExtendedWriter generalSbatchScript, SBATCHinfo sbatch , String inDir, String outDir,String timestamp, int count2){


		if(!IOTools.isDir(outDir))
			IOTools.mkDir(outDir);
		

		ArrayList <String> samples = IOTools.getDirectories(inDir);
		for(int i = 0; i < samples.size(); i++){
			DigiNormDir(generalSbatchScript,sbatch,inDir+"/"+samples.get(i) ,outDir+"/"+samples.get(i), timestamp,i);
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
				sbatch.printSBATCHinfo(EW,outDir,timestamp,count2,"DigiNorm_"+count2, time);
				EW.println("cd "+inDir);
				EW.println("PYTHONPATH=~/bin/screed/dist/screed-0.7-py2.6.egg:~/bin/khmer/python ");
				EW.println("export PYTHONPATH");
				
				// print all suffix files files to one file 
				String outFile = "all."+suffix;
				if(suffix.indexOf("gz") > -1)
					EW.println("zcat *."+suffix+" | gzip > "+outDir+"/all."+suffix);
				else{
					EW.println("cat *."+suffix+" | gzip > "+outDir+"/all."+suffix+".gz");
					outFile = "all."+suffix+".gz";
				}
				EW.println("cd "+outDir);
				EW.println(codeFile +" -C "+C+" -k "+k+" -N "+N+" -x "+x+" "+outFile);
				EW.println();
				
				EW.println("python ~/bin/khmer/sandbox/split-pe.py "+outFile+".keep" );
				EW.println("rm "+outFile);
				EW.println("rm "+outFile+".keep");
				EW.println("mv "+outFile+".keep.1 "+outFile+".keep.1.fa" );
				EW.println("mv "+outFile+".keep.2 "+outFile+".keep.2.fa" );
				
				
				
				EW.flush();
				EW.close();
				
			}catch(Exception E){E.printStackTrace();}
		}

	}


}
