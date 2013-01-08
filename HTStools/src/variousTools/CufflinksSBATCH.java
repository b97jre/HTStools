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

public class CufflinksSBATCH {

	String referenceFile;
	String time;
	String projectDir;
	String suffix;
	String split;
	String[] sep; 
	
	int innerDistance;
	int innerDistanceSTD;

	public CufflinksSBATCH(){
		this.referenceFile = "Database";

		
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
		CufflinksSBATCH shrimp = new CufflinksSBATCH();
		shrimp.run(T);
	}

	public void run(Hashtable<String,String> T){
		
		String inDir, outDir, logDir;
		inDir = outDir = logDir = null;
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
		
		if(T.containsKey("-o"))
			outDir= Functions.getValue(T, "-o", ".");
		else{
			System.out.println("must contain inDirectory -o");
			allPresent = false;
		}
		if(T.containsKey("-time"))
			time = Functions.getValue(T, "-time", ".");
		else{
			System.out.println("must contain likely time -time");
			allPresent = false;
		}
		suffix = Functions.getValue(T,"-suffix","sam");
		
		
		if(projectDir.compareTo("notPresent") != 0){
			inDir = projectDir+"/"+inDir;
			outDir = projectDir+"/"+outDir;
		}
		else{
			projectDir = inDir;
		}
		
		
		if(allPresent)
			cufflinks(sbatch, timeStamp, inDir, outDir);
		else
			System.out.println("\n\nAborting run because of missing arguments for cufflinks.");
	}


	public void cufflinks(SBATCHinfo sbatch, String timeStamp, String inDir, String outDir){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_cufflinks.sh"));
			cufflinksDir(EW,sbatch, timeStamp, inDir, outDir);
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	} 



	public void cufflinksDir( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String inDir, String outDir){

		
		ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir,suffix);
		
		if(fileNames.isEmpty()){
			System.out.println("No "+suffix+" files in folder :"+inDir);
			if(!IOTools.isDir(outDir))
				IOTools.mkDir(outDir);
		}
		else{
			if(!IOTools.isDir(outDir))
				IOTools.mkDir(outDir);
			if(!IOTools.isDir(outDir+"/reports"))
				IOTools.mkDir(outDir+"/reports");
			if(!IOTools.isDir(outDir+"/scripts"))
				IOTools.mkDir(outDir+"/scripts");
			try{
				
					
				for(int i = 0; i < fileNames.size(); i++){	
					String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_"+i+"_cufflinks.sbatch";
					generalSbatchScript.println("sbatch "+ sbatchFileName);
					
					ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
					String[] split1 = outDir.split("/");
					sbatch.printSBATCHinfoCore(EW,outDir,timestamp,i,"cufflinks_"+split1[split1.length-1], time);

					EW.println("cd "+inDir);
					EW.println("export LD_LIBRARY_PATH=/home/delhomme/lib/");

					EW.println();
					EW.println();
					String fileBase = fileNames.get(i).substring(0,fileNames.get(i).indexOf(suffix)-1);
					if(!IOTools.isDir(outDir+"/"+fileBase))
						IOTools.mkDir(outDir+"/"+fileBase);
					EW.println ("cufflinks  -o "+outDir+"/"+fileBase+" "+fileNames.get(i));
					EW.println();
					EW.println();
					EW.println("wait");
					EW.flush();
					EW.close();
				}
			}catch(Exception E){E.printStackTrace();}
		}
		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
		for(int i = 0; i < subDirs.size(); i++){
			cufflinksDir(generalSbatchScript,sbatch,timestamp,inDir+"/"+subDirs.get(i),outDir+"/"+subDirs.get(i));
		}
		
		

	}

	

}
