package variousTools;

import general.ExtendedWriter;
import general.Functions;
import general.IOTools;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Hashtable;

public class ShrimpSBATCH {

	String inDir;
	String referenceFile;
	String outDir; 
	String logDir;
	String time;
	String projectDir;
	String suffix;
	int length;
	int percentage;
	String programFile;

	public ShrimpSBATCH(){
		this.inDir ="reads";
		this.referenceFile = "Database";
		this.outDir = "results";
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
		ShrimpSBATCH shrimp = new ShrimpSBATCH();
		shrimp.run(T);
	}

	public void run(Hashtable<String,String> T){

		boolean allPresent = true;

		String timeStamp = Functions.getDateTime();
		SBATCHinfo sbatch = new SBATCHinfo();
		if(!sbatch.addSBATCHinfo(T)) allPresent = false;

		if(T.containsKey("-inDir"))
			inDir= Functions.getValue(T, "-inDir", ".");
		else{
			System.out.println("must contain inDirectory -inDir");
			allPresent = false;
		}
		
		projectDir= Functions.getValue(T, "-pDir", "notPresent");
		
		if(T.containsKey("-outDir"))
			outDir= Functions.getValue(T, "-outDir", ".");
		else{
			System.out.println("must contain inDirectory -outDir");
			allPresent = false;
		}

		if(T.containsKey("-referenceFile"))
			referenceFile= Functions.getValue(T, "-referenceFile", ".");
		else{
			System.out.println("must contain referenceFile -referenceFile");
			allPresent = false;
		}

		if(T.containsKey("-time"))
			time = Functions.getValue(T, "-time", ".");
		else{
			System.out.println("must contain likely time -time");
			allPresent = false;
		}
		logDir = Functions.getValue(T, "-logDir", outDir+"/log");
		length = Integer.parseInt(Functions.getValue(T,"-length","35"));
		programFile = Functions.getValue(T,"-programFile","gampper-cs");
		suffix = Functions.getValue(T,"-suffix","csfasta");
		
		
		if(projectDir.compareTo("notPresent") != 0){
			inDir = projectDir+"/"+inDir;
			outDir = projectDir+"/"+outDir;
			referenceFile = projectDir+"/"+referenceFile;
			logDir = projectDir+"/"+logDir; 
		}
		else{
			projectDir = inDir;
		}
		
		
		if(allPresent)
			Shrimp(sbatch, timeStamp);
		else
			System.out.println("\n\nAborting run because of missing arguments for shrimp.");
	}


	public void Shrimp(SBATCHinfo sbatch, String timeStamp){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_SHRIMP.sh"));
			shrimpDir(EW,sbatch, timeStamp);
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	} 




	public void shrimpDir(ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp){

		ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir,suffix);

		
		if(fileNames.isEmpty()){
			System.out.println("No "+suffix+" files in folder :"+inDir);
		}
		else{
			if(!IOTools.isDir(inDir+"/reports"))
				IOTools.mkDir(inDir+"/reports");
			if(!IOTools.isDir(inDir+"/scripts"))
				IOTools.mkDir(inDir+"/scripts");
			if(!IOTools.isDir(outDir))
				IOTools.mkDir(outDir);
			if(!IOTools.isDir(inDir+"/scripts"))
				IOTools.mkDir(inDir+"/scripts");
			int count = 0;

			try{
				
				String 	sbatchFileName = inDir+"/scripts/"+timestamp+"_"+count+".SHRIMP.sbatch";

				generalSbatchScript.println("sbatch "+ sbatchFileName);

				ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
				sbatch.printSBATCHinfo(EW,inDir,timestamp,count,"SHRIMP", time);

					
				EW.println("cd "+inDir);
				EW.println();
				EW.println();

				String [] temp =  referenceFile.split("/");
				String refName = referenceFile;
				if(temp.length > 1)
					refName = temp[temp.length-1];
					
				for(int i = 0; i < fileNames.size(); i++){	
					EW.println("gmapper-cs "+
								"  --strata -E -n 1 -N 1 --global -o 100 -h 85% "+
								inDir+"/"+fileNames.get(i)+" "+this.referenceFile+
								">"+outDir+"/"+fileNames.get(i)+"_"+refName+".sam 2>"+outDir+"/"+fileNames.get(i)+"_"+refName+".sam.log");
					if((i+1)%8 == 0  &&  fileNames.size()-i != 1){
						
						EW.println();
						EW.println();
						EW.println("wait");
						EW.flush();
						EW.close();
						count++;
						
						sbatchFileName = inDir+"/scripts/"+timestamp+"_"+count+".SHRIMP.sbatch";
						generalSbatchScript.println("sbatch "+ sbatchFileName);
						EW = new ExtendedWriter(new FileWriter(sbatchFileName));
						sbatch.printSBATCHinfo(EW,inDir,timestamp,count,"SHRIMP", time);
						EW.println("cd "+inDir);
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
			
			count++;
			
		}

	}


}
