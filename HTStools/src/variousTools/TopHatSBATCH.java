package variousTools;

import general.ExtendedWriter;
import general.Functions;
import general.IOTools;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Hashtable;

public class TopHatSBATCH {

	String referenceFile;
	String time;
	String projectDir;
	String suffix;
	String split;
	String[] sep; 

	int innerDistance;
	int innerDistanceSTD;

	public TopHatSBATCH(){
		this.referenceFile = "Database";
		this.innerDistance = 150;
		this.innerDistanceSTD = 60;
		this.split = ".";
		this.sep = new String[2];
		sep[0] = "1.fastq";
		sep[1] = "2.fastq";

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
		TopHatSBATCH shrimp = new TopHatSBATCH();
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

		if(T.containsKey("-refFile"))
			referenceFile= Functions.getValue(T, "-refFile", ".");
		else{
			System.out.println("must contain referenceFile -refFile");
			allPresent = false;
		}

		innerDistanceSTD= Integer.parseInt(Functions.getValue(T, "-rSTD", "70"));
		innerDistance= Integer.parseInt(Functions.getValue(T, "-r", "150"));
		if(!T.containsKey("-r")){
			System.out.println("must contain innerDistance -r");
			System.out.println("set to default 150");
		}
		if(T.containsKey("-time"))
			time = Functions.getValue(T, "-time", ".");
		else{
			System.out.println("must contain likely time -time");
			allPresent = false;
		}
		suffix = Functions.getValue(T,"-suffix","fastq");


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
			tophat(sbatch, timeStamp, inDir, outDir);
		else
			System.out.println("\n\nAborting run because of missing arguments for tophat.");
	}


	public void tophat(SBATCHinfo sbatch, String timeStamp, String inDir, String outDir){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_tophat.sh"));
			tophatDir(EW,sbatch, timeStamp, inDir, outDir);
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	} 



	public void tophatPair(ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String inDir, String outDir, int count, String right, String left){
		String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_"+count+"_tophat.sbatch";
		generalSbatchScript.println("sbatch "+ sbatchFileName);

		String [] temp =  referenceFile.split("/");
		String refName = referenceFile;
		if(temp.length > 1)
			refName = temp[temp.length-1];

		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
			String[] split1 = outDir.split("/");
			sbatch.printSBATCHinfoCore(EW,outDir,timestamp,count,"tophat_"+split1[split1.length-1]+"_"+count, time);

			
			if(!IOTools.isDir(outDir+"/pair_"+count))
				IOTools.mkDir(outDir);
			
			
			EW.println("cd "+inDir);
			EW.println("module load bioinfo-tools");
			EW.println("module load tophat/2.0.4");

			EW.println();
			EW.print("tophat -o "+outDir+"/pair_"+count+"  -r "+this.innerDistance+" --mate-std-dev "+this.innerDistanceSTD +" "+this.referenceFile);
			EW.print(" "+right);
			EW.print(" "+left);
			EW.println();
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}


	}

	public void tophatDir( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String inDir, String outDir){

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
			ArrayList <String[]> pairs = IOTools.findPairs(fileNames,sep);
			if(!pairs.isEmpty()){
				for(int i = 0; i < pairs.size(); i++){				
					tophatPair(generalSbatchScript,sbatch,timestamp,inDir,outDir,i,pairs.get(i)[0],pairs.get(i)[1]);
				}
			}
		}
		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
		for(int i = 0; i < subDirs.size(); i++){
			tophatDir(generalSbatchScript,sbatch,timestamp,inDir+"/"+subDirs.get(i),outDir+"/"+subDirs.get(i));
		}
	}
	
}
