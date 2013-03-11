package variousTools;

import general.ExtendedWriter;
import general.Functions;
import general.IOTools;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import Sequence.FastQSequences;
import Sequence.FastaSequences;
import Sequence.SequenceHandling;

public class deNovoExtension {

	String time;
	String projectDir;

	String suffix;
	String split;
	String[] sep; 

	public deNovoExtension(){
		this.split = ".";
		projectDir = time =  null;
	}




	public void run(Hashtable<String,String> T){

		boolean allPresent = true;

		SBATCHinfo sbatch = new SBATCHinfo();
		if(!sbatch.addSBATCHinfo(T)) allPresent = false;

		String timeStamp = Functions.getValue(T, "-TS",Functions.getDateTime());


		projectDir= Functions.getValue(T, "-pDir", IOTools.getCurrentPath());

		if(allPresent){
			ExtensionScript(sbatch,timeStamp, projectDir);
		}
		else
			System.out.println("\n\nAborting run because of missing arguments for trinity.");
	}



	public void ExtensionScript(SBATCHinfo sbatch, String timeStamp, String inDir){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_PE_trinity.sh"));
			EW.println("cd "+inDir);
			findSpecies(EW,sbatch, timeStamp, inDir);
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	} 

	public void findSpecies(ExtendedWriter EW, SBATCHinfo sbatch, String timeStamp, String inDir){

		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
		for(int i = 0; i < subDirs.size(); i++){
			findSamples(EW,sbatch,timeStamp,inDir+"/"+subDirs.get(i),subDirs.get(i));
		}
	}

	public void findSamples(ExtendedWriter EW, SBATCHinfo sbatch, String timeStamp, String inDir,String species){

		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
		for(int i = 0; i < subDirs.size(); i++){
			findMethod(EW,sbatch,timeStamp,inDir+"/"+subDirs.get(i),species, subDirs.get(i));
		}
	}


	public void findMethod(ExtendedWriter EW, SBATCHinfo sbatch, String timeStamp, String inDir,String species,String sample){

		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
		for(int i = 0; i < subDirs.size(); i++){
			if(subDirs.get(i).compareTo("trinity")==0){
				findProgress(EW,sbatch,timeStamp,inDir+"/"+subDirs.get(i),species, sample,"trinity");
			}
			if(subDirs.get(i).compareTo("oases")==0){
				findProgress(EW,sbatch,timeStamp,inDir+"/"+subDirs.get(i),species, sample,"oases");
			}
		}
	}


	public void findProgress(ExtendedWriter EW, SBATCHinfo sbatch, String timeStamp, String inDir,String species,String sample,String method){
		System.out.println("Checking progress of "+method+ "progress in sample "+ sample +" in species "+species);
		String round = "round";
		int count = 1;
		boolean lastFound = false;
		while(IOTools.fileExists(inDir+"/"+round+count+"/AllSequences_ends.fa")){count++;}
		System.out.println(count-1+" rounds are finished");
		System.out.println("Checking status of round"+count);
		ArrayList <String> samFile = IOTools.getSequenceFiles(inDir,".sam");
		if(samFile.size() == 1){
			if(IOTools.fileExists(inDir+"/"+round+count+"/Tend.txt")){
				System.out.println("Extraction of extensions are finished. Writing shellscript for merging extensions");
				if(count==1){
					EW.println("HTStools -p sbatch -script -node -i scripts/HTStools.merge.sh -time 4:00:00  -parameters "+species+"/"+sample+"/"+method+"/"+round+count+"/"+species+"."+sample+"."+method+".initial_"+species+"."+sample+"."+method+".initial.stict.sam initial");
				}
				else{
					EW.println("HTStools -p sbatch -script -node -i scripts/HTStools.merge.sh -time 4:00:00  -parameters "+species+"/"+sample+"/"+method+"/"+round+count+"/"+species+"."+sample+"."+method+".round"+(count-1)+"_"+species+"."+sample+"."+method+".round"+(count-1)+".stict.sam initial");
				}
			}
			else{
				System.out.println("Mapping of reads to transcripts are finished. Writing shellscript for extracting ends");
				if(count==1){
					EW.println("HTStools -p sbatch -script -node -i scripts/HTStools.merge.sh -time 4:00:00  -parameters "+species+"/"+sample+"/"+method+"/"+round+count+"/"+species+"."+sample+"."+method+".initial_"+species+"."+sample+"."+method+".initial.stict.sam initial");
				}
				else{
					EW.println("HTStools -p sbatch -script -node -i scripts/HTStools.merge.sh -time 4:00:00  -parameters "+species+"/"+sample+"/"+method+"/"+round+count+"/"+species+"."+sample+"."+method+".round"+(count-1)+"_"+species+"."+sample+"."+method+".round"+(count-1)+".stict.sam round"+(count-1));
				}
			}
		}else{
			System.out.println(" Writing shellscript to map reads to transcripts");
			if(count==1){
				EW.println("HTStools -p sbatch -script -node -i scripts/HTStools.map.sh -time 4:00:00  -parameters "+species+"/"+sample+"/"+method+"/"+round+count+" initial 30:00:00");
			}
			else{
				EW.println("HTStools -p sbatch -script -node -i scripts/HTStools.map.sh -time 4:00:00  -parameters "+species+"/"+sample+"/"+method+"/"+round+count+" round"+(count-1)+" 5:00:00");
			}
		}
	}

}

