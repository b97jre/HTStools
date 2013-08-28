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
	boolean extend;
	boolean blat;
	boolean Final;
	boolean clean;

	public deNovoExtension(){
		this.split = ".";
		projectDir = time =  null;
		blat = false;
		extend = true;
		Final = false;
		clean = false;
	}




	public void run(Hashtable<String,String> T){

		boolean allPresent = true;

		SBATCHinfo sbatch = new SBATCHinfo();
		if(!sbatch.addSBATCHinfo(T)) allPresent = false;

		String timeStamp = Functions.getValue(T, "-TS",Functions.getDateTime());
		projectDir= Functions.getValue(T, "-pDir", IOTools.getCurrentPath());

		if(T.containsKey("-blat")){
			extend = false;
			blat=true;
		}

		if(T.containsKey("-Final")){
			extend = false;
			Final=true;
		}

		if(T.containsKey("-clean")){
			extend = false;
			clean=true;
		}
	
		
		if(allPresent){
			ExtensionScript(sbatch,timeStamp, projectDir);
		}
		else
			System.out.println("\n\nAborting run because of missing arguments.");
	}



	public void ExtensionScript(SBATCHinfo sbatch, String timeStamp, String inDir){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_deNovoExtension.sh"));
			ExtendedWriter EW2 = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_deNovoExtension.info"));
			ExtendedWriter EW3 = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_deNovoExtension_RUNME.sh"));
			EW3.println("sh "+projectDir+"/scripts/"+timeStamp+"_deNovoExtension.sh");
			EW.println("cd "+inDir);
			findSpecies(EW, EW2, EW3, sbatch, timeStamp, inDir);
			EW.flush();
			EW.close();
			EW2.flush();
			EW2.close();
			EW3.flush();
			EW3.close();
		}catch(Exception E){E.printStackTrace();}
	} 

	public void findSpecies(ExtendedWriter EW,ExtendedWriter EW2,ExtendedWriter EW3, SBATCHinfo sbatch, String timeStamp, String inDir){
		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
		System.out.println("finding species");
		for(int i = 0; i < subDirs.size(); i++){
			findSamples(EW, EW2, EW3,sbatch,timeStamp+"_"+subDirs.get(i),inDir+"/"+subDirs.get(i),subDirs.get(i));
		}
	}

	public void findSamples(ExtendedWriter EW, ExtendedWriter EW2,ExtendedWriter EW3, SBATCHinfo sbatch, String timeStamp, String inDir,String species){
		System.out.println("finding samples");

		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
		for(int i = 0; i < subDirs.size(); i++){
			findMethod(EW, EW2, EW3,sbatch,timeStamp+"_"+subDirs.get(i),inDir+"/"+subDirs.get(i),species, subDirs.get(i));
		}
	}


	public void findMethod(ExtendedWriter EW, ExtendedWriter EW2,ExtendedWriter EW3, SBATCHinfo sbatch, String timeStamp, String inDir,String species,String sample){

		System.out.println("finding method "+timeStamp);
		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
		for(int i = 0; i < subDirs.size(); i++){
			if(subDirs.get(i).compareTo("trinity")==0 || subDirs.get(i).compareTo("oases")==0){
				if(extend){
					findProgress(EW, EW2, EW3,sbatch,timeStamp+"_"+subDirs.get(i),inDir+"/"+subDirs.get(i),species, sample,subDirs.get(i));
				}
				if(Final){
					mergeFinal(EW, EW2, EW3,sbatch,timeStamp+"_"+subDirs.get(i),inDir+"/"+subDirs.get(i),species, sample,subDirs.get(i));
				}
				if(blat){
					blat(EW, EW2, EW3,sbatch,timeStamp+"_"+subDirs.get(i),inDir+"/"+subDirs.get(i),species, sample,subDirs.get(i));
				}
				if(clean){
					clean(EW, EW2, EW3,sbatch,timeStamp+"_"+subDirs.get(i),inDir+"/"+subDirs.get(i),species, sample,subDirs.get(i));
				}
				
			}
		}
	}


	public void findProgress(ExtendedWriter EW, ExtendedWriter EW2,ExtendedWriter EW3, SBATCHinfo sbatch, String timeStamp, String inDir,String species,String sample,String method){
		EW2.println("Checking progress of "+method+ "progress in sample "+ sample +" in species "+species);
		System.out.println("Checking progress of "+method+ "progress in sample "+ sample +" in species "+species);

		String round = "round";
		int count = 1;
		while(IOTools.isDir(inDir+"/"+round+count) && IOTools.fileExists(inDir+"/"+round+count+"/AllSequences_ends.fa")){count++;}
		EW2.println(count-1+" rounds are finished");

		if(IOTools.isDir(inDir+"/"+round+count)){
			EW2.println("Checking status of round"+count);
			ArrayList <String> samFile = IOTools.getSequenceFiles(inDir+"/"+round+count,"sam");
			System.out.println(samFile.size());
			if(samFile.size() == 1){
				if(IOTools.isDir(inDir+"/"+round+count+"/Tend")){
					EW2.println("Extraction of extensions are finished. Writing shellscript for merging extensions");
					EW3.println("sbatch "+timeStamp+"_HTStools.merge.sh.sbatch");
					if(count==1){
						EW.println("java -Xmx10G -jar ~/bin/HTStools.jar -p sbatch -script -TS "+timeStamp+" -node -i scripts/HTStools.merge.sh -time 4:00:00  -parameters "+species+"/"+sample+"/"+method+"/"+round+count+"/"+samFile.get(0)+" initial");
					}
					else{
						EW.println("java -Xmx10G -jar ~/bin/HTStools.jar  -p sbatch -script -TS "+timeStamp+" -node -i scripts/HTStools.merge.sh -time 4:00:00  -parameters "+species+"/"+sample+"/"+method+"/"+round+count+"/"+samFile.get(0)+" round"+(count-1));
					}
				}
				else{
					EW2.println("Mapping of reads to transcripts are finished. Writing shellscript for extracting ends");
					EW3.println("sbatch "+timeStamp+"_HTStools.extract.sh.sbatch");
					if(count==1){
						EW.println("java -Xmx10G -jar ~/bin/HTStools.jar  -p sbatch -script -TS "+timeStamp+" -node -i scripts/HTStools.extract.sh -time 4:00:00  -parameters "+species+"/"+sample+"/"+method+"/"+round+count+"/"+samFile.get(0)+" initial");
					}
					else{
						EW.println("java -Xmx10G -jar ~/bin/HTStools.jar  -p sbatch -script -TS "+timeStamp+" -node -i scripts/HTStools.extract.sh -time 4:00:00  -parameters "+species+"/"+sample+"/"+method+"/"+round+count+"/"+samFile.get(0)+" round"+(count-1));
					}
				}
			}else{
				EW2.println(" Writing shellscript to map reads to transcripts");
				EW3.println("sbatch "+timeStamp+"_HTStools.map.sh.sbatch");
				if(count==1){
					EW.println("java -Xmx10G -jar ~/bin/HTStools.jar  -p sbatch -script -TS "+timeStamp+" -node -i scripts/HTStools.map.sh -time 4:00:00  -parameters "+species+"/"+sample+"/"+method+"/"+round+count+" initial 30:00:00");
				}
				else{
					EW.println("java -Xmx10G -jar ~/bin/HTStools.jar  -p sbatch -script -TS "+timeStamp+" -node -i scripts/HTStools.map.sh -time 4:00:00  -parameters "+species+"/"+sample+"/"+method+"/"+round+count+" round"+(count-1)+" 15:00:00");
				}
			}
		}else{
			EW2.println("All rounds according to folder structure are finished");
		}

		EW2.println();
		EW2.println();
		EW.println();
		EW.println();
	}


	public void clean(ExtendedWriter EW, ExtendedWriter EW2,ExtendedWriter EW3, SBATCHinfo sbatch, String timeStamp, String inDir,String species,String sample,String method){
		EW2.println("Checking progress of "+method+ "progress in sample "+ sample +" in species "+species);
		System.out.println("Checking progress of "+method+ "progress in sample "+ sample +" in species "+species);

		String round = "round";
		int count = 1;
		while(IOTools.isDir(inDir+"/"+round+count) && IOTools.fileExists(inDir+"/"+round+count+"/AllSequences_ends.fa")){
			ArrayList <String> samFile = IOTools.getSequenceFiles(inDir+"/"+round+count,"sam");
			//System.out.println(samFile.size());
			if(samFile.size() == 1){
				EW.println("java -Xmx10G -jar ~/bin/HTStools.jar  -p sbatch -script -TS "+timeStamp+"_round"+count+" -node -i scripts/HTStools.clean.sh -time 4:00:00 -modules samtools  -parameters "+species+"/"+sample+"/"+method+"/"+round+count+"/"+samFile.get(0));
				EW2.println("Writing script to remove all the  temporary files in "+inDir+"/"+round+count);
				EW3.println("sbatch "+timeStamp+"_round"+count+"_HTStools.clean.sh.sbatch");
			}
			count++;
		}
		EW2.println(count-1+" rounds are finished");


		EW2.println();
		EW2.println();
		EW.println();
		EW.println();
	}
	
	
	
	
	public void mergeFinal(ExtendedWriter EW, ExtendedWriter EW2,ExtendedWriter EW3, SBATCHinfo sbatch, String timeStamp, String inDir,String species,String sample,String method){
		EW2.println("Checking progress of "+method+ "progress in sample "+ sample +" in species "+species);
		System.out.println("Checking progress of "+method+ "progress in sample "+ sample +" in species "+species);

		String round = "round";
		int count = 1;
		while(IOTools.isDir(inDir+"/"+round+count) && IOTools.fileExists(inDir+"/"+round+count+"/AllSequences_ends.fa")){count++;}
		EW2.println(count-1+" rounds are finished");

		if(!IOTools.isDir(inDir+"/Final")){
			IOTools.mkDir(inDir+"/Final");
		}
		if(IOTools.isDir(inDir+"/Final")){
			EW2.println("Checking status of Final round ");
			ArrayList <String> samFile = IOTools.getSequenceFiles(inDir+"/Final","sam");
			//System.out.println(samFile.size());
			if(samFile.size() == 1){
//				if(!IOTools.fileExists(inDir+"/Final/AllSequences_merged.fa")){
					EW2.println("Extraction of extensions are finished. Writing shellscript for merging extensions");
					EW3.println("sbatch "+timeStamp+"_HTStools.mergeFinal.sh.sbatch");
					EW.println("java -Xmx10G -jar ~/bin/HTStools.jar  -p sbatch -script -TS "+timeStamp+" -node -i scripts/HTStools.mergeFinal.sh -time 12:00:00  -parameters "+species+"/"+sample+"/"+method+"/Final/"+samFile.get(0) +" round"+(count-1));
//				}else{
//					System.out.println("Final assembly is finished can be found at "+inDir+"/Final/AllSequences_merged.fa");
//				}
			}else{
				EW2.println(" Writing shellscript to map reads to transcripts");
				EW3.println("sbatch "+timeStamp+"_HTStools.mapFinal.sh.sbatch");
				EW.println("java -Xmx10G -jar ~/bin/HTStools.jar  -p sbatch -script -TS "+timeStamp+" -node -i scripts/HTStools.mapFinal.sh -time 48:00:00  -parameters "+species+"/"+sample+"/"+method+"/Final"+" round"+(count-1)+" initial 15:00:00");
			}
		}
		EW2.println();
		EW2.println();
		EW.println();
		EW.println();
	}
	
	
	
	
	public void blat(ExtendedWriter EW, ExtendedWriter EW2,ExtendedWriter EW3, SBATCHinfo sbatch, String timeStamp, String inDir,String species,String sample,String method){
		EW2.println("Blating files of "+method+ " in sample "+ sample +" in species "+species+" against reference transcriptome");
		if(!IOTools.isDir(projectDir+"/"+species+"/blatResults"))
			IOTools.mkDir(projectDir+"/"+species+"/blatResults");
		ArrayList <String> referenceFiles = IOTools.getSequenceFiles(projectDir+"/"+species+"/references","fa");
		ArrayList <String> assemblyFiles = IOTools.getSequenceFiles(inDir+"/references","fa");


		for(int i = 0; i < referenceFiles.size();i++){
			for(int j = 0; j < assemblyFiles.size();j++){
				Blat B = new Blat();
				B.blatScript(EW3, sbatch, timeStamp+"_"+assemblyFiles.get(j)+"_"+referenceFiles.get(i), inDir+"/references/"+assemblyFiles.get(j), projectDir+"/"+species+"/references/"+referenceFiles.get(i), projectDir+"/"+species+"/blatResults");
				B.blatScript(EW3, sbatch, timeStamp+"_"+referenceFiles.get(i)+"_"+assemblyFiles.get(j), projectDir+"/"+species+"/references/"+referenceFiles.get(i), inDir+"/references/"+assemblyFiles.get(j), projectDir+"/"+species+"/blatResults");
			}
		}

		EW2.println();
		EW2.println();
		EW.println();
		EW.println();
	}


}

