package variousTools;

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

public class Picard {

	String time;
	String projectDir;
	String suffix;
	String jarDir;
	int memory;

	boolean strandSpecifik;



	public Picard(){
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
		Picard Picard = new Picard();
		Picard.run(T);
	}

	public void run(Hashtable<String,String> T){

		String inDir, logDir;
		inDir = logDir = null;
		boolean allPresent = true;


		String timeStamp = Functions.getValue(T, "-TS", Functions.getDateTime());
		SBATCHinfo sbatch = new SBATCHinfo();
		if(!sbatch.addSBATCHinfo(T)){
			allPresent = false;
			return;
		}
		if(T.containsKey("-i"))
			inDir= Functions.getValue(T, "-i", ".");
		else{
			System.out.println("must contain inDirectory -i ");
			allPresent = false;
		}
		jarDir= Functions.getValue(T, "-jDir", "/bubo/home/h17/johanr/bin");
		time = Functions.getValue(T, "-time", "3:00:00");

		if(!T.containsKey("-jDir"))
			System.out.println("must contain directory where picard jar files are kept. -jDir  now set to default /bubo/home/h17/johanr/bin");
		if(T.containsKey("-time"))
			System.out.println("must contain likely time -time now set to default 3:00:00");

		memory = Functions.getInt(T, "-X", 3);


		suffix = Functions.getValue(T,"-suffix","bam");
		projectDir = Functions.getValue(T,"-projectDir",IOTools.getCurrentPath());
		if(suffix.indexOf('.')== 0)suffix = suffix.substring(1);

		if(allPresent){
			if(IOTools.isDir(projectDir+"/"+inDir)){
				inDir = projectDir+"/"+inDir;
			}
			else if(!IOTools.isDir(inDir)){
				System.out.println("Neither "+ inDir +" nor "+projectDir+"/"+inDir+ "was found");
				return;
			}
			PicardInitial(sbatch, timeStamp, inDir,T);
		}
		else
			System.out.println("\n\nAborting run because of missing arguments for Picard.");
	}

	public void PicardInitial(SBATCHinfo sbatch, String timeStamp, String inDir,Hashtable<String,String> T){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			String lastDir = inDir.substring(inDir.lastIndexOf("/")+1); 

			ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"."+lastDir+".Picard.sh"));
			PicardDir(EW,sbatch, timeStamp, inDir, T);

			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	} 




	public void PicardDir( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String inDir,Hashtable<String,String> T){

		ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir,suffix);
		if(!fileNames.isEmpty()){
			try{

				for(int i = 0; i < fileNames.size(); i++){
					System.out.println(fileNames.get(i));

					PicardFile(generalSbatchScript,sbatch,timestamp,inDir,fileNames.get(i),T,i);
				}
			}catch(Exception E){E.printStackTrace();}
		}
		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
		for(int i = 0; i < subDirs.size(); i++){
			PicardDir(generalSbatchScript,sbatch,timestamp,inDir+"/"+subDirs.get(i),T);
		}
	}

	
	
	public static String preGATK( ExtendedWriter EW, SBATCHinfo sbatch ,String timestamp,String outDir,  ArrayList <String> fileNames,Hashtable<String,String> T,String picardDir, int memory, String suffix){


		/*				
			for each sample (asumed to be in same dir)

		    lanes.bam <- merged lane.bams for sample
		    dedup.bam <- MarkDuplicates(lanes.bam) 
		    realigned.bam <- realign(dedup.bam) [with known sites included if available]
		    recal.bam <- recal(realigned.bam)
		    sample.bam <- recal.bam			
		 */

		EW.println();
		EW.println("cd "+ outDir);
		if(suffix.indexOf("sam")>-1){
			for(int i = 0; i < fileNames.size(); i++){
				String samFile = fileNames.get(i);
				String bamFile = samFile.substring(0,samFile.lastIndexOf("sam"))+"bam";
				fileNames.set(i, bamFile);
				SamtoolsSBATCH.sam2bam(EW, samFile, -1, -1, false, true, true, false, false);
			}
		}
		

		for(int i = 0; i < fileNames.size(); i++){
			System.out.println(fileNames.get(i));
			String bamFile = fileNames.get(i);
			String[] info =bamFile.split("_");
			String Name, Nucleotide, Barcode, Lane;
			Name = Nucleotide = Barcode = Lane = null;
			if(info.length >  3){
				Name = info[0];
				Nucleotide = info[1];
				Barcode = info[2];
				Lane = info[3];
			}

			//"Inter3-1_DNA_ATCACG_L005_R1_001fastq_Crubella_183.strict.sam"
			String RGLB = Functions.getValue(T, "-RGLB", "200");
			String RGPL = Functions.getValue(T, "-RGPL", "illumina");
			String RGPU = Functions.getValue(T, "-RGPU", Barcode);
			String RGSM = Functions.getValue(T, "-RGSM", Name);
			String RGID = Functions.getValue(T, "-RGID", RGSM+"_"+RGPU+"_"+RGLB+"_"+Nucleotide+"_"+Lane);

			bamFile = AddOrReplaceReadGroups(EW, bamFile, memory,picardDir,"bam",RGID,RGLB,RGPL,RGPU,RGSM,"coordinate")+"."+suffix;
			fileNames.set(i, bamFile);

		}


		String fileName = fileNames.get(0);
		String bamFile = fileName;
		if(fileNames.size()>1){
			fileName = Functions.getCommonPrefix(fileNames.get(0), fileNames.get(1));

			EW.println("# merging all bamFiles into one file named "+ fileName+".AddOrReplaceReadGroups.merged.bam");
			bamFile = mergeBamFiles(EW,memory,picardDir,fileNames,fileName+".AddOrReplaceReadGroups.merged.bam");
			EW.println();
			EW.println();
		}
		else{
			EW.println("# Only one file in folder so no merging of files");
		}

		EW.println("# Marking duplicates duplicates "+ bamFile);

		bamFile = markDuplicates(EW, bamFile, memory, picardDir,suffix);

		EW.println();
		EW.println();
		EW.println("wait");

		return bamFile;


	}

	public static String preGATKDNAsample( ExtendedWriter EW, SBATCHinfo sbatch ,String timestamp,String outDir,  String fileName,Hashtable<String,String> T,String picardDir, int memory, String suffix){


		/*				
			for each sample (asumed to be in same dir)

		    lanes.bam <- merged lane.bams for sample
		    dedup.bam <- MarkDuplicates(lanes.bam) 
		    realigned.bam <- realign(dedup.bam) [with known sites included if available]
		    recal.bam <- recal(realigned.bam)
		    sample.bam <- recal.bam			
		 */

		EW.println();
		EW.println("cd "+ outDir);
		if(suffix.indexOf("sam")>-1){
			String samFile = fileName;
			fileName = samFile.substring(0,samFile.lastIndexOf("sam"))+"bam";
		}

		String bamFile = fileName;
		String[] info =bamFile.split("_");
		String Name, Nucleotide, Barcode, Lane;
		Name = Nucleotide = Barcode = Lane = null;
		if(info.length >  3){
			Name = info[0];
			Nucleotide = info[1];
			Barcode = info[2];
			Lane = info[3];
		}

		//"Inter3-1_DNA_ATCACG_L005_R1_001fastq_Crubella_183.strict.sam"
		String RGLB = Functions.getValue(T, "-RGLB", Lane);
		String RGPL = Functions.getValue(T, "-RGPL", "illumina");
		String RGPU = Functions.getValue(T, "-RGPU", Barcode);
		String RGSM = Functions.getValue(T, "-RGSM", Name);
		String RGID = Functions.getValue(T, "-RGID", RGSM+"_"+RGPU+"_"+RGLB+"_"+Nucleotide+"_"+Lane);

		bamFile = AddOrReplaceReadGroups(EW, bamFile, memory,picardDir,suffix,RGID,RGLB,RGPL,RGPU,RGSM,"coordinate")+"."+suffix;
		EW.println("# Marking duplicates duplicates "+ bamFile);

		bamFile = markDuplicates(EW, bamFile, memory, picardDir,suffix);

		EW.println();
		EW.println();
		EW.println("wait");

		return bamFile;


	}




	public void PicardFile( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String outDir, String bamFile,Hashtable<String,String> T,int count){

		if(!IOTools.isDir(outDir+"/reports"))
			IOTools.mkDir(outDir+"/reports");
		if(!IOTools.isDir(outDir+"/scripts"))
			IOTools.mkDir(outDir+"/scripts");
		try{

			String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_"+count+"Picard.sbatch";
			generalSbatchScript.println("sbatch "+ sbatchFileName);

			ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
			sbatch.printSBATCHinfo(EW,outDir,timestamp,0,"Picard", time);

			EW.println();
			EW.println("cd "+ outDir);
			if(T.containsKey("-makrDuplicates")){
				bamFile = markDuplicates(EW, bamFile, memory, jarDir,suffix)+"."+suffix;
			}
			if(T.containsKey("-SortSam")){
				String sortOrder = Functions.getValue(T, "-SORT_ORDER", "coordinate");
				bamFile = SortSam(EW, bamFile, memory,jarDir,suffix, sortOrder)+"."+suffix;
			}
			if(T.containsKey("-preGATK")){
				System.out.println("Asumes that the sam filename has the following order RGSM_DNA_RGPU_RGLB_WHATEVER for DNA or");
				System.out.println("Intra8_2.Replicate.Tissue.RGPU.RGLB.Cr_reference.bam for RNAs");
				System.out.println(bamFile);
				if(T.containsKey("-RNA")){
					String[] info =bamFile.split("\\.");
					String Individual, Replicate, Tissue,Barcode, Lane;
					Individual = Replicate = Tissue = Barcode = Lane = null;
					if(info.length >  3){
						Individual = info[0];
						Replicate = info[1];
						Tissue = info[2];
						Barcode = info[3];
						Lane = info[4];
					}

					//"Intra8_2.1.L.AGTCAA.L007.Cr_reference.bam "
					String RGLB = Functions.getValue(T, "-RGLB", Lane);
					String RGPL = Functions.getValue(T, "-RGPL", "illumina");
					String RGPU = Functions.getValue(T, "-RGPU", Barcode);
					String RGSM = Functions.getValue(T, "-RGSM", Individual);
					String RGID = Functions.getValue(T, "-RGID", RGSM+"_"+Replicate+"_"+RGPU+"_"+RGLB+"_RNA_"+Lane+"_"+Tissue);


					prepareForGATK(EW, bamFile, memory, jarDir,  suffix, RGID, RGLB, RGPL, RGPU, RGSM, "coordinate");
				}else{
					
					String[] info =bamFile.split("_");
					String Name, Nucleotide, Barcode, Lane;
					Name = Nucleotide = Barcode = Lane = null;
					if(info.length >  3){
						Name = info[0];
						Nucleotide = info[1];
						Barcode = info[2];
						Lane = info[3];
					}

					//"Inter3-1_DNA_ATCACG_L005_R1_001fastq_Crubella_183.strict.sam"
					String RGLB = Functions.getValue(T, "-RGLB", Lane);
					String RGPL = Functions.getValue(T, "-RGPL", "illumina");
					String RGPU = Functions.getValue(T, "-RGPU", Barcode);
					String RGSM = Functions.getValue(T, "-RGSM", Name);
					String RGID = Functions.getValue(T, "-RGID", RGSM+"_"+RGPU+"_"+RGLB+"_"+Nucleotide+"_"+Lane);

				
					
					prepareForGATK(EW, bamFile, memory, jarDir, suffix, RGID, RGLB, RGPL, RGPU, RGSM, "coordinate");

				}


			}


			EW.println();
			EW.println();
			EW.println("wait");

			EW.flush();
			EW.close();

		}catch(Exception E){E.printStackTrace();}
	}




	public static String prepareForGATK( ExtendedWriter EW, String bamFile, int memory, String jarDir,String suffix, String RGID, String RGLB, String RGPL,String RGPU,String RGSM,String sortOrder){
		//		System.out.println("Asumes that the sam filename has the following order NAME_[RNA|DNA]_BARCODE_LANE_WHATEVER");
		//		System.out.println(bamFile);
		//"Inter3-1_DNA_ATCACG_L005_R1_001fastq_Crubella_183.strict.sam"
		
		if(suffix.indexOf("sam")>-1){
			bamFile = SamtoolsSBATCH.sam2bam(EW, bamFile,-1, -1, false, true, true, false, false);
		}
		
		bamFile = AddOrReplaceReadGroups(EW, bamFile, memory,jarDir,suffix,RGID,RGLB,RGPL,RGPU,RGSM,sortOrder)+"."+suffix;
		bamFile = markDuplicates(EW, bamFile, memory, jarDir,suffix);

		System.out.println("Name after picardSteps "+ bamFile);
		return bamFile;

	}



	public static String  AddOrReplaceReadGroups( ExtendedWriter EW, String inFile1, int memory, String jarDir,String suffix,
			String RGID,String RGLB,String RGPL,String RGPU,String RGSM,String sortOrder){
		
		String baseFile =inFile1.substring(0,inFile1.lastIndexOf(suffix)-1)+".AddOrReplaceReadGroups";
		EW.println();


		String bowtiecommand = "java -Xmx"+memory+"G -jar "+jarDir+"/AddOrReplaceReadGroups.jar INPUT="+inFile1+" OUTPUT="+baseFile+"."+suffix+
				" RGID="+RGID+" "+ //String	Read Group ID Default value: 1. This option can be set to 'null' to clear the default value.
				"RGLB="+RGLB+" "+ //String	Read Group Library Required. Insert size
				"RGPL="+RGPL+" "+ //String	Read Group platform (e.g. Illumina, solid) Required.
				"RGPU="+RGPU+" "+ //String	Read Group platform unit (eg. run barcode) Required.
				"RGSM="+RGSM+" "+ //String	Read Group sample name Required.
				"VALIDATION_STRINGENCY=LENIENT SORT_ORDER="+sortOrder; //Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT. Default value: null. Possible values: {unsorted, queryname, coordinate}



		System.out.println(bowtiecommand);
		EW.println(bowtiecommand);
		EW.println();
		return  baseFile;

	}


	public static String markDuplicates( ExtendedWriter EW, String inFile1, int memory, String jarDir,String suffix){

		String baseFile =inFile1.substring(0,inFile1.lastIndexOf(suffix)-1)+".markDuplicates";
		EW.println();
		String bowtiecommand = "java -Xmx"+memory+"G -jar "+jarDir+"/MarkDuplicates.jar INPUT="+inFile1+" OUTPUT="+baseFile+"."+suffix+" METRICS_FILE="+baseFile+".metrics VALIDATION_STRINGENCY=LENIENT";
		EW.println(bowtiecommand);
		EW.println();
		return baseFile+"."+suffix;
	}



	public static String SortSam( ExtendedWriter EW, String inFile1, int memory, String jarDir,String suffix, String sortOrder){

		String baseFile =inFile1.substring(0,inFile1.lastIndexOf(suffix)-1)+".SortSam";
		EW.println();
		String bowtiecommand = "java -Xmx"+memory+"G -jar "+jarDir+"/SortSam.jar INPUT="+inFile1+" OUTPUT="+baseFile+"."+suffix+" SORT_ORDER="+sortOrder +" VALIDATION_STRINGENCY=LENIENT ";
		EW.println(bowtiecommand);
		EW.println();
		return baseFile;

	}

	public static void CreateSequenceDictionary(ExtendedWriter EW, String jarDir, int memory,  String inFile1, String outFile){

		EW.println();
		//java -jar CreateSequenceDictionary.jar R= Homo_sapiens_assembly18.fasta O= Homo_sapiens_assembly18.dict - See more at: http://gatkforums.broadinstitute.org/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference#latest

		String bowtiecommand = "java -Xmx"+memory+"G -jar "+jarDir+"/CreateSequenceDictionary.jar  R= "+inFile1+" O= "+outFile;
		EW.println(bowtiecommand);
		EW.println();

	}


	public static String mergeBamFiles(ExtendedWriter EW, int memory, String jarDir, ArrayList <String> fileNames, String mergedFileName){
		String commandLine = "java -Xmx"+memory+"G -jar "+jarDir+"/MergeSamFiles.jar OUTPUT="+mergedFileName+" ";
		for(int i = 0; i < fileNames.size(); i++){
			commandLine += "INPUT="+fileNames.get(i)+" ";
		}
		commandLine += " VALIDATION_STRINGENCY=LENIENT ";

		EW.println(commandLine);	
		EW.println();

		EW.flush();
		return mergedFileName;


	}

}

