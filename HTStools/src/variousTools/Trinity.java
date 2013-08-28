package variousTools;

import general.ExtendedWriter;
import general.Functions;
import general.IOTools;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import Sequence.FastQSequences;

public class Trinity {

	String time;
	String projectDir;

	String suffix;
	String split;
	String[] sep; 
	
	String SS_lib_type;

	boolean pairedEnd;

	public Trinity(){
		this.split = ".";

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
		Trinity trinity = new Trinity();
		trinity.run(T);
	}

	public void run(Hashtable<String,String> T){

		if(T.containsKey("-deNovo")) return;
		
		
		String inDir, outDir;
		inDir = outDir =  null;
		boolean allPresent = true;
		this.SS_lib_type = null;
		String timeStamp = Functions.getValue(T, "-TS",Functions.getDateTime());
		SBATCHinfo sbatch = new SBATCHinfo();
		if(!sbatch.addSBATCHinfo(T)) allPresent = false;

		if(T.containsKey("-i"))
			inDir= Functions.getValue(T, "-i", ".");
		else{
			System.out.println("must contain inDirectory -i");
			allPresent = false;
		}

		boolean single = false;
		if(T.containsKey("-single")) single = true;
		if(T.containsKey("-SS_lib_type"))
			{this.SS_lib_type = Functions.getValue(T, "-SS_lib_type");
				System.out.println("--SS_lib_type works");
			}

		if(T.containsKey("-singleEnd")){
			this.pairedEnd = false;
		}
		else 
			this.pairedEnd = true;
		projectDir= Functions.getValue(T, "-pDir", IOTools.getCurrentPath());

		if(T.containsKey("-o"))
			outDir= Functions.getValue(T, "-o", ".");
		else{
			System.out.println("must contain outDirectory -o");
			allPresent = false;
		}

		time = Functions.getValue(T, "-time", "00.00.00");
		suffix = Functions.getValue(T,"-suffix","fastq");

		this.sep = new String[2];
		sep[0] = "1."+suffix;
		sep[1] = "2."+suffix;
		int insertSize = Integer.parseInt(Functions.getValue(T, "-insertSize", "400"));

		if(projectDir.compareTo("notPresent") != 0){
			inDir = projectDir+"/"+inDir;
			outDir = projectDir+"/"+outDir;
		}
		else{
			projectDir = inDir;
		}


		if(allPresent)
			if(!single)
				trinityScript(sbatch, timeStamp, inDir, outDir, pairedEnd, insertSize);
			else
				trinityScriptSingle(sbatch, timeStamp, inDir, outDir, pairedEnd, insertSize);
		else
			System.out.println("\n\nAborting run because of missing arguments for trinity.");
	}



	public void trinityScript(SBATCHinfo sbatch, String timeStamp, String inDir, String outDir, boolean PE, int instertSize){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			if(PE){
				ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_PE_trinity.sh"));
				trinityDirPE(EW,sbatch, timeStamp, inDir, outDir,instertSize);
				EW.flush();
				EW.close();

			}else{
				ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_SR_trinity.sh"));
				trinityDirSR(EW,sbatch, timeStamp, inDir, outDir);
				EW.flush();
				EW.close();
			}

		}catch(Exception E){E.printStackTrace();}
	} 

	public void trinityScriptSingle(SBATCHinfo sbatch, String timeStamp, String inDir, String outDir, boolean PE, int instertSize){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			if(PE){
				ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_Single_PE_trinity.sh"));
				trinityDirPE(EW,sbatch, timeStamp, inDir, outDir,instertSize);
				EW.flush();
				EW.close();
			}else{
				ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_Single_SR_trinity.sh"));
				trinityDirSRSingle(EW,sbatch, timeStamp, inDir, outDir);
				EW.flush();
				EW.close();
			}
		}catch(Exception E){E.printStackTrace();}
	} 



	public void trinityFilePE( ExtendedWriter EW , ArrayList <String[]> pairs, String outDir, int insertSize, String memory,int CPUs){
		// Trinity.pl --seqType fq --left 13_2/QC/13_2_1.fastq  --right 13_2/QC/13_2_2.fastq --output /bubo/proj/b2011098/private/trinity/13_2/output  --min_contig_length 200 --CPU 10 --bflyHeapSpace 10000M
		String catFiles = "cat ";
		for(int i = 0; i < pairs.size()-1;i++)
			catFiles +=pairs.get(i)[0]+" ";
		catFiles +=pairs.get(pairs.size()-1)[0];
		catFiles +=" >"+outDir+"/left.fq";
		EW.println(catFiles);
		catFiles = "cat ";
		for(int i = 0; i < pairs.size()-1;i++)
			catFiles +=pairs.get(i)[1]+" ";
		catFiles +=pairs.get(pairs.size()-1)[1];
		catFiles +=" >"+outDir+"/right.fq";
		EW.println(catFiles);
		EW.println("cd "+outDir);
		String trinitycommand = "Trinity.pl --seqType fq --left left.fq --right right.fq ";
		if(this.SS_lib_type != null){
			trinitycommand += " --SS_lib_type "+this.SS_lib_type+" ";
		}

		trinityFileEND(EW,trinitycommand,outDir, memory,CPUs);
	}

	public void trinityFileSR( ExtendedWriter EW , ArrayList <String> files, String outDir, String memory, int CPUs){
		// Trinity.pl --seqType fq --left 13_2/QC/13_2_1.fastq  --right 13_2/QC/13_2_2.fastq --output /bubo/proj/b2011098/private/trinity/13_2/output  --min_contig_length 200 --CPU 10 --bflyHeapSpace 10000M
		String catFiles = "cat ";
		for(int i = 0; i < files.size()-1;i++)
			catFiles +=files.get(i)+" ";
		catFiles +=files.get(files.size()-1);
		String trinitycommand = "Trinity.pl ";
		if(suffix.indexOf("fastq") == 0){
			catFiles +=" >"+outDir+"/single.fq";
			EW.println(catFiles);
			EW.println("cd "+outDir);
			trinitycommand += "--seqType fq --single single.fq";
			if(this.SS_lib_type != null){
				trinitycommand += " --SS_lib_type "+this.SS_lib_type+" ";
			}
			
		}
		else{
			catFiles +=" >"+outDir+"/single.fa";
			EW.println(catFiles);
			EW.println("cd "+outDir);
			trinitycommand += "--seqType fa --single single.fa";
		}

		trinityFileEND(EW,trinitycommand,outDir, memory,CPUs);
	}


	public void trinityFileSingle( ExtendedWriter EW , String file, String inDir, String outDir){
		// Trinity.pl --seqType fq --left 13_2/QC/13_2_1.fastq  --right 13_2/QC/13_2_2.fastq --output /bubo/proj/b2011098/private/trinity/13_2/output  --min_contig_length 200 --CPU 10 --bflyHeapSpace 10000M
		String trinitycommand = "Trinity.pl ";
		if(suffix.indexOf("fastq") == 0){
			trinitycommand += "--seqType fq --single "+inDir+"/"+file;
		}
		else{
			trinitycommand += "--seqType fa --single "+inDir+"/"+file;
		}
		if(this.SS_lib_type != null){
			trinitycommand += " --SS_lib_type "+this.SS_lib_type+" ";
		}
		
		trinityFileEndSingle(EW,trinitycommand,outDir+"/"+file.substring(0, file.indexOf(this.suffix)-1));
	}

	public static String getTime(int nrOfSequences){
		
		int nrofMB = nrOfSequences/1000000+1;
		int hours = nrofMB*3;
		if(hours > 168)hours = 168;	
		//System.out.println("Number of sequences are:"+ nrOfSequences);
		String newTime = hours+":00:00";
		return newTime;
		
	}
	
	public static String trinityFileStart(int nrOfSequences,String inDir, String outDir, String timestamp, int count, SBATCHinfo sbatch,ExtendedWriter EW , String time){

		int nrofMB = nrOfSequences/1000000+1;
		int hours = nrofMB*3;
		if(hours > 168)hours = 168;	
		System.out.println("Number of sequences are:"+ nrOfSequences);
		String newTime = "";
		if(time.compareTo("00.00.00") == 0){
			System.out.println("time not set will be set to "+hours +" hours");
			newTime =hours+":00:00";
		}
		else{
			newTime = time;
		}


		String[] split1 = outDir.split("/");
		String memory = "3G";
		if(nrofMB < 3){
			System.out.println("Memory allocated will be 3 GB");
			sbatch.printSBATCHinfoCore(EW,outDir,timestamp,count,"trinity_"+split1[split1.length-1], newTime);
			memory = "2G";
		}
		else if(nrofMB < 24){
			System.out.println("Memory allocated will be 24 GB");
			sbatch.printSBATCHinfo(EW,outDir,timestamp,count,"trinity_"+split1[split1.length-1], newTime);
			memory = "23G";
		}
		else if(nrofMB < 36){
			System.out.println("Memory allocated will be 36 GB");
			sbatch.printSBATCHinfoFat(EW,outDir,timestamp,count,"trinity_"+split1[split1.length-1], newTime);
			memory = "34G";
		}
		else if(nrofMB < 72){
			System.out.println("Memory allocated will be 72 GB");
			sbatch.printSBATCHinfo72GB(EW,outDir,timestamp,count,"trinity_"+split1[split1.length-1], newTime);
			memory = "70G";
		}				
		else{
			sbatch.printSBATCHinfohalvan(EW,outDir,timestamp,count,"trinity_"+split1[split1.length-1], newTime,nrofMB);
			memory = "72G";
		}

		EW.println("module load bioinfo-tools");
		EW.println("module load bowtie");
		EW.println("module load trinity/2012-04-27");

		EW.println("ulimit -s unlimited");

		EW.println("cd "+inDir);


		EW.println();
		EW.println();

		return memory;
	}
	
	
	
	public void trinityFileStartSingle(String inDir, String outDir, String timestamp, int count, SBATCHinfo sbatch,ExtendedWriter EW ){

		int hours = 20;
		String newTime = "";
		if(time.compareTo("00.00.00") == 0){
			System.out.println("time not set will be set to "+hours +" hours");
			newTime =hours+":00:00";
		}
		else{
			newTime = time;
		}

		sbatch.printSBATCHinfoCore(EW,outDir,timestamp,count,"trinity_"+count, newTime);
		
		EW.println("module load bioinfo-tools");
		EW.println("module load bowtie");
		EW.println("module load trinity/2012-04-27");

		EW.println("ulimit -s unlimited");

		EW.println("cd "+inDir);


		EW.println();
		EW.println();



	}


	public void trinityFileEND( ExtendedWriter EW ,String trinitycommand, String outDir,String memory,int CPUs){
		trinitycommand +=" --output "+outDir;
		trinitycommand +=" --CPU "+(CPUs-1)+" ";
		trinitycommand +=" --kmer_method jellyfish --max_memory "+memory;
		EW.println("echo START");
		EW.println();
		//EW.println("echo \""+trinitycommand+"\" 1>&2");
		EW.println(trinitycommand);
		EW.println();
		EW.println("echo DONE");
	}

	public void trinityFileEndSingle( ExtendedWriter EW ,String trinitycommand, String outDir){
		trinitycommand +=" --output "+outDir;
		trinitycommand +=" --min_contig_length 100 --CPU 1 --bflyHeapSpaceMax 2G";
		trinitycommand +=" --kmer_method jellyfish --max_memory 2G";
		EW.println("echo START");
		EW.println();
		//EW.println("echo \""+trinitycommand+"\" 1>&2");
		EW.println(trinitycommand);
		EW.println();
		EW.println("echo DONE");
	}

	public void trinityDirPE( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String inDir, String outDir, int insertSize ){

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
			int count = 0;
			try{
				ArrayList <String[]> pairs = IOTools.findPairs(fileNames,sep);
				int nrOfSequences = FastQSequences.countSequencesPE(inDir,pairs);
				String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_"+count+"_trinity.sbatch";
				generalSbatchScript.println("sbatch "+ sbatchFileName);

				ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));

				String memory = trinityFileStart(nrOfSequences,  inDir,  outDir,  timestamp,  count,  sbatch, EW,time);
				int CPUs = 8;
				trinityFilePE(EW, pairs,outDir,insertSize, memory,CPUs);

				EW.flush();
				EW.close();

			}catch(Exception E){E.printStackTrace();}

			count++;

		}
		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
		for(int i = 0; i < subDirs.size(); i++){
			trinityDirPE(generalSbatchScript,sbatch,timestamp,inDir+"/"+subDirs.get(i),outDir+"/"+subDirs.get(i), insertSize);
		}
	}


	public void trinityDirSR(ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String inDir, String outDir ){

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
			int count = 0;
			try{
				String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_"+count+"_trinity.sbatch";
				generalSbatchScript.println("sbatch "+ sbatchFileName);
				ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));

				int nrOfSequences = FastQSequences.countSequencesSR(inDir,fileNames);

				String memory = trinityFileStart(nrOfSequences,  inDir,  outDir,  timestamp,  count,  sbatch, EW,time );

				EW.println();
				EW.println();
				trinityFileSR(EW, fileNames,outDir, memory,8);

				EW.flush();
				EW.close();

			}catch(Exception E){E.printStackTrace();}

			count++;

		}
		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
		for(int i = 0; i < subDirs.size(); i++){
			trinityDirSR(generalSbatchScript,sbatch,timestamp,inDir+"/"+subDirs.get(i),outDir+"/"+subDirs.get(i));
		}
	}

	public void trinityDirSRSingle( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String inDir, String outDir ){

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
			int count = 0;
			int number = 0;
			int perRun = 100;
			while(count < fileNames.size()/perRun+1){ 
				try{
					String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_"+count+"_trinity.sbatch";
					generalSbatchScript.println("sbatch "+ sbatchFileName);
					
					ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));

					trinityFileStartSingle(inDir,  outDir,  timestamp,  count,  sbatch, EW);

					for(int i = 0; i < perRun;i++){
						if(number < fileNames.size()){
							EW.println();
							trinityFileSingle(EW, fileNames.get(number), inDir,  outDir);
							number++;
						}
					}
					EW.flush();
					EW.close();

				}catch(Exception E){E.printStackTrace();}

				count++;
			}
		}
	}
}
