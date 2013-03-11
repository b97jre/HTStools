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

public class deNovoAssembly {

	String time;
	String projectDir;

	String suffix;
	String split;
	String[] sep; 

	//
	int insertSize;
	int distribution;

	// trinity specific things
	String SS_lib_type;



	// oasesSpecificThings
	int maxKmer;
	int minKmer;
	int finalKmer;
	int step;

	boolean trinity;
	boolean oases;


	boolean pairedEnd;
	boolean strandSpecific;

	public deNovoAssembly(){
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
		deNovoAssembly assembly = new deNovoAssembly();
		assembly.run(T);
	}



	public void run(Hashtable<String,String> T){

		String inDir, trinityDir, oasesDir;
		inDir = trinityDir = oasesDir=  null;
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

		if(T.containsKey("-SS_lib_type")){
			this.SS_lib_type = Functions.getValue(T, "-SS_lib_type");
			System.out.println("--SS_lib_type works");
			this.strandSpecific = true;
		}
		if(T.containsKey("-strand_specific")){
			this.SS_lib_type = Functions.getValue(T, "-SS_lib_type", "forward");
			System.out.println("--SS_lib_type works");
			this.strandSpecific = true;
		}

		this.pairedEnd = true;
		if(T.containsKey("-singleEnd")){
			this.pairedEnd = false;
		}


		projectDir= Functions.getValue(T, "-pDir", IOTools.getCurrentPath());

		trinityDir= Functions.getValue(T, "-trinity", inDir+"_trinity");
		oasesDir= Functions.getValue(T, "-oases", inDir+"_oases");


		
		if(T.containsKey("-trinity"))this.trinity=true;
		else this.trinity = false;
		if(T.containsKey("-oases"))this.oases=true;
		else this.oases = false;

		time = Functions.getValue(T, "-time", "00.00.00");
		suffix = Functions.getValue(T,"-suffix","fastq");
		this.sep = new String[2];
		sep[0] = "1."+suffix;
		sep[1] = "2."+suffix;
		this.insertSize = Integer.parseInt(Functions.getValue(T, "-insertSize", "200"));
		
		this.finalKmer=Functions.getInt(T, "-g", 27);
		this.minKmer=Functions.getInt(T, "-m", 21);
		this.maxKmer=Functions.getInt(T, "-M", 29);
		this.step=Functions.getInt(T, "-s", 2);
		this.insertSize = Functions.getInt(T, "-insert_size", 200);

		if(allPresent){
			
		
			deNovoScript(sbatch, timeStamp, projectDir+"/"+inDir, projectDir+"/"+trinityDir,projectDir+"/"+oasesDir);
		}
		else
			System.out.println("\n\nAborting run because of missing arguments for trinity.");
	}



	public void deNovoScript(SBATCHinfo sbatch, String timeStamp, String inDir, String trinityDir, String oasesDir){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			if(this.pairedEnd){
				ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_PE_trinity.sh"));
				deNovoDirPE(EW,sbatch, timeStamp, inDir, trinityDir, oasesDir);
				EW.flush();
				EW.close();
			}else{
				ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_SR_trinity.sh"));
				deNovoDirSR(EW,sbatch, timeStamp, inDir, trinityDir, oasesDir);
				EW.flush();
				EW.close();
			}

		}catch(Exception E){E.printStackTrace();}
	} 






	public void deNovoDirPE( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String inDir, String trinityDir, String oasesDir){

		ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir,suffix);
		ArrayList <String[]> pairs = IOTools.findPairs(fileNames,sep);
		if(!pairs.isEmpty() && (oases || trinity)){
			System.out.println("creating de novo files for paired reads in "+inDir);
			int count = 0;
			int nrOfSequences = FastQSequences.countSequencesPE(inDir,pairs);
			if(trinity){
				try{
					if(!IOTools.isDir(trinityDir))
						IOTools.mkDirs(trinityDir);
					if(!IOTools.isDir(trinityDir+"/reports"))
						IOTools.mkDir(trinityDir+"/reports");
					if(!IOTools.isDir(trinityDir+"/scripts"))
						IOTools.mkDir(trinityDir+"/scripts");

					//trinitySpecific
					String 	sbatchFileName = trinityDir+"/scripts/"+timestamp+"_"+count+"_trinity.sbatch";
					generalSbatchScript.println("sbatch "+ sbatchFileName);
					ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
					String time = Trinity.getTime(nrOfSequences);
					String memory = Trinity.trinityFileStart(nrOfSequences,  inDir,  trinityDir,  timestamp,  count,  sbatch, EW,time);
					int CPUs = 8;
					Trinity run = new Trinity();
					run.suffix = this.suffix;
					
					run.trinityFilePE(EW, pairs,trinityDir,insertSize, memory,CPUs);
					EW.flush();
					EW.close();

				}catch(Exception E){E.printStackTrace();}
			}
			if(oases){
				try{
					if(!IOTools.isDir(oasesDir))
						IOTools.mkDirs(oasesDir);
					if(!IOTools.isDir(oasesDir+"/reports"))
						IOTools.mkDir(oasesDir+"/reports");
					if(!IOTools.isDir(oasesDir+"/scripts"))
						IOTools.mkDir(oasesDir+"/scripts");

					String 	sbatchFileName = oasesDir+"/scripts/"+timestamp+"_"+count+"_oases.sbatch";
					generalSbatchScript.println("sbatch "+ sbatchFileName);
					ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
					String time = Oases.getTime(nrOfSequences);
					String memory = Oases.FileStart(nrOfSequences,  inDir,  trinityDir,  timestamp,  count,  sbatch, EW,time);
					int CPUs = 8;
					Oases run = new Oases();
					run.setValues(maxKmer, minKmer, step, finalKmer, suffix, insertSize, projectDir, strandSpecific);
					
					run.oasesFilePE(EW, pairs, oasesDir, memory, CPUs);
					
					FastaSequences.getORFs(oasesDir+"/oasesPipelineMerged/transcripts.fa");
					
					
					EW.flush();
					EW.close();

				}catch(Exception E){E.printStackTrace();}
			}
			count++;
		}
		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
		for(int i = 0; i < subDirs.size(); i++){
			deNovoDirPE(generalSbatchScript,sbatch,timestamp,inDir+"/"+subDirs.get(i),trinityDir+"/"+subDirs.get(i),oasesDir+"/"+subDirs.get(i));
		}
	}


	public void deNovoDirSR(ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String inDir, String trinityDir, String oasesDir ){
		System.out.println("creating de novo files for single reads in "+inDir);

		ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir,suffix);
		if(!fileNames.isEmpty() && (oases || trinity)){
			int count = 0;
			int nrOfSequences = FastQSequences.countSequencesSR(inDir,fileNames);
			System.out.println("creating de novo files for single reads in "+inDir);
			if(trinity){
				try{
					if(!IOTools.isDir(trinityDir))
						IOTools.mkDirs(trinityDir);
					if(!IOTools.isDir(trinityDir+"/reports"))
						IOTools.mkDir(trinityDir+"/reports");
					if(!IOTools.isDir(trinityDir+"/scripts"))
						IOTools.mkDir(trinityDir+"/scripts");

					//trinitySpecific
					String 	sbatchFileName = trinityDir+"/scripts/"+timestamp+"_"+count+"_trinity.sbatch";
					generalSbatchScript.println("sbatch "+ sbatchFileName);
					ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
					String time = Trinity.getTime(nrOfSequences);
					String memory = Trinity.trinityFileStart(nrOfSequences,  inDir,  trinityDir,  timestamp,  count,  sbatch, EW,time);
					int CPUs = 8;
					Trinity run = new Trinity();
					run.suffix = this.suffix;
					run.trinityFileSR(EW, fileNames,trinityDir, memory,CPUs);
					EW.flush();
					EW.close();

				}catch(Exception E){E.printStackTrace();}
			}
			if(oases){
				try{
					if(!IOTools.isDir(oasesDir))
						IOTools.mkDirs(oasesDir);
					if(!IOTools.isDir(oasesDir+"/reports"))
						IOTools.mkDir(oasesDir+"/reports");
					if(!IOTools.isDir(oasesDir+"/scripts"))
						IOTools.mkDir(oasesDir+"/scripts");

					String 	sbatchFileName = trinityDir+"/scripts/"+timestamp+"_"+count+"_oases.sbatch";
					generalSbatchScript.println("sbatch "+ sbatchFileName);
					ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
					String time = Oases.getTime(nrOfSequences);
					String memory = Oases.FileStart(nrOfSequences,  inDir,  trinityDir,  timestamp,  count,  sbatch, EW,time);
					int CPUs = 8;
					Oases run = new Oases();
					run.setValues(maxKmer, minKmer, step, finalKmer, suffix, insertSize, projectDir, strandSpecific);
					run.oasesFileSR(EW, fileNames, oasesDir, memory, CPUs);
					EW.flush();
					EW.close();

				}catch(Exception E){E.printStackTrace();}
			}
			count++;
		}
		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
		for(int i = 0; i < subDirs.size(); i++){
			deNovoDirSR(generalSbatchScript,sbatch,timestamp,inDir+"/"+subDirs.get(i),trinityDir+"/"+subDirs.get(i),trinityDir+"/"+subDirs.get(i));
		}
	}

}
