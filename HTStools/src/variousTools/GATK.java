
// Asumes that the non-GATK steps have been made. For those steps please see Picard.java.
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
import java.util.Map.Entry;

import Sequence.FastaSequences;

public class GATK {

	String time;
	String projectDir;
	String suffix;
	String prefix;
	String picardDir;
	String GATKdir;
	String knownSNPVCF;
	String knownIndelVCF;
	String targetIntervals;
	int memory;
	int cutoff;
	int padding;
	String job;



	boolean PRIORITIZE;
	boolean HaploTypeCaller;
	boolean rerun;
	String Reference;




	public GATK(){
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
		GATK GATK = new GATK();
		GATK.run(T);
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


		if(T.containsKey("-R"))
			this.Reference= Functions.getValue(T, "-R");
		else{
			System.out.println("must contain ReferenceFile -R ");
			allPresent = false;
		}
		picardDir= Functions.getValue(T, "-picardDir", "/bubo/home/h17/johanr/bin");
		GATKdir= Functions.getValue(T, "-GATKDir", "/bubo/home/h17/johanr/bin");
		time = Functions.getValue(T, "-time", "3:00:00");
		cutoff = Functions.getInt(T, "-L", -1);


		if(!T.containsKey("-time"))
			System.out.println("must contain likely time -time now set to default 3:00:00");


		if(T.containsKey("-UnifiedGenotyper"))this.HaploTypeCaller=false;
		else this.HaploTypeCaller=true;


		knownSNPVCF = Functions.getValue(T,"-knownSNPs", null);
		if(IOTools.fileExists(IOTools.getCurrentPath()+"/"+knownSNPVCF))
			knownSNPVCF=IOTools.getCurrentPath()+"/"+knownSNPVCF;



		knownIndelVCF = Functions.getValue(T,"-knownIndels", null);
		targetIntervals = Functions.getValue(T,"-targetIntervals", null);
		
		
		if(IOTools.fileExists(IOTools.getCurrentPath()+"/"+knownIndelVCF))
			knownIndelVCF=IOTools.getCurrentPath()+"/"+knownIndelVCF;

		this.PRIORITIZE=true;
		if(T.containsKey("-UNIQUIFY"))PRIORITIZE=false;

		if(T.containsKey("-rerun"))rerun=true;

		memory = Functions.getInt(T, "-X", 20);
		if(memory < 3 )this.job="core";
		else if(memory< 24) this.job = "thin";
		else this.job = "fat";

		this.padding = 100;

		suffix = Functions.getValue(T,"-suffix","bam");
		prefix = Functions.getValue(T,"-prefix");
		if(suffix.indexOf('.')== 0)suffix = suffix.substring(1);

		projectDir = Functions.getValue(T,"-projectDir",IOTools.getCurrentPath());


		if(allPresent){
			if(IOTools.fileExists(projectDir+"/"+this.Reference)){
				this.Reference = projectDir+"/"+this.Reference;
			}
			if(IOTools.isDir(projectDir+"/"+inDir)){
				inDir = projectDir+"/"+inDir;
			}
			else if(!IOTools.isDir(inDir)){
				System.out.println("Neither "+ inDir +" nor "+projectDir+"/"+inDir+ "was found");
				return;
			}
			GATKInitial(sbatch, timeStamp, inDir,T);
		}
		else{
			System.out.println("\n\nAborting run because of missing arguments for GATK.");
			System.out.println("\n\nTypical line should be something like:");
			System.out.println("\n\nThis is the current :");
			System.out.println("HTStools -p sbatch -GATK -phase1 -ReassignOneMappingQuality  -ReAlign -BQSR -BQSRprint -ReduceReads -suffix markDuplicates.bam -i bamFiles -R /bubo/home/h17/johanr/capsellaDir/nobackup/Cr_reference/assembly/Crubella_183.fa -time 40:00:00 -knownSNPs /bubo/home/h17/johanr/capsellaDir/nobackup/knownSNPs/Cg/GATK_13samps.flt_excluded_sites.vcf -targetIntervals /bubo/home/h17/johanr/capsellaDir/nobackup/Cr_reference/annotation/Crubella_183_only_exons_unique.bed");
			System.out.println("HTStools -p sbatch -GATK -phase1  -i INDIR -R REFERENCEFILE -time 4:00:00 -knownSNPs knownSNPsFile  [-knownIndels knownIndelsFile] [-targetIntervals bedFile] -suffix bam [ -preGATK -RGPU Barcode] [-RGSM Name]] ");
			System.out.println("HTStools -p sbatch -GATK -phase2 -i INDIR -R REFERENCEFILE -time 4:00:00 -L 100000 [-rerun] -X 2 -suffix reduced.bam");
			System.out.println("HTStools -p sbatch -GATK -merge -i INDIR -R REFERENCEFILE -time 6-00:00:00 -X 23 -suffix vcf -prefix output.raw.snps.indels");
			System.out.println("HTStools -p sbatch -GATK --phaseSNPs -i INDIR -R REFERENCEFILE -time 6-00:00:00 -X 23 -knownSNPs knownSNPsFile -suffix real.BQSR.bam ");
		}
	}





	public void GATKInitial(SBATCHinfo sbatch, String timeStamp, String inDir,Hashtable<String,String> T){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			String lastDir = inDir.substring(inDir.lastIndexOf("/")+1); 
			if(lastDir.length() < 1){
				lastDir = inDir.substring(0,inDir.lastIndexOf("/")); 
				lastDir = lastDir.substring(lastDir.lastIndexOf("/")+1); 
			}
			if(T.containsKey("-phase1")){
				ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"."+lastDir+".GATK.phase1.sh"));
				GATKDirPhase1(EW,sbatch, timeStamp, inDir, T);
				EW.flush();
				EW.close();
			}
			if(T.containsKey("-phase2")){
				ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"."+lastDir+".GATK.phase2.sh"));
				GATKPhase2(EW,sbatch, timeStamp, inDir, T);
				EW.flush();
				EW.close();
			}
			if(T.containsKey("-merge")){

				ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"."+lastDir+".GATK.merge.sh"));
				GATKMerge2(EW,sbatch, timeStamp, inDir, T);
				EW.flush();
				EW.close();
			}

			
			if(T.containsKey("-phaseSNPs")){

				ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"."+lastDir+".GATK.phaseSNPs.sh"));
				GATKPhaseSNPs(EW,sbatch, timeStamp, inDir, T);
				EW.flush();
				EW.close();
			}
			
			if(T.containsKey("-Genotype")){

				ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"."+lastDir+".GATK.GATKGenotypeAndValidate.sh"));
				GATKGenotype(EW,sbatch, timeStamp, inDir, T);
				EW.flush();
				EW.close();
			}
			
			
			//			if(T.containsKey("-phase3")){
			//				GATKPhase2(sbatch, timeStamp, inDir, T);
			//			}


		}
		catch(Exception E){E.printStackTrace();}
	} 


	public void GATKPhaseSNPs( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String inDir,Hashtable<String,String> T){

		ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir,suffix);
		if(!fileNames.isEmpty()){
			try{
				GATKPhaseSNPsSample(generalSbatchScript,sbatch,timestamp,inDir,fileNames,T);
			}catch(Exception E){E.printStackTrace();}
		}
		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
		for(int i = 0; i < subDirs.size(); i++){
			GATKPhaseSNPs(generalSbatchScript,sbatch,timestamp,inDir+"/"+subDirs.get(i),T);
		}
	}

	public void GATKGenotype( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String inDir,Hashtable<String,String> T){

		ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir,suffix);
		if(!fileNames.isEmpty()){
			try{
				GATKGenotypeSample(generalSbatchScript,sbatch,timestamp,inDir,fileNames,T);
			}catch(Exception E){E.printStackTrace();}
		}
		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
		for(int i = 0; i < subDirs.size(); i++){
			GATKGenotype(generalSbatchScript,sbatch,timestamp,inDir+"/"+subDirs.get(i),T);
		}
	}
	
	

	public void GATKDirPhase1( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String inDir,Hashtable<String,String> T){

		ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir,suffix);
		if(!fileNames.isEmpty()){
			try{
				GATKPhase1Sample(generalSbatchScript,sbatch,timestamp,inDir,fileNames,T);
			}catch(Exception E){E.printStackTrace();}
		}
		ArrayList <String> subDirs = IOTools.getDirectories(inDir);
		for(int i = 0; i < subDirs.size(); i++){
			GATKDirPhase1(generalSbatchScript,sbatch,timestamp,inDir+"/"+subDirs.get(i),T);
		}
	}

	public ArrayList <String> getFilesWithSuffix(String dir){
		ArrayList <String> fileNames = new ArrayList <String>(); 
		ArrayList <String> subDirs = IOTools.getDirectories(dir);
		for(int i = 0; i < subDirs.size(); i++){
			ArrayList <String> newFileNames = getFilesWithSuffix(dir+"/"+subDirs.get(i));
			if(!newFileNames.isEmpty()){
				fileNames = Functions.mergeLists(fileNames, newFileNames);
			}
		}
		ArrayList <String> newFileNames = IOTools.getSequenceFilesFullPath(dir,suffix);
		if(!newFileNames.isEmpty()){
			fileNames = Functions.mergeLists(fileNames, newFileNames);
		}

		return fileNames;
	}


	public static void printSBATCHscript( ExtendedWriter MasterShellScriptFile, SBATCHinfo sbatch, String timestamp,int count, String outDir,String commandLine, String job, String outFile,String function,String time){
		try{
			String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_"+count+"_"+function+".sbatch";

			MasterShellScriptFile.println("sbatch "+ sbatchFileName);

			ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
			if(job.compareTo("thin") == 0)
				sbatch.printSBATCHinfo(EW,outDir,timestamp,count,function, time);
			else if(job.compareTo("core") == 0)
				sbatch.printSBATCHinfoCore(EW,outDir,timestamp,count,function, time);
			else if(job.compareTo("fat") == 0)
				sbatch.printSBATCHinfoFat(EW,outDir,timestamp,count,function, time);
			else
				sbatch.printSBATCHinfo(EW,outDir,timestamp,count,function, time);


			EW.println();
			EW.println();
			EW.println("cd "+outDir);
			EW.println();
			EW.println();
			EW.println(commandLine);
			EW.println();
			EW.println("mv "+outFile+".temp "+outFile);
			EW.println("mv "+outFile+".temp.idx "+outFile+".idx");
			EW.println();
			EW.println("wait");

			EW.flush();
			EW.close();

		}catch(Exception E){E.printStackTrace();}
	}


	public void GATKMerge2( ExtendedWriter MasterShellScriptFile, SBATCHinfo sbatch ,String timestamp,String outDir,Hashtable<String,String> T){

		/*				
		for each sample (asumed to be in same dir)

	    lanes.bam <- merged lane.bams for sample
	    dedup.bam <- MarkDuplicates(lanes.bam) 
	    realigned.bam <- realign(dedup.bam) [with known sites included if available]
	    recal.bam <- recal(realigned.bam)
	    sample.bam <- recal.bam			
		 */

		if(!IOTools.isDir(outDir+"/reports"))
			IOTools.mkDir(outDir+"/reports");
		if(!IOTools.isDir(outDir+"/scripts"))
			IOTools.mkDir(outDir+"/scripts");
		try{
			String prefix  = Functions.getValue(T,"-prefix");
			String suffix = Functions.getValue(T, "-suffix");
			
			
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+prefix+".Final.vcf"));
			ExtendedReader ER2 = new ExtendedReader(new FileReader(outDir+"/"+prefix+".part0."+suffix));
			while(ER2.lookAhead() == '#')EW.println(ER2.readLine());
			ER2.close();
			// Write header info once 

			ArrayList<Integer> lengths = FastaSequences.getLengths(Reference);
			ArrayList<String> names = FastaSequences.getNames(Reference);
			System.out.println();
			System.out.println();
			int count = 0;
			int pointer = 0;
			while(cutoff != -1 && pointer < names.size()){
				// ...

				int size = lengths.get(pointer);
				if(size > cutoff+padding){
					int start = 1;
					int stop = cutoff+padding;
					int first = 1;
					int last = cutoff;
					while(stop < size){
						
						ExtendedReader ER = new ExtendedReader(new FileReader(outDir+"/"+prefix+".part"+count+"."+suffix));
						while(ER.lookAhead() == '#')ER.readLine();
						while(ER.more()){
							String infoLine = ER.readLine(); 
							String[] info = infoLine.split("\t");
							int pos = Integer.parseInt(info[1]);
							if(pos >= first && pos <last ) EW.println(infoLine);
						}
						ER.close();
						start = stop-2*padding;
						stop = stop+cutoff;
						first=last;
						last=first+cutoff;
						count++;
						System.out.println(count);
					}
					stop = size;
					last=size+1;
					ExtendedReader ER = new ExtendedReader(new FileReader(outDir+"/"+prefix+".part"+count+"."+suffix));
					while(ER.lookAhead() == '#')ER.readLine();
					while(ER.more()){
						String infoLine = ER.readLine(); 
						String[] info = infoLine.split("\t");
						int pos = Integer.parseInt(info[1]);
						if(pos >= first && pos <last ) EW.println(infoLine);
					}
					ER.close();
					count++;
					System.out.println(count);
					pointer++;

				}else{
					pointer++;
					while(size < cutoff && pointer < names.size()){
						size += lengths.get(pointer);
						pointer++;
					}
					ExtendedReader ER = new ExtendedReader(new FileReader(outDir+"/"+prefix+".part"+count+"."+suffix));
					while(ER.lookAhead() == '#')ER.readLine();
					while(ER.more()){
						String infoLine = ER.readLine(); 
						EW.println(infoLine);
					}
					ER.close();
					count++;
					System.out.println(count);
				}

			}
			EW.flush();
			EW.close();
			ER2.close();

		}catch(Exception E){E.printStackTrace();}
	}



	public void GATKPhase2(  ExtendedWriter MasterShellScriptFile, SBATCHinfo sbatch ,String timestamp,String outDir,Hashtable<String,String> T){

		/*				
		for each sample (asumed to be in same dir)

	    lanes.bam <- merged lane.bams for sample
	    dedup.bam <- MarkDuplicates(lanes.bam) 
	    realigned.bam <- realign(dedup.bam) [with known sites included if available]
	    recal.bam <- recal(realigned.bam)
	    sample.bam <- recal.bam			
		 */

		if(!IOTools.isDir(outDir+"/reports"))
			IOTools.mkDir(outDir+"/reports");
		if(!IOTools.isDir(outDir+"/scripts"))
			IOTools.mkDir(outDir+"/scripts");
		try{
			ArrayList<Integer> lengths = FastaSequences.getLengths(Reference);
			ArrayList<String> names = FastaSequences.getNames(Reference);
			ArrayList <String> fileNames = getFilesWithSuffix(outDir);
			System.out.println("Files considered for final vcf file:");
			for(int i = 0; i < fileNames.size();i++){
				System.out.println(fileNames.get(i));
			}
			System.out.println();
			System.out.println();
			int count = 0;
			int pointer = 0;
			while(cutoff != -1 && pointer < names.size()){
				// ...

				int size = lengths.get(pointer);
				if(size > cutoff+padding){
					int start = 1;
					int stop = cutoff+padding;
					while(stop < size){
						ArrayList<String> SubNames = new ArrayList<String>();
						String interval = names.get(pointer)+":"+start+"-"+stop;
						SubNames.add(interval);
						String commandLine = "";
						String outFile = ""; 
						if(!HaploTypeCaller){
							outFile = "Unified.output.raw.snps.indels.part"+count+".vcf";
							commandLine = UnifiedGenotypeCaller(memory, GATKdir, fileNames, Reference, outFile+".temp", knownSNPVCF, SubNames);
						}
						else{
							outFile =  "Haplotype.output.raw.snps.indels.part"+count+".vcf";
							commandLine = HaplotypeCaller(memory, GATKdir, fileNames, Reference, outFile+".temp", knownSNPVCF, SubNames);

						}
						if(!rerun || !IOTools.fileExists(outDir+"/output.raw.snps.indels.part"+count+".vcf.idx"))
							printSBATCHscript(MasterShellScriptFile, sbatch, timestamp,count, outDir,commandLine,job,outFile, "GATK_phase2", time);
						start = stop-2*padding;
						stop = stop+cutoff;
						count++;
					}

					ArrayList<String> SubNames = new ArrayList<String>();
					stop = size;
					String interval = names.get(pointer)+":"+start+"-"+stop;
					SubNames.add(interval);
					String commandLine = "";
					String outFile = ""; 
					if(!HaploTypeCaller){
						outFile = "Unified.output.raw.snps.indels.part"+count+".vcf";
						commandLine = UnifiedGenotypeCaller(memory, GATKdir, fileNames, Reference, outFile+".temp", knownSNPVCF, SubNames);
					}
					else{
						outFile =  "Haplotype.output.raw.snps.indels.part"+count+".vcf";
						commandLine = HaplotypeCaller(memory, GATKdir, fileNames, Reference, outFile+".temp", knownSNPVCF, SubNames);

					}
					if(!rerun || !IOTools.fileExists(outDir+"/"+outFile))
						printSBATCHscript(MasterShellScriptFile, sbatch, timestamp,count, outDir,commandLine,job,outFile, "GATK_phase2", time);
					count++;
					pointer++;

				}else{
					ArrayList<String> SubNames = new ArrayList<String>();
					SubNames.add(names.get(pointer));
					pointer++;
					while(size < cutoff && pointer < names.size()){
						SubNames.add(names.get(pointer));
						size += lengths.get(pointer);
						pointer++;
					}
					String commandLine = "";
					String outFile = ""; 
					if(!HaploTypeCaller){
						outFile = "Unified.output.raw.snps.indels.part"+count+".vcf";
						commandLine = UnifiedGenotypeCaller(memory, GATKdir, fileNames, Reference, outFile+".temp", knownSNPVCF, SubNames);
					}
					else{
						outFile =  "Haplotype.output.raw.snps.indels.part"+count+".vcf";
						commandLine = HaplotypeCaller(memory, GATKdir, fileNames, Reference, outFile+".temp", knownSNPVCF, SubNames);

					}
					if(!rerun || !IOTools.fileExists(outDir+"/output.raw.snps.indels.part"+count+".vcf.idx"))
						printSBATCHscript(MasterShellScriptFile, sbatch, timestamp,count, outDir,commandLine,job,outFile, "GATK_phase2", time);
					count++;
				}

			}

			if(cutoff == -1){
				String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_"+count+"phase2_GATK.sbatch";
				System.out.println("sbatch "+ sbatchFileName);
				ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
				sbatch.printSBATCHinfo(EW,outDir,timestamp,0,"GATK_phase2", time);

				EW.println("cd "+outDir);

				if(HaploTypeCaller){
					String commandLine = HaplotypeCaller(memory, GATKdir, fileNames, Reference, "output.raw.snps.indels.vcf", knownSNPVCF, null);
					System.out.println(commandLine);
					EW.println(commandLine);
					count++;
				}
				else{
					System.out.println("I have not implemented UniGenotypeCaller"); 
				}
				EW.println();
				EW.println();
				EW.println("wait");

				EW.flush();
				EW.close();

			}


		}catch(Exception E){E.printStackTrace();}
	}

	public void GATKMerge(  ExtendedWriter MasterShellScriptFile, SBATCHinfo sbatch ,String timestamp,String outDir,Hashtable<String,String> T){

		/*				
		for each sample (asumed to be in same dir)

	    lanes.bam <- merged lane.bams for sample
	    dedup.bam <- MarkDuplicates(lanes.bam) 
	    realigned.bam <- realign(dedup.bam) [with known sites included if available]
	    recal.bam <- recal(realigned.bam)
	    sample.bam <- recal.bam			
		 */

		if(!IOTools.isDir(outDir+"/reports"))
			IOTools.mkDir(outDir+"/reports");
		if(!IOTools.isDir(outDir+"/scripts"))
			IOTools.mkDir(outDir+"/scripts");
		try{
			ArrayList <String> vcfFiles = IOTools.getSequenceFiles(outDir,suffix);
			vcfFiles = Functions.getStringsWithPrefix(vcfFiles,prefix);
			System.out.println("Files considered for final vcf file:");
			for(int i = 0; i < vcfFiles.size();i++){
				System.out.println(vcfFiles.get(i));
			}
			System.out.println();
			System.out.println();
			String commandLine = mergeVCFfiles(memory,GATKdir, vcfFiles, Reference, prefix+".final.vcf.temp",this.PRIORITIZE, 8 );
			printSBATCHscript(MasterShellScriptFile, sbatch, timestamp,0, outDir, commandLine,job,prefix+".final.vcf","GATK_mergeVCF",time);

		}catch(Exception E){E.printStackTrace();}
	}
	
	
	
	
	public void GATKPhase1Sample( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String outDir, ArrayList <String> fileNames,Hashtable<String,String> T){

		/*				
		for each sample (asumed to be in same dir)

	    lanes.bam <- merged lane.bams for sample
	    dedup.bam <- MarkDuplicates(lanes.bam) 
	    realigned.bam <- realign(dedup.bam) [with known sites included if available]
	    recal.bam <- recal(realigned.bam)
	    sample.bam <- recal.bam			
		 */

		if(!IOTools.isDir(outDir+"/reports"))
			IOTools.mkDir(outDir+"/reports");
		if(!IOTools.isDir(outDir+"/scripts"))
			IOTools.mkDir(outDir+"/scripts");
		try{


			boolean phase1 = true;
			if(T.containsKey("-ReAlign") || T.containsKey("-ReassignOneMappingQuality") || T.containsKey("-BQSR") || T.containsKey("-BQSRprint") || T.containsKey("-ReduceReads")){
				phase1 = false;
				for (int i = 0; i < fileNames.size();i++){ 
					String bamFile = fileNames.get(i);
					
					String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_"+i+"_phase1_GATK.sbatch";
					generalSbatchScript.println("sbatch "+ sbatchFileName);
					ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
					sbatch.printSBATCHinfo(EW,outDir,timestamp,0,"GATK_phase1", time);
					EW.println();
					EW.println("cd "+ outDir);
					
					// takes care of realining and recalibration
					boolean ReassignOneMappingQuality = false;
					if(!IOTools.fileExists(outDir+"/"+bamFile+".bai"))
						SamtoolsSBATCH.indexBamFile(bamFile, EW, true);

					if(T.containsKey("-ReassignOneMappingQuality")){
						bamFile = ReassignOneMappingQuality(EW, bamFile, memory, this.GATKdir,this.Reference, suffix);
					}
					if(T.containsKey("-ReAlign"))
						bamFile = ReAlign(EW, bamFile, memory, this.GATKdir,this.Reference, suffix,knownSNPVCF,knownIndelVCF,targetIntervals,ReassignOneMappingQuality);
					if(T.containsKey("-BQSR"))
						bamFile = BQSR(EW, bamFile, memory, this.GATKdir,this.Reference, suffix,knownSNPVCF,knownIndelVCF);
					if(T.containsKey("-BQSRprint"))
						bamFile = BQSRprint(EW, bamFile, memory, this.picardDir, this.GATKdir,this.Reference, suffix,knownSNPVCF,knownIndelVCF);
					if(T.containsKey("-ReduceReads"))
						ReduceReads(EW, bamFile, memory,  this.GATKdir,this.Reference);
					EW.println();
					EW.println();
					EW.flush();
					EW.close();
				}
			}
			if(phase1){ // run the entire pipeline
				if(T.containsKey("-preGATK")){
					String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_phase1_GATK.sbatch";
					generalSbatchScript.println("sbatch "+ sbatchFileName);
					ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
					sbatch.printSBATCHinfo(EW,outDir,timestamp,0,"GATK_phase1", time);
					EW.println();
					EW.println("cd "+ outDir);

					// takes care of merging, correct naming of bam files and marking dulicates
					String bamFile = Picard.preGATK(EW, sbatch, timestamp, outDir, fileNames, T, this.picardDir,this.memory ,this.suffix);
					
					// Changes suffix from sam to bam
					if(suffix.lastIndexOf("sam")>-1){
						this.suffix=this.suffix.substring(0,suffix.lastIndexOf("sam"))+"bam";
					}
				
					
					// takes care of realining and recalibration
					phase1(EW, bamFile, memory, this.picardDir, this.GATKdir,this.Reference, suffix,knownSNPVCF,knownIndelVCF,targetIntervals,true);
					EW.println();
					EW.println();
					EW.println("wait");

					EW.flush();
					EW.close();
				}else{
					// expects that all the files are in one bam file
					for (int i = 0; i < fileNames.size();i++){ 
						String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_"+i+"_phase1_GATK.sbatch";
						generalSbatchScript.println("sbatch "+ sbatchFileName);
						ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
						sbatch.printSBATCHinfo(EW,outDir,timestamp,0,"GATK_phase1", time);
						EW.println();
						EW.println("cd "+ outDir);

						String bamFile = fileNames.get(i);
						// takes care of realining and recalibration
						phase1(EW, bamFile, memory, this.picardDir, this.GATKdir,this.Reference, suffix,knownSNPVCF,knownIndelVCF,targetIntervals,true);
						EW.println();
						EW.println();
						EW.println("wait");

						EW.flush();
						EW.close();
					}
				}
			}


		}catch(Exception E){E.printStackTrace();}
	}

	public void GATKGenotypeSample( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String outDir, ArrayList <String> fileNames,Hashtable<String,String> T){
		System.out.println("running GATKGenotype");

		/*				
		for each sample (asumed to be in same dir)

	    lanes.bam <- merged lane.bams for sample
	    dedup.bam <- MarkDuplicates(lanes.bam) 
	    realigned.bam <- realign(dedup.bam) [with known sites included if available]
	    recal.bam <- recal(realigned.bam)
	    sample.bam <- recal.bam			
		 */

		if(!IOTools.isDir(outDir+"/reports"))
			IOTools.mkDir(outDir+"/reports");
		if(!IOTools.isDir(outDir+"/scripts"))
			IOTools.mkDir(outDir+"/scripts");
		try{


			for (int i = 0; i < fileNames.size();i++){ 
				String bamFile = fileNames.get(i);
					
				String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_"+i+"_Genotype_GATK.sbatch";
				generalSbatchScript.println("sbatch "+ sbatchFileName);
				ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
				sbatch.printSBATCHinfo(EW,outDir,timestamp,i,"GATK_Genotype", time);
				EW.println();
				EW.println("cd "+ outDir);
					
				String commandLine = Genotype(23, this.GATKdir, bamFile, this.Reference, bamFile+".Genotype.vcf", this.knownSNPVCF);
				EW.println(commandLine);
				System.out.println(commandLine);

				EW.println();
				EW.println();
				EW.flush();
				EW.close();
			}


		}catch(Exception E){E.printStackTrace();}
	}


	public void GATKPhaseSNPsSample( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String outDir, ArrayList <String> fileNames,Hashtable<String,String> T){
		System.out.println("running GATKPhaseSNPsSample");

		/*				
		for each sample (asumed to be in same dir)

	    lanes.bam <- merged lane.bams for sample
	    dedup.bam <- MarkDuplicates(lanes.bam) 
	    realigned.bam <- realign(dedup.bam) [with known sites included if available]
	    recal.bam <- recal(realigned.bam)
	    sample.bam <- recal.bam			
		 */

		if(!IOTools.isDir(outDir+"/reports"))
			IOTools.mkDir(outDir+"/reports");
		if(!IOTools.isDir(outDir+"/scripts"))
			IOTools.mkDir(outDir+"/scripts");
		try{


			for (int i = 0; i < fileNames.size();i++){ 
				String bamFile = fileNames.get(i);
					
				String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_"+i+"_phaseSNPs_GATK.sbatch";
				generalSbatchScript.println("sbatch "+ sbatchFileName);
				ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
				sbatch.printSBATCHinfo(EW,outDir,timestamp,i,"GATK_phase1", time);
				EW.println();
				EW.println("cd "+ outDir);
					
				String commandLine = ReadBackedPhasing(23, this.GATKdir, bamFile, this.Reference, bamFile+".phased.vcf", this.knownSNPVCF,10.0);
				EW.println(commandLine);
				System.out.println(commandLine);

				EW.println();
				EW.println();
				EW.flush();
				EW.close();
			}


		}catch(Exception E){E.printStackTrace();}
	}





	public static String phase1( ExtendedWriter EW, String bamFile, int memory, String picardDir, String GATKDir,String Reference, String suffix, String knownSNPVCF, String knownIndelsVCF,String targetIntervals, boolean ReassignOneMappingQuality){
		//		if(IOTools.fileExists(bamFile)){

		if(!IOTools.fileExists(bamFile+".bai")){
			EW.println("# index bam files  ");
			SamtoolsSBATCH.indexBamFile(bamFile, EW,true);
			EW.println();
			EW.println();
		}

		bamFile = ReAlign(EW, bamFile, memory, GATKDir, Reference, suffix, knownSNPVCF, knownIndelsVCF,targetIntervals,ReassignOneMappingQuality);
		bamFile = BQSR(EW, bamFile, memory, GATKDir, Reference, suffix, knownSNPVCF, knownIndelsVCF);
		bamFile = ReduceReads(EW, bamFile, memory, GATKDir, Reference);

		return bamFile;
		//		}
		//		return null;
	}

	public static String ReassignOneMappingQuality( ExtendedWriter EW, String bamFile, int memory, String GATKDir,String Reference, String suffix){
		
		String baseName = bamFile.substring(0,bamFile.indexOf(".bam"));			
		String commandLine = ReassignOneMappingQuality(memory,GATKDir, bamFile, Reference, baseName);
			
		EW.println("# ReassignOneMappingQuality  ");
		EW.println(commandLine);
		EW.println();
		EW.println("echo ReassignOneMappingQuality finished");
		EW.println();
	
		return baseName+".ReassignOneMappingQuality.bam";
	}



	public static String ReAlign( ExtendedWriter EW, String bamFile, int memory, String GATKDir,String Reference, String suffix, String knownSNPVCF, String knownIndelsVCF, String targetIntervals , boolean ReassignOneMappingQuality){
		
		String baseName = bamFile.substring(0,bamFile.lastIndexOf(suffix));
		if(targetIntervals == null){
			
			String commandLine = createIntervals(memory,GATKDir, bamFile, Reference, baseName,knownIndelsVCF);
			
			EW.println("# createIntervals  ");
			EW.println(commandLine);
			EW.println();
			EW.println("echo Create Intervals finished");
			EW.println();
			
			targetIntervals = baseName+".intervals";
			
			
		}

		
		
		
		String commandLine = LocalRealignment(memory,GATKDir, bamFile,Reference, baseName, knownIndelsVCF,targetIntervals);
		EW.println("# LocalRealignment  ");
		if(ReassignOneMappingQuality){
			commandLine = ReassignOneMappingQuality(commandLine);
			EW.println("# change mapping score for RNA seq data  ");
		}
		EW.println(commandLine);
		EW.println();

		EW.println("echo Local Realignment finished");
		EW.println();

		return baseName+".real.bam";
	}

	public static String BQSR( ExtendedWriter EW, String bamFile, int memory, String GATKDir,String Reference, String suffix, String knownSNPVCF, String knownIndelsVCF){
		String baseName = bamFile.substring(0,bamFile.indexOf(".bam"));
		EW.println("# BaseQualityRecalibration  ");

		String commandLine = BaseQualityReport(memory,GATKDir, baseName+".bam",Reference, baseName, knownSNPVCF,knownIndelsVCF);
		EW.println(commandLine);
		EW.println();
		EW.println("echo Base Quality Report finished");
		EW.println();

		String commandLine2 = BaseQualityRecalibration(memory,GATKDir, baseName+".bam",Reference, baseName, baseName+".grp");
		EW.println(commandLine2);
		EW.println();
		EW.println("echo Base Quality Recalibration finished");
		EW.println();

		return baseName+".BQSR.bam";
	}

	public static String BQSRprint( ExtendedWriter EW, String bamFile, int memory, String picardDir, String GATKDir,String Reference, String suffix, String knownSNPVCF, String knownIndelsVCF){
		String baseName = bamFile.substring(0,bamFile.indexOf(".bam"));
		EW.println("# BaseQualityRecalibration  ");

		String commandLine2 = BaseQualityRecalibration(memory,GATKDir, baseName+".bam",Reference, baseName, baseName+".grp");
		EW.println(commandLine2);
		EW.println();
		EW.println("echo Base Quality Recalibration finished");
		EW.println();

		return baseName+".BQSR.bam";
	}


	public static String ReduceReads( ExtendedWriter EW, String bamFile, int memory,  String GATKDir,String Reference){
		String baseName = bamFile.substring(0,bamFile.indexOf(".bam"));
		EW.println("# ReduceReads shell script  ");

		String commandLine2 = ReduceReads(memory,GATKDir, baseName,Reference);
		EW.println(commandLine2);
		EW.println();
		EW.println("echo ReduceReads finished");
		EW.println();

		return baseName+".reduced.bam";
	}

	public static String Genotype(int memory, String GATKDir, String bamFile, String reference, String baseName, String knownSNPVCF){

		/*		
		java -Xmx23g -jar /bubo/sw/apps/bioinfo/GATK/1.2.12/GenomeAnalysisTK.jar 
		-R /bubo/home/h12/pallol/glob/projects/b2010010/alignment/bwa/variation/assembly/fAlb13.masked.fa 
		-T UnifiedGenotyper 
		-nt 12 
		-I /bubo/home/h12/pallol/glob/projects/b2010010/alignment/bwa/variation/low_coverage/mult_genome_20101112/raw/vs_fAlb13rm/s_5.fastq.bwa.q33.sorted.real.md.recal.bam 
		-o COLL_34.UNIONSITES.ALLLIBS.SNP.vs_fAlb13rm.vcf 
		-alleles:masterAlleles /proj/b2010010/private/variation/vs_fAlb13rm/COLL-PIED.UNION.vcf 
		-gt_mode GENOTYPE_GIVEN_ALLELES 
		-out_mode EMIT_ALL_SITES 
		-BTI masterAlleles 
		-stand_call_conf 0.0 
		-G none
		
		
		
	  java
      -jar /GenomeAnalysisTK.jar
      -T  GenotypeAndValidate
      -R human_g1k_v37.fasta
      -I myNewTechReads.bam
      -alleles handAnnotatedVCF.vcf
      -L handAnnotatedVCF.vcf
*/
		String commandLine = "java -Xmx"+memory+"g -jar "+GATKDir+"/GenomeAnalysisTK.jar "+
				"-T UnifiedGenotyper "+
				"-nt 7 "+
				"-I "+bamFile+ " "+
				"-R "+reference+" "+
				"-o "+bamFile+"Genotype.vcf "+
				"-alleles  "+knownSNPVCF+" "+
				"-gt_mode GENOTYPE_GIVEN_ALLELES --output_mode EMIT_ALL_SITES "+				
				"-L  "+knownSNPVCF;

		return commandLine;

		
		
	}

	public static String BaseQualityReport(int memory, String GATKDir, String bamFile, String reference, String baseName, String knownSNPVCF, String knownIndelsVCF){

		/*
		java -Xmx4g -jar GenomeAnalysisTK.jar \
		   -T BaseRecalibrator \
		   -I my_reads.bam \
		   -R resources/Homo_sapiens_assembly18.fasta \
		   -knownSites bundle/hg18/dbsnp_132.hg18.vcf \
		   -knownSites another/optional/setOfSitesToMask.vcf \
		   -o recal_data.grp

		 */

		String commandLine = "java -Xmx"+memory+"g -jar "+GATKDir+"/GenomeAnalysisTK.jar "+
				"-T BaseRecalibrator "+
				"-I "+bamFile+ " "+
				"-R "+reference+" "+
				"-o "+baseName+".grp";
		if(knownSNPVCF!=null)
			commandLine+=" -knownSites "+knownSNPVCF;
		if(knownIndelsVCF!=null)
			commandLine+=" -knownSites "+knownIndelsVCF;

		return commandLine;

	}


	public static String ReduceReads(int memory, String GATKDir, String baseName, String reference){

		/*
 			java -Xmx4g -jar GenomeAnalysisTK.jar \
   			-R ref.fasta \
   -T ReduceReads \
   -I myData.bam \
   -o myData.reduced.bam

		 */

		String commandLine = "java -Xmx"+memory+"g -jar "+GATKDir+"/GenomeAnalysisTK.jar "+
				"-T ReduceReads "+
				"-I "+baseName+ ".bam "+
				"-R "+reference+" "+
				"-o "+baseName+".reduced.bam";

		return commandLine;

	}


	public static String BaseQualityRecalibration(int memory, String GATKDir, String bamFile, String reference, String baseName, String BSQRfile){

		/*
			java -jar GenomeAnalysisTK.jar \
			   -T PrintReads \
			   -R reference.fasta \
			   -I input.bam \
			   -BQSR recalibration_report.grp \
			   -o output.bam
			- See more at: http://gatkforums.broadinstitute.org/discussion/44/base-quality-score-recalibration-bqsr#sthash.8xVFUaHp.dpuf		 
		 */

		String commandLine = "java -Xmx"+memory+"g -jar "+GATKDir+"/GenomeAnalysisTK.jar "+
				"-T PrintReads "+
				"-I "+bamFile+ " "+
				"-R "+reference+" "+
				"-BQSR "+BSQRfile+" "+
				"-o "+baseName+".BQSR.bam";

		return commandLine;

	}

	public static String ReassignOneMappingQuality(int memory, String GATKDir, String bamFile, String reference, String baseName){

		/*
			java -jar GenomeAnalysisTK.jar \
				  -R ref.fasta \
				   -T PrintReads \
				   -o output.bam \
				   -I input1.bam \
				   -I input2.bam \
				   -rf ReassignOneMappingQuality
      			-RMQF 255
      			-RMQT 60
			- See more at: http://gatkforums.broadinstitute.org/discussion/44/base-quality-score-recalibration-bqsr#sthash.8xVFUaHp.dpuf		 
		 */

		String commandLine = "java -Xmx"+memory+"g -jar "+GATKDir+"/GenomeAnalysisTK.jar "+
				"-I "+bamFile+ " "+
				"-R "+reference+" "+
				"-T PrintReads "+
				"-o "+baseName+".ReassignOneMappingQuality.bam"+
				" -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 ";
		

		return commandLine;

	}



	public static String createIntervals(int memory, String GATKDir, String bamFile, String reference, String baseName, String knownIndelsVCF){
		/*		java -Xmx2g -jar GenomeAnalysisTK.jar \
		   -I input.bam \
		   -R ref.fasta \
		   -T RealignerTargetCreator \
		   -o forIndelRealigner.intervals \
		   [--known /path/to/indels.vcf]

		 */		
		String commandLine = "java -Xmx"+memory+"g -jar "+GATKDir+"/GenomeAnalysisTK.jar "+
				"-I "+bamFile+ " "+
				"-R "+reference+" "+
				"-T RealignerTargetCreator "+
				"-o "+baseName+".intervals";
		if(knownIndelsVCF!=null)
			commandLine+=" -known "+knownIndelsVCF;

		return commandLine;


	}
	
	public static String ReassignOneMappingQuality(String commandLine){
	/*
	
	java
    -jar GenomeAnalysisTK.jar
    -rf ReassignOneMappingQuality
    -RMQF 255
    -RMQT 60
    */
		
		commandLine += " -rf ReassignOneMappingQuality";
		return commandLine;
		
		
		
	}
	public static String LocalRealignment(int memory, String GATKDir, String bamFile, String reference, String baseName, String knownSNPVCF, String targetIntervals){


		/*			java -Xmx4g -jar GenomeAnalysisTK.jar \
			   -I input.bam \
			   -R ref.fasta \
			   -T IndelRealigner \
			   -targetIntervals intervalListFromRTC.intervals \
			   -o realignedBam.bam \
			   [-known /path/to/indels.vcf] \
			   [-compress 0]    (this argument recommended to speed up the process *if* this is only a temporary file; otherwise, use the default value)

		 */		


		String commandLine = "java -Xmx"+memory+"g -jar "+GATKDir+"/GenomeAnalysisTK.jar "+
				"-I "+bamFile+ " "+
				"-R "+reference+" "+
				"-T IndelRealigner "+
				"-targetIntervals "+targetIntervals+" "+
				"-o "+baseName+".real.bam";
		if(knownSNPVCF!=null)
			commandLine+=" -known "+knownSNPVCF;

		return commandLine;

	}

	public static String HaplotypeCaller(int memory, String GATKDir, ArrayList <String> bamFiles, String reference, String outFileVCF, String knownSNPVCF, ArrayList <String> contigs){

		/*  java
	     -jar GenomeAnalysisTK.jar
	     -T HaplotypeCaller
	     -R reference/human_g1k_v37.fasta
	     -I sample1.bam [-I sample2.bam ...] \
	     --dbsnp dbSNP.vcf \
	     -stand_call_conf [50.0] \
	     -stand_emit_conf 10.0 \
	     [-L targets.interval_list]
	     -o output.raw.snps.indels.vcf
		 */		




		String commandLine = "java -Xmx"+memory+"g -jar "+GATKDir+"/GenomeAnalysisTK.jar "+
				"-R "+reference+" "+
				"-T HaplotypeCaller "+
				"-o "+outFileVCF;
		for(int i = 0; i < bamFiles.size();i++){
			commandLine+=" -I "+bamFiles.get(i);
		}
		if(contigs != null){
			for(int i = 0; i < contigs.size();i++){
				commandLine+=" -L "+contigs.get(i);
			}
		}



		if(knownSNPVCF!=null)
			commandLine+=" --dbsnp "+knownSNPVCF;

		return commandLine;

	}


	public static String UnifiedGenotypeCaller(int memory, String GATKDir, ArrayList <String> bamFiles, String reference, String outFileVCF, String knownSNPVCF, ArrayList <String> contigs){

		/*  java
	     -jar GenomeAnalysisTK.jar
	     -T HaplotypeCaller
	     -R reference/human_g1k_v37.fasta
	     -I sample1.bam [-I sample2.bam ...] \
	     --dbsnp dbSNP.vcf \
	     -stand_call_conf [50.0] \
	     -stand_emit_conf 10.0 \
	     [-L targets.interval_list]
	     -o output.raw.snps.indels.vcf
		 */		




		String commandLine = "java -Xmx"+memory+"g -jar "+GATKDir+"/GenomeAnalysisTK.jar "+
				"-R "+reference+" "+
				"-T UnifiedGenotyper "+
				"-dcov 200 "+
				"-o "+outFileVCF;
		for(int i = 0; i < bamFiles.size();i++){
			commandLine+=" -I "+bamFiles.get(i);
		}
		if(contigs != null){
			for(int i = 0; i < contigs.size();i++){
				commandLine+=" -L "+contigs.get(i);
			}
		}



		if(knownSNPVCF!=null)
			commandLine+=" --dbsnp "+knownSNPVCF;

		return commandLine;

	}


	public static String ReadBackedPhasing(int memory, String GATKDir, String bamFile, String reference, String outFileVCF, String knownSNPVCF,double phaseQualityThresh){

		/*     java
      -jar GenomeAnalysisTK.jar
      -T ReadBackedPhasing
      -R reference.fasta
      -I reads.bam
      --variant SNPs.vcf
      -L SNPs.vcf
      -o phased_SNPs.vcf
      --phaseQualityThresh 20.0
 
		 */		




		String commandLine = "java -Xmx"+memory+"g -jar "+GATKDir+"/GenomeAnalysisTK.jar "+
				"-R "+reference+" "+
				"-T ReadBackedPhasing "+
				"--variant  "+knownSNPVCF+" "+
				"-L  "+knownSNPVCF+" "+
				"-o "+outFileVCF+
				" --phaseQualityThresh "+  phaseQualityThresh;
				commandLine+=" -I "+bamFile;

		return commandLine;

	}

	
	
	
	
	public static String mergeVCFfiles(int memory, String GATKDir, ArrayList <String> vcfFiles, String reference, String outFileVCF,boolean PRIORITIZE, int cores ){ 
		/*	 
 java -Xmx2g -jar GenomeAnalysisTK.jar \
   -R ref.fasta \
   -T CombineVariants \
   --variant input1.vcf \
   --variant input2.vcf \
   -o output.vcf \
   -genotypeMergeOptions UNIQUIFY

 java -Xmx2g -jar GenomeAnalysisTK.jar \
   -R ref.fasta \
   -T CombineVariants \
   --variant:foo input1.vcf \
   --variant:bar input2.vcf \
   -o output.vcf \
   -genotypeMergeOptions PRIORITIZE
   -priority foo,bar
		 */

		String commandLine = "java -Xmx"+memory+"g -jar "+GATKDir+"/GenomeAnalysisTK.jar "+
				"-R "+reference+" "+
				"-T CombineVariants "+
				"-o "+outFileVCF+
				" -nt "+cores;

		for(int i = 0; i < vcfFiles.size();i++){
			commandLine+=" --variant:f"+i+" "+vcfFiles.get(i);
		}
		if(PRIORITIZE){
			commandLine+=" -priority f0";
			for(int i = 1; i < vcfFiles.size();i++){
				commandLine+=",f"+i;
			}
		}
		else{
			commandLine+=" UNIQUIFY";
		}
		return commandLine;
	}

}

