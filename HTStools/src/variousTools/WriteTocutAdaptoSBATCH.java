package variousTools;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Hashtable;

import general.ExtendedReader;
import general.ExtendedWriter;
import general.Functions;
import general.IOTools;

public class WriteTocutAdaptoSBATCH {



	String projectNumber;
	String time;
	ArrayList <String> threeAdapters;
	ArrayList <String> otherAdapters;

	int cutoff;
	boolean cutadapt;
	boolean QC;
	String suffix;
	String[] sep;
	int length;
	int overlap;
	String timeStamp;
	boolean hiseq;
	boolean seqPrep;


	public static void main(String []args){



		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}
		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);

		WriteTocutAdaptoSBATCH sbatchscript = new WriteTocutAdaptoSBATCH();
		sbatchscript.run(T);


	}

	public WriteTocutAdaptoSBATCH(){
		projectNumber = time = null;
		cutadapt = true;
		QC = true;
		seqPrep = false;
		hiseq=false;
	}

	public void run(Hashtable<String,String> T){

		boolean allPresent = true;

		SBATCHinfo sbatch = new SBATCHinfo();
		if(!sbatch.addSBATCHinfo(T)) allPresent = false;
		if(!T.containsKey("-i")){
			System.out.println("must contain inDirectory -i");
			allPresent = false;
		}

		if(!T.containsKey("-seqPrep")){
			System.out.println();
			System.out.println("Does not contina seqPrep step for merging;");
			System.out.println("Add flag -seqPrep for merging step");
			System.out.println("Add flag -6 if the phred score are ascci from 31 (hiseq-instruments)");
		}else{
			seqPrep=true;
			if(T.containsKey("-seqPrep"))hiseq=true;
		}



		if(!T.containsKey("-t"))
			System.out.println("must contain likely time  (-t 1:00:00). Time now set to default 1:00:00");


		if(T.containsKey("-a") || T.containsKey("-b")){
			if(T.containsKey("-a")){
				String threePrimeAdaptersFile= Functions.getValue(T, "-a", ".");
				threeAdapters = getAdapters(threePrimeAdaptersFile);
			}
			if( T.containsKey("-b")){
				String allAdaptersFile= Functions.getValue(T, "-b", ".");
				otherAdapters = getAdapters(allAdaptersFile);

			}
		}
		else{
			System.out.println("must contain a file that contains -3Adapter or -Adapter now only QC of sequences will be carried out.");
			cutadapt = false;
		}

		String projectDir= Functions.getValue(T, "-pDir", IOTools.getCurrentPath());
		String inDir= Functions.getValue(T, "-i");
		String cutAdaptDir= Functions.getValue(T, "-cutadapt", inDir+"_cutadapt");
		String QCDir= Functions.getValue(T, "-QC", inDir+"_QC");
		String seqPrepDir= Functions.getValue(T, "-seqPrep", inDir+"_seqPrep");
		time = Functions.getValue(T, "-t", "1:00:00");
		timeStamp = Functions.getValue(T, "-TS", Functions.getDateTime());
		cutoff = Integer.parseInt(Functions.getValue(T, "-q", "-1"));
		overlap = Integer.parseInt(Functions.getValue(T, "-O", "-1"));
		length = Integer.parseInt(Functions.getValue(T, "-l", "30"));
		suffix = Functions.getValue(T,"-suffix","fastq");
		this.sep = new String[2];
		sep[0] = "1."+suffix;
		sep[1] = "2."+suffix;

		if(allPresent)
			cutAdapt(sbatch, projectDir,inDir, cutAdaptDir,QCDir,seqPrepDir);
		else
			System.out.println("\n\nAborting run because of missing arguments for cutadapt.");
	}


	public void cutAdapt(SBATCHinfo sbatch,String projectDir,String inDir, String cutAdaptDir, String QCDir, String seqPrepDir){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"cutadapt.sh"));
			
			
			cutAdaptSample(sbatch, EW,projectDir+"/"+inDir,projectDir+"/"+cutAdaptDir,projectDir+"/"+QCDir,projectDir+"/"+seqPrepDir);
			EW.println();
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	} 



	public  void cutAdaptSample(SBATCHinfo sbatch,ExtendedWriter generalSbatchScript, String inDir, String cutadaptDir,String QCDir,String SeqPrepDir){
		String finalDir = QCDir;
		try{
			ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir, suffix);
			if(fileNames != null && !fileNames.isEmpty()){
				if(cutadapt){
					if(!IOTools.isDir(cutadaptDir))
						IOTools.mkDirs(cutadaptDir);
					if(!IOTools.isDir(QCDir))
						IOTools.mkDirs(QCDir);
				}
				if(seqPrep){
					if(!IOTools.isDir(SeqPrepDir))
						IOTools.mkDirs(SeqPrepDir);
					finalDir = SeqPrepDir;
				}

				if(!IOTools.isDir(finalDir+"/reports"))
					IOTools.mkDir(finalDir+"/reports");
				if(!IOTools.isDir(finalDir+"/scripts"))
					IOTools.mkDir(finalDir+"/scripts");



				ArrayList <String[]> pairs = IOTools.findPairs(fileNames,sep);
				for(int i = 0; i < pairs.size(); i++){
					String sbatchFile = finalDir+"/scripts/"+timeStamp+"_"+i+"_cutadapt.sbatch";
					generalSbatchScript.println("sbatch "+ sbatchFile);
					ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFile));
					String[] split = SeqPrepDir.split("/");
					sbatch.printSBATCHinfoCore(EW,cutadaptDir,timeStamp,0,i+"_cutadapt_"+split[split.length-1], time);

					if(cutadapt){
						//cutAdapt step

						for(int k = 0; k < 2; k++){
							addCutAdaptStep(EW,inDir,cutadaptDir,pairs.get(i)[k]);
						}

						//QC step
						addQCstep(EW,cutadaptDir,QCDir,pairs.get(i)[0],pairs.get(i)[1]);
					}

					//SeqPrep step merging sequences
					if(seqPrep){
						SeqPrep SP = new SeqPrep();
						if(hiseq)SP.hiseq=true;
						SP.filter_fastqSample(EW, QCDir, SeqPrepDir, pairs.get(i)[0], pairs.get(i)[1], pairs.get(i)[2]);
					}
					FastQCSBATCH FQC = new FastQCSBATCH();
					//FastQCstep
					FQC.FastQCSample(EW, inDir, pairs.get(i)[0]);
					FQC.FastQCSample(EW, inDir, pairs.get(i)[1]);
					FQC.FastQCSample(EW, cutadaptDir, pairs.get(i)[0]);
					FQC.FastQCSample(EW, cutadaptDir, pairs.get(i)[1]);
					FQC.FastQCSample(EW, QCDir, pairs.get(i)[0]);
					FQC.FastQCSample(EW, QCDir, pairs.get(i)[1]);
					FQC.FastQCSample(EW, SeqPrepDir, pairs.get(i)[0]);
					FQC.FastQCSample(EW, SeqPrepDir, pairs.get(i)[1]);
					FQC.FastQCSample(EW, SeqPrepDir, pairs.get(i)[2]+"merged.fastq.gz");

					EW.println();
					EW.flush();
					EW.close();
				}
			}
			else{
				System.out.println("No .fastq files in folder :"+inDir);

			}
			ArrayList <String> samples = IOTools.getDirectories(inDir);
			if(samples != null){
				for(int i = 0; i < samples.size(); i++){
					cutAdaptSample(sbatch, generalSbatchScript,inDir+"/"+samples.get(i),cutadaptDir+"/"+samples.get(i),QCDir+"/"+samples.get(i),SeqPrepDir+"/"+samples.get(i));
				}
			}

		}catch(Exception E){E.printStackTrace();}
	}
	private void addCutAdaptStep(ExtendedWriter EW,String inDir, String outDir,String file){
		EW.println();
		EW.println();
		EW.println("#############################################################################################################");
		EW.println("Running cutadapt START");
		EW.println();
		EW.println();
		EW.println();
		EW.print("cutadapt");
		if(this.overlap > -1)EW.print(" --overlap="+this.overlap+" ");
		if(cutoff > -1)EW.print(" -q "+cutoff+" ");
		if(threeAdapters!= null){
			for (int j = 0; j < threeAdapters.size(); j++){
				EW.print(" -a "+threeAdapters.get(j));
			}
		}
		if(otherAdapters!=null){
			for (int j = 0; j < otherAdapters.size(); j++){
				EW.print(" -b "+otherAdapters.get(j));
			}
		}
		EW.print(" -o "+ outDir+"/"+file);
		EW.println(" "+inDir+"/"+file);
		EW.println();
		EW.println("running cutadapt DONE");
		EW.println("#############################################################################################################");
		EW.println();
		EW.println();

	}



	private void addQCstep(ExtendedWriter EW, String inDir, String outDir, String forward, String reverse){
		EW.println();
		EW.println();
		EW.println("#############################################################################################################");
		EW.println("Removing short reads START");
		EW.println();
		EW.println();
		EW.println("java -jar /bubo/home/h17/johanr/bin/HTStools.jar -p sequenceHandling QC -dir "+inDir+ " -o "+ outDir+" -l "+length +" -f1 "+forward+" -f2 "+reverse);
		EW.println();
		EW.println("Removing short reads DONE");
		EW.println("#############################################################################################################");
		EW.println();
		EW.println();

	}

	public  static ArrayList<String> getAdapters(String file){

		ArrayList<String> SequenceFiles = new ArrayList<String>();
		if(file != null){
			try{
				ExtendedReader ER = new ExtendedReader(new FileReader(file));
				while(ER.more()){
					SequenceFiles.add(ER.readLine());
				}
			}catch (Exception E){E.printStackTrace();}
		}
		return SequenceFiles;
	}
}
