package Sequence;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import general.ExtendedWriter;


import general.Functions;
import general.IOTools;
import general.RNAfunctions;

public class SequenceHandling {

	public static void main(String []args){
		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}
		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);
		SequenceHandling.run(T);
	}

	public static void run(Hashtable<String,String> T){

		String program = Functions.getValue(T, "-p", "h").toUpperCase();
		if(program.indexOf("UNIQUE") != -1){
			SequenceHandling.getOtherSequences(T);
		}
		else if(program.indexOf("MERGE") != -1){
			SequenceHandling.MergeSequences(T);
		}
		else if(program.indexOf("SAMTOOLS") != -1){
			SamSequences.run(T);
		}
		else if(program.indexOf("GENERATE") != -1){
			Generate.run(T);
		}
		else if(program.indexOf("INTERACTION") != -1){
		}
		else if(program.indexOf("LENGTH") != -1){
			SequenceHandling.extractAbove(T);
		}
		else if(program.indexOf("EXACT") != -1){
			SequenceHandling.extractSize(T);
		}
		else if(program.indexOf("ORFS") != -1){
			SequenceHandling.getORFs(T);
		}
		else if(program.indexOf("5END") != -1){
			SequenceHandling.print5Ends(T);
		}
		else if(program.indexOf("3END") != -1){
			SequenceHandling.print3Ends(T);
		}
		else if(program.indexOf("extractEnds".toUpperCase()) != -1){
			SequenceHandling.extractEnds(T);
		}
		else if(program.indexOf("removePrimer".toUpperCase()) != -1){

			SequenceHandling.removePrimers(T);
		}
		else if(program.indexOf("integrity".toUpperCase()) != -1){
			SequenceHandling.integrity(T);
		}
		else if(program.indexOf("check".toUpperCase()) != -1){
			SequenceHandling.check(T);
		}
		else if(program.indexOf("size".toUpperCase()) != -1){
			SequenceHandling.size(T);
		}
		else if(program.indexOf("fixSequenceNames".toUpperCase()) != -1){
			SequenceHandling.fixSequenceNames(T);
		}

		else if(program.indexOf("fastq2fasta".toUpperCase()) != -1){
			FastQSequences.fastq2fasta(T);
		}

		else if(program.indexOf("RNA2DNA".toUpperCase()) != -1){
			String infile = Functions.getValue(T, "-i");
			FastaSequences.RNA2DNA(infile);
		}


		else if(program.indexOf("mutationAnalysis".toUpperCase()) != -1){
			SequenceHandling.mutationalAnalysis(T);
		}

		else if(program.indexOf("extract".toUpperCase()) != -1){
			SequenceHandling.extract(T);
		}
		else if(T.containsKey("filter2")){
			CfastaSequences.run(T);
		}

		else if(T.containsKey("GCcount")){
			SequenceHandling.GCcount(T);
		}
		else if(program.indexOf("miRNAdistr".toUpperCase()) != -1){
			SequenceHandling.readsDistribution(T);
		}

		else if(program.indexOf("oneLine".toUpperCase()) != -1){
			SequenceHandling.oneLine(T);
		}

		else if(program.indexOf("split".toUpperCase()) != -1){
			SequenceHandling.split(T);
		}

		
		else if(program.indexOf("stratify".toUpperCase()) != -1){
			SequenceHandling.stratify(T);
		}
	
		else if(program.indexOf("one2two".toUpperCase()) != -1){
			SequenceHandling.one2two(T);
		}
	
		
		else if(program.indexOf("QC".toUpperCase()) != -1){
			SequenceHandling.QC(T);
		}

		else if(program.indexOf("mutual".toUpperCase()) != -1){
			SequenceHandling.mutual(T);
		}

		if(T.containsKey("-extra"))
			FastaSequences.run(T);
		if(T.containsKey("-subset"))
			FastaSequences.run(T);
	}

	public static FastaSequences getOtherSequences(Hashtable<String,String> T){
		String dir = Functions.getValue(T, "-d", "");
		String fastaFile1 = Functions.getValue(T, "-file1", "file1.fa");
		String fastaFile2 = Functions.getValue(T, "-file2", "file2.fa");

		FastaSequences FS1 = new FastaSequences(dir, fastaFile1);
		FastaSequences FS2 = new FastaSequences(dir, fastaFile2);

		FastaSequences UniqueNames1 = FastaSequences.getUniqueNames(FS1,FS2);
		FastaSequences UniqueNames2 = FastaSequences.getUniqueNames(FS2,FS1);

		FastaSequences Unique = FastaSequences.merge(UniqueNames1,UniqueNames2);
		String fastaFile3 = Functions.getValue(T, "-outFile", "outFile.fa");
		Unique.writeFastaFile(dir,fastaFile3);

		return Unique;
	}

	public static void removePrimers(Hashtable<String,String> T){
		boolean allRequired = true;
		String dir, outDir, SequenceFile, primerFile, suffix;
		dir = outDir = SequenceFile = primerFile  = suffix = null;
		if(T.containsKey("-id"))
			dir = Functions.getValue(T, "-id", "");
		else{
			allRequired = false;
			System.out.println("must have indirectory -id");

		}
		if(T.containsKey("-od"))
			outDir = Functions.getValue(T, "-od", "");
		else{
			allRequired = false;
			System.out.println("must have outdirectory -id");
		}
		if(T.contains("-i"))
			SequenceFile = Functions.getValue(T, "-i", "file1.fa");
		if(T.containsKey("-primerFile"))
			primerFile = Functions.getValue(T, "-primerFile", "primer.fa");
		else{
			allRequired = false;
			System.out.println("must have primerFile -primerFile");

		}
		suffix = Functions.getValue(T, "-suffix", "fastq");
		int lowCutoff = Integer.parseInt(Functions.getValue(T, "-lb", "17"));
		int highCutoff = Integer.parseInt(Functions.getValue(T, "-lb", "32"));

		if(allRequired){
			FastaSequences primers = new FastaSequences(primerFile);
			if(suffix.indexOf("fastq") != -1){}
			else if(suffix.indexOf("csfasta") != -1){
				CfastaSequences csPrimers = CfastaSequences.convertFasta2CSfasta(primers);
				if(T.containsKey("i")){
					try{
						ExtendedWriter EW3 = new ExtendedWriter(new FileWriter(outDir+"/info.txt"));
						CfastaSequences.removePrimers(EW3,dir,outDir, SequenceFile,csPrimers,lowCutoff,highCutoff);
					}catch(Exception E){E.printStackTrace();}

				}else{
					CfastaSequences.removePrimersDir(dir,outDir,csPrimers,suffix,lowCutoff,highCutoff);
				}
			}
			else if(suffix.indexOf("fa") != -1){}
			else{
				System.out.println("kind has to be either fastq, csfasta or fa");
			}
		}
	}


	public static void GCcount(Hashtable<String,String> T){
		boolean allRequired = true;
		String SequenceFile, suffix;
		SequenceFile =  suffix = null;
		SequenceFile = Functions.getValue(T, "-i", "file1.fa");
		suffix = Functions.getValue(T, "-suffix", "fasta");
		if(suffix.indexOf("fastq") != -1){}
		else if(suffix.indexOf("csfasta") != -1){}
		else if(suffix.indexOf("fasta") != -1){
			FastaSequences A = new FastaSequences();
			A.readFastaFile(SequenceFile);
			A.printGCcontent();
		}
		else{
			System.out.println("kind has to be either fastq, csfasta or fa");
		}
	}

	public static void readsDistribution(Hashtable<String,String> T){
		String SequenceFile, suffix, gmapperFile;
		SequenceFile =  suffix = null;
		SequenceFile = Functions.getValue(T, "-i", "file1.fa");
		suffix = Functions.getValue(T, "-suffix", "fasta");
		gmapperFile = Functions.getValue(T, "-r", "something.gmapper");
		int cutoff  = Integer.parseInt(Functions.getValue(T, "-x", "1"));
		if(suffix.indexOf("fastq") != -1){}
		else if(suffix.indexOf("csfasta") != -1){}
		else if(suffix.indexOf("fasta") != -1){
			FastaSequences A = new FastaSequences();
			A.readFastaFile(SequenceFile);
			A.printDistribution(gmapperFile, cutoff);
		}
		else{
			System.out.println("kind has to be either fastq, csfasta or fa");
		}
	}


	/*	public static void findDouble(Hashtable<String,String> T){
		boolean allRequired = true;
		String dir, outDir, SequenceFile, primerFile, suffix;
		dir = outDir = SequenceFile = primerFile  = suffix = null;
		if(T.containsKey("-id"))
			dir = Functions.getValue(T, "-id", "");
		else{
			allRequired = false;
			System.out.println("must have indirectory -id");

		}
		if(T.containsKey("-od"))
			outDir = Functions.getValue(T, "-od", "");
		else{
			allRequired = false;
			System.out.println("must have outdirectory -id");
		}
		if(T.contains("-i"))
			SequenceFile = Functions.getValue(T, "-i", "file1.fa");
		if(T.containsKey("-primerFile"))
			primerFile = Functions.getValue(T, "-primerFile", "primer.fa");
		else{
			allRequired = false;
			System.out.println("must have primerFile -primerFile");

		}
		suffix = Functions.getValue(T, "-suffix", "fastq");

		if(allRequired){
			FastaSequences primers = new FastaSequences(primerFile);
			if(suffix.indexOf("fastq") != -1){}
			else if(suffix.indexOf("csfasta") != -1){
				CfastaSequences csPrimers = CfastaSequences.convertFasta2CSfasta(primers);
				if(T.containsKey("i")){
					try{
						ExtendedWriter EW3 = new ExtendedWriter(new FileWriter(outDir+"/info.txt"));
						CfastaSequences.removePrimers(EW3,dir,outDir, SequenceFile,csPrimers,lowCutoff,highCutoff);
					}catch(Exception E){E.printStackTrace();}

				}else{
					CfastaSequences.removePrimersDir(dir,outDir,csPrimers,suffix,lowCutoff,highCutoff);
				}
			}
			else if(suffix.indexOf("fa") != -1){}
			else{
				System.out.println("kind has to be either fastq, csfasta or fa");
			}
		}
	}

	 */
	public static void integrity(Hashtable<String,String> T){
		System.out.println("Checking sequence integrity and saving the functional genes");
		boolean allRequired = true;
		String dir, outDir, SequenceFile, primerFile, suffix;
		dir = outDir = SequenceFile = primerFile  = suffix = null;
		dir = Functions.getValue(T, "-d", ".");
		if(!T.containsKey("-f1") && !T.containsKey("-f2"))
			allRequired = false;
		suffix = Functions.getValue(T, "-suffix", "fastq");

		if(allRequired){
			if(suffix.indexOf("fastq") != -1){
				if(T.containsKey("-f1") && T.containsKey("-f2")){
					String fileName1 = Functions.getValue(T, "-f1", ".");
					String fileName2 = Functions.getValue(T, "-f2", ".");
					FastQSequences.filterPairedFastQSequences(dir, fileName1, fileName2);
				}
				else if(T.containsKey("-f1")){
					//
				}
				//
			}
			else if(suffix.indexOf("csfasta") != -1){
				//CfastaSequences.removePrimersDir(dir);
			}
			else if(suffix.indexOf("fa") != -1){}
			else{
				System.out.println("kind has to be either fastq, csfasta or fa");
			}
		}
	}

	public static void check(Hashtable<String,String> T){
		System.out.println("checking sequencin integrity");
		boolean allRequired = true;
		String dir, outDir, SequenceFile, primerFile, suffix;
		dir = outDir = SequenceFile = primerFile  = suffix = null;
		dir = Functions.getValue(T, "-d", ".");
		if(!T.containsKey("-f1") && !T.containsKey("-f2"))
			allRequired = false;
		suffix = Functions.getValue(T, "-suffix", "fastaq");
		int length =  Integer.parseInt(Functions.getValue(T, "-l", "70"));


		if(allRequired){
			if(suffix.indexOf("fastq") != -1){
				if(T.containsKey("-f1") && T.containsKey("-f2")){
					String fileName1 = Functions.getValue(T, "-f1", ".");
					String fileName2 = Functions.getValue(T, "-f2", ".");
					FastQSequences.checkPairedFastQSequences(dir, fileName1, fileName2, length);
				}
				else if(T.containsKey("-f1")){
					//
				}
				//
			}
			else if(suffix.indexOf("csfasta") != -1){
				//CfastaSequences.removePrimersDir(dir);
			}
			else if(suffix.indexOf("fasta") != -1){

			}
			else{
				System.out.println("-suffix has to be either fastq, csfasta or fa");
			}
		}
	}

	public static void size(Hashtable<String,String> T){
		System.out.println("checking sequences size");
		String suffix = Functions.getValue(T, "-suffix", "fasta");

		if(suffix.indexOf("fastq") != -1){
			if(T.containsKey("-f1") && T.containsKey("-f2")){
				//					String fileName1 = Functions.getValue(T, "-f1", ".");
				//					String fileName2 = Functions.getValue(T, "-f2", ".");
				//					FastQSequences.checkPairedFastQSequences(dir, fileName1, fileName2);
			}
			else if(T.containsKey("-f1")){
				//
			}
			//
		}
		else if(suffix.indexOf("csfasta") != -1){
			String dir = Functions.getValue(T, "-d", ".");
			String fileName1 = Functions.getValue(T, "-f1", ".");
			int max = Integer.parseInt(Functions.getValue(T, "-max", "50"));
			CfastaSequences.printSizeDistribution(dir, fileName1, max);
		}
		else if(suffix.indexOf("fasta") != -1){
			String fileName1 = Functions.getValue(T, "-i", ".");
			System.out.println("Size is written to "+fileName1);
			FastaSequences.getLength(fileName1);
		}
		else{
			System.out.println("kind has to be either fastq, csfasta or fa");
		}
	}


	public static void extractAbove(Hashtable<String,String> T){
		System.out.println("checking sequences size");
		String suffix = Functions.getValue(T, "-suffix", "fasta");
		int length = Integer.parseInt(Functions.getValue(T, "-length", "0"));

		if(suffix.indexOf("fastq") != -1){
			if(T.containsKey("-f1") && T.containsKey("-f2")){
				//					String fileName1 = Functions.getValue(T, "-f1", ".");
				//					String fileName2 = Functions.getValue(T, "-f2", ".");
				//					FastQSequences.checkPairedFastQSequences(dir, fileName1, fileName2);
			}
			else if(T.containsKey("-f1")){
				//
			}
			//
		}
		else if(suffix.indexOf("csfasta") != -1){
			String dir = Functions.getValue(T, "-d", ".");
			String fileName1 = Functions.getValue(T, "-f1", ".");
			int max = Integer.parseInt(Functions.getValue(T, "-max", "50"));
			CfastaSequences.printSizeDistribution(dir, fileName1, max);
		}
		else if(suffix.indexOf("fasta") != -1){
			String inFile = Functions.getValue(T, "-i", ".");
			String outFile = Functions.getValue(T, "-o", inFile.substring(0,inFile.lastIndexOf("fa"))+length+".fa");
			String WD = Functions.getValue(T, "-d", IOTools.getCurrentPath());
			
			System.out.println("Extrackting sequence longer than "+length+" that is written to "+outFile);
			FastaSequences.extractSeqAboveLength(inFile,outFile,WD,length);
		}
		else{
			System.out.println("kind has to be either fastq, csfasta or fa");
		}
	}
	
	public static void extractSize(Hashtable<String,String> T){
		String suffix = Functions.getValue(T, "-suffix", "fa");
		int length = Integer.parseInt(Functions.getValue(T, "-l", "1"));
		System.out.println("Extracting "+ suffix +" sequences with size "+length);

		if(suffix.indexOf("fastq") != -1){
			if(T.containsKey("-f1") && T.containsKey("-f2")){
				//					String fileName1 = Functions.getValue(T, "-f1", ".");
				//					String fileName2 = Functions.getValue(T, "-f2", ".");
				//					FastQSequences.checkPairedFastQSequences(dir, fileName1, fileName2);
			}
			else if(T.containsKey("-i")){
				FastQSequences.fastqfixedLength(IOTools.getCurrentPath(), Functions.getValue(T, "-i"), length);
				//
			}
			//
		}
		else if(suffix.indexOf("csfasta") != -1){
			String dir = Functions.getValue(T, "-d", IOTools.getCurrentPath());
			String fileName1 = Functions.getValue(T, "-i", ".");
			int size = Integer.parseInt(Functions.getValue(T, "-max", "50"));
			CfastaSequences.printSizeDistribution(dir, fileName1, size);
		}
		else if(suffix.indexOf("fasta") != -1){
			String inFile = Functions.getValue(T, "-i", ".");
			String outFile = Functions.getValue(T, "-o", inFile.substring(0,inFile.lastIndexOf("fa"))+length+".fa");
			String WD = Functions.getValue(T, "-d", IOTools.getCurrentPath());
			
			System.out.println("Extrackting sequence longer than "+length+" that is written to "+outFile);
			FastaSequences.extractSeqAboveLength(inFile,outFile,WD,length);
		}
		else{
			System.out.println("kind has to be either fastq, csfasta or fa");
		}
	}
	
	public static void extractEnds(Hashtable<String,String> T){
		String suffix = Functions.getValue(T, "-suffix", "fasta");
		int length = Integer.parseInt(Functions.getValue(T, "-l", "300"));

		if(T.containsKey("-i")){
			String fileName1 = Functions.getValue(T, "-i", "test");

			if(suffix.indexOf("fastq") != -1){
				if(T.containsKey("-f1") && T.containsKey("-f2")){
					//					String fileName1 = Functions.getValue(T, "-f1", ".");
					//					String fileName2 = Functions.getValue(T, "-f2", ".");
					//					FastQSequences.checkPairedFastQSequences(dir, fileName1, fileName2);
				}
				else if(T.containsKey("-f1")){
					//
				}
				//
			}
			else if(suffix.indexOf("csfasta") != -1){
				//CfastaSequences.removePrimersDir(dir);
			}
			else if(suffix.indexOf("fasta") != -1){
				FastaSequences.extractEnds(fileName1, length);

			}
			else{				
				System.out.println("kind has to be either fastq, csfasta or fa");
			}
		}
		else{
			System.out.println("must contain a -i flag");
		}
	}

	public static void getORFs(Hashtable<String,String> T){
		String suffix = Functions.getValue(T, "-suffix", "fasta");
		String dir = Functions.getValue(T, "-d", IOTools.getCurrentPath());

		if(T.containsKey("-i")){
			String fileName1 = Functions.getValue(T, "-i", "test");

			if(suffix.indexOf("fastq") != -1){
				if(T.containsKey("-f1") && T.containsKey("-f2")){
					//					String fileName1 = Functions.getValue(T, "-f1", ".");
					//					String fileName2 = Functions.getValue(T, "-f2", ".");
					//					FastQSequences.checkPairedFastQSequences(dir, fileName1, fileName2);
				}
				else if(T.containsKey("-f1")){
					//
				}
				//
			}
			else if(suffix.indexOf("csfasta") != -1){
				//CfastaSequences.removePrimersDir(dir);
			}
			else if(suffix.indexOf("fasta") != -1){
				FastaSequences.getORFs(dir,fileName1);

			}
			else{				
				System.out.println("kind has to be either fastq, csfasta or fa");
			}
		}
		else{
			System.out.println("must contain a -i flag");
		}
	}

	
	
	public static void padEnds(Hashtable<String,String> T){
		String suffix = Functions.getValue(T, "-suffix", "fasta");
		String pad = Functions.getValue(T, "-pad", "ATGA");

		if(T.containsKey("-i")){
			String fileName1 = Functions.getValue(T, "-i", "test");

			if(suffix.indexOf("fastq") != -1){
				if(T.containsKey("-f1") && T.containsKey("-f2")){
					//					String fileName1 = Functions.getValue(T, "-f1", ".");
					//					String fileName2 = Functions.getValue(T, "-f2", ".");
					//					FastQSequences.checkPairedFastQSequences(dir, fileName1, fileName2);
				}
				else if(T.containsKey("-f1")){
					//
				}
				//
			}
			else if(suffix.indexOf("csfasta") != -1){
				//CfastaSequences.removePrimersDir(dir);
			}
			else if(suffix.indexOf("fasta") != -1){
				FastaSequences.extractEnds(fileName1, 2);

			}
			else{				
				System.out.println("kind has to be either fastq, csfasta or fa");
			}
		}
		else{
			System.out.println("must contain a -i flag");
		}
	}
	

	public static void extract(Hashtable<String,String> T){
		String suffix = Functions.getValue(T, "-suffix", "fasta");
		String name = Functions.getValue(T, "-name", "fasta");
		int start = Integer.parseInt(Functions.getValue(T, "-start", "0"));
		int stop = Integer.parseInt(Functions.getValue(T, "-stop", "1"));
		int surr = Integer.parseInt(Functions.getValue(T, "-surr", "0"));
		String extra = Functions.getValue(T, "-extra", "");


		if(suffix.indexOf("fastq") != -1){
			if(T.containsKey("-f1") && T.containsKey("-f2")){
				//					String fileName1 = Functions.getValue(T, "-f1", ".");
				//					String fileName2 = Functions.getValue(T, "-f2", ".");
				//					FastQSequences.checkPairedFastQSequences(dir, fileName1, fileName2);
			}
			else if(T.containsKey("-f1")){
				//
			}
			//
		}
		else if(suffix.indexOf("csfasta") != -1){
			//CfastaSequences.removePrimersDir(dir);
		}
		else if(suffix.indexOf("fasta") != -1){
			String fileName1 = Functions.getValue(T, "-f1", ".");
			FastaSequences A = new FastaSequences(fileName1);
			A.printSequenceSurr(name,start,stop, surr,extra);


		}
		else{
			System.out.println("kind has to be either fastq, csfasta or fa");
		}
	}


	public static void fixSequenceNames(Hashtable<String,String> T){
		System.out.println("Fixing sequencesize");
		String suffix = Functions.getValue(T, "-suffix", "fasta");

		if(suffix.indexOf("fastq") != -1){
			if(T.containsKey("-f1") && T.containsKey("-f2")){
				//					String fileName1 = Functions.getValue(T, "-f1", ".");
				//					String fileName2 = Functions.getValue(T, "-f2", ".");
				//					FastQSequences.checkPairedFastQSequences(dir, fileName1, fileName2);
			}
			else if(T.containsKey("-f1")){
				//
			}
			//
		}
		else if(suffix.indexOf("csfasta") != -1){
			//CfastaSequences.removePrimersDir(dir);
		}
		else if(suffix.indexOf("fasta") != -1){
			String fileName1 = Functions.getValue(T, "-i", ".");
			FastaSequences.fixFastaFile(fileName1);
		}
		else{
			System.out.println("kind has to be either fastq, csfasta or fa");
		}
	}

	public static void oneLine(Hashtable<String,String> T){
		System.out.println("Fixing sequencesize");
		String fileName1 = Functions.getValue(T, "-i", ".");
		FastaSequences.oneLine(fileName1);
	}

	public static void split(Hashtable<String,String> T){
		int n = Integer.parseInt(Functions.getValue(T, "-n", "1000000"));
		System.out.println("Splitting fastaSequences into parts with "+n +" sequences");
		String suffix = Functions.getValue(T, "-suffix", "fa");
		if(suffix.indexOf("fastq") > -1){
			String fileName1 = Functions.getValue(T, "-i", ".");
			if(IOTools.isDir(fileName1)){
				ArrayList <String> fileNames = IOTools.getSequenceFiles(fileName1,suffix);
				for(int i = 0; i < fileNames.size();i++){
					FastQSequences.split(fileName1+"/"+fileNames.get(i),n,suffix);
				}
			}
			else{
				FastQSequences.split(fileName1,n,suffix);
			}
		}
		else{
			String fileName1 = Functions.getValue(T, "-i", ".");
			if(IOTools.isDir(fileName1)){
				ArrayList <String> fileNames = IOTools.getSequenceFiles(fileName1,suffix);
				for(int i = 0; i < fileNames.size();i++){
					FastaSequences.split(fileName1+"/"+fileNames.get(i),n);
				}
			}
			else{
				FastaSequences.split(fileName1,n);
			}
		}
	}

	
	
	public static void one2two(Hashtable<String,String> T){
		if(T.containsKey("-i")){
			String fileName1 = Functions.getValue(T, "-i", ".");
			String suffix = Functions.getValue(T, "-suffix", "fastq");
			if(suffix.compareTo("fastq") == 0)
				FastQSequences.one2two(fileName1);
		}
	}



	public static void mutationalAnalysis(Hashtable<String,String> T){
		boolean allRequired = true;
		String dir,  muatantFile, mutantFile2 , Sequence, structure;
		dir = muatantFile = mutantFile2 = Sequence  = structure = null;
		if(T.containsKey("-d"))
			dir = Functions.getValue(T, "-d", "");
		else{
			allRequired = false;
		}
		String suffix = Functions.getValue(T, "-suffix", "fastq");
		int length =  Integer.parseInt(Functions.getValue(T, "-l", "70"));

		if(allRequired){
			if(suffix.indexOf("fastq") != -1){
				if(T.containsKey("-f1") && T.containsKey("-f2")){
					String fileName1 = Functions.getValue(T, "-f1", ".");
					String fileName2 = Functions.getValue(T, "-f2", ".");
					FastQSequences.checkPairedFastQSequences(dir, fileName1, fileName2,length);
				}
				else if(T.containsKey("-f1")){
					//
				}
				//
			}
			else if(suffix.indexOf("csfasta") != -1){
				//CfastaSequences.removePrimersDir(dir);
			}
			else if(suffix.indexOf("fa") != -1){}
			else{
				System.out.println("kind has to be either fastq, csfasta or fa");
			}
		}
	}

	public static void QC(Hashtable<String,String> T){
		System.out.println("Quality controll");
		boolean allRequired = true;
		String dir,  muatantFile, mutantFile2 , Sequence, structure;
		dir = muatantFile = mutantFile2 = Sequence  = structure = null;
		if(T.containsKey("-i"))
			dir = Functions.getValue(T, "-i", "");
		else if(T.containsKey("-dir"))
			dir = Functions.getValue(T, "-dir", "");
		else {
			dir = IOTools.getCurrentPath();
		}
		String suffix = Functions.getValue(T, "-suffix", "fastq");

		if(allRequired){
			if(suffix.indexOf("fastq") != -1){
				FastQSequences.QC(T);
			}
			else if(suffix.indexOf("csfasta") != -1){
				//CfastaSequences.removePrimersDir(dir);
			}
			else if(suffix.indexOf("fa") != -1){}
			else{
				System.out.println("kind has to be either fastq, csfasta or fa");
			}
		}
	}

	public static void stratify(Hashtable<String,String> T){
		System.out.println("Quality controll");
		boolean allRequired = true;
		String dir,  muatantFile, mutantFile2 , Sequence, structure;
		dir = muatantFile = mutantFile2 = Sequence  = structure = null;
		if(T.containsKey("-i"))
			dir = Functions.getValue(T, "-i", "");
		else if(!T.containsKey("-dir")){
			dir = IOTools.getCurrentPath();
		}
		String suffix = Functions.getValue(T, "-suffix", "fastq");

		if(allRequired){
			if(suffix.indexOf("fastq") != -1){
				FastQSequences.splitPair(T);
			}
			else if(suffix.indexOf("csfasta") != -1){
				//CfastaSequences.removePrimersDir(dir);
			}
//			else if(suffix.indexOf("fa") != -1){
//				FastaSequences.splitSequences(T);
//			}
			else{
				System.out.println("kind has to be either fastq, csfasta or fa");
			}
		}
	}
	
	
	
	public static void mutual(Hashtable<String,String> T){
		boolean allRequired = true;
		String seqFile1, seqFile2,outFile3;
		String dir = Functions.getValue(T, "-d", ".");
		seqFile1 = Functions.getValue(T, "-i1", "notNamned");
		seqFile2 = Functions.getValue(T, "-i2", "notNamned");
		outFile3 = Functions.getValue(T, "-o", "notNamned");

		if(seqFile1.compareTo(seqFile2) == 0) allRequired= false;
		if(seqFile1.compareTo(outFile3) == 0) allRequired= false;
		if(seqFile2.compareTo(outFile3) == 0) allRequired= false;

		String suffix = Functions.getValue(T, "-suffix", "fa");

		if(allRequired){
			if(suffix.indexOf("fastq") != -1){
				//				FastQSequences.QC(T);
			}
			else if(suffix.indexOf("csfasta") != -1){
				//CfastaSequences.removePrimersDir(dir);
			}
			else if(suffix.indexOf("fa") != -1){
				FastaSequences A = new FastaSequences(seqFile1);
				FastaSequences B = new FastaSequences(seqFile2);
				FastaSequences C = FastaSequences.getMutualTranscripts(A,B);
				C.printFastaRNA(dir, outFile3);
			}
			else{
				System.out.println("kind has to be either fastq, csfasta or fa");
			}
		}
		else{System.out.println("-i1, -i2 and -o must be specified");}
	}




	public static FastaSequences MergeSequences(Hashtable<String,String> T){
		String dir = Functions.getValue(T, "-d", "");
		String fastaFile1 = Functions.getValue(T, "-file1", "file1.fa");
		String fastaFile2 = Functions.getValue(T, "-file2", "file2.fa");

		FastaSequences FS1 = new FastaSequences(dir, fastaFile1);
		FastaSequences FS2 = new FastaSequences(dir, fastaFile2);

		FastaSequences Unique = FastaSequences.merge(FS1,FS2);
		String fastaFile3 = Functions.getValue(T, "-outFile", "outFile.fa");
		Unique.writeFastaFile(dir,fastaFile3);

		return Unique;
	}

	
	
	public static void getLengths(Hashtable<String,String> T){
		String dir = Functions.getValue(T, "-d", ".");
		String fastaFile1 = Functions.getValue(T, "-file1", "file1.fa");

		FastaSequences FS1 = new FastaSequences(dir, fastaFile1);

		int [] lengths = FS1.getLengths();
		for(int i = 0; i < lengths.length; i++)
			System.out.print(lengths[i]+",");
		System.out.println();
	}

	public static void print5Ends(Hashtable<String,String> T){
		String dir = Functions.getValue(T, "-d", ".");
		String fastaFile1 = Functions.getValue(T, "-file1", "file1.fa");

		FastaSequences FS1 = new FastaSequences(dir, fastaFile1);
		String outFile= Functions.getValue(T, "-o", dir+"/"+fastaFile1+ ".5end");
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outFile));
			FS1.print5end(EW);
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}

	public static void print3Ends(Hashtable<String,String> T){
		String dir = Functions.getValue(T, "-d", ".");
		String fastaFile1 = Functions.getValue(T, "-file1", "file1.fa");

		FastaSequences FS1 = new FastaSequences(dir, fastaFile1);
		String outFile= Functions.getValue(T, "-o", dir+"/"+fastaFile1+ ".3end");
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outFile));
			FS1.print3end(EW);
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}

	public static int matchSequence(int[] alignSeq, int[] RefSeq, FastaSequence seq ,  int errors){
		boolean found = false;
		int count = 0;
		int wrong  = 0;
		if(RefSeq.length != alignSeq.length) return -1;

		while(wrong < errors && count < RefSeq.length){
			if(alignSeq[count] == RefSeq[count])
				count++;
			else{
				wrong++;
				count++;
			}
		}
		if(wrong < errors){
			seq.Sequence = alignSeq;
			return wrong;
		}
		else
			return -1;

	}

	public static FastQSequence matchSequence(int[] RefSeq, FastQSequence seq ,  int errors){
		boolean found = false;

		int count2 = 0;
		int count = 0;
		int wrong  = 0;
		int[] alignSeq = RNAfunctions.RNAString2Int(seq.sequence);
		while (!found &&  RefSeq.length + count2 <= alignSeq.length){
			count = 0;
			while(wrong < errors && count < RefSeq.length){
				if(alignSeq[count + count2] == RefSeq[count])
					count++;
				else{
					wrong++;
					count++;
				}
			}
			if(wrong < errors && count == RefSeq.length ){
				seq.sequence = seq.sequence.substring(count2,count+count2);
				seq.quality = seq.quality.substring(count2, count+count2);
				found = true;
				return seq;
			}
			count2++;
			wrong = 0;
		}
		return matchReverseSequence(RefSeq, seq , errors);
	}

	public static FastQSequence matchReverseSequence(int[] RefSeq, FastQSequence seq ,  int errors){
		boolean found = false;

		int count2 = 0;
		int count = 0;
		int wrong  = 0;
		int[] alignSeq = RNAfunctions.getReverseComplement(RNAfunctions.RNAString2Int(seq.sequence));

		while (!found &&  RefSeq.length + count2 < alignSeq.length){
			count = 0;
			while(wrong < errors && count < RefSeq.length){
				if(alignSeq[count + count2] == RefSeq[count])
					count++;
				else{
					wrong++;
					count++;
				}
			}

			if(wrong < errors && count == RefSeq.length ){
				alignSeq = RNAfunctions.getSubsequence(alignSeq, count2, count+count2);
				seq.sequence = RNAfunctions.DNAInt2String(alignSeq);
				char[] array = seq.quality.toCharArray();
				char[] array2 = new char[alignSeq.length];
				for(int i = 0  ; i < array2.length; i++){
					array2[i] = array[array.length-1-count2-i];
				}
				seq.quality = new String(array2);
				found = true;
				return seq;
			}
			count2++;
			wrong = 0;
		}
		//		System.out.println("Seq Not Found");
		//		System.out.println(seq.sequence);
		//		System.out.println( RNAfunctions.DNAInt2String(RefSeq));
		//		System.out.println();
		return null;
	}

	public static int getNrOfErrors(int[] RefSeq, FastQSequence seq){
		boolean found = false;

		int count2 = 0;
		int count = 0;
		int wrong  = 0;
		int[] alignSeq = RNAfunctions.RNAString2Int(seq.sequence);
		for(int i = 0; i < RefSeq.length; i++){
			if(alignSeq[i] != RefSeq[i])
				wrong++;
		}
		return wrong;
	}



}
