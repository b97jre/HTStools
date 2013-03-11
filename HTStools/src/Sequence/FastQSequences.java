package Sequence;

import general.Functions;
import general.IOTools;
import general.RNAfunctions;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;

import general.ExtendedReader;
import general.ExtendedWriter;




public class FastQSequences implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//private String GenomeName;
	private String Name;
	ArrayList <FastQSequence> forward;
	ArrayList <FastQSequence> reverse;


	public static void main(String[] args) {
		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}
		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);
		run(T);
	}


	public static void run(Hashtable<String,String> T){
		String[] experiments = null;
		if(T.containsKey("-p"))
			experiments = T.get("-p").split(" ");
		for(int i = 0; i < experiments.length; i++){
			if(experiments[i].toUpperCase().compareTo("checkSeq".toUpperCase()) == 0){
				checkSeq(T);
			}
			if(experiments[i].toUpperCase().compareTo("filterSeq".toUpperCase()) == 0){
				checkSeq(T);
			}

			if(experiments[i].toUpperCase().compareTo("QC".toUpperCase()) == 0){
				QC(T);
			}

			if(experiments[i].toUpperCase().compareTo("unify".toUpperCase()) == 0){
				unify(T);
			}

			if(experiments[i].toUpperCase().compareTo("split".toUpperCase()) == 0){
				splitPair(T);
			}





		}
	}

	public static void split(String inFile, int nrOfSequences,String suffix){

		IOTools.mkDir("split");
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (inFile));

			String [] outFile  = inFile.split(suffix);
			int temp = nrOfSequences;
			int count = 0;
			int number = 0;
			while(ER.more()){
				ExtendedWriter EW = new ExtendedWriter(new FileWriter ("split/"+outFile[0]+"_"+number+"."+suffix));
				while(ER.more() && count < nrOfSequences){
					EW.println(ER.readLine());
					EW.println(ER.readLine());
					EW.println(ER.readLine());
					EW.println(ER.readLine());
					count++;
				}
				EW.flush();
				EW.close();
				nrOfSequences+= temp;
				number++;
			}
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}


	private static void checkSeq(Hashtable<String,String> T){
		String dir = Functions.getValue(T, "-d", ".");
		String fileName1, fileName2;
		fileName1 = fileName2 = null;
		int length =  Integer.parseInt(Functions.getValue(T, "-l", "70"));
		if(T.containsKey("-f1") && T.containsKey("-f2")){
			fileName1 = Functions.getValue(T, "-f1", ".");
			fileName2 = Functions.getValue(T, "-f2", ".");
			checkPairedFastQSequences(dir, fileName1, fileName2, length);

		}
		else if(T.containsKey("-f1")){
			fileName1 = Functions.getValue(T, "-f1", ".");
		}
	}

	public static void QC(Hashtable<String,String> T){
		boolean allRequired = true;
		String dir = null;
		int length =  Integer.parseInt(Functions.getValue(T, "-l", "70"));
		if(T.containsKey("-i"))
			dir = Functions.getValue(T, "-i", "");
		else if(T.containsKey("-dir"))
			dir = Functions.getValue(T, "-dir", "");
		else {
			dir = IOTools.getCurrentPath();
		}
		if(T.containsKey("-f1") && T.containsKey("-f2")){
			QC(dir, Functions.getValue(T, "-o", dir+"/QC"), Functions.getValue(T, "-f1"), Functions.getValue(T, "-f2"), length);
			return;
		}
		else
		{
			allRequired = false;
		}
		String outDir = Functions.getValue(T, "-o", dir+"/QC");

		String suffix = Functions.getValue(T, "-suffix", "fastq");

		if(allRequired){
			if(suffix.indexOf("fastq") != -1){
				System.out.println("Checking directory "+dir);
				ArrayList <String> fileNames = IOTools.getSequenceFiles(dir,suffix);
				String [] sep = new String[2];
				sep[0] = "1.fastq";
				sep[1] = "2.fastq";
				ArrayList <String[]> pairs = IOTools.findPairs(fileNames,sep);
				System.out.println("nr of pairs "+pairs.size());
				for(int i = 0; i < pairs.size(); i++){
					System.out.println(pairs.get(i)[0]+" "+pairs.get(i)[1]);
					QC(dir, outDir, pairs.get(i)[0], pairs.get(i)[1], length);
					System.out.println();
				}
			}
		}
		else{
			System.out.println("Not sufficient information to run program. You must either name a folder (-i) with pairs of fastq files that end with 1.fastq and 2.fastq or two files (-f1 and -f2) with the pairs that you want to QC");
			System.out.println("-i <folder> \t\tAll pairs withing folder will be QC. overrides  -f1 and -f2");
			System.out.println("-f1 <forward file> \tFirst file in fastq pair file");
			System.out.println("-f2 <reverse file> \tSecond file in fastq pair file");
			System.out.println("-o <folder> \t\tfolder were the checked fastq files end up. DEFAULT (\"current_folder/QC\").");
			System.out.println("-l <integer> \t\tLeast length needed to keep pair. DEFAULT (70).");

		}
	}


	public static void fixSequenceLengths(Hashtable<String,String> T){
		boolean allRequired = true;
		String dir = null;
		if(T.containsKey("-i"))
			dir = Functions.getValue(T, "-i", "");
		else {
			dir = IOTools.getCurrentPath();
		}
		String outDir = Functions.getValue(T, "-o", dir+"/QC");

		String suffix = Functions.getValue(T, "-suffix", "fastq");

		if(allRequired){
			if(suffix.indexOf("fastq") != -1){
				
				System.out.println("Checking directory "+dir);
				ArrayList <String> fileNames = IOTools.getSequenceFiles(dir,suffix);
				System.out.println("nr of Sequences "+fileNames.size());
				for(int i = 0; i < fileNames.size(); i++){
					System.out.println(fileNames.get(i));
					fixSequenceLengths(dir, outDir, fileNames.get(i));
					System.out.println();
				}
			}
		}
		else{
			System.out.println("Not sufficient information to run program. You must either name a folder (-i) with pairs of fastq files that end with 1.fastq and 2.fastq or two files (-f1 and -f2) with the pairs that you want to QC");
			System.out.println("-i <folder> \t\tAll pairs withing folder will be QC. overrides  -f1 and -f2");
			System.out.println("-f1 <forward file> \tFirst file in fastq pair file");
			System.out.println("-f2 <reverse file> \tSecond file in fastq pair file");
			System.out.println("-o <folder> \t\tfolder were the checked fastq files end up. DEFAULT (\"current_folder/QC\").");
			System.out.println("-l <integer> \t\tLeast length needed to keep pair. DEFAULT (70).");

		}
	}
	
	
	public static void splitPair(Hashtable<String,String> T){
		boolean allRequired = true;
		String dir = null;
		String [] sep = new String[2];
		sep[0] = Functions.getValue(T, "-s1", "1.fastq");;
		sep[1] = Functions.getValue(T, "-s2", "2.fastq");;
		String in =null; 
		if(T.containsKey("-i"))
			in = Functions.getValue(T, "-i", ".");
		else
		{
			System.out.println("must contain a infile or inDir (-i )");
			allRequired = false;
		}

		if(!T.containsKey("-s1") || !T.containsKey("-s2")){
			System.out.println("\"1.fastq\" and \"2.fastq\" are set as suffix for paired reds");
		}
		System.out.println(IOTools.getCurrentPath()+"/"+in+sep[0]);
		if(IOTools.fileExists(IOTools.getCurrentPath()+"/"+in+sep[0])){
			splitPair(IOTools.getCurrentPath(), Functions.getValue(T, "-o", IOTools.getCurrentPath()+"/split"), in, sep);
			return;
		}
		dir = in;
		String outDir = Functions.getValue(T, "-o", dir+"/QC");

		String suffix = Functions.getValue(T, "-suffix", "fastq");

		if(allRequired){
			if(suffix.indexOf("fastq") != -1){
				System.out.println("Checking directory "+dir);
				ArrayList <String> fileNames = IOTools.getSequenceFiles(dir,suffix);
				ArrayList <String[]> pairs = IOTools.findPairs(fileNames,sep);
				System.out.println("nr of pairs "+pairs.size());

				for(int i = 0; i < pairs.size(); i++){
					String[] fileName = pairs.get(i)[0].split(sep[0]);

					System.out.println(pairs.get(i)[0]+" "+pairs.get(i)[1]);
					splitPair(dir, outDir, fileName[0], sep);
					System.out.println();
				}
			}
		}
		else{
			System.out.println("Not sufficient information to run program. You must either name a folder (-i) with pairs of fastq files that end with 1.fastq and 2.fastq or two files (-f1 and -f2) with the pairs that you want to QC");
			System.out.println("-i <folder> \t\tAll pairs withing folder will be QC. overrides  -f1 and -f2");
			System.out.println("-f1 <forward file> \tFirst file in fastq pair file");
			System.out.println("-f2 <reverse file> \tSecond file in fastq pair file");
			System.out.println("-o <folder> \t\tfolder were the checked fastq files end up. DEFAULT (\"current_folder/QC\").");
			System.out.println("-l <integer> \t\tLeast length needed to keep pair. DEFAULT (70).");

		}
	}

	public static void unify(Hashtable<String,String> T){
		boolean allRequired = true;
		String dir = null;
		if(T.containsKey("-i"))
			dir = Functions.getValue(T, "-i", ".");
		else{
			allRequired = false;
		}
		String outDir = Functions.getValue(T, "-o", dir+"/UniSeq");
		int[] RefSeq =  RNAfunctions.RNAString2Int(Functions.getValue(T, "-seq", "AAAAAAAAAAA"));
		String suffix = Functions.getValue(T, "-suffix", "fastq");
		int errors = Integer.parseInt(Functions.getValue(T, "-e", "5"));

		if(allRequired){
			if(suffix.indexOf("fastq") != -1){
				ArrayList <String> fileNames = IOTools.getSequenceFiles(dir,suffix);
				for(int i = 0; i < fileNames.size(); i++){
					System.out.println(fileNames.get(i));
					Unify(dir, outDir, fileNames.get(i), RefSeq,errors);
				}
			}
		}
	}



	public static int countSequencesSR(String dir, ArrayList <String> pairs){
		int total = 0;
		for(int i = 0; i < pairs.size(); i++){
			total += countSequences(dir, pairs.get(i));
		}
		return total;
	}

	public static int countSequencesPE(String dir, ArrayList <String[]> pairs){
		int total = 0;
		for(int i = 0; i < pairs.size(); i++){
			total += countSequences(dir, pairs.get(i)[0]);
		}
		return total*2;
	}

	public static int countSequences(String dir, String file){
		try{
			return IOTools.countLines(dir+"/"+file);
		}
		catch(Exception E){
			E.printStackTrace(); 
			System.out.println("Could not read nr of lines on file "+dir+"/"+file);
		}
		return -1;
	}	


	public static void fastq2fasta(Hashtable<String,String> T){
		boolean allRequired = true;
		String dir = null;
		if(T.containsKey("-i"))
			dir = Functions.getValue(T, "-i", ".");
		else{
			allRequired = false;
		}
		String suffix = Functions.getValue(T, "-suffix", "fastq");

		if(allRequired){
			if(suffix.indexOf("fastq") != -1){
				ArrayList <String> fileNames = IOTools.getSequenceFiles(dir,suffix);
				for(int i = 0; i < fileNames.size(); i++){
					System.out.println(fileNames.get(i));
					fastq2fasta(dir, fileNames.get(i));
				}
			}
		}
	}


	public FastQSequences(){
		this.Name = "temp";
	}

	public FastQSequences(String Name){
		this.Name = Name;
	}


	public FastQSequences(String dir, String file){
		addFastQSequences(dir+"/"+file);
	}


	private void addForward(FastQSequence newSeq){
		this.forward.add(newSeq);
	}

	private void addReverse(FastQSequence newSeq){
		this.reverse.add(newSeq);
	}



	private void addFastQSequences(String fileName1){
		try{
			ExtendedReader ER1 = new ExtendedReader(new FileReader (fileName1));
			int nrOfGood = 0;
			while(ER1.more()){
				FastQSequence A = new FastQSequence();
				if(A.addInfo(ER1) ){
					nrOfGood++;
					this.forward.add(A);
				}
				else{
					System.out.println("this sequence is bad");
					System.out.println("seq 1 "+A.name );
				}
				if(forward.size() > nrOfGood){
					System.out.println("nr Of Good Sequences Reads " + nrOfGood);
					nrOfGood= nrOfGood + 1000000;
				}
			}
			ER1.close();
		}catch(Exception E){E.printStackTrace();}
	}

	private void addFastQSequences(String fileName1, String fileName2){
		try{
			ExtendedReader ER1 = new ExtendedReader(new FileReader (fileName1));
			ExtendedReader ER2 = new ExtendedReader(new FileReader (fileName2));
			int nrOfGood = 0;
			int count = 0;
			while(ER1.more()){
				count++;
				FastQSequence A = new FastQSequence();
				FastQSequence B = new FastQSequence();
				if(A.addInfo(ER1) && B.addInfo(ER2) && A.isPair(B)){
					nrOfGood++;
					this.forward.add(A);
					this.reverse.add(B);
				}
				else{
					System.out.println("this sequence is bad");
					if(!A.QC()){System.out.println("seq 1 "+A.name );}
					if(!B.QC()){System.out.println("seq 2 "+B.name );}
					if(!A.isPair(B)){System.out.println(B.name +" is not the pair sequence of "+A.name );}
				}
				if(forward.size() > nrOfGood){
					System.out.println("nr Of Sequences Read" + nrOfGood);
					nrOfGood= nrOfGood + 1000000;
				}
			}
			System.out.println("Total number of sequences = "+ count);
			System.out.println("Accepted number of sequences = "+ nrOfGood);
			System.out.println("Percentage number of sequences = "+ (nrOfGood*100/count) + "% ");


			ER1.close();
			ER2.close();
		}catch(Exception E){E.printStackTrace();}
	}

	public static void one2two(String fileName1){
		try{
			ExtendedReader ER1 = new ExtendedReader(new FileReader (fileName1));
			String base = fileName1.substring(0,fileName1.indexOf(".fastq"));
			ExtendedWriter EW1 = new ExtendedWriter(new FileWriter(base+".1.fastq"));
			ExtendedWriter EW2 = new ExtendedWriter(new FileWriter(base+".2.fastq"));
			int nrOfGood = 0;
			int total = 1000000;
			int count = 0;
			while(ER1.more()){
				count++;
				FastQSequence A = new FastQSequence();
				FastQSequence B = new FastQSequence();
				if(A.addInfo(ER1) && B.addInfo(ER1) && A.isPair(B)){
					nrOfGood++;
					A.printFastQ(EW1);
					B.printFastQ(EW2);
				}
				else{
					System.out.println("this sequence is bad");
					if(!A.QC()){System.out.println("seq 1 "+A.name );}
					if(!B.QC()){System.out.println("seq 2 "+B.name );}
					if(!A.isPair(B)){System.out.println(B.name +" is not the pair sequence of "+A.name );}
				}
				if(count > total){
					System.out.println("nr Of good Sequences Read" + nrOfGood);
					total= total + 1000000;
				}
			}
			System.out.println("Total number of sequences = "+ count);
			System.out.println("Accepted number of sequences = "+ nrOfGood);
			System.out.println("Percentage number of sequences = "+ (nrOfGood*100/count) + "% ");


			ER1.close();
			EW1.flush();
			EW1.close();
			EW2.flush();
			EW2.close();
		}catch(Exception E){E.printStackTrace();}
	}


	public void findPrimers(String primerSequence){
		int [] primer = RNAfunctions.RNAString2Int(primerSequence);
		int [] counter = new int[50];
		int start = 0;
		int count = 0;
		int maxLength = 15;
		int minMatches = 3;
		for(int i = 0; i < forward.size(); i++){
			int location = forward.get(i).findPrimer(2, start, maxLength, minMatches, primer);
			if(location>0){
				count++;
				counter[location]++;
			}
		}

		for(int i = start; i < counter.length; i++){
			System.out.println(i + "\t"+ counter[i]) ;
		}

		System.out.println(count);
		System.out.println(forward.size());
	}


	public int[] findPrimers(String[] primerSequence, int[] counter,FastQSequences notFound, int minLength){
		int [][] primers = RNAfunctions.fasta2CFasta(primerSequence);
		int start = 0;
		int count = 0;
		int maxLength = 15;
		int minMatches = 3;
		System.out.println("Sequences before :"+forward.size());
		for(int i = 0; i < forward.size(); i++){
			boolean found = false; 
			int j = 0;
			while(!found && j < primers.length){
				int location = forward.get(i).findPrimer(2, start, maxLength, minMatches, primers[j]);
				if(location < minLength){
					counter[location]++;
					found = true;
					forward.get(i).removePrimer(location);
				}
				j++;
			}
			if(!found){
				notFound.addForward(forward.get(i));
				forward.remove(i);
				i--;
				count++;
			}	
		}
		System.out.println("Sequences after +:"+forward.size());
		System.out.println("Sequences removed +:"+count);
		return counter;

	}



	public static boolean removePrimersDir(String inDir, String outDir, FastaSequences primers, String suffix){

		if(!IOTools.isDir(outDir))
			IOTools.mkDir(outDir);

		ArrayList <String> fileNames = IOTools.getSequenceFiles(inDir,suffix);

		if(fileNames.isEmpty()){
			System.out.println("No "+ suffix+"  files in folder :"+inDir);
			return false;
		}
		else{
			if(!IOTools.isDir(outDir+"/reports"))
				IOTools.mkDir(outDir+"/reports");
			if(!IOTools.isDir(outDir+"/scripts"))
				IOTools.mkDir(outDir+"/scripts");
			try{
				ExtendedWriter EW3 = new ExtendedWriter(new FileWriter(outDir+"/info.txt"));

				EW3.println("Sequences above length 15 and below length 30 is put i folder short");
				EW3.println("Sequences above length 29 is put i folder long");
				EW3.println();
				EW3.println("Distribution of sequences are:");

				EW3.print("cutoff SequenceFile ");
				for(int i = 0; i< 35; i++){
					EW3.print(i+" ");
				}
				EW3.println();


				for(int i = 0; i < fileNames.size(); i++){	
					System.out.println("removing primers in file "+fileNames.get(i));
					removePrimers(EW3,inDir, outDir, fileNames.get(i),primers);
				}
				EW3.flush();
				EW3.close();
			}
			catch(Exception E){E.printStackTrace();}


		}
		return true;
	}

	public static boolean removePrimers(ExtendedWriter EW3, String inDir, String outDir, String SequenceFile,FastaSequences primers){

		try{
			//System.out.println();
			double cutoff = 0.85;
			//			for(double cutoff = 0.0; cutoff < 1.0; cutoff = cutoff+0.1){
			ExtendedReader ER = new ExtendedReader(new FileReader (inDir+"/"+SequenceFile));
			if(!IOTools.isDir(outDir)) 	IOTools.mkDir(outDir);
			if(!IOTools.isDir(outDir+"/short")) IOTools.mkDir(outDir+"/short");
			if(!IOTools.isDir(outDir+"/long")) 	IOTools.mkDir(outDir+"/long");

			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/short/"+SequenceFile));
			ExtendedWriter EW2 = new ExtendedWriter(new FileWriter(outDir+"/long/"+SequenceFile));

			int[] dist = new int[36];
			while(ER.more()){
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					FastQSequence sequence = new FastQSequence();
					if(sequence.addInfo(ER)){
						dist[sequence.removePrimer(primers, cutoff)]++;
						sequence.printFastQ(EW,15,30);
						sequence.printFastQ(EW2,30);
					}
				}
			}
			int total = 0;
			EW3.print(cutoff+" "+SequenceFile+" ");
			for(int i = 0; i< dist.length; i++){
				EW3.print(dist[i]+" ");
				total = total + dist[i];
			}

			EW3.println(total);
			ER.close();
			EW.flush();
			EW.close();
			//			}

		}
		catch(Exception E){E.printStackTrace();}
		return false;

	}



	public static void checkPairedFastQSequences(String dir, String fileName1, String fileName2, int length){

		checkPairedFastQSequences(dir+"/"+fileName1,dir+"/"+fileName2, length );
	}




	private static void checkPairedFastQSequences(String fileName1, String fileName2, int length){
		try{
			ExtendedReader ER1 = new ExtendedReader(new FileReader (fileName1));
			ExtendedReader ER2 = new ExtendedReader(new FileReader (fileName2));
			int count = 0;
			int nrOfGood = 0;
			int checked = 0;
			while(ER1.more()){
				count++;
				FastQSequence A = new FastQSequence();
				FastQSequence B = new FastQSequence();
				boolean info1 = A.addInfo(ER1,length);
				boolean info2 = B.addInfo(ER2,length);
				if(info1 && info2 && A.isPair(B)){
					nrOfGood++;
				}
				//				else{
				//					System.out.println("this sequence is bad");
				//					if(!A.QC()){
				//						System.out.println("seq 1 "+A.name );
				//					}
				//					if(!B.QC()){
				//						System.out.println("seq 2 "+B.name );
				//					}
				//					if(!A.isPair(B)){System.out.println(B.name +" is not the pair sequence of "+A.name );}
				//				}
				if(count > checked){
					System.out.println("nr Of Sequences checked " + nrOfGood);
					checked= checked + 1000000;
				}

			}
			System.out.println("Total number of sequences = "+ count);
			System.out.println("Accepted number of sequences = "+ nrOfGood);
			System.out.println("Percentage of good sequences = "+ (nrOfGood*100/count) + "% ");
			ER1.close();
			ER2.close();
		}catch(Exception E){E.printStackTrace();}
	}

	public static void splitPair(String inDir, String outDir, String fileName,String[] extension){
		try{
			ExtendedReader ER1 = new ExtendedReader(new FileReader (inDir+ "/" + fileName+extension[0]));
			ExtendedReader ER2 = new ExtendedReader(new FileReader (inDir+ "/" + fileName+extension[1]));

			if(!IOTools.isDir(outDir) ){
				IOTools.mkDir(outDir);
			}
			ExtendedWriter ORF1 = new ExtendedWriter(new FileWriter (outDir+ "/"+ fileName+"ORF"+extension[0]));
			ExtendedWriter ORF2 = new ExtendedWriter(new FileWriter (outDir+ "/"+ fileName+"ORF"+extension[1]));
			ExtendedWriter FUTR1 = new ExtendedWriter(new FileWriter (outDir+ "/"+ fileName+"5UTR"+extension[0]));
			ExtendedWriter FUTR2 = new ExtendedWriter(new FileWriter (outDir+ "/"+ fileName+"5UTR"+extension[1]));
			ExtendedWriter TUTR1 = new ExtendedWriter(new FileWriter (outDir+ "/"+ fileName+"3UTR"+extension[0]));
			ExtendedWriter TUTR2 = new ExtendedWriter(new FileWriter (outDir+ "/"+ fileName+"3UTR"+extension[1]));
			ExtendedWriter NCRNA1 = new ExtendedWriter(new FileWriter (outDir+ "/"+ fileName+"NCRNA"+extension[0]));
			ExtendedWriter NCRNA2 = new ExtendedWriter(new FileWriter (outDir+ "/"+ fileName+"NCRNA"+extension[1]));
			ExtendedWriter WEIRD1 = new ExtendedWriter(new FileWriter (outDir+ "/"+ fileName+"WEIRD"+extension[0]));
			ExtendedWriter WEIRD2 = new ExtendedWriter(new FileWriter (outDir+ "/"+ fileName+"WEIRD"+extension[1]));

			int count = 0;
			int nrOfGood = 0;
			int nrOfORFs = 0 ;
			int nrOf5UTRs = 0;
			int nrOf3UTRs = 0;
			int nrOfncRNAs = 0;
			int nrOfWeirds = 0;
			int nrOfBad = 0;
			int checked = 0;
			int dot = 100000;
			int M = 1;
			
			while(ER1.more()){
				count++;
				FastQSequence A = new FastQSequence();
				FastQSequence B = new FastQSequence();
				boolean info1 = A.addInfo(ER1);
				boolean info2 = B.addInfo(ER2);
				if(info1 && info2 && A.isPair(B)){
					nrOfGood++;
				}
				int[][] Ainfo = A.getORFinfo();
				int[][] Binfo = B.getORFinfo();
				if(Ainfo[1][2] ==  0 && Binfo[1][2] == 0){
					A.printFastQ(ORF1);
					B.printFastQ(ORF2);
					nrOfORFs++;
				}
				else{ 
					A.printFastQ(NCRNA1);
					B.printFastQ(NCRNA2);
					nrOfncRNAs++;
				}
				if(count>dot){
					System.out.print(".");
					dot = dot+100000;
				}
				if(count>M*1000000){
					System.out.print(M+"M");
					M++;
				}
			}



			System.out.println("Total number of sequences = "+ count);
			System.out.println("Total number of good sequences = "+ count);
			System.out.println("Total number of ORF sequences"+ nrOfORFs);
			System.out.println("Total number of 5UTR sequences"+ nrOf5UTRs);
			System.out.println("Total number of 3UTR sequences"+ nrOf3UTRs);
			System.out.println("Total number of NCRNA sequences"+ nrOfncRNAs);
			System.out.println("Total number of WEIRD sequences"+ nrOfWeirds);
			System.out.println("Accepted number of sequences = "+ nrOfGood);
			System.out.println("Percentage of accepted sequences = "+ (nrOfGood*100/count) + "% ");
			System.out.println();

			ER1.close();
			ER2.close();
			ORF1.flush();
			ORF1.close();
			ORF2.flush();
			ORF2.close();
			FUTR1.flush();
			FUTR1.close();
			TUTR2.flush();
			TUTR2.close();
			NCRNA1.flush();
			NCRNA1.close();
			NCRNA2.flush();
			NCRNA2.close();
			WEIRD1.flush();
			WEIRD1.close();
			WEIRD2.flush();
			WEIRD2.close();
		}catch(Exception E){E.printStackTrace();}
	}

	public static void QC(String inDir, String outDir, String fileName1, String fileName2, int length){
		try{
			ExtendedReader ER1 = new ExtendedReader(new FileReader (inDir+ "/" + fileName1));
			ExtendedReader ER2 = new ExtendedReader(new FileReader (inDir+ "/" + fileName2));

			if(!IOTools.isDir(outDir) ){
				IOTools.mkDir(outDir);
			}
			ExtendedWriter EW1 = new ExtendedWriter(new FileWriter (outDir+ "/"+ fileName1));
			ExtendedWriter EW2 = new ExtendedWriter(new FileWriter (outDir+ "/"+ fileName2));

			int count = 0;
			int nrOfGood = 0;
			int nrOfBad = 0;
			int checked = 0;
			while(ER1.more()){
				count++;
				FastQSequence A = new FastQSequence();
				FastQSequence B = new FastQSequence();
				boolean info1 = A.addInfo(ER1, length);
				boolean info2 = B.addInfo(ER2, length);
				if(info1 && info2 && A.isPair(B)){
					nrOfGood++;
					A.printFastQ(EW1);
					B.printFastQ(EW2);
				}
				else{
					nrOfBad++;
				}

			}


			System.out.println("Total number of sequences = "+ count);
			System.out.println("Accepted number of sequences = "+ nrOfGood);
			System.out.println("Percentage of accepted sequences = "+ (double)(nrOfGood*100/count) + "% ");
			System.out.println();

			ER1.close();
			ER2.close();
			EW1.flush();
			EW1.close();
			EW2.flush();
			EW2.close();
		}catch(Exception E){E.printStackTrace();}
	}

	
	
	public static void QC(String inDir, String outDir, String fileName1, int length){
		try{
			ExtendedReader ER1 = new ExtendedReader(new FileReader (inDir+ "/" + fileName1));

			if(!IOTools.isDir(outDir) ){
				IOTools.mkDir(outDir);
			}
			ExtendedWriter EW1 = new ExtendedWriter(new FileWriter (outDir+ "/"+ fileName1));

			int count = 0;
			int nrOfGood = 0;
			while(ER1.more()){
				count++;
				FastQSequence A = new FastQSequence();
				boolean info1 = A.addInfo(ER1, length);
				if(info1){
					nrOfGood++;
					A.printFastQ(EW1);
				}
			}


			System.out.println("Total number of sequences = "+ count);
			System.out.println("Accepted number of sequences = "+ nrOfGood);
			System.out.println("Percentage of accepted sequences = "+ (double)(nrOfGood*100/count) + "% ");
			System.out.println();

			ER1.close();
			EW1.flush();
			EW1.close();
		}catch(Exception E){E.printStackTrace();}
	}
	
	
	public static void fixSequenceLengths(String inDir,String outDir, String fileName1){
		try{
			ExtendedReader ER1 = new ExtendedReader(new FileReader (inDir+ "/" + fileName1));
			if(!IOTools.isDir(outDir))IOTools.mkDir(outDir);
				ExtendedWriter EW1 = new ExtendedWriter(new FileWriter (outDir+ "/" + fileName1));

			int count = 0;
			int nrOfGood = 0;
			while(ER1.more()){
				count++;
				FastQSequence A = new FastQSequence();
				boolean info1 = A.addInfo2(ER1);
				if(info1){
					nrOfGood++;
					A.printFastQ(EW1);
				}

			}


			System.out.println("Total number of sequences = "+ count);
			System.out.println("Accepted number of sequences = "+ nrOfGood);
			System.out.println("Percentage of accepted sequences = "+ (double)(nrOfGood*100/count) + "% ");
			System.out.println();

			ER1.close();
			EW1.flush();
			EW1.close();
		}catch(Exception E){E.printStackTrace();}
	}
	
	
	private static void Unify(String inDir, String outDir, String fileName1, int[] RefSeq, int errors){
		try{
			ExtendedReader ER1 = new ExtendedReader(new FileReader (inDir+ "/" + fileName1));

			if(!IOTools.isDir(outDir) ){
				IOTools.mkDir(outDir);
			}
			ExtendedWriter EW1 = new ExtendedWriter(new FileWriter (outDir+ "/"+ fileName1));
			//ExtendedWriter EW2 = new ExtendedWriter(new FileWriter (outDir+ "/"+ fileName1+".perfect"));

			int[] errorDistribution = new int[errors];
			int count = 0;
			int nrOfGood = 0;
			int nrOfBad = 0;
			int checked = 0;
			while(ER1.more()){
				count++;
				FastQSequence A = new FastQSequence();
				boolean info1 = A.addInfo(ER1);
				A = SequenceHandling.matchSequence(RefSeq, A, errors);
				if(A != null){
					int info = SequenceHandling.getNrOfErrors(RefSeq, A);
					nrOfGood++;
					errorDistribution[info]++;
					if(info != 0){
						A.printFastQ(EW1);
					}
				}
				else{
					nrOfBad++;
				}
				if(count > checked){
					//System.out.println("nr Of Sequences checked " + nrOfGood);
					checked= checked + 1000000;
				}

			}


			System.out.println("Total number of sequences = "+ count);
			System.out.println("Accepted number of sequences = "+ nrOfGood);
			System.out.println("ErrorDistribution = "+ nrOfGood);
			for(int i = 0; i < errors;i++){
				System.out.println(i+"\t"+ errorDistribution[i]);
			}
			ER1.close();
			EW1.flush();
			EW1.close();
		}catch(Exception E){E.printStackTrace();}
	}

	private static void fastq2fasta(String inDir, String fileName1){
		try{
			ExtendedReader ER1 = new ExtendedReader(new FileReader (inDir+ "/" + fileName1));

			String outFile = fileName1.substring(0,fileName1.lastIndexOf('.'))+".fa";
			ExtendedWriter EW1 = new ExtendedWriter(new FileWriter (inDir+ "/"+ outFile));
			//ExtendedWriter EW2 = new ExtendedWriter(new FileWriter (outDir+ "/"+ fileName1+".perfect"));
			while(ER1.more()){
				FastQSequence A = new FastQSequence();
				boolean info1 = A.addInfo(ER1);

				A.printFasta(EW1);
			}
			ER1.close();
			EW1.flush();
			EW1.close();
		}catch(Exception E){E.printStackTrace();}
	}

	
	public static void fastqfixedLength(String inDir, String fileName1,int length){
		try{
			ExtendedReader ER1 = new ExtendedReader(new FileReader (inDir+ "/" + fileName1));

			String outFile = fileName1.substring(0,fileName1.lastIndexOf('.'))+"."+length+".fastq";
			ExtendedWriter EW1 = new ExtendedWriter(new FileWriter (inDir+ "/"+ outFile));
			//ExtendedWriter EW2 = new ExtendedWriter(new FileWriter (outDir+ "/"+ fileName1+".perfect"));
			while(ER1.more()){
				FastQSequence A = new FastQSequence();
				boolean info1 = A.addInfo(ER1,length);
				
				if(info1 && A.getIntSequence().length==length)
					A.printFastQ(EW1);
			}			
			ER1.close();
			EW1.flush();
			EW1.close();
		}catch(Exception E){E.printStackTrace();}
	}


	public static void filterPairedFastQSequences(String dir, String fileName1, String fileName2){
		filterPairedFastQSequences(dir+"/"+fileName1,dir+"/"+fileName2 );
	}
	private static void filterPairedFastQSequences(String fileName1, String fileName2){
		try{
			ExtendedReader ER1 = new ExtendedReader(new FileReader (fileName1));
			ExtendedReader ER2 = new ExtendedReader(new FileReader (fileName2));
			ExtendedWriter EW1 = new ExtendedWriter(new FileWriter (fileName1+".filtered.fastq"));
			ExtendedWriter EW2 = new ExtendedWriter(new FileWriter (fileName2+".filtered.fastq"));
			int count = 0;
			int nrOfGood = 0;
			int checked = 0;
			while(ER1.more()){
				count++;
				FastQSequence A = new FastQSequence();
				FastQSequence B = new FastQSequence();
				boolean info1 = A.addInfo(ER1);
				boolean info2 = B.addInfo(ER2);
				if(info1 && info2 && A.isPair(B)){
					nrOfGood++;
					A.printFastQ(EW1);
					B.printFastQ(EW2);
				}
				else{
					System.out.println("this sequence is bad");
					if(!A.QC()){
						System.out.println("seq 1 "+A.name );
					}
					if(!B.QC()){
						System.out.println("seq 2 "+B.name );
					}
					if(!A.isPair(B)){System.out.println(B.name +" is not the pair sequence of "+A.name );}
				}
				if(count > checked){
					System.out.println("nr Of Sequences checked " + nrOfGood);
					checked= checked + 1000000;
				}

			}

			System.out.println("Total number of sequences = "+ count);
			System.out.println("Accepted number of sequences = "+ nrOfGood);
			System.out.println("Percentage of good sequences = "+ (nrOfGood*100/count) + "% ");
			ER1.close();
			ER2.close();
			EW1.flush();
			EW1.close();
			EW2.flush();
			EW2.close();

		}catch(Exception E){E.printStackTrace();}
	}




	private void printSequences(ExtendedWriter EW){
		for(int i = 0; i < forward.size(); i++){
			forward.get(i).printFastQ(EW);
		}
	}

	private void printPairedSequences(ExtendedWriter EW, ExtendedWriter EW2){
		for(int i = 0; i < forward.size(); i++){
			forward.get(i).printFastQ(EW);
			reverse.get(i).printFastQ(EW2);
		}
	}

	private void printSequences(String outDir,String outFile){
		System.out.print("Printing sequences to file "+outDir+"/"+outFile+"  .......");
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile));
			printSequences(EW);
			EW.flush();
			EW.close();

		}catch(Exception E){E.printStackTrace();}
		System.out.println("finished");
	}

	private void printPairedSequences(String outDir,String outFile1, String outFile2 ){
		System.out.print("Printing sequences to file "+outDir+"/"+outFile1+"  .......");
		try{
			ExtendedWriter EW1 = new ExtendedWriter(new FileWriter(outDir+"/"+outFile1));
			ExtendedWriter EW2 = new ExtendedWriter(new FileWriter(outDir+"/"+outFile2));
			printPairedSequences(EW1, EW2);
			EW1.flush();
			EW1.close();
			EW2.flush();
			EW2.close();

		}catch(Exception E){E.printStackTrace();}
		System.out.println("finished");
	}

	public void setName(String name) {
		Name = name;
	}

	public String getName() {
		return Name;
	}

	public static void fixFastqFile(String inFile, String extender){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (inFile));
			ExtendedWriter EW = new ExtendedWriter(new FileWriter (inFile+".fixed"));

			while(ER.more()){

				if((char)ER.lookAhead() == '@'){
					String[] seqNames = ER.readLine().split(" ");
					String fixedSeqName = seqNames[0]+"/"+extender;
					for(int i = 1; i < seqNames.length;i++){
						 fixedSeqName += " "+seqNames[i];
					}
					EW.println(fixedSeqName);
				}
				else{
					EW.println(ER.readLine());
				}
			}
			EW.flush();
			EW.close();
			ER.close();
		}catch(Exception E){E.printStackTrace();}

	}



}




