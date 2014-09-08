package Sequence;

import extra.ExtendedReader;
import general.ExtendedWriter;
import general.Functions;
import general.GetFrequencies;

import java.io.File;
import java.io.FileWriter;
import java.util.Hashtable;

public class Generate {



	public static void run(Hashtable<String,String> T){

		String outFile, inFile, dir;
		int kmer, nrOfGeneratedSequences;
		outFile = inFile = dir = null;
		kmer = nrOfGeneratedSequences = 0;
		boolean allPresent = true;

		dir= Functions.getValue(T, "-d", ".");
		if(T.containsKey("-i"))
			inFile= Functions.getValue(T, "-i", ".");
		else{
			System.out.println("must contain inFile -i");
			allPresent = false;
		}

		kmer= Functions.getInt(T, "-k", 3);
		nrOfGeneratedSequences = Functions.getInt(T, "-n", 1);
		outFile= Functions.getValue(T, "-o", "generated");
		String[] freqArray = Functions.getValue(T, "-f","test").split(" ");
		int length = Integer.parseInt(Functions.getValue(T, "-l","0"));
	}

	private  static double[] getKmerFrequencies(String dir, String kmerFile){
		ExtendedReader ER = ExtendedReader.getFileReader(dir+"/"+kmerFile);
		ER.readLine();
		String[] freqArray = ER.readLine().split(" ");
		double[] freqs = null;
		if(freqArray.length > 3){
			freqs = new double[freqArray.length];
			for(int i = 0; i < freqArray.length; i++){
				freqs[i] = Double.parseDouble(freqArray[i]);
			}
		}
		return freqs;
	}


	public static void calculateKmerDistribution(String dir, String inFile, String kmerFile, int kmer){
		File infile = new File(dir+"/"+inFile);
		ExtendedWriter EW = ExtendedWriter.getFileWriter(dir+"/"+kmerFile);
		try{
			double[] freqs = GetFrequencies.getKmerFrequencies(GetFrequencies.parseBigFasta(infile),kmer);
			EW.println("Kmer freq used:");
			EW.print(freqs[0]);
			for(int i = 1; i < freqs.length; i++){
				EW.print(" "+freqs[i]);
			}
			EW.println();
		}catch(Exception E){E.printStackTrace();}
	}


	public static void generateSequenceFile(String dir, String seqFile, String kmerFile, String outFile, int kmers){
		FastaSequences A = new FastaSequences(dir,seqFile);
		double[] freqs = getKmerFrequencies(dir,kmerFile);
		try{
			ExtendedWriter EW = ExtendedWriter.getFileWriter(dir+"/"+outFile);
			for(int j = 0; j < A.size();j++){
				int length = A.get(j).Sequence.length;
				byte[] seq = GetFrequencies.generateSequence(length,kmers-1,freqs);
				char[] decodedSeq = GetFrequencies.decodeToChar(seq,"NACGTN");
				EW.println(A.getName()+"_generated");
				EW.println(decodedSeq);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}
	
	public static  void generateSequenceFile(String dir, String kmerFile, String outFile, int kmers, int start, double stop, int nrOfSequences){
		double[] freqs = getKmerFrequencies(dir,kmerFile);
		try{
			ExtendedWriter EW = ExtendedWriter.getFileWriter(dir+"/"+outFile);
			for(int j = 0; j < nrOfSequences;j++){
				int length = start+(int)(Math.random()*(stop-start));
				byte[] seq = GetFrequencies.generateSequence(length,kmers-1,freqs);
				char[] decodedSeq = GetFrequencies.decodeToChar(seq,"NACGTN");
				EW.println(">generated_"+j);
				EW.println(decodedSeq);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}
	
	


}
