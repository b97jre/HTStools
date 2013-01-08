package Sequence;

import general.ExtendedWriter;
import general.Functions;
import general.GetFrequencies;

import java.io.File;
import java.io.FileWriter;
import java.util.Hashtable;

public class Generate {

	String dir;
	String inFile;
	String outFile;
	int kmer;
	int nrOfGeneratedSequences;


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

		double[] freqs = null;
		if(freqArray.length > 3){
			freqs = new double[freqArray.length];
			for(int i = 0; i < freqArray.length; i++){
				freqs[i] = Double.parseDouble(freqArray[i]);
			}
		}



		FastaSequences A = new FastaSequences(dir,inFile);
		File infile = new File(dir+"/"+inFile);
		try{
			if(freqs == null)
				freqs = GetFrequencies.getKmerFrequencies(GetFrequencies.parseBigFasta(infile),kmer);
			System.out.println("Frequencies used:");
			System.out.print(freqs[0]);
			for(int i = 1; i < freqs.length; i++){
				System.out.print(" "+freqs[i]);
			}
			System.out.println();

			if(length == 0){
				for(int i = 0; i < nrOfGeneratedSequences;i++){
					try{
						ExtendedWriter EW = new ExtendedWriter(new FileWriter(dir+"/"+outFile+"_"+i+".fa"));
						for(int j = 0; j < A.size();j++){
							length = A.get(j).Sequence.length;
							byte[] seq = GetFrequencies.generateSequence(length,kmer-1,freqs);
							char[] decodedSeq = GetFrequencies.decodeToChar(seq,"NACGTN");
							EW.println(">generated_"+j);
							EW.println(decodedSeq);
						}
						EW.flush();
						EW.close();
					}catch(Exception E){E.printStackTrace();}
				}
			}
			else{
				try{
					ExtendedWriter EW = new ExtendedWriter(new FileWriter(dir+"/"+outFile+".fa"));
					for(int i = 0; i < nrOfGeneratedSequences;i++){
						byte[] seq = GetFrequencies.generateSequence(length,kmer-1,freqs);
						char[] decodedSeq = GetFrequencies.decodeToChar(seq,"NACGTN");
						EW.println(">generated_"+i);
						EW.println(decodedSeq);
					}
					EW.flush();
					EW.close();
				}catch(Exception E){E.printStackTrace();}
			}


		}catch(Exception E){E.printStackTrace();}
	}



}
