package Sequence;

import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;

import alignment.Gene;

import general.ExtendedReader;
import general.ExtendedWriter;
import general.RNAfunctions;

public class Illumina  {
	
	String Name;

	protected FastaSequences sequences;
	
	
	public Illumina(){
		this.Name = "default";
		this.sequences = new FastaSequences();
	}
	public Illumina(String Name){
		this.Name = Name;
		this.sequences = new FastaSequences();
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		// TODO Auto-generated method stub
		
		Illumina A = new Illumina();
		A.readIlluminaFile(args[0], args[1], args[2],args[3],args[4]);
//		A.printSequences2Fasta(args[0],args[3]);
	}
	
	
	
	public void readIlluminaFile(String dir, String fileName, String fileName2, String RefSeq,String outName){
		
		try{
		ExtendedReader ER = new ExtendedReader(new FileReader(dir+"/"+fileName));
		ExtendedReader ER2 = new ExtendedReader(new FileReader(dir+"/"+fileName2));
		int nrOfReads,nrOfReadsGood,nrOfReadsAligned,nrOfReadsAligned2;
		nrOfReads = nrOfReadsGood = nrOfReadsAligned = nrOfReadsAligned2 = 0;
		int checkpoint = 100000;
		
		int[] RefSeqInt = RNAfunctions.RNAString2Int(RefSeq);
		int [][] distribution = new int[RefSeqInt.length][6];
		int[] nrOfErrors = new int[4]; 
		ExtendedWriter[] EWs = new ExtendedWriter[nrOfErrors.length];
		for(int i = 0; i < EWs.length;i++){
			EWs[i] = new ExtendedWriter(new FileWriter(dir+"/"+outName+"_"+i+".fa"));
		}
		ExtendedWriter summary = new ExtendedWriter(new FileWriter(dir+"/"+outName+".summary"));
		
		summary.println("Name of sequenceFile1                 : "+ dir+"/"+fileName);
		summary.println("Name of sequenceFile2                 : "+ dir+"/"+fileName2);
		for(int i = 0; i < EWs.length;i++){
			summary.println("Name of outfile with "+i+" mutations:" + dir+"/"+outName+"_"+i+".fa");
		}

		
		while (ER.more()){
			nrOfReads++;
			String[] strings =  ER.readLine().split("\t");
			String[] strings2 =  ER2.readLine().split("\t");
			if(strings[9].indexOf("B") == -1 && strings2[9].indexOf("B") == -1){
				nrOfReadsGood++;
				String Name = strings[4]+"_"+strings[5];
				String Name2 = strings2[4]+"_"+strings2[5];
				if(Name.compareTo(Name2) == 0){
					int[] seq = RNAfunctions.RNAString2Int(strings[8]);
					int[] seq2 = RNAfunctions.getReverseComplement(RNAfunctions.RNAString2Int(strings2[8]));
					//System.out.println(RNAfunctions.DNAInt2String(seq));
					//System.out.println(RNAfunctions.DNAInt2String(seq2));
					
					int[] joinedSeq = joinSequences(seq, seq2,30);
					if(joinedSeq != null){ 
						
						nrOfReadsAligned++;
						FastaSequence A  = new FastaSequence();
						A.Name = Name;
						
						//System.out.println(RNAfunctions.RNAInt2String(joinedSeq));
						int errors = SequenceHandling.matchSequence(joinedSeq,RefSeqInt, A ,nrOfErrors.length);
						if(errors == -1 )
							errors = SequenceHandling. matchSequence(RNAfunctions.getReverseComplement(joinedSeq),RefSeqInt,A ,nrOfErrors.length);
						
						if(errors != -1){
							nrOfErrors[errors]++;
							A.printFasta(EWs[errors]);
						}
						else
							nrOfErrors[nrOfErrors.length-1]++;
						
						if(errors == 1){
							int[] sense = A.Sequence;
							for(int i = 0; i < sense.length; i++){
								distribution[i][sense[i]]++;
							}
						}
						
					}
					
					
/*					joinedSeq = joinSequences2(seq, seq2,30);
					if(joinedSeq != null){ 
						nrOfReadsAligned2++;
						this.addSequence(Name,joinedSeq);
					}
*/					
				}
				
			}
			if(nrOfReads > checkpoint){
				System.out.println(nrOfReads+"\t"+nrOfReadsGood+"\t"+nrOfReadsAligned);
				checkpoint = checkpoint + 100000;
				for(int i = 0; i < EWs.length;i++){
					EWs[i].flush();
				}
			}
		}
		summary.println("Reads\tQ1\tQ2\t0\t1\t2\t3");
		summary.print(nrOfReads+"\t"+nrOfReadsGood+"\t"+nrOfReadsAligned);
		for(int i = 0; i < nrOfErrors.length;i++){
			summary.print("\t"+nrOfErrors[i]);
		}
		summary.flush();
		summary.close();
		
		for(int i = 0; i < EWs.length;i++){
			EWs[i].flush();
			EWs[i].close();
		}
		
		
		for(int i = 0; i < distribution.length;i++){
			System.out.print(RNAfunctions.RNAInt2String(RefSeqInt[i])+"\t");
			int total = 0;
			int mutants = 0 ;
			for(int j = 0; j < distribution[i].length;j++){
				System.out.print(distribution[i][j]+"\t");
			}
			System.out.println();
		}
		}
		catch(Exception E ){E.printStackTrace();}
	}

	private int[] joinSequences(int[] seq, int[] seq2, int length){
		boolean found = false;
		int position = 64;
		while(position < seq.length-length && found ==  false){
			int count = 0;
			boolean same = true;
			while(same && position+count < seq.length){
				if(seq[position+count] == seq2[count])
					count++;
				else
					same = false;
			}
			if(same){
				int[] joinedSequences = new int[position+seq2.length];
				for(int i = 0; i < position;i++){
					joinedSequences[i] = seq[i];
				}	
				for(int i = 0; i < seq2.length;i++){
					joinedSequences[i+position] = seq2[i];
				}	
				return joinedSequences;
			}
			position++;
		}
		return null;
	}

	
	
	
}


