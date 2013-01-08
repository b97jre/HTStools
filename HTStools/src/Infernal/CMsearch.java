package Infernal;

import java.io.FileReader;
import java.util.ArrayList;
import java.util.Hashtable;

import Sequence.CfastaSequences;
import Sequence.FastaSequence;
import Sequence.FastaSequences;
import Sequence.Generate;
import Sequence.SequenceHandling;

import general.ExtendedReader;
import general.Functions;
import general.IOTools;
import general.RNAfunctions;

public class CMsearch extends ArrayList<CMresult>{

	
	
	public static void main(String[] args) {
		
		// TODO Auto-generated method stub
		
		CMsearch A = new CMsearch();
		A.readCMsearchResultTabFile(args[0], args[1]);
		
		FastaSequences FS = new FastaSequences(args[2],args[3]);
		
		System.out.println("Extracting sequences");
		System.out.println("printing sequences");
		
		FastaSequences FS3 = A.ExtractSequences(FS,150,0);
		
		System.out.println("printing sequences");
		FS3.writeFastaFile(args[2], args[4]+"_"+150+"nt_upstream.fa");
		
	}
	
	public CMsearch(){
		
	}
	public static void run(Hashtable<String,String> T){
		String cmFile  = Functions.getValue(T, "-i", "i");
		String refFile  = Functions.getValue(T, "-r", "r");
		String motifFile = Functions.getValue(T, "-m", "r");
		String dremeFile = Functions.getValue(T, "-d", "r");
		String pDir  = Functions.getValue(T, "-pDir", ".");
		int us = Integer.parseInt(Functions.getValue(T, "-us", "100")); 
		int ds = Integer.parseInt(Functions.getValue(T, "-ds", "100")); 
		
		double cutoff = Double.parseDouble(Functions.getValue(T, "-cutoff", "15.0")); 
	
		
		CMsearch A = new CMsearch();
		A.readCMsearchResultTabFile(pDir, cmFile);
		FastaSequences FS = new FastaSequences(pDir,refFile);

		
		if(T.containsKey("-extract")){
			
			System.out.println("Extracting sequences");
			
			FastaSequences FS2 = A.ExtractSequences(FS,cutoff);
			System.out.println("printing sequences");
			FS2.writeFastaFile(pDir, cmFile+"_"+cutoff+"_sequence.fa");
	//		FastaSequences FS3 = A.ExtractSequences(FS,us,ds);
			System.out.println("printing sequences");
	//		FS3.writeFastaFile(pDir, cmFile+"_"+us+"nt_upstream"+"_"+ds+"nt_upstream.fa");
			FastaSequences FS4 = A.ExtractUpstreamSequences(FS,us,cutoff);
			System.out.println("printing sequences");
			FS4.writeFastaFileDNA(pDir, cmFile+"_"+cutoff+"_"+us+"_nt_upstream.fa");
		}
		if(T.containsKey("-motif")){
			if(IOTools.fileExists(pDir, dremeFile))
				A.readDreme(pDir,dremeFile);
			A.readFimo(pDir, motifFile);
			A.printMotifs(refFile,FS,us);
		}
	}


	
	
	
	
	private void readDreme(String dir, String fileName){
		try{
		ExtendedReader ER = new ExtendedReader(new FileReader(dir+"/"+fileName));
		while (ER.more()){
			String info  = ER.readLine();
			if(info.indexOf("MOTIF") == 0){
				String[] info2 = info.split(" ");
				for(int i = 0; i < this.size();i++){
					Motif newMotif = new Motif();
					newMotif.setMotif(info2[1]);
					this.get(i).addDremeMotifs(newMotif);
				}
			}
		}
		
		}
		catch(Exception E){E.printStackTrace();}
	}

	private void readFimo(String dir, String fileName){
		try{
		ExtendedReader ER = new ExtendedReader(new FileReader(dir+"/"+fileName));
		while(ER.more()){
			while(ER.lookAhead() == '#'){ER.readLine();}
			readFimoGFF3Result(ER.readLine());
		}
		}
		catch(Exception E){E.printStackTrace();}
	}

	
	
	private void readCMsearchResultTabFile(String dir, String fileName){
		try{
		ExtendedReader ER = new ExtendedReader(new FileReader(dir+"/"+fileName));
		while(ER.more()){
			while(ER.lookAhead() == '#'){ER.readLine();}
			readTabResult(ER);
			while(ER.lookAhead() == '#'){ER.readLine();}
		}
		}
		catch(Exception E){E.printStackTrace();}
	}
	
	private void readTabResult(ExtendedReader ER){
		String info2 = ER.readLine();
		//System.out.println(info2);
		info2.trim();
		while(info2.indexOf("  ") != -1)
			info2 = info2.replaceAll("  ", " ");
		while(info2.indexOf(" ") == 0)
			info2 = info2.replaceFirst(" ", "");

		String[] info = info2.split(" ");
		String fastaFile = info[1];
		int targetStart = Integer.parseInt(info[2]);
		int targetStop = Integer.parseInt(info[3]);
		int queryStart = Integer.parseInt(info[4]);
		int queryStop = Integer.parseInt(info[5]);
		double bitScore = Double.parseDouble(info[6]);
		
		this.add(new CMresult(fastaFile, targetStart, targetStop, queryStart, queryStop, bitScore));
		
	}

	
	private void readFimoGFF3Result(String info2){
		if(info2.indexOf("#") != 0){
			boolean found = false;
			int count = 0;
			while(!found && count < this.size()){
				found = this.get(count).addMotif(info2);
				count++;	
			}
			if(!found)System.out.println("not found "+info2);
		}
		//MCCATAA DDB0232428_1892631_1892734_upstream_100_downstream_100_17.43    36      42      13.3007 6.51e-05                CCCATAA
		
	}
	
	private void printMotifs(String referenceFile, FastaSequences targetSequences, int upstream){
		System.out.println("Genome\tContig\tStart\tStop\tCMBitScore\tTotalScore\tDistance\tMotif1\tStart\tStop\tStrand\tbitScore\tSequence\tMotif2\tStart\tStop\tStrand\tbitScore\tSequence");
		for(int i = 0; i < this.size(); i++){
			//this.get(i).printMotifs();
			CMresult temp = this.get(i);
			int[] Sequence = targetSequences.getSequence(temp.fastaFile,temp.targetStart,temp.targetStop, upstream, 0);
			int[] upStream = Functions.getSubarray(Sequence, 0, upstream);
			String Sequence2 = RNAfunctions.DNAInt2String(upStream);
			temp.printSpecificMotifs(referenceFile,100,60,Sequence2);
		}
	}
	
	
	
	
	
	private FastaSequences ExtractSequences(FastaSequences targetSequences, double cutoff){
		
		FastaSequences RNAgenes = new FastaSequences();
		for(int i = 0; i < this.size();i++){
			if(this.get(i).bitScore>cutoff){
				CMresult temp = this.get(i);
				int[] Sequence = targetSequences.getSequence(temp.fastaFile,temp.targetStart,temp.targetStop);
			
				if(Sequence != null){
					RNAgenes.add(new FastaSequence(temp.fastaFile+"_"+temp.targetStart+"_"+temp.targetStop+"_"+temp.bitScore,Sequence));
				}
			}
			
		}
		return RNAgenes;
	}
	
	
	private FastaSequences ExtractSequences(FastaSequences targetSequences, int upstream, int downstream){
		FastaSequences RNAgenes = new FastaSequences();
		for(int i = 0; i < this.size();i++){
			CMresult temp = this.get(i);
			int[] Sequence = targetSequences.getSequence(temp.fastaFile,temp.targetStart,temp.targetStop, upstream, downstream);
			
			if(Sequence != null){
				RNAgenes.add(new FastaSequence(temp.fastaFile+"_"+temp.targetStart+"_"+temp.targetStop+"_upstream_"+upstream
						+"_downstream_"+downstream+"_"+temp.bitScore,Sequence));
			}
			
		}
		return RNAgenes;
	}

	private FastaSequences ExtractUpstreamSequences(FastaSequences targetSequences, int upstream, double cutoff){
		FastaSequences RNAgenes = new FastaSequences();
		for(int i = 0; i < this.size();i++){
			if(this.get(i).bitScore>cutoff){
			CMresult temp = this.get(i);
			int[] Sequence = targetSequences.getSequence(temp.fastaFile,temp.targetStart,temp.targetStop, upstream, 0);
			int[] upStream = Functions.getSubarray(Sequence, 0, upstream);
			if(Sequence != null){
				RNAgenes.add(new FastaSequence(temp.fastaFile+";"+temp.targetStart+";"+temp.targetStop+";"+temp.bitScore,upStream));
			}
			}
			
		}
		return RNAgenes;
	}

	
	
	

}

