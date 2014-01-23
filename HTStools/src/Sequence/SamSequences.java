package Sequence;

import general.Functions;
import general.IOTools;
import general.RNAfunctions;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;


import variousTools.Trinity;
import general.ExtendedReader;
import general.ExtendedWriter;




public class SamSequences extends Hashtable <String,FastaSequences> implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//private String GenomeName;
	private String Name;


	private int length;
	private int missmatches;
	private int inside;

	private Hashtable <String,Pair> pairs;
	boolean trinity;
	boolean oases;

	public static void main(String []args){
		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}
		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);
		SamSequences.run(T);
	}

	SamSequences(){
		length = 10;
		missmatches = 3;
		inside = 500;
		this.pairs = new Hashtable <String,Pair>();
		trinity = false;
		oases = false;
	}




	public static void run(Hashtable<String,String> T){
		String dir = Functions.getValue(T, "-d", IOTools.getCurrentPath());
		if(T.containsKey("-ref") ){
			String inFile = Functions.getValue(T, "-ref");


			if(T.containsKey("-extract") && T.containsKey("-i")){
				System.out.println("running samtools parsing....");
				String samFile = Functions.getValue(T, "-i");

				int cutoff = Functions.getInt(T, "-C", 10);

				Hashtable<String,String> flags = new Hashtable<String,String>();
				SamSequences A = new SamSequences();
				if(T.containsKey("-trinity")) A.trinity= true;
				if(T.containsKey("-oases") )A.oases=true;
				flags.put("101","101");
				flags.put("165","165");
				System.out.println();

				A.parseFile(dir+"/"+samFile,flags);
				A.removeLow(cutoff);
				System.out.println("Printing 5' ends");
				A.printFasta(dir+"/FivePrime");

				A = new SamSequences();
				flags.clear();
				flags.put("69","69");
				flags.put("133","133");
				Hashtable<String,Integer> sizes = A.getSizes(inFile,dir);
				A.parseReverseFile(dir+"/"+samFile,flags,sizes);
				A.removeLow(cutoff);
				System.out.println("Printing 3' ends");
				A.printFasta(dir+"/ThreePrime");

				System.out.println("Finnished running samtools parsing....");

			}

			else if(T.containsKey("-merge") ){
				System.out.println("merging new five and three prime ends....");
				SamSequences A = new SamSequences();
				if(T.containsKey("-trinity")) A.trinity= true;
				if(T.containsKey("-oases") )A.oases=true;
				System.out.println("Reading fasta files....");
				A.parseAllFasta(inFile,dir);
				if(T.containsKey("-samEnds") ){
					System.out.println("Merging ends....");

					if(T.containsKey("-R") ){
						System.out.println("Checking false discovery rate");
						A.mergeTrinityReadsFalse(dir);
					}
					else
						A.mergeTrinityReads(dir);
					A.printAllFasta(dir+"/AllSequences_ends.fa");

				} 
				if(T.containsKey("-samPairs") ){
					System.out.println("merging new five and three prime ends....");
					String samFile = Functions.getValue(T,"-i");
					Hashtable<String,Integer> sizes = A.getSizes(inFile,dir);
					Hashtable<String,String> flags = new Hashtable<String,String>();
					flags.clear();

					if(!T.containsKey("-direction")){
						flags.put("161","161");
						flags.put("97","97");
						flags.put("65","65");
						flags.put("129","129");
						flags.put("113","113");
						flags.put("177","177");
						flags.put("81","81");
						flags.put("145","145");
					}
					else{
						flags.put("161","161");
						flags.put("97","97");
					}


					System.out.println("Parsing pairs");
					A.parsePairs(samFile,flags,sizes);
					System.out.println("Reading fasta files....");
					A.printAllPairs(samFile+".merged");
					System.out.println("Merging ends....");
					A.mergePairs();
					System.out.println("Merging ends....");
					A.printAllFasta(dir+"/AllSequences_merged.fa");

				}

			}

		}
	}





	private void mergePairs() {
		try{
			for (Enumeration<String> e = pairs.keys(); e.hasMoreElements();){
				String hitName = e.nextElement();
				//System.out.println(hitName);
				Pair A = pairs.get(hitName);
				if(A.nrOfHits > 1){
					FastaSequences B = this.get(A.refName1);
					FastaSequences C = this.get(A.refName2);
					if(B != null && C != null){
						if(A.end1.compareTo("3end") == 0 && A.end2.compareTo("5end") == 0){
							//System.out.println("--->|--->  3|5");
							mergeSequences(B,C,false,false,A.nrOfHits);
						}
						else if(A.end1.compareTo("3end") == 0 && A.end2.compareTo("3end") == 0){
							//System.out.println("--->|<---	 3|3");
							mergeSequences(B,C,false,true,A.nrOfHits);
						}
						else if( A.end1.compareTo("5end") == 0 && A.end2.compareTo("5end") == 0){
							//System.out.println("<---|--->  5|5");
							mergeSequences(B,C,true,false,A.nrOfHits);
						}
						else if ( A.end1.compareTo("5end") == 0 && A.end2.compareTo("3end") == 0){
							//System.out.println("<---|<---  5|3"); 
							mergeSequences(B,C,true,true,A.nrOfHits);
						}
					}
				}

			}
		}catch(Exception E){E.printStackTrace();}
	}


	public void mergeSequencesORFs(FastaSequences temp, FastaSequences FE, boolean reverse, boolean reverse2, int nrOfHits){


		FastaSequences temp2 =  new FastaSequences();
		temp2.Name = temp.Name;
		FastaSequences FE2 =  new FastaSequences();
		FE2.Name = FE.Name;



		//int length = this.length;
		int length = 5;

		if(nrOfHits < 4)
			length = 10;
		else if(nrOfHits < 10)
			length = 5;
		int tempSize = temp.size();
		int FEsize = FE.size();
		//
		//		System.out.println();
		//		System.out.println("Before");
		//		System.out.println(FastaName+"\t"+this.get(FastaName).size()+"\t"+tempSize);
		//		System.out.println(OtherFastaName+"\t"+this.get(OtherFastaName).size()+"\t"+FEsize);
		//		System.out.println();
		//		
		//		
		for(int i = 0; i < tempSize;i++){
			FastaSequence tempSeq = temp.get(i);
			int [] seq2 = tempSeq.getSequence();
			if(reverse)seq2 = RNAfunctions.getReverseComplement(seq2);
			for(int j = 0; j < FEsize;j++){
				FastaSequence tempSeqFP = FE.get(j);
				int [] seqFP = tempSeqFP.getSequence();
				if(reverse2)seqFP = RNAfunctions.getReverseComplement(seqFP);
				int[] newSeq = RNAfunctions.merge3endIfORF(seqFP, seq2, length, missmatches, inside);
				int[] newSeqRev = RNAfunctions.getReverseComplement(newSeq);
				if(newSeq.length > seq2.length){
					String FastaName = tempSeq.getName().substring(1);
					String OtherFastaName = tempSeqFP.getName().substring(1);
					FastaSequence nSeq = null;
					if(reverse){
						
						nSeq = new FastaSequence(FastaName+"_Merged_"+OtherFastaName,newSeqRev); 
					}else{
						
						nSeq = new FastaSequence(FastaName+"_Merged_"+OtherFastaName,newSeq); 
					}
					temp2.add(nSeq);
					FastaSequence nSeq2 = null;
					if(reverse){
						nSeq2 = new FastaSequence(OtherFastaName+"_Merged_"+FastaName,newSeqRev); 
					}else{
						nSeq2 = new FastaSequence(OtherFastaName+"_Merged_"+FastaName,newSeq); 
					}
					FE2.add(nSeq2);
					System.out.println(FastaName+"_"+OtherFastaName+"_"+reverse+"_"+reverse2+"_"+nrOfHits);
				}
			}

		}
		for(int i = 0; i < tempSize;i++){temp2.add(temp.get(i));}
		for(int j = 0; j < FEsize;j++){FE2.add(FE.get(j));}
		//		if(temp2.size() > tempSize){
		//			for(int i = 0 ; i < temp2.size();i++){
		//				System.out.println(temp2.get(i).getName()+"\t"+temp2.get(i).getLength());
		//			}
		//		}
		//		
		//		if(FE2.size() > FEsize){
		//			for(int i = 0 ; i < FE2.size();i++){
		//				System.out.println(FE2.get(i).getName()+"\t"+FE2.get(i).getLength());
		//			}
		//		}
		//		
		String FastaName = temp.getName();
		String OtherFastaName = FE.getName();
		this.put(FastaName, temp2);
		this.put(OtherFastaName, FE2);
		//		
		//		System.out.println("After");
		//		System.out.println(FastaName+"\t"+this.get(FastaName).size());
		//		System.out.println(OtherFastaName+"\t"+this.get(OtherFastaName).size());
		//		System.out.println();
	}


	public void mergeSequences(FastaSequences temp, FastaSequences FE, boolean reverse, boolean reverse2, int nrOfHits){


		//int length = this.length;
		int length = 5;

		if(nrOfHits < 4)
			length = 25;
		else if(nrOfHits < 10)
			length = 10;
		int tempSize = temp.size();
		int FEsize = FE.size();
		//
		//		System.out.println();
		//		System.out.println("Before");
		//		System.out.println(FastaName+"\t"+this.get(FastaName).size()+"\t"+tempSize);
		//		System.out.println(OtherFastaName+"\t"+this.get(OtherFastaName).size()+"\t"+FEsize);
		//		System.out.println();
		//		
		//		
		for(int i = 0; i < tempSize;i++){
			FastaSequence tempSeq = temp.get(i);
			int [] seq2 = tempSeq.getSequence();
			if(reverse)seq2 = RNAfunctions.getReverseComplement(seq2);
			for(int j = 0; j < FEsize;j++){
				FastaSequence tempSeqFP = FE.get(j);
				int [] seqFP = tempSeqFP.getSequence();
				if(reverse2)seqFP = RNAfunctions.getReverseComplement(seqFP);
				int[] newSeq = RNAfunctions.merge3end(seqFP, seq2, length, missmatches, inside);
				int[] newSeqRev = RNAfunctions.getReverseComplement(newSeq);
				String FastaName = tempSeq.getName();
				String OtherFastaName = tempSeqFP.getName();

				
				if(newSeq.length > seq2.length){
					FastaSequences temp2 =  new FastaSequences();
					temp2.Name = FastaName+"_Merged_"+OtherFastaName;
					
					FastaSequence nSeq = null;
					if(reverse){
						nSeq = new FastaSequence(FastaName+"_Merged_"+OtherFastaName,newSeqRev); 
					}else{
						nSeq = new FastaSequence(FastaName+"_Merged_"+OtherFastaName,newSeq); 
					}
					temp2.add(nSeq);
					this.put(FastaName, temp2);
				}
			}

		}
	}


	
	
	public void removeRedundantSequences(){
		
	}
	
/*	public void mergeSequencesAgain(FastaSequences temp, FastaSequences FE){

		int tempSize = temp.size();
		int FEsize = FE.size();
		//
		//		System.out.println();
		//		System.out.println("Before");
		//		System.out.println(FastaName+"\t"+this.get(FastaName).size()+"\t"+tempSize);
		//		System.out.println(OtherFastaName+"\t"+this.get(OtherFastaName).size()+"\t"+FEsize);
		//		System.out.println();
		//		
		//		
		for(int i = 0; i < tempSize;i++){
			FastaSequence tempSeq = temp.get(i);
			int [] seq2 = tempSeq.getSequence();
			for(int j = 0; j < FEsize;j++){
				boolean found  = false;
				FastaSequence tempSeqFP = FE.get(j);
				int [] seqFP = tempSeqFP.getSequence();
				int[] newSeq = RNAfunctions.merge3end(seqFP, seq2, 100, 0, seq2.length-100);
				int[] newSeq1 = RNAfunctions.merge3end(RNAfunctions.getReverseComplement(seqFP), seq2, 100, 0, seq2.length-100);
				int[] newSeq2 = RNAfunctions.merge3end(seqFP, RNAfunctions.getReverseComplement(seq2), 100, 0, seq2.length-100);
				int[] newSeq3 = RNAfunctions.merge3end(RNAfunctions.getReverseComplement(seqFP), RNAfunctions.getReverseComplement(seq2), 100, 0, seq2.length-100);
					
					
				if(newSeq.length > seq2.length){
					FastaSequences temp2 =  new FastaSequences();
					temp2.Name = FastaName+"_Merged_"+OtherFastaName;
					
					FastaSequence nSeq = null;
					if(reverse){
						nSeq = new FastaSequence(FastaName+"_Merged_"+OtherFastaName,newSeqRev); 
					}else{
						nSeq = new FastaSequence(FastaName+"_Merged_"+OtherFastaName,newSeq); 
					}
					temp2.add(nSeq);
					this.put(FastaName, temp2);
				}
			}

		}
	}
	
	*/
	public void parseFile(String inFile, Hashtable<String,String> flags){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(inFile));
			int count = 1000;
			int lines = 1000000;
			int line = 0;
			int good = 0;
			while(ER.more()){
				line++;
				if((char)ER.lookAhead() == '@'){
					ER.skipLine();
				}
				else{
					good += addSamLineStart(ER, flags);
				}
				if(good > count){
					System.out.println("nr Of Sequences Read" + count);
					count= count + 1000;
				}
				if(line > lines){
					System.out.println("nr Of lines Read" + lines);
					lines= lines + 1000000;
				}
			}
			System.out.println("Finished total nr of sequences = " + this.size());
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}

	public void parsePairs(String inFile, Hashtable<String,String> flags,Hashtable<String,Integer> sizes){
		try{
			System.out.println("inFile");
			ExtendedReader ER = new ExtendedReader(new FileReader(inFile));
			//			int count = 1000;
			int line = 0;
			int lines = 1000000;
			int count  = 1000;
			while(ER.more() ){
				if((char)ER.lookAhead() == '@'){
					ER.skipLine();
				}
				else{
					addSamPair(ER,flags,sizes);
					line++;
				}
				if(pairs.size() > count){
					System.out.println("nr Of pairs found" + count);
					count= count + 1000;
				}
				if(line > lines){
					System.out.println("nr Of lines Read" + lines);
					lines= lines + 1000000;
				}
			}
			System.out.println("Finished total nr of sequences = " + pairs.size());
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}



	public void parseReverseFile(String inFile, Hashtable<String,String> flags,Hashtable<String,Integer> contigs){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(inFile));
			int count = 1000000;
			int good = 0;
			int pointer = 0;
			System.out.println("Parsing for three prime end reads.... ");
			while(ER.more()){
				if((char)ER.lookAhead() == '@'){
					ER.skipLine();
				}
				else{
					pointer++;
					good += addSamLineEnd(ER, flags,contigs);
				}
				if(pointer > count){
					System.out.println("Nr of lines read:" + count);
					System.out.println("Nr of sequences found" + this.size());
					System.out.println("Nr Of sequences with more than 10 reads" + good);
					System.out.println();

					count= count + 1000000;
				}
			}
			System.out.println("Finished total nr of sequences = " + this.size());
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}




	public void parseFasta(String inFile){
		try{
			FastaSequences.oneLine(inFile);
			ExtendedReader ER = new ExtendedReader(new FileReader(inFile+".oneLine"));
			int count = 1000;
			while(ER.more()){
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					addFasta(ER);
				}
				if(this.size() > count){
					System.out.println("nr Of Sequences Read" + count);
					count= count + 1000;
				}
			}
			System.out.println("Finished total nr of sequences = " + this.size());
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}

	public void parseAllFasta(String inFile, String WD){
		try{
			FastaSequences.oneLineFirst(inFile,WD+"/AllSequences.fa");
			ExtendedReader ER = new ExtendedReader(new FileReader(WD+"/AllSequences.fa"));
			int count = 1000;
			while(ER.more()){
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					addAllFasta(ER);
				}
				if(this.size() > count){
					System.out.println("nr Of Sequences Read" + count);
					count= count + 1000;
				}
			}
			System.out.println("Finished total nr of sequences = " + this.size());
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}


	public Hashtable<String,Integer> getSizes(String inFile, String WD){
		Hashtable<String,Integer> HT = new Hashtable<String,Integer>();
		try{
			FastaSequences.oneLine(inFile,WD+"/AllSequences.fa");
			System.out.println("Reading sizes....");
			ExtendedReader ER = new ExtendedReader(new FileReader(WD+"/AllSequences.fa"));
			int count = 1000;
			while(ER.more()){
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					getSizes(ER,HT);
				}
			}
			System.out.println("finished. Total nr of sequences = " + HT.size());
			ER.close();
		}catch(Exception E){E.printStackTrace();}
		return HT;
	}


	private void addSamLine(ExtendedReader EW){
		String[] info = EW.readLine().split("\t");
		//trinitySpecific
		String[] hitInfo = info[2].split("_");
		String hitName = hitInfo[0]+"_"+hitInfo[1]+"_"+hitInfo[2];
		FastaSequence A = new FastaSequence(info[0], info[9]);
		if(this.containsKey(hitName)){
			this.get(hitName).add(A);
		}
		else{
			FastaSequences FS = new FastaSequences();
			FS.add(A);
			this.put(hitName, FS);
		}
	}

	private int addSamLineStart(ExtendedReader EW, Hashtable<String,String> flags){
		String[] info = EW.readLine().split("\t");
		if(flags.containsKey(info[1])){
			//trinitySpecific
			String hitName = info[2];
			if(trinity){
				String[] hitInfo = info[2].split("_");
				hitName = hitInfo[0]+"_"+hitInfo[1]+"_"+hitInfo[2];
			}
			else if(oases){
				String[] hitInfo = info[2].split("/");
				hitName = hitInfo[0];
			}
			int location = Integer.parseInt(info[3]);
			if(location < inside ){
				FastaSequence A = new FastaSequence(info[0], info[9]);
				if(this.containsKey(hitName)){
					this.get(hitName).add(A);
					if(this.get(hitName).size() == 10) return 1;	
				}
				else{
					FastaSequences FS = new FastaSequences();
					FS.add(A);
					this.put(hitName, FS);
				}
			}
		}
		return 0;
	}



	private int addSamLineEnd(ExtendedReader EW, Hashtable<String,String> flags,Hashtable<String,Integer> sizes){
		String[] info = EW.readLine().split("\t");
		if(flags.containsKey(info[1])){
			//trinitySpecific
			String[] hitInfo = info[2].split("_");
			String hitName = hitInfo[0]+"_"+hitInfo[1]+"_"+hitInfo[2];
			int location = Integer.parseInt(info[3]);
			//			System.out.println(location+ " "+ hitName);
			//			System.out.println(sizes.get(hitName).intValue());
			if(sizes.containsKey(hitName) && location+inside > sizes.get(hitName).intValue()){
				FastaSequence A = new FastaSequence(info[0], RNAfunctions.getReverseComplement(info[9]));
				if(this.containsKey(hitName)){
					this.get(hitName).add(A);
					if(this.get(hitName).size() == 10) return 1;
				}
				else{
					FastaSequences FS = new FastaSequences();
					FS.add(A);
					this.put(hitName, FS);
				}
			}
		}
		return 0;
	}


	private void addSamPair(ExtendedReader EW, Hashtable<String,String> flags,Hashtable<String,Integer> sizes){
		String[] info = EW.readLine().split("\t");
		if(flags.containsKey(info[1])){
			int leftPos = Integer.parseInt(info[3]);
			int rightPos = Integer.parseInt(info[7]);
			//trinitySpecific
			String[] hitInfo = info[2].split("_");
			String hitName1 = hitInfo[0]+"_"+hitInfo[1]+"_"+hitInfo[2];
			String contig = hitInfo[0];
			hitInfo = info[6].split("_");
			if(hitInfo.length > 2){
				String hitName2 = hitInfo[0]+"_"+hitInfo[1]+"_"+hitInfo[2];
				int tag = Integer.parseInt(info[1]);
				if(hitName1.compareTo(hitName2) != 0 ){
					if((tag == 97) &&  leftPos+inside > sizes.get(hitName1).intValue() && rightPos < inside)
						addPair(hitName1,"3end",hitName2,"5end",tag);
					else if((tag == 65) &&  leftPos+inside > sizes.get(hitName1).intValue() && rightPos+inside > sizes.get(hitName1).intValue())
						addPair(hitName1,"3end",hitName2,"3end",tag);
					else if((tag == 113) && leftPos < inside && rightPos < inside)
						addPair(hitName1,"5end",hitName2,"5end",tag);
					else if ((tag == 81) && leftPos < inside && rightPos+inside > sizes.get(hitName1).intValue())
						addPair(hitName1,"5end",hitName2,"3end",tag);
				}
			}
		}



	}

	private void addPair(String hitName1, String end1, String hitName2, String end2, int code){
		if(this.pairs.containsKey(hitName1+"_"+hitName2+"_"+end1+"_"+end2))
			this.pairs.get(hitName1+"_"+hitName2+"_"+end1+"_"+end2).addHit();
		else{
			Pair newPair = new Pair(this.get(hitName1),end1,this.get(hitName2),end2,code);
			//System.out.println(hitName1+"_"+hitName2+"_"+end1+"_"+end2);
			this.pairs.put(hitName1+"_"+hitName2+"_"+end1+"_"+end2, newPair);
		}
	}



	private void mergeSequences(SamSequences other){
		for (Enumeration<String> e = other.keys(); e.hasMoreElements();){
			String hitName = e.nextElement();
			FastaSequences A = other.get(hitName);
			if(this.containsKey(hitName)){
				this.get(hitName).addAll(A);
			}
			else{
				this.put(hitName, A);
			}
		}
	}

	private void printAllFasta(String fileName){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(fileName));
			for (Enumeration<String> e = this.keys(); e.hasMoreElements();){
				String hitName = e.nextElement();
				FastaSequences A = this.get(hitName);
				if(A.size()>1){
					for(int i = 0; i < A.size();i++){
						A.get(i).Name =A.get(i).Name+"0"+i; 
						System.out.println(A.get(i).Name);
					}
				}
				A.printFasta(EW);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}

	private void printAllPairs(String fileName){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(fileName));
			for (Enumeration<String> e = pairs.keys(); e.hasMoreElements();){
				String hitName = e.nextElement();
				Pair A = pairs.get(hitName);
				A.printPair(EW);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}

	private void addFasta(ExtendedReader EW){
		//trinitySpecific
		String[] hitInfo = EW.readLine().split("_");
		String hitName = hitInfo[0].substring(1)+"_"+hitInfo[1]+"_"+hitInfo[2];
		String Sequence = EW.readLine();
		FastaSequence A = new FastaSequence(hitName, Sequence);
		if(this.containsKey(hitName)){
			this.get(hitName).add(A);
		}
	}

	private void addAllFasta(ExtendedReader EW){
		//trinitySpecific
		while(EW.more()){
			String infoLine = EW.readLine();
			//System.out.println(infoLine);
			String[] info = infoLine.split("\\ ");

			String hitName = info[0];
			if(trinity){
				String[] hitInfo = info[0].split("_");
				hitName = hitInfo[0]+"_"+hitInfo[1]+"_"+hitInfo[2];
			}
			else if(oases){
				String[] hitInfo = info[0].split("/");
				hitName = hitInfo[0];
			}

			String Sequence = EW.readLine();
			FastaSequence A = new FastaSequence(hitName, Sequence);
			FastaSequences FS = new FastaSequences();
			FS.add(A);
			hitName = hitName.substring(1);
			FS.Name = hitName;
			//System.out.println(hitName);
			this.put(hitName, FS);
		}
	}

	private void getSizes(ExtendedReader EW, Hashtable<String,Integer> HT){
		//trinitySpecific
		while(EW.more()){
			String infoLine = EW.readLine();
			//System.out.println(infoLine);
			String[] info = infoLine.split("\\ ");
			String[] hitInfo = info[0].split("_");
			//			int length = Integer.parseInt(info[1].split("=")[1]);

			String hitName = hitInfo[0].substring(1)+"_"+hitInfo[1]+"_"+hitInfo[2];
			String Sequence = EW.readLine();
			//			if(Sequence.length() == length)
			HT.put(hitName, new Integer(Sequence.length()));
			//			else
			//				System.out.println(hitName+"\t"+length+" "+Sequence.length());
		}
	}


	void mergeTrinityReads(String dir){
		try{
			System.out.print("Reading directory for 5' ends....");
			ArrayList<String> dirs = Functions.getStrings(dir+"/Fend.txt");
			System.out.println("Finished");
			System.out.println("Number of 5' ends to be tested :"+dirs.size());
			System.out.println("Merging 5' ends....");
			System.out.print("Percent finished..0%..");
			int percent = 10; 
			int point = 1;
			for(int i = 0; i < dirs.size();i++){
				if(IOTools.fileExists(dir+"/Fend/"+dirs.get(i)+"/Trinity.fasta")){
					//					System.out.println("checking " + dirs.get(i) );
					FastaSequences FE = new FastaSequences(dir+"/Fend/"+dirs.get(i)+"/Trinity.fasta");
					mergeSequencesFEnd(FE,dirs.get(i));
					//					System.out.println(i);
				}
				if((i*100)/dirs.size()>percent){
					System.out.print(""+percent+"%");
					percent +=10;
				}
				if((i*100)/dirs.size()>point){
					System.out.print(".");
					point +=1;
				}
			}



			System.out.println("..100%");
			System.out.print("Reading directory for 3' ends....");
			//dirs = IOTools.getDirectories(dir+"/Tend");
			dirs = Functions.getStrings(dir+"/Tend.txt");
			System.out.println("Finished");
			System.out.println("Merging 3' ends....");
			System.out.print("Percent finished..0%..");
			percent = 10;
			for(int i = 0; i < dirs.size();i++){
				if(IOTools.fileExists(dir+"/Tend/"+dirs.get(i)+"/Trinity.fasta")){
					FastaSequences TE = new FastaSequences(dir+"/Tend/"+dirs.get(i)+"/Trinity.fasta");
					mergeSequencesTEnd(TE,dirs.get(i));
				}
				if(i/dirs.size()>percent){
					System.out.print(".."+percent+"%..");
					percent +=10;
				}
			}
			System.out.println("..100%");
		}
		catch(Exception E){E.printStackTrace();}
	}

	void mergeTrinityReadsFalse(String dir){
		try{
			ArrayList<String> dirs = IOTools.getDirectories(dir+"/Fend");
			System.out.println("Merging 5' ends....");
			System.out.print("Percent finished..0%..");
			int percent = 10; 
			for(int i = 0; i < dirs.size();i++){
				if(IOTools.fileExists(dir+"/Fend/"+dirs.get(i)+"/Trinity.fasta")){
					FastaSequences FE = new FastaSequences(dir+"/Fend/"+dirs.get(i)+"/Trinity.fasta");
					mergeSequencesTEnd(FE,dirs.get(i));
				}
				if(i/dirs.size()>percent){
					System.out.print(".."+percent+"%..");
					percent +=10;
				}
			}
			System.out.println("..100%");
			dirs = IOTools.getDirectories(dir+"/Tend");
			//			ArrayList<String> dirs = IOTools.getDirectories(dir+"/Tend");
			System.out.println("Merging 3' ends....");
			System.out.print("Percent finished..0%..");
			percent = 10;
			for(int i = 0; i < dirs.size();i++){
				if(IOTools.fileExists(dir+"/Tend/"+dirs.get(i)+"/Trinity.fasta")){

					FastaSequences TE = new FastaSequences(dir+"/Tend/"+dirs.get(i)+"/Trinity.fasta");
					mergeSequencesFEnd(TE,dirs.get(i));
				}
				if(i/dirs.size()>percent){
					System.out.print(".."+percent+"%..");
					percent +=10;
				}
			}
			System.out.println("..100%");
		}
		catch(Exception E){E.printStackTrace();}
	}



	private void mergeSequencesFEnd(FastaSequences FE, String sequence){
		if(this.containsKey(sequence)){
			System.out.println(sequence);
			FastaSequences temp = this.get(sequence);
			FastaSequences nTemp = new FastaSequences();
			nTemp.Name = sequence;
			for(int i = 0; i < temp.size();i++){
				FastaSequence tempSeq = temp.get(i);
				int [] seq = tempSeq.getSequence();
				boolean better = false;
				for(int j = 0; j < FE.size();j++){
					//					System.out.println("test");
					FastaSequence tempSeqFP = FE.get(j);
					int [] seqFP = tempSeqFP.getSequence();
					int[] newSeq = RNAfunctions.merge5end(seqFP, seq, length, missmatches, inside);
					if(newSeq.length > seq.length){
						better = true; 
						FastaSequence nSeq = new FastaSequence(sequence,newSeq); 
						nTemp.add(nSeq);
					}
				}
				if(!better)
					nTemp.add(tempSeq);
			}
			this.put(sequence, nTemp);
		}
	}

	private void mergeSequencesTEnd(FastaSequences FE, String sequence){
		if(this.containsKey(sequence)){
			FastaSequences temp = this.get(sequence);
			FastaSequences nTemp = new FastaSequences();
			nTemp.Name= sequence;
			for(int i = 0; i < temp.size();i++){
				FastaSequence tempSeq = temp.get(i);
				int [] seq = tempSeq.getSequence();
				boolean better = false;
				for(int j = 0; j < FE.size();j++){
					FastaSequence tempSeqFP = FE.get(j);
					int [] seqFP = tempSeqFP.getSequence();
					int[] newSeq = RNAfunctions.merge3end(seqFP, seq, length, missmatches, inside);
					if(newSeq.length > seq.length){
						better = true; 
						FastaSequence nSeq = new FastaSequence(sequence,newSeq); 
						nTemp.add(nSeq);
					}
				}
				if(!better)
					nTemp.add(tempSeq);
			}
			this.put(sequence, nTemp);
		}
	}

	private void removeLow( int cutoff){
		System.out.println("removing ends with less than "+cutoff+ " reads");
		for (Enumeration<String> e = this.keys(); e.hasMoreElements();){
			String hitName = e.nextElement();
			if(this.get(hitName).size()<= cutoff){
				this.remove(hitName);
			}
		}
		System.out.println(this.size()+ " ends are left.");


	}

	private void printFasta(String dir){
		if(!IOTools.isDir(dir))IOTools.mkDir(dir);
		for (Enumeration<String> e = this.keys(); e.hasMoreElements();){
			String hitName = e.nextElement();
			FastaSequences temp = this.get(hitName);
			temp.printFasta(dir, hitName+".fa");
		}
	}


	private void callTrinity(String dir){
		Trinity newDir = new Trinity();

	}





}




