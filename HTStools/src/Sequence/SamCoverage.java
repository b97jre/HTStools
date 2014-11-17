package Sequence;

import general.Functions;
import general.IOTools;
import general.RNAfunctions;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;


import variousTools.Trinity;
import general.ExtendedReader;
import general.ExtendedWriter;




public class SamCoverage extends Hashtable <String,FastaSequence> implements Serializable{

	/**
	 * 
	 */


	public static void main(String []args){
		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}
		System.out.println();

		char[] Bitflag = Integer.toString(Integer.parseInt("0"), 2).toCharArray();

		System.out.println(Bitflag);

		Bitflag = Integer.toString(Integer.parseInt("352"), 2).toCharArray();
		System.out.println(Bitflag);

		Bitflag = Integer.toString(Integer.parseInt("1"), 2).toCharArray();
		System.out.println(Bitflag);

		Hashtable<String,String> T = Functions.parseCommandLine(args);



		SamCoverage.findCoverage(T);
	}

	int[] distribution; 

	private static final long serialVersionUID = 1L;
	//private String GenomeName;
	SamCoverage(int maxLength){
		this.distribution = new int[maxLength];
	}




	public static void findCoverage(Hashtable<String,String> T){
		String dir = Functions.getValue(T, "-d", IOTools.getCurrentPath());
		if(T.containsKey("-ref")){
			String inFile = Functions.getValue(T, "-ref");
			if(T.containsKey("-sam") || T.containsKey("-o")){
				String samFile = Functions.getValue(T, "-sam", null);
				String outbase = Functions.getValue(T, "-o", samFile);
				int maxInsertSize = Functions.getInt(T, "-insertLength",1000);
				SamCoverage A = new SamCoverage(maxInsertSize);
				

				System.out.println("Reading reference files....");
				A.parseFastaSize(inFile);
				System.out.println("Finding covergae of reads");			
				A.getCoverage(dir,samFile);
				System.out.println("Printing genomeCov files....");
				A.printGenomeCovFiles(dir, outbase);
			}
		}
	}



	public void printGenomeCovFiles(String outDir, String genomeCovFile){
		if(!IOTools.isDir(outDir+"/genomeCov")) IOTools.mkDir(outDir+"/genomeCov");
		if(!IOTools.isDir(outDir+"/genomeCov"))IOTools.mkDir(outDir+"/genomeCov/Fend");
		if(!IOTools.isDir(outDir+"/genomeCov/Tend"))IOTools.mkDir(outDir+"/genomeCov/Tend");
		if(!IOTools.isDir(outDir+"/genomeCov/Total"))IOTools.mkDir(outDir+"/genomeCov/Total");
		if(!IOTools.isDir(outDir+"/wig")) IOTools.mkDir(outDir+"/wig");
		if(!IOTools.isDir(outDir+"/wig/Fend"))IOTools.mkDir(outDir+"/wig/Fend");
		if(!IOTools.isDir(outDir+"/wig/Tend"))IOTools.mkDir(outDir+"/wig/Tend");
		if(!IOTools.isDir(outDir+"/wig/Total"))IOTools.mkDir(outDir+"/wig/Total");
		if(!IOTools.isDir(outDir+"/LengthDistribution")) IOTools.mkDir(outDir+"/LengthDistribution");
		
		
		
		for (Enumeration<String> e = this.keys(); e.hasMoreElements();){
			String contig = e.nextElement();
			FastaSequence A = this.get(contig);
			System.out.println("Print Five prime end coverage");
			A.printFends(outDir+"/genomeCov/Fend/"+genomeCovFile);
			A.printFendsWig(outDir+"/wig/Fend/"+genomeCovFile);
			System.out.println("Print Three prime end coverage");
			A.printTends(outDir+"/genomeCov/Tend/"+genomeCovFile);
			A.printTendsWig(outDir+"/wig/Tend/"+genomeCovFile);
			System.out.println("Print Total coverage");
			A.printTotal(outDir+"/genomeCov/Total/"+genomeCovFile);
			A.printTotalWig(outDir+"/wig/Total/"+genomeCovFile);
			System.out.println("Print length distribution");
			ExtendedWriter EW = ExtendedWriter.getFileWriter(outDir+"/LengthDistribution/"+genomeCovFile+".LengthDistribution");
			EW.println("length\tcount");
			
			int totalCount = 0;
			for(int i = 1; i < this.distribution.length;i++){
				totalCount += distribution[i];
				EW.println(i+"\t"+distribution[i]);
			}
			System.out.println("totalNumber of sequences analyzed :"+ totalCount);
			EW.flush();
			EW.close();
		}



	}

	public void getCoverage(String dir, String inFile){
		try{
			ExtendedReader ER;
			if(inFile != null){
				System.out.println(inFile);
				ER = new ExtendedReader(new FileReader(dir+"/"+inFile));
			}
			else
				ER = new ExtendedReader(new InputStreamReader(System.in));
			//			int count = 1000;
			int line = 1;
			int lines = 1000000;
			while((char)ER.lookAhead() == '@') 
				ER.skipLine();

			String[] info = (ER.readLine().split("\t"));
			boolean pairedEnd = isPairedEnd(info[1]);
			if(pairedEnd){
				addCoverageFilePairedEnd(info);
				System.out.println("sam file contains paired end reads alignment");
			}
			else{
				addCoverageFileSingleEnd(info);
				System.out.println("sam file contains single end reads alignment");
			}
			if(pairedEnd){
				while(ER.more() ){
					addCoverageFilePairedEnd(ER.readLine().split("\t"));
					line++;
					if(line > lines){
						System.out.println("Nr of alignments read" + line);
						lines= lines + 1000000;
					}
				}
			}else{
				while(ER.more() ){
					addCoverageFileSingleEnd(ER.readLine().split("\t"));
					line++;
					if(line > lines){
						System.out.println("Nr of alignments read" + lines);
						lines= lines + 1000000;
					}
				}
			}

			System.out.println("Finished");
			System.out.println("total nr of alignments = " + line);
			ER.close();

		}catch(Exception E){E.printStackTrace();}
	}

	
	private boolean isPairedEnd(String samBitflag){
		char[] Bitflag = Integer.toString(Integer.parseInt(samBitflag), 2).toCharArray();
		if(Bitflag[Bitflag.length-1] == '0') return false;
		return true;
	}



	private void addCoverageFilePairedEnd(String[] info){
		char[] Bitflag = Integer.toString(Integer.parseInt(info[1]), 2).toCharArray();
		if(Bitflag[Bitflag.length-7] == '0' || Bitflag[Bitflag.length-2] == '0') return;
		String contig  = info[2];
		int leftPos = Integer.parseInt(info[3]);
		int rightPos = Integer.parseInt(info[7]);
		int length =  Integer.parseInt(info[8]);
		if(Bitflag.length < 9 || Bitflag[Bitflag.length-9] == '0'){
			if(Math.abs(length)< this.distribution.length){
				this.distribution[Math.abs(length)]++;
			}
			else{
				this.distribution[this.distribution.length-1]++;
			}
		}

		if(Bitflag.length > 4 && Bitflag[Bitflag.length-5] == '1' ){
			//reverse
			//HWI-ST1299:92:H04D0ADXX:1:2102:5238:71962	83	gi|15674250|ref|NC_002737.1|	89426	255	68M	=	89396	-98	TTGTAAGTCGATGAGATCATAGTGCCCTTCCTTTGGTAGTTGTGTCAGACTAGTGTAATGGTCTTTGG	B@>?A;EA>ECIFECHEGG>GEGHBHFCICHHEIGGIIIJJJIJIIHBIJJJIIJIJJJJJHHHGHFE	NH:i:1	HI:i:1	AS:i:119	nM:i:3
			int fend = rightPos - length-1;
			int tend = rightPos;
			this.get(contig).addCoverage(fend,tend);
		}else{
			//forward
			int fend = leftPos;
			int tend = leftPos + length-1;
			this.get(contig).addCoverage(fend,tend);
		}
	}

	private void addCoverageFileSingleEnd(String[] info){

		char[] Bitflag = Integer.toString(Integer.parseInt(info[1]), 2).toCharArray();
		String contig  = info[2];
		int leftPos = Integer.parseInt(info[3]);
		int length = 0;
		try{
			length = Integer.parseInt(info[5].substring(0,info[5].length()-1));
		}
		catch(NumberFormatException e){
			return;
		}

		if(Bitflag.length < 9 || Bitflag[Bitflag.length-9] == '0'){
			if(Math.abs(length)< this.distribution.length){
				this.distribution[Math.abs(length)]++;
			}
			else{
				this.distribution[this.distribution.length-1]++;
			}
		}

		if(Bitflag.length > 4 && Bitflag[Bitflag.length-5] == '1' ){ // reverse
			// reverse
			int tend = leftPos;
			int fend = leftPos+length - 1;
			this.get(contig).addCoverage(fend,tend);

		}else{
			//forward
			int fend = leftPos;
			int tend = leftPos+length - 1;
			this.get(contig).addCoverage(fend,tend);

		}
	}




	public void parseFastaSize(String inFile){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(inFile));
			int count = 1000;
			while(ER.more()){
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					addFastaSequenceLength(ER);
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

	private void addFastaSequenceLength(ExtendedReader ER){
		String[] hitInfo = ER.readLine().split(" ");
		String hitName = hitInfo[0].substring(1);
		int length = 0;
		while(ER.more() && (char)ER.lookAhead() != '>'){
			length += ER.readLine().length();
		}

		FastaSequence A = new FastaSequence(hitName, length);
		if(!this.containsKey(hitName)){
			this.put(hitName, A);
		}else{
			System.out.println("This already existed in the fastaFile: "+hitName  );
		}  
	}




}




