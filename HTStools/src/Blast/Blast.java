package Blast;

import general.ExtendedReader;
import general.ExtendedWriter;
import general.Functions;
import general.IOTools;

import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;

import Sequence.FastaSequences;
import alignment.Gene;

public class Blast {


	public Hashtable<String,Gene> Genes;
	public Hashtable<String,Gene> DB;

	public static void main(String []args){
		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}
		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);
		Blast.run(T);
	}

	Blast(){}

	public static void run(Hashtable<String,String> T){
		String dir = Functions.getValue(T, "-d", IOTools.getCurrentPath());
		if(T.containsKey("-rd")){
			double cutoff = Double.parseDouble(Functions.getValue(T, "-rd", "0.75"));
			String inFile = Functions.getValue(T, "-i", ".");
			String outFile = Functions.getValue(T, "-o", inFile+".filtered");
			Blast blastDB = new Blast();
			System.out.println("parsing blast file");
			blastDB.parseInitialBlastFile(dir,inFile);
			System.out.println("finished");
			System.out.println("removing unsimilair matches");
			blastDB.parseAdditionalBlastFile(dir,inFile,cutoff);
			System.out.println("finished");
			if(T.containsKey("-print")){
				blastDB.printBlastHits(dir+"/"+outFile+".RD_"+Functions.getValue(T, "-rd", "0.75"));
			}
			System.out.println("Identifying very similair matches");
			System.out.println("removing no longer present blastHits");
			blastDB.removeNotPresentHits();
			System.out.println("finished");
			if(T.containsKey("-print")){
				System.out.println("printing matches");
				blastDB.printBlastHits(dir+"/"+outFile);
			}
			System.out.println("printing distribution");
			int bins = Integer.parseInt(Functions.getValue(T, "-distribution", "100"));
			blastDB.distribution(bins,dir+"/"+outFile+".RD_"+Functions.getValue(T, "-rd", "0.75")+".distribution");
			System.out.println("finished");
			if(T.containsKey("-cluster")){
				blastDB.checkDependencies();
				blastDB.printCluster(dir+"/"+outFile+".RD_"+Functions.getValue(T, "-rd", "0.75")+".clusters");
			}
		}
		if(T.containsKey("-remDup")){
			String inFile = Functions.getValue(T, "-i", ".");
			String outFile = Functions.getValue(T, "-o", inFile+".NoDup");
			Blast blastDB = new Blast();
			System.out.println("parsing blast file");
			blastDB.parseBlastFile(dir,inFile);
			System.out.println("finished");
			System.out.println("removing unsimilair matches");
			blastDB.removeWeakBlastHits(0.95);
			System.out.println("finished");
			System.out.println("removing duplicates and subsequences  A/B>0.9999");
			blastDB.removeDuplicates();
			System.out.println("finished");
			System.out.println("removing no longer present blastHits");
			blastDB.removeNotPresentHits();
			System.out.println("finished");
			System.out.println("printing matches to "+dir+"/"+outFile);
			blastDB.printBlastHits(dir+"/"+outFile);
			blastDB.printGenes(dir+"/"+outFile+".genes");
			System.out.println("finished");
		}



		if(T.containsKey("-distribution")){
			String inFile = Functions.getValue(T, "-i", ".");
			String outFile = Functions.getValue(T, "-o", inFile+".filtered");
			int bins = Integer.parseInt(Functions.getValue(T, "-distribution", "100"));
			Blast blastDB = new Blast();
			System.out.println("parsing blast file");
			blastDB.parseBlastFile(dir,inFile);
			System.out.println("finished");
			blastDB.distribution(bins,dir+"/"+outFile+".distribution");
		}
		if(T.containsKey("-merge")){
			String inFile = Functions.getValue(T, "-i", ".");
			String outFile = Functions.getValue(T, "-o", inFile+".merged.blast");
			int distance = Integer.parseInt(Functions.getValue(T, "-merge", "100"));
			int penalty = Integer.parseInt(Functions.getValue(T, "-penalty", "2"));
			Blast blastDB = new Blast();
			System.out.println("parsing blast file");
			blastDB.parseBlastFile(dir,inFile);
			System.out.println("finished");
			System.out.println("parsing blast file.....");
			blastDB.merge(distance,penalty);
			System.out.println("finished");
			System.out.println("print file....");
			blastDB.printBlastHits(dir+"/"+outFile);
			System.out.println("finished");
		}
		if(T.containsKey("-mergeTotal")){
			String inFile = Functions.getValue(T, "-i", ".");
			String outFile = Functions.getValue(T, "-o", inFile+".coverage");
			String fastaFile = Functions.getValue(T, "-fasta", inFile+".fa");
			Blast blastDB = new Blast();
			blastDB.getFastaNames(dir,fastaFile);
			System.out.println("parsing blast file");
			blastDB.parseBlastFileNoNew(dir,inFile);
			System.out.println("finished");
			System.out.println("Checking coverage.....");
			blastDB.mergeTotal(dir+"/"+outFile);
			System.out.println("finished");
		}


		if(T.containsKey("-bestHit")){
			String inFile = Functions.getValue(T, "-i", ".");
			String outFile = Functions.getValue(T, "-o", inFile+".bestHit");
			String fastaFile = Functions.getValue(T, "-fasta", inFile+".merged");
			Blast blastDB = new Blast();
			System.out.println("parsing fasta file");
			blastDB.getFastaNames(dir,fastaFile);
			System.out.println("parsing blast file");
			blastDB.parseBlastFile(dir,inFile);
			System.out.println("finished");
			blastDB.printBestHit(dir+"/"+outFile);
		}
		if(T.containsKey("-bestCoverage")){
			String inFile = Functions.getValue(T, "-i", ".");
			String outFile = Functions.getValue(T, "-o", inFile+".bestCoverage");
			String fastaFile = Functions.getValue(T, "-fasta", inFile+".fa");
			Blast blastDB = new Blast();
			System.out.println("parsing fasta file");
			blastDB.getFastaNames(dir,fastaFile);
			System.out.println("parsing blast file");
			blastDB.parseBlastFileNoNew(dir,inFile);
			System.out.println("finished");
			blastDB.printBestHitCoverage(dir+"/"+outFile);
		}


		if(T.containsKey("-cutoff")){
			double cutoff = Double.parseDouble(Functions.getValue(T, "-cutoff", "99"));
			String inFile = Functions.getValue(T, "-i", ".");
			String outFile = Functions.getValue(T, "-o", inFile+"."+cutoff+".blast");
			String fastaFile = Functions.getValue(T, "-fasta", inFile+".fa");
			Blast blastDB = new Blast();
			System.out.println("parsing fasta file");
			blastDB.getFastaNames(dir,fastaFile);
			System.out.println("parsing blast file");
			blastDB.parseBlastFileNoNew(dir,inFile,cutoff);
			System.out.println("finished");
			blastDB.printBlastHits(dir+"/"+outFile);
		}

		if(T.containsKey("-subset")){
			String inFile = Functions.getValue(T, "-i", ".");
			String outFile = Functions.getValue(T, "-o", inFile+".subset.blast");
			String fastaFile = Functions.getValue(T, "-fasta", inFile+".fa");
			String subset =  Functions.getValue(T, "-subset", inFile+".fa");
			Blast blastDB = new Blast();
			blastDB.parseNameFile(dir, subset);
			System.out.println("parsing fasta file");
			blastDB.getFastaNames(dir,fastaFile);
			System.out.println("parsing blast file");
			blastDB.parseBlastFileNoNewSubset(dir,inFile);
			System.out.println("finished");

			blastDB.printBlastHits(dir+"/"+outFile);
		}



	}


	private void getFastaNames(String dir, String fileName){
		FastaSequences A = new FastaSequences();
		A.getFastaFileNames(dir, fileName);
		this.Genes = Gene.convertToHashTable(A);
	}


	private void distribution(int bins,String file){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(file));
			int[] distribution = new int[bins+1]; 
			for (Enumeration<String> e = Genes.keys(); e.hasMoreElements();){
				int count = Genes.get(e.nextElement()).getNrOfBlastHits();
				if(count < bins) distribution[count]++;
				else distribution[bins]++;
			}

			for(int i = 0; i < distribution.length; i++){
				EW.println(i+"\t"+distribution[i]);
			}
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}
	}





	private void merge(int distance, int penalty){
		try{
			for (Enumeration<String> e = Genes.keys(); e.hasMoreElements();){
				Genes.get(e.nextElement()).mergeBlastHits(distance,penalty);
			}
		}
		catch(Exception E){E.printStackTrace();}



	}

	private void mergeTotal(String outFile){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outFile));
			for (Enumeration<String> e = Genes.keys(); e.hasMoreElements();){
				Genes.get(e.nextElement()).mergeBlastHitsTotal(EW);
			}
			EW.flush();
			EW.close();

		}
		catch(Exception E){E.printStackTrace();}



	}



	private void parseInitialBlastFile(String dir, String inFile ){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(dir+"/"+inFile));
			//trinity_4_1_1	trinity_4_1_1	100.00	79	0	0	1	79	1	79	1e-27	 118

			while(ER.more()){

				String Line = ER.readLine();
				String[] info = Line.split("\t");

				String	queryName = info[0];
				String hitName = info[1];
				if(queryName.compareTo(hitName) == 0){
					double similarity = Double.parseDouble(info[2]);
					int length = Integer.parseInt(info[3]);
					int missmatches = Integer.parseInt(info[4]);
					int gaps = Integer.parseInt(info[5]);
					int queryStart = Integer.parseInt(info[6]);
					int queryStop = Integer.parseInt(info[7]);
					int hitStart = Integer.parseInt(info[8]);
					int hitStop = Integer.parseInt(info[9]);
					double Evalue = Double.parseDouble(info[10]);
					double score = Double.parseDouble(info[11]);

					BlastHit B = new BlastHit(queryName,hitName,similarity,length,
							missmatches,gaps,queryStart,queryStop,hitStart,hitStop,Evalue, score);
					Gene A = null;
					if(Genes == null)Genes = new Hashtable<String, Gene>(1000000); 
					if(Genes.containsKey(queryName)){
						A = Genes.get(queryName);
						A.addBlastHit(B);
					}
					else{
						A = new Gene(queryName);
						A.addBlastHit(B);
						Genes.put(queryName, A);
					}
					Gene C = null;
					if(Genes.containsKey(hitName)){
						C = Genes.get(hitName);
						C.addBlastHit(B);
					}
					else{
						C = new Gene(hitName);
						C.addBlastHit(B);
						Genes.put(hitName, C);
					}

				}	
			}
			ER.close();
		}
		catch(Exception E){E.printStackTrace();}
	}	


	private void parseAdditionalBlastFile(String dir, String inFile, double cutoff){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(dir+"/"+inFile));
			//trinity_4_1_1	trinity_4_1_1	100.00	79	0	0	1	79	1	79	1e-27	 118

			while(ER.more()){

				String Line = ER.readLine();
				String[] info = Line.split("\t");

				String	queryName = info[0];
				String hitName = info[1];
				if(queryName.compareTo(hitName) != 0){
					double similarity = Double.parseDouble(info[2]);
					int length = Integer.parseInt(info[3]);
					int missmatches = Integer.parseInt(info[4]);
					int gaps = Integer.parseInt(info[5]);
					int queryStart = Integer.parseInt(info[6]);
					int queryStop = Integer.parseInt(info[7]);
					int hitStart = Integer.parseInt(info[8]);
					int hitStop = Integer.parseInt(info[9]);
					double Evalue = Double.parseDouble(info[10]);
					double score = Double.parseDouble(info[11]);

					BlastHit B = new BlastHit(queryName,hitName,similarity,length,
							missmatches,gaps,queryStart,queryStop,hitStart,hitStop,Evalue, score);
					Gene A = null;
					if(Genes == null)Genes = new Hashtable<String, Gene>(100000); 
					if(Genes.containsKey(queryName)){
						A = Genes.get(queryName);
						if(A.isAboveCutoff(score, cutoff)){
							A.addBlastHit(B);
						}
					}
					else{
						System.out.println("this should never happen :"+ queryName );
						A = new Gene(queryName);
						A.addBlastHit(B);
						Genes.put(queryName, A);
					}
					Gene C = null;
					if(Genes.containsKey(hitName)){
						C = Genes.get(hitName);
						C.addBlastHit(B);
					}
					else{
						C = new Gene(hitName);
						C.addBlastHit(B);
						Genes.put(hitName, C);
					}

				}	
			}
			ER.close();
		}
		catch(Exception E){E.printStackTrace();}
	}	

	private void parseNameFile(String dir, String inFile){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(dir+"/"+inFile));
			while(ER.more()){
				String Line = ER.readLine();
				if(DB == null)DB = new Hashtable<String, Gene>(1000000); 
				if(!DB.containsKey(Line)){
					Gene A = new Gene(Line);
					DB.put(Line, A);
				}
			}			
			ER.close();
	}
		catch(Exception E){E.printStackTrace();}
	}	



	private void parseBlastFile(String dir, String inFile ){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(dir+"/"+inFile));
			//trinity_4_1_1	trinity_4_1_1	100.00	79	0	0	1	79	1	79	1e-27	 118

			int count = 100000;
			int pointer = 0;
			System.out.println(pointer +" lines read");
			while(ER.more()){
				pointer++;
				if(pointer > count){
					System.out.println(count +" lines read");
					count = count+100000;
				}

				String Line = ER.readLine();
				String[] info = Line.split("\t");

				String	queryName = info[0];
				String hitName = info[1];
				double similarity = Double.parseDouble(info[2]);
				int length = Integer.parseInt(info[3]);
				int missmatches = Integer.parseInt(info[4]);
				int gaps = Integer.parseInt(info[5]);
				int queryStart = Integer.parseInt(info[6]);
				int queryStop = Integer.parseInt(info[7]);
				int hitStart = Integer.parseInt(info[8]);
				int hitStop = Integer.parseInt(info[9]);
				double Evalue = Double.parseDouble(info[10]);
				double score = Double.parseDouble(info[11]);

				BlastHit B = new BlastHit(queryName,hitName,similarity,length,
						missmatches,gaps,queryStart,queryStop,hitStart,hitStop,Evalue, score);
				Gene A = null;
				if(Genes == null)Genes = new Hashtable<String, Gene>(1000000); 
				if(Genes.containsKey(queryName)){
					A = Genes.get(queryName);
					A.addBlastHit(B);
				}
				else{
					A = new Gene(queryName);
					A.addBlastHit(B);
					Genes.put(queryName, A);
				}
				Gene C = null;
				if(Genes.containsKey(hitName)){
					B.specificGene = Genes.get(hitName);
				}
				else{
					C = new Gene(hitName);
					Genes.put(hitName, C);
					B.specificGene = Genes.get(hitName);
				}
			}			
			ER.close();
	}
		catch(Exception E){E.printStackTrace();}
	}	

	private void parseBlastFileNoNew(String dir, String inFile){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(dir+"/"+inFile));
			//trinity_4_1_1	trinity_4_1_1	100.00	79	0	0	1	79	1	79	1e-27	 118

			int count = 100000;
			int pointer = 0;
			System.out.println(pointer +" lines read");
			while(ER.more()){
				pointer++;
				if(pointer > count){
					System.out.println(count +" lines read");
					count = count+100000;
				}

				String Line = ER.readLine();
				String[] info = Line.split("\t");

				String	queryName = info[0];
				String hitName = info[1];
				double similarity = Double.parseDouble(info[2]);
				int length = Integer.parseInt(info[3]);
				int missmatches = Integer.parseInt(info[4]);
				int gaps = Integer.parseInt(info[5]);
				int queryStart = Integer.parseInt(info[6]);
				int queryStop = Integer.parseInt(info[7]);
				int hitStart = Integer.parseInt(info[8]);
				int hitStop = Integer.parseInt(info[9]);
				double Evalue = Double.parseDouble(info[10]);
				double score = Double.parseDouble(info[11]);
				if(similarity > 98.00){
					BlastHit B = new BlastHit(queryName,hitName,similarity,length,
							missmatches,gaps,queryStart,queryStop,hitStart,hitStop,Evalue, score);
					Gene A = null;
					if(Genes == null)Genes = new Hashtable<String, Gene>(1000000); 
					if(Genes.containsKey(queryName)){

						A = Genes.get(queryName);
						A.addBlastHit(B);
					}
					else{
						System.out.println(queryName +"not found");
					}
				}
			}			
			ER.close();
		}
		catch(Exception E){E.printStackTrace();}
	}	


	private void parseBlastFileNoNewSubset(String dir, String inFile){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(dir+"/"+inFile));
			//trinity_4_1_1	trinity_4_1_1	100.00	79	0	0	1	79	1	79	1e-27	 118

			int count = 100000;
			int pointer = 0;
			System.out.println(pointer +" lines read");
			while(ER.more()){
				pointer++;
				if(pointer > count){
					System.out.println(count +" lines read");
					count = count+100000;
				}

				String Line = ER.readLine();
				String[] info = Line.split("\t");

				String	queryName = info[0];
				String hitName = info[1];
				double similarity = Double.parseDouble(info[2]);
				int length = Integer.parseInt(info[3]);
				int missmatches = Integer.parseInt(info[4]);
				int gaps = Integer.parseInt(info[5]);
				int queryStart = Integer.parseInt(info[6]);
				int queryStop = Integer.parseInt(info[7]);
				int hitStart = Integer.parseInt(info[8]);
				int hitStop = Integer.parseInt(info[9]);
				double Evalue = Double.parseDouble(info[10]);
				double score = Double.parseDouble(info[11]);
				if(similarity > 98.00 && length > 100){
					BlastHit B = new BlastHit(queryName,hitName,similarity,length,
							missmatches,gaps,queryStart,queryStop,hitStart,hitStop,Evalue, score);
					Gene A = null;
					if(Genes == null)Genes = new Hashtable<String, Gene>(1000000); 
					if(Genes.containsKey(queryName)){
						if(DB.containsKey(queryName) || DB.containsKey(hitName)){
							A = Genes.get(queryName);
							A.addBlastHit(B);
						}
					}
					else{
						System.out.println(queryName +"not found");
					}
				}
			}			
			ER.close();
		}
		catch(Exception E){E.printStackTrace();}
	}	


	private void parseBlastFileNoNew(String dir, String inFile,double cutoff){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(dir+"/"+inFile));
			//trinity_4_1_1	trinity_4_1_1	100.00	79	0	0	1	79	1	79	1e-27	 118

			int count = 100000;
			int pointer = 0;
			System.out.println(pointer +" lines read");
			while(ER.more()){
				pointer++;
				if(pointer > count){
					System.out.println(count +" lines read");
					count = count+100000;
				}

				String Line = ER.readLine();
				String[] info = Line.split("\t");

				String	queryName = info[0];
				String hitName = info[1];
				double similarity = Double.parseDouble(info[2]);
				int length = Integer.parseInt(info[3]);
				int missmatches = Integer.parseInt(info[4]);
				int gaps = Integer.parseInt(info[5]);
				int queryStart = Integer.parseInt(info[6]);
				int queryStop = Integer.parseInt(info[7]);
				int hitStart = Integer.parseInt(info[8]);
				int hitStop = Integer.parseInt(info[9]);
				double Evalue = Double.parseDouble(info[10]);
				double score = Double.parseDouble(info[11]);
				if(similarity > cutoff){
					BlastHit B = new BlastHit(queryName,hitName,similarity,length,
							missmatches,gaps,queryStart,queryStop,hitStart,hitStop,Evalue, score);
					Gene A = null;
					if(Genes == null)Genes = new Hashtable<String, Gene>(1000000); 
					if(Genes.containsKey(queryName)){

						A = Genes.get(queryName);
						A.addBlastHit(B);
					}
					else{
						System.out.println(queryName +"not found");
					}
				}
			}			
			ER.close();
		}
		catch(Exception E){E.printStackTrace();}
	}	

	public static Hashtable<String, Gene> parseBlastFile(String dir, String inFile, Hashtable<String, Gene> Genes ){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(dir+"/"+inFile));
			//trinity_4_1_1	trinity_4_1_1	100.00	79	0	0	1	79	1	79	1e-27	 118

			int count = 100000;
			int pointer = 0;
			System.out.println(pointer +" lines read");
			while(ER.more()){
				pointer++;
				if(pointer > count){
					System.out.println(count +" lines read");
					count = count+100000;
				}

				String Line = ER.readLine();
				String[] info = Line.split("\t");

				String	queryName = info[0];
				String hitName = info[1];
				if(info[1].indexOf("|") > -1){
					String[] info2 = info[1].split("\\|");
					if(info2[0].indexOf(".")>-1)hitName = info2[0].substring(0,info2[0].indexOf("."));
					else hitName = info2[0];
				}

				double similarity = Double.parseDouble(info[2]);
				int length = Integer.parseInt(info[3]);
				int missmatches = Integer.parseInt(info[4]);
				int gaps = Integer.parseInt(info[5]);
				int queryStart = Integer.parseInt(info[6]);
				int queryStop = Integer.parseInt(info[7]);
				int hitStart = Integer.parseInt(info[8]);
				int hitStop = Integer.parseInt(info[9]);
				double Evalue = Double.parseDouble(info[10]);
				double score = Double.parseDouble(info[11]);

				BlastHit B = new BlastHit(queryName,hitName,similarity,length,
						missmatches,gaps,queryStart,queryStop,hitStart,hitStop,Evalue, score);
				Gene A = null;
				if(Genes == null)Genes = new Hashtable<String, Gene>(1000000); 
				if(Genes.containsKey(queryName)){
					A = Genes.get(queryName);
					A.addBlastHit(B);
				}
				else{
					System.out.println();				
				}
			}			
			ER.close();
		}
		catch(Exception E){E.printStackTrace();}
		return Genes;
	}	


//
//	private void addHitsBlastFile(String dir, String inFile ){
//		try{
//			ExtendedReader ER = new ExtendedReader(new FileReader(dir+"/"+inFile));
//			//trinity_4_1_1	trinity_4_1_1	100.00	79	0	0	1	79	1	79	1e-27	 118
//
//
//			while(ER.more()){
//
//				String Line = ER.readLine();
//				String[] info = Line.split("\t");
//
//				String	queryName = info[0];
//				String hitName = info[1];
//				double similarity = Double.parseDouble(info[2]);
//				int length = Integer.parseInt(info[3]);
//				int missmatches = Integer.parseInt(info[4]);
//				int gaps = Integer.parseInt(info[5]);
//				int queryStart = Integer.parseInt(info[6]);
//				int queryStop = Integer.parseInt(info[7]);
//				int hitStart = Integer.parseInt(info[8]);
//				int hitStop = Integer.parseInt(info[9]);
//				double Evalue = Double.parseDouble(info[10]);
//				double score = Double.parseDouble(info[11]);
//
//				BlastHit B = new BlastHit(queryName,hitName,similarity,length,
//						missmatches,gaps,queryStart,queryStop,hitStart,hitStop,Evalue, score);
//				Gene A = null;
//				if(Genes.containsKey(queryName)){
//					A = Genes.get(queryName);
//					A.addBlastHit(B);
//				}
//			}			
//		}
//		catch(Exception E){E.printStackTrace();}
//	}	

	public void removeAllhits(){
		for (Enumeration<String> e = Genes.keys(); e.hasMoreElements();){
			Genes.get(e.nextElement()).removeAllHits();
		}
	}



	public void removeSimilairHits( double cutoff){
		for (Enumeration<String> e = Genes.keys(); e.hasMoreElements();){
			String temp = Genes.get(e.nextElement()).removeSimilairHits(cutoff);
			if(temp.length() > 1){
				String[] Names = temp.split("\t");
				for(int i = 0; i < Names.length; i++){
					if(Genes.containsKey(Names[i]))
						Genes.remove(Names[i]);
				}
			}
		}
	}





	public void removeDuplicates(){
		for (Enumeration<String> e = Genes.keys(); e.hasMoreElements();){
			String Name = e.nextElement();
			String temp = Genes.get(Name).removeSimilairHits(0.99);
			if(temp.length() > 1){
				String[] Names = temp.split("\t");
				System.out.println(Genes.get(Name).getName()+"\t"+temp);
				for(int i = 0; i < Names.length; i++){
					if(Genes.containsKey(Names[i]))
						Genes.remove(Names[i]);
				}
			}
		}
	}


	
	public void keepBestHits(){
//		for (Enumeration<String> e = Genes.keys(); e.hasMoreElements();){
//			String Name = e.nextElement();
//			Genes.get(Name).findBesthit();
//			if(temp.length() > 1){
//				String[] Names = temp.split("\t");
//				System.out.println(Genes.get(Name).getName()+"\t"+temp);
//				for(int i = 0; i < Names.length; i++){
//					if(Genes.containsKey(Names[i]))
//						Genes.remove(Names[i]);
//				}
//			}
//		}
	}
	
	
	
	
	public void removeWeakBlastHits(double cutoff){
		for (Enumeration<String> e = Genes.keys(); e.hasMoreElements();){
			Genes.get(e.nextElement()).removeWeakHits(cutoff);
		}
	}
	public void removeNotPresentHits(){
		for (Enumeration<String> e = Genes.keys(); e.hasMoreElements();){
			Genes.get(e.nextElement()).removeNotPresentHits(Genes);
		}
	}


	public void printBlastHits(String outFile){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outFile));
			for (Enumeration<String> e = Genes.keys(); e.hasMoreElements();){
				Genes.get(e.nextElement()).printBlastHits(EW);
			}
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}
	}


	public void printBestHit(String outFile){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outFile));
			for (Enumeration<String> e = Genes.keys(); e.hasMoreElements();){
				Genes.get(e.nextElement()).printBestHit(EW,DB);
			}
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}
	}

	public void printBestHitCoverage(String outFile){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outFile));

			EW.println("Query\tHit\tPercent\tlength\tmissmatches\tgaps\tqueryStart\tqueryStop\thitStart\thitStop\tEvalue\tBitScore");
			for (Enumeration<String> e = Genes.keys(); e.hasMoreElements();){

				Genes.get(e.nextElement()).printBestCoverage(EW);
			}
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}
	}

	public void checkDependencies(){
		try{
			//ExtendedWriter EW = new ExtendedWriter(new FileWriter(outFile));
			for (Enumeration<String> e = Genes.keys(); e.hasMoreElements();){
				Genes.get(e.nextElement()).checkDependencies(Genes);
			}
			for (Enumeration<String> e = Genes.keys(); e.hasMoreElements();){
				Genes.get(e.nextElement()).getBlastSize();
			}

			//EW.flush();
			//EW.close();
		}
		catch(Exception E){E.printStackTrace();}
	}

	public void printCluster(String outFile){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outFile));
			ExtendedWriter EW2 = new ExtendedWriter(new FileWriter(outFile+".distribution"));
			int[] dist = new int[10000];
			int maxCluster = 0;
			int count = 0; 
			for (Enumeration<String> e = Genes.keys(); e.hasMoreElements();){
				Gene temp = Genes.get(e.nextElement());
				ArrayList <String> genes =temp.clusterHits(EW,Genes, count);
				count++;

				if(genes.size()> maxCluster){
					maxCluster = genes.size();
				}
				if(genes.size()> dist.length){
					dist[dist.length-1]++;
				}
				else dist[genes.size()]++;

				for(int i = 0; i < genes.size();i++){
					Genes.remove(genes.get(i));

				}
			}
			for(int i = 0; i < Math.min(maxCluster+1,dist.length);i++){
				EW2.println(i+"\t"+dist[i]);
			}
			EW.flush();
			EW.close();
			EW2.flush();
			EW2.close();
		}
		catch(Exception E){E.printStackTrace();}
	}





	public void printGenes(String outFile){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outFile));
			for (Enumeration<String> e = Genes.keys(); e.hasMoreElements();){
				EW.println(e.nextElement());
			}
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}
	}



}
