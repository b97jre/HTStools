package Ontologies;

import general.ExtendedReader;
import general.ExtendedWriter;
import general.Functions;
import general.IOTools;

import java.io.FileReader;
import java.io.FileWriter;
import java.util.Enumeration;
import java.util.Hashtable;

import alignment.Gene;

import Sequence.FastaSequences;

public class Panther {
	protected Hashtable <String,PantherClass> PantherClasses;
	protected Hashtable <String,PantherGene> PantherGenes;

	public Panther(){
		this.PantherClasses = new Hashtable <String,PantherClass>();
		this.PantherGenes = new Hashtable <String,PantherGene>(); 
	}

	public void run(Hashtable<String,String> T){
		System.out.println("running ontology panther...");
		if(T.containsKey("-subset")){
			String subset = Functions.getValue(T, "-subset", "subset.txt");
			String dir =  Functions.getValue(T, "-d", IOTools.getCurrentPath());
			String inFile =  Functions.getValue(T, "-i", "fileName");
			String outFile =  Functions.getValue(T, "-o", inFile+".subset.fa");
			System.out.println("writing sequences found in "+ inFile +" and "+subset+ "in file "+outFile);
			if(inFile.compareTo("fileName") != 0){
				Hashtable<String,String> HT = new Hashtable<String,String>(1000000);
				HT = Functions.getKeys(HT, subset);
				printPantherSubset(dir,inFile,outFile,HT);
			}
		}
		if(T.containsKey("-extractGO")){
			String dir =  Functions.getValue(T, "-d", IOTools.getCurrentPath());
			String hmm_classificationsFile = Functions.getValue(T, "-pantherRef", "PANTHER7.2_HMM_classifications");
			String inFile = Functions.getValue(T, "-i");
			String outFile = Functions.getValue(T, "-o",inFile+".pantherClasses");
			this.addClasses(dir+"/"+hmm_classificationsFile);
			this.addGenes(dir+"/"+inFile);
			this.linkGenes2Classes();
			this.printGenes(dir+"/"+outFile);
		}
		else{
			String dir =  Functions.getValue(T, "-d", IOTools.getCurrentPath());
			String hmm_classificationsFile = Functions.getValue(T, "-pantherRef", "PANTHER7.2_HMM_classifications");
			String inFile = Functions.getValue(T, "-pantherIn");
			String outFile = Functions.getValue(T, "-pantherOut",inFile+".pantherClasses");
			this.addClasses(dir+"/"+hmm_classificationsFile);
			this.addGenes(dir+"/"+inFile);
			this.linkGenes2Classes();
			this.printGenes(dir+"/"+outFile);
		}
		
	}


	public void printGenes(String outFile){
		try{
	
			ExtendedWriter EW= new ExtendedWriter(new FileWriter(outFile));
			for (Enumeration<String> e = this.PantherGenes.keys(); e.hasMoreElements();){
				this.PantherGenes.get(e.nextElement()).printClass(EW);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){
			E.printStackTrace();
		}


	}

	public void linkGenes2Classes(){
		for (Enumeration<String> e = this.PantherGenes.keys(); e.hasMoreElements();){
			String hitName = e.nextElement();
			PantherGene temp = this.PantherGenes.get(hitName);
			if(PantherClasses.containsKey(temp.PantherClass)){
				temp.PC = PantherClasses.get(temp.PantherClass);
			}
			this.PantherGenes.put(hitName, temp);
		}
	}

	public void linkGenes2PantherGenes(Hashtable <String,Gene> Genes){
		for (Enumeration<String> e = Genes.keys(); e.hasMoreElements();){
			String hitName = e.nextElement();
			if(PantherGenes.containsKey(hitName)){
				Genes.get(hitName).pantherGene = PantherGenes.get(hitName);
			}
		}
	}
	

	public void addClasses(String infile){
		System.out.print("adding classes");
		try{
			ExtendedReader ER= new ExtendedReader(new FileReader(infile));
			while(ER.more()){
				String line = ER.readLine();
				PantherClass temp = new PantherClass();
				if(temp.addInfo(line)){
					this.PantherClasses.put(temp.PANTHER_ID, temp);
				}else
					System.out.println(line);
			}
			ER.close();
		}catch(Exception E){
			E.printStackTrace();
		}
		System.out.println("... "+this.PantherClasses.size()+" classes added");
	}

	public void addGenes(String infile){
		System.out.print("adding genes");
		try{
			ExtendedReader ER= new ExtendedReader(new FileReader(infile));
			while(ER.more()){
				String line = ER.readLine();
				PantherGene temp = new PantherGene();
				if(temp.addInfo(line)){
					this.PantherGenes.put(temp.geneName, temp);
				}else
					System.out.println(line);
			}
			ER.close();
		}catch(Exception E){
			E.printStackTrace();
		}
		System.out.println("... "+this.PantherGenes.size()+" genes added");
	}




	public static void printPantherSubset(String dir, String infile, String outfile, Hashtable<String,String> HT){
		try{
			ExtendedWriter EW= new ExtendedWriter(new FileWriter(dir+"/"+outfile));
			ExtendedReader ER= new ExtendedReader(new FileReader(dir+"/"+infile));

			printPantherSubset(ER, HT, EW);
			ER.close();
			EW.flush();
			EW.close();
		}catch(Exception E){
			E.printStackTrace();
		}
	}

	private static void printPantherSubset(ExtendedReader ER,Hashtable<String,String> HT, ExtendedWriter EW){
		while(ER.more()){
			String InfoLine = ER.readLine();
			String[] info = InfoLine.split("\t");
			if(HT.containsKey(info[0])){
				EW.println(InfoLine);
			}
		}
	}

	public static void extractGOterms(String dir, String infile, String outfile, Hashtable<String,String> HT){
		try{
			ExtendedWriter EW= new ExtendedWriter(new FileWriter(dir+"/"+outfile));
			ExtendedReader ER= new ExtendedReader(new FileReader(dir+"/"+infile));

			extractGOterms(ER, HT, EW);
			ER.close();
			EW.flush();
			EW.close();
		}catch(Exception E){
			E.printStackTrace();
		}
	}
	private static void extractGOterms(ExtendedReader ER,Hashtable<String,String> HT, ExtendedWriter EW){
		while(ER.more()){
			String InfoLine = ER.readLine();
			String[] info = InfoLine.split("\t");
			if(HT.containsKey(info[0])){
				EW.println(InfoLine);
			}
		}
	}




}
