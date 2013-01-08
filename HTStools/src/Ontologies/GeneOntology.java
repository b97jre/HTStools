package Ontologies;

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

public class GeneOntology {
	protected Hashtable <String,GOClass> GOMolFun;
	protected Hashtable <String,GOClass> GOBioProc;
	protected Hashtable <String,GOClass> GOCellCom;
	protected Hashtable <String,GOGene> GOGenes;

	public GeneOntology(){
		this.GOMolFun = new Hashtable <String,GOClass>();
		this.GOBioProc = new Hashtable <String,GOClass>();
		this.GOCellCom = new Hashtable <String,GOClass>();

		this.GOGenes = new Hashtable <String,GOGene>(); 
	}

	public  void run(Hashtable<String,String> T){

		System.out.println("running ontology GO...");
			String dir =  Functions.getValue(T, "-d", IOTools.getCurrentPath());
			String inFile = Functions.getValue(T, "-GOin");
			String outFile = Functions.getValue(T, "-GOout",inFile+".GeneOntology");
			this.addClasses(dir+"/"+inFile);

			this.addGenes(dir+"/"+inFile);
			this.printGenes(dir+"/"+outFile);
	}

	public void printGenes(String outFile){
		try{

			ExtendedWriter EW= new ExtendedWriter(new FileWriter(outFile));
			for (Enumeration<String> e = this.GOGenes.keys(); e.hasMoreElements();){
				GOGene A = this.GOGenes.get(e.nextElement());
				ArrayList <GOClass> GOclasses = A.GOclasses;
				EW.print(A.geneName);
				EW.print("\t");
				int count = 0;
				for(int i = 0; i < GOclasses.size();i++){
					if(this.GOMolFun.containsKey(GOclasses.get(i).GO_ID)){
						if(count != 0)EW.print(", ");
						EW.print(GOclasses.get(i).term+"#"+GOclasses.get(i).GO_ID);
						count++;
					}
				}
				EW.print("\t");
				count = 0;
				for(int i = 0; i < GOclasses.size();i++){
					if(this.GOBioProc.containsKey(GOclasses.get(i).GO_ID)){
						if(count != 0)EW.print(", ");
						EW.print(GOclasses.get(i).term+"#"+GOclasses.get(i).GO_ID);
						count++;
					}
				}
				EW.print("\t");
				count = 0;
				for(int i = 0; i < GOclasses.size();i++){
					if(this.GOCellCom.containsKey(GOclasses.get(i).GO_ID)){
						if(count != 0)EW.print(", ");
						EW.print(GOclasses.get(i).term+"#"+GOclasses.get(i).GO_ID);
						count++;
					}
				}
			}
			EW.flush();
			EW.close();
		}catch(Exception E){
			E.printStackTrace();
		}
	}


	public void printGOinfo(GOGene A,ExtendedWriter EW){
		try{
			ArrayList <GOClass> GOclasses = A.GOclasses;
			int count = 0;
			for(int i = 0; i < GOclasses.size();i++){
				if(this.GOMolFun.containsKey(GOclasses.get(i).GO_ID)){
					if(count != 0)EW.print(", ");
					EW.print(GOclasses.get(i).term+"#"+GOclasses.get(i).GO_ID);
					count++;
				}
			}
				EW.print("\t");
			count = 0;
			for(int i = 0; i < GOclasses.size();i++){
				if(this.GOBioProc.containsKey(GOclasses.get(i).GO_ID)){
					if(count != 0)EW.print(", ");
					EW.print(GOclasses.get(i).term+"#"+GOclasses.get(i).GO_ID);
					count++;
				}
			}
				EW.print("\t");
			count = 0;
			for(int i = 0; i < GOclasses.size();i++){
				if(this.GOCellCom.containsKey(GOclasses.get(i).GO_ID)){
					if(count != 0)EW.print(", ");
					EW.print(GOclasses.get(i).term+"#"+GOclasses.get(i).GO_ID);
					count++;
				}
			}
		}catch(Exception E){
			E.printStackTrace();
		}
	}



	public void linkGenes2Classes(Hashtable <String,Gene> Genes){
		for (Enumeration<String> e =Genes.keys(); e.hasMoreElements();){
			String queryName = e.nextElement();
			String hitName = Genes.get(queryName).getBestHit();
			if(hitName != null && this.GOGenes.containsKey(hitName)){
				Genes.get(queryName).GOgene = this.GOGenes.get(hitName);
			}
		}
	}


	public void addClasses(String infile){
		System.out.print("adding classes");
		try{
			ExtendedReader ER= new ExtendedReader(new FileReader(infile));
			while(ER.more()){
				String line = ER.readLine();
				String[] info = line.split("\t");
				//System.out.println(info[0]+"\t"+info[5]);
				String GO_ID = info[5];
				String GOpart = info[7];
				String GOterm = info[4];
				
				boolean molFun,bioPro,cellCom;
				molFun=bioPro=cellCom=false;

				if(GOpart.compareTo("F") == 0){molFun=true;bioPro=false;cellCom=false;}
				else if(GOpart.compareTo("P") == 0){molFun=false;bioPro=true;cellCom=false;}
				else if(GOpart.compareTo("C") == 0){molFun=false;bioPro=false;cellCom=true;}

				GOClass temp = new GOClass(GO_ID,GOterm);

				if(molFun)
					this.GOMolFun.put(temp.GO_ID, temp);
				else if(cellCom)
					this.GOCellCom.put(temp.GO_ID, temp);
				else if(bioPro)
					this.GOBioProc.put(temp.GO_ID, temp);
				else
					System.out.println(line);
			}
			ER.close();
		}catch(Exception E){
			E.printStackTrace();
		}
	}

	public void addGenes(String infile){
		System.out.print("adding genes");
		try{
			ExtendedReader ER= new ExtendedReader(new FileReader(infile));
			while(ER.more()){
				String line = ER.readLine();
				String[] info = line.split("\t");
	//			System.out.println(info[0]+"\t"+info[5]);
				String GO_ID = info[5];
				String GOpart = info[7];
				String GOterm = info[4];
				
				if(!this.GOGenes.containsKey(info[0])){
					this.GOGenes.put(info[0], new GOGene(info[0]));
				}
				if(this.GOBioProc.containsKey(info[5]))
					this.GOGenes.get(info[0]).addGOterm(this.GOBioProc.get(info[5]));
				else if(this.GOCellCom.containsKey(info[5]))
					this.GOGenes.get(info[0]).addGOterm(this.GOCellCom.get(info[5]));
				else if(this.GOMolFun.containsKey(info[5]))
					this.GOGenes.get(info[0]).addGOterm(this.GOMolFun.get(info[5]));
				else
					System.out.println(info[0]+"\t"+info[5]+" NOT FOUND");

				
			}
			ER.close();
		}catch(Exception E){
			E.printStackTrace();
		}
		System.out.println("... "+this.GOGenes.size()+" genes added");
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
