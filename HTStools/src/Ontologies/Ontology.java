package Ontologies;

import general.ExtendedWriter;
import general.Functions;
import general.IOTools;

import java.io.FileWriter;
import java.util.Enumeration;
import java.util.Hashtable;

import Blast.Blast;
import Sequence.FastaSequences;
import alignment.Gene;

public class Ontology {

	public Hashtable <String,Gene> Genes;

	public static void main(String []args){
		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}
		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);
		Ontology.run(T);
	}

	public static void run(Hashtable<String,String> T){

		String program = Functions.getValue(T, "-p", "h").toUpperCase();

		String fastaFile = Functions.getValue(T, "-ref");
		String dir =  Functions.getValue(T, "-d", IOTools.getCurrentPath());
		String outFile = Functions.getValue(T, "-refOut",fastaFile+".annotation");
		Ontology O = new Ontology();
		O.getFastaNames(dir, fastaFile);
		if(T.containsKey("-blastFile")){
			O.Genes = Blast.parseBlastFile(dir, Functions.getValue(T, "-blastFile"), O.Genes);

		}
		if(program.indexOf("PANTHER") != -1){
			Panther A = new Panther();
			A.run(T);
			A.linkGenes2PantherGenes(O.Genes);
		}
		//		if(program.indexOf("PFAM") != -1){
		//SequenceHandling.MergeSequences(T);
		//		}
		GeneOntology A = new GeneOntology();
		if(program.indexOf("GO") != -1){
			A.run(T);
			A.linkGenes2Classes(O.Genes);
		}

		O.printGenesWithOntology(dir,outFile,A);

	}



	private void getFastaNames(String dir, String fileName){
		FastaSequences A = new FastaSequences();
		A.getFastaFileNames(dir, fileName);
		this.Genes = A.convertToHashTable();
	}


	private void printGenesWithOntology(String dir, String fileName, GeneOntology A){
		try{
			ExtendedWriter EW= new ExtendedWriter(new FileWriter(dir+"/"+fileName));

			EW.println("ORF\tlength\tPantherClass\tPantherName\tArabidopsis\tGO molecular functions\tGO Biological processes\tGO Cellular Compartment\tProtein Class\tPathway");
			for (Enumeration<String> e =Genes.keys(); e.hasMoreElements();){
				Genes.get(e.nextElement()).printAnnotation(EW,A);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}


}

