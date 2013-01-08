package MutationalAnalysis;

import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import structure.IntramolecularStructures;

import general.ExtendedReader;
import general.ExtendedWriter;
import general.Functions;
import general.IOTools;
import general.RNAfunctions;

import Ontologies.PantherClass;
import Sequence.FastaSequence;
import Sequence.FastaSequences;

public class MutationalAnalysis {

	private String Name;
	private int[] Sequence;
	private IntramolecularStructures Structure;

	ArrayList <Mutations> mutants;
	ArrayList <Mutations> compMutants;



	private Change[] singleMutation;
	private Change[][] doubleMutation;

	private int nrOfSequences;
	private int nrOfMutants;
	private int numberofSingleMutations;
	private int numberOfDoubleMutations;
	private double[][] SingleMutations;
	private double[][][] DoubleMutations;
	private double[][] SingleMutationsMean;
	private double[][] SingleMutationsStd;
	private double[][][] DoubleMutationsMean;
	private double[][][] DoubleMutationsStd;

	public static void main(String []args){
		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}

		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);
		MutationalAnalysis.run(T);
	}
	public static void run(Hashtable<String,String> T){
		String program = Functions.getValue(T, "-p", "HELP").toUpperCase();
		String sequence = Functions.getValue(T, "-seq", "AUUUAGGGCUGAUUUAUUACUACACACAGCAGUGCAACAUCUGUCAGUACUUCUGGUGCUUCUAUUUUAGAGGCAGCUGUCAGGUGUGCGAUCAAUAAAAAAAGCGGGGUUUCAUCAUGUUUAAUGAAGUCCAUAGUAUUCAUGGUCAUACAUUAUUG");

		if(T.containsKey("-analysis")){
			String structure = Functions.getValue(T, "-db", "HELP");

			String f1 = Functions.getValue(T, "-f1", "HELP");
			String f2 = Functions.getValue(T, "-f2", "HELP");
			String d = Functions.getValue(T, "-d", "HELP");
			MutationalAnalysis A = new MutationalAnalysis(sequence, d+"/"+f1);
			MutationalAnalysis B = new MutationalAnalysis(sequence, d+"/"+f2);

			IntramolecularStructures IS = IntramolecularStructures.addDotBracketAnnotation(structure.toCharArray());
			A.Structure = IS;
			A.findOverrepresentedMutants(B);

			System.out.println("overrepresented Mutants");
			A.sortByRealtiveFrequency();
			A.printMutants(5,A.numberofSingleMutations,A.numberOfDoubleMutations,B.numberofSingleMutations, B.numberOfDoubleMutations);
			System.out.println();
			System.out.println();
			System.out.println("compensatory Mutations");

			A.sortByDeviatingFreq();
			A.printCompensatoryMutants(5,A.numberofSingleMutations,A.numberOfDoubleMutations,B.numberofSingleMutations, B.numberOfDoubleMutations);
			System.out.println();
			System.out.println();

			A.printCompensatoryStructureMutants(3,A.numberofSingleMutations,A.numberOfDoubleMutations,B.numberofSingleMutations, B.numberOfDoubleMutations);
			//		A.printCompensatoryMutants(5,A.numberofSingleMutations,A.numberOfDoubleMutations,B.numberofSingleMutations, B.numberOfDoubleMutations);
			System.out.println();
			System.out.println();



			A.printKnownCompensatoryStructureMutations(5,A.numberofSingleMutations,A.numberOfDoubleMutations,B.numberofSingleMutations, B.numberOfDoubleMutations, 0.5);
			System.out.println();
			System.out.println();


			if(T.containsKey("-o")){
				String outFile = d+"/"+f1+"_"+f2;
				A.printMutantsSingle(outFile, A.numberofSingleMutations,A.numberOfDoubleMutations,B.numberofSingleMutations, B.numberOfDoubleMutations);
				A.printMutants(outFile, A.numberofSingleMutations,A.numberOfDoubleMutations,B.numberofSingleMutations, B.numberOfDoubleMutations);
				A.printMutantsStructureCompensatory(outFile, A.numberofSingleMutations,A.numberOfDoubleMutations,B.numberofSingleMutations, B.numberOfDoubleMutations);
				A.printKnownCompensatoryStructureMutations(outFile, A.numberofSingleMutations,A.numberOfDoubleMutations,B.numberofSingleMutations, B.numberOfDoubleMutations, 0.5);
			}
		}
		if(T.containsKey("-print")){
			String structure = Functions.getValue(T, "-db", "HELP");
			String suffix = Functions.getValue(T, "-suffix", "HELP");
			String dir = Functions.getValue(T, "-d", "HELP");

			ArrayList <String> fileNames = IOTools.getSequenceFiles(dir,suffix);
			MutationalAnalysis[] A = new MutationalAnalysis[fileNames.size()];
			for(int i = 0; i < fileNames.size(); i++){
				System.out.println(fileNames.get(i));
				A[i] = new MutationalAnalysis(sequence, dir+"/"+fileNames.get(i));
			}

			print2Matrix(A,dir,RNAfunctions.RNAString2Int(sequence));
		}

		if(T.containsKey("-change")){
			String structure = Functions.getValue(T, "-db", "HELP");
			String singleChange =Functions.getValue(T, "-f1", "HELP");
			String doubleChange =Functions.getValue(T, "-f2", "HELP");
			String outFile =Functions.getValue(T, "-o", "HELP");
			int length = 4*sequence.length();
			MutationalAnalysis A = new MutationalAnalysis();
			A.singleMutation = new Change[length];
			A.doubleMutation = new Change[length][length];

			A.addSingleChange(singleChange);
			A.addDoubleChange(doubleChange);
			A.findReversedDoubleMutants(0.1,117);
	//		A.findReversedDoubleMutantsSimple(0.01,114);

			A.printFoldChange(outFile);
		}
	}

	public void printFoldChange(String outFile){
		try{

			ExtendedWriter EW= new ExtendedWriter(new FileWriter(outFile));

			EW.print("NA");
			for(int i = 0; i < this.singleMutation.length; i++){
				if(this.singleMutation[i] != null &&(this.singleMutation[i].log2FoldChange > 1 || this.singleMutation[i].log2FoldChange < -1)){
					int loc = i/4+1;
					String seq = RNAfunctions.DNAInt2String(i%4+1);
					EW.print("\t\""+(loc-117)+"_"+seq+"\"");
				}
			}
			EW.println();
			EW.print("0");
			for(int i = 0; i < this.singleMutation.length; i++){
				if(this.singleMutation[i] != null &&(this.singleMutation[i].log2FoldChange > 1 || this.singleMutation[i].log2FoldChange < -1)){
					if(this.singleMutation[i] != null)
						EW.print("\t"+this.singleMutation[i].log2FoldChange);
					else
						EW.print("\t0");
				}
			}
			EW.println();
			for(int i = 0; i < this.singleMutation.length; i++){
				if(this.singleMutation[i] != null &&(this.singleMutation[i].log2FoldChange > 1 || this.singleMutation[i].log2FoldChange < -1)){
					if(this.singleMutation[i] != null)
						EW.print(this.singleMutation[i].log2FoldChange);
					else
						EW.print("0");
					for(int j = 0; j < this.singleMutation.length;j++){
						if(this.singleMutation[j] != null &&(this.singleMutation[j].log2FoldChange > 1 || this.singleMutation[j].log2FoldChange < -1)){
							if(this.doubleMutation[i][j] != null)
								EW.print("\t"+this.doubleMutation[i][j].log2FoldChange);
							else
								EW.print("\t0");
						}
					}
					EW.println();
				}
			}
			EW.println();
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}


	public void findReversedDoubleMutants(double cutoff, int offset){
		for(int i = 0; i < this.singleMutation.length; i++){
			for(int j = i+1; j < this.singleMutation.length;j++){
				if(this.singleMutation[i] != null && this.singleMutation[j] != null && this.doubleMutation[i][j] != null){
					if((this.singleMutation[i].log2FoldChange<0 && this.singleMutation[i].pval < cutoff) 
							&& (this.singleMutation[j].log2FoldChange<0 && this.singleMutation[j].pval < cutoff)
							&& (this.doubleMutation[i][j].log2FoldChange> 0 && this.doubleMutation[i][j].pval < cutoff )){
						System.out.println("Conflicting pairs");
						this.singleMutation[i].printInfo(offset);
						this.singleMutation[j].printInfo(offset);
						this.doubleMutation[i][j].printInfo(offset);
						System.out.println();
					}
					if((this.singleMutation[i].log2FoldChange>0 && this.singleMutation[i].pval < cutoff) 
							&& (this.singleMutation[j].log2FoldChange>0 && this.singleMutation[j].pval < cutoff)
							&& (this.doubleMutation[i][j].log2FoldChange< 0 && this.doubleMutation[i][j].pval < cutoff )){
						System.out.println("Conflicting pairs");
						this.singleMutation[i].printInfo(offset);
						this.singleMutation[j].printInfo(offset);
						this.doubleMutation[i][j].printInfo(offset);
						System.out.println();
					}
				}
			}
		}
	}

	public void findReversedDoubleMutantsSimple(double cutoff, int offset){
		for(int i = 0; i < this.singleMutation.length; i++){
			for(int j = i+1; j < this.singleMutation.length;j++){
				if(this.singleMutation[i] != null && this.singleMutation[j] != null && this.doubleMutation[i][j] != null){
					if(((this.singleMutation[i].log2FoldChange<0 && this.singleMutation[i].pval <cutoff) 
							|| (this.singleMutation[j].log2FoldChange<0 && this.singleMutation[j].pval < cutoff))
							&& (this.doubleMutation[i][j].log2FoldChange> 0 && this.doubleMutation[i][j].pval < cutoff )){
						System.out.println("Conflicting pairs simple");
						this.singleMutation[i].printInfo(offset);
						this.singleMutation[j].printInfo(offset);
						this.doubleMutation[i][j].printInfo(offset);
						System.out.println();
					}
					if(((this.singleMutation[i].log2FoldChange>0 && this.singleMutation[i].pval < cutoff) 
							|| (this.singleMutation[j].log2FoldChange>0 && this.singleMutation[j].pval < cutoff))
							&& (this.doubleMutation[i][j].log2FoldChange < 0 && this.doubleMutation[i][j].pval < cutoff )){
						System.out.println("Conflicting pairs simple");
						this.singleMutation[i].printInfo(offset);
						this.singleMutation[j].printInfo(offset);
						this.doubleMutation[i][j].printInfo(offset);
						System.out.println();
					}
				}
			}
		}
	}


	public void addSingleChange(String infile){
		System.out.print("adding classes");
		try{
			ExtendedReader ER= new ExtendedReader(new FileReader(infile));
			ER.skipLine();
			while(ER.more()){
				String line = ER.readLine();
				Change temp = new Change();
				temp.addInfo(line);
				String[] location = temp.id.split("\"");
				String[] subLocation = location[1].split("\\ ");
				int loc = Integer.parseInt(subLocation[0]);
				int[] seq = RNAfunctions.RNAString2Int(subLocation[1]);
				this.singleMutation[loc*4+seq[0]-5] = temp;
			}
			ER.close();
		}catch(Exception E){
			E.printStackTrace();
		}
	}

	public void addDoubleChange(String infile){
		System.out.print("adding classes");
		try{
			ExtendedReader ER= new ExtendedReader(new FileReader(infile));
			ER.skipLine();
			while(ER.more()){
				String line = ER.readLine();
				Change temp = new Change();
				temp.addInfo(line);
				String[] location = temp.id.split("\"");
				String[] subLocation = location[1].split("\\ ");
				int loc = Integer.parseInt(subLocation[0]);
				int[] seq = RNAfunctions.RNAString2Int(subLocation[1]);
				int loc2 = Integer.parseInt(subLocation[2]);
				int[] seq2 = RNAfunctions.RNAString2Int(subLocation[3]);
				this.doubleMutation[loc*4+seq[0]-5][loc2*4+seq2[0]-5] = temp;
				this.doubleMutation[loc2*4+seq2[0]-5][loc*4+seq[0]-5] = temp;

			}
			ER.close();
		}catch(Exception E){
			E.printStackTrace();
		}
	}




	MutationalAnalysis(){

	}


	MutationalAnalysis(String Sequence, String MutationFile){

		Name = MutationFile;
		this.Sequence = RNAfunctions.RNAString2Int(Sequence);
		this.mutants = new ArrayList <Mutations>();
		ArrayList<Sequence> mutants =  getMutations(MutationFile);
		getMutations(mutants);

	}

	public static void print2Matrix(MutationalAnalysis [] all,String dir,int[] Sequence){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(dir+".singleMutations.txt"));
			EW.print("Mutation");
			for(int i = 0; i < all.length;i++){
				EW.print("\t"+all[i].Name);
			}
			EW.println();
			//this.SingleMutations[mutants.get(i).get(0).location][mutants.get(i).get(0).nucleotide]++;

			for(int i = 0; i < Sequence.length ; i++){
				for(int j = 1; j < 5;j++){
					if(j != Sequence[i]){
						EW.print((i+1)+" "+RNAfunctions.RNAInt2String(j));
						for(int k = 0; k < all.length; k++){
							EW.print("\t"+(int)all[k].SingleMutations[i][j]);
						}
						EW.println();
					}
				}
			}

			EW.flush();
			EW.close();
			EW = new ExtendedWriter(new FileWriter(dir+".DoubleMutations.txt"));
			EW.print("Mutation");
			for(int i = 0; i < all.length;i++){
				EW.print("\t"+all[i].Name);
			}
			EW.println();
			//this.DoubleMutations[mutants.get(i).get(0).location][mutants.get(i).get(1).location][mutants.get(i).get(0).nucleotide*5+mutants.get(i).get(1).nucleotide] ++;
			for(int i = 0; i < Sequence.length ; i++){
				for(int j = 1; j < 5;j++){
					if(j != Sequence[i]){
						for(int k = i+1 ; k < Sequence.length ; k++){
							for(int l = 1; l < 5;l++){
								if(l != Sequence[k]){
									EW.print((i+1)+" "+RNAfunctions.RNAInt2String(j)+" "+(k+1)+" "+RNAfunctions.RNAInt2String(l));
									for(int m = 0; m < all.length; m++){
										EW.print("\t"+(int)all[m].DoubleMutations[i][k][j*5+l]);
									}
									EW.println();
								}
							}
						}
					}
				}
			}
			EW.flush();
			EW.close();



		}catch(Exception E){
			E.printStackTrace();
		}


	}

	private ArrayList <Sequence> getMutations(String sequenceFile){
		ArrayList <Sequence> mutants = new ArrayList<Sequence>();
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(sequenceFile));
			int count = 100000;
			while(ER.more()){
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					nrOfSequences++;
					FastaSequence A = FastaSequences.getFasta(ER);
					Sequence temp = new Sequence(this.Sequence, A.getSequence());
					if(temp.hasMutations())
						mutants.add(temp);
				}
				if(nrOfSequences > count){
					System.out.println("nr Of Sequences Read" + count);
					count= count + 100000;
				}
			}
			ER.close();
		}catch(Exception E){E.printStackTrace();}
		return mutants;
	}




	private void getMutations(ArrayList <Sequence> mutants){
		SingleMutations = new double[Sequence.length][5];
		this.DoubleMutations = new double[Sequence.length][Sequence.length][25];
		for(int i = 0; i < mutants.size();i++){
			int nrOfMutations = mutants.get(i).getNrOfMutations();

			if(nrOfMutations ==1 && mutants.get(i).get(0).nucleotide < 5){
				this.SingleMutations[mutants.get(i).get(0).location][mutants.get(i).get(0).nucleotide]++;
				nrOfMutants++;
				this.numberofSingleMutations++;
			}
			else if(nrOfMutations == 2 && mutants.get(i).get(0).nucleotide < 5 && mutants.get(i).get(1).nucleotide < 5){
				this.DoubleMutations[mutants.get(i).get(0).location][mutants.get(i).get(1).location][mutants.get(i).get(0).nucleotide*5+mutants.get(i).get(1).nucleotide] ++;
				nrOfMutants++;
				this.numberOfDoubleMutations++;
			}
			else if(nrOfMutations > 2)
				nrOfMutants++;
		}
		System.out.println(this.Name);
		System.out.println(this.numberofSingleMutations);
		System.out.println(this.numberOfDoubleMutations);
		System.out.println();

	}



	@SuppressWarnings("unchecked")
	private void findOverrepresentedMutants(MutationalAnalysis B){
		System.out.println(this.numberofSingleMutations);
		for( int i = 0; i < SingleMutations.length; i++){
			for(int j = 1; j < SingleMutations[i].length; j++){
				Mutations temp = new Mutations(i,j);
				temp.addFreq(this.SingleMutations[i][j],this.numberofSingleMutations, B.SingleMutations[i][j], B.numberofSingleMutations);
				this.mutants.add(temp);
			}
		}
		int count = 0;
		for( int i = 0; i < DoubleMutations.length; i++){
			for(int j = i+1; j < DoubleMutations[i].length; j++){
				for(int k = 0; k < DoubleMutations[i][j].length;k++){
					if(k%5 != 0){
						int[] loc = new int[2];
						loc[0] = i;
						loc[1] = j;
						int[] seq = new int[2];
						seq[0] = k/5;
						seq[1] = k%5;
						count++;
						Mutations temp = new Mutations(loc,seq);
						temp.addFreq(this.DoubleMutations[i][j][k],this.numberOfDoubleMutations, B.DoubleMutations[i][j][k], B.numberOfDoubleMutations);
						temp.addSuspectedFreq(this.SingleMutations[i][k/5], this.SingleMutations[j][k%5], B.SingleMutations[i][k/5], B.SingleMutations[j][k%5],this.numberofSingleMutations,B.numberofSingleMutations);
						this.mutants.add(temp);
					}
				}
			}
			System.out.println(this.mutants.size());
			System.out.println(count);
		}
	}


	private void sortByDeviatingFreq(){
		Collections.sort(mutants,new byDeviatingFrequency());
	}

	private void sortByRealtiveFrequency(){
		Collections.sort(mutants,new byRelativeFrequency());
	}


	private void printMutants(double cutoff, double AS, double AD,double BS, double BD){
		for(int i = 0; i < this.mutants.size();i++){
			mutants.get(i).printMutations(this.Sequence, cutoff, AS,AD,BS,BD);
		}
	}

	private void printMutants(String outFile, double AS, double AD,double BS, double BD){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outFile+".txt"));
			for(int i = 0; i < this.mutants.size();i++){
				mutants.get(i).printMutations(EW, this.Sequence,  AS,AD,BS,BD);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}


	private void printMutantsSingle(String outFile, double AS, double AD,double BS, double BD){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outFile+".single.txt"));
			for(int i = 0; i < this.mutants.size();i++){
				mutants.get(i).printMutationsSingle(EW, this.Sequence,  AS,AD,BS,BD);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}

	private void printMutantsStructureCompensatory(String outFile, double AS, double AD,double BS, double BD){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outFile+".compensatory.structure.txt"));
			for(int i = 0; i < this.mutants.size();i++){
				mutants.get(i).printMutationsStructureCompensatory(EW, this.Sequence,  AS,AD,BS,BD);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}

	private void printKnownCompensatoryStructureMutations(String outFile, double AS, double AD,double BS, double BD, double cutoff){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outFile+".compensatory.knownstructure.txt"));
			for(int i = 0; i < this.mutants.size();i++){
				mutants.get(i).printKnownCompensatoryStructureMutations(EW,this.Structure, this.Sequence,  AS,AD,BS,BD,cutoff);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}


	private void printCompensatoryMutants(double cutoff, double AS, double AD,double BS, double BD){
		for(int i = 0; i < this.mutants.size();i++){
			mutants.get(i).printCompensatoryMutations(this.Sequence, cutoff, AS,AD,BS,BD);
		}
	}


	private void printCompensatoryStructureMutants(double cutoff, double AS, double AD,double BS, double BD){
		for(int i = 0; i < this.mutants.size();i++){
			mutants.get(i).printCompensatoryStructureMutations(this.Sequence, cutoff, AS,AD,BS,BD);
		}
	}

	private void printKnownCompensatoryStructureMutations(double cutoff1,  double AS, double AD,double BS, double BD, double cutoff){
		System.out.println(this.mutants.size());
		for(int i = 0; i < this.mutants.size();i++){
			mutants.get(i).printKnownCompensatoryStructureMutations(cutoff1,this.Structure, this.Sequence,  AS,AD,BS,BD,cutoff);
		}
	}


	public class byRelativeFrequency implements java.util.Comparator {
		public int compare(Object A, Object B) {
			double sdif = ((Mutations)A).relFreq - ((Mutations)B).relFreq;
			if(sdif > 0)
				return  1;
			else
				return -1;
		}
	}

	public class byDeviatingFrequency implements java.util.Comparator {
		public int compare(Object A, Object B) {
			double sdif = ((Mutations)A).deviatingFreq - ((Mutations)B).deviatingFreq;
			if(sdif > 0)
				return  1;
			else
				return -1;
		}
	}

}






