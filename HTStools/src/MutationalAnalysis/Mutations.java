package MutationalAnalysis;

import structure.IntramolecularStructures;
import general.ExtendedWriter;
import general.Functions;
import general.RNAfunctions;

public class Mutations {
	int[] location;
	int[] nucleotide;
	double[] freq;
	double relFreq;
	double suspectedFreq;
	double deviatingFreq;
	double compensatoryMutations;

	Mutations(int loc, int nuc){
		this.location = new int[1];
		this.nucleotide = new int[1];
		this.location[0] = loc;
		this.nucleotide[0] = nuc;
	}

	Mutations(int[] loc, int []nucleotide){
		this.location = loc;
		this.nucleotide = nucleotide;
	}
	public void addFreq(double Afreq,int AnrOfSeq, double Bfreq,int BnrOfSeq){
		this.freq = new double[2];
		this.freq[0] = Afreq/(double)AnrOfSeq;
		this.freq[1] = Bfreq/(double)BnrOfSeq;
		if(this.freq[0] > 1)
			System.out.println("Something is wrong");
		if(this.freq[1] > 1)
			System.out.println("Something is wrong");
		if(Bfreq != 0 && Afreq != 0)
			this.relFreq = Math.log(freq[0]/freq[1]);
		else if(Bfreq == 0 && Afreq == 0)
			this.relFreq = 0;
		else if(Bfreq == 0)
			this.relFreq = Math.log(freq[0]/(1/(double)BnrOfSeq));
		else
			this.relFreq = Math.log((1/(double)AnrOfSeq)/freq[1]);
	}

	public void addSuspectedFreq(double A1, double A2, double B1,double B2, double NA, double NB){

		if(A1 == 0 && B1 == 0 && A2 == 0 && B2 == 0)
			this.suspectedFreq = 0;
		else{
			if(A1 == 0) A1 = 1;
			if(A2 == 0) A2 = 1;
			if(B1 == 0) B1 = 1;
			if(B2 == 0) B2 = 1;
			double A = A1*A2*2/(NA*NA);
			double B = B1*B2*2/(NB*NB);
			double C = A/B;
			this.suspectedFreq = Math.log(C);
			this.deviatingFreq = this.relFreq-this.suspectedFreq;
		}
	}

	public int sortFreq(Mutations B){
		if(this.relFreq < B.relFreq ) return -1;
		else return 1;
	}

	public void printMutations(int[] seq, double cutoff, double AS, double AD,double BS, double BD){
		if(Math.abs( relFreq) > cutoff){
			printMutants2Out(seq, AS,AD,BS,BD);
		}
	}

	public void printMutations(ExtendedWriter EW, int[] seq, double AS, double AD,double BS, double BD){
		printMutants2File(EW,seq, AS,AD,BS,BD);

	}

	
	public void printMutationsSingle(ExtendedWriter EW, int[] seq, double AS, double AD,double BS, double BD){
		if(this.location.length ==1){
			printMutants2File(EW,seq, AS,AD,BS,BD);

		}
	}

	public void printCompensatoryMutations(int[] seq, double cutoff, double AS, double AD,double BS, double BD){
		if(Math.abs(this.deviatingFreq) > cutoff && Math.abs( relFreq) > cutoff && this.location.length ==2){
			printMutants2Out(seq, AS,AD,BS,BD);
		}
	}
	
	
	public void printMutants2File(ExtendedWriter EW,int[] seq, double AS, double AD,double BS, double BD){
		EW.print(Functions.fixedLength(relFreq, 5)+"\t");
		if(this.location.length ==2){
			EW.print(Functions.fixedLength(suspectedFreq, 5)+"\t");
			EW.print(Functions.fixedLength(this.deviatingFreq, 5)+"\t");
			EW.print(Functions.fixedLength((Math.log(freq[0]*AD+1)/ Math.log(2)), 10)+"\t");
			EW.print(Functions.fixedLength((Math.log(freq[1]*BD+1)/ Math.log(2)), 10)+"\t");
		}
		else{
			EW.print(Functions.fixedLength("-", 5)+"\t");
			EW.print(Functions.fixedLength("-", 5)+"\t");
			EW.print(Functions.fixedLength((Math.log(freq[0]*AS+1)/ Math.log(2)), 25)+"\t");
			EW.print(Functions.fixedLength((Math.log(freq[1]*BS+1)/ Math.log(2)), 25)+"\t");
		}
		for(int i  = 0; i < location.length; i++){
			EW.print(location[i]+"\t"+RNAfunctions.DNAInt2String(seq[location[i]])+"->"+RNAfunctions.DNAInt2String(this.nucleotide[i])+"\t");
		}
		EW.println();
	}
	
	
	public void printMutants2Out(int[] seq, double AS, double AD,double BS, double BD){
		System.out.print(Functions.fixedLength(relFreq, 5)+"\t");
		if(this.location.length ==2){
			System.out.print(Functions.fixedLength(suspectedFreq, 5)+"\t");
			System.out.print(Functions.fixedLength(this.deviatingFreq, 5)+"\t");
			System.out.print(Functions.fixedLength((Math.log(freq[0]*AD+1)/ Math.log(2)), 10)+"\t");
			System.out.print(Functions.fixedLength((Math.log(freq[1]*BD+1) / Math.log(2)), 10)+"\t");
		}
		else{
			System.out.print(Functions.fixedLength("-", 5)+"\t");
			System.out.print(Functions.fixedLength("-", 5)+"\t");
			System.out.print(Functions.fixedLength((Math.log(freq[0]*AS+1)/ Math.log(2)), 25)+"\t");
			System.out.print(Functions.fixedLength((Math.log(freq[1]*BS+1)/ Math.log(2)), 25)+"\t");
		}
		for(int i  = 0; i < location.length; i++){
			System.out.print(location[i]+"\t"+RNAfunctions.DNAInt2String(seq[location[i]])+"->"+RNAfunctions.DNAInt2String(this.nucleotide[i])+"\t");
		}
		System.out.println();
	}
	
	public void printCompensatoryStructureMutations(int[] seq, double cutoff, double AS, double AD,double BS, double BD){
		if(Math.abs(this.deviatingFreq) > cutoff && Math.abs( relFreq) > cutoff && this.location.length ==2){
			if(RNAfunctions.isBasepair(seq[location[0]],seq[location[1]]) && RNAfunctions.isBasepair(this.nucleotide[0], this.nucleotide[1])){
				printMutants2Out(seq, AS,AD,BS,BD);
			}
		}
	}
	
	public void printMutationsStructureCompensatory(ExtendedWriter EW, int[] seq, double AS, double AD,double BS, double BD){
		if(this.location.length ==2 && RNAfunctions.isBasepair(seq[location[0]],seq[location[1]]) && RNAfunctions.isBasepair(this.nucleotide[0], this.nucleotide[1])){
				printMutants2File(EW,seq, AS,AD,BS,BD);
		}
	}
	
	public void printKnownCompensatoryStructureMutations(ExtendedWriter EW, IntramolecularStructures IS, int[] seq, double AS, double AD,double BS, double BD, double cutoff){
		if(this.location.length ==2 && IS.isBasepair(location[0], location[1], cutoff)){
				printMutants2File(EW,seq, AS,AD,BS,BD);

		}
	}
	
	
	
	public void printKnownCompensatoryStructureMutations(double cutoff1, IntramolecularStructures IS, int[] seq, double AS, double AD,double BS, double BD, double cutoff){
		if( this.location.length ==2 && IS.isBasepair(location[0], location[1], cutoff)){
			printMutants2Out(seq, AS,AD,BS,BD);
		}
	}
	
}
