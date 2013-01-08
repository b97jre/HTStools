package general;



import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Random;
//import javax.mail.*;
//import javax.mail.internet.*;
//import java.util.*;


import dbsr.jprogram.LocalProgram;
import general.ExtendedReader;
import general.ExtendedWriter;


public class RNAfunctions
{
	//////////Pattern Search Specific variables//////////////////////



	public static double[][] calculateDistances(char[][] alignSequences){
		int[][] Sequences = new int[alignSequences.length][0];
		for(int i = 0; i < alignSequences.length; i++){
			Sequences[i] = RNAchar2int(alignSequences[i]);
		}
		return calculateDistances(Sequences);		
	}

	public static double[][] calculateDistances(int[][] alignSequences){
		double [][] distance = new double[alignSequences.length][alignSequences.length];
		for(int j = 0; j < alignSequences.length; j++){
			for(int i = j+1; i < alignSequences.length; i++){

				distance[j][i] = getDistanceKimura(alignSequences[j],alignSequences[i]);
				distance[i][j] = distance[j][i];
			}
		}
		return distance;
	}

	public static double calculateDistanceKimura(char[] Seq1, char[] Seq2){
		int[] Sequences1 = RNAchar2int(Seq1);
		int[] Sequences2 = RNAchar2int(Seq2);
		return getDistanceKimura(Sequences1,Sequences2);
	}

	public static double calculateDistance(char[] RNA1,char[] RNA2){
		if(RNA1.length < 1 ||RNA2.length < 1) return 1;
		int nrOfNucleotides = 0;
		int nrOfMatchingNucleotides = 0;
		for(int i = 0; i < RNA1.length; i++){
			if(RNA1[i] != '-' && RNA1[i] != 'N' && RNA2[i] != '-' && RNA2[i] != 'N'){
				if(RNA1[i] == RNA2[i]){
					nrOfNucleotides++;
					nrOfMatchingNucleotides++;
				}
				else{
					nrOfNucleotides++;
				}

			}
		}
		return (double)(1.0-((double)nrOfMatchingNucleotides/(double)nrOfNucleotides));
	}		


	public static boolean theSame(int []A, int[] B){
		if(A.length != B.length) return false;
		for(int j = 0; j < A.length;j++){
			if(A[j] != B[j]) return false;
		}
		return true;

	}

	public static void printDistance(double[][] distance){
		for(int i = 0; i < distance.length; i++){
			for(int j = 0; j < distance[i].length;j++){
				System.out.print(distance[i][j]+"\t");
			}
			System.out.println();
		}
		System.out.println();
		System.out.println();
		System.out.println();
		System.out.println();

	}

	public static int countSequence(String A){
		char[] B = A.toCharArray();
		int Length = B.length;
		int count = 0;
		for(int i = 0; i < Length ; i++){
			if(B[i] != '-')count++;
		}
		return count;
	}

	public static String reverse(String A){

		char[] B = A.toCharArray();
		int Length = B.length;
		char[] C = new char[Length];
		for(int i = 0; i < Length ; i++){
			C[i] = B[Length-1-i];
		}
		return new String(C);
	}

	public static int[] reverse(int[] A){
		int Length = A.length;
		int[] C = new int[Length];
		for(int i = 0; i < Length ; i++){
			C[i] = A[Length-1-i];
		}
		return  C;
	}

	public static char[] reverse(char[] A){
		int Length = A.length;
		char[] C = new char[Length];
		for(int i = 0; i < Length ; i++){
			C[i] = A[Length-1-i];
		}
		return  C;
	}


	public static int getComplementary(int A){
		if(A == 1) return 4;
		else if(A == 2) return 3;
		else if(A == 3)return 2;
		else if(A == 4)return 1;
		else return A;
	}

	public static int[] getComplementary(int[] A){
		A = reverse(A);
		for(int i = 0; i < A.length; i ++){
			A[i] = getComplementary(A[i]);
		}
		return A;
	}

	public static int[] getReverseComplement(int[] A){
		A = reverse(A);
		for(int i = 0; i < A.length; i ++){
			A[i] = getComplementary(A[i]);
		}
		return A;
	}


	public static String getReverseComplement(String A){
		return RNAInt2String(getReverseComplement(RNAString2Int(A)));
	}

	public static int[] getAntisense(int[] A){
		int [] antisenseA = new int[A.length];
		for(int i = 0; i < A.length; i ++){
			antisenseA[i] = getComplementary(A[i]);
		}
		return antisenseA;
	}

	public static boolean isBasepair(char A, char B){
		return isBasepair(RNAchar2int2(A),RNAchar2int2(B));
	}

	public static boolean isBasepair(int A, int B){
		if(A == 1){
			if(B == 4)
				return true;
			return false;
		} 
		if(A == 2){
			if(B == 3)
				return true;
			return false;
		}
		else if(A == 3){
			if(B == 2)
				return true;
			if(B == 4)
				return true;
			return false;
		}
		else if(A == 4){
			if(B == 3)
				return true;
			if(B == 1)
				return true;
			return false;
		}
		return false;
	}


	public static boolean[][] getCanonicalBasepairsMM(){
		boolean[][] MM = new boolean[6][6];
		for(int i = 0; i< 6 ;i++){
			for(int j = 0; j< 6 ;j++){
				MM[i][j] = false;
			}
		}
		MM[1][4] = true;
		MM[2][3] = true;
		MM[3][2] = true;
		MM[4][1] = true;
		return MM;
	}

	public static boolean[][] getBasepairsMM(){
		boolean[][] MM = new boolean[6][6];
		for(int i = 0; i< 6 ;i++){
			for(int j = 0; j< 6 ;j++){
				MM[i][j] = false;
			}
		}
		MM[1][4] = true;
		MM[2][3] = true;
		MM[3][2] = true;
		MM[4][1] = true;
		MM[4][3] = true;
		MM[3][4] = true;
		return MM;
	}


	public static char getComplementary(char A){
		if(A == 'A' || A == 'a') return 'U';
		else if(A == 'C' || A == 'c')return 'G';
		else if(A == 'G' || A == 'g')return 'C';
		else if(A == 'T' || A == 't' || A == 'U' || A == 'u')return 'A';
		else return A;
	}

	public static char[] getComplementary(char[] A){
		A = reverse(A);
		for(int i = 0; i < A.length; i ++){
			A[i] = getComplementary(A[i]);
		}
		return A;
	}


	public static char[][] getComplementary(char[][] A){
		for(int i = 0; i < A.length; i++){
			if(A[i] != null)
				A[i] = getComplementary(A[i]);

		}
		return A;
	}


	public static int[][] splitSequence(int[] A){
		int [][] Sequences = new int [2][0];

		int Length = A.length;

		int[] C = new int[Length-40];
		for(int i = 0; i < Length-40 ; i++){
			C[i] = A[i];
		}
		int[] D = new int[20];
		for(int i = 0; i < 20; i++){
			D[i] = A[Length-6-i];
		}
		Sequences[0] = C;
		Sequences[1] = D;

		return  Sequences;
	}

	public static double getDistanceKimura(int[] Seq1 , int[] Seq2){
		int[][] SubstitutionMatrix = new int[4][4];
		int totalNrOfCompared = 0;
		for(int i = 0; i < Seq1.length; i++){
			if(Seq1[i] > 0 && Seq1[i] < 5 && Seq2[i] > 0 && Seq2[i] < 5){
				SubstitutionMatrix[Seq1[i]-1][Seq2[i]-1]++;
				totalNrOfCompared++;
			}
		}

		if(totalNrOfCompared < 30)
			return 1;


		int nrOfTransitions =   SubstitutionMatrix[2][0] + SubstitutionMatrix[0][2] + SubstitutionMatrix[3][1] + SubstitutionMatrix[1][3];
		int nrOfTransversions = SubstitutionMatrix[1][0] + SubstitutionMatrix[3][0] + SubstitutionMatrix[0][1] + SubstitutionMatrix[0][3] + 
		SubstitutionMatrix[2][1] + SubstitutionMatrix[1][2] + SubstitutionMatrix[2][3] + SubstitutionMatrix[3][2];

		double P = (double)nrOfTransitions/(double)totalNrOfCompared;
		double Q = (double)nrOfTransversions/(double)totalNrOfCompared;

		double distance = 1.0/2.0*java.lang.Math.log(1.0/(1.0-2.0*P-Q))+1.0/4.0*java.lang.Math.log(1.0/(1.0-2.0*Q));
		return distance;
	}

	public static String[] removeGapsInSpecificSequence(String[] Sequences, int position){
		char[][] tempSequences = new char[Sequences.length][0];
		for(int i = 0;i < Sequences.length;i++){
			tempSequences[i] = Sequences[i].toCharArray();
		}
		int count = 0;
		for(int i = 0 ; i < tempSequences[position].length ; i++){
			if(tempSequences[position][i] =='-')
				count++;
		}
		char[][] newSequences = new char[tempSequences.length][tempSequences[position].length-count];
		count = 0;
		for(int i = 0 ; i < tempSequences[position].length ; i++){
			if(tempSequences[position][i] =='-'){
				count++;	
			}
			else{
				for(int j = 0;j < Sequences.length;j++){
					newSequences[j][i-count] = tempSequences[j][i];
				}

			}

		}
		for(int j = 0;j < Sequences.length;j++){
			Sequences[j] = new String(newSequences[j]);
		}
		return Sequences;
	}	


	public static String removeGaps(String gapedSequence){
		return new String(removeGaps(gapedSequence.toCharArray()));
	}

	public static char[] removeGaps(char[] gapedSequence){
		int count = 0;
		for(int i = 0; i <  gapedSequence.length; i++){
			if(gapedSequence[i] != '-'){
				//				gapedSequence[count] = gapedSequence[i];
				count++;
			}
		}
		char[] Sequence = new char[count];
		count = 0;
		for(int i = 0; i <  gapedSequence.length; i++){
			if(gapedSequence[i] != '-'){
				Sequence[count] = gapedSequence[i];
				count++;
			}
		}
		return Sequence;
	}

	public static int[] removeGaps(int[] gapedSequence){
		int[] Sequence = null;
		for(int i = 0; i <  gapedSequence.length; i++){
			if(gapedSequence[i] != 0){
				Sequence = Functions.addInt(Sequence, gapedSequence[i]);
			}
		}
		return Sequence;
	}

	public static int[] mapGaps(int[] ungapedSequence,int[] gapedSequence ){
		if(ungapedSequence.length != gapedSequence.length)
			System.out.println("The sequences presented in mapGaps in RNAfunctions does not have the same length \n "+ ungapedSequence.length+" "+gapedSequence.length);
		for(int i = 0; i < ungapedSequence.length;i++){
			if(gapedSequence[i] == 0)
				ungapedSequence[i] = 0;
		}
		return ungapedSequence;
	}



	public static int[][] addIntSequence(int[][] SequenceArray, int[] newSequence){
		if(SequenceArray == null){
			SequenceArray = new int[1][0];
			SequenceArray[0] = newSequence;
			return SequenceArray;
		}

		int Length = SequenceArray.length;;
		int [][] newSequenceArray = new int[Length +1][0];
		for(int i = 0; i < Length;i++){
			newSequenceArray[i] = SequenceArray[i];
		}
		newSequenceArray[Length] = newSequence;
		return newSequenceArray;
	}

	public static int[] addIntSequence(int[] SequenceArray, int[] newSequence){
		if(SequenceArray == null){
			return newSequence;
		}

		int Length = SequenceArray.length;;
		int [] newSequenceArray = new int[Length+ newSequence.length];
		for(int i = 0; i < Length;i++){
			newSequenceArray[i] = SequenceArray[i];
		}
		for(int i = 0; i < newSequence.length;i++){
			newSequenceArray[Length+i] = newSequence[i];
		}
		return newSequenceArray;
	}



	public static int[] RNAString2Int(String A){
		char[] B = A.toCharArray();
		return RNAchar2int(B);
	}

	public static int[] RNAString2Int2(String A){
		char[] B = A.toCharArray();
		return RNAchar2int2(B);
	}

	public static int[] RNAchar2int2(char[] B){
		int[] IntArray = new int[B.length];
		for(int i = 0; i < B.length;i++){
			IntArray [i] = RNAchar2int2(B[i]);
		}
		return IntArray;
	}

	public static int[] RNAchar2int(char[] B){
		int[] IntArray = new int[B.length];
		for(int i = 0; i < B.length;i++){
			IntArray [i] = RNAchar2int(B[i]);
		}
		return IntArray;
	}

	public static int RNAchar2int(char B){
		if(B == 'A' ||B ==  'a'){return 1;}
		if(B == 'C'||B ==  'c'){return 2 ;}
		if(B == 'G'||B ==  'g'){return 3 ;}
		if(B == 'T'||B == 'U' || B == 't'|| B == 'u'){return 4;}
		// Ambigous code
		if(B == '-'){return 0;}
		if(B == 'N'||B ==  'n'){return 5 ;}//A+C+G+T

		if(B == 'V'||B ==  'v'){return 14 ;}//A+C+G
		if(B == 'H'||B ==  'h'){return 6 ;}//A+C+T
		if(B == 'D'||B ==  'd'){return 7 ;}//A+G+T

		if(B == 'M'||B ==  'm'){return 8 ;}//A+C
		if(B == 'R'||B ==  'r'){return 9 ;}//A+G
		if(B == 'W'||B ==  'w'){return 10;}//A+T
		if(B == 'S'||B ==  's'){return 11;}//C+G
		if(B == 'Y'||B ==  'y'){return 12;}//C+T
		if(B == 'K'||B ==  'k'){return 13;}//G+T

		return 5;
	}

	public static int RNAchar2int2(char B){
		if(B == 'A' ||B ==  'a'){return 1;}
		if(B == 'C'||B ==  'c'){return 2 ;}
		if(B == 'G'||B ==  'g'){return 3 ;}
		if(B == 'T'||B == 'U' || B == 't'|| B == 'u'){return 4;}
		// Ambigous code
		if(B == '-'){return 0;}
		if(B == 'N'||B ==  'n'){return 5 ;}//A+C+G+T

		return 5;
	}

	public static String RNAInt2String(int[] B){
		return new String(RNAInt2char(B));
	}

	public static String DNAInt2String(int[] B){
		return new String(DNAInt2char(B));
	}


	public static char[] RNAInt2char(int[] B){
		int Length = B.length;
		char[] charArray = new char[Length];
		for(int i = 0; i < Length;i++){
			charArray[i] = RNAInt2char(B[i]);
		}
		return charArray;
	}

	public static char[] DNAInt2char(int[] B){
		int Length = B.length;
		char[] charArray = new char[Length];
		for(int i = 0; i < Length;i++){
			charArray[i] = DNAInt2char(B[i]);
		}
		return charArray;
	}


	public static char RNAInt2char(int B){
		if(B == 1)return 'A';
		else if(B == 2)return 'C';
		else if(B == 3)return 'G';
		else if(B == 4)return 'U';

		else if(B == 5)return 'N' ;//A+C+G+T
		else if(B == 14)return 'V' ;//A+C+G
		else if(B == 6)return 'H' ;//A+C+T
		else if(B == 7)return 'D' ;//A+G+T

		else if(B == 8)return 'M' ;//A+C
		else if(B == 9)return 'R' ;//A+G
		else if(B == 10)return 'W';//A+T
		else if(B == 11)return 'S';//C+G
		else if(B == 12)return 'Y';//C+T
		else if(B == 13)return 'K';//G+T
		else if(B == 0)return '-';
		return 'N';

	}


	public static char DNAInt2char(int B){
		if(B == 1)return 'A';
		else if(B == 2)return 'C';
		else if(B == 3)return 'G';
		else if(B == 4)return 'T';

		else if(B == 5)return 'N' ;//A+C+G+T
		else if(B == 14)return 'V' ;//A+C+G
		else if(B == 6)return 'H' ;//A+C+T
		else if(B == 7)return 'D' ;//A+G+T

		else if(B == 8)return 'M' ;//A+C
		else if(B == 9)return 'R' ;//A+G
		else if(B == 10)return 'W';//A+T
		else if(B == 11)return 'S';//C+G
		else if(B == 12)return 'Y';//C+T
		else if(B == 13)return 'K';//G+T
		else if(B == 0)return '-';
		return 'N';

	}

	public static int countRestrictionSites(int[] primer, String[] restrictionEnzymeSequences){
		int count = 0;
		for(int i = 0; i < restrictionEnzymeSequences.length;i++ ){
			int[] RES = RNAString2Int(restrictionEnzymeSequences[i]);
			for(int j = 0; j < primer.length - RES.length; j++ ){
				boolean same = true;
				for(int k = 0; k < RES.length; k++){
					if(RES[k] != primer[j+k]){
						same = false;
						k = RES.length;
					}
				}
				if(same) count++;
			}
		}
		return count;
	}



	public static int[] printAnimalmiRNAprimer(String miRNASequence, String restrictionEnzymeSequence1, String restrictionEnzymeSequence2 ,double Afreq,double Cfreq, double Gfreq, double Ufreq, int totalLength){
		double[] freq = new double[] {Afreq,Cfreq,Gfreq,Ufreq};
		int[] miRNAsequence = RNAString2Int(miRNASequence);
		int[] RES1 = RNAString2Int(restrictionEnzymeSequence1);
		int[] RES2 = RNAString2Int(restrictionEnzymeSequence2);
		int[] primer = new int[totalLength];
		int[][] matrix = getAntisenseMatrix(1, 1, -1);
		int count = 0;
		while(count < 3){
			primer[count] = getRandomNucleotide(freq);
			count++;
		}
		int count2 = 0;
		while(count2 < RES1.length){
			primer[count] = RES1[count2];
			count++;
			count2++;
		}
		while(count < (totalLength-miRNAsequence.length -  RES1.length - RES2.length - 6)/2 +RES1.length +3 ){
			primer[count] = getRandomNucleotide(freq);
			count++;
		}
		//NNNNIIIIIIEEEIIIIIIIA
		count2 = miRNAsequence.length-1;
		while(count2 > -1){
			while(count2 > 16){
				primer[count] = getRandomNucleotide(freq);
				count++;
				count2--;
			}
			while(count2 > 10){
				primer[count] = getComplementary(miRNAsequence[count2]);
				count++;
				count2--;
			}
			while(count2 > 7){
				int seq = 0;
				while(seq == 0 || matrix[seq][miRNAsequence[count2]] > 0)
					seq = getRandomNucleotide(freq);
				primer[count] = seq;
				count++;
				count2--;
			}
			while(count2 > 0){
				primer[count] = getComplementary(miRNAsequence[count2]);
				count++;
				count2--;
			}
			primer[count] = 1;
			count++;
			count2--;

		}
		while(count < totalLength-RES2.length - 3 ){
			primer[count] = getRandomNucleotide(freq);
			count++;
		}
		count2 = 0;
		while(count2 < RES2.length){
			primer[count] = RES2[count2];
			count++;
			count2++;
		}
		while(count < primer.length){
			primer[count] = getRandomNucleotide(freq);
			count++;
		}
		//System.out.println();
		return primer;
	}


	public static int[] printPerfectmiRNAprimer(String miRNASequence, String restrictionEnzymeSequence1, String restrictionEnzymeSequence2 ,double Afreq,double Cfreq, double Gfreq, double Ufreq, int totalLength){
		double[] freq = new double[] {Afreq,Cfreq,Gfreq,Ufreq};
		int[] miRNAsequence = RNAString2Int(miRNASequence);
		int[] RES1 = RNAString2Int(restrictionEnzymeSequence1);
		int[] RES2 = RNAString2Int(restrictionEnzymeSequence2);
		int[] primer = new int[totalLength];
		int[][] matrix = getAntisenseMatrix(1, 1, -1);
		int count = 0;
		while(count < 3){
			primer[count] = getRandomNucleotide(freq);
			count++;
		}
		int count2 = 0;
		while(count2 < RES1.length){
			primer[count] = RES1[count2];
			count++;
			count2++;
		}
		while(count < (totalLength-miRNAsequence.length -  RES1.length - RES2.length - 6)/2 +RES1.length +3 ){
			primer[count] = getRandomNucleotide(freq);
			count++;
		}
		//NNNNIIIIIIEEEIIIIIIIA
		count2 = miRNAsequence.length-1;
		while(count2 > -1){
			primer[count] = getComplementary(miRNAsequence[count2]);
			count++;
			count2--;
		}
		while(count < totalLength-RES2.length - 3 ){
			primer[count] = getRandomNucleotide(freq);
			count++;
		}
		count2 = 0;
		while(count2 < RES2.length){
			primer[count] = RES2[count2];
			count++;
			count2++;
		}
		while(count < primer.length){
			primer[count] = getRandomNucleotide(freq);
			count++;
		}
		//System.out.println();
		return primer;
	}


	public static int getRandomNucleotide(double[] freqs) {
		double[][] kmers;
		kmers = getConditionalKmers(freqs,0);
		double[] cumFreqs;
		cumFreqs = cumsum(kmers[0]);
		double r = Math.random();
		//System.out.print(r+",");
		int j = 0;
		while (j < 4 && r >= cumFreqs[j])
			j++;
		return Math.min((j+1),4);
	}





	public static String RNAInt2String(int B){
		if(B == 1)return "A";
		else if(B == 2)return "C";
		else if(B == 3)return "G";
		else if(B == 4)return "U";

		else if(B == 5)return "N" ;//A+C+G+T
		else if(B == 14)return "V" ;//A+C+G
		else if(B == 6)return "H" ;//A+C+T
		else if(B == 7)return "D" ;//A+G+T

		else if(B == 8)return "M" ;//A+C
		else if(B == 9)return "R" ;//A+G
		else if(B == 10)return "W";//A+T
		else if(B == 11)return "S";//C+G
		else if(B == 12)return "Y";//C+T
		else if(B == 13)return "K";//G+T
		else if(B == 0)return "-";
		return "N";

	}

	public static String DNAInt2String(int B){
		if(B == 1)return "A";
		else if(B == 2)return "C";
		else if(B == 3)return "G";
		else if(B == 4)return "T";

		else if(B == 5)return "N" ;//A+C+G+T
		else if(B == 14)return "V" ;//A+C+G
		else if(B == 6)return "H" ;//A+C+T
		else if(B == 7)return "D" ;//A+G+T

		else if(B == 8)return "M" ;//A+C
		else if(B == 9)return "R" ;//A+G
		else if(B == 10)return "W";//A+T
		else if(B == 11)return "S";//C+G
		else if(B == 12)return "Y";//C+T
		else if(B == 13)return "K";//G+T
		else if(B == 0)return "-";
		return "N";

	}


	public static boolean generateAlignFile(String dir , String FileName,String Extension){

		int nrOfGenomes = getNrOfGenomesFasta(dir,FileName,Extension);
		if(nrOfGenomes == 1)
		{
			return convert2AlnFile(dir,FileName,Extension);
		}
		else
			return AlignSequences(dir,FileName,Extension);
	}

	public static boolean AlignSequences(String dir , String FileName,String Extension){
		return MAFFTAlign(dir,FileName,Extension); 
		//return muscleAlign(dir,FileName,Extension);
	}


	public static int[] RNA2CS(int[] RNAsequence){
		int[][] conversationMatrix = new  int[][]{
				//	  - ,A,	C, G, U, V, H, D, M, R, W, S, Y, K, N
				{-1, 1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// - 
				{-1, 0, 1, 2, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// A
				{-1, 1, 0, 3, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// C
				{-1, 2, 3, 0, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// G
				{-1, 3, 2, 1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// U
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// V
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// H
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// D
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// M
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// R
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// W
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// S
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// Y
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// K
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// N
		};
		int[] cs = new int[RNAsequence.length-1];
		for(int i = 1; i < RNAsequence.length; i++){
			cs[i-1] = conversationMatrix[RNAsequence[i-1]][RNAsequence[i]];
		}
		return cs;
	}

	public static int[][] fasta2CFasta(String[] FastaSequences){
		int[][] CFS = new int[FastaSequences.length][];
		for(int i = 0; i < FastaSequences.length; i++)
			CFS[i] = Fasta2CFasta(FastaSequences[i]);
		return CFS;

	}


	public static int[] Fasta2CFasta(String Fasta){
		System.out.println("Fasta sequence : "+Fasta);
		int[] RNAsequence = RNAString2Int(Fasta);
		int[] CS = RNA2CS(RNAsequence);
		System.out.print("CFasta sequence: ");
		System.out.print(RNAInt2char(RNAsequence[0]));
		for(int i = 0; i <CS.length; i++){
			System.out.print(CS[i]);
		}
		System.out.println();
		System.out.println();
		return CS;
	}







	public static int[] CS2RNA(int[] ColorSequence, int startCodon){
		int[][] conversationMatrix = new  int[][]{
				//	  - ,A,	C, G, U, V, H, D, M, R, W, S, Y, K, N
				{ -1,1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// - 
				{-1,0,1,2,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// -
				{-1,1,0,3,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// A
				{-1,2,3,0,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// C
				{-1,3,2,1,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// G
				{ -1,1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// U 
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// V
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// H
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// D
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// M
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// R
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// W
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// S
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// Y
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// K
				{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},// N
		};
		int[] RNAsequence = new int[ColorSequence.length+1];
		RNAsequence[0] = startCodon;
		for(int i = 1; i < RNAsequence.length; i++){
			int j = 1;
			while(RNAsequence[i]<1){
				int coloSequence = ColorSequence[i-1];
				int temp = conversationMatrix[RNAsequence[i-1]][j];
				if(temp == coloSequence)RNAsequence[i] = j;
				j++;
			}
		}
		return RNAsequence;
	}




	public static boolean clustalAlign(String dir , String FileName){
		try{
			StringBuffer stdout = new StringBuffer();
			StringBuffer stderr = new StringBuffer();

			LocalProgram prog = new LocalProgram("clustalw");
			prog.loadParameters();
			prog.loadResources();

			if(dir.length() > 0)
				prog.setParameter("infile",dir+"/"+ FileName+".fasta");
			else
				prog.setParameter("infile", FileName+".fasta");

			prog.setParameter("outorder","input");
			prog.exec(null, stdout, stderr);
			prog.clearResources();
			prog = null;


			//System.out.println(stdout);
			//System.out.println(stderr);
		}
		catch(java.io.IOException IE){
			IE.printStackTrace();
			return false;
		}
		catch(Exception E){
			E.printStackTrace();
			return false;
		}
		return true;
	}

	public static boolean muscleAlign(String dir , String FileName, String extension){
		try{
			StringBuffer stdout = new StringBuffer();
			StringBuffer stderr = new StringBuffer();

			LocalProgram prog = new LocalProgram("muscle");
			prog.loadParameters();
			prog.loadResources();

			if(dir.length() < 1)
				dir = ".";
			prog.setParameter("in",dir+"/"+ FileName+"."+extension);
			prog.setParameter("out",dir+"/"+ FileName+".aln");
			prog.setParameter("stable");
			prog.setParameter("clwstrict");
			prog.setParameter("quiet");
			prog.setParameter("seqtype","rna");
			prog.exec(null, stdout, stderr);
			prog.clearResources();
			prog = null;
			//System.out.println(stdout);
			//System.out.println(stderr);
		}
		catch(java.io.IOException IE){
			IE.printStackTrace();
			return false;
		}
		catch(Exception E){
			E.printStackTrace();
			return false;
		}
		return true;
	}



	public static boolean MAFFTAlign(String dir , String FileName, String extension){
		try{
			System.out.println();
			StringBuffer stdout = new StringBuffer();
			StringBuffer stderr = new StringBuffer();

			LocalProgram prog = new LocalProgram("MAFFT");
			prog.loadParameters();
			prog.loadResources();

			if(dir.length() < 1)
				dir = ".";
			prog.setParameter("input",dir+"/"+ FileName+"."+extension);
			prog.setParameter("kimura", "1");
			prog.setParameter("ep","0.0");
			prog.setParameter("quiet");
			prog.setParameter("clustalout");
			String[] cmdArray = prog.buildCmdArray();

			prog.exec(null, stdout, stderr);
			prog.clearResources();
			prog = null;
			ExtendedWriter B = null;
			if(dir.length() > 0)
				B = new ExtendedWriter(new FileWriter(dir+"/"+FileName + ".aln"));
			else
				B = new ExtendedWriter(new FileWriter(FileName + ".aln"));

			B.println(stdout);
			B.flush();
			B.close();
			//			System.out.println(stdout);
			//			System.out.println(stderr);
		}
		catch(java.io.IOException IE){
			IE.printStackTrace();
			return false;
		}
		catch(Exception E){
			E.printStackTrace();
			return false;
		}
		return true;
	}


	public static boolean alignGenomes(String dir , String FocalGenomeFile, String[] comparativeGenomesFiles, String outputName){
		try{

			System.out.println();
			StringBuffer stdout = new StringBuffer();
			StringBuffer stderr = new StringBuffer();

			LocalProgram prog = new LocalProgram("MAUVE");
			prog.loadParameters();
			prog.loadResources();

			String input = dir+"/"+FocalGenomeFile;
			for(int i = 0; i < comparativeGenomesFiles.length; i++){
				input = input+" "+dir+"/"+comparativeGenomesFiles[i];
			}

			if(dir.length() < 1)
				dir = ".";
			prog.setParameter("output-guide-tree",dir+"/"+ outputName+".treeFile");
			prog.setParameter("backbone-output",dir+"/"+ outputName+".backbone");
			prog.setParameter("output",dir+"/"+ outputName);
			prog.setParameter("inputFiles",input);

			String[] cmdArray = prog.buildCmdArray();
			for(int i = 0; i < cmdArray.length; i++){
				System.out.print(cmdArray[i]+" ");
			}
			System.out.println();


			prog.exec(null, stdout, stderr);
			prog.clearResources();
			prog = null;
			/*			ExtendedWriter B = null;

			if(dir.length() > 0)
				B = new ExtendedWriter(new FileWriter(dir+"/"+FileName + ".aln"));
			else
				B = new ExtendedWriter(new FileWriter(FileName + ".aln"));

			B.println(stdout);
			B.flush();
			B.close();
			 */
			System.out.println(stdout);
			System.out.println(stderr);
		}
		catch(java.io.IOException IE){
			IE.printStackTrace();
			return false;
		}
		catch(Exception E){
			E.printStackTrace();
			return false;
		}
		return true;
	}



	public static boolean stralAlign(String dir , String FileName){
		try{
			StringBuffer stdout = new StringBuffer();
			StringBuffer stderr = new StringBuffer();

			LocalProgram prog = new LocalProgram("stral");
			prog.loadParameters();
			prog.loadResources();

			prog.setParameter("n","FALSE");

			if(dir.length() > 0)
				prog.setParameter("File",dir+"/"+ FileName+".fasta");
			else
				prog.setParameter("File", FileName+".fasta");

			prog.exec(null, stdout, stderr);


			ExtendedWriter B = null;
			if(dir.length() > 0)
				B = new ExtendedWriter(new FileWriter(dir+"/"+FileName + ".stral"));
			else
				B = new ExtendedWriter(new FileWriter(FileName + ".stral"));

			B.println(stdout);
			B.flush();
			B.close();

		}
		catch(java.io.IOException IE){
			IE.printStackTrace();
			return false;
		}
		catch(Exception E){
			E.printStackTrace();
			return false;
		}
		return true;
	}

	public static boolean convertSingleSeqence(String dir,String FileName,String input,String outputFormat){
		if(outputFormat.indexOf("aln") != -1)
			return convert2AlnFile( dir, FileName, input);
		else
			return convert2FastaFile( dir, FileName, input);


	}

	public static boolean compareSequences(int [] Seq1, int[] Seq2 ){
		if(Seq1.length != Seq2.length) return false;
		for(int i  = 0; i < Seq1.length; i++){
			if(Seq1[i]!= Seq2[i])return false;
		}
		return true;
	}


	public static boolean convert2FastaFile(String Dir,String FileName,String input){
		int nrOfGenomes = 0;
		try{
			ExtendedReader br = new ExtendedReader(new FileReader(Dir+"/"+FileName+"."+input));
			ExtendedWriter fr = new ExtendedWriter(new FileWriter(Dir+"/"+FileName+".fasta"));
			while(br.more()){
				char A = (char)br.readChar();
				if(A == '>')
					nrOfGenomes++;
				String GenomeName = br.readLine();
				fr.println(">"+ GenomeName);
				while(br.more())
					fr.println(br.readLine());
			}				
			br.close();
			fr.flush();
			fr.close();
		}

		catch(Exception E){
			E.printStackTrace();
			return false;
		}
		return true;
	}

	
	
	public static int[] merge5end(int [] Fend, int [] main, int length, int percentErrors, int inside){
		if(Fend.length > length){
			boolean found = false;
			int errors = percentErrors*length/100;
			int pos = 0;
			inside = Math.min(inside, main.length-length);
			inside = Math.min(inside,  Fend.length-length);
			int error = 0;
			while(!found &&  pos < inside ){
				int foundErrors = 0;
				int count = 0;
				while (length-count > 0 && foundErrors <= errors){
					if(Fend[Fend.length-length+count] != main[pos+count]){
						foundErrors++;
					}
					count++;
				}
				if(foundErrors <= errors){
					int totalErrors = 0;
					for(int i = Fend.length - length - pos;i < Fend.length; i++ ){
						if(Fend[i] != main[i-( Fend.length - length - pos)]){
							totalErrors++;
						}
					}
					if(totalErrors <= (percentErrors*(Fend.length - length - pos)/100)){
						found = true;
						error = totalErrors;
					}
				}
				if(!found){
					pos++;
				}
			}
			if(found){
				if(length + pos < Fend.length){
					System.out.print("Extending the sequence with "+ (Fend.length - length - pos)+ " nt, ");
					System.out.println( (pos+ length) + " nt overlap nr of errors "+error+" ("+Fend.length+":"+main.length+")");
					int[] newSeq = new int[Fend.length+main.length-(length + pos)];
					for(int i = 0; i < Fend.length; i++){
						newSeq[i] = Fend[i];
					}
					for(int i = length + pos; i < main.length; i++ ){
						newSeq[Fend.length-(length + pos)+i] = main[i];
					}
					main = newSeq;
				}
			}
			
		}

		return main;
	}

	
	
	
	public static int[] merge3end(int [] Tend, int [] main, int length, int percentErrors, int inside){
		if(Tend.length > length){
			boolean found = false;
			int errors = percentErrors*length/100;
			inside = Math.min(inside, Tend.length -length );
			inside = Math.min(inside, main.length-length);
			int pos = main.length - inside-length;
			int error = 0;
			while(!found &&  pos+length < main.length ){
				int foundErrors = 0;
				int count = 0;
				while (length-count > 0 && foundErrors <= errors){
					if(Tend[count] != main[pos+count]){
						foundErrors++;
					}
					count++;
				}
				if(foundErrors <= errors){
					int totalErrors = 0;
					for(int i = 0; i < main.length-pos;i++){
						if(Tend[i] != main[pos+i]){
							totalErrors++;
						}
					}
					if(totalErrors <= (percentErrors*(main.length - pos)/100)){
						found = true;
						error = totalErrors;
					}
				}
				if(!found){
					pos++;
				}
			}
			if(found){
				if(main.length -  pos < Tend.length){
					System.out.print("Extending the sequence with "+ (Tend.length - (main.length -  pos))+ " nt, ");
					System.out.println( (main.length -  pos) + " nt overlap and "+ error +" errors ("+Tend.length+":"+main.length+")");
					int[] newSeq = new int[Tend.length+main.length-(main.length- pos)];
					for(int i = 0; i < pos; i++){
						newSeq[i] = main[i];
					}
					for(int i = 0; i < Tend.length; i++ ){
						newSeq[pos+i] = Tend[i];
					}
					main = newSeq;
				}
			}
		}

		return main;
	}
	
	

	public static boolean convert2AlnFile(String Dir,String FileName,String input){
		int nrOfGenomes = 0;
		try{
			ExtendedReader br = new ExtendedReader(new FileReader(Dir+"/"+FileName+"."+input));
			ExtendedWriter fr = new ExtendedWriter(new FileWriter(Dir+"/"+FileName+".aln"));
			fr.println("CLUSTAL W (1.83) multiple sequence alignment");
			fr.println();
			fr.println();

			char[] nucleotides = null;

			String GenomeName = null;
			while(br.more()){
				char A = (char)br.readChar();
				if(A == '>')
					nrOfGenomes++;
				GenomeName = br.readLine();
				while(br.more())
					nucleotides = Functions.addChar(nucleotides,(char)br.readChar());
			}				
			br.close();
			for (int i = 0; i <= nucleotides.length / 50 ; i++){
				int j = 0;
				if(i*50 != nucleotides.length){
					fr.print(Functions.fixedLength(GenomeName,30));
					fr.print(Functions.fixedLength(6));
					while(i*50 + j < nucleotides.length && j < 50){
						fr.print(nucleotides[i*50+j]);
						j++;
					}
					fr.println();
					fr.print(Functions.fixedLength(30));
					fr.print(Functions.fixedLength(6));
					fr.print(Functions.fixedLength(j,'*'));
					fr.println();
					fr.println();


				}
			}
			fr.flush();
			fr.close();
		}

		catch(Exception E){
			E.printStackTrace();
			return false;
		}
		return true;
	}

	public static void printPhyFile(String Dir,String FileName,String[] GenomeNames, int[][] genomeSequences){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(Dir+"/"+FileName+".phy"));

			EW.println("\t"+GenomeNames.length+"\t"+genomeSequences[0].length);
			for(int i = 0; i < GenomeNames.length; i++){
				EW.println(Functions.fixedLength(GenomeNames[i],10)+" "+DNAInt2String(genomeSequences[i]));
			}
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}
	}


	public static boolean printAlnFile(String Dir,String FileName,String[] GenomeNames, int[][] genomeSequences){
		try{
			ExtendedWriter fr = new ExtendedWriter(new FileWriter(Dir+"/"+FileName+".aln"));
			fr.println("CLUSTAL W (1.83) multiple sequence alignment");
			fr.println();
			fr.println();
			for (int i = 0; i <= genomeSequences[0].length / 50 ; i++){
				for(int k = 0; k < GenomeNames.length; k++){
					int j = 0;
					if(i*50 != genomeSequences[0].length){
						fr.print(Functions.fixedLength(GenomeNames[k],30));
						fr.print(Functions.fixedLength(6));
						while(i*50 + j < genomeSequences[0].length && j < 50){
							fr.print(RNAInt2char(genomeSequences[k][i*50+j]));
							j++;
						}
						fr.println();
					}
				}
				int j = 0;
				if(i*50 != genomeSequences[0].length){
					fr.print(Functions.fixedLength(30));
					fr.print(Functions.fixedLength(6));
					while(i*50 + j < genomeSequences[0].length && j < 50){
						boolean theSame = true;
						for(int k = 0; k < GenomeNames.length; k++){
							if(theSame && genomeSequences[k][i*50+j] != genomeSequences[0][i*50+j])
								theSame = false;
						}
						if(theSame)
							fr.print('*');
						else
							fr.print(' ');
						j++;

					}

					fr.println();
				}
				fr.println();
				fr.println();

			}
			fr.flush();
			fr.close();
		}
		catch(Exception E){
			E.printStackTrace();
			return false;
		}
		return true;
	}

	public static boolean printFastaFile(String Dir,String FileName,String[] GenomeNames, int[][] genomeSequences){
		try{
			ExtendedWriter fr = new ExtendedWriter(new FileWriter(Dir+"/"+FileName+".fa"));
			for(int k = 0; k < GenomeNames.length; k++){
				fr.println(">" + GenomeNames[k]);
				fr.println(RNAInt2char(genomeSequences[k]));
			}

			fr.flush();
			fr.close();
		}
		catch(Exception E){
			E.printStackTrace();
			return false;
		}
		return true;
	}



	public static boolean printStructureFile(String Dir,String FileName,char[] ConsensusSequence, char[] StructureSequence){
		try{
			ExtendedWriter fr = new ExtendedWriter(new FileWriter(Dir+"/"+FileName+".pattern"));
			fr.println(new String (ConsensusSequence));
			fr.println(new String (StructureSequence));
			fr.flush();
			fr.close();
		}
		catch(Exception E){
			E.printStackTrace();
			return false;
		}
		return true;
	}




	public static boolean convertFile(String dir , String FileName,String input, String outputFormat){
		if(!IOTools.fileExists(dir+"/"+ FileName+"."+outputFormat)){
			int nrOfGenomes = getNrOfGenomesFasta(dir,FileName,input);
			if(nrOfGenomes == 1)
				return convertSingleSeqence(dir,FileName,input,outputFormat);
			try{
				StringBuffer stdout = new StringBuffer();
				StringBuffer stderr = new StringBuffer();

				LocalProgram prog = new LocalProgram("clustalw");
				prog.loadParameters();
				prog.loadResources();

				if(dir.length() > 0)
					prog.setParameter("infile",dir+"/"+ FileName+"."+input);
				else
					prog.setParameter("infile", FileName+"."+input);

				prog.setParameter("outorder","input");
				prog.setParameter("convert");
				prog.setParameter("output",outputFormat);
				prog.exec(null, stdout, stderr);
				//				System.out.println(stdout);
				//				System.out.println(stderr);

			}
			catch(java.io.IOException IE){
				IE.printStackTrace();
				return false;
			}
			catch(Exception E){
				E.printStackTrace();
				return false;
			}
		}
		return true;

	}

	public static boolean convertFile(String dir , String FileName,String input){
		try{
			StringBuffer stdout = new StringBuffer();
			StringBuffer stderr = new StringBuffer();

			LocalProgram prog = new LocalProgram("clustalw");
			prog.loadParameters();
			prog.loadResources();

			if(dir.length() > 0)
				prog.setParameter("infile",dir+"/"+ FileName+"."+input);
			else
				prog.setParameter("infile", FileName+"."+input);

			prog.setParameter("outorder","input");
			prog.setParameter("convert");
			prog.exec(null, stdout, stderr);

			//System.out.println(stdout);
			//System.out.println(stderr);
		}
		catch(java.io.IOException IE){
			IE.printStackTrace();
			return false;
		}
		catch(Exception E){
			E.printStackTrace();
			return false;
		}
		return true;

	}

	public static boolean generatePhyFile(String dir , String FileName, String Format){
		try{
			StringBuffer stdout = new StringBuffer();
			StringBuffer stderr = new StringBuffer();

			LocalProgram prog = new LocalProgram("clustalw");
			prog.loadParameters();
			prog.loadResources();

			if(dir.length() > 0)
				prog.setParameter("infile",dir+"/"+ FileName+"."+Format);
			else
				prog.setParameter("infile", FileName+"."+Format);

			prog.setParameter("outorder","input");
			prog.setParameter("output","phylip");
			prog.exec(null, stdout, stderr);

			//System.out.println(stdout);
			//System.out.println(stderr);
		}
		catch(java.io.IOException IE){
			IE.printStackTrace();
			return false;
		}
		catch(Exception E){
			E.printStackTrace();
			return false;
		}
		return true;
	}

	/*

	public static boolean writeRootedPhymlFile(String dir,String fileName, String[] Sequences,String[] Genomes,tree rootedTree){
		try{
			System.out.println(dir+"/"+fileName+ "_rooted.phy is about to be printed");
			ExtendedWriter B = null;
			B = new ExtendedWriter(new FileWriter(dir+"/"+fileName+ "_rooted.phy"));
			B.println("\t"+Sequences.length+"\t"+Sequences[0].length());
			for(int i = 0; i < Sequences.length;i++){
				B.println(Functions.fixedLength(Genomes[i],10)+" "+RNA2DNA(Sequences[i]));
			}
			B.println("1");
			rootedTree.print(B);

			B.flush();
			B.close();


		}
		catch(java.io.IOException IE){
			IE.printStackTrace();
			return false;
		}
		catch(Exception E){
			E.printStackTrace();
			return false;
		}
		return true;

	}

	 */

	public static boolean generatePhymlScript(String dir , String FileName,String TrTv,String invar,String alpha,String Tree, String topology, String Branch){

		try{


			StringBuffer stdout = new StringBuffer();
			StringBuffer stderr = new StringBuffer();

			LocalProgram prog = new LocalProgram("phyml");
			prog.loadParameters();
			prog.loadResources();

			if(dir.length() > 0)
				prog.setParameter("infile",dir+"/"+ FileName+".phy");
			else
				prog.setParameter("infile", FileName+".phy");

			prog.setParameter("dataType","0");
			prog.setParameter("sequenceFormat","i");
			prog.setParameter("nrOfdatabases","1");
			prog.setParameter("nrOfbootstraps","0");
			prog.setParameter("SubstitutionModel","HKY");
			prog.setParameter("TrTVRation",TrTv);
			prog.setParameter("invar",invar);
			prog.setParameter("nb_categ","6");
			prog.setParameter("alpha",alpha);
			prog.setParameter("tree",Tree);
			prog.setParameter("topology",topology);
			prog.setParameter("branch",Branch);


			String[] cmdArray = prog.buildCmdArray();


			//			for(int i = 0; i < cmdArray.length; i++){
			//			System.out.println(cmdArray[i]);
			//			}
			//			System.out.println();

			for(int i = 0; i < cmdArray.length; i++){
				System.out.print(cmdArray[i]+" ");
			}
			System.out.println();



			//System.out.println("/mnt/E/phyml_v2.4.4/exe/phyml_linux "+dir+"/"+ FileName+".phy 0 i 1 0 HKY 2.4 0.11 1 e /mnt/e/microbology/COGsequences/targetRNAs/topology.tree.phy n y");

			prog.exec(null, stdout, stderr);

			//			System.out.println(stdout);
			//			System.out.println(stderr);
		}
		catch(java.io.IOException IE){
			IE.printStackTrace();
			return false;
		}
		catch(Exception E){
			E.printStackTrace();
			return false;
		}
		return true;
	}
	//	johan$ seq-gen -mhky -k 1 <ompA_TU0-4742.phy >ompA_seqGen.phy	
	public static boolean generateSeqGenScript(int numberOfGeneratedSequences,int length, String TrTv,String alpha,String freq, String TreeFile, String outFile,String fA, String fC, String fG, String fU,String invar,String categories,int root){
		try{
			StringBuffer stdout = new StringBuffer();
			StringBuffer stderr = new StringBuffer();

			LocalProgram prog = new LocalProgram("seqgen");
			prog.loadParameters();
			prog.loadResources();

			prog.setParameter("m","HKY");
			if(root > 0)
				prog.setParameter("k",Functions.Int2String(root));			
			prog.setParameter("n",Functions.Int2String(numberOfGeneratedSequences));
			prog.setParameter("l",Functions.Int2String(length));
			prog.setParameter("i",invar);
			prog.setParameter("f"); 
			prog.setParameter("fA",fA);
			prog.setParameter("fC",fC);
			prog.setParameter("fG",fG);
			prog.setParameter("fT",fU);

			prog.setParameter("g",categories);
			prog.setParameter("a",alpha);
			prog.setParameter("t",TrTv);
			prog.setParameter("q");
			prog.setParameter("o","r");
			prog.setParameter("treeFile",TreeFile);
			String[] cmdArray = prog.buildCmdArray();

			/*			

			for(int i = 0; i < cmdArray.length; i++){
				System.out.println(cmdArray[i]);
			}
			System.out.println();
			 */			
			for(int i = 0; i < cmdArray.length; i++){
				System.out.print(cmdArray[i]+" ");
			}
			System.out.println();



			//System.out.println("/mnt/E/phyml_v2.4.4/exe/phyml_linux "+dir+"/"+ FileName+".phy 0 i 1 0 HKY 2.4 0.11 1 e /mnt/e/microbology/COGsequences/targetRNAs/topology.tree.phy n y");

			prog.exec(null, stdout, stderr);

			ExtendedWriter B = new ExtendedWriter(new FileWriter(outFile));

			B.println(stdout);
			B.flush();
			B.close();

			//			System.out.println(stdout);
			System.out.println(stderr);
		}
		catch(java.io.IOException IE){
			IE.printStackTrace();
			return false;
		}
		catch(Exception E){
			E.printStackTrace();
			return false;
		}
		return true;
	}




	public static String getStructure(int[] RNAstructure){
		char[] RNAStruct = new char[RNAstructure.length];
		for(int j = 0; j < RNAstructure.length; j++){
			if(RNAstructure[j] > 0){
				if(RNAstructure[j]   == 1){
					RNAStruct[j] = ')';
				}
				else if(RNAstructure[j] == 2){
					RNAStruct[j] = '(';					
				}
				else 
					System.out.println(RNAStruct[j]);
			}
			else{
				if(RNAstructure[j] == -1)
					RNAStruct[j] = 'S';
				else if (RNAstructure[j] == -2)
					RNAStruct[j] = 'L';
				else 
					System.out.println(RNAStruct[j]);

			}
		}
		return new String(RNAStruct);
	}

	public static int[] fixRNA(int[] RNAseq){
		int length = RNAseq.length;
		int[] newRNAseq = new int[length];

		for(int i = 0; i < length; i++){
			if(RNAseq[i] > 5){
				newRNAseq[i] = 5;
			}
			else{
				newRNAseq[i] = RNAseq[i];
			}
		}
		return newRNAseq;


	}

	public static int getNrOfGenomesAln(String Dir,String FileName){
		String[] GenomeNames =null;
		try{
			ExtendedReader br = new ExtendedReader(new FileReader(Dir+"/"+FileName+".aln"));
			br.readLine();
			br.readLine();
			br.readLine();

			while(br.more()){
				String temp = br.readWord();
				if(temp.indexOf('*') ==  -1){
					boolean found = false;
					if(GenomeNames == null){
						GenomeNames = Functions.addString(GenomeNames, IOTools.fixFileName(temp));
					}
					else{
						for(int i = 0; i < GenomeNames.length; i++ ){
							if(GenomeNames[i].indexOf(temp) != -1){
								found = true;
							}
						}
						if(!found){
							GenomeNames = Functions.addString(GenomeNames, IOTools.fixFileName(temp));
						}
					}
				}
				br.skipLine();
			}
			br.close();
		}
		catch(Exception E){
			E.printStackTrace();
			return 0;
		}
		return GenomeNames.length;
	}

	public static int getNrOfGenomesPhy(String Dir,String FileName, String extension){
		int nrOfGenomes = 0;
		try{
			ExtendedReader br = new ExtendedReader(new FileReader(Dir+"/"+FileName+"."+extension));
			br.more();
			nrOfGenomes = br.readInt();
			br.close();
		}
		catch(Exception E){
			E.printStackTrace();
			return 0;
		}
		return nrOfGenomes;
	}


	/*	public static String[] readFasta(ExtendedReader ER){
		int nrOfGenomes = 0;
		try{
				while(br.more()){
					char A = (char)br.readChar();
					if(A == '>')
						nrOfGenomes++;
					br.skipLine();
				}				
				br.close();
		}

		catch(Exception E){
			E.printStackTrace();
			return 0;
		}
		return nrOfGenomes;
	}
	 */

	public static int getNrOfGenomesFasta(String Dir,String FileName, String extension){
		int nrOfGenomes = 0;
		try{
			ExtendedReader br = new ExtendedReader(new FileReader(Dir+"/"+FileName+"."+extension));
			while(br.more()){
				char A = (char)br.readChar();
				if(A == '>')
					nrOfGenomes++;
				br.skipLine();
			}				
			br.close();
		}

		catch(Exception E){
			E.printStackTrace();
			return 0;
		}
		return nrOfGenomes;
	}
	/*
	public static Gene loadStockholmRNA(String dir, String fileName){
		Gene RNA = new Gene();
		RNA.setName(fileName.substring(fileName.indexOf('.')));
		try{
			ExtendedReader br = new ExtendedReader(new FileReader(dir+"/"+fileName));
			br.skipLine();
			ArrayList<String> GenomeNames = new ArrayList<String>();
			ArrayList<String> Alignments = new ArrayList<String>();
			while(br.more()){
				if((char)br.lookAhead() != '/'){
					if((char)br.lookAhead() == '#'){
						br.readWord();
						String nxt = br.readWord();
						if(nxt.indexOf("STOCKHOLM")> -1) br.skipLine();
						else if(nxt.indexOf("SS_cons")> -1){
							String structure = br.readWord();
							IntramolecularStructures IS = IntramolecularStructures.addDotBracketAnnotation(structure.toCharArray());
							RNA.setConsensusStructure(IS);
						}
						else if(nxt.indexOf("RF")> -1){
							br.skipLine();
						}
					}
					else{
						String NameAndLocation = br.readWord();
						String Name = NameAndLocation.substring(0,NameAndLocation.indexOf('/'));
						GenomeNames.add(Name);
						Alignments.add(br.readWord());
					}
				}
				else
					br.skipLine();
			}
			br.close();
			for(int i = 0; i < GenomeNames.size(); i++){
				RNA.addGenomeAlignment(GenomeNames.get(i), Alignments.get(i)); 
			}

		}

		catch(Exception E){
			E.printStackTrace();
			return null;
		}
		return RNA;







	}


	public static sRNAGene loadStockholmsRNA(String dir, String fileName){
		sRNAGene RNA = new sRNAGene();
		RNA.setName(fileName.substring(0,fileName.indexOf('.')));
		try{
			ExtendedReader br = new ExtendedReader(new FileReader(dir+"/"+fileName));
			br.skipLine();
			ArrayList<String> GenomeNames = new ArrayList<String>();
			ArrayList<String> Alignments = new ArrayList<String>();
			while(br.more()){
				if((char)br.lookAhead() != '/'){
					if((char)br.lookAhead() == '#'){
						br.readWord();
						String nxt = br.readWord();
						if(nxt.indexOf("STOCKHOLM")> -1) br.skipLine();
						else if(nxt.indexOf("SS_cons")> -1){
							String structure = br.readWord();
							IntramolecularStructures IS = IntramolecularStructures.addDotBracketAnnotation(structure.toCharArray());
							RNA.setConsensusStructure(IS);
						}
						else if(nxt.indexOf("RF")> -1){
							br.skipLine();
						}
					}
					else{
						String NameAndLocation = br.readWord();
						String Name = NameAndLocation.substring(0,NameAndLocation.indexOf('/'));
						GenomeNames.add(Name);
						Alignments.add(br.readWord());
					}
				}
				else
					br.skipLine();
			}
			br.close();
			for(int i = 0; i < GenomeNames.size(); i++){
				RNA.addGenomeAlignment(GenomeNames.get(i), Alignments.get(i)); 
			}

		}

		catch(Exception E){
			E.printStackTrace();
			return null;
		}
		return RNA;
	}


	 */


	public static char[] getFastaSequence(String Dir,String FileName, String extension){
		char[] Sequence = null;
		try{
			ExtendedReader br = new ExtendedReader(new FileReader(Dir+"/"+FileName+"."+extension));
			while(br.more()){
				if((char)br.lookAhead() == '>'){
					br.skipLine();
				}
				else{
					Sequence = append(br.readLine().toCharArray(),Sequence);
				}
			}				
			br.close();
		}
		catch(Exception E){
			E.printStackTrace();
			return null;
		}
		return Sequence;
	}



	public static String[] getGenomes(String Dir,String FileName, String extension){
		try{
			ExtendedReader br = new ExtendedReader(new FileReader(Dir+"/"+FileName+"."+extension));
			int nrOfGenomes = br.readInt();
			br.skipLine();
			br.more();
			String[] Genomes = new String[nrOfGenomes];
			String[] Sequences = new String[nrOfGenomes];
			for(int i = 0;i < nrOfGenomes;i++){
				Genomes[i] = br.readWord();
				br.more();
				Sequences[i] = br.readLine();
			}
			br.close();
			return Genomes;
		}
		catch(Exception E){
			E.printStackTrace();
			return null;
		}
	}


	public static String[] getSequences(String Dir,String FileName, String extension){
		try{
			ExtendedReader br = new ExtendedReader(new FileReader(Dir+"/"+FileName+"."+extension));
			int nrOfGenomes = br.readInt();
			br.skipLine();
			br.more();
			String[] Genomes = new String[nrOfGenomes];
			String[] Sequences = new String[nrOfGenomes];
			for(int i = 0;i < nrOfGenomes;i++){
				Genomes[i] = br.readWord();
				br.more();
				Sequences[i] = br.readLine();
			}
			br.close();
			return Sequences;
		}
		catch(Exception E){
			E.printStackTrace();
			return null;
		}
	}


	public static int getLength(String Dir,String FileName, String extension){
		int length = 0;
		try{
			ExtendedReader br = new ExtendedReader(new FileReader(Dir+"/"+FileName+"."+extension));
			br.readInt();
			br.more();
			length = br.readInt();
			br.close();
		}
		catch(Exception E){
			E.printStackTrace();
			return 0;
		}
		return length;
	}



	public static int getNrOfGenomes(String Dir,String FileName, String extension){
		try{
			ExtendedReader br = new ExtendedReader(new FileReader(Dir+"/"+FileName+"."+extension));
			int nrOfGenomes = br.readInt();
			br.close();
			return nrOfGenomes;
		}
		catch(Exception E){
			E.printStackTrace();
			return 0;
		}

	}


	public static String[] getGenomeNamesFasta(String Dir,String FileName, String extension){
		String [] GenomeNames = null;
		try{
			ExtendedReader br = new ExtendedReader(new FileReader(Dir+"/"+FileName+"."+extension));
			while(br.more()){
				if((char)br.lookAhead() == '>')
					GenomeNames = Functions.addString(GenomeNames,br.readLine());
				else
					br.skipLine();
			}				
			br.close();
		}
		catch(Exception E){
			E.printStackTrace();
			return null;
		}
		return GenomeNames;

	}

	public static String[] getGenomeSequencesFasta(String Dir,String FileName, String extension){
		String [] GenomeSequences = null;
		String  GenomeSequence = "";
		try{
			ExtendedReader br = new ExtendedReader(new FileReader(Dir+"/"+FileName+"."+extension));
			while(br.more()){
				if((char)br.lookAhead() == '>'){
					br.skipLine();
					if(GenomeSequence.length() > 1){
						GenomeSequences = Functions.addString(GenomeSequences,GenomeSequence);
						GenomeSequence = "";
					}
				}
				else
					GenomeSequence = GenomeSequence + br.readLine();
			}				
			br.close();
			GenomeSequences = Functions.addString(GenomeSequences,GenomeSequence);
		}
		catch(Exception E){
			E.printStackTrace();
			return null;
		}
		return GenomeSequences;

	}

	public static boolean generatePatternFile(String dir,String FileName){
		try{

			StringBuffer stdout = new StringBuffer();
			StringBuffer stderr = new StringBuffer();

			LocalProgram prog = new LocalProgram("ALIfold");
			prog.loadParameters();
			prog.loadResources();

			if(dir.length() > 0)
				prog.setParameter("noLP", dir+"/"+ FileName+".aln");
			else
				prog.setParameter("noLP", FileName+".aln");

			prog.exec(null, stdout, stderr);


			ExtendedWriter B = null;
			if(dir.length() > 0)
				B = new ExtendedWriter(new FileWriter(dir+"/"+FileName + ".pattern"));
			else
				B = new ExtendedWriter(new FileWriter(FileName + ".pattern"));

			B.println(stdout);
			B.flush();
			B.close();
			//			System.out.println(stdout);
			//System.out.println(stderr);
		}
		catch(java.io.IOException IE){
			IE.printStackTrace();
			return false;
		}
		catch(Exception E){
			E.printStackTrace();
			return false;
		}
		return true;


	}





	public static boolean generateCONTRAfoldStructure(String dir,String FileName,String Extension,String structureFolder,double posteriorCutoff){
		if(!IOTools.isDir(structureFolder))
			IOTools.mkDir(structureFolder);
		int nrOfGenomes = getNrOfGenomesFasta(dir,FileName,Extension);
		if(nrOfGenomes == 1){
			if(!IOTools.fileExists(structureFolder,FileName+".fa"))
				try{
					StringBuffer stdout = new StringBuffer();
					StringBuffer stderr = new StringBuffer();

					LocalProgram prog = new LocalProgram("contrafold");
					prog.loadParameters();
					prog.loadResources();

					if(posteriorCutoff !=  -1){
						prog.setParameter("posteriors", Double.toString(posteriorCutoff));
						prog.setParameter("outputFile",structureFolder+"/"+FileName+".posteriors");
					}
					else{
						prog.setParameter("parens",structureFolder+"/"+FileName+".bpseq");
					}


					prog.setParameter("predict");
					//prog.setParameter("outputdir", structureFolder+"/");
					//prog.setParameter("parens", dir+"/"+FileName+".contrafold");
					prog.setParameter("input", dir+"/"+FileName+"."+ Extension);


					String[] cmdArray = prog.buildCmdArray();

					//					for(int i = 0; i < cmdArray.length; i++){
					//						System.out.print(cmdArray[i]+" ");
					//					}
					//					System.out.println();



					prog.exec(null, stdout, stderr);


					//			ExtendedWriter B = new ExtendedWriter(new FileWriter(structureFolder+"/"+FileName + ".structure"));
					//			B.println(stdout);
					//			B.flush();
					//			B.close();

					/*					System.out.println("testing");
					System.out.println(stdout);
					System.out.println(stderr);
					 */				}
			catch(java.io.IOException IE){
				IE.printStackTrace();
				return false;
			}
			catch(Exception E){
				E.printStackTrace();
				return false;
			}
			return true;
		}
		return true;

	}

	public static boolean generateCONTRAfoldStructure(String dir,String FileName,String Extension,double posteriorCutoff){
		int nrOfGenomes = getNrOfGenomesFasta(dir,FileName,Extension);
		if(nrOfGenomes == 1){
			if(!IOTools.fileExists(dir,FileName+".contrafold"))
				try{
					StringBuffer stdout = new StringBuffer();
					StringBuffer stderr = new StringBuffer();

					LocalProgram prog = new LocalProgram("contrafold");
					prog.loadParameters();
					prog.loadResources();

					if(posteriorCutoff !=  -1){
						prog.setParameter("posteriors", Double.toString(posteriorCutoff));
						//prog.setParameter("noncomplementary");

					}
					prog.setParameter("predict");
					//prog.setParameter("outputdir", structureFolder+"/");
					//prog.setParameter("parens", dir+"/"+FileName+".contrafold");
					prog.setParameter("input", dir+"/"+FileName+"."+ Extension);


					String[] cmdArray = prog.buildCmdArray();
					/*			
					for(int i = 0; i < cmdArray.length; i++){
						EW.print(cmdArray[i]);
					}
					EW.println();
					 */
					for(int i = 0; i < cmdArray.length; i++){
						System.out.println(cmdArray[i]);
					}
					System.out.println();


					prog.exec(null, stdout, stderr);


					ExtendedWriter B = new ExtendedWriter(new FileWriter(dir+"/"+FileName + ".contrafold"));
					B.println(stdout);
					B.flush();
					B.close();

					//					System.out.println("testing");
					//					System.out.println(stdout);
					//					System.out.println(stderr);
				}
			catch(java.io.IOException IE){
				IE.printStackTrace();
				return false;
			}
			catch(Exception E){
				E.printStackTrace();
				return false;
			}
			return true;
		}
		return true;

	}

	public static boolean generateRfoldProbabilities(String dir,String FileName,String Extension,String structureFolder,double posteriorCutoff,int distance){

		if(!IOTools.isDir(structureFolder))
			IOTools.mkDir(structureFolder);
		int nrOfGenomes = getNrOfGenomesFasta(dir,FileName,Extension);
		if(nrOfGenomes == 1){
			if(!IOTools.fileExists(structureFolder,FileName+"_"+distance+".rfold")){
				System.out.print(FileName+" is being folded ...........");
				try{
					StringBuffer stdout = new StringBuffer();
					StringBuffer stderr = new StringBuffer();

					LocalProgram prog = new LocalProgram("rfold");
					prog.loadParameters();
					prog.loadResources();

					prog.setParameter("command","COMPUTE_PROB");
					if(distance >  0)
						prog.setParameter("max_pair_dist",Integer.toString(distance));
					if(posteriorCutoff >  -1)
						prog.setParameter("print_prob_cutoff", Double.toString(posteriorCutoff));
					if(distance > 0)
						prog.setParameter("prob_file", structureFolder+"/"+FileName+"_"+distance+".rfold");
					else
						prog.setParameter("prob_file", structureFolder+"/"+FileName+".rfold");

					prog.setParameter("input", dir+"/"+FileName+".fa");

					String[] cmdArray = prog.buildCmdArray();


					//					for(int i = 0; i < cmdArray.length; i++){
					//					System.out.println(cmdArray[i]);
					//					}

					prog.exec(null, stdout, stderr);

					//				
					//				ExtendedWriter B = new ExtendedWriter(new FileWriter(structureFolder+"/"+FileName + ".structure"));
					//				B.println(stdout);
					//				B.flush();
					//				B.close();

					//					System.out.println(stdout);
					//					System.out.println(stderr);
					System.out.println("finished");

				}
				catch(java.io.IOException IE){
					IE.printStackTrace();
					return false;
				}
				catch(Exception E){
					E.printStackTrace();
					return false;
				}
				return true;
			}
		}
		return true;

	}

	public static boolean generateRfoldProbabilities(String dir,String FileName,String Extension,String structureFolder,double posteriorCutoff){
		return generateRfoldProbabilities(dir,FileName,Extension,structureFolder,posteriorCutoff,-1);
	}



	public static int[][] findLongestStretch(int length1, int length2, int nrOfGenes, int nrOfGenomes, double evolutionRate){
		int[][] stretches = new int[nrOfGenes][];
		for(int i = 0 ; i < nrOfGenes; i++){
			stretches[i] = findLongestStretch(length1,length2,nrOfGenomes,evolutionRate);
		}
		return stretches;
	}


	public static int[] findLongestStretch(int length1, int length2,int nrOfGenomes, double evolutionRate){
		int[] sRNA = generateIntSequence(length1);
		int[] mRNA = generateIntSequence(length2);
		int canoncialStretch = FindLongestConservedStretch(sRNA,mRNA,getCanonicalBasepairsMM(), nrOfGenomes, evolutionRate);
		int basepairStretch = FindLongestConservedStretch(sRNA,mRNA,getBasepairsMM() , nrOfGenomes, evolutionRate);
		int[] stretches = new int[2];
		stretches[0] = canoncialStretch;
		stretches[1] = basepairStretch;

		return stretches;
	}





	public static int FindLongestConservedStretch(int[] sRNA, int[] mRNA, boolean[][] MM, int nrOfGenomes, double evolutionRate){
		int sRNAstart,mRNAstart,mRNAstop,sRNAstop;
		sRNAstart = mRNAstart = 0;
		sRNAstop = sRNA.length;
		mRNAstop = mRNA.length;
		int longestStretch = 1;
		int[][] evolvedsRNASequences = new int[nrOfGenomes][sRNA.length];
		evolvedsRNASequences[0]= sRNA;
		int[][] evolvedmRNASequences = new int[nrOfGenomes][sRNA.length];
		evolvedmRNASequences[0] = mRNA;
		evolutionRate = evolutionRate*4/3;
		for(int k = 1; k < nrOfGenomes;k++){
			evolvedsRNASequences[k]= new int[sRNA.length];
			evolvedmRNASequences[k] = new int[mRNA.length];
			for(int i = mRNAstop-1; i > mRNAstart; i--){
				if(Math.random() > 1 - evolutionRate){
					evolvedmRNASequences[k][i]= RNAfunctions.generateNucleotide();
				}
				else
					evolvedmRNASequences[k][i]= mRNA[i];
			}
			for(int j = sRNAstop-1; j > sRNAstart; j--){
				if(Math.random() > 1 - evolutionRate)
					evolvedsRNASequences[k][j]= RNAfunctions.generateNucleotide();

				else
					evolvedsRNASequences[k][j]= sRNA[j];
			}
		}

		int[][] length = new int[mRNAstop- mRNAstart+1][sRNAstop-sRNAstart+1];
		int[][] Clength = new int[mRNAstop- mRNAstart+1][sRNAstop-sRNAstart+1];
		for(int i = mRNAstop-1; i > mRNAstart; i--){
			for(int j = sRNAstop-1; j > sRNAstart; j--){
				if(MM[sRNA[j]][mRNA[i]]){
					length[i][j] = length[i+1][j+1]+1;
					boolean cbp = true;
					for(int k = 1; k < nrOfGenomes;k++){
						if(!MM[evolvedmRNASequences[k][i]][evolvedsRNASequences[k][j]])
							cbp = false;
					}
					if(cbp)
						Clength[i][j] = Clength[i+1][j+1]+1;
					if(Clength[i][j] >= longestStretch  ){
						longestStretch = Clength[i][j];
						//System.out.println(longestStretch);

						//System.out.println(Clength[i][j] +" - "+length[i][j]);

						//}
					}
				}
			}
		}
		return longestStretch;
	}


	public int FindLongestConservedStretch(int[][] sRNA, int[][] mRNA, boolean[][] MM){
		int sRNAstart,mRNAstart,mRNAstop,sRNAstop;
		sRNAstart = mRNAstart = 0;
		sRNAstop = sRNA.length;
		mRNAstop = mRNA.length;
		int longestStretch = 1;
		int[][] length = new int[mRNAstop- mRNAstart+1][sRNAstop-sRNAstart+1];
		for(int i = mRNAstop-1; i > mRNAstart; i--){
			for(int j = sRNAstop-1; j > sRNAstart; j--){
				boolean conserved = true;
				int k = 0;
				while(conserved &&  k < Math.min(sRNA.length, mRNA.length)){
					if(!MM[mRNA[k][i]][sRNA[k][j]])
						conserved = false;
					k++;
				}
				if(conserved)
					length[i][j] = length[i+1][j+1]+1;
				if(length[i][j] > longestStretch )
					longestStretch = length[i][j];
			}
		}
		return longestStretch;
	}


	public static boolean generateRfoldStructures(String dir,String FileName,String Extension,String structureFolder,int distance){
		if(!IOTools.isDir(structureFolder))
			IOTools.mkDir(structureFolder);
		int nrOfGenomes = getNrOfGenomesFasta(dir,FileName,Extension);
		if(nrOfGenomes == 1){
			if(!IOTools.fileExists(structureFolder,FileName+"_"+distance+".rfold_structure"))
				try{
					StringBuffer stdout = new StringBuffer();
					StringBuffer stderr = new StringBuffer();

					LocalProgram prog = new LocalProgram("rfold");
					prog.loadParameters();
					prog.loadResources();

					prog.setParameter("command","COMPUTE_MEA_FOLD");
					if(distance >  0)
						prog.setParameter("max_pair_dist",Integer.toString(distance));


					prog.setParameter("outfile", structureFolder+"/"+FileName+"_"+distance+".rfold_structure");

					prog.setParameter("input", dir+"/"+FileName+".fa");

					String[] cmdArray = prog.buildCmdArray();


					for(int i = 0; i < cmdArray.length; i++){
						System.out.println(cmdArray[i]);
					}
					System.out.println();

					prog.exec(null, stdout, stderr);

					//				
					//				ExtendedWriter B = new ExtendedWriter(new FileWriter(structureFolder+"/"+FileName + ".structure"));
					//				B.println(stdout);
					//				B.flush();
					//				B.close();

					System.out.println(stdout);
					System.out.println(stderr);
				}
			catch(java.io.IOException IE){
				IE.printStackTrace();
				return false;
			}
			catch(Exception E){
				E.printStackTrace();
				return false;
			}
			return true;
		}
		return true;

	}

	public static double[] getFreq(int[] sequence){
		double[] freqs = new double[4];
		int count  = 0;
		for(int i = 0; i < sequence.length; i++){
			if(sequence[i] > 0 && sequence[i]< 5)
				freqs[sequence[i]-1]++;
			count++;
		}
		for(int i = 0; i < freqs.length; i++){
			freqs[i] = freqs[i]/sequence.length;
		}
		return freqs;
	}


	public static int[] generateSequence(int len, double[] freqs){
		int[] seq = new int[len];
		int i = 0;
		double p;
		freqs[1] += freqs[0];
		freqs[2] += freqs[1];
		freqs[3] += freqs[2];
		while (i < len)
		{
			p = Math.random();
			if (p < freqs[0])
				seq[i] = 1;
			else if (p < freqs[1])
				seq[i] = 2;
			else if (p < freqs[2])
				seq[i] = 3;
			else
				seq[i] = 4;
			i++;
		}
		return seq;
	}

	public static char findSynonym(char[] tempNT){
		if(tempNT == null){
			System.out.println("Errors in findSynonym: no sequence in");
			return 'N';
		}

		boolean[] newNT = new boolean[6];
		for(int i = 0; i < newNT.length; i++){
			newNT[i] = false;    
		}	    
		int[] intNT = RNAchar2int2(tempNT);
		for(int i = 1; i < tempNT.length; i++){
			newNT[intNT[i]] = true;
		}

		if(newNT[5])
			return 'N';// A+C+G+T
		if(newNT[1]){
			if(newNT[2]){
				if(newNT[3]){
					if(newNT[4])
						return 'N';// A+C+G+T
					else
						return 'V';//A+C+G
				}
				else if(newNT[4])
					return 'H';//A+C+T
				else
					return 'M' ;// A+C
			}
			else if(newNT[3]){// A+G
				if(newNT[4])//A+G+T
					return 'D';// A+G+T
				else
					return 'R';// A+G
			}
			else if(newNT[4])
				return 'W';//A+T
			else
				return 'A';
		}
		if(newNT[2]){
			if(newNT[3]){
				if(newNT[4])
					return 'B';// C+G+T
				else
					return 'S';//C+G
			}
			else if(newNT[4])
				return 'Y';//C+T
			else
				return 'C' ;// C
		}
		if(newNT[3]){
			if(newNT[4])
				return 'K';//G+T
			else
				return 'G';// G
		}
		else if(newNT[4])
			return 'U';//T
		else
			return 'N';
	}

	public static int findSynonym(int[] tempNT){
		return RNAchar2int(findCharSynonym(tempNT));
	}

	public static char findCharSynonym(int[] tempNT){
		if(tempNT == null){
			System.out.println("Errors in findSynonym: no sequence in");
			return 'N';
		}

		boolean[] newNT = new boolean[6];
		for(int i = 0; i < newNT.length; i++){
			newNT[i] = false;    
		}	    
		for(int i = 1; i < tempNT.length; i++){
			if(tempNT[i] < 6)
				newNT[tempNT[i]] = true;
			else
				newNT[5] = true;
		}

		if(newNT[5])
			return 'N';// A+C+G+T
		if(newNT[0])
			return 'N';
		if(newNT[1]){
			if(newNT[2]){
				if(newNT[3]){
					if(newNT[4])
						return 'N';// A+C+G+T
					else
						return 'V';//A+C+G
				}
				else if(newNT[4])
					return 'H';//A+C+T
				else
					return 'M' ;// A+C
			}
			else if(newNT[3]){// A+G
				if(newNT[4])//A+G+T
					return 'D';// A+G+T
				else
					return 'R';// A+G
			}
			else if(newNT[4])
				return 'W';//A+T
			else
				return 'A';
		}
		if(newNT[2]){
			if(newNT[3]){
				if(newNT[4])
					return 'B';// C+G+T
				else
					return 'S';//C+G
			}
			else if(newNT[4])
				return 'Y';//C+T
			else
				return 'C' ;// C
		}
		if(newNT[3]){
			if(newNT[4])
				return 'K';//G+T
			else
				return 'G';// G
		}
		else if(newNT[4])
			return 'U';//T
		else
			return 'N';
	}


	public static char DNA2RNA(char Nucleotide){
		if(Nucleotide == 'T'||Nucleotide == 't')
			return 'U';
		return  Nucleotide;

	}

	public static String RNA2DNA(String Nucleotides){
		return new String(RNA2DNA(Nucleotides.toCharArray()));
	}

	public static char[] RNA2DNA(char[] Nucleotides){
		for(int i = 0; i < Nucleotides.length; i++ )
			Nucleotides[i] = RNA2DNA(Nucleotides[i]);
		return  Nucleotides;

	}
	public static char RNA2DNA(char Nucleotide){
		if(Nucleotide == 'U'||Nucleotide == 'u')
			return 'T';
		return  Nucleotide;

	}


	public static int[][] getAntisenseMatrix(int M, int GU, int MM){
		int[][] AM = new  int[][]{
				{MM,MM,MM,MM,MM,MM},
				{MM,MM,MM,MM,M ,MM},
				{MM,MM,MM,M ,MM,MM},
				{MM,MM,M ,MM,GU,MM},
				{MM,M ,MM,GU,MM,MM},
				{MM,MM,MM,MM,MM,MM}
		};
		return AM;
	}

	public static int[][] getSenseMatrix(int M, int GU, int MM){
		int[][] AM = new  int[][]{
				{MM,MM,MM,MM,MM,MM},
				{MM,M ,MM,GU,MM,MM},
				{MM,MM,M ,MM,GU,MM},
				{MM,GU,MM,M ,MM,MM},
				{MM,MM,GU,MM,M ,MM},
				{MM,MM,MM,MM,MM,MM}
		};
		return AM;
	}


	public static void copyFastaFile(String fileName, String oldDir, String newDir,boolean addUnknown){
		try{
			String[] GenomeNames = RNAfunctions.getGenomeNamesFasta(oldDir,fileName,"fa");
			String[] Sequences = RNAfunctions.getGenomeSequencesFasta(oldDir,fileName,"fa");

			ExtendedWriter outFile = new ExtendedWriter(new FileWriter(newDir+"/"+fileName+".fa"));
			for(int i = 0; i < GenomeNames.length; i++){
				outFile.println(GenomeNames[i]);
				outFile.println(Sequences[i]);
			}
			if(GenomeNames.length == 2 && addUnknown){
				outFile.println(">unknown");
				outFile.println(Sequences[1]);
			}
			outFile.flush();
			outFile.close();
		}
		catch(Exception E){E.printStackTrace();}

	}


	public static String convertToScore(int[] B){
		int Length = B.length;
		String A = "";
		for(int i = 0; i < Length;i++){
			A = A+" "+String.valueOf(B[i]);
		}
		return A;
	}

	public static String convertToAllignment(int[] B){
		if(B == null){
			return "";
		}

		int Length = B.length;
		char[] IntArray = new char[Length];
		for(int i = 0; i < Length;i++){
			if(B[i] == 0){IntArray[i]  = ' ';}
			else if(B[i] == 1){IntArray[i]  = '|';}
			else if(B[i] == 2){IntArray[i]  = '!';}
			else{IntArray[i]  = (char)B[i];}

		}
		return new String(IntArray);
	}

	public static String convertToMatchPattern(int[] B){
		if(B == null){
			return "";
		}

		int Length = B.length;
		char[] IntArray = new char[Length];
		for(int i = 0; i < Length;i++){
			if(B[i] == 0){IntArray[i]  = ' ';}
			else if(B[i] == 1){IntArray[i]  = 'M';}
			else{IntArray[i]  = (char)B[i];}
		}
		return new String(IntArray);
	}

	public static int[][] splitRNA(int[] RNA,int start,int stop){
		int length = RNA.length;
		int[] leftArray = new int[start];
		int[] rightArray = new int[length-stop-1];

		for(int i = 0; i < start;i++){
			leftArray[i] = RNA[i];
		}
		int i = 0;
		while(stop+1+i < length){
			rightArray[i] = RNA[stop+i+1];
			i++;
		}

		int[][] intArrays = new int[3][0];
		intArrays[0] = leftArray;
		intArrays[1] = RNA;
		intArrays[2] = rightArray;

		return intArrays;
	}

	public static double[][] splitWeightVectors(double[] RNA,int start,int stop){
		int length = RNA.length;
		double[] leftArray = new double[start];
		double[] rightArray = new double[length-stop-1];

		for(int i = 0; i < start;i++){
			leftArray[i] = RNA[i];
		}
		int i = 0;
		while(stop+1+i < length){
			rightArray[i] = RNA[stop+i+1];
			i++;
		}

		double[][] intArrays = new double[3][0];
		intArrays[0] = leftArray;
		intArrays[1] = RNA;
		intArrays[2] = rightArray;

		return intArrays;
	}

	public static int[] getSubsequence(int[] RNA,int start,int stop){
		int[] subArray = new int[stop-start];

		for(int i = 0; i < stop-start;i++){
			subArray[i] = RNA[start+i];
		}

		return subArray;
	}

	public static String Randomize(String sequence){

		char[] Seq = sequence.toCharArray();

		int length = Seq.length;
		char altSeq;

		for(int i = length-1; i > -1 ; i--){
			int pick = (int)(Math.random() * (double)i);
			altSeq = Seq[pick];
			Seq[pick] = Seq[i];
			Seq[i] = altSeq;
		}

		String NewSeq = new String(Seq);

		return NewSeq;
	}

	public String generateStringSequence(int length){
		return new String(generateCharSequence(length));
	}

	private char[] generateCharSequence(int length){
		return RNAInt2char(generateIntSequence(length));
	}

	private static int[] generateIntSequence(int length){
		int[] Sequence = new int[length];
		for(int i = 0; i < length; i++){
			Sequence[i] = generateNucleotide(); 
		}
		return Sequence;

	}

	private static int generateNucleotide(){
		int nucleotide = (int)(Math.random() * (double)4)+ 1;
		return nucleotide;	
	}



	public static void ReadNrOfGenes(String FileName){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(FileName));
			int nrOfGenes = 0;			
			while(ER.more()){
				String newLine = ER.readLine();
				if(newLine.indexOf(">") != -1){
					nrOfGenes++;
					System.out.println(newLine);
				}
			}

			System.out.println(nrOfGenes);

			ER.close();


		}
		catch(Exception E){E.printStackTrace();}
	}

	public static double[] getKmerFrequencies(String[][] seqs, int k) {
		//double[] freqs;
		int n = 0;
		for (int i=1; i<=k; i++)
			n += (int) Math.pow(4,i);
		double[] totalFreqs = new double[n];
		byte[] seq = new byte[0];
		System.out.println("seqs.length = "+seqs.length);
		for (int i=0; i<seqs.length; i++) {
			if(seqs[i][1] != null)
				seq = append(addElement((byte) 0,encodeToByte(seqs[i][1])),seq);
		}
		totalFreqs = getKmerFrequencies(seq,k);
		return totalFreqs;
	}

	public static double[] getKmerFrequencies(byte[] seq, int k) {
		if (k <= 0)
			return new double[0];
		int noMotifs = 0;
		for (int j=1; j<=k; j++)
			noMotifs += (int) Math.pow(4,j);
		double[] n = new double[noMotifs];
		int[] count = new int[k];
		int index;
		int len;
		byte[] submotif;
		for (int i=0; i<(seq.length-k+1); i++) {
			submotif = subarray(i,i+k,seq);
			index = getKmerIndex(submotif);
			len = k;
			while (len > 0) {
				if (min(submotif) >= 1 && max(submotif) <= 4) {
					n[index]++;
					count[len-1]++;
					index = (index/4)-1;
				}
				else
					index = getKmerIndex(subarray(0,len-1,submotif));
				submotif[len-1] = 1;
				len--;
			}
		}
		for (int i=k-1; i>0; i--) {
			submotif = subarray(seq.length-i,seq.length,seq);
			index = getKmerIndex(submotif);
			len = i;
			while (len > 0) {
				if (min(submotif) >= 1 && max(submotif) <= 4) {
					n[index]++;
					count[len-1]++;
					index = (index/4)-1;
				}
				else
					index = getKmerIndex(subarray(0,len-1,submotif));
				submotif[len-1] = 1;
				len--;
			}
		}

		double[] freqs = new double[0];
		int lower;
		int upper;
		for (int i=1; i<=k; i++) {
			lower = getKmerIndex(newByteArray(i,(byte) 1));
			upper = getKmerIndex(newByteArray(i,(byte) 4));
			freqs = append(vectorDivide(subarray(lower,upper+1,n),count[i-1]),freqs);
		}

		return freqs;
	}	





	public static byte[] generateSequence(int len, int order, double[] freqs) {
		byte[] seq = new byte[len];
		double[][] kmers;
		double r;
		double[] cumFreqs;
		int next = 0;
		int j;
		int offset = -1;
		Random rand = new Random();
		while (next < order && next < len) {
			kmers = getConditionalKmers(freqs,next);
			r = rand.nextDouble();
			if (next == 0)
				cumFreqs = cumsum(kmers[0]);
			else
				cumFreqs = cumsum(kmers[getKmerIndex(subarray(0,next,seq))-offset]);
			j = 0;
			while (j < 4 && r >= cumFreqs[j])
				j++;
			seq[next] = (byte) Math.min((j+1),4);
			offset += (int) Math.pow(4,next);
			next++;
		}
		kmers = getConditionalKmers(freqs,order);
		for (int i=0; i<kmers.length; i++)
			kmers[i] = cumsum(kmers[i]);
		if (order == 0) {
			while (next < len) {
				j = 0;
				r = rand.nextDouble();
				while (j < 4 && r >= kmers[0][j])
					j++;
				seq[next] = (byte) Math.min((j+1),4);
				next++;
			}
		}
		else {
			byte[] prev = subarray(0,order,seq);
			int index;
			while (next < len) {
				index = getKmerIndex(prev)-offset;
				j = 0;
				r = rand.nextDouble();
				while (j < 4 && r >= kmers[index][j])
					j++;
				seq[next] = (byte) Math.min((j+1),4);
				prev = lshift(prev,(byte) Math.min((j+1),4));
				next++;
			}
		}
		return seq;
	}

	public static int getKmerIndex(byte[] seq) {
		int k = seq.length;
		int index = -1;
		for (int j=1; j<=k; j++)
			index += ((int) Math.pow(4,(j-1))) + ((int) Math.pow(4,(k-j)))*(seq[j-1]-1);
		return index;
	}  

	public static double[][] getConditionalKmers(double[] freqs, int order) {
		if (order == 0)
			return new double[][] {{freqs[0],freqs[1],freqs[2],freqs[3]}};
		double[][] kmers = new double[(int) Math.pow(4,order)][4];
		int lower = getKmerIndex(newByteArray(order+1,(byte) 1));
		int upper = getKmerIndex(newByteArray(order+1,(byte) 4));
		int k = lower;
		for (int i=0; k<upper; i++) {
			for (int j=0; j<4; j++) {
				if (freqs[k/4-1] != 0)
					kmers[i][j] = freqs[k]/freqs[k/4-1];
				else
					kmers[i][j] = 0;
				k++;
			}
		}
		return kmers;
	}  

	public static double[] cumsum(double[] arr) {
		double[] nArr = new double[arr.length];
		nArr[0] = arr[0];
		for (int i=1; i<arr.length; i++)
			nArr[i] = nArr[i-1]+arr[i];
		return nArr;
	}
	public static char[] decodeToChar(byte[] src, String cd) {
		char[] decoded = new char[src.length];
		for (int i=0; i<src.length; i++)
			decoded[i] = cd.charAt(src[i]);
		return decoded;
	}
	public static byte[] encodeToByte(String src) {
		String codeString = "NACGTN";
		src = src.toUpperCase();
		src = src.replace('T','U');
		codeString = codeString.replace('T','U');
		int sz = src.length();
		byte[] encoded = new byte[sz];
		for (int i=0; i<sz; i++)
			encoded[i] = (byte) Math.max(codeString.indexOf(src.charAt(i)),0);
		return encoded;
	}
	public static int indexOf(String obj, String[] arr) {
		if (arr == null || arr.length == 0)
			return -1;
		int sz = arr.length;
		for (int i=0; i<sz; i++)
			if (arr[i] != null && arr[i].compareTo(obj) == 0)
				return i;
		return -1;
	}

	public static byte[] lshift(byte[] arr, byte n) {
		for (int i=1; i<arr.length; i++)
			arr[i-1] = arr[i];
		arr[arr.length-1] = n;
		return arr;
	}
	public static String[][] parseBigFasta(File f) throws Exception	{
		int INITIAL_CAPACITY = 100;
		int INCREMENT = 50;
		int nextEntry = -1;
		String[][] entries = new String[INITIAL_CAPACITY][2];
		BufferedReader br = new BufferedReader(new FileReader(f));
		int chunkSize = 10*1024*1024;
		int rd;
		char[] cBuff;
		char[] entry = new char[chunkSize];
		String name;
		int j = 0;
		while (true) {
			cBuff = new char[chunkSize];
			rd = br.read(cBuff,0,chunkSize);
			if (rd < 0)
			{
				if (nextEntry >= 0)
					entries[nextEntry][1] = new String(subarray(0,j,entry));
				nextEntry++;
				break;
			}
			cBuff = subarray(0,rd,cBuff);
			for (int i=0; i<rd; i++)
			{
				if (j >= entry.length)
					entry = append(new char[chunkSize],entry);
				if (cBuff[i] != 10 && cBuff[i] != 13 && cBuff[i] != 62 && cBuff[i] != 32)
					entry[j] = cBuff[i];
				else if (cBuff[i] == 10 || cBuff[i] == 13 || cBuff[i] == 32)
					j--;
				else
				{
					if (nextEntry >= 0)
						entries[nextEntry][1] = new String(subarray(0,j,entry));
					entry = new char[chunkSize];
					j = -1;
					i++;
					name = new String();
					while (i < rd && cBuff[i] != 10 && cBuff[i] != 13)
					{
						name += cBuff[i];
						i++;
					}
					nextEntry++;
					if (nextEntry == entries.length)
						entries = append(new String[INCREMENT][2],entries);
					entries[nextEntry][0] = name;
				}
				j++;
			}
		}
		br.close();
		return subarray(0,nextEntry,entries);
	}


	public static byte[] addElement(byte n, byte[] arr) {
		if (arr == null || arr.length == 0)
			return new byte[] {n};
		int sz = arr.length;
		byte[] newArr = new byte[sz+1];
		for (int i=0; i<sz; i++)
			newArr[i] = arr[i];
		newArr[sz] = n;
		return newArr;
	}
	public static byte[] append(byte[] nE, byte[] src) {
		if (src == null)
			src = new byte[0];
		int sLt = src.length;
		int nLt = nE.length;
		byte[] newArr = new byte[sLt+nLt];
		for (int i=0; i<sLt; i++)
			newArr[i] = src[i];
		for (int i=0; i<nE.length; i++)
			newArr[sLt+i] = nE[i];
		return newArr;
	}
	public static char[] append(char[] nE, char[] src) {
		if (src == null)
			src = new char[0];
		int sLt = src.length;
		int nLt = nE.length;
		char[] newArr = new char[sLt+nLt];
		for (int i=0; i<sLt; i++)
			newArr[i] = src[i];
		for (int i=0; i<nE.length; i++)
			newArr[sLt+i] = nE[i];
		return newArr;
	}
	public static double[] append(double[] nE, double[] src) {
		if (src == null)
			src = new double[0];
		int sLt = src.length;
		int nLt = nE.length;
		double[] newArr = new double[sLt+nLt];
		for (int i=0; i<sLt; i++)
			newArr[i] = src[i];
		for (int i=0; i<nE.length; i++)
			newArr[sLt+i] = nE[i];
		return newArr;
	}
	public static String[][] append(String[][] nE, String[][] src) {
		if (src == null || src.length == 0)
			return nE;
		int sLt = src.length;
		int nLt = nE.length;
		String[][] newArr = new String[sLt+nLt][];
		for (int i=0; i<sLt; i++)
			newArr[i] = src[i];
		for (int i=0; i<nLt; i++)
			newArr[sLt+i] = nE[i];
		return newArr;
	}

	public static byte max(byte[] arr) {return (byte) max(vectorMultiply(arr,1.0));}	
	public static double max(double[] arr) {return maxElement(arr);}
	public static double maxElement(double[] arr) {
		int sz = arr.length;
		if (sz == 0)
			return 0;
		double max = arr[0];
		for (int i=1; i<sz; i++)
			if (arr[i] > max)
				max = arr[i];
		return max;
	}

	public static byte min(byte[] arr) {return (byte) min(vectorMultiply(arr,1.0));}
	public static double min(double[] arr) {return -1.0*max(vectorMultiply(arr,-1.0));}
	public static byte[] newByteArray(int len, byte fill) {
		byte[] c = new byte[len];
		for (int i=0; i<len; i++)
			c[i] = fill;
		return c;
	}
	public static double[] parseFreqs(File infile) throws Exception {
		String[] lines = readFile(infile).split(System.getProperty("line.separator"));
		double[] freqs = new double[lines.length-1];
		int next = 0;
		while (next < lines.length && lines[next+1].compareTo("# //") != 0) {
			freqs[next] = Double.parseDouble(lines[next+1]);
			next++;
		}
		freqs = subarray(0,next,freqs);
		return freqs;
	}
	public static final String readFile(File f) throws Exception {
		char[] buff = new char[(int) f.length()];
		FileReader fr = new FileReader(f);
		fr.read(buff,0,buff.length);
		return new String(buff);
	}
	public static byte[] subarray(int start, int stop, byte[] src) {
		if (start >= src.length || stop > src.length)
			return new byte[0];
		byte[] newArr = new byte[stop-start];
		for (int i=start; i<stop; i++)
			newArr[i-start] = src[i];
		return newArr;
	}	
	public static char[] subarray(int start, int stop, char[] src) {
		if (start >= src.length || stop > src.length)
			return new char[0];
		char[] newArr = new char[stop-start];
		for (int i=start; i<stop; i++)
			newArr[i-start] = src[i];
		return newArr;
	}	
	public static double[] subarray(int start, int stop, double[] src) {
		if (start >= src.length || stop > src.length)
			return new double[0];
		double[] newArr = new double[stop-start];
		for (int i=start; i<stop; i++)
			newArr[i-start] = src[i];
		return newArr;
	}
	public static String[][] subarray(int start, int stop, String[][] src)	{
		if (start >= src.length || stop > src.length)
			return new String[0][0];
		String[][] newArr = new String[stop-start][];
		for (int i=start; i<stop; i++)
			newArr[i-start] = src[i];
		return newArr;
	}
	public static double[] vectorDivide(double[] v1, double f) {return vectorMultiply(v1,1./f);}
	public static double[] vectorMultiply(byte[] v1, double f) {
		double[] arr = new double[v1.length];
		for (int i=0; i<v1.length; i++)
			arr[i] = (double) v1[i];
		return vectorMultiply(arr,f);
	}
	public static double[] vectorMultiply(double[] v1, double f)
	{
		int sz = v1.length;
		double[] nV = new double[sz];
		for (int i=0; i<sz; i++)
			nV[i] = v1[i] * f;
		return nV;
	}



	public static String cfasta2fasta(String cfasta){
		char[] colorSeq = cfasta.toCharArray();
		int startCodon = RNAchar2int(colorSeq[0]);
		int[] ColorSequence = new int[colorSeq.length-1];
		for(int i = 1; i < colorSeq.length;i++)
			ColorSequence[i-1] = colorSeq[i]-48;
		return DNAInt2String(CS2RNA(ColorSequence,startCodon)).substring(1);

	}


	public static void main(String[] args) {
		/*	int[] RNAsequence = RNAString2Int("TGCGGAACCGGCATGTCGTCGTTGCATCTCTTACTCCTTGGGCCCCGTC");
		System.out.println("TGCGGAACCGGCATGTCGTCGTTGCATCTCTTACTCCTTGGGCCCCGTC");
		int[] CS = RNA2CS(RNAsequence);
		System.out.print(RNAInt2char(RNAsequence[0]));
		for(int i = 0; i <CS.length; i++){
			System.out.print(CS[i]);
		}
		System.out.println();

		int[] sequence = CS2RNA(CS,RNAsequence[0]);

		System.out.println(DNAInt2String(sequence));




		String CSstring = "T011220303020201233232";
		System.out.println(">Consensus");
		System.out.println(CSstring);
		char[] colorSeq = CSstring.toCharArray();
		int startCodon = RNAchar2int(colorSeq[0]);
		int[] ColorSequence = new int[colorSeq.length-1];
		for(int i = 1; i < colorSeq.length;i++)
			ColorSequence[i-1] = colorSeq[i]-48;
		sequence = CS2RNA(ColorSequence,startCodon);
		System.out.println(DNAInt2String(sequence));

		 */


		double A = Double.parseDouble(args[3]);
		double C = Double.parseDouble(args[4]);
		double G = Double.parseDouble(args[5]);
		double U = Double.parseDouble(args[6]);
		int length = Integer.parseInt(args[7]);
		String[] restrictionSites = null;
		restrictionSites = Functions.addString(restrictionSites, args[1]);
		restrictionSites = Functions.addString(restrictionSites, args[2]);

		for(int i =8 ; i < args.length; i++){
			restrictionSites = Functions.addString(restrictionSites, args[i]);
		}

		System.out.println("miRNAsequence    \t : "+args[0]);
		System.out.println(getReverseComplement(args[0]));
		System.out.println("RE left sequence \t : "+args[1]);
		System.out.println("RE right sequence\t : "+args[2]);

		for(int i = 0; i < 10;i++){
			int[] primer = printAnimalmiRNAprimer(args[0],args[1],args[2],A,C,G,U,length);
			if(countRestrictionSites(primer,restrictionSites) == 2){
				System.out.println(">primer "+i);
				System.out.println(RNAInt2String(primer));
				System.out.println(">primer "+i+"_reverse complement");
				System.out.println(RNAInt2String(getReverseComplement(primer)));
			}
		}

		for(int i = 0; i < 10;i++){
			int[] primer = printPerfectmiRNAprimer(args[0],args[1],args[2],A,C,G,U,length);
			if(countRestrictionSites(primer,restrictionSites) == 2){
				System.out.println(">primer "+i+"_perfect");
				System.out.println(RNAInt2String(primer));
				System.out.println(">primer "+i+"_perfect_reverse complement");
				System.out.println(RNAInt2String(getReverseComplement(primer)));
			}
		}

	}


}






