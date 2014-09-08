package general;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Hashtable;

import alignment.StructuralVariation;

import general.ExtendedReader;


public abstract class Functions{

	public static void main(String []args){
		//to test the different subfunctions 
		//System.out.println(getFileWithoutSuffix("test.fastq","fastq"));
	}
	
	
	public static int countOccurrences(String haystack, String needleString){
		char[] needles = needleString.toCharArray(); 
	    int count = 0;
	    for (int i=0; i < haystack.length(); i++)
	    {
	        for(int j = 0; j < needles.length; j++){
	        	if (haystack.charAt(i) == needles[j]){
	             count++;
	             j = needles.length;
	        	}
	        }
	    }
	    return count;
	}
	

	
	public static int[] sum(int[] A, int [] B){
		if(A.length != B.length ){
			System.out.println("Error in Function.Sum");
			System.out.println("Two arrays of different lengths can not be summed");
			return null;
		}
		for(int i = 0; i < A.length; i++){
			A[i] = A[i]+B[i];
		}
		return A;
	}

	public static int[] subtract(int[] A, int [] B){
		if(A.length != B.length ){
			System.out.println("Error in Function.subtract");
			System.out.println("Two arrays of different lengths can not be subtracted");
			return null;
		}
		for(int i = 0; i < A.length; i++){
			A[i] = A[i]-B[i];
		}
		return A;
	}
	
	public static Hashtable<String, String> parseCommandLine(String[] args)
	{
		int len = args.length;
		Hashtable<String, String> commands = new Hashtable<String, String>(len);
		String key;
		String val;
		for (int i=0; i<len; i++)
		{
			if (args[i].charAt(0) == '-')
			{
				key = args[i];
				val = new String();
				if(i+1 < args.length && args[i+1].charAt(0) == '['){
					
					while (i+1 < args.length && (args[i+1].charAt(args[i+1].length()-1) != '[' || args[i+1].charAt(0) != '['))
					{
						val += args[i+1] + " ";
						i++;
					}
					val = val.substring(1,val.length()-2);
				}
				else{
					while (i+1 < args.length && args[i+1].charAt(0) != '-')
					{
						val += args[i+1] + " ";
						i++;
					}
				}
				val = val.trim();
				System.out.println(key+"\t"+val);
				commands.put(key,val);
			}
			else
				commands.put(args[i],new String());
		}
		return commands;
	}

	public static int[] copy(int[] original){
		int[] copy = new int[original.length];
		for(int i = 0; i < original.length;i++){
			copy[i] = original[i];
		}
		return copy;
	}

	public static double[][] copy(double[][] original){
		double[][] copy = new double[original.length][original[0].length];
		for(int i = 0; i < original.length;i++){
			for(int j = 0; j < original[i].length;j++){
				copy[i][j] = original[i][j];
			}
		}
		return copy;
	}

	public static String getDateTime() {
		DateFormat dateFormat = new SimpleDateFormat("yyyyMMdd_HH_mm_ss");
		Date date = new Date();
		return dateFormat.format(date);
	}


	public static int min(int a,int b){
		if(a < b) return a;
		return b;
	}

	public static int max(int a,int b){
		if(a > b) return a;
		return b;
	}

	public static int nearestInt(double score){
		int doubleScore =(int)(score * 10);
		int IntScore = ((int)score) * 10;
		if(doubleScore - IntScore > 0)
			return (int)(score+1);
		else
			return (int)(score);
	}

	public static int[] addInt(int[] intArray, int newInt){
		if(intArray == null){
			intArray = new int[1];

			intArray[0] = newInt;

			return intArray;
		}
		int Length = intArray.length;;
		int [] newIntArray = new int[Length +1];
		for(int i = 0; i < Length;i++){
			newIntArray[i] = intArray[i];
		}

		newIntArray[Length] = newInt;
		return newIntArray;
	}


	public static int[] addInt(int[] intArray, int[] newIntArray){
		if(intArray == null){
			return newIntArray;
		}
		int Length = intArray.length;;
		int [] totalIntArray = new int[Length +newIntArray.length];
		for(int i = 0; i < Length;i++){
			totalIntArray[i] = intArray[i];
		}
		for(int i = 0; i < newIntArray.length;i++){
			totalIntArray[Length+i] = newIntArray[i];
		}

		return totalIntArray;
	}


	public static boolean[] addBoolean(boolean[] booleanArray, boolean newBoolean){
		if(booleanArray == null){
			booleanArray = new boolean[1];

			booleanArray[0] = newBoolean;

			return booleanArray;
		}
		int Length = booleanArray.length;;
		boolean [] newBooleanArray = new boolean[Length +1];
		for(int i = 0; i < Length;i++){
			newBooleanArray[i] = booleanArray[i];
		}

		newBooleanArray[Length] = newBoolean;
		return newBooleanArray;
	}

	public static double[] addDouble(double[] doubleArray, double newDouble){
		if(doubleArray == null){
			doubleArray = new double[1];
			doubleArray[0] = newDouble;
			return doubleArray;
		}
		int Length = doubleArray.length;;
		double [] newDoubleArray = new double[Length +1];
		for(int i = 0; i < Length;i++){
			newDoubleArray[i] = doubleArray[i];
		}

		newDoubleArray[Length] = newDouble;
		return newDoubleArray;
	}

	public static double[] removeDouble(double[] doubleArray){
		if(doubleArray == null){
			return doubleArray;
		}
		if(doubleArray.length == 1)
			return null;

		int Length = doubleArray.length;;
		double [] newDoubleArray = new double[Length - 1];
		for(int i = 0; i < Length - 1;i++){
			newDoubleArray[i] = doubleArray[i];
		}

		return newDoubleArray;
	}

	public static int countOccurrences(String haystack, char[] needles)
	{
		int count = 0;
		for (int i=0; i < haystack.length(); i++){
			for(int j = 0; j < needles.length; j++){
				if (haystack.charAt(i) == needles[j]){
					//System.out.println("jsut checking");
					count++;
				}
			}
		}
		return count;
	}	


	public static Hashtable<String, String> getKeys(Hashtable<String,String> HT, String file){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(file));
			while(ER.more()){
				HT.put(ER.readLine(), " ");
			}
			ER.close();
		}
		catch(Exception E){

		}

		return HT;
	}

	public static ArrayList<String> getStrings(String file){
		ArrayList<String> HT = new ArrayList<String>();
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(file));
			while(ER.more()){
				HT.add(ER.readLine());
			}
			ER.close();
		}
		catch(Exception E){

		}

		return HT;
	}

	public static int[] removeInt(int[] intArray){
		if(intArray == null){
			return intArray;
		}
		if(intArray.length == 1)
			return null;

		int Length = intArray.length;
		int [] newIntArray = new int[Length - 1];
		for(int i = 0; i < Length - 1;i++){
			newIntArray[i] = intArray[i];
		}

		return newIntArray;
	}

	public static String[][] removeStrings(String[][] StringArrays, int pointer){
		if(StringArrays == null){
			return null;
		}
		if(StringArrays.length == 1)
			return null;

		int Length = StringArrays.length;
		String [][] newStringArrays = new String[Length - 1][];
		for(int i = 0; i < pointer;i++){
			newStringArrays[i] = StringArrays[i];
		}
		for(int i = pointer+1; i < Length;i++){
			newStringArrays[i-1] = StringArrays[i];
		}

		return newStringArrays;
	}



	public static boolean contains(int[] intArray,int query){
		if(intArray == null){
			return false;
		}
		for(int i = 0; i < intArray.length ;i++){
			if(intArray[i] == query)
				return true;
		}
		return false;
	}

	public static String getCommonPrefix(String seq1,String seq2){
		char[] CA1 = seq1.toCharArray();
		char[] CA2 = seq2.toCharArray();
		int pointer = 0;
		boolean common = true;
		while(pointer < CA1.length && common ){
			if(CA1[pointer]!= CA2[pointer]){
				common = false;
				pointer--;
			}
			pointer++;
		}
		return seq1.substring(0,pointer);
	}


	public static String[] addString(String[] StringArray, String newString){
		if(StringArray == null){
			StringArray = new String[1];
			StringArray[0] = newString;

			return StringArray;
		}
		int Length = StringArray.length;;
		String [] newStringArray = new String[Length +1];
		for(int i = 0; i < Length;i++){
			newStringArray[i] = StringArray[i];
		}

		newStringArray[Length] = newString;
		return newStringArray;
	}

	public static String[] addStrings(String[] StringArray, String[] newStrings){
		if(StringArray == null){
			return newStrings;
		}

		int Length = StringArray.length;
		int otherLength = newStrings.length;
		String [] newStringArray = new String[Length + otherLength];

		for(int i = 0; i < Length;i++){
			newStringArray[i] = StringArray[i];
		}
		for(int i = 0; i < otherLength;i++){
			newStringArray[Length+i] = newStrings[i];
		}
		return newStringArray;
	}



	public static char[] addChar(char[] charArray, char newchar){
		if(charArray == null){
			charArray = new char[1];
			charArray[0] = newchar;

			return charArray;
		}
		int Length = charArray.length;;
		char [] newcharArray = new char[Length +1];
		for(int i = 0; i < Length;i++){
			newcharArray[i] = charArray[i];
		}

		newcharArray[Length] = newchar;
		return newcharArray;
	}

	public static char[][] addCharArray(char[][] charArray, char[] newcharA){
		if(charArray == null){
			charArray = new char[1][0];
			charArray[0] = newcharA;

			return charArray;
		}
		int Length = charArray.length;;
		char [][] newcharArray = new char[Length +1][0];
		for(int i = 0; i < Length;i++){
			newcharArray[i] = charArray[i];
		}

		newcharArray[Length] = newcharA;
		return newcharArray;
	}


	public static int[][] addIntArray(int[][] intArray, int[] newintA){
		if(intArray == null){
			intArray = new int[1][0];
			intArray[0] = newintA;

			return intArray;
		}
		int Length = intArray.length;;
		int [][] newintArray = new int[Length +1][0];
		for(int i = 0; i < Length;i++){
			newintArray[i] = intArray[i];
		}

		newintArray[Length] = newintA;
		return newintArray;
	}
	public static String[][] addStringArray(String[][] intArray, String[] newintA){
		if(intArray == null){
			intArray = new String[1][0];
			intArray[0] = newintA;

			return intArray;
		}
		int Length = intArray.length;;
		String [][] newintArray = new String[Length +1][0];
		for(int i = 0; i < Length;i++){
			newintArray[i] = intArray[i];
		}

		newintArray[Length] = newintA;
		return newintArray;
	}

	public static String[][] addStringArray(String[][] intArray, String[][] newIntA){
		if(intArray == null){
			return newIntA;
		}
		if(newIntA == null)
			return intArray;
		int oldLength = intArray.length;
		int newLength = newIntA.length;
		String [][] newintArray = new String[oldLength+newLength][0];
		for(int i = 0; i < oldLength;i++){
			newintArray[i] = intArray[i];
		}
		for(int i = 0; i < newLength;i++){
			newintArray[oldLength+i] = newIntA[i];
		}
		return newintArray;
	}


	public static double[][] addDoubleArray(double[][] doubleArray, double[] newdoubleA){
		if(doubleArray == null){
			doubleArray = new double[1][0];
			doubleArray[0] = newdoubleA;

			return doubleArray;
		}
		int Length = doubleArray.length;;
		double [][] newdoubleArray = new double[Length +1][0];
		for(int i = 0; i < Length;i++){
			newdoubleArray[i] = doubleArray[i];
		}

		newdoubleArray[Length] = newdoubleA;
		return newdoubleArray;
	}



	public static double[] getReverse(double[] Array){
		double[] newArray = new double[Array.length];
		for(int i = 0; i < Array.length; i++){
			newArray[Array.length-i-1] = Array[i];
		}
		return newArray;
	}

	public static int[] getReverse(int[] Array){
		int[] newArray = new int[Array.length];
		for(int i = 0; i < Array.length; i++){
			newArray[Array.length-i-1] = Array[i];
		}
		return newArray;
	}

	public static char[] getReverse(char[] Array){
		char[] newArray = new char[Array.length];
		for(int i = 0; i < Array.length; i++){
			newArray[Array.length-i-1] = Array[i];
		}
		return newArray;
	}

	public static int[] getSubarray(int[] Array, int start, int stop){
		int[] newArray = new int[stop-start];
		for(int i = start; i < stop; i++){
			newArray[i-start] = Array[i];
		}
		return newArray;
	}

	public static String[] getSubarray(String[] Array, int start, int stop){
		String[] newArray = new String[stop-start];
		for(int i = start; i < stop; i++){
			newArray[i-start] = Array[i];
		}
		return newArray;
	}


	public static int[] getSubarray(int[] Array, int start){
		int stop = Array.length;
		int[] newArray = new int[stop-start];
		for(int i = start; i < stop; i++){
			newArray[i-start] = Array[i];
		}
		return newArray;
	}

	public static String[] getSubarray(String[] Array, int start){
		int stop = Array.length;
		if(stop <= start)
			return null;
		String[] newArray = new String[stop-start];
		for(int i = start; i < stop; i++){
			newArray[i-start] = Array[i];
		}
		return newArray;
	}


	public static char[] getSubarray(char[] Array, int start, int stop){
		char[] newArray = new char[stop-start];
		for(int i = start; i < stop; i++){
			newArray[i-start] = Array[i];
		}
		return newArray;
	}

	public static char[] getSubarray(char[] Array, int start){
		int stop = Array.length;
		char[] newArray = new char[stop-start];
		for(int i = start; i < stop; i++){
			newArray[i-start] = Array[i];
		}
		return newArray;
	}


	public static boolean[] getSubarray(boolean[] Array, int start, int stop){
		boolean[] newArray = new boolean[stop-start];
		for(int i = start; i < stop; i++){
			newArray[i-start] = Array[i];
		}
		return newArray;
	}

	public static boolean[] getSubarray(boolean[] Array, int start){
		int stop = Array.length;
		boolean[] newArray = new boolean[stop-start];
		for(int i = start; i < stop; i++){
			newArray[i-start] = Array[i];
		}
		return newArray;
	}




	public static double[] getSubarray(double[] Array, int start, int stop){
		double[] newArray = new double[stop-start];
		for(int i = start; i < stop; i++){
			newArray[i-start] = Array[i];
		}
		return newArray;
	}

	public static double[] getSubarray(double[] Array, int start){
		int stop = Array.length;
		double[] newArray = new double[stop-start];
		for(int i = start; i < stop; i++){
			newArray[i-start] = Array[i];
		}
		return newArray;
	}

	public static double[][] getSubmatrix(double[][] Matrix, int start1 ,int stop1,int start, int stop){
		double[][] newMatrix = new double[stop1-start1][stop-start];

		for(int i = start1; i < stop1; i++){
			for(int j = start; j < stop; j++){
				newMatrix[i-start1][j-start] = Matrix[i][j];
			}
		}
		return newMatrix;
	}

	public static int[][] getSubmatrix(int[][] Matrix, int start1 ,int stop1,int start, int stop){
		int[][] newMatrix = new int[stop1-start1][stop-start];

		for(int i = start1; i < stop1; i++){
			for(int j = start; j < stop; j++){
				newMatrix[i-start1][j-start] = Matrix[i][j];
			}
		}
		return newMatrix;
	}



	public static double[][] getSubmatrix(double[][] Matrix,int start1, int start){
		int stop1 = Matrix.length;
		int stop = Matrix[0].length;
		double[][] newMatrix = new double[stop1-start1][stop-start];

		for(int i = start1; i < stop1; i++){
			if(Matrix[i] != null){
				for(int j = start; j < stop; j++){
					if(Matrix[i].length > j){
						newMatrix[i-start1][j-start] = Matrix[i][j];
					}
				}
			}
		}
		return newMatrix;
	}



	public static char[][] getSubmatrix(char[][] Matrix, int start1 ,int stop1,int start, int stop){
		char[][] newMatrix = new char[stop1-start1][stop-start];

		for(int i = start1; i < stop1; i++){
			if(Matrix[i] != null){
				for(int j = start; j < stop; j++){
					if(Matrix[i].length > j){
						newMatrix[i-start1][j-start] = Matrix[i][j];
					}
				}
			}
			else
				newMatrix[i] = null;
		}
		return newMatrix;
	}




	public static char[][] getSubmatrix(char[][] Matrix,int start1, int start){
		int stop1 = Matrix.length;
		int stop = Matrix[0].length;
		char[][] newMatrix = new char[stop1-start1][stop-start];

		for(int i = start1; i < stop1; i++){
			if(Matrix[i] != null){
				for(int j = start; j < stop; j++){
					if(Matrix[i].length > j){
						newMatrix[i-start1][j-start] = Matrix[i][j];
					}
				}
			}
			else
				newMatrix[i] = null;
		}
		return newMatrix;
	}




	public static double[][][] addMatrix(double[][][] matrices, double[][] matrix){
		if(matrices == null){
			matrices = new double[1][0][0];
			matrices[0] = matrix;
			return matrices;
		}
		double [][][] newMatrices = new double[matrices.length +1][0][0];
		for(int i = 0; i < matrices.length; i++){
			newMatrices[i] = matrices[i];
		}

		newMatrices[matrices.length] = matrix;

		return newMatrices;
	}


	public static int findSequence(String[] Sequences, String Sequence){
		for(int i = 0; i < Sequences.length;i++){
			if(Sequences[i].compareTo(Sequence) == 0)
				return i;
		}
		return -1;
	}


	public static void printMatrix(double[][] matrix){
		for(int i = 0; i < matrix.length;i++){
			for(int j = 0; j < matrix[i].length; j++){
				System.out.print((int)matrix[i][j]+"\t");    
			}
			System.out.println();
		}
	}


	public static double[] initiateDoubleArray(int length, double weight){
		double[] array = new double[length];
		for(int i = 0; i < length; i++){
			array[i] = weight;
		}
		return array;
	}


	public static double[] changePre(double[] negativeScores){
		for(int i = 0; i < negativeScores.length;i++){
			negativeScores[i] = negativeScores[i]*-1;
		}
		return negativeScores;
	}

	public static int[] binScores(double[] unsortedScores, double min, double max, int nrOfBins){
		int[] sortedScores = new int[nrOfBins];
		double binLength = (max - min)/nrOfBins; 
		for(int i = 0; i < unsortedScores.length;i++){
			if(unsortedScores[i] < min)
				sortedScores[0]++;
			if(unsortedScores[i] >= max)
				sortedScores[nrOfBins-1]++;
			else{
				sortedScores[(int)((unsortedScores[i]- min)/binLength)]++;
			}
		}
		return sortedScores;
	}

	public static float findMaxFloat(ArrayList<Float> unsortedScores){
		float max = unsortedScores.get(0);
		for(int i = 1; i < unsortedScores.size();i++){
			if(unsortedScores.get(i) > max)
				max = unsortedScores.get(i);
		}
		return max;
	}
	public static int findMaxInt(ArrayList<Integer> unsortedScores){
		int max = unsortedScores.get(0);
		for(int i = 1; i < unsortedScores.size();i++){
			if(unsortedScores.get(i) > max)
				max = unsortedScores.get(i);
		}
		return max;
	}



	public static float findMinFLoat(ArrayList<Float> unsortedScores){
		float min = unsortedScores.get(0);
		for(int i = 1; i < unsortedScores.size();i++){
			if(unsortedScores.get(i) < min)
				min = unsortedScores.get(i);
		}
		return min;
	}

	public static int findMinInt(ArrayList<Integer> unsortedScores){
		int min = unsortedScores.get(0);
		for(int i = 1; i < unsortedScores.size();i++){
			if(unsortedScores.get(i) < min)
				min = unsortedScores.get(i);
		}
		return min;
	}


	public static int[] binScores(ArrayList<Float> unsortedScores,float min, float max,int nrOfBins){
		int[] sortedScores = new int[nrOfBins];

		double binLength = (max - min)/nrOfBins; 
		for(int i = 0; i < unsortedScores.size();i++){
			if(unsortedScores.get(i) < min)
				sortedScores[0]++;
			if(unsortedScores.get(i) >= max)
				sortedScores[nrOfBins-1]++;
			else{
				sortedScores[(int)((unsortedScores.get(i)- min)/binLength)]++;
			}
		}
		return sortedScores;
	}

	public static int[] binScores(ArrayList<Integer> unsortedScores,int min, int max,int nrOfBins){
		int[] sortedScores = new int[nrOfBins];
		double binLength = (max - min+1)/nrOfBins; 

		if(nrOfBins > max - min+1){
			sortedScores = new int[max - min+1];
			binLength = 1;
			nrOfBins = max -min+1;
		}

		for(int i = 0; i < unsortedScores.size();i++){
			if(unsortedScores.get(i) < min)
				sortedScores[0]++;
			if(unsortedScores.get(i) >= max)
				sortedScores[nrOfBins-1]++;
			else{
				sortedScores[(int)((unsortedScores.get(i)- min)/binLength)]++;
			}
		}
		return sortedScores;
	}



	public static int[] binScores(int[] unsortedScores, double min, double max, int nrOfBins){
		int[] sortedScores = new int[nrOfBins];
		double binLength = (max - min)/nrOfBins; 
		for(int i = 0; i < unsortedScores.length;i++){
			if(unsortedScores[i] < min)
				sortedScores[0]++;
			if(unsortedScores[i] >= max)
				sortedScores[nrOfBins-1]++;
			else{
				sortedScores[(int)((unsortedScores[i]- min)/binLength)]++;
			}
		}
		System.out.println();
		return sortedScores;
	}



	public static double[] sortScores(double[] unsortedScores){
		double[] sortedScores = null;
		for(int i = 0; i < unsortedScores.length;i++){
			sortedScores = Functions.sortScore(sortedScores, unsortedScores[i]);

		}
		return sortedScores;
	}




	public static double[] sortScoresAscending(double[] unsortedScores){
		double[] sortedScores = null;
		for(int i = 0; i < unsortedScores.length;i++){
			sortedScores = Functions.sortScoreAscending(sortedScores, unsortedScores[i]);

		}
		return sortedScores;
	}


	public static double[] sortScore(double[] oldArray,double newScore){
		if(oldArray == null){
			oldArray = new double[1];
			oldArray[0] = newScore;
			return oldArray;
		}
		double [] newArray = new double[oldArray.length + 1];
		int j = 0; 
		while(j < oldArray.length && newScore < oldArray[j]){
			newArray[j] = oldArray[j];
			j++;
		}
		newArray[j] = newScore;
		while(j < oldArray.length){
			newArray[j+1] = oldArray[j];
			j++;
		}
		return newArray;
	}

	public static int findLessThan(double[] Array,double Score){
		int j = 0;
		while(j< Array.length && Score <= Array[j]){
			j++;
		}
		return j;
	}


	public static int[] sortScore(int[] oldArray,int newScore){
		if(oldArray == null){
			oldArray = new int[1];
			oldArray[0] = newScore;
			return oldArray;
		}
		int [] newArray = new int[oldArray.length + 1];
		int j = 0; 
		while(j < oldArray.length && newScore < oldArray[j]){
			newArray[j] = oldArray[j];
			j++;
		}
		newArray[j] = newScore;
		while(j < oldArray.length){
			newArray[j+1] = oldArray[j];
			j++;
		}
		return newArray;
	}

	public static double[] sortScoreAscending(double[] oldArray,double newScore){
		if(oldArray == null){
			oldArray = new double[1];
			oldArray[0] = newScore;
			return oldArray;
		}
		double [] newArray = new double[oldArray.length + 1];
		int j = 0; 
		while(j < oldArray.length && newScore > oldArray[j]){
			newArray[j] = oldArray[j];
			j++;
		}
		newArray[j] = newScore;
		while(j < oldArray.length){
			newArray[j+1] = oldArray[j];
			j++;
		}
		return newArray;


	}





	public static int[] sortScoreAscending(int[] oldArray,int newScore){
		if(oldArray == null){
			oldArray = new int[1];
			oldArray[0] = newScore;
			return oldArray;
		}
		int [] newArray = new int[oldArray.length + 1];
		int j = 0; 
		while(j < oldArray.length && newScore > oldArray[j]){
			newArray[j] = oldArray[j];
			j++;
		}
		newArray[j] = newScore;
		while(j < oldArray.length){
			newArray[j+1] = oldArray[j];
			j++;
		}
		return newArray;
	}

	public static double getRank(double[] background,double score){

		if(background == null){
			return 1;
		}
		background = sortScores(background);

		int j = 0; 
		while(j < background.length && score < background[j]){
			j++;
		}
		return ((double)j/(double)background.length);


	}

	public static String fixedLength(float intIn, int length){
		return fixedLength((double)intIn,length);
	}

	public static String fixedLengthMax(double intIn, int length){
		String in = Double2String(intIn);
		if(in.length() > length)
			return in.substring(0,length);
		else 
			return in;
	}

	public static String fixedLength(double intIn, int length){
		String in = Double2String(intIn);
		if(in.length() > length)
			return in.substring(0,length);
		else if(in.length() == length)
			return in;
		else{
			char[] inChar = in.toCharArray();
			char[] outChar = new char[length];
			for(int i = 0; i < inChar.length; i++)
				outChar[i] = inChar[i];
			for(int i = inChar.length; i < outChar.length; i++)
				outChar[i] = ' ';
			return new String(outChar); 
		}
	}

	public static String fixedLength(int intIn, int length){
		String in = Int2String(intIn);
		if(in.length() > length)
			return in.substring(0,length);
		else if(in.length() == length)
			return in;
		else{
			char[] inChar = in.toCharArray();
			char[] outChar = new char[length];
			for(int i = 0; i < inChar.length; i++)
				outChar[i] = inChar[i];
			for(int i = inChar.length; i < outChar.length; i++)
				outChar[i] = ' ';
			return new String(outChar); 
		}
	}


	public static String fixedLength(String in, int length){
		if(in.length() > length)
			return in.substring(0,length);
		else if(in.length() == length)
			return in;
		else{
			char[] inChar = in.toCharArray();
			char[] outChar = new char[length];
			for(int i = 0; i < inChar.length; i++)
				outChar[i] = inChar[i];
			for(int i = inChar.length; i < outChar.length; i++)
				outChar[i] = ' ';
			return new String(outChar); 
		}
	}

	public static String fixedLengthCenter(String in, int length){
		if(in.length() > length)
			return in.substring(0,length);
		else if(in.length() == length)
			return in;
		else{
			char[] inChar = in.toCharArray();
			char[] outChar = new char[length];
			for(int i = 0; i < (length-inChar.length)/2; i++)
				outChar[i] = ' ';
			for(int i = 0; i < inChar.length; i++)
				outChar[(length-inChar.length)/2+i] = inChar[i];
			for(int i = inChar.length+(length-inChar.length)/2; i < outChar.length; i++)
				outChar[i] = ' ';
			return new String(outChar); 
		}
	}


	public static String fixedLength(int length){
		char[] outChar = new char[length];
		for(int i = 0; i < outChar.length; i++)
			outChar[i] = ' ';
		return new String(outChar); 
	}

	public static String fixedLength(int length,char specificChar){
		char[] outChar = new char[length];
		for(int i = 0; i < outChar.length; i++)
			outChar[i] = specificChar;
		return new String(outChar); 
	}

	public static String fixedLength(char specificChar,int length){
		char[] outChar = new char[length];
		outChar[0] = specificChar;
		for(int i = 1; i < outChar.length; i++)
			outChar[i] = specificChar;
		return new String(outChar); 
	}



	public static double[] getMeanArray(double[][] scores){
		double[] meanArray = new double[scores.length];
		for(int i = 0; i < scores.length; i++){
			meanArray[i] =getMean(scores[i]);
		}
		return meanArray;
	}

	public static double[] getMeanArray(int[][] scores){
		double[] meanArray = new double[scores.length];
		for(int i = 0; i < scores.length; i++){
			meanArray[i] =getMean(scores[i]);
		}
		return meanArray;
	}


	public static double  getMean(double[] scores){
		double mean = getSum(scores)/scores.length;

		return mean;
	}

	public static double  getMean(int[] scores){
		double mean = (double)getSum(scores)/scores.length;
		return mean;
	}



	public static double[][] transpose(double[][] Matrix){
		double[][] transpose = new double[Matrix[0].length][Matrix.length];
		for(int i = 0; i < Matrix.length; i++){
			for(int j = 0; j< Matrix[0].length; j++){
				transpose[j][i] = Matrix[i][j];
			}
		}
		return transpose;
	}


	public static int[][] transpose(int[][] Matrix){
		int[][] transpose = new int[Matrix[0].length][Matrix.length];
		for(int i = 0; i < Matrix.length; i++){
			for(int j = 0; j< Matrix[0].length; j++){
				transpose[j][i] = Matrix[i][j];
			}
		}
		return transpose;
	}


	public static int[][] transposeUnfilled(int[][] Matrix){

		ArrayList<ArrayList<Integer>> columns = new ArrayList<ArrayList<Integer>>();
		for(int i = 0; i < Matrix.length; i++){
			for(int j = 0; j< Matrix[i].length; j++){
				if(columns.size() <= j){
					ArrayList<Integer> row = new ArrayList<Integer>();
					columns.add(row);
				}
				if(Matrix[i].length > j)
					columns.get(j).add(Matrix[i][j]);
			}
		}

		int[][] transpose = new int[columns.size()][0];
		for(int j = 0; j< transpose.length; j++){
			transpose[j] = new int[columns.get(j).size()];
			for(int i = 0; i <columns.get(j).size(); i++)
				transpose[j][i] = columns.get(j).get(i);
		}
		return transpose;

	}




	public static double[] getSumArray(double[][] Matrix){
		double[] sumArray = new double[Matrix.length];
		for(int i = 0; i < Matrix.length; i++){
			sumArray[i] += getSum(Matrix[i]);
		}
		return sumArray;
	}


	public static double  getSum(double[] scores){
		double sum = 0;
		for(int i = 0; i < scores.length; i++){
			sum += scores[i];
		}
		return sum;
	}


	public static int  getSum(int[] scores){
		int sum = 0;
		for(int i = 0; i < scores.length; i++){
			sum += scores[i];
		}
		return sum;
	}


	public static double[] getSDArray(double[][] scores, double[] meanArray){
		double[] SDArray = new double[scores.length];
		for(int i = 0; i < scores.length; i++){
			SDArray[i] =getSD(scores[i], meanArray[i]);
		}
		return SDArray;
	}

	public static double[] getSDArray(int[][] scores, double[] meanArray){
		double[] SDArray = new double[scores.length];
		for(int i = 0; i < scores.length; i++){
			SDArray[i] =getSD(scores[i], meanArray[i]);
		}
		return SDArray;
	}



	public static double  getSD(double[] scores, double mean){
		double totalScore = 0;
		for(int i = 0; i < scores.length; i++){
			totalScore += java.lang.Math.pow((scores[i]-mean), 2);
		}
		totalScore = totalScore/scores.length;
		double SD = java.lang.Math.sqrt(totalScore);
		return SD;
	}

	public static double  getSD(int[] scores, double mean){
		double totalScore = 0;
		for(int i = 0; i < scores.length; i++){
			totalScore += java.lang.Math.pow(((double)scores[i]-mean), 2);
		}
		totalScore = totalScore/scores.length;
		double SD = java.lang.Math.sqrt(totalScore);
		return SD;
	}




	public static double getMeanFactor(double[] scores){
		double totalScore = 1;
		for(int i = 0; i < scores.length; i++){
			totalScore  = totalScore * scores[i];
		}
		//System.out.println(totalScore);
		double nrOfGenomes = (double)1/(double)scores.length;

		double mean = java.lang.Math.pow(totalScore,nrOfGenomes);

		return mean;
	}

	public static double getMinValue(double[] scores){
		double min = scores[0];
		for(int i = 1; i < scores.length; i++){
			if(scores[i] < min ) min = scores[i];
		}

		return min;
	}

	public static double getMinValue(double[][] scores){
		double min = getMinValue(scores[0]);
		for(int i = 1; i < scores.length; i++){
			double tempMin = getMinValue(scores[i]);
			if(tempMin < min ) 
				min = tempMin;
		}
		return min;
	}

	public static double getMaxValue(double[] scores){
		double max = scores[0];
		for(int i = 1; i < scores.length; i++){
			if(scores[i] > max ) max = scores[i];
		}
		return max;
	}

	public static double getMaxValue(double[][] scores){
		double max = getMaxValue(scores[0]);
		for(int i = 1; i < scores.length; i++){
			if(getMaxValue(scores[i]) > max ) max = getMaxValue(scores[i]);
		}
		return max;
	}

	public static int String2Int(String number){
		//		System.out.println(number);
		Integer value = new Integer(number);
		return value.intValue();

	}

	public static String Int2String(int number){
		return String.valueOf(number);
	}

	public static String Double2String(double number){
		return String.valueOf(number);
	}

	public static String[] toArray(String[][] matrix){
		int count = 0;
		for(int i = 0;i < matrix.length; i++){
			for(int j = 0; j< matrix[i].length; j++){
				count++;
			}
		}
		String[] array = new String[count];
		count = 0;
		for(int i = 0; i<matrix.length; i++){
			for(int j = 0; j< matrix[i].length; j++){
				array[count]= matrix[i][j];
				count++;
			}
		}
		return array;

	}


	public static double String2Double(String number){
		Double value = new Double(number);
		return value.doubleValue();

	}


	public static double[] getStdArray(double[][] scores){
		double[] StdArray = new double[scores.length];
		for(int i = 0; i < scores.length; i++){
			StdArray[i] =getStd(scores[i]);
		}
		return StdArray;
	}

	public  static double getStd(double[] scores){
		double sumOfScores = 0;
		double sumOfSquaredScores = 0;
		for(int i = 0; i < scores.length; i++){
			sumOfScores += scores[i];
			sumOfSquaredScores += scores[i]*scores[i];
		}
		double std = java.lang.Math.sqrt((sumOfSquaredScores-(sumOfScores*sumOfScores/scores.length))/(scores.length-1));
		return std;
	}	




	public static char String2Char(String SingleChar){
		char[] A = SingleChar.toCharArray();
		if(A.length > 1)
			System.out.println("String is longer than one char : "+ SingleChar);
		char B = A[0];
		return B;

	}


	public static String[] findUniqueStrings(String[] fastaFiles,String[] alnFiles){

		String[] unalignedRNAs = null;
		for(int i = 0; i < fastaFiles.length; i++){

			int j = 0;
			boolean found = false;
			if(alnFiles != null){
				while(j < alnFiles.length && !found){
					if(fastaFiles[i].compareTo(alnFiles[j]) == 0)
						found = true;
					j++;
				}
			}
			if(!found)
				unalignedRNAs = Functions.addString(unalignedRNAs,fastaFiles[i]);
		}
		return unalignedRNAs;
	}


	public static String[] splitUp(String concatenated){
		String[] newStrings = null;
		concatenated = IOTools.fixFileName(concatenated);
		newStrings = concatenated.split("_");
		return newStrings;

	} 


	public static boolean contains(ArrayList<String> StringArray,String newString){
		if(StringArray == null)
			return false;
		newString = IOTools.fixFileName(newString);
		for(int i = 0; i < StringArray.size(); i++){
			if(StringArray.get(i).compareTo(newString) == 0) return true;
		}
		return false;
	}


	public static boolean contains(String[] StringArray,String newString){
		if(StringArray == null)
			return false;
		newString = IOTools.fixFileName(newString);
		int length = newString.length();
		for(int i = 0; i < StringArray.length; i++){
			if(StringArray[i].length() >= length && StringArray[i].indexOf(newString) > -1)
				return true;
			else if(StringArray[i].length() < length && newString.indexOf(StringArray[i]) > -1)
				return true;
		}
		return false;
	}


	public static int getNrOfCommonStrings(String[] stringArray1, String[] stringArray2){
		int common = 0;
		for(int i = 0; i < stringArray1.length; i++){
			if(contains(stringArray2,stringArray1[i]))
				common++;
		}
		return common;

	}


	public static int getPosition(String[] StringArray,String newString){
		newString = IOTools.fixFileName(newString);
		int length = newString.length();
		for(int i = 0; i < StringArray.length; i++){
			if(StringArray[i].length() >= length && StringArray[i].indexOf(newString) > -1)
				return i;
			else if(StringArray[i].length() < length && newString.indexOf(StringArray[i]) > -1)
				return i;
		}
		return -1;
	}


	public static String[] parseString(String mRNAString){
		System.out.println(mRNAString);
		char[] mRNAArray = mRNAString.toCharArray();
		String [] mRNAs = null;
		char[] mRNA = null;
		for(int i = 0; i < mRNAArray.length; i++){
			if(mRNAArray[i] != ' ')
				mRNA = Functions.addChar(mRNA,mRNAArray[i]);
			else{
				mRNAs = Functions.addString(mRNAs,new String(mRNA));
				mRNA = null;
			}
		}
		mRNAs = Functions.addString(mRNAs,new String(mRNA));
		if(mRNAs != null){
			for(int i = 0; i < mRNAs.length; i++){
				System.out.println(mRNAs[i]);
			}
		}
		return mRNAs;

	}	


	public static ArrayList<String> mergeLists(ArrayList<String> listOne, ArrayList<String> listTwo){
		ArrayList<String> newList = new ArrayList<String>();
		newList.addAll(listOne);
		newList.addAll(listTwo);
		return newList;
	}

	public static ArrayList<String> getStringsWithPrefix(ArrayList<String> listOne, String prefix){
		ArrayList<String> newList = new ArrayList<String>();
		for(int i = 0; i< listOne.size(); i++){
			if(listOne.get(i).startsWith(prefix))
				newList.add(listOne.get(i));
		}
		return newList;
	}



	public static Object[] addObject(Object[] objects, Object newObject ){
		if(objects == null){
			objects = new Object[1];
			objects[0] = newObject;
			return objects;
		}
		Object[] newObjects = new Object[objects.length+1];
		for(int i = 0; i < objects.length;i++){
			newObjects[i] = objects[i];
		}
		newObjects[objects.length] = newObject;
		return newObjects;
	}


	public  static double getGumbellPValue(double score,double Xi , double Theta){
		return (1-(java.lang.Math.exp(-java.lang.Math.exp(-(score-Xi)/Theta))));
	}

	public static String getValue(Hashtable<String,String> T,String key, String defaultValue){
		if(T.containsKey(key)){
			String value =  T.get(key);
			if (value.length() > 0 ) 
				return value;
		}
		return defaultValue;
	}

	public static int getInt(Hashtable<String,String> T,String key, int defaultValue){
		if(T.containsKey(key)){
			String value =  T.get(key);
			if (value.length() > 0 ){
				return Integer.parseInt(value);
			}
		}
		return defaultValue;
	}

	public static double getDouble(Hashtable<String,String> T,String key, double defaultValue){
		if(T.containsKey(key)){
			String value =  T.get(key);
			if (value.length() > 0 ){
				return Double.parseDouble(value);
			}
		}
		return defaultValue;
	}

	public static String getFileWithoutSuffix(String filename, String suffix){
		if(filename.indexOf(suffix) > -1 && filename.lastIndexOf(suffix) + suffix.length() == filename.length()){
			String fileNameWithoutSuffix = filename.substring(0,filename.lastIndexOf(suffix));
			fileNameWithoutSuffix.trim();
			if(fileNameWithoutSuffix.lastIndexOf(".") == (fileNameWithoutSuffix.length()-1))
				fileNameWithoutSuffix = fileNameWithoutSuffix.substring(0,fileNameWithoutSuffix.length()-1);
			return fileNameWithoutSuffix;
		}
		System.out.println("Suffix not found");
		return filename;
	}

	public static String getValue(Hashtable<String,String> T,String key){
		if(T.containsKey(key)){
			String value =  T.get(key);
			if (value.length() > 0 ) 
				return value;
		}
		return null;
	}


}




