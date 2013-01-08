package miRNA;

import java.io.FileReader;

import org.apache.commons.math.stat.descriptive.moment.Mean;
import org.apache.commons.math.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math.stat.descriptive.summary.Sum;

import general.ExtendedReader;

public class RNAevalParser {
	int[] sequence;
	int miRNAstart;
	int miRNAstarStart;
	
	
	public static void main(String[] args) {
		RNAevalParser A = new RNAevalParser(args[0]);
		
	}
	
	
	RNAevalParser(String inFile){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(inFile));
			parseFile(ER);
		}
		catch(Exception E){E.printStackTrace();}
	}

	
	private void parseFile(ExtendedReader ER){
		double[] stats = new double[20000];
		int pointer = 0;
		while (ER.more()){
			while(ER.lookAhead() != '>')
				ER.skipLine();
			
			String[] info = ER.readLine().split(" ");
			String name = info[0];
			miRNAstart = Integer.parseInt(info[1].split("=")[1]);
			if(info.length > 2)
				miRNAstarStart = Integer.parseInt(info[2].split("=")[1]);
			else 
				miRNAstarStart = 0;
			double[] localDeltaG = new double[300];
			int leftStart = 1000;
			int leftStop =1;
			int rightStart = 1;
			int rightStop = 1;
		
			while(ER.more() && ER.lookAhead() != '>'){
			
			String Line = ER.readLine();
			if(Line.indexOf("Interior loop") == 0){
				info = Line.split("\\(");
				String from = info[1].split("\\)")[0];
				String to = info[2].split("\\)")[0];
				String[] info1 = from.split(",");
				String[] info2 = to.split(",");
				int fromLeft = Integer.parseInt(info1[0].trim());
				int fromRight= Integer.parseInt(info1[1].trim());
				int toLeft = Integer.parseInt(info2[0].trim());
				int toRight= Integer.parseInt(info2[1].trim());
				int DeltaG = Integer.parseInt(info[2].split(":")[1].trim());
				if(fromLeft < leftStart){
					leftStart = fromLeft;
					rightStart = fromRight;
				}
				if(toLeft > leftStop){
					leftStop = toLeft;
					rightStop = toRight;
					
				}
				for(int i = fromLeft; i <=toLeft;i++){
					localDeltaG[i] +=	(double)DeltaG/(double)(toLeft-fromLeft+1);
				}
				for(int i = fromRight; i >= toRight; i--){
					localDeltaG[i] += (double)DeltaG/(double)(fromRight-toRight+1);
				}
			}
		}
		
		for(int i = leftStart+1; i < leftStop-1; i++ ){
			stats[pointer] = localDeltaG[i];
			pointer++;
		}
		for(int i =rightStop+1; i < rightStart-1;i++){
			stats[pointer] = localDeltaG[i];
			pointer++;
		}
		
		
		double[] trimmedStats = new double[pointer];
		for(int i = 0; i < pointer; i++){
			trimmedStats[i] = stats[i];
//			System.out.println(stats[i]);
		}
		StandardDeviation STD = new StandardDeviation();
		Mean Mean = new Mean();
		Sum Sum = new Sum();
		double std = STD.evaluate(trimmedStats);
		double mean = Mean.evaluate(trimmedStats);
		double sum = Sum.evaluate(trimmedStats);
		
//		for(int i = leftStart; i < rightStart; i++){System.out.println(i +"\t"+ localDeltaG[i]);}
		
		System.out.print(name+"\t"+localDeltaG[miRNAstart-3]+"\t"+localDeltaG[miRNAstart-2]+"\t"+localDeltaG[miRNAstart-1]+"\t"+localDeltaG[miRNAstart]+"\t"+localDeltaG[miRNAstart+1]+"\t"+localDeltaG[miRNAstart+2]+"\t"+localDeltaG[miRNAstart+3]+"\t");
		if(miRNAstarStart > 0)
			System.out.print(localDeltaG[miRNAstarStart-3]+"\t"+localDeltaG[miRNAstarStart-2]+"\t"+localDeltaG[miRNAstarStart-1]+"\t"+localDeltaG[miRNAstarStart]+"\t"+localDeltaG[miRNAstarStart+1]+"\t"+localDeltaG[miRNAstarStart+2]+"\t"+localDeltaG[miRNAstart+3]+"\t");
		else
			System.out.print(" \t \t \t \t \t");
			
		System.out.print(sum+"\t");
		System.out.print(mean+"\t");
		System.out.println(std);
		}
		
	}
}
//>ddi-mir-1176 MI0006244 miRNAstart=21 miRNAstarStart=83
//External loop                           :  -110
//Interior loop (  1,119) CG; (  2,118) CG:  -330
//Interior loop (  2,118) CG; (  3,117) CG:  -330
//Interior loop (  3,117) CG; (  4,116) AU:  -210
//Interior loop (  4,116) AU; (  7,115) GC:   330
//Interior loop (  7,115) GC; (  8,114) UA:  -220
//Interior loop (  8,114) UA; (  9,113) CG:  -240
//Interior loop (  9,113) CG; ( 10,112) GC:  -240
//Interior loop ( 10,112) GC; ( 11,111) UA:  -220
//Interior loop ( 11,111) UA; ( 12,110) AU:  -130
//Interior loop ( 12,110) AU; ( 13,109) UA:  -110
//Interior loop ( 13,109) UA; ( 14,108) CG:  -240
//Interior loop ( 14,108) CG; ( 15,107) AU:  -210
//Interior loop ( 15,107) AU; ( 16,106) GC:  -210
//Interior loop ( 16,106) GC; ( 17,105) GC:  -330
//Interior loop ( 17,105) GC; ( 18,104) UA:  -220
//Interior loop ( 18,104) UA; ( 19,103) GC:  -210
//Interior loop ( 19,103) GC; ( 20,102) GC:  -330
//Interior loop ( 20,102) GC; ( 22,100) CG:    40
//Interior loop ( 22,100) CG; ( 23, 99) AU:  -210
//Interior loop ( 23, 99) AU; ( 24, 98) AU:   -90
//Interior loop ( 24, 98) AU; ( 25, 97) UA:  -110
//Interior loop ( 25, 97) UA; ( 26, 96) UA:   -90
//Interior loop ( 26, 96) UA; ( 27, 95) UA:   -90
//Interior loop ( 27, 95) UA; ( 28, 94) UA:   -90
//Interior loop ( 28, 94) UA; ( 29, 93) UA:   -90
//Interior loop ( 29, 93) UA; ( 31, 91) UA:   170
//Interior loop ( 31, 91) UA; ( 32, 90) CG:  -240
//Interior loop ( 32, 90) CG; ( 33, 89) AU:  -210
//Interior loop ( 33, 89) AU; ( 34, 88) AU:   -90
//Interior loop ( 34, 88) AU; ( 35, 87) GC:  -210
//Interior loop ( 35, 87) GC; ( 36, 86) GC:  -330
//Interior loop ( 36, 86) GC; ( 37, 85) AU:  -240
//Interior loop ( 37, 85) AU; ( 38, 84) AU:   -90
//Interior loop ( 38, 84) AU; ( 40, 82) GC:   110
//Interior loop ( 40, 82) GC; ( 41, 81) CG:  -340
//Interior loop ( 41, 81) CG; ( 42, 80) UA:  -210
//Interior loop ( 42, 80) UA; ( 43, 79) GC:  -210
//Interior loop ( 43, 79) GC; ( 45, 77) AU:   110
//Interior loop ( 45, 77) AU; ( 46, 76) UA:  -110
//Interior loop ( 46, 76) UA; ( 47, 75) CG:  -240
//Interior loop ( 47, 75) CG; ( 48, 74) AU:  -210
//Interior loop ( 48, 74) AU; ( 49, 73) UA:  -110
//Interior loop ( 49, 73) UA; ( 50, 72) CG:  -240
//Interior loop ( 50, 72) CG; ( 51, 71) AU:  -210
//Interior loop ( 51, 71) AU; ( 52, 70) AU:   -90
//Interior loop ( 52, 70) AU; ( 53, 69) GC:  -210
//Interior loop ( 53, 69) GC; ( 55, 66) GU:   210
//Interior loop ( 55, 66) GU; ( 56, 65) CG:  -250
//Hairpin  loop ( 56, 65) CG              :   480
//CCCAAAGUCGUAUCAGGUGGCCAAUUUUUAUCAAGGAAAGCUGUAUCAUCAAGUGCCUCCCUCUGUCUCUUGAUGAUUCAGCCUUCCUUGACAAAAAUUGCCCACCUGAUACGACUGGGA
//((((..((((((((((((((.((((((((.((((((((.((((.(((((((((.((........))..))))))))).)))).)))))))).)))))))).)))))))))))))))))). (-70.50)
