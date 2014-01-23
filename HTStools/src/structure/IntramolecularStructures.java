package structure;

import general.Functions;
import general.RNAfunctions;
import general.IOTools;
import general.energy;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Stack;

import general.ExtendedReader;
import general.ExtendedWriter;


public class IntramolecularStructures  implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private Structure []structures;
	
	
	public double[][] getAllStructures(){
		double[][] allSubStructures = new double[structures.length][];
		int[][] locations = new int[structures.length][];
		double[][] IASs = new double[structures.length][];
		
		for(int j = 0 ; j < structures.length; j++)
			allSubStructures[j] = new double[structures.length-j];
		
		for(int j = 0 ; j < structures.length; j++){
			if(structures[j]!= null){
				locations[j] = structures[j].getOtherLocations(locations[j]);
				IASs[j] = structures[j].getIAS(IASs[j]);
				double sum = 0;
				for(int k = 0; k < locations[j].length; k++){
					if(locations[j][k] > j){
						sum += IASs[j][k];
					}
				}
				for(int i = j ; i < structures.length; i++)
					allSubStructures[0][i] += sum;
			}
		}
		
		for(int j = 1 ; j < structures.length; j++){
			for(int k = 0; k < allSubStructures[j].length; k++){
				allSubStructures[j][k] = allSubStructures[j-1][k+1];
			}
			double sum = 0;
			
			if(locations[j-1]!= null){
				for(int k = 0; k < locations[j-1].length; k++){
					if(locations[j-1][k] >= j){
						for(int l = 0; l < locations[j-1][k]-j; l++){
							allSubStructures[j][l] = allSubStructures[j][l] - IASs[j-1][k];
						}
					}
					else{
						sum += IASs[j-1][k];
					}
				}
				if(sum > 0){
					for(int l = 0; l < allSubStructures[j].length; l++){
						allSubStructures[j][l] = allSubStructures[j][l] - sum;
					}
				}
				
			}
		}
		
		return allSubStructures;
	}
	
	
	
	public IntramolecularStructures getSpecificStructure(int[] sequence){
		IntramolecularStructures focalStructure = this.copy();
		int pointer = 0;
		for(int i = 0; i < this.structures.length; i++){
			if(sequence[i] == 0){
				focalStructure.removeStructureTotaly(pointer);
			}
			else
				pointer++;
		}
		return focalStructure;
	}

	public static void main(String []args){
		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
		}
		Hashtable<String,String> T = Functions.parseCommandLine(args);
		String Dir = ".";
		if(T.containsKey("-d"))
			Dir = (String)T.get("-d");
		
		String FileName = (String)T.get("-f");
		IntramolecularStructures IS = null;
		
		String outputFile = FileName + "_prepencity.txt";
		if(T.containsKey("-o"))
			outputFile = (String)T.get("-o");
		
		
		switch(((String)T.get("-p")).charAt(0) ) {
		case 'R':
		case 'r':
			@SuppressWarnings("unused")
			int width = 100;

			if( T.containsKey("-w"))
				width = Functions.String2Int((String)T.get("-w"));

			IS = getRfoldProbabilities(FileName,Dir);
			break;
		case 'c':
		case 'C':

			width = 100;

			if( T.containsKey("-w"))
				width = Functions.String2Int((String)T.get("-w"));

			IS = getRfoldStructure(FileName,Dir);
			break;
		case 'D':
		case 'd':
			IS = getDotBracketStructure(FileName,Dir);

		}
		
		if(IS == null)
			return;
		
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(Dir+"/"+outputFile));
			IS.printPrepencity(EW);
		}
		catch(Exception E){
				E.printStackTrace();
		}
		
			
	}
		
	
	public IntramolecularStructures(){}
	
	public IntramolecularStructures(int length){
		structures = new Structure[length];
		for(int i = 0; i < length;i++){
				this.structures[i] = new Structure();
		 }

		
	}
	
	public IntramolecularStructures(Structure[] structures){
		this.structures = structures;
		
	}
	
	
	public IntramolecularStructures copy(){
		Structure[] newStructures = new Structure[structures.length];
		for(int i = 0; i < structures.length; i ++ ){
			if(structures[i] != null){
				newStructures[i]= structures[i].copy();
			}
		}
		
		return new IntramolecularStructures(newStructures);
	}

	

	public IntramolecularStructures getRNAplfoldStructure(String FileName,String Dir, int genomeSequenceLength){
		try{
			ExtendedReader PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName+".ps"));
			//skipFAline
			String Info = PTTFileReader.readLine();
			System.out.println(Info);
			
			//> 1:232-33953 + E:\align\microAlignment\EColik12.gbk
			int nrOfRows = 0;	
			while(Info.indexOf("/winSize") < 0 ){
				PTTFileReader.skipLine();
				Info = PTTFileReader.readWord();
				//System.out.println(Info);
			}
			int winLength = PTTFileReader.readInt();

			while(PTTFileReader.readLine().compareTo("drawgrid_turn") != 0){
				nrOfRows = 0;
			}
			
			IntramolecularStructures tempStructure  = new IntramolecularStructures(genomeSequenceLength);
			while(PTTFileReader.lookAhead() != 's'){
				int location1 = PTTFileReader.readInt();
				int location2 = PTTFileReader.readInt();
				float prob = (float)PTTFileReader.readDouble();
				PTTFileReader.skipLine();

				addStructure(location1,location2,prob);
			}
			PTTFileReader.close();
			return tempStructure;
		}
		catch(Exception E){
			E.printStackTrace();
		}
		System.out.println("Something wrong with extracting the structure");
		return null;
		
	}

	public static IntramolecularStructures getRfoldProbabilities(String FileName,String Dir){
		int sequenceLength = getRfoldProbabilityLength(FileName, Dir);
		return getRfoldProbabilities(FileName,Dir,sequenceLength);
	}
	
	
	public static IntramolecularStructures getRfoldProbabilities(String FileName,String Dir, int sequenceLength){
		try{
			ExtendedReader PTTFileReader = null;
			if(FileName.endsWith(".rfold"))
				PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName));
			else
				 PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName+".rfold"));
			//skipFAline
			PTTFileReader.more();
			PTTFileReader.skipLine();

			//	String Info = PTTFileReader.readLine();
			//System.out.println(Info);
			
			IntramolecularStructures tempStructure  = new IntramolecularStructures(sequenceLength);
			while(PTTFileReader.more())
				tempStructure.addRfoldProbability(PTTFileReader);			
			PTTFileReader.close();
			Structure[] tempStructures = tempStructure.getStructures();
			
			tempStructure.setStructures(tempStructures);
			return tempStructure;
		}
		catch(Exception E){
			E.printStackTrace();
		}
		System.out.println("Something wrong with extracting the structure");
		return null;
	}
	
	
	
	public void addRfoldProbability( ExtendedReader PTTFileReader){
		int location1 = PTTFileReader.readInt();
		if(location1 < structures.length){
			int location2 = PTTFileReader.readInt();
			if(location1 != location2){
				float prob = (float)PTTFileReader.readDouble();
				PTTFileReader.skipLine();
				addStructure(location1-1, location2-1, prob);
			}
			else
				PTTFileReader.skipLine();
		}
		else
			PTTFileReader.skipLine();
	}
	
	public void setWeight(energy A){
		for(int i = 0; i < structures.length; i++){
			if(structures[i] != null && structures[i].getKindOfStructure() > 0){
				structures[i].weightStructure(i,A);
			}
		}
	}
	
	public static IntramolecularStructures generateContraFoldStructure(String structureDir, String sequenceDir,String fileName){
		if(IOTools.fileExists(structureDir, fileName+".CONTRAfold.bpseq"))
			IOTools.deleteFile(structureDir, fileName+".CONTRAfold.bpseq");
			RNAfunctions.generateCONTRAfoldStructure(sequenceDir, fileName, "fa", structureDir,-1);
		return 	getContraFoldStructure(structureDir,fileName+".CONTRAfold.bpseq",-1);
	}
	
	public static IntramolecularStructures generateContraFoldStructure(String structureDir, String sequenceDir,String fileName, IntramolecularStructures A){
		if(IOTools.fileExists(structureDir, fileName+".CONTRAfold.bpseq"))
			IOTools.deleteFile(structureDir, fileName+".CONTRAfold.bpseq");
			RNAfunctions.generateCONTRAfoldStructure(sequenceDir, fileName, "fa", structureDir,-1);
		return 	getContraFoldStructure(structureDir,fileName+".CONTRAfold.bpseq",-1);
	}

	
	
	public static IntramolecularStructures generateRfoldStructure(String structureDir, String sequenceDir,String fileName,int width){
		if(!IOTools.fileExists(structureDir, fileName+"_"+width+".rfold_structure"))
			RNAfunctions.generateRfoldStructures(sequenceDir, fileName, "fa", structureDir,width);
		return 	getRfoldStructure(structureDir,fileName+"_"+width);
	}
	
	public static IntramolecularStructures getRfoldStructure(String Dir,String FileName){
		int sequenceLength = getRfoldStructureLength(FileName, Dir);
		return getRfoldStructure(FileName,Dir,sequenceLength);
	}
	

	public static IntramolecularStructures getRfoldStructure(String FileName,String Dir, int sequenceLength){
		
		try{
			ExtendedReader PTTFileReader = null;
			if(FileName.endsWith(".rfold_structure"))
				PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName));
			else
				 PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName+".rfold_structure"));
			//skipFAline
			PTTFileReader.more();
			PTTFileReader.skipLine();

//				String Info = PTTFileReader.readLine();
//			System.out.println(Info);
			
			IntramolecularStructures tempStructure  = new IntramolecularStructures(sequenceLength);
			int structureStart = 0;
			while(PTTFileReader.more() )
				structureStart = tempStructure.addRfoldStructure(PTTFileReader);			
			PTTFileReader.close();
			Structure[] tempStructures =tempStructure.getStructures();
			 for(int i = 0; i < sequenceLength;i++){
				 if(tempStructures[i] == null){
					 tempStructures[i] = new Structure();
					 tempStructures[i].setBackbone();
				 }
			 }
			tempStructure.setStructures(tempStructures);
			
			return tempStructure;
		}
		catch(Exception E){
			E.printStackTrace();
		}
		
		System.out.println("Something wrong with extracting the structure");
		return null;
	}


	public static int getRfoldStructureLength(String FileName,String Dir){
		
		try{
			ExtendedReader PTTFileReader = null;
			if(FileName.endsWith(".rfold_structure"))
				PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName));
			else
				 PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName+".rfold_structure"));
			//skipFAline
			PTTFileReader.more();
			PTTFileReader.skipLine();

			//	String Info = PTTFileReader.readLine();
			//System.out.println(Info);
			
			int structureEnd = 0;
			while(PTTFileReader.more() )
				structureEnd = getRfoldStructureEnd(PTTFileReader);			
			PTTFileReader.close();
			
			return structureEnd;
		}
		catch(Exception E){
			E.printStackTrace();
		}
		
		System.out.println("Something wrong with extracting the structure");
		return 0;
	}
	
	
	public IntramolecularStructures getRfoldSubStructure(int start, int stop){

		Structure[] tempStructures = this.getStructures();
		
		int structureStart = start,pointer = start;
		
		while(pointer > -1 && tempStructures.length > pointer){
			if(tempStructures[pointer].isBackbone()){
				structureStart = pointer;
				pointer = -1;
			}
			pointer--;
		}
		int structureEnd = pointer = stop;
		
		while(tempStructures.length > pointer){
			if(tempStructures[pointer].isBackbone() ){
				structureEnd = pointer;
				pointer = tempStructures.length;
			}
			pointer++;
		}
		
		Structure[] newStructures = new Structure[structureEnd-structureStart];
		for(int i = 0; i < newStructures.length; i++){
			newStructures[i] = tempStructures[structureStart+i];
			newStructures[i].setLocation2(newStructures[i].getLocation2() - structureStart);
		}
		IntramolecularStructures newStructure = new IntramolecularStructures();
		newStructure.setStructures(newStructures);
		return newStructure;
	}

	public IntramolecularStructures getSubStructure(int start, int stop){

		Structure[] tempStructures = this.getStructures();
		
		
		Structure[] newStructures = new Structure[stop-start];
		for(int i = 0; i < newStructures.length; i++){
			newStructures[i] = tempStructures[start+i];
			//newStructures[i].setNewLocation(start);
		}
		IntramolecularStructures newStructure = new IntramolecularStructures();
		newStructure.setStructures(newStructures);
		return newStructure;
	}
	
	
	public int getRfoldSubStructureStart(int start){

		Structure[] tempStructures = this.getStructures();
		
		int pointer = start;
		while(pointer > -1 && tempStructures.length > pointer){
			if(tempStructures[pointer].isBackbone()){
				return pointer;
			}
			pointer--;
		}
		return 0;
	}

	public void addIntramolecularStructure(IntramolecularStructures newStructure){
		int thisLength = this.getLength();
		int otherLength = newStructure.getLength();
		
		Structure[] totalStructure = new Structure[thisLength+otherLength];
		for(int i = 0; i < thisLength; i++){
			totalStructure[i] = this.structures[i];
		}
		for(int i = 0; i < otherLength; i++){
			totalStructure[i+thisLength] = newStructure.structures[i].copy();
			totalStructure[i+thisLength].shift(thisLength);
		}
		this.structures = totalStructure;
		
	}
	
	public void addSingleBackbone(){
		int thisLength = this.getLength();
		
		Structure[] totalStructure = new Structure[thisLength+1];
		for(int i = 0; i < thisLength; i++){
			totalStructure[i] = this.structures[i];
		}
		Structure sb = new Structure();
		sb.setBackbone();
		totalStructure[thisLength] = sb;
		this.structures = totalStructure;
	}

	public void addSingleBasepair(){
		int thisLength = this.getLength();
		
		Structure[] totalStructure = new Structure[thisLength+2];
		for(int i = 0; i < thisLength; i++){
			totalStructure[i] = this.structures[i];
		}
		totalStructure[thisLength] = new Structure(thisLength+1,1,1);;
		totalStructure[thisLength+1] = new Structure(thisLength,1,1);
		
		this.structures = totalStructure;
	}
	
	public void addSingleStrandedRegion(int length){
		int thisLength = this.getLength();
		
		Structure[] totalStructure = new Structure[thisLength+length];
		for(int i = 0; i < thisLength; i++){
			totalStructure[i] = this.structures[i];
		}
		for(int i = 0; i< length; i++){
			totalStructure[i+thisLength] = new Structure();
		}
		
		this.structures = totalStructure;
	}
	
	
	public int addRfoldStructure( ExtendedReader PTTFileReader){
		int start = PTTFileReader.readInt();
		PTTFileReader.more();
		int end = PTTFileReader.readInt();
		
		PTTFileReader.more();
		PTTFileReader.readDouble();
		PTTFileReader.more();
		char[] structure = PTTFileReader.readWord().toCharArray();
		Stack<Integer> stack = new Stack<Integer>();
		for(int i = 0; i < structure.length; i++){
			int pos = i+start-1;
			if(structure[i] == '<')
				stack.push(pos);
			else if(structure[i] == '>'){
				int pos1 = stack.pop();
				int pos2 = pos;
				this.addStructure(pos1,pos2,1);
			}
			else{
				structures[pos] = new Structure();
				structures[pos].setLoop();
			}
		}
		PTTFileReader.more();
		return end;
	}
		
		
	public static int getRfoldStructureEnd( ExtendedReader PTTFileReader){
		PTTFileReader.readInt();
		PTTFileReader.more();
		int end = PTTFileReader.readInt();
		
		PTTFileReader.skipLine();
		return end;
	}
	
	private static int getRfoldProbabilityLength(String FileName,String Dir){	
		try{
			ExtendedReader PTTFileReader = null;
			if(FileName.endsWith(".rfold"))
				PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName));
			else
				 PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName+".rfold"));
				//skipFAline
			PTTFileReader.more();
//			String Info = PTTFileReader.readLine();
			PTTFileReader.readLine();
			int max = 0;
			while(PTTFileReader.more() ){
				PTTFileReader.readInt();
//				int location1 = PTTFileReader.readInt();
				int location2 = PTTFileReader.readInt();
				PTTFileReader.skipLine();
				if(location2 > max)
					max = location2;
			}
			return max;
		}
		catch(Exception E){
			E.printStackTrace();
		}
		return 0;
	}
	
	
	public void  changeDirection(){
		Structure [] newStructures = new Structure[structures.length];
		for(int i = 0; i < structures.length; i++){
			if(structures[i]!= null){
//S				System.out.println(i+"\t"+structures[i].printLocations(""));
				structures[i].changeDir(structures.length);
				newStructures[structures.length-i-1] = structures[i];
			}
		}
/*		for(int i = 0; i < newStructures.length; i++){
			if(newStructures[i]!= null){
				System.out.println(i+"\t"+newStructures[i].printLocations(""));
			}
		}
	*/	
		
		structures = newStructures;
	}
	
	
	public static IntramolecularStructures getContraFoldStructure(String FileName,String Dir, int sequenceLength,double cutoff){
		//System.out.println(Dir);
		try{
			ExtendedReader PTTFileReader = null;
			if(cutoff > -1){
				if(IOTools.fileExists(Dir+"/"+FileName+".CONTRAfold.posteriors"))
					PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName+".CONTRAfold.posteriors"));
				else if(IOTools.fileExists(Dir+"/"+FileName+".fa"))
					PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName+".fa"));
				else
					PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName));
					
			}
			else{
				if(IOTools.fileExists(Dir+"/"+FileName+".bpseq"))
					PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName+".bpseq"));
				else if(IOTools.fileExists(Dir+"/"+FileName+".fa"))
					PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName+".fa"));
				else
					PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName));
			}
			
			//skipFAline
			PTTFileReader.more();
			IntramolecularStructures tempStructure  = new IntramolecularStructures(sequenceLength);
			while(PTTFileReader.more())
				tempStructure.addContraFoldStructure(PTTFileReader);			
			PTTFileReader.close();
			return tempStructure;

		}
		catch(Exception E){
			E.printStackTrace();
		}
		
		System.out.println("Something wrong with extracting the structure");
		return null;
	}
	
	public static IntramolecularStructures getContraFoldStructure(String FileName,String Dir){
		int SequenceLength = getContraFoldSequenceLength(FileName,Dir);
		return getContraFoldStructure(FileName,Dir, SequenceLength,0.05);
	}

	public static IntramolecularStructures getContraFoldStructure(String Dir,String FileName, double cutoff){
		int SequenceLength = getContraFoldSequenceLength(FileName,Dir);
		return getContraFoldStructure(FileName,Dir, SequenceLength,cutoff);
	}
	
	
	private static int getContraFoldSequenceLength(String FileName, String Dir){
		try{
			ExtendedReader PTTFileReader = null;
			PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName));
			//skipFAline
			PTTFileReader.more();
			int max = 0;
			while(PTTFileReader.more()){
				int location1 = PTTFileReader.readInt();
				if(location1 > max) 
					max = location1;
				PTTFileReader.readWord();
				//String nucleotide = PTTFileReader.readWord();
				String info = PTTFileReader.readLine();
				info = info.trim();

				if(info.length() != 0){
					//System.out.println(location1+"   "+nucleotide+"  "+info);
					if(info.indexOf(":") == -1){
						int location2 = Integer.parseInt(info);
						if(location2 > max )
							max = location2;
					}
					else{
						while(info.indexOf(":") != -1){
							info= info.trim();
							String loc = info.substring(0,info.indexOf(":"));
							info = info.substring(info.indexOf(":")+1);
//							String pr = null;
//							if(info.indexOf(":") == -1)
//								pr = info;
							if(info.indexOf(":") != -1){
								//pr = info.substring(0,info.indexOf(" "));
								info = info.substring(info.indexOf(" ")+1);
							}
							int location2 = Integer.parseInt(loc);
							if(location2 > max )
								max = location2;
						}
					}	
				}

			}
			PTTFileReader.close();
			return max;

		}
		catch(Exception E){
			E.printStackTrace();
		}

		System.out.println("Something wrong with extracting the structure");
		return 0;
	}

	
	public void addContraFoldStructure(ExtendedReader PTTFileReader){
		int location1 = PTTFileReader.readInt();
		if(location1 < structures.length){
			//String nucleotide = PTTFileReader.readWord();
			PTTFileReader.readWord();
			String info = PTTFileReader.readLine();
			info = info.trim();
			
			if(info.length() != 0){
				//System.out.println(location1+"   "+nucleotide+"  "+info);
				if(info.indexOf(":") == -1){
					int location2 = Integer.parseInt(info);
					//System.out.println(location2);
					float prob = 1;
					if(location2 > 0)
						addStructure(location1-1, location2-1, prob);			
				}
				else{
					while(info.indexOf(":") != -1){
						info= info.trim();
						String loc = info.substring(0,info.indexOf(":"));
						info = info.substring(info.indexOf(":")+1);
						String pr = null;
						if(info.indexOf(":") == -1)
							pr = info;
							else{
								pr = info.substring(0,info.indexOf(" "));
								info = info.substring(info.indexOf(" ")+1);
							}float prob = (float)Double.parseDouble(pr);
							int location2 = Integer.parseInt(loc);
							//System.out.println(loc+"   "+pr+"testing");
							//System.out.println(location2+"   "+prob+"testing");
							addStructure(location1-1, location2-1, prob);
					}
				}	
			}
		}
		else
			PTTFileReader.skipLine();
	}
	
	
	
	
	public static IntramolecularStructures getDotBracketStructure(String FileName,String Dir, int sequenceLength){
		try{
			ExtendedReader PTTFileReader = null;
			PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName));
			//skipFAline
			while(PTTFileReader.more()){
				String word = PTTFileReader.readWord();
				char[] array = word.toCharArray();
				if(array.length == sequenceLength){
					if(analyseDotBracket(array))
						return addDotBracketAnnotation(array);
				}
			}
			PTTFileReader.close();
			

		}
		catch(Exception E){
			E.printStackTrace();
		}
		
		System.out.println("Something wrong with extracting the structure");
		return null;
	}

	public static IntramolecularStructures getDotBracketStructure(String FileName,String Dir){
		ExtendedReader PTTFileReader = null;
		try{
			PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName));
			//skipFAline
			while(PTTFileReader.more()){
				String word = PTTFileReader.readWord();
				char[] array = word.toCharArray();
				if(analyseDotBracket(array)){
					PTTFileReader.close();
					return addDotBracketAnnotation(array);
				}
			}
			PTTFileReader.close();
		}
		catch(Exception E){
			E.printStackTrace();
		}
		finally {
			try{
				PTTFileReader.close();
			}catch(Exception e){e.printStackTrace();}
		}

		System.out.println("Something wrong with extracting the structure");
		return null;
	}

	
	
	
	
	public static IntramolecularStructures addDotBracketAnnotation(char[] Pattern){
		IntramolecularStructures tempStructure  = new IntramolecularStructures(Pattern.length);
		//int[] RNAStructureVector = null;
		Stack<Integer> stack = new Stack<Integer>();
		for(int i = 0;i < Pattern.length; i++){
			if(Pattern[i] == ')' || Pattern[i] == '>'){
				try{
					int location1 = stack.pop();
					int location2 = i;
					float prob = 1;
					tempStructure.addStructure(location1, location2, prob);
				}
				catch(Exception E){
					System.out.println("Something wrong with the structure at location:" + i);
				}
			}
			else if(Pattern[i] == '(' || Pattern[i] == '<'){
					stack.push(i);
			}
		}
		return tempStructure;
	}
	
	private static boolean analyseDotBracket(char[] array){
		int count = 0;
		for(int i = 0; i < array.length; i++){
			if(array[i] == '.')
				count++;
			else if(array[i] == '(')
				count++;
			else if(array[i] == ')')
				count++;
			else if(array[i] == ']')
				count++;
			else if(array[i] == '}')
				count++;
			else if(array[i] == '[')
				count++;
			else if(array[i] == '{')
				count++;
		}
		if(count == array.length)
			return true;
		else
			return false;
	}

	
	
	public void convertStructure (double cutoff){
		for(int i = 0; i < structures.length;i++){
			if(structures[i] != null){
				if(structures[i].getProbability() < cutoff)
					structures[i] = structures[i].removeStructure((float)cutoff);
			}
		}
		
	}

	public boolean isBasepair(int location1, int location2, double cutoff){
		if(structures[location1] == null || structures[location2] == null) return false;
		else return structures[location1].isBasepairing(location2,cutoff);
	}

	public void addStructure(int location1, int location2, float probability){
		try{
		Structure S = new Structure(location2 , probability);
		
		if(structures[location1] == null|| structures[location1].getKindOfStructure() < 0)
			structures[location1] = S;
		else{
			structures[location1].addStructure(S);
		}
		Structure S2 = new Structure(location1 , probability);
		if(location2 > -1 && location2 < structures.length){ 
			if(location1 != location2){
				if(structures[location2] == null || structures[location2].getKindOfStructure() < 0)
					structures[location2] = S2;
				else
					structures[location2].addStructure(S2);
			
			}
		}
		}catch(Exception E){E.printStackTrace();}
	}

	public void addStructure(int location1, int location2, float probability, int kindOfStructure){
		try{
		Structure S = new Structure(location2 , probability,kindOfStructure);
		
		if(structures[location1] == null|| structures[location1].getKindOfStructure() < 0)
			structures[location1] = S;
		else{
			structures[location1].addStructure(S);
		}
		Structure S2 = new Structure(location1 , probability,kindOfStructure);
		if(location2 > -1 && location2 < structures.length){ 
			if(location1 != location2){
				if(structures[location2] == null || structures[location2].getKindOfStructure() < 0)
					structures[location2] = S2;
				else
					structures[location2].addStructure(S2);
			
			}
		}
		}catch(Exception E){E.printStackTrace();}
	}
	
	
	
	

	public void addKnot(int location1, int location2,int knotNR){
		try{
		Structure S = new Structure(location2 , 1,Structure.knot, knotNR);
		
		if(structures[location1] == null|| structures[location1].getKindOfStructure() < 0)
			structures[location1] = S;
		else{
			structures[location1].addStructure(S);
		}
		Structure S2 = new Structure(location1 , 1,Structure.knot,knotNR);
		if(location2 > -1 && location2 < structures.length){ 
			if(location1 != location2){
				if(structures[location2] == null || structures[location2].getKindOfStructure() < 0)
					structures[location2] = S2;
				else
					structures[location2].addStructure(S2);
			
			}
		}
		}catch(Exception E){E.printStackTrace();}
	}
	
	
	
	
	
	public boolean isBackbone(int location1, int location2){
		int[] locations = null;
		locations = structures[location1].getInteractions(locations,location1);
		if(locations!= null)
			return false;
		locations = structures[location2].getInteractions(locations,location2);
		if(locations!= null)
			return false;
		
		
		for(int i = location1+1; i < location2 ; i++){
			if(structures[i] != null){
				if(structures[i].getKindOfStructure() ==1 ){
					int temp = structures[i].getLocation2();
					if(temp < location1 || temp > location2)
						return false;
					else 
						i= temp;
				}	
				else if(structures[i].getKindOfStructure() != -1){
					locations = null;
					locations = structures[i].getInteractions(locations,i);
					if(locations != null){
						for(int j = 0; j < locations.length; j++){
							if(locations[j] <= location1 || locations[j] >=location2)
								return false;
						}
					}
				}
			}
		}
		return true;

	}

	
	public void changeStructure(int location1, int location2, float probability){
		this.structures[location1] = new Structure(location2 , probability);
		this.structures[location2] = new Structure(location1 , probability);

	}
	
	
	public void setStructure(int location1, int location2, float probability){
		this.addStructure(location1, location2, probability);
	}
	
	public void setStructure(int location1, int location2, float probability,int kindOfStructure){
		this.addStructure(location1, location2, probability,kindOfStructure);
	}

	
	
	public boolean isKnot(int location1, int location2){
		int[] locations = null;
		if(structures[location1].getInteractions(locations,location1).length > 1 || structures[location1].getKindOfStructure() == 2)
			return true;
		if(structures[location2].getInteractions(locations, location2).length > 1 || structures[location2].getKindOfStructure() ==2)
			return true;
		return false;		
	}
		
		

	
	

	public char[] getDotBracketAnnotation(){
		char[] dotBracketStructures = new char[structures.length];
		for(int i = 0; i < structures.length; i++){
			if(structures[i]!= null){
				dotBracketStructures[i] = structures[i].getDotBracketStructure(i);
				if(dotBracketStructures[i] == (char)'!')
					System.out.println(i);
			}
			else
				System.out.println(i);
		}
		return dotBracketStructures;
		
	}
	

	
	public int[][] getIntermolecularStructures(String Sequence){
		int length = Sequence.indexOf("&");
		int[] mRNAstructure = new int[length];
		int[] sRNAstructure = new int[this.structures.length- length - 2];
		for(int i = 0; i < structures.length; i++){
			if(i < length){
				if(structures[i].getLocation2() > length ){
					mRNAstructure[i] = 0;
				}
				else
					mRNAstructure[i] = -1;

			}
			else if(i > length + 1){
				if(structures[i].getLocation2() < length )
					sRNAstructure[i - length-2 ] = 0;
				else
					sRNAstructure[i - length-2 ] = -1;

			}
		}
/*		for(int i = 0; i < sRNAstructure.length; i++){
			System.out.println((i+1)+"\t"+Sequence.charAt(i+2+length)+"\t"+sRNAstructure[i]);
		}
		for(int i = 0; i < mRNAstructure.length; i++){
			System.out.println((i+1)+"\t"+Sequence.charAt(i)+"\t"+mRNAstructure[i]);
		}
*/
		
		int[][] structures = new int[2][];
		structures[0]= mRNAstructure;
		structures[1] = sRNAstructure;
		return structures;
		

	}
	
	
	public void printStructure(ExtendedWriter EW, char[] sequence){
		for(int i = 0; i < structures.length; i ++){
			EW.print((i+1)+"  "+sequence[i]);
			if(structures[i]!= null){
				structures[i].printIntraStructure(EW, i, sequence);
			}
			EW.println();
		}
	}

	public void printStructure(ExtendedWriter EW){
		for(int i = 0; i < structures.length; i ++){
			if(structures[i]!= null){
				structures[i].printIntraStructure(EW, i);
			}
		}
	}

	
	public void printStructures(int location){
			if(structures[location]!= null){
				structures[location].printIntraStructure(location);
			}
	}

	public void printPrepencity(int location){
			if(structures[location]!= null){
				structures[location].printIntraStructurePrepencity(location);
			}
	}
	
	public void printPrepencity(ExtendedWriter EW,int location){
		if(structures[location]!= null){
			float totalProb = structures[location].getStructureProbabilities(location);
			EW.println("   "+totalProb);
		}
		else
			EW.println("   "+0);
		
	}

	public void printPrepencity(ExtendedWriter EW){
		for(int i = 0; i < structures.length; i ++){
			if(structures[i]!= null){
				float totalProb = structures[i].getStructureProbabilities(i);
				EW.println(totalProb);
			}
			else
				EW.println(0);
		}
		EW.flush();
	}
	
	public void getPrepencity(ArrayList <Float> prepencity){
		for(int i = 0; i < structures.length; i ++){
			if(structures[i]!= null){
				prepencity.add(structures[i].getStructureProbabilities(i));
			}
			else
				prepencity.add((float)0);
		}
	}
	
	public double[] getPrepencity(){
		double[] prepencities = new double[structures.length];  
		for(int i = 0; i < structures.length; i ++){
			if(structures[i]!= null){
				prepencities[i] = structures[i].getStructureProbabilities(i);
			}
		}
		return prepencities;
	}
	
	public void printXMLPrepencity(ExtendedWriter EW){
		EW.println("<structures>");
		for(int i = 0; i < structures.length; i ++){
			EW.println("<structure>");
			EW.println("<location>"+ (i+1) + "</location>");
			if(structures[i]!= null){
				float totalProb = structures[i].getStructureProbabilities(i);
				EW.println("<prepencity>"+totalProb+"</prepencity>");
			}
			else
				EW.println("<prepencity>"+0+"</prepencity>");
		}
		EW.println("</structures>");
		EW.flush();
	}

	public void checkStructures(){
		int above = 0;
		for(int i = 0; i < structures.length; i ++){
			
			if(structures[i]!= null){
				
				float totalProb = structures[i].getStructureProbabilities(i);
				if(totalProb > 1.0 )
					above++;
				
			}
		}
		System.out.println("nr above 1: "+ above+ "of total: " +structures.length );
	}
	

	public boolean checkStructures(int[] sequence){
		try{
		int above = 0;
		for(int i = 0; i < structures.length; i ++){
			if(structures[i]!= null){
				int [] locations = null;
				locations = structures[i].getInteractions(locations,i);
				if(locations[0]> 0 && locations[0] < sequence.length && !RNAfunctions.isBasepair(sequence[i],sequence[locations[0]])){
					System.out.println("noncanonical basepair : "+ i +" - "+ locations[0]+  "     "+ sequence[i]+ "- "+sequence[locations[0]]);
					above++;
				
				}
			}
		}
		if(above > 0)	
			System.out.print("noncanonical basepairs: "+ above+ "of total: " +structures.length );
		else
			return true;
		return false;
		}
				catch(Exception E){
			E.printStackTrace();
		}
		return false;
		
	}


	
	public float getStructurePenalty(int location,int startLocation){
		if(structures[location] != null)
			return structures[location].getStructureProbabilities(startLocation,location,0);
		else
			return 0;
	}

	public float getStructurePenalty(int gapStart,int gapWidth,int initiationStart){
		float structurePenalty = 0; 
		for(int i = gapStart+gapWidth; i > gapStart-1; i--){
			if(structures[i] != null){
				structurePenalty += getStructurePenalty(i, initiationStart);
			}
		}
		return structurePenalty;
	}

	
	public void removeStructure(int RNAstart, int RNAstop){
		if(RNAstart >RNAstop){
			for(int i = RNAstart-1 ; i >= RNAstop; i--)
				removeStructure(i);
		}
		else if(RNAstart <RNAstop){
			for(int i = RNAstart+1 ; i <= RNAstop; i++)
				removeStructure(i);
		}
		
		
	}
	
	public void removeStructures(int RNAstart,int RNAstop){
		for(int i = RNAstart ; i <= RNAstop; i++){
			removeStructure(i);
		}
	}

	public void removeStructure(int location){
		if(structures[location] != null){
			int[] locations = null;
			locations = structures[location].getInteractions(locations, location);
			structures[location] = new Structure();
			reasign(location);
			if(locations != null){
				for(int j = 0 ; j < locations.length; j++){
					if(locations[j] >-1 && locations[j] < structures.length){
						structures[locations[j]] = structures[locations[j]].removeStructure(location);
						if(structures[locations[j]] == null){
							structures[locations[j]] = new Structure();
							reasign(locations[j]);
						}
					}
				}
			}
		}
	}

	public void removeStructureTotaly(int location){
		if(structures[location] != null){
			int[] locations = null;
			locations = structures[location].getInteractions(locations, location);
			structures[location] = new Structure();
			reasign(location);
			if(locations != null){
			for(int j = 0 ; j < locations.length; j++){
				if(locations[j] >-1 && locations[j] < structures.length){
					structures[locations[j]] = structures[locations[j]].removeStructure(location);
					if(structures[locations[j]] == null){
						structures[locations[j]] = new Structure();
						reasign(locations[j]);
					}
				}
			}
			}
		}
		Structure[] newStructures = new Structure[structures.length-1];
		for(int i = 0; i < structures.length; i++){
			if(i < location){
				newStructures[i] = structures[i];
				newStructures[i].changePointer(location);
			}
			if(i > location){
				newStructures[i-1] = structures[i];
				newStructures[i-1].changePointer(location);
			}
		}
		this.structures = newStructures;
	}
	
	
	
	public void reasign(int location){
			if(location != 0){
				if(this.structures[location-1].getKindOfStructure() == -1){
					this.structures[location].setBackbone();
					if(location+1 < structures.length && structures[location+1].getKindOfStructure() == -2)
						reasign(location+1);
					return;
				}
			}
			if(location != structures.length-1){
				if(this.structures[location+1].getKindOfStructure() == -1){
					this.structures[location].setBackbone();
					if(location-1 > -1 && structures[location-1].getKindOfStructure() == -2)
						reasign(location-1);
					return;
				}
				
			}
			structures[location].setLoop();
	}
	
	public void reasign(){
		for(int i = 0; i < this.structures.length-1;i++){
			if(this.structures[i].getKindOfStructure() > 0 && i != this.structures.length-1)
				i= setLoops(i+1,structures[i].getLocation2());
			else if(this.structures[i].getKindOfStructure() == -2)
				if(this.structures[i-1].getKindOfStructure() == -1 || this.structures[i+1].getKindOfStructure() == -1 )
					this.structures[i].setKindOfStructure(-1);
		}
	
	}

	
	
	private int setLoops(int start, int stop){
		for(int i = start;i < stop; i++){
			if(i> this.structures.length-1){
				System.out.println("fel i loops start"+ start+ " stop: "+stop );
			}
			if(this.structures[i].getKindOfStructure() > 0){
				if(structures[i].getLocation2() > stop){
					stop = structures[i].getLocation2();
				
				}
				if(structures[i].getLocation2() < start){
					i = setLoops(structures[i].getLocation2(),stop);
				}
			}
			else
				structures[i].setKindOfStructure(-2);
		}
		return stop;
	}
	
	public IntramolecularStructures extractStructures(int start, int stop){
		Structure[] newStructure = new Structure[stop-start];
		for(int i = start ; i < stop ;i ++ ){
			if(structures[i] != null && structures[i].getKindOfStructure() >0)
				newStructure[i-start] = structures[i].extractStructures(start);
			else
				newStructure[i-start] = structures[i];
		}
		return new IntramolecularStructures(newStructure);
	}

	
	public boolean isOpen(int location, double cutoff){
		if(structures[location] != null &&!structures[location].isOpen(location, cutoff))
			return false;
		return true;
	}


	/**
	 * @return the structures
	 */
	public Structure[] getStructures() {
		return structures;
	}


	/**
	 * @param structures the structures to set
	 */
	public void setStructures(Structure[] structures) {
		this.structures = structures;
	}
	

	/**
	 * @return the structures
	 */
	public Structure getStructure(int i) {
		if(structures[i]!= null)
			return structures[i];
		return null;
	}


	/**
	 * @param structures the structures to set
	 */
	public void setStructures(Structure structure, int i) {
		this.structures[i] = structure;
	}
	
	
	public int getLength(){
		if(this.structures != null)
			return structures.length;
		else
			return 0;
	}
	
	
}
