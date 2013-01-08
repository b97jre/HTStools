package structure;

import general.Functions;

import java.io.FileWriter;
import java.util.Hashtable;

import general.ExtendedReader;
import general.ExtendedWriter;




public class GenomeStructure {
	
	private IntramolecularStructures forwardStructure;
	private IntramolecularStructures reverseStructure;
	private int length;
	
	
	
	public static void main(String []args){
		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
		}
		Hashtable<String,String> T = Functions.parseCommandLine(args);
		if(T.containsKey("-f")) {
			

		String FileName = T.get("-f");
		
		String dir = getDir(T);
		
		GenomeStructure GS = new GenomeStructure(FileName,dir,4639676,"rfold");
		
		IntramolecularStructures newS = GS.extractStructures(986205+100,986205-100);
		//IntramolecularStructures newS2 = GS.extractStructures(80,70);
		try{
		ExtendedWriter EW = new ExtendedWriter(new FileWriter("testing.ps"));
//		newS.printStructure(EW);
	//	EW.flush();
		//EW.close();
		//EW = new ExtendedWriter(new FileWriter("testing2.ps"));
		newS.printPrepencity(EW);
		EW.flush();
		EW.close();
		}
		catch(Exception E){
			E.printStackTrace();
		}                                                                           
		}
		else{
			System.out.println("You have to give NCnumber alignmentFilename ");
			System.out.println("-n \t\t NCnumber of the focal genome");
			System.out.println("-f \t\t Name of alignmentfiles");
			System.out.println("-ptt \t\t if you want to use a ptt file for mRNA location");
		}
		
	}
	
	
	public GenomeStructure(String NCNumber, String Dir,int length, String Algorithm){
		
		this.length = length;
		System.out.println("Getting forward structure");
		if(Algorithm.indexOf("RNAplfold") == -1)
			forwardStructure = IntramolecularStructures.getRfoldProbabilities(NCNumber,Dir,length);
		System.out.println("Done");
		forwardStructure.checkStructures();
		System.out.println("Getting reverse structure");
		reverseStructure = IntramolecularStructures.getRfoldProbabilities(NCNumber+"_complementary",Dir,length);
		reverseStructure.checkStructures();
		System.out.println("Done");
		//reverseStructure = getStructure(NCNumber,Dir,length);
			
	}
	

	private static String getDir(Hashtable<String,String> T){
		if(T.containsKey("-d"))
			return T.get("-d");
		else
			return ".";
	}
	
	

	public GenomeStructure(String NCNumber, String Dir,int length, String Algorithm,boolean forward){
		
		this.length = length;
		if(forward){
		
			System.out.println("Getting forward structure");
			if(reverseStructure != null)
				reverseStructure = null;
			if(Algorithm.indexOf("RNAplfold") == -1)
				forwardStructure = IntramolecularStructures.getRfoldProbabilities(NCNumber,Dir,length);
			forwardStructure.checkStructures();
			System.out.println("Done");
		}
		else{
		System.out.println("Getting reverse structure");
			if(forwardStructure != null)
				forwardStructure = null;
			reverseStructure = IntramolecularStructures.getRfoldProbabilities(NCNumber+"_complementary",Dir,length);
			reverseStructure.checkStructures();
			System.out.println("Done");
		//reverseStructure = getStructure(NCNumber,Dir,length);
		}
	}

	
	
	
	public IntramolecularStructures extractStructures(int start, int stop){
		if(start < stop){
			if(forwardStructure != null)
				return forwardStructure.extractStructures(start,stop);
			return null;
		}
		else if(reverseStructure != null){
			
			start = length-start+1;
			stop = length- stop;
			//System.out.println(start+" "+stop);
			return reverseStructure.extractStructures(start,stop);
		}
		return null;
		
	}
	
	
	
}





