package Ontologies;

import general.ExtendedWriter;

public class PantherGene {
	
	String geneName;
	String PantherClass;
	double pValue;
	double score;
	String location;
	
	PantherClass PC;

//trinity_24533_1_1       PTHR23155               5.3e-10 44.1    3-145,

	PantherGene(){}
	
	public boolean addInfo(String info){
		String[] tabInfo = info.split("\t");
		if(tabInfo.length != 6) {
			System.out.print(tabInfo.length +"\t");
			return false;
		}
		this.geneName = tabInfo[0];
		this.PantherClass  = tabInfo[1];
		this.pValue = Double.parseDouble(tabInfo[3]);
		this.score = Double.parseDouble(tabInfo[4]);
		this.location  = tabInfo[5];
		return true;
	}
	
	public void printClass(ExtendedWriter EW){
		EW.print(this.geneName+"\t");
		PC.printInfo(EW);
	}

	public void printInfo(ExtendedWriter EW){
		PC.printInfo(EW);
	}
	public void printPantherClassInfo(ExtendedWriter EW){
		PC.printClassInfo(EW);
	}
	
	public boolean printPantherGOInfo(ExtendedWriter EW){
		if(PC.hasGOinfo()){
			PC.printGOInfo(EW);
			return true;
		}return false;
	}
	
	public boolean printOtherInfo(ExtendedWriter EW){
		if(PC.hasOtherInfo()){
			PC.printOtherInfo(EW);
			return true;
		}return false;
	}
	
	
}
