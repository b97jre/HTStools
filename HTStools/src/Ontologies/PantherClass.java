package Ontologies;

import sun.tools.tree.ThisExpression;
import general.ExtendedWriter;



//1) PANTHER ID:  for example, PTHR11258 or PTHR12213:SF6.  ":SF" indicates the subfamily ID
//2) Name:  The annotation assigned by curators to the PANTHER family or subfamily
//3) Molecular function*:  PANTHER GO slim molecular function terms assigned to families and subfamilies
//4) Biological process*:  PANTHER GO slim biological process terms assigned to families and subfamilies
//5) Cellular components*:  PANTHER GO slim cellular component terms assigned to families and subfamilies
//6) Protein class*  PANTHER protein class terms assigned to families and subfamilies
//7) Pathway***: PANTHER pathways have been assigned to families and subfamilies.  


public class PantherClass {
	String PANTHER_ID;
	String Name;
	String Mol_Fun;
	String Bio_Fun;
	String Cel_Fun;
	String Prot_Class;
	String Pathway;


	public PantherClass(){}
	
	public boolean addInfo(String info){
		String[] tabInfo = info.split("\t",7);
		if(tabInfo.length != 7){ 
			System.out.print(tabInfo.length+"\t");
			return false;
		}
		this.PANTHER_ID = tabInfo[0];
		this.Name = tabInfo[1];
		this.Mol_Fun = tabInfo[2];
		this.Bio_Fun = tabInfo[3];
		this.Cel_Fun = tabInfo[4];
		this.Prot_Class = tabInfo[5];
		this.Pathway = tabInfo[6];
		return true;
	}
	
	boolean hasGOinfo(){
		int totalLength = this.Mol_Fun.length()+this.Bio_Fun.length()+this.Cel_Fun.length();
		if(totalLength == 0) return false;
		return true;
	}
	
	boolean hasOtherInfo(){
		int totalLength = this.Prot_Class.length() + this.Pathway.length();
		if(totalLength == 0) return false;
		return true;
	}
	
	public void printInfo(ExtendedWriter EW){
		EW.println(this.PANTHER_ID+"\t"+this.Name+"\t-\t"+this.Mol_Fun+"\t"+this.Bio_Fun+"\t"+this.Cel_Fun+"\t"+this.Prot_Class+"\t"+this.Pathway);
	}

	public void printClassInfo(ExtendedWriter EW){
		EW.print(this.PANTHER_ID+"\t"+this.Name);
	}
	
	public void printGOInfo(ExtendedWriter EW){
		EW.print(this.Mol_Fun+"\t"+this.Bio_Fun+"\t"+this.Cel_Fun);
		
	}
	
	public void printOtherInfo(ExtendedWriter EW){
		EW.print(this.Prot_Class+"\t"+this.Pathway);
		
	}
	
	
	
	
}
