package Ontologies;

import java.util.ArrayList;

import general.ExtendedWriter;

public class GOGene{
	
	public String geneName;
	
	ArrayList <GOClass> GOclasses;

//trinity_24533_1_1       PTHR23155               5.3e-10 44.1    3-145,

	GOGene(){
		this.GOclasses = new ArrayList <GOClass>();
	}

	GOGene(String Name){
		this.geneName = Name;
		this.GOclasses = new ArrayList <GOClass>();
	}
	
	public boolean addGOterm(GOClass GOC){
		this.GOclasses.add(GOC);
		return true;
	}
	

}
