package Ontologies;

import general.ExtendedWriter;



//1) PANTHER ID:  for example, PTHR11258 or PTHR12213:SF6.  ":SF" indicates the subfamily ID
//2) Name:  The annotation assigned by curators to the PANTHER family or subfamily
//3) Molecular function*:  PANTHER GO slim molecular function terms assigned to families and subfamilies
//4) Biological process*:  PANTHER GO slim biological process terms assigned to families and subfamilies
//5) Cellular components*:  PANTHER GO slim cellular component terms assigned to families and subfamilies
//6) Protein class*  PANTHER protein class terms assigned to families and subfamilies
//7) Pathway***: PANTHER pathways have been assigned to families and subfamilies.  


public class GOClass {
	String GO_ID;
	String term;
	String Description;


	public GOClass(){}
	
	public GOClass(String ID,String term){
		this.GO_ID=ID;
		this.term=term;
	}
	
	public boolean addInfo(String info){
		
		System.out.println(info);
		String[] tabInfo = info.split("\t");
		int offset = 0;
		if (tabInfo.length == 2)offset=1;
		String[] info2 = tabInfo[1-offset].split("\\(");
		this.term = info2[0].trim();
		this.GO_ID = info2[1].substring(0,info2[1].indexOf(")"));
		this.Description= tabInfo[2-offset];
		return true;
	}
	
	public void printInfo(ExtendedWriter EW){
		EW.println(this.term+"("+this.GO_ID+")"+"\t"+this.Description);
	}
	
	public void printInfoShort(ExtendedWriter EW){
		EW.print(this.term+"("+this.GO_ID+")");
	}
	
	
	
}
