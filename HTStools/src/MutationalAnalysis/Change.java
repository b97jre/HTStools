package MutationalAnalysis;

import general.RNAfunctions;

public class Change {
	public String id; 
	public double baseMean;
	public double baseMeanA;
	public double baseMeanB;
	public double foldChange;
	public double log2FoldChange;
	public double pval;
	public double padj;
	
	
	public Change(){}
	
	public void addInfo(String info){
		String[] tabInfo = info.split("\t");
		this.id = tabInfo[1];
		this.baseMean = Double.parseDouble(tabInfo[2]);
		this.baseMeanA = Double.parseDouble(tabInfo[3]);
		this.baseMeanB = Double.parseDouble(tabInfo[4]);
		if(tabInfo[5].compareTo("Inf")== 0 || tabInfo[6].compareTo("Inf")== 0){
			if(this.baseMean > 10){
				this.foldChange = 4;
				this.log2FoldChange = 4;
			}else{
				this.foldChange = 0;
				this.log2FoldChange = 0;
			}
				
				this.pval = Double.parseDouble(tabInfo[7]);
			this.padj = Double.parseDouble(tabInfo[8]);
		}
		else if(tabInfo[5].compareTo("-Inf")== 0  || tabInfo[6].compareTo("-Inf")== 0 ){
			if(this.baseMean > 10){
				this.foldChange =-4;
				this.log2FoldChange = -4;
			}else{
				this.foldChange = 0;
				this.log2FoldChange = 0;
			}
			this.pval = Double.parseDouble(tabInfo[7]);
			this.padj = Double.parseDouble(tabInfo[8]);
		}
		else{
			if(tabInfo[6].compareTo("NA") != 0 &&tabInfo[5].compareTo("NA") != 0  ){
				this.foldChange = Double.parseDouble(tabInfo[5]);
				this.log2FoldChange = Double.parseDouble(tabInfo[6]);
				this.pval = Double.parseDouble(tabInfo[7]);
				this.padj = Double.parseDouble(tabInfo[8]);
			}
			else{
				this.foldChange = 0;
				this.log2FoldChange = 0;
				this.pval = 1;
				this.padj = 1;
			}
		}
	}

	
	
	public void printInfo(int offset){
		String[] location = this.id.split("\"");
		if(location.length > 1){
		String[] subLocation = location[1].split("\\ ");
		if(subLocation.length ==2){
			int loc = Integer.parseInt(subLocation[0]);
			int[] seq = RNAfunctions.RNAString2Int(subLocation[1]);
			String DNA = RNAfunctions.DNAInt2String(seq[0]);
			this.id = (loc-114)+"_"+DNA;
		}
		else if(subLocation.length ==4){
			int loc = Integer.parseInt(subLocation[0]);
			int[] seq = RNAfunctions.RNAString2Int(subLocation[1]);
			int loc2 = Integer.parseInt(subLocation[2]);
			int[] seq2 = RNAfunctions.RNAString2Int(subLocation[3]);
			String DNA = RNAfunctions.DNAInt2String(seq[0]);
			String DNA2 = RNAfunctions.DNAInt2String(seq2[0]);
			this.id = (loc-114)+"_"+DNA+"_"+(loc2-114)+"_"+DNA2;
		}
		}
		
		System.out.println(this.id+"\t"+this.baseMean+"\t"+this.baseMeanA+"\t"+this.baseMeanB+"\t"+this.foldChange+"\t"+this.log2FoldChange+"\t"+this.pval+"\t"+this.padj);
	}
}
