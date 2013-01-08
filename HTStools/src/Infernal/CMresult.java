package Infernal;

import java.util.ArrayList;

import general.ExtendedReader;
import general.ExtendedWriter;

public class CMresult {
	protected String fastaFile;
	protected int targetStart; 
	protected int targetStop;
	protected int queryStart;
	protected int queryStop;
	protected double bitScore;
	ArrayList <Motif> motifs;
	
	public CMresult(String fastaFile, int targetStart, int targetStop, int queryStart, int queryStop, double bitScore){
		this.fastaFile = fastaFile;
		this.targetStart = targetStart;
		this.targetStop = targetStop;
		this.queryStart = queryStart;
		this.queryStop = queryStop;
		this.bitScore = bitScore;
		this.motifs = new ArrayList <Motif>();
	}
	
	public void addDremeMotifs(Motif newMotif){
		this.motifs.add(newMotif);
		
	}
	
	public boolean addMotif(String info){
		String tab[] = info.split("\t");
		String info2[] = tab[0].split(";");
		if(info2[0].indexOf(this.fastaFile)  == 0 && this.targetStart == Integer.parseInt(info2[1])){
			boolean found = false;
			int count = 0;
			while(!found && count < motifs.size()){
				found = motifs.get(count).addNewMotif(info);
				count++;
			}
			if(!found)motifs.add(new Motif(info));
			return true;
		}
		else
			return false;
	}
	
	public boolean findGene(String info){
		return true;
	}


	public void printMotifs() {
		System.out.print(this.fastaFile+"\t"+this.targetStart+"\t"+this.bitScore+"\t");
		if(this.motifs.size() > 0){
			for(int i =0 ; i < this.motifs.size(); i++){
				this.motifs.get(i).printMotif();
				System.out.print("\t");
			}
			
		}
		System.out.println();
		// TODO Auto-generated method stub
		
	}
	
	public void printSpecificMotifs(String referenceFile, int length, int distance, String Sequence) {
		if(this.motifs.size() > 0 && this.motifs.size() < 3){
			System.out.print(referenceFile +"\t"+this.fastaFile+"\t"+this.targetStart+"\t"+this.targetStop+"\t"+this.bitScore+"\t");
			if(this.motifs.size() == 1){
				Motif TSS = new Motif();
				TSS.motif = "TSS";
				TSS.addTSS(length);
				this.motifs.add(TSS);
			}
			if(this.motifs.size() == 2){
				this.motifs.get(0).printSpecificMotif(100,60,4,this.motifs.get(1),this.bitScore,Sequence);
			}
			System.out.println();
		}
		else printMotifs();
		// TODO Auto-generated method stub
		
	}
	


}
