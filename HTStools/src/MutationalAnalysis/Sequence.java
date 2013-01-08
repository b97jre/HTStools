package MutationalAnalysis;

import java.util.ArrayList;

public class Sequence extends ArrayList<Mutation> {
		
	Sequence(int[] RefSeq, int[] seq){
		findMutations(RefSeq,seq);
	}
	
	public void findMutations(int[] RefSeq, int[] seq){
		for(int i = 0; i < RefSeq.length; i++){
			if(RefSeq[i] != seq[i])
				this.add(new Mutation(i,seq[i]));
		}
	}
	
	public boolean hasMutations(){
		if(this.size() > 0 ) return true;
		return false;
	}
	
	
	public int getNrOfMutations(){return this.size();}
	

}
