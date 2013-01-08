package Infernal;

import java.util.ArrayList;

public class Motif {

	String motif;
	ArrayList <MotifScore> hits;
	//DDB0232431;959644;959587;34.69
	//fimo	
	//polypeptide_motif	
	//90	
	//96	
	//1.84877
	//+
	//.
	//motif_name=CCCAYAA;pvalue=0.00216;qvalue=;sequence=CACATAA

	
	public Motif(){
		this.hits = new ArrayList<MotifScore>();
	}
	public Motif(String info){
		this.addMotif(info);
	}
	
	
	public void addMotif(String info){
		String tab[] = info.split("\t");
		MotifScore A = new MotifScore();
		A.start = Integer.parseInt(tab[3]);
		A.stop = Integer.parseInt(tab[4]);
		if(tab[6].compareTo("+") == 0)
			A.clockwise = true;
		else
			A.clockwise = false;
		
		A.bitScore = Double.parseDouble(tab[5]);
		
		String tab2[] = tab[8].split(";");
		for(int i = 0; i < tab2.length; i++){
			String[] inf = tab2[i].split("=");
			if(inf[0].compareTo("motif_name") == 0) motif = inf[1];
			if(inf[0].compareTo("sequence") == 0) A.sequence =inf[1];
			if(inf[0].compareTo("pvalue") == 0) A.pValue = Double.parseDouble(inf[1]);
		}
		this.hits.add(A);
	}
	
	public void addTSS(int length){
		MotifScore A = new MotifScore();
		A.start = length;
		A.stop = length;
		A.clockwise = true;
		A.bitScore = 0;
		A.sequence ="TSS";
		A.pValue = 1;
		this.hits.add(A);
	}

	public void addNegative(double score, int length){
		MotifScore A = new MotifScore();
		A.start = length;
		A.stop = length;
		A.clockwise = true;
		A.bitScore = -10;
		A.sequence ="NoHit";
		A.pValue = 1;
		this.hits.add(A);
	}
	
	
	public void switchMotif(String info){
		String tab[] = info.split("\t");
		MotifScore A = new MotifScore();
		A.start = Integer.parseInt(tab[3]);
		A.stop = Integer.parseInt(tab[4]);
		if(tab[6].compareTo("+") == 0)
			A.clockwise = true;
		else
			A.clockwise = false;
		
		A.bitScore = Double.parseDouble(tab[5]);
		
		String tab2[] = tab[8].split(";");
		for(int i = 0; i < tab2.length; i++){
			String[] inf = tab2[i].split("=");
			if(inf[0].compareTo("motif_name") == 0) motif = inf[1];
			if(inf[0].compareTo("sequence") == 0) A.sequence =inf[1];
			if(inf[0].compareTo("pvalue") == 0) A.pValue = Double.parseDouble(inf[1]);
		}
		if(this.hits.size() > 0)
			this.hits.removeAll(hits);
		this.hits.add(A);
	}
	
	
	public boolean changeMotif(String info){
		boolean change = false;
		boolean found = false;
		
		if(this.motif == null) change = true;
		
		String tab[] = info.split("\t");
		String tab2[] = tab[8].split(";");
		if(!change){
			for(int i = 0; i < tab2.length; i++){
				String[] inf = tab2[i].split("=");
				if(inf[0].compareTo("motif_name") == 0 && this.motif.compareTo(inf[1]) == 0){
					found = true;
					if(this.hits.size() == 0 )//|| this.hits.get(0).bitScore < Double.parseDouble(tab[5]))
						change = true;
				}
			}
		}
		if(change){
			switchMotif(info);
		}
		return found;
	}
	
	public boolean addNewMotif(String info){
		boolean add = false;
		boolean found = false;
		boolean change = false;
		
		if(this.motif == null) {
			add = true;
			found =true;
		}
		
		String tab[] = info.split("\t");
		String tab2[] = tab[8].split(";");
		if(!add){
			for(int i = 0; i < tab2.length; i++){
				String[] inf = tab2[i].split("=");
				if(inf[0].compareTo("motif_name") == 0 && this.motif.compareTo(inf[1]) == 0){
					found = true;
					if(this.hits.size() == 0) add = true;
					else if (this.hits.get(0).bitScore == Double.parseDouble(tab[5]))
						add = true;
					else if(this.hits.get(0).bitScore < Double.parseDouble(tab[5])) 
						change = true;
				}
			}
		}
		if(add){
			addMotif(info);
		}
		else if (change)
			switchMotif(info);
		return found;
	}
	
//	public void printDistance(Motif otherMotif){
//		if(this.clockwise != otherMotif.clockwise ){
//			System.out.print("*");
//			System.out.print(this.sequence+"\t"+this.start+"\t"+this.stop);
//			System.out.print(-100);
//			System.out.print(otherMotif.sequence+"\t"+otherMotif.start+"\t"+otherMotif.stop);
//			
//		}
//		else{
//			if(this.start < otherMotif.start){
//				if(this.clockwise)System.out.print("+\t");
//				else System.out.print("-\t");
//				System.out.print(this.sequence+"\t"+this.start+"\t"+this.stop);
//				System.out.print(otherMotif.start-this.stop);
//				System.out.print(otherMotif.sequence+"\t"+otherMotif.start+"\t"+otherMotif.stop);
//				System.out.print("\t"+this.pValue+"\t"+otherMotif.pValue);
//			}
//			else{
//				if(this.clockwise)System.out.print("+\t");
//				else System.out.print("-\t");
//				System.out.print(otherMotif.sequence+"\t"+otherMotif.start+"\t"+otherMotif.stop);
//				System.out.print(this.start-otherMotif.stop);
//				System.out.print(this.sequence+"\t"+this.start+"\t"+this.stop);
//				System.out.print("\t"+otherMotif.pValue+"\t"+this.pValue);
//			}
//		}
//	}


	public void printMotif() {
		System.out.print(this.motif+"\t");
		if(this.hits.size() > 0 ){
		for(int i = 0; i < this.hits.size()-1;i++)System.out.print(this.hits.get(i).start+";");
		System.out.print(this.hits.get(this.hits.size()-1).start+"\t");
		for(int i = 0; i < this.hits.size()-1;i++)System.out.print(this.hits.get(i).stop+";");
		System.out.print(this.hits.get(this.hits.size()-1).stop+"\t");
		for(int i = 0; i < this.hits.size()-1;i++)System.out.print(this.hits.get(i).clockwise+";");
		System.out.print(this.hits.get(this.hits.size()-1).clockwise+"\t");
		for(int i = 0; i < this.hits.size()-1;i++)System.out.print(this.hits.get(i).pValue+";");
		System.out.print(this.hits.get(this.hits.size()-1).pValue+"\t");
		for(int i = 0; i < this.hits.size()-1;i++)System.out.print(this.hits.get(i).sequence+";");
		System.out.print(this.hits.get(this.hits.size()-1).sequence+"\t");
		}
		else 
			System.out.print("-1\t-1\t-1\t1000\tNotFound");
	}

	public void printSpecificMotif(int length, int distance) {
		int pointer = 0;
		int closestDist = length;
		for(int i = 0; i < this.hits.size();i++){
			if(this.hits.get(i).clockwise){
				int dist = length - (this.hits.get(i).stop-(this.hits.get(i).stop - this.hits.get(i).start)/2);
				if(Math.abs(dist - distance )< closestDist) pointer = i;
			}
		}
		System.out.print(this.motif+"\t");
		if(this.hits.size() > 0 ){
		System.out.print(this.hits.get(pointer).start+"\t");
		System.out.print(this.hits.get(pointer).stop+"\t");
		System.out.print(this.hits.get(pointer).clockwise+"\t");
		System.out.print(this.hits.get(pointer).bitScore+"\t");
		System.out.print(this.hits.get(pointer).sequence+"\t");
		}
		else 
			System.out.print("-1\t-1\t-1\t0\tNotFound");
	}
	
	public void printSpecificMotif(int length, int distance, int variation,  Motif otherMotif, double CMbitScore,String Sequence) {
		this.addNegative(-10,length);
		otherMotif.addTSS(length);
		int pointer = 0;
		int otherPointer = 0;
		int closestDist = length;
		boolean Dist = false;;
		double bestScore = -16;
		int bestDist = length;
		for(int i = 0; i < this.hits.size();i++){
			double tempBitScore =  this.hits.get(i).bitScore;
			int loc = length - (this.hits.get(i).stop-(this.hits.get(i).stop - this.hits.get(i).start)/2);
			if(!this.hits.get(i).clockwise){
				tempBitScore = -10;
				loc = -1;
			}
			for(int j = 0; j< otherMotif.hits.size();j++){
				double otherBitScore =  otherMotif.hits.get(j).bitScore;
				int loc2 = length - (otherMotif.hits.get(j).stop-(otherMotif.hits.get(j).stop - otherMotif.hits.get(j).start)/2);
				if(!otherMotif.hits.get(j).clockwise){
					otherBitScore = -5;
					loc2 = -1;
				}
				int dist = loc -loc2;
				
				if(distance-variation <= dist && dist < distance+variation && tempBitScore+otherBitScore+5 > bestScore){
//					System.out.println();
//					System.out.println(dist+"\t"+loc+"\t"+loc2);
					pointer = i;
					otherPointer = j;
					bestScore = tempBitScore+otherBitScore+5;
					Dist = true;
					bestDist = dist;
				}
				else if(tempBitScore+otherBitScore-5 > bestScore){
					pointer = i;
					otherPointer = j;
					bestScore = tempBitScore+otherBitScore-5;
					Dist = false;
					bestDist = dist;
				}
			}
		}
		System.out.print(CMbitScore+bestScore+"\t");
		System.out.print(bestDist+"\t");
		System.out.print(this.motif+"\t");
		System.out.print(this.hits.get(pointer).start-length+"\t");
		System.out.print(this.hits.get(pointer).stop-length+"\t");
		System.out.print(this.hits.get(pointer).clockwise+"\t");
		System.out.print(this.hits.get(pointer).bitScore+"\t");
		System.out.print(this.hits.get(pointer).sequence+"\t");
		if(this.hits.get(pointer).clockwise)
			System.out.print(Sequence.substring(Math.max(this.hits.get(pointer).start-5,0),Math.min(this.hits.get(pointer).stop+5,length))+"\t");
		else
			System.out.print("NNNNN\t");
			
		System.out.print(otherMotif.motif+"\t");
		System.out.print(otherMotif.hits.get(otherPointer).start-length+"\t");
		System.out.print(otherMotif.hits.get(otherPointer).stop-length+"\t");
		System.out.print(otherMotif.hits.get(otherPointer).clockwise+"\t");
		System.out.print(otherMotif.hits.get(otherPointer).bitScore+"\t");
		System.out.print(otherMotif.hits.get(otherPointer).sequence+"\t");
		if(this.hits.get(pointer).clockwise)
			System.out.print(Sequence.substring(Math.max(otherMotif.hits.get(otherPointer).start-5,0),Math.min(otherMotif.hits.get(otherPointer).stop+5,length))+"\t");
		else
			System.out.print("NNNNN\t");
			
		
		if(Math.abs(closestDist)<6)System.out.print("1");
		else System.out.print("0");
			
	}
	
	
	public String getMotif() {
		return motif;
	}
	public void setMotif(String motif) {
		this.motif = motif;
	}
	public class MotifScore{
		public String sequence;
		public double bitScore;
		public int  start;
		public int stop;
		public boolean clockwise;
		public double pValue;
		
	}

}
