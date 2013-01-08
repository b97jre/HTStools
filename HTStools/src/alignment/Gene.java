package alignment;
import general.Functions;
import general.ExtendedReader;
import general.ExtendedWriter;


import java.io.FileWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;

import Blast.BlastHit;
import Ontologies.GOGene;
import Ontologies.GeneOntology;
import Ontologies.PantherGene;
import Sequence.FastaSequence;


public class Gene extends Object implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	protected String ID;
	protected String Name;
	protected int left;
	protected int right;
	protected String description;
	protected boolean plusStrand;
	protected int kind;

	protected final int mRNA = 1;
	protected final int ncRNA = 2;
	protected final int intergenic = 3;
	protected final int antisense = 4;
	protected final int repeat = 5;

	public FastaSequence fastaSeq;
	public int length;


	protected ArrayList <Double> fpkm_values;

	protected ArrayList <Hit> hits;
	protected ArrayList <BlastHit> blastHits;
	public GOGene GOgene;
	public PantherGene pantherGene;

	Gene(){
	}

	
	
	public Gene(String Name){
		this.Name = Name;
	}

	
	
	Gene(int left, int right,boolean plusStrand,String ID,String Name,String description){
		this.ID = ID;
		this.Name = Name;
		this.plusStrand = plusStrand;
		this.left = left;
		this.right = right;
		this.description = description;
	}

	
	public void addFPKM(double newFPKM){
		if(this.fpkm_values == null)
			this.fpkm_values = new ArrayList<Double>();
		this.fpkm_values.add(newFPKM);
	}
	
	public void setInfo(Gene otherGene){
		this.ID = otherGene.ID;
		this.Name = otherGene.Name;
		this.plusStrand = otherGene.plusStrand;
		this.left = otherGene.left;
		this.right = otherGene.right;
		this.description = otherGene.description;
	}

	
	Gene(int left, int right,String ID,String Name,String description){
		this.ID = ID;
		this.Name = Name;
		this.left = left;
		this.right = right;
		this.description = description;
	}


	public void addHit(Hit newHit){
		if(hits == null)
			hits = new ArrayList<Hit>();
		hits.add(newHit);
	}

	public void addBlastHit(BlastHit newHit){
		if(blastHits == null)
			blastHits = new ArrayList<BlastHit>();
		blastHits.add(newHit);
	}

	public int getNrOfBlastHits(){
		return this.blastHits.size();
	}

	public void removeWeakHits(double cutoff){
		double score = this.getHighestScore();
		for(int i = 0; i < blastHits.size(); i++){
			if(blastHits.get(i).getScore() / score < cutoff){
				blastHits.remove(i);
				i--;
			}
		}
	}

	public void mergeBlastHits(int distance, int penalty){
		System.out.println(this.Name);

		if(this.blastHits != null){
			for(int i = 0; i < blastHits.size(); i++){
				int before = -1;
				int beforeDist = distance+1;
				int afterDist = distance+1;
				int after = -1;
				for(int j = i+1; j <  blastHits.size(); j++){
					if(blastHits.get(i).getHitName().compareTo(blastHits.get(j).getHitName()) == 0){
						if(blastHits.get(i).isForward() && blastHits.get(j).isForward()){
							if(blastHits.get(i).getHitStart()-blastHits.get(j).getHitStop() < beforeDist && blastHits.get(i).getHitStart()-blastHits.get(j).getHitStop() > 0){
								beforeDist = blastHits.get(i).getHitStart()-blastHits.get(j).getHitStop();
								before = j;
							}
							if(blastHits.get(j).getHitStart()-blastHits.get(i).getHitStop() < beforeDist && blastHits.get(j).getHitStart()-blastHits.get(i).getHitStop() > 0){
								afterDist = blastHits.get(j).getHitStart()-blastHits.get(i).getHitStop();
								after = j;
							}
						}
						else if(!blastHits.get(i).isForward() && !blastHits.get(j).isForward()){
							if(blastHits.get(i).getHitStop()-blastHits.get(j).getHitStart() < beforeDist && blastHits.get(i).getHitStop()-blastHits.get(j).getHitStart() > 0){
								afterDist = blastHits.get(i).getHitStop()-blastHits.get(j).getHitStart();
								after = j;
							}
							if(blastHits.get(j).getHitStop()-blastHits.get(i).getHitStart() < beforeDist && blastHits.get(j).getHitStop()-blastHits.get(i).getHitStart() > 0){
								beforeDist = blastHits.get(j).getHitStop()-blastHits.get(i).getHitStart();
								before = j;
							}

						}
					}
				}
				if(before > -1 && blastHits.get(i).getScore()+ blastHits.get(before).getScore()-beforeDist*penalty > Math.max(blastHits.get(i).getScore(), blastHits.get(before).getScore()))
					blastHits.get(i).merge(blastHits.get(before),2,beforeDist);
				else
					before = -1;
				if(after > -1 && blastHits.get(i).getScore()+ blastHits.get(after).getScore()-afterDist*penalty > Math.max(blastHits.get(i).getScore(), blastHits.get(after).getScore()))
					blastHits.get(i).merge(blastHits.get(after),2, afterDist);
				else 
					after = -1;
				if(after != -1) {
					blastHits.remove(after);
					if(after < before) before--;
				}
				if(before != -1) blastHits.remove(before);
				if(before != -1 || after != -1) i--;

			}
		}
	}



	public boolean isAboveCutoff(double otherScore, double cutoff){
		double score = blastHits.get(0).getScore();
		if(otherScore / score < cutoff)return false;
		return true;
	}

	public double getHighestScore(){
		double highestScore = 0;
		for(int i = 0; i < blastHits.size(); i++){
			if(blastHits.get(i).getScore() > highestScore){
				highestScore = blastHits.get(i).getScore();
			}
		}
		return highestScore;
	}

	
	
	public String getBestHit(){
		double highestScore = 0;
		int pointer = -1;
		if(this.blastHits != null){
			for(int i = 0; i < blastHits.size(); i++){
				if(blastHits.get(i).getScore() > highestScore){
					highestScore = blastHits.get(i).getScore();
					pointer = i;
				}
			}
			return blastHits.get(pointer).getHitName();
		}
		return null;
	}

	
	public int getBestHitPointer(){
		double highestScore = 0;
		int pointer = -1;
		if(this.blastHits != null){
			for(int i = 0; i < blastHits.size(); i++){
				if(blastHits.get(i).getScore() > highestScore){
					highestScore = blastHits.get(i).getScore();
					pointer = i;
				}
			}
			return pointer;
		}
		return -1;
	}
	

	public int getLongestAssembly(){
		double highestScore = 0;
		int longestAssembly = 0;
		for(int i = 0; i < blastHits.size(); i++){
			if(blastHits.get(i).getScore() > highestScore){
				highestScore = blastHits.get(i).getScore();
				longestAssembly = blastHits.get(i).getLength();
			}
		}
		return longestAssembly;
	}



	public String removeSimilairHits(double cutoff){
		double score = getHighestScore();
		String removedGenes = "";
		for(int i = 0; i < blastHits.size(); i++){
			if(blastHits.get(i).getScore() / score > cutoff && blastHits.get(i).getHitName().compareTo(this.Name)!= 0){

				Gene A = blastHits.get(i).getSpecificGene();
				if(A != null){
					double otherHighestScore = A.getHighestScore();
					if(otherHighestScore <= score  ){
						String otherName = blastHits.get(i).getHitName();
						System.out.println(this.Name+"\t"+ score+"\t"+ otherName+"\t"+otherHighestScore );
						blastHits.remove(i);
						i--;
						removedGenes += otherName+"\t";
					}
				}else{
					blastHits.remove(i);
					i--;
				}
			}
		}
		removedGenes.trim();
		return removedGenes;
	}



	public void removeAllHits(){
		while(this.blastHits.size() > 1)
			this.blastHits.remove(0);
	}

	public void removeNotPresentHits(Hashtable<String,Gene> Genes){
		for(int i = 0; i < blastHits.size(); i++){
			if(!Genes.containsKey(blastHits.get(i).getHitName())){
				blastHits.remove(i);
				i--;
			}
		}
	}

	public void getBlastSize(){
		for(int i = 0; i < blastHits.size(); i++){
			if(blastHits.get(i).getHitName().compareTo(blastHits.get(i).getQueryName())== 0){
				this.length = blastHits.get(i).getQueryStop() -  blastHits.get(i).getQueryStart() +1;
				i = blastHits.size();
			}
		}
	}

	public void checkDependencies(Hashtable<String,Gene> Genes){
		for(int i = 0; i < blastHits.size(); i++){
			if(Genes.containsKey(blastHits.get(i).getHitName())){
				Gene temp = Genes.get(blastHits.get(i).getHitName());
				if(!temp.containsQuery(blastHits.get(i))){
					//					System.out.println("1\t"+Genes.get(blastHits.get(i).getHitName()).blastHits.size());
					Genes.get(blastHits.get(i).getHitName()).blastHits.add(BlastHit.reversed(blastHits.get(i)));
					//					System.out.println("2\t"+Genes.get(blastHits.get(i).getHitName()).blastHits.size());
				}
			}
		}
	}



	public ArrayList <String> clusterHits(ExtendedWriter EW, Hashtable<String,Gene> Genes, int number){
		for(int i = 0; i < blastHits.size(); i++){
			if(Genes.containsKey(blastHits.get(i).getHitName())){
				Gene temp = Genes.get(blastHits.get(i).getHitName());
				for(int j = 0; i < temp.blastHits.size(); i++){
					if(!this.contains(temp.blastHits.get(i))){
						this.blastHits.add(temp.blastHits.get(i));
					}
				}
			}
		}
		ArrayList <String> genes = new ArrayList <String>();
		int length = 0;
		int count = 0;
		String seq =  "";
		for(int i = 0; i < blastHits.size(); i++){
			if(Genes.containsKey(blastHits.get(i).getHitName())){
				length += Genes.get(blastHits.get(i).getHitName()).length;
				count++;
				genes.add(blastHits.get(i).getHitName());
			}

		}
		length = length / count;
		EW.print("cluster_"+number+"\t"+count+"\t"+length+"\t");
		for(int i = 0; i < genes.size()-1;i++){
			EW.print(genes.get(i)+",");
		}
		EW.println(genes.get(genes.size()-1));
		return genes;
	}


	public ArrayList <BlastHit> printcluster(Hashtable<String,Gene> Genes){
		int count = 0;
		int length = 0;
		for(int i = 0; i < blastHits.size(); i++){
			if(Genes.containsKey(blastHits.get(i).getHitName())){
				Gene temp = Genes.get(blastHits.get(i).getHitName());

			}
		}
		for(int i = 0; i < blastHits.size(); i++){
			if(Genes.containsKey(blastHits.get(i).getHitName())){
				Genes.get(blastHits.get(i).getHitName()).blastHits = this.blastHits;
			}
		}
		return blastHits;
	}


	public boolean contains(BlastHit BH){
		for(int i = 0; i < blastHits.size(); i++){
			if(blastHits.get(i).getHitName().compareTo(BH.getHitName()) == 0)
				return true;
		}
		return false;
	}

	public boolean containsQuery(BlastHit BH){
		for(int i = 0; i < blastHits.size(); i++){
			if(blastHits.get(i).getHitName().compareTo(BH.getQueryName()) == 0)
				return true;
		}
		return false;
	}

	public void printBlastHits(ExtendedWriter EW){

		if(this.blastHits!= null){
			for(int i = 0; i < blastHits.size(); i++){
				blastHits.get(i).printHit(EW);
			}
		}
	}


	public void printAnnotation(ExtendedWriter EW,GeneOntology A){
		String Name = this.Name;
		Name = Name.replace("trinity", "Acr");
		Name = Name.substring(0,Name.length()-2);
		EW.print(Name+"\t"+this.length*3+"\t");

		if(this.pantherGene!= null){
			this.pantherGene.printPantherClassInfo(EW);
		}
		else{
			EW.print("\t");			
		}
		EW.print("\t");	
		if(this.GOgene!= null){
			EW.print(this.GOgene.geneName);
		}
		EW.print("\t");	
		if(this.pantherGene == null || !this.pantherGene.printPantherGOInfo(EW) ){
			if(this.GOgene!= null){
				A.printGOinfo(this.GOgene, EW);
			}
			else
				EW.print("\t\t");
				
		}
		EW.print("\t");	
		if(this.pantherGene == null || !this.pantherGene.printOtherInfo(EW)){
			EW.print("\t");
		}
		EW.println();
	}



	public void printBestHit(ExtendedWriter EW, Hashtable HT ){
		if(this.blastHits!= null){
			double score = getHighestScore();
			int pointer = 0;
			boolean found = false;
			BlastHit temp = this.blastHits.get(0);
			while(!found && blastHits.size() != 0){
				double nextBestScore = 0;
				if(this.blastHits!= null){
					for(int i = 0; i < blastHits.size(); i++){
						if(blastHits.get(i).getScore() > nextBestScore){
							nextBestScore = blastHits.get(i).getScore();
							pointer = i;
						}
					}
					if(!HT.containsKey(blastHits.get(pointer).getHitName())){
						blastHits.get(pointer).printBestHit(EW,score);
						found = true;
					}  
					else{
						blastHits.remove(pointer);
					}
				}
			}
			if(!found)
				temp.printNoHit(EW,score);

		}
	}

	
	public void printBestCoverage(ExtendedWriter EW){
		if(this.blastHits!= null){
			int pointer = getBestHitPointer();
			if(pointer != -1 ){
				BlastHit temp = this.blastHits.get(pointer);
				temp.printHit(EW);
//				int length = this.length;
//				int blastLength = temp.getQueryStop() - temp.getQueryStart();
//				double start = (double)temp.getQueryStart()/(double)length;
//				double stop = (double)temp.getQueryStop()/(double)length;
//				double fraction = (double)blastLength / (double)length;
			
//				EW.println(this.Name+"\t"+temp.getHitName()+"\t"+fraction+"\t"+start+"\t"+stop+"\t"+blastLength+"\t"+length+"\t"+temp.getHitStart()+"\t"+temp.getHitStop());
			}
		}
		else
			EW.println(this.Name+"\tNaN\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0");
	}
	
	
	public void printBestLength(ExtendedWriter EW, Hashtable HT ){
		if(this.blastHits!= null){
			double length = getLongestAssembly();
			int pointer = 0;
			boolean found = false;
			BlastHit temp = this.blastHits.get(0);
			while(!found && blastHits.size() != 0){
				double nextBestScore = 0;
				if(this.blastHits!= null){
					for(int i = 0; i < blastHits.size(); i++){
						if(blastHits.get(i).getScore() > nextBestScore){
							nextBestScore = blastHits.get(i).getScore();
							pointer = i;
						}
					}
					if(!HT.containsKey(blastHits.get(pointer).getHitName())){
						blastHits.get(pointer).printBestLength(EW,length);
						found = true;

					}  
					else{
						blastHits.remove(pointer);
					}
				}
			}
			if(!found)
				temp.printNoHit(EW,length);

		}
	}



	public int getNrOfHits(){
		if(hits != null)
			return hits.size();
		return 0;
	}

	public void printNrOfHits(ExtendedWriter EW, String chromosomeName, int chromosomeLength){
		if(hits != null){
			EW.println(Functions.fixedLength(this.Name, 20)+"\t"+chromosomeName+"\t+\t"+this.kind+"\t"+this.hits.size()+"\t"+this.left+"\t"+this.right+"\t"+(this.left-chromosomeLength)+"\t"+(this.right-chromosomeLength) +"\t"+chromosomeName+":"+this.left+".."+this.right);
		}
	}


	public int size(){
		return this.right- this.left;
	}

	public boolean isOverlapping(Gene otherGene){
		if(this.left <= otherGene.left && this.right >= otherGene.left)return true;
		if(this.left <= otherGene.right && this.right >= otherGene.right)return true;
		return false;
	}


	public void printNrOfDirectedHits(ExtendedWriter EW, String chromosomeName, int chromosomeLength){
		if(hits != null){
			int sense = 0; 
			int antisense = 0;
			double weightedSense = 0;
			double weightedAntisense = 0;
			for(int i =0 ; i < hits.size(); i++){
				if(this.plusStrand != hits.get(i).plusStrand){
					antisense++;
					weightedAntisense += hits.get(i).getWeightedHit();
				}
				else{
					sense++;
					weightedSense += hits.get(i).getWeightedHit();
				}
			}
			double w1 = sense;
			double w2 = antisense;
			if(weightedSense < 2) w1 = 2;
			if(weightedAntisense < 2) w2 = 2;

			double difference = java.lang.Math.abs(java.lang.Math.log10(w1/w2));
			if(isAboveCutoff(0.7,40) )
				EW.println(Functions.fixedLength(this.Name, 20)+"\t"+chromosomeName+"\t"+this.kind
						+"\t"+(sense+antisense)+"\t"+Functions.fixedLength((weightedSense+weightedAntisense), 5)
						+"\t"+sense+"\t"+Functions.fixedLength(weightedSense, 5)
						+"\t"+antisense+"\t"+Functions.fixedLength(weightedAntisense, 5)
						+"\t"+Functions.fixedLength(difference, 5)+"\t"+this.left+"\t"+this.right+"\t"
						+"\t"+chromosomeName+":"+this.left+".."+this.right);
		}
	}

	public boolean isAboveCutoff(double diffCutOff, int countCutoff ){
		this.sortHits();
		//		int maxSense = this.getHighestPlusStrandHit();
		//		int maxAntisense = this.getHighestMinusStrandHit();
		//		if(maxSense == 0) maxSense = 1;
		//		if(maxAntisense == 0) maxAntisense = 1;
		//		double difference = java.lang.Math.abs(java.lang.Math.log10(maxSense/maxAntisense));
		//		if(difference >diffCutOff  || ( maxSense > countCutoff || maxAntisense > countCutoff))
		return true;
		//		return false;
	}




	public boolean printpremiRNAstructures(String Dir, String file){
		if(hits != null && isAboveCutoff(1.0,5)){

			try{
				ExtendedWriter EW = null;
				if(this.plusStrand == true)
					EW = new ExtendedWriter(new FileWriter(Dir+"/"+this.Name+"_+."+file+".hits"));
				else
					EW = new ExtendedWriter(new FileWriter(Dir+"/"+this.Name+"_-."+file+".hits"));
				for(int i = 0; i < this.hits.size();i++){
					this.hits.get(i).printHit(EW);
				}
				EW.flush();
				EW.close();
				removeDuplicates();
				if(this.hits.size() < 30)
					return true;
			}
			catch(Exception E){
				E.printStackTrace(); 
			}
		}

		return false;
	}

	public void printSurrounding(Chromosome C, ExtendedWriter EW, int surrounding){
		for(int i = 0; i < this.hits.size();i++){
			//double difference = java.lang.Math.abs(java.lang.Math.log10((double)this.hits.get(i).nrOfReads/(double)min));
			//if(difference > 0.9 )
			this.hits.get(i).printSurrounding(C, 5, surrounding, this.Name, EW);
		}
	}

	public void print3UTR(Chromosome C, ExtendedWriter EW, int length){
		if(this.plusStrand)
			C.print3UTR(this.right, length, this.plusStrand, this.Name, EW);
		else
			C.print3UTR(this.left, length, this.plusStrand, this.Name, EW);

	}


	public void print3UTRtab(Chromosome C, ExtendedWriter EW, int length){
		if(this.plusStrand)
			C.print3UTRtab(this.right, length, this.plusStrand, this.Name, EW);
		else
			C.print3UTRtab(this.left, length, this.plusStrand, this.Name, EW);

	}



	private void sortHits(){
		for(int i = 1; i < this.hits.size();i++){
			int location = findHit(i, this.hits.get(i));
			if(location != -1){
				this.hits.get(location).nrOfReads++;
			}
		}
		this.hits.trimToSize();
	}

	private void removeDuplicates(){
		for(int i = 1; i < this.hits.size();i++){
			int location = findHit(i, this.hits.get(i));
			if(location != -1){
				this.hits.remove(i);
				i--;
			}
		}
		this.hits.trimToSize();
	}




	public int findHit(int stop, Hit newHit){
		for(int i = 0; i < stop;i++){
			if(this.hits.get(i).hasSameLocation(newHit))
				return i;
		}
		return -1;
	}

	public int getHighestPlusStrandHit(){
		int max = 0;
		for(int i = 0; i < this.hits.size();i++){
			if(this.hits.get(i).plusStrand && this.hits.get(i).nrOfReads > max)
				max = this.hits.get(i).nrOfReads;
		}
		return max;
	}

	public int getHighestMinusStrandHit(){
		int max = 0;
		for(int i = 0; i < this.hits.size();i++){
			if(!this.hits.get(i).plusStrand && this.hits.get(i).nrOfReads > max)
				max = this.hits.get(i).nrOfReads;
		}
		return max;
	}

	public int getPlusStrandHits(){
		int sum = 0;
		for(int i = 0; i < this.hits.size();i++){
			if(this.hits.get(i).plusStrand )
				sum += this.hits.get(i).nrOfReads;
		}
		return sum;
	}

	public int getMinusStrandHits(){
		int sum = 0;
		for(int i = 0; i < this.hits.size();i++){
			if(!this.hits.get(i).plusStrand )
				sum += this.hits.get(i).nrOfReads;
		}
		return sum;
	}





	public void printHits(String Dir, String file){
		if(hits != null){
			int sense = 0; 
			int antisense = 0;
			double weightedSense = 0;
			double weightedAntisense = 0;
			for(int i =0 ; i < hits.size(); i++){
				if(this.plusStrand != hits.get(i).plusStrand){
					antisense++;
					weightedAntisense += hits.get(i).getWeightedHit();
				}
				else{
					sense++;
					weightedSense += hits.get(i).getWeightedHit();
				}
			}
			double w1 = weightedSense;
			double w2 = weightedAntisense;
			if(weightedSense < 2) w1 = 2;
			if(weightedAntisense < 2) w2 = 2;
			double difference = java.lang.Math.abs(java.lang.Math.log10(w1/w2));
			if(difference > 0.5  && ( weightedSense > 50 || weightedAntisense > 50)){
				try{
					ExtendedWriter EW = new ExtendedWriter(new FileWriter(Dir+"/"+this.Name+"."+file+".hits"));
					for(int i = 0; i < this.hits.size();i++){
						this.hits.get(i).printHit(EW);
					}
					EW.flush();
					EW.close();

				}
				catch(Exception E){
					E.printStackTrace(); 
				}
			}
		}
	}





	public void compareNROfHits(String chromosomeName, double otherNrOfHits, ExtendedWriter ER){
		double nrOfHits = this.getNrOfHits();

		double difference = 0;
		if(nrOfHits > 10 || otherNrOfHits > 10){
			if(nrOfHits < 10)
				nrOfHits = 10;
			else if(otherNrOfHits < 10)
				otherNrOfHits = 10;
			else
				difference = Math.log(otherNrOfHits/nrOfHits);
		}
		if(difference != 0)
			ER.println(this.Name+","+chromosomeName +","+difference+","+nrOfHits+","+otherNrOfHits+","+chromosomeName+":"+this.left+".."+this.right);



	}

	public int getKind(boolean plusStrand){
		if(kind == intergenic)
			return intergenic;

		else if(kind == ncRNA){
			if(this.plusStrand == plusStrand)
				return ncRNA;
			else
				return antisense;
		}
		else if(kind == mRNA){
			if(this.plusStrand == plusStrand)
				return mRNA;
			else
				return antisense;
		}
		System.out.println("Should not end up here in getKind(bool) in Gene");
		return -1;
	}



	public void print(){
		if(plusStrand){
			System.out.println(ID+"\t"+left+"\t"+right+"\t+\t"+description);
		}
		else
			System.out.println(ID+"\t"+left+"\t"+right+"\t-\t"+description);
	}

	public boolean isParent(mRNA newmRNA){
		if(newmRNA.parent.compareTo(ID) == 0){
			return true;
		}
		return false;
	}

	public String getName() {
		return Name;
	}

	public void setName(String name) {
		Name = name;
	}

	public String getDescription() {
		return description;
	}

	public void setDescription(String description) {
		this.description = description;
	}

	public FastaSequence getFastaSeq() {
		return fastaSeq;
	}

	public void setFastaSeq(FastaSequence fastaSeq) {
		this.fastaSeq = fastaSeq;
	}


}
