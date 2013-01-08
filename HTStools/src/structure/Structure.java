package structure;



import java.io.Serializable;

import general.Functions;
import general.RNAfunctions;
import general.energy;

import general.ExtendedWriter;




public class Structure  implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private int location2;
	private float probability;
	private float weight;
	private int knotNr = 0;
	
	private int kindOfStructure;
	
	public static final int backbone = -1;
	public static final int loop = -2;
	public static final int basepair = 1;
	public static final int knot = 2;

	
	
	private Structure nextStructure;
	

	Structure(){
		kindOfStructure = backbone;
	}

	
	Structure(int location2, float probability, float weight){
		this.kindOfStructure = basepair;
		this.location2 = location2;
		this.probability = probability;		
		this.weight = weight;
		
		this.nextStructure = null;
	}
	
	
	Structure(int location2, float probability){
		this.kindOfStructure = basepair;
		this.location2 = location2;
		this.probability = probability;		
		this.weight = 1;
		
		this.nextStructure = null;
	}
	
	Structure(int location2, float probability, int kindOfStructure){
		this.kindOfStructure = kindOfStructure;
		this.location2 = location2;
		this.probability = probability;		
		this.weight = 1;
		
		this.nextStructure = null;
	}

	
	Structure(int location2, float probability, int kindOfStructure, int knotNr){
		this.kindOfStructure = kindOfStructure;
		this.location2 = location2;
		this.probability = probability;		
		this.weight = 1;
		this.knotNr = knotNr;
		
		this.nextStructure = null;
	}
	

	public Structure removeStructure(float cutoff){
		if(this.probability < cutoff){
			if(this.nextStructure == null)
				return null;
			else
				return this.nextStructure.removeStructure(cutoff);
		}
		if(this.nextStructure != null)
			this.nextStructure = this.nextStructure.removeStructure(cutoff);
		return this; 
	}
	
	public boolean isBasepairing(int location2, double cutoff){
		if(this.probability >= cutoff && this.location2 == location2) return true;
		if(this.nextStructure == null)return false;
		return this.nextStructure.isBasepairing(location2, cutoff);
	}
	
	
	public Structure copy(){
		Structure newStructure = new Structure();
		newStructure.setLocation2(this.getLocation2());
		newStructure.setProbability(this.getProbability());
		newStructure.setKindOfStructure(this.getKindOfStructure());
		newStructure.setWeight(this.getWeight());
		if(this.nextStructure != null)
			newStructure.setNextStructure(this.nextStructure.copy());
		return newStructure;
	}
	
	public void shift(int shift){
		this.location2 = this.location2+shift;
		if(this.nextStructure != null)
			this.nextStructure.shift(shift);
	}

	
	public void changePointer(int location){
		if(this.location2 > location)
			location2--;
		if(this.nextStructure != null)
			this.nextStructure.changePointer(location);
		
	}
	public void setBackbone(){
		kindOfStructure = backbone;
	}
	
	public boolean isBackbone(){
		if(kindOfStructure == backbone)
			return true;
		return false;
	}

	public boolean isUnstructured(){
		if(kindOfStructure < 0)
			return true;
		return false;
	}
	
	
	public void setLoop(){
		kindOfStructure = loop;
	}
	
	public char getDotBracketStructure(int thisLocation){
		if(kindOfStructure == backbone)
			return '.';
		else if(kindOfStructure == loop)
			return ':';
		else if(kindOfStructure == basepair){
			if(this.getLocation2() < thisLocation)
				return ')';
			else
				return '(';
		}
		else if(kindOfStructure == knot){
			if(knotNr == 0){
				if(this.getLocation2() < thisLocation)
					return ']';
				else
					return '[';
			}
			else{
				if(this.getLocation2() < thisLocation)
					return '}';
				else
					return '{';
			}
		}
		System.out.println("Something is wrong structure is not preset");
		return '!';
	}

	
	public void addStructure(int location2, float probability){
		if(nextStructure == null)
			nextStructure = new Structure(location2,probability, 1);
		else
			nextStructure.addStructure(location2,probability, 1);
	}
	public void addStructure(Structure S){
		if(nextStructure == null)
			nextStructure = S;
		else
			nextStructure.addStructure(S);
	}
	
	
	public void addStructure(int location2, float probability,float weight){
		if(nextStructure == null)
			nextStructure = new Structure( location2,probability, weight);
		else
			nextStructure.addStructure(location2,probability,weight);
	}
	
	
	public void printStructure(int location){
		System.out.println(location+"\t"+(location+location2)+"\t"+probability);
		if(nextStructure != null)
			nextStructure.printStructure(location);
	}
	
	public void printIntraStructure(ExtendedWriter EW, int location1){
		if(location1 <= location2)
			EW.println(location1+"\t"+location2+"\t"+probability);
		else if(location2 < 0)
			EW.println(location1+"\t"+location2+"\t"+probability);
		if(nextStructure != null)
			nextStructure.printIntraStructure(EW, location1);
	}

	public void printIntraStructure(ExtendedWriter EW, int location1, char[] nucleotides){
		if(location1 <= location2){
			if(location2 < nucleotides.length){
				EW.print(" "+(location2+1)+":"+probability);
				if(!RNAfunctions.isBasepair(nucleotides[location1], nucleotides[location2]))
					System.out.println("Something wrong with basepair  :"+nucleotides[location1]+"-"+ nucleotides[location2]);
			}
			else
				EW.print(" "+(location2+1)+":"+probability);
		}
		else if(location2 < 0)
			EW.print(" "+(location2+1)+":"+probability);
		if(nextStructure != null)
			nextStructure.printIntraStructure(EW, location1,nucleotides);
	}

	
	public void printIntraStructure(int location1){
		boolean printed= false;
		if(location1 < location2){
			System.out.print(location2+":"+probability);
			printed = true;
		}
		else if(location2 < 0){
			System.out.print(location2+":"+probability);
			printed = true;
		}
		if(nextStructure != null){
			if(printed)
				System.out.print("   ");
			nextStructure.printIntraStructure(location1);
		}
	}

	public void printIntraStructurePrepencity(int location){
		float prob = getStructureProbabilities(location);
		if(prob > 0)
			System.out.print(prob);;
	}

	public float getStructureProbabilities(int location){
		float prob = 0;
		if(location != location2)
			prob = probability;
		if(nextStructure == null)
			return prob;
		else 
			return prob + nextStructure.getStructureProbabilities(location);
	
	}
	
	public boolean isOpen(int location, double cutoff){
		if(location != location2 && cutoff < probability)
			return false;
		if(nextStructure == null)
			return true;
		else 
			return nextStructure.isOpen(location, cutoff);
	}
	

/*	public float getStructureProbabilities(int location){
		if (location2 == location)
			return probability;
		
		if(nextStructure == null){
			System.out.println("There was no structure to "+location);
			return 0;
		}
		else
			return nextStructure.getStructureProbabilities(location);
	}
*/
	public float getStructureProbabilities(int startLocation, int thisLocation ,float prob){
		if(this.isUnstructured())
			return 0;
		
		if (this.location2 < thisLocation  || this.location2 > startLocation)
			prob =  this.probability*this.weight*(-1) + prob;
//		if(startLocation != thisLocation)
//		System.out.println(startLocation+"\t"+location2+"\t"+thisLocation);
		if(nextStructure == null){
			return prob;
		}
		else
			return nextStructure.getStructureProbabilities(startLocation,thisLocation,prob);
	}
	
	public int getOtherLocation(){
		return location2;
		
	}
	
	public int[] getOtherLocations(int[] otherLocations){
		otherLocations = Functions.addInt(otherLocations, this.location2);
		if(this.nextStructure != null)
			return this.nextStructure.getOtherLocations(otherLocations);
		return otherLocations;
	}
	
	public double[] getIAS(double[] IAS){
		IAS = Functions.addDouble(IAS, this.probability*this.weight*(-1));
		if(this.nextStructure != null)
			return this.nextStructure.getIAS(IAS);
		return IAS;
	}
	
	
	public Structure extractStructures(int subtraction){
		Structure newStructure = new Structure(location2 - subtraction, this.probability,this.weight);
		if(nextStructure != null)
			newStructure.nextStructure = nextStructure.extractStructures(subtraction);
		return newStructure;
	}
	
	
	public int[] getInteractions(int[] locations, int thisLocation){
		if(nextStructure != null)
			locations = nextStructure.getInteractions(locations,thisLocation);
		if(kindOfStructure > 0 ){
				locations = Functions.addInt(locations,location2);
		}else
			locations =  null;
		return locations;	
	}
	
	public Structure removeStructure(int location){
		if(location == location2){
			if(nextStructure != null)
				return nextStructure;
			else
				return null;
		}
		else{
			if(nextStructure != null)
				nextStructure = nextStructure.removeStructure(location);
			else
				System.out.println("There was no structure to "+location);
		}
		return this;
	}
	
	public void changeDir(int length){
		if(kindOfStructure > 0){
			location2 = length-1-location2;
			if(nextStructure != null)
				nextStructure.changeDir(length);
		}
	}
	
	public String printLocations(String locations){
		locations ="\t"+Functions.Int2String(location2);
		if(nextStructure != null)
			locations += nextStructure.printLocations(locations);
		return locations;
	}
	
	
	
	
	/**
	 * Returns the value of location2.
	 */
	public int getLocation2()
	{
		return location2;
	}

	/**
	 * Sets the value of location2.
	 * @param location2 The value to assign location2.
	 */
	public void setLocation2(int location2)
	{
		this.location2 = location2;
	}

	/**
	 * Returns the value of probability.
	 */
	public float getProbability()
	{
		return probability;
	}

	/**
	 * Sets the value of probability.
	 * @param probability The value to assign probability.
	 */
	public void setProbability(float probability)
	{
		this.probability = probability;
	}

	/**
	 * Returns the value of nextStructure.
	 */
	public Structure getNextStructure()
	{
		return nextStructure;
	}

	
	/**
	 * Returns if it has a nextStructure.
	 */
	public boolean hasNextStructure()
	{
		if(nextStructure != null)
			return true; 
		return false;
	}
	
	
	/**
	 * Sets the value of nextStructure.
	 * @param nextStructure The value to assign nextStructure.
	 */
	public void setNextStructure(Structure nextStructure)
	{
		this.nextStructure = nextStructure;
	}


	/**
	 * @return the kindOfStructure
	 */
	public int getKindOfStructure() {
		return kindOfStructure;
	}


	/**
	 * @param kindOfStructure the kindOfStructure to set
	 */
	public void setKindOfStructure(int kindOfStructure) {
		this.kindOfStructure = kindOfStructure;
	}
	
	
	public void setNewLocation(int shift){
		if(kindOfStructure > 0)
			this.location2 = this.location2-shift;
		if(nextStructure != null)
			nextStructure.setNewLocation(shift);
	}
	
	public void weightStructure(int thisLocation, energy A){
		if(this.location2 < thisLocation){
			this.weight = (float)A.intraEnergy(location2, thisLocation);
		}
		else{
			this.weight = (float)A.intraEnergy(thisLocation, location2);
		}
		
		if(this.nextStructure != null)
			nextStructure.weightStructure(thisLocation, A);
	}


	/**
	 * @return the weight
	 */
	public float getWeight() {
		return weight;
	}


	/**
	 * @param weight the weight to set
	 */
	public void setWeight(float weight) {
		this.weight = weight;
	}



}
