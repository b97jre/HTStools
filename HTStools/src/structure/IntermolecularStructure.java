package structure;

import java.io.Serializable;

import general.Functions;

import general.ExtendedWriter;




public class IntermolecularStructure  implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private int mRNApointer;
	private int sRNApointer;
	private int initiationLength;
	private float intermolecularScore;
	private float mRNAstructurePrepencity;
	private float sRNAstructurePrepencity;
	
	private float consensusScore;
	
	private  float consensusWeight;
	private float intramolecularWeight;
	private float intermoleucularWeight;

	private char[] mRNAsequence;
	private char[] sRNAsequence;

	private float previousScore;
	private float danglingScore;
	
	private boolean initiationStructure;
	private boolean stableStructure;


	private int mRNAstart;
	private int sRNAstart;
	

	public boolean isInitiationPoint(int mRNApointer, int sRNApointer){
		if(mRNApointer == this.getmRNApointer() && sRNApointer == this.getsRNApointer())
			return true;
		return false;
	}
	
	public boolean isInitiationRegion(){
		return this.initiationStructure;
	}
	
	public void printInteraction(ExtendedWriter A , int sRNApos,int mRNApos){
			if(mRNAsequence != null )
				A.println(Functions.fixedLength("  "+this.mRNAstructurePrepencity,6)+ Functions.fixedLength(4)
						+Functions.fixedLength(mRNApos,6)+" "+mRNAsequence[0]+"-"+sRNAsequence[0]
						+Functions.fixedLength(3)+Functions.fixedLength(sRNApos+1,6)+Functions.fixedLength(3)
						+ Functions.fixedLength(this.sRNAstructurePrepencity,6)+Functions.fixedLength(6)
						+Functions.fixedLength(getConsensusScore(),6)+Functions.fixedLength(6)
						+Functions.fixedLength(getIntermolecularScore(),6)+Functions.fixedLength(6)
						+Functions.fixedLength(getTotalScore(),6)
						
				);
			else 
				A.println(Functions.fixedLength(6)+ Functions.fixedLength(3)
						+Functions.fixedLength(6)+" x-x "						
				);
//			A.print(Functions.fixedLength(getTotalScore(),6)+"\t");
//			A.print(Functions.fixedLength(getIntermolecularScore(),6)+"\t");
//			A.print(Functions.fixedLength(this.sRNAstructurePrepencity,6)+"\t");				
//			A.print(Functions.fixedLength(this.mRNAstructurePrepencity,6)+"\t");				
//			A.println(Functions.fixedLength(getConsensusScore(),6)+"\t");
/*			A.print(Functions.fixedLength(getPreviousScore(),6)+"\t");
			A.print(Functions.fixedLength(getsRNApointer(),6)+"\t");
			A.print(Functions.fixedLength(getmRNApointer(),6)+"\t");
			A.print(Functions.fixedLength(getsRNAstart(),6)+"\t");
			A.print(Functions.fixedLength(getmRNAstart(),6)+"\t");
			A.println(Functions.fixedLength(getInitiationLength(),6)+"\t");
*/			if(mRNAsequence != null){
				for(int i = 1; i < Functions.max(mRNAsequence.length, sRNAsequence.length);i++){
					A.print(Functions.fixedLength(16));
					if(mRNAsequence.length> i)
						A.print(mRNAsequence[i]);
					else
						A.print("|");
					A.print("   ");
					if(sRNAsequence.length> i)
						A.println(sRNAsequence[i]);
					else
						A.println("|");
				}
			}
	}
	
	
	
	public String printInteraction(String A , int sRNApos,int mRNApos){
		if(mRNAsequence != null )
			A+=Functions.fixedLength("  "+this.mRNAstructurePrepencity,6)+ Functions.fixedLength(4)
					+Functions.fixedLength(mRNApos,6)+" "+mRNAsequence[0]+"-"+sRNAsequence[0]
					+Functions.fixedLength(3)+Functions.fixedLength(sRNApos,6)+Functions.fixedLength(3)
					+ Functions.fixedLength(this.sRNAstructurePrepencity,6)+Functions.fixedLength(6)
					+Functions.fixedLength(getConsensusScore(),6)+Functions.fixedLength(6)
					+Functions.fixedLength(getIntermolecularScore(),6)+Functions.fixedLength(6)
					+Functions.fixedLength(getTotalScore(),6)
					+"\n";
		else 
			A+=Functions.fixedLength(6)+ Functions.fixedLength(3)
					+Functions.fixedLength(6)+" x-x "						
					+"\n";
//		A.print(Functions.fixedLength(getTotalScore(),6)+"\t");
//		A.print(Functions.fixedLength(getIntermolecularScore(),6)+"\t");
//		A.print(Functions.fixedLength(this.sRNAstructurePrepencity,6)+"\t");				
//		A.print(Functions.fixedLength(this.mRNAstructurePrepencity,6)+"\t");				
//		A.println(Functions.fixedLength(getConsensusScore(),6)+"\t");
/*			A.print(Functions.fixedLength(getPreviousScore(),6)+"\t");
		A.print(Functions.fixedLength(getsRNApointer(),6)+"\t");
		A.print(Functions.fixedLength(getmRNApointer(),6)+"\t");
		A.print(Functions.fixedLength(getsRNAstart(),6)+"\t");
		A.print(Functions.fixedLength(getmRNAstart(),6)+"\t");
		A.println(Functions.fixedLength(getInitiationLength(),6)+"\t");
*/			if(mRNAsequence != null){
			for(int i = 1; i < Functions.max(mRNAsequence.length, sRNAsequence.length);i++){
				A += Functions.fixedLength(16);
;
				if(mRNAsequence.length> i)
					A += mRNAsequence[i];
				else
					A +="|";
				A +="   ";
				if(sRNAsequence.length> i)
					A+=sRNAsequence[i]+"\n";
				else
					A+="|"+"\n";
			}
		}
	return A;
	}
	

	public void printFinalScore(ExtendedWriter A , int sRNApos,int mRNApos){
		A.println(Functions.fixedLength("|",6)+ Functions.fixedLength(3)
				+Functions.fixedLength(6)+"   |  "
				+Functions.fixedLength(3)+Functions.fixedLength(6)+Functions.fixedLength(3)
				+ Functions.fixedLength(7)+Functions.fixedLength("|",4)
				+Functions.fixedLength(6)+Functions.fixedLength(6)
				+Functions.fixedLength(getDanglingScore(),6)+Functions.fixedLength(6)
				+Functions.fixedLength(getFinalScore(),6)
		);
	}

	public void printFinalScore( int sRNApos,int mRNApos){
		System.out.println(Functions.fixedLength("|",6)+ Functions.fixedLength(3)
				+Functions.fixedLength(6)+"   |  "
				+Functions.fixedLength(3)+Functions.fixedLength(6)+Functions.fixedLength(3)
				+ Functions.fixedLength(7)+Functions.fixedLength("|",4)
				+Functions.fixedLength(6)+Functions.fixedLength(6)
				+Functions.fixedLength(getDanglingScore(),6)+Functions.fixedLength(6)
				+Functions.fixedLength(getFinalScore(),6)
		);
	}

	public String printFinalScore( int sRNApos,int mRNApos,String A){
		A+= Functions.fixedLength("|",6)+ Functions.fixedLength(3)
				+Functions.fixedLength(6)+"   |  "
				+Functions.fixedLength(3)+Functions.fixedLength(6)+Functions.fixedLength(3)
				+ Functions.fixedLength(7)+Functions.fixedLength("|",4)
				+Functions.fixedLength(6)+Functions.fixedLength(6)
				+Functions.fixedLength(getDanglingScore(),6)+Functions.fixedLength(6)
				+Functions.fixedLength(getFinalScore(),6)
		+"\n";
		return A;
	}
	
	public void printInteraction( int sRNApos,int mRNApos){
		if(mRNAsequence != null )
			System.out.println(Functions.fixedLength("  "+this.mRNAstructurePrepencity,6)+ Functions.fixedLength(4)
					+Functions.fixedLength(mRNApos,6)+" "+mRNAsequence[0]+"-"+sRNAsequence[0]
					+Functions.fixedLength(3)+Functions.fixedLength(sRNApos,6)+Functions.fixedLength(3)
					+ Functions.fixedLength(this.sRNAstructurePrepencity,6)+Functions.fixedLength(6)
					+Functions.fixedLength(getConsensusScore(),6)+Functions.fixedLength(6)
					+Functions.fixedLength(getIntermolecularScore(),6)+Functions.fixedLength(6)
					+Functions.fixedLength(getTotalScore(),6)
					
			);
		else 
			System.out.println(Functions.fixedLength(6)+ Functions.fixedLength(3)
					+Functions.fixedLength(6)+" x-x "						
			);
//		A.print(Functions.fixedLength(getTotalScore(),6)+"\t");
//		A.print(Functions.fixedLength(getIntermolecularScore(),6)+"\t");
//		A.print(Functions.fixedLength(this.sRNAstructurePrepencity,6)+"\t");				
//		A.print(Functions.fixedLength(this.mRNAstructurePrepencity,6)+"\t");				
//		A.println(Functions.fixedLength(getConsensusScore(),6)+"\t");
/*			A.print(Functions.fixedLength(getPreviousScore(),6)+"\t");
		A.print(Functions.fixedLength(getsRNApointer(),6)+"\t");
		A.print(Functions.fixedLength(getmRNApointer(),6)+"\t");
		A.print(Functions.fixedLength(getsRNAstart(),6)+"\t");
		A.print(Functions.fixedLength(getmRNAstart(),6)+"\t");
		A.println(Functions.fixedLength(getInitiationLength(),6)+"\t");
*/			if(mRNAsequence != null){
			for(int i = 1; i < Functions.max(mRNAsequence.length, sRNAsequence.length);i++){
				System.out.print(Functions.fixedLength(16));
				if(mRNAsequence.length> i)
					System.out.print(mRNAsequence[i]);
				else
					System.out.print("|");
				System.out.print("   ");
				if(sRNAsequence.length> i)
					System.out.println(sRNAsequence[i]);
				else
					System.out.println("|");
			}
		}
	}
	public void printDanglingEnd(){
		System.out.print("          \t");
		System.out.print(Functions.fixedLength(getFinalScore(),6)+"\t");
		System.out.println(Functions.fixedLength(getDanglingScore(),6)+"\t");
}

	



	public void changeDir(int sRNAlength, int mRNAlength){
		int tempmRNAstart = mRNAlength - 1 - sRNAstart;
		int tempmRNApointer = mRNAlength - 1 - sRNApointer;
		this.sRNAstart = sRNAlength- this.mRNAstart-1;
		this.sRNApointer = sRNAlength - this.mRNApointer-1;
		
		
		this.mRNAstart = tempmRNAstart;
		this.mRNApointer = tempmRNApointer;
		
		char[] tempSeq = this.sRNAsequence;
		this.sRNAsequence = this.mRNAsequence;
		this.mRNAsequence = tempSeq;

		float tempStruct = this.mRNAstructurePrepencity;
		this.mRNAstructurePrepencity =  this.sRNAstructurePrepencity;
		this.sRNAstructurePrepencity = tempStruct;
		
	}	
	
	

	public IntermolecularStructure(){}
	
	public void setPosition(int mRNApointer, int sRNApointer){
		this.mRNApointer = mRNApointer;
		this.sRNApointer = sRNApointer;
	}
	public void setIntramolecularPrepencity(float mRNAstructurePrepencity, float sRNAstructurePrepencity){
		this.mRNAstructurePrepencity = mRNAstructurePrepencity;
		this.sRNAstructurePrepencity = sRNAstructurePrepencity;
	}
	
	public void setWeights(float intermolecularWeight, float intramolecularWeight, float consensusWeight){
		this.intermoleucularWeight = intermolecularWeight;
		this.intramolecularWeight = intramolecularWeight;
		this.consensusWeight = consensusWeight;
	}
	
	
	public IntermolecularStructure copy(){
		IntermolecularStructure copy = new IntermolecularStructure();
		copy.setPosition(this.mRNApointer, this.sRNApointer);
		copy.setInitiationLength(this.initiationLength);
		copy.setIntermolecularScore(this.intermolecularScore);
		copy.setIntramolecularPrepencity(this.mRNAstructurePrepencity, this.sRNAstructurePrepencity);
		copy.setConsensusScore(this.consensusScore);
		copy.setWeights(this.intermoleucularWeight,this.intramolecularWeight,this.consensusWeight);
		copy.setDanglingScore(this.getDanglingScore());
		copy.setmRNAstart(this.getmRNAstart());
		copy.setsRNAstart(this.getsRNAstart());
		
		return copy;
	}
	
	public void setStructure(int mRNApointer, int sRNApointer,float intermolecularScore,
			float mRNAstructurePrepencity, float sRNAstructurePrepencity, float consensusScore){
			this.mRNApointer = mRNApointer;
			this.sRNApointer = sRNApointer;
			this.intermolecularScore = intermolecularScore;
			this.mRNAstructurePrepencity = mRNAstructurePrepencity;
			this.sRNAstructurePrepencity = sRNAstructurePrepencity;
			this.consensusScore = consensusScore;
	}
	
	
	
	
	
	public static IntermolecularStructure[] addIntermolecularStructure(IntermolecularStructure[] objects, IntermolecularStructure newIntermolecularStructure ){
		if(objects == null){
			objects = new IntermolecularStructure[1];
			objects[0] = newIntermolecularStructure;
			return objects;
		}
		IntermolecularStructure[] newIntermolecularStructures = new IntermolecularStructure[objects.length+1];
		for(int i = 0; i < objects.length;i++){
			newIntermolecularStructures[i] = objects[i];
		}
		newIntermolecularStructures[objects.length] = newIntermolecularStructure;
		return newIntermolecularStructures;
	}

	/**
	 * Returns the value of totalScore.
	 */
	public float getTotalScore()
	{
		return intermolecularScore*this.intermoleucularWeight 
				+ this.intramolecularWeight *(this.sRNAstructurePrepencity+ this.mRNAstructurePrepencity)
				+ this.consensusScore*this.consensusWeight
				+ this.previousScore;
		
	}
	
	public float getDeltaScore(){
		return intermolecularScore*this.intermoleucularWeight 
		+ this.intramolecularWeight *(this.sRNAstructurePrepencity+ this.mRNAstructurePrepencity)
		+ this.consensusScore*this.consensusWeight;
		
	}
	
	
	public float getIntermolecularScore(){
		return intermolecularScore*this.intermoleucularWeight;
	}
	
	public float getIntramolecularScore(){
		return this.intramolecularWeight *(this.sRNAstructurePrepencity+ this.mRNAstructurePrepencity);
	}
	
	public float getConsensusScore(){
		return this.consensusScore*this.consensusWeight;
	}
	
	public float getFinalScore(){
		return getTotalScore()+ this.danglingScore*this.intermoleucularWeight;	
	}
	
	
	/**
	 * Returns the value of mRNApointer.
	 */
	public int getmRNApointer()
	{
		return mRNApointer;
	}

	/**
	 * Sets the value of mRNApointer.
	 * @param mRNApointer The value to assign mRNApointer.
	 */
	public void setmRNApointer(int mRNApointer)
	{
		this.mRNApointer = mRNApointer;
	}

	/**
	 * Returns the value of sRNApointer.
	 */
	public int getsRNApointer()
	{
		return sRNApointer;
	}

	/**
	 * Sets the value of sRNApointer.
	 * @param sRNApointer The value to assign sRNApointer.
	 */
	public void setsRNApointer(int sRNApointer)
	{
		this.sRNApointer = sRNApointer;
	}

	/**
	 * Returns the value of initiationLength.
	 */
	public int getInitiationLength()
	{
		return initiationLength;
	}

	/**
	 * Sets the value of initiationLength.
	 * @param initiationLength The value to assign initiationLength.
	 */
	public void setInitiationLength(int initiationLength)
	{
		this.initiationLength = initiationLength;
	}



	/**
	 * Sets the value of intermolecularScore.
	 * @param intermolecularScore The value to assign intermolecularScore.
	 */
	public void setIntermolecularScore(float intermolecularScore)
	{
		this.intermolecularScore = intermolecularScore;
	}
/**

	/**
	 * Sets the value of consensusScore.
	 * @param consensusScore The value to assign consensusScore.
	 */
	public void setConsensusScore(float consensusScore)
	{
		this.consensusScore = consensusScore;
	}

	
	/**
	 * Returns the value of previousScore.
	 */
	public float getPreviousScore()
	{
		return previousScore;
	}

	/**
	 * Sets the value of previousScore.
	 * @param previousScore The value to assign previousScore.
	 */
	public void setPreviousScore(float previousScore)
	{
		this.previousScore = previousScore;
	}


	/**
	 * Returns the value of mRNAstart.
	 */
	public int getmRNAstart()
	{
		return mRNAstart;
	}

	/**
	 * Sets the value of mRNAstart.
	 * @param mRNAstart The value to assign mRNAstart.
	 */
	public void setmRNAstart(int mRNAstart)
	{
		this.mRNAstart = mRNAstart;
	}

	/**
	 * Returns the value of sRNAstart.
	 */
	public int getsRNAstart()
	{
		return sRNAstart;
	}

	/**
	 * Sets the value of sRNAstart.
	 * @param sRNAstart The value to assign sRNAstart.
	 */
	public void setsRNAstart(int sRNAstart)
	{
		this.sRNAstart = sRNAstart;
	}

	/**
	 * Returns the value of danglingScore.
	 */
	public float getDanglingScore()
	{
		return danglingScore;
	}

	/**
	 * Sets the value of danglingScore.
	 * @param danglingScore The value to assign danglingScore.
	 */
	public void setDanglingScore(float danglingScore)
	{
		this.danglingScore = danglingScore;
	}
		/**
	 * Returns the value of mRNAsequence.
	 */
	public char[] getmRNAsequence()
	{
		return mRNAsequence;
	}

	/**
	 * Sets the value of mRNAsequence.
	 * @param mRNAsequence The value to assign mRNAsequence.
	 */
	public void setmRNAsequence(char[] mRNAsequence)
	{
		this.mRNAsequence = mRNAsequence;
	}

	/**
	 * Returns the value of sRNAsequence.
	 */
	public char[] getsRNAsequence()
	{
		return sRNAsequence;
	}

	/**
	 * Sets the value of sRNAsequence.
	 * @param sRNAsequence The value to assign sRNAsequence.
	 */
	public void setsRNAsequence(char[] sRNAsequence)
	{
		this.sRNAsequence = sRNAsequence;
	}

	
	public boolean isInitiationStructure() {
		return initiationStructure;
	}



	public void setInitiationStructure(boolean initiationStructure) {
		this.initiationStructure = initiationStructure;
	}



	public boolean isStableStructure() {
		return stableStructure;
	}



	public void setStableStructure(boolean stableStructure) {
		this.stableStructure = stableStructure;
	}
	
	
	
	
}
