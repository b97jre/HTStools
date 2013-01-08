package alignment;

import java.io.Serializable;


public class Repeat extends Gene implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	

	Repeat(){}
	
	Repeat(int left, int right,boolean plusStrand,String ID){
		this.ID = ID;
		this.Name = ID;
		this.plusStrand = plusStrand;
		this.left = left;
		this.right = right;
		kind = this.repeat;
		}

}