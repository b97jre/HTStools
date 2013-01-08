package alignment;
import general.Functions;
import general.ExtendedReader;
import general.ExtendedWriter;

import general.IOTools;
import general.RNAfunctions;
import general.energy;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;


public class mRNA extends Gene implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	protected String parent;
	protected ArrayList <Exon> exons;

	
	mRNA(int left, int right,boolean plusStrand,String ID,String Name, String parent, String description){
		this.ID = ID;
		this.Name = Name;
		this.plusStrand = plusStrand;
		this.left = left;
		this.right = right;
		this.parent = parent;
		this.description = description;
	}
	
	public boolean addExon(Exon newExon){
		if(newExon.parent.compareTo(ID) == 0){
			if(this.exons == null)
				this.exons = new ArrayList<Exon>();
			this.exons.add(newExon);
			return true;
		}
		return false;
	}
	
	
	

}
