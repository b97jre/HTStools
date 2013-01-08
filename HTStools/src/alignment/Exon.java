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


public class Exon extends Object implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	protected int left;
	protected int right;
	protected String parent;
	
	Exon(int left, int right,String parent){
		this.left = left;
		this.right = right;
		this.parent = parent;
	}
}
