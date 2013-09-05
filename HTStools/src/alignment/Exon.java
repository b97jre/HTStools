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


public class Exon extends Object implements Comparable<Exon>{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public int left;
	public int right;
	protected String parent;
	protected String dir;
	protected String name;
	protected String kind;


	Exon(int left, int right,String parent){
		this.left = left;
		this.right = right;
		this.parent = parent;
	}


	Exon(String[] bedInfo){
		//	0			1		2     3		  4			5		6		7		8		9			 
		// #CHROM(0)  Start     Stop     Name     .     dir    INFO  type    Something  XTR 	AInfo   
		// scaffold_1      767     2124    PAC:20891551.exon.3     .       -       phytozome8_0    exon    .       ID=PAC:20891551.exon.3;Parent=PAC:20891551;pacid=20891551
		this.left = Integer.parseInt(bedInfo[1]);
		this.right = Integer.parseInt(bedInfo[2]);
		this.name = bedInfo[3];
		if(bedInfo.length > 4){
		this.dir = bedInfo[5];
		this.kind = bedInfo[7];
		this.parent = bedInfo[9];
		}
	}

	
	public boolean join(Exon otherExon){

		if(this.left > otherExon.right || this.right < otherExon.left) {
			return false;
		}
		if(this.left > otherExon.left) this.left = otherExon.left;
		if(this.right < otherExon.right) this.right = otherExon.right;
		if(this.name.compareTo(otherExon.name)!= 0)
			this.name=this.name+"_"+otherExon.name;
		if(this.dir != null){
		if(this.dir.compareTo(otherExon.dir)!= 0)
			this.dir=this.dir+""+otherExon.dir;
		if(this.kind.compareTo(otherExon.kind)!= 0)
			this.kind=this.kind+""+otherExon.kind;
		if(this.parent.compareTo(otherExon.parent)!= 0)
			this.parent=this.parent+""+otherExon.parent;
		}
		return true;
	}

	public int compareTo(Exon n) {
		int lastCmp = left-n.left;
		return (lastCmp != 0 ? lastCmp : right- n.right);
	}

	public void printBED(String chrom){
		System.out.println(chrom+"\t"+this.left+"\t"+this.right+"\t"+this.name+"\t.\t"+this.dir+"\t.\t"+this.kind+"\t.\t"+this.parent);
	}

}
