package alignment;

import general.ExtendedWriter;

import java.util.ArrayList;



public class FUTR extends Gene implements Comparable<FUTR>{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	protected String parent;
	protected String dir;
	protected String name;
	protected String kind;


	
	FUTR(int left, int right){
		this.left = left;
		this.right = right;
	}
	
	FUTR(int left, int right,String parent){
		this.left = left;
		this.right = right;
		this.parent = parent;
	}

	FUTR(int left, int right,String parent,String Name){
		this.left = left;
		this.right = right;
		this.parent = parent;
	}
	
	FUTR(int left, int right,String parent,String Name,  String dir, String kind){
		this.left = left;
		this.right = right;
		this.parent = parent;
		this.name = Name;
		this.dir = dir;
		this.kind = kind;
	}
	


	FUTR(String[] bedInfo){
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

	
	public boolean join(FUTR otherExon){

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

	public boolean overlaps(FUTR otherExon){

		if(this.left > otherExon.right || this.right < otherExon.left) {
			return false;
		}return true;
	}
	
	
	
	public ArrayList<FUTR> split(FUTR otherExon){
		
		if(this.left > otherExon.right || this.right < otherExon.left) {
			System.out.println("this should not happen 2");
			return null;
		}
		ArrayList<FUTR> newExons = new ArrayList<FUTR>();
		if(this.left > otherExon.left && this.right > otherExon.right){
			int newLeft = this.left;
			int newRight = otherExon.right;
			otherExon.right = newLeft-1;
			this.left = newRight+1;
			FUTR newExon = new FUTR(newLeft,newRight, this.parent+"_overlap_"+otherExon.parent, this.name+"_overlap_"+otherExon.name, this.dir+"_overlap_"+otherExon.dir, this.kind+"_overlap_"+otherExon.kind);
			newExons.add(this);
			newExons.add(newExon);
			newExons.add(otherExon);
		} 
		else if(this.left > otherExon.left &&this.right < otherExon.right){
			int newLeft = this.left;
			int newRight = this.right;
			FUTR newExon = new FUTR(otherExon.left,newLeft-1,otherExon.parent, otherExon.name+"_left",otherExon.dir, otherExon.kind);
			FUTR newExon2 = new FUTR(newLeft,newRight, this.parent+"_overlap_"+otherExon.parent, this.name+"_overlap_"+otherExon.name, this.dir+"_overlap_"+otherExon.dir, this.kind+"_overlap_"+otherExon.kind);
			FUTR newExon3 = new FUTR(newRight+1,otherExon.right,otherExon.name+"_right", otherExon.name+"_right",otherExon.dir, otherExon.kind);
			newExons.add(newExon);
			newExons.add(newExon2);
			newExons.add(newExon3);
		}
		else if(this.left < otherExon.left && this.right < otherExon.right){
			FUTR newExon = new FUTR(this.left,otherExon.left-1,this.parent, this.name,this.dir,this.kind);
			FUTR newExon2 = new FUTR(otherExon.left,this.right, this.parent+"_overlap_"+otherExon.parent, this.name+"_overlap_"+otherExon.name, this.dir+"_overlap_"+otherExon.dir, this.kind+"_overlap_"+otherExon.kind);
			FUTR newExon3 = new FUTR(this.right+1,otherExon.right,otherExon.parent, otherExon.name,otherExon.dir,otherExon.kind);
			newExons.add(newExon);
			newExons.add(newExon2);
			newExons.add(newExon3);
		}else if(this.left < otherExon.left && this.right > otherExon.right){
			int newLeft = otherExon.left;
			int newRight = otherExon.right;
			FUTR newExon = new FUTR(this.left,newLeft-1,this.parent, this.name+"_left",this.dir, this.kind);
			FUTR newExon2 = new FUTR(newLeft,newRight,this.parent+"_overlap_"+otherExon.parent, this.name+"_overlap_"+otherExon.name, this.dir+"_overlap_"+otherExon.dir, this.kind+"_overlap_"+otherExon.kind);
			FUTR newExon3 = new FUTR(newRight+1,this.right,this.parent, this.name+"_right",this.dir, this.kind);
			newExons.add(newExon);
			newExons.add(newExon2);
			newExons.add(newExon3);
		}else if(this.left == otherExon.left && this.right > otherExon.right){
			newExons.add(this);
		}else if(this.left == otherExon.left && this.right <= otherExon.right){
			newExons.add(otherExon);
		}
		else if(this.left < otherExon.left && this.right == otherExon.right){
			newExons.add(this);
		}else if(this.left >= otherExon.left && this.right == otherExon.right){
			newExons.add(otherExon);
		}
		else{
			System.out.println("this should not happen");
			return null;
		}
		
		//System.out.println(this.right);
		return newExons;
	}
		
	public int compareTo(FUTR n) {
		int lastCmp = left-n.left;
		return (lastCmp != 0 ? lastCmp : right- n.right);
	}

	public void printBED(String chrom){
		System.out.println(chrom+"\t"+this.left+"\t"+this.right+"\t"+this.name+"\t.\t"+this.dir+"\t.\t"+this.kind+"\t.\t"+this.parent);
	}

	public void printBED(String chrom, ExtendedWriter EW){
		EW.println(chrom+"\t"+this.left+"\t"+this.right+"\t"+this.name+"\t.\t"+this.dir+"\t.\t"+this.kind+"\t.\t"+this.parent);
	}
	
	
}
