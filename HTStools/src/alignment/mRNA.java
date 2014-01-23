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
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;


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

	public void printPremRNA( ExtendedWriter EW, int[] sequence){

		if(plusStrand){
			int[] premRNAsequence = Functions.getSubarray(sequence, left, right);
			EW.println(">premRNA,"+this.Name+","+left+"->"+right+",+");
			EW.println(RNAfunctions.RNAInt2String(premRNAsequence));
			//System.out.println(">"+this.Name+","+geneName+","+nrOfReads+","+start+","+length+",+,"+width+"_nt_upstreamAndDownstream");

		}
		else{
			int[] premRNAsequence = RNAfunctions.getComplementary(
					Functions.getSubarray(sequence, left, right));
			EW.println(">premRNA,"+this.Name+","+right+"->"+left+",-");
			EW.println(RNAfunctions.RNAInt2String(premRNAsequence));
		}

	}

	public void printmRNA( ExtendedWriter EW, int[] sequence,String ContigName){

		if(plusStrand){
			EW.println(this.Name+","+left+"->"+right+",+");
			System.out.println(Functions.fixedLength(ContigName, 20)+Functions.fixedLength(this.Name, 20)+Functions.fixedLength(left,10 )+Functions.fixedLength(right,10 )+Functions.fixedLength("+",5 ));
			Collections.sort(exons,new Comparator<Exon>() {
				public int compare(Exon exon1, Exon exon2) {
					return  exon1.left - exon2.left ;
				}
			});
			for(int i = 0;i < this.exons.size();i++){
				System.out.println(Functions.fixedLength(exons.get(i).left,10 )+Functions.fixedLength(exons.get(i).right,10 ));
				EW.print(RNAfunctions.DNAInt2String(Functions.getSubarray(sequence, exons.get(i).left-1, exons.get(i).right)));
			}
			EW.println();
			//System.out.println(">"+this.Name+","+geneName+","+nrOfReads+","+start+","+length+",+,"+width+"_nt_upstreamAndDownstream");

		}
		else{
			EW.println(this.Name+","+right+"->"+left+",-");
			Collections.sort(exons,new Comparator<Exon>() {
				public int compare(Exon exon1, Exon exon2) {
					return  exon2.left - exon1.left ;
				}
			});
			System.out.println(Functions.fixedLength(ContigName, 20)+Functions.fixedLength(this.Name, 20)+Functions.fixedLength(right,10 )+Functions.fixedLength(left,10 )+Functions.fixedLength("-",5 ));

			for(int i = 0;i < this.exons.size();i++){
				System.out.println(Functions.fixedLength(exons.get(i).right,10 )+Functions.fixedLength(exons.get(i).left,10 ));
				EW.print(RNAfunctions.DNAInt2String(
						RNAfunctions.getComplementary(
								Functions.getSubarray(sequence, exons.get(i).left-1, exons.get(i).right)
								)
						));
			}
			EW.println();
		}
	}
	public void printmRNA( ExtendedWriter EW, int[] sequence,String ContigName,Hashtable <Integer,StructuralVariation> SVs, String SampleName){


		if(plusStrand){
			EW.println(this.Name+","+left+"->"+right+",+");
			System.out.println(Functions.fixedLength(ContigName, 20)+Functions.fixedLength(this.Name, 20)+Functions.fixedLength(left,10 )+Functions.fixedLength(right,10 )+Functions.fixedLength("+",5 ));
			Collections.sort(exons,new Comparator<Exon>() {
				public int compare(Exon exon1, Exon exon2) {
					return  exon1.left - exon2.left ;
				}
			});
			for(int i = 0;i < this.exons.size();i++){
				System.out.println(Functions.fixedLength(exons.get(i).left,10 )+Functions.fixedLength(exons.get(i).right,10 ));
				EW.print(RNAfunctions.DNAInt2String(Functions.getSubarray(sequence, exons.get(i).left-1, exons.get(i).right)));
			}
			EW.println();
			//System.out.println(">"+this.Name+","+geneName+","+nrOfReads+","+start+","+length+",+,"+width+"_nt_upstreamAndDownstream");

		}
		else{
			EW.println(this.Name+","+right+"->"+left+",-");
			Collections.sort(exons,new Comparator<Exon>() {
				public int compare(Exon exon1, Exon exon2) {
					return  exon2.left - exon1.left ;
				}
			});
			System.out.println(Functions.fixedLength(ContigName, 20)+Functions.fixedLength(this.Name, 20)+Functions.fixedLength(right,10 )+Functions.fixedLength(left,10 )+Functions.fixedLength("-",5 ));

			for(int i = 0;i < this.exons.size();i++){
				System.out.println(Functions.fixedLength(exons.get(i).right,10 )+Functions.fixedLength(exons.get(i).left,10 ));
				EW.print(RNAfunctions.DNAInt2String(
						RNAfunctions.getComplementary(
								Functions.getSubarray(sequence, exons.get(i).left-1, exons.get(i).right)
								)
						));
			}
			EW.println();
		}


	}




	public void printPersonalmRNA( ExtendedWriter[] EWs, int[] sequence,String ContigName,Hashtable <Integer,
			StructuralVariation> SVs, ArrayList<Integer> SVorder,ArrayList<String> samples){

		int startPointer = 0;
		while(startPointer < SVorder.size() && SVorder.get(startPointer)< left) startPointer++;
		int stopPointer = startPointer;
		while(stopPointer < SVorder.size() && SVorder.get(stopPointer)<= right) stopPointer++;

		int[] startPhase = SVs.get(SVorder.get(startPointer)).getPhases(samples);
		int[] stopPhase = SVs.get(SVorder.get(stopPointer)).getPhases(samples);
		int[] heterozygousSites = new int[samples.size()];  

		if(plusStrand){
			Collections.sort(exons,new Comparator<Exon>() {
				public int compare(Exon exon1, Exon exon2) {
					return  exon1.left - exon2.left ;
				}
			});
			for(int i = 0;i < this.exons.size();i++){
				heterozygousSites = Functions.sum(heterozygousSites,
						exons.get(i).getNumberOfHeterozygousSites(SVs, SVorder,samples,startPointer));
				System.out.println(Functions.fixedLength(exons.get(i).left,10 )+Functions.fixedLength(exons.get(i).right,10 ));
			}

		}
		else{
			Collections.sort(exons,new Comparator<Exon>() {
				public int compare(Exon exon1, Exon exon2) {
					return  exon2.left - exon1.left ;
				}
			});

			for(int i = 0;i < this.exons.size();i++){
				System.out.println(Functions.fixedLength(exons.get(i).right,10 )+Functions.fixedLength(exons.get(i).left,10 ));
				heterozygousSites = Functions.sum(heterozygousSites,
						exons.get(i).getNumberOfHeterozygousSites(SVs, SVorder,samples,startPointer));
			}
		}


	}




	public void printPersonalmRNAInfo( String ContigName,Hashtable <Integer,
			StructuralVariation> SVs, ArrayList<Integer> SVorder, ArrayList<String> samples,ExtendedWriter info){

		int startPointer = 0;
		while(startPointer < SVorder.size() && SVorder.get(startPointer)< left) startPointer++;
		int stopPointer = startPointer;
		while(stopPointer < SVorder.size() && SVorder.get(stopPointer)<= right) stopPointer++;
		stopPointer--;
		if(stopPointer < 0 ) stopPointer = 0;
		if(startPointer == SVorder.size()){

			//			System.out.print(Functions.fixedLength(ContigName, 20)+Functions.fixedLength(this.Name, 20)+Functions.fixedLength(left,10 )+Functions.fixedLength(right,10 ));
			//			for(int i = 0; i < samples.size();i++){
			//				System.out.print(Functions.fixedLength(0, 5)+Functions.fixedLength(0, 5));
			//			}
			//			System.out.println();

			return;
		}
		int[] startPhase = SVs.get(SVorder.get(startPointer)).getPhases(samples);
		int[] stopPhase = SVs.get(SVorder.get(stopPointer)).getPhases(samples);
		int[] heterozygousSites = new int[samples.size()];  


		if(plusStrand){
			//EW.println(this.Name+","+left+"->"+right+",+");
			Collections.sort(exons,new Comparator<Exon>() {
				public int compare(Exon exon1, Exon exon2) {
					return  exon1.left - exon2.left ;
				}
			});
			for(int i = 0;i < this.exons.size();i++){
				heterozygousSites = Functions.sum(heterozygousSites,
						exons.get(i).getNumberOfHeterozygousSites(SVs, SVorder,samples,startPointer));
				//System.out.println(Functions.fixedLength(exons.get(i).left,10 )+Functions.fixedLength(exons.get(i).right,10 ));
				//EW.print(RNAfunctions.DNAInt2String(Functions.getSubarray(sequence, exons.get(i).left-1, exons.get(i).right)));
			}
			//EW.println();
			//System.out.println(">"+this.Name+","+geneName+","+nrOfReads+","+start+","+length+",+,"+width+"_nt_upstreamAndDownstream");

		}
		else{
			//EW.println(this.Name+","+right+"->"+left+",-");
			Collections.sort(exons,new Comparator<Exon>() {
				public int compare(Exon exon1, Exon exon2) {
					return  exon2.left - exon1.left ;
				}
			});
			//System.out.println(Functions.fixedLength(ContigName, 20)+Functions.fixedLength(this.Name, 20)+Functions.fixedLength(right,10 )+Functions.fixedLength(left,10 )+Functions.fixedLength("-",5 ));

			for(int i = 0;i < this.exons.size();i++){
				//System.out.println(Functions.fixedLength(exons.get(i).right,10 )+Functions.fixedLength(exons.get(i).left,10 ));
				heterozygousSites = Functions.sum(heterozygousSites,
						exons.get(i).getNumberOfHeterozygousSites(SVs, SVorder,samples,startPointer));

				//EW.print(RNAfunctions.DNAInt2String(
				//						RNAfunctions.getComplementary(
				//								Functions.getSubarray(sequence, exons.get(i).left-1, exons.get(i).right)
				//								)
				//						));
			}


			//EW.println();


		}

		info.print(Functions.fixedLength(ContigName, 20)+Functions.fixedLength(this.Name, 20)+Functions.fixedLength(left,10 )+Functions.fixedLength(right,10 ));
		for(int i = 0; i < samples.size();i++){
			info.print(Functions.fixedLength(stopPhase[i]- startPhase[i], 5)+Functions.fixedLength(heterozygousSites[i], 5));
		}
		info.println();


	}

	public void printPersonalmRNAInfo( String ContigName,Hashtable <Integer,
			StructuralVariation> SVs, ArrayList<Integer> SVorder, String Sample,ExtendedWriter info){

		int startPointer = 0;
		while(startPointer < SVorder.size() && SVorder.get(startPointer)< left) startPointer++;
		int stopPointer = startPointer;
		while(stopPointer < SVorder.size() && SVorder.get(stopPointer)<= right) stopPointer++;
		stopPointer--;
		if(stopPointer < 0 ) stopPointer = 0;

		if(startPointer == SVorder.size()){

			//			System.out.print(Functions.fixedLength(ContigName, 20)+Functions.fixedLength(this.Name, 20)+Functions.fixedLength(left,10 )+Functions.fixedLength(right,10 ));
			//			for(int i = 0; i < samples.size();i++){
			//				System.out.print(Functions.fixedLength(0, 5)+Functions.fixedLength(0, 5));
			//			}
			//			System.out.println();

			return;
		}
		int startPhase = SVs.get(SVorder.get(startPointer)).getPhase(Sample);
		int stopPhase = SVs.get(SVorder.get(stopPointer)).getPhase(Sample);
		int heterozygousSites = 0;  


		//EW.println(this.Name+","+left+"->"+right+",+");
		Collections.sort(exons,new Comparator<Exon>() {
			public int compare(Exon exon1, Exon exon2) {
				return  exon1.left - exon2.left ;
			}
		});

		int phase = startPhase;
		int pointer = startPointer;
		int motherCount = 0; 
		int fatherCount = 0;
		int leftStart = this.left;
		int nrOfHeterozygousSites = 0;
		for(int i = 0;i < this.exons.size();i++){
			int left = this.exons.get(i).left;
			int right = this.exons.get(i).right;

			while(pointer< SVorder.size() && SVorder.get(pointer)< left){ 
				pointer++;
			}
			while(pointer< SVorder.size() && SVorder.get(pointer)<= right){
				if(SVs.get(SVorder.get(pointer)).getPhase(Sample) != phase){
					info.println(Functions.fixedLength(ContigName, 20)+Functions.fixedLength(this.Name, 20)+Functions.fixedLength(phase, 20)+Functions.fixedLength(leftStart,10 )+
							Functions.fixedLength(SVorder.get(pointer-1), 10) +Functions.fixedLength(nrOfHeterozygousSites, 5)+Functions.fixedLength(motherCount, 7)+Functions.fixedLength(fatherCount, 7));
					phase = SVs.get(SVorder.get(pointer)).getPhase(Sample);
					motherCount = fatherCount = nrOfHeterozygousSites = 0;
					leftStart = SVorder.get(pointer);
				}
				if(SVs.get(SVorder.get(pointer)).isHeterozygous(Sample)){
					nrOfHeterozygousSites++;
					motherCount += SVs.get(SVorder.get(pointer)).getMotherCount(Sample);
					fatherCount += SVs.get(SVorder.get(pointer)).getFatherCount(Sample);
				}
				pointer++;
			}
			
		}
		info.println(Functions.fixedLength(ContigName, 20)+Functions.fixedLength(this.Name, 20)+Functions.fixedLength(phase, 20)+Functions.fixedLength(leftStart,10 )+
				Functions.fixedLength(this.right,10 )+Functions.fixedLength(nrOfHeterozygousSites, 5)+Functions.fixedLength(motherCount, 7)+Functions.fixedLength(fatherCount, 7));

	}


}
