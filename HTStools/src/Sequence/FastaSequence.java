package Sequence;

import java.io.FileReader;
import java.io.Serializable;
import java.util.ArrayList;

//import sun.security.pkcs11.wrapper.Functions;

import general.ExtendedReader;
import general.ExtendedWriter;
import general.GeneticCode;
import general.RNAfunctions;
import general.Functions;



public class FastaSequence extends Object implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	String Name;
	int[] Sequence;
	public double[] nrOfHits;
	public String structure;
	public float covered;
	public int length;
	
	
	
	FastaSequence(String Name, String Sequence){
		this.Name = Name;
		this.Sequence = RNAfunctions.RNAString2Int2(Sequence);
		this.nrOfHits = new double[1];
		this.covered = 0;
		
	}
	
	FastaSequence(String Name){
		this.Name = Name;
		this.Sequence = null;
		this.nrOfHits = null;
		this.covered = 0;
		
	}
	
	
	boolean isStop(int A,int B, int C){
		if(A==4 && B==1 && C ==3) return true;//UAG
		if(A==4 && B==1 && C ==1) return true;//UAA
		if(A==4 && B==3 && C ==1) return true;//UGA
		return false;		
	}
	boolean isStart(int A,int B, int C){
		if(A==1 && B==4 && C ==3) return true;//AUG
		return false;		
	}
	
	public void findLongestORF(double cutoff,ExtendedWriter info,ExtendedWriter ORFs,ExtendedWriter F2PSeq, ExtendedWriter proteinSeq,GeneticCode GC){
		int[][] allInfo =  findLongestORFinfo();
		
		info.println(this.Name+"\t"+allInfo[0][3]+"\t"+allInfo[1].length+"\t"+this.Sequence.length+"\t"+allInfo[0][0]+"\t"+allInfo[0][1]+"\t"+allInfo[0][2]);
		ORFs.println(this.Name);
		ORFs.println(RNAfunctions.RNAInt2String(allInfo[1]));
		if(allInfo[0][3]==1){
			F2PSeq.println(this.Name);
			F2PSeq.println(RNAfunctions.RNAInt2String(this.Sequence));
			proteinSeq.println(this.Name);
			if(allInfo[0][0] != 0)
				proteinSeq.println(GC.TranslateRNAseq(allInfo[0][0], this.Sequence));
			else
				proteinSeq.println(GC.TranslateRNAseq(allInfo[0][4], this.Sequence));
		}
		else if(allInfo[0][3]==-1){
			F2PSeq.println(this.Name);
			int[] seq = RNAfunctions.getReverseComplement(this.Sequence);
			F2PSeq.println(RNAfunctions.RNAInt2String(seq));
			proteinSeq.println(this.Name);
			if(allInfo[0][0] != 0)
				proteinSeq.println(GC.TranslateRNAseq(allInfo[0][0], seq));
			else
				proteinSeq.println(GC.TranslateRNAseq(allInfo[0][4], seq));
			
			
		}
		else{
			F2PSeq.println(this.Name+"_forward");
			F2PSeq.println(RNAfunctions.RNAInt2String(this.Sequence));
			proteinSeq.println(this.Name+"_forward");
			if(allInfo[0][0] != 0)
				proteinSeq.println(GC.TranslateRNAseq(allInfo[0][0], this.Sequence));
			else
				proteinSeq.println(GC.TranslateRNAseq(allInfo[0][4]/10, this.Sequence));
			
			F2PSeq.println(this.Name+"_reverse");
			int[] seq = RNAfunctions.getReverseComplement(this.Sequence);
			F2PSeq.println(RNAfunctions.RNAInt2String(seq));
			proteinSeq.println(this.Name+"_reverse");
			if(allInfo[0][0] != 0)
				proteinSeq.println(GC.TranslateRNAseq(allInfo[0][0], seq));
			else
				proteinSeq.println(GC.TranslateRNAseq((allInfo[0][4]-(((int)allInfo[0][4]/10)*10)), seq));
		}
		
	}
	
	private int[][] findLongestORFinfo(){
		int[] forward = findLongestORFStrand(this.Sequence);
		int[] reverse = findLongestORFStrand(RNAfunctions.getReverseComplement(this.Sequence));
		int[] finalInfo = new int[5];
		int[] finalSeq = null;
		if(forward.length == reverse.length  && reverse.length == this.Sequence.length+4){
			finalInfo[0] = forward[forward.length-3];//start
			finalInfo[1] = forward[forward.length-2];//stop
			finalInfo[2] = forward[forward.length-4];//kind
			finalInfo[3] = 0;//direction
			finalInfo[4] = forward[forward.length-1]*10 + reverse[reverse.length-1];//frame
			finalSeq = removeInfo(forward);
		}
		else if(forward.length>reverse.length){
			finalInfo[0] = forward[forward.length-3];//start
			finalInfo[1] = forward[forward.length-2];//stop
			finalInfo[2] = forward[forward.length-4];//kind
			finalInfo[3] = 1;//direction
			finalInfo[4] = forward[forward.length-1];//frame
			finalSeq = removeInfo(forward);
			
		}
		else{
			finalInfo[0] = reverse[reverse.length-3];//start
			finalInfo[1] = reverse[reverse.length-2];//stop
			finalInfo[2] = reverse[reverse.length-4];//kind
			finalInfo[3] = -1;//direction
			finalInfo[4] = 0;//frame
			finalSeq = removeInfo(reverse);
		}
		int[][] allInfo = new int[2][0];
		allInfo[0] = finalInfo;
		allInfo[1] = finalSeq;
		return allInfo;
	}

	public int[] removeInfo(int[] seq){
		int[] newSeq = new int[seq.length-4];
		for(int i = 0; i < newSeq.length; i++){
			newSeq[i] = seq[i];
		}
		return newSeq;
	}
	
	public static int[][] findLongestORFinfo(String sequence){
		FastaSequence A = new FastaSequence("fool", sequence);
		
		return A.findLongestORFinfo();
	}
	
	
	
 public int[] findLongestORFStrand(int[] Sequence){
		int length = 0;
		int start = 0;
		int frame = 0;

		for(int i = 0; i < 3;i++){
			int j = 0+i;
			int tempStart = 0;
			boolean started = true;
			while(j < Sequence.length-2){
				if(isStop(Sequence[j],Sequence[j+1],Sequence[j+2]) && started){
					if(j - tempStart > length){
						start = tempStart;
						length = j+2 - tempStart;
						frame =i;
					}
					started = false;
				}
				if(!started && isStart(Sequence[j],Sequence[j+1],Sequence[j+2])){
					tempStart = j;
					started=true;
				}
				j = j+3;
			}
			if(started){
				if(Sequence.length - tempStart > length){
					start = tempStart;
					length = Sequence.length - tempStart;
					frame =i;
				}
				started = false;
			}
		}
		int[] ORF =new int[length+4];
		for(int i = start; i < start+length; i++){
			ORF[i-start] = Sequence[i];
		}
		ORF[length+1]=start;
		ORF[length+2]=start+length;
		ORF[length+3]=frame;
		if(start == 0 && length == Sequence.length){
			ORF[length]=0;
		}
		else if(start+length == Sequence.length){
			ORF[length]=1;
		}
		else if(start == 0){
			ORF[length]=2;
		}
		else
			ORF[length]=3;
		return ORF;
		
	}
	

	public FastaSequence(){
		this.Name = null;
		this.Sequence = null;
		this.nrOfHits = null;
		this.covered = 0;
		
	}
	
	public int getLength(){
		return Sequence.length;
	}
	
	public boolean hasTwoMutations(int[][][] doubleMuationsFreq, int[] RefSeq){
		int[] mutLoc = new int[3];
		int[] mutSeq = new int[3];
		int nrOfMut = 0;
		int pointer = 0;
		int total = 0;
		while(nrOfMut < 3 && pointer < RefSeq.length){
			if(Sequence[pointer] != RefSeq[pointer]){
				mutLoc[nrOfMut] = pointer;
				mutSeq[nrOfMut] = Sequence[pointer]-1;
				nrOfMut++;
			}
			pointer++;
		}
		if(nrOfMut == 2){
			doubleMuationsFreq[mutLoc[0]][mutLoc[1]][mutSeq[0]*4+mutSeq[1]]++;
			return true;
		}
		return false;
		
	}
	
	public void padSequenceBothEnds(String Sequence){
		int[] pad  = RNAfunctions.RNAString2Int(Sequence);
		int[] reversePad = RNAfunctions.getReverseComplement(pad);
		int[] newString = new int[this.Sequence.length+pad.length*2];
		int pointer = 0;
		for(int i = 0; i < pad.length;i++){
			newString[pointer] = pad[i];
			pointer++;
		}
		for(int i = 0; i < this.Sequence.length;i++){
			newString[pointer] = this.Sequence[i];
			pointer++;
		}
		for(int i = 0; i < pad.length;i++){
			newString[pointer] = reversePad[i];
			pointer++;
		}
	}
	
	public boolean hasSingleMutation(int[][] singleMuationsFreq, int[] RefSeq){
		int[] mutLoc = new int[2];
		int[] mutSeq = new int[2];
		int nrOfMut = 0;
		int pointer = 0;
		int total = 0;
		while(nrOfMut < 2 && pointer < RefSeq.length){
			if(Sequence[pointer] != RefSeq[pointer]){
				mutLoc[nrOfMut] = pointer;
				mutSeq[nrOfMut] = Sequence[pointer]-1;
				nrOfMut++;
			}
			pointer++;
		}
		if(nrOfMut == 1){
			singleMuationsFreq[mutLoc[0]][mutSeq[0]]++;
			return true;
		}
		return false;
		
	}
	
	
	public FastaSequence(String Name, int[] Sequence){
		this.Name = Name;
		this.Sequence = Sequence;
		this.nrOfHits = new double[1];
		
	}

	
	
	FastaSequence(String Name, String Sequence, int exp){
		this.Name = Name;
		this.Sequence = RNAfunctions.RNAString2Int2(Sequence);
		this.nrOfHits = new double[exp];
	}
	
	public void printFasta(ExtendedWriter EW, String ExtraName){
		EW.println(">"+ExtraName+" "+this.Name.substring(1));
		EW.println(RNAfunctions.DNAInt2String(Sequence));
	}
	
	public void printFasta(ExtendedWriter EW){
		String seqName = this.Name;
		char[] Narray = seqName.toCharArray();
		int pointer = 0;
		while(Narray[pointer] == '>')pointer++;
		seqName = seqName.substring(pointer);
		EW.println(">"+seqName);

		EW.println(RNAfunctions.DNAInt2String(this.Sequence));
	}

	public void printFastaRNA(ExtendedWriter EW){
		EW.println(">"+this.Name);
		EW.println(RNAfunctions.RNAInt2String(this.Sequence));
	}

	public void printName(ExtendedWriter EW){
		String[] temp = this.Name.split(" ");
		if(temp[0].indexOf(">")> -1)
			EW.println(temp[0].substring(1));
		else
			EW.println(temp[0]);
			
	}

	
	
	public void printNrOfHits(ExtendedWriter EW){
		EW.print(Name +"\t"+Sequence.length);
		for(int i =0 ; i < nrOfHits.length; i++)	
			EW.print("\t"+nrOfHits[i]);
		EW.println();
	}
	
	public void addHit(double nrOfHits){
		this.nrOfHits[0] += nrOfHits;
	}
	
	public void addHit(double nrOfHits,int exp){
		this.nrOfHits[exp] += nrOfHits;
	}
	
	public int[] getSequence() {
		return Sequence;
	}

	public String getSequence(int start, int stop) {
		int[] newSeq = new int[Math.abs(stop-start)+1];
		if(stop > start){
			for(int i = 0; i < newSeq.length; i++){
				newSeq[i] = this.Sequence[start-1+i];
			}
		}else{
			for(int i = 0; i < newSeq.length; i++){
				newSeq[i] = RNAfunctions.getComplementary(this.Sequence[start-1-i]);
			}
		}
		return RNAfunctions.RNAInt2String(newSeq);
	}
	
	public String getSequenceSurr(int start, int stop, int surr) {
		int[] newSeq = new int[Math.abs(stop-start)+1+surr+surr];
		if(stop > start){
			for(int i = 0; i < newSeq.length; i++){
				newSeq[i] = this.Sequence[start-1+i-surr];
			}
		}else{
			for(int i = 0; i < newSeq.length; i++){
				newSeq[i] = RNAfunctions.getComplementary(this.Sequence[start-1-i+surr]);
			}
		}
		return RNAfunctions.RNAInt2String(newSeq);
	}
	
	public float getGCcontent(){
		float count = 0;
		for(int i = 0; i< this.Sequence.length; i++){
			if(Sequence[i] == 2 || Sequence[i] == 3)count++;
		}
		return  count/(float)Sequence.length;
	}
	
	
	public int[][] getDistr(String gmapperFile){
		int[][] dist  = new int[this.Sequence.length][this.Sequence.length];
		try{
			ExtendedReader ER= new ExtendedReader(new FileReader(gmapperFile));
			while(ER.more()){
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				String[] info = ER.readLine().split("\t");
				String readname = info[0];    //>855_1396_1408_F3
				String contigname = info[1];  // adaptor + bc1 + P2 
				String strand = info[2];      //		+ 
				int contigstart = Integer.parseInt(info[3]);    // 		2
				int contigend = Integer.parseInt(info[4]);      //		25
				int readstart  = Integer.parseInt(info[5]);     // 		24
				int readend  = Integer.parseInt(info[6]);     // 		24
				int readlength = Integer.parseInt(info[7]);     //      	50
				int score = Integer.parseInt(info[8]);          //      	184
				String editstring = info[9];  //     	x10x2x8x4
				String readsequence = null;
				if(this.Name.indexOf(contigname) == 1){
					dist[contigstart-1][readlength]++;
				}
			}
		}catch(Exception E){E.printStackTrace();}
		return dist;
	}
	
	public int[] getSizeDistribution(int[] sizeDistribution,int[][] dist){
		for(int i = 0; i < dist.length;i++){
			for(int j = 0; j < sizeDistribution.length; j++){
				sizeDistribution[j] += dist[i][j];
			}
		}
		return sizeDistribution;
		
		
	}
	
	
	
	public int printDist(int[][] dist, int cutoff){
		int count = 0;
		for(int i = 0; i < this.Sequence.length;i++){
			for(int j = 0; j < this.Sequence.length; j++){
				count += dist[i][j];
				if(dist[i][j] >= cutoff){
					for(int k = 0; k < i;k++)System.out.print(".");
					for(int k = 0; k < j;k++)if(i+k < this.Sequence.length)System.out.print(RNAfunctions.RNAInt2char(this.Sequence[i+k]));
					for(int k = i+j; k < this.Sequence.length;k++)System.out.print(".");
					System.out.println("\t"+dist[i][j]);
				}
			}
		}
		return count;
	}
	
	public String getStringSequence() {
		return RNAfunctions.RNAInt2String(this.Sequence);
	}
	public void setSequence(int[] sequence) {
		Sequence = sequence;
	}
	
	public String getName() {
		return Name;
	}

	public void setName(String name) {
		Name = name;
	}
	
	public int findLocation(int[] sequence){
		int loc  = -1;
		for(int i = 0; i < this.Sequence.length - sequence.length + 1; i++ ){
			int pointer = 0;
			boolean same = true;
			while(pointer < sequence.length && same){
				if(this.Sequence[i+pointer]  != sequence[pointer]) same = false;
				pointer++;
			}
			if(same) return i;
			
		}
		return -1;
	}
	
	
	public boolean FixSurroundingSequenceLength(int length, int[] sequence){
		int location = findLocation(sequence);
		if(location > -1){
			double arg;
			double[] freqs = RNAfunctions.getFreq(this.Sequence);
			int upstream = length - location;
			int downstream = length - (this.Sequence.length-(sequence.length+location));
			if(upstream > 0)
				addSequenceUpstream(upstream, freqs);
			else
				removeUpstreamSequence(-upstream);
			if(downstream > 0)
				addSequenceDownstream(downstream,freqs);
			else
				removeDownstreamSequence(-downstream);
			return true;
		}
		else{
			System.out.println("could not find subsequence"+this.Name);
			System.out.println(RNAfunctions.RNAInt2String(this.Sequence));
			System.out.println(RNAfunctions.RNAInt2String(sequence));
			System.out.println();
		}
		return false;
		
	}
	
	public void removeUpstreamSequence(int length){
		int [] newSequence = new int[this.Sequence.length - length];
		for(int i = length ; i < this.Sequence.length; i++){
			newSequence[i-length] = this.Sequence[i];
		}
		this.Sequence = newSequence;
	}

	public void removeDownstreamSequence(int length){
		int [] newSequence = new int[this.Sequence.length - length];
		for(int i = 0 ; i < this.Sequence.length-length; i++){
			newSequence[i] = this.Sequence[i];
		}
		this.Sequence = newSequence;
	}
	
	
	public void addSequenceUpstream(int length, double freqs[]){
		int[] seq = RNAfunctions.generateSequence(length, freqs);
		int[] newSeq = new int[seq.length+this.Sequence.length];
		for(int i = 0; i < seq.length; i++) newSeq[i] = seq[i];
		for(int i = 0; i < this.Sequence.length; i++)newSeq[seq.length+i] = this.Sequence[i];
		this.Sequence = newSeq;
	}

	public void addSequenceDownstream(int length, double freqs[]){
		int[] seq = RNAfunctions.generateSequence(length, freqs);
		int[] newSeq = new int[seq.length+this.Sequence.length];
		for(int i = 0; i < this.Sequence.length; i++)newSeq[i] = this.Sequence[i];
		for(int i = 0; i < seq.length; i++) newSeq[this.Sequence.length+i] = seq[i];
		this.Sequence = newSeq;
	}
	
	public boolean sameSequence(int[] otherSequence){
		return RNAfunctions.theSame(this.Sequence, otherSequence);
	}
	

}
