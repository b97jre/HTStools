package alignment;

import general.ExtendedWriter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;


public class StructuralVariation {

	String ID;
	char[] major;
	char[] minor;
	double QUAL;
	String FILTER;
	String INFO;
	String FORMAT;



	Hashtable <String, StructuralVariationSample> Samples;


	public 	StructuralVariation(ArrayList<String> sampleNames, String [] StructuralVariationLine){
		// 0			1		2     3		  4			5		6		7		8		9				10			11			12				13				14				15				16		 
		// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
		// scaffold_1      5       .       A       T       56.41   .       AC=3;AF=0.188;AN=16;BaseQRankSum=-1.804;DP=72;Dels=0.00;FS=0.000;HaplotypeScore=2.4206;MLEAC=2;MLEAF=0.125;MQ=26.05;MQ0=11;MQRankSum=0.618;QD=2.69;ReadPosRankSum=1.053 GT:AD:DP:GQ:PL  0|0:14,0:14:30:0,30,268 0|0:1,0:1:3:0,3,29      0/1:13,5:18:72:72,0,247 0|0:14,0:14:39:0,39,316 0|0:13,0:13:39:0,39,328 0|0:4,0:4:3:0,3,25      1|1:2,1:3:3:23,3,0      0|0:5,0:5:3:0,3,25

		this.ID = StructuralVariationLine[2];
		this.major = StructuralVariationLine[3].toCharArray();
		this.minor = StructuralVariationLine[4].toCharArray();
		this.QUAL = Double.parseDouble(StructuralVariationLine[5]);
		this.FILTER = StructuralVariationLine[6];
		this.INFO = StructuralVariationLine[7];
		this.FORMAT = StructuralVariationLine[8];
		//TODO	
		this.Samples = new Hashtable<String, StructuralVariationSample>();

		for(int i = 0; i < sampleNames.size();i++){
			//System.out.println(sampleNames.get(i)+"\t"+StructuralVariationLine[9+i]);
			Samples.put(sampleNames.get(i), new StructuralVariationSample(StructuralVariationLine[9+i]) );
		}
	}

	
	public void addSamples(ArrayList<String> sampleNames, String [] StructuralVariationLine){
		// 0			1		2     3		  4			5		6		7		8		9				10			11			12				13				14				15				16		 
		// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
		// scaffold_1      5       .       A       T       56.41   .       AC=3;AF=0.188;AN=16;BaseQRankSum=-1.804;DP=72;Dels=0.00;FS=0.000;HaplotypeScore=2.4206;MLEAC=2;MLEAF=0.125;MQ=26.05;MQ0=11;MQRankSum=0.618;QD=2.69;ReadPosRankSum=1.053 GT:AD:DP:GQ:PL  0|0:14,0:14:30:0,30,268 0|0:1,0:1:3:0,3,29      0/1:13,5:18:72:72,0,247 0|0:14,0:14:39:0,39,316 0|0:13,0:13:39:0,39,328 0|0:4,0:4:3:0,3,25      1|1:2,1:3:3:23,3,0      0|0:5,0:5:3:0,3,25


		for(int i = 0; i < sampleNames.size();i++){
			//System.out.println(sampleNames.get(i)+"\t"+StructuralVariationLine[9+i]);
			Samples.put(sampleNames.get(i), new StructuralVariationSample(StructuralVariationLine[9+i]) );
		}
	}
	
	

	public void printSamples(ExtendedWriter EW, ArrayList<String> samples){
		EW.print(ID+"\t"+new String(this.major)+"\t"+new String(this.minor)+"\t"+this.QUAL+"\t"+this.FILTER+"\t"+this.INFO+"\t"+this.FORMAT);
		for(int i = 0; i < samples.size();i++)
		{
			EW.print("\t");
			this.Samples.get(samples.get(i)).print(EW);
		}

	}

	
	public void printSamplesRfriendly(ExtendedWriter EW, ArrayList<String> samples){
		EW.print(ID+"\t"+new String(this.major)+"\t"+new String(this.minor)+"\t"+this.QUAL+"\t"+this.FILTER+"\t"+this.INFO+"\t"+this.FORMAT);
		for(int i = 0; i < samples.size();i++)
		{
			EW.print("\t");
			this.Samples.get(samples.get(i)).printRfriendly(EW);
		}

	}
	
	
	
	public boolean isHomozygous(ArrayList<String> samples){
		for(int i = 0; i < samples.size();i++)
		{	
			if(this.Samples.containsKey(samples.get(i))){
				if(!this.Samples.get(samples.get(i)).isHomozygous()) return false;
			}else{
				System.out.println(samples.get(i)+" not found");
			}
		}
		return true;

	}






}
