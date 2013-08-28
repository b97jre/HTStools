package alignment;

import general.ExtendedWriter;

public class StructuralVariationSample {
	
	char chr1;
	char chr2;
	
	int chr1Count;
	int chr2Count;
	
	boolean notPresent;
	boolean phased;
	String info;
	
	int call;
	//	0|0 = 0
	// 	0|1 = 1
	//	1|0 = 2
	//	1|1 = 3
	//	./. = 4
	
	
	public StructuralVariationSample(String info) {
		// 0/1:36,36:72:99:656,0,661 			(unphased)
		// 0|1:23,13:36:99:181,0,446:692.21 	(phased)
		// ./. 									(Not present)
		this.info = info;
		if(info.compareTo("./.")== 0){
			notPresent = true;
		}
		else{
			notPresent = false;
			String[] infoArray = info.split(":");
			if(infoArray[0].contains("|"))phased=true;
			else phased = false;
			char[] info0 = infoArray[0].toCharArray();
			this.chr1=info0[0];
			this.chr2=info0[2];

			String[] info1 = infoArray[1].split(",");
			this.chr1Count = Integer.parseInt(info1[0]);
			this.chr2Count = Integer.parseInt(info1[1]);
		}
		// TODO Auto-generated constructor stub
	}
	
	
	public boolean isHomozygous(){
		if(chr1 == chr2) return true;
		return false;
	}

	
/*	compareTwoDistributions(SampleName1, SampleName2){
		
		
		
	}
	*/
	
	
	public void print(ExtendedWriter EW){
		EW.print(info);
//		if(notPresent)EW.print("./.");
//		else{
//			if(phased){
//				
//			}
//			else{
//				
//			}
//			
//		}
	}

	
	public void printRfriendly(ExtendedWriter EW){
		//EW.print("\t"+this.samples.get(i)+"_Call\t"+this.samples.get(i)+"_Count1\t"+this.samples.get(i)+"_Count2\t"+this.samples.get(i)+"_Total\t"+this.samples.get(i)+"_fraction");

		if(this.notPresent)
			EW.print("./.\t0\t0\t0\t0");
		else if(phased){
			EW.print(this.chr1+"|"+this.chr2+"\t"+this.chr1Count+"\t"+this.chr2Count+"\t"+(this.chr1Count+this.chr2Count)+"\t"+(((double)this.chr1Count-(double)this.chr2Count)/((double)this.chr1Count+(double)this.chr2Count)));
		}else{
			EW.print(this.chr1+"/"+this.chr2+"\t"+this.chr1Count+"\t"+this.chr2Count+"\t"+(this.chr1Count+this.chr2Count)+"\t"+(((double)this.chr1Count-(double)this.chr2Count)/((double)this.chr1Count+(double)this.chr2Count)));
		}
//		if(notPresent)EW.print("./.");
//		else{
//			if(phased){
//				
//			}
//			else{
//				
//			}
//			
//		}
	}
	
	
}
