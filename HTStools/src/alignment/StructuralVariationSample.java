package alignment;

import general.ExtendedWriter;

public class StructuralVariationSample {
	
	char chr1;
	char chr2;
	
	int majorCount;
	int minorCount;
	
	boolean notPresent;
	boolean phased;
	String info;
	int phaseNr;
	
	//int call;
	//	0|0 = 0
	// 	0|1 = 1
	//	1|0 = 2
	//	1|1 = 3
	//	./. = 4
	

	
	public int[]  changeCount(StructuralVariationSample A,int[] info){
		if(A.notPresent){
			this.majorCount = 0;
			this.minorCount = 0;
			info[4]++;
		}
		else if(this.isHomozygous() && A.isHomozygous())info[0]++;
		else if(!this.isHomozygous() && !A.isHomozygous())info[1]++;
		else if(!this.isHomozygous() && A.isHomozygous())info[2]++;
		else if(this.isHomozygous() && !A.isHomozygous())info[3]++;
		this.majorCount = A.majorCount;
		this.minorCount = A.minorCount;
		return info;
	}
	
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
			this.majorCount = Integer.parseInt(info1[0]);
			this.minorCount = Integer.parseInt(info1[1]);
			this.info = "";
			for(int i = 2; i < infoArray.length;i++){
				this.info += ":"+infoArray[i];
			}
		}
		// TODO Auto-generated constructor stub
		phaseNr =0;
	}

	public boolean comparePhasing(StructuralVariationSample other) {
		if(!this.isHomozygous() && !other.isHomozygous()){
			if(this.chr1 == other.chr1 && this.chr2 == other.chr2) return true;
		}
		return false;
	}

	
	public boolean isSame(StructuralVariationSample SVS){
		if(this.chr1 == SVS.chr1 && this.chr2 == SVS.chr2) return true;
		return false;
	}
	
	
	public boolean isHomozygous(){
		if(chr1 == chr2) return true;
		return false;
	}
	
	public int getMotherCount(){
		if(chr1 == '0')	
			return this.majorCount;
		else
			return this.minorCount;
			
	}

	public int getFatherCount(){
		if(chr2 == '0')	
			return this.majorCount;
		else
			return this.minorCount;
	}

	public int addPhase(int phaseNr, boolean readPhased){
		if(!readPhased )
			this.phaseNr = phaseNr;
		else if(this.phased){
			this.phaseNr = phaseNr;
		}
		else
			this.phaseNr = phaseNr+1;

		return this.phaseNr;
	}

	
	public int addInitialPhase(int phaseNr){
			this.phaseNr = phaseNr;

		return this.phaseNr;
	}
	
/*	compareTwoDistributions(SampleName1, SampleName2){
		
		
		
	}
	*/
	
	
	public void print(ExtendedWriter EW){
//		EW.print(info);
		if(notPresent){
			if(phased){
				EW.print(this.chr1+"|"+this.chr2);
				EW.print(":"+this.majorCount+","+this.minorCount+this.info);
			}
			EW.print("./.");
		}
		else{
			if(phased){
				EW.print(this.chr1+"|"+this.chr2);
			}
			else{
				EW.print(this.chr1+"/"+this.chr2);
			}
			EW.print(":"+this.majorCount+","+this.minorCount+this.info);
		}
		
	}

	
	public void printRfriendly(ExtendedWriter EW){
		//EW.print("\t"+this.samples.get(i)+"_Call\t"+this.samples.get(i)+"_Count1\t"+this.samples.get(i)+"_Count2\t"+this.samples.get(i)+"_Total\t"+this.samples.get(i)+"_fraction");

		if(this.notPresent)
			EW.print("./.\t0\t0\t0\t0");
		else if(phased){
			EW.print(this.chr1+"|"+this.chr2+"\t"+this.majorCount+"\t"+this.minorCount+"\t"+(this.majorCount+this.minorCount)+"\t"+(((double)this.majorCount-(double)this.minorCount)/((double)this.majorCount+(double)this.minorCount)));
		}else{
			EW.print(this.chr1+"/"+this.chr2+"\t"+this.majorCount+"\t"+this.minorCount+"\t"+(this.majorCount+this.minorCount)+"\t"+(((double)this.majorCount-(double)this.minorCount)/((double)this.majorCount+(double)this.minorCount)));
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
