package alignment;

import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Hashtable;


import general.ExtendedReader;
import general.ExtendedWriter;
import general.Functions;

public class mirfold {


	int nrOfHits;
	String Name;
	String chromosome;
	int location;
	int penalty;
	int length; 
	double dG;
	String sequence;
	String structure;
	int start;
	int stop;

	mirfold(){

	}

	public static void main(String[] args) {
		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}
		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);
		String Dir = Functions.getValue(T, "-dir", ".");
		String File = Functions.getValue(T, "-file", "file");
		ArrayList<mirfold> miRNAs = mirfold.parseFile(Dir, File);
		System.out.println(miRNAs.size());
		miRNAs = mirfold.JoinmiRNAs(miRNAs);
		System.out.println(miRNAs.size());

//		if(T.containsKey("-solidFile")){
//			String CDir = Functions.getValue(T, "-solidDir", Dir);
//			String CFile = Functions.getValue(T, "-solidFile", "you have to specify -solidFile");
//			CfastaSequences C1 = new CfastaSequences();
//			C1.addSolidSequences(CDir,CFile);
//			miRNAs = mirfold.JoinmiRNAs(miRNAs, C1);
//		}
//		
		System.out.println(miRNAs.size());
		String outDir = Functions.getValue(T, "-outDir", Dir);
		String outFile = Functions.getValue(T, "-outFile", File+".best");
		mirfold.writeMirfoldHits(outDir, outFile, miRNAs);

	}	



	public boolean theSame(mirfold otherHit){
		if(this.chromosome.compareTo(otherHit.chromosome) == 0 && this.location == otherHit.location){
			this.nrOfHits++;
			return true;
		}
		return false;
	}

	public void parseHit(ExtendedReader ER){
		String info = ER.readLine();
		System.out.println(info);
		int index = info.indexOf("no structure found");
		if(info.indexOf("no structure found") < 0){
			//1,thug-S 1627 2192,36,4915329,21,-,150_nt_upstreamAndDownstream_up l=155nt penalty=0 dG=-94.66 (-0.611/nt) miR@118-138
			String[] info2 = info.split(",");
			this.chromosome = info2[0];
			this.Name = info2[1];
			this.nrOfHits = Integer.parseInt(info2[2]);
			this.location = Integer.parseInt(info2[3]);
			String rest = info2[6].substring(info2[6].indexOf(" "));
			index = rest.indexOf("l=")+2;
			int index2 = rest.indexOf("nt");
			this.length = Integer.parseInt(rest.substring(index,index2));
			rest = rest.substring(rest.indexOf("penalty"));

			this.penalty = Integer.parseInt(rest.substring(rest.indexOf("penalty=")+8,rest.indexOf(" dG")));
			rest = rest.substring(rest.indexOf("dG"));

			this.dG = Double.parseDouble(rest.substring(rest.indexOf("dG")+3,rest.indexOf(" (")));
			rest = rest.substring(rest.indexOf("@")+1);
			this.start = Integer.parseInt(rest.substring(0, rest.indexOf("-")));
			this.stop = Integer.parseInt(rest.substring(rest.indexOf("-")+1));
			this.sequence = ER.readLine();
			this.structure = ER.readLine();
			this.nrOfHits++;
		}
	}

	
	public void parseHitSparse(ExtendedReader ER){
		String info = ER.readLine();
		System.out.println(info);
		if(info.indexOf("no structure found") < 0 &&  info.indexOf("dG=  0.00") == -1){
			//1,thug-S 1627 2192,36,4915329,21,-,150_nt_upstreamAndDownstream_up l=155nt penalty=0 dG=-94.66 (-0.611/nt) miR@118-138
			
			String[] info2 = null;
			if(info.indexOf("_up")>0){
				this.Name = info.substring(0,info.indexOf("_up"));
				info2 = info.substring(info.indexOf("_up")).split(" ");
			}
			else if(info2[0].indexOf("_dwn")>0){
				this.Name = info.substring(0,info.indexOf("_dwn"));
				info2 = info.substring(info.indexOf("_dwn")).split(" ");
			}
			else
				System.out.println(info2[0]);
			this.chromosome = "N";
			this.nrOfHits = 1;
			this.location = 1;
		
			int index2 = info2[1].indexOf("nt");
			this.length = Integer.parseInt(info2[1].substring(2,index2));

			this.penalty = Integer.parseInt(info2[2].substring(info2[2].indexOf("penalty=")+8));
				this.dG = Double.parseDouble(info2[3].substring(info2[3].indexOf("dG")+3).trim());
			
			String rest = info2[5].substring(info2[5].indexOf("@")+1);
			this.start = Integer.parseInt(rest.substring(0, rest.indexOf("-")));
			this.stop = Integer.parseInt(rest.substring(rest.indexOf("-")+1));
			this.sequence = ER.readLine();
			this.structure = ER.readLine();
			this.nrOfHits++;
		}
	}
	
	public static ArrayList<mirfold> parseFileSparse(String Dir, String file){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(Dir+"/"+file));
			ArrayList<mirfold> miRNAs = new ArrayList<mirfold>();

			while(ER.more()){
				mirfold hit = new mirfold();
				hit.parseHitSparse(ER);
				if( hit.dG < 0  )
					miRNAs.add(hit);
			}
			
			ER.close();
			miRNAs = selectBestmiRNA(miRNAs);
			return miRNAs;

		}catch(Exception E){E.printStackTrace();}
		return null;
	}

	
	
	public static ArrayList<mirfold> parseFile(String Dir, String file){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(Dir+"/"+file));
			ArrayList<mirfold> miRNAs = new ArrayList<mirfold>();

			while(ER.more()){
				mirfold hit = new mirfold();
				hit.parseHit(ER);
				if(hit.penalty < 5 && hit.dG < -25 && hit.dG / (double)hit.length < -0.2 && hit.length > 55 )
					miRNAs.add(hit);

			}
			ER.close();
			return miRNAs;

		}catch(Exception E){E.printStackTrace();}
		return null;
	}
	
	
	
	

	public static ArrayList<mirfold> JoinmiRNAs(ArrayList<mirfold> miRNAs){
		ArrayList<mirfold> miRNAs2 = new ArrayList<mirfold>();
		for(int i = 0; i< miRNAs.size();i++){
			int pointer = 0;
			boolean found = false;
			while(miRNAs2.size()> pointer && !found){
				found = miRNAs2.get(pointer).theSame(miRNAs.get(i));
				pointer++;
			}
			if(!found)
				miRNAs2.add(miRNAs.get(i));
		}
		return miRNAs2;
	}

//	public static ArrayList<mirfold> JoinmiRNAs(ArrayList<mirfold> miRNAs, CfastaSequences miRNASequences){
//		ArrayList<mirfold> miRNAs2 = new ArrayList<mirfold>();
//		for(int i = 0; i< miRNAs.size();i++){
//			int count = 0; 
//			for(int j = 0; j < miRNASequences.size(); j++){
//				if(miRNASequences.get(j).sameLocation(miRNAs.get(i).chromosome, miRNAs.get(i).location))
//					count++;
//			}
//			miRNAs.get(i).nrOfHits = count;
//			if(count != 0)
//				miRNAs2.add(miRNAs.get(i));
//		}
//		return miRNAs2;
//	}
	
	public static ArrayList<mirfold> selectBestmiRNA(ArrayList<mirfold> miRNAs){
		for(int i = 0; i< miRNAs.size();i++){
			int count = 0; 
			for(int j = i+1; j < miRNAs.size(); j++){
				if(miRNAs.get(j).Name.compareTo(miRNAs.get(i).Name) == 0)
					if(miRNAs.get(i).penalty<=miRNAs.get(j).penalty){
						System.out.println(miRNAs.get(j).Name);
						miRNAs.remove(j);
						j = miRNAs.size();
					}
					else{
						System.out.println(miRNAs.get(i).Name);
						miRNAs.remove(i);
						i--;
						j = miRNAs.size();
					}
					
			}
		}
		return miRNAs;
	}
	
	

	public static void writeMirfoldHits(String Dir, String file, ArrayList<mirfold> miRNAs){
		try{	
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(Dir+"/"+file));

			for(int i = 0; i< miRNAs.size();i++){
				miRNAs.get(i).printHit(EW);
			}

			EW.flush();
			EW.close();

			EW = new ExtendedWriter(new FileWriter(Dir+"/"+file+".parsed"));

			for(int i = 0; i< miRNAs.size();i++){
				miRNAs.get(i).printHitParsed(EW);
			}

			EW.flush();
			EW.close();

			ArrayList<mirfold> miRNAs2 = new ArrayList<mirfold>();
			for(int i = 0; i< miRNAs.size();i++){
				if(miRNAs.get(i).penalty < 5 && miRNAs.get(i).dG < -25 && miRNAs.get(i).dG / (double)miRNAs.get(i).length < -0.4 && miRNAs.get(i).length > 55 )
					miRNAs2.add(miRNAs.get(i));
			}
			
			EW = new ExtendedWriter(new FileWriter(Dir+"/"+file+".parsed.passed"));

			for(int i = 0; i< miRNAs2.size();i++){
				miRNAs2.get(i).printHitParsed(EW);
			}

			EW.flush();
			EW.close();
			

			EW.flush();
			EW.close();


		}catch(Exception E){E.printStackTrace();}
	}


	private void printHit(ExtendedWriter EW){
		EW.println(
				this.Name+","+
				this.chromosome+","+
				this.location+","+
				this.nrOfHits+","+
				" l="+this.length+
				" penalty="+this.penalty+
				" dG="+this.dG+
				" ("+(float)(this.dG/(double)this.length)+"/nt)"+
				" miRNAstart="+this.start);
		EW.println(this.sequence);
		EW.println(this.structure);
		EW.println();


	}
	private void printHitParsed(ExtendedWriter EW){
		EW.println(
				this.Name+"\t"+
				this.chromosome+"\t"+
				this.location+"\t"+
				this.length+"\t"+
				this.penalty+"\t"+
				this.dG+"\t"+
				(float)(this.dG/(double)this.length)+"\t"+
				this.start+"\t"+
				this.stop+"\t"+
				this.sequence.subSequence(this.start-1, this.stop)+"\t"+
				this.nrOfHits);		
	}

}
