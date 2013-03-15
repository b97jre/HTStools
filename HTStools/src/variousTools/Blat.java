package variousTools;

import general.ExtendedWriter;
import general.Functions;
import general.IOTools;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import Sequence.FastQSequences;

public class Blat {

	String time;
	String projectDir;

	String suffix;
	String split;
	String[] sep; 


	boolean pairedEnd;
	int insertLength;

	int maxKmer;
	int minKmer;
	int step;
	int mergeKmer;
	boolean strandSpecific;

	public Blat(){
		this.split = ".";
		projectDir = time =  null;

	}

	public void blatScript( ExtendedWriter generalSbatchScript, SBATCHinfo sbatch ,String timestamp,String inFile, String referenceFile, String outDir ){

		if(!IOTools.isDir(outDir))
			IOTools.mkDir(outDir);
		if(!IOTools.isDir(outDir+"/reports"))
			IOTools.mkDir(outDir+"/reports");
		if(!IOTools.isDir(outDir+"/scripts"))
			IOTools.mkDir(outDir+"/scripts");
		int count = 0;
		try{
			String 	sbatchFileName = outDir+"/scripts/"+timestamp+"_blat.sbatch";
			generalSbatchScript.println("sbatch "+ sbatchFileName);

			ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
			sbatch.printSBATCHinfoCore(EW,outDir,timestamp,count,"blat", time);
			EW.println("module load bioinfo-tools");
			EW.println("module load blat/34");


			String [] temp =  referenceFile.split("/");
			String refName = referenceFile;
			if(temp.length > 1)
				refName = temp[temp.length-1];

			temp =  inFile.split("/");
			String inName = inFile;
			if(temp.length > 1)
				inName = temp[temp.length-1];

			EW.println("blat -out=blast8 -minScore=100 "+referenceFile+" "+inFile+" "+outDir+"/"+refName+"_"+inFile+".blat");

			EW.flush();
			EW.close();

		}catch(Exception E){E.printStackTrace();}

		count++;

	}
}

