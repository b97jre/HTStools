package alignment;

import general.Functions;
import general.IOTools;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.Serializable;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import sun.org.mozilla.javascript.internal.EcmaError;

import com.sun.org.apache.xpath.internal.functions.Function;

import Sequence.CfastaSequences;
import Sequence.Solid;

import general.ExtendedReader;
import general.ExtendedWriter;




public class Databases implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//private String GenomeName;
	private String Name;
	protected Hashtable <String,Database> Databases2;
	protected String[] DatabaseOrder;
	protected ArrayList <String> samples;
	protected CfastaSequences solidSequences;

	public static void main(String[] args) {
		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}

		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);
		Databases.run(T);

	}




	Databases(String Name,String DatabasesDir,String DatabasesFile ){
		this.Name = Name;
		this.Databases2 = new Hashtable<String,Database>();
		//System.out.print("Reads database files .......");
		this.loadDatabasesFile(DatabasesDir, DatabasesFile);
		//System.out.println("finished");
	}	

	Databases(String DatabasesFile ){
		this.Name = "temp";
		this.Databases2 = new Hashtable<String,Database>();
		//System.out.print("Reads database files .......");
		this.loadDatabasesFile(DatabasesFile);
		//System.out.println("finished");
	}	

	public Databases(){
		this.Databases2 = new Hashtable<String,Database>();
	}


	public void addSolidSequences(String solidDir, String solidFile){
		this.loadSolidSequences(solidDir, solidFile);
		this.mapSolidSequences();
	}

	public void addRmapperSequences(String rmapperDir, String rmapperFile){
		this.loadSolidSequences(rmapperDir, rmapperFile);
		this.mapSolidSequences();
	}





	public static void run(Hashtable<String,String> T){

		String dir = Functions.getValue(T, "-d", IOTools.getCurrentPath());
		String resultDir = Functions.getValue(T, "-resultDir",IOTools.getCurrentPath());
		String resultFile = Functions.getValue(T, "-resultFile", "test.out.vcf");
		String DatabasesDir = Functions.getValue(T, "-DatabasesDir", resultDir);
		String DatabasesFile = Functions.getValue(T, "-DatabasesFile", resultDir);
		String outDir = Functions.getValue(T, "-outDir", resultDir);


		if(T.containsKey("-loadVCFinfo")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);
			String VCF = Functions.getValue(T, "-VCF");
			A.loadVCFFile(VCF);
			A.printVCFSample(outDir, VCF+".Sample.vcf",Functions.getValue(T, "-samples"),VCF);
		}


		if(T.containsKey("-getPhasedInterVCFinfo")){
			boolean readPhased = true;
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();

			A.getDatabasesSizes(dir+"/"+file);
			String VCF = Functions.getValue(T, "-VCF");
			A.addPhase(readPhased);
			String sample = Functions.getValue(T, "-sample");
			String outFile = Functions.getValue(T, "-o",sample+".heterozygous."+VCF );
			A.loadVCFFile(dir+"/"+VCF);

			A.printVCFHeterozygousSample(dir,outFile,sample,VCF);
		}
		

		if(T.containsKey("-printTransmissionPhasedVCFfile")){
			boolean readPhased = false;
			String file = Functions.getValue(T, "-R", ".");
			String VCF2 = Functions.getValue(T, "-transmissionPhasedVCF");
			String sample = Functions.getValue(T, "-sample");
			Databases  A= new Databases();

			String referenceName = new File(file).getName();
			String vcfFileName = new File(VCF2).getName();
			if(IOTools.fileExists(file))
				A.loadDatabasesFile(file);
			else
				A.loadDatabasesFile(dir+"/"+file);

			if(IOTools.fileExists(VCF2))
				A.loadVCFFile(VCF2);
			else
				A.loadVCFFile(dir+"/"+VCF2);
			A.addPhase(!readPhased);	

			
			try{
				ExtendedWriter EW =  ExtendedWriter.getFileWriter(dir+"/"+sample+"."+vcfFileName+"_Mother.vcf");
				String extraInfo = "_"+sample+"_phased_mother.fa";

				A.printPhasedVCFSample(EW,sample,VCF2,extraInfo);
				EW.flush();
				EW.close();

				EW =  ExtendedWriter.getFileWriter(dir+"/"+sample+"."+vcfFileName+"_Father.vcf");
				extraInfo = "_"+sample+"_phased_father.fa";

				A.printPhasedVCFSample(EW,sample,VCF2,extraInfo);
				EW.flush();
				EW.close();

				EW =  ExtendedWriter.getFileWriter(dir+"/"+sample+"."+vcfFileName+".printTransmissionPhasedVCFfile.info");
				EW.println("Files generated using java -Xmx20G -jar /glob/johanr/bin/HTStools.jar ");
				EW.println("Flags to get this file");
				for (Map.Entry<String,String> entry : T.entrySet()) {
					String key = entry.getKey();
					String value = entry.getValue();
					EW.println(key+" "+value);
					// do stuff
				}
				EW.println();
				EW.println();
				EW.println("Reference :"+ new File(file).getAbsolutePath());
				EW.println("ÊSNPfile :"+ new File(VCF2).getAbsolutePath());
				EW.println();
				EW.println("File with left phased information and non-phased heterozygous sites are replaced with N");
				EW.println(dir+"/"+sample+"."+VCF2+"_Mother.vcf");
				EW.println();
				EW.println("File with right phased information and non-phased heterozygous sites are replaced with N");
				EW.println(dir+"/"+sample+"."+VCF2+"_Father.vcf");
				EW.println();
				EW.println();
				
				EW.println("Date created: "+Functions.getDateTime());
				EW.flush();
				EW.close();


			}catch(Exception E){
				E.printStackTrace();
			}
		}
	
		
		if(T.containsKey("-printTransmissionPhasedGenome")){
			boolean readPhased = false;
			String file = Functions.getValue(T, "-R", ".");
			String VCF2 = Functions.getValue(T, "-transmissionPhasedVCF");
			Databases  A= new Databases();

			String referenceName = new File(file).getName();
			if(IOTools.fileExists(file))
				A.loadDatabasesFile(file);
			else
				A.loadDatabasesFile(dir+"/"+file);

			if(IOTools.fileExists(VCF2))
				A.loadVCFFile(VCF2);
			else
				A.loadVCFFile(dir+"/"+VCF2);
			A.addPhase(!readPhased);	

			String sample = Functions.getValue(T, "-sample");
			try{
				ExtendedWriter EW =  ExtendedWriter.getFileWriter(dir+"/"+sample+"."+referenceName+".phaseByTransmission_Mother.fa");
				ExtendedWriter EW2 =  ExtendedWriter.getFileWriter(dir+"/"+sample+"."+referenceName+".phaseByTransmission_Mother.changes");
				A.printHaploGenome(EW,sample,false,EW2);
				EW.flush();
				EW.close();
				EW2.flush();
				EW2.close();


				EW =  ExtendedWriter.getFileWriter(dir+"/"+sample+"."+referenceName+".phaseByTransmission_Father.fa");
				EW2 =  ExtendedWriter.getFileWriter(dir+"/"+sample+"."+referenceName+".phaseByTransmission_Father.changes");
				A.printHaploGenome(EW,sample,true,EW2);
				EW.flush();
				EW.close();
				EW2.flush();
				EW2.close();
				
				EW =  ExtendedWriter.getFileWriter(dir+"/"+sample+"."+referenceName);
				A.printHomozygousGenome(EW,sample);
				EW.flush();
				EW.close();
			
				EW =  ExtendedWriter.getFileWriter(dir+"/"+sample+"."+referenceName+".phaseByTransmission.info");
				EW.println("Files generated using java -Xmx20G -jar /glob/johanr/bin/HTStools.jar ");
				EW.println("Flags to get this file");
				for (Map.Entry<String,String> entry : T.entrySet()) {
					String key = entry.getKey();
					String value = entry.getValue();
					EW.println(key+" "+value);
					// do stuff
				}
				EW.println();
				EW.println();
				EW.println("Reference :"+ new File(file).getAbsolutePath());
				EW.println("SNPfile :"+ new File(VCF2).getAbsolutePath());
				EW.println();
				EW.println("File with left phased information and non-phased heterozygous sites are replaced with N");
				EW.println(dir+"/"+sample+"."+referenceName+".phaseByTransmission_Mother.fa");
				EW.println();
				EW.println("File with right phased information and non-phased heterozygous sites are replaced with N");
				EW.println(dir+"/"+sample+"."+referenceName+".phaseByTransmission_Father.fa");
				EW.println();
				EW.println();
				EW.println("File with where all heterozygous sites are replaced with N");
				EW.println(dir+"/"+sample+"."+referenceName);
				EW.println();
				EW.println();
				EW.println();
				EW.println();
				
				EW.println("Date created: "+Functions.getDateTime());
				EW.flush();
				EW.close();


			}catch(Exception E){
				E.printStackTrace();
			}
		}
		
		
		if(T.containsKey("-parseMpileUpToVCF")){
			boolean readPhased = false;
			String file = Functions.getValue(T, "-R", ".");
			String VCF2 = Functions.getValue(T, "-transmissionPhasedVCF");
			Databases  A= new Databases();

			String referenceName = new File(VCF2).getName();
			System.out.println("loading reference");
			if(IOTools.fileExists(file))
				A.loadDatabasesFile(file);
			else
				A.loadDatabasesFile(dir+"/"+file);

			System.out.println("loading VCFfile");
			if(IOTools.fileExists(VCF2))
				A.loadVCFFile(VCF2);
			else
				A.loadVCFFile(dir+"/"+VCF2);
			A.addPhase(!readPhased);	
			

			String fileName = Functions.getValue(T, "-mpileupFile");
			String mpileUpFileName = new File(fileName).getName();
			String sample = Functions.getValue(T, "-sample");
			String motherSep =  Functions.getValue(T, "-mother");
			String fatherSep =  Functions.getValue(T, "-father");
			String whatCounts = Functions.getValue(T, "-whatCounts",".,");

			System.out.println("genotype "+ VCF2 +" sample "+sample+" using mpileup data from "+fileName);

			
			A.genotypemPileUpFile(fileName, sample, fatherSep, motherSep, whatCounts);
			A.printVCFSample(outDir, mpileUpFileName+"_"+referenceName, sample, VCF2);
			
			
		}
		
		
		

		if(T.containsKey("-getVCFinfo")){
			boolean readPhased = false;
			String file = Functions.getValue(T, "-R", ".");
			String VCF2 = Functions.getValue(T, "-VCF");
			Databases  A= new Databases();

			String referenceName = new File(file).getName();
			if(IOTools.fileExists(file))
				A.loadDatabasesFile(file);
			else
				A.loadDatabasesFile(dir+"/"+file);

			if(IOTools.fileExists(VCF2))
				A.loadVCFFile(VCF2);
			else
				A.loadVCFFile(dir+"/"+VCF2);
			if(T.containsKey("-readPhased"))
				A.addPhase(readPhased);	
			else if(T.containsKey("-transmissionPhased"))
				A.addPhase(!readPhased);	
			
			
			try{
				ExtendedWriter EW =  ExtendedWriter.getFileWriter(dir+"/"+VCF2+".info");
				EW.println("Contig\tHomozygousMajor\tHomozygousMinor\tAmbigousHomozygousPhased\tHeterozygous\tHeterozygousPhased\tSample");
				A.printVCFinfo(EW);
				EW.flush();
				EW.close();


			}catch(Exception E){
				E.printStackTrace();
			}
		}




		if(T.containsKey("-printHomozygousGenome")){
			boolean readPhased = false;
			String file = Functions.getValue(T, "-R", ".");
			String VCF2 = Functions.getValue(T, "-transmissionPhasedVCF");
			Databases  A= new Databases();

			String referenceName = new File(file).getName();
			if(IOTools.fileExists(file))
				A.loadDatabasesFile(file);
			else
				A.loadDatabasesFile(dir+"/"+file);

			if(IOTools.fileExists(VCF2))
				A.loadVCFFile(VCF2);
			else
				A.loadVCFFile(dir+"/"+VCF2);
			A.addPhase(!readPhased);

			String sample = Functions.getValue(T, "-sample");
			try{
				ExtendedWriter EW =  ExtendedWriter.getFileWriter(dir+"/"+sample+"."+referenceName);
				A.printHomozygousGenome(EW,sample);
				EW.flush();
				EW.close();
				EW =  ExtendedWriter.getFileWriter(dir+"/"+sample+"."+referenceName+".info");
				EW.println("File generated using java -Xmx20G -jar /glob/johanr/bin/HTStools.jar ");
				EW.println("Flags to get this file");
				for (Map.Entry<String,String> entry : T.entrySet()) {
					String key = entry.getKey();
					String value = entry.getValue();
					EW.println(key+" "+value);
					// do stuff
				}
				EW.println();
				EW.println();
				EW.println("Reference :"+ new File(file).getAbsolutePath());
				EW.println("SNPfile :"+ new File(VCF2).getAbsolutePath());
				EW.println();
				EW.println();
				EW.println("Date created: "+Functions.getDateTime());
				EW.flush();
				EW.close();


			}catch(Exception E){
				E.printStackTrace();
			}
		}



		if(T.containsKey("-comparePhasedVCFinfo")){
			boolean readPhased = true;
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();

			A.getDatabasesSizes(dir+"/"+file);

			String VCF = Functions.getValue(T, "-readPhasedVCF");
			A.loadVCFFile(dir+"/"+VCF);
			A.addPhase(readPhased);

			Databases  B= new Databases();
			B.getDatabasesSizes(dir+"/"+file);
			String VCF2 = Functions.getValue(T, "-transmissionPhasedVCF");
			B.loadVCFFile(dir+"/"+VCF2);
			B.addPhase(!readPhased);

			String sample = Functions.getValue(T, "-sample");
			try{
				ExtendedWriter EW =  ExtendedWriter.getFileWriter(dir+"/"+sample+".phased.diff");
				A.comparePhasedVCFFiles(B,sample,EW);
			}catch(Exception E){
				E.printStackTrace();
			}
		}

		if(T.containsKey("-intersectVCFinfo")){
			boolean readPhased = true;
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();

			if(IOTools.fileExists(dir+"/"+file))
				A.getDatabasesSizes(dir+"/"+file);
			else
				A.getDatabasesSizes(file);

			String VCF = Functions.getValue(T, "-readPhasedVCF");
			A.loadVCFFile(dir+"/"+VCF);
			A.addPhase(readPhased);

			Databases  B= new Databases();
			if(IOTools.fileExists(dir+"/"+file))
				B.getDatabasesSizes(dir+"/"+file);
			else
				B.getDatabasesSizes(file);
			String VCF2 = Functions.getValue(T, "-unphasedVCF");
			B.loadVCFFile(dir+"/"+VCF2);
			B.addPhase(!readPhased);

			String sample = Functions.getValue(T, "-sample");
			try{
				A.intersectTwoVCFFiles(B,sample);
				A.printVCFSample(outDir, VCF.substring(0,VCF.lastIndexOf('.'))+"."+VCF2.substring(0,VCF2.lastIndexOf('.'))+".intersect.vcf",Functions.getValue(T, "-samples"),VCF);
			}catch(Exception E){
				E.printStackTrace();
			}
		}



		if(T.containsKey("-compareVCFinfo")){
			boolean readPhased = true;
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();

			A.getDatabasesSizes(dir+"/"+file);

			String VCF = Functions.getValue(T, "-VCF1");
			A.loadVCFFile(dir+"/"+VCF);

			Databases  B= new Databases();
			B.getDatabasesSizes(dir+"/"+file);
			String VCF2 = Functions.getValue(T, "-VCF2");
			B.loadVCFFile(dir+"/"+VCF2);

			String VCF1unique = Functions.getValue(T, "-o1",VCF.substring(0,VCF.length()-3)+"unique.vcf");
			String VCF2unique = Functions.getValue(T, "-o2",VCF2.substring(0,VCF2.length()-3)+"unique.vcf");


			try{
				ExtendedWriter EW1 =  ExtendedWriter.getFileWriter(dir+"/"+VCF1unique);
				ExtendedWriter EW2 =  ExtendedWriter.getFileWriter(dir+"/"+VCF2unique);
				A.compareVCFFiles(B,EW1,EW2);
			}catch(Exception E){
				E.printStackTrace();
			}
		}




		/*		
		if(T.containsKey("-getTranscripts")){
			boolean readPhased = false;
			if(T.containsKey("-readPhased"))readPhased = true;
			if(T.containsKey("-transmissionPhased"))readPhased = false;
			if(!T.containsKey("-readPhased") && !T.containsKey("-transmissionPhased")){
				System.out.println("must contain flag if it is phased by reads or transmission");
				System.out.println("-readPhased or -transmissionPhased");				
				return;
			}
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			String annotationFile = null;
			A.loadDatabasesFile(dir+"/"+file);
			if(T.containsKey("-gff3")){
				A.initiateAnnotation();
				annotationFile = Functions.getValue(T, "-gff3", ".");
				A.addGFF3Annotation(dir,annotationFile);
			}
			if(T.containsKey("-bed")){
				System.out.println("bed annotation not written yet");
				//				String annotationFile = Functions.getValue(T, "-bed", ".");
				//				A.addBEDAnnotation(annotationFile);
			}
			if(T.containsKey("-gtf")){
				System.out.println("GTF annotation not written yet");

				//				String annotationFile = Functions.getValue(T, "-gtf", ".");
				//				A.addGTFAnnotation(annotationFile);
			}

			String VCF = Functions.getValue(T, "-VCF");
			if(!IOTools.fileExists(VCF)){
				if(IOTools.fileExists(dir+"/"+VCF)) VCF = dir+"/"+VCF;
			}


			A.loadVCFFile(VCF);
			if(load)
			A.addPhase(readPhased);

			A.printPeronsalmRNAsInfo(dir, annotationFile+".mRNA");

			A.printpreMRNAs(dir, annotationFile+".preMRNAs.fa");
			A.printmRNAs(dir, annotationFile+".mRNA.fa");
		}
		 */

		if(T.containsKey("-phaseVCFfile")){
			boolean readPhased = false;
			if(T.containsKey("-readPhased"))readPhased = true;
			if(T.containsKey("-transmissionPhased"))readPhased = false;
			if(!T.containsKey("-readPhased") && !T.containsKey("-transmissionPhased")){
				System.out.println("must contain flag if it is phased by reads or transmission");
				System.out.println("-readPhased or -transmissionPhased");				
				return;
			}
			String file = Functions.getValue(T, "-R", ".");
			if(!IOTools.fileExists(file)){
				if(IOTools.fileExists(dir+"/"+file)) file = dir+"/"+file;
				else {System.out.println(dir+"/"+file+" not found"); return;}
			}
			Databases  A= new Databases();
			String annotationFile = null;
			A.loadDatabasesFile(file);

			String VCF = Functions.getValue(T, "-VCF");
			if(!IOTools.fileExists(VCF)){
				if(IOTools.fileExists(dir+"/"+VCF)) VCF = dir+"/"+VCF;
				else {System.out.println(dir+"/"+VCF+" not found"); return;}
			}
			String sample = Functions.getValue(T, "-sample",null);



			A.loadVCFFile(VCF);
			A.addPhase(readPhased);
			System.out.println();
			String[] unphased = Functions.getValue(T, "-unphased").split(",");
			for(int i = 0; i < unphased.length;i++){
				String filePath = unphased[i];
				System.out.println(filePath);
				if(!IOTools.fileExists(filePath)){
					if(IOTools.fileExists(dir+"/"+unphased)) filePath = dir+"/"+filePath;
				}
				System.out.println(filePath);
				String fileName = new File(filePath).getName();
				String UnphasedDir = new File(filePath).getParent();
				A.phaseVCFfile(UnphasedDir,fileName,sample);
			}
			A.printPeronsalmRNAsInfo(dir, annotationFile+".mRNA");
			//			
			//			A.printpreMRNAs(dir, annotationFile+".preMRNAs.fa");
			//			A.printmRNAs(dir, annotationFile+".mRNA.fa");
		}


		if(T.containsKey("-printVCFinfoRfriendly")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);

			String VCF = Functions.getValue(T, "-VCF");
			A.loadVCFFile(VCF);



			A.printVCFSampleRfriendly(outDir, VCF+".Rfriendly",Functions.getValue(T, "-samples"));

		}



		if(T.containsKey("-addVCFinfo")){

			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();

			file = dir+"/"+file;
			A.getDatabasesSizes(file);

			String VCF = Functions.getValue(T, "-VCF");
			A.loadVCFFile(VCF);

			String addVCFFiles= Functions.getValue(T, "-addVCFfiles");

			String[] files = addVCFFiles.split(",");

			for(int i = 0; i < files.length;i++){
				String[] sampleInfo = files[i].substring(files[i].lastIndexOf("/")+1).split("\\.");
				String sampleName = sampleInfo[0]+"_"+sampleInfo[1]+"_"+sampleInfo[2];
				A.addVCFInfo(files[i],sampleName);
			}


			A.printVCFSample(outDir, VCF+".Sample.vcf",Functions.getValue(T, "-samples"),VCF);

		}



		if(T.containsKey("-filterVCFinfo")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);

			String VCF = Functions.getValue(T, "-VCF");
			A.loadVCFFile(VCF);

			String BED = Functions.getValue(T, "-BED");
			A.loadFilterBedFile(BED);
			A.mergeFilters();


			if(!T.containsKey("-complement"))
				A.filterVCFfiles();
			else
				A.filterVCFfilesComplement();

			A.printVCFSample(outDir, VCF.substring(0,VCF.lastIndexOf('.'))+"."+BED.substring((BED.lastIndexOf("/")+1),BED.lastIndexOf('.'))+".vcf",Functions.getValue(T, "-samples"),VCF);

		}


		if(T.containsKey("-gatherGeneinfo")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);

			String VCF = Functions.getValue(T, "-VCF");
			A.loadVCFFile(VCF);

			String BED = Functions.getValue(T, "-BED");
			A.loadFilterBedFile(BED);
			A.mergeFilters();


			if(!T.containsKey("-complement"))
				A.filterVCFfiles();
			else
				A.filterVCFfilesComplement();

			A.printVCFSample(outDir, VCF.substring(0,VCF.lastIndexOf('.'))+"."+BED.substring((BED.lastIndexOf("/")+1),BED.lastIndexOf('.'))+".vcf",Functions.getValue(T, "-samples"),VCF);

		}




		if(T.containsKey("-binBed")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);

			String BED = Functions.getValue(T, "-BED");
			A.loadFilterBedFile(BED);
			A.mergeFilters();
			String binFile = BED+".bin";
			int binSize = Functions.getInt(T, "-binBed", 50000);
			String newBedFile = binFile+".bed";
			double fraction = Functions.getDouble(T, "-fraction", 0.3);
			A.countBEDfilters(binSize,fraction,binFile, newBedFile);


		}

		if(T.containsKey("-binBedCount")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);

			String BED = Functions.getValue(T, "-BED");
			A.loadFilterBedFile(BED);
			A.mergeFilters();
			String binFile = BED+".bin";
			int binSize = Functions.getInt(T, "-binBed", 50000);
			String newBedFile = binFile+".bed";
			double fraction = Functions.getDouble(T, "-fraction", 0.3);
			A.countBEDfilters(binSize,fraction,binFile, newBedFile);
		}


		if(T.containsKey("-annotateSNPs")){
			String file = Functions.getValue(T, "-R", ".");
			String FullPathfile = dir+"/"+file;

			Databases  A= new Databases();
			A.getDatabasesSizes(FullPathfile);

			String VCF = Functions.getValue(T, "-VCF");
			String FullPathVCF = dir+"/"+VCF;
			A.loadVCFFile(FullPathVCF);

			String BED = Functions.getValue(T, "-BED");
			if(!IOTools.fileExists(BED))
				BED = dir+"/"+BED;

			A.loadFilterBedFile(BED);
			A.mergeFilters();	

			A.annotateVCFfiles(VCF.substring(Math.max(VCF.lastIndexOf('/')+1, 0),VCF.lastIndexOf('.'))+"."+BED.substring(Math.max(BED.lastIndexOf('/')+1, 0),BED.lastIndexOf('.'))+".VCFannotaion");

			//A.printVCFSample(outDir, VCF.substring(0,VCF.lastIndexOf('.'))+"."+BED.substring(0,BED.lastIndexOf('.'))+".vcf",Functions.getValue(T, "-samples"));

		}

		if(T.containsKey("-SkellyFormat")){

			boolean readPhased = false;
			if(T.containsKey("-readPhased"))readPhased = true;

			String file = Functions.getValue(T, "-R", ".");
			if(!IOTools.fileExists(file)){
				if(IOTools.fileExists(dir+"/"+file)) file = dir+"/"+file;
				else {System.out.println(dir+"/"+file+" not found"); return;}
			}

			Databases  A= new Databases();
			A.getDatabasesSizes(file);
			String BED = Functions.getValue(T, "-BED");
			if(!IOTools.fileExists(BED)){
				if(IOTools.fileExists(dir+"/"+BED)) BED = dir+"/"+BED;
				else {System.out.println(dir+"/"+BED+" not found"); return;}
			}
			if(!IOTools.fileExists(BED))
				BED = dir+"/"+BED;

			String VCF = Functions.getValue(T, "-VCF");
			if(!IOTools.fileExists(VCF)){
				if(IOTools.fileExists(dir+"/"+VCF)) VCF = dir+"/"+VCF;
				else {System.out.println(dir+"/"+VCF+" not found"); return;}
			}




			A.loadVCFFile(VCF);
			A.addPhase(readPhased);

			A.loadFilterBedFile(BED);
			A.mergeFilters();	



			A.SkellyFormat(dir+"/"+VCF.substring(Math.max(VCF.lastIndexOf('/')+1, 0),VCF.lastIndexOf('.'))+"."+BED.substring(Math.max(BED.lastIndexOf('/')+1, 0),BED.lastIndexOf('.')));

			//A.printVCFSample(outDir, VCF.substring(0,VCF.lastIndexOf('.'))+"."+BED.substring(0,BED.lastIndexOf('.'))+".vcf",Functions.getValue(T, "-samples"));

		}





		if(T.containsKey("-onlyHeterozygous")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);

			String VCF = Functions.getValue(T, "-VCF");
			A.loadVCFFile(VCF);

			A.removeHomozygous(Functions.getValue(T, "-sample"));
			A.printVCFSample(outDir,VCF.substring(0,VCF.lastIndexOf('.'))+".heterozygous.vcf",Functions.getValue(T, "-sample"),VCF);

		}

		if(T.containsKey("-onlyHomozygous")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);

			String VCF = Functions.getValue(T, "-VCF");
			A.loadVCFFile(VCF);

			A.removeHeterozygous(Functions.getValue(T, "-sample"));
			A.printVCFSample(outDir, VCF.substring(0,VCF.lastIndexOf('.'))+".homozygous.vcf",null,VCF);

		}



		if(T.containsKey("-sortBed")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);
			String BED = Functions.getValue(T, "-BED");
			A.loadFilterBedFile(BED);
			A.printFilters();
			A.sortFilters();
			A.printFilters();

			A= new Databases();
			A.getDatabasesSizes(file);
			BED = Functions.getValue(T, "-BED");
			A.loadFilterBedFile(BED);
			A.printFilters();
			A.mergeFilters();
			A.printFilters();
		}

		if(T.containsKey("-splitBed")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(dir+"/"+file);
			String BED = Functions.getValue(T, "-BED");
			A.loadFilterBedFile(dir+"/"+BED);
			A.splitFilters();
			String outBed = BED.substring(0,BED.lastIndexOf("."))+"split.bed"; 
			A.printFilters(dir+"/"+outBed);

		}




		if(T.containsKey("-loadGFF")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);

			String GTF = Functions.getValue(T, "-GFF", ".");
			A.loadGTFFile(GTF);
			A.printCoverage();
		}

		if(T.containsKey("-coverage")){
			String file = Functions.getValue(T, "-R", ".");
			Databases  A= new Databases();
			A.getDatabasesSizes(file);

			String GTF = Functions.getValue(T, "-GFF", ".");
			A.loadGTFFile(GTF);
			A.printCoverage();

			String VCF = Functions.getValue(T, "-VCF");
			A.loadVCFFile(VCF);
		}


		if(T.containsKey("-printDistribution")){
			Databases B = new Databases("temp", DatabasesDir, DatabasesFile);

			int cutoff = Integer.parseInt(Functions.getValue(T, "-x", "1"));
			B.loadRmapperSequences(resultDir, resultFile);
			B.mapSolidSequences();
			B.printDistribution(cutoff, outDir,resultFile);

		}

		System.out.println("End of databases");

		//runPrograms(DatabasesDir, DatabasesFile,solidDir, solidFile,resultDir, resultFile,outDir,
		//			"temp",0,0,T);
	}


	public void addGFF3Annotation(String dir, String fileName){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));

			while(ER.more()){
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					readInfo(ER.readLine());
				}
			}
		}catch(Exception E){E.printStackTrace();}
	}

	private void readInfo(String GFF3line){
		String[] columns = GFF3line.split("\t");
		String DatabaseID = columns[0];

		if(this.Databases2.containsKey(DatabaseID))
			this.Databases2.get(DatabaseID).addGFF3LineInfo(columns);
		else
			System.out.println(GFF3line);
	}





	public static void runPrograms(String DatabasesDir, String DatabasesFile, 
			String solidDir, String solidFile, String resultDir, String resultFile, String outDir,
			String experiment, int nucleotideLength,int nrOfHits,
			Hashtable<String,String> T){

		int cutoff = Integer.parseInt(Functions.getValue(T, "-cutoff", "100"));

		Databases B = new Databases(experiment, DatabasesDir, DatabasesFile);
		if(T.containsKey("-solidDir")){
			B.loadSolidSequences(solidDir, solidFile);
		}
		if(T.containsKey("-resultDir")){
			if(!experiment.contains("temp")){
				resultDir = resultDir+"/"+experiment;
				if(nucleotideLength > 0)
					resultDir = resultDir+"/"+nucleotideLength;
				if(nrOfHits > 0){
					resultFile = experiment+"."+DatabasesFile+".rmapper."+nrOfHits;
					if(T.containsKey("-extra"))
						resultFile += "."+ Functions.getValue(T, "-extra", "rmapper");
				}
			}
			B.loadRmapperSequences(resultDir, resultFile);
		}

		B.mapSolidSequences();

		if(T.containsKey("-countGroups")){
			B.solidSequences.countGroups();
		}

		if(T.containsKey("-countlocations")){
			int total = B.countLocations();
			System.out.println(experiment+"\t"+nucleotideLength+"\t"+nrOfHits+"\t"+B.getNrOfSequences()+"\t"+B.getNrOfHits2()+"\t"+total);
		}

		if(T.containsKey("-printDistribution")){
			B.printDistribution(cutoff, outDir,resultFile);
		}

		if(T.containsKey("-printMaxSequences")){
			int length = Integer.parseInt(Functions.getValue(T, "-length", "21"));
			B.printMaxSequences(outDir,resultFile,length, cutoff);
		}

	}


	private void loadSolidSequences(String dir, String file){
		this.solidSequences = new CfastaSequences();
		this.solidSequences.addSolidSequences(dir,file);
	}

	private void loadRmapperSequences(String dir, String file){
		//System.out.print("reading rmapper files.......");
		this.solidSequences = new CfastaSequences();
		this.solidSequences.addRmapperSequences(dir,file);
		//System.out.println("finished");
		//this.solidSequences.countHits();
	}

	private void removeNonredundantRmapperSequences(String dir, String file){
		this.solidSequences = new CfastaSequences();
		this.solidSequences.addRmapperSequences(dir,file);




	}


	public static int[] countGroups(String experiment, int SequenceLength, int nrOfHits, Hashtable<String,String> T){

		String resultDir = Functions.getValue(T, "-resultDir", ".");
		String resultFile = Functions.getValue(T, "-resultFile", ".");
		String DatabasesDir = Functions.getValue(T, "-DatabasesDir", resultDir);
		String DatabasesFile = Functions.getValue(T, "-DatabasesFile", resultDir);

		String[] chromosomes = null; 
		String finalDir = resultDir+"/"+experiment+"/"+SequenceLength+"/";

		String finalFile = experiment+"."+DatabasesFile+".rmapper."+nrOfHits; 
		String extra = Functions.getValue(T, "-extra", "");

		if(extra.length() > 1)
			finalFile+="."+extra;
		Databases C = new Databases("temp",DatabasesDir, DatabasesFile);
		C.loadRmapperSequences(finalDir, finalFile);
		C.mapSolidSequences();
		return C.solidSequences.countGroups();
	}




	private void mapSolidSequences(){
		for(int  i = 0 ; i < this.solidSequences.size(); i++){
			Solid hit = this.solidSequences.get(i);
			for(int j = 0; j < this.Databases2.size(); j++){
				this.Databases2.get(j).mapSolidSequence2Database(hit);
			}
		}
	}



	private int getNrOfHits(){
		int total = 0;
		for(int  i = 0; i < this.Databases2.size(); i++){
			total += Databases2.get(this.DatabaseOrder[i]).getNrOfHits();
		}
		return total;
	}


	private int getNrOfHits2(){
		int total = 0;
		return this.solidSequences.getNrOfHits();
	}

	private int getNrOfSequences(){
		return this.solidSequences.size();
	}


	private int countLocations(){
		int total = 0;
		for(int  i = 0; i < this.Databases2.size(); i++){
			total += Databases2.get(this.DatabaseOrder[i]).countLocations();
		}
		return total;
	}



	private void printDistribution(int cutoff,String outDir,String resultFile){

		int total = 0;
		try{
			if(!IOTools.isDir(outDir))
				IOTools.mkDir(outDir);
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+resultFile+".distribution.overview"));
			for(int  i = 0; i < this.Databases2.size(); i++){
				Databases2.get(this.DatabaseOrder[i]).printOverallDistribution(EW,cutoff);
			}
			EW.flush();
			EW.close();

			for(int  i = 0; i < this.Databases2.size(); i++){
				String Name = IOTools.fixFileName(Databases2.get(this.DatabaseOrder[i]).getName());
				EW = new ExtendedWriter(new FileWriter(outDir+"/"+Name+"."+resultFile+".distribution"));
				Databases2.get(this.DatabaseOrder[i]).printDistribution(EW,cutoff);
				int length = Databases2.get(this.DatabaseOrder[i]).getLength();
				if(Name.indexOf("miRNAstart=") > -1){
					//int A = Integer.parseInt(Name.substring(Name.indexOf("miRNAstart=")+11));
					if(Name.indexOf(",") != -1)
						System.out.println("plot.miRNAdist(dir,experiments,\""+ Name+".distribution\",21,0,"+length+",dir,\""+Name.substring(0,Name.indexOf(","))+"\")");
					else{
						System.out.println("plot.miRNAdist(dir,experiments,\""+ Name+".distribution\",21,0,"+length+",dir,\""+Name.substring(0,Name.indexOf("_"))+"\")");
					}

				}
				else
					System.out.println("plot.Distribution.infestans(dir,experiments,\""+ Name+"\","+length+",dir,legendInfo)");


				EW.flush();
				EW.close();
			}	
		}
		catch(Exception E){E.printStackTrace();}

	}



	private void printVCFSampleRfriendly(String outDir,String resultFile,String sampleString){
		int total = 0;
		try{
			if(!IOTools.isDir(outDir))
				IOTools.mkDir(outDir);
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+resultFile));
			// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
			EW.print("CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
			if(sampleString != null){
				this.samples = new ArrayList<String>();
				String[] sampleArray = sampleString.split(","); 
				for(int i = 0; i < sampleArray.length; i++){
					samples.add(sampleArray[i]);
				}
			}
			for(int i  = 0 ; i < this.samples.size(); i++){
				EW.print("\t"+this.samples.get(i)+"_Call\t"+this.samples.get(i)+"_Count1\t"+this.samples.get(i)+"_Count2\t"+this.samples.get(i)+"_Total\t"+this.samples.get(i)+"_fraction");
			}	
			EW.println();
			for(int  i = 0; i < this.Databases2.size(); i++){
				Databases2.get(this.DatabaseOrder[i]).printVCFinfoSamplesRfriendly(this.samples, EW);
			}
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}
	}


	private void printHeader(ExtendedWriter EW, String inFile){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(inFile));
			while(ER.more() && ER.lookAhead() == '#'){
				EW.println(ER.readLine());
			}
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}

	private void printHeaderPhased(ExtendedWriter EW, String inFile,String extraInfo){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(inFile));
			
			while(ER.more() && ER.lookAhead() == '#'){
				String inLine = ER.readLine();
				if(inLine.indexOf("##contig=<ID=")== -1 )EW.println(inLine);
				else{
					String[] info =inLine.split(",");
					
					EW.println(info[0]+extraInfo+","+info[1]);
				}
			}
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}

	
	
	private void printVCFHeterozygousSample(String outDir,String resultFile,String sample, String inFile){

		System.out.print("Writing the VCF sample to the 2"+resultFile+"........");
		try{
			if(!IOTools.isDir(outDir))
				IOTools.mkDir(outDir);
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+resultFile));
			printHeader(EW,inFile);

			// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
			EW.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+sample);
			for(int  i = 0; i < this.Databases2.size(); i++){
				Databases2.get(this.DatabaseOrder[i]).printVCFHetinfoSamples(sample, EW);
			}
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}
		System.out.println("Finished");
	}

	private void printVCFSample(String outDir,String resultFile,String sampleString, String inFile){

		System.out.print("Writing the VCF sample to the 2"+resultFile+"........");
		try{
			if(!IOTools.isDir(outDir))
				IOTools.mkDir(outDir);
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+resultFile));
			printHeader(EW,inFile);

			// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
			EW.print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
			if(sampleString != null){
				this.samples = new ArrayList<String>();
				String[] sampleArray = sampleString.split(","); 
				for(int i = 0; i < sampleArray.length; i++){
					samples.add(sampleArray[i]);
				}
			}
			for(int i  = 0 ; i < this.samples.size(); i++){
				System.out.println(this.samples.get(i));
				EW.print("\t"+this.samples.get(i));
			}
			EW.println();
			for(int  i = 0; i < this.Databases2.size(); i++){
				Databases2.get(this.DatabaseOrder[i]).printVCFinfoSamples(this.samples, EW);
			}
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}
		System.out.println("Finished");
	}

	
	private void printPhasedVCFSample(ExtendedWriter EW,String sample, String inFile,String extraInfo){

		printHeaderPhased(EW,inFile,extraInfo);
			
			// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
		EW.print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
		this.samples = new ArrayList<String>();
		samples.add(sample);
		for(int i  = 0 ; i < this.samples.size(); i++){
			System.out.println(this.samples.get(i));
			EW.print("\t"+this.samples.get(i));
		}
		EW.println();
		for(int  i = 0; i < this.Databases2.size(); i++){
			Databases2.get(this.DatabaseOrder[i]).printVCFinfoSamples(this.samples, EW,extraInfo);
		}
		EW.flush();
		EW.close();
	}



	private void printFilters(){
		for(int  i = 0; i < this.Databases2.size(); i++){
			Databases2.get(this.DatabaseOrder[i]).printFilters();
		}
	}

	private void printFilters(String file){
		try{
			ExtendedWriter EW = ExtendedWriter.getFileWriter(file);
			for(int  i = 0; i < this.Databases2.size(); i++){
				Databases2.get(this.DatabaseOrder[i]).printFilters(EW);
			}
		}catch(Exception E){
			E.printStackTrace();
		}
	}


	private void sortFilters(){
		int total = 0;
		for(int  i = 0; i < this.Databases2.size(); i++){
			Databases2.get(this.DatabaseOrder[i]).sortBEDfilters();
		}

	}

	private void mergeFilters(){
		int total = 0;
		for(int  i = 0; i < this.Databases2.size(); i++){
			Databases2.get(this.DatabaseOrder[i]).mergeBEDfilters();
		}
	}

	private void splitFilters(){
		int total = 0;
		for(int  i = 0; i < this.Databases2.size(); i++){
			Databases2.get(this.DatabaseOrder[i]).splitBEDfilters();
		}
	}



	private void countBEDfilters(int size, double fraction, String binFile, String bedFile){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(binFile)); 
			ExtendedWriter BedEW = new ExtendedWriter(new FileWriter(bedFile));
			for(int  i = 0; i < this.Databases2.size(); i++){
				Databases2.get(this.DatabaseOrder[i]).countBEDfilters(size,fraction, EW,BedEW);
			}
			EW.flush();
			EW.close();
			BedEW.flush();
			BedEW.close();
		}catch(Exception E){E.printStackTrace();}
	}


	private void filterVCFfiles(){
		System.out.print("Keeping reads that are within the regioins reported in the bed file......");
		for(int  i = 0; i < this.Databases2.size(); i++){
			Databases2.get(this.DatabaseOrder[i]).filterVCFinfoOutside();
		}
		System.out.println("Finished");
	}


	private void filterVCFfilesComplement(){
		System.out.print("Keeping reads that are not within the regioins reported in the bed file......");
		for(int  i = 0; i < this.Databases2.size(); i++){
			Databases2.get(this.DatabaseOrder[i]).filterVCFinfoInside();
		}
		System.out.println("Finished");
	}


	private void annotateVCFfiles(String VCFfile){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(VCFfile));
			for(int  i = 0; i < this.Databases2.size(); i++){
				Databases2.get(this.DatabaseOrder[i]).annotateVCFinfo(EW);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){
			E.printStackTrace();
		}
	}

	private void SkellyFormat(String SkellyFile){
		try{
			ExtendedWriter[] EW = new ExtendedWriter[samples.size()];
			for(int j = 0; j < samples.size();j++){
				EW[j] = ExtendedWriter.getFileWriter(SkellyFile+"_"+samples.get(j)+".Skelly");
				for(int  i = 0; i < this.Databases2.size(); i++){
					Databases2.get(this.DatabaseOrder[i]).SkellyFormat(EW[j],samples.get(j));
				}
				EW[j].flush();
				EW[j].close();
			}
		}catch(Exception E){
			E.printStackTrace();
		}
	}



	private void removeHomozygous(String sampleString){
		ArrayList<String> samples2 =null;
		if(sampleString != null){
			samples2 = new ArrayList<String>();
			String[] sampleArray = sampleString.split(","); 
			for(int i = 0; i < sampleArray.length; i++){
				samples2.add(sampleArray[i]);
			}
		}else{
			samples2 = this.samples;
		}
		for(int  i = 0; i < this.Databases2.size(); i++){
			Databases2.get(this.DatabaseOrder[i]).removeHomozygous(samples2);
		}

	}

	private void removeHeterozygous(String sampleString){
		ArrayList<String> samples2 =null;
		if(sampleString != null){
			samples2 = new ArrayList<String>();
			String[] sampleArray = sampleString.split(","); 
			for(int i = 0; i < sampleArray.length; i++){
				samples2.add(sampleArray[i]);
			}
		}else{
			samples2 = this.samples;
		}
		for(int  i = 0; i < this.Databases2.size(); i++){
			Databases2.get(this.DatabaseOrder[i]).removeHeterozygous(samples2);
		}

	}



	//
	//	private void printVCFSampleExons(String outDir,String resultFile,String sample,){
	//
	//		int total = 0;
	//		try{
	//			if(!IOTools.isDir(outDir))
	//				IOTools.mkDir(outDir);
	//			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+resultFile));
	//			for(int  i = 0; i < this.Databases2.size(); i++){
	//				Databases2.get(this.DatabaseOrder[i]).printVCFinfoSample(sample, EW);
	//			}
	//			EW.flush();
	//			EW.close();
	//		}
	//		catch(Exception E){E.printStackTrace();}
	//
	//	}
	//


	private void printVCFDistribuionSample(String outDir,String resultFile,String sample){

		int total = 0;
		try{
			if(!IOTools.isDir(outDir))
				IOTools.mkDir(outDir);
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+resultFile));
			ArrayList<String> samples = new ArrayList<String>();
			for(int  i = 0; i < this.Databases2.size(); i++){
				Databases2.get(this.DatabaseOrder[i]).printVCFinfoSamples(samples, EW);
			}
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}

	}



	private void printMaxSequences(String outDir,String resultFile, int length, int cutoff){

		try{
			if(!IOTools.isDir(outDir))
				IOTools.mkDir(outDir);
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+resultFile+".maxSequences.overview"));
			for(int  i = 0; i < this.Databases2.size(); i++){
				Databases2.get(this.DatabaseOrder[i]).printMaxSequence(EW, length, cutoff);
			}
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}

	}


	public void compareDifferentRuns(Databases otherRun, ExtendedWriter ER,String exp1, String exp2){
		int nrOfHits = this.getNrOfHits();
		System.out.println("Number of hits in first genome : "+nrOfHits);
		int otherNrOfHits = otherRun.getNrOfHits();
		System.out.println("Number of hits in second genome: "+otherNrOfHits);

		double ratio = (double)nrOfHits/(double)otherNrOfHits;
		System.out.println("The ratio is: "+ratio);
		ER.println("Name,chromosome,ratio,"+exp1+","+exp2+",location,nrOfHits in "+ exp1+" = "+nrOfHits+",nrOfHits in "+ exp2+" = "+otherNrOfHits);

		for(int  i = 0; i < this.Databases2.size(); i++){
			Databases2.get(this.DatabaseOrder[i]).compareDistribution(otherRun.Databases2.get(i),ER);
		}
	}

	public void loadDatabasesFile(String dir, String fileName){
		loadDatabasesFile(dir+"/"+fileName);
	}


	public void loadDatabasesFile( String fileName){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(fileName));
			int count = 0; 
			int progress = 1000;
			while(ER.more()){
				while(ER.lookAhead() != '>')ER.skipLine();
				String DatabaseName = ER.readLine().substring(1);
				this.DatabaseOrder = Functions.addString(DatabaseOrder, DatabaseName);
				this.Databases2.put(DatabaseName,new Database(DatabaseName,ER));
				count++;
				if(count > progress){
					System.out.println("sequences read :"+ progress);
					progress += 1000;
				}
			}
		}catch(Exception E){E.printStackTrace();}
	}

	public void getDatabasesSizes(String fileName){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(fileName));
			int count = 0; 
			int progress = 1000;
			while(ER.more()){
				while(ER.lookAhead() != '>')ER.skipLine();
				String DatabaseName = ER.readLine().substring(1);
				Database A = new Database(DatabaseName);
				A.getChromosomeSequenceSize(ER);
				this.Databases2.put(DatabaseName,A);
				this.DatabaseOrder = Functions.addString(this.DatabaseOrder, DatabaseName);
				count++;
				if(count > progress){
					System.out.println("sequences read :"+ progress);
					progress += 1000;
				}
			}
			System.out.println("Reading fastafile  :"+ fileName);
			System.out.println("Number of contigs  :"+ this.Databases2.size());


		}catch(Exception E){E.printStackTrace();}
	}


	public void loadGTFFile(String dir, String fileName){
		loadGTFFile(dir+"/"+fileName);
	}

	public void loadGTFFile(String fileName){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(fileName));
			int count = 0;
			int progress = 1000;
			while(ER.more()){
				count++;
				readGTFInfo(ER);
				if(count > progress ){
					progress += 1000;
					System.out.println("gtf transcripts read: "+progress );
				}
			}
		}catch(Exception E){E.printStackTrace();}
	}


	private void readGTFInfo(ExtendedReader ER){
		String GFF3line = ER.readLine();
		String[] columns = GFF3line.split("\t");
		if(columns.length > 2){
			if(columns[2].indexOf("transcript") == 0){
				addCoverage(columns);
			}
		}
		else
			System.out.println(GFF3line);
	}


	public void loadVCFFile(String fileName){
		//	0			1		2     3		  4			5		6		7		8		9				10			11			12				13				14				15				16		 
		// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
		// scaffold_1      5       .       A       T       56.41   .       AC=3;AF=0.188;AN=16;BaseQRankSum=-1.804;DP=72;Dels=0.00;FS=0.000;HaplotypeScore=2.4206;MLEAC=2;MLEAF=0.125;MQ=26.05;MQ0=11;MQRankSum=0.618;QD=2.69;ReadPosRankSum=1.053 GT:AD:DP:GQ:PL  0|0:14,0:14:30:0,30,268 0|0:1,0:1:3:0,3,29      0/1:13,5:18:72:72,0,247 0|0:14,0:14:39:0,39,316 0|0:13,0:13:39:0,39,328 0|0:4,0:4:3:0,3,25      1|1:2,1:3:3:23,3,0      0|0:5,0:5:3:0,3,25

		String info = null;
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(fileName));
			int count = 0;
			int DatabasePointer = 0;
			int progress = 1000;
			info = ER.readLine();
			while(ER.more() && ER.lookAhead() == '#'){
				info = ER.readLine();
			}
			String [] infoArray = info.split("\t");
			this.samples = new ArrayList<String>();
			for(int i = 9 ; i < infoArray.length;i++){
				samples.add(infoArray[i]);
				System.out.println(infoArray[i]);
			}

			System.out.print("Now gathering vcf data from :"+this.DatabaseOrder[DatabasePointer]);
			while(ER.more()){
				info = ER.readLine();
				String[] LineInfo = info.split("\t");
				while(DatabasePointer < this.Databases2.size() && this.DatabaseOrder[DatabasePointer].compareTo(LineInfo[0])!=0){
					System.out.println("Finished");
					DatabasePointer++;
					if(DatabasePointer < this.Databases2.size())
						System.out.print("Now gathering vcf data for :"+this.DatabaseOrder[DatabasePointer]+".");
				}
				if(DatabasePointer > this.Databases2.size())
					System.out.println(LineInfo[0]+"\t"+LineInfo[1]);
				else
					this.Databases2.get(this.DatabaseOrder[DatabasePointer]).addVCFinfo(samples, LineInfo);
				//				}
				count++;
				if(count > progress ){
					progress += 10000;
					System.out.print(".");
				}
			}
			System.out.println("Finished");
		}catch(Exception E){
			E.printStackTrace();
			System.out.println(info);
		}
	}

	public void genotypemPileUpFile(String fileName, String sample, String father,String mother, String WhatCounts){
		//	0			1		2     3		  4			5		6		7		8		9				10			11			12				13				14				15				16		 
		String info = null;
			// remove counts on all positions
			for(int i = 0; i <this.DatabaseOrder.length; i++){	
				this.Databases2.get(this.DatabaseOrder[i]).removeCounts(sample);
			}
			// Add count for father
			genotypemPileUpFile(fileName,sample,mother,true,WhatCounts);
			// Add count for mother
			genotypemPileUpFile(fileName,sample,father,false,WhatCounts);
	}

	public void genotypemPileUpFile(String fileName, String sample, String sep , boolean mother,String WhatCounts){
		ExtendedReader ER = ExtendedReader.getFileReader(fileName);
		int DatabasePointer = 0;
		System.out.print("Now gathering mPileUpData data from :"+this.DatabaseOrder[DatabasePointer] +" with suffix "+sep);			

		int pointerCheck = 1000000;
		int pointer = 0;
		while(ER.more()){
			String info = ER.readLine();
			String[] LineInfo = info.split("\t");
			if(LineInfo[0].indexOf(sep) > 0){
				//scaffold_1_Inter4-1_phased_mother.fa	773	T	28	...................,,,,....,	JGJDIFDHJJHJFIJJIJJDDD@FFF@@
				pointer ++;
				if(pointer >pointerCheck){
					System.out.print(".");
					pointerCheck += 100000;
				}
					
				String DatabaseName = LineInfo[0].substring(0,LineInfo[0].indexOf(sep));
				int count = Functions.countOccurrences(LineInfo[4], WhatCounts);
				int location = Integer.decode(LineInfo[1]);
				while(DatabasePointer < this.Databases2.size() && this.DatabaseOrder[DatabasePointer].compareTo(DatabaseName)!=0 ){
					System.out.println("....Finished");
					DatabasePointer++;
					while(DatabasePointer < this.Databases2.size() && this.Databases2.get(this.DatabaseOrder[DatabasePointer]).SVs== null){
						System.out.println("Now gathering mPileUpData data from :"+this.DatabaseOrder[DatabasePointer] +" with suffix "+sep+"....Finished");
						DatabasePointer++;
					}
					if(DatabasePointer < this.Databases2.size()){
						System.out.print("Now gathering mPileUpData data from :"+this.DatabaseOrder[DatabasePointer] +" with suffix "+sep);
					}
					
				}

				if(DatabasePointer < this.Databases2.size())
					this.Databases2.get(this.DatabaseOrder[DatabasePointer]).addCounts(location,sample,count,mother);
			}
		}
		System.out.println("Finished");
	}
	
	
	


	public void compareVCFFiles(Databases OtherDataBase,  ExtendedWriter EW, ExtendedWriter EW2){
		for(int i = 0; i < this.DatabaseOrder.length ; i++){
			if(OtherDataBase.Databases2.containsKey(DatabaseOrder[i])){
				this.Databases2.get(DatabaseOrder[i]).compareVCFinfo(OtherDataBase.Databases2.get(DatabaseOrder[i]),this.samples, EW,EW2);
			}
		}

	}


	public void comparePhasedVCFFiles(Databases OtherDataBase, String sample, ExtendedWriter EW){
		EW.println("Name\tcurrentPhase\tstart\tstop\tnrOfSNPs\tnrOfIdentical");
		for(int i = 0; i < this.DatabaseOrder.length ; i++){
			if(OtherDataBase.Databases2.containsKey(DatabaseOrder[i])){
				this.Databases2.get(DatabaseOrder[i]).comparePhasedVCFinfo(OtherDataBase.Databases2.get(DatabaseOrder[i]),sample, EW);
			}
		}

	}

	public void intersectTwoVCFFiles(Databases OtherDataBase, String sample){
		for(int i = 0; i < this.DatabaseOrder.length ; i++){
			if(OtherDataBase.Databases2.containsKey(DatabaseOrder[i])){
				this.Databases2.get(DatabaseOrder[i]).trimPhasedVCFinfo(OtherDataBase.Databases2.get(DatabaseOrder[i]),sample);
			}
		}

	}


	public void phaseVCFfile(String dir, String fileName,String sample){
		//	0			1		2     3		  4			5		6		7		8		9				10			11			12				13				14				15				16		 
		// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
		// scaffold_1      5       .       A       T       56.41   .       AC=3;AF=0.188;AN=16;BaseQRankSum=-1.804;DP=72;Dels=0.00;FS=0.000;HaplotypeScore=2.4206;MLEAC=2;MLEAF=0.125;MQ=26.05;MQ0=11;MQRankSum=0.618;QD=2.69;ReadPosRankSum=1.053 GT:AD:DP:GQ:PL  0|0:14,0:14:30:0,30,268 0|0:1,0:1:3:0,3,29      0/1:13,5:18:72:72,0,247 0|0:14,0:14:39:0,39,316 0|0:13,0:13:39:0,39,328 0|0:4,0:4:3:0,3,25      1|1:2,1:3:3:23,3,0      0|0:5,0:5:3:0,3,25


		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(dir+"/"+fileName));
			ExtendedWriter EW = ExtendedWriter.getFileWriter(fileName.substring(0,fileName.length()-3)+"phased.vcf");
			int count = 0;
			int DatabasePointer = 0;
			int progress = 1000;
			String info = ER.readLine();
			while(ER.more() && ER.lookAhead() == '#'){
				info = ER.readLine();
				EW.println(info);
			}
			EW.println(info);
			String [] infoArray = info.split("\t");
			if(sample == null)
				sample = infoArray[9];
			int[] phaseInfo = new int[5]; 
			System.out.println("Now phasing vcf data for sample "+sample +"from file "+fileName);
			System.out.println();
			System.out.println ("Now phasing vcf data for :"+this.DatabaseOrder[DatabasePointer]);
			while(ER.more()){
				String[] LineInfo = ER.readLine().split("\t");
				while(DatabasePointer < this.Databases2.size() && this.DatabaseOrder[DatabasePointer].compareTo(LineInfo[0])!=0){
					System.out.println("Finished");
					DatabasePointer++;
					if(DatabasePointer < this.Databases2.size())
						System.out.print("Now phasing vcf data for :"+this.DatabaseOrder[DatabasePointer]);
				}
				//				if(DatabasePointer > this.Databases2.size())System.out.println(LineInfo[0]+"\t"+LineInfo[1]);
				//				else{
				this.Databases2.get(this.DatabaseOrder[DatabasePointer]).addVCFphase(sample, LineInfo, EW, phaseInfo);
				//				}
				count++;
				if(count > progress ){
					progress += 10000;
					System.out.print(".");
				}
			}
			System.out.println("Finished");

			System.out.println("PhaseInfo");
			System.out.println(phaseInfo[0]+"\t"+phaseInfo[1]+"\t"+phaseInfo[2]+"\t"+phaseInfo[3]+"\t"+phaseInfo[4]);
			ER.close();
			EW.flush();
			EW.close();


		}catch(Exception E){E.printStackTrace();}
	}


	public void markSubset(String fileName){
		//	0			1		2     3		  4			5		6		7		8		9				10			11			12				13				14				15				16		 
		// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
		// scaffold_1      5       .       A       T       56.41   .       AC=3;AF=0.188;AN=16;BaseQRankSum=-1.804;DP=72;Dels=0.00;FS=0.000;HaplotypeScore=2.4206;MLEAC=2;MLEAF=0.125;MQ=26.05;MQ0=11;MQRankSum=0.618;QD=2.69;ReadPosRankSum=1.053 GT:AD:DP:GQ:PL  0|0:14,0:14:30:0,30,268 0|0:1,0:1:3:0,3,29      0/1:13,5:18:72:72,0,247 0|0:14,0:14:39:0,39,316 0|0:13,0:13:39:0,39,328 0|0:4,0:4:3:0,3,25      1|1:2,1:3:3:23,3,0      0|0:5,0:5:3:0,3,25

		int count = 0; 
		int progress = 10000;

		try{
			ExtendedReader ER = ExtendedReader.getFileReader(fileName);
			while(ER.more() && ER.lookAhead() == '#'){
				ER.readLine();
			}
			int DatabasePointer = 0;

			while(ER.more()){
				String[] LineInfo = ER.readLine().split("\t");
				while(DatabasePointer < this.Databases2.size() && this.DatabaseOrder[DatabasePointer].compareTo(LineInfo[0])!=0){
					System.out.println("Finished");
					DatabasePointer++;
					if(DatabasePointer < this.Databases2.size())
						System.out.print("Now phasing vcf data for :"+this.DatabaseOrder[DatabasePointer]);
				}
				this.Databases2.get(this.DatabaseOrder[DatabasePointer]).markSubset(LineInfo);
				count++;
				if(count > progress ){
					progress += 10000;
					System.out.print(".");
				}
			}
			System.out.println("Finished");

		}catch(Exception E){E.printStackTrace();}
	}


	public void addVCFInfo(String fileName,String sampleName){
		//	0			1		2     3		  4			5		6		7		8		9				10			11			12				13				14				15				16		 
		// #CHROM(0)  POS     ID     REF     ALT(0)     QUAL    FILTER  INFO    FORMAT  Cr1GR1-2-KS3    Cr_39_1 Inter3-1        Inter4-1        Inter5-1        Intra6-3        Intra7-2        Intra8-2
		// scaffold_1      5       .       A       T       56.41   .       AC=3;AF=0.188;AN=16;BaseQRankSum=-1.804;DP=72;Dels=0.00;FS=0.000;HaplotypeScore=2.4206;MLEAC=2;MLEAF=0.125;MQ=26.05;MQ0=11;MQRankSum=0.618;QD=2.69;ReadPosRankSum=1.053 GT:AD:DP:GQ:PL  0|0:14,0:14:30:0,30,268 0|0:1,0:1:3:0,3,29      0/1:13,5:18:72:72,0,247 0|0:14,0:14:39:0,39,316 0|0:13,0:13:39:0,39,328 0|0:4,0:4:3:0,3,25      1|1:2,1:3:3:23,3,0      0|0:5,0:5:3:0,3,25


		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(fileName));
			int count = 0;
			int DatabasePointer = 0;
			int progress = 1000;
			String info = ER.readLine();
			while(ER.more() && ER.lookAhead() == '#'){
				info = ER.readLine();
			}
			String [] infoArray = info.split("\t");

			ArrayList<String> newSamples = new ArrayList<String>();
			for(int i = 9 ; i < infoArray.length;i++){
				while(Functions.contains(samples, sampleName)) sampleName = sampleName+"_1";
				samples.add(sampleName);
				newSamples.add(sampleName);
				System.out.println(sampleName);
			}
			infoArray = ER.readLine().split("\t");
			while(ER.more()){

				System.out.println("Now gathering vcf data from :"+this.Databases2.get(DatabasePointer).getName());
				infoArray = this.Databases2.get(DatabasePointer).addVCFSamples(ER, newSamples,infoArray);
				DatabasePointer++;
				//				}
				count++;
				if(count > progress ){
					progress += 10000;
					System.out.print(".");
				}
			}
			System.out.println("Finished");
		}catch(Exception E){E.printStackTrace();}
	}


	public int findPointer(String name){
		int DatabasePointer = 0;
		while(DatabasePointer < this.Databases2.size() && this.Databases2.get(DatabasePointer).getName().compareTo(name) != 0){
			DatabasePointer++;
		}
		if(DatabasePointer == this.Databases2.size()){
			System.out.println("Could not find "+name);
			return -1;
		}
		return DatabasePointer;



	} 

	public void loadFilterBedFile(String fileName){
		//	0			1		2     3		  4			5		6		7		8		9			 
		// #CHROM(0)  Start     Stop     Name     .     .    INFO  type    Something  XTR 	AInfo   
		// scaffold_1      767     2124    PAC:20891551.exon.3     .       -       phytozome8_0    exon    .       ID=PAC:20891551.exon.3;Parent=PAC:20891551;pacid=20891551
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(fileName));
			String currentDatabase = "";
			int count = 0;
			int DatabasePointer = 0;
			int progress = 1000;
			while(ER.more() && ER.lookAhead() == '#'){
				ER.readLine();
			}

			System.out.print("loading info from "+fileName+"................");

			while(ER.more()){
				String[] LineInfo = ER.readLine().split("\t");
				//				if(currentDatabase.compareTo(LineInfo[0])!=0){
				//					System.out.println("Finished");
				//					currentDatabase = LineInfo[0];
				//
				//					DatabasePointer = findPointer(currentDatabase);
				//					System.out.println(currentDatabase+"\t=  "+ this.Databases2.get(DatabasePointer).getName());
				//					System.out.print("Now gathering filter data for :"+this.Databases2.get(DatabasePointer).getName());						
				//				}
				//				//				if(DatabasePointer > this.Databases2.size())System.out.println(LineInfo[0]+"\t"+LineInfo[1]);
				//				//				else{
				if(this.Databases2.containsKey(LineInfo[0]))
					this.Databases2.get(LineInfo[0]).addBEDfilterInfo(LineInfo);
				else{
					System.out.print("No data for  :"+LineInfo[0]);
				}
				//				}
				count++;
				if(count > progress ){
					progress += 10000;
					System.out.print(".");
				}
			}
			//
			//			if(currentDatabase.compareTo(LineInfo[0])!=0){
			//				System.out.println("Finished");
			//				currentDatabase = LineInfo[0];
			//				DatabasePointer = findPointer(currentDatabase);
			//				System.out.print("Now gathering filter data for :"+this.Databases2.get(DatabasePointer).getName());						
			//			}
			//
			//			this.Databases2.get(DatabasePointer).addBEDfilterInfo(LineInfo);
			System.out.println("Finished");
		}catch(Exception E){E.printStackTrace();}
	}




	private void addCoverage(String[] columns){
		int databaseNr = findDatabase(columns[0]);
		if(databaseNr > -1)this.Databases2.get(databaseNr).addCoverage(Integer.parseInt(columns[3]), Integer.parseInt(columns[4]), columns[8]);
		else System.out.println(columns[0] +" not found in Database");

	}
	private int findDatabase(String name){
		int pointer = -1;
		for(int i = 0; i < this.Databases2.size(); i++){
			if(Databases2.get(this.DatabaseOrder[i]).isDatabase(name)){
				pointer = i;
				i = this.Databases2.size();
			}
		}
		return pointer;
	}

	private void printCoverage(){
		int pointer = -1;
		for(int i = 0; i < this.Databases2.size(); i++){
			Databases2.get(this.DatabaseOrder[i]).printCoverage();
		}
	}


	public void addPhase(boolean readPhased){	
		Enumeration<String> enumKey = Databases2.keys();
		while(enumKey.hasMoreElements()) {
			String key = enumKey.nextElement();
			this.Databases2.get(key).addPhase(readPhased);
		}
	}


	private void printpreMRNAs(String dir, String outFile){
		try{
			ExtendedWriter EW = ExtendedWriter.getFileWriter(dir+"/"+outFile);
			Enumeration<String> enumKey = Databases2.keys();
			while(enumKey.hasMoreElements()) {
				String key = enumKey.nextElement();
				this.Databases2.get(key).printpreMRNAs(EW);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}



	public void printPeronsalmRNAs(String dir, String outFile, String[] Samples){
		try{
			ExtendedWriter[] EWs =  new ExtendedWriter[Samples.length];
			for(int i = 0 ; i < Samples.length;i++){
				EWs[i] = ExtendedWriter.getFileWriter(dir+"/"+outFile+"_"+Samples[i]);
			}
			Enumeration<String> enumKey = Databases2.keys();
			while(enumKey.hasMoreElements()) {
				String key = enumKey.nextElement();
				this.Databases2.get(key).printPeronsalmRNAs(EWs,Samples);
			}
			for(int i = 0 ; i < Samples.length;i++){
				EWs[i].flush();
				EWs[i].close();
			}
		}catch(Exception E){E.printStackTrace();}
	}

	public void printPeronsalmRNAsInfo(String dir, String outFile){
		try{
			ExtendedWriter info = ExtendedWriter.getFileWriter(dir+"/"+outFile+".info");
			Enumeration<String> enumKey = Databases2.keys();
			while(enumKey.hasMoreElements()) {
				String key = enumKey.nextElement();
				this.Databases2.get(key).printPersonalmRNAsInfo(samples,info);
			}
			info.flush();
			info.close();
			System.out.println();
			for(int i = 0; i < samples.size();i++){
				enumKey = null;
				info = ExtendedWriter.getFileWriter(dir+"/"+outFile+"."+samples.get(i)+".info");
				enumKey = Databases2.keys();
				while(enumKey.hasMoreElements()) {
					String key = enumKey.nextElement();
					this.Databases2.get(key).printPersonalmRNAsInfo(samples.get(i),info);
				}
				info.flush();
				info.close();
			}
		}catch(Exception E){E.printStackTrace();}

	}	


	public void printVCFinfo(ExtendedWriter EW){
		try{
			for(int j = 0 ; j < this.samples.size();j++){
				for(int i = 0; i< this.Databases2.size();i++){
					this.Databases2.get(this.DatabaseOrder[i]).printVCFinfo(EW,this.samples.get(j));
				}
			}
		}catch(Exception E){E.printStackTrace();}

	}	
	
	
	public void printHaploGenome(ExtendedWriter EW, String sample, boolean father,ExtendedWriter SNPinfo ){
		try{
			for(int i = 0; i< this.Databases2.size();i++){
				this.Databases2.get(this.DatabaseOrder[i]).printHaploContig(EW,sample,father,SNPinfo);
				EW.println();
			}
		}catch(Exception E){E.printStackTrace();}

	}	

	public void printHomozygousGenome(ExtendedWriter EW, String sample){
		try{
			for(int i = 0; i< this.Databases2.size();i++){
				this.Databases2.get(this.DatabaseOrder[i]).printNNContig(EW,sample);
				EW.println();
			}
		}catch(Exception E){E.printStackTrace();}

	}	



	private void printmRNAs(String dir, String outFile){
		try{
			ExtendedWriter EW = ExtendedWriter.getFileWriter(dir+"/"+outFile);
			Enumeration<String> enumKey = Databases2.keys();
			while(enumKey.hasMoreElements()) {
				String key = enumKey.nextElement();
				this.Databases2.get(key).printmRNAs(EW);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}



	private void initiateAnnotation(){
		Enumeration<String> enumKey = Databases2.keys();
		while(enumKey.hasMoreElements()) {
			String key = enumKey.nextElement();
			Databases2.get(key).initiateAnnotation();
		}

	}

}




