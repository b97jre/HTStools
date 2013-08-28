package variousTools;

import general.ExtendedWriter;
import general.Functions;
import general.IOTools;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Hashtable;

public class SamtoolsSBATCH {

	String inDir;
	String projectDir;
	String outDir;
	String suffix;
	String time;
	int F;
	int f;

	boolean merge;
	boolean sort;
	boolean header;
	boolean view;
	boolean index;
	boolean flagstat;

	public SamtoolsSBATCH(){
		inDir = projectDir= suffix = time = null;
		merge = sort = header = view = index = flagstat = false;
	}

	public static void main(String []args){

		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}
		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);
		FilterFastqSBATCH filter = new FilterFastqSBATCH();
		filter.run(T);
	}

	public void run(Hashtable<String,String> T){

		boolean allPresent = true;

		String timeStamp = Functions.getDateTime();
		SBATCHinfo sbatch = new SBATCHinfo();
		if(!sbatch.addSBATCHinfo(T)) allPresent = false;

		if(T.containsKey("-i"))
			inDir= Functions.getValue(T, "-i", ".");
		else{
			System.out.println("must contain inDirectory -i");
			allPresent = false;
		}

		outDir= Functions.getValue(T, "-o", inDir);

		projectDir= Functions.getValue(T, "-pDir", IOTools.getCurrentPath());
		if(!T.containsKey("-time"))
			System.out.println("must contain time e.g. -time 01:00:00 . Now set to default 30 minutes");
		time = Functions.getValue(T, "-time", "30:00");
		if(IOTools.isDir(inDir))
			suffix = Functions.getValue(T, "-suffix", "sam");
		else
			suffix = inDir.substring(inDir.lastIndexOf(".")+1);
		if(!T.containsKey("-suffix"))
			System.out.println("must contain a suffix e.g. -suffix sam/bam now set to "+suffix);
		F = Functions.getInt(T, "-F", 0);
		f = Functions.getInt(T, "-f", 0);

		if(T.containsKey("-view")) this.view = true;
		if(T.containsKey("-h")) this.header = true;
		if(T.containsKey("-sort")) this.sort = true;
		if(T.containsKey("-flagstat")) this.flagstat = true;
		if(T.containsKey("-index")) this.index = true;
		if(T.containsKey("-All")) {sort = header = view = index = flagstat = true;}



		String program = Functions.getValue(T, "-samtools", "all");

		if(T.containsKey("-merge")) this.merge = true;


		if(allPresent)
			samtoolsTop(sbatch, timeStamp,suffix, program);
		else
			System.out.println("\n\nAborting run because of missing arguments for samtools.");
	}


	public void samtoolsTop(SBATCHinfo sbatch, String timeStamp,String suffix,String program){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			if(IOTools.fileExists(projectDir+"/"+inDir) && !IOTools.isDir(projectDir+"/"+inDir)){
				samtoolsFile(sbatch,projectDir+"/"+inDir,projectDir, timeStamp,suffix,0);
			}else if(IOTools.isDir(projectDir+"/"+inDir)){
				ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_samtools.sh"));
				samtoolsDir(EW,sbatch,projectDir+"/"+inDir,projectDir+"/"+outDir, timeStamp,suffix,0);
				EW.flush();
				EW.close();
			}else
				System.out.println(projectDir+"/"+inDir+" was not found");
		}catch(Exception E){E.printStackTrace();}
	} 


	public void samtoolsDir(ExtendedWriter generalSbatchScript, SBATCHinfo sbatch , String finalInDir, String finalOutDir, String timestamp,String suffix, int count){

		if(!IOTools.isDir(finalOutDir))
			IOTools.mkDir(finalOutDir);


		ArrayList <String> fileNames = IOTools.getSequenceFiles(finalInDir,suffix);

		if(fileNames.isEmpty()){

			System.out.println("No "+suffix+" files in folder :"+finalInDir);

		}
		else{
			if(!IOTools.isDir(finalInDir+"/reports"))
				IOTools.mkDir(finalInDir+"/reports");
			if(!IOTools.isDir(finalInDir+"/scripts"))
				IOTools.mkDir(finalInDir+"/scripts");

			try{
				String[] dirs = finalInDir.split("/");
				String lastDir = dirs[dirs.length-1];

				int start = 0;
				//convert sam to bam
				for(int i = 0; i < fileNames.size(); i++){	
					System.out.println(fileNames.get(i) +" found");
					String 	sbatchFileName = finalInDir+"/scripts/"+timestamp+"_"+i+"_"+lastDir+".samtools.sbatch";

					generalSbatchScript.println("sbatch "+ sbatchFileName);

					ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
					sbatch.printSBATCHinfoCore(EW,finalInDir,timestamp,count,"samtools"+"_"+lastDir, time);


					EW.println("module load bioinfo-tools");
					EW.println();

					EW.println("module load samtools");
					EW.println();

					EW.println("cd "+finalInDir);
					EW.println();


					String nameBase = fileNames.get(i).substring(0,fileNames.get(i).indexOf(suffix));
					if(suffix.indexOf("sam")>-1)
						sam2bam(EW,fileNames.get(i),f, F, sort, header,view,flagstat,index);
					else
						bam2bam(EW,fileNames.get(i),f, F, sort, header,view,flagstat,index);


					EW.flush();
					EW.close();
				}
			}catch(Exception E){E.printStackTrace();}

		}
		ArrayList <String> subDirs = IOTools.getDirectories(finalInDir);
		for(int i = 0; i < subDirs.size(); i++){
			samtoolsDir(generalSbatchScript,sbatch,finalInDir+"/"+subDirs.get(i), finalOutDir+"/"+subDirs.get(i), timestamp,suffix,count);
		}
	}

	public void samtoolsFile(SBATCHinfo sbatch , String fileName, String finalOutDir, String timestamp,String suffix, int count){

		if(!IOTools.isDir(finalOutDir+"/reports"))
			IOTools.mkDir(finalOutDir+"/reports");
		if(!IOTools.isDir(finalOutDir+"/scripts"))
			IOTools.mkDir(finalOutDir+"/scripts");

		try{
			String[] dirs = finalOutDir.split("/");
			String lastDir = dirs[dirs.length-1];

			int start = 0;
			//convert sam to bam
			String 	sbatchFileName = finalOutDir+"/scripts/"+timestamp+"_"+lastDir+".samtools.sbatch";


			ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
			sbatch.printSBATCHinfoCore(EW,finalOutDir,timestamp,count,"samtools"+"_"+lastDir, time);


			EW.println("cd "+finalOutDir);
			EW.println();

			String nameBase = fileName.substring(0,fileName.indexOf(suffix));
			if(suffix.indexOf("sam")>0)
				sam2bam(EW,fileName,f, F, sort, header,view,flagstat,index );
			else
				bam2bam(EW,fileName,f, F, sort, header,view,flagstat,index);


			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}

	}



	public static String sam2bam(ExtendedWriter EW, String samFile, int f, int F, boolean sort, boolean header, boolean view, boolean flagstat, boolean index){
		String suffix = ".sam";
		String nameBase = samFile;
		if(samFile.endsWith(".sam"))
			nameBase = samFile.substring(0,samFile.lastIndexOf(suffix));
		EW.println();
		EW.println("module load bioinfo-tools");
		EW.println("module load samtools");
		EW.println();
		if(view){
			if(!header)
				EW.print("samtools view -bS");
			else
				EW.print("samtools view -bSh");
			if(F > 0 || f > 0 ){
				if(F > 0)
					EW.print(" -F "+F);
				if(f > 0)
					EW.print(" -f "+f);
				EW.println(" -o "+nameBase+".F_"+F+".f_"+f+".bam "+samFile );
				nameBase = nameBase+".F_"+F+".f_"+f;
			}else{
				EW.println(" -o "+nameBase+".bam "+samFile );
			}
			if(sort){
				EW.println("samtools sort "+nameBase+".bam "+nameBase+".sorted ");
				nameBase = nameBase+".sorted";
			}
			if(index){
				EW.println("samtools index "+nameBase+".bam");
			}
			if(flagstat)
				EW.println("samtools flagstat "+nameBase+".bam >"+nameBase+".bam.flagstat");
			
			return nameBase+".bam";
			
		}
		else{
			System.out.println(samFile +" not converted. Add flag '-view' to convert to bam File");
		}
		return null;
	}

	public static void bam2bam(ExtendedWriter EW, String bamFile, int f, int F, boolean sort, boolean header, boolean view, boolean flagstat,boolean index){
		String nameBase = bamFile;
		if(bamFile.endsWith(".bam"))
			nameBase = bamFile.substring(0,bamFile.lastIndexOf(".bam"));
		EW.println();
		EW.println("module load bioinfo-tools");
		EW.println("module load samtools");
		EW.println();
		if(view && (F > 0 || f > 0 )){
			if(!header)
				EW.print("samtools view -b");
			else
				EW.print("samtools view -bh");
			if(F > 0)
				EW.print(" -F "+F);
			if(f > 0)
				EW.print(" -f "+f);
			EW.println(" -o "+nameBase+".F_"+F+".f_"+f+".bam "+bamFile );
			nameBase = nameBase+".F_"+F+".f_"+f;
		}
		if(sort){
			EW.println("samtools sort "+nameBase+".bam "+nameBase+".sorted ");
			nameBase = nameBase+".sorted";
		}
		if(index){
			EW.println("samtools index "+nameBase+".bam");
		}
		if(flagstat)
			EW.println("samtools flagstat "+nameBase+".bam >"+nameBase+".bam.flagstat");
	}





	public void samtoolsSortReads(ExtendedWriter generalSbatchScript, SBATCHinfo sbatch , String finalInDir ,String finalOutDir,String timestamp,String suffix, int count){
		// asuming that it is in bam format
		if(!IOTools.isDir(finalOutDir))
			IOTools.mkDir(finalOutDir);

		ArrayList <String> fileNames = IOTools.getSequenceFiles(finalInDir,suffix);
		if(fileNames.isEmpty()){
			System.out.println("No "+suffix+" files in folder :"+finalInDir);
		}
		else{
			if(!IOTools.isDir(finalInDir+"/reports"))
				IOTools.mkDir(finalInDir+"/reports");
			if(!IOTools.isDir(finalInDir+"/scripts"))
				IOTools.mkDir(finalInDir+"/scripts");

			try{
				String[] dirs = finalInDir.split("/");
				String lastDir = dirs[dirs.length-1];


				for(int i = 0; i < fileNames.size(); i++){	
					String nameBase = fileNames.get(i).substring(0,fileNames.get(i).indexOf(suffix));
					String 	sbatchFileName = finalInDir+"/scripts/"+timestamp+"_"+count+"_"+lastDir+".samtools.sbatch";
					generalSbatchScript.println("sbatch "+ sbatchFileName);
					ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
					sbatch.printSBATCHinfoCore(EW,finalInDir,timestamp,count,"samtools"+"_"+lastDir, time);

					EW.println("module load bioinfo-tools");
					EW.println();
					EW.println("module load samtools");
					EW.println();

					EW.println("cd "+finalInDir);
					EW.println();
					EW.println();
					EW.println("samtools sort  "+fileNames.get(i)+ " "+finalOutDir+"/"+nameBase+"sorted &");
					EW.println();
					EW.println();
					EW.println("wait");

					EW.flush();
					EW.close();
					count++;

				}

			}catch(Exception E){E.printStackTrace();}

		}
		ArrayList <String> subDirs = IOTools.getDirectories(finalInDir);
		for(int i = 0; i < subDirs.size(); i++){
			samtoolsSortReads(generalSbatchScript,sbatch,finalInDir+"/"+subDirs.get(i),finalOutDir+"/"+subDirs.get(i), timestamp,suffix,count);
		}


	}





	public void samtoolsSortName(ExtendedWriter generalSbatchScript, SBATCHinfo sbatch , String finalInDir, String finalOutDir, String timestamp,String suffix, int count){


		ArrayList <String> fileNames = IOTools.getSequenceFiles(finalInDir,suffix);


		if(fileNames.isEmpty()){
			System.out.println("No "+suffix+" files in folder :"+finalInDir);
		}
		else{
			if(!IOTools.isDir(finalInDir+"/reports"))
				IOTools.mkDir(finalInDir+"/reports");
			if(!IOTools.isDir(finalInDir+"/scripts"))
				IOTools.mkDir(finalInDir+"/scripts");

			try{
				String[] dirs = finalInDir.split("/");
				String lastDir = dirs[dirs.length-1];

				String 	sbatchFileName = finalInDir+"/scripts/"+timestamp+"_"+count+"_"+lastDir+".samtools.sbatch";

				generalSbatchScript.println("sbatch "+ sbatchFileName);

				ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
				sbatch.printSBATCHinfo(EW,finalInDir,timestamp,count,"samtools"+"_"+lastDir, time);


				EW.println("module load bioinfo-tools");
				EW.println("module load samtools");

				EW.println("cd "+finalInDir);
				EW.println();
				EW.println();

				int start = 0;
				//convert sam to bam
				for(int i = 0; i < fileNames.size(); i++){	
					String[] nameParts = fileNames.get(i).split("\\.");
					String nameBase = "";
					for(int k = 0; k < nameParts.length-1; k++){
						nameBase += nameParts[k]+".";
					}
					if(F > 0)
						EW.println("samtools view -bS -F "+F+" -o "+finalInDir+"/"+nameBase+"bam "+fileNames.get(i)+" &");
					else
						EW.println("samtools view -bS  -o "+finalInDir+"/"+nameBase+"bam "+fileNames.get(i)+" &");

					if((i+1)%8 == 0  &&  fileNames.size()-i != 1){


						EW.println();
						EW.println();
						EW.println("wait");

						for(int j = start; j < i+1; j++){
							nameParts = fileNames.get(j).split("\\.");
							nameBase = "";
							for(int k = 0; k < nameParts.length-1; k++){
								nameBase += nameParts[k]+".";
							}
							EW.println("samtools sort "+finalInDir+"/"+nameBase+"bam "+ finalInDir+"/"+nameBase+"sorted &");
						}


						EW.println();
						EW.println();
						EW.println("wait");

						for(int j = start; j < i+1; j++){
							nameParts = fileNames.get(j).split("\\.");
							nameBase = "";
							for(int k = 0; k < nameParts.length-1; k++){
								nameBase += nameParts[k]+".";
							}
							EW.println("samtools index "+finalInDir+"/"+nameBase+"sorted.bam &");
						}


						EW.println();
						EW.println();
						EW.println("wait");

						EW.flush();
						EW.close();
						count++;

						sbatchFileName = finalInDir+"/scripts/"+timestamp+"_"+count+"_"+lastDir+".samtools.sbatch";
						generalSbatchScript.println("sbatch "+ sbatchFileName);
						EW = new ExtendedWriter(new FileWriter(sbatchFileName));
						sbatch.printSBATCHinfo(EW,finalInDir,timestamp,count,"samtools"+"_"+lastDir, time);

						EW.println("module load bioinfo-tools");
						EW.println("module load samtools");

						EW.println("cd "+finalInDir);
						EW.println();
						EW.println();
						start = i+1;
					}
				}
				EW.println();
				EW.println();
				EW.println("wait");



				for(int j = start; j < fileNames.size(); j++){
					String[] nameParts = fileNames.get(j).split("\\.");
					String nameBase = "";
					for(int k = 0; k < nameParts.length-1; k++){
						nameBase += nameParts[k]+".";
					}
					EW.println("samtools sort "+finalInDir+"/"+nameBase+"bam "+ finalInDir+"/"+nameBase+"sorted &");
				}
				EW.println();
				EW.println();
				EW.println("wait");


				for(int j = start; j < fileNames.size(); j++){
					String[] nameParts = fileNames.get(j).split("\\.");
					String nameBase = "";
					for(int k = 0; k < nameParts.length-1; k++){
						nameBase += nameParts[k]+".";
					}
					EW.println("samtools index "+finalInDir+"/"+nameBase+"sorted.bam &");
				}


				EW.println();
				EW.println();
				EW.println("wait");



				EW.flush();
				EW.close();



			}catch(Exception E){E.printStackTrace();}

			count++;
		}
		ArrayList <String> subDirs = IOTools.getDirectories(finalInDir);
		for(int i = 0; i < subDirs.size(); i++){
			samtoolsSortName(generalSbatchScript,sbatch,finalInDir+"/"+subDirs.get(i), finalOutDir+"/"+subDirs.get(i), timestamp,suffix,count);
		}


	}


	public static void indexBamFile(String bamFile, ExtendedWriter EW, boolean loadModules){

		if(loadModules){
			EW.println("module load bioinfo-tools");
			EW.println();
			EW.println("module load samtools");
			EW.println();
		}

		EW.println("samtools index  "+bamFile);
		EW.println();

		EW.flush();


	}

	public static void flagstat(String bamFile, ExtendedWriter EW, boolean loadModules){

		if(loadModules){
			EW.println("module load bioinfo-tools");
			EW.println();
			EW.println("module load samtools");
			EW.println();
		}

		EW.println("samtools flagstat  "+bamFile);
		EW.println();

		EW.flush();


	}



	public static boolean indexFastaFile(String fastaFile, ExtendedWriter EW, boolean loadModules){

		if(loadModules){
			EW.println("module load bioinfo-tools");
			EW.println();
			EW.println("module load samtools");
			EW.println();
		}

		EW.println("samtools faidx  "+fastaFile);
		EW.println();

		EW.flush();
		return true;


	}




}
