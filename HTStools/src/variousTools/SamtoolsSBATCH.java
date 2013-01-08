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
	String F;
	String f;
	
	boolean merge;

	public SamtoolsSBATCH(){
		inDir = projectDir= suffix = time = null;
		merge = false;
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
		if(T.containsKey("-suffix"))
			suffix = Functions.getValue(T, "-suffix", "sam");
		else{
			System.out.println("must contain a suffix e.g. -suffix sam/bam now set to sam");
		}
		F = Functions.getValue(T, "-F", "-1");
		f = Functions.getValue(T, "-f", "-1");
		
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
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_samtools.sh"));
			if(program.compareTo("sort") == 0)
				samtoolsSortReads(EW,sbatch,projectDir+"/"+inDir, projectDir+"/"+outDir,timeStamp,suffix,0);
			else
				samtoolsDir(EW,sbatch,projectDir+"/"+inDir,projectDir+"/"+outDir, timeStamp,suffix,0);
			EW.flush();
			EW.close();
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

				String 	sbatchFileName = finalInDir+"/scripts/"+timestamp+"_"+count+"_"+lastDir+".samtools.sbatch";

				generalSbatchScript.println("sbatch "+ sbatchFileName);

				ExtendedWriter EW = new ExtendedWriter(new FileWriter(sbatchFileName));
				sbatch.printSBATCHinfo(EW,finalInDir,timestamp,count,"samtools"+"_"+lastDir, time);


				EW.println("module load bioinfo-tools");
				EW.println();
				
				EW.println("module load samtools");
				EW.println();

				EW.println("cd "+finalInDir);
				EW.println();
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
						EW.print("samtools view -bS");
						if(Integer.parseInt(F) > 0)
							EW.print(" -F "+F);
						if(Integer.parseInt(f) > -1)
							EW.print(" -f "+f);
						EW.println(" -o "+finalOutDir+"/"+nameBase+"bam "+fileNames.get(i)+" &");

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
							EW.println("samtools sort "+finalOutDir+"/"+nameBase+"bam "+ finalOutDir+"/"+nameBase+"sorted &");
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
							EW.println("samtools index "+finalOutDir+"/"+nameBase+"sorted.bam &");
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
					EW.println("samtools sort "+finalOutDir+"/"+nameBase+"bam "+ finalOutDir+"/"+nameBase+"sorted &");
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
					EW.println("samtools index "+finalOutDir+"/"+nameBase+"sorted.bam &");
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
			samtoolsDir(generalSbatchScript,sbatch,finalInDir+"/"+subDirs.get(i), finalOutDir+"/"+subDirs.get(i), timestamp,suffix,count);
		}


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
					if(Integer.parseInt(F) > 0)
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




}
