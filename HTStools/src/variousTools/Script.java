package variousTools;

import general.ExtendedReader;
import general.ExtendedWriter;
import general.Functions;
import general.IOTools;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Hashtable;



public class Script {

	String inFile;
	String projectDir;
	String time;
	boolean core;
	String parameters;

	public Script(){
		inFile = projectDir= time = null;
		core = true;
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

		
		if(T.containsKey("-node"))core = false;
		
		if(T.containsKey("-i"))
			inFile= Functions.getValue(T, "-i");
		else{
			System.out.println("must contain inDirectory -i");
			allPresent = false;
		}
		projectDir= Functions.getValue(T, "-d", IOTools.getCurrentPath());
		parameters= Functions.getValue(T, "-parameters", "");
		
		if(T.containsKey("-time"))
			time = Functions.getValue(T, "-time");
		else{
			System.out.println("time is not set. Now set to default 15:00");
			time = Functions.getValue(T, "-time", "15:00");
		}
		if(allPresent)
			shellScript(sbatch, timeStamp);
		else
			System.out.println("\n\nAborting run because of missing arguments for script.");
	}


	public void shellScript(SBATCHinfo sbatch, String timeStamp){
		try{
			if(!IOTools.isDir(projectDir+"/scripts"))
				IOTools.mkDir(projectDir+"/scripts");
			if(!IOTools.isDir(projectDir+"/reports"))
				IOTools.mkDir(projectDir+"/reports");
			String info1 = null;
			if(inFile.indexOf("/")>-1){
				String[] info = inFile.split("/");
				info1 = info[info.length-1];
			}else{
				info1= inFile;
			}
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(projectDir+"/scripts/"+timeStamp+"_"+info1+".sbatch"));
			shellScript(EW,sbatch,timeStamp);
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	} 


	public void shellScript(ExtendedWriter EW, SBATCHinfo sbatch, String timestamp){
		String info1 = null;
		if(inFile.indexOf("/")>-1){
			String[] info = inFile.split("/");
			info1 = info[info.length-1];
		}else{
			info1= inFile;
		}
		if(!core)sbatch.printSBATCHinfo(EW,projectDir,timestamp,0,info1, time);
		else sbatch.printSBATCHinfoCore(EW,projectDir,timestamp,0,info1, time);
		
		EW.println("cd "+this.projectDir);
		EW.println();
		EW.println();
		EW.println("sh "+inFile +" "+parameters);
		EW.println("echo \" Shellscript can be found here: "+projectDir+"/reports/"+0+"_inFile_"+timestamp+".shellScript\"");
		EW.println("echo \" Parameters added were: "+parameters+"\"");
		try{
			IOTools.copy(inFile,projectDir+"/reports/"+0+"_inFile_"+timestamp+".shellScript");
		}catch(Exception E){E.printStackTrace();}
	}

	

}
