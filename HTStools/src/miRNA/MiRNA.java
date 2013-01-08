package miRNA;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import alignment.mirfold;

import general.ExtendedWriter;
import general.Functions;
import Sequence.FastaSequences;

public class MiRNA {

	public static void main(String []args){
		fixLength(200,args[0], args[1], args[2],args[3]);
		parsemirfoldData(args[0],args[3]+".mirfold");
	}
	
	public static void run(Hashtable<String,String> T){
		String dir = Functions.getValue(T, "-d", ".");
		String matureFile = Functions.getValue(T, "-m", "-m not found");
		String hairpinFile = Functions.getValue(T, "-h", "-h not found");
		String outFile = Functions.getValue(T, "-o", "-o not found");
		int length = Integer.parseInt(Functions.getValue(T, "-l", "0"));
		if(length > 0){
			fixLength(length, dir, matureFile,hairpinFile, outFile);
		}
		parsemirfoldData(dir,outFile+".mirfold");
		
		
	}
	public static void fixLength(int length, String dir, String matureFile, String hairpinFile,String outFile){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(dir+"/"+outFile));
			FastaSequences matures = new FastaSequences(dir,matureFile);
			FastaSequences hairpins = new FastaSequences(dir,hairpinFile);

			for(int i = 0; i < hairpins.size(); i++){
				String[] parts = hairpins.get(i).getName().split("_");
				System.out.println(parts[0]);
				String[] Name = parts[0].split("-");
				String mir = "";
				for(int j = 0; j < Math.min(Name.length,3);j++){
					mir += Name[j];
				}
				boolean found = false;
				int pointer = 0;
				while (!found && pointer < matures.size()){
					parts = matures.get(pointer).getName().split(" ");
					Name = parts[0].split("-");
					String mature = "";
					for(int j = 0; j < Math.min(Name.length,3);j++){
						mature += Name[j];
					}
					if(mir.toLowerCase().compareTo(mature.toLowerCase())== 0){
						if(hairpins.get(i).FixSurroundingSequenceLength(200, matures.get(pointer).getSequence()))
							hairpins.get(i).printFasta(EW);
						found = true;
					}
					pointer	++;
				} 
			}
			EW.flush();
			EW.close();
			
		}catch(Exception E){E.printStackTrace();}
	}
	
	
	public static void parsemirfoldData(String dir, String mirfoldFile){
		
		ArrayList<mirfold> miRNAs = mirfold.parseFileSparse(dir, mirfoldFile);
		mirfold.writeMirfoldHits(dir,mirfoldFile+".best",miRNAs);
		
	}
}
