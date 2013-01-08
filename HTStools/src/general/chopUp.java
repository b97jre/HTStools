package general;

import java.io.FileReader;
import java.io.FileWriter;
import general.ExtendedReader;
import general.ExtendedWriter;




public class chopUp{

	public static void chopUpSequence(String fileName,String outfile){
		try{
			
			ExtendedReader PTTFileReader = new ExtendedReader(new FileReader(fileName));
			ExtendedWriter newFile = new ExtendedWriter(new FileWriter(outfile));
			PTTFileReader.more();
			PTTFileReader.skipLine();
//			String Name = Line.substring(1);
			PTTFileReader.more();
			int count = 0;
			int nr = 1;
			String Sequence = "";
			while(PTTFileReader.more()){
				if(count != 4){
					Sequence +=  PTTFileReader.readLine();
					count++;
				}
				else{
					if(Math.random() < 0.33){
						newFile.println(">"+outfile+"_"+nr);
						newFile.println(Sequence);
					}
					nr++;
					count = 0;
					Sequence = "";
				}
				if(nr > 10000){
					PTTFileReader.close();
					newFile.flush();
					newFile.close();
					break;
				}
			}
		
		}
		catch(Exception E){E.printStackTrace();}
		
		
	}
	
	public static void main(String []args){
//		chopUp.chopUpSequence(args[0],args[1]);
		chopUp.parseMirFile(args[0]);
//		chopUp.removeN(args[0],args[1]);
	}

	
	public static void parseMirFile(String fileName){
		try{
			//generated_seq_freqs_generated_seq_6_dwn l=91nt penalty=6 dG=-20.97 (-0.230/nt) miR@20-41
			int nrOfstructuresNotFound = 0;
			ExtendedReader PTTFileReader = new ExtendedReader(new FileReader(fileName));
			while(PTTFileReader.more()){
				String Sequence =  PTTFileReader.readLine();
				if(Sequence.indexOf("no structure found")>0){
					nrOfstructuresNotFound++;
				}
				else if(Sequence.indexOf("l=") > 0){
					//System.out.println(Sequence);	
					String L = Sequence.substring(Sequence.indexOf("l=")+2,Sequence.indexOf("nt"));
					String p = Sequence.substring(Sequence.indexOf("y=")+2,Sequence.indexOf(" dG"));
					String DG = Sequence.substring(Sequence.indexOf("dG=")+3,Sequence.indexOf(" ("));
						String Name = Sequence.substring(0,Sequence.indexOf(" "));
					System.out.println(Name+"\t"+L+"\t"+p+"\t"+DG);	
				}
			}
		
		}
		catch(Exception E){E.printStackTrace();}
		
		
	}
	

	public static void removeN(String fileName,String outfile){
		try{
			
			ExtendedReader PTTFileReader = new ExtendedReader(new FileReader(fileName));
			ExtendedWriter newFile = new ExtendedWriter(new FileWriter(outfile));
			PTTFileReader.more();
			String Line = PTTFileReader.readLine();
			newFile.println(Line);
			int count = 0;
			while(PTTFileReader.more()){
				char nt = (char)PTTFileReader.readChar();
				
				if(nt != 'N'){
					newFile.print(nt);
					count++;
				}
				if(count == 70){
					count = 0;
					newFile.println();
				
				}
					
			}
			newFile.flush();
			newFile.close();
		
		}
		catch(Exception E){E.printStackTrace();}
		
		
	}


	
	//parse mirfoldFile

	chopUp(){}
	
	
	
	
}





