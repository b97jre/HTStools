package alignment;
import general.Functions;
import general.RNAfunctions;

import java.io.FileReader;
import java.util.Hashtable;
import general.ExtendedReader;

public class Alignments {

	private String[] GenomeNames;
	private Alignment alignments;	

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
		}
		Hashtable<String,String> T = Functions.parseCommandLine(args);
		if(T.containsKey("-n") && T.containsKey("-f")) {
			
		String dir = T.get("-d"); 	
		//String NCnumber = (String)T.get("-n");
		String FileName = T.get("-f");
		System.out.println("kommer jag hit?");
		Alignments A = new Alignments(dir,FileName);
		System.out.println("kommer jag hit?");
		char[][] Align = A.getAlignment(2812823,2812823 ,2812897,"Coli");
		System.out.println("kommer jag hit igen??");
		for(int i = 0 ; i < Align.length; i++){
			System.out.println(new String(Align[i]));
		}

		
		Align = A.getAlignment(2812823,2812897 ,2812823 ,"Coli");
		for(int i = 0 ; i < Align.length; i++){
			System.out.println(new String(Align[i]));
		}


		}

	}
	
	
	Alignments(String dir,String FileName){
		getSMLalignments(FileName,dir);
		
		
	}
	
	public static String[] getGenomes(String FileName,String Dir){
		String[] GenomeNames = null;
		try{
			
			ExtendedReader PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName));
			while(PTTFileReader.lookAhead() != '>'){
				
				String Line = PTTFileReader.readLine();
				if(Line.indexOf("#Sequence") != -1 && Line.indexOf("File") != -1){
					//System.out.println(Line);
					String tempName = Line.substring(Line.lastIndexOf("/")+1,Line.lastIndexOf("."));
					GenomeNames = Functions.addString(GenomeNames,tempName);
				}
			}
			
			
			for(int i = 0; i < GenomeNames.length;i++){
				System.out.println(GenomeNames[i]);
			}
			PTTFileReader.close();
		}
		catch(Exception E){E.printStackTrace();}

			return GenomeNames;
			
	
	}
	
	private void getSMLalignments(String FileName,String Dir){
		
		try{
			
			ExtendedReader PTTFileReader = new ExtendedReader(new FileReader(Dir+"/"+FileName));
			while(PTTFileReader.lookAhead() != '>'){
				
				String Line = PTTFileReader.readLine();
				if(Line.indexOf("#Sequence") != -1 && Line.indexOf("File") != -1){
					//System.out.println(Line);
					String tempName = Line.substring(Line.lastIndexOf("/")+1,Line.lastIndexOf("."));
					GenomeNames = Functions.addString(GenomeNames,tempName);
				}
			}
			
			for(int i = 0; i < GenomeNames.length;i++){
				System.out.println(GenomeNames[i]);
			}
			
			int[] start = new int[GenomeNames.length];
			int[] stop = new int[GenomeNames.length];
			for(int i = 0; i < GenomeNames.length;i++)
				start[i] = stop[i] = -1;
			char[][] alignedSequences = new char[GenomeNames.length][0];
			while(PTTFileReader.more()){
				if(PTTFileReader.lookAhead() == '='){
					for(int i = 0; i < GenomeNames.length;i++){
						if(start[i] ==  -1)
							alignedSequences[i] = null;	
					}
										
					addAlignment(start,stop,alignedSequences);
					start = new int[GenomeNames.length];
					stop = new int[GenomeNames.length];
					for(int i = 0; i < GenomeNames.length;i++)
						start[i] = stop[i] = -1;
					alignedSequences = new char[GenomeNames.length][0];
					PTTFileReader.skipLine();
				}
				else if(PTTFileReader.lookAhead() != '>'){
					PTTFileReader.skipLine();
				}
				else{
					String Info = PTTFileReader.readLine();
					//System.out.println(Info);
					//> 1:232-33953 + E:\align\microAlignment\EColik12.gbk
				
					int Genome = Integer.parseInt(Info.substring(2,Info.indexOf(':')));
					start[Genome-1] = Integer.parseInt(Info.substring(Info.indexOf(':')+1,Info.indexOf('-')));
					stop[Genome-1] = Integer.parseInt(Info.substring(Info.indexOf('-')+1,Info.indexOf(' ',Info.indexOf('-'))));
					String direction = Info.substring(Info.indexOf(' ',Info.indexOf('-'))+1, Info.indexOf(' ',Info.indexOf('-'))+2);
					if(direction.indexOf("+") == -1){
						int temp = start[Genome-1];
						start[Genome-1] = stop[Genome-1];
						stop[Genome-1] = temp;
					}
						
					//System.out.println(Genome+" "+start[Genome-1]+" "+stop[Genome-1]);
									
					int nrOfRows = 0;
					String [] Alignment = new String[300000];
					while(PTTFileReader.lookAhead() != '>' && PTTFileReader.lookAhead() != '='){
						//PTTFileReader.skipLine();
						
						Alignment[nrOfRows] = PTTFileReader.readLine();
						//Alignment[nrOfRows] = Alignment[nrOfRows].toUpperCase();
						nrOfRows++;
					}
					int AlignmentLength = 0;
					if(nrOfRows > 1)
						AlignmentLength = (nrOfRows-1) * Alignment[0].length() + Alignment[nrOfRows-1].length(); 
					else
						AlignmentLength = Alignment[0].length();
					char[] Alignment2 = new char[AlignmentLength];
					for(int i = 0; i < nrOfRows;i++){
						char[] temp = Alignment[i].toCharArray();
						for(int j = 0; j < temp.length; j++){
							Alignment2[i*80+j] = RNAfunctions.DNA2RNA(temp[j]);
						}
					}
					
					if(start[Genome-1] - stop[Genome-1] != 0)
						alignedSequences[Genome-1] = Alignment2;
					
					
				}
			}
			PTTFileReader.close();
			
		}
		catch(Exception E){E.printStackTrace();}
	}
	
	private void addAlignment(int[] start,int[] stop, char[][] charAlignments){
		Alignment newAlignment = new Alignment(start,stop,charAlignments);
		
		if(alignments == null)
			alignments = newAlignment;
		else
			alignments.addAlignment(newAlignment);
	}
	
	public char[][] getAlignment(int location, int upstream, int downstream, String focalGenome){
		int Genome = 0;
		while(Genome < GenomeNames.length && focalGenome.compareTo(GenomeNames[Genome]) != 0)
			Genome++;
		if(Genome == GenomeNames.length){
			System.out.println(focalGenome+" does not exist");
			return null;
			}
		else
			return alignments.getAlignment(location,upstream,downstream,Genome);
	}
	
	/**
	 * Returns the value of GenomeNames.
	 */
	public String[] getGenomeNames()
	{
		return GenomeNames;
	}

	/**
	 * Sets the value of GenomeNames.
	 * @param GenomeNames The value to assign GenomeNames.
	 */
	public void setGenomeNames(String[] GenomeNames)
	{
		this.GenomeNames = GenomeNames;
	}
	
	

	
	
}
