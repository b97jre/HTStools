package alignment;

import general.ExtendedReader;
import general.Functions;
import general.IOTools;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Set;

import Sequence.CfastaSequences;
import Sequence.FastaSequence;
import Sequence.FastaSequences;
import Sequence.Solid;

import general.ExtendedWriter;




public class GFF3info extends Hashtable <String,Contig> implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//private String GenomeName;
	private String Name;
	private FastaSequences Contigs;


	public static void main(String[] args) {
		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}
		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);

		String dir = Functions.getValue(T, "-d",IOTools.getCurrentPath());
		String gff3File = Functions.getValue(T, "-gff3");
		String fastaFile = Functions.getValue(T,"-fasta");
		GFF3info test = new GFF3info(fastaFile,gff3File, dir);
		

		//GFF3info.run(T);

	}


	
	public GFF3info(String fastaFile , String gff3File, String WD){
		System.out.println("Parsing Sequence file....");
		this.Contigs = FastaSequences.getFastaSequences(fastaFile,WD);

		ArrayList<String> ContigNames = Contigs.getAllSequenceNames();

		System.out.println("Parsing GFF3 file....");
		for(int i = 0; i < ContigNames.size();i++){
			this.put(ContigNames.get(i), new Contig(ContigNames.get(i), Contigs.get(i)));
		}

		System.out.println("Parsing GFF3 file....");

		addGFF3info(WD, gff3File);
		printUpstream_5UTR_FirstIntron_Sequence(WD, gff3File+".upstream.1000.fa",1000);

	}


	public void addGFF3info(String dir, String fileName){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));
			while(ER.more()){
				if((char)ER.lookAhead() == '#'){
					String GFF3line = ER.readLine();
				}
				else{
					readInfo(ER);
				}
			}
		}catch(Exception E){E.printStackTrace();}
	}


	private void readInfo(ExtendedReader ER){
		String GFF3line = ER.readLine();
		String[] columns = GFF3line.split("\t");
		if(this.containsKey(columns[0])){
			this.get(columns[0]).addInfo(columns);
		}else{
			System.out.println(GFF3line);
		}
	}



/*	private void printCodingSequence(String dir,String file){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(dir+"/"+file+".fa"));
			// Get a set of all the entries (key - value pairs) contained in the Hashtable
			Set<Entry<String, Contig>> entrySet = this.entrySet();
			// Obtain an Iterator for the entries Set
			Iterator<Entry<String, Contig>> it = entrySet.iterator();
			// Iterate through Hashtable entries
			while(it.hasNext()){
				this.get(it.next()).printCodingSequence(EW);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}
*/

	private void printUpstream_5UTR_FirstIntron_Sequence(String dir,String file, int upstremLength){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(dir+"/"+file+".fa"));
			// Get a set of all the entries (key - value pairs) contained in the Hashtable
			for(Enumeration<String> e = this.keys();
					e.hasMoreElements();
					this.get(e.nextElement()).printUpstream_5UTR_FirstIntron_Sequence(EW, upstremLength)
				);
			EW.flush();
			EW.close();
			
		}
		catch(Exception E){E.printStackTrace();}
	}


	public static String getExtraInfo(String Column9Info, String info){
		String [] extra = Column9Info.split(";");
		String outInfo = null;
		if(extra != null){
			for(int i = 0; i < extra.length; i++){
				if(extra[i].indexOf(info+"=") == 0){
					String[] IDs = extra[i].split("=");
					outInfo = IDs[1];
				}
			}
		}
		return outInfo;

	}




}




