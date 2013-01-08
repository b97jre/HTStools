package variousTools;

import general.ExtendedReader;
import general.ExtendedWriter;
import general.Functions;

import java.io.FileReader;
import java.io.FileWriter;
import java.util.Hashtable;

public class TabBlast2OrthoMCLdbo {


	public static void run(Hashtable<String,String> T){

		String inFile = Functions.getValue(T, "-i", ".");
		String outFile = Functions.getValue(T, "-o", ".");
		if(inFile.compareTo(".") == 0 || outFile.compareTo(".") == 0){
			System.out.println("blast2OrthoMCL must contain a in file (-i) and a outfile (-o)");
			return;
		}
		//int  A = checkTabularBlastFile(inFile);
		
		covertTabularBlastFile2OrthoMCLformat(inFile, outFile);
	}
	public static int checkTabularBlastFile(String TabulatedBlastFile){
		int A = 1;
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(TabulatedBlastFile));
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(TabulatedBlastFile+".wrong"));

			while(ER.more()){
				String Line = ER.readLine();
				String[] info = Line.split("\t");
				if(info.length != 12){
					EW.println(Line);
					A = -1;
				}else{
					int B = checkBasicInfo(info[2]);
					if(B == -1){
						System.out.println(Line);
						A=-1;
					}
				}
			}
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}
		return A;
	}


	public static int checkBasicInfo(String PercentMatch){
		try{
			PercentMatch.substring(0,PercentMatch.indexOf('.'));
		}
		catch(Exception E){
			return -1;
		}
		return 1;
	}

	public static void printBasicInfo(ExtendedWriter EW,int DBID, String queryName, String targetName, String Evalue, String PercentMatch){
		try{
			String match = PercentMatch.substring(0,PercentMatch.indexOf('.'));
			EW.print(DBID+";"+queryName+";0;"+targetName+";0;"+Evalue+";"+match+";");
		}
		catch(Exception E){
			System.out.println("something wrong with the combination queryName and targetName");
		}
	}
	public static void printHSPinfo(ExtendedWriter EW, int HSPID, String queryStart, String queryEnd, String targetStart, String targetEnd){
		EW.print(HSPID+":"+queryStart+"-"+queryEnd+":"+targetStart+"-"+targetEnd+".");
	}
	/*****************
		 convertTabularBlastFile2OrthoMCLformat	
		 Takes a blastTabularFile and converts it to the OrthoMCL format that is needed in for module 4
		Blast tabular format is 
		1	Query
		2 	Subject
		3	% id
		4	alignment length
		5	mistmatches
		6	gap openings
		7	q.start
		8	q.end
		9	s.start
		10	s.end
		11	e-value
		12	bit score

		orthoMCLformat is
		FORMAT of all.bpo or "usr_bpo_file"
		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		1;At1g01190;535;At1g01190;535;0.0;97;1:1-535:1-535.
		2;At1g01190;535;At1g01280;510;2e-56;29;1:69-499:28-474.
		3;At1g01190;535;At1g11600;510;1e-45;27;1:59-531:21-509.




		Each line represents each query-subject similarity relation. And all the info is
		separated by ";", which are, in order, similarity id, query id, query length, 
		subject id, subject length, BLAST E-value, percent identity, HSP info (each HSP
		is in the format of HSP_id:query_start-query_end:subject_start-subject_end. 
		different HSP info are seperated by "." )
		IMPORTANT: 1. Similarity ID represents BPO file line id, so it should start 
		              from 1 for the first line, and be consecutive for the whole file.
		           2. BPO file is a parsing result from BLAST, so for each query gene
		              id, its hits can't be scattered in the file, but should be listed 
		              in ajacent lines.
		           3. For BLAST m8 format (i.e. $BLAST_FORMAT="compact"), sequence length
		              information is not stored. So when running OrthoMCL in mode 3,
		              the corresponding columns (i.e. query length, and subject length)
		              will be 0. Please do not use Percent match cutoff in this case.

		Blast output format 
		LOC_Os12g44380.3	jgi|Araly1|484195|fgenesh2_kg.5__187__AT2G02860.1	49.27	205	96	4	116	314	348	550	6e-47	 187
		LOC_Os12g44380.3	jgi|Araly1|484195|fgenesh2_kg.5__187__AT2G02860.1	48.15	54	28	0	25	78	205	258	1e-08	60.5
		is converted to
		71817815;LOC_Os12g44380.3;0;jgi|Araly1|484195|fgenesh2_kg.5__187__AT2G02860.1;0;6e-47;50;1:116-314:348-550.2:25-78:205-258.


	 */
	
	public static void covertTabularBlastFile2OrthoMCLformat(String TabulatedBlastFile, String usr_bpo_file){
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(usr_bpo_file));
			ExtendedReader ER = new ExtendedReader(new FileReader(TabulatedBlastFile));

			String targetName = "";
			String queryName = "";
			int  ID =1;
			int HSPID=0;
			String Line = ER.readLine();
			String[] info = Line.split("\t");
			if(info.length == 12){
				queryName = info[0];
				targetName = info[1];
				printBasicInfo(EW,ID,info[0],info[1],info[10],info[2]);
				ID++;
				HSPID = 1;
				printHSPinfo(EW,HSPID,info[6],info[7], info[8], info[9]);
				HSPID++;
			}
			else{
				queryName = info[0];
				targetName = info[2];
				printBasicInfo(EW,ID,info[0],info[2],info[11],info[3]);
				ID++;
				HSPID = 1;
				printHSPinfo(EW,HSPID,info[7],info[8], info[9], info[10]);
				HSPID++;
			}
			while(ER.more()){
				Line = ER.readLine();
				info = Line.split("\t");
				if(info.length == 12){
					if(info[0].compareTo(queryName) != 0){ // if new query print names and info
						queryName = info[0];
						targetName = info[1];
						EW.println();
						printBasicInfo(EW,ID,info[0],info[1],info[10],info[2]);
						ID++;
						HSPID = 1;
					}
					else if(info[1].compareTo(targetName) != 0){ // if new target print names and info
						targetName = info[1];
						EW.println();
						printBasicInfo(EW,ID,info[0],info[1],info[10],info[2]);
						ID++;
						HSPID = 1;
					}
					printHSPinfo(EW,HSPID,info[6],info[7], info[8], info[9]);
					HSPID++;
				}
				else{
					if(info[0].compareTo(queryName) != 0){ // if new query print names and info
						queryName = info[0];
						targetName = info[2];
						EW.println();
						printBasicInfo(EW,ID,info[0],info[2],info[11],info[3]);
						ID++;
						HSPID = 1;
					}
					else if(info[2].compareTo(targetName) != 0){ // if new target print names and info
						targetName = info[2];
						EW.println();
						printBasicInfo(EW,ID,info[0],info[2],info[11],info[3]);
						ID++;
						HSPID = 1;
					}
					printHSPinfo(EW,HSPID,info[7],info[8], info[9], info[10]);
					HSPID++;
				}
			}
			EW.println();
			EW.flush();
			EW.close();
		}
		catch(Exception E){E.printStackTrace();}
	}



}
