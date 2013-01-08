package general;

import java.io.File;
import java.io.FileReader;

public class ClusterBlast
{
	
	public static void main(String[] args)
	{
		try{
			File blastFile = new File(args[0]);
			String[] cutoffs = new String[17];
			if (args.length > 1)
				cutoffs = subarray(1,args);
			int[][] groups = getGroups(blastFile,cutoffs);
			for (int i=0; i<groups.length; i++)
			{
				if(groups[i].length > 100){
				for (int j=0; j<groups[i].length; j++)
					System.out.print(String.valueOf(groups[i][j])+"\t");
				System.out.println();
			}
			}
		} catch (Exception e) {e.printStackTrace();}
	}
	
	// blastOutput is a tab separated BLAST result file (tabular output)
	// cutoffs is an array with cutoff values according to:
	//	[0] -> Maximum E-value
	//	[1] -> Minimum bit score
	//	[2] -> Minimum raw score
	//	[3] -> Minimum total aligned length
	//	[4] -> Minimum number of identical residues
	//	[5] -> Maximum number of mismatches
	//	[6] -> Minimum percent identity
	//	[7] -> Minimum percent positive scores
	//	[8] -> Maximum number of gaps in query
	//	[9] -> Maximum number of gaps in subject
	//	[10] -> Query frame (without '+'-sign)
	//	[11] -> Nothing
	//	[12] -> Nothing
	// 	[13] -> Subject frame (without '+'-sign)
	//	[14] -> Nothing
	//	[15] -> Nothing
	//      [16] -> "O" Means that query and subject should have the opposite orientation. "S" means that query and subject should have the same orientation 
	public static int[][] getGroups(File blastOutput, String[] cutoffs) throws Exception
	{
		String contents = readFile(blastOutput);
		String[] lines = contents.split(System.getProperty("line.separator"));
		int sz = lines.length;
		
		int regId, targetId, alignlength, identical, mismatches, querygaps, subjectgaps;
//		int querystart, queryend, subjectstart, subjectend;
		double evalue, bitscore, rawscore, percentidentity, percentpositive;
		String queryframe, subjectframe;
		boolean passed;
		String[] fields;
		
		int capacity = 1000;
		int next = 0;
		int n = 0;
		int[][] groups = new int[capacity][capacity];
		
		for (int i=0; i<sz; i++)
		{
			fields = lines[i].split("\\t+");
			if (fields.length >= 22)
			{
				regId = Integer.parseInt(fields[0]);
				targetId = Integer.parseInt(fields[1]);
				
				evalue = Double.parseDouble(fields[2]);
				bitscore = Double.parseDouble(fields[4]);
				rawscore = Double.parseDouble(fields[5]);
				alignlength = Integer.parseInt(fields[6]);
				identical = Integer.parseInt(fields[7]);
				mismatches = Integer.parseInt(fields[9]);
				percentidentity = Double.parseDouble(fields[10]);
				percentpositive = Double.parseDouble(fields[11]);
				querygaps = Integer.parseInt(fields[13]);
				subjectgaps = Integer.parseInt(fields[15]);
				queryframe = fields[16];
				if (queryframe.charAt(0) == '+')
					queryframe = queryframe.substring(1);
//				querystart = Integer.parseInt(fields[17]);
//				queryend = Integer.parseInt(fields[18]);
				subjectframe = fields[19];
				if (subjectframe.charAt(0) == '+')
					subjectframe = subjectframe.substring(1);
//				subjectstart = Integer.parseInt(fields[20]);
//				subjectend = Integer.parseInt(fields[21]);
				passed = true;
	
				for (int j=0; j<cutoffs.length && passed; j++)
				{
					if (cutoffs[j] != null)
					{
						switch (j)
						{
							case 0:
								if (evalue > Double.parseDouble(cutoffs[j]))
									passed = false;
								break;
							case 1:
								if (bitscore < Double.parseDouble(cutoffs[j]))
									passed = false;
								break;
							case 2:
								if (rawscore < Double.parseDouble(cutoffs[j]))
									passed = false;
								break;
							case 3:
								if (alignlength < Integer.parseInt(cutoffs[j]))
									passed = false;
								break;
							case 4:
								if (identical < Integer.parseInt(cutoffs[j]))
									passed = false;
								break;
							case 5:
								if (mismatches > Integer.parseInt(cutoffs[j]))
									passed = false;
								break;
							case 6:
								if (percentidentity < Double.parseDouble(cutoffs[j]))
									passed = false;
								break;
							case 7:
								if (percentpositive < Double.parseDouble(cutoffs[j]))
									passed = false;
								break;
							case 8:
								if (querygaps > Integer.parseInt(cutoffs[j]))
									passed = false;
								break;
							case 9:
								if (subjectgaps > Integer.parseInt(cutoffs[j]))
									passed = false;
								break;
							case 10:
								if (queryframe.compareTo(cutoffs[j]) != 0)
									passed = false;
								break;
							case 13:
								if (subjectframe.compareTo(cutoffs[j]) != 0)
									passed = false;
								break;
							case 16:
								int sign = Integer.parseInt(queryframe)*Integer.parseInt(subjectframe);
								if (sign > 0 && cutoffs[j].compareTo("O") == 0)
									passed = false;
								else if (sign < 0 && cutoffs[j].compareTo("S") == 0)
									passed = false;
								break;
						}
					}
				}
				if (passed)
				{
					if (regId != groups[next][0])
					{
						if (next > 0 || groups[next][0] > 0)
						{
							groups[next] = subarray(0,n,groups[next]);
							next++;
						}
						n = 0;
						if (next == groups.length)
							groups = append(new int[capacity][capacity],groups);
						groups[next][n] = regId;
						n++;
						groups[next][n] = Integer.parseInt(queryframe);
						n++;
					}
					if (n > groups[next].length-2)
						groups[next] = append(new int[capacity],groups[next]);
					groups[next][n] = targetId;
					n++;
					groups[next][n] = Integer.parseInt(subjectframe)*Integer.parseInt(subjectframe)*Integer.parseInt(queryframe);
					n++;
				}
			}
		}
		groups[next] = subarray(0,n,groups[next]);
		groups = subarray(0,next,groups);
		return groups;
	}
	
	public static int[][] groupGroups(int[][] groups)
	{
		int id, dir;
		int[][] occurrencies;
		int[][] grouped = new int[groups.length][];
		int[][][] result;
		int next = 0;
		int index;
		for (int i=0; i<groups.length; i++)
		{
			if (groups[i].length > 0)
			{
				id = groups[i][0];
				dir = groups[i][1];
				grouped[next] = new int[] {id,dir};
				for (int j=2; j<groups[i].length; j+=2)
				{
					if ((index = indexOf(groups[i][j],grouped[next])) < 0 || index%2 != 0)
						grouped[next] = append(new int[] {groups[i][j],groups[i][j+1]},grouped[next]);
				}
				for (int j=2; j<grouped[next].length; j+=2)
				{
					result = getOccurrencies(grouped[next][j],groups);
					occurrencies = result[0];
					groups = result[1];
					for (int k=0; k<occurrencies.length; k++)
					{
						for (int l=2; l<occurrencies[k].length; l+=2)
						{
							if ((index = indexOf(occurrencies[k][l],grouped[next])) < 0 || index%2 != 0)
								grouped[next] = append(new int[] {occurrencies[k][l],occurrencies[k][l+1]*occurrencies[k][1]*grouped[next][j+1]},grouped[next]);
						}
					}
				}
				next++;
			}
		}
		grouped = subarray(0,next,grouped);
		return grouped;
	}

	public static int[][][] getOccurrencies(int id, int[][] groups)
	{
		int[][] occurrencies = new int[0][];
		int index, t1, t2;
	
		for (int i=0; i<groups.length; i++)
		{
			if ((index = indexOf(id,groups[i])) >= 0 && index%2 == 0)
			{
				occurrencies = addElement(groups[i],occurrencies);
				if (index > 0)
				{
					t1 = occurrencies[occurrencies.length-1][0];
					t2 = occurrencies[occurrencies.length-1][1];
					occurrencies[occurrencies.length-1][0] = occurrencies[occurrencies.length-1][index];
					occurrencies[occurrencies.length-1][1] = occurrencies[occurrencies.length-1][index+1];
					occurrencies[occurrencies.length-1][index] = t1;
					occurrencies[occurrencies.length-1][index+1] = t2;
				}
				groups[i] = new int[0];
			}
		}
		return new int[][][] {occurrencies,groups};
	}
					
	public static final String readFile(File f) throws Exception
	{
		char[] buff = new char[(int) f.length()];
		FileReader fr = new FileReader(f);
		fr.read(buff,0,buff.length);
		return new String(buff);
	}

	public static int[][] addElement(int[] nE, int[][] arr)
	{
		if (arr == null || arr.length == 0)
			return new int[][] {nE};
		int sz = arr.length;
		int[][] newArr = new int[sz+1][];
		for (int i=0; i<sz; i++)
			newArr[i] = arr[i];
		newArr[sz] = nE;
		return newArr;
	}

	public static int[] append(int[] nE, int[] src)	{return append(0,nE.length,nE,src);}
	public static int[] append(int start, int stop, int[] nE, int[] src)
	{
		if (src == null)
			src = new int[0];
		int sLt = src.length;
		int nLt = stop-start;
		int[] newArr = new int[sLt+nLt];
		for (int i=0; i<sLt; i++)
			newArr[i] = src[i];
		for (int i=start; i<stop; i++)
			newArr[sLt+i] = nE[i];
		return newArr;
	}	
	public static int[][] append(int[][] nE, int[][] src)	{return append(0,nE.length,nE,src);}
	public static int[][] append(int start, int stop, int[][] nE, int[][] src)
	{
		if (src == null)
			src = new int[0][0];
		int sLt = src.length;
		int nLt = stop-start;
		int[][] newArr = new int[sLt+nLt][];
		for (int i=0; i<sLt; i++)
			newArr[i] = src[i];
		for (int i=start; i<stop; i++)
			newArr[sLt+i] = nE[i];
		return newArr;
	}

	public static int indexOf(int obj, int[] arr)
	{
		if (arr == null || arr.length == 0)
			return -1;
		int sz = arr.length;
		for (int i=0; i<sz; i++)
			if (arr[i] == obj)
				return i;
		return -1;
	}
	
	public static int[][] subarray(int start, int[][] src) {return subarray(start,src.length,src);}
	public static int[][] subarray(int start, int stop, int[][] src)
	{
		if (start >= src.length || stop > src.length)
			return new int[0][0];
		int[][] newArr = new int[stop-start][];
		for (int i=start; i<stop; i++)
			newArr[i-start] = src[i];
		return newArr;
	}
	
	public static int[] subarray(int start, int[] src) {return subarray(start,src.length,src);}
	public static int[] subarray(int start, int stop, int[] src)
	{
		if (start >= src.length || stop > src.length)
			return new int[0];
		int[] newArr = new int[stop-start];
		for (int i=start; i<stop; i++)
			newArr[i-start] = src[i];
		return newArr;
	}
	public static String[] subarray(int start, String[] src) {return subarray(start,src.length,src);}
	public static String[] subarray(int start, int stop, String[] src)
	{
		if (start >= src.length || stop > src.length)
			return new String[0];
		String[] newArr = new String[stop-start];
		for (int i=start; i<stop; i++)
			newArr[i-start] = src[i];
		return newArr;
	}
}
