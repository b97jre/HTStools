package general;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Random;


public class GetFrequencies {

	public static void main(String[] args) {
		int index;
		int k = 1;
		int len = 0;
		File infile = null;
		if ((index = indexOf("-k",args)) >= 0)
			k = Integer.parseInt(args[index+1]);
		if ((index = indexOf("-l",args)) >= 0)
			len = Integer.parseInt(args[index+1]);
		if ((index = indexOf("-i",args)) >= 0)
			infile = new File(args[index+1]);
		try {
			if (indexOf("-f",args) >= 0 && infile != null) {
				double[] freqs = getKmerFrequencies(parseBigFasta(infile),k);
				System.out.println("# "+String.valueOf(k));
				for (int i=0; i<freqs.length; i++)
				    System.out.println(String.valueOf(freqs[i]));
				System.out.println("# //");
				System.exit(0);
			}
			if (indexOf("-g",args) >= 0 && infile != null) {
				double[] freqs = parseFreqs(infile);
				byte[] seq = generateSequence(len,k-1,freqs);
				char[] decodedSeq = decodeToChar(seq,"NACGTN");
				int wrap = 60;
				int laps = seq.length/60;
				for (int i=0; i<laps; i++)
					System.out.println(new String(subarray(i*wrap,(i+1)*wrap,decodedSeq)));
				System.out.println(new String(subarray(laps*wrap,seq.length,decodedSeq)));
				System.exit(0);
			}
		} catch (Exception e) {e.printStackTrace();}
	}		

    public static double[] getKmerFrequencies(String[][] seqs, int k) {
	int n = 0;
	for (int i=1; i<=k; i++)
	    n += (int) Math.pow(4,i);
	double[] totalFreqs = new double[n];
	byte[] seq = new byte[0];
	for (int i=0; i<seqs.length; i++) {
		seq = append(addElement((byte) 0,encodeToByte(seqs[i][1])),seq);
	}
	totalFreqs = getKmerFrequencies(seq,k);
	return totalFreqs;
    }

    public static double[] getKmerFrequencies(byte[] seq, int k) {
	if (k <= 0)
	    return new double[0];
	int noMotifs = 0;
	for (int j=1; j<=k; j++)
	    noMotifs += (int) Math.pow(4,j);
	double[] n = new double[noMotifs];
	int[] count = new int[k];
	int index;
	int len;
	byte[] submotif;
	for (int i=0; i<(seq.length-k+1); i++) {
	    submotif = subarray(i,i+k,seq);
	    index = getKmerIndex(submotif);
	    len = k;
	    while (len > 0) {
		if (min(submotif) >= 1 && max(submotif) <= 4) {
		    n[index]++;
		    count[len-1]++;
		    index = (index/4)-1;
		}
		else
		    index = getKmerIndex(subarray(0,len-1,submotif));
		submotif[len-1] = 1;
		len--;
	    }
	}
	for (int i=k-1; i>0; i--) {
	    submotif = subarray(seq.length-i,seq.length,seq);
	    index = getKmerIndex(submotif);
	    len = i;
	    while (len > 0) {
		if (min(submotif) >= 1 && max(submotif) <= 4) {
		    n[index]++;
		    count[len-1]++;
		    index = (index/4)-1;
		}
		else
		    index = getKmerIndex(subarray(0,len-1,submotif));
		submotif[len-1] = 1;
		len--;
	    }
	}
	
	double[] freqs = new double[0];
	int lower;
	int upper;
	for (int i=1; i<=k; i++) {
	    lower = getKmerIndex(newByteArray(i,(byte) 1));
	    upper = getKmerIndex(newByteArray(i,(byte) 4));
	    freqs = append(vectorDivide(subarray(lower,upper+1,n),count[i-1]),freqs);
	}
	    
	return freqs;
    }		
  
    public static byte[] generateSequence(int len, int order, double[] freqs) {
		byte[] seq = new byte[len];
		double[][] kmers;
		double r;
		double[] cumFreqs;
		int next = 0;
		int j;
		int offset = -1;
		Random rand = new Random();
		while (next < order && next < len) {
			kmers = getConditionalKmers(freqs,next);
			r = rand.nextDouble();
			if (next == 0)
				cumFreqs = cumsum(kmers[0]);
			else
				cumFreqs = cumsum(kmers[getKmerIndex(subarray(0,next,seq))-offset]);
			j = 0;
			while (j < 4 && r >= cumFreqs[j])
				j++;
			seq[next] = (byte) Math.min((j+1),4);
			offset += (int) Math.pow(4,next);
			next++;
		}
		kmers = getConditionalKmers(freqs,order);
		for (int i=0; i<kmers.length; i++)
			kmers[i] = cumsum(kmers[i]);
			if (order == 0) {
				while (next < len) {
					j = 0;
					r = rand.nextDouble();
					while (j < 4 && r >= kmers[0][j])
						j++;
					seq[next] = (byte) Math.min((j+1),4);
					next++;
				}
			}
			else {
				byte[] prev = subarray(0,order,seq);
				int index;
				while (next < len) {
					index = getKmerIndex(prev)-offset;
					j = 0;
					r = rand.nextDouble();
					while (j < 4 && r >= kmers[index][j])
						j++;
					seq[next] = (byte) Math.min((j+1),4);
					prev = lshift(prev,(byte) Math.min((j+1),4));
					next++;
				}
			}
			return seq;
    }

    public static byte[] generateWeightedSequence(int len, int order, double[] freqs,double weightFreq, byte[] antisenseSeq,double[][] antisenseFreq){
		byte[] seq = new byte[len];
		double[][] kmers;
		double r;
		double[] cumFreqs;
		int next = 0;
		int j;
		int offset = -1;
		Random rand = new Random();
		while (next < order && next < len) {
			r = rand.nextDouble();
			if(r < weightFreq){
				kmers = getConditionalKmers(freqs,next);
				r = rand.nextDouble();
				if (next == 0)
					cumFreqs = cumsum(kmers[0]);
				else
					cumFreqs = cumsum(kmers[getKmerIndex(subarray(0,next,seq))-offset]);
				j = 0;
				while (j < 4 && r >= cumFreqs[j])
					j++;
				seq[next] = (byte) Math.min((j+1),4);
				offset += (int) Math.pow(4,next);
			}
			else{
				r = rand.nextDouble();
				seq[next] = getAntinsenseNucleotide(seq,antisenseSeq,next,antisenseFreq,r);
			}
			next++;
		}
		kmers = getConditionalKmers(freqs,order);
		for (int i=0; i<kmers.length; i++)
			kmers[i] = cumsum(kmers[i]);
				if (order == 0) {
					while (next < len) {
						r = rand.nextDouble();
						if(r < weightFreq){
					
							j = 0;
							r = rand.nextDouble();
							while (j < 4 && r >= kmers[0][j])
								j++;
							seq[next] = (byte) Math.min((j+1),4);
							next++;
						}
						else{
							
						}
				}
			}
			else {
				byte[] prev = subarray(0,order,seq);
				int index;
				while (next < len) {
					index = getKmerIndex(prev)-offset;
					j = 0;
					r = rand.nextDouble();
					while (j < 4 && r >= kmers[index][j])
						j++;
					seq[next] = (byte) Math.min((j+1),4);
					prev = lshift(prev,(byte) Math.min((j+1),4));
					next++;
				}
			}
		return seq;
    }
    
	
	public static byte getAntinsenseNucleotide(byte[] seq,byte[] antisenseSeq,int next,double[][] antisenseFreq, double r){
		int j = 0;
		while (j < 4 && r >= antisenseFreq[antisenseSeq[next]][j])
			j++;		
		return (byte)j;

	}
	

    
    public static int getKmerIndex(byte[] seq) {
	int k = seq.length;
	int index = -1;
	for (int j=1; j<=k; j++)
	    index += ((int) Math.pow(4,(j-1))) + ((int) Math.pow(4,(k-j)))*(seq[j-1]-1);
	return index;
    }  
  
    public static double[][] getConditionalKmers(double[] freqs, int order) {
	if (order == 0)
	    return new double[][] {{freqs[0],freqs[1],freqs[2],freqs[3]}};
	double[][] kmers = new double[(int) Math.pow(4,order)][4];
	int lower = getKmerIndex(newByteArray(order+1,(byte) 1));
	int upper = getKmerIndex(newByteArray(order+1,(byte) 4));
	int k = lower;
	for (int i=0; k<upper; i++) {
	    for (int j=0; j<4; j++) {
		if (freqs[k/4-1] != 0)
		    kmers[i][j] = freqs[k]/freqs[k/4-1];
		else
		    kmers[i][j] = 0;
  		k++;
	    }
	}
	return kmers;
    }  
    
	public static double[] cumsum(double[] arr) {
		double[] nArr = new double[arr.length];
		nArr[0] = arr[0];
		for (int i=1; i<arr.length; i++)
			nArr[i] = nArr[i-1]+arr[i];
		return nArr;
	}
	public static char[] decodeToChar(byte[] src, String cd) {
		char[] decoded = new char[src.length];
		for (int i=0; i<src.length; i++)
			decoded[i] = cd.charAt(src[i]);
		return decoded;
	}
	public static byte[] encodeToByte(String src) {
		String codeString = "NACGTN";
		src = src.toUpperCase();
		src = src.replace('T','U');
		codeString = codeString.replace('T','U');
		int sz = src.length();
		byte[] encoded = new byte[sz];
		for (int i=0; i<sz; i++)
			encoded[i] = (byte) Math.max(codeString.indexOf(src.charAt(i)),0);
		return encoded;
	}
	public static int indexOf(String obj, String[] arr) {
		if (arr == null || arr.length == 0)
			return -1;
		int sz = arr.length;
		for (int i=0; i<sz; i++)
			if (arr[i] != null && arr[i].compareTo(obj) == 0)
				return i;
		return -1;
	}

    public static byte[] lshift(byte[] arr, byte n) {
	for (int i=1; i<arr.length; i++)
	    arr[i-1] = arr[i];
	arr[arr.length-1] = n;
	return arr;
    }
	public static String[][] parseBigFasta(File f) throws Exception	{
		int INITIAL_CAPACITY = 100;
		int INCREMENT = 50;
		int nextEntry = -1;
		String[][] entries = new String[INITIAL_CAPACITY][2];
		BufferedReader br = new BufferedReader(new FileReader(f));
		int chunkSize = 10*1024*1024;
		int rd;
		char[] cBuff;
		char[] entry = new char[chunkSize];
		String name;
		int j = 0;
		while (true) {
			cBuff = new char[chunkSize];
			rd = br.read(cBuff,0,chunkSize);
			if (rd < 0)
			{
				if (nextEntry >= 0)
					entries[nextEntry][1] = new String(subarray(0,j,entry));
				nextEntry++;
				break;
			}
			cBuff = subarray(0,rd,cBuff);
			for (int i=0; i<rd; i++)
			{
				if (j >= entry.length)
					entry = append(new char[chunkSize],entry);
				if (cBuff[i] != 10 && cBuff[i] != 13 && cBuff[i] != 62 && cBuff[i] != 32)
					entry[j] = cBuff[i];
				else if (cBuff[i] == 10 || cBuff[i] == 13 || cBuff[i] == 32)
					j--;
				else
				{
					if (nextEntry >= 0)
						entries[nextEntry][1] = new String(subarray(0,j,entry));
					entry = new char[chunkSize];
					j = -1;
					i++;
					name = new String();
					while (i < rd && cBuff[i] != 10 && cBuff[i] != 13)
					{
						name += cBuff[i];
						i++;
					}
					nextEntry++;
					if (nextEntry == entries.length)
						entries = append(new String[INCREMENT][2],entries);
					entries[nextEntry][0] = name;
				}
				j++;
			}
		}
		br.close();
		return subarray(0,nextEntry,entries);
	}
	

	public static byte[] addElement(byte n, byte[] arr) {
		if (arr == null || arr.length == 0)
			return new byte[] {n};
		int sz = arr.length;
		byte[] newArr = new byte[sz+1];
		for (int i=0; i<sz; i++)
			newArr[i] = arr[i];
		newArr[sz] = n;
		return newArr;
	}
	public static byte[] append(byte[] nE, byte[] src) {
		if (src == null)
			src = new byte[0];
		int sLt = src.length;
		int nLt = nE.length;
		byte[] newArr = new byte[sLt+nLt];
		for (int i=0; i<sLt; i++)
			newArr[i] = src[i];
		for (int i=0; i<nE.length; i++)
			newArr[sLt+i] = nE[i];
		return newArr;
	}
	public static char[] append(char[] nE, char[] src) {
		if (src == null)
			src = new char[0];
		int sLt = src.length;
		int nLt = nE.length;
		char[] newArr = new char[sLt+nLt];
		for (int i=0; i<sLt; i++)
			newArr[i] = src[i];
		for (int i=0; i<nE.length; i++)
			newArr[sLt+i] = nE[i];
		return newArr;
	}
	public static double[] append(double[] nE, double[] src) {
		if (src == null)
			src = new double[0];
		int sLt = src.length;
		int nLt = nE.length;
		double[] newArr = new double[sLt+nLt];
		for (int i=0; i<sLt; i++)
			newArr[i] = src[i];
		for (int i=0; i<nE.length; i++)
			newArr[sLt+i] = nE[i];
		return newArr;
	}
	public static String[][] append(String[][] nE, String[][] src) {
		if (src == null || src.length == 0)
			return nE;
		int sLt = src.length;
		int nLt = nE.length;
		String[][] newArr = new String[sLt+nLt][];
		for (int i=0; i<sLt; i++)
			newArr[i] = src[i];
		for (int i=0; i<nLt; i++)
			newArr[sLt+i] = nE[i];
		return newArr;
	}

	public static byte max(byte[] arr) {return (byte) max(vectorMultiply(arr,1.0));}	
	public static double max(double[] arr) {return maxElement(arr);}
	public static double maxElement(double[] arr) {
		int sz = arr.length;
		if (sz == 0)
			return 0;
		double max = arr[0];
		for (int i=1; i<sz; i++)
			if (arr[i] > max)
				max = arr[i];
		return max;
	}

	public static byte min(byte[] arr) {return (byte) min(vectorMultiply(arr,1.0));}
	public static double min(double[] arr) {return -1.0*max(vectorMultiply(arr,-1.0));}
	public static byte[] newByteArray(int len, byte fill) {
		byte[] c = new byte[len];
		for (int i=0; i<len; i++)
			c[i] = fill;
		return c;
	}
	public static double[] parseFreqs(File infile) throws Exception {
		String[] lines = readFile(infile).split(System.getProperty("line.separator"));
		double[] freqs = new double[lines.length-1];
		int next = 0;
		while (next < lines.length && lines[next+1].compareTo("# //") != 0) {
		    freqs[next] = Double.parseDouble(lines[next+1]);
		    next++;
		}
		freqs = subarray(0,next,freqs);
		return freqs;
	}
	public static final String readFile(File f) throws Exception {
		char[] buff = new char[(int) f.length()];
		FileReader fr = new FileReader(f);
		fr.read(buff,0,buff.length);
		return new String(buff);
	}
	public static byte[] subarray(int start, int stop, byte[] src) {
		if (start >= src.length || stop > src.length)
			return new byte[0];
		byte[] newArr = new byte[stop-start];
		for (int i=start; i<stop; i++)
			newArr[i-start] = src[i];
		return newArr;
	}	
	public static char[] subarray(int start, int stop, char[] src) {
		if (start >= src.length || stop > src.length)
			return new char[0];
		char[] newArr = new char[stop-start];
		for (int i=start; i<stop; i++)
			newArr[i-start] = src[i];
		return newArr;
	}	
	public static double[] subarray(int start, int stop, double[] src) {
		if (start >= src.length || stop > src.length)
			return new double[0];
		double[] newArr = new double[stop-start];
		for (int i=start; i<stop; i++)
			newArr[i-start] = src[i];
		return newArr;
	}
	public static String[][] subarray(int start, int stop, String[][] src)	{
		if (start >= src.length || stop > src.length)
			return new String[0][0];
		String[][] newArr = new String[stop-start][];
		for (int i=start; i<stop; i++)
			newArr[i-start] = src[i];
		return newArr;
	}
	public static double[] vectorDivide(double[] v1, double f) {return vectorMultiply(v1,1./f);}
	public static double[] vectorMultiply(byte[] v1, double f) {
		double[] arr = new double[v1.length];
		for (int i=0; i<v1.length; i++)
			arr[i] = v1[i];
		return vectorMultiply(arr,f);
	}
	public static double[] vectorMultiply(double[] v1, double f)
	{
		int sz = v1.length;
		double[] nV = new double[sz];
		for (int i=0; i<sz; i++)
			nV[i] = v1[i] * f;
		return nV;
	}
}
