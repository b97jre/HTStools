package general;


import java.awt.FileDialog;
import java.awt.Frame;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.HashSet;


import javax.swing.JFileChooser;

public abstract class IOTools
{
	public static File[] addElement(File nO, File[] arr)
	{
		if (arr == null || arr.length == 0)
			return new File[] {nO};
		int sz = arr.length;
		File[] newArr = new File[sz+1];
		for (int i=0; i<sz; i++)
			newArr[i] = arr[i];
		newArr[sz] = nO;
		return newArr;
	}


	public static String getCurrentPath(){
		try{
			return new java.io.File(".").getCanonicalPath();
		}catch(Exception E){E.printStackTrace();}


		return ".";
	}

	
	
	public static int countLines(String filename) throws IOException {
		System.out.println("starting reading sequences in "+filename);
		long startTime = System.nanoTime();
		InputStream is = new BufferedInputStream(new FileInputStream(filename));
		try {
			byte[] c = new byte[1024];
			int count = 0;
			int readChars = 0;
			while ((readChars = is.read(c)) != -1) {
				for (int i = 0; i < readChars; ++i) {
					if (c[i] == '\n')
						++count;
				}
			}
			long endTime = System.nanoTime();
			long duration = endTime - startTime;
			System.out.println("finished. Took  "+duration/1000000000.0+" seconds");

			return count;
		} finally {
			is.close();
		}


	}

	public static void concatFile(File src, File concat) throws Exception
	{
		src.createNewFile();
		BufferedWriter bw = new BufferedWriter(new FileWriter(src,true));
		BufferedReader br = new BufferedReader(new FileReader(concat));
		String buff;
		while ((buff = br.readLine()) != null)
		{
			bw.write(buff);
			bw.newLine();
		}
		bw.close();
		br.close();
	}

	public static void deleteFiles(File[] tF)
	{
		int sz;
		if (tF == null || (sz = tF.length) == 0)
			return;
		for (int i=sz; i>0; i--)
			if (tF[i-1] != null)
				tF[i-1].delete();
	}


	public static String fixFileName(String name)
	{
		name = name.replace('\\','-');
		name = name.replace('/','-');
		name = name.replace(':','-');
		name = name.replace('*','-');
		name = name.replace('?','-');
		name = name.replace('"','-');
		name = name.replace('<','-');
		name = name.replace('>','-');
		name = name.replace('|','-');
		name = name.replace(' ','_');
		name = name.replace('\t', '_');
		name = name.replace('\n', '_');
		return name;
	}

	public static String fixFastaNames(String name)
	{
		name = name.replace('\\','-');
		name = name.replace('/','-');
		name = name.replace(':','-');
		name = name.replace('*','-');
		name = name.replace('?','-');
		name = name.replace('"','-');
		name = name.replace('<','-');
		name = name.replace('|','-');
		name = name.replace(' ','_');
		name = name.replace('\t', '_');
		name = name.replace('\n', '_');
		return name;
	}


	public static String fixFileName(String name, char ex)
	{
		name = name.replace('\\',ex);
		name = name.replace('/',ex);
		name = name.replace(':',ex);
		name = name.replace('*',ex);
		name = name.replace('?',ex);
		name = name.replace('"',ex);
		name = name.replace('<',ex);
		name = name.replace('>',ex);
		name = name.replace('|',ex);
		name = name.replace(' ',ex);
		name = name.replace('\t', ex);
		name = name.replace('\n', ex);
		return name;
	}



	public static String getContents(String s) throws Exception {return getContents(new File(s));}
	public static final String getContents(File f) throws Exception {return readFile(f);}

	public static final String getExtension(File f)
	{
		String path = f.getAbsolutePath();
		int index = path.lastIndexOf('.');
		if (index < 0)
			return new String();
		return path.substring(index+1);
	}

	public static final String getName(File f)
	{
		String path = f.getAbsolutePath();
		int index = path.lastIndexOf('.');
		if (index < 0)
			index = path.length();
		return path.substring(0,index);
	}

	public static final String getFileName(File f)
	{
		String name = f.getName();
		int index = name.lastIndexOf('.');
		if (index < 0)
			index = name.length();
		return name.substring(0,index);
	}

	public static final byte[] readByteFile(File f) throws Exception
	{
		if (f == null || !f.exists())
			return new byte[0];
		FileInputStream br = new FileInputStream(f);
		int len = (int) f.length();
		byte[] arr = new byte[len];
		br.read(arr,0,len);
		br.close();
		return arr;
	}

	public static final String readFile(String s) throws Exception {return readFile(new File(s));}
	public static final String readFile(File f) throws Exception
	{
		char[] buff = new char[(int) f.length()];
		FileReader fr = new FileReader(f);
		fr.read(buff,0,buff.length);
		return new String(buff);
	}

	public static final Object readObject() throws Exception
	{
		File f = selectFile();
		if (f == null || !f.exists())
			return null;
		return readObject(f);
	}
	public static final Object readObject(File f) throws Exception
	{
		ObjectInputStream iS = new ObjectInputStream(new FileInputStream(f));
		Object obj = iS.readObject();
		iS.close();
		return obj;
	}

	public static final File selectFile() {return selectFile(new Frame());}
	public static final File selectFile(Frame parent)
	{
		FileDialog fD = new FileDialog(parent,"Select a file");
		fD.setVisible(true);;
		String file = fD.getFile();
		String dir = fD.getDirectory();
		if (file == null || dir == null)
			return null;
		return new File(dir+file);
	}
	public static final File[] selectFiles() {return selectFiles(new Frame());}
	public static final File[] selectFiles(boolean useSwing) {return selectFiles(new Frame(), useSwing);}
	public static final File[] selectFiles(Frame parent) {return selectFiles(parent,false);}
	public static final File[] selectFiles(Frame parent, boolean useSwing)
	{
		File[] tF = new File[0];
		if (!useSwing)
		{
			FileDialog fD = new FileDialog(parent,"Select a file");
			String file,dir;
			while (true)
			{
				fD.setVisible(true);
				file = fD.getFile();
				dir = fD.getDirectory();
				if (file == null || dir == null)
					break;
				tF = addElement(new File(dir+file),tF);
			}
		}
		else
		{
			JFileChooser fileChooser = new JFileChooser("");
			fileChooser.setDialogTitle("Select files");
			fileChooser.setDialogType(JFileChooser.OPEN_DIALOG);
			fileChooser.setMultiSelectionEnabled(true);
			if (fileChooser.showDialog(null,"Ok") != JFileChooser.APPROVE_OPTION)
				return null;
			tF = fileChooser.getSelectedFiles();
		}
		return tF;
	}


	public static String setExtension(File f, String ext) {return setExtension(f.getPath(),ext);}
	public static String setExtension(String name, String ext)
	{
		int index = name.lastIndexOf('.');
		if (index < 0)
			return name + '.' + ext;
		return name.substring(0,index+1) + ext;
	}

	public static String setSuffix(File f, String s)
	{
		String name = f.getPath();
		int index = name.lastIndexOf('.');
		if (index < 0)
			name += '_' + s;
		else
			name = name.substring(0,index) + '_' + s + name.substring(index);
		return name;
	}

	public static final File[] splitFile(File src, int sz) throws Exception
	{
		sz *= 1024*1024;
		int remSize = (int) src.length();
		if (remSize <= sz)
			return new File[] {src};
		FileReader fr = new FileReader(src);
		FileWriter fw;
		File[] splittedFiles = new File[1];
		String name = src.getName();
		int index = name.lastIndexOf('.');
		if (index < 0)
			index = name.length();
		String extension = name.substring(index,name.length());
		name = name.substring(0,index);
		File tF = new File(src.getParent()+"\\"+name+"\\");
		tF.mkdir();
		name = tF.getPath()+"\\"+name;
		splittedFiles[0] = tF;
		int fileCounter = 1,rd;
		char[] buff;
		while (remSize > 0)
		{
			tF = new File(name+"_"+String.valueOf(fileCounter)+extension);
			fw = new FileWriter(tF);
			buff = new char[sz];
			rd = fr.read(buff,0,sz);
			remSize -= rd;
			fw.write(buff,0,rd);
			if (rd > 0)
			{
				while ((rd = fr.read()) != '\n' && rd != -1)
				{
					fw.write(rd);
					remSize--;
				}
				if (rd != -1)
					fw.write(rd);
				remSize--;
			}
			fw.close();
			splittedFiles = addElement(tF,splittedFiles);
			fileCounter++;
		}
		fr.close();
		return splittedFiles;
	}

	public static final void writeByteFile(byte[] arr, File f) throws Exception
	{
		if (f == null)
			return;
		FileOutputStream bw = new FileOutputStream(f);
		bw.write(arr,0,arr.length);
		bw.close();
	}

	public static final void writeObject(Object obj) throws Exception
	{
		File f = selectFile();
		if (f == null)
			return;
		writeObject(obj,f,false);
	}
	public static final void writeObject(Object obj, File f) throws Exception {writeObject(obj,f,false);}
	public static final void writeObject(Object obj, File f, boolean append) throws Exception
	{
		ObjectOutputStream oS = new ObjectOutputStream(new FileOutputStream(f,append));
		oS.writeObject(obj);
		oS.flush();
		oS.close();
	}

	public static final void writeToFile(String contents) throws Exception {writeToFile(contents,selectFile(new Frame()));}
	public static final void writeToFile(String contents, File destFile) throws Exception
	{
		if (contents == null || destFile == null)
			return;
		BufferedWriter bw = new BufferedWriter(new FileWriter(destFile));
		bw.write(contents);
		bw.flush();
		bw.close();
	}
	public static final void writeToFile(ArrayList<Object> arr, File destFile) throws Exception
	{
		if (arr == null || destFile == null)
			return;
		BufferedWriter bw = new BufferedWriter(new FileWriter(destFile));
		int sz = arr.size();
		for (int i=0; i<sz; i++)
		{
			bw.write((arr.get(i)).toString());
			bw.newLine();
		}
		bw.close();
	}

	public static void copy(String src, String dst) throws IOException {
		InputStream in = new FileInputStream(src);
		OutputStream out = new FileOutputStream(dst);

		// Transfer bytes from in to out
		byte[] buf = new byte[1024];
		int len;
		while ((len = in.read(buf)) > 0) {
			out.write(buf, 0, len);
		}
		in.close();
		out.close();
	}


	public static void copy(File src, File dst) throws IOException {
		InputStream in = new FileInputStream(src);
		OutputStream out = new FileOutputStream(dst);

		// Transfer bytes from in to out
		byte[] buf = new byte[1024];
		int len;
		while ((len = in.read(buf)) > 0) {
			out.write(buf, 0, len);
		}
		in.close();
		out.close();
	}

	public static void mkDirs(String directoryName){
		String[] subDirectories = directoryName.split("/");
		String Dir="";
		for(int i = 0; i < subDirectories.length; i++){
			Dir+=subDirectories[i];
			if(!isDir(Dir)){
				mkDir(Dir);
			}
			Dir+="/";
		}
		boolean success = (new File(directoryName)).mkdir();
		if(!success){
			// Directory creation failed
		}
	}

	public static void mkDir(String directoryName){
		boolean success = (new File(directoryName)).mkdir();
		if (!success) {
			// Directory creation failed
		}
	}


	public static boolean isDir(String directoryName){
		boolean success = (new File(directoryName)).isDirectory();
		return success;
	}



	public static boolean fileExists(String dir, String fileName){
		return (new File(dir+"/"+fileName)).exists();
	}

	public static boolean deleteFile(String dir, String fileName){
		return (new File(dir+"/"+fileName)).delete();
	}

	public static boolean fileExists(String fileName){
		return (new File(fileName)).exists();
	}



	public static ArrayList<String> getDirectories(String Dir){
		ArrayList<String> dirs = new ArrayList<String>();
		File dir = new File(Dir);

		String[] children = dir.list();
		File[] files = dir.listFiles();

		if (children == null) {
			return null;
			// Either dir does not exist or is not a directory
		} else {
			for (int i=0; i<children.length; i++) {
				// Get filename of file or directory
				if(files[i].isDirectory()){
					dirs.add(children[i]);
				}
			}
		}
		return dirs;
	}


	public static ArrayList<String> getSequenceFiles(String Dir){
		ArrayList<String> SequenceFiles = new ArrayList<String>();
		File dir = new File(Dir);

		String[] children = dir.list();
		File[] files = dir.listFiles();
		if (children == null) {
			return null;
			// Either dir does not exist or is not a directory
		} else {
			for (int i=0; i<children.length; i++) {
				// Get filename of file or directory
				String filename = children[i];

				if(!files[i].isDirectory()){
					SequenceFiles.add(filename);
				}
			}
		}
		return SequenceFiles;
	}



	public static String longestCommonPrefix(ArrayList<String> strs) {
		String prefix = new String();
		if(strs.size() > 0)
			prefix = strs.get(0);
		for(int i = 1; i < strs.size(); ++i) {
			String s = strs.get(i);
			int j = 0;
			for(; j < Math.min(prefix.length(), s.length()); ++j) {
				if(prefix.charAt(j) != s.charAt(j)) {
					break;
				}
			}
			prefix = prefix.substring(0, j);
		}
		return prefix;
	}

	public static String longestCommonPrefix(String A, String B) {
		String prefix = A;
		int j = 0;   
		for(; j < Math.min(prefix.length(), B.length()); ++j) {
			if(prefix.charAt(j) != B.charAt(j)) {
				break;
			}
		}
		prefix = prefix.substring(0, j);
		return prefix;
	}
	
	public static String longestCommonSuffix(String A, String B) {
		String suffix = A;
		int pointer = A.length()-1;
		int pointer2 = B.length()-1;
		int j = 0;
		for(; j < Math.min(suffix.length(), B.length()); ++j) {
			if(suffix.charAt(pointer-j) != B.charAt(pointer2-j)) {
				break;
			}
		}
		suffix = suffix.substring(pointer-j+1);
		return suffix;
	}
	
	public static String longestCommonSuffix(ArrayList<String> strs) {
		String suffix = strs.get(0);
		for(int j = 1; j < strs.size(); ++j) {
			suffix = longestCommonSuffix(suffix,strs.get(j));
		}
		return suffix;
	}


	public static String removeLastDot(String fileName){
		String inFileBase = fileName;
		while(inFileBase.lastIndexOf('.') == inFileBase.length()-1){
			inFileBase = inFileBase.substring(0,inFileBase.length()-1);
		}
		return inFileBase;
	}

	public static String removeFirstDot(String fileName){
		String inFileBase = fileName;
		while(inFileBase.indexOf('.') == 0){
			inFileBase = inFileBase.substring(1);
		}
		return inFileBase;
	}




	public static String getFileBase(String fileName, String suffix){
		if(suffix != null){
			String inFileBase = fileName.substring(0,fileName.lastIndexOf(suffix));
			while(inFileBase.lastIndexOf('.') == inFileBase.length()-1){
				inFileBase = inFileBase.substring(0,inFileBase.length()-1);
			}	
			return inFileBase;
		}else{
			String inFileBase = fileName.substring(0,fileName.lastIndexOf('.'));
			while(inFileBase.lastIndexOf('.') == inFileBase.length()-1){
				inFileBase = inFileBase.substring(0,inFileBase.length()-1);
			}	
			return inFileBase;
		}
	}
	
	


	public static ArrayList<String> getSequenceFilesFullPath(String Dir, String suffix){
		ArrayList<String> SequenceFiles = new ArrayList<String>();
		File dir = new File(Dir);

		String[] children = dir.list();
		if (children == null) {
			return null;
			// Either dir does not exist or is not a directory
		} else {
			for (int i=0; i<children.length; i++) {
				// Get filename of file or directory
				String filename = children[i];

				if(filename.indexOf(suffix) > -1 && filename.lastIndexOf(suffix) + suffix.length() == filename.length()){
					SequenceFiles.add(Dir+"/"+filename);
				}
			}
		}
		return SequenceFiles;
	}




	public static ArrayList<String> getFilesSuffix(String Dir, String suffix){
		ArrayList<String> SequenceFiles = new ArrayList<String>();
		File dir = new File(Dir);

		String[] children = dir.list();
		if (children == null) {
			return null;
			// Either dir does not exist or is not a directory
		} else {
			for (int i=0; i<children.length; i++) {
				// Get filename of file or directory
				String filename = children[i];
				if(suffix == null)
					SequenceFiles.add(filename);
				else if(filename.indexOf(suffix) > -1 && filename.lastIndexOf(suffix) + suffix.length() == filename.length()){
					SequenceFiles.add(filename);
				}
			}
		}
		return SequenceFiles;
	}


	public static ArrayList<String> getFilesPrefix(String Dir, String prefix){
		ArrayList<String> SequenceFiles = new ArrayList<String>();
		File dir = new File(Dir);

		String[] children = dir.list();
		if (children == null) {
			return null;
			// Either dir does not exist or is not a directory
		} else {
			for (int i=0; i<children.length; i++) {
				// Get filename of file or directory
				String filename = children[i];
				if(prefix == null)
					SequenceFiles.add(filename);
				else if(filename.indexOf(prefix) == 0 ){
					SequenceFiles.add(filename);
				}
			}
		}
		return SequenceFiles;
	}



	public static <T> List<T> union(List<T> list1, List<T> list2) {
		Set<T> set = new HashSet<T>();

		set.addAll(list1);
		set.addAll(list2);

		return new ArrayList<T>(set);
	}

	public static  <T> List<T> intersection(List<T> list1, List<T> list2) {
		List<T> list = new ArrayList<T>();

		for (T t : list1) {
			if(list2.contains(t)) {
				list.add(t);
			}
		}

		return list;
	}



	public static ArrayList<String> getSequenceFiles(String Dir, String suffix){
		ArrayList<String> SequenceFiles = new ArrayList<String>();
		File dir = new File(Dir);

		String[] children = dir.list();
		if (children == null) {
			return null;
			// Either dir does not exist or is not a directory
		} else {
			for (int i=0; i<children.length; i++) {
				// Get filename of file or directory
				String filename = children[i];

				if(filename.indexOf(suffix) > -1 && filename.lastIndexOf(suffix) + suffix.length() == filename.length()){
					SequenceFiles.add(filename);
				}
			}
		}
		return SequenceFiles;
	}




	public static ArrayList<String> getSequenceFilesPrefix(String Dir, String prefix){
		ArrayList<String> SequenceFiles = new ArrayList<String>();
		File dir = new File(Dir);

		String[] children = dir.list();
		if (children == null) {
			return null;
			// Either dir does not exist or is not a directory
		} else {
			for (int i=0; i<children.length; i++) {
				// Get filename of file or directory
				String filename = children[i];

				if(filename.indexOf(prefix) == 0 ){
					SequenceFiles.add(filename);
				}
			}
		}
		return SequenceFiles;
	}



	public static ArrayList <String[]> findPairs(ArrayList <String> fileNames,  String [] sep){

		ArrayList <String[]> pairs = new ArrayList <String[]>();
		for(int i = 0; i < fileNames.size(); i++){	
			String fName = fileNames.get(i);
			if(fName.contains(sep[0])){
				String[] fileName = fName.split(sep[0]);
				for(int j = 0; j < fileNames.size();j++){
					String fName2 = fileNames.get(j);
					if(fName2.contains(sep[1])){
						String[] fileName2 = fName2.split(sep[1]);
						int count = 0;
						for(int k = 0; k < fileName.length; k++){
							//							System.out.println(fileName[k] +"\t"+fileName2[k]);
							if(fileName[k].compareTo(fileName2[k]) != 0){
								count++;
							}
						}
						if(count ==0 ){
							String[] pair = new String[3];
							pair[0] = fileNames.get(i);
							pair[1] = fileNames.get(j);
							pair[2] = fileName[0];


							pairs.add(pair);
							if(i < j){
								fileNames.remove(j);
								fileNames.remove(i);
							}else{
								fileNames.remove(i);
								fileNames.remove(j);
								i--;
							}
							j = fileNames.size();
							i--;
						}
					}
				}
			}
		}
		return pairs;

	}

	public static void copyFile(File source, File dest) throws IOException {
		FileChannel sourceChannel = null;
		FileChannel destChannel = null;
		try {
			sourceChannel = new FileInputStream(source).getChannel();
			destChannel = new FileOutputStream(dest).getChannel();
			destChannel.transferFrom(sourceChannel, 0, sourceChannel.size());
		}finally{
			sourceChannel.close();
			destChannel.close();
		}
	}

}
