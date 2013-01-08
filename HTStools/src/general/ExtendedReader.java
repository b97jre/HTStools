package general;
import java.io.*;
import java.text.*;
import java.util.*;

public class ExtendedReader extends BufferedReader {
  // texts for error messages
  private String illegalInt;
  private String illegalNum;
  private String eofText;
 
  private NumberFormat nf; 

  // constructors
  public ExtendedReader (Reader  r) {
    super(r);
    setFormat(NumberFormat.getInstance());
    setLang(Locale.getDefault().getLanguage());
  }

  public ExtendedReader (Reader  r, int bufferSize) {
    super(r, bufferSize);
    setFormat(NumberFormat.getInstance());
    setLang(Locale.getDefault().getLanguage());  }

  // private instance methods
  private void flush() {
    if (this == Std.in)
      Std.out.flush();
  }

  private static void abort(String message) {
    Std.err.println(); Std.err.println();
    Std.err.println(message);
    System.exit(1);
  }

  // public static methods
  public static ExtendedReader getFileReader(String fileName) {
    try {
      return new ExtendedReader(new FileReader(fileName));
    }
    catch (FileNotFoundException e) {
      return null;
    }
  }  
  
  // public instance methods
  public int readChar() {
    try {
      flush();
      return read();  // returns -1 at end of file
    }  
    catch (IOException e) {
      e.printStackTrace();
      System.exit(1);
    }
    return 0; // will never get here
  }

  public int lookAhead() {
    try {
      flush();
      mark(2);
      int i = read();
      if (i != -1) 
        reset();  
      return i;  // returns -1 at end of file
    }
    catch (IOException e) {
      e.printStackTrace();
      System.exit(1);
    }
    return 0;  // will never get here
  }

  public long skip(long n) {  // returns the number of skipped characters
    flush();
    try {
      return super.skip(n);
    }
    catch (IOException e) {
      e.printStackTrace();
      System.exit(1);
    }
    return 0;  // will never get here
  }

  public void skipLine() {
    flush();
    try {  
      int c;
      while ((c = read()) != -1 && c != '\r' && c != '\n')
        ;
      if (c == '\r' && lookAhead() == '\n')  // MS-DOS 
        skip(1);
    }
    catch (IOException e) {
      e.printStackTrace();
      System.exit(1);
    }
  }

  public boolean more() { // skips blanks, gives false at end of file
    int c;
    while ((c = lookAhead()) == ' ' || 
           c == '\t' || c == '\r' || c == '\n')
      skip(1);
    return c != -1;
  }

  public String readLine() {
    try {  
      flush();
      return super.readLine();  // returns null at end of file
    }
    catch (IOException e) {
      e.printStackTrace();
      System.exit(1);
    }
    return null;  // will never get here
  }

  public String readWord() {
    if (more()) {
      String s = "";
      int c;
      while ((c = lookAhead()) != -1 && c != ' ' && 
              c != '\t' && c != '\r' && c != '\n')
        s = s + (char) readChar();
       return s;
    }
     return null;      
  }
                    
  public Number readNumber(String errMess) { 
    // returns null at end of file
    while (true) {         
      String s = readWord();
      if (s == null)   // end of file
        return null;
      ParsePosition pos = new ParsePosition(0);
      Number n = nf.parse(s, pos);
      if (n != null && pos.getIndex() == s.length())
        return n;
      else {
        skipLine();
        Std.err.println(errMess);
      }
    }
  }                       

  public int readInt() {
    nf.setParseIntegerOnly(true);
    Number n = readNumber(illegalInt);
    nf.setParseIntegerOnly(false);
    if (n == null)      
      abort(eofText+"readInt");
    return n.intValue();
  }

  public long readLong() {
    nf.setParseIntegerOnly(true);
    Number n = readNumber(illegalInt);
    nf.setParseIntegerOnly(false);
    if (n == null)      
      abort(eofText+"readLong");
    return n.longValue();
  }
    
  public double readDouble() {
    Number n = readNumber(illegalNum);
    if (n == null)      
      abort(eofText+"readDouble");
    return n.doubleValue();  
  }

  public void close() {
    try {
      super.close();
    }
    catch (IOException e) {
      e.printStackTrace();
      System.exit(1);
    }
  }

  public void setFormat(NumberFormat f) {
    nf = (NumberFormat) f.clone();
    nf.setGroupingUsed(false);
  }
    
  public NumberFormat getFormat() {
    return nf;
  }

  public void setLang(String lang) {
    // initilize texts for error messages      
    if (lang.equals("sv")) { // Swedish
      illegalInt = "Felaktigt heltal. Fšrsšk igen";
      illegalNum = "Felaktigt tal. Fšrsšk igen";
      eofText    = "Filslut i metoden ";
    }
    // or whatever language you wish ...
    else { // English default
      illegalInt = "Illegal whole number. Try again";
      illegalNum = "Illegal number. Try again";
      eofText    = "End of file in ";
    }
  } 
}