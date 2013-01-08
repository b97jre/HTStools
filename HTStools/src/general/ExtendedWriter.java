package general;

import java.io.*;
import java.text.*;

public class ExtendedWriter extends PrintWriter {
  private char fillChar = ' ';
  private NumberFormat nf;
  private static NumberFormat snf;

  static {
    setStaticFormat(NumberFormat.getInstance());
  }

  // constructors
  public ExtendedWriter(Writer w) {
    super(w);
    setFormat(NumberFormat.getInstance());
  }

  public ExtendedWriter(Writer w, boolean autoFlush) {
    super(w, autoFlush);
    setFormat(NumberFormat.getInstance());
  }

  public ExtendedWriter(OutputStream out) {
    super(out);
    setFormat(NumberFormat.getInstance());
  }

  public ExtendedWriter(OutputStream out, boolean autoFlush) {
    super(out, autoFlush);
    setFormat(NumberFormat.getInstance());
  }

 // private class methods

  private static String formatFixed(double num, int noPos, int noDecMax, 
                           int noDecMin, char fill, NumberFormat r) {
    int max = r.getMaximumFractionDigits();
    int min = r.getMinimumFractionDigits();
    r.setMaximumFractionDigits(noDecMax);
    r.setMinimumFractionDigits(noDecMin);
    String s = adjustRight(r.format(num), noPos, fill);
    r.setMaximumFractionDigits(max);
    r.setMinimumFractionDigits(min);
    return s;
  }

  private static String formatExponent(double num, int noPos, 
                              int noDec, char fill, NumberFormat r) {
    String s;
    if (Math.abs(num) < 1.0e-3 || Math.abs(num) > 1.0e7) // these values are in exponent form
      s = Double.toString(num);
    else if (Math.abs(num) < 1.0e0) {  // multiply by 1.0e-10 to get exponent form
      s = Double.toString(num*1.0e-10);
      s = s.substring(0, s.indexOf("E")) + "E-" + s.charAt(s.length()-1);
    }
    else {       // multiply by 1.0e10 to get exponent form
      s = Double.toString(num*1.0e10);
      s = s.substring(0, s.indexOf("E")) + "E" + s.charAt(s.length()-1);
    }
    if (s.indexOf("E")<0)
      s = s + "E0";

    // arrange the correct number of decimals
    String digits;
    if (num >=0 )
      digits = s.substring(0, s.indexOf("E"));
    else
      digits = s.substring(1, s.indexOf("E"));
    if (noDec < 0)
      noDec = 0;
    if (noDec < digits.length()-2) {  // round to noDec decimals
      long l = Math.round(Double.valueOf(digits).doubleValue()*Math.pow(10,noDec));
      digits = Double.toString(l*Math.pow(10,-noDec));
      if (digits.charAt(2) == '.')  { // two digits to the left of the decimal point
        digits = digits.charAt(0) + "." + digits.charAt(1) + digits.substring(3); 
        int exp = Integer.parseInt(s.substring(s.indexOf("E")+1));
        s = "E" + (exp+1);
      }
    }
    if (noDec > digits.length()-2)
      digits = adjustLeft(digits, noDec+2, '0');
    else if (noDec>0) 
      digits = digits.substring(0, noDec+2);
    else // noDec==0
      digits = digits.substring(0, 1);
    if (num >= 0)
      s = digits + s.substring(s.indexOf("E"));
    else 
      s = '-' + digits + s.substring(s.indexOf("E"));
    s = adjustRight(s, noPos, fill);
    if (r instanceof DecimalFormat) 
      // change decimal point to local decimal symbol (i.e. decimal comma)  
      s = s.replace('.', ((DecimalFormat) r).getDecimalFormatSymbols().getDecimalSeparator());
    return s;
  }

  private static String formatNumber(float num, int noPos, int noDecMax, 
                               int noDecMin, char fill, NumberFormat r) {
    if (useExponent(num, noDecMax))
      return formatExponent((double) num, noPos, 6, fill, r);
    else
      return formatFixed((double) num, noPos, 
                         noDecMax, noDecMin, fill, r); 
  }

  private static String formatNumber(double num, int noPos, int noDecMax, 
                              int noDecMin, char fill, NumberFormat r) {
    if (useExponent(num, noDecMax))
      return formatExponent(num, noPos, 14, fill, r);
    else
      return formatFixed(num, noPos, noDecMax, noDecMin, fill, r);
  }
 

 // public class methods

  public static ExtendedWriter getFileWriter(String fileName) {
    try {
      return new ExtendedWriter(new BufferedWriter(
                            new FileWriter(fileName)), true);
    }
    catch (IOException e) {
      return null;
    }
  }
    
  public static ExtendedWriter getFileWriter(String fileName, boolean append) {
    try {
      return new ExtendedWriter(new BufferedWriter(
                            new FileWriter(fileName, append)), true);
    }
    catch (IOException e) {
      return null;
    }
  }

  public static void setStaticFormat(NumberFormat f) {
    snf = (NumberFormat) f.clone();
    snf.setGroupingUsed(false);
  }
    
  public static NumberFormat getStaticFormat() {
    return snf;
  } 
	
  public static String adjustRight(String s, int noPos, char fill) {
   if (s.length() < noPos) {
     char[] blanks = new char[noPos-s.length()];
     for (int i=0; i < blanks.length; i++)
       blanks[i] = fill;
     return new String(blanks) + s;
   }
   else
     return s;
  }

  public static String adjustLeft(String s, int noPos, char fill) {
    if (s.length() < noPos) {
      char[] blanks = new char[noPos-s.length()];
      for (int i=0; i<blanks.length; i++)
        blanks[i] = fill;
      return s + new String(blanks);
    }
    else
      return s;
  }

  public static String toFixedLength(String s, int length) {
    return adjustLeft(s, length, ' ').substring(0,length);
  }

  public static boolean useExponent(double num, int noDec) {
    // The method format in java.text.NumberFormat cannot
    // format very small or very large numbers.
    // This method tests if the number must be printed in exponent form
    return Math.abs(num*Math.pow(10,noDec)) >= Long.MAX_VALUE;
  }

  public static String formatNum(long num) {
    return snf.format(num);
  }

  public static String formatNum(long num, int noPos) {
    return adjustRight(snf.format(num), noPos, ' ');
  }

  public static String formatNum(long num, int noPos, char fill) {
    return adjustRight(snf.format(num), noPos, fill);
  }


  public static String formatNum(double num) {
    return formatNumber(num, 1, 6, 0, ' ', snf);
  }

  public static String formatNum(double num, int noPos) {
    return formatNumber(num, noPos, 6, 0, ' ', snf);
  }

  public static String formatNum(double num, int noPos, char fill) {
    return formatNumber(num, noPos, 6, 0, fill, snf);
  }

  public static String formatNum(double num, int noPos, int noDec) {
    return formatNumber(num, noPos, noDec, noDec, ' ', snf);
  }

  public static String formatNum(double num, int noPos, int noDec, char fill) {
    return formatNumber(num, noPos, noDec, noDec, fill, snf);
  }

  public static String formatExp(double num) {
    return formatExponent(num, 1, 14, ' ', snf);
  }

  public static String formatExp(double num, int noPos) {
    return formatExponent(num, noPos, 14, ' ', snf);
  }

  public static String formatExp(double num, int noPos, char fill) {
    return formatExponent(num, noPos, 14, fill, snf);
  }

  public static String formatExp(double num, int noPos, int noDec) {
    return formatExponent(num, noPos, noDec, ' ', snf);
  }

  public static String formatExp(double num, int noPos, int noDec, char fill) {
    return formatExponent(num, noPos, noDec, fill, snf);
  }

  public static String formatExp(float num) {
    return formatExponent(num, 1, 6, ' ', snf);
  }

  public static String formatExp(float num, int noPos) {
    return formatExponent(num, noPos, 6, ' ', snf);
  }

  public static String formatExp(float num, int noPos, char fill) {
    return formatExponent(num, noPos, 6, fill, snf);
  }


  // instance methods

  public void setFormat(NumberFormat f) {
    nf = (NumberFormat) f.clone();
    nf.setGroupingUsed(false);
  }
    
  public NumberFormat getFormat() {
    return nf;
  }  

  public void setFillChar(char c) {
    fillChar = c;
  }
		 		
  public void print(String s, int noPos) {
    print(adjustLeft(s, noPos, fillChar));
  }

  public void print(int num, int noPos) {
    print((long) num, noPos);
  }

  public void print(int num) {
    print((long) num);
  }
	
  public void print(long num, int noPos) {
    print(adjustRight(nf.format(num), noPos, fillChar));
  }

  public void print(long num) {
    print(nf.format(num));
  }

  public void print(float num, int noPos, int noDec) {
    print(formatNumber(num, noPos, noDec, noDec, fillChar, nf));
  }

  public void print(float num, int noPos) {
    print(formatNumber(num, noPos, 6, 0, fillChar, nf));
  }

  public void print(float num) {
    print(formatNumber(num, 1, 6, 0, fillChar, nf));
  }

  public void print(double num, int noPos, int noDec) {
    print(formatNumber(num, noPos, noDec, noDec, fillChar, nf));
  }

  public void print(double num, int noPos) {
    print(formatNumber(num, noPos, 6, 0, fillChar, nf));
  }

  public void print(double num) {
    print(formatNumber(num, 1, 6, 0, fillChar, nf));
  }
   
  public void printExp(double num) {
    print(formatExponent(num, 1, 14, fillChar, nf));
  }

  public void printExp(double num, int noPos) {
    print(formatExponent(num, noPos, 14, fillChar, nf));
  }

  public void printExp(double num, int noPos, int noDec) {
    print(formatExponent(num, noPos, noDec, fillChar, nf));
  }  
 
  public void printExp(float num) {
    print(formatExponent((double) num, 1, 6, fillChar, nf));
  }

  public void printExp(float num, int noPos) {
    print(formatExponent((double) num, noPos, 6, fillChar, nf));
  } 

  public void printExp(float num, int noPos, int noDec) {
    print(formatExponent((double) num, noPos, noDec, fillChar, nf));
  }  

  public void print(boolean b, int noPos) {
    print(String.valueOf(b), noPos);
  }

  public void print(char c) {
    super.print(c);
  }

  public void print(char c, int noPos) {
    print(adjustRight("" + c, noPos, fillChar));
  }

  public void print(char[] a, int noPos) {
    print(new String(a), noPos);
  }

  public void print(Object obj, int noPos) {
    print(obj.toString(), noPos);
  }

  public void println(String s, int noPos) {
    print(s,noPos); println();  
  }

  public void println(int num, int noPos) {   
    print(num, noPos); println();  
  }

  public void println(int num) {
    print(num); println();
  }
	
  public void println(long num, int noPos) {  
    print(num, noPos); println();  
  }

  public void println(long num) {
    print(num); println();
  }

  public void println(float num, int noPos, int noDec) {
    print(num, noPos, noDec); println();   
  }

  public void println(float num, int noPos) {
    print(num, noPos); println();   
  }

  public void println(float num) { 
    print(num); println();   
  }
	
  public void println(double num, int noPos, int noDec) {
    print(num, noPos, noDec); println();   
  }

  public void println(double num, int noPos) {
    print(num, noPos); println();   
  }

  public void println(double num) { 
    print(num); println();   
  }

  public void printlnExp(float num, int noPos, int noDec) {
    printExp(num, noPos, noDec); println();   
  }

  public void printlnExp(float num, int noPos) {
    printExp(num, noPos); println();   
  }

  public void printlnExp(float num) { 
    printExp(num); println();   
  }
	
  public void printlnExp(double num, int noPos, int noDec) {
    printExp(num, noPos, noDec); println();   
  }

  public void printlnExp(double num, int noPos) {
    printExp(num, noPos); println();   
  }

  public void printlnExp(double num) { 
    printExp(num); println();   
  }

  public void println(boolean b, int noPos) {
    print(b, noPos); println();   
  }

  public void println(char c) {
    super.println(c);
  }

  public void println(char c, int noPos) {
    print(c, noPos); println();
  }

  public void println(char[] a, int noPos) {
    print(a, noPos); println();
  }

  public void println(Object obj, int noPos) {
    print(obj, noPos); println();
  }       
}


 