package general;
import java.io.*;
import java.text.*;
import java.util.*;

public class Std {
  // class variables
  public static ExtendedReader in;
  public static ExtendedWriter out;
  public static ExtendedWriter err;
  private static boolean dos;

  // initialization of the streams in, out and err
  private static void initStreams() {
    String code;
    if (dos) {
      String lang = System.getProperty("user.language"); 
      if (lang.equals("en"))
        code = "Cp437";  // English
      else if (lang.equals("ru"))
        code = "Cp866";  // Russian
      else if (lang.equals("pt"))
        code = "Cp860";  // Portuguese
      else if (lang.equals("ar"))
        code = "Cp864";  // Arabic
      else if (lang.equals("tr"))
        code = "Cp857";  // Turkish
      else if (lang.equals("el"))
        code = "Cp737";  // Greek
      else if (lang.equals("he"))
        code = "Cp862";  // Hebrew
      else if (lang.equals("is"))
        code = "Cp861";  // Icelandic
      else if (lang.equals("pl") || lang.equals("cz"))
        code = "Cp852";  // Slavic languages 
      else if (lang.equals("da") || lang.equals("no"))
        code = "Cp865";  // Danish or Norwegian
      else
        code = "Cp850";   // multilingual (Swedish)      
    }
    else // not dos
      code = System.getProperty("file.encoding");  // default-value
    try {
      in  = new ExtendedReader(new InputStreamReader (System.in,  code), 1);
      out = new ExtendedWriter(new OutputStreamWriter(System.out, code), true); 
      err = new ExtendedWriter(new OutputStreamWriter(System.err, code), true);
    }
    catch (UnsupportedEncodingException e) {
      in  = new ExtendedReader(new InputStreamReader (System.in),  1);
      out = new ExtendedWriter(new OutputStreamWriter(System.out), true); 
      err = new ExtendedWriter(new OutputStreamWriter(System.err), true);
      err.println("Std: Warning. The code " + code + " is not supported"); 
    }
  }

  // initialization of static variables
  static {
    dos = System.getProperty("os.name").startsWith("Windows");
    initStreams();     
  }

  public static void setDos(boolean b) {
    if (b != dos) {
      dos = b;
      initStreams();
    }
  }

  public static boolean getDos() {
    return dos;
  }

  public  static void reinit() {  
  // should be called after a change of default Locale
    in.setLang(Locale.getDefault().getLanguage());
    in.setFormat(NumberFormat.getInstance());
    out.setFormat(NumberFormat.getInstance());
    err.setFormat(NumberFormat.getInstance());
    ExtendedWriter.setStaticFormat(NumberFormat.getInstance());
  }       
 
  // corrections of erronous Swedish Locale
  public static synchronized Collator getCollator(Locale loc) {
    if (loc.getLanguage().equals("sv")) {
      // Create an en_US Collator object
      RuleBasedCollator col = (RuleBasedCollator)
                  Collator.getInstance(new Locale("en", "US", ""));
      // Create a new Collator object with additional Swedish rules
      String addRules = "& Z < Œ, < Š,€ < š,…";
      try { 
        return new RuleBasedCollator(col.getRules()+addRules);
      }
      catch (ParseException e) {return Collator.getInstance(loc);}
    }        
    else 
      return Collator.getInstance(loc);
  } 

  public static synchronized Collator getCollator() {
    return getCollator(Locale.getDefault());
  }

}