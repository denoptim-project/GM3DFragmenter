/**
 * Licence
 */

import java.io.IOException;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.File;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.HashMap;

import java.util.Date;

/**
 * Prepare files and scripts to run a list of independent fragmentation jobs
 */

public class ParallelizeGM3DF
{
    //List of options 
    private static Map<String,String> options = new HashMap<String,String>();
    private static Map<String,String> optionsHelp = new HashMap<String,String>();

    //number of processes
    private static int num;

    //Names
    private static String uniqueID;
    private static String scriptName;
    private static String inFile;
    private static List<String> outFiles = new ArrayList<String>();

    //Path of the working directory
    private static String path;

    //newline
    private static String sep = System.getProperty("line.separator");

//------------------------------------------------------------------------------
   public static void main(String[] args) throws Exception
   {
	//Initial message
        System.out.println("\n    Parallelize GM3DFragmenter jobs");
        System.out.println("    ===============================\n");
	
	//Initialize options
	setDefault();

	//Usage message
	if (args.length == 0)
	    fatal();

	//Read options
	for (String s : args)
	{
	    for (String key : options.keySet())
	    {
		
		if (s.startsWith(key))
		{
		    String val = s.substring(key.length());
		    options.put(key,val);
		}
	    }
	}

	//Set all
        Date date = new Date();
        uniqueID = String.valueOf(date.getTime());
	num = Integer.parseInt(options.get("-n"));
	path = options.get("-p");
	if (!path.equals("") && !path.endsWith("/"))
	    path = path+"/";
        inFile = args[0];
	for (int i=0; i<num; i++)
	{
	    String outName = new String();
	    outName = path+"LibFrgJob_"+i+".sdf";
	    outFiles.add(outName);
	}
	String java = options.get("-jp");
	String gm3df = options.get("-GM3DFp");
	String inp = options.get("-i");
	int wt =  Integer.parseInt(options.get("-wt"));
	int nn =  Integer.parseInt(options.get("-a"));
	int np =  Integer.parseInt(options.get("-b"));
	int mem =  Integer.parseInt(options.get("-mem"));
	String sub = options.get("-queue");

//Prepare portions of the library
	//Get number of objects in the SDF file
	int totObj = 0;
	try {
 	    BufferedReader reader = new BufferedReader(new FileReader(inFile));
	    String line;
            System.out.print(" Reading the file...");
	    while ((line = reader.readLine()) != null) {
	        if (line.startsWith("$$$$"))
		    totObj++;
	    }
	    reader.close();
	} catch (Throwable t) {
            System.out.println("ERROR in reading the input file!");
	    t.printStackTrace();
	    System.exit(-1);
	}
	System.out.println(" Found "+totObj+" objects.");

	//Generate smaller files
	int objPerFile = totObj / num;
	System.out.println(" Splitting in "+num+" groups ("+objPerFile+" objects each)");
	try {
            BufferedReader reader = new BufferedReader(new FileReader(inFile));

	    //generate writers
	    List<FileWriter> writers = new ArrayList<FileWriter>();
	    try {
		for (String fname : outFiles)
		{
		    FileWriter w = new FileWriter(fname,true);
		    writers.add(w);
		}
	    } catch (Throwable t2) {
                System.out.println("ERROR in generating writers!");
	 	t2.printStackTrace();
                System.exit(-1);
	    }

            String line;
            int locObj = 0;
	    int loc = 0;
	    System.out.println(" Output: "+outFiles.get(loc));
            while ((line = reader.readLine()) != null) {
		writers.get(loc).write(line+sep);
//		if (line.startsWith("$$$$"))
// 		    System.out.println(line+" loc:"+loc+" locObj:"+locObj);
                if (line.startsWith("$$$$"))
                    locObj++;
		if (locObj == objPerFile) 
		{
		    if ((loc + 1)  < outFiles.size())
		    { 
		        loc++;
			locObj = 0;
			System.out.println(" Output: "+outFiles.get(loc));
			writers.get(loc-1).close();
		    }
		}
            }
            reader.close();
	    writers.get(loc).close();
        } catch (Throwable t) {
            System.out.println("ERROR in writing output file!");
            t.printStackTrace();
            System.exit(-1);	   
        }
	
//Prepare script
        scriptName = "RunAllGM3DFragmentation_"+uniqueID+".sh";
	writeTXTAppend(scriptName,"#! /bin/bash"+sep,true);
        writeTXTAppend(scriptName,"# Run a list of independent fragmentation jobs"+sep,true);
        writeTXTAppend(scriptName,"cd "+path+sep,true);

	for (int i=0; i<num; i++)
	{
	    String jobfile = path+"RunFrgJob_"+i+".sh";
	    String locLog = path+"logFrgJob_"+i;
            writeTXTAppend(scriptName,"chmod 754 "+jobfile+sep,true);
	    //select submission protocoll
	    if (sub.equals("FIMM"))
	    {
                writeTXTAppend(scriptName,"qsub -A kjemisk -l nodes="+nn+":ppn="+np+",walltime="+wt+":10:00,pmem="+mem+"mb "+jobfile+sep,true);
	    } else if (sub.equals("187")) {
		if (path.equals(""))
		    writeTXTAppend(scriptName,"./"+jobfile+" "+sep,true);
		else
		    writeTXTAppend(scriptName,jobfile+" "+sep,true);
	    } else if (sub.equals("STALLO")) {
                writeTXTAppend(scriptName,"qsub -A nn2506k -l nodes="+nn+":ppn="+np+",walltime="+wt+":10:00,pmem="+mem+"mb "+jobfile+sep,true);
	    } else {

		System.out.println("\n ERROR! Submission protocol NOT FOUND!");
                System.out.println(" Check option '-queue' = "+sub+"\n");
		System.exit(-1);
	    }
	    //preprare parameters for locale fragmentation job
	    String locInp = path+"InpFrgJob_"+i+".par";
	    updateParametersFile(inp,locInp,outFiles.get(i));
	    //write jobfile
	    writeTXTAppend(jobfile,"#! /bin/bash"+sep,true);
            writeTXTAppend(jobfile,"cd "+path+sep,true);
	    writeTXTAppend(jobfile,java+" -jar "+gm3df+" "+locInp+" > "+locLog+" "+sep,true); 
	}
	
	//Change permissions
	File sc = new File(scriptName);
	sc.setExecutable(true,false);
	

	//Final message
	System.out.println("\n Normal termination! :-)");
        System.out.println("\n Submit all the independent jobs with the command: ");
	System.out.println("\t ./"+scriptName+"\n");

   }

//------------------------------------------------------------------------------

/**
 */

    public static void updateParametersFile(String in, String out, String newStrFile)
    {
        try {
            BufferedReader reader = new BufferedReader(new FileReader(in));
            FileWriter writer = new FileWriter(out,true);
            String line;
            while ((line = reader.readLine()) != null) {
		if (line.startsWith("STRUCTURESFILE"))
		{
		    writer.write("STRUCTURESFILE "+newStrFile+sep);
		} else {
		    writer.write(line+sep);
		}
            }
            reader.close();
            writer.close();
        } catch (FileNotFoundException fnf) {
            System.err.println("\nFile Not Found: " + in);
            System.err.println("Check option '-i'\n");
            System.exit(-1);
        } catch (Throwable t) {
            System.out.println("\nERROR in writing output file!");
            t.printStackTrace();
            System.exit(-1);
        }
    }

//------------------------------------------------------------------------------

/**
 * Write a line on a TXT file
 * @param text the line to be written
 * @param file name of the target file
 */
    public static void writeTXTAppend(String filename, String txt,
                                      boolean append)
    {
        FileWriter writer = null;
        try
        {
            writer = new FileWriter(filename,append);
            writer.write(txt);
        } catch (Throwable t) {
            System.err.println("\nFailure in writing TXT: " + t);
            System.exit(-1);
        } finally {
             try {
                 if(writer != null)
                     writer.close();
             } catch (IOException ioe) {
                 System.err.println("\nError in writing: " + ioe);
                 System.exit(-1);
             }
        }
    }

//------------------------------------------------------------------------------
    private static void setDefault()
    {
	//Default value for the options
	String optKey = "-n";
	String optHlp = "define the number of independet jobs to create";
	String optDef = "2";
	options.put(optKey,optDef);
        optionsHelp.put(optKey,optHlp);

        String optKey1 = "-p";
        String optDef1 = "";
        String optHlp1 = "path to working directory [default: "+optDef1+"]";
        options.put(optKey1,optDef1);
        optionsHelp.put(optKey1,optHlp1);

        String optKey2 = "-jp";
        String optDef2 = "/usr/bin/java";
        String optHlp2 = "path to java executable [default: "+optDef2+"]";
        options.put(optKey2,optDef2);
        optionsHelp.put(optKey2,optHlp2);

        String optKey3 = "-GM3DFp";
        String optDef3 = "./GM3DFragmenter.jar";
        String optHlp3 = "path to GM3DFragmenter.jar [default: "+optDef3+"] ";
        options.put(optKey3,optDef3);
        optionsHelp.put(optKey3,optHlp3);

        String optKey5 = "-a";
        String optDef5 = "1";
        String optHlp5 = "number of nodes per each job (FIMM qsub format) [default:"+optDef5+" ]";
        options.put(optKey5,optDef5);
        optionsHelp.put(optKey5,optHlp5);

        String optKey6 = "-i";
        String optDef6 = "input.par";
        String optHlp6 = "path to parameters' file for GM3DFragmenter [default:"+optDef6+" ]";
        options.put(optKey6,optDef6);
        optionsHelp.put(optKey6,optHlp6);

        String optKey7 = "-b";
        String optDef7 = "2";
        String optHlp7 = "number of processors per node (FIMM qsub format) [default:"+optDef7+" ]";
        options.put(optKey7,optDef7);
        optionsHelp.put(optKey7,optHlp7);

        String optKey8 = "-mem";
        String optDef8 = "500";
        String optHlp8 = "memory requirements (FIMM qsub format) [default:"+optDef8+"mb ]";
        options.put(optKey8, optDef8);
        optionsHelp.put(optKey8, optHlp8);

        String optKeyA = "-wt";
        String optDefA = "1";
        String optHlpA = "walltime (FIMM qsub format) [default:"+optDefA+" h]";
        options.put(optKeyA,optDefA);
        optionsHelp.put(optKeyA,optHlpA);

        String optKey9 = "-queue";
        String optDef9 = "187";
        String optHlp9 = "selects type of submission protocol: 187 or FIMM [default:"+optDef9+" ]";
        options.put(optKey9,optDef9);
        optionsHelp.put(optKey9,optHlp9);

/*
        String optKey = "";
        String optDef = "";
        String optHlp = "[default:"+optDef+" ]";
        options.put(optKey,optDef);
        optionsHelp.put(optKey,optHlp);
*/

    }
//------------------------------------------------------------------------------
    private static void fatal()
    {
            printUsage();
            System.exit(0);
    }
//------------------------------------------------------------------------------
    private static void fatal(String msg)
    {
            System.err.println("ERROR! "+msg);
	    fatal();
    }
//------------------------------------------------------------------------------
    private static void printUsage()
    {
            System.out.println(" Usage:\n\t java -jar ParallelizeGM3DF.jar <libname> <parameters> ");
            System.out.println(" where:\n\t <libname> = name of the file containing the molecules to chop");
            System.out.println("\t <parameters> = ");
            for (String key : options.keySet())
                System.out.println("\t "+key+" => "+optionsHelp.get(key));
            System.out.println("\n ");
            System.out.println(" Example:\n java -jar ParallelizeGM3DF.jar library.sdf -n4\n");
    }
}
