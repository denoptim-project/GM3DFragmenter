
import java.io.File;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.HashMap;

import java.io.IOException;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.File;

import java.util.Date;

public class SplitAndMerge
{
    //List of options 
    private static Map<String,String> options = new HashMap<String,String>();
    private static Map<String,String> optionsHelp = new HashMap<String,String>();

    private static String uniqueID;
    private static String path;
    private static String java;
    private static String gm3df;
    private static String sep = System.getProperty("line.separator");
    private static String scriptName;
    private static String sub;
    private static int binSize = 5;

    //Identify input libraries
    private static String headSplit = "Fragments_";
    private static String headMerge = "MWBin_";

//------------------------------------------------------------------------------
    public static void main(String[] args) throws Exception
    {

	//initialize options
	setDefault();

        if (args.length < 1)
	{
	    printUsage();
            System.exit(-1);
	} 

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

        //Set all parameters
        Date date = new Date();
        uniqueID = String.valueOf(date.getTime());	
	path = options.get("-p");
	File folder = new File(path);
	if (!folder.exists())
	{
                System.err.println("\n ERROR! Working directory does not exist!");
                System.err.println(" Check option '-p' = "+path+"\n");
                System.exit(-1);
	}

	sep = System.getProperty("line.separator");
	java = options.get("-jp");
        gm3df = options.get("-GM3DFp");
	headSplit = options.get("-sws");
	headMerge = options.get("-swm");
	sub = options.get("-queue");

        if (args[0].equals("split"))
 	    prepareScriptForSplitting();
        else if (args[0].equals("merge"))
            prepareScriptForMerging();
        else {  
            System.out.println("\n ERROR! Type of job not specified. I got "+args[0]+" instead of 'split/merge'");
            printUsage();
            System.exit(-1);
        }

        //Goodbye
        System.out.println("\n Normal termination! :-)");
        System.out.println("\n Submit all the independent jobs with the command: ");
        System.out.println("\t ./"+scriptName+"\n");

    }

//------------------------------------------------------------------------------

    private static void prepareScriptForSplitting()
    {
        //Get list of Fragments Libraries to be splitted
        File folder = new File(path);
        File[] listOfFiles = folder.listFiles();
        List<File> libFiles = new ArrayList<File>();
        for (int i=0; i<listOfFiles.length; i++)
        {
            File f = listOfFiles[i];
            if (f.getName().startsWith(headSplit))
                libFiles.add(f);
        }

        int wt =  Integer.parseInt(options.get("-wt"));
        int nn =  Integer.parseInt(options.get("-nn"));
        int np =  Integer.parseInt(options.get("-np"));
        int mem =  Integer.parseInt(options.get("-mem"));

	//Prepare sciptForSplitting
	scriptName = "RunAllSPLITTING_"+uniqueID+".sh";
        writeTXTAppend(scriptName,"#! /bin/bash"+sep,true);
        writeTXTAppend(scriptName,"#Run all jobs splitting Libraries (frags MW range)"+sep,true);
        writeTXTAppend(scriptName,"cd "+path+sep,true);
        for (File lib : libFiles)
        {
	    //Get reference of fragmentation job
	    String name = lib.getName();
	    int init = headSplit.length();
	    int end = name.lastIndexOf(".sdf");
	    String frgRef = name.substring(init,end);
	    String[] p = frgRef.split("_");
	    String libRef = p[1]+"_"+p[2];
            String jobName = "RunSpltJob_"+libRef+".sh";
            String GM3DparFile = "InpSpltJob_"+libRef+".par";
	    String logName = "logSpltJob_"+libRef;
	    if (!(path.equals("") || path.endsWith("/")))
	    {
		jobName = path+"/"+jobName;
		logName = path+"/"+logName;
		GM3DparFile = path+"/"+GM3DparFile;
	    }

            writeTXTAppend(scriptName,"echo \"NOFRAGMENTATION\" > "+GM3DparFile+sep,true);
            writeTXTAppend(scriptName,"echo \"REPORT 1\" >> "+GM3DparFile+sep,true);
            writeTXTAppend(scriptName,"echo \"MWSPLITTING "+binSize+"\" >> "+GM3DparFile+sep,true);
            writeTXTAppend(scriptName,"echo \"STRUCTURESFILE "+lib+"\" >> "+GM3DparFile+sep,true);
	    //Prepare job cript
            writeTXTAppend(jobName,"#! /bin/bash"+sep,true);
            writeTXTAppend(jobName,"# Split library in small range according to MW"+sep,true);
            writeTXTAppend(jobName,"cd "+path+sep,true);
            writeTXTAppend(jobName,java+" -jar "+gm3df+" "+GM3DparFile+" > "+logName+sep,true);
            //add to script submitting all jobs
            writeTXTAppend(scriptName,"chmod 754 "+jobName+sep,true);
            if (sub.equals("FIMM"))
            {
                writeTXTAppend(scriptName,"qsub -A kjemisk -l nodes="+nn+":ppn="+np+",walltime="+wt+":10:00,pmem="+mem+"mb "+jobName+sep,true);
            } else if (sub.equals("187")) {
                if (path.equals(""))
                    writeTXTAppend(scriptName,"./"+jobName+" "+sep,true);
                else
                    writeTXTAppend(scriptName,jobName+" "+sep,true);
            } else {
                System.err.println("\n ERROR! Submission protocol NOT FOUND!");
                System.err.println(" Check option '-queue' = "+sub+"\n");
                System.exit(-1);
            }
	}

        //Change permissions
        File sc = new File(scriptName);
        sc.setExecutable(true,false);

    }

//------------------------------------------------------------------------------

    private static void prepareScriptForMerging()
    {
	//Get list of Fragments Libraries to be merged
	File folder = new File(path);
	File[] listOfFiles = folder.listFiles();
	List<File> libFiles = new ArrayList<File>();
	for (int i=0; i<listOfFiles.length; i++) 
  	{
	    File f = listOfFiles[i];
	    String fname =f.getName();
	    if (fname.startsWith(headMerge) && !fname.contains("AllFrg"))
		libFiles.add(f);
	}

	System.out.println("\n Found "+libFiles.size()+" sub libraries of fragments");

	//Get files referring to same MW range
	Map<String,ArrayList<File>> mwRanges = new HashMap<String,ArrayList<File>>();
        List<String> jobNames = new ArrayList<String>();
        for (int i=0; i<libFiles.size(); i++)
        {
            String[] parts = libFiles.get(i).getName().split("_");
            //Store bin name
            String bin = parts[1];
            if (mwRanges.keySet().contains(bin))
            {
                mwRanges.get(bin).add(libFiles.get(i));
            } else {
                ArrayList<File> l = new ArrayList<File>();
                l.add(libFiles.get(i));
                mwRanges.put(bin,l);
            }
	    //Store job name
	    String jobName = "";
	    String jobId = parts[parts.length-1];
            if (jobId.contains("."))
                jobId = jobId.substring(0,jobId.indexOf("."));
	    jobName = parts[parts.length-2]+"_"+jobId;
	    if (!jobNames.contains(jobName))
                jobNames.add(jobName);  
        }

	//Check for files with isomer count
        for (String jbn : jobNames)
	{
	    String ismFileName = "IsomerCountsJob_"+jbn+".dat";
	    File isomerFile = new File(ismFileName);
	    if (!isomerFile.exists())
	    {
		System.out.println("\n ERROR! File "+ismFileName+" NOT FOUND!\n");
		System.exit(0);
	    }
	}

        int wt =  Integer.parseInt(options.get("-wt"));
        int nn =  Integer.parseInt(options.get("-nn"));
        int np =  Integer.parseInt(options.get("-np"));
        int mem =  Integer.parseInt(options.get("-mem"));

	//Prepare scripts to merge all portions of the same range
	scriptName = "RunAllMERGING_"+uniqueID+".sh";
	writeTXTAppend(scriptName,"#! /bin/bash"+sep,true);
        writeTXTAppend(scriptName,"#Run all jobs merging MW-Libraries"+sep,true);
        writeTXTAppend(scriptName,"cd "+path+sep,true);

        for (String bin : mwRanges.keySet())
	{
	    //Prepare names
	    String GM3DparFile = "InpMrgBin_"+bin+".par";
	    String jobName = "RunMrgBin_"+bin+".sh";
	    String logName = "logMrgBin_"+bin;
            if (!(path.equals("") || path.endsWith("/")))
	    {
                GM3DparFile = path+"/"+GM3DparFile;
		jobName = path+"/"+jobName;
		logName = path+"/"+logName;
 	    }
	    
	    //Prepare scripts
	    writeTXTAppend(jobName,"#! /bin/bash"+sep,true);
	    writeTXTAppend(jobName,"# Merge fragment libraries deriving from MW-based splitting"+sep,true);
	    writeTXTAppend(jobName,"cd "+path+sep,true);
	    writeTXTAppend(scriptName,"echo \"NOFRAGMENTATION\" > "+GM3DparFile+sep,true);
            writeTXTAppend(scriptName,"echo \"MWREORDER\" >> "+GM3DparFile+sep,true);
            writeTXTAppend(scriptName,"echo \"REPORT 1\" >> "+GM3DparFile+sep,true);
	    writeTXTAppend(scriptName,"echo \"MWMERGE ",true);
	    for (File f : mwRanges.get(bin))
	        writeTXTAppend(scriptName," "+f,true);
            writeTXTAppend(scriptName,"\" >> "+GM3DparFile+sep,true);
	    //submit job
	    writeTXTAppend(jobName,java+" -jar "+gm3df+" "+GM3DparFile+" > "+logName+sep,true);

	    //add to script submitting all jobs
            writeTXTAppend(scriptName,"chmod 754 "+jobName+sep,true);
            if (sub.equals("FIMM"))
            {
                writeTXTAppend(scriptName,"qsub -A kjemisk -l nodes="+nn+":ppn="+np+",walltime="+wt+":10:00,pmem="+mem+"mb "+jobName+sep,true);
            } else if (sub.equals("187")) {
                if (path.equals(""))
                    writeTXTAppend(scriptName,"./"+jobName+" "+sep,true);
                else
                    writeTXTAppend(scriptName,jobName+" "+sep,true);
            } else {
                System.out.println("\n ERROR! Submission protocol NOT FOUND!");
                System.out.println(" Check option '-queue' = "+sub+"\n");
                System.exit(-1);
            }
	}

        //Change permissions
        File sc = new File(scriptName);
        sc.setExecutable(true,false);

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
            System.err.println("Failure in writing TXT: " + t);
            System.exit(-1);
        } finally {
             try {
                 if(writer != null)
                     writer.close();
             } catch (IOException ioe) {
                 System.err.println("Error in writing: " + ioe);
                 System.exit(-1);
             }
        }
    }

//------------------------------------------------------------------------------
    private static void setDefault()
    {
        //Default value for the options
        String optKey = "-p";
        String optHlp = "path to working directory";
        String optDef = ".";
        options.put(optKey,optDef);
        optionsHelp.put(optKey,optHlp);

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

        String optKey4 = "-wt";
        String optDef4 = "2";
        String optHlp4 = "walltime for each job [default:"+optDef4+" ]";
        options.put(optKey4,optDef4);
        optionsHelp.put(optKey4,optHlp4);

        String optKey5 = "-sws";
        String optDef5 = "Fragments_";
        String optHlp5 = "first characters of fragment libraries [default:"+optDef5+" ]";
        options.put(optKey5,optDef5);
        optionsHelp.put(optKey5,optHlp5);

        String optKey6 = "-nn";
        String optDef6 = "1";
        String optHlp6 = "number of nodes per each job (FIMM qsub format) [default:"+optDef6+" ]";
        options.put(optKey6,optDef6);
        optionsHelp.put(optKey6,optHlp6);

        String optKey7 = "-np";
        String optDef7 = "2";
        String optHlp7 = "number of processors per node (FIMM qsub format) [default:"+optDef7+" ]";
        options.put(optKey7,optDef7);
        optionsHelp.put(optKey7,optHlp7);

        String optKey8 = "-mem";
        String optDef8 = "500";
        String optHlp8 = "memory requirements (FIMM qsub format) [default:"+optDef8+"mb ]";
        options.put(optKey8, optDef8);
        optionsHelp.put(optKey8, optHlp8);


        String optKey9 = "-queue";
        String optDef9 = "187";
        String optHlp9 = "selects type of submission protocol: 187 or FIMM [default:"+optDef9+" ]";
        options.put(optKey9,optDef9);
        optionsHelp.put(optKey9,optHlp9);

        String optKeyA = "-swm";
        String optDefA = "MWBin_";
        String optHlpA = "first characters of MW range libraries to merge [default:"+optDefA+" ]";
        options.put(optKeyA,optDefA);
        optionsHelp.put(optKeyA,optHlpA);


/*
        String optKey = "";
        String optDef = "";
        String optHlp = "[default:"+optDef+" ]";
        options.put(optKey,optDef);
        optionsHelp.put(optKey,optHlp);
*/

    }

//------------------------------------------------------------------------------
    private static void printUsage()
    {
            System.out.println("\n Usage: \n java -jar SplitAndMerge.jar <jobType> <options>\n");
            System.out.println(" Types of job:");
            System.out.println(" -> split "+
            "\n  Search for all SDF files starting with \""+headSplit+"\" and assumes that"+
            "\n  these are output from parallel fragmentation jobs designed by 'ParallelizeGM3DF'. "+
            "\n  These libraries are then splitted in sub-libraries containing only fragments with"+
            "\n  a specified range of molecular weight (BIN).");
            System.out.println(" -> merge "+
            "\n  Search for all SDF files starting with \""+headMerge+"<Num>\" and produces a library by "+
            "\n  including all the unique fragments found in the list of \"SubLib<Num>\" referring"+
            "\n  to the same BIN <Num> (range of molecular weight).");
//TODO avoid usage of the shell script or run them directly

            System.out.println("\n List of options:");
            for (String key : options.keySet())
                System.out.println(key+"\t=> \t"+optionsHelp.get(key));
            System.out.println("\n ");
    }

//------------------------------------------------------------------------------
}
