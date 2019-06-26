/*    
 *   GM3DFragmenter
 *   Copyright (C) 2019 Marco Foscato <marco.foscato@uib.no>
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Affero General Public License as published
 *   by the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Affero General Public License for more details.
 *
 *   You should have received a copy of the GNU Affero General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import java.util.Set;
import java.util.HashSet;
import java.util.Map;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.SortedMap;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Date;
import java.io.File;
import java.text.DateFormat;
import java.text.SimpleDateFormat;

import org.openscience.cdk.interfaces.IChemObject;

/**
 * Setting and storage of parameters and user specified options 
 * 
 * @author Marco Foscato (University of Bergen)
 */

public class Parameters
{
//Job details
    private static String inputName;
    private static String jobName;
    private static String dateForm;

//FLAGS defining the type of jobb
    //Flag: pre-filtering of molecules
    public static boolean preFiltering;

    //Flag: single output file that can be used directly by other routines
    public static boolean singleOutput;

    //Flag: add dummy atoms to avoid vanishing torsions in internal coordinates
    public static boolean addDummyOnLinear;

    //Flag: add index to APclass name - this is to make sure each pair of
    //      APs as an unique APclass and allow re-building the same molecule
    //      by combination of its own fragments
    public static boolean addIDtoAPClass;
    //Flag: read the cutting rules
    public static boolean readRules;

    //Flag: check agreement between formula and SD structure
    public static boolean chkFormula;

    //Flag: analyze and fix input structures
    public static boolean analyzeAndFix;

    //Flag: filtering of a library
    public static boolean onlyFiltering;

    //Flag: fragmentation
    public static boolean chopMols;

    //Flag: remove duplicates
    public static boolean rmDuplicates;

    //Flag: regroup rotamers
    public static boolean groupingRotamers;

    //Flag: merge with other library
    public static boolean mergeLibraries;

    //Flag: change format of Fragment library
    public static boolean covertLibrary;

    //Flag: reorder according to MW 
    public static boolean orderMW;
    public static boolean MWascending;

    //Flag: split library according to fragmenr MW
    public static boolean MWsplitting;
    public static int massBinSize = 2;

    //Flag: merge libraries  deriving from MW splitting
    public static boolean MWMerge;

    //Flag: delete intermediate files
    public static boolean removeIntermFiles;

    //Extraction of fragments from libraries
    public static boolean extractClass;
    public static boolean extractSMARTS;

    //Flag: just create CPMatrix from cutting rules
    public static boolean makeCPMap;

    //Flag: ignore fragments that are already known
    public static boolean ignoreKnownFrags;

    //Flag: collect only fragments that match one of the targets
    public static boolean lookForTargets;

    //Amount of information printed on screen 
    //(0=minimum, 1=intermediate,  2=maximum, 3=only_for_debug)
    public static int report = 0;

//FILE NAMES
    //KeyFile - User defined parameters
    public static String input;

    //SDF input file contining molecules to be processed
    private static String sdfInFile;

    //Updated pointer to SDFfile
    private static String upToDateSDFile;

    //SDF intermediate  files
    private static Set<String> interFiles = new HashSet<String>();

    //TXT input file containing chemical information for CCDC entries
    public static String txtFile;

    //TXT input file defining the cutting rules
    public static String rulesFile;

    //Existing SDF file containing GM3DFragments
    public static String previousLibrary;
    public static String uniqueLibrary;

    //SDF with known fragments
    public static String ignorableFragsFile;
    public static String ignorableFragsFileFormat;

    //SDF with fragments to be matched during the fragmentation
    public static String targetFragsFile;
    public static String targetFragsFileFormat;

    //List of MWLibraries to merge and Isomer counts
    public static Set<String> MWMpath = new HashSet<String>();

    // pre-filtering of molecules with SMARTS
    public static Set<String> preFilterSMARTS = new HashSet<String>();

//COMMON STUFF STORAGE
    //Storage for cutting rules
    public static Map<String,GM3DCuttingRule> rules = new HashMap<String,GM3DCuttingRule>();
    public static SortedMap<Integer,String> sortedRules = new TreeMap();
    public static List<String> anyAtm = new ArrayList<String>();

    //Storage for found CLASSes (<rule_name>+<SubClass_flag>)
    public static ArrayList<String> classes = new ArrayList<String>();

//FRAGMENTATION
    //Maximum number of classes
    public static int maxNumClasses = 1000;

    //CLASS string format
    public static String moreAtmsSeparator = " ";
    public static String moreAPSeparator = ",";
    public static String atmSeparator = "#";
    public static String detailsSeparator = ":";
    public static String coordSeparator = "%";

    //Rejection rules. Which fragment will be rejected?
    // -> fragments containing AP with these CLASS+SubCLASS
    public static Set<String> rejClasses = new HashSet<String>();
    // -> fragments containing any combination of APs with these CLASS or portion of CLASS' name
    // NB: only the combinations are detected. Fragments containing one AP of one of these these classes but not combined with nothet one are not rejected
    public static Set<List<String>> rejClassCombination = new HashSet<List<String>>();
    // -> fragments matching these SMARTS
    public static Set<String> rejSMARTS = new HashSet<String>();
    // -> fragments NOT matching these SMARTS
    public static Set<String> rejNotSMARTS = new HashSet<String>();
    // -> fragments containing these elements
    public static Set<String> rejElements = new HashSet<String>();
    // -> fragments containing more than the number of atoms specified for each element
    public static Set<List<String>> rejFormulaMax = new HashSet<List<String>>();
    // -> fragments containing less than the number of atoms specified for each elements
    public static Set<List<String>> rejFormulaMin = new HashSet<List<String>>(); 
    // -> fragments havimg more/less atoms than this limit
    public static int fragmentMaxSize = 30;
    public static int fragmentMinSize = 0;
    // -> fragments containing uncommon isotopes
    public static boolean rejIsotopes = true;

//FORMAT of FRAGMENT LIBRARIES
    private static String libFormat = "DENOPTIM";
    private static String oldFormat = "DENOPTIM"; //for reading
    private static String newFormat = "DENOPTIM"; //for writing

//PERIODIC TABLE
/* 
   [alkali] 
"Li", "Na", "K", "Rb", "Cs", "Fr"
   [alkaline earth] 
"Be", "Mg", "Ca", "Sr", "Ba", "Ra"
   [poor] 
"Al", "Ga", "Ge", "In", "Sn", "Sb", "Tl", "Pb", "Bi", "Po"
   [TM] 
"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Ac", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"
   [1st row] 
"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"
   [2nd row]
"Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd"
   [3rd row]
"Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"
   [Lanthanides]
"La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"
   [Actinides]
"Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"
*/
   public static ArrayList<String> metals = new ArrayList<String>(
	Arrays.asList("Li", "Na", "K", "Rb", "Cs", "Fr",
	"Be", "Mg", "Ca", "Sr", "Ba", "Ra",
	"Al", "Ga", "Ge", "In", "Sn", "Sb", "Tl", "Pb", "Bi", "Po",
	"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
	"Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
	"Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
	"La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
	"Tm", "Yb", "Lu",
	"Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", 
	"Md", "No", "Lr"));

//DUMMY ATOM SYMBOL
   public static final String duSymbol = "Du";

//CLOSE-TO-LINEAR BEND THRESHOLD
   public static double linearBendThld = 170.0;

//GM3DF parameters for CDK's SMARTSQueryTool
   public static int maxRingSizeMF = 7;
   public static long maxTimeAllRingFinder = 10000;


//------------------------------------------------------------------------------

    public static void setAll(String[] CLOpts)
    {
        //Set language settings
        Locale.setDefault(Locale.ENGLISH);

        //Set unique date-and-time dependent string
        Date date = new Date();
        String dts = String.valueOf(date.getTime());

        DateFormat df = new SimpleDateFormat("dd/MM/yyyy");
        dateForm = df.format(new Date());

	preFiltering = false;
        addDummyOnLinear = false;
        addIDtoAPClass = false;
        singleOutput = true;
        inputName = "input";
        input = inputName+".par";
        chkFormula = false;
        chopMols = true;
	analyzeAndFix = false;
        rmDuplicates = true;
	lookForTargets = false;
        groupingRotamers = false;
        orderMW = false;
        mergeLibraries = false;
        covertLibrary = false;        
        MWascending = true;
        MWsplitting = false;
        MWMerge = false;
        extractClass = false;
        extractSMARTS = false;
        removeIntermFiles = false;
        sdfInFile = "structures.sdf";        
        rulesFile = "cuttingrules.rul";
        readRules = true;
	makeCPMap = false;
	onlyFiltering = false;
        report = 0;

        //Read User's params from input file
        if (CLOpts.length > 0)
        {
            input = CLOpts[0];
            File inputFile = new File(input);
            if (!inputFile.exists())
            {
                System.err.println("ERROR! Input file "+input+" not found!");
                System.exit(0);
            }
            inputName = inputFile.getName();
            if (inputName.contains("."))
                inputName = inputName.substring(0,inputName.lastIndexOf("."));
        }
        readParametersFile(input);

        //Setting standard name of all the files generated by this job
        jobName = inputName + "_" + dts;
        upToDateSDFile = sdfInFile;

        // Import cutting rules
        if (readRules)
            setCuttingRules();

        //Reporting
        if (report >=1)
            printParameters();

        if (chopMols)
            if (report >= 1)
                printCuttingRules();
        
    }
//------------------------------------------------------------------------------

    public static String getJobName()
    {
        return jobName;
    }

//------------------------------------------------------------------------------

    public static void readParametersFile(String filename)
    {
        ArrayList<String> lines = new ArrayList<String>();
        lines = IOtools.readTXT(filename);
        for (int i = 0; i< lines.size(); i++)
        {
            String line = lines.get(i);
            if (line.startsWith("#"))
                continue;
            String[] words = line.split("\\s+");

            if (words[0].equals("CHECKFORMULA"))
            {
                chkFormula = true;
            } else if (words[0].equals("ANALYSEANDFIX"))
            {
                analyzeAndFix = true;
            } else if (words[0].equals("NOFRAGMENTATION"))
            {
                chopMols = false;
            } else if (words[0].equals("KEEPDUPLICATES"))
            {
                rmDuplicates = false;
            } else if (words[0].equals("IGNOREKNOWNFRAGS"))
            {
                ignoreKnownFrags = true;
                if (words.length != 3)
                {
                    killDueToParams("The keyword 'IGNOREKNOWNFRAGS' requires 2 arguments: <filename> <fragment's format>. Check the input!");
                }
                ignorableFragsFile = words[1];
                ignorableFragsFileFormat = words[2];
            } else if (words[0].equals("LOOKFORTARGETFRAGS"))
            {
		lookForTargets = true;
                rmDuplicates = false;
		if (words.length != 3)
		{
		    killDueToParams("The keyword 'LOOKFORTARGETFRAGS' requires 2 arguments: <filename> <fragment's format>. Check the input!");
		}
		targetFragsFile = words[1];
		targetFragsFileFormat = words[2];
            } else if (words[0].equals("MWREORDER"))
            {
                orderMW = true;
            } else if (words[0].equals("GROUPROTAMERS"))
            {
                groupingRotamers = true;
            } else if (words[0].equals("ONLYFILTER"))
            {
                onlyFiltering = true;
		chopMols = false;
		readRules = false;
                singleOutput = false;
            } else if (words[0].equals("MWSPLITTING"))
            {
//TODO make it working for the fragmentation and all other steps dealing with libraries
//                killDueToParams("MWSPLITTING NOT IMPLEMENTED YET");
                MWsplitting = true;
                readRules = false;
                singleOutput = false;
                if (words.length == 2)
                    massBinSize = Integer.parseInt(words[1]);
            } else if (words[0].equals("MWMERGE"))
            {
//TODO make it as a continuation of a MW-splitting job
                MWMerge = true;
                readRules = false;
                for (int j=0; j<words.length; j++)
                {
                    if (j == 0)
                        continue;
                    File subLibFile = new File(words[j]);
                    if (!subLibFile.exists())
                        killDueToParams("Library "+subLibFile+" NOT FOUND!");
                    MWMpath.add(words[j]);
                }

            } else if (words[0].equals("MERGELIBRARY"))
            {
                mergeLibraries = true;
                previousLibrary = words[1];
                if (words.length > 2)
                    uniqueLibrary = words[2];
            } else if (words[0].equals("CONVERTLIBRARY"))
            {
                covertLibrary = true;
                if (words.length < 3)
                    killDueToParams("Missing library format for conversion");
                else {
                    oldFormat = words[2];
                    newFormat = words[3];
                }
            } else if (words[0].equals("REVERSEMWORDER"))
            {
                MWascending = false;
            } else if (words[0].equals("REMOVEINTERFILES"))
            {
                removeIntermFiles = true;
            } else if (words[0].equals("FORMULATXTFILE"))
            {
                txtFile = words[1];
            } else if (words[0].equals("STRUCTURESFILE"))
            {
                sdfInFile = words[1];
            } else if (words[0].equals("REPORT"))
            {
                report = Integer.parseInt(words[1]);
            } else if (words[0].equals("RULESFILE"))
            {
                rulesFile = words[1];
            } else if (words[0].equals("REJECTCLASS"))
            {
                rejClasses.add(words[1]);
            } else if (words[0].equals("REJECTCLASSCOMBINATION"))
            {
		List<String> pair = new ArrayList<String>();
		if (words.length < 3)
		{
                    killDueToParams("Keyword 'REJECTCLASSCOMBINATION' requires 2 strings as argument");
		}
		pair.add(words[1]);
		pair.add(words[2]);
		rejClassCombination.add(pair);
            } else if (words[0].equals("MAXFRAGSIZE"))
            {
                fragmentMaxSize = Integer.parseInt(words[1]);
            } else if (words[0].equals("MINFRAGSIZE"))
            {
                fragmentMinSize = Integer.parseInt(words[1]);
            } else if (words[0].equals("REJECTSMARTS"))
            {
                rejSMARTS.add(words[1]);
            } else if (words[0].equals("PREFILTERSMARTS"))
            {
		preFiltering = true;
		preFilterSMARTS.add(words[1]);
            } else if (words[0].equals("REJECTNOTSMARTS"))
            {
		rejNotSMARTS.add(words[1]);
            } else if (words[0].equals("REJELEMENT"))
            {
                rejElements.add(words[1]);
            } else if (words[0].equals("REJFORMULAMORETHAN"))
            {
		List<String> formulaLim = new ArrayList<String>();
		for (int j=1; j<words.length; j++)
                {
		    formulaLim.add(words[j]);
		}
                rejFormulaMax.add(formulaLim);
            } else if (words[0].equals("REJFORMULALESSTHAN"))
            {
                List<String> formulaLim = new ArrayList<String>();
                for (int j=1; j<words.length; j++)
                {
                    formulaLim.add(words[j]);
                }
                rejFormulaMin.add(formulaLim);

            } else if (words[0].equals("KEEPISOTOPES"))
            {
                rejIsotopes = false;
            } else if (words[0].equals("LINEARBONDTHLD"))
            {
                linearBendThld = Double.parseDouble(words[1]);
            } else if (words[0].equals("ADDDUMMYONLINEAR"))
            {
                addDummyOnLinear = true;
            } else if (words[0].equals("EXTRACTCLASS"))
            {
                txtFile = words[1];
                extractClass = true;

                readRules = false;
                chkFormula = false;
                chopMols = false;
                rmDuplicates = false;
                groupingRotamers = false;
                mergeLibraries = false;
                covertLibrary = false;
                orderMW = false;
                MWascending = false;
                MWsplitting = false;
                MWMerge = false;
                removeIntermFiles = false;
                extractSMARTS = false;
		singleOutput = false;

            } else if (words[0].equals("EXTRACTSMARTS"))
            {
                txtFile = words[1];
                extractSMARTS = true;

                readRules = false;
                chkFormula = false;
                chopMols = false;
                rmDuplicates = false;
                groupingRotamers = false;
                mergeLibraries = false;
                covertLibrary = false;
                orderMW = false;
                MWascending = false;
                MWsplitting = false;
                MWMerge = false; 
                removeIntermFiles = false;
                extractClass = false;
		singleOutput = false;
            } else if (words[0].equals("MAKECPMFROMRULES"))
            {
		makeCPMap = true;
            } else if (words[0].equals("ADDIDTOAPCLASS"))
            {
                addIDtoAPClass = true;

/*
            } else if (words[0].equals(""))
            {
                 = words[1];

*/
            } //end of keyword identification
        } //end of loop over lines

        //Check Compatibility between parameters
	if (!makeCPMap && !chopMols)
	    readRules = false;
        if (fragmentMaxSize < fragmentMinSize)
            killDueToParams("Check Fragments Size Requirements");
        if (rmDuplicates & groupingRotamers)
            killDueToParams("Check keywords: I cannot goup the rotamers if I have to delete the duplicates fragments");
    }
//------------------------------------------------------------------------------

    private static void killDueToParams(String message)
    {
        System.err.println("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        System.err.println("      Uncompatible keywords or settings!");
        System.err.println("\n"+message);
        System.err.println("\nExecution terminates now!");
        System.exit(0);
    }

//------------------------------------------------------------------------------

    public static void setCuttingRules()
    {

        //Store all-atoms rules: the rule's components used to match any element
        ArrayList<String> keyLinesA = new ArrayList<String>();
        keyLinesA = IOtools.readTXTKeyword(rulesFile,"ANY");
        if (keyLinesA.size() == 0)
        {
            anyAtm.add("[$([*;!#1])]");
            anyAtm.add("[$(*)]");
        } else {
            for (int i = 0; i<keyLinesA.size(); i++)
            {
                try
                {
                    String[] words = keyLinesA.get(i).split("\\s+");
                    String anyAtmRule = words[1]; //rule's component matching any atom
                    anyAtm.add(anyAtmRule);
                } catch (Throwable t) {
                    System.out.println("ERROR in getting 'any-atom' labels!");
                    System.out.println("Chech KEYWORD-labelled line of file \n "+rulesFile+"\n > "+keyLinesA.get(i));
                    System.out.println("Program will terminate now."+t);
                    t.printStackTrace();
                    System.exit(-1);
                }
            }
	}

        //Now get the list of cutting rules
        ArrayList<String> keyLines = new ArrayList<String>();
        keyLines = IOtools.readTXTKeyword(rulesFile,"CTR");
        for (int i = 0; i<keyLines.size(); i++)
        {
            try
            {
                String[] words = keyLines.get(i).split("\\s+");
                String name = words[1]; //name of the rule
                if (words.length < 6)
                {
                    System.err.println("Check Rule! Less than 6 fields");
                    for (String s : words)
                        System.err.println(" I got "+s);
                    System.err.println("Make sure field are separed by spaces");
                }

                // further details in map of options
		ArrayList<String> opts = new ArrayList<String>();
                if (words.length >= 7)
                {
                    for (int wi=6; wi<words.length; wi++)
                    {
                        opts.add(words[wi]);
                    }
                }

		int priority = Integer.parseInt(words[2]);

		GM3DCuttingRule rule = new GM3DCuttingRule(name,
						words[3],  //atom1
						words[4],  //atom2
						words[5],  //bond between 1 and 2
						priority,  
						opts);     

		rules.put(name,rule);

		// ordered list of rules
                if (!sortedRules.containsKey(priority))
                    sortedRules.put(priority,name);
                else {
                    System.out.println("ERROR in getting the cutting rules!");
                    System.out.println("Check priority of rule "+name);
                    System.exit(-1);
                }
            } catch (Throwable t) {
                System.out.println("ERROR in getting the cutting rules!");
                System.out.println("Chech KEYWORD-labelled line number "+(i+1)+" of file \n "+rulesFile+"\n > "+keyLines.get(i));
                System.out.println("Program will terminate now."+t);
                t.printStackTrace();
                System.exit(-1);
            }
        }
    }

//------------------------------------------------------------------------------

    public static Set<String> getIntermediateFiles()
    {
        return interFiles;
    }

//------------------------------------------------------------------------------

    public static String getCurrentSDFile()
    {
        return upToDateSDFile;
    }

//------------------------------------------------------------------------------
//TODO prepared for library conversion
/*    public static String getOldFormat()
    {
        return oldFormat;
    }
//------------------------------------------------------------------------------

    public static String getNewFormat()
    {
        return newFormat;
    }*/
//------------------------------------------------------------------------------

    public static String getLibFormat()
    {
        return libFormat;
    }


//------------------------------------------------------------------------------

    public static void setLibFormat(String format)
    {
        libFormat = format;
    }

//------------------------------------------------------------------------------

    public static void updateStructureFilePointer(String newName)
    {
        //move old file to the list of intermediate files
        if (!upToDateSDFile.equals(sdfInFile))
            interFiles.add(upToDateSDFile);

        //update file name
        upToDateSDFile = newName;
    }

//------------------------------------------------------------------------------

    public static void printParameters()
    {
//TODO update
        System.out.println("\nParameters: ");
        System.out.println(" User's defined parameters from file "+input);
        System.out.println("\n # Type of job #");
        System.out.println(" - check input (with FORMULA): "+chkFormula);
        System.out.println(" - fragmentation:              "+chopMols);
        System.out.println(" - remove duplicates:          "+rmDuplicates);
        System.out.println(" - reorder according to MW:    "+orderMW);
        System.out.println(" - remove intermediate files:  "+removeIntermFiles);
        System.out.println(" - reporting on screen level:  "+report);
//        System.out.println(" -   :"+);
        System.out.println("\n # Input files #");
        System.out.println(" 3D structures from file -> "+sdfInFile);
        if (chkFormula)
            System.out.println(" Formulae from  file -> "+txtFile);
        if (chopMols)
            System.out.println(" Cutting rules from file -> "+rulesFile);
   }

//------------------------------------------------------------------------------

    public static void printCuttingRules()
    {
        System.out.println("\nCutting Rules priority: ");
        for (int i : sortedRules.keySet())
        {
            System.out.println(" Rule "+i+" => "+sortedRules.get(i));
        }
   }

//------------------------------------------------------------------------------

    public static String getUniqueNum()
    {
        //Set unique date-and-time dependent string
        Date date = new Date();
        String dts = String.valueOf(date.getTime());
        return dts;
   }

//------------------------------------------------------------------------------

    public static String getFormattedDate()
    {
	return dateForm;
    }

//------------------------------------------------------------------------------

}
