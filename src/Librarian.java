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

import java.util.Map;
import java.util.HashMap;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.List;
import java.util.Set;
import java.util.ArrayList;
import java.util.Collections;
import java.io.FileNotFoundException;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import java.io.FileInputStream;
import java.io.File;


/**
 * Librarian groups a list of useful methods dealing with 
 * libraries of fragments
 *
 * @author Marco Foscato (University of Bergen)
 */

class Librarian
{

    //Input file
    private String inFile;
    private String inFormat; //prepared for conversion

    //Output
    private String outFile;
    private String outFormat; //prepared for conversion

    //Job name
    private String jobName;

    //Report on screen
    private int repOnScreen;
    private String preStr = "Librarian";

//-------------------------------------------------------------------

    /**
     * Creates a <code>Librarian</code> setting variables according 
     * an <code>Parameters</code>
     */
    public Librarian()
    {
        repOnScreen = Parameters.report;
    }

//-------------------------------------------------------------------

    /**
     * Extract fragments that NOT match the rejection criteria
     * in the same way as for the fragmentation
     */

    public void filter()
    {
        repOnScreen = Parameters.report;
        inFile = Parameters.getCurrentSDFile();
        inFormat = Parameters.getLibFormat();
        String txtFile = Parameters.txtFile;
        String keepFile = "keptMols.sdf";
        String rejFile = "rejectedMols.sdf";

        if (repOnScreen >= 0)
            System.out.println("\n==== Filtering Library of Fragments ====");

        int i = 0;
        int rejNum = 0;
        int keptNum = 0;
        try {
            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(inFile), DefaultChemObjectBuilder.getInstance());
            while (reader.hasNext())
            {
                i++;
                IAtomContainer mol = reader.next();
                if (repOnScreen >= 1)
                    System.out.println("Evaluating entry "+i+": '"+MolecularUtils.getNameOrID(mol)+"'");

//TODO: add automatic detection of fragment's format
                GM3DFragment frag = new GM3DFragment(mol,"DENOPTIM");

                if (FragmentFilter.keepFragment(frag))
                {
                    keptNum++;
                    IOtools.writeSDFAppend(keepFile,mol,true);
                } else {
                    rejNum++;
                    IOtools.writeSDFAppend(rejFile,mol,true);
                }

//Only for debug
//IOtools.pause();
            }
        } catch (FileNotFoundException fnf) {
            System.err.println("File Not Found: " + inFile);
            System.err.println(fnf.getMessage());
            System.exit(-1);
        } catch (Throwable t) {
            System.err.println("\nERROR in reading the file with IteratorMDLReader (extractSMARTS). "+t);
            t.printStackTrace();
            System.out.println("\nERROR: cannot iterate through MDL file. This "
                                + "might be becouse the MDL format is prior to "
                                + "V2000 or because there are bond with "
                                + "type='any' (type=8 in connectivity matrix). "                                + "Please, try to provide a V2000 format and "
                                + "convert avoid bond type=8.");
            System.exit(0);
        }

        if (repOnScreen >= 0)
        {
            System.out.println("\n=========== Filtering - DONE ===========");
            System.out.println("I red "+i+" molecules: \n"+
                                keptNum+" molecules kept (see file "+keepFile+");\n"+
                                rejNum+" hit molecules one of the rejection criteria (see file "+rejFile+")");
        }

    }

//-------------------------------------------------------------------

    /**
     * Extract fragments matching the SMARTS 
     * (target.sdf = output mols matching the criteria)
     * (others.sdf = mols not matching the criteria)
     * @param  
     */

    public void extractSMARTS()
    {
        repOnScreen = Parameters.report;
        inFile = Parameters.getCurrentSDFile();
        inFormat = Parameters.getLibFormat();
        String txtFile = Parameters.txtFile;
        String trgtFile = "targets.sdf";
        String othersFile = "others.sdf";

        if (repOnScreen >= 0)
            System.out.println("\n=== Filtering Library (using SMARTS) ====");


        //get Smarts
        Map<String,String> smarts = new HashMap<String,String>();
        ArrayList<String> listOfQueries = IOtools.readTXT(txtFile);
        for (int i=0; i<listOfQueries.size(); i++)
            smarts.put(Integer.toString(i),listOfQueries.get(i));

        int i = 0;
        int mi = 0;
        int nmi = 0;
        try {
            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(inFile), DefaultChemObjectBuilder.getInstance());
            while (reader.hasNext())
            {
                i++;
                IAtomContainer mol = reader.next();
                if (repOnScreen >= 1)
                    System.out.println("Evaluating entry "+i+": '"+MolecularUtils.getNameOrID(mol)+"'");

                ManySMARTSQuery msq = new ManySMARTSQuery(mol,smarts);
                int totMatches = msq.getTotalMatches();
                if (totMatches != 0)
                {
                    if (repOnScreen >= 1)
                        System.out.println(" -> Matches per query ID: "+msq.getNumMatchesMap());
                    IOtools.writeSDFAppend(trgtFile,mol,true);
                    mi++;
                } else {
                    if (repOnScreen >= 2)
                        System.out.println(" -> No match");
                    nmi++;
                    IOtools.writeSDFAppend(othersFile,mol,true);
                }
//Only for debug
//IOtools.pause();
            }
        } catch (FileNotFoundException fnf) {
            System.err.println("File Not Found: " + inFile);
            System.err.println(fnf.getMessage());
            System.exit(-1);
        } catch (Throwable t) {
            System.err.println("\nERROR in reading the file with IteratorMDLReader (extractSMARTS). "+t);
            t.printStackTrace();
            System.out.println("\nERROR: cannot iterate through MDL file. This "
                                + "might be becouse the MDL format is prior to "
                                + "V2000 or because there are bond with "
                                + "type='any' (type=8 in connectivity matrix). "                                + "Please, try to provide a V2000 format and "
                                + "convert avoid bond type=8.");
            System.exit(0);
        }

        if (repOnScreen >= 0)
        {
            System.out.println("\n===== SMARTS-bases extraction DONE =====");
            System.out.println("I red "+i+" molecules: \n"+
                                mi+" matched the SMARTS query (see "+trgtFile+");\n"+
                                nmi+" didn't (see "+othersFile+")");
	}
    }

//-------------------------------------------------------------------

    /**
     * Extract fragments from the library according to the CLASS
     * of their attachment points.
     * (target.sdf = output mols matching the criteria)
     * (others.sdf = mols not matching the criteria)
     * @param  
     */

    public void extractClass()
    {
        repOnScreen = Parameters.report;
        inFile = Parameters.getCurrentSDFile();
        inFormat = Parameters.getLibFormat();
        String txtFile = Parameters.txtFile;
        String trgtFile = "targets.sdf";
        String othersFile = "others.sdf";

        if (repOnScreen >= 0)
            System.out.println("\n==== Filtering Library (using CLASS) ====");

        ArrayList<String> listOfTargets = IOtools.readTXT(txtFile);
        int i = 0;
        int mi = 0;
        int nmi = 0;
        try {
            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(inFile), DefaultChemObjectBuilder.getInstance());
            while (reader.hasNext())
            {
                i++;
                IAtomContainer mol = reader.next();

                if (repOnScreen >= 1)
                    System.out.println("Evaluating entry "+i+": '"+MolecularUtils.getNameOrID(mol)+"'");

                GM3DFragment frag = new GM3DFragment(mol,inFormat);

                if (frag.getNumberOfAttachmentPoints() == 0)
                {
                    if (repOnScreen >= 1)
                        System.out.println("Molecule "+i+" has no attachment point! Discharged!");
                    continue;
                }

                ArrayList<GM3DAttachmentPoint> aps = frag.getAllAPs();

                boolean found = false;
                for (GM3DAttachmentPoint ap : aps)
                {
                    String apclass = ap.getAPClass();
                    for (String target : listOfTargets)
                     {
//System.out.println("Comparing: "+apclass+" and "+target);
                        if (apclass.contains(target))
                        {
                            IOtools.writeSDFAppend(trgtFile,mol,true);
                            mi++;
                            if (repOnScreen >= 1)
                                System.out.println(" -> Matches target '"+target+"'");
                            found = true;
                            break;
                        }
                    }

                    if (found)
                        break;
                }

                if (!found)
                {
                    nmi++;
                    IOtools.writeSDFAppend(othersFile,mol,true);
                }
            }
        } catch (FileNotFoundException fnf) {
            System.err.println("File Not Found: " + inFile);
            System.err.println(fnf.getMessage());
            System.exit(-1);
        } catch (Throwable t) {
            System.err.println("\nERROR in reading the file with IteratorMDLReader. "+t);
            t.printStackTrace();
            System.out.println("\nERROR: cannot iterate through MDL file. This "
                                + "might be becouse the MDL format is prior to "
                                + "V2000 or because there are bond with "
                                + "type='any' (type=8 in connectivity matrix). "                                + "Please, try to provide a V2000 format and "
                                + "convert avoid bond type=8.");
            System.exit(0);
        }

        if (repOnScreen >= 0)
	{
            System.out.println("\n===== Filtering using CLASS - DONE ======");
            System.out.println("I red "+i+" molecules: \n"+
                                mi+" contained target Classes (see "+trgtFile+");\n"+
                                nmi+" didn't (see "+othersFile+").");
	}
    }

//-------------------------------------------------------------------

    /**
     *
     * @param  
     */

    public void mergeLibraries()
    {        
        preStr = preStr+"MERGE-";
        if (repOnScreen >= 1)
            System.out.println("\n======= Merging libraries of fragments =======");

        //Set formats
        inFormat = Parameters.getLibFormat();
        outFormat = Parameters.getLibFormat();

        //Set Filenames
        inFile = Parameters.getCurrentSDFile();
        jobName = Parameters.getJobName();
        String oldFile = Parameters.previousLibrary; //already existing library
        String uniqueLib = Parameters.uniqueLibrary; //list of unique fragments
        String outUniqueLib = "Update-"+jobName+"_"+uniqueLib;
        outFile = "Combined_"+jobName+"_"+oldFile;
        
        //Decide the type of job accordingly to the existence of a previous list of unique fragments
        if (uniqueLib == null)
        {
            if (repOnScreen >= 1)
            {
                System.out.println("No uniqueFragments Library provided =>");
                System.out.println("File: "+inFile+" is just appended to a copy of file "+oldFile);
            }
            //just append to outFile
            IOtools.copyFile(oldFile,outFile,false);
            IOtools.copyFile(inFile,outFile,true);
        } else {
            if (repOnScreen >= 1)
                System.out.println("Using Library "+uniqueLib+" to group the rotamers");
            //First copy the existing library and uniqueLib (Avoid to compromize these files)
            IOtools.copyFile(oldFile,outFile,false);
            IOtools.copyFile(uniqueLib,outUniqueLib,false);
            Parameters.uniqueLibrary = outUniqueLib;
            //Then add the new fragments generating the 'ISOMER' field
            groupRotamers(inFile,outFile,outUniqueLib);
        }

        //Print a summary
        if (repOnScreen >= 1)
        {
            System.out.println("\n================ Merging DONE ================");
            System.out.println("Combined library of all fragments   : "+outFile);
            System.out.println("Combined library of unique fragments: "+outUniqueLib);
        }

        //Redirect name of input SDF file
        Parameters.updateStructureFilePointer(outFile);
    }

//-------------------------------------------------------------------

    /**
     *
     * @param  
     */

    public void mwMerge()
    {
        //Get parameters
        Set<String> files = Parameters.MWMpath;
        String thisJob = Parameters.getJobName(); // thisJob = bin+"_"+jobUnq
        Map<String,String> binAndJob = getMWBinAndJob(thisJob);
        String bin = binAndJob.get("bin");
        int numMWBin = files.size();
        outFile = "AllMWBin_"+thisJob+".sdf";
        String isoCount = "IsomerCounts-AllMWBin_"+thisJob+".txt";
        String frgIsoKeyConverter ="RunAllFrgKeyConversion.sh";

        if (repOnScreen >= 1)
        {
            System.out.println("\n========= Merging MW-range libraries =========");
            System.out.println("Number of sub libraries to merge: "+numMWBin);
        }

        //Check for isomer count files
        for (String f : files)
        {
            File binFile = new File(f);
            Map<String,String> locBinAndJob = getMWBinAndJob(f);
            String isoCountFileName = "IsomerCountsJob_"+locBinAndJob.get("job")+"_"+locBinAndJob.get("jobUnq")+".dat";
            File isoCountFile = new File(isoCountFileName);
            if (!isoCountFile.exists())
            {
                System.err.println("\nERROR! File "+isoCountFile+" NOT FOUND! It's required by 'mwMErge()' from "+f);
                System.exit(-1);
            }
        }

        //Process all the libraries contributing to this bin
//        int totUnqFrgs = 0;
        boolean first = true;
        String propName = "ISOMER";
        Map<String,Integer> counts = new HashMap<String,Integer>();
//        Map<String,ArrayList<String>> listIDequivalence = new HashMap<String,ArrayList<String>>();
        for (String subLibName : files)
        {
            if (repOnScreen >= 1)
                System.out.println("Merging unique fragments from "+subLibName);

            //Set name of the piece of library containins all fragments
            Map<String,String> locBinAndJob = getMWBinAndJob(subLibName);
            String allFrgInBinFile = "MWBin_"+bin+"_AllFrg_Job_"+locBinAndJob.get("job")+"_"+locBinAndJob.get("jobUnq")+".sdf";

            if (first)
            {
                //Append ALL Frags from the first library
                try {
                    IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(subLibName), DefaultChemObjectBuilder.getInstance());
                    while (reader.hasNext())
                    {
                        IAtomContainer frag = reader.next();

                        //TODO add automated detection of the fragment format
                        String fragFormat = "DENOPTIM";

// WARNING!
// WARNING! Assuming that these are Unique fragments and that there is a file with counting
// WARNING!

                        String oldKey = frag.getProperty(propName).toString();
//         System.err.println("oldKey = "+oldKey);
                        int localCount = getLocalCount(oldKey);
//                        ArrayList<String> equivKeys = new ArrayList<String>();
//                        equivKeys.add(oldKey);
//                        listIDequivalence.put(oldKey,equivKeys);
                        counts.put(oldKey,localCount);
//                        frag.setProperty(propName,groupID);
//         System.err.println("oldKey = "+oldKey+" groupID= "+groupID);
                        IOtools.writeSDFAppend(outFile,frag,true);
                    }
                    reader.close();
	        } catch (FileNotFoundException fnf) {
	            System.err.println("File Not Found: " + inFile);
	            System.err.println(fnf.getMessage());
	            System.exit(-1);
                } catch (Throwable t) {
                    System.err.println("\nERROR in reading the file with IteratorMDLReader."+t);
                    t.printStackTrace();
                    System.out.println("\nERROR: cannot iterate through MDL "
				+ "file. This "
                                + "might be becouse the MDL format is prior to "
                                + "V2000 or because there are bond with "
                                + "type='any' (type=8 in connectivity matrix). "                                + "Please, try to provide a V2000 format and "
                                + "convert avoid bond type=8.");
                    System.exit(0);
                }
                first=false;
            } else {
                try {
                    IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(subLibName), DefaultChemObjectBuilder.getInstance());
                    while (reader.hasNext())
                    {
                        IAtomContainer mol = reader.next();
                        boolean isUnique = true;
                        String oldKey = mol.getProperty(propName).toString();
                        int localCount = getLocalCount(oldKey);

                        //TODO add automated detection of the fragment format
                        String fragFormat = "DENOPTIM";
                        GM3DFragment frag = new GM3DFragment(mol,fragFormat);                        

                        IteratingMDLReader readUnq = new IteratingMDLReader(new FileInputStream(outFile), DefaultChemObjectBuilder.getInstance());
                        while (readUnq.hasNext())
                        {
                            IAtomContainer molUnq = readUnq.next();
                            GM3DFragment fragUnq = new GM3DFragment(molUnq,fragFormat);
                            if (frag.sameFragOf(fragUnq))
                            {
                                //Get the ID of the uniqueFragment
                                String unqKey = fragUnq.getProperty(propName).toString();
                                //Update count of group members
                                int oldCount = counts.get(unqKey);
                                int newCount = oldCount + localCount;
                                counts.put(unqKey,newCount);
                                //Set equivalence key
//                                if (!listIDequivalence.get(unqKey).contains(oldKey))
//                                    listIDequivalence.get(unqKey).add(oldKey);
                                //Add Isomer Key conversion
                                IOtools.writeTXTAppend(frgIsoKeyConverter,"sed -i -e 's/^"+oldKey+"$/"+unqKey+"/g' "+allFrgInBinFile+" \n",true);
                                //Set flag to avoid reporting it as unique
                                isUnique = false;
                                break;
                            }
                        }
                        readUnq.close();
                        //Report as unique
                        if (isUnique)
                        {
//                            ArrayList<String> equivKeys = new ArrayList<String>();
//                            listIDequivalence.put(oldKey,equivKeys);
                            counts.put(oldKey,localCount);
                            IOtools.writeSDFAppend(outFile,mol,true);
                        } 
                    }
                    reader.close();
	        } catch (FileNotFoundException fnf) {
	            System.err.println("File Not Found: " + inFile);
	            System.err.println(fnf.getMessage());
	            System.exit(-1);
                } catch (Throwable t) {
                    System.err.println("\nERROR in reading the file with IteratorMDLReader. "+t);
                    t.printStackTrace();
                    System.out.println("\nERROR: cannot iterate through MDL "
				+ "file. This "
                                + "might be becouse the MDL format is prior to "
                                + "V2000 or because there are bond with "
                                + "type='any' (type=8 in connectivity matrix). "                                + "Please, try to provide a V2000 format and "
                                + "convert avoid bond type=8.");
                    System.exit(0);
                }
            }
        }

        //Report isomer counting
        int tot = 0;
        int totUnq = 0;
        for (String key : counts.keySet())
        {
            String line ="Isomer#"+thisJob+"#"+key+"# entries "+counts.get(key)+"\n";
            totUnq = totUnq + 1;
            tot = tot + counts.get(key);
            IOtools.writeTXTAppend(isoCount,line,true);
        }
        String ending = "In BIN "+bin+" Total Fragments = "+tot
                        +" Total Unique Fragments = "+totUnq+" \n";
        IOtools.writeTXTAppend(isoCount,ending,true);
        String allBINs = "all_BINs_Counting.dat";
        IOtools.writeTXTAppend(allBINs,ending,true);
/*
        //Prepare script for converting frag ISOMER property
        String frgIsoKeyConverter ="RunAllFrgKeyConversion.sh";
String allFrgBin = "MWBin_"+subLibIdx+"_AllFrg_Job_"+thisJob+".sdf";
        String allFragmentsFile = "allFragments.sdf";
        IOtools.writeTXTAppend(frgIsoKeyConverter,"sed -i ",true);
        for (String unqID : listIDequivalence.keySet())
        {
            ArrayList<String> allEquiv= listIDequivalence.get(unqID);
            for (int i=0; i< allEquiv.size(); i++)
            {
                System.out.println(" -e 's/"+allEquiv.get(i)+"/"+unqID+"/g' ");
            }
        }
        IOtools.writeTXTAppend(frgIsoKeyConverter," "+allFragmentsFile,true);

*/

        //Goodbye
        if (repOnScreen >= 1)
        {
//            int totUnqFrgs = listIDequivalence.keySet().size();
            int totUnqFrgs = counts.keySet().size();
            System.out.println("\n================ Merging DONE ================");
            System.out.println("Number of Unique fragments found in BIN "+bin+": "+totUnqFrgs);
            System.out.println("Combined library of all fragments   : "+outFile);
            System.out.println("Isomers counting                    : "+isoCount);
            System.out.println("See also the global file            : "+allBINs);

        }

        //Redirect name of input SDF file
        Parameters.updateStructureFilePointer(outFile);
    }
//------------------------------------------------------------------------------
    private static Map<String,String> getMWBinAndJob(String file)
    {
        Map<String,String> binAndJob = new HashMap<String,String>();
        File f = new File(file);
        String name = f.getName();
        //get rid of extension
        if (name.contains("."))
            name = name.substring(0,name.indexOf("."));
        String[] parts = name.split("_");
        //In case of MW-range fragments library
        if (name.startsWith("MWBin_"))
        {
            binAndJob.put("bin",parts[1]);
            binAndJob.put("job",parts[3]);
            binAndJob.put("jobUnq",parts[4]);
        } else if (name.startsWith("IsomerCountsJob_"))
        {
            binAndJob.put("bin",null);
            binAndJob.put("job",parts[1]);
            binAndJob.put("jobUnq",parts[2]);
        } else if (name.startsWith("Fragments_"))
        {
            binAndJob.put("bin",null);
            binAndJob.put("job",parts[2]);
            binAndJob.put("jobUnq",parts[3]);
        } else if (name.startsWith("InpMrgBin_"))
        {
            binAndJob.put("bin",parts[1]);
            binAndJob.put("job",null);
            String[] p2 = parts[2].split("#");
            binAndJob.put("jobUnq",p2[0]);
/*
        } else if (name.startsWith(""))
        {
            binAndJob.put("bin",parts[]);
            binAndJob.put("job",parts[]);
            binAndJob.put("jobUnq",parts[]);
*/
        } else {
            System.err.println("\nERROR! Problem in gettin BIN and JOB details from "+file+"\n");
            System.exit(-1);
        }
        
        return binAndJob;
    }

//------------------------------------------------------------------------------
    /**
     * Looks for the number of times a fragment has been found in a local
     * library generated by MW-based splitting. This method assumes that text 
     * files generated by GM3DFragmenter with method <code>mwSplitLib</code>
     */
    private static int getLocalCount(String isomerID)
    {
        int num = 0;
/*
        Map<String,String> binAndJob = getMWBinAndJob(libraryName);
        String bin = binAndJob.get("bin");
        String job = binAndJob.get("job");
        String jobUnq = binAndJob.get("jobUnq");

/*        File subLibFile = new File(libraryName);
        
        String[] parts = subLibFile.getName().split("_");
        //parts[0] is eqbinAndJobual to MWBin<BIN_INDEX>
        String bin = parts[0].substring(6); 
        String[] parts2 = parts[1].split("\\.");
        Long jbn = Long.parseLong(parts2[0]);
        //As for 'mwSplitLib' the name of the txt file is the following
*/
        String[] p = isomerID.split("-");
        String groupID = p[1];
        String[] jp = p[0].split("_");
        String job = jp[0];
        String jobUnq = jp[1];

        //look for the file with the previous counting
        String ismCntFileName = "IsomerCountsJob_"+job+"_"+jobUnq+".dat"; 
        String key = "Isomer:"+groupID;
        ArrayList<String> lines = IOtools.readTXTKeyword(ismCntFileName,key);
        if (lines.size() != 1)
        {
            System.err.println("ERROR! Unable to get isomer counting for "+ismCntFileName+" key:_"+key+"_\n");
            System.exit(-1);
        }
        String[] words = lines.get(0).split("\\s+");
        String[] entries = words[1].split(":");
        num = Integer.parseInt(entries[1]);        
        return num;
    }

//-------------------------------------------------------------------
    /**
     *
     * @param  
     */

    public void groupRotamers()
    {
        if (repOnScreen >= 1)
            System.out.println("\n========= Search of Rotamers: STARTS =========");

        inFile = Parameters.getCurrentSDFile();
        jobName = Parameters.getJobName();
        outFile = "RotoGrouped-"+inFile;
        String uniqueFile = "UniqueFrags-"+jobName+".sdf";

        inFormat = Parameters.getLibFormat();
        outFormat = Parameters.getLibFormat();

        //reorder
        groupRotamers(inFile,outFile,uniqueFile);

        //Redirect name of input SDF file
        Parameters.updateStructureFilePointer(outFile);
    }

//-------------------------------------------------------------------

    /**
     *
     * @param  
     */

    public void groupRotamers(String inFile, String outFile, String uniqueFile)
    {
        //Prepare the storage for unique fragments and infos
        int totFrags = 0;
        String propName = "ISOMER";
        SortedMap<Integer,Integer> counts = new TreeMap<Integer,Integer>();
        int groupID = -1;
        //get previously existing count of unique fragments, if any
        File unqFile = new File(uniqueFile);
        if (unqFile.exists())
        {
            try {
                IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(unqFile), DefaultChemObjectBuilder.getInstance());
                while (reader.hasNext())
                {
                        String prop = reader.next().getProperty(propName).toString();
                        int oldGroupID = Integer.parseInt(prop);
                    if (!counts.keySet().contains(oldGroupID))
                        counts.put(oldGroupID,1);
                    else {
                        System.err.println("\nERROR! File "+uniqueFile+" contains dublicated "+propName+"-ID"); 
                        System.exit(0);
                    }
                }
                reader.close();
            } catch (FileNotFoundException fnf) {
                System.err.println("File Not Found: " + inFile);
                System.err.println(fnf.getMessage());
                System.exit(-1);
            } catch (Throwable t) {
                System.err.println("\nERROR in reading the file with IteratorMDLReader. "+t);
                t.printStackTrace();
                System.out.println("\nERROR: cannot iterate through MDL file. This "
                                + "might be becouse the MDL format is prior to "
                                + "V2000 or because there are bond with "
                                + "type='any' (type=8 in connectivity matrix). "                                + "Please, try to provide a V2000 format and "
                                + "convert avoid bond type=8.");
                System.exit(0);
            }
        }

        //Loop over input fragments
        try {
            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(inFile), DefaultChemObjectBuilder.getInstance());
            while (reader.hasNext())
            {
                totFrags++;
                if (repOnScreen >= 2)
                    System.out.println("Check fragment "+totFrags);
                IAtomContainer mol = reader.next();
                GM3DFragment frag = new GM3DFragment(mol,inFormat);
                boolean isUnique = true;
                if (unqFile.exists())
                {
                //Compare the fragment with the unique fragments
                IteratingMDLReader readUnq = new IteratingMDLReader(new FileInputStream(unqFile), DefaultChemObjectBuilder.getInstance());
                while (readUnq.hasNext())
                {
                    //Compare the fragment with the unique fragments
                    IAtomContainer molUnq = readUnq.next();

//TODO does this really speed up the execution?
//TODO It was designed before having all this file reading/writing
                    //Rough compatison to speed up the loop
                    if (mol.getAtomCount() != molUnq.getAtomCount())
                        continue;                        

                    //Deep comparison
                    GM3DFragment fragUnq = new GM3DFragment(molUnq,inFormat);
                    if (frag.sameFragOf(fragUnq))
                    {
                        //Get the ID of the uniqueFragment
                        String prop = fragUnq.getProperty(propName).toString();
                        int existingGroupID = Integer.parseInt(prop);
                        //Update count of group members
                        int oldCount = counts.get(existingGroupID);
                        counts.put(existingGroupID,oldCount + 1);
                        //Set the flag of the new fragment
                        mol.setProperty(propName,existingGroupID);
                        //Store the new fragment
                        IOtools.writeSDFAppend(outFile,mol,true);
                        isUnique = false;
                        break;
                    }
                } 
                readUnq.close();
                }

                //Add a new unique fragment
                if (isUnique)
                {
                    groupID++;
                    frag.setProperty(propName,groupID);
                    mol.setProperty(propName,groupID);
                    IOtools.writeSDFAppend(outFile,mol,true);
                    IOtools.writeSDFAppend(uniqueFile,frag.toIAtomContainer(outFormat),true);
                    counts.put(groupID,1);
                }
            }
            reader.close();
        } catch (FileNotFoundException fnf) {
            System.err.println("File Not Found: " + inFile);
            System.err.println(fnf.getMessage());
            System.exit(-1);
        } catch (Throwable t) {
            System.err.println("\nERROR in reading the file with IteratorMDLReader. "+t);
            t.printStackTrace();
            System.out.println("\nERROR: cannot iterate through MDL file. This "
                                + "might be becouse the MDL format is prior to "
                                + "V2000 or because there are bond with "
                                + "type='any' (type=8 in connectivity matrix). "                                + "Please, try to provide a V2000 format and "
                                + "convert avoid bond type=8.");
            System.exit(0);
        }

        //Print a summary
        if (repOnScreen >= 1)
        {
            System.out.println("Total number of fragments added       : "+totFrags);
            System.out.println("Total number of Unique Fragments added: "+counts.size());
        }
    }

//-------------------------------------------------------------------

    /**
     *
     * @param  
     */

    public void reorderMW(boolean ascending)
    {

        if (repOnScreen >= 1)
        {
            String ord = "ascending";
            if (!ascending)
                ord = "decending";
            System.out.println("\n========== Reordering fragments (MW) =========");
            System.out.println("Reorder according to "+ord+" Molecular Weight");
        }

        inFile = Parameters.getCurrentSDFile();
        outFile = "MWsorted-"+inFile;

        //reorder
        reorderMW(inFile,outFile,ascending);

        //Redirect name of input SDF file
        Parameters.updateStructureFilePointer(outFile);
    }

//-------------------------------------------------------------------

    /**
     * Reorder a file containing molecular objects according
     * to their molecular weight
     * @param fileIn SDF file with imput list
     * @param fileOut SDF file output file
     * @param ascending set <code>true</code> for ascending order
     */
    public void reorderMW(String fileIn, String fileOut, boolean ascending)
    {

        //get input
        ArrayList<IAtomContainer> listIn = IOtools.readSDF(fileIn);

        //sort the list
        Collections.sort(listIn, new MWComparator());

        if (!ascending)
            Collections.reverse(listIn);

        //report output
        for (int i = 0; i<listIn.size(); i++)
        {
            IOtools.writeSDFAppend(outFile,listIn.get(i),true);
        }
    }

//-------------------------------------------------------------------

    /**
     * Initializer for splitting the library of fragments in sub 
     * libraries of given MW wideness. 
     */
    public void mwSplitting()
    {
        //for initializing parameters
        //TODO

        //Split existing library
        //TODO controll keep clones
        mwSplitLib(Parameters.getCurrentSDFile(),Parameters.massBinSize,false);
    }

//-------------------------------------------------------------------

    /**
     * Split a library of fragments into sub libraries that
     * consider only a selected range on molecular weights
     * @param inLib name of the input file - original library
     * @param binSize Molecular Weight wideness of each sub library
     * @param keepClones if <code>true</code> duplicates entries are
     * ignored
     */

    public void mwSplitLib(String inLib, int binSize, boolean keepClones)
    {
        //Starting message
        if (repOnScreen >= 1)
        {
            System.out.println("\n======== Splitting Library - MW-based ========");
            System.out.println("Library: "+inLib+"\nMW-wideness: "+binSize);
        }

        //Counters and storage of infos
        int totFrags = 0;
        int groupID = -1;
        String propName = "ISOMER";
        SortedMap<Integer,Integer> counts = new TreeMap<Integer,Integer>();
        Map<Integer,String> isomerProp = new HashMap<Integer,String>();
        //Get job details
        Map<String,String> binAndJob = getMWBinAndJob(inLib);
        String jobUnq = binAndJob.get("jobUnq");
//        String bin = binAndJob.get("bin");
      String job = binAndJob.get("job");
        String thisJob = job+"_"+jobUnq;

/*        File in = new File(inLib);
        String nameOfSource = in.getName();
        String[] pNameOfSource = nameOfSource.split("_");
        String jobId = pNameOfSource[pNameOfSource.length-1];
        //Unique jobId
        String thisJob = "";
        if (jobId.contains("."))
            jobId = jobId.substring(0,jobId.indexOf("."));
        if (jobId.length() == 13)
            thisJob = jobId;
        else
            thisJob = Parameters.getUniqueNum();
        //Local job ID
        if (pNameOfSource.length > 2)
             thisJob = pNameOfSource[pNameOfSource.length-2] + "_" + thisJob;
*/

        //Read the library and generate sub libraries covering a range of MW
        try {
            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(inLib), DefaultChemObjectBuilder.getInstance());
            while (reader.hasNext())
            {
                IAtomContainer mol = reader.next();

                //TODO add automated detection of the fragment format
                String fragFormat = "DENOPTIM";
                double mw = -1;
                try {

                    MolecularFormulaManipulator mf = new MolecularFormulaManipulator();
                    mw = mf.getNaturalExactMass(mf.getMolecularFormula(mol));
                } catch (Throwable t) {
                    System.out.println("ERROR! Exception while calculating Molecular Weight (LIB)");
                    System.exit(0);
                }
                double ratio = mw / (double) binSize;
                int subLibIdx = (int) ratio;
                String subLibName = "MWBin_"+subLibIdx+"_Job_"+thisJob+".sdf";
                String allFrgBin = "MWBin_"+subLibIdx+"_AllFrg_Job_"+thisJob+".sdf";
                if (!keepClones)
                {
                    File f = new File(subLibName);
                    if (f.exists())
                    {
                        boolean isUnique = true;
                        GM3DFragment frag = new GM3DFragment(mol,fragFormat);
                        IteratingMDLReader readUnq = new IteratingMDLReader(new FileInputStream(subLibName), DefaultChemObjectBuilder.getInstance());

                        String frgUnqID = "";
                        while (readUnq.hasNext())
                        {
                            IAtomContainer molUnq = readUnq.next();
                            GM3DFragment fragUnq = new GM3DFragment(molUnq,fragFormat);
                            if (frag.sameFragOf(fragUnq))
                                {
                                //Get the ID of the uniqueFragment
                                frgUnqID = fragUnq.getProperty(propName).toString();
                                String[] pw = frgUnqID.split("-");
                                int existingGroupID = Integer.parseInt(pw[1]);
                                //Update count of group members
                                int oldCount = counts.get(existingGroupID);
                                counts.put(existingGroupID,oldCount + 1);
                                isUnique = false;
                                break;
                            }
                        }
                        readUnq.close();

                        if (isUnique)
                        {
                            groupID++;
                            counts.put(groupID,1);
                            String isomerName = thisJob+"-"+groupID;
                            String props = " Ref:"+isomerName+" BIN:"+subLibIdx+" MW:"+mw;
                            isomerProp.put(groupID,props);
                            mol.setProperty(propName,isomerName);
                            IOtools.writeSDFAppend(subLibName,mol,true);
                            IOtools.writeSDFAppend(allFrgBin,mol,true);
                        } else {
                            mol.setProperty(propName,frgUnqID);
                            IOtools.writeSDFAppend(allFrgBin,mol,true);                            
                        }
                    } else {
                        groupID++;
                        counts.put(groupID,1);
                        String isomerName = thisJob+"-"+groupID;
                        String props = " Ref:"+isomerName+" BIN:"+subLibIdx+" MW:"+mw;
                        isomerProp.put(groupID,props);
                        mol.setProperty(propName,isomerName);
                        IOtools.writeSDFAppend(subLibName,mol,true);
                        IOtools.writeSDFAppend(allFrgBin,mol,true);
                    }
                } else {
                    //in case of keepClones==true
                    IOtools.writeSDFAppend(subLibName,mol,true);
                }
            }
            reader.close();
        } catch (FileNotFoundException fnf) {
            System.err.println("File Not Found: " + inFile);
            System.err.println(fnf.getMessage());
            System.exit(-1);
        } catch (Throwable t) {
            System.err.println("\nERROR in reading the file with IteratorMDLReader. "+t);
            t.printStackTrace();
            System.out.println("\nERROR: cannot iterate through MDL file. This "
                                + "might be becouse the MDL format is prior to "
                                + "V2000 or because there are bond with "
                                + "type='any' (type=8 in connectivity matrix). "                                + "Please, try to provide a V2000 format and "
                                + "convert avoid bond type=8.");
            System.exit(0);
        }

        //Report counting
        String countsReport = "IsomerCountsJob_"+thisJob+".dat";
        for (Integer i : counts.keySet())
        {
            String line ="Isomer:"+i+" Entries:"+counts.get(i)+isomerProp.get(i)+"\n";
            IOtools.writeTXTAppend(countsReport,line,true);
        }
    }

//-------------------------------------------------------------------

//-------------------------------------------------------------------

    /**
     * TODO
     * @param  
     */
/*
    public 
    {
    }
*/
//-------------------------------------------------------------------
}
