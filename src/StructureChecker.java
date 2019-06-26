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
import java.util.Set;
import java.util.HashSet;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import org.openscience.cdk.ChemFile;

import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.DefaultChemObjectBuilder;

/**
 * Toolbox checking chemical structures.
 * 
 * @author Marco Foscato (University of Bergen) Discussions: Giovanni Occhipinti and Vidar R. Jensen (University of Bergen)
 */

public class StructureChecker
{

    //Amount of information written on screen
    private int repOnScreen;

    //Formulae file
    private String formulaeFile;

    //Structure file
    private String struct3dFile;

    //identification of the jobfiles
    private String jobName;

//------------------------------------------------------------------------------

    public StructureChecker()
    {
        formulaeFile = Parameters.txtFile;
        struct3dFile = Parameters.getCurrentSDFile();
        jobName = Parameters.getJobName();
        repOnScreen = Parameters.report;
    }

//------------------------------------------------------------------------------

    /**
     * Creates a StructureChecker specifying the input files
     * @param struct3dFile SD/SDF file containing molecular structures to check
     * @param formulaeFile TXT file containing molecular formulae
     */
    public StructureChecker(String struct3dFile, String formulaeFile)
    {
        this.struct3dFile = struct3dFile;
        this.formulaeFile = formulaeFile;
        this.jobName = struct3dFile.substring(struct3dFile.length()-4,
                                                struct3dFile.length());
        this.repOnScreen = 1;
    }

//------------------------------------------------------------------------------
    /**
     * Change the pointers so the input file of the next operation 
     * is the output of the one from which this method is called
     * @param freshOut output file from the method calling
     * updateFileName
     */
    private void updateFileName(String freshOut)
    {
        this.struct3dFile = freshOut;
        Parameters.updateStructureFilePointer(freshOut);        
    }

//------------------------------------------------------------------------------
    /**
     * Prefiltering of structures
     */

    public void preFilter()
    {

        if (repOnScreen >= 0)
        {
            String header = "\n========= pre-Filter of structures ==========";
            System.out.println(header);
        }

        int totMols = 0;
        int checkMols = 0;
        String checkfile = "check-preFilter_"+jobName+".sdf";
        String outSDfile = "3DreadyFixedFiltered_"+jobName+".sdf";

        //Iterate over structures
        try {
            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(struct3dFile), DefaultChemObjectBuilder.getInstance());
            //Analysis for each molecule
            while (reader.hasNext())
            {
                totMols++;
                IAtomContainer mol = reader.next();
                if (repOnScreen >= 1)
                    System.out.println("Working on mol "+totMols+" -> "+MolecularUtils.getNameOrID(mol));

		Map<String,String> smarts = new HashMap<String,String>();
		int is = 0;
		for (String s : Parameters.preFilterSMARTS)
		{
		    is++;
		    smarts.put(Integer.toString(is),s);
		}
		ManySMARTSQuery msq = new ManySMARTSQuery(mol,smarts);
	        if (msq.hasProblems())
        	{
	            String cause = msq.getMessage();
		    checkMols++;
		    if (repOnScreen >= 1)
	                System.out.println("ERROR in detecting SMARTS for pre-filtering. Rejecting molecule ("+cause+")");
		    continue;
		} else if (msq.getTotalMatches() > 0) 
		{
                    checkMols++;
                    if (repOnScreen >= 1)
                        System.out.println("Molecule rejected during pre-filtering of structures");
                    IOtools.writeSDFAppend(checkfile,mol,true);
                    continue;
		}

                //Write out the surviving molecule
                IOtools.writeSDFAppend(outSDfile,mol,true);
            }
        } catch (FileNotFoundException fnf) {
            System.err.println("File Not Found: " + struct3dFile);
            System.err.println(fnf.getMessage());
            System.exit(-1);
        } catch (Throwable t) {
            System.err.println("\nERROR in reading the file with IteratorMDLReader (in pre-Filter). "+t);
            t.printStackTrace();
            System.out.println("\nERROR: cannot iterate through MDL file. This "
                                + "might be becouse the MDL format is prior to "
                                + "V2000 or because there are bond with "
                                + "type='any' (type=8 in connectivity matrix). "                                + "Please, try to provide a V2000 format and "
                                + "convert avoid bond type=8.");
            System.exit(0);
        }

        //Redirect name of input SDF file
        updateFileName(outSDfile);

        //Print out summary of results
        System.out.println("\n========== Results - pre-Filtering ==========");
        System.out.println("Total number of molecules: "+totMols);
        System.out.println("Output file: "+outSDfile);
        if (checkMols > 0)
            System.out.println("\nCheck "+checkMols+" problematic cases in "+checkfile);
    }

//------------------------------------------------------------------------------
    /**
     * Analyze the molecules and trys to fix chemical representation problems
     * like aromatic form (kekulization)
     */

    public void fixChemicalRepresentation()
    {
        if (repOnScreen >= 0)
        {
            String header = "\n======== Analysing and fixing inputs ========";
            System.out.println(header);
        }

        int totMols = 0;
        int fixedMols = 0;
        int notModMols = 0;
        int checkMols = 0;
        String checkfile = "check-AnalyzeAndFix_"+jobName+".sdf";
        String outSDfile = "3DreadyFixed_"+jobName+".sdf";

        //Iterate over structures
        try {
            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(struct3dFile), DefaultChemObjectBuilder.getInstance());
            //Analysis for each molecule
            while (reader.hasNext())
            {
                totMols++;
                IAtomContainer mol = reader.next();
                if (repOnScreen >= 1)
                    System.out.println("Working on mol "+totMols+" -> "+MolecularUtils.getNameOrID(mol));

                //Deal with aromatic bond and aromaticity
                String aromErr = MolecularUtils.missmatchingAromaticity(mol);
                if (!aromErr.equals(""))
                {
                    //Try to fix it
                    fixedMols++;
                    mol = MolecularUtils.fixSDAromaticity(mol);

                    //Double check modified molecule
                    String checkArom = MolecularUtils.missmatchingAromaticity(mol);
                    if (!checkArom.equals(""))
                    {
                        checkMols++;
                        if (repOnScreen >= 1)
                            System.out.println("Check this molecule! "+checkArom);

                        mol.setProperty("REJECTED",checkArom);
                        IOtools.writeSDFAppend(checkfile,mol,true);
                        continue;
                    }
                } else {
                    notModMols++;
                }

                //Write out the surviving molecule
                IOtools.writeSDFAppend(outSDfile,mol,true);
            }
        } catch (FileNotFoundException fnf) {
            System.err.println("File Not Found: " + struct3dFile);
            System.err.println(fnf.getMessage());
            System.exit(-1);
        } catch (Throwable t) {
            System.err.println("\nERROR in reading the file with IteratorMDLReader (in fixChemicalRepresentation). "+t);
            t.printStackTrace();
            System.out.println("\nERROR: cannot iterate through MDL file. This "
                                + "might be becouse the MDL format is prior to "
                                + "V2000 or because there are bond with "
                                + "type='any' (type=8 in connectivity matrix). "                                + "Please, try to provide a V2000 format and "
                                + "convert avoid bond type=8.");
            System.exit(0);
        }

        //Redirect name of input SDF file
        updateFileName(outSDfile);

        //Print out summary of results
        System.out.println("\n========== Results - Analyle & Fix ==========");
        System.out.println("Total number of molecules: "+totMols);
        System.out.println("Unmodified molecules:      "+notModMols);
        System.out.println("Modified molecules:        "+(fixedMols-checkMols));

        System.out.println("Output file: "+outSDfile);
        if (checkMols > 0)
            System.out.println("\nCheck "+checkMols+" problematic cases in "+checkfile);
    }

//------------------------------------------------------------------------------
    /**
     * Checks the agreement between the SDF molecular representation 
     * (as SDF file) and the molecular formula declared in a txt file
     * Compares the number of atoms per each element in the two 
     * representations identifying possible disagreements,
     * particularry suitable to identify missing atoms in SDF atoms.
     */

    public void check3DvsFormula()
    {
        if (repOnScreen >= 0)
        {
            String header = "\n========== check3DvsFormula starts ==========";
            System.out.println(header);
        }


        //Read and store chemical formula
        Map<String,String> formulae = null;
        formulae = getFormulae(formulaeFile);
        int numForms = formulae.keySet().size();
        if (repOnScreen >= 1)
            System.out.println("Found "+numForms+" chemical formulae in file "+formulaeFile);
 
        //Iterate over structures
        int num3ds = 0;
        int okMols = 0;
        boolean badStructures = false;
        String checkfile = "check-3DvsFORMULA_"+jobName+".sdf";
        String outSDfile = "3Dready_"+jobName+".sdf";
        Set<String> bad3d = new HashSet<String>();
        List<Boolean> copyTXTMol = new ArrayList<Boolean>();
        try {
            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(struct3dFile), DefaultChemObjectBuilder.getInstance());
            //Element analysis for each molecule
            while (reader.hasNext())
            {
                num3ds++;
                IAtomContainer mol = reader.next();

                //Skip if Formula is missing
                String refcode = mol.getProperty("cdk:Title").toString();
                if (!formulae.containsKey(refcode))
		{
		    if (repOnScreen >= 0)
			System.out.println("\nSkipping object "+num3ds+": missing formula with jey '"+refcode+"'");
                    continue;
		}

                String oneForm = formulae.get(refcode);
                if (repOnScreen >= 1)
                    System.out.println("\n===== CCDC object "+num3ds+" => "+refcode+" =====");

                //Analysis and comparison of SD and Formula
                ElementalAnalyser elAl = new ElementalAnalyser(oneForm,mol);

                //Return problematic cases to the user with detailed information
                if (!elAl.getAllElementsAgree())
                {
                    badStructures = true;
                    mol.setProperty("Formula",oneForm);
                    mol.setProperty("ElementAnalysisOverview",elAl.getSDFormulaAgreeOnElement());
                    mol.setProperty("fromFormula",elAl.getElemAnalFromFormula());
                    mol.setProperty("fromSD",elAl.getElemAnalFromSD());
                    IOtools.writeSDFAppend(checkfile,mol,true);
                    copyTXTMol.add(false);
                } else {
//                    mol.setProperty("Formula",oneForm);
//                    IOtools.writeSDFAppend(outSDfile,mol,true);
                    copyTXTMol.add(true);
                    okMols++;
                }
            }

            reader.close();

            //Make a copy of input filtering the molecules
            //Here we work with TXT to avoid CDK limitations on dealing with bond order = 4
            if (repOnScreen >= 1)
                System.out.println("\nCopying SD file as TXT. Please wait...");
            IOtools.filterSDasTXT(struct3dFile,outSDfile,copyTXTMol);
        } catch (FileNotFoundException fnf) {
            System.err.println("File Not Found: " + struct3dFile);
            System.err.println(fnf.getMessage());
            System.exit(-1);
        } catch (Throwable t) {
            System.err.println("\nERROR in reading the file with IteratorMDLReader. "+t);
            t.printStackTrace();
	    System.out.println("\nERROR: cannot iterate through MDL file. This "
	    			+ "might be becouse the MDL format is prior to "
				+ "V2000 or because there are bond with "
				+ "type='any' (type=8 in connectivity matrix). "
				+ "Please, try to provide a V2000 format and "
				+ "convert avoid bond type=8.");
            System.exit(0);
        }

        //Redirect name of input SDF file
        updateFileName(outSDfile);

        //Print out summary of results
        System.out.println("\n======== Results - check3DvsFormula =========");
        System.out.println("Total number of processed molecules:           "+num3ds);
        System.out.println("Molecules with all atoms of molecular formula: "+okMols);
        if (numForms != num3ds)
            System.out.println("\nWARNING!\nInput files do not contain the same number of chemical objects.\nThere are "+Math.abs(numForms-num3ds)+" objects that do not have structural information (in the SD file) or formula (in the TXT file).\n");
        if (badStructures)
           System.out.println("\nCheck the problematic cases in file "+checkfile);
    }

//------------------------------------------------------------------------------

    /**
     * Extract reference code and chemical formula for all the entries
     * in the txt file
     * @param intxt txt file with chemical informations generated by 
     * ConQuest-1.14 (CCDC)
     * @return the list of Reference Codes with the declared formula
     */

    private Map<String,String> getFormulae(String filename)
    {
        Map<String,String> allFormulae = new HashMap<String,String>();
        BufferedReader buffRead = null;
        try {
            //Read the file line by line
            buffRead = new BufferedReader(new FileReader(filename));
            String lineAll = null;
            String refcode = "";
            String formula = "";
            while ((lineAll = buffRead.readLine()) != null) 
            {
                String[] lineArgs = lineAll.split(":");
                //Get the name
                if (lineArgs[0].equals("REFCODE")) 
                 refcode = lineArgs[1].trim();

                //Get the formula
                if (lineArgs[0].equals("  Formula")) 
                {
                    formula = lineArgs[1].trim();
                    //Store formula
                    allFormulae.put(refcode,formula);
                    //Clean fields
                    refcode = "";
                    formula = "";
                }
            }
        } catch (FileNotFoundException fnf) {
            System.err.println("File Not Found: " + filename);
            System.err.println(fnf.getMessage());
            System.exit(-1);
        } catch (IOException ioex) {
            System.err.println(ioex.getMessage());
            System.exit(-1);           
        } finally {
            try {
                if (buffRead != null)
                    buffRead.close();
            } catch (IOException ioex2) {
                System.err.println(ioex2.getMessage());
                System.exit(-1);
            }
        }

        //return the map with refcodes and formulae
        return allFormulae;
   }

//------------------------------------------------------------------------------
}
