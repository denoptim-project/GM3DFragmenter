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
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.SortedMap;
import java.util.Iterator;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.smsd.Isomorphism;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.graph.PathTools;

import org.openscience.cdk.io.iterator.IteratingMDLReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;

import org.openscience.cdk.interfaces.ISingleElectron;
import org.openscience.cdk.interfaces.IElectronContainer;
import org.openscience.cdk.interfaces.ILonePair;
import org.openscience.cdk.config.IsotopeFactory;

/**
 * CLASS-based fragmentation of 3D structured 
 * chemical objects. Input and parameters are taken from 
 * <code>Parameters</code>
 * 
 * @author Marco Foscato (University of Bergen)
 */

public class Fragmenter
{
    //Job ID
    private String jobName;

    //Input file: molecular structures to chop
    private String inFile;

    //Output file: fragments
    private String outFile;

    //Output file: rejected molecules
    private String checkfile;
    private boolean someMolRejected;
    private int numRejected;
    private boolean thereAre2D;
    private int num2D;

    //Output file: number of molecules over number of fragments 
    private String MlFrRatioFile;

    //Fragments library formats
    private String inFormat;
    private String outFormat;

    //Compatibility matrix
    private CompatibilityMatrix compMat = new CompatibilityMatrix();

    //Output file for Compatibility matrix
    private String compMatFile;

    //Bond order
    private Map<String,Integer> classBndOrd = new HashMap<String,Integer>();

    //Storage of cutting rules 
    private Map<String,GM3DCuttingRule> cutRules = new HashMap<String,GM3DCuttingRule>();
    // Sorted list of cutting rules - Sorting based on the 'priority' filed of cutRules
    private SortedMap<Integer,String> sortedCutRules = new TreeMap();
    private List<String> any = new ArrayList<String>();

    //Rejection rules 
    // -> fragments containing AP with these CLASS+SubCLASS
    private Set<String> rejClasses = new HashSet<String>();
    // -> fragments containing any combination of APs with these CLASS or portion of CLASS' name
    // NB: only the combinations are detected. Fragments containing one AP of one of these these classes but not combined with nothet one are not rejected
    public static Set<List<String>> rejClassCombination = new HashSet<List<String>>();
    // -> fragments matching these SMARTS
    private static Set<String> rejSMARTS = new HashSet<String>();
    // -> fragments containing these elements
    private static Set<String> rejElements = new HashSet<String>();
    // -> fragments having more/less atoms than this
    private static int fragmentMaxSize = 30;
    private static int fragmentMinSize = 0;
    // -> fragments containing uncommon isotopes
    private static boolean rejIsotopes = true;

    //Remove duplicates
    private boolean removeDuplicates = true;

    //Ignore well known fragments
    private boolean ignoreKnownFrags = false;
    private String ignorableFile;
    private String ignorableFormat;

    //Collect only target fragment 
    private boolean lookForTargets = false;
    private String targetFile;
    private String targetFormat;
    private String fragCollectingDir;

    //Add dummy aton on linear system
    private boolean addDuOnLinear = false;
    private String duSymbol;

    //Execution mode - HIGH MEMORY MAY BE REQUIRED!
    private boolean highMemProfile = true;
//TODO the high memory profile is not available anymore!
//TODO remove the flags to run high Memomry profile

    //Reporting flag
    private int repOnScreen;

    //Error messages
    private String errMess="";

    //Recursion flag
    private int recNum = 1;

    //preString for reporting on screen
    private String pre = "Fragmenter: ";

//----------------------------------------------------------------------------------------
    /**
     * Creates a Fragmenter importing all input parameters from 
     * a dedicated object <code>Parameters</code>
     */
    public Fragmenter()
    {
        //Set all the specifications
        // -> CUTTING RULES
        cutRules = Parameters.rules;
        sortedCutRules = Parameters.sortedRules;
        any = Parameters.anyAtm; 
        // -> amount of writing on screen
        repOnScreen = Parameters.report;
        // -> rejection rules 
        rejClasses =  Parameters.rejClasses;
        rejClassCombination =  Parameters.rejClassCombination;
        rejSMARTS = Parameters.rejSMARTS;
        rejElements = Parameters.rejElements;
        fragmentMaxSize = Parameters.fragmentMaxSize;
        fragmentMinSize = Parameters.fragmentMinSize;
        rejIsotopes = Parameters.rejIsotopes;
        // -> remove duplicates
        removeDuplicates = Parameters.rmDuplicates;
	// -> ignore known fragments
	ignoreKnownFrags = Parameters.ignoreKnownFrags;
	ignorableFile = Parameters.ignorableFragsFile;
        ignorableFormat = Parameters.ignorableFragsFileFormat;
	// -> collect only target fragment 
	lookForTargets = Parameters.lookForTargets;
	targetFile = Parameters.targetFragsFile;
	targetFormat = Parameters.targetFragsFileFormat;

        // -> add Du on linear
        addDuOnLinear = Parameters.addDummyOnLinear;
	duSymbol = Parameters.duSymbol;

        inFormat = Parameters.getLibFormat();
        outFormat = Parameters.getLibFormat();

        //deal with file names
	jobName = Parameters.getJobName();
        inFile = Parameters.getCurrentSDFile();
        outFile = "Fragments_"+jobName+".sdf";
        compMatFile = "CPMap_"+jobName+".par";
        checkfile = "check-Fragmenter_"+jobName+".sdf";
        MlFrRatioFile = "MolFrag-ratio_"+jobName+".dat";
	fragCollectingDir = "fragsCollected_"+jobName;

	//Preparation of folder tree
        if (lookForTargets)
        {
            FileUtils.makeDir(fragCollectingDir);
        }


        someMolRejected = false;
        numRejected = 0;
    }

//----------------------------------------------------------------------------------------
    /**
     * Fragmentation of molecules for the production of 
     * three-dimensionally structured molecular fragments. 
     * Identifis and cuts chemical bonds according to a list of 
     * cutting rules given as SMARTS.
     * (see: http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html).
     * Input parameters, structures and cutting rules are defined via
     * <code>Parameters</code>
     */

    public void chopMolecules()
    {
        if (repOnScreen >= 0)
            System.out.println("\n============ Fragmentation Starts ============");

        // loop over molecules
        int molnum = 0;
        int numTotFrag = 0;
        int num2D = 0;
        try {
            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(inFile), DefaultChemObjectBuilder.getInstance());
//            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(inFile),NoNotificationChemObjectBuilder.getInstance()); // returns ERROR! AtomContainer does'n looks like a fragment in 'DENOPTIM' format!
            while (reader.hasNext()) 
            {
                //Start working with the molecule
                IAtomContainer mol = reader.next();
                String name = MolecularUtils.getNameOrID(mol);        
                if (repOnScreen >= 1)
                    System.out.println("\nChopping Molecule "+name+" - "+mol.getAtomCount()+" atoms");

                //Check dimensionality of the objgct
                int dimensions = MolecularUtils.getDimensions(mol);
                boolean is2D = false;
                if (dimensions == 2)
                {
                    is2D = true;
                } else if (dimensions != 3) {
                    rejectMol(mol,"Unable to get coordinates (at least 2D) for some atom");
                    continue;
                }


                // Get list of matching atoms with classes
                Map<String, ArrayList<GM3DTargetBond>> matchingbonds = getMatchingBondsAllInOne(mol);

                if (matchingbonds.keySet().size() == 0)
                {
		    if (errMess.equals(""))
		    {
		        errMess = "WARNING! No match found for any of the cutting rules.";
		    }
                    rejectMol(mol,errMess);
		    errMess = "";
                    continue;
                }

                // Cut previously selected bonds
                for (int i : sortedCutRules.keySet()) 
                {
                    String ruleName = sortedCutRules.get(i);
                    GM3DCuttingRule rule = cutRules.get(ruleName);

                    // Skip unmatched rules
                    if (!matchingbonds.keySet().contains(ruleName))
                        continue;

                    if (repOnScreen >= 1)
                        System.out.println("Cutting bonds matching rule "+ ruleName);

                    for (GM3DTargetBond tb: matchingbonds.get(ruleName)) 
                    {
                        IAtom atmA = tb.getAtmSubClass0();
                        IAtom atmB = tb.getAtmSubClass1();

                        //ignore if bond already broken
                        if (!mol.getConnectedAtomsList(atmA).contains(atmB))
                        { 
                            continue;
                        }

                        //treatment of n-hapto ligands
                        if (rule.isHAPTO())
                        {
                            if (repOnScreen >= 1)
                                System.out.println("Attempt to generate ligand with hapticity > 1");

                            //Get central atom index for the selected pair
                            //As a convention the central atom has subclass '0'
                            IAtom centralAtm = atmA;

                            //Get list of candidates for hapto: same cutting Rule and central metal
                            ArrayList<Integer> candidatesForHapto = new ArrayList<Integer>();
                            for (GM3DTargetBond tbForHapto : matchingbonds.get(ruleName))
                            {
                                //Consider only bond invorving same metal
                                if (tbForHapto.getAtmSubClass0() == centralAtm)
                                    candidatesForHapto.add(tbForHapto.getIDSubClass1());
                            }

                            //Get atms in n-hapto system: contiguous neigbours with same 
                            // type of bond with the same central atom
                            ArrayList<Boolean> flags = new ArrayList<Boolean>();
                            for (int iatm=0; iatm<mol.getAtomCount(); iatm++)
                                flags.add(false);
                            Set<Integer> atmsInHapto = new HashSet<Integer>();
                            atmsInHapto.add(tb.getIDSubClass1());
                            atmsInHapto = exploreHapticity(tb.getIDSubClass1(),candidatesForHapto,atmsInHapto,mol,flags);
                            if (atmsInHapto.size() == 1)
                            {
                                if (repOnScreen > 2)
                                    System.out.println("Unable to find other bonds involved in high-hapticity ligand! Bond ignored.");
                                continue;
                            }

                            if (repOnScreen >= 1)
                                System.out.println("Hapticity: "+atmsInHapto.size()+" - Atoms involved: "+atmsInHapto);

                            //Check existence of all bonds involved in multihapto system
                            boolean isSystemIntact = true;
                            for (int ligIdx : atmsInHapto)
                            {
                                IAtom ligAtm = mol.getAtom(ligIdx);
                                List<IAtom> nbrsOfLigAtm = mol.getConnectedAtomsList(ligAtm);

//TODO remove
//System.out.println("Controlling bond: "+ligAtm.getSymbol()+mol.getAtomNumber(ligAtm)+" "+centralAtm.getSymbol()+mol.getAtomNumber(centralAtm)+" = "+(!nbrsOfLigAtm.contains(centralAtm)));

                                if (!nbrsOfLigAtm.contains(centralAtm))
                                {
                                    isSystemIntact = false;
                                    break;
                                }
                            } 

                            //If not, it means that another rule already acted on the system 
                            //thus kill this attempt without generating dummyatom
                            if (!isSystemIntact)
                                continue;

                            //A dummy atom will be used to define attachment point of
                            // ligand with high hapticity
                            Point3d dummyP3d = new Point3d(); //Used also for reporting 2D
                            int pointsIn3D = 0;
                            for (int ligIdx : atmsInHapto)
                            {
                                Point3d ligP3d = MolecularUtils.getCoords3d(mol.getAtom(ligIdx));
                                dummyP3d.x = dummyP3d.x + ligP3d.x;
                                dummyP3d.y = dummyP3d.y + ligP3d.y;
                                dummyP3d.z = dummyP3d.z + ligP3d.z;
                            }

                            dummyP3d.x = dummyP3d.x / (double) atmsInHapto.size();
                            dummyP3d.y = dummyP3d.y / (double) atmsInHapto.size();
                            dummyP3d.z = dummyP3d.z / (double) atmsInHapto.size();

                            //Add Dummy atom to molecular object
                            //if no other Du is already in the same position
                            boolean addDu = true;
                            int oldDuIdx = -1;
                            for (IAtom oldDu : mol.atoms())
                            {
                                if (oldDu.getSymbol() == duSymbol)
                                {
                                    Point3d oldDuP3d = oldDu.getPoint3d();
                                    if (oldDuP3d.distance(dummyP3d) < 0.002)
                                    {
                                        addDu = false;
                                        oldDuIdx = mol.getAtomNumber(oldDu);
                                        break;
                                    }
                                } 
                            }
                        
                            IAtom dummyAtm;
                            if (addDu)
                            {
//TODO move Du to variable
                                dummyAtm = new Atom(duSymbol);
                                dummyAtm.setPoint3d(dummyP3d);
                                mol.addAtom(dummyAtm);
//TODO: remove                        MolecularUtils.addDummyOnCentralAtom(dummyAtm,mol,true);
                            } else {
                                dummyAtm = mol.getAtom(oldDuIdx);
                            }

                            //Modify connectivity of atoms involved in high-hapticity coordination
                            //creation of Du-to-ATM bonds
                            IBond.Order border = IBond.Order.valueOf("SINGLE");

                            // Translate bond order to int
                            int intborder = MolecularUtils.bondorderToint(border);
                            
                            for (int ligIdx : atmsInHapto)
                            {
                                IAtom ligAtm = mol.getAtom(ligIdx);
                                //Check for existing bond
                                List<IAtom> nbrsOfDu = mol.getConnectedAtomsList(dummyAtm);
//System.out.println("\nnbrsOfDu size: "+nbrsOfDu.size()+" "+ligIdx);
//for (IAtom aaa : nbrsOfDu)
//    System.out.println("nbrsOfDu: "+mol.getAtomNumber(aaa));
//System.out.println("T/F: "+nbrsOfDu.contains(ligAtm));

                                if (!nbrsOfDu.contains(ligAtm))
                                {
                                    //Add bond with dummy
                                    Bond bnd = new Bond(dummyAtm,ligAtm,border);
                                    mol.addBond(bnd);
                                }
                                //Remove bonds between coordinating atoms and central atom
                                IBond oldBnd = mol.getBond(centralAtm,ligAtm);
                                mol.removeBond(oldBnd);
                            }

                            if (repOnScreen >= 1)
                                System.out.println("Multihapto ligand created");

			    String locRuleName = ruleName;
                            String locCSC0 = tb.getClassSubClass0();
                            String locCSC1 = tb.getClassSubClass1();
                            if (Parameters.addIDtoAPClass)
			    {
				locRuleName = ruleName 
			        + Integer.toString(UniqueIndex.getUnqInt()) 
				+ "unq";
				locCSC0.replaceAll(ruleName,
							     locRuleName);
				locCSC1.replaceAll(ruleName,
							     locRuleName);
			    }

                            // Add rule to compatibility matrix
                            compMat.addTrueEntry(locCSC0,locCSC1);
                            compMat.addTrueEntry(locCSC1,locCSC0);

                            // Report CLASS-BondOrder key
                            if (!classBndOrd.keySet().contains(locRuleName))
                                classBndOrd.put(locRuleName,intborder);

                            //Set attachment points on central atom 
                            addAttachmentPoint(centralAtm,dummyAtm,locRuleName,
                                            getSubClass(ruleName,centralAtm),
                                            intborder,mol);
                            addAttachmentPoint(dummyAtm,centralAtm,locRuleName,
                                            getSubClass(ruleName,atmB),
                                            intborder,mol);
                        } else {
                            //tratment of mono-hapto ligands
                            if (repOnScreen >= 1)
                                System.out.println("Cutting bond between atoms "+tb.getIDSubClass0()+"-"+tb.getIDSubClass1());

                            //Identify object Bond
                            IBond bnd = mol.getBond(atmA,atmB);

                            // Remember bond order
//TODO make this a parameter that the user can change
                            IBond.Order border = IBond.Order.valueOf("SINGLE");
                            if (!rule.getSMARTSBnd().equals("~"))
                                border = bnd.getOrder();

                            // Translate bond order to int
                            int intborder = MolecularUtils.bondorderToint(border); 
                            String locRuleName = ruleName;
                            String locCSC0 = tb.getClassSubClass0();
                            String locCSC1 = tb.getClassSubClass1();
                            if (Parameters.addIDtoAPClass)
                            {
                                locRuleName = ruleName
                                + Integer.toString(UniqueIndex.getUnqInt())
                                + "unq";
                                locCSC0.replaceAll(ruleName,
                                                             locRuleName);
                                locCSC1.replaceAll(ruleName,
                                                             locRuleName);
                            }
                            
                            // Add rule to compatibility matrix
                            compMat.addTrueEntry(locCSC0,locCSC1);
                            compMat.addTrueEntry(locCSC1,locCSC0);

                            // Report CLASS-BondOrder key
                            if (!classBndOrd.keySet().contains(locRuleName))
                                classBndOrd.put(locRuleName,intborder);

                            // now cut the bond
                            mol.removeBond(bnd);

                            //Set attachment points on central atom 
                            addAttachmentPoint(atmA,atmB,locRuleName,
                                            getSubClass(ruleName,atmA),
                                            intborder,mol);
                            addAttachmentPoint(atmB,atmA,locRuleName,
                                            getSubClass(ruleName,atmB),
                                            intborder,mol);
                        } //end of if (hapticity>1)
                    } //end of loop over matching bonds
                } //end of loop over rules

                //Report eventual 2D
                if (is2D)
                {
                    thereAre2D = true;
                    num2D++;
                }

                //Isolation and analysis of the fragments
                try {
                    // Split the broken molecule producing all fragments
                    AtomContainerSet frags = isolateFrags(mol);

                    // Analyze the fragments
                    int i = 0;
                    for (IAtomContainer atCont : frags.atomContainers()) 
                    {
                        i++;
                        if (repOnScreen >= 3)
                             System.out.println("Analysing fragment "+ i);

                        GM3DFragment frag = new GM3DFragment(atCont);

                        //Get rid of fragments having no attachment point
                        //they cannot become GM3DFragments
                        if (frag.getNumberOfAttachmentPoints() == 0)
                        {
                            if (repOnScreen >= 3)
                                System.out.println("Molecule "+i+" has no attachment point! Discharged!");
                            continue;
                        }

                        //Set title to the fragment
                        frag.setProperty("cdk:Title","From_"+name+"_"+i);
                        //Change remarks
                        if (is2D)
                            frag.setProperty("cdk:Remark","From GM3DFragmenter - 2D");
                        else
                            frag.setProperty("cdk:Remark","From GM3DFragmenter");

                        //Add fragment to the output library
                        if (repOnScreen >= 1)
                            System.out.print("Attempt to add a new fragment... ");

                        //Add dummy atoms for use of internal coordinates
                        if ((!is2D) && addDuOnLinear)
                            MolecularUtils.addDummiesOnLinearities(frag);

                        //Check this fragments for rejection criteria
                        if (FragmentFilter.keepFragment(frag))
                        {
			    //Compare with list of frags to ignore
			    if (ignoreKnownFrags)
			    {
				if (hitIgnorableFragment(frag,ignorableFile,ignorableFormat))
                                {
                                    if (repOnScreen >= 1)
                                        System.out.println("Ignorable fragment.");
				    continue;
                                }
			    }
                            if (removeDuplicates)
                            {
                                //Compare frag with the alreagy generated frags
                                if (newFragment(frag,outFile))
                                {
                                    if (repOnScreen >= 1)
                                        System.out.println("NEW Fragment added!");
                                    IAtomContainer ac = frag.toIAtomContainer(outFormat);
                                    IOtools.writeSDFAppend(outFile, ac, true);
                                    numTotFrag++;
                                } else {
                                    if (repOnScreen >= 1)
                                        System.out.println("Not a new fragment.");
                                }
                            } else
                            {
				if (lookForTargets)
				{
				    //Compare frag with the library of targets
				    if (hitTargetFragment(frag,targetFile,targetFormat))
				    {
					String hit = frag.getProperty("TARGETHIT").toString();
					String fragFile = fragCollectingDir+"/"+"hittingTarget_"+hit+".sdf";
					IAtomContainer ac = frag.toIAtomContainer(outFormat);
					IOtools.writeSDFAppend(fragFile, ac, true);
					numTotFrag++;
				    }
				} else {
                                    if (repOnScreen >= 1)
                                        System.out.println("KEEP-FRAGMENTS MODE: Fragment added to the ouput list");
                                    IAtomContainer ac = frag.toIAtomContainer(outFormat);
                                    IOtools.writeSDFAppend(outFile, ac, true);
                                    numTotFrag++;
                                }
			    }
                        }
                    } //end loop over fragments

                    molnum++;

                    if (removeDuplicates)
                           reportMolFragRatio(molnum,numTotFrag);

                } catch (Throwable t) {
                    System.err.println("\nWARNING! Brutal exit for molecule "+name);
                    System.err.println(" EXCEPTION: "+t);
                    t.printStackTrace();
                }
            } //end of loop over molecules
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

        if (repOnScreen >= 0)
        {
            System.out.println("\n============= Fragmentation DONE =============");
            System.out.println("Total number of fragmented molecules: "+molnum);
            System.out.println("Total number of stored fragments:     "+numTotFrag);
            if (someMolRejected)
                System.out.println("\nCheck "+numRejected+" rejected molecules in "+checkfile);
            if (thereAre2D)
                System.out.println("\nFound "+num2D+" molecules in 2D. Check fragments labeled with '2D'");
        }

        //In case of no fragments generated
        if (numTotFrag == 0)
        {
            System.out.println("\nNO FRAGMENT was generated!\n");
            System.exit(0);
        }

        //Report Compatibility matrix for this fragmentation
        compMat.writeCPMapFile(compMatFile,classBndOrd);

        //Redirect input file of the next step
        Parameters.updateStructureFilePointer(outFile);
    }

//-----------------------------------------------------------------------------

    /**
     * Create a new attachment point to be reported as property of the source
     * atom in the molecular object. The information reported contains useful
     * description of the torget bond (type, direction, bond order)
     * @param fromAtm source atom and tail of the attachment point vector
     * @param toAtm tartget atom and head of the attachment point vector
     * @param rule name of the cutting rule used
     * @param subClass code of the sub class (0 or 1) identifying which of the 
     * two extremities of the target bond is represented by this attachment point
     * @param bndOrd bond order of the target bond
     * @param mol the whole molecular object 
     */

    private void addAttachmentPoint(IAtom fromAtm, IAtom toAtm, String rule, int subClass, int bndOrder, IAtomContainer mol)
    {
        //Identify atom
        int atmid = mol.getAtomNumber(fromAtm);

        //Define direction vector of AP
        ArrayList<Double> pointer = new ArrayList<Double>();
        Point3d p3d = MolecularUtils.getCoords3d(toAtm);
        pointer.add(p3d.x);
        pointer.add(p3d.y);
        pointer.add(p3d.z);

        //Create object attachment point
        GM3DAttachmentPoint ap = new GM3DAttachmentPoint(atmid,
                                                          rule,
                                                      subClass,
                                                      bndOrder,
                                                       pointer);

        //Add attachment point to the atom in properties
        try {
            ArrayList<GM3DAttachmentPoint> oldAPs = (ArrayList<GM3DAttachmentPoint>) fromAtm.getProperty("AttachmentPoints");
            oldAPs.add(ap);
            fromAtm.setProperty("AttachmentPoints",oldAPs);
        } catch (Throwable t ) {
            ArrayList<GM3DAttachmentPoint> aps = new ArrayList<GM3DAttachmentPoint>();
            aps.add(ap);
            fromAtm.setProperty("AttachmentPoints",aps);
        }
    }

//-----------------------------------------------------------------------------

    /**
     * Creates the list of pairs (rulename;SMARTSquery) from the list of 
     * cutting rules. The output list includes the whole cutting rule
     * [<atm0>]-[<atm1>], the first query [<atm0>], and the 
     * second query [<atm1>]
     */
    private Map<String,String> getAllSMARTSQueries()
    {
        Map<String,String> smarts = new HashMap<String,String>();

        for (int ir : sortedCutRules.keySet())
        {
            String ruleName = sortedCutRules.get(ir);
            GM3DCuttingRule rule = cutRules.get(ruleName);

            // add the whole cutting rule to the map
            smarts.put(ruleName,rule.getWholeSMARTSRule());

            // add the first atom query of this rule
            smarts.put(rule.getSubClassName0(),rule.getSMARTSSubClass0());

            // add the second atom query of this rule
            smarts.put(rule.getSubClassName1(),rule.getSMARTSSubClass1());
        }

        return smarts;
    }

//-----------------------------------------------------------------------------

    /**
     * Identification of the bonds matching a list of SMARTS rules
     * @param mol chemical system to be analyzed
     * @return list or couple of atoms (as integer idexes) per each rule name
     */

    private Map<String, ArrayList<GM3DTargetBond>> getMatchingBondsAllInOne(IAtomContainer mol)
    {
	Map<String, ArrayList<GM3DTargetBond>> matchingBonds = new HashMap<String, ArrayList<GM3DTargetBond>>();

        // Get all SMARTS queries
        Map<String,String> smarts = getAllSMARTSQueries();

        // Get all the matches
        ManySMARTSQuery msq = new ManySMARTSQuery(mol,smarts);
        if (msq.hasProblems())
        {
            String cause = msq.getMessage();
	    errMess = cause;
	    return matchingBonds;
        }

        // Loop over cutting rules tryingto match bonds
        for (int ir : sortedCutRules.keySet())
        {
            String ruleName = sortedCutRules.get(ir);
            GM3DCuttingRule rule = cutRules.get(ruleName);

            if (msq.getNumMatchesOfQuery(ruleName) == 0)
            {
                continue;
            }
            if (repOnScreen >= 1)
                System.out.println("Rule '"+ruleName+"'\n - Mathces: " + msq.getNumMatchesOfQuery(ruleName));

            // Apply further options of cutting rule
            List<List<Integer>> purgedPairs = filterListOfMatches(msq.getMatchesOfSMARTS(ruleName),rule,mol);
            if (repOnScreen >= 1)
                System.out.println(" - Mathces (post-filtering): " + purgedPairs.size());

            if (purgedPairs.size() == 0)
            {
                continue;
            }

            // Evaluate subclass membership and eventually store target bonds
            ArrayList<GM3DTargetBond> ruledBonds = new ArrayList<GM3DTargetBond>();
            for (int i = 0; i < purgedPairs.size(); i++)
            {
                Map<Integer,List<Boolean>> subClassMembership = defineSubClasses(
                                                purgedPairs.get(i).get(0),
                                                purgedPairs.get(i).get(1),
                                                msq,
                                                rule);
		if (repOnScreen >= 2)
		    System.out.println(" - Evaluation of SubClass: "+subClassMembership);

                // Finally, if the bond matches the rule unambiguously, store it
                GM3DTargetBond tb = new GM3DTargetBond();
                if (evaluateSubClass(rule,mol,purgedPairs.get(i),subClassMembership,tb))
                {
                    // write CLASS and SubCLASS in atoms' properties (as preCLASS)
                    storePreClassOnAtoms(mol,tb);
                    // store target bond in the output list
                    ruledBonds.add(tb);
                }
            }

            if (!ruledBonds.isEmpty())
                matchingBonds.put(ruleName,ruledBonds);
        }

        return matchingBonds;
    }

//-----------------------------------------------------------------------------

    /**
     * Apply cutting rule options
     */
    private List<List<Integer>> filterListOfMatches(List<List<Integer>> inList, GM3DCuttingRule rule, IAtomContainer mol)
    {
        // temporary storage
        List<List<Integer>> workList = new ArrayList<List<Integer>>();
        for (List<Integer> l : inList)        
        {
            workList.add(l);
        }

        if (rule.hasOptions())
        {
	    //Count Du that may contribute to ring size
	    int numDu = 0;
	    for (IAtom a : mol.atoms())
	    {
		if (a.getSymbol().equals(Parameters.duSymbol) && 
		    mol.getConnectedAtomsCount(a)>1)
		{
		    numDu++;
		}
	    }

            ArrayList<String> opts = rule.getOptions();
            for (int wi=0; wi<opts.size(); wi++)
            {
                if (inList.size() == 0)
                {
                    break;
                }

                String opt = opts.get(wi);
                //Deal with extra ring size requirements
                if (opt.startsWith("RING>"))
                {
                    String[] parts = opt.split(">");
                    int sizeLim = Integer.parseInt(parts[1]);
                    if (repOnScreen >= 2)
                        System.out.println(" - Got smallest ring-size requirement (SIZE > "+sizeLim+")");
                    List<List<Integer>> purgedList = new ArrayList<List<Integer>>();
                    for (int pi=0; pi<workList.size(); pi++)
                    {
                        List<Integer> pair = workList.get(pi);
                        int ringSize = MolecularUtils.getSmallestRing(
                                        mol.getAtom(pair.get(0)),
                                        mol.getAtom(pair.get(1)),mol);
                        if ((ringSize == -1) || (ringSize > sizeLim))
                            purgedList.add(pair);                               
                    }
                    workList = purgedList;
                }

		if (opt.startsWith("OMRING>"))
                {
                    String[] parts = opt.split(">");
                    int sizeLim = Integer.parseInt(parts[1]);
                    if (repOnScreen >= 2)
                        System.out.println(" - Got smallest metal-involving ring-size requirement (SIZE > "+sizeLim+")");
                    List<List<Integer>> purgedList = new ArrayList<List<Integer>>();
                    for (int pi=0; pi<workList.size(); pi++)
                    {
                        List<Integer> pair = workList.get(pi);
                        int ringSize = MolecularUtils.getSmallestOMRing(
	                                        mol.getAtom(pair.get(0)),
	                                        mol.getAtom(pair.get(1)),
						mol,
						sizeLim + numDu);
                        if ((ringSize == -1) || (ringSize > sizeLim))
                            purgedList.add(pair);
                    }
                    workList = purgedList;
                }

                // Add other requirements here
//                if (opt.startsWith("???"))
//                {
//                }

            }
        }

        return workList;
    }
 
//-----------------------------------------------------------------------------

    /**
     * Identification of the bonds matching a list of SMARTS rules
     * @param mol chemical system to be analyzed
     * @return list or couple of atoms (as integer idexes) per each rule name
     */

/*    private Map<String, ArrayList<List<Integer>>> getMatchingBonds(IAtomContainer mol)
    {
        Map<String, ArrayList<List<Integer>>> matchingBonds = new HashMap<String, ArrayList<List<Integer>>>();

        // Loop over the cutting rules
        for (int ir : sortedCutRules.keySet()) 
        {
            String rule = sortedCutRules.get(ir);
            ArrayList<String> ruleTerms = cutRules.get(rule);
            String smirule = ruleTerms.get(0) + ruleTerms.get(2) + ruleTerms.get(1);
            if (repOnScreen >= 1)
            {
                System.out.print(" - Query: "+rule+" => "+smirule+" ");
                if (ruleTerms.size() > 4)
                {
                    for (int wi=4; wi<ruleTerms.size(); wi++)
                    {
                        System.out.print(ruleTerms.get(wi)+" ");
                    }
                }
                System.out.print("\n");
            }

            boolean found = false;
            List<List<Integer>> listOfCouples = new ArrayList<List<Integer>>();
            ArrayList<List<Integer>> ruledBonds = new ArrayList<List<Integer>>();
            try {
                SMARTSQueryTool query = new SMARTSQueryTool(smirule);
                if (query.matches(mol))
                    found = true;
                listOfCouples = query.getUniqueMatchingAtoms();
            } catch (CDKException cdkEx) {
                String cause = cdkEx.getCause().getMessage(); 
                if (cause.equals("Timeout for AllringsFinder exceeded"))
                {
                    if (repOnScreen >= 1)
                        System.out.println("\nWARNING! "+cause);
                    matchingBonds.put(rule,ruledBonds);
                    rejectMol(mol,cause);
                    break;
                } else {
                    if (repOnScreen >= 1)
                    {
                        System.out.println("\nWARNING! Not able to search SMARTS for rule "+rule);
                        System.out.println("Cause: "+cause);
                    }
                    matchingBonds.put(rule,ruledBonds);
                    rejectMol(mol,cause);
                    break;
*/
/*
// better not to go on: you rinsk other rules can match the bond leading to a messy library

                    continue; 
*/
/*                }
            } catch (Throwable t) {
                if (repOnScreen >= 1)
                {
                    System.out.println("\nWARNING! Not able to search SMARTS for rule "+rule);
                    System.out.println("Using SMARTSQueryTool: "+t);
                    t.printStackTrace();
                }
                matchingBonds.put(rule,ruledBonds);
                String cause = "Not able to search SMARTS for rule "+rule;
                rejectMol(mol,cause);
                break;
/*
// better not to go on: you rinsk other rules can match the bond leading to a messy library

                    continue; 
*/
/*            }

            // Apply further conditions to filter list matches
            if (found) 
            {
                if (cutRules.get(rule).size() > 4)
                {
                    for (int wi=4; wi<ruleTerms.size(); wi++)
                    {
                        if (listOfCouples.size() == 0)
                        {
                            found = false;
                            break;
                        }

                        String opt = cutRules.get(rule).get(wi);
                        //Deal with extra ring size requirements
                        if (opt.startsWith("RING>"))
                        {
                            String[] parts = opt.split(">");
                            int sizeLim = Integer.parseInt(parts[1]);
                            if (repOnScreen >= 2)
                                System.out.println("Got smallest ring-size requirement (SIZE > "+sizeLim+")");
                            List<List<Integer>> purgedList = new ArrayList<List<Integer>>();
                            for (int pi=0; pi<listOfCouples.size(); pi++)
                            {
                                List<Integer> pair = listOfCouples.get(pi);
                                int ringSize = MolecularUtils.getSmallestRing(
                                                mol.getAtom(pair.get(0)),
                                                mol.getAtom(pair.get(1)),mol);
                                if ((ringSize == -1) || (ringSize > sizeLim))
                                    purgedList.add(pair);                                
                            }
                            listOfCouples = purgedList;
                        }

                        // Add other requirements here
//                        if (opt.startsWith("???"))
//                        {
//                        }

                    }
                }
            }

            // Evaluate and store the remaining matches
            if (found)
            {
                int nmatch = listOfCouples.size();
                if (repOnScreen >= 1)
                    System.out.println(" - Mathces: " + nmatch);
                for (int i = 0; i < nmatch; i++) 
                {
                    if (repOnScreen >= 1)
                        System.out.println(" - - > "+i+" "+listOfCouples.get(i));

                    // Evaluate and eventually store (sub)CLASSes as property of the two Atoms
                    if (evaluateSubClass(mol,rule,listOfCouples.get(i)))
                        ruledBonds.add(listOfCouples.get(i));
                }
            } 

            // Store bonds per each class
            matchingBonds.put(rule,ruledBonds);
        }

        return matchingBonds;
   }
*/
//-----------------------------------------------------------------------------

    /**
     *TODO fix documentation
     * Definition of the SubClass of both atoms matching a cutting 
     * rule. Since a cutting rule is written as a couple of two
     * atoms having peculiar properties (hybridation,chemical 
     * environment,charge,bonds, etc.) each rule defines two 
     * SubClasses: one per each atom used to write the cutting rule.
     * This method identifies (if possible) to which SubClass each
     * of the two atoms belongs and store this information as atom'
     * propety 'preCLASS'. If it is not possible to define the 
     * SubClass properly, which means that we cannot identify 
     * 'which is which' in the couple of atoms, no SubClass is 
     * provided. 
     * @param mol molecular system containing the matching atoms 
     * @param rule name of the cutting rule for which the 
     * SubClass must be defined 
     * @param atimid couple of integer indexes identifying the atoms
     * matching the cutting rule 
     * @return <code>true</code> if the definition of the SubClass
     * succeeds
     */

    private boolean evaluateSubClass(GM3DCuttingRule rule, IAtomContainer mol, List<Integer> atmid, Map<Integer,List<Boolean>> subClassMembership, GM3DTargetBond tb)
    {
        //We need to take into account a list of possibilities 
        //{consider tha boolean vector [True/False for subrule0,True/False for subrule1]}
        //An atom may respect the criteria for 
        //              --> only one of the subclasses {[1,0] or [0,1]}
        //              --> both subclasses [1,1]
        //              --> none of the subclasses [0,0]
        //So having two atoms we need to take into account 16 possibilities
        //
        //                     a           b           c           d
        //            Atom0: [0,0]       [1,0]       [0,1]       [1,1]
        //      Atom1
        //  a   [0,0]     [0,0][0,0]  [0,0][1,0]  [0,0][0,1]  [0,0][1,1]  
        //  b   [1,0]     [1,0][0,0]  [1,0][1,0]  [1,0][0,1]  [1,0][1,1]  
        //  c   [0,1]     [0,1][0,0]  [0,1][1,0]  [0,1][0,1]  [0,1][1,1]  
        //  d   [1,1]     [1,1][0,0]  [1,1][1,0]  [1,1][0,1]  [1,1][1,1]  
        //
        // No possibility to define the CLASS for 6 out of 16 cases: aa, bb, cc, dd, da,ad
        // The very easy cases are (on atom respect ONLY one subclass): cb, bc
        // Other casess that, for some reason, lead to an atom not being described by any class 
        // (i.e. ester O-R by RECAP) are: ab, ac, ba, ca
        // For these (ab, ac, ba, ca) the the definition of the CLASS is possible since one
        // of the atoms belongs to one and only one subclass.
        // When a rule may be defined by the "any atom" definition for one of the atoms
        // identifining the bond to be broken (i.e. SMARTS = [...whatever...]!@[$(*)]), both
        // atoms will fall in this subclass. Like in the remaining cases: bd, db (subrule0 is 
        // always 1) and cd, dc (subrule1 is always 1). We wanna deal with these 4 cases only
        // if it is rally true that this is due to the definition of an "any-atom" rule.
        // A similar behaviour was detected in case of vicinal groups. This is treated as
        // a critical case below.
        //

        String smirule0 = rule.getSMARTSSubClass0();
        String smirule1 = rule.getSMARTSSubClass1();

        String sbcl0 = rule.getSubClassName0();
        String sbcl1 = rule.getSubClassName1();

        //'class0' referts to the class of atom of which the index is in atmid.get(0)
        String class0 = "";
        //'class2' referts to the class of atom of which the index is in atmid.get(1)
        String class1 = "";

        // Check for case 'aa' = [0,0][0,0]
        boolean allFalse = true;
        for (Integer i : subClassMembership.keySet())
        {
            for (Boolean entry : subClassMembership.get(i))
            {
                if (entry)
                    allFalse = false;
            }
        }
        // run away in case of no possibility to define subclasses (case 'aa' = [0,0][0,0])
        if (allFalse)
        {
            if (repOnScreen >= 3)
                printInfosDEBUG(subClassMembership);
            return false;        
        }

        // Let's see if we can work out the subclass membership
        boolean ignore = false;
        // first deal with any-atom matching rules (cases like da,db,dc,dd,ad,bd,cd)
        // 
        if (any.contains(smirule0)) {
            if (subClassMembership.get(1).get(1) && !subClassMembership.get(1).get(0)) {
                class0 = sbcl0;
                class1 = sbcl1;
                if (repOnScreen >= 3)
                    printInfosDEBUG(subClassMembership);
            } else if (!subClassMembership.get(1).get(1) && subClassMembership.get(1).get(0)) {
                class0 = sbcl1;
                class1 = sbcl0;
                if (repOnScreen >= 3)
                    printInfosDEBUG(subClassMembership);
            } else {
                if (repOnScreen >= 3)
                    System.out.println("GEN not found - "+subClassMembership);
                ignore = true;
            }
        } else if (any.contains(smirule1)) {
            if (subClassMembership.get(0).get(0) && !subClassMembership.get(0).get(1)) {
                class0 = sbcl0;
                class1 = sbcl1;
                if (repOnScreen >= 3)
                    printInfosDEBUG(subClassMembership);
            } else if (!subClassMembership.get(0).get(0) && subClassMembership.get(0).get(1)) {
                class0 = sbcl1;
                class1 = sbcl0;
                if (repOnScreen >= 3)
                    printInfosDEBUG(subClassMembership);
            } else {
                if (repOnScreen >= 3)
                    System.out.println("GEN not found - "+subClassMembership);
                ignore = true;
            }
        } else {
        // Deal with other cases not containing any-atom matching rules
            if ((subClassMembership.get(0).get(0) && !subClassMembership.get(1).get(0)) && 
                (subClassMembership.get(1).get(1) && !subClassMembership.get(0).get(1))) {
                class0 = sbcl0;
                class1 = sbcl1;
                if (repOnScreen >= 3)
                    printInfosDEBUG(subClassMembership);
            } else if ((subClassMembership.get(0).get(1) && !subClassMembership.get(1).get(1)) && 
                (subClassMembership.get(1).get(0) && !subClassMembership.get(0).get(0))) {
                class0 = sbcl1;
                class1 = sbcl0;
                if (repOnScreen >= 3)
                    printInfosDEBUG(subClassMembership);
            } else if (!subClassMembership.get(0).get(0) && !subClassMembership.get(1).get(0)) {
        //NB: Here we use the subclass from only atom '1' to decide upon subslasses
                if (subClassMembership.get(1).get(1) && !subClassMembership.get(0).get(1)) {
                    class0 = sbcl0;
                    class1 = sbcl1;
                    if (repOnScreen >= 3)
                        printInfosDEBUG(subClassMembership);
                } else if (subClassMembership.get(0).get(1) && !subClassMembership.get(1).get(1)) {
                    class0 = sbcl1;
                    class1 = sbcl0;
                    if (repOnScreen >= 3)
                        printInfosDEBUG(subClassMembership);
                }
            } else if (!subClassMembership.get(0).get(1) && !subClassMembership.get(1).get(1)) {
        //NB: Here we use the subclass from only atom '0' to decide upon subslasses
                if (subClassMembership.get(1).get(0) && !subClassMembership.get(0).get(0)) {
                    class0 = sbcl1;
                    class1 = sbcl0;
                    if (repOnScreen >= 3)
                        printInfosDEBUG(subClassMembership);
                } else if (subClassMembership.get(0).get(0) && !subClassMembership.get(1).get(0)) {
                    class0 = sbcl0;
                    class1 = sbcl1;
                    if (repOnScreen >= 3)
                        printInfosDEBUG(subClassMembership);
                }
            } else {
                ignore = true;
                if (repOnScreen >= 3) {
                    System.out.println("Case NOT COVERED for "+rule);
                    printInfosDEBUG(subClassMembership);
                }
            }
        }

        //Treatment of equal atom queries
        if (rule.isSymmetric()) {
            ignore = false;
            class0 = sbcl0;
            class1 = sbcl0;
        }

        //Try to recover the critical cases
        if (ignore)
        {
            int trueEntries = howManyTrueEntries(subClassMembership);

            if ((trueEntries == 0) || (trueEntries == 2) || (trueEntries == 4))
            {
                // Here we got an even number of true entries in subClassMembership

                //Treatment of what looks like a bug in CDK... 
                // the whole SMARTS query matches the bond but one
                // of the two 1-atom SMARTS doesn't match any atom

                if (repOnScreen >= 3) 
                {
                    System.out.println("Analysing critical case: ");
                    printInfosDEBUG(subClassMembership);
                }

                //Identify which one of the SMARTS to use

                //TODO now it's hard coded: first 0 and then 1: define criteria to choose order
                // so to avoid one of the two steps and speed up the execution

                //TODO check fo $(...) and !$(...)
                
                int res0 = matchSimplifiedSMARTS(mol,atmid.get(0),atmid.get(1),smirule0);
                if (res0 != -1)
                {
//System.out.println("matchSimplifiedSMARTS-0: "+res0);
                    if (res0 == atmid.get(0))
                    {
                        ignore = false;
                        class0 = sbcl0;
                        class1 = sbcl1;
                    } else if (res0 == atmid.get(1))
                    {
                        ignore = false;
                        class0 = sbcl1;
                        class1 = sbcl0;
                    } 
                } else {
                    int res1 = matchSimplifiedSMARTS(mol,atmid.get(0),atmid.get(1),smirule1);
//System.out.println("matchSimplifiedSMARTS-1: "+res1);
                    if (res1 != -1)
                    {
                        if (res1 == atmid.get(0))
                        {
                            ignore = false;
                            class0 = sbcl1;
                            class1 = sbcl0;
                        } else if (res1 == atmid.get(1))
                        {
                            ignore = false;
                            class0 = sbcl0;
                            class1 = sbcl1;
                        }
                    }
                }
            } else if ((trueEntries == 1) || (trueEntries == 3))
            {
                // Here we got an odd number of true entries in subClassMembership 
                // We can still get the subclass membership without further matching

                if (repOnScreen >= 3)
                {
                    System.out.println("Recovering odd case: ");
                    printInfosDEBUG(subClassMembership);
                }

                int atmMatched = -1;
                for (Integer i : subClassMembership.keySet())
                {
                    int numTrue = 0;
                    for (Boolean entry : subClassMembership.get(i))
                    {
                        if (entry)
                            numTrue++;
                    }

                    if (numTrue ==1)
                    {
                        atmMatched = i;
                        ignore = false;
                        break;
                    }
                }

                String recoveredClassOfMatched = "";
                String recoveredClassOfNotMatched = "";
                if (subClassMembership.get(atmMatched).get(1) && !subClassMembership.get(atmMatched).get(0)) 
                { 
                    // for 'atmMatched' we have [F,T]
                    recoveredClassOfMatched = sbcl1;
                    recoveredClassOfNotMatched = sbcl0;
  
                } else if (!subClassMembership.get(atmMatched).get(1) && subClassMembership.get(atmMatched).get(0)) 
                {
                    // for 'atmMatched' we have [T,F]
                    recoveredClassOfMatched = sbcl0;
                    recoveredClassOfNotMatched = sbcl1;
                }

                if (!ignore)
                {
                    if (atmMatched == 0)
                    {
                        class0 = recoveredClassOfMatched;
                        class1 = recoveredClassOfNotMatched;
                    } else if (atmMatched == 1)
                    {
                        class0 = recoveredClassOfNotMatched;
                        class1 = recoveredClassOfMatched;
                    }
                }
            }
        }
/*
//TODO: this is used for design of cutting rules. Make it dependen on input keyword

//to check which side of the SMARTS string is not working properly
System.out.println("class0: "+class0+" "+atmid.get(0));
System.out.println("class1: "+class1+" "+atmid.get(1));
printInfosDEBUG(subClassMembership);
System.out.println("HERE");
IOtools.pause();

//printInfosDEBUG(subClassMembership);
//endTODO
*/

        // Make the target bond
        boolean breakbond = false;
        if (!ignore) 
        {
            if (class0.equals(sbcl0) && class1.equals(sbcl1))
            {
                tb.setAtmSubClass0(mol,atmid.get(0));
                tb.setAtmSubClass1(mol,atmid.get(1));
                tb.setSubClass0(sbcl0);
                tb.setSubClass1(sbcl1);
                tb.setRuleName(rule.getName());
                tb.setSymmetricSubClass(rule.isSymmetric());            
            } else if (class0.equals(sbcl1) && class1.equals(sbcl0))
            {
                tb.setAtmSubClass0(mol,atmid.get(1));
                tb.setAtmSubClass1(mol,atmid.get(0));
                tb.setSubClass0(sbcl0);
                tb.setSubClass1(sbcl1);
                tb.setRuleName(rule.getName());
                tb.setSymmetricSubClass(rule.isSymmetric());    
            } else if (class0.equals(sbcl0) && class1.equals(sbcl0))
            {
                tb.setAtmSubClass0(mol,atmid.get(0));
                tb.setAtmSubClass1(mol,atmid.get(1));
                tb.setSubClass0(sbcl0);
                tb.setSubClass1(sbcl0);
                tb.setRuleName(rule.getName());
                tb.setSymmetricSubClass(rule.isSymmetric());    
            }
            breakbond = true;
        }

        return breakbond;
    }

//-----------------------------------------------------------------------------

    /**
     * write CLASS and SubCLASS in atoms' properties (as preCLASS)
     */
    private void storePreClassOnAtoms(IAtomContainer mol, GM3DTargetBond bnd)
    {
        //do it for atom in <subclass>0
        IAtom atm0 = mol.getAtom(bnd.getIDSubClass0());
        String old = "";
        try {
            old = atm0.getProperty("preCLASS").toString()+"_";
        } catch (Throwable t0) {
            old = "";
        }
        atm0.setProperty("preCLASS",old+bnd.getClassSubClass0());

        //do it for the other one in <subclass>1
        String old1 = "";
        IAtom atm1 = mol.getAtom(bnd.getIDSubClass1());
        try {
            old1 = atm1.getProperty("preCLASS").toString()+"_";
        } catch (Throwable t1) {
            old1 = "";
        }
        // 'bnd' already containg the proper subclass in case of symmetric cutting rule
        atm1.setProperty("preCLASS",old1+bnd.getClassSubClass1());

        if (repOnScreen >= 3) 
        {
            System.out.println("preCLASS-atm0 "+old+" "+Integer.toString(mol.getAtomNumber(atm0))+":"+bnd.getClassSubClass0());
            System.out.println("preCLASS-atm1 "+old1+" "+Integer.toString(mol.getAtomNumber(atm1))+":"+bnd.getClassSubClass1());
        }
    }

//-----------------------------------------------------------------------------

    /**
     * Create a 2x2 matrix-like map containing sublass membership of the two 
     * atoms given
     */
    private Map<Integer,List<Boolean>> defineSubClasses(int atm0, int atm1, ManySMARTSQuery msq, GM3DCuttingRule rule)
    {
        Map<Integer,List<Boolean>> subclassesVSatms = new HashMap<Integer,List<Boolean>>();

        String subclass0Name = rule.getSubClassName0(); 
        String subclass1Name = rule.getSubClassName1();

        boolean found0 = msq.hasMatches(subclass0Name);
        int numMatches0 = msq.getNumMatchesOfQuery(subclass0Name);

        boolean found1 = msq.hasMatches(subclass1Name);
        int numMatches1 = msq.getNumMatchesOfQuery(subclass1Name);

        List<List<Integer>> listOfCouples0 = new ArrayList<List<Integer>>();
        List<List<Integer>> listOfCouples1 = new ArrayList<List<Integer>>();

        // Which atom matches subclass 0?
        List<Boolean> subclass0 = new ArrayList<Boolean>(Arrays.asList(false,false));
        if (found0)
        {
            listOfCouples0 = msq.getMatchesOfSMARTS(subclass0Name);
            boolean doneAtm0 = false;
            boolean doneAtm1 = false;
            // for each group of atoms matching the SMARTS query 
            for (int im = 0; im < numMatches0; im++)
            {
                List<Integer> sublist = listOfCouples0.get(im);
                // There should be only one atom per each group but this is more general
                // For each index in the group
                for (int id : sublist)
                {
                    if (id == atm0)
                    {
                        subclass0.set(0,true);
                        doneAtm0 = true;
                    }
                    if (id == atm1)
                    {
                        subclass0.set(1,true);
                        doneAtm1 = true;
                    }
                }
                // just to speed up
                if (doneAtm0 & doneAtm1)
                    break;
            }
        }
        subclassesVSatms.put(0,subclass0);

        // Which atom matches subclass 1?
        List<Boolean> subclass1 = new ArrayList<Boolean>(Arrays.asList(false,false));
        if (found1)
        {
            listOfCouples1 = msq.getMatchesOfSMARTS(subclass1Name);
            boolean doneAtm0 = false;
            boolean doneAtm1 = false;
            // For each group of indeces
            for (int im = 0; im < numMatches1; im++)
            {
                List<Integer> sublist = listOfCouples1.get(im);

                // For each index in the group
                for (int id : sublist)
                {
                    if (id == atm0)
                    {
                        subclass1.set(0,true);
                        doneAtm0 = true;
                    }
                    if (id == atm1)
                    {
                        subclass1.set(1,true);
                        doneAtm1 = true;
                    }
                }
                if (doneAtm0 & doneAtm1)
                    break;
            }
        }
        subclassesVSatms.put(1,subclass1);
        
        return subclassesVSatms;
    }

//-----------------------------------------------------------------------------

/*    private Map<Integer,List<Boolean>> defineSubClasses(IAtomContainer mol, int atm0, int atm1, String smarts0, String smarts1)
    {
        //create the 4x4 matrix of subclass membership 
        //Key = subclass; Value = couple of booleans for Atom0 and Atom1
        Map<Integer,List<Boolean>> subclassesVSatms = new HashMap<Integer,List<Boolean>>();

        boolean found0 = false;
        int numMatches0 = -1;
        boolean found1 = false;
        int numMatches1 = -1;
        List<List<Integer>> listOfCouples0 = new ArrayList<List<Integer>>();
        List<List<Integer>> listOfCouples1 = new ArrayList<List<Integer>>();

        //get a smaller molecule for better performances
        int radius = 5;
        IAtomContainer smallerMol = new AtomContainer();

        try {
            smallerMol = (IAtomContainer) mol.clone();
            for (int i = 0; i < smallerMol.getAtomCount(); i++)
            {
                IAtom testAtm = smallerMol.getAtom(i);
                testAtm.setProperty("ORIGINAL_NUM",i);
                visitAtomsInRadius(smallerMol,atm0,radius,10);
                visitAtomsInRadius(smallerMol,atm1,radius,10);
            }
            //delete atoms that are too far
            for (IAtom delatm : smallerMol.atoms())
            {
                if (!delatm.getFlag(10))
                {
                    for (IBond delbond : smallerMol.getConnectedBondsList(delatm))
                        smallerMol.removeBond(delbond);
                    smallerMol.removeAtom(delatm);
                }
            }
        } catch (Throwable t) {
            if (repOnScreen >= 2)
                System.out.println("WARNING! Cloning failed. I'll try to go on anyway (but it will take longer)...");
            //try to proceed anyway... the more time consuming approach may stil works
            smallerMol = mol;
        }

        //Evaluate who belongs to which subclass 
        try {
            SMARTSQueryTool query0 = new SMARTSQueryTool(smarts0);
            SMARTSQueryTool query1 = new SMARTSQueryTool(smarts1);
            if (query0.matches(smallerMol)) {
                found0 = true;
                listOfCouples0 = query0.getUniqueMatchingAtoms();
                numMatches0 = listOfCouples0.size();
            }
            if (query1.matches(smallerMol)) {
                found1 = true;
                listOfCouples1 = query1.getUniqueMatchingAtoms();
                numMatches1 = listOfCouples1.size();
            }

//TODO remove
//System.out.println("0: "+query0.matches(smallerMol));
//System.out.println("1: "+query1.matches(smallerMol));

        } catch (Throwable t) {
            if (repOnScreen >= 1)
                System.out.println("ERROR! Troubles in getting the subclasses. Bond between "+atm0+" and "+atm1+" will ot be broken. "+t.getMessage());
            return subclassesVSatms;
        }
//TODO remove
//System.out.println("0: "+listOfCouples0);
//System.out.println("1: "+listOfCouples1);

        List<Boolean> subclass0 = new ArrayList<Boolean>(Arrays.asList(false,false));
        //fill the 4x4 matrix (Map) of subclass membership
        //subclass 0
        if (found0)
        {
            boolean doneAtm0 = false;
            boolean doneAtm1 = false;
            for (int im = 0; im < numMatches0; im++)
            {
                List<Integer> sublist = listOfCouples0.get(im);
                for (int newi : sublist)
                {
                    String ondiStr = smallerMol.getAtom(newi).getProperty("ORIGINAL_NUM").toString();
                    int oldi = Integer.parseInt(ondiStr);
                    if (oldi == atm0)
                    {
                        subclass0.set(0,true);
                        doneAtm0 = true;
                    }
                    if (oldi == atm1)
                    {
                        subclass0.set(1,true);
                        doneAtm1 = true;
                    }
                }
                if (doneAtm0 & doneAtm1)
                    break;
            }
        }
        subclassesVSatms.put(0,subclass0);
        // subclass 1
        List<Boolean> subclass1 = new ArrayList<Boolean>(Arrays.asList(false,false));
        if (found1)
        {
            boolean doneAtm0 = false;
            boolean doneAtm1 = false;
            for (int im = 0; im < numMatches1; im++)
            {
                List<Integer> sublist = listOfCouples1.get(im);
                for (int newi : sublist)
                {
                    String ondiStr = smallerMol.getAtom(newi).getProperty("ORIGINAL_NUM").toString();
                    int oldi = Integer.parseInt(ondiStr);
                    if (oldi == atm0)
                    {
                        subclass1.set(0,true);
                        doneAtm0 = true;
                    }
                    if (oldi == atm1)
                    {
                        subclass1.set(1,true);
                        doneAtm1 = true;
                    }
                }
                if (doneAtm0 & doneAtm1)
                    break;
            }
        }
        subclassesVSatms.put(1,subclass1);
        return subclassesVSatms;
    }

*/
//-----------------------------------------------------------------------------

    /**
     * Try to define subClass membership by recursive simplification of the
     * SMARTS string
     * @param mol the molecule
     * @param atm0 atom in first position in pair of matches
     * @param atm1 atom in second position in pair of matches
     * @param smarts SMARTS string to be simplified
     * @return the index of the atom matched  or -1 in case of no match
     */

    private int matchSimplifiedSMARTS(IAtomContainer mol, int atm0, int atm1, String smarts)
    {

        GM3DSMARTS s = new GM3DSMARTS(smarts);
        if (repOnScreen >= 2)        
        System.out.println("Trying to semplify GM3DSMARTS: "+smarts.toString());

        int res = -1;
        Map<String,String> simplerSmarts = new HashMap<String,String>();
        for (int i=0; i<s.getMaxNumSemplification(); i++)
        {
            String simpler = s.getSimplerSMARTS();
            simplerSmarts.put(Integer.toString(i),simpler);
        }

        // Get all the matches
        ManySMARTSQuery msq = new ManySMARTSQuery(mol,simplerSmarts);

        if (msq.hasProblems())
        {
            String cause = msq.getMessage();
            rejectMol(mol,cause);
        }

        for (int i=0; i<s.getMaxNumSemplification(); i++)
        {
            // Get matches for the simplitied queri corresponding to 'i'
            String ref = Integer.toString(i);
            if(!msq.hasMatches(ref))
                continue;

            boolean matches0 = false;
            boolean matches1 = false;
            List<List<Integer>> listOfMatches = msq.getMatchesOfSMARTS(ref);
            for (int im = 0; im < msq.getNumMatchesOfQuery(ref); im++)
            {
                List<Integer> sublist = listOfMatches.get(im);
                if (sublist.contains(atm0))
                    matches0 = true;
                if (sublist.contains(atm1))
                    matches1 = true;
            }

            if (matches0)
            {
                if (matches1) // V,V => give up!
                {
                    System.out.println("Case: "+matches0+" "+matches1);
                    break;
                } else { //V,F => ok, found it!
                    System.out.println("Case: "+matches0+" "+matches1);
                    res = atm0;
                    break;
                }
            } else {
                if (matches1) // F,V => ok, found it
                {
                    System.out.println("Case: "+matches0+" "+matches1);
                    res = atm1;
                    break;
                } else { //F,F => keep searching
                    System.out.println("Case: "+matches0+" "+matches1);
                    continue;
                }
            }
        }

        if (repOnScreen >= 3)
            System.out.println("Simplifies SMARTS search returns: "+res);

        return res;
    }

//-----------------------------------------------------------------------------

    /**
     * Analyze the matrix of booleans resulting from the subClass identification
     * and detects how many entries are true. Assuming that the general case in
     * which we have 1 true entry per each key of the map was detected before,
     * all cases having 0, 2 or 4 true entries in the whole matrix  are critical.
     * @param matrix coordinates used to identify the subCLASS
     * @return <code>true</code> if there are 0,2 or 4 true entries
     */
/* OLD to be deleted
    private boolean criticalCaseForSubClassDetection(Map<Integer,List<Boolean>> matrix)
    {
        int numTrue = 0;
        for (Integer i : matrix.keySet())
        {
            for (Boolean entry : matrix.get(i))
            {
                if (entry)
                    numTrue++;
            }
        }

        boolean res = false;
        if ((numTrue == 0) || (numTrue == 2) || (numTrue == 4))
            res = true;

        return res;
    }
*/
//-----------------------------------------------------------------------------

    /**
     * Analyze the matrix of booleans resulting from the subClass identification
     * counting all true entries.
     * @param matrix coordinates used to identify the subCLASS
     * @return the number of true entries 
     */
    private int howManyTrueEntries(Map<Integer,List<Boolean>> matrix)
    {
        int numTrue = 0;
        for (Integer i : matrix.keySet())
        {
            for (Boolean entry : matrix.get(i))
            {
                if (entry)
                    numTrue++;
            }
        }

        return numTrue;
    }

//-----------------------------------------------------------------------------

    /**
     * Analyze a list of booleans counting all true entries.
     * @param list of booleans to be evaluated
     * @return the number of true entries 
     */
    private int howManyTrueEntries(List<Boolean> list)
    {
        int numTrue = 0;
        for (Boolean entry : list)
        {
            if (entry)
                numTrue++;
        }

        return numTrue;
    }

//-----------------------------------------------------------------------------

    /**
     * Identify which of the SMARTS queries returns only false
     */
/*
//TODO maybe this is not necessary but leave it here for future developments
    private int whichIsCriticalSubClassDetection(Map<Integer,List<Boolean>> matrix)
    {
        int res = -1;
        ArrayList<Integer> numFalse = new ArrayList<Integer>();
        for (Integer i : matrix.keySet())
        {
            numFalse.add(0);
            int ii = numFalse.size() - 1;
            for (Boolean entry : matrix.get(i))
            {
                if (!entry)
                {
                    int old = numFalse.get(ii);
                    numFalse.add(ii, old +1);
                }
            }
        }
System.out.println("numFalse: "+numFalse);

    }
*/

//-----------------------------------------------------------------------------

    /**
     * Explore the molecule following existing bonds.
     * @param mol the molecule
     * @param seed index of the atom from which the exploration starts
     * @param ragius radius of the sphere of neighbours (in number of bonds) 
     * that will be visited in a reqursive fashion
     * @param flagID index of CDK flag that is used to mark visited atoms
     */

    private void visitAtomsInRadius(IAtomContainer mol, int seed, int radius, int flagID)
    {
        mol.getAtom(seed).setFlag(flagID,true);
        if (radius > 0)
        {
            for (IAtom atm : mol.getConnectedAtomsList(mol.getAtom(seed)))
                if (!atm.getFlag(flagID))
                    visitAtomsInRadius(mol,mol.getAtomNumber(atm),radius-1,flagID);
        }
        
    }

//-----------------------------------------------------------------------------
    /**
     * Analysis of an atom looking for the preCLASS property
     * which derives from a defined cutting rule and the 
     * subclass flag. The subclass flag is a numerical flag 
     * identifying which one of the two SMARTS
     * of the cutting rule, is satisfied by the atom.
     * @param rule name of the cutting rule to be searched
     * @param atm atom to analyze
     * @return the subCLASS flag as integer
     */

    private int getSubClass(String rule, IAtom atm)
    {
        int subClass = -1;
        String preClass = "";
        try {
            preClass = atm.getProperty("preCLASS").toString();
            String[] st = preClass.split("_");

            for (String pc : st) 
            {
                //compare only class without subclass flag 
                //(which is an integer as last character)
                String ppc = pc.subSequence(0,pc.length()-2).toString();
                int sc = Integer.parseInt(pc.substring(pc.length()-1));
                if (ppc.equals(rule))
		{
		   subClass = sc;
                   break;
		}
            }
        } catch (Throwable t) {
            System.err.println("ERROR in getting SubCLASS! "+t);
            System.err.println("Error occurred in 'Fragmenter.java' while dealing with "+rule);
            System.exit(0);
        }

        return subClass;
   }

//-----------------------------------------------------------------------------
    /**
     * Identifies and isolates disconnected portions of 
     * a chemical system. If a portion of the system is 
     * completely disconnected to the rest of the 
     * chemical system, it becomes a new fragment. 
     * @param mol atom container to analyse looking from
     *  disconnected fragments
     * @return the Set containing all the fragments as 
     * isolated AtomContainers
     */

    private AtomContainerSet isolateFrags(IAtomContainer mol)
    {
        AtomContainerSet frags = new AtomContainerSet();
        int num = 0;
        for (IAtom a : mol.atoms()) 
        {
            num++;

            if (repOnScreen >= 3)
                System.out.println("Isolating fragments - Attempt "+num);

            if (!a.getFlag(0)) {
                try {
                    // Clone the molecule to generate a new fragment
                    IAtomContainer frag = (IAtomContainer) mol.clone();

                    // As an oil drop, we will identify the fragment by esploring
                    // the still connected neighbourhood of a firstly picked atom.
                    // This atom is the root of the current exploration
                    IAtom af = frag.getAtom(mol.getAtomNumber(a));
                    af.setFlag(10,true);

                    // Debuging
                    if (repOnScreen >= 3) 
                    {
                        System.out.println(" DEBUG - FRAG before exploring the nighbouthood of atoms");
                        fragDEBUG(frag);
                    }

                    // Set label by exploring the neighbour atoms 
                    exploreNeighbour(af,frag,10);

                    // Debuging After
                    if (repOnScreen >= 3) 
                    { 
                        System.out.println(" DEBUG - FRAG after exploring the nighbouthood of atoms");
                        fragDEBUG(frag);
                    }

                    // Check for possible incompatibilities
                    if (mol.getAtomCount() != frag.getAtomCount()) 
                    {
                        System.out.println("ERROR! Unexpected mismatch between MOL and FRAG. Please send this error to the author.");
                        System.exit(0);
                    }

                    // Transfer info on fragment membership to mol
                    for (int i = 0; i < mol.getAtomCount(); i++) 
                    {
                        if (!mol.getAtom(i).getFlag(0))
                            mol.getAtom(i).setFlag(0,frag.getAtom(i).getFlag(10));
                    }

/*
System.out.println("BEFORE \nbonds: "+frag.getBondCount()+
                          "\natoms: "+frag.getAtomCount()+
                          "\nLP : "+frag.getLonePairCount()+
                          "\nEC : "+frag.getElectronContainerCount()+
                          "\nSEC: "+frag.getSingleElectronCount());
*/
                    // Remove all atoms and related objecte that do NOT belong to the fragment
                    boolean goon = true;
                    while (goon) 
                    {
                        for (IAtom delAtm : frag.atoms()) 
                        {
                            if (!delAtm.getFlag(10)) 
                            {
                                //bonds
                                for (IBond delBond : frag.getConnectedBondsList(delAtm))
                                    frag.removeBond(delBond);
                                //electroncontainers
                                for (IElectronContainer  delEC : frag.getConnectedElectronContainersList(delAtm))
                                    frag.removeElectronContainer(delEC);
                                //Lone Pairs
                                for (ILonePair delLP : frag.getConnectedLonePairsList(delAtm))
                                    frag.removeLonePair(delLP);
                                //SingleElectronContainers
                                for (ISingleElectron delSE : frag.getConnectedSingleElectronsList(delAtm))
                                    frag.removeSingleElectron(delSE);
                                frag.removeAtom(delAtm);
                            }
                        }
                        if (!MolecularUtils.containsFalseFlag(frag,10))
                            goon = false;
                    }
/*
System.out.println("AFTER \nbonds: "+frag.getBondCount()+
                          "\natoms: "+frag.getAtomCount()+
                          "\nLP : "+frag.getLonePairCount()+
                          "\nEC : "+frag.getElectronContainerCount()+
                          "\nSEC: "+frag.getSingleElectronCount()); 
*/

                    // Store fragment
                    frags.addAtomContainer(frag);

                    if (repOnScreen >= 3) 
                    {
                        System.out.println("\n NEW FRAGMET counts "+frag.getAtomCount()+" atoms and Bonds: "+frag.getBondCount());
                        fragDEBUG(frag);
                    }
                } catch (Throwable t) {
                    System.err.println("CLONING molecule failed! "+t);
                    System.exit(0);
                }
            }
        }
        if (repOnScreen >= 3)
            System.out.println("Got "+frags.getAtomContainerCount()+" fragments!");
        return frags;
    }

//---------------------------------------------------------------------------------------
    /**
     * Identifies atoms involved in the same n-hapto ligand of 
     * the seed atom
     * @param seed atom index of the starting point
     * @param candidates list of atoms that may belong to the 
     * multihapto ligand
     * @param atmsInHapto list of atoms already recognized to be
     * members of the same multihapto ligand
     * @param mol the molecule
     * @param flags
     * @return an updated list of members of the same multihapto 
     * ligand of the seed atom
     */

    private Set<Integer> exploreHapticity(int seed, ArrayList<Integer> candidates, Set<Integer> atmsInHapto, IAtomContainer mol, ArrayList<Boolean> flags)
    {
        for (IAtom connectedAtm : mol.getConnectedAtomsList(mol.getAtom(seed)))
        {
            int indxConnectedAtm = mol.getAtomNumber(connectedAtm);

            if (flags.get(indxConnectedAtm))
                continue;
            else 
                flags.set(indxConnectedAtm,true);

            if (candidates.contains(indxConnectedAtm))
            {
                atmsInHapto.add(indxConnectedAtm);
                atmsInHapto = exploreHapticity(indxConnectedAtm,candidates,atmsInHapto,mol,flags);
            }
        }
        return atmsInHapto;
    }

//---------------------------------------------------------------------------------------

    /**
     * Identifies all the atoms that can be reached by moving only 
     * along the existing connections (any bond in the connectivity 
     * matrix) starting from atom <code>root</code>.
     * For each reachable atom sets the flag number <code>flag</code>
     * of the container <code>frag</code> to <code>true</code>.
     * @param root atom used to start the exploration
     * @param frag the molecular container
     * @param flag identifier of the flag to be used to keep track 
     * of the exploration
     */

    private void exploreNeighbour(IAtom root, IAtomContainer frag, int flag)
    {
        if (repOnScreen >= 3)
            System.out.println("EXPLORATION on atom: "+frag.getAtomNumber(root)+" is "+root.getSymbol());

        for (IAtom connectedAtom : frag.getConnectedAtomsList(root))
        {
            String recFlag = "";
            for (int ri = 0; ri < recNum; ri++)
                recFlag = recFlag+"-";

            if (repOnScreen >= 3)
                System.out.println(recFlag+"> connected atom: "+frag.getAtomNumber(connectedAtom)+" is "+connectedAtom.getSymbol()+" Branches: "+frag.getConnectedAtomsCount(connectedAtom)+" Flag: "+connectedAtom.getFlag(flag));

            // in case of rings, avoid to bite your tail!
            if (connectedAtom.getFlag(flag))
                continue;
            connectedAtom.setFlag(flag,true);

            // move to the next shell of atoms
            if (frag.getConnectedAtomsCount(connectedAtom) > 1)
            {
                recNum++;
                if (repOnScreen >= 3)
                    System.out.println(recFlag+"> resursion on atom "+frag.getAtomNumber(connectedAtom)+" which is "+connectedAtom.getSymbol());
                exploreNeighbour(connectedAtom,frag,flag);
                recNum--;
            }
        }
    }

//-----------------------------------------------------------------------------
/*
//TODO delete: moved to other file
    /**
     * Check if a molecular fragment should NOT be rejected 
     * according to the list of recjection criteria
     * @param frag fragment to check
     * @return <code>true</code> if all criteria are satisfied
     */
/*
    private boolean keepFragment(GM3DFragment frag)
    {
        //Get all attachment points
        ArrayList<GM3DAttachmentPoint> allAPs = frag.getAllAPs();

        //get rid of fragments containing radicals
//TODO: make it optional when new CDK release is available
        if (frag.getSingleElectronCount() != 0)
        {
            if (repOnScreen >= 1)
                System.out.println("Fragment is a radical => rejected!");
            return false;
        }
//TODO: RADICALS are not written to output file due to lack in CDK release
// The current release is not able to handle radicals so we get rid of them here.
// Later, with the new release of CDK (1.5.x), this check must be moven in 'keepFragment'
//endTODO

        //loop over CLASS-based rejection criteria to identify unwanted CLASS
        for (String rejCls : rejClasses)
        {
            //loop over attachment points
            for (int i = 0; i < allAPs.size(); i++)
            {
                GM3DAttachmentPoint ap = allAPs.get(i);
                String apClass = ap.getAPClass();
        
                if (apClass.equals(rejCls))
                {
                    if (repOnScreen >= 1)
                        System.out.println("Fragment containing AP of class "+rejCls+" => rejected!");
                    return false;
                }
            }
        }

	//loop over rules for rejecting combinations of classes
	for (List l : rejClassCombination)
	{
	    String apcA = (String) l.get(0);
	    String apcB = (String) l.get(1);
            //loop over attachment points
            for (int i = 0; i < allAPs.size(); i++)
            {
                GM3DAttachmentPoint ap = allAPs.get(i);
                String apClass = ap.getAPClass();

                if (apClass.startsWith(apcA))
                {
		    for (int ib = 0; ib < allAPs.size(); ib++)
            	    {
			if (i == ib)
			{
			    continue;
			} else {
        	            GM3DAttachmentPoint apB = allAPs.get(ib);
	                    String apClassB = apB.getAPClass();
			    if (apClassB.startsWith(apcB))
			    {
	                        if (repOnScreen >= 1)
            		            System.out.println("Fragment containing combination of AP of classes '"+apClass+"'-'"+apClassB+"' => rejected!");
	                        return false;
			    }
			}
		    }
                }
            }
	}

        //Ckech the size requirements of the fragment
        int totHeavyAtm = 0;
        for (IAtom atm : frag.atoms())
        {
            if ((atm.getSymbol() != "H") && (atm.getSymbol() != Parameters.duSymbol))
                totHeavyAtm++;
        }
//        if (frag.getAtomCount() > fragmentMaxSize)
//TODO: update documentation this counts only heavy atoms
        if (totHeavyAtm > fragmentMaxSize)
        {
            if (repOnScreen >= 1)
                System.out.println("Fragment exceeds the maximun number of atoms => rejected!");
            return false;
        } else if (frag.getAtomCount() < fragmentMinSize)
        {
            if (repOnScreen >= 1)
                System.out.println("Fragment is too small => rejected!");
            return false;
        }

        //reject dummy or R-group that where contained in the initial SDF input
        for (IAtom atm : frag.atoms())
        {
            String smb = atm.getSymbol();
            if (smb.equals("R")) 
            {
                if (repOnScreen >= 1)
                    System.out.println("Fragment contains R-group => rejected!");
                return false;
            } else if (smb.equals("*"))
            {
                if (repOnScreen >= 1)
                    System.out.println("Fragment contains R-group => rejected!");
                return false;
            }
        }

        //loop over Elements' Symbol to reject 
        for (String rejSymb : rejElements)
        {
            for (IAtom atm : frag.atoms())
            {
                String smb = atm.getSymbol();
                if (smb.equals(rejSymb))
                {
                    if (repOnScreen >= 1)
                        System.out.println("Fragment contains element '"+rejSymb+"' => rejected!");
                    return false;
                }
            }
        }

        //isotopes
        if (rejIsotopes)
        {
            for (IAtom atm : frag.atoms())
            {
                String symb = atm.getSymbol();
                int a = -1;
                if (symb != duSymbol)
                {
                    a = atm.getMassNumber();
                    try {
                        IsotopeFactory factory = IsotopeFactory.getInstance(DefaultChemObjectBuilder.getInstance());
                        IIsotope major = factory.getMajorIsotope(symb);
                        if (a != major.getMassNumber())
                        {
                            if (repOnScreen >= 1)
                                System.out.println("Fragment contains isotope "+symb+a+" => rejected!");
                            return false;
                        }
                    } catch (Throwable t) {
                        if (repOnScreen >= 1)
                            System.out.println("\nWARNING! Not able to handle IsotopeFactory.");
                            System.out.println("ISOTOPE-based rejection criteria will be ignored!\n");
                        return false;
                    }
                }
            }
        }

        //Incomplete fragmentation: when an atom (or Du) has the same coords of an AP.
        //NB: min 2 bonds must be brocken to produce 2+ frags from a cyclic system
        for (int i=0; i<frag.getNumberOfAttachmentPoints(); i++)
        {
            GM3DAttachmentPoint ap = allAPs.get(i);
            ArrayList<Double> vector = ap.getAPVector();
            Point3d ap3d = new Point3d(vector.get(0),vector.get(1),vector.get(2));
            for (IAtom atm : frag.atoms())
            {
                Point3d atm3d = MolecularUtils.getCoords3d(atm);
                double dist = ap3d.distance(atm3d);
                if (dist < 0.0002)
                {
                    if (repOnScreen >= 1)
                       System.out.println("Attachment point "+i+" and atom "
                                +atm.getSymbol()+frag.getAtomNumber(atm)
                                +" coincide => rejected!");
                    return false;
                }
            }
        }

//TODO
//TODO //Add here other rejection critera
//TODO

        //loop over SMARTS-based rejection criteria
        for (String smrt : rejSMARTS)
        {
            try {
                SMARTSQueryTool query = new SMARTSQueryTool(smrt);
                if (query.matches(frag.toIAtomContainer()))
                {
                    if (repOnScreen >= 1)
                        System.out.println("Fragment matches SMARTS-based rejection criteria \""+smrt+"\" => rejected!");
                    return false;
                }
            } catch (Throwable t0) {
                if (repOnScreen >= 1)
                {
                    System.out.println("\nWARNING! Unable to check fragment's SMARTS.");
                    System.out.println("SMARTS-based rejection criteria will be ignored!\n");
                }
            }
        }

        //If not yet rejected, accept this fragment
        return true;
    }
*/
//-----------------------------------------------------------------------------

    /**
     * Compare fragment <code>frag</code> to the list of fragments 
     * in the library <code>libFile</code>
     * @param frag candidate new <code>GM3DFragment</code>
     * @param libFile name of SDF file containing the existing library 
     * <code>GM3DFragment</code>
     * @return <code>true</code> if <code>frag</code> is not found in
     * <code>list</code>, so it's a new molecular entity
     */

    private boolean newFragment(GM3DFragment frag, String libFile)
    {
        try {
            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(libFile), DefaultChemObjectBuilder.getInstance());
//            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(inFile),NoNotificationChemObjectBuilder.getInstance()); // returns ERROR! AtomContainer does'n looks like a fragment in 'DENOPTIM' format!
            while (reader.hasNext())
            {
                GM3DFragment oldFrag = new GM3DFragment(reader.next(),inFormat);
                if (oldFrag.sameFragOf(frag))
		{
                    return false;
		}
            }
        } catch (FileNotFoundException fnf) {
            if (repOnScreen >= 1)
                System.out.println("NO previous library found - Create NEW file: "+libFile);
            return true;
        } catch (Throwable t) {
	    System.err.println("\nPROBLEM in reading library "+libFile+" with Iterator. "+t);
            t.printStackTrace();
            System.out.println("\nERROR: cannot iterate through MDL file. This "
                                + "might be becouse the MDL format is prior to "
                                + "V2000 or because there are bond with "
                                + "type='any' (type=8 in connectivity matrix). "                                + "Please, try to provide a V2000 format and "
                                + "convert avoid bond type=8.");
            System.exit(0);
        }

        return true;
    }

//-----------------------------------------------------------------------------

    /**
     * Compare a <code>GM3DFragment</code> with the fragments the user wants
     * to ignore, which are provided as SDF library.
     * @param frag candidate <code>GM3DFragment</code> 
     * @param libIngorable name of SDF file containing all the fragment that 
     * can be ignored
     * @param libFormat format of fragments in <code>libIgnorable</code>
     * @return <code>true</code> if <code>frag</code> matches a fragment in
     * <code>libIgnorable</code>
     */

    private boolean hitIgnorableFragment(GM3DFragment frag, String libIgnorable, String libFormat) 
    {
        try {
            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(libIgnorable), DefaultChemObjectBuilder.getInstance());
            while (reader.hasNext())
            {
                GM3DFragment oldFrag = new GM3DFragment(reader.next(),libFormat);
                if (oldFrag.sameFragOf(frag))
                {
                    String hit =  MolecularUtils.getNameOrID(oldFrag);
                    frag.setProperty("IGNORABLE",hit);
                    return true;
                }
            }
        } catch (FileNotFoundException fnf) {
            System.out.println("\nLibrary of ignorable fragsmets '"+libIgnorable+"' NOT FOUND!");
            fnf.printStackTrace();
            System.exit(0);
        } catch (Throwable t) {
            System.err.println("\nPROBLEM in reading library "+libIgnorable+" with Iterator. "+t);
            t.printStackTrace();
            System.out.println("\nERROR: cannot iterate through MDL file. This "
                                + "might be becouse the MDL format is prior to "
                                + "V2000 or because there are bond with "
                                + "type='any' (type=8 in connectivity matrix). "                                + "Please, try to provide a V2000 format and "
                                + "convert avoid bond type=8.");
            System.exit(0);
        }

        return false;
    }

//-----------------------------------------------------------------------------

    /**
     * Compare a <code>GM3DFragment</code> to the target fragments
     * in the SDF library <code>libTargets</code> in which fragments are
     * reported with the format <code>format</code>.  
     * A dedicated property will be modified/created in <code>frag</code>
     * to record the ID of the target fragment matched.
     * @param frag candidate <code>GM3DFragment</code>
     * @param libTargets name of SDF file containing the target fragments
     * @param libFormat format of fragments in <code>libFile</code>
     * @return <code>true</code> if <code>frag</code> matches a fragment in
     * <code>libFile</code>
     */

    private boolean hitTargetFragment(GM3DFragment frag, String libTargets, String libFormat)
    {
        try {
            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(libTargets), DefaultChemObjectBuilder.getInstance());
            while (reader.hasNext())
            {
                GM3DFragment oldFrag = new GM3DFragment(reader.next(),libFormat);
                if (oldFrag.sameFragOf(frag))
                {
		    String hit =  MolecularUtils.getNameOrID(oldFrag);
		    frag.setProperty("TARGETHIT",hit);
                    return true;
                }
            }
        } catch (FileNotFoundException fnf) {
            System.out.println("\nLibrary of target fragsmets '"+libTargets+"' NOT FOUND!");
            fnf.printStackTrace();
            System.exit(0);
        } catch (Throwable t) {
            System.err.println("\nPROBLEM in reading library "+libTargets+" with Iterator. "+t);
            t.printStackTrace();
            System.out.println("\nERROR: cannot iterate through MDL file. This "
                                + "might be becouse the MDL format is prior to "
                                + "V2000 or because there are bond with "
                                + "type='any' (type=8 in connectivity matrix). "                                + "Please, try to provide a V2000 format and "
                                + "convert avoid bond type=8.");
            System.exit(0);
        }

        return false;
    }

//-----------------------------------------------------------------------------

    /**
     * Store molecules that will not be treated becouse of troubles
     * during Fragmenter execution. 
     * @param mol input molecular/chemical object
     * @param reason brief explamention of the reason for the rejection
     */
    private void rejectMol(IAtomContainer mol, String reason)
    {
        someMolRejected = true;
        numRejected++;
        mol.setProperty("REJECTED",reason);
        IOtools.writeSDFAppend(checkfile,mol,true);
        if (repOnScreen >= 1)
            System.out.println("Molecule Rejected! "+reason);
    }

//-----------------------------------------------------------------------------

    /**
     * Print the number of processed molecules and the total number of 
     * fragments stored at the moment when this method is called.
     * @param molnum total number of processed objects
     * @param fragNum number of fragment produced
     * @param file name of the output file
     */

    private void reportMolFragRatio(int molNum, int fragNum)
    {
        String line = " "+molNum+" "+fragNum+" \n";
        IOtools.writeTXTAppend(MlFrRatioFile,line,true);
    }

//-----------------------------------------------------------------------------

    /**
     * Prints the boolean list of 4 coordinates used to identify the subCLASS
     *  membership of a couple of atoms
     * @param mol molecule
     * @param rule cutting rule name for which check the atoms
     * @param atmid atom number of the two atoms to be checked
     */

    private void printInfosDEBUG(Map<Integer,List<Boolean>> matrix)
    {
        System.out.println("Case detected "+matrix);
    }

//-----------------------------------------------------------------------------

    /**
     * Prints useful infos on a fragment. For debuggigng purposes
     * @param frag fragment to be analyzed
     */
     private void fragDEBUG(IAtomContainer frag)
     {
        for (int j = 0; j < frag.getAtomCount(); j++) 
        {          
             try {
                  System.out.println("FRAG-atom "+j+" is "+frag.getAtom(j).getSymbol()+
                                     " FLAG-0: "+frag.getAtom(j).getFlag(0)+
                                     " FLAG-10: "+frag.getAtom(j).getFlag(10)+
                                     " -> '"+frag.getAtom(j).getProperty("CLASS").toString()+"'");
             } catch (Throwable t) {
                 System.out.println("FRAG-atom "+j+" is "+frag.getAtom(j).getSymbol()+ 
                                    " FLAG-0: "+frag.getAtom(j).getFlag(0)+
                                    " FLAG-10: "+frag.getAtom(j).getFlag(10)+
                                    " NOCLASS");
             }
        }
    }
//-----------------------------------------------------------------------------
}
