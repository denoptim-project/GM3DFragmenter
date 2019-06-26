import java.io.*;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import javax.vecmath.Point3d;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;

import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Atom;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.tools.FormatStringBuffer;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.geometry.alignment.KabschAlignment;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.tools.FormatStringBuffer;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.smsd.Isomorphism;
import org.openscience.cdk.smsd.interfaces.Algorithm;
import org.openscience.cdk.io.XYZWriter;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.config.IsotopeFactory;




public class ExtractMostCommonConformer
{

    //Verbosity level
    private static int verbosity = 1;

    //Name of input file
    private static String inputFileName;

    //List of options 
    private static Map<String,String> options = new HashMap<String,String>();
    private static Map<String,String> optionsHelp = new HashMap<String,String>();

    //Map of isomers
    private static Map<String,String> isomersFileName = new HashMap<String,String>();
    private static Map<String,Integer> isomersCount = new HashMap<String,Integer>();

//------------------------------------------------------------------------------

    public static void main(String[] args)
    {
        //initialize options and environment
        if (verbosity > 0)
            printInit();

        setDefault();

        //Update options and environment
        readOpts(args);
        
        ExtractMostCommonConformer grmsd = new ExtractMostCommonConformer();
        
        //Define file names
        File f = new File(inputFileName);
        String onlyName = f.getName();
        String[] partsOfName = onlyName.split("[.]");
        String rootName = "";
        if (partsOfName.length > 1)
        {
            for (int i=0; i<(partsOfName.length - 1); i++)
            {
                rootName = rootName + partsOfName[i];
            }
        } else {
            rootName = onlyName;
        }
        String outFile = rootName + "_mostRepMol.sdf";

        boolean useRMSDIntDist = false;
        if (options.get("-useRMSdevIntDist").equals("NO"))
            useRMSDIntDist = false;
        else if (options.get("-useRMSdevIntDist").equals("YES"))
            useRMSDIntDist = true;
        else
            terminator("Ambiguous value of '-useRMSdevIntDist'. Use only 'YES' or 'NO'");

//String tmpFileVedi = "combinationsAligned.sdf";
//boolean noCenterOfMass = true;
        boolean bondSensitive = true;
        boolean ringmatch = false;
        boolean stereoMatch = false;
        boolean fragmentMinimization = false;
        boolean energyMinimization = false;

        //Split the library in files containing only the first N rotamers
        int maxStorage = 100; // => this is N
        getGroupsOfRotamers(inputFileName,maxStorage);

        //Loop over each group of rotamers
	int numGroupsDone = 0;
        for (String isomerID : isomersFileName.keySet())
        {
	    numGroupsDone++;
	    String localInputFile = isomersFileName.get(isomerID);
            System.out.println( numGroupsDone + "/" + isomersFileName.size() 
				+ " Working on " + localInputFile
                                + " => " + isomersCount.get(isomerID));

            //Get molecules
            ArrayList<IAtomContainer> mols = grmsd.readSDFFile(localInputFile);
       
            //Initialize matrix of RMSD 
            int size = mols.size();
            DistanceMatrix dm = new DistanceMatrix(size);

            //Calculate and store RMSD
            try
            {
                for (int i=0; i<mols.size()-1; i++)            
                {
                    //Get the reference
                    IAtomContainer refmol = mols.get(i);
		    MolecularUtils.ensure3d(refmol);

                    //Make it a fragment to decode AP-related information
                    GM3DFragment refFrag = new GM3DFragment(refmol,"DENOPTIM");

                    IAtomContainer refApAsAtms = refFrag.getContainerWithAPasAtoms();

                    //Collect the AP classes
                    ArrayList<GM3DAttachmentPoint> allAPs = refFrag.getAllAPs();
                    Set<String> allClasses = new HashSet<String>();
                    for (GM3DAttachmentPoint ap : allAPs)
                    {
                        allClasses.add(ap.getAPClass());
                    }
                    Map<String,String> classToElement = refFrag.getClassToElement(allClasses,refApAsAtms);
                    //Change pseudo atoms 'AP' and 'Du' to real elements to allow alignment
                    for (IAtom atm : refApAsAtms.atoms())
                    {
                        if (atm.getSymbol().equals("AP"))
                        {
                            GM3DAttachmentPoint ap = (GM3DAttachmentPoint) atm.getProperty("CLASS");
                            String apClass = ap.getAPClass();
                            atm.setSymbol(classToElement.get(apClass));
                        } else if (atm.getSymbol().equals("Du")) {
                            atm.setSymbol("Tc");
                        }
                    }

                    IAtomContainer genRef = (IAtomContainer) refApAsAtms.clone();

//IAtomContainer genRef = (IAtomContainer) refmol.clone();

                    //Loop on molecules to calculate distance with the reference
                    for (int j=i+1; j<mols.size(); j++)
                    {
//System.out.println("Combination "+i+" - "+j);
                        IAtomContainer query = mols.get(j);
			MolecularUtils.ensure3d(query);

                        GM3DFragment queryFrag = new GM3DFragment(query,"DENOPTIM");
                        IAtomContainer queryApAsAtms = queryFrag.getContainerWithAPasAtoms();
/*
queryApAsAtms.setProperty("cdk:Title","AP");
IOtools.writeSDFAppend("query_on_the_way.sdf",queryApAsAtms,true);
*/
                        //Assuming that the two fragments are two copies of the same GM3DFragment,
                        // use the class-to-element crrespondence defined for the reference

                        //Change pseudo atoms 'AP' and 'Du' to real elements to allow alignment
                        for (IAtom atm : queryApAsAtms.atoms())
                        {
                            if (atm.getSymbol().equals("AP"))
                            {    
                                GM3DAttachmentPoint ap = (GM3DAttachmentPoint) atm.getProperty("CLASS");
                                String apClass = ap.getAPClass();
                                atm.setSymbol(classToElement.get(apClass));
                            } else if (atm.getSymbol().equals("Du")) {
                                atm.setSymbol("Tc");
                            }

                        }
/*
queryApAsAtms.setProperty("cdk:Title","AP-to-El");
IOtools.writeSDFAppend("query_on_the_way.sdf",queryApAsAtms,true);
*/
                        IAtomContainer genQuery = (IAtomContainer) queryApAsAtms.clone();
/*
genQuery.setProperty("cdk:Title","GENUINE");
IOtools.writeSDFAppend("query_on_the_way.sdf",genQuery,true);
*/
//IAtomContainer genQuery = (IAtomContainer) query.clone();


//System.out.println("Molecule J: "+queryApAsAtms);
//IOtools.pause();

                        Isomorphism comparison = new Isomorphism(Algorithm.DEFAULT, true);
                        // set molecules, 'remove hydrogens=false', 'clean and configure molecule'=false
                        comparison.init(genRef, genQuery, false, false);
                        comparison.setChemFilters(stereoMatch, fragmentMinimization, energyMinimization);

                        double minValue = Double.MAX_VALUE;
                        int minValueMapping = -1;
                        int mapIdx = 0;

                        //Get RMSD alignment and RMS deviantion of interatomic distances
                        if (comparison.getFirstAtomMapping() != null)
                        {
                            IAtomContainer Mol1 = comparison.getReactantMolecule();
                            IAtomContainer Mol2 = comparison.getProductMolecule();

                            for (Map<Integer,Integer> mapping : comparison.getAllMapping())
                            {
                                mapIdx++;
                                double rmsdAlign = -1.0;
                                double rmsdIntDist = -1.0;
  
                                if (verbosity > 2)
                                    System.out.println(" Mapping "+mapIdx+": ");
 
                                //Prepare vectors of atoms
                                int sz = mapping.entrySet().size();
				if (sz != Mol1.getAtomCount())
				{
				    if (verbosity > 2)
                                        System.out.println(" skipping uncomplete mapping");

				    continue;
				}
//System.out.println(" Tot. atoms: " + Mol1.getAtomCount() + " - " + Mol2.getAtomCount());
//System.out.println(" Size mapping: " + sz);
                                Atom[] a1 = null;
                                Atom[] a2 = null;
                                a1 = new Atom[sz];
                                a2 = new Atom[sz];

                                //Fill vector of atoms
/*NOTE: weight doesnt seem to take effect
                                double[] weightsAlignment = new double[sz];
                            double maxWt = 0.0;
*/
                                int k = 0;
                                for (Map.Entry map : mapping.entrySet())
                                {
                                    int eIdx = (int) map.getKey();
                                    int pIdx = (int) map.getValue();
                                    IAtom eAtom = Mol1.getAtom(eIdx);
                                    IAtom pAtom = Mol2.getAtom(pIdx);

/*NOTE: weight doesnt seem to take effect
                                //Get weight for alignment
                                double wt = 1.0;
                                if (noCenterOfMass)
                                {
                                    double eMass = getAtomicMass(eAtom);
                                    double pMass = getAtomicMass(pAtom);
                                    wt = eMass + pMass;
                                    wt = wt / (double) 2.0;

                                    //Keep track of the max mass
                                    if (wt > maxWt)
                                        maxWt = wt;

                                    wt = (double) 1.0 / wt;
                                }
                                weightsAlignment[k] = wt;
*/

                                    if (verbosity > 2) 
                                    {
                                        String eAtmRef = MolecularUtils.getAtomRef(eAtom,Mol1);
                                        String pAtmRef = MolecularUtils.getAtomRef(pAtom,Mol2);
                                        System.out.println(" - "+eAtmRef+" -> "+pAtmRef);
                                    }

                                    //Finally put the atom in the vector
                                    a1[k] = (Atom) eAtom;
                                    a2[k] = (Atom) pAtom;
                                    k++;
                                }

/*
/NOTE: weight doesnt seem to take effect
                            //scale weight so to have 1 for the smallest
                            for (int iw=0; iw<sz; iw++) 
                            {
                                weightsAlignment[iw] = weightsAlignment[iw] * maxWt;
                                if (iw == 1)
                                    weightsAlignment[iw] = 10000000.0;
                            }
*/

/*
System.out.print("Weights: ");
for (int iii=0; iii<sz;iii++)
    System.out.print(weightsAlignment[iii]+" ");
System.out.println(" ");
*/

                                //Align molecules
                                KabschAlignment sa = null;
                                try
                                {
				    if (verbosity > 2)
                                    {
					System.out.println(" KabschAlignment");
				    }
//NOTE: weight doesnt seem to take effect
//                                sa = new KabschAlignment(a1, a2, weightsAlignment);
                                    sa = new KabschAlignment(a1, a2);
                                    sa.align();
                                }
                                catch (CDKException cdke)
                                {
                                    System.out.println("ERROR in using KabschAlignment!");
//TODO deal with the error
                                    System.exit(-1);
                                }

                                //Rototranslation of molecule 1
                                Point3d cm1 = sa.getCenterOfMass();
                                for (int ia = 0; ia < Mol1.getAtomCount(); ia++)
                                {
                                    Atom a = (Atom) Mol1.getAtom(ia);
                                    Point3d newPlace = new Point3d(a.getPoint3d().x - cm1.x, a.getPoint3d().y - cm1.y, a.getPoint3d().z - cm1.z);
                                    a.setPoint3d(newPlace);
                                }

                                //Rototranslation of molecule 2
                                sa.rotateAtomContainer(Mol2);
/*
Mol1.setProperty("cdk:Title","Combination_"+i+"-"+j+"_mol-"+i);
IOtools.writeSDFAppend(tmpFileVedi,Mol1,true);
Mol2.setProperty("cdk:Title","Combination_"+i+"-"+j+"_mol-"+j);
IOtools.writeSDFAppend(tmpFileVedi,Mol2,true);
*/

                                //Get RMSD for structure superposition
                                rmsdAlign = sa.getRMSD();

                                //Get RMS deviation of Itramolecular Distances
                                rmsdIntDist = getRMSDevIntramolecularDistances(a1,a2);

                                //Chose the value to be minimized over the mappings
                                double localValue = 0.0d;
                                if (useRMSDIntDist)
                                {
                                    localValue = rmsdIntDist;        
                                } else {
                                    localValue = rmsdAlign;
                                }

                                //Compare with previous mappings
                                if (localValue < minValue)
                                {
                                    minValue = localValue;
                                    minValueMapping = mapIdx;
                                }

                                //Report results
                                if (verbosity > 2)
                                {
                                    System.out.println("RSMD from the alignmen: "+rmsdAlign);
                                    System.out.println("RSM deviation of intramolecular distances: "+rmsdIntDist);
                                }

                            } //End loop over atom mapping

                        } else //if no mapping was foung betwee the molecules
                        {
                            if (verbosity > 0)
                            {
                                System.err.println("WARNING! No valid mapping of atoms between structures.");
                                System.err.println("Setting entry in distance matrix to "+minValue+".");
                            }
                        }

                        //Store value in distance matrix
                        dm.setElement(i,j,minValue);
                        dm.setElement(j,i,minValue);

                        if (verbosity > 1)
                        {
                            System.out.println("Setting element "+i+";"+j+" = "+minValue);
                            System.out.println("Setting element "+j+";"+i+" = "+minValue);
                        }

                    } //End loop over j-molecules

                } //End loop over i-molecules

            }
            catch (Exception ex)
            {
                grmsd.printExceptionChain(ex);
                System.exit(-1);
            }

            //Prune the distance matrix up to the closest 10% of the entries
            for (int i=0; i < (size-1); i++)
            {
                dm.pruneMostDistant();
//              System.out.println(dm.toString());
//              IOtools.pause();
            }

            //Identify the resulting molecule 
            int mostRepMolID = dm.getLabels().get(0);

            //Store the output
            IAtomContainer mostRepMol = mols.get(mostRepMolID);
	    mostRepMol.setProperty("CONFORMATIONS_SAMPLE_SIZE",size);
            IOtools.writeSDFAppend(outFile, mostRepMol, true);
            if (verbosity > 0)
                System.out.println("Most representative mol is: "+mostRepMolID);

	} //end of loop over groups of isomers

        //Set exit value and... ... bye, bye!
        if (verbosity > 0)
            printFinal();

        System.exit(0);
    }

//------------------------------------------------------------------------------

    private static double getRMSDevIntramolecularDistances(Atom[] a, Atom[] b)
    {
        double result = 0.0d;

        //Check if conditions for running this calculation are satisfied
        if (a.length != b.length)
        {
            System.err.println("ERROR! Attempt to calculate RMS deviation of "+
                                "structures with different number of atoms");
            System.exit(-1);
        }

        for (int i=0; i<a.length; i++)
        {
            //Get first point for both structures
            Point3d pAi = a[i].getPoint3d();
            Point3d pBi = b[i].getPoint3d();

            for (int j=i+1; j<a.length; j++)
            {
                //Get second point for both structures
                Point3d pAj = a[j].getPoint3d();
                Point3d pBj = b[j].getPoint3d();

                //Get distances
                double dAij = pAi.distance(pAj);
                double dBij = pBi.distance(pBj);
/*
System.out.println("Points A: "+pAi+" "+pAj);
System.out.println("Points B: "+pBi+" "+pBj);
System.out.println("Distanced: "+dAij+" "+dBij);
*/
                //deviation
                double dev = dAij - dBij;
                double devSquare = (dev) * (dev);

                //sum to others
                result = result + devSquare;
            }
        }

//System.out.println("summ: "+result);

        //Complete calucation with prefactor...
        int n = a.length;
        double denominator = n * (n - 1);
        double preFactor = 2.0 / denominator; 
//System.out.println("Prefactor: "+preFactor);
        result = preFactor * result;
        // ... and square root
        result = Math.sqrt(result);

        return result;
    }

//------------------------------------------------------------------------------

    /**
     * Execute an external process as a command
     * @param command
     * @param waitForResponse
     * @return value <code>0</code> if the command was executed successfully
     */

    private static int executeCommand(String command, boolean waitForResponse)
    {
        int exitVal = 0;

        ProcessBuilder pb = new ProcessBuilder("sh", "-c", command);
        pb.redirectErrorStream(true);

        try
        {
            final Process proc = pb.start();
            Runtime rt = Runtime.getRuntime();
            rt.addShutdownHook(new Thread()
            {
                @Override
                public void run()
                {
                    //System.out.println("Running Shutdown Hook");
                    proc.destroy();
                }
            });

            if (waitForResponse)
            {
                // To capture output from the shell
                InputStream shellIn = proc.getInputStream();

                // Wait for the shell to finish and get the return code
                exitVal = proc.waitFor();

                shellIn.close();
            }
        }
        catch (IOException | InterruptedException ioe)
        {
            String mess = ioe.getMessage();
            System.out.println("ERROR! Interrupted command: "+mess);
            IOtools.pause();
        }

        return  exitVal;
    }

//------------------------------------------------------------------------------

    private ArrayList<IAtomContainer> readSDFFile(String filename)
    {
        MDLV2000Reader mdlreader = null;
        ArrayList<IAtomContainer> lstContainers = new ArrayList<IAtomContainer>();

        try
        {
            mdlreader = new MDLV2000Reader(new FileReader(new File(filename)));
                      ChemFile chemFile = (ChemFile) mdlreader.read((ChemObject) new ChemFile());
            lstContainers.addAll(
                            ChemFileManipulator.getAllAtomContainers(chemFile));
        }
        catch (CDKException cdke)
        {
            System.err.println(cdke.getMessage());
            System.exit(-1);
        }
        catch (IOException ioe)
                {
            System.err.println(ioe.getMessage());
            System.exit(-1);
                }
                finally
                {
                   try
                        {
                                if(mdlreader != null)
                                        mdlreader.close();
                        }
                        catch (IOException ioe)
                        {
                            ioe.printStackTrace();
                            System.exit(-1);
                        }
                }

                if (lstContainers.isEmpty())
                {
                    System.err.println("No data found in " + filename);
            System.exit(-1);
                }

                return lstContainers;
        //return lstContainers.toArray(new IAtomContainer[lstContainers.size()]);
    }    
    
//------------------------------------------------------------------------------
 
    // this is the javaboutique version
    private void printExceptionChain(Throwable thr)
    {
        
        StackTraceElement elements[] = thr.getStackTrace();
        for (int i = 0, n = elements.length; i < n; i++) 
        {
            System.err.println(elements[i].getFileName() +  ":" 
                            + elements[i].getLineNumber() + ">> " 
                            + elements[i].getMethodName() + "()"
                            );
        }
    
    
        StackTraceElement[] steArr;
        Throwable cause = thr;
        while (cause != null) 
        {
            System.err.println("-------------------------------");
            steArr = cause.getStackTrace();
            StackTraceElement s0 = steArr[0];
            System.err.println(cause.toString());
            System.err.println("  at " + s0.toString());
            if (cause.getCause() == null) 
            {
                System.err.println("-------------------------------");
                cause.printStackTrace();
            }
            cause = cause.getCause();
        }
    }

//------------------------------------------------------------------------------

    private static double getAtomicMass(IAtom a) 
    {
        IsotopeFactory factory = null;
        try {
            factory = IsotopeFactory.getInstance(a.getBuilder());
        } catch (Exception e) {
            System.err.println("Error while instantiating the isotope factory: "+e.getMessage());
        }

        assert factory != null;
        double am = factory.getMajorIsotope( a.getSymbol() ).getExactMass();

        return am;
    }

    
//------------------------------------------------------------------------------

    private double alignMols(IAtomContainer ac1, IAtomContainer ac2, Atom[] a1, Atom[] a2)
    {
        KabschAlignmentMF sa = null;
 
        try 
        {
            sa = new KabschAlignmentMF(a1, a2);
            sa.align();
        } 
        catch (CDKException cdke)
        {
                System.out.println("ERROR in using KabschAlignment!");
                System.exit(-1);
        }


        Point3d cm1 = sa.getCenterOfMass();
        for (int i = 0; i < ac1.getAtomCount(); i++) 
        {
           Atom a = (Atom) ac1.getAtom(i);
           Point3d newPlace = new Point3d(a.getPoint3d().x - cm1.x, a.getPoint3d().y - cm1.y, a.getPoint3d().z - cm1.z);
           a.setPoint3d(newPlace);
        }

        sa.rotateAtomContainer(ac2);

System.out.println("SEE Aligned");
IOtools.writeSDFAppend("seeMeAlign.sdf",ac1,false);
IOtools.writeSDFAppend("seeMeAlign.sdf",ac2,true);
IOtools.pause();
        
        return sa.getRMSD();
    }

//------------------------------------------------------------------------------

    /**
     * Split the initial file in a list of files: one file per each isomer
     */
    private static void getGroupsOfRotamers(String filename, int maxStorage)
    {
        // loop over molecules
        int molnum = 0;
        try {
            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(filename), DefaultChemObjectBuilder.getInstance());
            while (reader.hasNext())
            {
                molnum++;
                IAtomContainer mol = reader.next();
                String isomerID = mol.getProperty("ISOMER").toString();

                if (isomersFileName.keySet().contains(isomerID))
                {
                    String isomerStorageFile = isomersFileName.get(isomerID);
                    int numIsomersRead = isomersCount.get(isomerID);
                    if (numIsomersRead < maxStorage)
                    {
                        IOtools.writeSDFAppend(isomerStorageFile,mol,true);
                        isomersCount.put(isomerID, numIsomersRead + 1);
                    }
                } else {
                    String isomerStorageFile = "first" + maxStorage + "_"
                                                + isomerID + ".sdf";
                    IOtools.writeSDFAppend(isomerStorageFile,mol,true);
                    isomersFileName.put(isomerID,isomerStorageFile);
                    isomersCount.put(isomerID,1);
                }
            }
            reader.close();
        } catch (Throwable t) {
            System.err.println("\nERROR in reading the file with IteratorMDLReader. "+t);
            t.printStackTrace();
            System.exit(0);
        }
    }

//------------------------------------------------------------------------------

    /**
     * Define all default variables and options
     */
    private static void setDefault()
    {
        //Default value for the options
        String          optKey = "-useRMSdevIntDist";
        String          optDef = "NO";
        String          optHlp = "Set to 'YES' to use RMS deviation of interatomic distances instead of RMSD of superposition";
        options.put(    optKey ,
                        optDef );
        optionsHelp.put(optKey ,
                        optHlp );

/*
        String          optKey = "-";
        String          optDef = "";
        String          optHlp = " [default: "+
                        optDef +"]";
        options.put(    optKey ,
                        optDef );
        optionsHelp.put(optKey ,
                        optHlp );
*/

    }

//------------------------------------------------------------------------------
    /**
     * Read options
     */
    private static void readOpts(String[] args)
    {
        if (args.length < 1)
        {
            printUsage();
            terminator("Check command line arguments");
        }

        boolean inpIsFound = false;
        for (String s : args)
        {
            boolean optIndetified = false;
            for (String key : options.keySet())
            {
                if (s.startsWith(key))
                {
                    String val = s.substring(key.length());
                    options.put(key,val);
                    optIndetified = true;
                }
            }
            if (!optIndetified)
            {
                if (inpIsFound)
                    terminator("Check option "+s);
                inpIsFound = true;
                inputFileName = s;
            }
        }
    }

//------------------------------------------------------------------------------
    /**
     * Terminate execution with error message
     */
    private static void  terminator(String msg)
    {
        System.out.println("\n Unexpected termination!"
                +" Cause or hint to solve the problem: \n"+msg+"\n");
        System.exit(0);
    }

//------------------------------------------------------------------------------
    /**
     * Write the initial message
     */
    private static void printInit()
    {
        System.out.println("\n\n*************************************"
                +"\n  Extract Representative Fragment"
                +"\n*************************************\n");
    }

//------------------------------------------------------------------------------
    /**
     * Write the final message
     */
    private static void printFinal()
    {
        System.out.println("\n Normal termination\n");
    }

//------------------------------------------------------------------------------
    /**
     * Write the manual/usage information
     */
    private static void printUsage()
    {
        System.out.println("\n Usage: \n java -jar ExtractMostCommonConformer.jar <options> <input.sdf>\n");
        System.out.println("\n\n Options:");
        for (String key : options.keySet())
            System.out.println(key+"\t=> \t"+optionsHelp.get(key));
        System.out.println("\n");
    }
    
//------------------------------------------------------------------------------            
    
}
