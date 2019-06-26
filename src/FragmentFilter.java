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
 * Toolbox for filtering GM3DFragments
 * 
 * @author Marco Foscato (University of Bergen)
 */


public class FragmentFilter
{
    //Rejection rules 
    // -> fragments containing AP with these CLASS+SubCLASS
    private static Set<String> rejClasses = Parameters.rejClasses;
    // -> fragments containing any combination of APs with these CLASS or portion of CLASS' name
    // NB: only the combinations are detected. Fragments containing one AP of one of these these classes but not combined with nothet one are not rejected
    public static Set<List<String>> rejClassCombination = Parameters.rejClassCombination;
    // -> fragments matching these SMARTS
    private static Set<String> rejSMARTS = Parameters.rejSMARTS;
    // -> fragments NOT matching these SMARTS
    private static Set<String> rejNotSMARTS = Parameters.rejNotSMARTS;
    // -> fragments containing these elements
    private static Set<String> rejElements = Parameters.rejElements;
    // -> fragments containing more than the number of atoms specified for each element
    public static Set<List<String>> rejFormulaMax = Parameters.rejFormulaMax;
    // -> fragments containing less than the number of atoms specified for each elements
    public static Set<List<String>> rejFormulaMin = Parameters.rejFormulaMin;
    // -> fragments hanimg more atoms than this limit
    private static int fragmentMaxSize = Parameters.fragmentMaxSize;
    private static int fragmentMinSize = Parameters.fragmentMinSize;
    // -> fragments containing uncommon isotopes
    private static boolean rejIsotopes = Parameters.rejIsotopes;

    private static String duSymbol = Parameters.duSymbol;

    //Reporting flag
    private static int repOnScreen = Parameters.report;

//----------------------------------------------------------------------------------------
    /**
     * Empty constructor
     */
    public FragmentFilter()
    {
    }

//-----------------------------------------------------------------------------
    /**
     * Check if a molecular fragment should NOT be rejected 
     * according to the list of recjection criteria
     * @param frag fragment to check
     * @return <code>true</code> if all criteria are satisfied
     */

    public static boolean keepFragment(GM3DFragment frag)
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

        
                if (apClass.startsWith(rejCls))
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

//TODO create option to define which alternative to use
//                if (apClass.contains(apcA))
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
//TODO create option to define which alternative to use
//                            if (apClassB.contains(apcB))
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
	    String symb = atm.getSymbol();
            if ((!symb.equals("H")) && (!symb.equals(Parameters.duSymbol)))
                totHeavyAtm++;
        }
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

        //reject dummy (not those from GM3DFragment) or R-group
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

	//reject fragments by molecular formula
	List<String> formCriteria = new ArrayList<String>();
	for (List<String> formulaLim : rejFormulaMax)
	{
	    formCriteria.add(formulaLim.toString());
	    List<Boolean> results = new ArrayList<Boolean>();
	    for (String elLim : formulaLim)
	    {
		String[] parts = elLim.split("(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)");
		String limitedElSymb = parts[0];
		int maxAtms = Integer.parseInt(parts[1]);
		int actualNum = 0;
		boolean violation = false;
		for (IAtom atm : frag.atoms())
		{
		    if (atm.getSymbol().equals(limitedElSymb))
		    {
			actualNum++;
			if (actualNum > maxAtms)
			{
			    violation = true;
			    break;
			}
		    }
		}
		if (violation)
		{
		    results.add(true);
		} else {
		    results.add(false);
		}
	    }

	    //are all conditions specified in the formula violated?
	    boolean allViolated = true;
	    int limID = 0;
	    for (int il = 0; il<results.size(); il++)
	    {
	        if (!results.get(il))
		{
		    allViolated = false;
		    limID = il;
		    break;
		}
	    }
	    if (allViolated)
	    {
                if (repOnScreen >= 1)
                    System.out.println("Fragment violates the limitations on molecular formula: MAX "+formCriteria.get(limID)+"  => rejected!");
		return false;
	    }
	}

        //reject fragments by molecular formula (MIN)
        List<String> formCriteriaMin = new ArrayList<String>();
        for (List<String> formulaLim : rejFormulaMin)
        {
            formCriteriaMin.add(formulaLim.toString());
            List<Boolean> results = new ArrayList<Boolean>();
            for (String elLim : formulaLim)
            {
                String[] parts = elLim.split("(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)");
                String limitedElSymb = parts[0];
                int minAtms = Integer.parseInt(parts[1]);
                int actualNum = 0;
                for (IAtom atm : frag.atoms())
                {
                    if (atm.getSymbol().equals(limitedElSymb))
                    {
                        actualNum++;
                    }
                }
                if (actualNum < minAtms)
                {
                    results.add(true);
                } else {
                    results.add(false);
                }
            }
            //are all conditions specified in the formula violated?
            boolean allViolated = true;
            int limID = 0;
            for (int il = 0; il<results.size(); il++)
            {
                if (!results.get(il))
                {
                    allViolated = false;
                    limID = il;
                    break;
                }
            }
            if (allViolated)
            {
                if (repOnScreen >= 1)
                    System.out.println("Fragment violates the limitations on molecular formula: MIN "+formCriteriaMin.get(limID)+"  => rejected!");
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
                if (!symb.equals(duSymbol))
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
        //NB: min 2 bonds must be broken to produce 2+ frags from a cyclic system
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
        // Check if 2 AP point to the same atom 
        for (int i=0; i<frag.getNumberOfAttachmentPoints(); i++)
        {
            GM3DAttachmentPoint apA = allAPs.get(i);
	    ArrayList<Double> vecHeadA = apA.getAPVector();
	    Point3d apA3d = new Point3d(vecHeadA.get(0),vecHeadA.get(1),vecHeadA.get(2));

	    for (int j=i+1; j<frag.getNumberOfAttachmentPoints(); j++)
            {
		GM3DAttachmentPoint apB = allAPs.get(j);
                ArrayList<Double> vecHeadB = apB.getAPVector();
		Point3d apB3d = new Point3d(vecHeadB.get(0),vecHeadB.get(1),vecHeadB.get(2));
		double dist = apB3d.distance(apA3d);
		if (dist < 0.0002)
		{
                    if (repOnScreen >= 1)
                       System.out.println("AP "+i+" and AP "+j
                                +" point to the same atom. Cyclic system fragmented by chopping two bonds of the same atom. This makes the building process more difficult => rejected!");
		    return false;
		}
	    }
	}

//TODO
//TODO //Add here other rejection critera
//TODO
/*
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
*/
        //SMARTS-based rejection criteria
        int rsi = 0;
        Map<String,String> rejSMARTSAllInOne = new HashMap<String,String>();
        for (String smrt : rejSMARTS)
        {
            rsi++;
            rejSMARTSAllInOne.put("rej"+rsi,smrt);
        }
        if (rsi>0)
        {
            ManySMARTSQuery msq = new ManySMARTSQuery(frag.toIAtomContainer(),rejSMARTSAllInOne);
            if (msq.hasProblems())
            {
                String cause = msq.getMessage();
                if (repOnScreen >= 1)
                {
                    System.out.println("\nWARNING! Problems in evaluating SMARTS-based rejection criteria. "+cause);
                }
            }
	    
	    Map<String,Integer> mapMatches = msq.getNumMatchesMap();
	    for (String s : mapMatches.keySet())
	    {
		if (repOnScreen >= 1)
                    System.out.println("Fragment matches SMARTS-based rejection criteria \""+rejSMARTSAllInOne.get(s)+"\" => rejected!");
		return false;
	    }
	}

        //SMARTS-based rejection criteria (NEGATION)
        int rnsi = 0;
        Map<String,String> rejNotSMARTSAllInOne = new HashMap<String,String>();
        for (String smrt : rejNotSMARTS)
        {
            rnsi++;
            rejNotSMARTSAllInOne.put("rej"+rnsi,smrt);
        }
        if (rnsi>0)
        {
            ManySMARTSQuery msq = new ManySMARTSQuery(frag.toIAtomContainer(),rejNotSMARTSAllInOne);
            if (msq.hasProblems())
            {
                String cause = msq.getMessage();
                if (repOnScreen >= 1)
                {
                    System.out.println("\nWARNING! Problems in evaluating SMARTS-based rejection criteria (with NOT operator). "+cause);
                }
            }

	    Map<String,Integer> mapMatches = msq.getNumMatchesMap();
System.out.println("matches: "+mapMatches);

	    for (String s : rejNotSMARTSAllInOne.keySet())
            {
	        if (!mapMatches.keySet().contains(s))
		{
                    if (repOnScreen >= 1)
                        System.out.println("Fragment matches NOT(SMARTS) rejection criteria \""+rejNotSMARTSAllInOne.get(s)+"\" => rejected!");
                    return false;
		}
            }
        }


        //If not yet rejected, accept this fragment
        return true;
    }

//-----------------------------------------------------------------------------
}
