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

/**
 * Toolbox for molecular objects
 * 
 * @author Marco Foscato (University of Bergen)
 */

import java.util.Set;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.Bond;
import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import org.openscience.cdk.smiles.FixBondOrdersTool;


/**
 * Toolbox for molecular objects
 */
public class MolecularUtils
{

    //Reporting flag
    private static int repOnScreen = Parameters.report;

//------------------------------------------------------------------------------

/**
 * Counts the number of atoms for each element in the system and returns
 * a string with the molecular formula (elements ordered in alphabetic order)
 */

    public static String getMolecularFormulaAsStringForAtomContainer(IAtomContainer mol)
    {
        Map<String,Integer> formulaMap = new HashMap<String,Integer>();
        for (IAtom atm : mol.atoms())
        {
            String elSymbol = atm.getSymbol();
            if (formulaMap.keySet().contains(elSymbol))
            {
                int num = formulaMap.get(elSymbol) + 1;
                formulaMap.put(elSymbol,num);
            } else {
                formulaMap.put(elSymbol,1);
            }
        }

        String formula = "";

        Set<String> elSet = formulaMap.keySet();
        ArrayList<String> els = new ArrayList<String>();
        for (String el : elSet)
            els.add(el);

        Collections.sort(els);

        for (String el : els)
        {
            formula = formula + el + formulaMap.get(el);
        }

	return formula;
    }

//------------------------------------------------------------------------------

/**
 * Counts the number of atoms for each element in the system
 */

    private Map<String,Integer> getMolecularFormulaOfAtomContainer(IAtomContainer mol)
    {
	Map<String,Integer> formula = new HashMap<String,Integer>();
        for (IAtom atm : mol.atoms())
        {
            String elSymbol = atm.getSymbol();
            if (formula.keySet().contains(elSymbol))
            {
                int num = formula.get(elSymbol) + 1;
                formula.put(elSymbol,num);
            } else {
                formula.put(elSymbol,1);
            }
        }
	return formula;
    }

//------------------------------------------------------------------------------

    /**
     * Returns a string with the element symbol and the atom number (1-n)
     * of the given atom
     */
    public static String getAtomRef(IAtom atm, IAtomContainer mol)
    {
	String ref = atm.getSymbol() + (mol.getAtomNumber(atm) +1);
	return ref;
    }

//------------------------------------------------------------------------------

    /**
     * Tool for detecting the smallest ring involving two atoms
     * @param atmS first atom involved (seed of the search)
     * @param atmT second atom involved 
     * @param mol the molecular object
     * @return the size of the ring or -1 if no ring has been found
     */
    public static int getSmallestRing(IAtom atmS, IAtom atmT, IAtomContainer mol)
    {
        if (repOnScreen >= 3)
            System.out.println("Looking for ring involving "
                        +getAtomRef(atmS, mol)
                        +" and "
                        +getAtomRef(atmT, mol));
        int size=-1;

        boolean moreTrials = true;
        boolean first = true;
        int level = 0;

        Map<IAtom,Integer> visited = new HashMap<IAtom,Integer>();
        visited.put(atmS,level);
        while (moreTrials)
        {
            // Get atoms from previous level
            List<IAtom> inLevel = new ArrayList<IAtom>();
            for (IAtom a : visited.keySet())
                if (visited.get(a) == level)
                    inLevel.add(a);

            // Get atoms in next level
            level++;
            List<IAtom> nextLevel = new ArrayList<IAtom>();
            for (IAtom sfpl : inLevel)  //sfpl: seed from previous level
            {
                List<IAtom> neighbours = mol.getConnectedAtomsList(sfpl);
                //in case it's a dummy ingore and move on
                for (IAtom nbr : neighbours)
                {
                    if (nbr.getSymbol().equals(Parameters.duSymbol))
                    {
                        List<IAtom> neighboursOfDu = mol.getConnectedAtomsList(nbr);
                        visited.put(nbr,-1);
                        nextLevel.remove(nbr);
                        nextLevel.removeAll(neighboursOfDu); //remove possible duplicates
			for (IAtom n : neighboursOfDu)
			{
			    if (!Parameters.metals.contains(n.getSymbol()))
			    {
                                nextLevel.add(n); //add new candidate src
			    }
			}
                    }
                }
                nextLevel.removeAll(neighbours); //remove possible duplicates
		for (IAtom n : neighbours)
                {
                    if (!Parameters.metals.contains(n.getSymbol()))
                    {
                        nextLevel.add(n); //add new candidate src
                    }
                }
            }
            nextLevel.removeAll(visited.keySet());
            if (nextLevel.size() == 0 )
            {
                moreTrials = false;
            } else {
                for (IAtom nbr : nextLevel)
                {
                    if (nbr == atmT)
                    {
                        if (!first)
                        {
                            size = level + 1;
                            moreTrials = false;
                            if (repOnScreen >= 3)
                                System.out.println("\nFound target: "+getAtomRef(nbr,mol)+" - Ring size: "+size);
                        } else {
                            first = false;
                        }
                    } else {
                        if (repOnScreen >= 3)
                            System.out.print((mol.getAtomNumber(nbr)+1)+" ");
                        visited.put(nbr,level);
                    }
                }
                if (repOnScreen >= 3)
                    System.out.print("\n");
            }
        }

        return size;
    }

//------------------------------------------------------------------------------

    /**
     * Tool for detecting the smallest ring involving two given atoms and 
     * any metal
     * @param atmS first atom involved (seed of the search)
     * @param atmT second atom involved
     * @param mol the molecular object
     * @param maxSz maximum size
     * @return the size of the ring or -1 if no ring has been found within the
     * given size
     */
    public static int getSmallestOMRing(IAtom atmS, IAtom atmT, IAtomContainer mol, int maxSz)
    {
        if (repOnScreen >= 3)
            System.out.println("Looking for OM-ring involving "
                        +getAtomRef(atmS, mol)
                        +" and "
                        +getAtomRef(atmT, mol));

        int size=-1;

	// Find metals
	ArrayList<IAtom> metals = new ArrayList<IAtom>();
	for (IAtom a : mol.atoms())
	{
	    if (Parameters.metals.contains(a.getSymbol()))
	    {
		metals.add(a);
	    }
	}
	if (metals.size() < 1)
	{
	    return size;
	}

	// Find smallest om-rings
        AllRingsFinder arf = new AllRingsFinder();
        IRingSet allRings;
        try {
            arf.setTimeout(Parameters.maxTimeAllRingFinder);
            allRings = arf.findAllRings(mol,maxSz);
        } catch (CDKException e) {
            System.err.println("Unable to identify OM-ring");
            e.printStackTrace();
	    return Integer.MAX_VALUE;
        }
	IBond bnd = mol.getBond(atmS,atmT);
	IRingSet ringsOnBnd = allRings.getRings(bnd);
	for (IAtomContainer ring : ringsOnBnd.atomContainers())
	{
//TODO del
//System.out.println(" Ring size on bond "+atmS.getSymbol()+mol.getAtomNumber(atmS)+"-"+atmT.getSymbol()+mol.getAtomNumber(atmT)+": "+ ring.getAtomCount());
	    for (IAtom metal : metals)
	    {
	        if (ring.contains(metal) && ring.getAtomCount()>size)
	        {
		    int numDu = 0;
		    for (IAtom a : ring.atoms())
		    {
			if (a.getSymbol().equals(Parameters.duSymbol))
			    numDu++;
		    }
		    size = ring.getAtomCount() - numDu;  
		    break;
		}
	    }
	}
//TODO del
//System.out.println("Smallest: "+size);
	return size;	
    }

//------------------------------------------------------------------------------
    /**
     * Evaluates the need of changing the bond orders to match aromaticity
     */
    public static String missmatchingAromaticity(IAtomContainer mol)
    {
	String cause = "";
        for (IAtom atm : mol.atoms())
        {
            //Check of carbons with or without aromatic flags
	    if (atm.getSymbol().equals("C")) 
	    {
		if (atm.getFormalCharge() == 0)
		{
		    if (mol.getConnectedAtomsCount(atm) == 3)
		    {
			if (atm.getFlag(CDKConstants.ISAROMATIC))
			{
			    int n = numOfBondsWithBO(atm,mol,"DOUBLE");
		            if (n == 0)
		            {
				cause = "Aromatic atom "+getAtomRef(atm,mol)+" has 3 connected atoms but no double bonds";
/*
		                if (repOnScreen >= 1)
		                    System.out.println("Aromatic atom "+getAtomRef(atm,mol)+" has 3 connected atoms but no double bonds");
*/
		                return cause;
		            }
			} else {
			    for (IAtom nbr : mol.getConnectedAtomsList(atm))
			    {
				if (nbr.getSymbol().equals("C"))
				{
				    if (nbr.getFormalCharge() == 0)
				    {
					if (mol.getConnectedAtomsCount(nbr) == 3)
					{
					    int nNbr = numOfBondsWithBO(nbr,mol,"SINGLE");
					    int nAtm = numOfBondsWithBO(atm,mol,"SINGLE");
					    if ((nNbr == 3) && (nAtm == 3))
					    {
						cause = "Connected atoms "+getAtomRef(atm,mol)+" "+getAtomRef(nbr,mol)+" have 3 connected atoms but no double bond. They are likely to be aromatic but no aromaticity has been reported";
						return cause;
/*
                		                return true;
*/
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
/*
            if (!atm.getFlag(CDKConstants.ISAROMATIC))
	        continue;
	    if (!atm.getSymbol().equals("C"))
		continue;
            if (atm.getFormalCharge() != 0)
                continue;
	    if (mol.getConnectedAtomsCount(atm) != 3)
		continue;
	    int n = 0;
	    for (IBond bnd : mol.getConnectedBondsList(atm))
	    {
		if (bnd.getOrder() == IBond.Order.valueOf("DOUBLE"))
		    n++;
	    }
	    if (n == 0)
	    {
		if (repOnScreen >= 1)
	 	    System.out.println("Aromatic atom "+atm.getSymbol()+mol.getAtomNumber(atm)+" has 3 connected atoms but no double bonds");
	        return true;
	    }
*/
	}

        return cause;
    }

//------------------------------------------------------------------------------
   
    /**
     * Returns the number of bonds, with a certain bond order, surrounding an atom
     */
    public static int numOfBondsWithBO(IAtom atm, IAtomContainer mol, String ord)
    {
	int res = -1;
        int n = 0;
        for (IBond bnd : mol.getConnectedBondsList(atm))
        {
            if (bnd.getOrder() == IBond.Order.valueOf(ord))
                n++;
        }

	res = n;
	return res;
    }

//------------------------------------------------------------------------------
    /**
     * Fix aromaticity for SD files where bond order=4 is used
     */
    public static IAtomContainer fixSDAromaticity(IAtomContainer mol)
    {
	IAtomContainer newMol = new AtomContainer();
	try {
            newMol = (IAtomContainer) mol.clone();

	    AtomContainerManipulator.percieveAtomTypesAndConfigureUnsetProperties(newMol);
/*
	    for (IAtom atm : newMol.atoms())
	    {
//System.out.print("ATOM "+atm.getSymbol()+newMol.getAtomNumber(atm)+" (Ar: "+atm.getFlag(CDKConstants.ISAROMATIC)+") ");
//if (atm.getFlag(CDKConstants.VISITED))
//    continue;
//System.out.print("Hybr. '"+atm.getHybridization()+"' ");
		if (atm.getFlag(CDKConstants.ISAROMATIC))
		{
		    atm.setHybridization(org.openscience.cdk.interfaces.IAtomType.Hybridization.valueOf("SP2"));
//		    exploreAromaticSystem(atm,mol);
		}		    
//System.out.println(" becomes '"+atm.getHybridization()+"'");
	    }
*/
        } catch (Throwable t) {
            System.err.println("CLONING molecule failed! "+t);
	    t.printStackTrace();
            System.exit(0);
        }

	FixBondOrdersTool fbt = new FixBondOrdersTool();

	try {
	    newMol = fbt.kekuliseAromaticRings(newMol);
	} catch (Throwable t) {
            System.err.println("KEKULIZATION FAILED! "+t);
	    t.printStackTrace();
            System.exit(0);
        }

/*	
	for (IBond b : newMol.bonds())
	{
	    System.out.println("Bond "+b.getOrder());
	}
*/
	
	return newMol;
    }

//------------------------------------------------------------------------------
    /**
     * Append dummy atoms so to allow the use of internal coordinates.
     * Dummy atoms are connected to the central atom of a linear 
     * (or close-to-linear) bend.
     * @param mol the molecule to be modified
     */
    public static void addDummiesOnLinearities(IAtomContainer mol)
    {
        double angLim = Parameters.linearBendThld;
        for (IAtom atm : mol.atoms())
        {
            boolean fastWayOut = false;
            List<IAtom> nbrs = mol.getConnectedAtomsList(atm);
            for(int i=0; i<nbrs.size(); i++)
            {
                IAtom atmL = nbrs.get(i);
                for (int j=i+1; j<nbrs.size(); j++)
                {
                    IAtom atmR = nbrs.get(j);
                    double angle = calculateBondAngle(atmL,atm,atmR);
                    if (angle > angLim)
                    {
                        if (repOnScreen >= 3)
                            System.out.println("+Adding dummy atom on atom "
                                    +getAtomRef(atm,mol));
    
                        // Create dummy atom
                        IAtom dummyAtm = getDummyInSafeDirection(atm,atmR,mol);
                        mol.addAtom(dummyAtm);
    
                        // Connect dummy atom
                        IBond dummyBnd = new Bond(dummyAtm,atm);
                        mol.addBond(dummyBnd);
    
                        fastWayOut = true;
                    }
    
                    if (fastWayOut)
                        break;
                }
                if (fastWayOut)
                    break;
            }
        }
    }

//------------------------------------------------------------------------------
    /**
     * Append dummy atoms so to allow the use of internal coordinates.
     * Dummy atoms are connected to the central atom of a linear 
     * (or close-to-linear) bend. Attachment points are also taken into account
     * @param frag the fragment to be modified
     */
    public static void addDummiesOnLinearities(GM3DFragment frag)
    {
        double angLim = Parameters.linearBendThld;
	//Use the 'apOnMol' which is the fragment with extra atoms in lieu 
	//of attachment points.
	IAtomContainer apOnMol = frag.apOnMol;
        for (int atmIdx=0; atmIdx<apOnMol.getAtomCount(); atmIdx++)
        {
	    //Atom on 'apOnMol' (includes Du atoms from multihapto systems)
	    IAtom atm = apOnMol.getAtom(atmIdx);

            boolean fastWayOut = false;
            List<IAtom> nbrs = apOnMol.getConnectedAtomsList(atm);
            for(int i=0; i<nbrs.size(); i++)
            {
                IAtom atmL = nbrs.get(i);
                for (int j=i+1; j<nbrs.size(); j++)
                {
                    IAtom atmR = nbrs.get(j);
                    double angle = calculateBondAngle(atmL,atm,atmR);
                    if (angle > angLim)
                    {
                        if (repOnScreen >= 3)
                            System.out.println("Adding dummy atom on atom "
                                    +getAtomRef(atm,apOnMol)
				    +" (GM3DFragment.apOnMol atom count)");

                        // Create dummy atom
                        IAtom dummyAtm = getDummyInSafeDirection(atm,atmR,apOnMol);
			//WARNING! we work 'apOnMol' but we add dummy on the fragment!
                        frag.addAtom(dummyAtm);

                        // Connect dummy atom
			//WARNING! we work 'apOnMol' but dummy must be on 'frag'!
			//Therefore make sure you connect dummy to 'frag'
			//atoms "AP" are appended so atom numbereng is the same
                        IBond dummyBnd = new Bond(dummyAtm,frag.getAtom(atmIdx));
                        frag.addBond(dummyBnd);
			//Add also the bond on 'apOnMol'
			IBond dummyBndApOnMol = new Bond(dummyAtm,atm);
			apOnMol.addBond(dummyBndApOnMol);

                        fastWayOut = true;
                    }

                    if (fastWayOut)
                        break;
                }
                if (fastWayOut)
                    break;
            }
        }
    }

//------------------------------------------------------------------------------

    /**
     * Generates a dummy atom 0.1 nm from atom A in a place that is safe for
     * a dummy atom. The be safe the position of the dummy has to respect 2 
     * criteria: 
     * 1) bo not too close to an existing atom connected to A,
     * 2) bo not too close to be opposite to an atom connected to A, which
     * means, do not create linearities.
     * The method try to use 90 and 45 deg angles with respect to vector AB.
     * @param atmA first atom used to place the dummy
     * @param atmB second atom used to place the dummy
     * @param mol the molecular system
     * @return the dummy atom
     */
    public static IAtom getDummyInSafeDirection(IAtom atmA, IAtom atmB, IAtomContainer mol)
    {
        // Get coordinates: A central atom B the surely linear bond
        Point3d pA = getCoords3d(atmA);
        Point3d pB = getCoords3d(atmB);

        // Get the forbidden areas: those in proximity of
        // atoms connected to A, or opposite (180 deg) to them
        ArrayList<Vector3d> allBusyDirections =new ArrayList();
        List<IAtom> nbrs = mol.getConnectedAtomsList(atmA);
        for (IAtom nbr : nbrs)
        {
            Point3d pLoc = new Point3d();
            pLoc = getCoords3d(nbr);
            Vector3d vLoc = new Vector3d();
            vLoc = CartesianSpaceUtils.getVectorFromTo(pA,pLoc);
            vLoc.normalize();
            allBusyDirections.add(vLoc);

//Delete. it was for checking the code
/*
System.out.println("vLoc: "+vLoc);
Vector3d vLoc2  = new Vector3d(vLoc.x,vLoc.y,vLoc.z);
CartesianSpaceUtils.translateOrigin(vLoc2,pA);
IAtom dummyB = new Atom("H",new Point3d(vLoc2.x,vLoc2.y,vLoc2.z));
mol.addAtom(dummyB);
*/
            Vector3d vLocOpposite = new Vector3d(vLoc.x * -1.0,vLoc.y * -1.0,vLoc.z * -1.0);
            vLocOpposite.normalize();
            allBusyDirections.add(vLocOpposite);
/*
//Delete. it was for checking the code
System.out.println("vLocOpposite: "+vLocOpposite);
Vector3d vLocOpposite2 = new Vector3d(vLocOpposite.x,vLocOpposite.y,vLocOpposite.z);
CartesianSpaceUtils.translateOrigin(vLocOpposite2,pA);
IAtom dummyb = new Atom("H",new Point3d(vLocOpposite2.x,vLocOpposite2.y,vLocOpposite2.z));
mol.addAtom(dummyb);
*/
        }
        
        //Get vector perpendicular (normal) to AB
        Vector3d vAB = CartesianSpaceUtils.getVectorFromTo(pA,pB);
        vAB.normalize();
        Vector3d vNorm = CartesianSpaceUtils.getNormalDirection(vAB);
        vNorm.normalize();

        // Rotate the normal vector to check different angles
        ArrayList<Vector3d> allAttempts = new ArrayList<Vector3d>();
        ArrayList<Double> allAttemptsMinVal = new ArrayList<Double>();
        double angStep = 22.0;
        double maxStepD =  360.0 / angStep;
        int maxStep = (int) maxStepD;
        double forbiddenRadius = 0.2;
        for(int j=1; j<3; j++)
        {
//System.out.println("\n Using j: "+j+" "+allAttempts);
            for(int step = 0; step<maxStep; step++)
            {
		Vector3d vADuTry = new Vector3d();
		if (j==1){
		    vADuTry = new Vector3d(vNorm.x,vNorm.y,vNorm.z);
		} else if (j==2) {
                    vADuTry = CartesianSpaceUtils.getSumOfVector(vAB,
				new Vector3d(vNorm.x*(Math.sqrt(2.0)),
                                                vNorm.y*(Math.sqrt(2.0)),
						vNorm.z*(Math.sqrt(2.0))));
		}
                vADuTry.normalize();
	
                // Get the new candidate position by rotation
		double ang = angStep * step;
                CartesianSpaceUtils.rotatedVectorWAxisAngle(vADuTry,vAB,ang);
//System.out.println("Attempt "+step+" angle: "+angStep*step);
            
                // Check if the candidate position is too close to a forbidden place
                // Use distance since they are all originating from A and normalized
                boolean skip = false;
                double min = 10.0;
                for (int i=0; i<allBusyDirections.size(); i++)
                {
                    Vector3d busyDir = allBusyDirections.get(i);
                    Vector3d diffVec = CartesianSpaceUtils.getDiffOfVector(vADuTry,busyDir);
                    double l = diffVec.length();
//System.out.println("diffLength :"+diffVec.length()+" step:"+step+" i:"+i);
                    if (l < forbiddenRadius)
                    {
                        skip = true;
                        break;
                    }
                    if (l<min)
                        min = l;
                }

                if (skip)
                    continue;

                // Store the surviving ones:
                // those that are not too close to forbidden regions;
                allAttempts.add(vADuTry);
                allAttemptsMinVal.add(min);
//System.out.println("ATTEMPT "+allAttemptsMinVal.indexOf(min)+" Min: "+min+"\n\t\tVEC: "+vADuTry);
            } //end of loop over angles around AB (torsion of AB)
        } //end of loop over angle with AB

        // Find the best candidate: the most distant from forbidden areas
        double max = Collections.max(allAttemptsMinVal);
        int best = allAttemptsMinVal.indexOf(max);
        Vector3d vADu = new Vector3d();
        vADu = allAttempts.get(best);

        // Create the dummy in the best position found
        CartesianSpaceUtils.translateOrigin(vADu,pA);
        Point3d duP3dB = new Point3d(vADu.x, vADu.y, vADu.z);
        IAtom dummyAtm = new Atom(Parameters.duSymbol,duP3dB);
        mol.addAtom(dummyAtm);

//IOtools.writeSDFAppend("checkthis.sdf",mol,true);
//System.out.println("END of SafeDirection");
//IOtools.pause();

        return dummyAtm;

    }

//------------------------------------------------------------------------------

    /**
     * Generates a dummy atom 0.1 nm from atom A and forming a bond angle of 
     * 45.0 degree with atom B
     * @param atmA first atom used to place the dummy
     * @param atmB second atom used to place the dummy
     * @return the dummy atom
     */
    public static IAtom getDummyIn45DegreeDirection(IAtom atmA, IAtom atmB)
    {
        Point3d pA = getCoords3d(atmA);
        Point3d pB = getCoords3d(atmB);

        // Vector from A to B
        Vector3d vAB = CartesianSpaceUtils.getVectorFromTo(pA,pB);

        // Get position of Du
        Vector3d norm = CartesianSpaceUtils.getNormalDirection(vAB);
        vAB.normalize();
        Vector3d vADu = CartesianSpaceUtils.getSumOfVector(vAB,norm);
        vADu.normalize();

        // Translate 
        vADu.x = vADu.x + pA.x;
        vADu.y = vADu.y + pA.y;
        vADu.z = vADu.z + pA.z;

        Point3d duP3d = new Point3d(vADu.x, vADu.y, vADu.z);

        //Create new dummy atom
        IAtom dummy = new Atom(Parameters.duSymbol,duP3d);

        return dummy;
    }

//------------------------------------------------------------------------------

    /**
     * Generates a dummy atom 0.1 nm from atom A and forming a bond angle of 
     * 90.0 degree with atom B
     * @param atmA first atom used to place the dummy
     * @param atmB second atom place the dummy
     * @return the dummy atom
     */
    public static IAtom getDummyInPerpendicularDirection(IAtom atmA, IAtom atmB)
    {
        Point3d pA = getCoords3d(atmA);
        Point3d pB = getCoords3d(atmB);

        // Vector from A to B
        Vector3d vAB = CartesianSpaceUtils.getVectorFromTo(pA,pB);

        // Get position of Du
        Vector3d vADu = CartesianSpaceUtils.getNormalDirection(vAB);

        // Translate 
        vADu.x = vADu.x + pA.x;
        vADu.y = vADu.y + pA.y;
        vADu.z = vADu.z + pA.z;

        Point3d duP3d = new Point3d(vADu.x, vADu.y, vADu.z);

        //Create new dummy atom
        IAtom dummy = new Atom(Parameters.duSymbol,duP3d);

        return dummy;
    }

//------------------------------------------------------------------------------

    /**
     * Get Cartesian coordinates in 3D space of an atom (even if this is in 2D)
     */
    public static Point3d getCoords3d(IAtom atom)
    {
        Point3d p3d = new Point3d();
        try {
            Point2d atp2d = new Point2d();
            atp2d = atom.getPoint2d();
            p3d.x = atp2d.x;
            p3d.y = atp2d.y;
            p3d.z = 0.0000;
        } catch (Throwable t) {
            Point3d atp3d = new Point3d();
            atp3d = atom.getPoint3d();
            p3d.x = atp3d.x;
            p3d.y = atp3d.y;
            p3d.z = atp3d.z;            
        }
        return p3d;
    }

//-----------------------------------------------------------------------------

    /**
     * Determines the dimensionality of the chemical object submitted
     * @param mol input molecular/chemical object
     * @return dimensionality of this object (2 or 3) or -1
     */
    public static int getDimensions(IAtomContainer mol)
    {
        final int is2D = 2;
        final int is3D = 3;
        final int not2or3D = -1;

        int numOf2D = 0;
        int numOf3D = 0;

        for (IAtom atm : mol.atoms())
        {
            Point2d p2d = new Point2d();
            Point3d p3d = new Point3d();
            p2d = atm.getPoint2d();
            boolean have2D = true;
            if (p2d == null)
            {
                have2D = false;
                p3d = atm.getPoint3d();
                if (p3d == null)
                {
                    return not2or3D;
                }
            }
            ArrayList<Double> pointer = new ArrayList<Double>();
            try {
                if (have2D)
                {
                      pointer.add(p2d.x);
                    pointer.add(p2d.y);
                    numOf2D++;
                } else {
                    pointer.add(p3d.x);
                    pointer.add(p3d.y);
                    pointer.add(p3d.z);
                    numOf3D++;
                }
            } catch (Throwable t) {
                return not2or3D;
            }
        }

        if (numOf2D == mol.getAtomCount())
            return is2D;
        else if (numOf3D == mol.getAtomCount())
            return is3D;
        else
            return not2or3D;
    }

//-----------------------------------------------------------------------------

    /**
     * Looks for a proprty referring to the name or ID of the molecule.
     * Recognized cdk:Title,ChEBI ID,TTD DRUGID
     * @param mol molecule
     * @return the name or ID if any. Otherwise 'noname'
     */

    public static String getNameOrID(IAtomContainer mol)
    {
        String name = "noname";

        //ChEBI
        try {
            name = mol.getProperty("ChEBI ID").toString();
        } catch (Throwable t1) {
            //TTD
            try {
                name = mol.getProperty("DRUGID").toString();
            } catch (Throwable t2) {
                //CDK
                try {
                    name = mol.getProperty("cdk:Title").toString();
                } catch (Throwable t3) {
                    //General case using title
                    try {
                        name = mol.getProperty(CDKConstants.TITLE).toString();
                    } catch (Throwable t) {
                        if (repOnScreen >= 3)
                            System.out.println("Molecule name not found. Set to "+name);
                    }
                }
            }
        }

        return name;
    }

//------------------------------------------------------------------------------

    /**
     * Calculate the bond angle given the coordinates
     * @param atomLeft
     * @param atomCentre
     * @param atomRight
     * @return the bond angle
     */

    public static double calculateBondAngle(IAtom atomLeft, IAtom atomCentre, IAtom atomRight)
    {
        double angle = 0.0;
        Point3d l3d = getCoords3d(atomLeft);
        Point3d c3d = getCoords3d(atomCentre);
        Point3d r3d = getCoords3d(atomRight);

        double xab = l3d.x - c3d.x;
        double yab = l3d.y - c3d.y;
        double zab = l3d.z - c3d.z;
        double xcb = r3d.x - c3d.x;
        double ycb = r3d.y - c3d.y;
        double zcb = r3d.z - c3d.z;
        double rab2 = xab*xab + yab*yab + zab*zab;
        double rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
        double rabc = Math.sqrt(rab2 * rcb2);
        if (rabc != 0.0)
        {
            double cosine = (xab*xcb + yab*ycb + zab*zcb) / rabc;
            cosine = Math.min(1.0, Math.max(-1.0, cosine));
            angle = 57.29577951308232088 * Math.acos(cosine);
        }
        return angle;
    }

//-----------------------------------------------------------------------------

    /**
     * Translate integer to a CDK bond order
     * @param bndOrder number defining the bond order
     * @return <code>IBond.Order</code> corresponding to integer bond order
     */
    public static Order intToBondOrder(int bndOrder)
    {
        if (bndOrder == 1)
            return IBond.Order.valueOf("SINGLE");
        else if (bndOrder == 2)
            return IBond.Order.valueOf("DOUBLE");
        else if (bndOrder == 3)
            return IBond.Order.valueOf("TRIPLE");
        else if (bndOrder == 4)
            return IBond.Order.valueOf("QUADRUPLE");
        else if (bndOrder == 0)
            return IBond.Order.valueOf("UNSET");
        else {
            System.err.println("Trying to get a not defined bond order!");
            System.exit(0);
            return IBond.Order.valueOf("UNSET");
        }
    }

//-----------------------------------------------------------------------------

    /**
     * Translate a CDK bond order to an integer number
     * @param bndOrder string-like CDK bond order
     * @return integer corresponding to CDK bond order
     */
    public static int bondorderToint(IBond.Order bndOrder)
    {
        if (bndOrder.equals(IBond.Order.valueOf("SINGLE")))
            return 1;
        else if (bndOrder.equals(IBond.Order.valueOf("DOUBLE")))
            return 2;
        else if (bndOrder.equals(IBond.Order.valueOf("TRIPLE")))
            return 3;
        else if (bndOrder.equals(IBond.Order.valueOf("QUADRUPLE")))
            return 4;
        else if (bndOrder.equals(IBond.Order.valueOf("UNSET")))
            return 0;
        else {
            System.err.println("Trying to get a not defined bond order!");
            System.exit(0);
            return -1;
        }
   }

//-----------------------------------------------------------------------------

    /**
     * Analysis of CDK Flags over the whole AtomContainer
     * @param mol molecular system to check
     * @param flagid number identifying the flag to check
     * @return <code>true</code> if at least one flag in <code>false</code>
     */

    public static boolean containsFalseFlag(IAtomContainer mol, int flagid)
    {
        for (IAtom a : mol.atoms())
        {
            if (!a.getFlag(flagid))
                return true;
        }
        return false;
   }

//------------------------------------------------------------------------------
}
