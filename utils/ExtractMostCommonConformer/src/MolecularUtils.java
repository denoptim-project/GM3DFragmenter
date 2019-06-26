/**
 * Toolbox for molecular objects
 * 
 * @author Marco Foscato (University of Bergen)
 */

import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import org.openscience.cdk.tools.periodictable.PeriodicTable;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.Bond;
import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class MolecularUtils
{

    //Reporting flag
    private static int repOnScreen = 0;
    private static String duSymbol = "Du";
    private static double linearBendThld = 175.0;

//------------------------------------------------------------------------------

/**
 * Counts the number of atoms for each element in the system
 */

    public static  Map<String,Integer> getMolecularFormulaOfAtomContainer(IAtomContainer mol)
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
                    if (nbr.getSymbol().equals(duSymbol))
                    {
                        List<IAtom> neighboursOfDu = mol.getConnectedAtomsList(nbr);
                        visited.put(nbr,-1);
                        nextLevel.remove(nbr);
                        nextLevel.removeAll(neighboursOfDu); //remove possible duplicates
                        nextLevel.addAll(neighboursOfDu); //add new candidates
                    }
                }
                nextLevel.removeAll(neighbours); //remove possible duplicates
                nextLevel.addAll(neighbours); //add new candidates
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

//------------------------------------------------------------------------------
    /**
     * Append dummy atoms so to allow the use of internal coordinates.
     * Dummy atoms are connected to the central atom of a linear 
     * (or close-to-linear) bend.
     * @param mol the molecule to be modified
     */
    public static void addDummiesOnLinearities(IAtomContainer mol)
    {
        double angLim = linearBendThld;
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
/*
// TODO remove
System.out.println("CHECK 0");
IOtools.writeSDFAppend("checkthis.sdf",mol,true);
IOtools.pause();
*/
    }

//------------------------------------------------------------------------------

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

//TODO delete. it was for checking the code
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
//TODO delete. it was for checking the code
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
        IAtom dummyAtm = new Atom(duSymbol,duP3dB);
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
        IAtom dummy = new Atom(duSymbol,duP3d);

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
        IAtom dummy = new Atom(duSymbol,duP3d);

        return dummy;
    }

//------------------------------------------------------------------------------

    /**
     * Makes sure that all atoms have a <code>Point3d</code> representation.
     * CDK assumes that z-coords=0.0000 means 2d rathar than 3d coordinates
     * thus we apply a small distorsion to by-pass this issue.
     */
    public static void ensure3d(IAtomContainer mol)
    {
	for (IAtom a : mol.atoms())
	{
	    Point3d p3d = getCoords3d(a);
	    a.setPoint3d(p3d);
	}
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
     **/

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
     * Add connections between atoms if distance is below summ of vdW radii
     * Does not remove existing bonds.
     */

    public static void addConnectionsByVDWRadius(IAtomContainer mol, String el, double tolerance)
    {
	PeriodicTable pt = new PeriodicTable();
        for (IAtom atmA : mol.atoms())
        {
            if (!atmA.getSymbol().equals(el))
                continue;

	    //get van der Waals radius
	    double rA = 1.5;
	    try {
	        rA = pt.getVdwRadius(el);
	    } catch (Throwable t) {
		System.out.println("WARNING! Element "+el+" not in Periodic Table. Using "+rA+" as van der Waals radius.");
	    }

	    //Already connected atoms
	    List<IAtom> nbrs = mol.getConnectedAtomsList(atmA);

            for (IAtom atmB : mol.atoms())
            {
		if (atmB.equals(atmA))
		    continue;

		//get van der Waals radius
		double rB = 1.5;
		String sB = atmB.getSymbol();
                try {
                    rB = pt.getVdwRadius(sB);
                } catch (Throwable t) {
                    System.out.println("WARNING! Element "+sB+" not in Periodic Table. Using "+rB+" as van der Waals radius.");
                }
		//Evaluate possibility of generatin a new bond
                if (!nbrs.contains(atmB))
                {
                    double dist = MolecularUtils.calculateInteratomicDistance(atmA,atmB);
		    double refDist = rA + rB;
		    refDist = refDist - (refDist * tolerance);
//System.out.println("dist: "+dist+" ref:"+refDist+" vdWSum:"+(rA + rB));
                    if (dist < refDist)
                    {
                        IBond b = new Bond(atmA,atmB,IBond.Order.valueOf("SINGLE"));
			mol.addBond(b);
//System.out.println("Added new bond: "+getAtomRef(atmA,mol)+"-"+getAtomRef(atmB,mol));
//System.out.println("CASE dist: "+dist+" ref:"+refDist+" vdWSum:"+(rA + rB));
//IOtools.pause(); 
                    }
                }
            }
        }
    }

//------------------------------------------------------------------------------

    /**
     * Calculate the bond angle given the coordinates
     * @param atmA
     * @param atmB
     * @return distance between atmA and atmA
     */

    public static double calculateInteratomicDistance(IAtom atmA, IAtom atmB)
    {
        Point3d pa = getCoords3d(atmA);
        Point3d pb = getCoords3d(atmB);
	double dx = pa.x - pb.x;
        double dy = pa.y - pb.y;
        double dz = pa.z - pb.z;
	double dist = Math.sqrt((Math.pow(dx,2.0D)) +
				(Math.pow(dy,2.0D)) +
                                (Math.pow(dz,2.0D)));
	return dist;
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

//------------------------------------------------------------------------------

    /**
     * Calculate the dihedral angle A-B-C-D given the atomd
     *
     * @param A atom in position A
     * @param B atom in position B
     * @param C atom in position C
     * @param D atom in position D
     * @return the dihedral angle
     */
    public static double calculateTorsionAngle(IAtom atmA, IAtom atmB, IAtom atmC, IAtom atmD)
    {
	Point3d pA = getCoords3d(atmA);
        Point3d pB = getCoords3d(atmB);
        Point3d pC = getCoords3d(atmC);
        Point3d pD = getCoords3d(atmD);

	double[] A = new double[] {pA.x,pA.y,pA.z};
	double[] B = new double[] {pB.x,pB.y,pB.z};
	double[] C = new double[] {pC.x,pC.y,pC.z};
	double[] D = new double[] {pD.x,pD.y,pD.z};
        double angle = 0.0D;
        double xba = B[0] - A[0];
        double yba = B[1] - A[1];
        double zba = B[2] - A[2];
        double xcb = C[0] - B[0];
        double ycb = C[1] - B[1];
        double zcb = C[2] - B[2];
        double xdc = D[0] - C[0];
        double ydc = D[1] - C[1];
        double zdc = D[2] - C[2];
        double xt = yba * zcb - ycb * zba;
        double yt = xcb * zba - xba * zcb;
        double zt = xba * ycb - xcb * yba;
        double xu = ycb * zdc - ydc * zcb;
        double yu = xdc * zcb - xcb * zdc;
        double zu = xcb * ydc - xdc * ycb;
        double rt2 = xt * xt + yt * yt + zt * zt;
        double ru2 = xu * xu + yu * yu + zu * zu;
        double rtru = Math.sqrt(rt2 * ru2);
        if (rtru != 0.0)
        {
            double cosine = (xt * xu + yt * yu + zt * zu) / rtru;
            cosine = Math.min(1.0, Math.max(-1.0, cosine));
            angle = 57.29577951308232088 * Math.acos(cosine);
            double sign = xba * xu + yba * yu + zba * zu;
            if (sign < 0.0)
            {
                angle = -angle;
            }
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
     **/

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
