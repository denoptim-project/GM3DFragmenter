import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import java.util.Collections;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @author Marco Foscato (University of Bergen) Discussions: Giovanni Occhipinti and Vidar R. Jensen (University of Bergen)
 */

public class GM3DLigand
{
    //Molecule
    private IAtomContainer mol = new AtomContainer();

    //Root of the Ligand
    private IAtom root;

    //First atom of the ligand
    private IAtom seed;

//
//   ...         <...part of the ligand...>
//    /         /
//  <root>---<seed>---<...part of the ligand...>
//    \         \ 
//    ...        <...part of the ligand...>
//

    //Order of atoms
    private Map<Integer,List<IAtom>> orderAtoms = new HashMap<Integer,List<IAtom>>();

    //Flag for reporting
    private int repOnScreen = 0;
    private String pre = "GM3DLigand: ";
    //Recursion flag
    static int recNum = 1;

    //Flags per atoms
    private ArrayList<ArrayList<Boolean>> flags = new ArrayList<ArrayList<Boolean>>();


//------------------------------------------------------------------------------

/**
 * TODO
 */

    public GM3DLigand(IAtom seed, IAtom root, IAtomContainer mol)
    {
	this.seed = seed;
        this.root = root;
        this.mol = mol;
	
	exploreLigand();
    }

//------------------------------------------------------------------------------
    public int getFreeAtomsFlag()
    {
        int freeFlag = -1;
        //create a vector with false entries
        int atoms = this.mol.getAtomCount();
        ArrayList<Boolean> flg = new ArrayList<Boolean>();
        for (int i = 0; i<atoms; i++)
                flg.add(false);

        //add the vector to the list of flags
        flags.add(flg);
        freeFlag = flags.indexOf(flg);

        return freeFlag;
    }

//------------------------------------------------------------------------------
    public void setAtomsFlag(int flagID, int atm, boolean value)
    {
        this.flags.get(flagID).set(atm,value);
    }

//------------------------------------------------------------------------------
    public boolean getAtomsFlag(int flagID, int atm)
    {
        return this.flags.get(flagID).get(atm);
    }
//------------------------------------------------------------------------------
    public boolean getAtomsFlag(int flagID, IAtom atm)
    {
	int atmindex = this.mol.getAtomNumber(atm);
        return this.flags.get(flagID).get(atmindex);
    }

//------------------------------------------------------------------------------
    public void deleteAtomsFlag(int flagID)
    {
        this.flags.remove(flagID);
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */
    private void exploreLigand()
    {
        //set string for reporting and debugging
        String recFlag = "";
        for (int ri = 0; ri < recNum; ri++)
             recFlag = recFlag+"-";

	if (repOnScreen >= 3)
	    System.out.println(pre+recFlag+" Exploring GM3DLigand");

	int flag = getFreeAtomsFlag();

	//root is already visited by definition of GM3DLigand
	setAtomsFlag(flag,mol.getAtomNumber(root),true);
	addAtomToOrder(seed, flag, "");

        //explore the rest of the connected system if any
        int conAtms = this.mol.getConnectedAtomsCount(seed);
        if (conAtms > 1)
	{
	    recNum++;
	    exploreSubLigand(seed,flag);
	    recNum--;
	}

        if (repOnScreen >= 3)
            System.out.println(pre+recFlag+" Number of reported levels: "+orderAtoms.keySet().size());
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */
    private void exploreSubLigand(IAtom seedTree, int flag)
    {
        //set string for reporting and debugging
        String recFlag = "";
        for (int ri = 0; ri < recNum; ri++)
             recFlag = recFlag+"-";

	//get list of neghbours
	List<IAtom> neigbourAtoms = this.mol.getConnectedAtomsList(seedTree);

	//Add this shell of atoms
	List<IAtom> todoNeigbours = new ArrayList<IAtom>();
        for (IAtom connectedAtom : neigbourAtoms)
	    if (!getAtomsFlag(flag, connectedAtom))
	    {
		todoNeigbours.add(connectedAtom);
	        addAtomToOrder(connectedAtom, flag, recFlag);
	    }

        for (IAtom connectedAtom : todoNeigbours)
        {
            // move to the next shell of atoms
            if (this.mol.getConnectedAtomsCount(connectedAtom) > 1)
            {
		recNum++;
		if (repOnScreen >= 3)
                    System.out.println(pre+recFlag+"> resursion on atom "+this.mol.getAtomNumber(connectedAtom)+connectedAtom.getSymbol());
		exploreSubLigand(connectedAtom,flag);
                recNum--;
            }
	}
    }
//------------------------------------------------------------------------------
    private void addAtomToOrder(IAtom atm, int flagID, String recFlag)
    {
	int atmidx = this.mol.getAtomNumber(atm);
	setAtomsFlag(flagID,atmidx,true);
	int key = recNum-1;
	if (!orderAtoms.containsKey(key))
	{
	    List<IAtom> lst = new ArrayList<IAtom>();
	    lst.add(atm);
	    orderAtoms.put(key,lst);
	} else 
	{
	    orderAtoms.get(key).add(atm);

            //Ordering according to the following criteria
            // MASS NUMBER in DESCENDING ORDER (low index == high mass number)
            int atNumIn;
            if (atm.getSymbol().equals("AP") || atm.getSymbol().equals("Du"))
                atNumIn = 0;
            else
                atNumIn = atm.getMassNumber();

            // NUMBER OF CONNECTIONS in DESCENDING ORDER (low index == many connected atoms)
            int conIn = this.mol.getConnectedAtomsCount(atm);

            //Now compare entries one-by-one and, in case, swap them
            int size = orderAtoms.get(key).size();
            size = size -1;
            for (int i = size; i >= 0; i--)
            {
                //Identify comparing atom
                IAtom aout = orderAtoms.get(key).get(i);

//TODO remove: for debug only
//System.out.println("For atom: "+aout.getSymbol()+this.mol.getAtomNumber(aout)+"res = "+aout.getSymbol().equals("AP")+aout.getSymbol().equals("Du"));

                int atNumOut;
                if (aout.getSymbol().equals("AP") || aout.getSymbol().equals("Du"))
                    atNumOut = 0;
                else {

//TODO remove: for debug only
//System.out.println("For atom: "+aout.getSymbol()+this.mol.getAtomNumber(aout));

                    atNumOut = aout.getMassNumber();
		}

                int conOut = this.mol.getConnectedAtomsCount(aout);
                //One-by-one comparison
                if (atNumIn > atNumOut)
                {
                    //Swap based on mass number
                    int aid = orderAtoms.get(key).indexOf(atm);
                    Collections.swap(orderAtoms.get(key),aid,i);
                } else if (atNumIn < atNumOut)
                {
                    break;
                } else
                {
                    if (conIn > conOut)
                    {
                        //Swap based on number of connections
                        int aid = orderAtoms.get(key).indexOf(atm);
                        Collections.swap(orderAtoms.get(key),aid,i);
                    } else if (conIn < conOut)
                    {
                    break;
                    } else
                    {
			if (repOnScreen >= 1)
			{
                            System.out.println(pre+recFlag+" WARNING! EQUAL atoms or Atom-Atom comparison NOT COVERED");
                            System.out.println(pre+recFlag+" Atom numbers: "+this.mol.getAtomNumber(atm)+atm.getSymbol()+" and "+this.mol.getAtomNumber(aout)+aout.getSymbol());
			}
                    }
                }
            }
	}
	if (repOnScreen >= 3)
	   System.out.println(pre+recFlag+"Adding Ordered Atom "+atmidx+atm.getSymbol()+" to level "+key);
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */
    public IAtom getSeed()
    {
	return seed;
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */
    public IAtomContainer getMol()
    {
        return mol;
    }

//------------------------------------------------------------------------------
/**
 * TODO
 */
    public Map<Integer,List<IAtom>> getOrderedList()
    {
        return orderAtoms;
    }
//------------------------------------------------------------------------------
}
