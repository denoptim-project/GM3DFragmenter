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
import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import java.util.Collections;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import net.sf.jniinchi.INCHI_RET;

import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smsd.Isomorphism;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIGenerator;

/**
 * GM3DFragment is a portion of a molecular object containing information
 * on the position and type of connections towards other fragments.
 * 
 * @author Marco Foscato (University of Bergen) Discussions: Giovanni Occhipinti and Vidar R. Jensen (University of Bergen)
 */

public class GM3DFragment extends AtomContainer implements IAtomContainer
{ 
    //List of attachment points
    private ArrayList<GM3DAttachmentPoint> allAPs = new ArrayList<GM3DAttachmentPoint>();

    //AtomContainer repreentation of fragment with attachment points as atoms
    /**
     * Alternative representation of the fragment in which the
     * attachment points are represented by dummy atoms
     */
    public IAtomContainer apOnMol;

    //Alternative atom order
    private Map<GM3DAttachmentPoint,ArrayList<Integer>> atomOrders = new HashMap<GM3DAttachmentPoint,ArrayList<Integer>>();

    //Flags per atoms
    private ArrayList<ArrayList<Boolean>> flags = new ArrayList<ArrayList<Boolean>>();

    //Dimensionality: 2 or 3 (-1 if error)
    private int dimensions;

    //Utilities
    //Level or information printed on screen
    private int repOnScreen = 0;
    String pre = "GM3DFragment: ";
    private static String checkfile = "check-GM3DFragment_"+Parameters.getJobName()+".sdf";
    //Recursion flag for reporting infos
    static int recNum = 1;

//------------------------------------------------------------------------------

    public GM3DFragment()
    {
        super();
        repOnScreen = Parameters.report;
    }
//------------------------------------------------------------------------------

    public GM3DFragment(IAtomContainer mol)
    {
        super(mol);
        repOnScreen = Parameters.report;

	dimensions = getDimensions(mol);

        //Re-elaboration of information on attachment points
        //update atom number
        updateAPs();
        //generate a list of objects GM3DAttachmentPoint
        extractAllAPs();

        //make the AtomContainer representation having 
        //attachment points as atoms with symbol "AP"
        apOnMol = new AtomContainer(mol);
        makeAPPseudoAtoms();
    }

//------------------------------------------------------------------------------

    public GM3DFragment(IAtomContainer mol, String format)
    {
        super(mol);
        repOnScreen = Parameters.report;

        dimensions = getDimensions(mol);

	//copy paste properties
	this.setProperties(mol.getProperties());

        //Re-elaboration of information on attachment points
	movePropertyToAP(format);

        //make the AtomContainer representation having 
        //attachment points as atoms with symbol "AP"
        apOnMol = new AtomContainer(mol);
        makeAPPseudoAtoms();
    }

//------------------------------------------------------------------------------

    public int getNumberOfAttachmentPoints()
    {
        return allAPs.size();
    }

//------------------------------------------------------------------------------

    public int getDimensions()
    {
        return dimensions;
    }

//------------------------------------------------------------------------------

/**
 * Compare to another <code>GM3DFragment</code>. 
 * <strong>WARNING! Needs fixed aromaticity to work properly</strong>.
 *
 * @param other the <code>GM3DFragment</code> to be compared with this one
 * @return <code>true</code> if <code>this</code> correspond to the same 
 * fragment of <code>other</code>
 */

    public boolean sameFragOf(GM3DFragment other)
    {
        boolean areEqual = false;
        //Set preCompString for reporting messages on screen
        String nameThis = "";
//TODO update to use method from MolecularUtils
        try { 
            nameThis = this.getProperty("cdk:Title").toString();
        } catch (Throwable t) {
            nameThis = "noname";
        }
        String nameOther = "";
        try { 
            nameOther = other.getProperty("cdk:Title").toString();
        } catch (Throwable t) {
            nameOther = "noname";
        }
        String preComp = "Comparing GM3DFragments: "+nameThis+" "+nameOther+" => ";

        //Are both fragments 'real fragments' having at leas one attchment point?
        if (!(this.getNumberOfAttachmentPoints() > 0) || !(other.getNumberOfAttachmentPoints() > 0))
        {
            System.err.println("\nERROR! Attempt to compare a GM3Dfragment with a molecule without attachment points!\n");
            System.exit(-1);
        }

        //Number of Atoms
        if (this.getAtomCount() != other.getAtomCount())
        {
            if (repOnScreen >= 2)
                System.out.println(preComp+"Different number of atoms.");
            return false;
        }

        //Number of Bonds
        if (this.getBondCount() != other.getBondCount())
        {
            if (repOnScreen >= 2)
                System.out.println(preComp+"Different number of bonds.");
            return false;
        }

        //Compare molecular structure and AP location

        //Couples of similar AP on 'this' and 'other'
        ArrayList<ArrayList<GM3DAttachmentPoint>> likeAPs = new ArrayList<ArrayList<GM3DAttachmentPoint>>();
        //Prepare list of classes (used for the next comparison approach using InChi)
        Set<String> allClasses = new HashSet<String>();

        //Identify similar Attachment Points
        for (int i = 0; i<allAPs.size(); i++)
        {
            GM3DAttachmentPoint ap = allAPs.get(i);
            for (int j = 0; j<other.getAllAPs().size(); j++)
            {
                GM3DAttachmentPoint otherap = other.getAllAPs().get(j);
                if (ap.getAPClass().equals(otherap.getAPClass()))
                {
                    ArrayList<GM3DAttachmentPoint> couple = new ArrayList<GM3DAttachmentPoint>();
                    couple.add(ap);
                    couple.add(otherap);
                    likeAPs.add(couple);
		    allClasses.add(ap.getAPClass());
                }
            }
        }

        //Compare number of similar AP on the two fragments
        if (likeAPs.size() == 0)
        {
            if (repOnScreen >= 2)
                System.out.println(preComp+"No similar AP on the two fragments");
            return false;
        } else {
            if (likeAPs.size() < allAPs.size() || likeAPs.size() < other.getAllAPs().size())
            {
                if (repOnScreen >= 2)
                    System.out.println(preComp+"Not all the APs have their counterpart on the other fragment!");
                return false;
            }
            else
            {
                if (repOnScreen >= 3)
                    System.out.println(preComp+"Number of possibilities for similar AP: "+likeAPs.size());
                //Do all the APs have their counterpart on the other fragment?
                //Check APs from this fragment
                for (int i = 0; i<allAPs.size(); i++)
                {
                    boolean apHasCounterpart = false;
                    GM3DAttachmentPoint ap = allAPs.get(i);
                    for (int j = 0; j<likeAPs.size(); j++)
                    {
                        ArrayList<GM3DAttachmentPoint> twoAps = likeAPs.get(j);
                        if (twoAps.contains(ap))
                        {
                            apHasCounterpart = true;
                            break;
                        }
                    }
                    if (!apHasCounterpart)
                    {
                        if (repOnScreen >= 2)
                            System.out.println(preComp+"AP "+ap+" doesn't have a counterpart on the other fragment.");
                        return false;
                    }
                }
                //Check APs from other fragment
                for (int i = 0; i<other.getAllAPs().size(); i++)
                {
                    boolean apHasCounterpart = false;
                    GM3DAttachmentPoint ap = other.getAllAPs().get(i);
                    for (int j = 0; j<likeAPs.size(); j++)
                    {
                        ArrayList<GM3DAttachmentPoint> twoAps = likeAPs.get(j);
                        if (twoAps.contains(ap))
                        {
                            apHasCounterpart = true;
                            break;
                        }
                    }
                    if (!apHasCounterpart)
                    {
                        if (repOnScreen >= 2)
                            System.out.println(preComp+"AP "+ap+" doesn't have a counterpart on the this fragment.");
                        return false;
                    }
                }
            }
        }

//TODO do we need this?

        //Try common atom orther and compare atom by atom
//        for (int i = 0; i<likeAPs.size(); i++)
	if (likeAPs.size() == 1)
        {
	    int i = 0;
            if (repOnScreen >= 3)
                System.out.println(preComp+"Comparing atom order of two fragments:");

            ArrayList<GM3DAttachmentPoint> twoAps = likeAPs.get(i);
            ArrayList<Integer> Ord = this.getReorderAtoms(twoAps.get(0));
            ArrayList<Integer> otherOrd = other.getReorderAtoms(twoAps.get(1));

            boolean equalOrder = true;
            int fail = -1;
            //WARNING! looping on atoms AND AP (using apOnMol)
            for (int j = 0; j<this.apOnMol.getAtomCount(); j++)
            {
                if (repOnScreen >= 3)
                {
                    System.out.println(preComp+"Comparing entry "+j+": "+
                    Ord.get(j)+this.apOnMol.getAtom(Ord.get(j)).getSymbol()+" "+
                    otherOrd.get(j)+other.apOnMol.getAtom(otherOrd.get(j)).getSymbol());
                }

                IAtom a = this.apOnMol.getAtom(Ord.get(j));
                IAtom aOther = other.apOnMol.getAtom(otherOrd.get(j));

                //Compare atoms
                GM3DAtomComparator comAtms = new GM3DAtomComparator();
                int res = comAtms.compareAtoms(a,this.apOnMol,aOther,other.apOnMol);

                //Compare atom number in the order of atoms
                if (res != 0)
                {
                    equalOrder = false;
                    fail = j;
                    break;
                }
            }

            if (!equalOrder)
            {
                if (repOnScreen >= 2)
                {
                    System.out.println(preComp+"Different atom: "+
                Ord.get(fail)+this.apOnMol.getAtom(Ord.get(fail)).getSymbol()+" "+
                otherOrd.get(fail)+other.apOnMol.getAtom(otherOrd.get(fail)).getSymbol());
                }
                return false;
            }
        }

        //Use InChi on 'apOnMol' to compare the two fragments
        String inchiThis = "";
        String inchiOther = "";
        boolean inchiDone = false;
        Map<String,String> classToElement = getClassToElement(allClasses,this,other);
        inchiThis = this.getInChiForFragmentWithAP(classToElement);
        inchiOther = other.getInChiForFragmentWithAP(classToElement);
	if (inchiThis.equals("") || inchiOther.equals(""))
	{
	    if (repOnScreen >= 1)
                System.out.println(preComp+"INCHI representation FAILED. Fragments will be sees as equal! (see "+checkfile+") ");
		rejectMol(this.apOnMol,"INCHI FAILED");
		inchiDone = false;
	} else {
	    inchiDone = true;
	}
        if (inchiDone & !(inchiThis.equals(inchiOther)))
        {
            if (repOnScreen >= 2)
                System.out.println(preComp+"Not same InChi representation - Fragments differ."); 
            return false;
        } 

        //All previous attempts to distinguish the two fragments failed: fragments do not differ!
        if (repOnScreen >= 2)
            System.out.println(preComp+" All criteria satisfied! The two GM3DFragments are equal.");

        return true;

    }
//------------------------------------------------------------------------------

    public ArrayList<Integer> getReorderAtoms(GM3DAttachmentPoint ap)
    {
        String name = this.getProperty("cdk:Title").toString();
        if (atomOrders.keySet().contains(ap))
        {
            if (repOnScreen >= 3)
                System.out.println("Atom Order of "+name+" => "+this.atomOrders.get(ap));
            return atomOrders.get(ap);
        }

        if (repOnScreen >= 3)
            System.out.println(pre+"Reordering atoms of "+name);

        int flag = getFreeAtomsFlag();

        ArrayList<Integer> reorderedAtms = new ArrayList<Integer>();
        this.atomOrders.put(ap,reorderedAtms);

        //explore the molecule starting from the AP
        int seed = getAtomNumberOfAP(this.apOnMol,ap);
        recNum = 1;
        exploreMolecule(seed,flag,ap);
        deleteAtomsFlag(flag);

        if (repOnScreen >= 3)
            System.out.println(pre+"Atom Order of "+name+" => "+this.atomOrders.get(ap));

        return this.atomOrders.get(ap);
    }

//------------------------------------------------------------------------------
    public int getFreeAtomsFlag()
    {
        int freeFlag = -1;
        //create a vector with false entries
        int atoms = this.apOnMol.getAtomCount();
        ArrayList<Boolean> flg = new ArrayList<Boolean>();
        for (int i = 0; i<atoms; i++)
                flg.add(false);

        //add the vector to the list of flags
        flags.add(flg);
//TODO remove
//System.out.println("SIZE FLAG-IN2: "+flags.size());

//        freeFlag = flags.indexOf(flg);
	freeFlag = flags.size() - 1;

//TODO remove
//System.out.println("RETURN: "+freeFlag);

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
    public void deleteAtomsFlag(int flagID)
    {
        this.flags.remove(flagID);
    }
//------------------------------------------------------------------------------

    /**
     * @return the molecular weight of the fragment
     */

    public Double getFragMW()
    {
        Double mw = 0.0;
        try {
            MolecularFormulaManipulator mf = new MolecularFormulaManipulator();
            mw = mf.getNaturalExactMass(mf.getMolecularFormula(this));
        } catch (Throwable t) {
            System.out.println("ERROR! Exception while calculating Molecular Weight");
            System.exit(0);
        }
	return mw;
    }
//------------------------------------------------------------------------------

    private void exploreMolecule(int seed, int flag, GM3DAttachmentPoint ap)
    {
        //set string for reporting and debugging
        String recFlag = "";
        for (int ri = 0; ri < recNum; ri++)
             recFlag = recFlag+"-";

        //deal with the starting point of this search
        IAtom seedAtm = this.apOnMol.getAtom(seed);
        addAtomToOrder(seedAtm,ap,flag,recFlag);

        //get the list of neighbours (connected atoms)
        List<IAtom> neighbourAtoms = this.apOnMol.getConnectedAtomsList(seedAtm);

        //through away already done
        List<IAtom> purgedList = new ArrayList<IAtom>();
        for (IAtom connectedAtom : neighbourAtoms)
        {
            int atmidx = this.apOnMol.getAtomNumber(connectedAtom);
            if (!getAtomsFlag(flag, atmidx))
                purgedList.add(connectedAtom);
        }

        //In case nothing alse to do
        if (purgedList.size() == 0)
            return;

        //reorder according to priority
        if (purgedList.size() > 1)
            purgedList = prioritizeAtomList(purgedList, seedAtm);

        for (IAtom connectedAtom : purgedList)
        {
            //identify atom
            int atmidx = this.apOnMol.getAtomNumber(connectedAtom);

            if (repOnScreen >= 3)
                System.out.println(pre+recFlag+"> connected atom: "+ atmidx+" is "+connectedAtom.getSymbol()+" - Branches: "+this.apOnMol.getConnectedAtomsCount(connectedAtom)+" Flag: "+connectedAtom.getFlag(flag));

            // Use flag to avoid giving a bite on your own tail!
            if (getAtomsFlag(flag, atmidx))
                continue;

            // move to the next shell of atoms
            if (this.apOnMol.getConnectedAtomsCount(connectedAtom) > 1)
            {
                recNum++;
                if (repOnScreen >= 3)
                    System.out.println(pre+recFlag+"> recursion on atom "+atmidx+" which is "+connectedAtom.getSymbol());
                exploreMolecule(atmidx,flag,ap);
                if (repOnScreen >= 3)
                    System.out.println(pre+recFlag+"> recursion ALL DONE!");
                recNum--;
            } else {
                connectedAtom.setFlag(flag,true);
                addAtomToOrder(connectedAtom,ap,flag,recFlag);
            }
        }
    }

//------------------------------------------------------------------------------

    private void addAtomToOrder(IAtom atm, GM3DAttachmentPoint ap, int visitedFlag, String recFlag)
    {
        int atmidx = this.apOnMol.getAtomNumber(atm);

        //set visited flag
        setAtomsFlag(visitedFlag,atmidx,true);

        this.atomOrders.get(ap).add(atmidx);
        if (repOnScreen >= 3)
            System.out.println(pre+recFlag+"> Adding to "+ap.getAPClass()+": "+atmidx+" "+atm.getSymbol());        
    }

//------------------------------------------------------------------------------

/**
 * root is the atom connected to all the ones in the inList 
 * (before it was called seed, now it becomes root of the GM3DLigands)
 *
 */
    private List<IAtom> prioritizeAtomList(List<IAtom> inList, IAtom root)
    {

        //set string for reporting and debugging
        String recFlag = "";
        for (int ri = 0; ri < recNum; ri++)
            recFlag = recFlag+"-";

        List<IAtom> outList = new ArrayList<IAtom>();
        List<GM3DLigand> ligands = new ArrayList<GM3DLigand>();
        if (repOnScreen >= 3)
            System.out.println(pre+recFlag+"PRIORITAZING LIST with root: "+this.apOnMol.getAtomNumber(root));

        for (IAtom seed : inList)
        {
            if (repOnScreen >= 3)
                System.out.println(pre+recFlag+"Creating GM3DLigand on atom "+this.apOnMol.getAtomNumber(seed)+seed.getSymbol());
            GM3DLigand lig = new GM3DLigand(seed,root,this.apOnMol);
            ligands.add(lig);
        }

        if (repOnScreen >= 3)
        {
            System.out.println(pre+recFlag+" LIGANDS: ");
            for (int i = 0; i<ligands.size(); i++)
                System.out.println(pre+recFlag+"     "+this.apOnMol.getAtomNumber(ligands.get(i).getSeed()));
        }
        
        Collections.sort(ligands, new GM3DLigandComparator()); 


        if (repOnScreen >= 3)   
        {       
            System.out.println(pre+recFlag+" SORTED LIGANDS: ");
            for (int i = 0; i<ligands.size(); i++)
                System.out.println(pre+recFlag+"  (S) "+this.apOnMol.getAtomNumber(ligands.get(i).getSeed()));
        }

        for (int i = 0; i<ligands.size(); i++)
        {
            IAtom a = ligands.get(i).getSeed();
            outList.add(a);
        }
        
        return outList;
    }

//------------------------------------------------------------------------------
    /**
     * Generated the INCHI code for a this fragments using pseuso atoms
     * as attachment points. 
     */

    public String getInChiForFragmentWithAP(Map<String,String> classToElement)
    {
        //Use 'weird' elements instead of attachment points and dummies
        int apFlag = this.getFreeAtomsFlag();
        int dummyFlag = this.getFreeAtomsFlag();
        for (int i = 0; i<this.apOnMol.getAtomCount(); i++)
        {
            IAtom atm = this.apOnMol.getAtom(i);
            if (atm.getSymbol().equals("AP"))
            {
                setAtomsFlag(apFlag,i,true);
                GM3DAttachmentPoint ap = (GM3DAttachmentPoint) atm.getProperty("CLASS");
                String apClass = ap.getAPClass();
                atm.setSymbol(classToElement.get(apClass));
                if (repOnScreen >= 2)
                    System.out.println(pre+" Setting Symbol of class "+apClass+" to "+classToElement.get(apClass)+" for atom "+i);
            } else if (atm.getSymbol().equals("Du")) {
		setAtomsFlag(dummyFlag,i,true);
		atm.setSymbol("Tc");
                if (repOnScreen >= 2)
                    System.out.println(pre+" Setting Symbol 'Tc' of atom "+i);
	    }
        }

//TODO remove
//IOtools.writeSDFAppend("forInchi.sdf",this.apOnMol,true);

        //Get the InChi for the fragment having AP represented by uncommon elements
        String inchi = "";
	try {
	    //Get non standard INCHI
	    List<net.sf.jniinchi.INCHI_OPTION> opts = new ArrayList<net.sf.jniinchi.INCHI_OPTION>();
	    opts.add(net.sf.jniinchi.INCHI_OPTION.valueOf("RecMet"));
            opts.add(net.sf.jniinchi.INCHI_OPTION.valueOf("DoNotAddH"));
            opts.add(net.sf.jniinchi.INCHI_OPTION.valueOf("SPXYZ"));
            InChIGeneratorFactory inchiFac = InChIGeneratorFactory.getInstance();
            InChIGenerator inchiGen = inchiFac.getInChIGenerator(this.apOnMol);
            INCHI_RET retStat = inchiGen.getReturnStatus();
            inchi = inchiGen.getInchiKey();
	} catch (Throwable t) {
	    inchi = "";
//TODO: remove. This was for debugging
/*
System.out.println("INCHI FAILED: ");
IOtools.writeSDFAppend("checkINCHI.sdf",this.apOnMol,false);
t.printStackTrace();
IOtools.pause();
*/
	}

        if (repOnScreen >= 3)
            System.out.println(pre+"InChI is: "+inchi);

//TODO remove
//System.out.println("SIZE FLAG-4: "+flags.size());

//TODO remove
//System.out.println(pre+"FLAGS AP: "+apFlag);
//System.out.println(pre+"FLAGS Du: "+dummyFlag);

        //Reset back atom symbols
        for (int i = 0; i<this.apOnMol.getAtomCount(); i++)
        {
            IAtom atm = this.apOnMol.getAtom(i);
            if (getAtomsFlag(apFlag,i))
                atm.setSymbol("AP");                
        }
        for (int i = 0; i<this.apOnMol.getAtomCount(); i++)
        {
            IAtom atm = this.apOnMol.getAtom(i);
            if (getAtomsFlag(dummyFlag,i))
                atm.setSymbol("Du");
        }

//TODO remove
//IOtools.writeSDFAppend("forInchi.sdf",this.apOnMol,true);

        //Clean up
        deleteAtomsFlag(dummyFlag);
        deleteAtomsFlag(apFlag);

        return inchi;
    }

//------------------------------------------------------------------------------

    private Map<String,String> getClassToElement(Set<String> classSet, IAtomContainer molA, IAtomContainer molB) 
    {
        Map<String,String> c2e = new HashMap<String,String>();
        Set<String> oneCandEls = new HashSet<String>(Arrays.asList("D","T","He","Ne","Ar","Kr","Xe","Rn","Cs","Fr"));

        for (String apClass : classSet)
        {
            boolean foundOneFree = false;
            for (String el : oneCandEls)
            {
                //avoid to use twice the same element
                if (c2e.containsValue(el))
                    continue;

                //Check if the molecules already contain this element
                boolean useIt = true;
                for (IAtom a : molA.atoms())
                    if (a.getSymbol().equals(el))
                        useIt = false;

                for (IAtom b : molB.atoms())
                    if (b.getSymbol().equals(el))
                        useIt = false;
                
                //OK, now you can use it!
                if (useIt)
                {
                    if (repOnScreen >= 3)
                        System.out.println(pre+"Using element "+el+" for class "+apClass);
                    foundOneFree = true;
                    c2e.put(apClass,el);
                    break;
                }
            }
            //Well, if you have to deal with very uncommon chemistry this may require something more
            if (!foundOneFree)
                throw new Error("ERROR! Free Element not found (generating InChi). Please report this error to the author. He will fix the code for you!");
        }
        return c2e;
    }

//------------------------------------------------------------------------------

    public void resetFlag(IAtomContainer frag, int flag)
    {
        for (IAtom atm : frag.atoms())
        {
            atm.setFlag(flag,false);
        }
    }
//------------------------------------------------------------------------------

    private int getAtomNumberOfAP(IAtomContainer frag, GM3DAttachmentPoint ap)
    {
        int num = -1;
        int atmid = ap.getAPAtm();
        IAtom atm = frag.getAtom(atmid);
        List<IAtom> lst = frag.getConnectedAtomsList(atm);
        for (IAtom a : lst)
        {
            if (a.getSymbol().equals("AP"))
                if (a.getProperty("CLASS").equals(ap))
                   num = frag.getAtomNumber(a);
        }
        return num;
    }

//------------------------------------------------------------------------------
    public void updateAPs()
    {
        for (int atomid = 0; atomid<this.getAtomCount(); atomid++)
        {
            IAtom atm = this.getAtom(atomid);
            try {
                ArrayList<GM3DAttachmentPoint> APs = (ArrayList<GM3DAttachmentPoint>)atm.getProperty("AttachmentPoints");
                for (int i = 0; i < APs.size(); i++)
                {
                    GM3DAttachmentPoint ap = APs.get(i);
                    ap.setAPAtm(atomid);
                }
            } catch (Throwable t ) {
                continue;
            }
        }
    }

//------------------------------------------------------------------------------
    public ArrayList<GM3DAttachmentPoint> getAllAPs()
    {
        return allAPs;
    }

//------------------------------------------------------------------------------
    private void extractAllAPs()
    {
        //extract info from atomic properties
        for (IAtom atm : this.atoms())
        {
            try {
                ArrayList<GM3DAttachmentPoint> APs = (ArrayList<GM3DAttachmentPoint>)atm.getProperty("AttachmentPoints");
                for (int i = 0; i < APs.size(); i++)
                {
                    allAPs.add(APs.get(i));
                }
            } catch (Throwable t ) {
                continue;
            }
        }

        //Reorder according to GM3DAttachmentPoint priority
        Collections.sort(allAPs, new GM3DAttachmentPointComparator());
    }

//------------------------------------------------------------------------------
    private void movePropertyToAP(String format)
    {

	if (format.equals("DENOPTIM"))
	{
            //extract info from molecular properties
	    String allAtomsProp = "";
	    try {
                allAtomsProp = this.getProperty("CLASS").toString();
            } catch (Throwable t ) {
                System.out.println("ERROR! AtomContainer does'n looks like a fragment in '"+format+"' format!");
		System.out.println("Execution killed!"+t);
		System.exit(-1);
            }

	    String[] atomsPRop = allAtomsProp.split(Parameters.moreAtmsSeparator);
	    for (int i = 0; i< atomsPRop.length; i++)
	    {
		String onThisAtm = atomsPRop[i];
		if (onThisAtm.contains(Parameters.moreAPSeparator))
		{
		    String[] moreAPonThisAtm = onThisAtm.split(Parameters.moreAPSeparator);
		    GM3DAttachmentPoint ap = new GM3DAttachmentPoint(moreAPonThisAtm[0],format);
		    int atmID = ap.getAPAtm();
		    //DENOPTIM's format used [1 to n+1] instead of [0 to n]
		    atmID = atmID-1;
		    ap.setAPAtm(atmID);
                    allAPs.add(ap);
		    for (int j = 1; j<moreAPonThisAtm.length; j++ )
		    {
			GM3DAttachmentPoint apMany = new GM3DAttachmentPoint(atmID,moreAPonThisAtm[j],format);
			allAPs.add(apMany);
		    }
		} else {
		    GM3DAttachmentPoint ap = new GM3DAttachmentPoint(onThisAtm,format);
		    int atmID = ap.getAPAtm();
                    //DENOPTIM's format used [1 to n+1] instead of [0 to n]
                    atmID = atmID-1;
		    ap.setAPAtm(atmID);
		    allAPs.add(ap);
		}
	    }

            //Reorder according to GM3DAttachmentPoint priority
            Collections.sort(allAPs, new GM3DAttachmentPointComparator());

	    //Write attachment points also on the atoms
            for (int i = 0; i < allAPs.size(); i++)
            {
                GM3DAttachmentPoint ap = allAPs.get(i);
                int atmID = ap.getAPAtm();
                IAtom atm = this.getAtom(atmID);
                try {
                    ArrayList<GM3DAttachmentPoint> oldAPs = (ArrayList<GM3DAttachmentPoint>)atm.getProperty("AttachmentPoints");
                    oldAPs.add(ap);
                    atm.setProperty("AttachmentPoints",oldAPs);
                } catch (Throwable t ) {
                    ArrayList<GM3DAttachmentPoint> aps = new ArrayList<GM3DAttachmentPoint>();
                    aps.add(ap);
                    atm.setProperty("AttachmentPoints",aps);
                }
	    } 


//TODO add other fragments here
/*
	} else if (format.equals("");
	    .... other format ...
*/
//TODO
	} else {
	    System.out.println("ERROR! Wrong fromat for improting fragments");
	    System.exit(-1);
	}
    }

//-------------------------------------------------------------------

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

//------------------------------------------------------------------------------
    private void makeAPPseudoAtoms()
    {
        for (int i = 0; i < allAPs.size(); i++)
        {
            GM3DAttachmentPoint ap = allAPs.get(i);

            //Create the pseudo-atom
            ArrayList<Double> vec = ap.getAPVector();
            Point3d pt = new Point3d(vec.get(0),vec.get(1),vec.get(2));
// FIRST alternative: independent from Fragmentation
            IAtom atm = new Atom("AP",pt);
            atm.setProperty("CLASS",ap);
// SECOND alternative: integer defining class uniquely using the list of all classe
// to be implemented
//
            apOnMol.addAtom(atm);

            //Define connection between AP and molecule
//MF: debugging the remotion of duplicate fragments containing double bonds to AP
//    The APClass takes already into account the bond order so in 'DENOPTIM' format
//    we do not read the fiels <ATTACMENT_POINT> that contains the bond order.
//    Since the format 'DENOPTIM' does not report the bond order in the same property
//    it's rinsky to read in paralled the two fields assuming that the order is always the same.
//    HERE it's better to use only singol bond for pseudo atoms representing AP

//            Order bndOrd = intToBondOrder(ap.getAPBondOrder());
            Order bndOrd = intToBondOrder(1);

//MF:debug end
            IBond bnd = new Bond(apOnMol.getAtom(ap.getAPAtm()),atm,bndOrd);
            apOnMol.addBond(bnd);
        }
    }

//------------------------------------------------------------------------------
    public void moveAPsTOProperties(String format)
    {

	if (format.equals("DENOPTIM"))
	{
            //WARNING! Here we use enumeration 1-to-n instead of 0-to-(n-1)
            //         to produce a file readable by DENOPTIM

	    //Look for AP on the same atom
	    Map<Integer,Set<Integer>> allAPOnOneAtm = new HashMap<Integer,Set<Integer>>();
	    for (int i = 0; i < allAPs.size(); i++)
            {
		GM3DAttachmentPoint ap = new GM3DAttachmentPoint();
                ap = allAPs.get(i);
		int atmID = ap.getAPAtm();
		
		if (allAPOnOneAtm.keySet().contains(atmID))
		{
		    Set<Integer> oldAP = allAPOnOneAtm.get(atmID);
		    oldAP.add(i);
		    allAPOnOneAtm.put(atmID,oldAP);
		} else {
		    Set<Integer> newAP = new HashSet<Integer>();
		    newAP.add(i);
		    allAPOnOneAtm.put(atmID,newAP);
		}
	    }
            
	    //Translate information from GM3DAttachmentPoint into Strings
            String apClass = "";
            String apBond = "";

	    for (int atmID : allAPOnOneAtm.keySet())
	    {
		boolean firstCL = true;
		for (int APref : allAPOnOneAtm.get(atmID))
		{
		    GM3DAttachmentPoint ap = new GM3DAttachmentPoint();
		    ap = allAPs.get(APref);
		    String stingAPP = ""; //String Attachment Point Property

		    //Build SDF property "CLASS"
		    if (firstCL)
		    {
			firstCL = false;
			stingAPP = ap.toString(format);
			if (apClass.equals(""))
	                    stingAPP = stingAPP.substring(1,stingAPP.length());
		    } else {
			stingAPP = ap.toStringFurther(format);
		    }
		    apClass = apClass + stingAPP;

		    //Build SDF property "ATTACHMENT_POINT"
		    int BndOrd = ap.getAPBondOrder();
		    String sBO = Integer.toString(BndOrd);
		    String stBnd = " " + Integer.toString(atmID+1)+":"+sBO;
		    if (apBond.equals(""))
                        stBnd = stBnd.substring(1,stBnd.length());
		    apBond = apBond + stBnd;
		}
	    }

            //Set properties
            this.setProperty("CLASS",apClass);
            this.setProperty("ATTACHMENT_POINT",apBond);

//TODO add other formats here
/*
        } else if (format.equals("");
            .... other format ...
*/
//TODO

        } else {
            System.out.println("ERROR! Not recognized format requested for storing fragments");
            System.exit(-1);
        }
    }

//-----------------------------------------------------------------------------

    /**
     * Determines the dimensionality of the chemical object submitted
     * @param mol input molecular/chemical object
     * @return dimensionality of this object (2 or 3) or -1
     */
    private int getDimensions(IAtomContainer mol)
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

//------------------------------------------------------------------------------

    public IAtomContainer toIAtomContainer(String format)
    {
        IAtomContainer mol = new AtomContainer();
        moveAPsTOProperties(format);
        mol = this;

        return mol;
    }

//------------------------------------------------------------------------------

    public IAtomContainer toIAtomContainer()
    {
        IAtomContainer mol = new AtomContainer();
        mol = this;

        return mol;
    }

//------------------------------------------------------------------------------
    /**
     * Store molecules that will not be treated becouse of troubles
     * during execution of task.
     * @param mol input molecular/chemical object
     * @param reason brief explamention of the reason for the rejection
     */
    private void rejectMol(IAtomContainer mol, String reason)
    {
        mol.setProperty("REJECTED",reason);
        IOtools.writeSDFAppend(checkfile,mol,true);
    }

//-----------------------------------------------------------------------------
}
