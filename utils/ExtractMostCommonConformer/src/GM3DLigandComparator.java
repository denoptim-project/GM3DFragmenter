import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

import java.util.Comparator;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;


/**
 *
 * @author Marco Foscato (University of Bergen) Discussions: Giovanni Occhipinti and Vidar R. Jensen (University of Bergen)
 */

public class GM3DLigandComparator implements Comparator<GM3DLigand>
{

    //Flag for reporting
    private int repOnScreen = 0;
    private String pre = "GM3DLigandComparator: ";

    @Override
    public int compare(GM3DLigand a, GM3DLigand b)
    {
	final int FIRST = 1;
	final int EQUAL = 0;
        final int LAST = -1;
	String frs = "First GM3DLigand has priority! ";
	String sec = "Second GM3DLigand has priority! ";

	Map<Integer,List<IAtom>> listA = new HashMap<Integer,List<IAtom>>();
	listA = a.getOrderedList();
        Map<Integer,List<IAtom>> listB = new HashMap<Integer,List<IAtom>>();
        listB = b.getOrderedList();

        //compare number of levels
        int levelsA = listA.keySet().size();
        int levelsB = listB.keySet().size();

        if (repOnScreen >= 2)
	    System.out.println(pre+"Compare list on atoms of GM3DLigands having seed "+
				a.getMol().getAtomNumber(a.getSeed())+a.getSeed().getSymbol()+
				" and "+
				b.getMol().getAtomNumber(b.getSeed())+b.getSeed().getSymbol()+":");

	int min = 0;
	boolean useA = true;
	if (levelsA > levelsB)
	{
	    useA = false;
	    min = levelsB;
	}
	else
	    min = levelsA;

        for (int lev = 0; lev<min; lev++)
        {
	    if (repOnScreen >= 2)
	        System.out.println(pre+"Entries level "+lev);

	    //identify maximun number of atoms to compare (as the minimum between sizes)
	    int minAtms = 0;
	    if (listA.get(lev).size() > listB.get(lev).size())
		minAtms = listB.get(lev).size();
	    else
                minAtms = listA.get(lev).size();

	    for (int iatm = 0; iatm<minAtms; iatm++)
            {
	        //Get the elements to compate
		IAtom atmA = listA.get(lev).get(iatm);
                IAtom atmB = listB.get(lev).get(iatm);

                if (repOnScreen >= 2)
                    System.out.println(pre+" -> \t"+a.getMol().getAtomNumber(listA.get(lev).get(iatm))+atmA.getSymbol() + " \t "+ b.getMol().getAtomNumber(listB.get(lev).get(iatm)) + atmB.getSymbol());

		//Compare atoms
		GM3DAtomComparator comAtms = new GM3DAtomComparator();
		int res = comAtms.compareAtoms(atmA,a.getMol(),atmB,b.getMol());
		if (res > 0)
		    return FIRST;
		else if (res < 0)
		    return LAST;
	    }
	}

        if (repOnScreen >= 2)
            System.out.println(pre+" Ligands are EQUAL!");

	return EQUAL;
    }
}
