import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;


/**
 *
 * @author Marco Foscato (University of Bergen) Discussions: Giovanni Occhipinti and Vidar R. Jensen (University of Bergen)
 */

public class GM3DAtomComparator
{
    //Flag for reporting
    private int repOnScreen = 0;
    private String pre = "GM3DAtomComparator: ";
    private String frs = "First atom has priority! ";
    private String sec = "Second atom has priority! ";

//------------------------------------------------------------------------------

    public GM3DAtomComparator()
    {
    }

//------------------------------------------------------------------------------

    public int compareAtoms(IAtom atmA, IAtomContainer molA, IAtom atmB, IAtomContainer molB)
    {
        final int FIRST = 1;
        final int EQUAL = 0;
        final int LAST = -1;

        int atNumA;
        if (atmA.getSymbol().equals("AP") || atmA.getSymbol().equals("Du"))
            atNumA = 0;
        else
            atNumA = atmA.getMassNumber();

        int atNumB;
        if (atmB.getSymbol().equals("AP") || atmB.getSymbol().equals("Du"))
            atNumB = 0;
        else
            atNumB = atmB.getMassNumber();

        int numConAtmA = molA.getConnectedAtomsCount(atmA);
        int numConAtmB = molB.getConnectedAtomsCount(atmB);

        if (atNumA > atNumB)
        {
            if (repOnScreen >= 2)
                System.out.println(pre+frs+"Due to Mass number");
            return FIRST;
        } else if (atNumA < atNumB)
        {
            if (repOnScreen >= 2)
                System.out.println(pre+sec+"Due to Mass number");
            return LAST;
        }

        //Until here we have atoms with the same number of mass
        //Check number of connections
        if (numConAtmA > numConAtmB)
        {
            if (repOnScreen >= 2)
                System.out.println(pre+frs+"Due to Number of Connections");
            return FIRST;
        } else if (numConAtmA < numConAtmB)
        {
            if (repOnScreen >= 2)
                System.out.println(pre+sec+"Due to Number of Connections");
            return LAST;
        }

//TODO
                //Here we have the same atoms, with the same number of connections
                //add other criteria here!
//TODO
	
	return EQUAL;
    }

//------------------------------------------------------------------------------
}
