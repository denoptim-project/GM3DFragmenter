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
 * SMARTS string with capacity of getting simpler by increasing removal of
 * structural details. In practice, atoms' definition are moved to 'any atom'
 * at each simplification of the SMARTS string.
 *
 * @author Marco Foscato (University of Bergen) Discussions: Giovanni Occhipinti and Vidar R. Jensen (University of Bergen)
 */

public class GM3DSMARTS
{ 
    //String
    private String smarts;
    //Simplified SMARTS
    private String simplerSmarts;
    //Current simplification level
    private int simpLev = 0;
    //Maximum number of semplifications
    private int maxSemp = 1;
    //List of atoms suitable for simplification
    private ArrayList<List<Integer>> extremesOfAtomsToSimplify = new ArrayList<List<Integer>>();


//------------------------------------------------------------------------------

    /**
     * Construct an empty object <code>GM3DSMARTS</code>
     */
    public GM3DSMARTS()
    {
    }

//------------------------------------------------------------------------------

    /**
     * Construct a <code>GM3DSMARTS</code> defining the SMARTS string
     */
    public GM3DSMARTS(String smarts)
    {
        //Set the string
        this.smarts = smarts;
        //identify possible simplifications
        this.findAtomsInSquares();
    }
//------------------------------------------------------------------------------

    public String toString()
    {
        return smarts;
    }

//-----------------------------------------------------------------------------

    private void findAtomsInSquares()
    {
        String oSP = "[";
        String cSP = "]";
        int size = smarts.length();
        int init = 0;
        boolean goon = true;
        ArrayList<Integer> nonMostInternal = new ArrayList<Integer>();
        while (goon)
        {
            if (!smarts.substring(init).contains(oSP))
                break;

            int a = smarts.indexOf(oSP,init);
            int c = smarts.indexOf(cSP,a);

            int next = smarts.indexOf(oSP, (a + 1) );

            if ((next > c) || (next == -1))
            {
                List<Integer> pair = new ArrayList<Integer>();
                pair.add(a);
                pair.add(c);
                extremesOfAtomsToSimplify.add(pair);
                init = c;
            } else {
                init = a +1;
            }
        }
	this.maxSemp = extremesOfAtomsToSimplify.size();
    }

//-----------------------------------------------------------------------------

    public int getMaxNumSemplification()
    {
	return maxSemp;
    }

//-----------------------------------------------------------------------------

    /**
     * Returns a simpler SMARTS meaning that the last atom non equal to [*]
     * is changed into [*] and the new SMARTS is returned. Everytime this
     * method runs the simplified version of the smile gets simpler.
     * <strong>WARNING! Once the method is called older simplifications of this 
     * SMARTS cannot be recovered!</strong>
     */
    public String getSimplerSMARTS()
    {
	updateSimpler();
        return simplerSmarts;
    }

//-----------------------------------------------------------------------------

    private void updateSimpler()
    {
        if (simpLev < extremesOfAtomsToSimplify.size())
            simpLev++;
        String newSmarts = "";
        newSmarts = smarts.substring(0,extremesOfAtomsToSimplify.get(0).get(0));
        for (int j=0; j<extremesOfAtomsToSimplify.size(); j++)
        {
            int init = extremesOfAtomsToSimplify.get(j).get(0);
            int end = extremesOfAtomsToSimplify.get(j).get(1);
//System.out.println("IN :"+init+" "+end);
            int whereToChange = extremesOfAtomsToSimplify.size() - simpLev;
            if (j < whereToChange)
            {
                newSmarts = newSmarts + smarts.substring(init,end + 1);
            } else {
                newSmarts = newSmarts + "[*]";
            }

            if (end < smarts.length())
            {
                if ((j+1) < extremesOfAtomsToSimplify.size())
                {
                    int i2 = end + 1;
                    int e2 = extremesOfAtomsToSimplify.get(j+1).get(0);
//System.out.println("IN :"+i2+" "+e2);
                    newSmarts = newSmarts + smarts.substring(i2,e2);
                } else {
                    newSmarts = newSmarts + smarts.substring(end + 1);
                }
            }
//System.out.println("IN newSmarts = "+newSmarts);
        }
	this.simplerSmarts = newSmarts;
    }
//------------------------------------------------------------------------------
}
