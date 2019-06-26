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

import java.util.Comparator;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

/**
 * Comparator for GMFragments based on Molecular Weight
 *
 * @author Marco Foscato (University of Bergen) Discussions: Giovanni Occhipinti and Vidar R. Jensen (University of Bergen)
 */

public class MWComparator implements Comparator<IAtomContainer>
{
//-------------------------------------------------------------------
    @Override
    public int compare(IAtomContainer a, IAtomContainer b)
    {
        Double mwa = 0.0;
        Double mwb = 0.0;

        try {
            MolecularFormulaManipulator mf = new MolecularFormulaManipulator();
            mwa = mf.getNaturalExactMass(mf.getMolecularFormula(a));
            mwb = mf.getNaturalExactMass(mf.getMolecularFormula(b));
        } catch (Throwable t) {
            System.out.println("ERROR! Exception while calculating Molecular Weight");
            System.exit(0);
        }

        return mwa.compareTo(mwb);
    }
//-------------------------------------------------------------------
}


