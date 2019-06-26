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

/**
 * Comparator for attachment point classes
 * 
 * @author Marco Foscato (University of Bergen) Discussions: Giovanni Occhipinti and Vidar R. Jensen (University of Bergen)
 */

public class GM3DAPClassComparator implements Comparator<GM3DAPClass>
{
   @Override
   public int compare(GM3DAPClass apc1, GM3DAPClass apc2)
    {
	Boolean hasMet1 = apc1.involvesMetal();
        Boolean hasMet2 = apc2.involvesMetal();
	int res = hasMet1.compareTo(hasMet2);
	if (res != 0)
	{
	    return res;
	}

        Boolean hpt1 = apc1.isMultiHapto();
        Boolean hpt2 = apc2.isMultiHapto();
        res = hpt1.compareTo(hpt2);
        if (res != 0)
        {
            return res;
        }

        Boolean s1 = apc1.isSymmetric();
        Boolean s2 = apc2.isSymmetric();
        res = s1.compareTo(s2);
        if (res != 0)
        {
            return res;
        }

	//WARNING! there is a minus here!
	Integer bo1 = apc1.getBondOrder();
        Integer bo2 = apc2.getBondOrder();
        res = -bo1.compareTo(bo2);
        if (res != 0)
        {
            return res;
        }

	String n1 = apc1.getAPClass();
	String n2 = apc2.getAPClass();
	res = n1.compareTo(n2);

        return res;
    }

}

