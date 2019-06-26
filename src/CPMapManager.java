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

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;

import java.util.Collections;

/**
 * CPMapManager is a list of tools for dealing with the 
 * compatibility matrix
 *
 * @author Marco Foscato (University of Bergen)
 */

class CPMapManager
{

    //Output name
    private String compMatFile;

    //Report on screen
    private int repOnScreen;
    private String preStr = "CPMapManager";

//-------------------------------------------------------------------

    /**
     * Creates a <code>CPMapManager</code> setting variables 
     * according 
     * to <code>Parameters</code>
     */
    public CPMapManager()
    {
	compMatFile = "CPMap_"+Parameters.getJobName()+".par";
        repOnScreen = Parameters.report;
    }

//-------------------------------------------------------------------

    /**
     * Creates the compatibility matrix on the basis of the list of
     * cutting rules
     */

    public void makeFromCuttingRules()
    {
	//Get cutting rules
	Map<String,GM3DCuttingRule> cutRules = Parameters.rules;

	//Make list of AP classes
        List<GM3DAPClass> apClasses = new ArrayList<GM3DAPClass>();
	for (String ruleName : cutRules.keySet())
	{
	    GM3DCuttingRule rule = cutRules.get(ruleName);
	    if (rule.isSymmetric())
	    {
	        apClasses.add(rule.getAPClass0());
	    } else {
                apClasses.add(rule.getAPClass0());
                apClasses.add(rule.getAPClass1());
	    }
	}
/*
for (GM3DAPClass apc : apClasses)
    System.out.println("PRE: "+
                        " M:"+apc.involvesMetal()+
                        " BO:"+apc.getBondOrder()+
                        " Hp:"+apc.isMultiHapto()+
                        " Sym:"+apc.isSymmetric()+" "+apc.getAPClass());
*/

	//Reorganize list of cutting rules
	GM3DAPClassComparator APClassComp = new GM3DAPClassComparator();
	Collections.sort(apClasses,APClassComp);
	Collections.reverse(apClasses);

/*
for (GM3DAPClass apc : apClasses)
    System.out.println("POST: "+
                        " M:"+apc.involvesMetal()+
                        " BO:"+apc.getBondOrder()+
                        " Hp:"+apc.isMultiHapto()+
                        " Sym:"+apc.isSymmetric()+" "+apc.getAPClass());
*/

	// create the matrix with this properties:
	// - only complementary AP classes are set to compatible
	//
	CompatibilityMatrix cpm = new CompatibilityMatrix();
	for (GM3DAPClass apc : apClasses)
	{
	    GM3DCuttingRule rule = cutRules.get(apc.getRuleName());
//TODO tome nodified after making CompatibilityMatrix.java wirking with APClass objects
	    cpm.addTrueEntry(apc.getAPClass(),rule.getComplementaryAPClass(apc).getAPClass());
	}


        // class-to-bond order section
	Map<String,Integer> classBndOrd = new HashMap<String,Integer>();
        for (GM3DAPClass apc : apClasses)
        {
	    classBndOrd.put(apc.getRuleName(),apc.getBondOrder());
	    if (apc.getBondOrder() < 0)
		System.err.println(preStr+" WARNING! Check bond order of: "+apc.getRuleName());
        }

	// capping
        Map<String,String> cap = new HashMap<String,String>();
//TODO get a proper capping group on the dasis of APClass
	String capAPClass = "hyd:1";
        for (GM3DAPClass apc : apClasses)
        {
            cap.put(apc.getAPClass(),capAPClass);
        }

	// write CPMap on file
	cpm.writeCPMapFile(compMatFile,classBndOrd,cap);	

    }

//------------------------------------------------------------------------------

    /**
     * Writes the compatibility matrix on a text file
     * @param filename name of the output file
     */
/*
    public static void writeCPMapFile(CompatibilityMatrix mp, String filename)
    {
        String date = Parameters.getDate();

        // Print head of the file
        String head = "#\n# CompatibilityMatrix for Class Based Builders\n#";
        head = head + "\n# Created by GM3DFragmenter - [d/m/y] " + date;
        head = head + "\n# Format";
        head = head + "\n# parentClass compatibleClass1:subClass ";
        head = head + "compatibleClass2:subClass [etc.]\n#";

        IOtools.writeTXTAppend(head,filename,true);

        // Print compatibility matrix
        StringBuffer allCl = new StringBuffer();
        allCl.append("# List of all classes: ");
        Set<String> keys = mp.keySet();
        List<String> keylist = new ArrayList<String>(keys);
        Collections.sort(keylist);
        boolean one = true;
        for (String parentclass : keylist)
        {
            StringBuffer line = new StringBuffer();
            line.append(mpKey+" ");
            line.append(parentclass);
            if (one)
            {
                allCl.append(parentclass);
                one = false;
            } else {
                allCl.append(" "+parentclass);
            }

            for (String childclass : mp.get(parentclass))
                line.append(" "+childclass);

            IOtools.writeTXTAppend(line.toString(),filename,true);
        }
        IOtools.writeTXTAppend(allCl.toString(),filename,true);
    }
*/
//-------------------------------------------------------------------
}
