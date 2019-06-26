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

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;


/**
 * Container for list of atoms matching a list of SMARTS
 *
 * @author Marco Foscato (University of Bergen)
 */

class ManySMARTSQuery
{
    //Container
    private Map<String,List<List<Integer>>> allMatches = new HashMap<String,List<List<Integer>>>();
    //Counts
    private int totNum;
    private Map<String,Integer> numMatches = new HashMap<String,Integer>();

    //Problems
    private boolean problems = false;
    private String message = "";

    //Level of information printed on screen
    private int repOnScreen;

//------------------------------------------------------------------------------

    public ManySMARTSQuery()
    {
	super();
	repOnScreen = Parameters.report;
    }

//------------------------------------------------------------------------------

    public ManySMARTSQuery(IAtomContainer mol, Map<String,String> smarts)
    {
        super();
        repOnScreen = Parameters.report;
	totNum = 0;
	String blankSmarts = "[*]";

	String err="";

	try {
                SMARTSQueryTool query = new SMARTSQueryTool(blankSmarts);
		query.setAllRingsFinderTimeout(Parameters.maxTimeAllRingFinder);
                query.setMaxRingSize(Parameters.maxRingSizeMF);
		for (String smartsRef : smarts.keySet())
		{
		    //get the new query
                    String oneSmarts = smarts.get(smartsRef);
		    err = smartsRef;

                    if (repOnScreen >= 3)
                    {
                        System.out.println("Attempt to match query "+smartsRef);
                        System.out.println("SMARTS: "+oneSmarts);
                    }

                    //Update the query tool
		    query.setSmarts(oneSmarts);

		    
		    if (query.matches(mol))
		    {
			//Store matches
			List<List<Integer>> listOfIds = new ArrayList<List<Integer>>();
			listOfIds = query.getUniqueMatchingAtoms();
			allMatches.put(smartsRef,listOfIds);
			//Store number
//CDK BUG here! this number is somehow wrong
//			int num = query.countMatches();
			int num = listOfIds.size();
			numMatches.put(smartsRef,num);
			totNum = totNum + num;
			if (repOnScreen >= 2)
			{
			    System.out.println("Matches for query '"+smartsRef+"': "+num+" => Atoms: "+listOfIds);
			}
 		    }
		}
        } catch (CDKException cdkEx) {
                String cause = cdkEx.getCause().getMessage();
		err = "\nWARNING! For query " + err + " => " + cause;
//TODO del
//cdkEx.printStackTrace();
		problems = true;
		message = err;
	} catch (Throwable t) {
		java.lang.StackTraceElement[] stes = t.getStackTrace();
		String cause = "";
		int s = stes.length;
		if (s >= 1)
		{
		    java.lang.StackTraceElement ste = stes[0];
		    cause = ste.getClassName();
		} else {
		    cause = "'unknown' (try to process this molecule alone to get more infos)";
		}
                err = "\nWARNING! For query " + err + " => Exception returned by "+cause;
                problems = true;
                message = err;
//TODO del
//t.printStackTrace();
	}
    }

//------------------------------------------------------------------------------

    public boolean hasProblems()
    {
        return problems;
    }

//------------------------------------------------------------------------------

    public String getMessage()
    {
	return message;
    }

//------------------------------------------------------------------------------

    public int getTotalMatches()
    {
	return totNum;
    }

//------------------------------------------------------------------------------

    public Map<String,Integer> getNumMatchesMap()
    {
        return numMatches;
    }

//------------------------------------------------------------------------------

    public int getNumMatchesOfQuery(String query)
    {
	if (numMatches.keySet().contains(query))
	    return numMatches.get(query);
	else
	    return 0;
    }

//------------------------------------------------------------------------------

    public boolean hasMatches(String query)
    {
        if (numMatches.keySet().contains(query))
            return true;
        else
            return false;
    }

//------------------------------------------------------------------------------

    public List<List<Integer>> getMatchesOfSMARTS(String ref)
    {
        return allMatches.get(ref);
    }

//------------------------------------------------------------------------------
}
