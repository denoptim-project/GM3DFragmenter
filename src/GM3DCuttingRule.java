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

import java.util.ArrayList;

/**
 * A cutting rule with three SMARTS queries (atom 1, bond, atom2) and options
 *
 * @author Marco Foscato (University of Bergen) Discussions: Giovanni Occhipinti and Vidar R. Jensen (University of Bergen)
 */

public class GM3DCuttingRule
{
    // Rule name
    private String ruleName;
    private String subClass0 = "0";
    private String subClass1 = "1";

    // AP classes deriving from this rule
    private GM3DAPClass apc0;
    private GM3DAPClass apc1; //null for symmetric rules

    // SMARTS queries 
    private String smartsAtm0;
    private String smartsAtm1;
    private String smartsBnd;

    // Priority index
    private int priority;

    // Options
    private ArrayList<String> opts;


//------------------------------------------------------------------------------

    public GM3DCuttingRule(String ruleName, String smartsAtm0, String smartsAtm1, String smartsBnd, int priority, ArrayList<String> opts)
    {
        this.ruleName = ruleName;
        this.smartsAtm0 = smartsAtm0;
        this.smartsAtm1 = smartsAtm1;
        this.smartsBnd = smartsBnd;
        this.priority = priority;
        this.opts = opts;
	this.apc0 = new GM3DAPClass(ruleName,subClass0,this.getBondOrder(),this.involvesMetal(),this.isHAPTO(),this.isSymmetric());
	if (!this.isSymmetric())
	{
	    this.apc1 = new GM3DAPClass(ruleName,subClass1,this.getBondOrder(),this.involvesMetal(),this.isHAPTO(),this.isSymmetric());
	} else {
	    apc1 = null;
	}
    }

//------------------------------------------------------------------------------

/**
 * Returns the name of the cutting rule
 */
    public String getName()
    {
        return ruleName;
    }

//------------------------------------------------------------------------------

/**
 * Returns the name of SubClass0
 */
    public String getSubClassName0()
    {
	return apc0.getAPClass();
    }

//------------------------------------------------------------------------------

/**
 * Returns the name of SubClass1
 */
    public String getSubClassName1()
    {
        if (this.isSymmetric())
            return apc0.getAPClass();
	else
	    return apc1.getAPClass();
    }

//------------------------------------------------------------------------------

/**
 * Get the AP class with sub class 0
 */
    public GM3DAPClass getAPClass0()
    {
        return apc0;
    }

//------------------------------------------------------------------------------

/**
 * Get the AP class with sub class 1
 */
    public GM3DAPClass getAPClass1()
    {
	if (this.isSymmetric())
	    return apc0;
	else
            return apc1;
    }

//------------------------------------------------------------------------------

/**
 * Get complementary class
 */
    public GM3DAPClass getComplementaryAPClass(GM3DAPClass apc)
    {
	if (this.isSymmetric() && apc.equals(apc0))
	{
	    return apc0;
	} else if (!this.isSymmetric() && apc.equals(apc0)) {
	    return apc1;
        } else if (!this.isSymmetric() && apc.equals(apc1)) {
            return apc0;
	} else {
	    System.err.println("ERROR! Attempt to use SubClass1 from a symmetric cutting rule!");
	    return null;
	}
    }

//------------------------------------------------------------------------------

/**
 * Returns the SMARTS query of the whole rule
 */
    public String getWholeSMARTSRule()
    {
        return smartsAtm0+smartsBnd+smartsAtm1;
    }

//------------------------------------------------------------------------------

/**
 * Get the SMARTS query of the first atom (SubClass 0)
 */
    public String getSMARTSSubClass0()
    {
	return smartsAtm0;
    }

//------------------------------------------------------------------------------

/**
 * Get the SMARTS query of the second atom (SubClass 1)
 */
    public String getSMARTSSubClass1()
    {
        return smartsAtm1;
    }

//------------------------------------------------------------------------------

/**
 * Get the SMARTS query of the bond
 */
    public String getSMARTSBnd()
    {
        return smartsBnd;
    }

//------------------------------------------------------------------------------

/**
 * @return true if this rule has further options
 */
    public boolean hasOptions()
    {
	if (opts.size() > 0)
            return true;
	else
	    return false;
    }

//------------------------------------------------------------------------------

/**
 * @return true if this rule mathcer multihapto ligands
 */
    public boolean isHAPTO()
    {
        if (opts.contains("HAPTO"))
            return true;
        else
            return false;
    }

//------------------------------------------------------------------------------

/**
 * @return true if this is a symmetric rule (atom SMARTS coincide)
 */
    public boolean isSymmetric()
    {
        if (smartsAtm0.equals(smartsAtm1))
            return true;
        else
            return false;
    }

//------------------------------------------------------------------------------

/**
 * Returns the list of options
 */
    public ArrayList<String> getOptions()
    {
        return opts;
    }


//------------------------------------------------------------------------------

/**
 * Identify the bond order of the matched bond
 */
    public int getBondOrder()
    {
	int res = -1;

	//Easy case
	String s = smartsBnd;
	if (s.contains("!@"))
	{
	    s = s.substring(0,s.indexOf("!@"));
	}
	if (s.equals("-"))
	    res = 1;
	else if (s.equals("="))
	    res = 2;
	else if (s.equals("#"))
            res = 3;

        return res;

	// OK, it's not so easy...
//TODO make it more general
	    
    }

//------------------------------------------------------------------------------

/**
 * @return <code>true</code> if the rule involves metals
 */
    public boolean involvesMetal()
    {
	ArrayList<String> metals = Parameters.metals;
	//analyze atom 0
	for (String el : metals)
	{
            if (smartsAtm0.contains(el))
	    {
		return true;
	    }
	}

        //analyze atom 1
        for (String el : metals)
        {
            if (smartsAtm1.contains(el))
            {
                return true;
            }
        }

	return false;
    }

//------------------------------------------------------------------------------

/**
 * Returns the string representing this rule
 */
    public String toString()
    {
	String str = ruleName+"_"+smartsAtm0+smartsBnd+smartsAtm1+
				"_priority:"+priority+
				"_opts:"+opts;
        return str;
    }

//------------------------------------------------------------------------------

}
