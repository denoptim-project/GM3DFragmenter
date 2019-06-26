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

/**
 * Class of an attachment point: it represents the cutting rule by which 
 * action the attachment point is created
 *
 * @author Marco Foscato (University of Bergen) Discussions: Giovanni Occhipinti and Vidar R. Jensen (University of Bergen)
 */

public class GM3DAPClass
{
    // String names
    private String apClass;
    private String ruleName;
    private String subClassID;

    // Bond order
    private int bndOrd;

    // Rule involves metal
    private boolean metal;

    // MultiHapto
    private boolean multiHapto;

    // For symmetric cutting rules complementary AP class does not exist
    private boolean symmetric;

    // Separator
    private String sep = ":";


//------------------------------------------------------------------------------

    public GM3DAPClass()
    {
	this.ruleName = "";
        this.subClassID = "";
        this.apClass = "";
        this.bndOrd = -1;
        this.metal = false;
        this.multiHapto = false;
        this.symmetric = false;
    }

//------------------------------------------------------------------------------

    public GM3DAPClass(String ruleName, String subClassID)
    {
        this.ruleName = ruleName;
        this.subClassID = subClassID;
        this.apClass = ruleName+sep+subClassID;
        this.bndOrd = -1;
        this.metal = false;
        this.multiHapto = false;
	this.symmetric = false;
    }

//------------------------------------------------------------------------------

    public GM3DAPClass(String ruleName, String subClassID, int bndOrd, boolean metal, boolean multiHapto, boolean symmetric)
    {
        this.ruleName = ruleName;
        this.subClassID = subClassID;
	this.apClass = ruleName+sep+subClassID;
        this.bndOrd = bndOrd;
        this.metal = metal;
        this.multiHapto = multiHapto;
	this.symmetric = symmetric;
    }

//------------------------------------------------------------------------------

/**
 * Returns AP class [cutting rule][separator][subclass index]
 */ 
    public String getAPClass()
    {
	return apClass;
    }

//------------------------------------------------------------------------------

/**
 * Returns the name of the cutting rule
 */
    public String getRuleName()
    {
        return ruleName;
    }

//------------------------------------------------------------------------------

/**
 * Returns the subclass index
 */
    public String getSubClass()
    {
        return subClassID;
    }

//------------------------------------------------------------------------------

/**
 * @return <code>true</code> if the AP class corresponds to a multihapto ligand
 */
    public boolean isMultiHapto()
    {
        return multiHapto;
    }

//------------------------------------------------------------------------------

/**
 * @return the bond order as an integer or -1 if not set
 */
    public int getBondOrder()
    {
        return bndOrd;
    }

//------------------------------------------------------------------------------

/**
 * @return <code>true</code> if the AP class derives from a cutting rule involving metals
 */
    public boolean involvesMetal()
    {
        return metal;
    }

//------------------------------------------------------------------------------

/**
 * @return <code>true</code> if AP class doesn't have complementary class
 */
    public boolean isSymmetric()
    {
        return symmetric;
    }

//------------------------------------------------------------------------------

/**
 * Returns the string representing this rule
 */
    public String toString()
    {
	String str = apClass+"_BO:"+bndOrd+"_met:"+metal+"_hpt:"+multiHapto+"_symm:"+symmetric;
        return str;
    }

//------------------------------------------------------------------------------

}
