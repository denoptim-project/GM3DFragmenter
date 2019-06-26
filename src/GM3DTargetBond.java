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

import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * The bond matching a cutting rule
 *
 * @author Marco Foscato (University of Bergen) Discussions: Giovanni Occhipinti and Vidar R. Jensen (University of Bergen)
 */

public class GM3DTargetBond
{
    //Atoms 
    private IAtom atm0;  //atom matching subclass 0
    private IAtom atm1;  //atom matching subclass 1

    //Indeces of atoms forming the bond
    private int idAtm0;  //atom matching subclass 0
    private int idAtm1;  //atom matching subclass 1

    //Class and Subclass strings
    private String classSubClass0;
    private String classSubClass1;

    //Flag for Symmetric subclass
    private boolean isSymm;

    //Cutting rule
    private String cutRul;


//------------------------------------------------------------------------------

    public GM3DTargetBond()
    {
        this.idAtm0 = -1;
        this.idAtm1 = -1;
	this.cutRul = "";
	this.isSymm = false;
    }

//------------------------------------------------------------------------------
/*
    public GM3DTargetBond(int idAtm0, int idAtm1, String cutRul, boolean isSymm)
    {
        this.idAtm0 = idAtm0;
        this.idAtm1 = idAtm1;
        this.isSymm = isSymm;

        this.cutRul = cutRul;
    }
*/
//------------------------------------------------------------------------------
/*
    public GM3DTargetBond(int idAtm0, int idAtm1, String cutRul, IAtomContainer mol)
    {
//TODO 
System.out.println("You are trying to use a part of the code that is under development");
System.out.println("To avoiud misleading results, I kill the execution")
System.exit(0);

	this.idAtm0 = idAtm0;
        this.idAtm1 = idAtm1;

	this.atm0 = mol.getAtom(idAtm0);
	this.atm1 = mol.getAtom(idAtm1);

	this.bnd = mol.getBond(atm0,atm1);

	this.cutRul = cutRul;
    }
*/
//------------------------------------------------------------------------------

/**
 * Returns the name of the cutting rule this bond matches
 */
    public String getRule()
    {
        return cutRul;
    }

//------------------------------------------------------------------------------

/**
 * Returns the atom number of the atom matching subclass '1'
 */
    public int getIDSubClass1()
    {
	return idAtm1;
    }

//------------------------------------------------------------------------------

/**
 * Returns the atom number of the atom matching subclass '0'
 */
    public int getIDSubClass0()
    {
        return idAtm0;
    }

//------------------------------------------------------------------------------

/**
 * Returns the atom mathcing subclass '1'
 */
    public IAtom getAtmSubClass1()
    {
        return atm1;
    }

//------------------------------------------------------------------------------

/**
 * Returns the atom matching subclass '0'
 */
    public IAtom getAtmSubClass0()
    {
        return atm0;
    }

//------------------------------------------------------------------------------

/**
 * Returns the string encoding class and subclass-0 membership
 */
    public String getClassSubClass0()
    {
        return classSubClass0;
    }

//------------------------------------------------------------------------------

/**
 * Returns the string encoding class and subclass-1 membership
 */
    public String getClassSubClass1()
    {
        return classSubClass1;
    }

//------------------------------------------------------------------------------

/**
 * @return true if this mond had the same subclass on both atoms
 */
    public boolean hasSymmetricSubClass()
    {
        return isSymm;
    }

//------------------------------------------------------------------------------

/**
 * Set the atom matching subclass '0'
 */
    public void setAtmSubClass0(IAtomContainer mol, int idAtm0)
    {
        this.idAtm0 = idAtm0;
        this.atm0 = mol.getAtom(idAtm0);
    }

//------------------------------------------------------------------------------

/**
 * Set the atom matching subclass '1'
 */
    public void setAtmSubClass1(IAtomContainer mol, int idAtm1)
    {
        this.idAtm1 = idAtm1;
        this.atm1 = mol.getAtom(idAtm1);
    }

//------------------------------------------------------------------------------

/**
 * Set the string encoding class and subclass-0 membership
 */
    public void setSubClass0(String str)
    {
        this.classSubClass0 = str;
    }

//------------------------------------------------------------------------------

/**
 * Set the string encoding class and subclass-1 membership
 */
    public void setSubClass1(String str)
    {
        this.classSubClass1 = str;
    }
       
//------------------------------------------------------------------------------

/**
 * Set the name of the cutting rule that is matched by this bond
 */
    public void setRuleName(String name)
    {
        this.cutRul = name;
    }

//------------------------------------------------------------------------------

/**
 * Set the flag for both atoms with the same class:subclass
 */
    public void setSymmetricSubClass(boolean flg)
    {
        this.isSymm = flg;
    }

//------------------------------------------------------------------------------

/**
 * Returns the string representing this target bond
 */
    public String toString()
    {
	String str = idAtm0+"-"+classSubClass0+"_"+idAtm1+"-"+classSubClass1;

        return str;
    }

//------------------------------------------------------------------------------

}
