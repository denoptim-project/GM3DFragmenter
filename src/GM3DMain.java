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

import java.io.IOException;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.File;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;

import java.util.Date;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import org.openscience.cdk.ChemFile;


/**
 * GM3DMain is the entering door to access GM3DFragmenter and related tools.
 * The suite is meant for generation and managment of molecular fragments 
 * that are generated from molecular structures (2D or 3D)
 * according to a given list of cutting rules.
 *
 * @author Marco Foscato (University of Bergen) Discussions: Giovanni Occhipinti and Vidar R. Jensen (University of Bergen)
 */

public class GM3DMain
{
   public static void main(String[] args) throws Exception
   {
        //Welcome message
        writeHeader();
	//Set up all parameters
	Parameters.setAll(args);

	//Check the quality of the 3D structures from CCDC
	StructureChecker stuChk = new StructureChecker();
        if (Parameters.chkFormula)
        {
            stuChk.check3DvsFormula();
        }
	if (Parameters.analyzeAndFix) 
	{
            stuChk.fixChemicalRepresentation();
        }

	//Prefilter
	if (Parameters.preFiltering)
	{
	    stuChk.preFilter();
	}

	//Now chop the molecules generating the fragments
	if (Parameters.chopMols)
	{
	    Fragmenter fr = new Fragmenter();
	    fr.chopMolecules();
	}

	Librarian lib = new Librarian();
        //Applay rejection rules to library of fragments
        if (Parameters.onlyFiltering)
            lib.filter();

	//Split library in small MW-range pieces
	if (Parameters.MWsplitting)
	    lib.mwSplitting();

	//Merge MW-splitted libraries
	if (Parameters.MWMerge)
	    lib.mwMerge();

        //Merge libraries
        if (Parameters.mergeLibraries)
            lib.mergeLibraries();

        //Reorder according to molecular weight
        if (Parameters.orderMW)
            lib.reorderMW(Parameters.MWascending);

	//Group rotamers
        if (Parameters.groupingRotamers & !Parameters.mergeLibraries)
            lib.groupRotamers();

	//Extract Classes
	if (Parameters.extractClass)
	    lib.extractClass();
	if (Parameters.extractSMARTS)
            lib.extractSMARTS();

	//Managment of compatibility matrix
	CPMapManager cpmm = new CPMapManager();
	if (Parameters.makeCPMap)
	    cpmm.makeFromCuttingRules();

	//Delete temporary files
	if (Parameters.removeIntermFiles)
	{
            if (Parameters.report >= 1)
		System.out.println("\nRemoving intermediare files");
	    Set<String> list = Parameters.getIntermediateFiles();
	    for (String fName : list)
		IOtools.deleteFile(fName);
	}

        //Goodbye messages
	printTail();
   }

//------------------------------------------------------------------------------

/**
 * Writes the header with version information and welcome message
 */

   private static void writeHeader()
   {
        String GM3DFVersion="0.1.0";
        while (GM3DFVersion.length() < 8)
           GM3DFVersion=GM3DFVersion+" ";
         
        System.out.println("\n");
        System.out.println("#####################################################################");
        System.out.println("##                                                                 ##");
        System.out.println("##                         GM3DFragmenter                          ##");
        System.out.println("##                                                                 ##");
        System.out.println("##         Generation of 3-Dimensional Molecular Fragments         ##");
        System.out.println("##               from crystallographic data (CCDC).                ##");
        System.out.println("##                                                                 ##");
        System.out.println("##                         Version  "+GM3DFVersion+"                       ##");
        System.out.println("##                                                                 ##");
        System.out.println("##                     Author: Marco Foscato                       ##");
        System.out.println("##      Discussions: Giovanni Occhipinti and Vidar R. Jensen       ##");
        System.out.println("##                 University of Bergen - Norway                   ##");
        System.out.println("##                                                                 ##");
        System.out.println("#####################################################################\n");
   }

//------------------------------------------------------------------------------

    private static void printTail()
    {
        System.out.println("\n=============== Task completed ===============");
	if (Parameters.singleOutput)
            System.out.println("Final output file: "+Parameters.getCurrentSDFile());
        System.out.println("\nThanks for using Foscato's sripts and softwares. \nMandi! :) \n");
    }

}
