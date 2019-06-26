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

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.formula.MolecularFormula;
import org.openscience.cdk.Isotope;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;


/**
 * ElementalAnalyser compare the number of atoms reported in an 
 * <code>AtomContainer</code>
 * with the number reported in the molecular formula to check whether
 * the system has missing atoms or misstyped entries. 
 * ElementalAnalyser is ment to work with information retrived from 
 * The Cambridge Crystallographic Data Centre (CCDC) 
 * [see: www.ccdc.cam.ac.uk].
 * 
 * @author Marco Foscato (University of Bergen)
 */

public class ElementalAnalyser
{
    //Presence of tunable stoichiometric factors
    private boolean doTuning;

    //Define the maximum 'n' for tunable stoichiometric factors
    private int limSF = 6;

    //Contains only the largest molecule
    private boolean largest;
    //MAtching problems
    private boolean matching;

    //Formula
    private String formula;

    //Chemical system == Molecular informations
    private IAtomContainer mol;

    //Elemetal analysis
    private Map<String,ArrayList<Double>> elemAnalFormula = new HashMap<String,ArrayList<Double>>();
    private Map<String,Double> elemAnalMolInfo = new HashMap<String,Double>();

    //Comparison
    private Map<String,Boolean> SDFormulaAgreeOnElement = new HashMap<String,Boolean>();
    private boolean allElementsAgree;

    //Output file: rejected molecules
    private String checkfile;
    private boolean someMolRejected;
    private int numRejected;

    //Reporting flag
    private int repOnScreen;

//------------------------------------------------------------------------------

/**
 * Constructor
 *
 * @param formula is the chemical formula as reported in CCDC under the 'Chemical' 
 * type of information. Get this info from a TXT output of ConQuest-1.14
 * @param mol molecular system. It contains all the information on the chemical 
 * system as reported in SD output format from ConQuest-1.14.
 */

    public ElementalAnalyser(String formula, IAtomContainer mol)
    {
        repOnScreen = Parameters.report;
        this.mol = mol;
        this.formula = formula;
        this.doTuning = false;
        this.largest = onlyLargest(mol);
	this.matching = matchingProblems(mol);
	checkfile = "check-ElementalAnalyser_"+Parameters.getJobName()+".sdf";

	if (repOnScreen >= 1)
	    System.out.println("Elemental Analysis of "+formula);

	try {
	    //Count atoms in SD and Formula
	    AnalysisOfMolInfo();
	    AnalysisOfFormula();
            //Compare entries
            compareCounts();
	} catch (Throwable t) {
	    rejectMol(mol,"Faiulure in Elemental Anaysis: "+t);
	    allElementsAgree = false;
	}
    }

//------------------------------------------------------------------------------

/**
 * @return <code>true</code> if the molecular formula reflects the atom composition 
 * for every element in the system
 */

    public boolean getAllElementsAgree()
    {
	return allElementsAgree;
    }

//------------------------------------------------------------------------------

/**
 * @return the overview over all atom on the agreement between the number of
 * atoms defined in the molecular formula and the number of atoms in the coordinate 
 * section.
 */

    public Map<String,Boolean> getSDFormulaAgreeOnElement()
    {
        return SDFormulaAgreeOnElement;
    }

//------------------------------------------------------------------------------

/**
 * Check if for an element the numer of atoms in the coordinate section correspond
 * to the number of atoms expected from the molecular formula.
 * @param element is the symbol of the element to query.
 * @return <code>false</code> if SD and Formula disagree or the element is not in the list
 */

    public boolean getSDFormulaAgreementElement(String element)
    {
	if (!SDFormulaAgreeOnElement.containsKey(element))
	   return false;
        return SDFormulaAgreeOnElement.get(element);
    }

//------------------------------------------------------------------------------

/**
 * @return the result of elemental analysis on the formula including the attempt of
 * tuning the variable stoichiometric factor, if any [I.e. n(H2O),2n(CO2)]
 */

    public Map<String,ArrayList<Double>> getElemAnalFromFormula()
    {
        return elemAnalFormula;
    }

//------------------------------------------------------------------------------

/**
 * @return  the result of elemental analysis on the molecular system from the SD file
 */

    public Map<String,Double> getElemAnalFromSD()
    {
        return elemAnalMolInfo;
    }

//------------------------------------------------------------------------------

/**
 * Define if the SD file provides only the largest molecule according to ConQuest
 * @param molecule chemical object as Iatomcontainer
 * @return <code>true</code> if the atoms list consider only the largest molecule
 */

    public boolean onlyLargest(IAtomContainer molecule)
    {
        //Read remarks 
	//suitable for ConQuest-1.14
        String remarks = molecule.getProperty("cdk:Remark").toString();

        //Decide
        boolean large = false;
        if (remarks.contains("Largest molecule"))
            large = true;

        if (repOnScreen >= 2) 
            System.out.println("Largest molecule: "+large);

        return large;
    }

//------------------------------------------------------------------------------

/**
 * Define if the SD file is labelled for matching problems by ConQuest-1.14
 * @param molecule chemical system as provided by ConQuest-1.14 (SD output file)
 * @return <code>true</code> if the atoms list has matching problems
 */

    public boolean matchingProblems(IAtomContainer molecule)
    {
        //Read remarks 
        String remarks = molecule.getProperty("cdk:Remark").toString();

	boolean found = false;
        if (remarks.contains("Matching problem"))
            found  = true;

        if (repOnScreen >= 2)
            System.out.println("Matching problems: "+found);

        return found;
    }

//------------------------------------------------------------------------------

/**
 * Compare the atom counting for every element in the system
 */

    private void compareCounts()
    {
        if (repOnScreen >= 3)
            System.out.println("Comparing the atom number per each element");

	allElementsAgree = true;

//System.out.println("Results of atom number comparison\n"+SDFormulaAgreeOnElement);


	Map<String,ArrayList<Integer>> whichOneAgreed = new HashMap<String,ArrayList<Integer>>();
        //loop over elements
        for (String chkEl : elemAnalFormula.keySet())
        {
	    if (repOnScreen >= 3)
                System.out.println("Comparing entries for element "+chkEl);

	    SDFormulaAgreeOnElement.put(chkEl,false);
	    ArrayList<Integer> idexOfFound = new ArrayList<Integer>();

/*	    //really bad case: no such atom in coordinates section!
	    if (!elemAnalMolInfo.containsKey(chkEl))
	    {
		allElementsAgree = false;
		continue;
	    }
*/	 
	    //No such atom in coordinates file 
            if (!elemAnalMolInfo.containsKey(chkEl))
		elemAnalMolInfo.put(chkEl,0.0);
 
//System.out.println("HERE: "+elemAnalFormula.get(chkEl)+" "+elemAnalMolInfo.get(chkEl));
	    //Check for agreement over the list of possibilities
	    if (elemAnalFormula.get(chkEl).contains(elemAnalMolInfo.get(chkEl)))
	    {
		//which entry is in agreement?
		for (int i = 0; i< elemAnalFormula.get(chkEl).size(); i++)
		{
		    if (elemAnalFormula.get(chkEl).get(i).equals(elemAnalMolInfo.get(chkEl)))
		    {
			idexOfFound.add(i);
		        SDFormulaAgreeOnElement.put(chkEl,true);
		        if (repOnScreen >= 3) 
        	            System.out.println("Found agreement between Formula and SD in "+i);
		    }
		}
	    }
/*
	    } else {
                for (int x = 1; x < limSF; x++)
                {
		    double xd = x;
                    if ((elemAnalMolInfo.get(chkEl) * xd) == elemAnalFormula.get(chkEl).get(0))
                    {
			idexOfFound.add(0);
                        SDFormulaAgreeOnElement.put(chkEl,true);
//                        if (repOnScreen >= 3)
                            System.out.println("Found agreement between Formula and SD using Stoich. pre-Fact. 1/"+x);
			break;
                    }
                }
	    }
*/
	    whichOneAgreed.put(chkEl,idexOfFound);
	    if (!SDFormulaAgreeOnElement.get(chkEl))
		allElementsAgree = false;
	}

	//In case of all elements show agreement with one of the candidates
	//check if the agreement is found for the same candidate for all the elements
	if (allElementsAgree)
	{
	    if (repOnScreen >= 3)
	        System.out.println("Check for common agreement between elements");

	    //get the shortest list of agreements
	    String shortest = "";
	    int shortestSize = -1;
	    for (String key : whichOneAgreed.keySet())
	    {
		if ((whichOneAgreed.get(key).size() < shortestSize) | (shortestSize == -1))
		{
		    shortestSize = whichOneAgreed.get(key).size();
		    shortest = key;
		}
	    }
	    if (repOnScreen >= 3)
	        System.out.println("Shortest list of agreements found for "+shortest+" ("+shortestSize+" cases)");
	    //Look for common entries through elements
	    boolean hasCommon = true;
	    int idOfCommon = -1;
	    for (int idx : whichOneAgreed.get(shortest))
	    {
	        if (repOnScreen >= 3)
                    System.out.println("Cheking list "+idx);
		hasCommon = true;
		for (String key : whichOneAgreed.keySet())
		{
		     if (repOnScreen >= 3)
 	                 System.out.println("Does "+key+" has agreement in list "+idx+"? "+whichOneAgreed.get(key).contains(idx));
		    if (!whichOneAgreed.get(key).contains(idx))
		    {
			hasCommon = false;
			break;
		    }
		}
		if (hasCommon)
		{
		    idOfCommon = idx;
		    if (repOnScreen >= 3)
                        System.out.println("Found agreement for list "+idx);
		    break;
		}
	    }
	    if (!hasCommon)
	    {
		if (repOnScreen >= 2)
		    System.out.println("Other elements do not mach list of "+shortest);
		allElementsAgree = false;
	    }
	}

	if (repOnScreen >= 1)
	{
            System.out.println("Agreement found for these elements (in at least one configuration)\n"+SDFormulaAgreeOnElement);
            System.out.println("Full agreement (SD<=>Formula): "+allElementsAgree);
	}
    }

//------------------------------------------------------------------------------

/**
 * Counts the number of atoms for each element in the system
 */

    private void AnalysisOfMolInfo()
    {
	if (repOnScreen >= 2)
            System.out.print("\nAnalysing Molecular Information from SD... ");

	boolean init = true;
        for (IAtom atm : mol.atoms())
        {
	    String elSymbol = atm.getSymbol();
	    //Deal with deuterium symbol
	    if (atm.getMassNumber() == 2)
		elSymbol = "D";
	    if (elemAnalMolInfo.keySet().contains(elSymbol))
  	    {
	        double num = elemAnalMolInfo.get(elSymbol) + 1.0;
	        elemAnalMolInfo.put(elSymbol,num);
	    } else {
	        elemAnalMolInfo.put(elSymbol,1.0);
	    }
        }

        if (repOnScreen >= 2)
            System.out.println("done.");
        if (repOnScreen >= 1)
            System.out.println("Atom Counts from from SD/SDF file: "+elemAnalMolInfo);
    }

//------------------------------------------------------------------------------

/**
 * Counts the number of atoms for each element in the molecular formula
 */

    private void AnalysisOfFormula() throws Throwable
    {
        if (repOnScreen >= 2)
           System.out.println("\nAnalysis of Formula "+formula);

        //Create container for counting atoms over the list of molecules
        ArrayList<Map<String,Integer>> allAtmCounts = new ArrayList<Map<String,Integer>>();

        //identify all isolated molecules
        String[] mols = formula.split(",");

        //Deal with stoichiometric factors
        List<Double> stocFact = new ArrayList<Double>();
        List<Boolean> tuneStFact = new ArrayList<Boolean>();
        List<Boolean> tuneStFactToInt = new ArrayList<Boolean>();
        for (int i=0; i < mols.length; i++)
        {
            if (repOnScreen >= 2)
                System.out.println("Decripting isolated contribution: "+mols[i]);
            stocFact.add(i,1.0);
            tuneStFact.add(i,false);
	    tuneStFactToInt.add(i,false);
            String locForm = mols[i];
            boolean found = false;
            //Check for stoichiometric factor
            if (mols[i].contains("("))
            {  //There is a stoichiometric factor
                if (repOnScreen >= 2)
                    System.out.println("Stoich. Fact. => There is a stoichiometric factor");
                if (mols[i].lastIndexOf(")") == (mols[i].length()-1))
                {  // the Stoc. Fact. is at the beginning => nothing at the end!
                    if (repOnScreen >= 2)
                        System.out.println("Stoich. Fact. => Stoc. Fact. is at the beginning");
                    String[] spt = mols[i].split("[()]");
                    if (spt[0].contains("n"))
                    {  //Tunable Stoic. Fact.
                        if (repOnScreen >= 2)
                            System.out.println("Stoich. Fact. => Tunable Stoic. Fact. (n)");
                        if (spt[0].lastIndexOf("n") == 0)
                        {  // I see only N
                            if (repOnScreen >= 2)
                                System.out.println("Stoich. Fact. => I see only N in position 0");
                            //n-case 0: __n(formula)__
                            tuneStFact.set(i,true);
                            locForm = spt[1];
                            found = true;
                        } else if (spt[0].lastIndexOf("n") == (spt[0].length()-1))
                        {  //NUM times n
                            if (repOnScreen >= 2)
                                System.out.println("Stoich. Fact. => NUM times n ");
                            //n-case 1: __NUMn(formula)__
                            stocFact.set(i,Double.parseDouble(spt[0].substring(0,spt[0].length() -1)));
                            tuneStFact.set(i,true);
                            locForm = spt[1];
                            found = true;
                        }
                    } else if (spt[0].contains("x"))
                    {  //Tunable Stoic. Fact. - other version
                        if (repOnScreen >= 2)
                            System.out.println("Stoich. Fact. => Tunable Stoic. Fact. (x)");
                        if (spt[0].lastIndexOf("x") == 0)
                        {  // I see only X
                            if (repOnScreen >= 2)
                                System.out.println("Stoich. Fact. => I see only X in position 0");
                            //n-case 0: __x(formula)__
                            tuneStFact.set(i,true);
                            locForm = spt[1];
                            found = true;
                        } else if (spt[0].lastIndexOf("x") == (spt[0].length()-1))
                        {  //NUM times x
                            if (repOnScreen >= 2)
                                System.out.println("Stoich. Fact. => NUM times x ");
                            //n-case 1: __NUMx(formula)__
                            stocFact.set(i,Double.parseDouble(spt[0].substring(0,spt[0].length() -1)));
                            tuneStFact.set(i,true);
                            locForm = spt[1];
                            found = true;
                        }
                    } else {
                        if (repOnScreen >= 2)
                            System.out.println("Stoich. Fact. Regular Case => NUM times (formula)");
                        //Regular case: __NUM(formula)__
                        stocFact.set(i,Double.parseDouble(spt[0]));
                        locForm = spt[1];
                        found = true;
                    }
                } else if (mols[i].lastIndexOf("(") == 0)
                {  //Stoic. Fact. at the end
                    if (repOnScreen >= 2)
                        System.out.println("Stoich. Fact. => Stoic. Fact. is at the end");
                    String[] sptR = mols[i].split("[()]");
                    //Take care! With this regular expression you get an empty string as first entry
//          TODO    System.out.println("Splitted in __"+sptR[0]+"_and_"+sptR[1]+"__"+sptR[2]+"_");
                    locForm = sptR[1];
 
                    //n-case 2: __(formula)n__
                    if (sptR[2].length() == 1)
                    {
                        if (repOnScreen >= 2)
                            System.out.println("Stoich. Fact. => Tunable Stoic. Fact.");
                        tuneStFact.set(i,true);
                        found = true;
                    } else if (sptR[2].lastIndexOf("x") == (sptR[0].length()-1)) {
			//NUM times x
  		       if (repOnScreen >= 2)
                            System.out.println("Stoich. Fact. => NUM times x ");
                        //n-case 1: __(formula)NUMx__
                        stocFact.set(i,Double.parseDouble(sptR[2].substring(0,sptR[2].length() -1)));
                        tuneStFact.set(i,true);
                        locForm = sptR[1];
                        found = true;
                    }

                }
 
                if (!found)
                {
                    System.out.println("\nERROR! Formula contains a syntax not recognized!!!"+mols[i]);
                    System.out.println("Please report this issue to the authors: marco.foscato@kj.uib.no\n");
                    System.exit(-1);
                }

		//Deal with fractionary Stoich. Fact.
		double mod = stocFact.get(i) % 1.0;
		if (mod != 0.0)
		    tuneStFactToInt.set(i,true);

                //Store extracted formula
                mols[i] = locForm;
 
            } //end of stochiometric factor analysis
 
            //Identify elements to be counted
            Map<String,Integer> locCount = new HashMap<String,Integer>();
            String [] lmnts = locForm.split("\\s+");
            for (int l = 0; l < lmnts.length; l++)
            {
                //Get rid of the carge, if any
                if (lmnts[l].endsWith("+") || lmnts[l].endsWith("-"))
                    continue;
 
                //Count of a single element
                String elSymbol = "";
                int elCount = 0;
                if (repOnScreen >= 3)
                    System.out.println("Decripting element "+lmnts[l]);
 
                int fd = 0;
                while (Character.isLetter(lmnts[l].charAt(fd)))
                    fd++;
                elSymbol = lmnts[l].substring(0,fd);
                elCount = Integer.parseInt(lmnts[l].substring(fd));
 
                //add this element and its count to the map
                locCount.put(elSymbol,elCount);
            }
            //Move the atoms count to the general storage
            allAtmCounts.add(i,locCount);
 
        } //end of loop over mols
 
        //Report results in analysing the formula
        if (repOnScreen >= 2)
        {
            System.out.println("Results from Formula "+formula);
            System.out.println("-> Molecules:  "+mols.length);
            System.out.println("-> tuneStFact: "+tuneStFact);
            System.out.println("-> tuneStFactToInt: "+tuneStFactToInt);
            System.out.println("-> stocFact:   "+stocFact);
            for (int i = 0; i < mols.length; i++)
                System.out.println("-> Formula "+i+":   "+mols[i]);
            System.out.println("-> Map of Counts: "+allAtmCounts);
        }

        //Find largest molecule (from formula weight)
        int largestMol = 0;
	double largestMass = 0.0;
        for (int i = 0; i < mols.length; i++)
        {
	    //calculater formula weigh
	    MolecularFormula molForm = new MolecularFormula();
            for (String el : allAtmCounts.get(i).keySet())
	    {
		IIsotope is = new Isotope(el);
		molForm.addIsotope(is,allAtmCounts.get(i).get(el));
	    }
	    MolecularFormulaManipulator mfm = new MolecularFormulaManipulator();
	    double mass = mfm.getNaturalExactMass(molForm);
	    if (mass > largestMass)
	    {
		largestMass = mass;
		largestMol = i;
	    }
        }
        if (repOnScreen >= 2)
             System.out.println("-> The largest molecule is "+largestMol+" with mass: "+largestMass);

        //In case of reporting only the largest molecule
/*        if (largest)
        {
            for (String el : allAtmCounts.get(largestMol).keySet())
            {
                ArrayList<Double> num = new ArrayList<Double>();
                if (tuneStFact.get(largestMol))
                {
                    for (int pf = 1; pf < limSF ; pf ++)
                    {
                        double pfd = pf;
			if (tuneStFactToInt.get(largestMol))
			{
			    pfd = pfd * (1.0 / stocFact.get(largestMol));
			}			    
                        num.add(allAtmCounts.get(largestMol).get(el) * pfd * stocFact.get(largestMol));
                    }
                } else {
		    double pfd = 1.0;
		    if (tuneStFactToInt.get(largestMol))
                    {
                        pfd = pfd * (1.0 / stocFact.get(largestMol));
                    }
                    num.add(allAtmCounts.get(largestMol).get(el) * pfd * stocFact.get(largestMol));
                }
                elemAnalFormula.put(el,num);
            }
        } else if (matching) //in case of matching problems one or more molecules may be missing
	{
*/
            //Get all the elements 
            Set<String> allEl = new HashSet<String>();
            for (int i = 0; i < mols.length; i++)
                for (String el : allAtmCounts.get(i).keySet())
                    allEl.add(el);

            //report counting for all SINGLE molecules
            for (int i = 0; i < mols.length; i++)
            {
                for (String el : allEl)
                {
                    double num = 0.0;
                    if (allAtmCounts.get(i).containsKey(el))
                    {
                        num = allAtmCounts.get(i).get(el);
                    }
                    if (!elemAnalFormula.containsKey(el))
                        elemAnalFormula.put(el,new ArrayList<Double>(Arrays.asList(num)));
                    else
                        elemAnalFormula.get(el).add(num);
                }
            }

	    //report counting for all SINGLE unit [NUM(molecule)]
	    for (int i = 0; i < mols.length; i++)
	    {
		for (String el : allEl)
		{
		    double num = 0.0;
		    if (allAtmCounts.get(i).containsKey(el))
		    {
			double pf = 1.0;
			if (tuneStFactToInt.get(i))
			    pf = pf * (1.0 / stocFact.get(i));
		        num = allAtmCounts.get(i).get(el) * pf * stocFact.get(i);
		    }
		    if (!elemAnalFormula.containsKey(el))
 		        elemAnalFormula.put(el,new ArrayList<Double>(Arrays.asList(num)));
		    else
			elemAnalFormula.get(el).add(num);
		}
	    }
	    //report counting for sums of the first n molecules up to all molecules
	    if (mols.length > 1)
	    {
	        for (int n = 2; n <= mols.length; n++)
	        {
		    //check for tunable stoich. factors in all molecules 
		    //(assuming ONLY ONE factor has to be tuned)
		    doTuning = true;
                    for (int i = 0; i < mols.length; i++)
                        if (!tuneStFact.get(i))
                        {
                            doTuning = false;
                            if (repOnScreen >= 2)
                                System.out.println("No tuning of stoichionetric factors according to mol "+i);
                        }

	            //Calculate number of atoms of per each element
		    for (String el : allEl)
                    {
			double thisElCount = 0.0;
			//sums of the first n molecules Ignoring stoichiometric factors
		        for (int i = 0; i < n; i++)
            	        {
                    	    if (allAtmCounts.get(i).containsKey(el))
			    {
				double pf = 1.0;
				if (tuneStFactToInt.get(i))
				    pf = pf * (1.0 / stocFact.get(i));
                        	thisElCount = thisElCount + allAtmCounts.get(i).get(el) * pf * stocFact.get(i);
			    }
            	        }
			elemAnalFormula.get(el).add(thisElCount);

	                //tune stoichiometric factors if required
        	        if (doTuning)
			{
			    //per each prefactor within the limits of the tuning procedure
			    for (int pf = 2; pf < limSF ; pf ++)
                    	    {
				double pfd = pf;
			        double thisElCountTune = 0.0;
				//sums of the first n molecules with stoichiometric factors
			        for (int i = 0; i < n; i++)
                                {
				    if (allAtmCounts.get(i).containsKey(el))
                                    {   
                                	thisElCountTune = thisElCountTune + allAtmCounts.get(i).get(el) * pfd * stocFact.get(i);
                            	    }
			        }
				elemAnalFormula.get(el).add(thisElCountTune);
			    }
			}
		    }
	        }
	    }
/*
        } else { //In case of ALL molecules
            //Get all the elements 
            Set<String> allEl = new HashSet<String>();
            for (int i = 0; i < mols.length; i++)
                for (String el : allAtmCounts.get(i).keySet())
                    allEl.add(el);
            //Sum over all molecules for each element
            for (String el : allEl)
            {
                //check for tunable stoich. factors
	        //If a S.F. for a molecule has to be tunned then all the others should be tuned according to the same factor
                doTuning = true;
                for (int i = 0; i < mols.length; i++)
                    if (!tuneStFact.get(i))
		    {
                        doTuning = false;
		        if (repOnScreen >= 2)
		        System.out.println("No tuning of stoichionetric factors according to mol "+i);
		    }

                ArrayList<Double> num = new ArrayList<Double>();
                //Calculate number of atoms for this element
                double fromAllMols = 0.0;
                for (int i = 0; i < mols.length; i++)
                    if (allAtmCounts.get(i).keySet().contains(el))
                        fromAllMols = fromAllMols + (allAtmCounts.get(i).get(el) * stocFact.get(i));
                num.add(fromAllMols);

                //tune siochiometric factors if required
                if (doTuning)
                {
                    for (int pf = 2; pf < limSF ; pf ++)
                    {
                         double pfd = pf;
                         num.add(num.get(0) * pfd);
                    }
                }

                //store the summ of atoms of this element
                elemAnalFormula.put(el,num);
            }
        }
*/

        if (repOnScreen >= 2)
        {
            System.out.println("-> RESULTS (with largest = "+largest+" - MatchingProblems = "+matching+")");
            for (String elm : elemAnalFormula.keySet())
                System.out.println("   Element: "+elm+" => "+elemAnalFormula.get(elm));
        }	
    }

//-----------------------------------------------------------------------------

    /**
     * Store molecules that will not be treated for some reason
     * @param mol input molecular/chemical object
     * @param reason brief explamention of the reason for the rejection
     */
    private void rejectMol(IAtomContainer mol, String reason)
    {
        someMolRejected = true;
        numRejected++;
        mol.setProperty("REJECTED",reason);
        IOtools.writeSDFAppend(checkfile,mol,true);
    }

//------------------------------------------------------------------------------
} 
