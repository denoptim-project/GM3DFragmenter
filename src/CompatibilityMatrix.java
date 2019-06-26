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

import java.io.*;
import java.lang.*;
import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;
import java.util.Date;
import java.text.DateFormat;
import java.text.SimpleDateFormat;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.ChemObject;

/**
 * Compatibility Matrix defines the compatibility 
 * between attachment point classes. It also defines the relation 
 * between attachment point class and bond order.
 *
 * @author Marco Foscato (University of Bergen)
 */

public class CompatibilityMatrix
{

    /** 
     * The map of TRUE entries (TRUE = classes are compatible)
     */
    protected static Map<String, ArrayList<String>> mp = new HashMap<String, ArrayList<String>>();

    //Keyword for compatibility matrix to text tile convertion
    private static String mpKey = "RCN";

    //Keyword for CLASS-to-bond order relation
    private static String boKey = "RBO";

    //Keyword for capping rules
    private static String capKey = "CAP";

//------------------------------------------------------------------------------
    /**
     * Creates an empty compatibility matrix
     */

    public CompatibilityMatrix()
    {
        this.mp.clear();
    }

//------------------------------------------------------------------------------

    /**
     * Add a signle true entry to the compatibility matrix
     * @param parentClass class in parent field (row index of the matrix)
     * @param childClass class in the child field (column index of the matrix)
     */

    public void addTrueEntry(String parentClass, String childClass)
    {
        if (!mp.keySet().contains(parentClass)) {
            ArrayList<String> entry = new ArrayList<String>();
            entry.add(childClass);
            mp.put(parentClass,entry);
        } else {
            if (!mp.get(parentClass).contains(childClass))
               mp.get(parentClass).add(childClass);
        }
    }

//------------------------------------------------------------------------------

    /**
     * Check the compatibility between two classes
     * @return <code>true</code> if the two classes are compatible
     */

    public boolean isCompatible(String parentclass, String childclass)
    {
        ArrayList<String> compclasses = mp.get(parentclass);
        boolean isthere = compclasses.contains(childclass);
        return isthere;
    }

//------------------------------------------------------------------------------

    /**
     * @param key check compatibility with this class
     * @return the list of classes which ar compatible with key
     */

    public ArrayList<String> getCompClasses(String key)
    {
        if (!mp.containsKey(key))
            return null;
        return mp.get(key);
    }        

//------------------------------------------------------------------------------

    /**
     * Read the compatibility matrix from file
     * @param parfile name of the file to read
     */
    public void fillCPMap(String parfile)
    {
        ArrayList<String> table = readKEYFile(parfile,mpKey);
        for (String line : table)
        {
            String[] words = line.split(" ");
	    String parentName = words[1];
            ArrayList<String> allvalue = new ArrayList<String>();
            if (words.length > 2)
            {
		for (int i = 2; i<words.length; i++)
                    allvalue.add(words[i]);
            } 
            mp.put(parentName,allvalue);
        }
    }

//------------------------------------------------------------------------------

    /**
     * Print the compatibility matrix on standard output
     */
    public void printCPMap()
    {
        Set<String> keys = mp.keySet();
        List<String> keylist = new ArrayList<String>(keys);
        Collections.sort(keylist);
        for(String s : keylist)
             System.out.println("Comp-Classes for <"+s+"> are: "+mp.get(s));
    }

//------------------------------------------------------------------------------

    /**
     * Get the CLASS-to-BondOrder conversion key from file
     * @param filename file to read 
     */
    public static Map<String,Integer> getCL2BO(String filename)
    {
        Map<String,Integer> mp = new HashMap<String,Integer>();
        ArrayList<String> table = readKEYFile(filename,boKey);
        for (String s : table) 
	{
            String[] ent = s.split("\\s");
            mp.put(ent[1],Integer.parseInt(ent[2]));
        }
        return mp;
    }

//------------------------------------------------------------------------------

    /**
     * Get capping rules from text file 
     * @param filename file to read
     */
    public static Map<String,Set<String>> getCAPP(String filename)
    {
        Map<String,Set<String>> mp = new HashMap<String,Set<String>>();
        ArrayList<String> table = readKEYFile(filename,capKey);
        for (String s : table) 
	{
            String[] ent = s.split("\\s");
            Set<String> val = new HashSet<String>();
            val.addAll(Arrays.asList(ent[2].split(",")));
            mp.put(ent[1],val);
        }
        return mp;
    }

//------------------------------------------------------------------------------

    /**
     * Reads the lines starting with the keyword
     * @param filename name of the file to read
     * @param key first word in the line is the keyword
     * @return list of lines from file
     */
    public static ArrayList<String> readKEYFile(String filename, String key)
    {
        ArrayList<String> table = new ArrayList<String>();
        try {
            BufferedReader fileinp =  new BufferedReader(new FileReader(filename));
            String line = null;
            while ((line = fileinp.readLine()) != null) 
            {
                int sp = key.length();
                if (line.length() > sp && sp > 0) 
                    if (line.substring(0,sp).equals(key))
                        table.add(line);
                    else
                       continue;
            }
        } catch (Throwable t) {
            System.out.println("While reading KEY file: Exception occours! " + t.getMessage());
        }
        return table;
    }

//------------------------------------------------------------------------------

    /**
     * Writes the compatibility matrix on a text file
     * @param filename name of the output file
     */
    public static void writeCPMapFile(String filename)
    {
        DateFormat df = new SimpleDateFormat("dd/MM/yyyy");
        String date = df.format(new Date());

        // Print head of the file
        String head = "#\n# CompatibilityMatrix for Class Based Builders\n#";
	head = head + "\n# Created by GM3DFragmenter - [d/m/y] " + date; 
        head = head + "\n# Format";
        head = head + "\n# parentClass compatibleClass1:subClass ";
        head = head + "compatibleClass2:subClass [etc.]\n#";

        writeToFile(head,filename,true);

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

            writeToFile(line.toString(),filename,true);
        }
        writeToFile(allCl.toString(),filename,true);
    }

//------------------------------------------------------------------------------

    /**
     * Writes the compatibility matrix on a text file and appends the 
     * CLASS-BondOrder correspondence 
     * @param filename name of the output file
     * @param classBndOrd CLASS-BondOrder conversion key
     */
    public static void writeCPMapFile(String filename, Map<String,Integer> classBndOrd)
    {
        // Compatibility matrix
        writeCPMapFile(filename);

        // CLASS-BondOrder conversion key
        writeToFile("# CLASS-to-BondOrder conversion",filename,true);
        Set<String> keys = classBndOrd.keySet();
        List<String> keylist = new ArrayList<String>(keys);
        Collections.sort(keylist);
        for (String cl : keylist) 
	{
            String line = " "+Integer.toString(classBndOrd.get(cl));
            writeToFile(boKey+" "+cl+line,filename,true);
        }
    }

//------------------------------------------------------------------------------

    /**
     * Writes the compatibility matrix on a text file and appends the 
     * CLASS-BondOrder correspondence 
     * @param filename name of the output file
     * @param classBndOrd CLASS-BondOrder conversion key
     * @param cap section
     */
    public static void writeCPMapFile(String filename, Map<String,Integer> classBndOrd, Map<String,String> cap)
    {
        // Compatibility matrix
        writeCPMapFile(filename);

        // CLASS-BondOrder conversion key
        writeToFile("# CLASS-to-BondOrder conversion",filename,true);
        Set<String> keys = classBndOrd.keySet();
        List<String> keylist = new ArrayList<String>(keys);
        Collections.sort(keylist);
        for (String cl : keylist)
        {
            String line = " "+Integer.toString(classBndOrd.get(cl));
            writeToFile(boKey+" "+cl+line,filename,true);
        }

	// CAPPING
        writeToFile("# Capping groups",filename,true);
        Set<String> ck = cap.keySet();
        List<String> ckl = new ArrayList<String>(ck);
        Collections.sort(ckl);
        for (String cl : ckl)
        {
	     String line = " "+cap.get(cl);
            writeToFile(capKey+" "+cl+line,filename,true);
        }
    }

//------------------------------------------------------------------------------
    /**
     * Write a string to a text tile
     * @param s string to write
     * @param filename text file to be written
     * @param append if <code>true</code> string will be appended to existing text 
     */
    private static void writeToFile(String s, String filename, boolean append) 
    {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(filename), append));
            bw.write(s);
            bw.newLine();
            bw.close();
        } catch (Exception e) {
            System.out.println("Problems while writing on file "+filename);
        }
    }

//------------------------------------------------------------------------------
}
