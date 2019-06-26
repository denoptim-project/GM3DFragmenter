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

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.FileReader;
import java.io.File;
import java.io.IOException;
import java.io.FileNotFoundException;

import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;

import org.openscience.cdk.ChemObject;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.exception.CDKException;

/**
 * Tools for Input/Output
 * 
 * @author Marco Foscato (University of Bergen)
 */

public class IOtools
{
//------------------------------------------------------------------------------

/**
 * Reads TXT files (Suitable for small files - do NOT use this for huge files!)
 * @param filename file to be read
 * @return all the lines into an <code>ArrayList</code>
 */
    public static ArrayList<String> readTXT(String filename)
    {
	ArrayList<String> allLines = new ArrayList<String>();
        BufferedReader buffRead = null;
        try {
            buffRead = new BufferedReader(new FileReader(filename));
            String line = null;
            while ((line = buffRead.readLine()) != null) 
	        allLines.add(line);
        } catch (FileNotFoundException fnf) {
            System.err.println("File Not Found: " + filename);
            System.err.println(fnf.getMessage());
            System.exit(-1);
        } catch (IOException ioex) {
            System.err.println(ioex.getMessage());
            System.exit(-1);
        } finally {
            try {
                if (buffRead != null)
                    buffRead.close();
            } catch (IOException ioex2) {
                System.err.println(ioex2.getMessage());
                System.exit(-1);
            }
        }
	return allLines;
    }
//------------------------------------------------------------------------------

/**
 * Reads keyword-labelled lines in TXT files (Suitable for small files - do NOT use this for huge files!)
 * @param filename file to be read
 * @param keyword labedel identifying the wanted lines
 * @return all the lines biginning with the <code>keyword</code> into an <code>ArrayList</code>
 */
    public static ArrayList<String> readTXTKeyword(String filename,String keyword)
    {
        ArrayList<String> allLines = new ArrayList<String>();
        BufferedReader buffRead = null;
        try {
            buffRead = new BufferedReader(new FileReader(filename));
            String line = null;
            while ((line = buffRead.readLine()) != null)
	    {
		String[] words = line.split("\\s+");
		if (words[0].equals(keyword))
                    allLines.add(line);
	    }
        } catch (FileNotFoundException fnf) {
            System.err.println("File Not Found: " + filename);
            System.err.println(fnf.getMessage());
            System.exit(-1);
        } catch (IOException ioex) {
            System.err.println(ioex.getMessage());
            System.exit(-1);
        } finally {
            try {
                if (buffRead != null)
                    buffRead.close();
            } catch (IOException ioex2) {
                System.err.println(ioex2.getMessage());
                System.exit(-1);
            }
        }
        return allLines;
    }

//------------------------------------------------------------------------------

/**
 * Reads SDF files (chemical format)
 * @param filename SDF file to be read
 * @return all the chemical objects into an <code>ArrayList</code>
 */

    public static ArrayList<IAtomContainer> readSDF(String filename)
    {
        MDLV2000Reader mdlreader = null;
        ArrayList<IAtomContainer> lstContainers = new ArrayList<IAtomContainer>();
        try {
            mdlreader = new MDLV2000Reader(new FileReader(new File(filename)));
            ChemFile chemFile = (ChemFile) mdlreader.read((ChemObject) new ChemFile());
            lstContainers.addAll(ChemFileManipulator.getAllAtomContainers(chemFile));
        } catch (Throwable t2) {
            System.err.println("Failure in reading SDF: " + t2);
            System.exit(-1);
        }

        return lstContainers;
    }

//------------------------------------------------------------------------------

/**
 * Writes on SDF a new file or appends to an existing file
 * @param filename target SDF file (new or existing)
 * @param mol atom container to be written on the SDF file
 * @param append <code>true</code> to append to existing file
 */

    public static void writeSDFAppend(String filename, IAtomContainer mol,
                                      boolean append)
    {
        SDFWriter sdfWriter = null;
        try {
            sdfWriter = new SDFWriter(new FileWriter(new File(filename), append));
            sdfWriter.write(mol);
	} catch (CDKException e) {
	    if (e.getMessage().contains("For input string: \"#\""))
	    {
		System.err.println("CDK unable to write MDL file for "+MolecularUtils.getNameOrID(mol));
	    }
	    
        } catch (Throwable t2) {
            System.err.println("Failure in writing SDF: " + t2);
            System.exit(-1);
        } finally {
             try {
                 if(sdfWriter != null)
                     sdfWriter.close();
             } catch (IOException ioe) {
                 System.err.println("Error in writing: " + ioe);
                 System.exit(-1);
             }
        }
    }

//------------------------------------------------------------------------------

/**
 * Writes on SDF a new file or appends to an existing file
 * @param filename target SDF file (new or existing)
 * @param mols Set of atom containers to be written on the SDF file
 * @param append <code>true</code> to append to existing file
 */

    public static void writeSDFAppendSet(String filename, IAtomContainerSet mols,
                                      boolean append)
    {
        SDFWriter sdfWriter = null;
        try {
            sdfWriter = new SDFWriter(new FileWriter(new File(filename), append));
            sdfWriter.write(mols);
        } catch (Throwable t2) {
            System.err.println("Failure in writing SDF: " + t2);
            System.exit(-1);
        } finally {
             try {
                 if(sdfWriter != null)
                     sdfWriter.close();
             } catch (IOException ioe) {
                 System.err.println("Error in writing: " + ioe);
                 System.exit(-1);
             }
        }
    }

//------------------------------------------------------------------------------

/**
 * Write a line on a TXT file
 * @param filename name of the target file
 * @param txt the line to be written
 * @param append <code>true</code> to append to existing file
 */
    public static void writeTXTAppend(String filename, String txt,
                                      boolean append)
    {
	FileWriter writer = null;
        try
        {
            writer = new FileWriter(filename,append); 
            writer.write(txt);
        } catch (Throwable t) {
            System.err.println("Failure in writing TXT: " + t);
            System.exit(-1);
        } finally {
             try {
                 if(writer != null)
                     writer.close();
             } catch (IOException ioe) {
                 System.err.println("Error in writing: " + ioe);
                 System.exit(-1);
             }
        }
    }

//------------------------------------------------------------------------------

/**
 * Copy the content of text files into other files
 * @param inFile path and name of the source file
 * @param outFile path and name of the target file
 */

    public static void copyFile(String inFile, String outFile, boolean append)
    {
	try {
        InputStream inStr = new FileInputStream(inFile);
	OutputStream outStr = new FileOutputStream(outFile,append);
	byte[] buf = new byte[1024];
	int len;
	while ((len = inStr.read(buf)) > 0)
	{
            outStr.write(buf, 0, len);
        }
	inStr.close();
	outStr.close();
        } catch (Throwable t) {
            System.err.println("\nERROR in copying file "+inFile+": "+t);
            System.exit(0);
        } 
    }

//------------------------------------------------------------------------------

/**
 * Delete a file 
 * @param fileName path and name of the file to remove
 */

    public static void deleteFile(String fileName)
    {
	File f = new File(fileName);
	if (f.exists())
	    if (f.canWrite())
		f.delete();
    }

//------------------------------------------------------------------------------

/**
 * Read an SD/SDF file as text and make a copy including only selected molecules
 */
    public static void filterSDasTXT(String inFile, String outFile, List<Boolean> filter)
    {
        BufferedReader buffRead = null;
        try {
            buffRead = new BufferedReader(new FileReader(inFile));
            String line = null;
	    int i = 0;
            while ((line = buffRead.readLine()) != null)
            {
		if (line.startsWith("$$$$"))
		{
		    if (filter.get(i))
		        writeTXTAppend(outFile,line+"\n",true);
		    i++;
		    
		} else if (filter.get(i)) {
		    writeTXTAppend(outFile,line+"\n",true);
		}
            }
        } catch (FileNotFoundException fnf) {
            System.err.println("File Not Found: " + inFile);
            System.err.println(fnf.getMessage());
            System.exit(0);
        } catch (IOException ioex) {
            System.err.println(ioex.getMessage());
            System.exit(0);
        } finally {
            try {
                if (buffRead != null)
                    buffRead.close();
            } catch (IOException ioex2) {
                System.err.println(ioex2.getMessage());
                System.exit(0);
            }
        }

    }

//------------------------------------------------------------------------------

/**
 * Stop execution, Ask the user to pres RETURN key and move on
 */

    public static void pause()
    {
        System.err.println("Press <RETURN> to continue");
        try
        {
            int inchar = System.in.read();
        }
        catch (IOException e)
        {
            System.err.println("Error reading from user");
        }
    }

//------------------------------------------------------------------------------
}
