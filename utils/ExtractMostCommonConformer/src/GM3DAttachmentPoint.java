/**
 * Licence
 **/

import java.util.List;
import java.util.ArrayList;
import java.text.DecimalFormat;

/**
 * TODO
 *
 * @author Marco Foscato (University of Bergen)
 */

class GM3DAttachmentPoint
{
    //AtomID
    private int apAtm;

    //CLASS = RULE+SubCLASS
    private String apClass;

    //RULE
    private String apRule;

    //subCLASS
    private int apSubClass;

    //BondOrder
    private int apBndOrd;

    //Direction Vector (0 = X, 1 = Y, 2 = Z)
    private ArrayList<Double> apVector = new ArrayList<Double>(3);

    //Torsion
//    private 

    //AP string detailsSeparators
    private String defSeparator = ":";
    private String valueSeparator = "=";
    //Separators pro-format
    private String moreAtmsSeparator;
    private String moreAPSeparator;
    private String atmSeparator;
    private String detailsSeparator;
    private String coordSeparator;


//------------------------------------------------------------------------------

/**
 * TODO
 */

    public GM3DAttachmentPoint()
    {
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */

    public GM3DAttachmentPoint(int atm, String rule, int subCls, int bndOrd)
    {
        this.apAtm = atm;
	this.apRule = rule;
        this.apClass = rule + Integer.toString(subCls);
        this.apSubClass = subCls;
        this.apBndOrd = bndOrd;
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */

    public GM3DAttachmentPoint(int atm, String rule, int subCls, int bndOrd, ArrayList<Double> vector)
    {
        this.apAtm = atm;
        this.apRule = rule;
        this.apClass = rule + Integer.toString(subCls);
        this.apSubClass = subCls;
        this.apBndOrd = bndOrd;
        this.apVector = vector;
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */

    public GM3DAttachmentPoint(String apString)
    {
        String[] details = apString.split(defSeparator);
        for (String s : details)
        {
            String[] sParts = s.split(valueSeparator);
            String key = sParts[0];
            if (key.equals("ATM"))
                this.apAtm = Integer.parseInt(sParts[1]);
            else if (key.equals("RUL"))
                this.apRule = sParts[1];
            else if (key.equals("SCL"))
                this.apSubClass = Integer.parseInt(sParts[1]); 
            else if (key.equals("BNO"))
                this.apBndOrd = Integer.parseInt(sParts[1]);
            else if (key.equals("DX"))
                this.apVector.set(0,Double.parseDouble(sParts[1]));
            else if (key.equals("DY"))
                this.apVector.set(1,Double.parseDouble(sParts[1]));
            else if (key.equals("DZ"))
                this.apVector.set(2,Double.parseDouble(sParts[1]));
	    this.apClass = this.apRule + Integer.toString(this.apSubClass);
        }
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */

    public GM3DAttachmentPoint(String apStrFor, String format)
    {
        setSeparators();

        if (format.equals("DENOPTIM"))
        {
            try {
                String[] parts = apStrFor.split(atmSeparator);
                this.apAtm = Integer.parseInt(parts[0]);

                String[] details = parts[1].split(detailsSeparator);
                this.apRule = details[0];
                this.apSubClass = Integer.parseInt(details[1]);
		this.apClass = this.apRule + Integer.toString(this.apSubClass);

                String[] coord = details[2].split(coordSeparator);
                ArrayList<Double> pointer = new ArrayList<Double>();
                pointer.add(Double.parseDouble(coord[0]));
                pointer.add(Double.parseDouble(coord[1]));
                pointer.add(Double.parseDouble(coord[2]));
		this.apVector = pointer;
                //DENOPTIM format doesn't containg bond order so set to 1
                this.apBndOrd = 1;
	    } catch (Throwable t) {
		System.err.println("ERROR in getting AP details with format: "+format);
		t.printStackTrace();
                System.exit(0);
	    }
        }
//Add other formats here
//TODO
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */

    public GM3DAttachmentPoint(int id, String apStrFor, String format)
    {
        setSeparators();

        if (format.equals("DENOPTIM"))
        {
            try {
                this.apAtm = id;

                String[] details = apStrFor.split(detailsSeparator);
                this.apRule = details[0];
                this.apSubClass = Integer.parseInt(details[1]);
		this.apClass = this.apRule + Integer.toString(this.apSubClass);

                String[] coord = details[2].split(coordSeparator);
                ArrayList<Double> pointer = new ArrayList<Double>();
                pointer.add(Double.parseDouble(coord[0]));
                pointer.add(Double.parseDouble(coord[1]));
                pointer.add(Double.parseDouble(coord[2]));
		this.apVector = pointer;
		//DENOPTIM format doesn't containg bond order so set to 1
		this.apBndOrd = 1;
            } catch (Throwable t) {
                System.err.println("ERROR in getting AP details with format: "+format);
		t.printStackTrace();
                System.exit(0);
            }
        }
//Add other formats here
//TODO
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */

    private void setSeparators()
    {
//For now HARD CODED
this.moreAtmsSeparator = moreAtmsSeparator = " ";
this.moreAPSeparator = moreAPSeparator = ",";
this.atmSeparator = atmSeparator = "#";
this.detailsSeparator = detailsSeparator = ":";
this.coordSeparator = coordSeparator = "%";
/*
        this.moreAtmsSeparator = Parameters.moreAtmsSeparator;
        this.moreAPSeparator = Parameters.moreAPSeparator;
        this.atmSeparator = Parameters.atmSeparator;
        this.detailsSeparator = Parameters.detailsSeparator;
        this.coordSeparator = Parameters.coordSeparator;
*/
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */

    public void setAPAtm(int num)
    {
	apAtm = num;
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */

    public String getAPClass()
    {
        return apClass;
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */

    public String getAPRule()
    {
        return apRule;
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */

    public int getAPSubClass()
    {
        return apSubClass;
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */

    public int getAPAtm()
    {
        return apAtm;        
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */

    public int getAPBondOrder()
    {
        return apBndOrd;
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */

    public ArrayList<Double> getAPVector()
    {
        return apVector;
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */

    public boolean sameAs(GM3DAttachmentPoint other)
    {
        //Compare CLASS name
        String cl = this.getAPClass();
        String othercl = other.getAPClass();
        if (!cl.equals(othercl))
            return false;

        //Compare Bond Order that NOW should be equal by the definition of the CLASS
        int bo = this.getAPBondOrder();
        int otherbo = other.getAPBondOrder();
        if (bo != otherbo)
        {
            System.err.println("Unexpected difference in Bond Order while CLASS is equal!");
            System.exit(-1);
        }

        //Compare AtomID if CLASS is equal
        int atm = this.getAPAtm();
        int otheratm = other.getAPAtm();
        if (atm != otheratm)
            return false;

	return true;
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */

    public int compareTo(GM3DAttachmentPoint other)
    {
        final int BEFORE = -1;
        final int EQUAL = 0;
        final int AFTER = 1;
	
	if (this == other)
	    return EQUAL;

	//Compare CLASS name
	String cl = this.getAPClass();
	String othercl = other.getAPClass();
	int res = cl.compareTo(othercl);
	if (res < 0)
	    return BEFORE;
        else if (res > 0)
            return AFTER;

	//Since CLASS is equal
        //apRule and apSubClass must be equals by definition of GM3DAttachmentPoint
	//Also apBndOrd should be equals but let us check it

        //Compare Bond Order that NOW should be equal by the definition of the CLASS
        int bo = this.getAPBondOrder();
        int otherbo = other.getAPBondOrder();
	if (bo != otherbo)
	{
	    System.err.println("Unexpected difference in Bond Order while CLASS is equal!");
	    System.exit(-1);
	}

	//Compare AtomID if CLASS is equal
	int atm = this.getAPAtm();
	int otheratm = other.getAPAtm(); 
	if (atm < otheratm)
            return BEFORE;
        if (atm > otheratm)
            return AFTER;

        //Compare Direction Vector if AtomID is equal
        ArrayList<Double> vec = this.getAPVector();
        ArrayList<Double> othervec = other.getAPVector();
	//Compare X
	if (vec.get(0) < othervec.get(0))
            return BEFORE;
	else if (vec.get(0) > othervec.get(0))
            return AFTER;
        //Compare Y
        if (vec.get(1) < othervec.get(1))
            return BEFORE;
        else if (vec.get(1) > othervec.get(1))
            return AFTER;
        //Compare Z
        if (vec.get(2) < othervec.get(2))
            return BEFORE;
        else if (vec.get(2) > othervec.get(2))
            return AFTER;

	assert this.equals(other) : "GM3DAttachmentPoint.compareTo inconsistent with equals.";
	return EQUAL;
    }
//------------------------------------------------------------------------------

/**
 * TODO
 */

    public String toString()
    {
        String apString = "";
        apString = "ATM" + valueSeparator + Integer.toString(apAtm);
        apString = apString + defSeparator + "RUL" + valueSeparator + apClass;
        apString = apString + defSeparator + "SCL" + valueSeparator + apSubClass;
        apString = apString + defSeparator + "BNO" + valueSeparator + apBndOrd;
        apString = apString + defSeparator + "DX" + valueSeparator + apVector.get(0);
        apString = apString + defSeparator + "DY" + valueSeparator + apVector.get(1);
        apString = apString + defSeparator + "DZ" + valueSeparator + apVector.get(2);
//add other like Torsion
        return apString;
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */

    public String toString(String format)
    {
	setSeparators();

        String apString = "";
        if (format.equals("DENOPTIM"))
        {
	    apString = moreAtmsSeparator + Integer.toString(apAtm+1);
            apString = apString + atmSeparator + apRule;
            apString = apString + detailsSeparator + apSubClass;
            apString = apString + detailsSeparator + getAPVecStrg("x") + coordSeparator + getAPVecStrg("y") + coordSeparator + getAPVecStrg("z");
        }
//add other formats 
//TODO
        return apString;
    }

//------------------------------------------------------------------------------

/**
 * TODO
 */

    public String toStringFurther(String format)
    {
        setSeparators();

        String apString = "";
        if (format.equals("DENOPTIM"))
        {
            apString = moreAPSeparator + apRule;
            apString = apString + detailsSeparator + apSubClass;
            apString = apString + detailsSeparator + getAPVecStrg("x") + coordSeparator + getAPVecStrg("y") + coordSeparator + getAPVecStrg("z");
        }
//add other formats 
//TODO
        return apString;
    }

//------------------------------------------------------------------------------

    /**
     * 
     * @param coord direction for which the coordinate has to be returned
     * @return a string repredenting the coordinate in the requested axes
     */

    private String getAPVecStrg(String dir)
    {
	Double dblVal = 0.0;
	if (dir.equals("x"))
	    dblVal = apVector.get(0);
	else if (dir.equals("y"))
	    dblVal = apVector.get(1);
	else if (dir.equals("z"))
	    dblVal = apVector.get(2);
        else {
	    System.err.println("ERROR! Unable to write the coordinates for AP: "+this);
	    System.exit(-1);
	}

        DecimalFormat digits = new DecimalFormat("###.####");
	digits.setMinimumFractionDigits(4);
        String strVal = digits.format(dblVal);

	return strVal;
    }

//------------------------------------------------------------------------------
/**
 * TODO
 */
/*
    public
    {
    }
*/

//------------------------------------------------------------------------------
}
