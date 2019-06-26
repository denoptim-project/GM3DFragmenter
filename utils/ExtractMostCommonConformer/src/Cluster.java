import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;

/**
 *
 */

public class Cluster
{
    //Cluster ID
    private String id;
 
    //Number of members included in this cluster
    private int numMembers;

    //List of members as String 'V1, V22, V43, ...'
    private ArrayList<String> membersStr = new ArrayList<String>();

    //List of members as Integer '1, 22, 43, ...'
    private ArrayList<Integer> membersInt = new ArrayList<Integer>();

    //Labeb
    private final String label = "V";

//------------------------------------------------------------------------------

    public Cluster()
    {
    }

//------------------------------------------------------------------------------

    public Cluster(String id)
    {
	this.id = id;
    }

//------------------------------------------------------------------------------

    public Cluster(String id, int n, ArrayList<String> membersStr, ArrayList<Integer> membersInt)
    {
	this.id = id;
	this.numMembers = n;
	this.membersStr = membersStr;
	this.membersInt= membersInt;
    }

//------------------------------------------------------------------------------

    public void addMember(String str)
    {
	membersStr.add(str);
	int num = Integer.parseInt(str.substring(label.length(),str.length()));
	membersInt.add(num);
	numMembers = membersStr.size();
    }

//------------------------------------------------------------------------------

    public void addMember(int num)
    {
	membersInt.add(num);
	String str = label + Integer.toString(num);
	membersStr.add(str);
        numMembers = membersStr.size();
    }

//------------------------------------------------------------------------------

    public int getSize()
    {
	return numMembers;
    }

//------------------------------------------------------------------------------

    public ArrayList<Integer> getMembersIDInt()
    {
	return membersInt;
    }

//------------------------------------------------------------------------------

    public ArrayList<String> getMembersIDStr()
    {
        return membersStr;
    }

//------------------------------------------------------------------------------

    public int getMemberIntIDInPlace(int i)
    {
	if (i > numMembers)
	{
	    System.out.println("ERROR! Trying to get "+i+"-th member of the cluster that hs only "+numMembers+" members!");
	     System.exit(0);
	}

        return membersInt.get(i);
    }

//------------------------------------------------------------------------------

    public String  getMemberStrIDInPlace(int i)
    {
        if (i > numMembers)
        {
            System.out.println("ERROR! Trying to get "+i+"-th member of the cluster that hs only "+numMembers+" members!");
             System.exit(0);
        }

        return membersStr.get(i);
    }


//------------------------------------------------------------------------------
    public String toString()
    {
	String s = "[Cluster " + id +
		 " NumEntries: " + numMembers +
		 " EntriesStr: " + membersStr + 
		 " EntriesInt: " + membersInt + "]";
	return s;
    }
}
