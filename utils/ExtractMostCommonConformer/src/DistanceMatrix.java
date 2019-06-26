import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;

/**
 *
 */

public class DistanceMatrix
{

    //List of elements' reperences
    private ArrayList<Integer> labels;

    //Matrix elements
    private double[][] mtx;

    //Size
    private int size;

//------------------------------------------------------------------------------

    public DistanceMatrix(int size)
    {
	this.labels = new ArrayList<Integer>();
	this.mtx = new double[size][size];
        for (int i=0; i < size; i++)
	{
            this.mtx[i][i] = 0;
	    this.labels.add(i);
	}
	this.size = size;
    }

//------------------------------------------------------------------------------

    public DistanceMatrix(ArrayList<Integer> labels)
    {
        this.labels = labels;
        this.size = labels.size();
        this.mtx = new double[labels.size()][labels.size()];
    }

//------------------------------------------------------------------------------

    public DistanceMatrix(ArrayList<Integer> labels, double[][] mtx)
    {
	this.mtx = mtx;
	this.labels = labels;
	int sizeLab = labels.size();
	int sizeArr = mtx.length;
	if (sizeLab != sizeArr)
	{
	    System.err.println("ERROR! Attempt to create a distance matrix with non-same number of labels and entries!");
	    System.exit(0);
	} else {
	    this.size = sizeLab;
	}
    }

//------------------------------------------------------------------------------

    public void setElement(int i, int j, double value)
    {
	this.mtx[i][j] = value;
    }

//------------------------------------------------------------------------------

    public double getElement(int i, int j)
    {
        return mtx[i][j];
    }

//------------------------------------------------------------------------------

    public double getMaxValue()
    {
	double max = Double.MIN_VALUE;
	for (int i=0; i<size; i++)
	{
	    for (int j=0; j<size; j++)
	    {
		double value = mtx[i][j];
		if (value > max)
		{
		    max = value;
		} 
	    }
	}
	return max;
    }

//------------------------------------------------------------------------------

    public double getMinValue()
    {
        double min = Double.MAX_VALUE;
        for (int i=0; i<size; i++)
        {
            for (int j=0; j<size; j++)
            {
                double value = mtx[i][j];
                if (value < min)
                {
                    min = value;
                }
            }
        }
        return min;
    }

//------------------------------------------------------------------------------

    public int[] getMaxValueEntry()
    {
        double max = Double.MIN_VALUE;
	int[] maxIdx = new int[2];
        for (int i=0; i<size; i++)
        {
            for (int j=0; j<size; j++)
            {
                double value = mtx[i][j];
                if (value > max)
                {
                    max = value;
		    maxIdx[0] = i;
		    maxIdx[1] = j;
                }
            }
        }
        return maxIdx;
    }

//------------------------------------------------------------------------------

    public int[] getMinValueEntry()
    {
        double min = Double.MAX_VALUE;
        int[] minIdx = new int[2];
        for (int i=0; i<size; i++)
        {
            for (int j=0; j<size; j++)
            {
                double value = mtx[i][j];
                if (value < min)
                {
                    min = value;
                    minIdx[0] = i;
                    minIdx[1] = j;
                }
            }
        }
        return minIdx;
    }

//------------------------------------------------------------------------------

    /**
     * Returns the number (0-n) of the most distant vector to all the 
     * others, or 0 if all verctors are equally distant to each other.
     */

    public int getMaxSumColumnVector()
    {
        double max = - Double.MIN_VALUE;
        int maxSumIdx = 0;
        for (int i=0; i<size; i++)
        {
            double summ = - Double.MIN_VALUE;
            for (int j=0; j<size; j++)
            {
                double value = mtx[i][j];
                summ = summ + value;
            }

            if (summ > max)
            {
                max = summ;
                maxSumIdx = i;
            }
        }
        return maxSumIdx;
    }

//------------------------------------------------------------------------------

    /**
     * Returns the number (0-n) of the vector closer to all the 
     * others, or 0 if all verctors are equally distant to each other.
     */

    public int getMinSumColumnVector()
    {
        double min = Double.MIN_VALUE;
	int minSumIdx = -1;
        for (int i=0; i<size; i++)
        {
	    double summ = 0.000;
            for (int j=0; j<size; j++)
            {
                double value = mtx[i][j];
                summ = summ + value;
            }

	    if (summ < min)
	    {
		min = summ;
		minSumIdx = i;
	    }
        }
        return minSumIdx;
    }

//------------------------------------------------------------------------------

    public int getSize()
    {
	return size;
    }

//------------------------------------------------------------------------------

    public ArrayList<Integer> getLabels()
    {
	return labels;
    }

//------------------------------------------------------------------------------

    public int getLabel(int n)
    {
        return labels.get(n);
    }

//------------------------------------------------------------------------------

    public void pruneMostDistant()
    {
	//identify the "bad guy"
	int mostdistant = getMaxSumColumnVector();

	//prepare the new matrix
	int newSize = this.size - 1;
	double[][] newMtx = new double[newSize][newSize];

	//copy all entries but the most distance
        for (int i=0; i<this.size; i++)
        {
            for (int j=0; j<this.size; j++)
            {
		if ((i < mostdistant) && (j < mostdistant))
		{
                    newMtx[i][j] = mtx[i][j];
		} else if ((i > mostdistant) && (j > mostdistant)) 
		{
		    newMtx[i-1][j-1] = mtx[i][j];
		} else if ((i > mostdistant) && (j < mostdistant))
                {
                    newMtx[i-1][j] = mtx[i][j];
                } else if ((i < mostdistant) && (j > mostdistant))
                {
                    newMtx[i][j-1] = mtx[i][j];
                } else if ((i == mostdistant) || (j == mostdistant))
		{
		    continue;
		} else {
		    System.out.println("ERROR! Case i: "+i+" j: "+j+" not covered!");
                    System.out.println("Check code! This shouldn't happen.");
		    System.exit(0);
		}
            }
        }

	//Update list of labels
	ArrayList<Integer> newLabels = new ArrayList<Integer>();
	for (int i=0; i < this.size; i++)
        {
	    if (i != mostdistant)
		newLabels.add(labels.get(i));
	}

        //Update matrix
        this.mtx = newMtx;

        //Update labels
        this.labels = newLabels;

	//Update size
	this.size = this.size - 1;
    }

//------------------------------------------------------------------------------

    public String toStringUpperTriangle()
    {
        String s = "";
        for (int i=0; i<size; i++)
        {
            String line = "";
            for (int j=0; j<size; j++)
            {
		if (j>j)
                    line = line + String.format("%.2f ", this.mtx[i][j]);
		else
		    line = line + "0.00 ";
            }
            s = s + line + "\n";
        }

        return s;
    }

//------------------------------------------------------------------------------

    public String toString()
    {
	String s = "[DistanceMatrix " + 
		 " size: " + size +
		 " labels: " + labels + 
		 " mtx: \n";
        for (int i=0; i<size; i++)
        {
            String line = "";
            for (int j=0; j<size; j++)
            {
                line = line + String.format("%.2f ", this.mtx[i][j]);
            }
            s = s + line + "\n";
        }
	
	return s;
    }

//------------------------------------------------------------------------------

}
