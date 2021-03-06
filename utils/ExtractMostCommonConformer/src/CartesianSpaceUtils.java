/**
 * Toolbox for Operations in Cartesian space
 * 
 * @author Marco Foscato (University of Bergen)
 */

import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.Arrays;
import java.util.ArrayList;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;
import javax.vecmath.Matrix3d;
import javax.vecmath.Quat4d;
import javax.vecmath.AxisAngle4d;

public class CartesianSpaceUtils
{

    //Reporting flag
    private static int repOnScreen = 0;

//------------------------------------------------------------------------------

    public static void translateOrigin(Vector3d v, Point3d newOrigin)
    {
	v.x = v.x + newOrigin.x;
        v.y = v.y + newOrigin.y;
        v.z = v.z + newOrigin.z;
    }

//------------------------------------------------------------------------------

    public static Quat4d getQuaternion(Vector3d v1, Vector3d v2)
    {
        Quat4d q = new Quat4d();
        System.out.println("getQuaternion!!!");
//TODO tobe written!
IOtools.pause();

        return q;
    }

//------------------------------------------------------------------------------

    public static Vector3d getVectorFromTo(Point3d a, Point3d b)
    {
        double x = b.x - a.x;
        double y = b.y - a.y;
        double z = b.z - a.z;

        Vector3d v = new Vector3d(x, y, z);

        return v;        
    }

//------------------------------------------------------------------------------

    /**
     * Get sum of vector A and B
     * @param A vector A
     * @param B vector B
     * @return the vector sum
     */
    public static Vector3d getSumOfVector(Vector3d A, Vector3d B)
    {
	Vector3d resVec = new Vector3d((A.x + B.x),
				   (A.y + B.y),
				   (A.z + B.z));
	return resVec;
    }

//------------------------------------------------------------------------------

    /**
     * Calculates the vector difference of vectors A and B
     * @param A vector A
     * @param B vector B
     * @return the vector difference
     */
    public static Vector3d getDiffOfVector(Vector3d A, Vector3d B)
    {
        Vector3d resVec = new Vector3d((A.x - B.x),
                                   (A.y - B.y),
                                   (A.z - B.z));
        return resVec;
    }

//------------------------------------------------------------------------------

    /**
     * Generate a vector that is perpedicular to the given one
     * No control on which perpendicular direction will be chosen
     * @param dir input direction
     * @return a perpendicular/normal direction
     */

    public static Vector3d getNormalDirection(Vector3d dir)
    {
        Vector3d normalDir = new Vector3d();

        Vector3d dirX = new Vector3d(1.0, 0.0, 0.0);
        Vector3d dirY = new Vector3d(0.0, 1.0, 0.0);
        Vector3d dirZ = new Vector3d(0.0, 0.0, 1.0);
        List<Vector3d> candidates = new ArrayList<Vector3d>();
        candidates.add(dirX);
        candidates.add(dirY);
        candidates.add(dirZ);

        // Check for the lucky case... one of the candidates IS the solution
        List<Double> dotProds = new ArrayList<Double>();
        boolean found = false;
        double max = 0.0;
        for (int i=0; i<candidates.size(); i++)
        {
            double res = dir.dot(candidates.get(i));
	    double absRes = Math.abs(res);
//TODO max or min? possible bug change to <
            if (absRes > max)
                max = absRes;

            if (res == 0.0)
            {
                normalDir = candidates.get(i);
                found = true;
                break;
            } else {
                dotProds.add(absRes);
            }
        }

        // So, since you are not that lucky use the cross-product to get a
        // normal direction using the most divergent of the previous candidates
        if (!found)
        {
            int mostDivergent = dotProds.indexOf(max);
            normalDir.cross(dir,candidates.get(mostDivergent));
            normalDir.normalize();
	}

        return normalDir;
    }

//------------------------------------------------------------------------------

    /**
     * Generate a vector that is perpedicular to the given one
     * and with a given angle to the most divergent reference axis
     * @param dir input direction
     * @param ang ange with the most divergent axis
     * @return the perpendicular/normal direction
     */
/*
TODO DELETE never completed
    public static Vector3d getNormalDirection(Vector3d dir, double ang)
    {
        Vector3d normalDir = new Vector3d();

	dir.normalize();

	// Candidates
        Vector3d dirX = new Vector3d(1.0, 0.0, 0.0);
        Vector3d dirY = new Vector3d(0.0, 1.0, 0.0);
        Vector3d dirZ = new Vector3d(0.0, 0.0, 1.0);
        List<Vector3d> candidates = new ArrayList<Vector3d>();
        candidates.add(dirX);
        candidates.add(dirY);
        candidates.add(dirZ);

	// Find the less and most divergent candidates
	double min = 1.0;
        double max = 0.0;
	int lessDivergent = -1;
	int mostDivergent = -1;
        for (int i=0; i<candidates.size(); i++)
        {
            double res = dir.dot(candidates.get(i));
            double absRes = Math.abs(res);

            if (absRes > max)
	    {
                max = absRes;
		lessDivergent = i;
	    }
            if (absRes < min)
            {
                min = absRes;
                mostDivergent = i;
            }
        }

	// Get rotated vector
	int planeA = mostDivergent;
	int planeB = -1;
	for (int tryOther=0; tryOther<2 ; tryOther++)
	if ((lessDivergent != tryOther) && (mostDivergent != tryOther))
	    planeB = tryOther;
	Vector3d rotVec = getRotatedVector(candidates.get(mostDivergent),
						candidates.get(planeA),
						candidates.get(planeB),
						,ang);
	
	// Get Normal direction from the rotated vector
	normalDir.cross(dir,rotVec);
        normalDir.normalize();

        return normalDir;
    }
*/
//------------------------------------------------------------------------------

    /**
     * Rotate a vector according to axis and angle
     */
    public static void rotatedVectorWAxisAngle(Vector3d v, Vector3d axis, double ang)
    {
	axis.normalize();
	double rad = Math.toRadians(ang);
//System.out.println("Ang: "+ang+" toRad: "+rad);
	AxisAngle4d aa = new AxisAngle4d(axis.x,axis.y,axis.z,rad);
	Matrix3d rotMatrix = new Matrix3d();
	rotMatrix.set(aa);
	rotMatrix.transform(v);	
    }
//------------------------------------------------------------------------------

}
