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

/**
 * Toolbox for operations in Cartesian space.
 * 
 * @author Marco Foscato (University of Bergen)
 */

public class CartesianSpaceUtils
{

    //Reporting flag
    private static int repOnScreen = Parameters.report;

//------------------------------------------------------------------------------

    /**
     * Changes the origin of a vector.
     * @param v the original vector
     * @param newOrigin the new origin as a point in Cartesian space
     */
    public static void translateOrigin(Vector3d v, Point3d newOrigin)
    {
	v.x = v.x + newOrigin.x;
        v.y = v.y + newOrigin.y;
        v.z = v.z + newOrigin.z;
    }

//------------------------------------------------------------------------------

    /**
     * Creates an object <code>Vector3d</code> that originates from point
     * <code>a</code> and goes to point <code>b</code>.
     * @param a starting point
     * @param b termination point
     * @return the vector that goes from <code>a</code> to <code>b</code>
     */
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
     * Generate a vector that is perpedicular to the given one.
     * No control over which perpendicular direction will be chosen.
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
     * Rotate a vector according to a given rotation axis and angle.
     * @param v original vector to be rotated
     * @param axis rotation axis
     * @param ang rotation angle
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
