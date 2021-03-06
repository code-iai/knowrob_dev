/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Stefan Profanter - initial API and implementation, Year: 2012
 ******************************************************************************/
package org.knowrob.vis.model.util.algorithm;

import java.util.ArrayList;

import javax.vecmath.Point3f;
import javax.vecmath.Vector3f;

import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

/**
 * Best fits a line to given points in 3D.
 * 
 * @author Stefan Profanter
 * 
 */
public class BestFitLine3D {

	/**
	 * Best fits a line to given points <tt>points</tt>. New line is represented by <tt>center</tt>
	 * and <tt>direction</tt>
	 * 
	 * @param points
	 *            List of points to fit line
	 * @param dir
	 *            normalized direction of best fit line
	 * @param center
	 *            center of best fit line
	 */
	public static void getBestFitLine(ArrayList<Point3f> points, Vector3f dir, Point3f center) {
		getBestFitLine(points, dir, null, center);
	}

	/**
	 * Best fits a line to given points <tt>points</tt>. New line is represented by <tt>center</tt>
	 * and <tt>direction</tt>. Additionally you get perpendicular direction vector for perpendicular
	 * line.
	 * 
	 * If the number of points is 0, the function returns immediately without changing anything.
	 * 
	 * @param points
	 *            List of points to fit line
	 * @param dir
	 *            normalized direction of best fit line
	 * @param dirPerp
	 *            normalized perpendicular direction to best fit line *
	 * @param center
	 *            center of best fit line
	 */
	public static void getBestFitLine(ArrayList<Point3f> points, Vector3f dir, Vector3f dirPerp,
			Point3f center) {
		// Try to best fit a line through the intersection points:

		if (points.size() == 0) {
			return;
		}

		// Number of maximum rows for SVD. If value is too large SVD will throw a Out of Heap
		// Exception
		int maxRows = 1000;

		// Calculate the increment value if number of intersections is larger than maxRows.
		int increment = (points.size() / maxRows);
		if (increment == 0)
			increment = 1;
		int numberOfPoints = points.size() / increment;
		if (points.size() % maxRows > 0 && increment != 1)
			numberOfPoints++;

		SimpleMatrix A = new SimpleMatrix(
				numberOfPoints == 2 ? numberOfPoints + 1 : numberOfPoints, 3);

		int row = 0;

		center.x = 0;
		center.y = 0;
		center.z = 0;

		for (int i = 0; i < points.size(); i += increment) {
			center.add(points.get(i));
		}
		center.scale(1f / points.size());

		for (int i = 0; i < points.size(); i += increment) {
			Point3f p = points.get(i);
			A.setRow(row++, 0, (p.x - center.x), (p.y - center.y), (p.z - center.z));
		}

		if (numberOfPoints == 2) {
			A.setRow(row++, 0, center.x, center.y, center.z);
		}

		@SuppressWarnings("rawtypes")
		SimpleSVD svd = A.svd();

		if (dirPerp != null) {
			dirPerp.x = (float) svd.getV().get(0, 2);
			dirPerp.y = (float) svd.getV().get(1, 2);
			dirPerp.z = (float) svd.getV().get(2, 2);
		}
		// Take the column with biggest singular value (it is the first one).
		dir.x = (float) svd.getV().get(0, 0);
		dir.y = (float) svd.getV().get(1, 0);
		dir.z = (float) svd.getV().get(2, 0);

	}
}
