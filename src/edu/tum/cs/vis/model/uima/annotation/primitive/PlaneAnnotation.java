/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Stefan Profanter - initial API and implementation, Year: 2012
 ******************************************************************************/
package edu.tum.cs.vis.model.uima.annotation.primitive;

import java.awt.Color;
import java.util.HashMap;
import java.util.Map.Entry;

import javax.vecmath.Vector3f;

import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

import processing.core.PGraphics;
import edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation;
import edu.tum.cs.vis.model.util.Vertex;

/**
 * @author Stefan Profanter
 * 
 */
public class PlaneAnnotation extends PrimitiveAnnotation {

	/**
	 * 
	 */
	private static final long	serialVersionUID	= 7758656289829843165L;

	private final Vector3f		planeNormal			= new Vector3f();

	private Vector3f			centroid			= new Vector3f();

	/**
	 * @param annotationColor
	 */
	public PlaneAnnotation() {
		super(new Color(250, 200, 255));
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation#drawAnnotation(processing.core.PGraphics)
	 */
	@Override
	public void drawPrimitiveAnnotation(PGraphics g) {
		g.fill(255, 0, 0);
		g.stroke(255, 255, 255);
		g.strokeWeight(2);

		g.line(centroid.x, centroid.y, centroid.z, centroid.x + planeNormal.x, centroid.y
				+ planeNormal.y, centroid.z + planeNormal.z);

	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation#fitAnnotation()
	 */
	@Override
	public void fitAnnotation() {
		/*
		 * Best fitting plane.
		 * 
		 * Formula for plane: a*x + b*y + c*z + d = 0
		 * d = -(a*x + b*y + c*z)
		 * 
		 * Centroid (average values): x0,y0,z0
		 * 
		 * To minimize: Sum [|a(xi-x0) + b(yi-y0) + c(zi-z0)|^2]
		 * 
		 * 
		 *   v^T  = [a  b  c]
		 *
		 *        --                             --
		 *        | x1 - x0    y1 - y0    z1 - z0 |
		 *        | x2 - x0    y2 - y0    z2 - z0 |
		 *   M =  |    .          .          .    |
		 *        |    .          .          .    |
		 *        |    .          .          .    |
		 *        | xn - x0    yn - y0    zn - z0 |
		 *        --                             --
		 * 
		 */

		HashMap<Vertex, Float> vertices = new HashMap<Vertex, Float>();
		centroid = getVerticesWithWeight(vertices);

		int numberOfPoints = vertices.size();

		SimpleMatrix A = new SimpleMatrix(
				numberOfPoints == 3 ? numberOfPoints + 1 : numberOfPoints, 4);

		int row = 0;

		// Weighted SVD according to
		// http://www.mathworks.com/matlabcentral/newsreader/view_thread/262996
		for (Entry<Vertex, Float> e : vertices.entrySet()) {
			float weight = (float) Math.sqrt(e.getValue());
			Vertex v = e.getKey();

			A.setRow(row++, 0, (v.x - centroid.x) * weight, (v.y - centroid.y) * weight,
					(v.z - centroid.z) * weight, 1);
		}

		if (numberOfPoints == 3) {
			A.setRow(row++, 0, 0, 0, 0);
		}

		@SuppressWarnings("rawtypes")
		SimpleSVD svd = A.svd();

		SimpleMatrix U = svd.getU();
		SimpleMatrix W = svd.getW();
		SimpleMatrix V = svd.getV();

		planeNormal.x = (float) svd.getV().get(0, 3);
		planeNormal.y = (float) svd.getV().get(1, 3);
		planeNormal.z = (float) svd.getV().get(2, 3);
		if (planeNormal.length() == 0) {
			planeNormal.get(mesh.getTriangles().get(0).getNormalVector());
		} else
			planeNormal.normalize();
	}

	/**
	 * @return the centroid
	 */
	public Vector3f getCentroid() {
		return centroid;
	}

	/**
	 * @return the planeNormal
	 */
	public Vector3f getPlaneNormal() {
		return planeNormal;
	}

}
