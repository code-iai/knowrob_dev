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

import javax.vecmath.Matrix4f;
import javax.vecmath.Tuple3f;
import javax.vecmath.Vector3f;

import processing.core.PGraphics;
import edu.tum.cs.vis.model.Model;
import edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation;
import edu.tum.cs.vis.model.util.Curvature;
import edu.tum.cs.vis.model.util.Vertex;

/**
 * @author Stefan Profanter
 * 
 */
public final class PlaneAnnotation extends PrimitiveAnnotation<PlaneAnnotation> {

	/**
	 * auto generated
	 */
	private static final long	serialVersionUID	= 7758656289829843165L;

	/**
	 * Plane representing the plane annotation
	 */
	private final Plane			plane;

	/**
	 * Creates a new plane annotation.
	 * 
	 * @param curvatures
	 *            Map of curvatures for vertices
	 * @param model
	 *            parent model
	 */
	public PlaneAnnotation(HashMap<Vertex, Curvature> curvatures, Model model) {
		super(PlaneAnnotation.class, curvatures, model, new Color(250, 200, 255));
		plane = new Plane();
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation#drawAnnotation(processing.core.PGraphics)
	 */
	@Override
	public void drawPrimitiveAnnotation(PGraphics g, Color color) {
		plane.draw(g, color == null ? getDrawColor() : color);

	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation#fitAnnotation()
	 */
	@Override
	public boolean fitAnnotation() {
		HashMap<Vertex, Float> vertices = new HashMap<Vertex, Float>();
		Vector3f centroid = getVerticesWithWeight(vertices);

		return plane.fit(centroid, vertices.keySet(), vertices, mesh.getTriangles());
	}

	/**
	 * Get centroid of plane at unscaled position
	 * 
	 * @return the centroid
	 */
	public Tuple3f getCentroid() {
		return model.getUnscaled(plane.getCentroid());
	}

	/**
	 * Get unsaled plane corners
	 * 
	 * @return the corner
	 */
	public Tuple3f[] getCorner() {

		return model.getUnscaled(plane.getCorner());
	}

	/**
	 * Get unscaled vector (length, direction) of long side.
	 * 
	 * @return the longSide
	 */
	public Vector3f getLongSide() {
		return model.getUnscaled(plane.getLongSide());
	}

	/**
	 * Get the unerlying plane.
	 * 
	 * WARNING: Sphere parameters are scaled. To get unscaled values put them through
	 * model.getUnscaled().
	 * 
	 * @return fitted plane for annotation
	 */
	public Plane getPlane() {
		return plane;
	}

	/**
	 * Get plane normal
	 * 
	 * @return the planeNormal
	 */
	public Vector3f getPlaneNormal() {
		return plane.getPlaneNormal();
	}

	/**
	 * Get pose matrix for plane.
	 * 
	 * @return 4x4 pose matrix of the plane relative to the object centroid
	 */
	@Override
	public Matrix4f getPoseMatrix() {
		return plane.getPoseMatrix();
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation#getPrimitiveArea()
	 */
	@Override
	public float getPrimitiveArea() {
		return model.getUnscaled(plane.getArea());
	}

	/**
	 * Get unscaled vector (length, direction) of long side.
	 * 
	 * @return the shortSide
	 */
	public Vector3f getShortSide() {
		return model.getUnscaled(plane.getShortSide());
	}

}