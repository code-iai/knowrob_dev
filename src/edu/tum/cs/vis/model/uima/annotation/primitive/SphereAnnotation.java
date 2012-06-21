/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Stefan Profanter - initial API and implementation, Year: 2012
 ******************************************************************************/
package edu.tum.cs.vis.model.uima.annotation.primitive;

import java.awt.Color;
import java.util.HashSet;

import javax.vecmath.Vector3f;

import processing.core.PGraphics;
import edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation;
import edu.tum.cs.vis.model.util.Vertex;

/**
 * @author Stefan Profanter
 * 
 */
public class SphereAnnotation extends PrimitiveAnnotation {

	/**
	 * 
	 */
	private static final long		serialVersionUID	= 4870579150170536881L;

	private final boolean			concav;

	private final HashSet<Vertex>	vertices			= new HashSet<Vertex>();

	private final Vector3f			center				= new Vector3f();

	private float					radius				= 0;

	public SphereAnnotation(boolean concav) {
		super(concav ? new Color(0, 255, 0) : new Color(255, 0, 0));
		this.concav = concav;
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation#drawAnnotation(processing.core.PGraphics)
	 */
	@Override
	public void drawPrimitiveAnnotation(PGraphics g) {
		g.fill(getDrawColor().getRed(), getDrawColor().getGreen(), getDrawColor().getBlue(), 50);
		g.translate(center.x, center.y, center.z);
		g.sphere(radius);
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation#fitAnnotation()
	 */
	@Override
	public void fitAnnotation() {
		/*
		 * Fitting sphere iteratively according to http://www.geometrictools.com/Documentation/LeastSquaresFitting.pdf
		 */

		Vector3f centroid = getVertices(vertices);

		float a = centroid.x;
		float b = centroid.y;
		float c = centroid.z;

		float a2 = Float.MAX_VALUE, b2 = Float.MAX_VALUE, c2 = Float.MAX_VALUE;

		int iterations = 0;

		while (Math.abs(a - a2) + Math.abs(b - b2) + Math.abs(c - c2) != 0 && iterations < 500) {

			a2 = a;
			b2 = b;
			c2 = c;

			float L = getL(a2, b2, c2);

			a = centroid.x + L * getLa(a2, b2, c2);
			b = centroid.y + L * getLb(a2, b2, c2);
			c = centroid.z + L * getLc(a2, b2, c2);

			iterations++;

		}

		center.x = a;
		center.y = b;
		center.z = c;

		radius = getL(a, b, c);

		// System.out.println("Center: " + center.toString() + " Radius: " + radius);

	}

	/**
	 * @return the center
	 */
	public Vector3f getCenter() {
		return center;
	}

	private float getL(float a, float b, float c) {
		float sum = 0;
		for (Vertex v : vertices) {
			sum += getLi(v, a, b, c);
		}
		return sum / vertices.size();
	}

	private float getLa(float a, float b, float c) {
		float sum = 0;
		for (Vertex v : vertices) {
			sum += (a - v.x) / getLi(v, a, b, c);
		}
		return sum / vertices.size();
	}

	private float getLb(float a, float b, float c) {
		float sum = 0;
		for (Vertex v : vertices) {
			sum += (b - v.y) / getLi(v, a, b, c);
		}
		return sum / vertices.size();
	}

	private float getLc(float a, float b, float c) {
		float sum = 0;
		for (Vertex v : vertices) {
			sum += (c - v.z) / getLi(v, a, b, c);
		}
		return sum / vertices.size();
	}

	private float getLi(Vertex v, float a, float b, float c) {
		return (float) Math
				.sqrt(Math.pow(v.x - a, 2) + Math.pow(v.y - b, 2) + Math.pow(v.z - c, 2));
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation#getPrimitiveArea()
	 */
	@Override
	public float getPrimitiveArea() {
		return (float) (4f * Math.PI * (radius * radius));
	}

	/**
	 * @return the radius
	 */
	public float getRadius() {
		return radius;
	}

	private float getResiduum(float a, float b, float c) {
		float sum = 0;
		float r = getL(a, b, c);
		for (Vertex v : vertices) {
			sum += getLi(v, a, b, c) - r;
		}
		return Math.abs(sum / vertices.size());
	}

	/**
	 * @return the concav
	 */
	public boolean isConcav() {
		return concav;
	}

}
