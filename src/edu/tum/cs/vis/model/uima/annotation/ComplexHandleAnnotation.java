/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Stefan Profanter - initial API and implementation, Year: 2012
 ******************************************************************************/
package edu.tum.cs.vis.model.uima.annotation;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;

import javax.vecmath.Matrix4f;
import javax.vecmath.Tuple3f;
import javax.vecmath.Vector3f;

import processing.core.PGraphics;
import edu.tum.cs.vis.model.Model;
import edu.tum.cs.vis.model.uima.annotation.primitive.Cone;
import edu.tum.cs.vis.model.uima.annotation.primitive.PlaneAnnotation;
import edu.tum.cs.vis.model.util.DrawSettings;
import edu.tum.cs.vis.model.util.Triangle;
import edu.tum.cs.vis.model.util.Vertex;

/**
 * Annotation for a complex handle.
 * 
 * A complex handle is a convex part in the model.
 * 
 * 
 * @author Stefan Profanter
 * 
 */
public class ComplexHandleAnnotation extends DrawableAnnotation implements HandleAnnotation,
		PrologAnnotationInterface {

	/**
	 * generated
	 */
	private static final long				serialVersionUID		= 1983829446921229660L;

	/**
	 * Set of primitive annotations which represent the complex handle.
	 */
	@SuppressWarnings("rawtypes")
	private final Set<PrimitiveAnnotation>	primitiveAnnotations	= new HashSet<PrimitiveAnnotation>();

	/**
	 * Fitted cone of complex handle
	 */
	private final Cone						cone;

	/**
	 * The parent model of complex handle.
	 */
	private final Model						model;
	/**
	 * allowed tolerance between normal vectors. Tolerance is indicated as the result of the dot
	 * product. The resulting angle is PI-acos(DOT_NORMAL_TOLERANCE), which is for 0.1 approx. 6
	 * degree.
	 */
	public final static float				DOT_NORMAL_TOLERANCE	= 0.1f;

	/**
	 * Create new complex handle annotation.
	 * 
	 * @param model
	 *            parent model
	 */
	public ComplexHandleAnnotation(Model model) {
		super(new Color(255, 41, 255, 200));
		this.model = model;
		cone = new Cone(false);
	}

	/**
	 * Add given annotation to the complex handle. After all annotations are added, you have to call
	 * fit afterwards to refit cone.
	 * 
	 * @param pa
	 *            annotation to add.
	 */
	public void addAnnotation(@SuppressWarnings("rawtypes") PrimitiveAnnotation pa) {
		primitiveAnnotations.add(pa);

	}

	/**
	 * Check if the neighborHandle violates constraints for complex handle.
	 * 
	 * @param neighborHandle
	 *            handle to check
	 * @return true if neighborHandle combined with this handle doesn't violate handle constraint
	 */
	@SuppressWarnings("rawtypes")
	public boolean allowMerge(ComplexHandleAnnotation neighborHandle) {
		for (PrimitiveAnnotation pa : primitiveAnnotations) {
			if (!(pa instanceof PlaneAnnotation))
				continue;
			PlaneAnnotation plane1 = (PlaneAnnotation) pa;
			for (PrimitiveAnnotation pn : neighborHandle.primitiveAnnotations) {
				if (!(pn instanceof PlaneAnnotation))
					continue;
				PlaneAnnotation plane2 = (PlaneAnnotation) pn;

				// we have two planes. Now check if they violate the angle constraint. (the angle
				// between two planes must be bigger or equal to 180 degree).
				float dot = plane1.getPlaneNormal().dot(plane2.getPlaneNormal());
				// allow a small error
				if (dot > DOT_NORMAL_TOLERANCE)
					return false;
			}
		}
		return true;

	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.DrawableAnnotation#containsTriangle(edu.tum.cs.vis.model.util.Triangle)
	 */
	@Override
	public boolean containsTriangle(Triangle t) {
		for (@SuppressWarnings("rawtypes")
		PrimitiveAnnotation pa : primitiveAnnotations) {
			if (pa.containsTriangle(t))
				return true;
		}
		return false;
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.DrawableAnnotation#drawAnnotation(processing.core.PGraphics, edu.tum.cs.vis.model.util.DrawSettings)
	 */
	@Override
	protected void drawAnnotation(PGraphics g, DrawSettings drawSettings) {
		DrawSettings tmp = (DrawSettings) drawSettings.clone();
		if (tmp.getOverrideColor() == null)
			tmp = drawSettings.getTemporaryOverride(getDrawColor());
		tmp.forceDraw = true;
		for (@SuppressWarnings("rawtypes")
		PrimitiveAnnotation pa : primitiveAnnotations) {
			pa.draw(g, tmp);
		}
	}

	/**
	 * Draw fitted cone of complex handle.
	 * 
	 * @param g
	 *            Graphics context
	 * @param color
	 *            color to draw. If null, default handle color is used.
	 */
	public void drawPrimitiveAnnotation(PGraphics g, Color color) {

		// Creating new color to set alpha to 255
		Color c = color == null ? getDrawColor() : color;
		cone.draw(g, new Color(c.getRed(), c.getGreen(), c.getBlue(), 255));

	}

	/**
	 * Fit cone to annotation.
	 * 
	 * @return true if successfully fit.
	 */
	@SuppressWarnings("unchecked")
	public boolean fit() {

		LinkedHashMap<Vertex, Float> vertices = new LinkedHashMap<Vertex, Float>();
		Vector3f centroid = new Vector3f();
		Set<Triangle> triangles = new HashSet<Triangle>();

		for (@SuppressWarnings("rawtypes")
		PrimitiveAnnotation pa : primitiveAnnotations) {
			centroid.add(pa.getVerticesWithWeight(vertices));
			triangles.addAll(pa.mesh.getTriangles());
		}
		centroid.scale(1f / primitiveAnnotations.size());
		return cone.fit(centroid, vertices.keySet(), vertices, triangles);
	}

	/**
	 * Gets the area of all triangles.
	 * 
	 * @return Area of cone.
	 */
	@Override
	public float getArea() {
		float trianglesArea = 0;
		for (@SuppressWarnings("rawtypes")
		PrimitiveAnnotation pa : primitiveAnnotations) {
			trianglesArea += pa.getArea();
		}
		return model.getUnscaled(trianglesArea);
	}

	/**
	 * Get value between > 0 for area coverage which indicates how good primitive annotation is fit
	 * into mesh. 1 indicates perfect fit, because area of triangles is exactly the same as area of
	 * primitive annotation.
	 * 
	 * @return value > 0
	 */
	@Override
	public float getAreaCoverage() {
		return getPrimitiveArea() / getArea();
	}

	/**
	 * get centroid of cone at unscaled position
	 * 
	 * @return the centroid
	 */
	public Tuple3f getCentroid() {
		return model.getUnscaled(cone.getCentroid());
	}

	@Override
	public Cone getCone() {
		return cone;
	}

	/**
	 * 
	 * get direction (unscaled) of cone. Direction is aligned with generating line and shows from
	 * centroid to small radius cap. Length of direction is half height of the cone (center to one
	 * end).
	 * 
	 * @return the direction
	 */
	public Vector3f getDirection() {
		return new Vector3f(model.getUnscaled(cone.getDirection()));
	}

	/**
	 * Get total unscaled height of cone from bottom cap to top cap which is 2*directionUnscaled.
	 * 
	 * @return unscaled height
	 */
	public float getHeight() {
		return getDirection().length() * 2;
	}

	/**
	 * Get parent model of annotation.
	 * 
	 * @return the model
	 */
	public Model getModel() {
		return model;
	}

	/**
	 * Get pose matrix for cone.
	 * 
	 * @return 4x4 pose matrix of the plane relative to the object centroid
	 */
	@Override
	public Matrix4f getPoseMatrix() {

		return cone.getPoseMatrix();
	}

	/**
	 * Get set of primitive annotations which are part of this complex handle annotation.
	 * 
	 * @return set of primitive annotations representing this handle
	 */
	@SuppressWarnings("rawtypes")
	public Set<PrimitiveAnnotation> getPrimitiveAnnotations() {
		return primitiveAnnotations;
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrologBaseAnnotation#getPrimitiveArea()
	 */
	@Override
	public float getPrimitiveArea() {
		return model.getUnscaled(cone.getArea());
	}

	/**
	 * Get average radius (unscaled) of cone which is the average between small and large radius
	 * 
	 * @return average radius unscaled
	 */
	public float getRadiusAvg() {
		return model.getUnscaled(cone.getRadiusAvg());
	}

	/**
	 * Get large radius (unscaled), which is at the bottom of cone.
	 * 
	 * @return the radiusLarge
	 */
	public float getRadiusLarge() {
		return model.getUnscaled(cone.getRadiusLarge());
	}

	/**
	 * Get small radius (unscaled), which is at the bottom of cone.
	 * 
	 * @return the radiusSmall
	 */
	public float getRadiusSmall() {
		return model.getUnscaled(cone.getRadiusSmall());
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrologBaseAnnotation#getTriangles()
	 */
	@Override
	public Triangle[] getTriangles() {
		List<Triangle> t = new ArrayList<Triangle>();
		for (@SuppressWarnings("rawtypes")
		PrimitiveAnnotation a : primitiveAnnotations)
			t.addAll(Arrays.asList(a.getTriangles()));
		return t.toArray(new Triangle[0]);
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrologBaseAnnotation#getVertices()
	 */
	@Override
	public Vertex[] getVertices() {
		Set<Vertex> v = new HashSet<Vertex>();
		for (@SuppressWarnings("rawtypes")
		PrimitiveAnnotation a : primitiveAnnotations)
			v.addAll(Arrays.asList(a.getVertices()));
		return v.toArray(new Vertex[0]);
	}

	/**
	 * Get unscaled volume of cone.
	 * 
	 * @return the volume unscaled.
	 */
	public float getVolume() {

		return model.getUnscaled(cone.getVolume());
	}

	/**
	 * Is cone concave or convex?
	 * 
	 * @return the concave
	 */
	public boolean isConcave() {
		return cone.isConcave();
	}

	/**
	 * Merge given complex handle into this one. You have to call fit to refit the cone.
	 * 
	 * @param an
	 *            annotation to merge into this.
	 */
	public void merge(ComplexHandleAnnotation an) {
		primitiveAnnotations.addAll(an.primitiveAnnotations);
	}
}
