/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Stefan Profanter - initial API and implementation, Year: 2012
 ******************************************************************************/
package edu.tum.cs.vis.model.uima.annotation;

import java.awt.Color;
import java.util.HashSet;

import processing.core.PGraphics;
import edu.tum.cs.vis.model.Model;
import edu.tum.cs.vis.model.uima.cas.MeshCas;
import edu.tum.cs.vis.model.util.Mesh;
import edu.tum.cs.vis.model.util.Triangle;

/**
 * Base class for all mesh annotations. A mesh annotation is an annotation over multiple triangles
 * or lines.
 * 
 * @author Stefan Profanter
 * @param <S>
 *            Type of mesh annotation
 * 
 */
@SuppressWarnings("rawtypes")
public abstract class MeshAnnotation<S extends MeshAnnotation> extends DrawableAnnotation {

	/**
	 * Class of derived annotation.
	 */
	private final Class<S>		clazz;

	/**
	 * auto generated
	 */
	private static final long	serialVersionUID	= 1222063742441634463L;

	/**
	 * The annotation color for this type of annotation.
	 */
	protected Color				annotationColor;
	/**
	 * A random annotation color. Each annotation gets also a random color.
	 */
	private final Color			randomAnnotationColor;

	/**
	 * Use random color for drawing the annotation
	 */
	protected boolean			useRandomColor		= false;

	/**
	 * Mesh which contains the referenced Triangles for which this annotation stands.
	 */
	protected Mesh				mesh				= new Mesh();

	/**
	 * Parent model of this annotation
	 */
	protected Model				model;

	/**
	 * Default constructor. Sets the annotation color. Each type of annotation should have a
	 * different color.
	 * 
	 * 
	 * @param clazz2
	 *            Class of derived annotation.
	 * @param model
	 *            parent model for this annotation
	 * 
	 * @param annotationColor
	 *            The annotation color for this type of annotation
	 */
	public MeshAnnotation(Class<S> clazz2, Model model, final Color annotationColor) {
		super();
		this.clazz = clazz2;
		randomAnnotationColor = new Color((int) (Math.random() * 255), (int) (Math.random() * 255),
				(int) (Math.random() * 255));
		this.annotationColor = annotationColor;
		this.model = model;
	}

	@Override
	protected void drawAnnotation(PGraphics g) {
		mesh.drawLines(g, getDrawColor());
		mesh.drawTriangles(g, getDrawColor());
	}

	/**
	 * Returns the color for drawing this annotation. If useRandomColor is true, the random color
	 * will be returned.
	 * 
	 * @return randomAnnotationColor if useRandomColor. Otherwise: annotationColor
	 */
	public Color getDrawColor() {
		if (useRandomColor)
			return randomAnnotationColor;

		return annotationColor;
	}

	/**
	 * Get mesh of annotation
	 * 
	 * @return the mesh
	 */
	public Mesh getMesh() {
		return mesh;
	}

	/**
	 * Get all annotations of same type as this which are direct neighbors of this annotation by
	 * getting all annotations where direct neighbor triangles are a member of.
	 * 
	 * @param cas
	 *            main mesh cas
	 * @return Set of found annotations which are the same type of this annotation
	 */
	public HashSet<S> getNeighborAnnotations(MeshCas cas) {
		return getNeighborAnnotations(cas, clazz);
	}

	/**
	 * Get all annotations of given type which are direct neighbors of this annotation by getting
	 * all annotations where direct neighbor triangles are a member of.
	 * 
	 * @param cas
	 *            main cas
	 * @param parClazz
	 *            type of annotations you want
	 * @return Set of all direct neighbor annotations of given type
	 */
	public <T extends MeshAnnotation> HashSet<T> getNeighborAnnotations(MeshCas cas,
			Class<T> parClazz) {
		HashSet<T> annotations = new HashSet<T>();
		for (Triangle t : getMesh().getTriangles()) {
			annotations.addAll(getNeighborAnnotationsForTriangle(cas, parClazz, t));
		}
		if (annotations.contains(this))
			annotations.remove(this);
		return annotations;
	}

	/**
	 * Get all annotations of given type which are direct neighbors of given triangle by getting all
	 * annotations where direct neighbor triangles are a member of.
	 * 
	 * @param cas
	 *            main cas
	 * @param parClazz
	 *            type of annotations you want
	 * @param t
	 *            Get annotations of all direct neighbors of this triangle.
	 * @return Set of found annotations
	 */
	public <T extends MeshAnnotation> HashSet<T> getNeighborAnnotationsForTriangle(MeshCas cas,
			Class<T> parClazz, Triangle t) {

		HashSet<T> annotations = new HashSet<T>();
		/* 
		 * Triangle has max 3 neighbors
		 * First iteration set a1 to annotation of neighbor 1, second set a2 to annotation of neighbor 2.
		 * On third iteration 
		 */

		// Check all neighbors of the triangle which annotation they have
		for (Triangle neig : t.getNeighbors()) {
			// neighbor is in same annotation, skip
			if (getMesh().getTriangles().contains(neig))
				continue;

			// Get annotation of triangle
			T ma = cas.findAnnotation(parClazz, neig);
			if (ma != null)
				annotations.add(ma);
		}
		return annotations;
	}

	/**
	 * Indicates if random color is used to draw this annotation
	 * 
	 * @return the useRandomColor
	 */
	public boolean isUseRandomColor() {
		return useRandomColor;
	}

	/**
	 * Checks if this annotation includes the triangle <code>p</code>.
	 * 
	 * @param p
	 *            triangle to check for
	 * @return true if annotation includes triangle
	 */
	public boolean meshContainsTriangle(final Triangle p) {
		return mesh.getTriangles().contains(p);
	}

	/**
	 * Base mesh of this annotation.
	 * 
	 * @param mesh
	 *            the mesh to set
	 */
	public void setMesh(Mesh mesh) {
		this.mesh = mesh;
	}

	/**
	 * Set to true if random color should be used to draw annotation instead of predefined one.
	 * 
	 * @param useRandomColor
	 *            the useRandomColor to set
	 */
	public void setUseRandomColor(boolean useRandomColor) {
		this.useRandomColor = useRandomColor;
	}

}
