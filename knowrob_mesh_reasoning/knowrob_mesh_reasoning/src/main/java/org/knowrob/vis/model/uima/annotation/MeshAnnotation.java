/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: 
 * 				Stefan Profanter - initial API and implementation, Year: 2012
 * 				Andrei Stoica - refactored implementation during Google Summer of Code 2014
 ******************************************************************************/
package org.knowrob.vis.model.uima.annotation;

import java.awt.Color;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import processing.core.PGraphics;

import com.google.common.collect.HashMultimap;

import org.knowrob.vis.model.Model;
import org.knowrob.vis.model.uima.cas.MeshCas;
import org.knowrob.vis.model.util.DrawSettings;
import org.knowrob.vis.model.util.Edge;
import org.knowrob.vis.model.util.Mesh;
import org.knowrob.vis.model.util.Triangle;
import org.knowrob.vis.model.util.Vertex;

/**
 * Base class for all mesh annotations. A mesh annotation is an annotation over multiple triangles
 * or lines.
 * 
 * @author Stefan Profanter
 * @author Andrei Stoica (refactor on neighboring annotations getters)
 * 
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
	 * Mesh which contains the referenced triangles for which this annotation stands.
	 */
	protected Mesh				mesh				= new Mesh();

	/**
	 * Parent model of this annotation
	 */
	protected Model				model;

	/**
	 * Default constructor. 
	 * Sets the annotation color. Each type of annotation should have a different color.
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
		super(annotationColor);
		this.clazz = clazz2;
		this.model = model;
	}

	/**
	 * Checks if this annotation includes the triangle {@code p}.
	 * 
	 * @param p
	 *            triangle to check for
	 * @return true if annotation includes triangle
	 */
	@Override
	public boolean containsTriangle(final Triangle p) {
		return mesh.getTriangles().contains(p);
	}

	/**
	 * Draws the selected annotation as called from the visualization
	 * thread.
	 */
	@Override
	protected void drawAnnotation(PGraphics g, DrawSettings drawSettings) {
		DrawSettings tmpSet;
		if (drawSettings.getOverrideColor() != null)
			tmpSet = drawSettings;
		else
			tmpSet = drawSettings.getTemporaryOverride(getDrawColor());
		mesh.drawLines(g, tmpSet);
		mesh.drawTriangles(g, tmpSet);
	}

	/**
	 * Gets the mesh which the annotation belongs to
	 * 
	 * @return the mesh
	 */
	public Mesh getMesh() {
		return mesh;
	}

	/**
	 * Gets the parent model of annotation.
	 * 
	 * @return the model
	 */
	public Model getModel() {
		return model;
	}

	/**
	 * Get all annotations of same type as this which are direct neighbors of this annotation by
	 * getting all annotations where direct neighbor triangles are a member of.
	 * 
	 * @param cas
	 *            main cas
	 * @return Set of found annotations which are the same type of this annotation
	 */
	public Set<S> getNeighborAnnotations(MeshCas cas) {
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
	public <T extends MeshAnnotation> Set<T> getNeighborAnnotations(MeshCas cas, Class<T> parClazz) {
		Set<T> annotations = new HashSet<T>();
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
	public <T extends MeshAnnotation> Set<T> getNeighborAnnotationsForTriangle(MeshCas cas,
			Class<T> parClazz, Triangle t) {

		HashSet<T> annotations = new HashSet<T>();
		/* 
		 * Triangle has max 3 neighbors
		 * First iteration set a1 to annotation of neighbor 1, second set a2 to annotation of neighbor 2.
		 * On third iteration 
		 */

		// Check all neighbors of the triangle which annotation they have
		Edge[] edges = t.getEdges();
		for (int i = 0 ; i < edges.length ; ++i) {
			List<Triangle> neighbors = t.getNeighborsOfEdge(edges[i]);
			for (int j = 0 ; j < neighbors.size() ; ++j) {
				// neighbor not in same annotation, then process
				if (!mesh.getTriangles().contains(neighbors.get(j))) {
					// get annotation of triangle
					T neighAnnotation = cas.findAnnotation(parClazz, neighbors.get(j));
					if (neighAnnotation != null) {
						annotations.add(neighAnnotation);
					}
				}
			}
		}
		return annotations;
	}

	/**
	 * Finds all vertices and triangles which are on the edge between this annotation and the given
	 * neighbor annotation.
	 * 
	 * @param cas
	 *            MeshCas
	 * @param neighbor
	 *            The neighbor annotation between which the edge should be found
	 * @param edgeVertices
	 *            All vertices found on the edge will be added to this set. It should (but doesn't
	 *            have to) be empty.
	 * @param edgeTriangles
	 *            All triangles found on the edge will be added to this map. It should (but doesn't
	 *            have to) be empty.
	 */
	public <T extends MeshAnnotation> void getNeighborEdge(MeshCas cas, T neighbor,
			Set<Vertex> edgeVertices, HashMultimap<Triangle, Triangle> edgeTriangles) {

		for (Triangle t : mesh.getTriangles()) {
			for (Triangle neigTriangle : t.getNeighbors()) {
				if (neighbor.containsTriangle(neigTriangle)) {
					// arriving here, we found two triangles which represent the edge
					Edge commEdge = t.getCommonEdge(neigTriangle);
					if (commEdge != null) {
						edgeVertices.add(commEdge.getVerticesOfEdge()[0]);
						edgeVertices.add(commEdge.getVerticesOfEdge()[1]);
						edgeTriangles.put(t, neigTriangle);
					}
				}
			}
		}	
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

}
