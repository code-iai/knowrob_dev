/*******************************************************************************
 * Copyright (c) 2014 Andrei Stoica. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Andrei Stoica - initial API and implementation during
 * 								 Google Summer of Code 2014
 ******************************************************************************/

package edu.tum.cs.vis.model.util;

import edu.tum.cs.vis.model.util.Curvature;
import edu.tum.cs.vis.model.util.Vertex;

import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;

/**
 * Class which implements a cluster of vertices with similar curvature properties.
 * The cluster contains an identifying label and the vertices that describe it.
 * 
 * @author Andrei Stoica
 */
public class Cluster {

	/**
	 * Cluster numeric identifier
	 */
	private final int				id;

	/**
	 * Array of floating point values containing the minimum curvature kMin,
	 * maximum curvature kMax and the kMinkMax curvature values of the cluster
	 * of 3-length data points associated with the vertices.
	 */
	private final float[]			centroid		= new float[3];

	/**
	 * List of all vertices that are part of this cluster instance.
	 */
	private final List<Vertex>	clusterVertices	= new ArrayList<Vertex>();

	/**
	 * Constructor for {@code Cluster} class which 
	 * assigns the id and sets the centroid to {0, 0, 0}.
	 * 
	 * @param newId 
	 * 			id assigned to the cluster
	 */
	public Cluster(final int newId) {
		this.id = newId;
		this.centroid[0] = this.centroid[1] = this.centroid[2] = 0.0f;
	}

	/**
	 * Gets the label integer identifier of the cluster instance
	 * 
	 * @return integer id of cluster
	 */
	public int getLabelId() {
		return this.id;
	}

	/**
	 * Gets the centroid of the cluster. 
	 * 
	 * @return array of 3 curvature properties
	 * 			pos[0] = kMin
	 * 			pos[1] = kMax
	 * 			pos[2] = kMinkMax
	 */
	public float[] getCentroid() {
		return this.centroid;
	}

	/**
	 * Gets the list of vertices that are part of the cluster
	 * 
	 * @return list of cluster vertices
	 */
	public List<Vertex> getVertices() {
		return this.clusterVertices;
	}

	/**
	 * Adds vertex to the cluster if this is not already present in it
	 * 
	 * @param v
	 * 			vertex to be added to cluster
	 */
	public void addVertexToCluster(Vertex v) {
		if (v != null && !this.clusterVertices.contains(v)) {
			this.clusterVertices.add(v);
		}
	}

	/**
	 * Removes vertex from the cluster if this is present in it
	 * 
	 * @param v
	 * 			vertex to be removed from cluster
	 */
	public void removeVertexFromCluster(Vertex v) {
		if (v != null && this.clusterVertices.contains(v)) {
			this.clusterVertices.remove(v);
		}
	}

	/**
	 * Recalculates the centroid of the cluster based on the current available vertices
	 * and their corresponding estimated curvatures.
	 * 
	 * @param curvatures
	 * 			curvature map containing all the model vertices as keys 
	 * 			and their estimated curvatures as the corresponding values
	 */
	public void updateCentroid(HashMap<Vertex, Curvature> curvatures) {
		if (curvatures == null || this.clusterVertices.size() == 0) {
			// exit if no vertices or invalid curvatures
			return;
		}
		float numValidVertices = 0f;
		this.centroid[0] = this.centroid[1] = this.centroid[2] = 0.0f;
		for (Vertex v : this.clusterVertices) {
			if (curvatures.get(v) != null) {
				this.centroid[0] += curvatures.get(v).getCurvatureMin();
				this.centroid[1] += curvatures.get(v).getCurvatureMax();
				this.centroid[2] += curvatures.get(v).getCurvatureMinMax();
				numValidVertices += 1f;
			}
		}
		this.centroid[0] /= numValidVertices;
		this.centroid[1] /= numValidVertices;
		this.centroid[2] /= numValidVertices;
	}
	
	@Override
	public String toString() {
		String print = "Cluster ID " + this.id + "\n";
		print = print + "Centroid: (" + this.centroid[0] + ", " + this.centroid[1] + ", " + this.centroid[2] + ")\n";
		print = print + "Vertices: " + this.clusterVertices.size() + "\n";
		return print;
	}
}