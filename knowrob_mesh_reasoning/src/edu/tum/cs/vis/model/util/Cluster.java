/*******************************************************************************
 * Copyright (c) 2014 Andrei Stoica. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Andrei Stoica - initial API and implementation, Year: 2014
 ******************************************************************************/

package edu.tum.cs.vis.model.util;

import edu.tum.cs.vis.model.util.Curvature;
import edu.tum.cs.vis.model.util.Vertex;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Class that implements API of a cluster object to be used in the KMeans classification algorithm
 * 
 * @author Andrei Stoica
 */
public class Cluster {

	/**
	 * Cluster numeric id
	 */
	private final int				id;

	/**
	 * Centroid of the cluster containing Kmin on first position and Kmax on second position
	 */
	private final float[]			centroid		= new float[2];

	/**
	 * Array list of all vertices that are classified to this cluster instance
	 */
	private final ArrayList<Vertex>	clusterVertices	= new ArrayList<Vertex>();

	/**
	 * Constructor for class Cluster
	 */
	public Cluster(final int newId) {
		this.id = newId;
		this.centroid[0] = this.centroid[1] = 0.0f;
	}

	/**
	 * Getter for the label id
	 */
	public int getLabelId() {
		return id;
	}

	/**
	 * Getter for the centroid
	 */
	public float[] getCentroid() {
		return centroid;
	}

	/**
	 * Getter for the ArrayList
	 */
	public ArrayList<Vertex> getVertices() {
		return clusterVertices;
	}

	/**
	 * Adds vertex to the cluster if this is not already present in it
	 */
	public void addVertexToCluster(Vertex v) {
		if (!this.clusterVertices.contains(v)) {
			this.clusterVertices.add(v);
		}
	}

	/**
	 * Removes vertex from the cluster if this is present in it
	 */
	public void removeVertexFromCluster(Vertex v) {
		if (this.clusterVertices.contains(v)) {
			this.clusterVertices.remove(v);
		}
	}

	/**
	 * Recalculates centroid of cluster with the available vertices
	 */
	public void updateCentroid(HashMap<Vertex, Curvature> curvatures) {
		if (clusterVertices.size() == 0) {
			// exit if no vertices
			return;
		}
		centroid[0] = centroid[1] = 0.0f;
		for (Vertex v : clusterVertices) {
			centroid[0] += curvatures.get(v).getCurvatureMin();
			centroid[1] += curvatures.get(v).getCurvatureMax();
		}
		centroid[0] /= clusterVertices.size();
		centroid[1] /= clusterVertices.size();
	}
}