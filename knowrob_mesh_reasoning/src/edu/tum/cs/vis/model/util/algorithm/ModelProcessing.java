/*******************************************************************************
 * Copyright (c) 2014 Andrei Stoica. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Andrei Stoica - initial API and implementation, Year: 2014
 ******************************************************************************/
package edu.tum.cs.vis.model.util.algorithm;

import java.lang.Math;

import java.util.LinkedList;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.concurrent.Callable;

import javax.vecmath.Vector3f;

import edu.tum.cs.vis.model.Model;
import edu.tum.cs.vis.model.util.Vertex;
import edu.tum.cs.vis.model.util.Triangle;
import edu.tum.cs.ias.knowrob.utils.ThreadPool;

public class ModelProcessing{

	/**
	 * Model to be processed
	 */
	protected Model				model;
	
	/**
	 * Flag that shows if sharp edge detection has been done
	 */
	private boolean 			modelSharpEdgeDetectionCheck = false;
	
	/**
	 * Flag that shows if the classification of the curvatures has been done
	 */
	private boolean				modelKMeansCurvatureClassification = false;
	
	/**
	 * Default constructor of the ModelProcessing class
	 */
	public ModelProcessing() {
		this.model = null;
	}
	
	/**
	 * Constructor of the ModelProcessing class
	 * 
	 * @param newModel 
	 * 			model to be analyzed
	 */
	public ModelProcessing(Model newModel) {
		this.model = newModel;
	}
	
	/**
	 * Constructor of the ModelProcessing class
	 * 
	 * @param newModel
	 * 			model to be analyzed
	 * @param isSharpEdgeDetected
	 * 			sets the flag of the sharp edge detection feature
	 * @param isKMeansCLassified
	 * 			sets the flag of the KMeansClassification process
	 */
	public ModelProcessing(Model newModel, boolean isSharpEdgeDetected, boolean isKMeansClassified) {
		this.model = newModel;
		this.modelSharpEdgeDetectionCheck = isSharpEdgeDetected;
		this.modelKMeansCurvatureClassification = isKMeansClassified;
	}
	
	/**
	 * Getter for the model of an instance
	 */
	public Model getModel() {
		return model;
	}
	
	/**
	 * Setter for the model of an instance
	 */
	public void setModel(Model newModel) {
		this.model = newModel;
	}
	
	/**
	 * Getter for the sharp edge detection flag
	 */
	public boolean isSharpEdgeDetectionChecked() {
		return modelSharpEdgeDetectionCheck;
	}
	
	/**
	 * Setter for the sharp edge detection flag
	 */
	public void SharpEdgeDetectionChecked(final boolean value) {
		this.modelSharpEdgeDetectionCheck = value;
	}
	
	/**
	 * Getter for the K Means classification of curvatures
	 */
	public boolean isKMeansCurvatureClassified() {
		return modelKMeansCurvatureClassification;
	}
	
	/**
	 * Setter for the K Means classification flag
	 */
	public void KMeansCurvatureClassified(final boolean value) {
		this.modelKMeansCurvatureClassification = value;
	}
	
	/**
	 * Detects the model sharp edges, marks them and adds additional 
	 * points for the "sharp triangles" in order to correct
	 * for problematic curvature computation
	 */
	public void sharpEdgeDetection() {
		this.removeCollinearTriangles();
		
		List<Callable<Void>> threads = new LinkedList<Callable<Void>>();

		final int interval = 500;

		for (int start = 0; start < model.getTriangles().size(); start += interval) {
			final int st = start;
			threads.add(new Callable<Void>() {

				@Override
				public Void call() throws Exception {
					int end = Math.min(st + interval, model.getTriangles().size());
					for (int i = st; i < end; i++) {
						sharpEdgeDetectionForTriangle(model.getTriangles().get(i));
					}
					return null;
				}

			});
		}
		ThreadPool.executeInPool(threads);
		
		for (int i = 0 ; i < model.getTriangles().size() ; ++i) {
			model.getTriangles().get(i).checkIsSharpTriangle();
			if (model.getTriangles().get(i).isSharpTriangle()) {
				addTrianglesToModel(model.getTriangles().get(i));
				model.getTriangles().remove(i);
			}
		}
		this.removeCollinearTriangles();
		this.modelSharpEdgeDetectionCheck = true;
	}
	
	/**
	 * Performs the sharp detection at the triangle level
	 */
	private void sharpEdgeDetectionForTriangle(Triangle t) {
		Iterator<Triangle> it = t.getNeighbors().iterator();
		while (it.hasNext()) {
			Triangle n = it.next();
			float angleOfNormals = (float)Math.toDegrees(t.getNormalVector().angle(n.getNormalVector()));
			if ((angleOfNormals >= 80.0) && (angleOfNormals <= 110.0)) {
				List<Vertex> vShared = findSharedVertices(t,n);
				Vector3f edge = new Vector3f();
				edge.sub(vShared.get(0),vShared.get(1));
				Vector3f edgeInv = new Vector3f();
				edgeInv.negate(edge);
				vShared.get(0).isSharpVertex(true);
				vShared.get(1).isSharpVertex(true);
				for (Vector3f edge1 : t.getEdges()) {
					if (edge1.equals(edge) || edge1.equals(edgeInv)) {
						t.addSharpEdge(edge1);
						break;
					}
				}
				for (Vector3f edge2 : n.getEdges()) {
					if (edge2.equals(edge) || edge2.equals(edgeInv)) {
						n.addSharpEdge(edge2);
						break;
					}
				}
			}
		}
	}
	
	/**
	 * Finds the two common points of two neighboring trinagles
	 */
	private List<Vertex> findSharedVertices(Triangle t, Triangle n) {
		List<Vertex> v = new ArrayList<Vertex>(2);
		for (int i = 0 ; i < t.getPosition().length ; ++i) {
			for (int j = 0 ; j < n.getPosition().length ; ++j) {
				if (t.getPosition()[i].sameCoordinates(n.getPosition()[j])) {
					v.add(t.getPosition()[i]);
					break;
				}
			}
		}
		return v;
	}
	
	/**
	 * Adds triangles to model using the centroid of the triangle
	 * 
	 * @param t 
	 * 			triangle decomposed in 3 smaller triangles
	 */
	private void addTrianglesToModel(Triangle t) {
		Vertex newVertex = new Vertex(t.getCentroid().x, t.getCentroid().y, t.getCentroid().z);
		model.getVertices().add(newVertex);
		Triangle[] newTriangle = new Triangle[3];
		for (int i = 0 ; i < 3 ; ++i) {
			newTriangle[i] = new Triangle(t.getPosition()[i],newVertex,t.getPosition()[(i+1)%3]);
			newTriangle[i].setAppearance(t.getAppearance());
			newTriangle[i].setNormalVector(t.getNormalVector());
		}
		
		// add neighbors inside the big triangle
		newTriangle[0].addNeighbor(newTriangle[1]);
		newTriangle[0].addNeighbor(newTriangle[2]);
		newTriangle[1].addNeighbor(newTriangle[2]);
		
		// add triangle neighbors outside the original triangle 
		List<Triangle> neighbors = new ArrayList<Triangle>();
		neighbors.addAll(t.getNeighbors());
		for (int i = 0 ; i < 3 ; ++i) {
			for (int j = 0 ; j < neighbors.size() ; ++j) {
				Triangle n = neighbors.get(j);
				int cont = 0;
				for (Vertex v : n.getPosition()) {
					if ((v.sameCoordinates(t.getPosition()[i])) || (v.sameCoordinates(t.getPosition()[(i+1)%3]))) {
						cont++;
					}
				}
				// if the neighboring triangle (exact 2 common vertices)
				if (cont == 2) {
					t.removeNeighbor(n);
					neighbors.remove(n);
					newTriangle[i].addNeighbor(n);
				}
			}
		}
		
		// add sharp edges if any to the new 3 created triangles
		for (Vector3f sharpEdge : t.getSharpEdges()) {
			for (int i = 0 ; i < t.getEdges().length ; ++i) {
				Vector3f sharpEdgeInv = new Vector3f();
				sharpEdgeInv.negate(sharpEdge);
				if (sharpEdge.equals(t.getEdges()[i]) || sharpEdgeInv.equals(t.getEdges()[i])) {
					newTriangle[(i+2) % t.getEdges().length].addSharpEdge(sharpEdge);
				}
			}
		}
		
		// add new vertex as direct neighbor and old vertices as neighbors for new one
		// compute centroids of new triangles and add them to the model
		for (int i = 0 ; i < 3 ; ++i) {
			t.getPosition()[i].addNeighbor(newVertex);
			newVertex.addNeighbor(t.getPosition()[i]);
			newTriangle[i].updateCentroid();
			model.getTriangles().add(newTriangle[i]);
		}
	}
	
	/**
	 * Removes "colinear triangles", i.e. "triangles" which have 3 "colinear" vertices
	 */
	private void removeCollinearTriangles() {
		List<Triangle> allTriangles = new ArrayList<Triangle>();
		allTriangles.addAll(model.getTriangles());
		for (int i = 0 ; i < allTriangles.size() ; ++i) {
			Vector3f crossProd = new Vector3f();
			crossProd.cross(allTriangles.get(i).getEdges()[0], allTriangles.get(i).getEdges()[1]);
			if (crossProd.length() == 0.0) {
				System.out.println("removing" + allTriangles.get(i));
				List<Triangle> tn = new ArrayList<Triangle>();
				tn.addAll(allTriangles.get(i).getNeighbors());
				for (int j = 0 ; j < tn.size() ; ++j) {
					tn.get(j).removeNeighbor(allTriangles.get(i));
				}
				model.getTriangles().remove(allTriangles.get(i));
				allTriangles.remove(allTriangles.get(i));
			}
		}
	}
}

