/*******************************************************************************
 * Copyright (c) 2014 Andrei Stoica. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Andrei Stoica - initial API and implementation, Year: 2014
 ******************************************************************************/

package edu.tum.cs.vis.model.util.algorithm;

import java.lang.Math;

import org.apache.log4j.Logger;

import java.awt.Color;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.concurrent.Callable;

import javax.vecmath.Vector3f;

import edu.tum.cs.util.PrintUtil;
import edu.tum.cs.vis.model.Model;
import edu.tum.cs.vis.model.util.Appearance;
import edu.tum.cs.vis.model.util.Cluster;
import edu.tum.cs.vis.model.util.Curvature;
import edu.tum.cs.vis.model.util.Edge;
import edu.tum.cs.vis.model.util.Region;
import edu.tum.cs.vis.model.util.Triangle;
import edu.tum.cs.vis.model.util.Vertex;
import edu.tum.cs.ias.knowrob.utils.ThreadPool;


/**
 * Class that implements the processing blocks used to
 * render better segmentation results of the models 
 * analyzed. Blocks can be either pre- or post-processing
 * units in the computation flow applied to the model
 * 
 * @author Andrei Stoica
 */
public class ModelProcessing{

	private static Logger		logger = Logger.getLogger(ModelProcessing.class);
	/**
	 * Model to be processed
	 */
	protected Model						model;
	
//	/**
//	 * Flag that shows if sharp edge detection has been done
//	 */
//	private boolean 					modelSharpEdgeDetectionCheck = false;
	
	/**
	 * Flag that shows if the classification of the curvatures has been done
	 */
	private boolean						modelKMeansCurvatureClassification = false;
	
	/**
	 * Defines the number of clusters used for classifying the vertices
	 */	
	private static int					NUM_OF_CLUSTERS = 30;
	
	/**
	 * Defines upper iteration limit for the KMeans algorithm
	 */
	private final static int			UPPER_ITERATION_LIMIT = 100;
	
	/**
	 * Defines upper iteration limit for the region merging
	 */
	private static int			UPPER_ITERATION_LIMIT_GROWING = 180;
	
	/**
	 * Defines the minimal area for small regions merging step
	 */
	private final static float			AREA_MIN_LIMIT = 5e-2f;
	
	/**
	 * Define minimum distance threshold for the region merging which stops the merging
	 */
	private final static float			MIN_DISTANCE_THRESHOLD = 5e-12f;
	
	/**
	 * Defines the value of the area weighting on the distance calculation
	 * of the merging process
	 */
	private final static float			EPSILON = 1e-5f;
	
	/**
	 * Defines the maximum number for the 32-bit floating point precision (+Inf)
	 */
	private final static float			PINF = Float.MAX_VALUE;
	
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
	 * Setter for the model of an instance
	 */
	public void setModel(Model newModel) {
		this.model = newModel;
	}
	
	/**
	 * Getter for the model of an instance
	 */
	public Model getModel() {
		return model;
	}
	
//	/**
//	 * Getter for the number of the added triangles
//	 */
//	public int getNumAddedTriangles() {
//		return numAddedTriangles;
//	}
	
	/**
	 * Gets number of clusters defined
	 */
	public int getNumOfClusters() {
		return NUM_OF_CLUSTERS;
	}
	
//	/**
//	 * Getter for the sharp edge detection flag
//	 */
//	public boolean isSharpEdgeDetectionChecked() {
//		return modelSharpEdgeDetectionCheck;
//	}
	
	/**
	 * Getter for the K Means classification of curvatures
	 */
	public boolean isKMeansCurvatureClassified() {
		return modelKMeansCurvatureClassification;
	}
	
//	/**
//	 * Detects the model sharp edges, marks them and adds additional 
//	 * points for the "sharp triangles" in order to correct
//	 * for problematic curvature computation
//	 */
//	public void sharpEdgeDetection() {
//		// remove any wrong tesselated triangles in the model
//		this.removeCollinearTriangles();
//		
//		this.numAddedTriangles = this.model.getTriangles().size();
//		
//		// perform the sharp edge detection for individual triangles (multi-threaded)
//		List<Callable<Void>> threads = new LinkedList<Callable<Void>>();
//
//		final int interval = 500;
//
//		for (int start = 0; start < model.getTriangles().size(); start += interval) {
//			final int st = start;
//			threads.add(new Callable<Void>() {
//
//				@Override
//				public Void call() throws Exception {
//					int end = Math.min(st + interval, model.getTriangles().size());
//					for (int i = st; i < end; i++) {
//						sharpEdgeDetectionForTriangle(model.getTriangles().get(i));
//					}
//					return null;
//				}
//
//			});
//		}
//		ThreadPool.executeInPool(threads);
//
////		for (int i = 0 ; i < model.getTriangles().size() ; ++i) {
////			Triangle t = model.getTriangles().get(i);
////			System.out.println(i + " " + t.getSharpEdges());
////		}
//		
//		// for sharp triangles add new point and introduce 3 new triangles
//		List<Triangle> toRemove = new ArrayList<Triangle>();
//		for (int i = 0 ; i < model.getTriangles().size() ; ++i) {
//			Triangle t = model.getTriangles().get(i);
//			if (t.checkIsSharpTriangle()) {
//				addTrianglesToModel(t);
//				toRemove.add(t);
//			}
//		}
//		model.getTriangles().removeAll(toRemove);
//		
////		this.removeCollinearTriangles();
//	
//		this.numAddedTriangles = this.model.getTriangles().size() - this.numAddedTriangles;
//		if (this.numAddedTriangles != 0) {
//			logger.debug("Added " + this.numAddedTriangles + " triangles to the model");
//		}
//		else {
//			logger.debug("No triangles added to the model");
//		}
//		//
//		this.modelSharpEdgeDetectionCheck = true;
//		
//		for (int i = 0 ; i < model.getTriangles().size() ; ++i) {
//			Triangle t = model.getTriangles().get(i);
//			System.out.println(i + " " + model.getTriangles().get(i) + "sharp edges:\n" + t.getSharpEdges() + "\n" + t.getNeighbors() + "\n");
//		}
//	}
//	
//	/**
//	 * Performs the sharp detection at the triangle level
//	 */
//	private void sharpEdgeDetectionForTriangle(Triangle t) {
//		synchronized (t) {
//		Iterator<Triangle> it = t.getNeighbors().iterator();
//		while (it.hasNext()) {
//			Triangle n = it.next();
//			synchronized (n) {
//			float angleOfNormals = (float)Math.toDegrees(t.getNormalVector().angle(n.getNormalVector()));
//			if ((angleOfNormals >= 80.0) && (angleOfNormals <= 110.0)) {
//				List<Vertex> vShared = findSharedVertices(t,n);
//				Edge edge = new Edge(vShared.get(0), vShared.get(2));
//				vShared.get(0).isSharpVertex(true);
//				vShared.get(1).isSharpVertex(true);
//				vShared.get(2).isSharpVertex(true);
//				vShared.get(3).isSharpVertex(true);
//				t.addSharpEdge(edge);
//				n.addSharpEdge(edge);
//			}
//			}
//		}
//		if (t.getNeighbors().size() < 3) {
//			Edge[] edges = t.getEdges();
//			for (int i = 0 ; i < edges.length ; ++i) {
//				Triangle n = t.getNeighborOfEdge(edges[i]);
//				if (n == null) {
//					// mark vertices that define the edge as being sharp
//					t.getPosition()[(i+1) % edges.length].isSharpVertex(true);
//					t.getPosition()[(i+2) % edges.length].isSharpVertex(true);
//					t.addSharpEdge(edges[i]);
//				}
//				
//			}
//		}
//		}
//	}
//	
//	/**
//	 * Finds the two common points of two neighboring trinagles
//	 */
//	private List<Vertex> findSharedVertices(Triangle t, Triangle n) {
//		List<Vertex> v = new ArrayList<Vertex>(4);
//		for (int i = 0 ; i < t.getPosition().length ; ++i) {
//			for (int j = 0 ; j < n.getPosition().length ; ++j) {
//				if (t.getPosition()[i].sameCoordinates(n.getPosition()[j])) {
//					v.add(t.getPosition()[i]);
//					v.add(n.getPosition()[j]);
//					break;
//				}
//			}
//		}
//		return v;
//	}
//	
//	/**
//	 * Adds triangles to model using the centroid of the triangle
//	 * 
//	 * @param t 
//	 * 			triangle decomposed in 3 smaller triangles
//	 */
//	private void addTrianglesToModel(Triangle t) {
//		Vertex newVertex = new Vertex(t.getCentroid().x, t.getCentroid().y, t.getCentroid().z);
//		model.getVertices().add(newVertex);
//		Triangle[] newTriangle = new Triangle[3];
//		for (int i = 0 ; i < 3 ; ++i) {
//			newTriangle[i] = new Triangle(t.getPosition()[i],t.getPosition()[(i+1)%3],newVertex);
//			newTriangle[i].setAppearance(t.getAppearance());
//			newTriangle[i].setNormalVector(t.getNormalVector());
//		}
//
//		// add neighbors inside the big triangle
//		newTriangle[0].addNeighbor(newTriangle[1]);
//		newTriangle[0].addNeighbor(newTriangle[2]);
//		newTriangle[1].addNeighbor(newTriangle[2]);
//		
//		// add triangle neighbors outside the original triangle 
//		Edge[] edges = t.getEdges();
//		for (int i = 0 ; i < edges.length ; ++i) {
//			Triangle n = t.getNeighborOfEdge(edges[i]);
//			if (n == null) {
//				continue;
//			}
//			for (int j = 0 ; j < 3 ; ++j) {
//				if (newTriangle[j].containsEdge(edges[i])) {
//					n.removeNeighbor(t);
//					n.addNeighbor(newTriangle[j]);
//					newTriangle[j].addNeighbor(n);
//					break;
//				}
//			}
//		}
//		
////		
////		List<Triangle> neighbors = new ArrayList<Triangle>();
////		neighbors.addAll(t.getNeighbors());
////		for (int i = 0 ; i < 3 ; ++i) {
////			for (int j = 0 ; j < neighbors.size() ; ++j) {
////				logger.debug("j = " + j);
////				Triangle n = neighbors.get(j);
////				int cont = 0;
////				for (Vertex v : n.getPosition()) {
////					if ((v.sameCoordinates(newTriangle[i].getPosition()[0])) || (v.sameCoordinates(newTriangle[i].getPosition()[1]))) {
////						cont++;
////					}
////				}
////				// if the neighboring triangle (exact 2 common vertices)
////				if (cont == 2) {
////					t.removeNeighbor(n);
////					neighbors.remove(n);
////					newTriangle[i].addNeighbor(n);
////					n.addNeighbor(newTriangle[i]);
////					break;
////				}
////			}
////		}
//		
//		// add sharp edges if any to the new 3 created triangles
//		for (Edge sharpEdge : t.getSharpEdges()) {
//			for (int i = 0 ; i < 3 ; ++i) {
//				newTriangle[i].addSharpEdge(sharpEdge);
//			}
//		}
//		
//		// add new vertex as direct neighbor and old vertices as neighbors for new one
//		// compute centroids of new triangles and add them to the model
//		for (int i = 0 ; i < 3 ; ++i) {
//			t.getPosition()[i].addNeighbor(newVertex);
//			newVertex.addNeighbor(t.getPosition()[i]);
//			newTriangle[i].updateCentroid();
//			model.getTriangles().add(newTriangle[i]);
//		}
//	}
//	
//	/**
//	 * Removes "colinear triangles", i.e. "triangles" which have 3 "colinear" vertices
//	 */
//	private void removeCollinearTriangles() {
//		List<Triangle> allTriangles = new ArrayList<Triangle>();
//		allTriangles.addAll(model.getTriangles());
//		int rmTriangles = 0;
//		for (int i = 0 ; i < allTriangles.size() ; ++i) {
//			Edge[] edges = allTriangles.get(i).getEdges();
//			Vector3f crossProd = new Vector3f();
//			crossProd.cross(edges[0].getEdgeValue(), edges[1].getEdgeValue());
//			if (crossProd.length() == 0.0 || edges[0].getEdgeValue().length() == 0.0 || 
//					edges[1].getEdgeValue().length() == 0.0 || edges[2].getEdgeValue().length() == 0.0) {
//				logger.debug("removing" + allTriangles.get(i));
//				List<Triangle> tn = new ArrayList<Triangle>();
//				tn.addAll(allTriangles.get(i).getNeighbors());
//				for (int j = 0 ; j < tn.size() ; ++j) {
//					tn.get(j).removeNeighbor(allTriangles.get(i));
//				}
//				model.getTriangles().remove(allTriangles.get(i));
//				allTriangles.remove(allTriangles.get(i));
//				rmTriangles++;
//			}
//		}
//		if (rmTriangles > 0) {
//			logger.debug("Removed " + rmTriangles + " triangles");
//		}
//	}
	
	/**
	 * KMeans algorithm implementation for vertex
	 * curvature classification
	 */
	public void KMeansVCClassification(HashMap<Vertex,Curvature> curvatures) {
		if (NUM_OF_CLUSTERS >= model.getVertices().size() / 10) {
			logger.debug("Number of vertices in the model smaller than number of clusters chosen");
			int temp = NUM_OF_CLUSTERS;
			NUM_OF_CLUSTERS = temp / 5;
			logger.debug("Number of clusters has been reduced from " +
			temp + " to " + NUM_OF_CLUSTERS);
		}
		Cluster[] clusters = new Cluster[NUM_OF_CLUSTERS];
		
		// randomly initialize clusters with one element each
		List<Integer> dummyRnd = new ArrayList<Integer>();
		final List<Integer> pickedData = new ArrayList<Integer>();
		List<Vertex> verticesData = new ArrayList<Vertex>(model.getVertices());
		for (int i = 0 ; i < verticesData.size() ; ++i) {
			dummyRnd.add(i);
		}
		Collections.shuffle(dummyRnd);
		pickedData.addAll(dummyRnd.subList(0, NUM_OF_CLUSTERS));
		for (int i = 0 ; i < NUM_OF_CLUSTERS ; ++i) {
			clusters[i] = new Cluster(i);
			// System.out.println(clusters[i].getLabelId());
			Collections.shuffle(dummyRnd);
			// store random picked element to cluster
			Vertex v = verticesData.get(pickedData.get(i));
			clusters[i].addVertexToCluster(v);
			clusters[i].updateCentroid(curvatures);
			v.setClusterLabel(clusters[i].getLabelId());
			v.setClusterCurvatureVal(clusters[i].getCentroid()[0],clusters[i].getCentroid()[1]);
		}
		
		// remove randomly picked elements
		for (int i = 0 ; i < NUM_OF_CLUSTERS ; ++i) {
			verticesData.remove(clusters[i].getVertices().get(0));
		}
		
		// add all remaining vertices to closest cluster according to Euclidean
		// distance applied on the curvature space
		while (!verticesData.isEmpty()) {
			Vertex v = verticesData.remove(0);
			float distMin = this.distEuclid(curvatures.get(v), clusters[0]);
			int clusterIndex = 0;
			for (int i = 1 ; i < NUM_OF_CLUSTERS ; ++i) {
				float distTmp = this.distEuclid(curvatures.get(v), clusters[i]);
				if (distTmp < distMin) {
					distMin = distTmp;
					clusterIndex = i;
				}
			}
			clusters[clusterIndex].addVertexToCluster(v);
			clusters[clusterIndex].updateCentroid(curvatures);
			v.setClusterLabel(clusters[clusterIndex].getLabelId());
			v.setClusterCurvatureVal(clusters[clusterIndex].getCentroid()[0],clusters[clusterIndex].getCentroid()[1]);
		}
		
//		for (int i = 0 ; i < NUM_OF_CLUSTERS ; ++i) {
//			System.out.println("Cluster Id " + clusters[i].getLabelId());
//			System.out.println("Vertices: " + clusters[i].getVertices());
//		}
		
		// iterate until no significant changes are visible or iteration
		// limit is exceeded
		boolean isRunning = true;
		int iteration = 0;
		while (isRunning && iteration < UPPER_ITERATION_LIMIT) {
			isRunning = false;
			for (int i = 0 ; i < model.getVertices().size() ; ++i) {
				float distMin = this.distEuclid(curvatures.get(model.getVertices().get(i)), clusters[0]);
				int clusterIndex = 0;
				for (int j = 1 ; j < NUM_OF_CLUSTERS ; ++j) {
					float distTmp = this.distEuclid(curvatures.get(model.getVertices().get(i)), clusters[j]);
					if (distTmp < distMin) {
						distMin = distTmp;
						clusterIndex = j;
					}
				}
				//System.out.println(model.getVertices().get(i));
				if (model.getVertices().get(i).getClusterLabel() != clusterIndex) {
					//System.out.println(model.getVertices().get(i).getClusterLabel());
					clusters[model.getVertices().get(i).getClusterLabel()].removeVertexFromCluster(model.getVertices().get(i));
					clusters[model.getVertices().get(i).getClusterLabel()].updateCentroid(curvatures);
					clusters[clusterIndex].addVertexToCluster(model.getVertices().get(i));
					clusters[clusterIndex].updateCentroid(curvatures);
					model.getVertices().get(i).setClusterLabel(clusterIndex);
					model.getVertices().get(i).setClusterCurvatureVal(clusters[clusterIndex].getCentroid()[0], clusters[clusterIndex].getCentroid()[1]);
					isRunning = true;
				}
			}
			iteration++;
		}
		this.modelKMeansCurvatureClassification = true;
		
		int classifiedVertices = 0;
		for (int i = 0 ; i < NUM_OF_CLUSTERS ; ++i) {
			classifiedVertices += clusters[i].getVertices().size();
			for (int j = 0 ; j < clusters[i].getVertices().size() ; ++j) {
				Curvature c = curvatures.get(clusters[i].getVertices().get(j));
				clusters[i].getVertices().get(j).setClusterCurvatureVal(clusters[i].getCentroid()[0], clusters[i].getCentroid()[1]);
				c.setCurvatureMin(clusters[i].getCentroid()[0]);
				c.setCurvatureMax(clusters[i].getCentroid()[1]);
			}
		}
		
		// compute new coloring used for primitive fitting
		CurvatureCalculation.setCurvatureHueSaturation(curvatures, model, 1f);
		
//		for (int i = 0 ; i < NUM_OF_CLUSTERS ; ++i) {
//			System.out.println(clusters[i]);
//		}
		
		logger.debug("Classified " + classifiedVertices + " vertices out of " + model.getVertices().size() + " into " + NUM_OF_CLUSTERS + " clusters");
	}
	
	/**
	 * Computes the Euclidean distance between the centroid of a cluster
	 * and the curvature of a vertex
	 *  
	 * @param c
	 * 			curvature associated with selected vertex
	 * @param cluster
	 * 			cluster considered
	 * @return
	 * 			Euclidean distance between c and cluster.centroid
	 */
	private float distEuclid(Curvature c, Cluster cluster) {
		return (float)Math.sqrt(Math.pow(c.getCurvatureMin() - cluster.getCentroid()[0] , 2) 
				+ Math.pow(c.getCurvatureMax() - cluster.getCentroid()[1] , 2));
	}
	
	/**
	 * Method to colorize the traingles of the mesh randomly
	 * according to their regions
	 * 
	 * @param regions 
	 * 			to be colorized
	 */
	private void colorizeRegions(List<Region> regions) {
		int min = 0;
		int max = 255;
		for (int i = 0 ; i < regions.size() ; ++i) {
			Region r = regions.get(i);
			Appearance a = new Appearance();
			int R,G,B;
			R = min + (int)(Math.random() * ((max - min) + 1));
			G = min + (int)(Math.random() * ((max - min) + 1));
			B = min + (int)(Math.random() * ((max - min) + 1));
			Color c = new Color(R,G,B);
			for (Triangle t : r.getTriangles()) {
				a.setColorFill(c);
//				a.setColorLine(c);
				t.setAppearance(a);
			}
		}
	}
	
	private void computeAdjacencyDistances(final float[][] adjacencyMatrix, final Region r) {
		List<Region> regions = model.getRegions();
		for (int j = 0 ; j < regions.size() ; ++j) {
			Region n = regions.get(j);
			if (!r.isNeighbor(n)) {
				continue;
			}
			float[] curvatureDistanceBoundary = {0.0f, 0.0f};
			float curvatureDistance = 0.0f;
			float commonBorder = 0.0f;
			float neighborhoodDistance = 0.0f;
			float sizeDistance = 0.0f;
			float distance = 0.0f;
			List<Edge> commonEdges = r.getCommonEdges(n);
			if (commonEdges.size() == 1) {
				// if common edge is sharp cannot merge the two regions
				if (commonEdges.get(0).getIsSharpEdge()) {
					adjacencyMatrix[r.getRegionId()][n.getRegionId()] = PINF;
					adjacencyMatrix[n.getRegionId()][r.getRegionId()] = PINF;
					continue;
				}
				curvatureDistanceBoundary[0] = (commonEdges.get(0).getVerticesOfEdge()[0].getClusterCurvatureVal()[0] 
						+ commonEdges.get(0).getVerticesOfEdge()[1].getClusterCurvatureVal()[0]) / 2.0f;
				curvatureDistanceBoundary[1] = (commonEdges.get(0).getVerticesOfEdge()[0].getClusterCurvatureVal()[1] 
						+ commonEdges.get(0).getVerticesOfEdge()[1].getClusterCurvatureVal()[1]) / 2.0f;
				commonBorder += commonEdges.get(0).getEdgeValue().length();
			}
			else {
				HashMap<Vertex,Integer> edgesVerticesAdjacency = new HashMap<Vertex,Integer>();
				int cont = 0;
				for (int k = 0 ; k < commonEdges.size() && distance != PINF ; ++k) {
					Edge e = commonEdges.get(k);
					if (e.getIsSharpEdge()) {
						distance = PINF;
						adjacencyMatrix[r.getRegionId()][n.getRegionId()] = distance;
						adjacencyMatrix[n.getRegionId()][r.getRegionId()] = distance;
						break;
					}
					if (!edgesVerticesAdjacency.containsKey(e.getVerticesOfEdge()[0])) {
						edgesVerticesAdjacency.put(e.getVerticesOfEdge()[0], 1);
					}
					else {
						edgesVerticesAdjacency.put(e.getVerticesOfEdge()[0], edgesVerticesAdjacency.get(e.getVerticesOfEdge()[0]) + 1);
					}
					if (!edgesVerticesAdjacency.containsKey(e.getVerticesOfEdge()[1])) {
						edgesVerticesAdjacency.put(e.getVerticesOfEdge()[1], 1);
					}
					else {
						edgesVerticesAdjacency.put(e.getVerticesOfEdge()[1], edgesVerticesAdjacency.get(e.getVerticesOfEdge()[1]) + 1);
					}
					commonBorder += e.getEdgeValue().length();
				}
				if (distance == PINF) {
					continue;
				}
				for (Vertex v : edgesVerticesAdjacency.keySet()) {
					if (edgesVerticesAdjacency.get(v) == 2) {
						curvatureDistanceBoundary[0] += v.getClusterCurvatureVal()[0];
						curvatureDistanceBoundary[1] += v.getClusterCurvatureVal()[1];
						cont++;
					}
				}
				curvatureDistanceBoundary[0] /= (float)cont;
				curvatureDistanceBoundary[1] /= (float)cont;
			}
			// compute curvature distance: cD = ||c_reg_r - cDB|| + ||c_reg_n - cDB||
			curvatureDistance = (float)((Math.sqrt(Math.pow(r.getCurvatureMinMaxOfRegion()[0] - curvatureDistanceBoundary[0], 2) 
					+ Math.pow(r.getCurvatureMinMaxOfRegion()[1] - curvatureDistanceBoundary[1], 2))) + 
					(Math.sqrt(Math.pow(n.getCurvatureMinMaxOfRegion()[0] - curvatureDistanceBoundary[0], 2) 
					+ Math.pow(n.getCurvatureMinMaxOfRegion()[1] - curvatureDistanceBoundary[1], 2))) + 
					(Math.sqrt(Math.pow(n.getCurvatureMinMaxOfRegion()[0] - r.getCurvatureMinMaxOfRegion()[0], 2) 
					+ Math.pow(n.getCurvatureMinMaxOfRegion()[1] - r.getCurvatureMinMaxOfRegion()[1], 2))));
			if (curvatureDistance == 0.0f) {
				curvatureDistance = EPSILON;
			}
			// compute neighborhood distance: nD = min(per_r,per_n) / per_common;
			neighborhoodDistance = Math.min(r.getPerimeterOfRegion(), n.getPerimeterOfRegion()) / commonBorder;
			// size distance: sD = 
			if (r.getAreaOfRegion() < AREA_MIN_LIMIT || n.getAreaOfRegion() < AREA_MIN_LIMIT) {
				sizeDistance = EPSILON;
			}
			else {
				sizeDistance = 1.0f;
			}
			distance = curvatureDistance * neighborhoodDistance * sizeDistance;
			// add distance to region adjacency matrix (which is symmetric)
			adjacencyMatrix[r.getRegionId()][n.getRegionId()] = distance;
			adjacencyMatrix[n.getRegionId()][r.getRegionId()] = distance;
		}
	}
	
	private void computeAdjacencyDistances(final float[][] adjacencyMatrix) {
		List<Region> regions = model.getRegions();
		for (int i = 0 ; i < regions.size() ; ++i) {
			Region r = regions.get(i);
			for (int j = 0 ; j <= i ; ++j) {
				Region n = regions.get(j);
				if (!r.isNeighbor(n)) {
					continue;
				}
				float[] curvatureDistanceBoundary = {0.0f, 0.0f};
				float curvatureDistance = 0.0f;
				float commonBorder = 0.0f;
				float neighborhoodDistance = 0.0f;
				float sizeDistance = 0.0f;
				float distance = 0.0f;
				List<Edge> commonEdges = r.getCommonEdges(n);
				if (commonEdges.size() == 1) {
					// if common edge is sharp cannot merge the two regions
					if (commonEdges.get(0).getIsSharpEdge()) {
						adjacencyMatrix[r.getRegionId()][n.getRegionId()] = PINF;
						adjacencyMatrix[n.getRegionId()][r.getRegionId()] = PINF;
						continue;
					}
					curvatureDistanceBoundary[0] = (commonEdges.get(0).getVerticesOfEdge()[0].getClusterCurvatureVal()[0] 
							+ commonEdges.get(0).getVerticesOfEdge()[1].getClusterCurvatureVal()[0]) / 2.0f;
					curvatureDistanceBoundary[1] = (commonEdges.get(0).getVerticesOfEdge()[0].getClusterCurvatureVal()[1] 
							+ commonEdges.get(0).getVerticesOfEdge()[1].getClusterCurvatureVal()[1]) / 2.0f;
					commonBorder += commonEdges.get(0).getEdgeValue().length();
				}
				else {
					HashMap<Vertex,Integer> edgesVerticesAdjacency = new HashMap<Vertex,Integer>();
					int cont = 0;
					for (int k = 0 ; k < commonEdges.size() && distance != PINF ; ++k) {
						Edge e = commonEdges.get(k);
						if (e.getIsSharpEdge()) {
							distance = PINF;
							adjacencyMatrix[r.getRegionId()][n.getRegionId()] = distance;
							adjacencyMatrix[n.getRegionId()][r.getRegionId()] = distance;
							break;
						}
						if (!edgesVerticesAdjacency.containsKey(e.getVerticesOfEdge()[0])) {
							edgesVerticesAdjacency.put(e.getVerticesOfEdge()[0], 1);
						}
						else {
							edgesVerticesAdjacency.put(e.getVerticesOfEdge()[0], edgesVerticesAdjacency.get(e.getVerticesOfEdge()[0]) + 1);
						}
						if (!edgesVerticesAdjacency.containsKey(e.getVerticesOfEdge()[1])) {
							edgesVerticesAdjacency.put(e.getVerticesOfEdge()[1], 1);
						}
						else {
							edgesVerticesAdjacency.put(e.getVerticesOfEdge()[1], edgesVerticesAdjacency.get(e.getVerticesOfEdge()[1]) + 1);
						}
						commonBorder += e.getEdgeValue().length();
					}
					if (distance == PINF) {
						continue;
					}
					for (Vertex v : edgesVerticesAdjacency.keySet()) {
						if (edgesVerticesAdjacency.get(v) == 2) {
							curvatureDistanceBoundary[0] += v.getClusterCurvatureVal()[0];
							curvatureDistanceBoundary[1] += v.getClusterCurvatureVal()[1];
							cont++;
						}
					}
					curvatureDistanceBoundary[0] /= (float)cont;
					curvatureDistanceBoundary[1] /= (float)cont;
				}
				// compute curvature distance: cD = ||c_reg_r - cDB|| + ||c_reg_n - cDB||
				curvatureDistance = (float)((Math.sqrt(Math.pow(r.getCurvatureMinMaxOfRegion()[0] - curvatureDistanceBoundary[0], 2) 
						+ Math.pow(r.getCurvatureMinMaxOfRegion()[1] - curvatureDistanceBoundary[1], 2))) + 
						(Math.sqrt(Math.pow(n.getCurvatureMinMaxOfRegion()[0] - curvatureDistanceBoundary[0], 2) 
						+ Math.pow(n.getCurvatureMinMaxOfRegion()[1] - curvatureDistanceBoundary[1], 2))) + 
						(Math.sqrt(Math.pow(n.getCurvatureMinMaxOfRegion()[0] - r.getCurvatureMinMaxOfRegion()[0], 2) 
						+ Math.pow(n.getCurvatureMinMaxOfRegion()[1] - r.getCurvatureMinMaxOfRegion()[1], 2))));
				if (curvatureDistance == 0.0f) {
					curvatureDistance = EPSILON;
				}
				// compute neighborhood distance: nD = min(per_r,per_n) / per_common;
				neighborhoodDistance = Math.min(r.getPerimeterOfRegion(), n.getPerimeterOfRegion()) / commonBorder;
				// size distance: sD = 
				if (r.getAreaOfRegion() < AREA_MIN_LIMIT || n.getAreaOfRegion() < AREA_MIN_LIMIT) {
					sizeDistance = EPSILON;
				}
				else {
					sizeDistance = 1.0f;
				}
				distance = curvatureDistance * neighborhoodDistance * sizeDistance;
				// add distance to region adjacency matrix (which is symmetric)
				adjacencyMatrix[r.getRegionId()][n.getRegionId()] = distance;
				adjacencyMatrix[n.getRegionId()][r.getRegionId()] = distance;
			}
		}
	}
	
	/**
	 * Merges regions r1 and r2 by integrating r2 into r1 and
	 * recalculating all properties of the r1 instance, while
	 * deleting the r2 instance from the list of regions existent 
	 * in the model
	 * 
	 * @param r1
	 * @param r2
	 */
	private void mergeRegions(Region r1, Region r2) {
		for (Triangle t : r2.getTriangles()) {
			// unset region label and add triangle to r1
			t.setRegionLabel(-1);
			r1.addTriangleToRegion(t);
		}
		
		// update boundary
		r1.updateRegionBoundary();
		
		// update curvature values weighted by area
		float newKMin = 0.0f, newKMax = 0.0f;
		newKMin = (r1.getAreaOfRegion() * r1.getCurvatureMinMaxOfRegion()[0] + r2.getAreaOfRegion() * r2.getCurvatureMinMaxOfRegion()[0]) 
				/ (r1.getAreaOfRegion() + r2.getAreaOfRegion());
		newKMax = (r1.getAreaOfRegion() * r1.getCurvatureMinMaxOfRegion()[1] + r2.getAreaOfRegion() * r2.getCurvatureMinMaxOfRegion()[1]) 
				/ (r1.getAreaOfRegion() + r2.getAreaOfRegion());
		r1.setCurvatureMin(newKMin);
		r1.setCurvatureMax(newKMax);
		
		// remove r2 from the regions list
		model.getRegions().remove(r2);
		
		// re-determine the neighboring relations for updated r1
		r1.updateRegionNeighbors(model.getRegions());
	}
	
	/**
	 * Reverts the coloring of the triangles to default
	 * 
	 * @param regions
	 */
	public void resetRegionsColor(List<Region> regions) {
		for (int i = 0 ; i < model.getTriangles().size() ; ++i) {
			model.getTriangles().get(i).setAppearance(null);
		}
	}
	
	/**
	 * Implements the region growing block for the triangles
	 * inside the model. Should be called after the classification
	 * process has ended.
	 */
	public void processRegionGrowing() {
		logger.debug("Building up regions ...");
		long duration = System.currentTimeMillis();
		int trianglesToClassifyToRegions = model.getTriangles().size();
		int unclassifiedNum = 0;
		List<Region> regions = new ArrayList<Region>();
		int contId = 0;
		
		logger.debug("Triangles to classify: " + trianglesToClassifyToRegions);
		logger.debug("Starting with seed triangles ...");
		for (int i = 0 ; i < model.getTriangles().size() ; ++i) {
			Triangle tr = model.getTriangles().get(i);
			if (tr.updateIsSeedTriangle() && tr.getRegionLabel() == -1) {
				Region newRegion = new Region(contId,tr);
//				if (newRegion.getCurvatureMinMaxOfRegion()[0] == 0.0 && newRegion.getCurvatureMinMaxOfRegion()[1] == 0.0) {
//					System.out.println("0.0 curv @" + tr);
//				}
				newRegion.buildUpRegion();
				regions.add(newRegion);
				contId++;
			}
		}
		
		for (int i = 0 ; i < model.getTriangles().size() ; ++i) {
			if (model.getTriangles().get(i).getRegionLabel() == -1) {
				unclassifiedNum++;
			}
		}
		logger.debug("Classified: " + (trianglesToClassifyToRegions - unclassifiedNum) + " out of: " + trianglesToClassifyToRegions);

		while (unclassifiedNum != 0) {
			
			Collections.sort(regions, new Comparator<Region>() {
				@Override
				public int compare(Region r1, Region r2) {
					Float area1 = r1.getAreaOfRegion();
					Float area2 = r2.getAreaOfRegion();
					int value = area1.compareTo(area2);
					return value * (-1);
				}
			});

			for (int i = 0 ; i < regions.size() ; ++i) {
				List<Triangle> neighborsAtRegion = regions.get(i).getOutsideBoundaryUnlabelled();
				for (Triangle t : neighborsAtRegion) {
					Edge edge = regions.get(i).getCommonEdge(t);
					// if common edge is not sharp add triangle to region
					if (edge == null) {
						continue;
					}
					if (!edge.getIsSharpEdge()) {
						regions.get(i).addTriangleToRegion(t);
					}
				}
//				regions.get(i).addTrianglesToRegion(neighborsAtRegion);
			}
			
			unclassifiedNum = 0;
			for (int i = 0 ; i < model.getTriangles().size() ; ++i) {
				if (model.getTriangles().get(i).getRegionLabel() == -1) {
					unclassifiedNum++;
				}
			}
		}
		
		Collections.sort(regions, new Comparator<Region>() {
			@Override
			public int compare(Region r1, Region r2) {
				Integer id1 = r1.getRegionId();
				Integer id2 = r2.getRegionId();
				return id1.compareTo(id2);
			}
		});
		
//		int contcom = 0;
//		for (int i = 0 ; i < model.getTriangles().size() ; ++i) {
//			Triangle t = model.getTriangles().get(i);
//			int cont = 0;
//			for (int j = 0 ; j < regions.size() ; ++j) {
//				if (regions.get(j).getTriangles().contains(t)) {
//					cont++;
//				}
//			}
//			if (cont != 1) {
//				contcom++;
//			}
////			logger.debug("cont is:" + cont);
//		}
//		
//		logger.debug("common triangles " + contcom);
//		logger.debug("Analyzing unvisited triangles ...");
//		if (regions.size() != 0) {
//			int conti = 0;
//			for (Triangle tr : model.getTriangles()) {
//				if (tr.getRegionLabel() == -1 && tr.isVisited == false) {
//					fillRegionCrackForTriangle(regions, tr);
//					conti++;
//				}
//			}
//			logger.debug("Cases visited: " + conti);
//		}
////		else {
////			for (Triangle tr : model.getTriangles()) {
////				if (tr.getRegionLabel() == -1) {
////					Region newRegion = new Region(contId,tr);
////					newRegion.buildUpRegion();
////					regions.add(newRegion);
////					contId++;
////				}
////			}
////		}
//		duration = System.currentTimeMillis() - duration;
//		logger.debug("Ended. Took: " + PrintUtil.prettyMillis(duration));
		for (int i = 0 ; i < regions.size() ; ++i) {
			regions.get(i).updateRegionNeighbors(regions);
//			System.out.println(regions.get(i));
		}
		
		// mark regions by coloring
		this.colorizeRegions(regions);
		
		// add regions processed to the CAD model
		model.setRegions(regions);
		
		logger.debug("Ended. Took: " + PrintUtil.prettyMillis(System.currentTimeMillis() - duration));
		logger.debug("Triangles classified: " + (trianglesToClassifyToRegions - unclassifiedNum) + " out of: " + trianglesToClassifyToRegions);
		logger.debug("CAD model has: " + regions.size() +" classified regions");
	}
	
	public void processRegionMerging() {
		logger.debug("Merging existent regions ...");
		long duration = System.currentTimeMillis();
		int regionsToMerge = model.getRegions().size();
		
		logger.debug("Regions to merge: " + regionsToMerge);
		HashMap<Integer, Region> regions = model.getRegionsMap();
		float[][] adjacencyMatrix = new float[regionsToMerge][regionsToMerge];

		
		// initialize adjacency matrix entries with maximum distance values
		for (int i = 0 ; i < adjacencyMatrix.length ; ++i) {
			for (int j = 0 ; j <= i ; ++j) {
				adjacencyMatrix[i][j] = PINF;
				adjacencyMatrix[j][i] = PINF;
			}
		}
		
		// compute for the first time adjacency distances
		computeAdjacencyDistances(adjacencyMatrix);
		
		float min = PINF;
		int rI = 0, rJ = 0, iteration = 0;
		if (UPPER_ITERATION_LIMIT_GROWING > model.getRegions().size()) {
			UPPER_ITERATION_LIMIT_GROWING = model.getRegions().size();
		}
		while (min > MIN_DISTANCE_THRESHOLD && iteration < UPPER_ITERATION_LIMIT_GROWING) {
			min = PINF;
			for (int i = 0 ; i < adjacencyMatrix.length ; ++i) {
				for (int j = 0 ; j <= i ; ++j) {
					if (min > adjacencyMatrix[i][j]) {
						min = adjacencyMatrix[i][j];
						rI = i;
						rJ = j;
					}
				}
			}
			if (min == PINF) {
				logger.debug("All regions are distinct. Stopping ...");
				break;
			}
			if (model.getRegions().size() == 1) {
				logger.debug("Only one regions left. Stopping ...");
				break;
			}
			if (regions.get(rJ) == null) {
				logger.debug("Stopping ...");
				break;
			}
			
			// merge regions rI and rJ (rJ into rI)
			mergeRegions(regions.get(rI), regions.get(rJ));
			iteration++;
			
			// disconnect merged region rJ from any other regions
			for (int i = 0 ; i < adjacencyMatrix.length ; ++i) {
				adjacencyMatrix[i][rJ] = PINF;
				adjacencyMatrix[rJ][i] = PINF;
			}
			
			// update neighbors of neighbors of merged regions
			for (Region nr : regions.get(rJ).getNeighbors()) {
				nr.updateRegionNeighbors(model.getRegions());				
			}
			
			// remove deprecated merged region rJ
			regions.remove(rJ);
			
			// re-compute necessary entries in the adjacency matrix (only for the existent merged region)
			computeAdjacencyDistances(adjacencyMatrix, regions.get(rI));
		}
		logger.debug("Exiting min distance value: " + min);
		logger.debug("Merged " + regionsToMerge + " regions into " + model.getRegions().size() + " regions in " + iteration + " iterations");
		logger.debug("Marking regions by coloring ...");
		colorizeRegions(model.getRegions());
		logger.debug("Ended. Took: " + PrintUtil.prettyMillis(System.currentTimeMillis() - duration));
	}
	
	public void updateCurvaturesBasedOnRegions(HashMap<Vertex,Curvature> curvatures) {
		HashMap<Integer, Region> regions = model.getRegionsMap();
		for (int i = 0 ; i < model.getTriangles().size() ; ++i) {
			Triangle t = model.getTriangles().get(i);
			for (int j = 0 ; j < t.getPosition().length ; ++j) {
				t.getPosition()[j].setClusterCurvatureVal(regions.get(t.getRegionLabel()).getCurvatureMinMaxOfRegion()[0], 
						regions.get(t.getRegionLabel()).getCurvatureMinMaxOfRegion()[1]);
				Curvature c = curvatures.get(t.getPosition()[j]);
				c.setCurvatureMax(t.getPosition()[j].getClusterCurvatureVal()[0]);
				c.setCurvatureMin(t.getPosition()[j].getClusterCurvatureVal()[1]);
				c.setCurvatureMinMax(t.getPosition()[j].getClusterCurvatureVal()[0] * t.getPosition()[j].getClusterCurvatureVal()[0]);
			}
		}
		CurvatureCalculation.setCurvatureHueSaturation(curvatures, model, 0.5f);
	}
	
//	/**
//	 * Crack filling recursive algorithm that assigns unlabelled 
//	 * triangles to the maximally spread neighboring region of itself 
//	 * or of its neighbors
//	 * 
//	 * @param regions
//	 * 			built regions of the model mesh
//	 * 
//	 * @param tr
//	 * 			triangle to be filled in one of the regions
//	 */
//	public void fillRegionCrackForTriangle(List<Region> regions, Triangle tr) {
//		float maxArea = 0.0f;
//		int regionLabelId = -1;
//		tr.isVisited = true;
//		for (Triangle neighbor : tr.getNeighbors()) {
//			if (neighbor.getRegionLabel() == -1 && neighbor.isVisited == false) {
//				fillRegionCrackForTriangle(regions, neighbor);
//			} 
//			if (neighbor.getRegionLabel() != -1) {
//				if (maxArea < regions.get(neighbor.getRegionLabel()).getAreaOfRegion()) {
//					maxArea = regions.get(neighbor.getRegionLabel()).getAreaOfRegion();
//					regionLabelId = regions.get(neighbor.getRegionLabel()).getRegionId();
//				}
//			}
//		}
//		if (regionLabelId == -1) {
//			
//			Region newRegion = new Region(regions.size(),tr);
//			regions.add(newRegion);
//		}
//		else {
//
//			regions.get(regionLabelId).addTriangleToRegion(tr);
//		}
//	}
}

