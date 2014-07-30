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
import java.util.List;
import java.util.ArrayList;

import edu.tum.cs.util.PrintUtil;
import edu.tum.cs.vis.model.Model;
import edu.tum.cs.vis.model.util.Appearance;
import edu.tum.cs.vis.model.util.Cluster;
import edu.tum.cs.vis.model.util.Curvature;
import edu.tum.cs.vis.model.util.Edge;
import edu.tum.cs.vis.model.util.Region;
import edu.tum.cs.vis.model.util.Triangle;
import edu.tum.cs.vis.model.util.UtilityValues;
import edu.tum.cs.vis.model.util.Vertex;


/**
 * Class that implements the processing blocks used to
 * render better segmentation results of the models 
 * analyzed. Blocks can be either pre- or post-processing
 * units in the computation flow applied to the model
 * 
 * @author Andrei Stoica
 */
public class ModelProcessing{

	private static Logger				logger = Logger.getLogger(ModelProcessing.class);
	/**
	 * Model to be processed
	 */
	protected Model						model;
	
	/**
	 * Flag that shows if the classification of the curvatures has been done
	 */
	private boolean						modelKMeansCurvatureClassification = false;
	
	/**
	 * Defines the number of clusters used for classifying the vertices
	 */	
	private static int					NUM_CLUSTERS = UtilityValues.NUM_CLUSTERS;
	
	/**
	 * Defines upper iteration limit for the KMeans algorithm
	 */
	private final static int			ITERATION_LIMIT = UtilityValues.ITERATIONS_LIM;
	
	/**
	 * Defines upper iteration limit for the region merging
	 */
	private static int					ITERATION_LIMIT_GROWING = UtilityValues.ITERATIONS_LIM_REGION_MERGING;
	
	/**
	 * Defines the minimal area for small regions merging step
	 */
	private final static float			AREA_MIN_LIMIT = UtilityValues.AREA_LIM;
	
	/**
	 * Define minimum distance threshold for the region merging which stops the merging
	 */
	private final static float			MIN_DISTANCE_THRESHOLD = UtilityValues.MIN_DISTANCE_TOL;
	
	/**
	 * Defines the value of the area weighting on the distance calculation
	 * of the merging process
	 */
	private final static float			EPSILON = UtilityValues.EPSILON;
	
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
	
	/**
	 * Gets number of clusters defined
	 */
	public int getNumOfClusters() {
		return NUM_CLUSTERS;
	}
	
	/**
	 * Getter for the K Means classification of curvatures
	 */
	public boolean isKMeansCurvatureClassified() {
		return modelKMeansCurvatureClassification;
	}
	
	/**
	 * KMeans algorithm implementation for vertex
	 * curvature classification
	 */
	public void KMeansVCClassification(HashMap<Vertex,Curvature> curvatures) {
		if (NUM_CLUSTERS >= model.getVertices().size() / 10) {
			logger.debug("Number of vertices in the model smaller than number of clusters chosen");
			int temp = NUM_CLUSTERS;
			NUM_CLUSTERS = temp / 5;
			logger.debug("Number of clusters has been reduced from " +
			temp + " to " + NUM_CLUSTERS);
		}
		Cluster[] clusters = new Cluster[NUM_CLUSTERS];
		
		// randomly initialize clusters with one element each
		List<Integer> dummyRnd = new ArrayList<Integer>();
		final List<Integer> pickedData = new ArrayList<Integer>();
		List<Vertex> verticesData = new ArrayList<Vertex>(model.getVertices());
		for (int i = 0 ; i < verticesData.size() ; ++i) {
			dummyRnd.add(i);
		}
		Collections.shuffle(dummyRnd);
		pickedData.addAll(dummyRnd.subList(0, NUM_CLUSTERS));
		for (int i = 0 ; i < NUM_CLUSTERS ; ++i) {
			clusters[i] = new Cluster(i);
			// System.out.println(clusters[i].getLabelId());
			Collections.shuffle(dummyRnd);
			// store random picked element to cluster
			Vertex v = verticesData.get(pickedData.get(i));
			clusters[i].addVertexToCluster(v);
			clusters[i].updateCentroid(curvatures);
			v.setClusterLabel(clusters[i].getLabelId());
			v.setClusterCurvatureVal(clusters[i].getCentroid()[0],clusters[i].getCentroid()[1],clusters[i].getCentroid()[2]);
		}
		
		// remove randomly picked elements
		for (int i = 0 ; i < NUM_CLUSTERS ; ++i) {
			verticesData.remove(clusters[i].getVertices().get(0));
		}
		
		// add all remaining vertices to closest cluster according to Euclidean
		// distance applied on the curvature space
		while (!verticesData.isEmpty()) {
			Vertex v = verticesData.remove(0);
			float distMin = this.distEuclid(curvatures.get(v), clusters[0]);
			int clusterIndex = 0;
			for (int i = 1 ; i < NUM_CLUSTERS ; ++i) {
				float distTmp = this.distEuclid(curvatures.get(v), clusters[i]);
				if (distTmp < distMin) {
					distMin = distTmp;
					clusterIndex = i;
				}
			}
			clusters[clusterIndex].addVertexToCluster(v);
			clusters[clusterIndex].updateCentroid(curvatures);
			v.setClusterLabel(clusters[clusterIndex].getLabelId());
			v.setClusterCurvatureVal(clusters[clusterIndex].getCentroid()[0],clusters[clusterIndex].getCentroid()[1],clusters[clusterIndex].getCentroid()[2]);
		}
		
//		for (int i = 0 ; i < NUM_CLUSTERS ; ++i) {
//			System.out.println("Cluster Id " + clusters[i].getLabelId());
//			System.out.println("Vertices: " + clusters[i].getVertices());
//		}
		
		// iterate until no significant changes are visible or iteration
		// limit is exceeded
		boolean isRunning = true;
		int iteration = 0;
		while (isRunning && iteration < ITERATION_LIMIT) {
			isRunning = false;
			for (int i = 0 ; i < model.getVertices().size() ; ++i) {
				float distMin = this.distEuclid(curvatures.get(model.getVertices().get(i)), clusters[0]);
				int clusterIndex = 0;
				for (int j = 1 ; j < NUM_CLUSTERS ; ++j) {
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
					model.getVertices().get(i).setClusterCurvatureVal(clusters[clusterIndex].getCentroid()[0], clusters[clusterIndex].getCentroid()[1], clusters[clusterIndex].getCentroid()[2]);
					isRunning = true;
				}
			}
			iteration++;
		}
		this.modelKMeansCurvatureClassification = true;
		
		int classifiedVertices = 0;
		for (int i = 0 ; i < NUM_CLUSTERS ; ++i) {
			classifiedVertices += clusters[i].getVertices().size();
			for (int j = 0 ; j < clusters[i].getVertices().size() ; ++j) {
				Curvature c = curvatures.get(clusters[i].getVertices().get(j));
				clusters[i].getVertices().get(j).setClusterCurvatureVal(clusters[i].getCentroid()[0], clusters[i].getCentroid()[1], clusters[i].getCentroid()[2]);
				c.setCurvatureMin(clusters[i].getCentroid()[0]);
				c.setCurvatureMax(clusters[i].getCentroid()[1]);
				c.setCurvatureMinMax(clusters[i].getCentroid()[2]);
			}
		}
		
		// compute new coloring used for primitive fitting
		CurvatureCalculation.setCurvatureHueSaturation(curvatures, model, 1f);
		
//		for (int i = 0 ; i < NUM_CLUSTERS ; ++i) {
//			System.out.println(clusters[i]);
//		}
		
		logger.debug("Classified " + classifiedVertices + " vertices out of " + model.getVertices().size() + " into " + NUM_CLUSTERS + " clusters");
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
				+ Math.pow(c.getCurvatureMax() - cluster.getCentroid()[1] , 2) 
				+ Math.pow(c.getCurvatureMinMax() - cluster.getCentroid()[2], 2));
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
				boolean detachedEdges = true;
				float[] curvatureDistanceBoundaryLinked = new float[2];
				curvatureDistanceBoundaryLinked[0] = curvatureDistanceBoundary[0];
				curvatureDistanceBoundaryLinked[1] = curvatureDistanceBoundary[1];
				int contLinked = cont;
				for (Vertex v : edgesVerticesAdjacency.keySet()) {
					if (edgesVerticesAdjacency.get(v) == 1) {
						curvatureDistanceBoundary[0] += v.getClusterCurvatureVal()[0];
						curvatureDistanceBoundary[1] += v.getClusterCurvatureVal()[1];
						cont++;
					}
					if (edgesVerticesAdjacency.get(v) == 2) {
						curvatureDistanceBoundaryLinked[0] += v.getClusterCurvatureVal()[0];
						curvatureDistanceBoundaryLinked[1] += v.getClusterCurvatureVal()[1];
						contLinked++;
						detachedEdges=false;
					}
				}
				if (detachedEdges) {
					curvatureDistanceBoundary[0] /= (float)cont;
					curvatureDistanceBoundary[1] /= (float)cont;
				}
				else {
					curvatureDistanceBoundary[0] = curvatureDistanceBoundaryLinked[0] / (float)contLinked;
					curvatureDistanceBoundary[1] = curvatureDistanceBoundaryLinked[1] / (float)contLinked;
				}
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
					boolean detachedEdges = true;
					float[] curvatureDistanceBoundaryLinked = new float[2];
					curvatureDistanceBoundaryLinked[0] = curvatureDistanceBoundary[0];
					curvatureDistanceBoundaryLinked[1] = curvatureDistanceBoundary[1];
					int contLinked = cont;
					for (Vertex v : edgesVerticesAdjacency.keySet()) {
						if (edgesVerticesAdjacency.get(v) == 1) {
							curvatureDistanceBoundary[0] += v.getClusterCurvatureVal()[0];
							curvatureDistanceBoundary[1] += v.getClusterCurvatureVal()[1];
							cont++;
						}
						if (edgesVerticesAdjacency.get(v) == 2) {
							curvatureDistanceBoundaryLinked[0] += v.getClusterCurvatureVal()[0];
							curvatureDistanceBoundaryLinked[1] += v.getClusterCurvatureVal()[1];
							contLinked++;
							detachedEdges=false;
						}
					}
					if (detachedEdges) {
						curvatureDistanceBoundary[0] /= (float)cont;
						curvatureDistanceBoundary[1] /= (float)cont;
					}
					else {
						curvatureDistanceBoundary[0] = curvatureDistanceBoundaryLinked[0] / (float)contLinked;
						curvatureDistanceBoundary[1] = curvatureDistanceBoundaryLinked[1] / (float)contLinked;
					}
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

		int cachedUnclassifiedNum;
		while (unclassifiedNum != 0) {
			
			cachedUnclassifiedNum = unclassifiedNum;
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
				Triangle t = model.getTriangles().get(i);
				if (t.getRegionLabel() == -1) {
					unclassifiedNum++;
				}
			}
			logger.debug("Remained triangles to classify: " + unclassifiedNum);
			
			// if not all triangles have been successfully classified into regions
			// continue and remedy this by looking at their immediate neighbors that
			// belong to a region; this solution works as the number of unclassified triangles
			// (if any at all) is very low (a couple of triangles)
			if (unclassifiedNum == cachedUnclassifiedNum) {
				for (int i = 0 ; i < model.getTriangles().size() ; ++i) {
					Triangle t = model.getTriangles().get(i);
					if (t.getRegionLabel() == -1) {
						float maxArea = 0.0f;
						int index = -1;
						for (Triangle n : t.getNeighbors()) {
							if (n.getRegionLabel() != -1 && maxArea < n.getArea()) {
								maxArea = n.getArea();
								index = n.getRegionLabel();
							}
						}
						regions.get(index).addTriangleToRegion(t);
					}
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
		
//		boolean[] isStillRegion = new boolean[regionsToMerge];
//		for (int i = 0 ; i < isStillRegion.length ; ++i) {
//			isStillRegion[i] = true;
//		}

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
		if (ITERATION_LIMIT_GROWING > model.getRegions().size()) {
			ITERATION_LIMIT_GROWING = model.getRegions().size();
		}
		while (min > MIN_DISTANCE_THRESHOLD && iteration < ITERATION_LIMIT_GROWING) {
			min = PINF;
			for (int i = 0 ; i < adjacencyMatrix.length ; ++i) {
//				if (isStillRegion[i]) {
					for (int j = 0 ; j <= i ; ++j) {
						if (/*isStillRegion[j] && */min > adjacencyMatrix[i][j]) {
							min = adjacencyMatrix[i][j];
							rI = i;
							rJ = j;
						}
					}
//				}
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
//			isStillRegion[rJ] = false;
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
						regions.get(t.getRegionLabel()).getCurvatureMinMaxOfRegion()[1], regions.get(t.getRegionLabel()).getCurvatureMinMaxOfRegion()[2]);
				Curvature c = curvatures.get(t.getPosition()[j]);
				c.setCurvatureMax(t.getPosition()[j].getClusterCurvatureVal()[0]);
				c.setCurvatureMin(t.getPosition()[j].getClusterCurvatureVal()[1]);
				c.setCurvatureMinMax(t.getPosition()[j].getClusterCurvatureVal()[2]);
			}
		}
		CurvatureCalculation.setCurvatureHueSaturation(curvatures, model, 1.0f);
	}
}

