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

import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;

import edu.tum.cs.util.PrintUtil;
import edu.tum.cs.vis.model.uima.cas.MeshCas;
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
public class CasProcessing{

	/**
	 * Log4J logger
	 */
	private static Logger				logger = Logger.getLogger(CasProcessing.class);
	
	/**
	 * Mesh UIMA CAS to be processed
	 */
	private MeshCas						cas;
	
	/**
	 * Flag that shows if the classification of the curvatures has been done
	 */
	private boolean						KMeansCurvatureClassified = false;
	
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
	 * Default constructor of the class
	 */
	public CasProcessing() {
		this.cas = null;
	}
	
	/**
	 * Constructor of the CasProcessing class with the desired MeshCas
	 * 
	 * @param cas 
	 * 			cas to be analyzed
	 */
	public CasProcessing(MeshCas cas) {
		this.cas = cas;
	}
	
	/**
	 * Setter for the cas of an object instance
	 * 
	 * @param cas
	 * 			cas to be set
	 */
	public void setModel(MeshCas cas) {
		this.cas = cas;
	}
	
	/**
	 * Getter for the MeshCas of an object instance
	 * 
	 * @return cas
	 * 			MeshCas of the object instance
	 */
	public MeshCas getCas() {
		return this.cas;
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
		return KMeansCurvatureClassified;
	}
	
	/**
	 * KMeans algorithm implementation for vertex
	 * curvature classification
	 */
	public void KMeansVCClassification() {
		if (NUM_CLUSTERS >= cas.getModel().getVertices().size()* 0.025) {
			logger.debug("Number of vertices in the model smaller than number of clusters chosen");
			int temp = NUM_CLUSTERS;
			NUM_CLUSTERS = (int)(cas.getModel().getVertices().size() * 0.025);
			if (NUM_CLUSTERS < 5) {
				NUM_CLUSTERS = 5;
			}
			logger.debug("Number of clusters has been changed from " + temp + " to " + NUM_CLUSTERS);
		}
		Cluster[] clusters = new Cluster[NUM_CLUSTERS];
		
		// randomly initialize clusters with one element each
		List<Integer> dummyRnd = new ArrayList<Integer>();
		final List<Integer> pickedData = new ArrayList<Integer>();
		List<Vertex> verticesData = new ArrayList<Vertex>(cas.getModel().getVertices());
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
			clusters[i].updateCentroid(cas.getCurvatures());
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
			float distMin = this.distEuclid(cas.getCurvature(v), clusters[0]);
			int clusterIndex = 0;
			for (int i = 1 ; i < NUM_CLUSTERS ; ++i) {
				float distTmp = this.distEuclid(cas.getCurvature(v), clusters[i]);
				if (distTmp < distMin) {
					distMin = distTmp;
					clusterIndex = i;
				}
			}
			clusters[clusterIndex].addVertexToCluster(v);
			clusters[clusterIndex].updateCentroid(cas.getCurvatures());
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
			for (int i = 0 ; i < cas.getModel().getVertices().size() ; ++i) {
				float distMin = this.distEuclid(cas.getCurvature(cas.getModel().getVertices().get(i)), clusters[0]);
				int clusterIndex = 0;
				for (int j = 1 ; j < NUM_CLUSTERS ; ++j) {
					float distTmp = this.distEuclid(cas.getCurvature(cas.getModel().getVertices().get(i)), clusters[j]);
					if (distTmp < distMin) {
						distMin = distTmp;
						clusterIndex = j;
					}
				}
				//System.out.println(model.getVertices().get(i));
				if (cas.getModel().getVertices().get(i).getClusterLabel() != clusterIndex) {
					//System.out.println(model.getVertices().get(i).getClusterLabel());
					clusters[cas.getModel().getVertices().get(i).getClusterLabel()].removeVertexFromCluster(cas.getModel().getVertices().get(i));
					clusters[cas.getModel().getVertices().get(i).getClusterLabel()].updateCentroid(cas.getCurvatures());
					clusters[clusterIndex].addVertexToCluster(cas.getModel().getVertices().get(i));
					clusters[clusterIndex].updateCentroid(cas.getCurvatures());
					cas.getModel().getVertices().get(i).setClusterLabel(clusterIndex);
					cas.getModel().getVertices().get(i).setClusterCurvatureVal(clusters[clusterIndex].getCentroid()[0], clusters[clusterIndex].getCentroid()[1], clusters[clusterIndex].getCentroid()[2]);
					isRunning = true;
				}
			}
			iteration++;
		}
		this.KMeansCurvatureClassified = true;
		
		int classifiedVertices = 0;
		for (int i = 0 ; i < NUM_CLUSTERS ; ++i) {
			classifiedVertices += clusters[i].getVertices().size();
			for (int j = 0 ; j < clusters[i].getVertices().size() ; ++j) {
				Curvature c = cas.getCurvature(clusters[i].getVertices().get(j));
				clusters[i].getVertices().get(j).setClusterCurvatureVal(clusters[i].getCentroid()[0], clusters[i].getCentroid()[1], clusters[i].getCentroid()[2]);
				c.setCurvatureMin(clusters[i].getCentroid()[0]);
				c.setCurvatureMax(clusters[i].getCentroid()[1]);
				c.setCurvatureMinMax(clusters[i].getCentroid()[2]);
			}
		}
		
		// compute new coloring used for primitive fitting
		CurvatureCalculation.setCurvatureHueSaturation(cas.getCurvatures(), cas.getModel(), UtilityValues.CURV_SMOOTHING);
		
//		for (int i = 0 ; i < NUM_CLUSTERS ; ++i) {
//			System.out.println(clusters[i]);
//		}
		
		logger.debug("Classified " + classifiedVertices + " vertices out of " + cas.getModel().getVertices().size() + " into " + NUM_CLUSTERS + " clusters");
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
	
	private void computeAdjacencyDistancesOfRegions(final float[][] adjacencyMatrix, final Region r, final Region n) {
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
				return;
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
			float sharpEdgeBorder = 0.0f;
			for (int k = 0 ; k < commonEdges.size() ; ++k) {
				Edge e = commonEdges.get(k);
				if (e.getIsSharpEdge()) {
					sharpEdgeBorder += e.getEdgeValue().length();
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
			if (sharpEdgeBorder >= (0.4 * commonBorder)) {
				adjacencyMatrix[r.getRegionId()][n.getRegionId()] = PINF;
				adjacencyMatrix[n.getRegionId()][r.getRegionId()] = PINF;
				return;
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
	
	private void computeAdjacencyDistances(final float[][] adjacencyMatrix, final Region r) {
		List<Region> regions = cas.getModel().getRegions();
		for (int j = 0 ; j < regions.size() ; ++j) {
			Region n = regions.get(j);
			synchronized (n) {
				if (r.isNeighbor(n)) {
					computeAdjacencyDistancesOfRegions(adjacencyMatrix, r, n);;
				}
			}
		}
	}
	
	private void computeAdjacencyDistances(final float[][] adjacencyMatrix) {
		List<Region> regions = cas.getModel().getRegions();
		for (int i = 0 ; i < regions.size() ; ++i) {
			Region r = regions.get(i);
			for (int j = 0 ; j <= i ; ++j) {
				Region n = regions.get(j);
				synchronized (n) {
					if (r.isNeighbor(n)) {
						computeAdjacencyDistancesOfRegions(adjacencyMatrix, r, n);;
					}
				}
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
		// update curvature values weighted by area
		float newKMin = 0.0f, newKMax = 0.0f, newKMinKMax = 0.0f;
		newKMin = (r1.getAreaOfRegion() * r1.getCurvatureMinMaxOfRegion()[0] + r2.getAreaOfRegion() * r2.getCurvatureMinMaxOfRegion()[0]) 
				/ (r1.getAreaOfRegion() + r2.getAreaOfRegion());
		newKMax = (r1.getAreaOfRegion() * r1.getCurvatureMinMaxOfRegion()[1] + r2.getAreaOfRegion() * r2.getCurvatureMinMaxOfRegion()[1]) 
				/ (r1.getAreaOfRegion() + r2.getAreaOfRegion());
		newKMinKMax = (r1.getAreaOfRegion() * r1.getCurvatureMinMaxOfRegion()[2] + r2.getAreaOfRegion() * r2.getCurvatureMinMaxOfRegion()[2])
				/ (r1.getAreaOfRegion() + r2.getAreaOfRegion());
		r1.setCurvatureMin(newKMin);
		r1.setCurvatureMax(newKMax);
		r1.setCurvatureMinMax(newKMinKMax);
		
		r1.addTrianglesToRegion(r2.getTriangles());
		
		// update boundary
		r1.updateRegionBoundary();
		
		// remove r2 from the regions list
		cas.getModel().getRegions().remove(r2);
		
		// re-determine the neighboring relations for updated r1 and for deprecated r2
		for (Region reg : r2.getNeighbors()) {
			reg.updateRegionNeighbors(cas.getModel().getRegions());				
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
		int trianglesToClassifyToRegions = cas.getModel().getTriangles().size();
		List<Triangle> toClassify = new ArrayList<Triangle>();
		List<Region> regions = new ArrayList<Region>();
		int contId = 0;
		
		logger.debug("Triangles to classify: " + trianglesToClassifyToRegions);
		logger.debug("Starting with seed triangles ...");
		for (int i = 0 ; i < cas.getModel().getTriangles().size() ; ++i) {
			Triangle tr = cas.getModel().getTriangles().get(i);
			if (tr.updateIsSeedTriangle() && tr.getRegionLabel() == -1) {
				Region newRegion = new Region(contId,tr);
				newRegion.buildUpRegion();
				regions.add(newRegion);
				contId++;
			}
		}
		
		for (int i = 0 ; i < cas.getModel().getTriangles().size() ; ++i) {
			if (cas.getModel().getTriangles().get(i).getRegionLabel() == -1) {
				toClassify.add(cas.getModel().getTriangles().get(i));
			}
		}
		logger.debug("Classified: " + (trianglesToClassifyToRegions - toClassify.size()) + " out of: " + trianglesToClassifyToRegions);

		int cachedUnclassifiedNum;
		while (toClassify.size() != 0) {
			
			cachedUnclassifiedNum = toClassify.size();
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
						toClassify.remove(t);
					}
				}
			}
			
			logger.debug("Remained to classify: " + toClassify.size() + " triangles");
			
			// if not all triangles have been successfully classified into regions
			// continue and remedy this by looking at their immediate neighbors that
			// belong to a region; this solution works as the number of unclassified triangles
			// (if any at all) is very low (a couple of triangles)
			if (toClassify.size() == cachedUnclassifiedNum) {
				for (int i = 0 ; i < cas.getModel().getTriangles().size() ; ++i) {
					Triangle t = cas.getModel().getTriangles().get(i);
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
						toClassify.remove(t);
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
			regions.get(i).updateRegionBoundary();
		}
		
		for (int i = 0 ; i < regions.size() ; ++i) {
			regions.get(i).updateRegionNeighbors(regions);
//			System.out.println(regions.get(i));
		}
		
		// mark regions by coloring
//		this.colorizeRegions(regions);
		
		// add regions processed to the CAD model
		cas.getModel().setRegions(regions);
		
		logger.debug("Ended. Took: " + PrintUtil.prettyMillis(System.currentTimeMillis() - duration));
		logger.debug("Triangles classified: " + (trianglesToClassifyToRegions - toClassify.size()) + " out of: " + trianglesToClassifyToRegions);
		logger.debug("CAD model has: " + regions.size() +" classified regions");
	}
	
	public void processRegionMerging() {
		logger.debug("Merging existent regions ...");
		long duration = System.currentTimeMillis();
		int regionsToMerge = cas.getModel().getRegions().size();
		
		logger.debug("Regions to merge: " + regionsToMerge);
		HashMap<Integer, Region> regions = cas.getModel().getRegionsMap();
		float[][] adjacencyMatrix = new float[regionsToMerge][regionsToMerge];
		
		boolean[] isStillRegion = new boolean[regionsToMerge];
		for (int i = 0 ; i < isStillRegion.length ; ++i) {
			isStillRegion[i] = true;
		}

		// initialize adjacency matrix entries with maximum distance values
		for (int i = 0 ; i < adjacencyMatrix.length ; ++i) {
			for (int j = 0 ; j <= i ; ++j) {
				adjacencyMatrix[i][j] = PINF;
				adjacencyMatrix[j][i] = PINF;
			}
		}
		
		// compute for the first time adjacency distances
		computeAdjacencyDistances(adjacencyMatrix);
		
		float min = 0.0f;
		int rI = 0, rJ = 0, iteration = 0;
		ITERATION_LIMIT_GROWING = cas.getModel().getRegions().size();
//		if (ITERATION_LIMIT_GROWING > cas.getModel().getRegions().size()) {
//			ITERATION_LIMIT_GROWING = cas.getModel().getRegions().size();
//		}
		while (iteration < ITERATION_LIMIT_GROWING) {
			min = PINF;
			for (int i = 0 ; i < adjacencyMatrix.length ; ++i) {
				if (isStillRegion[i]) {
					for (int j = 0 ; j <= i ; ++j) {
						if (isStillRegion[j] && min > adjacencyMatrix[i][j]) {
							min = adjacencyMatrix[i][j];
							rI = i;
							rJ = j;
						}
					}
				}
			}
			if (min == PINF) {
				logger.debug("All regions are distinct. Stopping ...");
				break;
			}
			if (min > MIN_DISTANCE_THRESHOLD) {
				logger.debug("Minimum adjacency distance is bigger than threshold ("+ MIN_DISTANCE_THRESHOLD +"). Stopping ...");
				break;
			}
			if (regions.get(rJ) == null) {
				logger.debug("Stopping ...");
				break;
			}
			
			// merge regions rI and rJ: smaller region in bigger region
			if (regions.get(rI).getAreaOfRegion() < regions.get(rJ).getAreaOfRegion()) {
				int tmp;
				tmp = rI;
				rI = rJ;
				rJ = tmp;
			}
		    
			// merge and update neighboring relations for all neigbors of reg rJ
			mergeRegions(regions.get(rI), regions.get(rJ));
			iteration++;
			
			// disconnect merged region rJ from any other regions
			for (int i = 0 ; i < adjacencyMatrix.length ; ++i) {
				adjacencyMatrix[i][rJ] = PINF;
				adjacencyMatrix[rJ][i] = PINF;
			}

			// remove deprecated merged region rJ
			isStillRegion[rJ] = false;
			regions.remove(rJ);
			
			// re-compute necessary entries in the adjacency matrix (only for the existent merged region)
			computeAdjacencyDistances(adjacencyMatrix, regions.get(rI));
		}
		logger.debug("Exiting min distance value: " + min);
		logger.debug("Merged " + regionsToMerge + " regions into " + cas.getModel().getRegions().size() + " regions in " + iteration + " iterations");
		logger.debug("Ended. Took: " + PrintUtil.prettyMillis(System.currentTimeMillis() - duration));
		
		duration = System.currentTimeMillis();
		logger.debug("Updating vertices and triangles curvatures based on regions ...");
		this.updateCurvaturesBasedOnRegions();
		logger.debug("Curvatures values for the CAD model updated. Took: " + (PrintUtil.prettyMillis(System.currentTimeMillis() - duration)));
	}
	
	private void updateCurvaturesBasedOnRegions() {
		HashMap<Integer, Region> regions = cas.getModel().getRegionsMap();
		for (int i = 0 ; i < cas.getModel().getTriangles().size() ; ++i) {
			Triangle t = cas.getModel().getTriangles().get(i);
			for (int j = 0 ; j < t.getPosition().length ; ++j) {
				t.getPosition()[j].setClusterCurvatureVal(regions.get(t.getRegionLabel()).getCurvatureMinMaxOfRegion()[0], 
						regions.get(t.getRegionLabel()).getCurvatureMinMaxOfRegion()[1], regions.get(t.getRegionLabel()).getCurvatureMinMaxOfRegion()[2]);
				Curvature c = cas.getCurvature(t.getPosition()[j]);
				c.setCurvatureMax(t.getPosition()[j].getClusterCurvatureVal()[0]);
				c.setCurvatureMin(t.getPosition()[j].getClusterCurvatureVal()[1]);
				c.setCurvatureMinMax(t.getPosition()[j].getClusterCurvatureVal()[2]);
			}
			t.setCurvatureLevels(regions.get(t.getRegionLabel()).getCurvatureMinMaxOfRegion()[0], 
					regions.get(t.getRegionLabel()).getCurvatureMinMaxOfRegion()[1], 
					regions.get(t.getRegionLabel()).getCurvatureMinMaxOfRegion()[2]);
		}
		CurvatureCalculation.setCurvatureHueSaturation(cas.getCurvatures(), cas.getModel(), UtilityValues.CURV_SMOOTHING);
	}
}

