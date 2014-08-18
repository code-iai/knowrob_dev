/*******************************************************************************
 * Copyright (c) 2014 Andrei Stoica. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Andrei Stoica - initial API and implementation during Google Summer of Code 2014
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
import edu.tum.cs.vis.model.util.ThresholdsReasoning;
import edu.tum.cs.vis.model.util.Vertex;


/**
 * Class that implements the processing blocks used to filter the estimated curvatures.
 * The functionality implemented is directed to as pre-processing block before the primitive
 * fitting analyser ({@link src.edu.tum.cs.vis.model.uima.analyser.PrimitiveAnalyser}.
 * 
 * The filtering steps implemented are:
 * 		- KMeans clustering of the KMin, KMax and KMinMax values estimated by the curvature
 * calculation step
 * 		- Region grouping of triangles based on the clustered curvature values
 * 		- Region merging combining locally connected regions with similar curvature properties
 * into bigger regions
 * 		- Updater of the curvature properties based on the grown and merged regions
 * 
 * @author Andrei Stoica
 */
public class CasProcessing{

	/**
	 * Log4J logger
	 */
	private static final Logger		LOGGER = Logger.getLogger(CasProcessing.class);
	
	/**
	 * Mesh UIMA CAS to be processed
	 */
	private MeshCas					cas;
	
	/**
	 * Defines the number of clusters used for classifying the vertices
	 */	
	private int						NUM_CLUSTERS = ThresholdsReasoning.NUM_CLUSTERS;
		
	/**
	 * Defines the maximum number for the 32-bit floating point precision (+Inf)
	 */
	private final static float		PINF = Float.MAX_VALUE;
	
	/**
	 * Default constructor of the class. Initializes the MeshCas to null.
	 */
	public CasProcessing() {
		this.cas = null;
	}
	
	/**
	 * Constructor of the CasProcessing class based on the MeshCas to be analyzed
	 * 
	 * @param cas 
	 * 			MeshCas to be analyzed
	 */
	public CasProcessing(MeshCas cas) {
		this.cas = cas;
	}
	
	/**
	 * Sets the MeshCas of the instance if this was not previously set at
	 * initialization time
	 * 
	 * @param cas
	 * 			MeshCas to be set
	 */
	public void setModel(MeshCas cas) {
		this.cas = (this.cas != null && cas != null) ? cas : this.cas;
	}
	
	/**
	 * Gets the MeshCas of the instance
	 * 
	 * @return MeshCas of the instance
	 */
	public MeshCas getCas() {
		return this.cas;
	}
	
	/**
	 * Gets the number of clusters used to group similar curvature properties
	 * The defined default number of clusters is obtained only before the {@code process()}
	 * method is called on the instance, while the computed number of clusters
	 * can be obtained after the processing is done.
	 * 
	 * @return number of clusters used to group similar curvatures
	 */
	public int getNumOfClusters() {
		return this.NUM_CLUSTERS;
	}
	
	/**
	 * Groups vertices in clusters according to the similarity of their curvature
	 * properties. This is done based on the K-Means clustering algorithm and on
	 * the euclidean distance applied to the 3-length feature vector formed of the
	 * minimum curvature KMin, maximum curvature KMax and respectively KMinMax curvature
	 * property of the individual vertices of the model contained by the MeshCas.
	 * 
	 * The grouping creates clusters which assign new averaged values of the curvature
	 * properties to the similar vertices within them.
	 */
	private void kMeansCurvClassify() {
		if (this.cas == null) {
			return;
		}
		LOGGER.debug("Classifying vertices of CAD model by curvature using KMeans ...");
		long KMeansStartTime = System.currentTimeMillis();
		if (NUM_CLUSTERS >= cas.getModel().getVertices().size()* 0.025) {
			LOGGER.debug("Number of vertices in the model smaller than number of clusters chosen");
			int temp = NUM_CLUSTERS;
			NUM_CLUSTERS = (int)(cas.getModel().getVertices().size() * 0.025);
			if (NUM_CLUSTERS < 5) {
				NUM_CLUSTERS = 5;
			}
			LOGGER.debug("Number of clusters has been changed from " + temp + " to " + NUM_CLUSTERS);
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
			Collections.shuffle(dummyRnd);
			// store random picked element to cluster
			Vertex v = verticesData.get(pickedData.get(i));
			clusters[i].addVertexToCluster(v);
			clusters[i].updateCentroid(cas.getCurvatures());
			v.setClusterLabel(clusters[i].getLabelId());
			v.setClusterCurvatureVal(clusters[i].getCentroid()[0],clusters[i].getCentroid()[1],
					clusters[i].getCentroid()[2]);
		}
		
		// remove randomly picked elements
		for (int i = 0 ; i < NUM_CLUSTERS ; ++i) {
			verticesData.remove(clusters[i].getVertices().get(0));
		}
		
		// add all remaining vertices to closest cluster according to Euclidean
		// distance applied on the curvature space
		while (!verticesData.isEmpty()) {
			Vertex v = verticesData.remove(0);
			float distMin = distEuclid(cas.getCurvature(v), clusters[0]);
			int clusterIndex = 0;
			for (int i = 1 ; i < NUM_CLUSTERS ; ++i) {
				float distTmp = distEuclid(cas.getCurvature(v), clusters[i]);
				if (distTmp < distMin) {
					distMin = distTmp;
					clusterIndex = i;
				}
			}
			clusters[clusterIndex].addVertexToCluster(v);
			clusters[clusterIndex].updateCentroid(cas.getCurvatures());
			v.setClusterLabel(clusters[clusterIndex].getLabelId());
			v.setClusterCurvatureVal(clusters[clusterIndex].getCentroid()[0],clusters[clusterIndex].getCentroid()[1],
					clusters[clusterIndex].getCentroid()[2]);
		}
		
		// iterate until no significant changes are visible or iteration
		// limit is exceeded
		boolean isRunning = true;
		int iteration = 0;
		while (isRunning && iteration < ThresholdsReasoning.ITERATIONS_LIM) {
			isRunning = false;
			for (int i = 0 ; i < cas.getModel().getVertices().size() ; ++i) {
				float distMin = distEuclid(cas.getCurvature(cas.getModel().getVertices().get(i)), clusters[0]);
				int clusterIndex = 0;
				for (int j = 1 ; j < NUM_CLUSTERS ; ++j) {
					float distTmp = distEuclid(cas.getCurvature(cas.getModel().getVertices().get(i)), clusters[j]);
					if (distTmp < distMin) {
						distMin = distTmp;
						clusterIndex = j;
					}
				}
				if (cas.getModel().getVertices().get(i).getClusterLabel() != clusterIndex) {
					clusters[cas.getModel().getVertices().get(i).getClusterLabel()].removeVertexFromCluster(cas.getModel().getVertices().get(i));
					clusters[cas.getModel().getVertices().get(i).getClusterLabel()].updateCentroid(cas.getCurvatures());
					clusters[clusterIndex].addVertexToCluster(cas.getModel().getVertices().get(i));
					clusters[clusterIndex].updateCentroid(cas.getCurvatures());
					cas.getModel().getVertices().get(i).setClusterLabel(clusterIndex);
					cas.getModel().getVertices().get(i).setClusterCurvatureVal(clusters[clusterIndex].getCentroid()[0], 
							clusters[clusterIndex].getCentroid()[1], clusters[clusterIndex].getCentroid()[2]);
					isRunning = true;
				}
			}
			iteration++;
		}
		
		int classifiedVertices = 0;
		for (int i = 0 ; i < NUM_CLUSTERS ; ++i) {
			classifiedVertices += clusters[i].getVertices().size();
			for (int j = 0 ; j < clusters[i].getVertices().size() ; ++j) {
				Curvature c = cas.getCurvature(clusters[i].getVertices().get(j));
				clusters[i].getVertices().get(j).setClusterCurvatureVal(clusters[i].getCentroid()[0], 
						clusters[i].getCentroid()[1], clusters[i].getCentroid()[2]);
				c.setCurvatureMin(clusters[i].getCentroid()[0]);
				c.setCurvatureMax(clusters[i].getCentroid()[1]);
				c.setCurvatureMinMax(clusters[i].getCentroid()[2]);
			}
		}
		
		// compute new coloring used for primitive fitting
		CurvatureCalculation.setCurvatureHueSaturation(cas.getCurvatures(), cas.getModel(), ThresholdsReasoning.CURV_SMOOTHING);
		
		LOGGER.debug("Classified " + classifiedVertices + " vertices out of "+
		cas.getModel().getVertices().size() + " into " + NUM_CLUSTERS + " clusters");
		long KMeansDuration = System.currentTimeMillis() - KMeansStartTime;
		LOGGER.debug("Ended. Took: " + PrintUtil.prettyMillis(KMeansDuration));
	}
	
	/**
	 * Computes the Euclidean distance between the center point of a cluster
	 * and the curvature of a single vertex based on the KMin, KMax and KMinKMax
	 * properties.
	 *  
	 * @param c
	 * 			curvature associated with selected vertex
	 * @param cluster
	 * 			cluster considered
	 * @return Euclidean distance between <b>c</b> and <b>cluster</b> centroid
	 */
	private static float distEuclid(final Curvature c, final Cluster cluster) {
		return (float)Math.sqrt(Math.pow(c.getCurvatureMin() - cluster.getCentroid()[0] , 2) 
				+ Math.pow(c.getCurvatureMax() - cluster.getCentroid()[1] , 2) 
				+ Math.pow(c.getCurvatureMinMax() - cluster.getCentroid()[2], 2));
	}
	
	/**
	 * Computes the adjacency distance between two neighboring regions. This is done based on the
	 * paper "A new CAD mesh segmentation method, based on curvature tensor analysis",
	 * Guillaume Lavoue, Florent Dupont, Atilla Baskurt, Computer-Aided Design 37 (2005) 
	 * 975–987, and it takes into account only the KMin and KMax curvature values.The result 
	 * is stored into the corresponding entry of the {@code adjacencyMatrix} associated.
	 * 
	 * Two similar curvature-wise regions tend to have a low distance associated, while two distinct
	 * ones tend to have a high distance associated. Two regions that are bounded at some boundary
	 * by sharp edges are completely distinct and they score {@link java.lang.Float.Float.MAX_VALUE}
	 * 
	 * @param adjacencyMatrix
	 * 			symmetric matrix that contains all distances between regions
	 * @param r
	 * 			first region
	 * @param n
	 * 			second region
	 * @see <a href="http://dl.acm.org/citation.cfm?id=1649921">"A new CAD mesh segmentation method, based on curvature tensor analysis"</a>
	 */
	private static void computeAdjacencyDistancesOfRegions(final float[][] adjacencyMatrix, final Region r, final Region n) {
		float[] curvatureDistanceBoundary = {0.0f, 0.0f};
		float curvatureDistance = 0.0f;
		float commonBorder = 0.0f;
		float neighborhoodDistance = 0.0f;
		float sizeDistance = 0.0f;
		float distance = 0.0f;
		List<Edge> commonEdges = r.getCommonEdges(n);
		if (commonEdges.size() == 1) {
			// if common edge is sharp cannot merge the two regions
			if (commonEdges.get(0).isSharpEdge()) {
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
				if (e.isSharpEdge()) {
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
		// compute curvature distance: cD = ||c_reg_r - cDB|| + ||c_reg_n - cDB|| + ||c_reg_r - c_reg_n||
		curvatureDistance = (float)((Math.sqrt(Math.pow(r.getCurvatureMinMaxOfRegion()[0] - curvatureDistanceBoundary[0], 2) 
				+ Math.pow(r.getCurvatureMinMaxOfRegion()[1] - curvatureDistanceBoundary[1], 2))) + 
				(Math.sqrt(Math.pow(n.getCurvatureMinMaxOfRegion()[0] - curvatureDistanceBoundary[0], 2) 
				+ Math.pow(n.getCurvatureMinMaxOfRegion()[1] - curvatureDistanceBoundary[1], 2))) + 
				(Math.sqrt(Math.pow(n.getCurvatureMinMaxOfRegion()[0] - r.getCurvatureMinMaxOfRegion()[0], 2) 
				+ Math.pow(n.getCurvatureMinMaxOfRegion()[1] - r.getCurvatureMinMaxOfRegion()[1], 2))));
		if (curvatureDistance == 0.0f) {
			curvatureDistance = ThresholdsReasoning.EPSILON;
		}
		// compute neighborhood distance: nD = min(per_r,per_n) / per_common;
		neighborhoodDistance = Math.min(r.getPerimeterOfRegion(), n.getPerimeterOfRegion()) / commonBorder;
		// size distance: sD = 
		if (r.getAreaOfRegion() < ThresholdsReasoning.MIN_AREA || n.getAreaOfRegion() < ThresholdsReasoning.MIN_AREA) {
			sizeDistance = ThresholdsReasoning.EPSILON;
		}
		else {
			sizeDistance = 1.0f;
		}
		distance = curvatureDistance * neighborhoodDistance * sizeDistance;
		// add distance to region adjacency matrix (which is symmetric)
		adjacencyMatrix[r.getRegionId()][n.getRegionId()] = distance;
		adjacencyMatrix[n.getRegionId()][r.getRegionId()] = distance;
	}
	
	/**
	 * Computes the adjacency distances of a region with respect to all 
	 * its available neighbors according to {@link computeAdjacencyDistancesOfRegions(float, Region, Region)}.
	 * The results are stored in the associated places of the matrix {@code adjacencyMatrix} according to the
	 * region pairs under calculation.
	 * 
	 * @param adjacencyMatrix
	 * 				matrix of adjacency distances over all the regions
	 * @param r
	 * 				region to compute the adjacency distances for
	 */
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
	
	/**
	 * Computes the adjacency distances of all the regions existent in the CAD model with respect to their 
	 * subsequent neighbors according to {@link computeAdjacencyDistancesOfRegions(float, Region, Region)}.
	 * The results are stored in the associated places of the matrix {@code adjacencyMatrix} according to the
	 * region pairs under calculation.
	 * 
	 * @param adjacencyMatrix
	 */
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
	 * Merges regions {@code r1} and {@code r2} by integrating {@code r2} into {@code r1} and
	 * recalculating all curvature properties and neighbors of the {@code r1} instance, while
	 * still keeping the {@code r2} instance in the list of regions existent in the model.
	 * The deletion of the 
	 * 
	 * @param r1
	 * 				region to merge into
	 * @param r2
	 * 				region to merge
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
	 * Grows regions (connected triangles) on the surface of the CAD mesh based on the
	 * clustered vertices of the model. The region growing mechanism is based on the
	 * paper "A new CAD mesh segmentation method, based on curvature tensor analysis",
	 * Guillaume Lavoue, Florent Dupont, Atilla Baskurt, Computer-Aided Design 37 (2005),
	 * 975–987. The regions have therefore triangles with the same curvature properties,
	 * grouping connected components together on the mesh. The grown regions are added 
	 * to a list of regions which is then added to the MeshCas instance.
	 * 
	 * @see <a href="http://dl.acm.org/citation.cfm?id=1649921">"A new CAD mesh segmentation method, based on curvature tensor analysis"</a>
	 */
	private void growRegions() {
		if (this.cas == null) {
			// if no cas to analyze, exit
			return;
		}
		LOGGER.debug("Building up regions ...");
		long duration = System.currentTimeMillis();
		int trianglesToClassifyToRegions = cas.getModel().getTriangles().size();
		List<Triangle> toClassify = new ArrayList<Triangle>();
		List<Region> regions = new ArrayList<Region>();
		int contId = 0;
		
		LOGGER.debug("Triangles to classify: " + trianglesToClassifyToRegions);
		LOGGER.debug("Starting with seed triangles ...");
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
		LOGGER.debug("Classified: " + (trianglesToClassifyToRegions - toClassify.size()) + " out of: " + trianglesToClassifyToRegions);

		// crack filling algorithm needed to add not grouped triangles to
		// already existent regions
		int cachedUnclassifiedNum;
		while (toClassify.size() != 0) {
			
			cachedUnclassifiedNum = toClassify.size();
			// weight filling by area of regions
			Collections.sort(regions, new Comparator<Region>() {
				@Override
				public int compare(Region r1, Region r2) {
					Float area1 = r1.getAreaOfRegion();
					Float area2 = r2.getAreaOfRegion();
					int value = area1.compareTo(area2);
					return value * (-1);
				}
			});
			
			// fill the cracks in the region segmentation of the model
			for (int i = 0 ; i < regions.size() ; ++i) {
				List<Triangle> neighborsAtRegion = regions.get(i).getOutsideBoundaryUnlabelled();
				for (Triangle t : neighborsAtRegion) {
					Edge edge = regions.get(i).getCommonEdge(t);
					// if common edge is not sharp add triangle to region
					if (edge == null) {
						continue;
					}
					if (!edge.isSharpEdge()) {
						regions.get(i).addTriangleToRegion(t);
						toClassify.remove(t);
					}
				}
			}
			
			LOGGER.debug("Remained to classify: " + toClassify.size() + " triangles");
			
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
		// sort regions back based on index
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
		}
		
		// add regions processed to the CAD model
		cas.getModel().setRegions(regions);
		
		LOGGER.debug("Ended. Took: " + PrintUtil.prettyMillis(System.currentTimeMillis() - duration));
		LOGGER.debug("Triangles classified: " + (trianglesToClassifyToRegions - toClassify.size()) + " out of: " + trianglesToClassifyToRegions);
		LOGGER.debug("CAD model has: " + regions.size() +" classified regions");
	}
	
	/**
	 * Merges together already existent regions which part the mesh surface of the
	 * CAD model. The merging process is implemented based on the paper "A new CAD mesh 
	 * segmentation method, based on curvature tensor analysis", Guillaume Lavoue, 
	 * Florent Dupont, Atilla Baskurt, Computer-Aided Design 37 (2005) 975–987. The region
	 * content of the MeshCas is modified as the deprecated merged regions are removed, while
	 * the regions which are merged into are updated. The merging is continued until 
	 * either all remaining regions are completely distinct (all entries are at {@link java.lang.Float.Float.MAX_VALUE}
	 * and the regions are therefore separated by sharp edges only) or a maximum disimilarity distance threshold 
	 * is reached (threshold set in {@link src.edu.tum.cs.util.TresholdReasoning}). After the merging is completed
	 * the updated curvature information is used in order to compute the curvature properties' values needed for the
	 * primitive fitting part of the segmentation.
	 * 
	 * @see <a href="http://dl.acm.org/citation.cfm?id=1649921">"A new CAD mesh segmentation method, based on curvature tensor analysis"</a>
	 */
	private void mergeRegions() {
		if (this.cas == null) {
			// if no cas to analyze exit
			return;
		}
		LOGGER.debug("Merging existent regions ...");
		long duration = System.currentTimeMillis();
		int regionsToMerge = cas.getModel().getRegions().size();
		int iterationLimit = 0;
		
		LOGGER.debug("Regions to merge: " + regionsToMerge);
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
		iterationLimit = cas.getModel().getRegions().size();
		while (iteration < iterationLimit) {
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
				LOGGER.debug("All regions are distinct. Stopping ...");
				break;
			}
			if (min > ThresholdsReasoning.MIN_DISTANCE_TOL) {
				LOGGER.debug("Minimum adjacency distance is bigger than threshold ("+ ThresholdsReasoning.MIN_DISTANCE_TOL +"). Stopping ...");
				break;
			}
			if (regions.get(rJ) == null) {
				LOGGER.debug("Stopping ...");
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
		LOGGER.debug("Exiting min distance value: " + min);
		LOGGER.debug("Merged " + regionsToMerge + " regions into " + cas.getModel().getRegions().size() + " regions in " + iteration + " iterations");
		LOGGER.debug("Ended. Took: " + PrintUtil.prettyMillis(System.currentTimeMillis() - duration));
	}
	
	/**
	 * Updates the curvature properties' values based on the existent regions
	 * inside the MeshCas. The calculated values are used into the primitive 
	 * classification part of the segmentation and they are stored as properties
	 * of individual vertices.
	 */
	private void updateCurvaturesProperties() {
		if (this.cas == null || this.cas.getTriangles() == null) {
			return;
		}
		long duration = System.currentTimeMillis();
		LOGGER.debug("Updating vertices and triangles curvatures based on regions ...");
		HashMap<Integer, Region> regions = cas.getModel().getRegionsMap();
		for (int i = 0 ; i < cas.getModel().getTriangles().size() ; ++i) {
			if (cas.getModel().getTriangles().get(i) == null) {
				continue;
			}
			Triangle t = cas.getModel().getTriangles().get(i);
			for (int j = 0 ; j < t.getPosition().length ; ++j) {
				t.getPosition()[j].setClusterCurvatureVal(regions.get(t.getRegionLabel()).getCurvatureMinMaxOfRegion()[0], 
						regions.get(t.getRegionLabel()).getCurvatureMinMaxOfRegion()[1], 
						regions.get(t.getRegionLabel()).getCurvatureMinMaxOfRegion()[2]);
				if (cas.getCurvature(t.getPosition()[j]) == null) {
					continue;
				}
				Curvature c = cas.getCurvature(t.getPosition()[j]);
				c.setCurvatureMax(t.getPosition()[j].getClusterCurvatureVal()[0]);
				c.setCurvatureMin(t.getPosition()[j].getClusterCurvatureVal()[1]);
				c.setCurvatureMinMax(t.getPosition()[j].getClusterCurvatureVal()[2]);
			}
			t.setCurvatureLevels(regions.get(t.getRegionLabel()).getCurvatureMinMaxOfRegion()[0], 
					regions.get(t.getRegionLabel()).getCurvatureMinMaxOfRegion()[1], 
					regions.get(t.getRegionLabel()).getCurvatureMinMaxOfRegion()[2]);
		}
		CurvatureCalculation.setCurvatureHueSaturation(cas.getCurvatures(), cas.getModel(), ThresholdsReasoning.CURV_SMOOTHING);
		LOGGER.debug("Curvatures values for the CAD model updated. Took: " + 
		(PrintUtil.prettyMillis(System.currentTimeMillis() - duration)));
	}
	
	/**
	 * Runs all the processing filters on the MeshCas model of the instance.
	 * The changes made by filtering touch the curvature of vertices and triangles
	 * of the model as well as the hue and saturation properties of a curvature object.
	 * In addition, the clustering information of vertices is updated and regions of
	 * triangles with similar curvature properties are created. Triangles or vertices 
	 * are not added or removed along the processing stages. All the changes are stored
	 * in the utilities Vertex, Triangle, Curvature and then stored upstream in the 
	 * model object of the MeshCas.
	 */
	public void process() {
		if (this.cas == null) {
			LOGGER.debug("No MeshCas to process. Received MeshCas: null");
			return;
		}
		// perform the KMeans clustering for the vertices curvature properties
		this.kMeansCurvClassify();
		
		// create regions based on the clustering and then merge them
		this.growRegions();
		this.mergeRegions();
		
		// update curvature properties now
		this.updateCurvaturesProperties();
	}
}

