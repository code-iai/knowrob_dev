/*******************************************************************************
 * Copyright (c) 2014 Andrei Stoica. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Andrei Stoica - initial API and implementation during
 * 								 Google Summer of Code 2014
 ******************************************************************************/
package edu.tum.cs.vis.model.uima.analyser;

import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.Callable;

import org.apache.log4j.Logger;

import edu.tum.cs.util.PrintUtil;
import edu.tum.cs.vis.model.util.Appearance;
import edu.tum.cs.vis.model.util.Edge;
import edu.tum.cs.vis.model.util.Triangle;
import edu.tum.cs.vis.model.util.ThresholdsReasoning;
import edu.tum.cs.vis.model.util.Vertex;
import edu.tum.cs.vis.model.uima.cas.MeshCas;
import edu.tum.cs.ias.knowrob.utils.ThreadPool;

/**
 * Analyzer for detecting the sharp edges of the CAD model. A sharp edge is
 * an edge which separates two surfaces whose normals form an angle bigger
 * than a modifiable threshold. The sharp edges play a role in building up 
 * regions of triangles with similar curvature properties. Such regions are
 * bounded by the sharp edges, since at the level of these edges the curvature
 * levels are not well defined
 * 
 * The analyser can add / remove triangles and vertices to the models. It is 
 * an implementation based on the paper "A new CAD mesh segmentation method, based on
 * curvature tensor analysis", Guillaume Lavoue, Florent Dupont, Atilla Baskurt, 
 * Computer-Aided Design 37 (2005) 975–987.
 * 
 * @author Andrei Stoica 
 */
public class EdgeAnalyser extends MeshAnalyser {
	
	/**
	 * Log4j logger
	 */
	private static final Logger		LOGGER = Logger.getLogger(EdgeAnalyser.class);
	
	/**
	 * Number of currently processed triangles. Used for progress status
	 */
	private final AtomicInteger		trianglesProcessed = new AtomicInteger(0);
	
	/**
	 * MeshCas containing the model under analysis
	 */
	private MeshCas					cas;
	
	/**
	 * List of all vertices of the CAD model
	 */
	private List<Vertex> 			allVertices;
	
	/**
	 * List of all triangles of the CAD model
	 */
	private List<Triangle>			allTriangles;
	
	/**
	 * Number of added triangles to the model
	 */
	private int 					numAddedTriangles;
	
	/**
	 * Gets the number of added triangles to the model
	 * 
	 * @return number of added triangles to the model
	 */
	public int getNumAddedTriangles() {
		return numAddedTriangles;
	}
	
	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.analyser.MeshAnalyser#getLogger()
	 */
	@Override
	public Logger getLogger() {
		return LOGGER;
	}
	
	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.analyser.MeshAnalyser#getName()
	 */
	@Override
	public String getName() {
		return "Edge";
	}
	
	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.analyser.MeshAnalyser#processStart(edu.tum.cs.vis.model.uima.cas.MeshCas)
	 */
	@Override
	public void processStart(MeshCas newCas) {
		if (newCas == null || newCas.getModel() == null || 
				newCas.getModel().getTriangles() == null || 
				newCas.getModel().getVertices() == null) {
			return;
		}
		cas = newCas;
		allTriangles = cas.getModel().getTriangles();
		allVertices = cas.getModel().getVertices();
		trianglesProcessed.set(0);
		
		this.numAddedTriangles = allTriangles.size();
		
		// perform the sharp edge detection for individual triangles
		List<Callable<Void>> threads = new LinkedList<Callable<Void>>();

		final int interval = 500;

		for (int start = 0; start < allTriangles.size(); start += interval) {
			final int st = start;
			threads.add(new Callable<Void>() {

				@Override
				public Void call() throws Exception {
					int end = Math.min(st + interval, allTriangles.size());
					for (int i = st; i < end; i++) {
						sharpEdgeDetectionForTriangle(allTriangles.get(i));
					}
					return null;
				}

			});
		}
		ThreadPool.executeInPool(threads);
		
		// after sharp edges are detected build up 3 triangles
		// in place of the sharp triangles and remove sharp
		// triangles from the model
		trianglesProcessed.set(0);
		for (int i = 0 ; i < allTriangles.size() ; ++i) {
			Triangle t = allTriangles.get(i);
			trianglesProcessed.incrementAndGet();
			if (t.checkIsSharpTriangle()) {
				addTrianglesToModel(t);
				trianglesProcessed.decrementAndGet();
			}
		}
		
		// reload the list of vertices
		cas.getModel().reloadVertexList();
	
		this.numAddedTriangles = allTriangles.size() - this.numAddedTriangles;
		if (this.numAddedTriangles != 0) {
			LOGGER.debug("Added " + this.numAddedTriangles + " triangles to the model");
		}
		else {
			LOGGER.debug("No triangles added to the model");
		}
		if (numAddedTriangles != 0) {
			// if there were triangles added to the model re-calculate normals
			long start = System.currentTimeMillis();
			LOGGER.debug("Re-calculating vertex normals ...");
			cas.getModel().updateVertexNormals();
			LOGGER.debug("Took: "
					+ PrintUtil.prettyMillis(System.currentTimeMillis() - start) + " (Vertices: "
					+ cas.getModel().getVertices().size() + ", Lines: " + cas.getModel().getLines().size()
					+ ", Triangles: " + cas.getModel().getTriangles().size() + ")");
		}
	}
	
	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.analyser.MeshAnalyser#updateProgress()
	 */
	@Override
	public void updateProgress() {
		if (allTriangles != null)
			setProgress((float) trianglesProcessed.get() / (float) allTriangles.size() * 100.0f);
	}
	
	/**
	 * Detects the sharp edges and sharp vertices of a triangle by comparing
	 * the angle of its normal vector with the normal vector of all its neighbors.
	 * The edge shared between two triangles with such an angle over the threshold
	 * is marked down as sharp edge for both of them and the edge's vertices are
	 * marked down as sharp vertices as well. 
	 * 
	 * In addition an edge can be sharp if the triangle processed does not have
	 * any neighbors associated with that edge.
	 * 
	 * @param t
	 * 			triangle to be analyzed
	 */
	private void sharpEdgeDetectionForTriangle(Triangle t) {
		synchronized (t) {
			Edge[] edges = t.getEdges();
			for (int i = 0 ; i < edges.length ; ++i) {
				List<Triangle> neighbors = t.getNeighborsOfEdge(edges[i]);
				// if no neighbors at the edge, then triangle is at the boundary 
				// of an open part of the mesh, so mark boundary as sharp edge
				if (neighbors.size() == 0) {
					edges[i].getVerticesOfEdge()[0].isSharpVertex(true);
					edges[i].getVerticesOfEdge()[1].isSharpVertex(true);
					t.addSharpEdge(edges[i]);
				}
				else if (neighbors.size() == 1 && !neighbors.get(0).isVisited()) {
					float angleOfNormals = (float) (Math.toDegrees(t.getNormalVector().angle(neighbors.get(0).getNormalVector())));
					if (angleOfNormals >= ThresholdsReasoning.SHARP_EDGE_ANGLE_TOL) {
						// get shared edge from the neighboring triangle
						Edge neighborCommonEdge = neighbors.get(0).getCommonEdge(t);
						edges[i].getVerticesOfEdge()[0].isSharpVertex(true);
						edges[i].getVerticesOfEdge()[1].isSharpVertex(true);
						neighborCommonEdge.getVerticesOfEdge()[0].isSharpVertex(true);
						neighborCommonEdge.getVerticesOfEdge()[1].isSharpVertex(true);
						t.addSharpEdge(edges[i]);
						neighbors.get(0).addSharpEdge(neighborCommonEdge);
					}
				}
				else {
					// if multiple neighbors then decided on individual basis if 
					// edges of neighbors are edge and weight on this the decision for
					// the triangle edge
					
					// perimeter of edge of sharp and non sharp triangle neighbors
					float perNonSharpNeighbors = 0.0f;
					for (int j = 0 ; j < neighbors.size() ; ++j) {
						float angleOfNormals = (float) (Math.toDegrees(t.getNormalVector().angle(neighbors.get(j).getNormalVector())));
						Edge neighborCommonEdge = neighbors.get(j).getCommonEdge(t);
						if (angleOfNormals >= ThresholdsReasoning.SHARP_EDGE_ANGLE_TOL && !neighbors.get(j).isVisited()) {
							neighborCommonEdge.getVerticesOfEdge()[0].isSharpVertex(true);
							neighborCommonEdge.getVerticesOfEdge()[1].isSharpVertex(true);
							neighbors.get(j).addSharpEdge(neighborCommonEdge);
							// also mark vertices on the triangle t to be sharp if same coordinates with the the ones from the neighbor
							if (edges[i].getVerticesOfEdge()[0].sameCoordinates(neighborCommonEdge.getVerticesOfEdge()[0])
									|| edges[i].getVerticesOfEdge()[0].sameCoordinates(neighborCommonEdge.getVerticesOfEdge()[1])) {
								edges[i].getVerticesOfEdge()[0].isSharpVertex(true);
							}
							else if (edges[i].getVerticesOfEdge()[1].sameCoordinates(neighborCommonEdge.getVerticesOfEdge()[0])
									|| edges[i].getVerticesOfEdge()[1].sameCoordinates(neighborCommonEdge.getVerticesOfEdge()[1])) {
								edges[i].getVerticesOfEdge()[1].isSharpVertex(true);
							}
						}
						else {
							perNonSharpNeighbors += neighborCommonEdge.getEdgeValue().length();
						}
					}
					// if perimeter of non sharp edges is 0, then also edge is sharp
					if (perNonSharpNeighbors == 0.0f) {
						t.addSharpEdge(edges[i]);
					}
				}
			}
			t.setIsVisited(true);
		}
	}
	
	/**
	 * Adds 3 triangles to the model if the passed triangle is a 
	 * sharp triangle as defined in "A new CAD mesh segmentation method, based on
	 * curvature tensor analysis", Guillaume Lavoue, Florent Dupont, Atilla Baskurt, 
	 * Computer-Aided Design 37 (2005) 975–987, and previously checked based
	 * on the processing done inside of {@link edu.tum.cs.vis.model.uima.analyser.EdgeAnalyser#sharpEdgeDetectionForTriangle(Triangle)}
	 * 
	 * This method builds up 3 triangles starting from the given one by adding the
	 * centroid of it as a new vertex. The new triangles neighboring relations are 
	 * handled based on the edge neighbors of the parent triangle, while the parent
	 * is removed from the neighbors list of its neighbors and also removed from the
	 * group of triangles of the CAD model 
	 * 
	 * @param t 
	 * 			triangle decomposed in 3 smaller triangles and then removed
	 */
	private void addTrianglesToModel(Triangle t) {
		Vertex newVertex = new Vertex(t.getCentroid().x, t.getCentroid().y, t.getCentroid().z);
		allVertices.add(newVertex);
		Triangle[] newTriangle = new Triangle[3];
		for (int i = 0 ; i < 3 ; ++i) {
			newTriangle[i] = new Triangle(t.getPosition()[i],t.getPosition()[(i+1)%3],newVertex);
			if (t.getAppearance() != null) {
				Appearance appearance = new Appearance(t.getAppearance());
				newTriangle[i].setAppearance(appearance);
			}
			// update now its edges and set normal vector and appearance
			newTriangle[i].updateEdges();
			newTriangle[i].setTexPosition(t.getTexPosition());
			newTriangle[i].setNormalVector(t.getNormalVector());
		}

		// add neighbors inside the big triangle
		newTriangle[0].addNeighbor(newTriangle[1]);
		newTriangle[0].addNeighbor(newTriangle[2]);
		newTriangle[1].addNeighbor(newTriangle[2]);
		
		// add triangle neighbors outside the original triangle 
		Edge[] edges = t.getEdges();
		for (int i = 0 ; i < edges.length ; ++i) {
			List<Triangle> toAdd = t.getNeighborsOfEdge(edges[i]);
			if (toAdd.size() == 0) {
				continue;
			}
			for (Triangle n : toAdd) {
				for (int j = 0 ; j < 3 ; ++j) {
					if (newTriangle[j].isDirectNeighbor(n, ThresholdsReasoning.FAST_NEIGHBOR_DETECTION)) {
						n.removeNeighbor(t);
						n.addNeighbor(newTriangle[j]);
						newTriangle[j].addNeighbor(n);
						break;
					}
				}
			}
		}
		
		// add sharp edges if any to the new 3 created triangles
		for (Edge sharpEdge : t.getSharpEdges()) {
			for (int i = 0 ; i < 3 ; ++i) {
				newTriangle[i].addSharpEdge(sharpEdge);
			}
		}
		
		// add new vertex as direct neighbor and old vertices as neighbors for new one
		// compute centroids of new triangles and add them to the model
		for (int i = 0 ; i < 3 ; ++i) {
			t.getPosition()[i].addNeighbor(newVertex);
			newVertex.addNeighbor(t.getPosition()[i]);
			newTriangle[i].updateCentroid();
			cas.getModel().getGroup().addTriangle(newTriangle[i]);
		}
		
		// now remove parent triangle from the group of the model
		cas.getModel().getGroup().removeTriangle(t);
	}
}