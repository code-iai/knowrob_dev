/*******************************************************************************
 * Copyright (c) 2014 Andrei Stoica. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Andrei Stoica - initial API and implementation, Year: 2014
 ******************************************************************************/
package edu.tum.cs.vis.model.uima.analyser;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.Callable;

import javax.vecmath.Vector3f;

import org.apache.log4j.Logger;

import edu.tum.cs.vis.model.util.Edge;
import edu.tum.cs.vis.model.util.Triangle;
import edu.tum.cs.vis.model.util.Vertex;
import edu.tum.cs.vis.model.uima.cas.MeshCas;
import edu.tum.cs.ias.knowrob.utils.ThreadPool;

public class EdgeAnalyser extends MeshAnalyser {
	
	/**
	 * Log4j logger
	 */
	private static Logger		logger = Logger.getLogger(EdgeAnalyser.class);
	
	/**
	 * Number of currently processed triangles. Used for progress status
	 */
	final AtomicInteger			trianglesProcessed = new AtomicInteger(0);
	
	/**
	 * List of all vertices of the CAD model
	 */
	List<Vertex> 				allVertices;
	
	/**
	 * List of all triangles of the CAD model
	 */
	List<Triangle>				allTriangles;
	
	/**
	 * Number of added triangles to the model
	 */
	private int 				numAddedTriangles;
	
	private MeshCas				cas;
	
	/**
	 * Getter for the number of the added triangles
	 */
	public int getNumAddedTriangles() {
		return numAddedTriangles;
	}
	
	@Override
	public Logger getLogger() {
		return logger;
	}
	
	@Override
	public String getName() {
		return "Edge";
	}
	
	@Override
	public void processStart(MeshCas newCas) {
		cas = newCas;
		allTriangles = cas.getModel().getTriangles();
		allVertices = cas.getModel().getVertices();
		trianglesProcessed.set(0);
		
		// remove any wrong tesselated triangles in the model
//		this.removeCollinearTriangles();
		
		this.numAddedTriangles = allTriangles.size();
		
		// perform the sharp edge detection for individual triangles (multi-threaded)
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
					trianglesProcessed.addAndGet(end - st);
					return null;
				}

			});
		}
		ThreadPool.executeInPool(threads);
		
		trianglesProcessed.set(0);
		List<Triangle> toRemove = new ArrayList<Triangle>();
		for (int i = 0 ; i < allTriangles.size() ; ++i) {
			Triangle t = allTriangles.get(i);
			trianglesProcessed.incrementAndGet();
			if (t.checkIsSharpTriangle()) {
				addTrianglesToModel(t);
				toRemove.add(t);
				trianglesProcessed.decrementAndGet();
			}
		}
		cas.getModel().getGroup().removeTriangle(toRemove);
		cas.getModel().reloadVertexList();
//		allTriangles.remove(toRemove);
		
//		this.removeCollinearTriangles();
	
		this.numAddedTriangles = allTriangles.size() - this.numAddedTriangles;
		if (this.numAddedTriangles != 0) {
			logger.debug("Added " + this.numAddedTriangles + " triangles to the model");
		}
		else {
			logger.debug("No triangles added to the model");
		}
	}
	
	@Override
	public void updateProgress() {
		if (allTriangles != null)
			setProgress((float) trianglesProcessed.get() / (float) allTriangles.size() * 100.0f);
	}
	
//	/**
//	 * Removes "colinear triangles", i.e. "triangles" which have 3 "colinear" vertices
//	 */
//	private void removeCollinearTriangles() {
//		List<Triangle> allTrianglesToProcess = new ArrayList<Triangle>();
//		allTrianglesToProcess.addAll(allTriangles);
//		int rmTriangles = 0;
//		for (int i = 0 ; i < allTrianglesToProcess.size() ; ++i) {
//			Triangle tr = allTrianglesToProcess.get(i);
////			Edge[] edges = tr.getEdges();
//			if (tr.getNeighbors().size() == 0) {
//				logger.debug("Removing " + tr);
//				for (int j = 0 ; j < allTrianglesToProcess.size() ; ++j) {
//					Triangle tn = allTrianglesToProcess.get(j);
//					if (tn.getNeighbors().contains(tr)) {
//						tn.removeNeighbor(tr);
//					}
//				}
//				cas.getModel().getGroup().removeTriangle(tr);
//				allTriangles.remove(tr);
//				allTrianglesToProcess.remove(tr);
//				rmTriangles++;
//			}
//		}
//		if (rmTriangles > 0) {
//			logger.debug("Removed " + rmTriangles + " triangles");
//		}
//	}
	
	/**
	 * Performs the sharp detection at the triangle level
	 */
	private void sharpEdgeDetectionForTriangle(Triangle t) {
		synchronized (t) {
		Iterator<Triangle> it = t.getNeighbors().iterator();
		while (it.hasNext()) {
			Triangle n = it.next();
			//synchronized (n) {
			float angleOfNormals = (float)Math.toDegrees(t.getNormalVector().angle(n.getNormalVector()));
			if (angleOfNormals >= 85.0f) {
				List<Vertex> vShared = findSharedVertices(t,n);
				Edge edge = new Edge(vShared.get(0), vShared.get(2));
				vShared.get(0).isSharpVertex(true);
				vShared.get(1).isSharpVertex(true);
				vShared.get(2).isSharpVertex(true);
				vShared.get(3).isSharpVertex(true);
				t.addSharpEdge(edge);
				n.addSharpEdge(edge);
			}
			//}
		}
//		Edge[] edges = t.getEdges();
//		for (int i = 0 ; i < edges.length ; ++i) {
//			List<Triangle> neighbors = t.getNeighborsOfEdge(edges[i]);
//			if (neighbors.size() == 0) {
//				edges[i].getVerticesOfEdge()[0].isSharpVertex(true);
//				edges[i].getVerticesOfEdge()[1].isSharpVertex(true);
//				t.addSharpEdge(edges[i]);
//			}
//		}
		}
	}
	
	/**
	 * Finds the two common points of two neighboring triangles
	 */
	private List<Vertex> findSharedVertices(Triangle t, Triangle n) {
		List<Vertex> v = new ArrayList<Vertex>(4);
		Edge[] tEdges = t.getEdges();
		Edge[] nEdges = n.getEdges();
		for (int i = 0 ; i < tEdges.length ; ++i) {
			for (int j = 0 ; j < nEdges.length ; ++j) {
				if (tEdges[i].isDirectNeighbor(nEdges[j])) {
					v.add(tEdges[i].getVerticesOfEdge()[0]);
					v.add(nEdges[j].getVerticesOfEdge()[0]);
					v.add(tEdges[i].getVerticesOfEdge()[1]);
					v.add(nEdges[j].getVerticesOfEdge()[1]);
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
		allVertices.add(newVertex);
		Triangle[] newTriangle = new Triangle[3];
		for (int i = 0 ; i < 3 ; ++i) {
			newTriangle[i] = new Triangle(t.getPosition()[i],t.getPosition()[(i+1)%3],newVertex);
			// update now its edges and set normal vector and appearance
			newTriangle[i].updateEdges();
			newTriangle[i].setAppearance(t.getAppearance());
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
					if (newTriangle[j].isDirectNeighbor(n)) {
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
	}
}