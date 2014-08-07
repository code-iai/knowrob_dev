/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Stefan Profanter - initial API and implementation, Year: 2012
 ******************************************************************************/
package edu.tum.cs.vis.model.uima.analyser;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import org.apache.log4j.Logger;

import edu.tum.cs.ias.knowrob.utils.ThreadPool;
import edu.tum.cs.vis.model.uima.cas.MeshCas;
import edu.tum.cs.vis.model.util.Edge;
import edu.tum.cs.vis.model.util.Triangle;
import edu.tum.cs.vis.model.util.Vertex;

/**
 * Analyzer for a mesh which sets direct neighbors of a triangle.
 * 
 * The neighbor information is used in other analyzers for better performance.
 * 
 * @author Stefan Profanter
 * 
 */
public class NeighborAnalyser extends MeshAnalyser {

	/**
	 * Log4j logger
	 */
	private static Logger	logger				= Logger.getLogger(NeighborAnalyser.class);

	/**
	 * Number of currently processed triangles. Used for progress status.
	 */
	final AtomicInteger		trianglesElaborated	= new AtomicInteger(0);

	/**
	 * When calling <code>process</code> all triangles of the group and its children are collected
	 * in this list to process them afterwards.
	 */
	List<Triangle>			allTriangles;
	
	/**
	 * Stores initial number of triangles
	 */
	int 					initialNumOfTriangles;

	@Override
	public Logger getLogger() {
		return logger;
	}

	@Override
	public String getName() {
		return "Neighbor";
	}

	@Override
	public void processStart(MeshCas cas) {

		trianglesElaborated.set(0);
		allTriangles = cas.getModel().getTriangles();

		removeCollinearTriangles(allTriangles);
		
		initialNumOfTriangles = allTriangles.size();
		
		final int interval = 100;

		final Lock lock = new ReentrantLock();

		List<Callable<Void>> threads = new LinkedList<Callable<Void>>();

		for (int start = 0; start < allTriangles.size(); start += interval) {
			final int st = start;
			threads.add(new Callable<Void>() {

				@Override
				public Void call() throws Exception {
					int end = Math.min(st + interval, allTriangles.size());
					for (int i = st; i < end; i++) {
						Triangle tr = allTriangles.get(i);
						for (Vertex v : tr.getPosition()) {
							// vertices on same triangle are always direct neighbors
							v.addNeighbor(tr.getPosition());
						}
						for (int j = i + 1; j < allTriangles.size(); j++) {
							Triangle n = allTriangles.get(j);
							// check and add triangle as neighbor
							if (n.addNeighbor(tr, lock)) {
								// they are neighbors, update vertex neighbors for vertices which
								// are on both triangles
								for (Vertex v : tr.getPosition()) {
									for (Vertex vn : n.getPosition()) {
										if (v != vn && v.sameCoordinates(vn)) {
											// not same vertex reference, but same position ->
											// combine neighbor vertices
											v.addNeighbor(vn.getNeighbors());
											vn.addNeighbor(v.getNeighbors());
										}
									}
								}
							}
						}
					}
					trianglesElaborated.addAndGet(end - st);
					return null;
				}

			});
		}
		ThreadPool.executeInPool(threads);
		updateProgress();
	}

	@Override
	public void updateProgress() {
		if (allTriangles != null)
			setProgress((float) trianglesElaborated.get() / (float) initialNumOfTriangles * 100.0f);
	}
	
	/**
	 * Checks all the triangles associated to the UIMA CAS of the CAD model in order to see
	 * if they collapse to trivial triangle annotations ("collinear" triangles, where two vertices
	 * are the same, or where two edges are parallel)
	 * 
	 * @param allTriangles
	 * 				list of all triangles in the model to be checked
	 */
	private void removeCollinearTriangles(final List<Triangle> allTriangles) {
		List<Triangle> toRemove = new ArrayList<Triangle>();
		for (Triangle triangle : allTriangles) {
			Vertex v0 = triangle.getPosition()[0];
			Vertex v1 = triangle.getPosition()[1];
			Vertex v2 = triangle.getPosition()[2];
			if ((v0.x == v1.x && v0.y == v1.y && v0.z == v1.z) || (v0.x == v2.x && v0.y ==v2.y && v0.z == v2.z) || (v1.x == v2.x && v1.y == v2.y && v1.z == v2.z)) {
				toRemove.add(triangle);
				continue;
			}
			Edge[] edge = triangle.getEdges();
			float angle = (float)(Math.toDegrees(edge[0].getEdgeValue().angle(edge[1].getEdgeValue())));
			if (angle == 0.0f || angle == 180.0f) {
				toRemove.add(triangle);
				continue;
			}
			angle = (float)(Math.toDegrees(edge[0].getEdgeValue().angle(edge[2].getEdgeValue())));
			if (angle == 0.0f || angle == 180.0f) {
				toRemove.add(triangle);
				continue;
			}
			angle = (float)(Math.toDegrees(edge[1].getEdgeValue().angle(edge[2].getEdgeValue())));
			if (angle == 0.0f || angle == 180.0f) {
				toRemove.add(triangle);
				continue;
			}	
		}
		logger.debug("Removed " + toRemove.size() + " collinear triangles from the model.");
		allTriangles.removeAll(toRemove);
	}
}
