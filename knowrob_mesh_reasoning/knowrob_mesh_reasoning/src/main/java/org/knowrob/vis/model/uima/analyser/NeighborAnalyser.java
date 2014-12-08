/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: 
 * 				Stefan Profanter - initial API and implementation, Year: 2012
 * 				Andrei Stoica - refactored implementation during Google Summer of Code 2014
 ******************************************************************************/
package org.knowrob.vis.model.uima.analyser;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import javax.vecmath.Vector3f;

import org.apache.log4j.Logger;

import org.knowrob.utils.ThreadPool;
import org.knowrob.vis.model.uima.cas.MeshCas;
import org.knowrob.vis.model.util.Triangle;
import org.knowrob.vis.model.util.Vertex;
import org.knowrob.vis.util.PrintUtil;
import org.knowrob.vis.model.util.Edge;
import org.knowrob.utils.ThreadPool;


/**
 * Analyzer for a mesh which sets the direct neighbors of the triangles and
 * vertices in the CAD model parsed. The neighboring information is used
 * throughout the entire implementation to describe the CAD model
 * 
 * @author Stefan Profanter
 * @author Andrei Stoica (added collinear removal and refactored neighboring
 * to also make up for badly tesselated meshes with artifacts, e.g. gaps, not
 * connected vertices, floating edges etc.)
 */
public class NeighborAnalyser extends MeshAnalyser {

	/**
	 * Log4j logger
	 */
	private static final Logger	LOGGER				= Logger.getLogger(NeighborAnalyser.class);

	/**
	 * Number of currently processed triangles. Used for progress status.
	 */
	private final AtomicInteger			trianglesElaborated	= new AtomicInteger(0);

	/**
	 * List of all the triangles of the associated MeshCas of the CAD model
	 */
	private List<Triangle>				allTriangles;
	
	/**
	 * Initial number of triangles
	 */
	private int 						initialNumOfTriangles;

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
		return "Neighbor";
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.analyser.MeshAnalyser#processStart(edu.tum.cs.vis.model.uima.cas.MeshCas)
	 */
	@Override
	public void processStart(MeshCas cas) {
		if (cas == null || cas.getModel() == null || cas.getModel().getTriangles() == null) {
			return;
		}
		// in ply (and also collada) files there may be double sided triangles,
		// so remove them if they exist in the current model
		long startTime = System.currentTimeMillis();
		LOGGER.debug("Checking for double sided triangles ...");
		LOGGER.debug("Removed " + cas.getModel().removeDoubleSidedTriangles() + " triangles. Took: "
				+ PrintUtil.prettyMillis(System.currentTimeMillis() - startTime));
		
		trianglesElaborated.set(0);
		allTriangles = cas.getModel().getTriangles();
		
		// now remove collinear triangles
		removeCollinearTriangles(cas);
		
		// update number of valid triangles in the model
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
		
		// update the vertex sharing and neighboring information
		startTime = System.currentTimeMillis();
		LOGGER.debug("Checking for vertex sharing and calculating vertex normals ...");
		cas.getModel().updateVertexSharing();

		cas.getModel().updateVertexNormals();
		LOGGER.debug("Model initialized. Took: "
				+ PrintUtil.prettyMillis(System.currentTimeMillis() - startTime) + " (Vertices: "
				+ cas.getModel().getVertices().size() + ", Lines: " + cas.getModel().getLines().size()
				+ ", Triangles: " + cas.getModel().getTriangles().size() + ")");
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.analyser.MeshAnalyser#updateProgress()
	 */
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
	private void removeCollinearTriangles(MeshCas cas) {
		List<Triangle> toRemove = new ArrayList<Triangle>();
		for (Triangle triangle : allTriangles) {
			Vertex v0 = triangle.getPosition()[0];
			Vertex v1 = triangle.getPosition()[1];
			Vertex v2 = triangle.getPosition()[2];
			// if two vertices share the same coordinates than remove triangle as collinear
			if ((v0.x == v1.x && v0.y == v1.y && v0.z == v1.z) || 
					(v0.x == v2.x && v0.y ==v2.y && v0.z == v2.z) || 
					(v1.x == v2.x && v1.y == v2.y && v1.z == v2.z)) {
				toRemove.add(triangle);
				continue;
			}
			// if any two edges are parallel (cross product has absolute value 0) then remove triangle as collinear
			Edge[] edge = triangle.getEdges();
			Vector3f crossProd = new Vector3f();
			crossProd.cross(edge[0].getEdgeValue(), edge[1].getEdgeValue());
			if (crossProd.length() == 0.0f) {
				toRemove.add(triangle);
				continue;
			}
		}
		LOGGER.debug("Removed " + toRemove.size() + " collinear triangles from the model.");
		allTriangles.removeAll(toRemove);
		// remove triangles also from the group
		cas.getModel().getGroup().removeTriangle(toRemove);
	}
}
