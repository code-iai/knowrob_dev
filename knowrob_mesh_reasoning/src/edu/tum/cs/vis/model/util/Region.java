/*******************************************************************************
 * Copyright (c) 2014 Andrei Stoica. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Andrei Stoica - initial API and implementation, Year: 2014
 ******************************************************************************/

package edu.tum.cs.vis.model.util;

import java.util.HashSet;
import java.util.Set;

import edu.tum.cs.vis.model.util.Triangle;

/**
 * Class that implements the data storage and functionality
 * of a region as a connected area of triangles on the mesh
 * surface that share common curvature properties. The triangles
 * are situated on relatively contained and localized position
 * on the mesh.
 */
public class Region {
	
	/**
	 * Stores the regions label id
	 */
	private final int 				id;
	
	/**
	 * The area of the entire region
	 */
	private float					area;
	
	/**
	 * The curvature parameters of the boundary.
	 * The first parameter stores the minimum value
	 * of the curvature and the second one the maximum
	 * value of the parameter
	 */
	private float[]					curvatureMinMax = new float[2];
	
	/**
	 * A list of all the triangles inside the region
	 */
	private Set<Triangle>			triangles = new HashSet<Triangle>();
	
	/**
	 * A list of all the triangles on the boundary 
	 * of the region
	 */
	private Set<Triangle>			boundaryTriangles = new HashSet<Triangle>();
	
	/**
	 * Constructor for a region that takes only the region id
	 */
	public Region(final int newId) {
		this.id = newId;
		this.area = 0;
		this.curvatureMinMax[0] = this.curvatureMinMax[1] = 0.0f;
	}
	
	/**
	 * Constructor for a region that takes the region
	 * id and the curvature values for a region.
	 */
	public Region(final int newId, final float[] newCurvatureValues) {
		this.id = newId;
		this.area = 0;
		this.curvatureMinMax[0] = newCurvatureValues[0];
		this.curvatureMinMax[1] = newCurvatureValues[1];
	}
	
	/**
	 * Constructor for a region that takes the region id
	 * and the initial seed triangle from which the region
	 * build-up starts
	 */
	public Region(final int newId, final Triangle initSeedTr) {
		this.id = newId;
		this.area = initSeedTr.getArea();
		this.curvatureMinMax[0] = initSeedTr.getCurvatureValues()[0];
		this.curvatureMinMax[1] = initSeedTr.getCurvatureValues()[1];
	}
	
	/**
	 * Getter for the region's id
	 */
	public int getRegionId() {
		return this.id;
	}
	
	/**
	 * Getter for area property
	 */
	public float getAreaOfRegion() {
		return this.area;
	}
	
	/**
	 * Getter for a region's triangles
	 */
	public Set<Triangle> getTriangles() {
		return this.triangles;
	}
	
	/**
	 * Getter for a region's triangle boundary
	 */
	public Set<Triangle> getBoundaryTriangles() {
		return this.boundaryTriangles;
	}
	
	/**
	 * Calculates the region surface area based on
	 * the individual areas of the triangles included
	 * in the region
	 */
	public void updateAreaOfRegion() {
		area = 0;
		for (Triangle tr : triangles) {
			area += tr.getArea();
		}
	}
	
	/**
	 * Adds up the area of a specified triangle
	 * to the area of the region
	 */
	public void updateAreaOfRegion(final Triangle tr) {
		area += tr.getArea();
	}
	
	/**
	 * Adds triangle to the list of triangles 
	 * included in the region. The addition is
	 * performed only if the triangle was not
	 * previously inside the region. If the triangle
	 * was already in the region, this function does
	 * nothing.
	 */
	public void addTriangleToRegion(final Triangle newTr) {
		if (triangles.add(newTr)) {
			updateAreaOfRegion(newTr);
			addTriangleToRegionBoundary(newTr);
		}
	}
	
	/**
	 * Checks if a triangle is at a region's boundary
	 * and if it is it adds it to the list and then
	 * removes any triangles from the list that are not
	 * longer at the boundary. If the triangle is not at
	 * the boundary, this function does nothing. The 
	 * check is based on the neighbors of the triangle
	 * compared against the triangles already in the 
	 * region.
	 */
	public void addTriangleToRegionBoundary(final Triangle newTr) {
		if (!triangles.containsAll(newTr.getNeighbors())) {
			Set<Triangle> toRemove = new HashSet<Triangle>();
			// check for old triangles that are not on boundary
			// of current updated region
			for (Triangle tr : boundaryTriangles) {
				if (triangles.containsAll(tr.getNeighbors())) {
					toRemove.add(tr);
				}
			}
			// remove triangles that are not on the boundary anymore
			boundaryTriangles.removeAll(toRemove);
			// add current Triangle object located on boundary
			boundaryTriangles.add(newTr);
		}
	}
	
	/**
	 * Removes triangle from region and updates
	 */
}