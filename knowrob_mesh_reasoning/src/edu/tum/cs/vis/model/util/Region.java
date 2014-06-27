/*******************************************************************************
 * Copyright (c) 2014 Andrei Stoica. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Andrei Stoica - initial API and implementation, Year: 2014
 ******************************************************************************/

package edu.tum.cs.vis.model.util;

import java.util.ArrayList;
import java.util.List;

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
	private List<Triangle>			triangles = new ArrayList<Triangle>();
	
	/**
	 * A list of all the triangles on the boundary 
	 * of the region
	 */
	private List<Triangle>			boundaryTriangles = new ArrayList<Triangle>();
	
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
		this.addTriangleToRegion(initSeedTr);
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
	public List<Triangle> getTriangles() {
		return this.triangles;
	}
	
	/**
	 * Getter for a region's triangle boundary
	 */
	public List<Triangle> getBoundaryTriangles() {
		return this.boundaryTriangles;
	}
	
	/**
	 * Getter for the curvature values associated with the region
	 */
	public float[] getCurvatureMinMaxOfRegion() {
		return this.curvatureMinMax;
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
	 * Subtracts the area of a specified triangle
	 * from the area of the region
	 */
	public void updateAreaOfRegion(final Triangle tr, final boolean toRemove) {
		area -= tr.getArea();
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
		if (!triangles.contains(newTr)) {
			triangles.add(newTr);
			newTr.setRegionLabel(id);
			this.updateAreaOfRegion(newTr);
			this.addTriangleToRegionBoundary(newTr);
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
			List<Triangle> toRemove = new ArrayList<Triangle>();
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
	 * Removes triangle from region and updates the list 
	 * of triangles and the boundary ones
	 */
	public void removeTriangleFromRegion(final Triangle toRemove) {
		if (triangles.remove(toRemove)) {
			updateAreaOfRegion(toRemove,true);
			if (boundaryTriangles.remove(toRemove)) {
				System.out.println("Triangle " + toRemove.toString() + " also removed from the boundary region");
				this.updateRegionBoundary();
			}
		}
		// unset region label
		toRemove.setRegionLabel(-1);
	}
	
	/**
	 * Updates the boundary region based on the existing
	 * triangles inside the region itself
	 */
	public void updateRegionBoundary() {
		for (Triangle tr : triangles) {
			if (!triangles.containsAll(tr.getNeighbors())) {
				boundaryTriangles.add(tr);
			} else {
				boundaryTriangles.remove(tr);
			}
		}
	}
	
	/**
	 * Builds the region starting from the initial 
	 * triangle seed used to initialize the region
	 */
	public void buildUpRegion() {
		for (int i = 0 ; i < triangles.size() ; ++i) {
			Triangle tr = triangles.get(i);
			Edge[] nonSharpEdges = tr.getNonSharpEdges();
			for (int j = 0 ; j < nonSharpEdges.length ; ++j) {
				Triangle neighbor = tr.getNeighborOfEdge(nonSharpEdges[j]);
				if (neighbor == null) {
//					System.out.println(tr.getNeighbors().size());
//					System.out.println("null @: " + i + " " + j);
					continue;
				}
				Vertex oppositeVertex = neighbor.getOppositeVertexFromEdge(nonSharpEdges[j]);
				if ((oppositeVertex.isSharpVertex()) || ((oppositeVertex.getClusterCurvatureVal()[0] == curvatureMinMax[0]) 
						&& (oppositeVertex.getClusterCurvatureVal()[1] == curvatureMinMax[1]))) {
					this.addTriangleToRegion(neighbor);
				}
			}
		}
	}
	
	@Override
	public String toString() {
		String print = "Region ID " + id + "\n";
		print = print + "Curvatures: KMin = " + curvatureMinMax[0] + ", KMax = " + curvatureMinMax[1] + "\n";
		print = print + "Triangles: " + triangles.size() + "\n";
		print = print + "Boundary Triangles: " + boundaryTriangles.size() + "\n";
		return print;
	}
}