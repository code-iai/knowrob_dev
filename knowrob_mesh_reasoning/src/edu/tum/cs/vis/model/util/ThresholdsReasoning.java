/*******************************************************************************
 * Copyright (c) 2014 Andrei Stoica. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Andrei Stoica - initial API and implementation during Google Summer of Code 2014
 ******************************************************************************/

package edu.tum.cs.vis.model.util;

/**
 * Class that implements the main thresholds and constants used in the mesh reasoning
 * process. This is meant to be a utility class allowing fast access to some of the
 * knobs that allow some performance increase of the segmentation from model to model
 * 
 * The main switches to car for are:
 * 			-{@code NUM_CLUSTERS}
 * 			-{@code MIN_DISTANCE_TOL} 
 * 			-{@code SHARP_EDGE_ANGLE_TOL}
 * 
 * @author Andrei Stoica
 */
public class ThresholdsReasoning extends Thresholds {

	/**
	* Parameter used as the minimum value of the distance measure calculations
	* for the region adjacency graph of the CAD model.
	* 
	* {@link edu.tum.cs.vis.model.util.algorithm.CasProcessing#process()}
	*/
	public static final float EPSILON 						= 1e-5f;
	
	/**
	* Default number of clusters used to classify the curvature properties of
	* the vertices of the model. 
	* 
	* {@link edu.tum.cs.vis.model.util.algorithm.CasProcessing#process()}
	* 
	* This value will not be changed by {@code CodeProcessing} as long as it
	* does not exceed 2.5% of the number of model vertices. If this is the case
	* it is replaced by this value.
	*/
	public static final int NUM_CLUSTERS 					= 15;
	
	/**
	* Limit of iterations that should be performed on the KMeans clustering algorithm
	* after the initial set up of the clusters.
	* 
	* {@link edu.tum.cs.vis.model.util.algorithm.CasProcessing#process()}
	*/
	public static final int ITERATIONS_LIM 					= 200;
	
	/**
	* Threshold for the minimum area of a region. This parameter is used for faster
	* integration of small regions into bigger regions during the region growing
	* algorithm of the MeshCas processing.
	* 
	* The regions that have an area smaller than this threshold have the area distance 
	* set to {@code ThresholdsReasoning.EPSILON} during the adjacency score computation
	* 
	* {@link edu.tum.cs.vis.model.util.algorithm.CasProcessing#process()}
	*/
	public static final float MIN_AREA						= 0.05f;
	
	/**
	* Minimum distance threshold for the region adjacency distances. This parameter limits
	* the region merging stopping the merging process when the regions are too disimilar 
	* (the smallest distance between to regions is bigger than this threshold)
	* 
	* {@link edu.tum.cs.vis.model.util.algorithm.CasProcessing#process()}
	*/
	public static final float MIN_DISTANCE_TOL 				= 5f;
	
	/**
	* Weighting parameter used in computing the curvature hue and saturation properties.
	* 
	* {@link edu.tum.cs.vis.model.util.algorithm.CurvatureCalculation#setCurvatureHueSaturation(java.util.HashMap, edu.tum.cs.vis.model.Model, float)}
	* 
	* NOTE: since the class {@link edu.tum.cs.vis.model.util.algorithm.CasProcessing} performs an averaging filtering of the curvature properties, 
	* this parameter should be kept to low values.
	*/
	public static final float CURV_SMOOTHING 				= 0.01f;
	
	/**
	* Minimum threshold angle used in deciding whether or not an edge shared by two neighboring
	* triangles is sharp or not. The angle between the normals of two triangle neighbors is compared
	* to this value and if it is above it, then the shared edge is considered to be sharp for both
	* triangle instances.
	* 
	* {@link edu.tum.cs.vis.model.util.algorithm.CasProcessing#process()}
	*/
	public static final float SHARP_EDGE_ANGLE_TOL 			= 65.0f;
	
	/**
	 * PrimitiveAnalyser maximal threshold for the angle between
	 * the normals of two neighboring plane annotations in order
	 * to combine them in one larger plane annotation.
	 * 
	 *  {@link edu.tum.cs.vis.model.uima.analyser.PrimitiveAnalyser}
	 */
	public static final float PA_PLANE_COMBINE_DEGREE		= 5.0f;
	
	/**
	 * PrimitiveAnalyser maximal threshold for the angle between 
	 * the normals of two neighboring triangles to add them to the 
	 * same plane annotation.
	 * 
	 * {@link edu.tum.cs.vis.model.uima.analyser.PrimitiveAnalyser}
	 */
	public static final float PA_PLANE_TOLERANCE			= 2.0f;
	
	/**
	* Default constructor to avoid creation of new instances
	*/
	private ThresholdsReasoning() {
		return;
	}
}
