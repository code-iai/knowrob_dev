/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: 
 * 		Stefan Profanter - initial API and implementation, Year: 2012
 * 		Andrei Stoica - refactored implementation during Google Summer of Code 2014
 ******************************************************************************/
package edu.tum.cs.vis.model.uima.analyser;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;

import javax.vecmath.Vector3f;

import org.apache.log4j.Logger;

import com.google.common.collect.HashMultimap;

import edu.tum.cs.ias.knowrob.utils.ThreadPool;
import edu.tum.cs.uima.Annotation;
import edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation;
import edu.tum.cs.vis.model.uima.annotation.primitive.ConeAnnotation;
import edu.tum.cs.vis.model.uima.annotation.primitive.PlaneAnnotation;
import edu.tum.cs.vis.model.uima.annotation.primitive.PrimitiveType;
import edu.tum.cs.vis.model.uima.annotation.primitive.SphereAnnotation;
import edu.tum.cs.vis.model.uima.cas.MeshCas;
import edu.tum.cs.vis.model.util.Curvature;
import edu.tum.cs.vis.model.util.Edge;
import edu.tum.cs.vis.model.util.Region;
import edu.tum.cs.vis.model.util.Triangle;
import edu.tum.cs.vis.model.util.ThresholdsReasoning;
import edu.tum.cs.vis.model.util.Vertex;

/**
 * Mesh analyzer which assigns each triangle to its corresponding primitive type 
 * according to its curvature estimated and filtered properties.
 * 
 * @author Stefan Profanter
 * @author Andrei Stoica (refactor of the processing block)
 */
public class PrimitiveAnalyser extends MeshAnalyser {

	/**
	 * Log4J Logger
	 */
	private static Logger		logger					= Logger.getLogger(PrimitiveAnalyser.class);

	/**
	 * Threshold degree between normal vertices for allowing combination of neighboring 
	 * plane annotations
	 */
	private static final float	PLANE_COMBINE_DEGREE	= ThresholdsReasoning.PA_PLANE_COMBINE_DEGREE;

	/**
	 * Tolerance in degrees between two surface normals of neighboring triangles in order
	 * to add them to the same plane annotation
	 */
	private static double		PLANE_TOLERANCE			= ThresholdsReasoning.PA_PLANE_TOLERANCE;
	
	/**
	 * Map which maps a triangle to a certain primitive type: PLANE, SPHERE.CONVEX, 
	 * SPHERE.CONCAVE, CONE.CONVEX, CONE.CONCAVE
	 */
	private final HashMap<Triangle, PrimitiveType>	trianglePrimitiveTypeMap	= new HashMap<Triangle, PrimitiveType>();

	/**
	 * List of all the vertices of the CAD model under analysis
	 */
	List<Vertex>									allVertices;

	/**
	 * List of all the triangles of the CAD model under analysis
	 */
	List<Triangle>									allTriangles;
	
	/**
	 * List of all the regions of the CAD model under analysis
	 */
	List<Region>									allRegions;

	/**
	 * The number of triangles already elaborated / processed by the analyser. 
	 * Used for indicating current progress in the GUI
	 */
	final AtomicInteger								itemsElaborated	= new AtomicInteger(0);

	/**
	 * The number of triangles remained to be processed by the analyser.
	 * Used for indicating current progress in the GUI
	 */
	private int										itemsToElaborate = 0;

	/**
	 * Analyze given triangle and try to grow annotations (PLANE, SPHERE, CONE) based
	 * on the processed vertices primitives types and on the same properties of its
	 * connected neighbors. PLANE, SPHERE.CONVEX, SPHERE.CONCAVE, CONE.CONVEX, CONE.CONCAVE
	 * annotations can be grown by this routine. The triangles added to the annotations
	 * are marked under the {@code alreadyInAnnotation} set of already evaluated triangles
	 * 
	 * @param cas
	 *            MeshCas
	 * @param triangle
	 *            triangle to analyze
	 * @param alreadyInAnnotation
	 *            set of already analyzed triangles
	 */
	protected void analyseTriangle(MeshCas cas, Triangle triangle, Set<Triangle> alreadyInAnnotation) {
		if (alreadyInAnnotation.contains(triangle))
			return;

		@SuppressWarnings("rawtypes")
		PrimitiveAnnotation annotation;
		PrimitiveType type = getTrianglePrimitiveType(cas.getCurvatures(), triangle);

		if (type == PrimitiveType.PLANE)
			annotation = new PlaneAnnotation(cas.getCurvatures(), cas.getModel());
		else if (type == PrimitiveType.SPHERE_CONCAVE || type == PrimitiveType.SPHERE_CONVEX)
			annotation = new SphereAnnotation(cas.getCurvatures(), cas.getModel(), type == PrimitiveType.SPHERE_CONCAVE);
		else
			annotation = new ConeAnnotation(cas.getCurvatures(), cas.getModel(), type == PrimitiveType.CONE_CONCAVE);

		synchronized (annotation.getMesh().getTriangles()) {
			annotation.getMesh().getTriangles().add(triangle);
		}
		alreadyInAnnotation.add(triangle);

		synchronized (cas.getAnnotations()) {
			cas.addAnnotation(annotation);
		}

		// List of already visited triangles for BFS
		HashSet<Triangle> visited = new HashSet<Triangle>();
		visited.add(triangle);

		// FIFO queue for triangles to visit for BFS
		LinkedList<Triangle> queue = new LinkedList<Triangle>();
		Edge[] edges = triangle.getEdges();
		for (int i = 0 ; i < edges.length ; ++i) {
			if (!edges[i].getIsSharpEdge()) {
				List<Triangle> neighbors = triangle.getNeighborsOfEdge(edges[i]);
				for (int j = 0 ; j < neighbors.size() ; ++j) {
					queue.add(neighbors.get(j));
				}
			}
		}

		while (!queue.isEmpty()) {
			Triangle currNeighbor = queue.pop();
			visited.add(currNeighbor);
			if (alreadyInAnnotation.contains(currNeighbor))
				continue;

			boolean isSameType = (type == getTrianglePrimitiveType(cas.getCurvatures(),currNeighbor));

			boolean isSameNormal = planeAngleWithinTolerance(triangle.getNormalVector(),currNeighbor.getNormalVector());

			if (isSameType && type == PrimitiveType.PLANE)
				isSameType = isSameNormal; // only combine triangles which lie on the same plane

			if (isSameType) {
				synchronized (annotation.getMesh().getTriangles()) {
					annotation.getMesh().getTriangles().add(currNeighbor);
				}
				alreadyInAnnotation.add(currNeighbor);

				// Add all neighbors of current triangle to queue
				Edge[] currEdges = currNeighbor.getEdges();
				for (int i = 0 ; i < currEdges.length ; ++i) {
					if (!currEdges[i].getIsSharpEdge()) {
						List<Triangle> currNeighbors = currNeighbor.getNeighborsOfEdge(currEdges[i]);
						for (int j = 0 ; j < currNeighbors.size() ; ++j) {
							synchronized (annotation.getMesh()) {
								synchronized (annotation.getMesh().getTriangles()) {
									if (visited.contains(currNeighbors.get(j)) || annotation.getMesh().getTriangles().contains(currNeighbors.get(j))) {
										continue;
									}
								}
							}
							queue.add(currNeighbors.get(j));
						}
					}
				}
			}
		}
	}
	
	/**
	 * Analyze given triangle and try to grow a plane annotation. Only PLANE annotations are generated
	 * in this method. Curvature calculation has some problems with simple planes where a circle is
	 * cut out. Thus we try to fix this problem by first searching planes regardless the curvature
	 * value by only comparing their surface normal.
	 * 
	 * @param cas
	 *            MeshCas
	 * @param triangle
	 *            triangle from where to start region growing
	 * @param alreadyInAnnotation
	 *            list of triangles which are already in an annotation
	 * @return plane annotation found or null
	 */
	protected static PlaneAnnotation analyseTrianglePlane(MeshCas cas, Triangle triangle,
			Set<Triangle> alreadyInAnnotation) {
		// List of already visited triangles for BFS
		Set<Triangle> visited = new HashSet<Triangle>();
		visited.add(triangle);

		// FIFO queue for triangles to visit for BFS
		LinkedList<Triangle> queue = new LinkedList<Triangle>();
		Edge[] edges = triangle.getEdges();
		for (int i = 0 ; i < edges.length ; ++i) {
			if (!edges[i].getIsSharpEdge()) {
				queue.addAll(triangle.getNeighborsOfEdge(edges[i]));
			}
		}

		PlaneAnnotation annotation = null;

		while (!queue.isEmpty()) {
			Triangle currNeighbor = queue.pop();
			visited.add(currNeighbor);
			if (alreadyInAnnotation.contains(currNeighbor))
				continue;

			// check if triangles belong to same plane annotation (normals have same direction with some tolerance)
			boolean isSameNormal = planeAngleWithinTolerance(triangle.getNormalVector(),currNeighbor.getNormalVector());

			if (isSameNormal) {
				// found two triangles representing a plane
				if (annotation == null) {
					annotation = new PlaneAnnotation(cas.getCurvatures(), cas.getModel());

					synchronized (annotation.getMesh().getTriangles()) {
						annotation.getMesh().getTriangles().add(triangle);
					}
					alreadyInAnnotation.add(triangle);
				}

				synchronized (annotation.getMesh().getTriangles()) {
					annotation.getMesh().getTriangles().add(currNeighbor);
				}
				alreadyInAnnotation.add(currNeighbor);

				Edge[] nEdges = currNeighbor.getEdges();
				for (int i = 0 ; i < nEdges.length ; ++i) {
					if (!edges[i].getIsSharpEdge()) {
						List<Triangle> neighbors = currNeighbor.getNeighborsOfEdge(nEdges[i]);
						for (int j = 0 ; j < neighbors.size() ; ++j) {
							synchronized (annotation.getMesh()) {
								synchronized (annotation.getMesh().getTriangles()) {
									if (visited.contains(neighbors.get(j)) || 
											annotation.getMesh().getTriangles().contains(neighbors.get(j))) {
										continue;
									}
								}
							}
							queue.add(neighbors.get(j));
						}
					}
				}
			}
		}
		return annotation;
	}

	/**
	 * Sets the primitive type of a vertex
	 * 
	 * @param curvatures
	 *            list of all curvature values for the vertices of the CAD model
	 * @param v
	 *            vertex to analyze
	 */
	static void analyseVertex(HashMap<Vertex, Curvature> curvatures, Vertex v) {

		Curvature c = curvatures.get(v);
		c.setPrimitiveType(getPrimitiveType(curvatures, v));
	}

	/**
	 * Combines similar neighboring annotations into one annotation. PLANE, SPHERE.CONVEX,
	 * SPHERE.CONCAVE, CONE.CONVEX, CONE.CONCAVE annotations can be merged by this routine.
	 * Additional checking is done on merging of two PLANE annotations based on the angle
	 * of their normal vectors and the threshold {@code PLANE_COMBINE_DEGREE}. The method
	 * also returns in the parameter set {@code toRefit} the list of triangle annotations
	 * needed to be refit because they have been changed (merged).
	 * 
	 * @param cas
	 *            MeshCas
	 * @param annotations
	 *            list of all annotations to check
	 * @param toRefit
	 *            annotations to refit due to merging.
	 */
	@SuppressWarnings("rawtypes")
	private static void combineSameNeighboringAnnotations(MeshCas cas,
			List<Annotation> annotations, Set<Annotation> toRefit) {
		// Combine neighboring annotations which were previously divided by smaller annotations and
		// are now neighbors
		Set<Annotation> toRemove = new HashSet<Annotation>();
		for (Iterator<Annotation> it = annotations.iterator(); it.hasNext();) {
			Annotation a = it.next();
			if (a instanceof PrimitiveAnnotation) {
				PrimitiveAnnotation pa = (PrimitiveAnnotation) a;

				@SuppressWarnings("unchecked")
				Set<PrimitiveAnnotation> neighborAnnotations = pa.getNeighborAnnotations(cas, pa.getClass());
				for (PrimitiveAnnotation a1 : neighborAnnotations) {
					if (toRemove.contains(a1))
						continue;
					if (pa instanceof ConeAnnotation 
							&& ((ConeAnnotation) pa).isConcave() != ((ConeAnnotation) a1).isConcave()) {
						continue;
					} else if (pa instanceof SphereAnnotation
							&& ((SphereAnnotation) pa).isConcave() != ((SphereAnnotation) a1).isConcave()) {
						continue;
					} else if (pa instanceof PlaneAnnotation
							&& Math.toDegrees(((PlaneAnnotation) pa).getPlaneNormal().angle(((PlaneAnnotation) a1).getPlaneNormal())) > PLANE_COMBINE_DEGREE) {
						// Two planes, but angle bigger than PLANE_COMBINE_DEGREE degree
						continue;
					}

					// Found same annotations
					toRemove.add(pa);

					synchronized (a1.getMesh().getTriangles()) {
						a1.getMesh().getTriangles().addAll(pa.getMesh().getTriangles());
					}
					if (toRefit != null)
						toRefit.add(a1);
					break;
				}
			}
		}
		synchronized (annotations) {
			annotations.removeAll(toRemove);
		}
	}

	/**
	 * Combines small annotations with their surrounding larger ones and performs
	 * a filtered cleaning of possible unwanted artifacts. The annotations types 
	 * that can be merged into other annotations are:
	 * 		small PLANE -> big SPHERE / CONE
	 * 		small CONE.CONVEX -> big CONE.CONCAVE (artifact)
	 * 		small CONE.CONCAVE -> big CONE.CONVEX (artifact)
	 * 		small CONE -> big SPHERE (possible artifact) 
	 * 		reset very small CONE / SPHERE annotations (artifact)
	 * 
	 * @param cas
	 *            MeshCas
	 */
	@SuppressWarnings({ "rawtypes", "unchecked" })
	private static void combineSmallAnnotations(MeshCas cas) {
		Set<Annotation> toRemove = new HashSet<Annotation>();
		for (Iterator<Annotation> it = cas.getAnnotations().iterator(); it.hasNext();) {
			Annotation a = it.next();
			if (a instanceof PrimitiveAnnotation) {
				PrimitiveAnnotation pa = (PrimitiveAnnotation) a;

				// zero area, remove
				if (pa.getArea() == 0) {
					toRemove.add(a);
					continue;
				}

				Set<PrimitiveAnnotation> neighborAnnotations = pa.getNeighborAnnotations(cas,
						PrimitiveAnnotation.class);

				// check for cone or sphere annotation which has less than 5 triangles
				if (neighborAnnotations.size() > 0 && (a instanceof ConeAnnotation || a instanceof SphereAnnotation)
						&& ((PrimitiveAnnotation) a).getMesh().getTriangles().size() < 5) {
					toRemove.add(a);

					PrimitiveAnnotation a1 = neighborAnnotations.iterator().next();

					synchronized (a1.getMesh().getTriangles()) {
						a1.getMesh().getTriangles().addAll(pa.getMesh().getTriangles());
					}
					continue;
				}

				// merge small annotations here
				for (PrimitiveAnnotation a1 : neighborAnnotations) {

					boolean merge = false;
					if (pa instanceof ConeAnnotation
							&& a1 instanceof ConeAnnotation
							&& ((ConeAnnotation) pa).isConcave() != ((ConeAnnotation) a1).isConcave()) {
						// merge pa into a1 if area of pa is significantly smaller than a1, because
						// neighboring cones, one convex and one concave with different size are
						// very unlikely.
						if (pa.getArea() < a1.getArea() * 0.30f)
							merge = true;
					}

					if (!merge && pa instanceof PlaneAnnotation && a1 instanceof ConeAnnotation) {
						if (pa.getArea() < a1.getArea() * 0.05f) {
							// Found a very small plane annotation neighboring a cone annotation. Check if
							// plane annotation should be part of cone annotation by checking if plane
							// normal is perpendicular to the generating axis.
							// Primitives are not fitted yet, therefore we need to compare triangles on
							// the edge for approximately the same normal vector
							HashMultimap<Triangle, Triangle> edgeTriangles = HashMultimap.create(a1.getMesh().getTriangles().size(), 2);
							Set<Vertex> edgeVertices = new HashSet<Vertex>();
							pa.getNeighborEdge(cas, a1, edgeVertices, edgeTriangles);

							// check for a triangle pair where the angle between triangle normals is
							// bigger or equal to 180 degree with a tolerance of twice the plane combine degree.
							for (Triangle t : edgeTriangles.keySet()) {
								Set<Triangle> partnerSet = edgeTriangles.get(t);
								for (Triangle partner : partnerSet) {
									float angle = (float)(Math.toDegrees(t.getNormalVector().angle(partner.getNormalVector())));
									if (angle <= (2 * PLANE_COMBINE_DEGREE) || angle >= 180f - (2 * PLANE_COMBINE_DEGREE)) {
										merge = true;
										break;
									}
								}
							}
						}
						else if (pa.getArea() >= a1.getArea() * 0.05f && pa.getArea() < a1.getArea() * 0.30f) {
							// found a small plane annotation neighboring cone annotation. Check if
							// plane annotation should be part of cone annotation by checking if plane
							// normal is perpendicular to the generating axis.
							// Primitives are not fitted yet, therefore we need to compare triangles on
							// the edge for approximately the same normal vector
							HashMultimap<Triangle, Triangle> edgeTriangles = HashMultimap.create(a1.getMesh().getTriangles().size(), 2);
							Set<Vertex> edgeVertices = new HashSet<Vertex>();
							pa.getNeighborEdge(cas, a1, edgeVertices, edgeTriangles);

							// check for a triangle pair where the angle between triangle normals is
							// bigger or equal to 180 degree with a tolerance of twice the plane combine degree.
							for (Triangle t : edgeTriangles.keySet()) {
								Set<Triangle> partnerSet = edgeTriangles.get(t);
								for (Triangle partner : partnerSet) {
									float angle = (float)(Math.toDegrees(t.getNormalVector().angle(partner.getNormalVector())));
									if (angle <= PLANE_TOLERANCE || angle >= 180f - PLANE_TOLERANCE) {
										merge = true;
										break;
									}
								}
							}
						}
					}
					
					if (!merge && pa instanceof PlaneAnnotation && (a1 instanceof SphereAnnotation || a1 instanceof ConeAnnotation)
							&& pa.getArea() < a1.getArea() * 0.02f) {
						// found very small plane (covering less than 2 percent of neighboring sphere or cone annotation)
						// so merge it into the big annotation as these results are just unwanted artifacts after the region growing
						merge = true;
					}
					
					if (!merge && pa instanceof ConeAnnotation && a1 instanceof SphereAnnotation 
							&& ((ConeAnnotation)pa).isConcave() == ((SphereAnnotation)a1).isConcave()
							&& pa.getArea() < a1.getArea() * 0.05f) {
						// found a small cone annotation next to a sphere annotation of the same convexity, so merge them together and 
						// this would be reverted at the sphere checking time if the fit error is too big
						merge = true;
					}
					
					if (merge) {
						toRemove.add(a);
						synchronized (a1.getMesh().getTriangles()) {
							a1.getMesh().getTriangles().addAll(pa.getMesh().getTriangles());
						}
						break;
					}
				}
			}
		}
		synchronized (cas.getAnnotations()) {
			cas.getAnnotations().removeAll(toRemove);
		}
	}

	/**
	 * Performs a geometric fit of all the {@code annotations} passed by as
	 * an argument. The fitted geometric shapes are:
	 * 		PLANE			: purple-white
	 * 		SPHERE CONVEX	: red
	 * 		SPHERE CONCAVE	: green
	 * 		CONE CONVEX		: yellow
	 * 		CONE CONCAVE	: turquoise	  
	 * 
	 * @param annotations
	 *            list of annotations to search through
	 * @return subset of annotations list of the given annotation type
	 */
	@SuppressWarnings("rawtypes")
	private static List<? extends PrimitiveAnnotation> fitAnnotations(
			final Collection<Annotation> annotations) {

		List<Callable<Void>> threads = new LinkedList<Callable<Void>>();
		final List<PrimitiveAnnotation> failedFittings = new LinkedList<PrimitiveAnnotation>();
		for (final Annotation a : annotations) {
			if (a instanceof PrimitiveAnnotation) {
				threads.add(new Callable<Void>() {
					@Override
					public Void call() throws Exception {
						if (!((PrimitiveAnnotation) a).fit()) {
							synchronized (failedFittings) {
								failedFittings.add((PrimitiveAnnotation) a);
							}
						} else
							((PrimitiveAnnotation) a).updateAnnotationArea();
						return null;
					}
				});
			}
		}

		ThreadPool.executeInPool(threads);
		return failedFittings;
	}

	/**
	 * Determines the primitive type of a given vertex by checking its curvature properties.
	 * The primitive type returned can be PLANE, SPHERE.CONVEX, SPHERE.CONCAVE, CONE.CONVEX,
	 * CONE.CONCAVE.
	 * 
	 * @param curvatures
	 *            list of curvatures for vertices
	 * @param v
	 *            vertex to analyze
	 * @return primitive type of vertex
	 */
	private static PrimitiveType getPrimitiveType(HashMap<Vertex, Curvature> curvatures, Vertex v) {

		Curvature c = curvatures.get(v);
		
		if (Math.abs(c.getSaturation()) < 0.20)
			return PrimitiveType.PLANE;

		float hue = c.getHue();

		if (hue < 15 * Math.PI / 180)
			return PrimitiveType.SPHERE_CONVEX;
		else if (hue >= 15 * Math.PI / 180 && hue < 75 * Math.PI / 180)
			return PrimitiveType.CONE_CONVEX;
		else if (hue >= 75 * Math.PI / 180 && hue < 150 * Math.PI / 180
				|| hue >= 230 * Math.PI / 180)
			return PrimitiveType.SPHERE_CONCAVE;
		else
			return PrimitiveType.CONE_CONCAVE;
	}

	/**
	 * Determines the primitive type given the property counts of each primitive type. 
	 * The biggest number indicates the primitive type.
	 * 
	 * @param planeCnt
	 *            number of plane vertices for triangle
	 * @param sphereConvexCnt
	 *            number of convex sphere vertices for triangle
	 * @param sphereConcaveCnt
	 *            number of concave sphere vertices for triangle
	 * @param coneConvexCnt
	 *            number of convex cone vertices for triangle
	 * @param coneConcaveCnt
	 *            number of concave cone vertices for triangle
	 * @return primitive type indicated by biggest number value
	 */
	private static PrimitiveType getTypeForCounts(int planeCnt, int sphereConvexCnt,
			int sphereConcaveCnt, int coneConvexCnt, int coneConcaveCnt) {
		int max = Math.max(
				planeCnt,
				Math.max(sphereConvexCnt,
						Math.max(sphereConcaveCnt, Math.max(coneConvexCnt, coneConcaveCnt))));

		if (max == planeCnt) {
			return PrimitiveType.PLANE;
		} else if (max == sphereConvexCnt) {
			return PrimitiveType.SPHERE_CONVEX;
		} else if (max == sphereConcaveCnt) {
			return PrimitiveType.SPHERE_CONCAVE;
		} else if (max == coneConvexCnt) {
			return PrimitiveType.CONE_CONVEX;
		} else
			return PrimitiveType.CONE_CONCAVE;
	}

	/**
	 * Checks if two plane annotations represent the same plane by comparing the 
	 * angle between their surface angles to the threshold tolerance value from
	 * {@code PLANE_TOLERANCE}
	 * 
	 * @param a1
	 *            plane 1
	 * @param a2
	 *            plane 2
	 * @return true if plane annotations should be combined into one
	 */
	@SuppressWarnings({ "rawtypes", "unused" })
	private static boolean isSamePlane(PrimitiveAnnotation a1, PrimitiveAnnotation a2) {
		if (!(a1 instanceof PlaneAnnotation && a2 instanceof PlaneAnnotation))
			return true;
		return planeAngleWithinTolerance(((PlaneAnnotation) a1).getPlaneNormal(),
				((PlaneAnnotation) a2).getPlaneNormal());
	}

	/**
	 * Merges given set of annotations with neighboring annotations. The neighbor annotation to merge
	 * into is chosen by it's size. If an annotation of the same type is one of the neighbors, this
	 * one is preferred over the bigger one.
	 * 
	 * @param cas
	 *            MeshCas
	 * @param toMerge
	 *            set of annotations to merge
	 * @param toRefit
	 *            returns a set of annotations which need to be refitted because their 
	 *            content has changed
	 */
	@SuppressWarnings("rawtypes")
	private static void mergeWithNeighbors(MeshCas cas, Set<PrimitiveAnnotation> toMerge,
			Set<Annotation> toRefit) {

		for (PrimitiveAnnotation pa : toMerge) {

			@SuppressWarnings("unchecked")
			Set<PrimitiveAnnotation> neighborAnnotations = pa.getNeighborAnnotations(cas,
					PrimitiveAnnotation.class);
			float maxArea = 0;
			PrimitiveAnnotation bestNeighbor = null;
			boolean bestIsSameInstance = false;
			// find the biggest neighbor annotation or one of the same type
			for (PrimitiveAnnotation a1 : neighborAnnotations) {

				if (toMerge.contains(a1))
					continue;

				if (pa instanceof PlaneAnnotation
						&& a1 instanceof PlaneAnnotation
						&& Math.acos(((PlaneAnnotation) pa).getPlaneNormal().dot(
								((PlaneAnnotation) a1).getPlaneNormal())) > PLANE_COMBINE_DEGREE
								* Math.PI / 180.0) {
					// Two planes, but angle bigger than XX degree
					continue;
				}

				boolean sameConvexity = pa.getClass() == a1.getClass();
				if (sameConvexity && pa instanceof ConeAnnotation && a1 instanceof ConeAnnotation) {
					sameConvexity = ((ConeAnnotation) pa).isConcave() == ((ConeAnnotation) a1)
							.isConcave();
				}
				if (sameConvexity && pa instanceof SphereAnnotation
						&& a1 instanceof SphereAnnotation) {
					sameConvexity = ((SphereAnnotation) pa).isConcave() == ((SphereAnnotation) a1)
							.isConcave();
				}

				if (sameConvexity) {
					bestIsSameInstance = true;
					if (a1.getArea() > maxArea) {
						bestNeighbor = a1;
						maxArea = a1.getArea();
					}
				} else if (bestIsSameInstance)
					continue;

				if (a1.getArea() > maxArea) {
					bestNeighbor = a1;
					maxArea = a1.getArea();
				}
			}
			if (bestNeighbor != null) {
				// merge
				synchronized (bestNeighbor.getMesh().getTriangles()) {
					bestNeighbor.getMesh().getTriangles().addAll(pa.getMesh().getTriangles());
				}
				toRefit.add(bestNeighbor);
			} else {
				// no neighbor found, add as new annotation
				synchronized (cas.getAnnotations()) {
					cas.getAnnotations().add(pa);
				}
				toRefit.add(pa);
			}
		}
	}

	/**
	 * Checks if angle between the given surface normals is within <tt>PLANE_TOLERANCE</tt>.
	 * 
	 * @param norm1
	 *            surface normal 1
	 * @param norm2
	 *            surface normal 2
	 * @return true if angle is within tolerance
	 */
	private static boolean planeAngleWithinTolerance(Vector3f norm1, Vector3f norm2) {
		float angle = (float)(Math.toDegrees(norm1.angle(norm2)));
		return (angle <= PLANE_TOLERANCE || angle >= 180 - PLANE_TOLERANCE);
	}

	/**
	 * Gets the primitive type of a triangle
	 * 
	 * @param curvatures
	 *            curvature property for each vertex
	 * @param triangle
	 *            triangle to analyze
	 * @return primitive type of triangle
	 */
	private PrimitiveType getTrianglePrimitiveType(HashMap<Vertex, Curvature> curvatures,
			Triangle triangle) {
		return getTrianglePrimitiveType(curvatures, triangle, true);
	}

	/**
	 * Gets primitive type of triangle by optionally averaging over neighboring triangles. First
	 * determines primitive type of triangle and then checks if neighboring triangles are of the
	 * same type. If triangle is totally different then type of neighboring triangles is returned.
	 * 
	 * @param curvatures
	 *            curvature property for each vertex
	 * @param triangle
	 *            triangle to analyze
	 * @param checkNeighbors
	 *            set to true if smoothing by neighbor triangles should be enabled
	 * @return primitive type for triangle
	 */
	private PrimitiveType getTrianglePrimitiveType(HashMap<Vertex, Curvature> curvatures,
			Triangle triangle, boolean checkNeighbors) {
		if (!checkNeighbors && trianglePrimitiveTypeMap.containsKey(triangle))
			return trianglePrimitiveTypeMap.get(triangle);

		// determine type of triangle
		int planeCnt = 0;
		int sphereConvexCnt = 0;
		int sphereConcavCnt = 0;
		int coneConvexCnt = 0;
		int coneConcavCnt = 0;
		for (Vertex v : triangle.getPosition()) {
			if (!v.isSharpVertex()) {
				Curvature c = curvatures.get(v);
				if (c.getPrimitiveType() == PrimitiveType.PLANE)
					planeCnt++;
				else if (c.getPrimitiveType() == PrimitiveType.SPHERE_CONVEX)
					sphereConvexCnt++;
				else if (c.getPrimitiveType() == PrimitiveType.SPHERE_CONCAVE)
					sphereConcavCnt++;
				else if (c.getPrimitiveType() == PrimitiveType.CONE_CONVEX)
					coneConvexCnt++;
				else if (c.getPrimitiveType() == PrimitiveType.CONE_CONCAVE)
					coneConcavCnt++;
			}
		}

		if (checkNeighbors) {
			// smooth type by neighbors
			Edge[] edges = triangle.getEdges();
			for (int i = 0 ; i < edges.length ; ++i) {
				if (!edges[i].getIsSharpEdge()) {
					List<Triangle> neighbors = triangle.getNeighborsOfEdge(edges[i]);
					for (int j = 0 ; j < neighbors.size() ; ++j) {
						PrimitiveType type = getTrianglePrimitiveType(curvatures, neighbors.get(j), false);
						if (type == PrimitiveType.PLANE)
							planeCnt += 1;
						else if (type == PrimitiveType.SPHERE_CONVEX)
							sphereConvexCnt += 1;
						else if (type == PrimitiveType.SPHERE_CONCAVE)
							sphereConcavCnt += 1;
						else if (type == PrimitiveType.CONE_CONVEX)
							coneConvexCnt += 1;
						else if (type == PrimitiveType.CONE_CONCAVE)
							coneConcavCnt += 1;
					}
				}
			}
		}

		trianglePrimitiveTypeMap.put(triangle, getTypeForCounts(planeCnt, sphereConvexCnt, sphereConcavCnt, coneConvexCnt, coneConcavCnt));
		return trianglePrimitiveTypeMap.get(triangle);
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.analyser.MeshAnalyser#processStart(edu.tum.cs.vis.model.uima.cas.MeshCas)
	 */
	@SuppressWarnings("rawtypes")
	@Override
	public void processStart(final MeshCas cas) {
		allVertices = cas.getModel().getVertices();
		allTriangles = cas.getModel().getTriangles();
		allRegions = cas.getModel().getRegions();
		float totalSurfaceArea = 0.0f;
		
		// compute the total surface area of the CAD model
		for (int i = 0 ; i < allTriangles.size() ; ++i) {
			totalSurfaceArea += allTriangles.get(i).getArea();
		}
		// unscale total surface area
		totalSurfaceArea /= cas.getModel().getScale();
		
		// set primitive type for all vertices
		List<Callable<Void>> threads = new LinkedList<Callable<Void>>();

		final int interval = 500;

		// Add approx 10% of triangles because we don't know yet how much work we have to do after
		// we found the primitives and have to merge them
		itemsToElaborate = (int) (allVertices.size() + allTriangles.size() + (allTriangles.size() * 0.1));

		for (int start = 0; start < allVertices.size(); start += interval) {
			final int st = start;
			threads.add(new Callable<Void>() {

				@Override
				public Void call() throws Exception {
					int end = Math.min(st + interval, allVertices.size());
					for (int i = st; i < end; i++) {
						analyseVertex(cas.getCurvatures(), allVertices.get(i));
						itemsElaborated.incrementAndGet();
					}
					return null;
				}

			});
		}

		ThreadPool.executeInPool(threads);

		// set primitive type for all triangles
		final Set<Triangle> alreadyInAnnotation = new HashSet<Triangle>();

		final Set<Triangle> toElaborate = new HashSet<Triangle>();
		for (Region r : allRegions) {
			float areaRegionUnscaled = r.getAreaOfRegion() / cas.getModel().getScale();
			for (Triangle t : r.getTriangles()) {
				if (alreadyInAnnotation.contains(t))
					continue;
				PlaneAnnotation found = analyseTrianglePlane(cas, t, alreadyInAnnotation);
				if (found != null) {
					// check if the plane is big enough to be stored as a valid annotation
					if (found.getArea() > (0.05 * totalSurfaceArea) || found.getArea() > (0.3 * areaRegionUnscaled) 
							|| found.getMesh().getTriangles().size() > 20) {
						synchronized (cas.getAnnotations()) {
							cas.addAnnotation(found);
						}
						itemsElaborated.incrementAndGet();
					}
					else {
						// plane too small, maybe part of another annotation
						// give another chance in analyseTriangle
						 toElaborate.addAll(found.getMesh().getTriangles());
						 alreadyInAnnotation.removeAll(toElaborate);
					}
				} else {
					toElaborate.add(t);
				}
			}
		}
		
		for (PlaneAnnotation pa : cas.findAnnotations(PlaneAnnotation.class)) {
			// there may be some triangles added to toElaborate but afterwards added into a
			// annotation because of region grow
			toElaborate.removeAll(pa.getMesh().getTriangles());
			alreadyInAnnotation.addAll(pa.getMesh().getTriangles());
		}
		
		// Split number of triangles into alreadyAdded and to elaborate
		itemsToElaborate = (int) (allVertices.size() + alreadyInAnnotation.size() + toElaborate.size() + (allTriangles.size() * 0.1));

		alreadyInAnnotation.removeAll(toElaborate);
		
		for (Region r : allRegions) {
			for (Triangle t : r.getTriangles()) {
				if (toElaborate.contains(t)) {
					analyseTriangle(cas, t, alreadyInAnnotation);
					itemsElaborated.incrementAndGet();
				}
			}
		}
		
		fitAnnotations(cas.getAnnotations());
		
		logger.debug("Pre-combining neighbors ...");
		// now merge small annotations into bigger ones
		combineSmallAnnotations(cas);
		// and combine neighbors of same type
		combineSameNeighboringAnnotations(cas, cas.getAnnotations(), null);

		final Set<PrimitiveAnnotation> failedFittings = new HashSet<PrimitiveAnnotation>();
		failedFittings.addAll(fitAnnotations(cas.getAnnotations()));

		Set<Annotation> toRefit = new HashSet<Annotation>();
		if (failedFittings.size() > 0) {
			logger.debug("Merging failed fittings into neighbors (" + failedFittings.size() + ")");
			mergeWithNeighbors(cas, failedFittings, toRefit);
			combineSameNeighboringAnnotations(cas, cas.getAnnotations(), toRefit);
		}
		fitAnnotations(toRefit);
		toRefit.clear();

		// now check if all sphere annotations are spheres or if they should be 
		// cones by evaluating the fit error
		// fit error should be smaller than 0.005 for good fitted spheres

		logger.debug("Checking spheres ...");

		Set<Annotation> toAdd = new HashSet<Annotation>();
		Set<Annotation> toRemove = new HashSet<Annotation>();

		// Merge spheres with small area coverage into neighbor
		Set<PrimitiveAnnotation> smallSpheres = new HashSet<PrimitiveAnnotation>();

		for (Iterator<Annotation> it = cas.getAnnotations().iterator(); it.hasNext();) {
			Annotation a = it.next();
			if (a instanceof SphereAnnotation) {
				SphereAnnotation sa = (SphereAnnotation) a;

				if (sa.getAreaCoverage() < 0.01) {
					smallSpheres.add(sa);
					continue;
				}

				if (sa.getSphere().getFitError() < 0.005)
					continue;

				// create a temporary cone annotation and let it fit:
				ConeAnnotation tmp = new ConeAnnotation(cas.getCurvatures(), cas.getModel(), sa.isConcave());
				synchronized (tmp.getMesh().getTriangles()) {
					tmp.getMesh().getTriangles().addAll(sa.getMesh().getTriangles());
				}
				if (!tmp.fit())
					continue;
				tmp.updateAnnotationArea();
				if (tmp.getCone().getFitError() < sa.getSphere().getFitError()) {
					toRemove.add(a);
					toAdd.add(tmp);
					// logger.debug("Changing sphere annotation to cone because fit error is smaller: "
					// + tmp.getCone().getFitError() + "<" + sa.getSphere().getFitError());
				}
			}
		}

		synchronized (cas.getAnnotations()) {
			cas.getAnnotations().removeAll(smallSpheres);
			cas.getAnnotations().removeAll(toRemove);
			cas.getAnnotations().addAll(toAdd);
		}
		mergeWithNeighbors(cas, smallSpheres, toRefit);
		
		logger.debug("Combining neighbors ...");
		// combine neighboring annotations which were previously different types before changing sphere to cone
		combineSameNeighboringAnnotations(cas, cas.getAnnotations(), toRefit);
		
		fitAnnotations(toRefit);

		// everything has been processed at this point
		itemsElaborated.set(itemsToElaborate);
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.analyser.MeshAnalyser#updateProgress()
	 */
	@Override
	public void updateProgress() {
		if (allVertices != null && allVertices.size() > 0 && itemsToElaborate != 0) {
			setProgress((float) itemsElaborated.get() / (float) itemsToElaborate * 100.0f);
		}

	}
	
	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.analyser.MeshAnalyser#getLogger()
	 */
	@Override
	public Logger getLogger() {
		return logger;
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.analyser.MeshAnalyser#getName()
	 */
	@Override
	public String getName() {
		return "Primitive";
	}
}
