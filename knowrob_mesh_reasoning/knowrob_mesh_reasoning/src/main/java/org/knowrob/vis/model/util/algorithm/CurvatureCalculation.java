/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: 
 * 				Stefan Profanter - initial API and implementation, Year: 2012
 * 				Andrei Stoica - refactored original implementation 
 * 								during Google Summer of Code 2014
 ******************************************************************************/
package org.knowrob.vis.model.util.algorithm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicLong;

import javax.vecmath.Point3f;
import javax.vecmath.Vector3f;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.log4j.Logger;

import org.knowrob.utils.ThreadPool;
import org.knowrob.vis.model.Model;
import org.knowrob.vis.model.util.Curvature;
import org.knowrob.vis.model.util.Triangle;
import org.knowrob.vis.model.util.Vertex;
import org.knowrob.vis.model.util.Edge;
import org.knowrob.vis.model.util.ThresholdsReasoning;


/**
 * Accumulator class. Ported from trimesh2 (2.12).
 * 
 * @see <a href="https://github.com/fcole/qrtsc/tree/master/trimesh2">trimesh2</a>
 * 
 * @author Stefan Profanter
 */
abstract class ACCUM {
	/**
	 * The accumulator functor
	 */
	@SuppressWarnings("javadoc")
	public abstract void a(final Model m, HashMap<Vertex, Curvature> curvatures, Vertex v0,
			Vertex c, float w, Vertex v);
}

/**
 * Accumulator class for curvature accumulation. Ported from trimesh2 (2.12).
 * 
 * @see <a href="https://github.com/fcole/qrtsc/tree/master/trimesh2">trimesh2</a>
 * 
 * @author Stefan Profanter
 */
final class AccumCurv extends ACCUM {

	@Override
	public void a(final Model m, HashMap<Vertex, Curvature> curvatures, Vertex v0, Vertex c,
			float w, Vertex v) {
		Curvature curv = curvatures.get(v);
		Curvature curv0 = curvatures.get(v0);
		if (curv == null || curv0 == null)
			return;
		float ncurv[] = CurvatureCalculation.proj_curv(curv.getPrincipleDirectionMax(),
				curv.getPrincipleDirectionMin(), curv.getCurvatureMax(), 0, curv.getCurvatureMin(),
				curv0.getPrincipleDirectionMax(), curv0.getPrincipleDirectionMin());
		Vector3f tmp = new Vector3f(ncurv);
		tmp.scale(w);
		c.add(tmp);
	}
}

/**
 * Class which implements methods to estimate the curvature of the vertices
 * of a CAD model's mesh. The algorithms computing the curvatures within this
 * class are mainly ported from the open source trimesh2 (2.12) 
 * (Szymon Rusinkiewicz, Princeton University) library. The curvature computation
 * is based on the normals variation of the vertices and makes use of the second
 * normal form and its least square solution computed with the LU decomposition.
 * 
 * @see <a href="https://github.com/fcole/qrtsc/tree/master/trimesh2">trimesh2</a>
 * @see <a href="http://gfx.cs.princeton.edu/pubs/_2004_ECA/curvpaper.pdf">"Estimating Curvatures and Their Derivatives on Triangle Meshes", S. Rusinkiewicz</a>
 * 
 * @author Stefan Profanter
 * @author Andrei Stoica - minor refactor on existing implementation plus addition of new methods and functionality
 */
public class CurvatureCalculation {

	/**
	 * Log4J Logger
	 */
	private static final Logger LOGGER = Logger.getLogger(CurvatureCalculation.class);
	
	/**
	 * Calculates curvatures for each vertex of a model
	 * 
	 * @param curvatures
	 *            resulting vertex curvature mapping
	 * @param m
	 *            model to calculate curvature for
	 */
	public static void calculateCurvatures(HashMap<Vertex, Curvature> curvatures, final Model m) {
		if (curvatures == null || m == null || m.getVertices() == null || m.getVertices().size() == 0) {
			return;
		}
		calculateVoronoiArea(m);
		calculateCurvature(curvatures, m);
		setCurvatureHueSaturation(curvatures, m, ThresholdsReasoning.CURV_SMOOTHING);
		setCurvatureBoundaryValues(curvatures, m);
	}
	
	/**
	 * Calculate hue and saturation for curvature properties.
	 * 
	 * @param curvatures
	 *            maps curvatures to vertices
	 * @param smoothSigma
	 *            sigma value for smoothing the curvature. Set to 0 to disable smoothing.
	 * @param m
	 *            model needed to calculate hue saturation scale
	 */
	public static void setCurvatureHueSaturation(HashMap<Vertex, Curvature> curvatures, Model m,
			float smoothSigma) {
		if (smoothSigma > 0.0f) {
			float scaledSigma = smoothSigma * m.feature_size();
			diffuse_curv(m, curvatures, scaledSigma);
		}

		float cscale = 120.0f * typical_scale(curvatures, m);
		cscale = cscale * cscale;
		int nv = m.getVertices().size();
		for (int i = 0; i < nv; i++) {
			Curvature c = curvatures.get(m.getVertices().get(i));
			float H = 0.5f * (c.getCurvatureMax() + c.getCurvatureMin()); 
			// mean curvature -> don't set the mean curvature here as we want only the initial one to be stored and this method is public
			float K = c.getCurvatureMax() * c.getCurvatureMin(); 
			// Gaussian curvature -> don't set the Gaussian curvature here as we want only the initial one to be stored and this method is public
			float h = (float) (4.0f / 3.0f * Math
					.abs(Math.atan2(H * H - K, H * H * Math.signum(H))));
			float s = (float) ((2 / Math.PI) * Math.atan((2.0f * H * H - K) * cscale));
			c.setHue(h);
			c.setSaturation(s);
		}
	}
	
	/**
	 * Calculate curvature for each vertex of model
	 * 
	 * @param curvatures
	 *            resulting vertex curvature mapping
	 * @param m
	 *            model to calculate curvature for
	 */
	private static void calculateCurvature(final HashMap<Vertex, Curvature> curvatures,
			final Model m) {

		// Set up an initial coordinate system with min and max curvature as u and v per vertex
		for (int i = 0; i < m.getTriangles().size(); i++) {
			Triangle t = m.getTriangles().get(i);
			for (int j = 0; j < 3; j++) {
				Vertex v = t.getPosition()[j];
				Curvature c = new Curvature();
				curvatures.put(v, c);
				c.setPrincipleDirectionMax(new Vector3f(t.getPosition()[(j + 1) % 3]));
				c.getPrincipleDirectionMax().sub(v);
			}
		}
		
		for (Iterator<Vertex> it = m.getVertices().iterator(); it.hasNext();) {
			Vertex v = it.next();

			Curvature c = curvatures.get(v);
			if (c == null) {
				// vertex isn't referenced in any triangle
				it.remove();
				continue;
			}
			
			Vector3f tmp = new Vector3f();
			tmp.cross(c.getPrincipleDirectionMax(), v.getNormalVector());
			tmp.normalize();
			c.setPrincipleDirectionMax(tmp);

			tmp = new Vector3f();
			tmp.cross(v.getNormalVector(), c.getPrincipleDirectionMax());
			tmp.normalize();
			c.setPrincipleDirectionMin(tmp);

		}

		// Compute curvature per-face

		List<Callable<Void>> threads = new LinkedList<Callable<Void>>();

		final int interval = 500;

		for (int start = 0; start < m.getTriangles().size(); start += interval) {
			final int st = start;
			threads.add(new Callable<Void>() {

				@Override
				public Void call() throws Exception {
					int end = Math.min(st + interval, m.getTriangles().size());
					for (int i = st; i < end; i++) {
						calculateCurvatureForTriangle(curvatures, m.getTriangles().get(i));
					}
					return null;
				}

			});
		}

		ThreadPool.executeInPool(threads);
		
		// Compute principal directions and curvatures at each vertex
		for (int i = 0; i < m.getVertices().size(); i++) {
			Curvature c = curvatures.get(m.getVertices().get(i));
			Vector3f pdirRet[] = new Vector3f[2];
			float kRet[] = new float[2];

			diagonalize_curv(c.getPrincipleDirectionMax(), c.getPrincipleDirectionMin(),
					c.getCurvatureMax(), c.getCurvatureMinMax(), c.getCurvatureMin(), m
							.getVertices().get(i).getNormalVector(), pdirRet, kRet);
			c.setPrincipleDirectionMax(pdirRet[1]);
			c.setPrincipleDirectionMin(pdirRet[0]);
			c.setCurvatureMax(kRet[0]);
			c.setCurvatureMin(kRet[1]);
		}
	}

	/**
	 * Calculates curvature for the vertices of a triangle. Based on
	 * the trimesh2 (2.12) (Szymon Rusinkiewicz Princeton University)
	 * curvature calculation algorithm.
	 * 
	 * @see <a href="https://github.com/fcole/qrtsc/tree/master/trimesh2">trimesh2</a>
	 * 
	 * @param curvatures
	 *            vertex curvature mapping
	 * @param tri
	 *            triangle to calculate curvature for
	 */
	private static void calculateCurvatureForTriangle(HashMap<Vertex, Curvature> curvatures, Triangle tri) {
		// Edges
		Edge e[] = tri.getEdges();

		// N-T-B coordinate system per face
		Vector3f t = new Vector3f(e[0].getEdgeValue());
		t.normalize();
		Vector3f n = new Vector3f();
		n.cross(e[0].getEdgeValue(), e[1].getEdgeValue());
		//n.normalize();
		Vector3f b = new Vector3f();
		b.cross(n, t);
		b.normalize();

		// Estimate curvature based on variation of normals
		// along edges
		float m[] = { 0.0f, 0.0f, 0.0f };
		double w[][] = { { 0.0d, 0.0d, 0.0d }, { 0.0d, 0.0d, 0.0d }, { 0.0d, 0.0d, 0.0d } };
		for (int j = 0; j < 3; j++) {

			float u = e[j].getEdgeValue().dot(t);
			float v = e[j].getEdgeValue().dot(b);
			w[0][0] += u * u;
			w[0][1] += u * v;
			// w[1][1] += v*v + u*u;
			// w[1][2] += u*v;
			w[2][2] += v * v;
			Vector3f dn = new Vector3f(tri.getPosition()[(j + 2) % 3].getNormalVector());
			Vector3f dnn = new Vector3f(tri.getPosition()[(j + 1) % 3].getNormalVector());
			dn.sub(dnn);
			float dnu = dn.dot(t);
			float dnv = dn.dot(b);
			m[0] += dnu * u;
			m[1] += dnu * v + dnv * u;
			m[2] += dnv * v;
		}
		w[1][1] = w[0][0] + w[2][2];
		w[1][2] = w[0][1];

		RealMatrix coefficients = new Array2DRowRealMatrix(w, false);
		DecompositionSolver solver = new LUDecomposition(coefficients).getSolver();
		if (!solver.isNonSingular()) {
			return;
		}

		RealVector constants = new ArrayRealVector(new double[] { m[0], m[1], m[2] }, false);
		RealVector solution = solver.solve(constants);

		m[0] = (float) solution.getEntry(0);
		m[1] = (float) solution.getEntry(1);
		m[2] = (float) solution.getEntry(2);

		// flag for ill-solutioned system (NaN = no valid info on the curvatures => then consider 
		// the curvatures to be the ones of a simple plane)
		boolean okSolution = true;
		if (isNaN(m[0]) || isNaN(m[1]) || isNaN(m[2])) {
			okSolution = false;
		}
		
		// Push it back out to the vertices
		for (int j = 0; j < tri.getPosition().length; j++) {
			Vertex vj = tri.getPosition()[j];
			float c1, c12, c2;
			float[] ret;

			// guard against ill-conditioned systems (NaN curvatures are viewed as a simple planar
			// corresponding curves		
			if (okSolution) {
				ret = proj_curv(t, b, m[0], m[1], m[2], curvatures.get(vj).getPrincipleDirectionMax(), curvatures.get(vj).getPrincipleDirectionMin());
				if (isNaN(ret[0]) || isNaN(ret[1]) || isNaN(ret[2])) {
					c1 = 0.0f;
					c12 = 0.0f;
					c2 = 0.0f;
				}
				else {
					c1 = ret[0];
					c12 = ret[1];
					c2 = ret[2];
				}
			}
			else {
				c1 = 0.0f;
				c12 = 0.0f;
				c2 = 0.0f;
			}
			
			Curvature c = curvatures.get(vj);

			double wt = 0.0d;
			if (vj.getPointarea() != 0.0f) {
				if (j == 0)
					wt = tri.getCornerarea().x / vj.getPointarea();
				else if (j == 1)
					wt = tri.getCornerarea().y / vj.getPointarea();
				else
					wt = tri.getCornerarea().z / vj.getPointarea();
			}
			
			synchronized (c) {
				c.setCurvatureMax((float) (c.getCurvatureMax() + wt * c1));
				c.setCurvatureMinMax((float) (c.getCurvatureMinMax() + wt * c12));
				c.setCurvatureMin((float) (c.getCurvatureMin() + wt * c2));
			}
		}
	}

	/**
	 * Ported from trimesh2 (2.12) (Szymon Rusinkiewicz Princeton University)
	 * 
	 * Compute the area "belonging" to each vertex or each corner of a triangle (defined as Voronoi
	 * area restricted to the 1-ring of a vertex, or to the triangle).
	 * 
	 * @see <a href="https://github.com/fcole/qrtsc/tree/master/trimesh2">trimesh2</a>
	 * 
	 * @param m
	 *            Model to calculate Voronoi area for
	 */
	private static void calculateVoronoiArea(final Model m) {
		// Compute from faces

		List<Callable<Void>> threads = new LinkedList<Callable<Void>>();

		final int interval = 500;

		for (Vertex v : m.getVertices())
			v.setPointarea(0);

		/*for (int start = 0; start < m.getTriangles().size(); start++) {
			calculateVoronoiAreaForTriangle(m.getTriangles().get(start));
		};*/

		for (int start = 0; start < m.getTriangles().size(); start += interval) {
			final int st = start;
			threads.add(new Callable<Void>() {

				@Override
				public Void call() throws Exception {
					int end = Math.min(st + interval, m.getTriangles().size());
					for (int i = st; i < end; i++) {
						calculateVoronoiAreaForTriangle(m.getTriangles().get(i));
					}
					return null;
				}

			});
		}

		ThreadPool.executeInPool(threads);
	}
	
	/**
	 * Calculate Voronoi area for triangle and its vertices
	 * 
	 * @param t
	 *            triangle
	 */
	private static void calculateVoronoiAreaForTriangle(Triangle t) {
		// Edges
		Edge e[] = t.getEdges();

		// Compute corner weights
		Vector3f tmp = new Vector3f();
		tmp.cross(e[0].getEdgeValue(), e[1].getEdgeValue());
		float area = 0.5f * tmp.length();
		float l2[] = { e[0].getEdgeValue().lengthSquared(), e[1].getEdgeValue().lengthSquared(), e[2].getEdgeValue().lengthSquared() };
		float ew[] = { l2[0] * (l2[1] + l2[2] - l2[0]), l2[1] * (l2[2] + l2[0] - l2[1]),
				l2[2] * (l2[0] + l2[1] - l2[2]) };
		t.setCornerarea(new Vector3f());
		boolean ok = false;
		if (ew[0] <= 0.0f) {
			float d1 = e[0].getEdgeValue().dot(e[2].getEdgeValue());
			float d2 = e[0].getEdgeValue().dot(e[1].getEdgeValue());
			if (d1 != 0 && d2 != 0) {
				t.getCornerarea().y = -0.25f * l2[2] * area / d1;
				t.getCornerarea().z = -0.25f * l2[1] * area / d2;
				t.getCornerarea().x = area - t.getCornerarea().y - t.getCornerarea().z;
				ok = true;
			}
		} else if (ew[1] <= 0.0f) {
			float d1 = e[1].getEdgeValue().dot(e[0].getEdgeValue());
			float d2 = e[1].getEdgeValue().dot(e[2].getEdgeValue());
			if (d1 != 0 && d2 != 0) {
				t.getCornerarea().z = -0.25f * l2[0] * area / d1;
				t.getCornerarea().x = -0.25f * l2[2] * area / d2;
				t.getCornerarea().y = area - t.getCornerarea().z - t.getCornerarea().x;
				ok = true;
			}
		} else if (ew[2] <= 0.0f) {

			float d1 = e[2].getEdgeValue().dot(e[1].getEdgeValue());
			float d2 = e[2].getEdgeValue().dot(e[0].getEdgeValue());
			if (d1 != 0 && d2 != 0) {
				t.getCornerarea().x = -0.25f * l2[1] * area / d1;
				t.getCornerarea().y = -0.25f * l2[0] * area / d2;
				t.getCornerarea().z = area - t.getCornerarea().x - t.getCornerarea().y;
				ok = true;
			}
		}
		if (!ok) {
			float ewscale = 0.5f * area / (ew[0] + ew[1] + ew[2]);
			t.getCornerarea().x = ewscale * (ew[1] + ew[2]);
			t.getCornerarea().y = ewscale * (ew[2] + ew[0]);
			t.getCornerarea().z = ewscale * (ew[0] + ew[1]);
		}
		synchronized (t.getPosition()[0]) {
			t.getPosition()[0]
					.setPointarea(t.getPosition()[0].getPointarea() + t.getCornerarea().x);
		}
		synchronized (t.getPosition()[1]) {
			t.getPosition()[1]
					.setPointarea(t.getPosition()[1].getPointarea() + t.getCornerarea().y);
		}
		synchronized (t.getPosition()[2]) {
			t.getPosition()[2]
					.setPointarea(t.getPosition()[2].getPointarea() + t.getCornerarea().z);
		}
	}

	/**
	 * Given a curvature tensor, find its principal directions and curvatures. Makes sure that pdir1 and
	 * pdir2 are perpendicular to the normal vector of the triangle. The principle directions
	 * are stored in the array of 3D flaoting point vectors passed as parameter {@code pdir}. 
	 * The code is ported from trimesh2 (2.12) (Szymon Rusinkiewicz Princeton University).
	 * 
	 * @see <a href="https://github.com/fcole/qrtsc/tree/master/trimesh2">trimesh2</a> 
	 * 
	 * @param old_u
	 * 		initial direction associated with the maximum curvature
	 * @param old_v
	 * 		initial direction associated with the minimum curvature
	 * @param ku
	 * 		maximum curvature value - kMax
	 * @param kuv
	 * 		min-max curvature property - kMinMax
	 * @param kv
	 * 		minimum curvature value - kMin
	 * @param new_norm
	 * 		parameter that contains the initial normal vector of the
	 * 		vertex of interest and returns the direction after diagonalization
	 * 		of the curvature directional matrix
	 * @param pdir
	 * 		array of 3D vectors of the principle directions of the curvature
	 * @param k
	 * 		array of the minimum and maximum principle directions of the curvature
	 */
	private static void diagonalize_curv(final Vector3f old_u, final Vector3f old_v, float ku,
			float kuv, float kv, final Vector3f new_norm, Vector3f pdir[], float k[]) {
		Vector3f r_old_u = new Vector3f(), r_old_v = new Vector3f();
		rot_coord_sys(old_u, old_v, new_norm, r_old_u, r_old_v);

		float c = 1, s = 0, tt = 0;
		if (kuv != 0.0f) {
			// Jacobi rotation to diagonalize
			float h = 0.5f * (kv - ku) / kuv;
			tt = (float) ((h < 0.0f) ? 1.0f / (h - Math.sqrt(1.0f + h * h)) : 1.0f / (h + Math
					.sqrt(1.0f + h * h)));
			c = (float) (1.0f / Math.sqrt(1.0f + tt * tt));
			s = tt * c;
		}

		k[0] = ku - tt * kuv;
		k[1] = kv + tt * kuv;

		if (Math.abs(k[0]) >= Math.abs(k[1])) {
			r_old_u.scale(c);
			r_old_v.scale(s);
			r_old_u.sub(r_old_v);
			pdir[0] = r_old_u;
		} else {
			float kt = k[0];
			k[0] = k[1];
			k[1] = kt;

			r_old_u.scale(s);
			r_old_v.scale(c);
			r_old_u.add(r_old_v);
			pdir[0] = r_old_u;
		}
		pdir[1] = new Vector3f();
		pdir[1].cross(new_norm, pdir[0]);
	}

	/**
	 * Diffuses the curvatures across the mesh of the CAD model.
	 * Ported from trimesh2 (2.12) (Szymon Rusinkiewicz Princeton University)
	 * 
	 * @see <a href="https://github.com/fcole/qrtsc/tree/master/trimesh2">trimesh2</a>
	 * 
	 * @param m
	 *            the main model
	 * @param curvatures
	 *            map assigning a curvature to each vertex of the model
	 * @param sigma
	 *            smoothing sigma which determines the spread of the diffusion
	 */
	private static void diffuse_curv(Model m, HashMap<Vertex, Curvature> curvatures, float sigma) {
		int nv = m.getVertices().size();

		float invsigma2 = (float) (1.0f / Math.pow(sigma, 2));

		Vertex[] cflt = new Vertex[nv];
		// TODO #pragma omp parallel
		{
			// Thread-local flags
			Map<Vertex, Long> flags = new HashMap<Vertex, Long>(nv);
			AtomicLong flag_curr = new AtomicLong(0);

			// TODO #pragma omp for
			ACCUM accumCurv = new AccumCurv();
			for (int i = 0; i < nv; i++) {
				cflt[i] = new Vertex(0, 0, 0);
				diffuse_vert_field(m, curvatures, flags, flag_curr, accumCurv, i, invsigma2,
						cflt[i]);
			}

			// TODO #pragma omp for
			for (int i = 0; i < nv; i++) {
				Vertex v = m.getVertices().get(i);
				Curvature c = curvatures.get(v);
				Vector3f pdir[] = new Vector3f[] { c.getPrincipleDirectionMax(),
						c.getPrincipleDirectionMin() };
				float k[] = new float[] { c.getCurvatureMax(), c.getCurvatureMin() };
				diagonalize_curv(c.getPrincipleDirectionMax(), c.getPrincipleDirectionMin(),
						cflt[i].x, cflt[i].y, cflt[i].z, v.getNormalVector(), pdir, k);
				c.setPrincipleDirectionMax(pdir[0]);
				c.setPrincipleDirectionMin(pdir[1]);
				c.setCurvatureMax(k[0]);
				c.setCurvatureMin(k[1]);
			}
		} // #pragma omp parallel
	}

	/**
	 * Diffuses a vector field around one vertex, weighted by a Gaussian of width 
	 * {@code 1/sqrt(invsigma2)}. Ported from trimesh2 (2.12) (Szymon Rusinkiewicz Princeton University)
	 * 
	 * @see <a href="https://github.com/fcole/qrtsc/tree/master/trimesh2">trimesh2</a>
	 * 
	 * @param m
	 * 			CAD model analyzed
	 * @param curvatures
	 * 			curvatures HashMap with model vertices as keys and their cuvatures as values
	 * @param flags
	 * 			map of vertex to long value (default initialized with 0f and of the same length with the list of vertices of the model)
	 * @param flag_curr
	 * 			atomic long value
	 * @param accum
	 * 			accumulator
	 * @param invsigma
	 * 			spread of the Gaussian used in weighting
	 * @param flt
	 * 			3D vector field diffused based on the curvature
	 */
	@SuppressWarnings("javadoc")
	private static void diffuse_vert_field(final Model m, HashMap<Vertex, Curvature> curvatures,
			Map<Vertex, Long> flags, AtomicLong flag_curr, final ACCUM accum, int v,
			float invsigma2, Vertex flt) {
		Vertex vert = m.getVertices().get(v);
		if (vert.getNeighbors().size() == 0) {
			// flt.set(0, 0, 0);
			accum.a(m, curvatures, vert, flt, .5f, vert);
			return;
		}

		// flt.set(0, 0, 0);
		accum.a(m, curvatures, vert, flt, vert.getPointarea(), vert);
		float sum_w = vert.getPointarea();
		final Vector3f nv = vert.getNormalVector();

		long flag_curr_val = flag_curr.incrementAndGet();
		flags.put(vert, flag_curr_val);
		LinkedList<Vertex> boundary = new LinkedList<Vertex>();
		boundary.addAll(vert.getNeighbors());
		while (boundary.size() > 0) {
			Vertex n = boundary.pop();
			if (flags.get(n) != null && flags.get(n) == flag_curr_val)
				continue;
			flags.put(n, flag_curr_val);
			if (nv.dot(n.getNormalVector()) <= 0.0f)
				continue;
			// Gaussian weight
			float w = wt(n, vert, invsigma2);
			if (w == 0.0f)
				continue;
			// Downweight things pointing in different directions
			w *= nv.dot(n.getNormalVector());
			// Surface area "belonging" to each point
			w *= n.getPointarea();
			// Accumulate weight times field at neighbor
			accum.a(m, curvatures, vert, flt, w, n);
			sum_w += w;
			for (Vertex nn : n.getNeighbors()) {
				if (flags.get(nn) != null && flags.get(nn) == flag_curr_val)
					continue;
				boundary.push(nn);
			}
		}
		flt.scale(1 / sum_w);
	}

	/**
	 * Reprojects a curvature tensor from the basis spanned by old_u and old_v (which are assumed to
	 * be unit-length and perpendicular) to the new_u, new_v basis. Implementation
	 * ported from trimesh2 (2.12) (Szymon Rusinkiewicz Princeton University).
	 * The newly projected directions of the minimum and maxium curvatures are stored in
	 * {@code new_u}, and respectively, {@code new_v}.
	 * 
	 * @see <a href="https://github.com/fcole/qrtsc/tree/master/trimesh2">trimesh2</a>
	 * 
	 * @param old_u
	 * 		old principle curvature direction for maximum curvature
	 * @param old_v 
	 * 		old principle curvature direction for minimum curvature
	 * @param old_ku
	 * 		old maximum curvature value - kMax
	 * @param old_kuv
	 * 		old min-max curvature property - kMinMax
	 * @param old_kv
	 * 		old minimum curvature value - kMin
	 * @param new_u
	 * 		newly projected maximum curvature principle direction
	 * @param new_v 
	 * 		newly projected maximum curvature principle direction
	 * @return ret
	 * 		array of new updated values of the curvatures:
	 * 		ret[0] = maximum curvature kMax;
	 * 		ret[1] = min-max curvature kMinMax;
	 * 		ret[2] = minimum curvature kMin;
	 */
	protected static float[] proj_curv(final Vector3f old_u, final Vector3f old_v, float old_ku,
			float old_kuv, float old_kv, final Vector3f new_u, final Vector3f new_v) {
		Vector3f r_new_u = new Vector3f(), r_new_v = new Vector3f();
		Vector3f tmp = new Vector3f();
		tmp.cross(old_u, old_v);
		rot_coord_sys(new_u, new_v, tmp, r_new_u, r_new_v);

		float u1 = r_new_u.dot(old_u);
		float v1 = r_new_u.dot(old_v);
		float u2 = r_new_v.dot(old_u);
		float v2 = r_new_v.dot(old_v);
		float ret[] = new float[3];
		ret[0] = old_ku * u1 * u1 + old_kuv * (2.0f * u1 * v1) + old_kv * v1 * v1;
		ret[1] = old_ku * u1 * u2 + old_kuv * (u1 * v2 + u2 * v1) + old_kv * v1 * v2;
		ret[2] = old_ku * u2 * u2 + old_kuv * (2.0f * u2 * v2) + old_kv * v2 * v2;
		return ret;
	}

	/**
	 * Rotates a coordinate system to be perpendicular to the normal vector given
	 * as a parameter. Implementation ported from trimesh2 (2.12) (Szymon Rusinkiewicz Princeton University)
	 * 
	 * @see <a href="https://github.com/fcole/qrtsc/tree/master/trimesh2">trimesh2</a>
	 * 
	 * @param old_u
	 * 			old principle curvature direction associated with the maximum curvature
	 * @param old_v
	 * 			old principle curvature direction associated with the minimum curvature
	 * @param new_norm
	 * 			newly rotated norm of the vertex
	 * @param new_u
	 * 			newly rotated principle curvature direction associated with the maximum curvature
	 * @param new_v
	 * 			newly rotated principle curvature direction associated with the minimum curvature
	 */
	private static void rot_coord_sys(final Vector3f old_u, final Vector3f old_v,
			final Vector3f new_norm, Vector3f new_u, Vector3f new_v) {
		old_u.get(new_u);
		old_v.get(new_v);

		Vector3f old_norm = new Vector3f();
		old_norm.cross(old_u, old_v);
		float ndot = old_norm.dot(new_norm);
		if (ndot <= -1.0f) {
			new_u.scale(-1);
			new_v.scale(-1);
			return;
		}
		Vector3f perp_old = new Vector3f(new_norm);
		Vector3f tmp = new Vector3f(old_norm);
		tmp.scale(ndot);
		perp_old.sub(tmp);

		Vector3f dperp = new Vector3f(old_norm);
		dperp.add(new_norm);
		dperp.scale(1.0f / (1 + ndot));

		tmp = new Vector3f(dperp);
		tmp.scale(new_u.dot(perp_old));
		new_u.sub(tmp);

		tmp = new Vector3f(dperp);
		tmp.scale(new_v.dot(perp_old));
		new_v.sub(tmp);
	}

	/**
	 * Computes a "typical scale" for the mesh as 1% of the reciprocal of the 10-th
	 * percentile of the curvature
	 * 
	 * @param curvatures
	 *            curvature values for each vertex of the model
	 * @param m
	 *            model analyzed
	 * @return the typical scale value for the mesh
	 */
	private static float typical_scale(HashMap<Vertex, Curvature> curvatures, Model m) {
		float frac = 0.1f;
		float mult = 0.01f;

		int nv = m.getVertices().size();
		int nsamp = Math.min(nv, 500);
		if (nsamp == 0)
			return 0;

		float samples[] = new float[nsamp * 2];
		int idx = 0;
		for (int i = 0; i < nsamp; i++) {
			int ind = (int) (Math.random() * nv);

			samples[idx++] = Math.abs(curvatures.get(m.getVertices().get(ind)).getCurvatureMax());
			samples[idx++] = Math.abs(curvatures.get(m.getVertices().get(ind)).getCurvatureMin());
		}

		int which = (int) (frac * samples.length);
		Arrays.sort(samples);

		float f = 0;
		if (samples[which] == 0.0f || Float.isNaN(samples[which])) {
			f = mult * m.getBoundingSphere().getR();
			// LOGGER.warn("Couldn't determine typical scale. Using bsphere value: " + f + ".");
		} else {
			f = mult / samples[which];
		}
		return f;
	}

	/**
	 * Approximation to Gaussian scaling used in filtering the curvature information 
	 * between two vertices on the mesh. This is based on the distance between the two
	 * vertices Ported from trimesh2 (2.12) (Szymon Rusinkiewicz Princeton University)
	 * 
	 * @see <a href="https://github.com/fcole/qrtsc/tree/master/trimesh2">trimesh2</a>
	 * 
	 * @param p1
	 * 			position of first vertex
	 * @param p2
	 * 			position of second vertex
	 * @param invsigma2
	 * 			spread of the Gaussian filtering
	 * @return the Gaussian approximation
	 */
	private static float wt(Point3f p1, Point3f p2, float invsigma2) {
		float d2 = invsigma2 * p1.distanceSquared(p2);
		return (float) ((d2 >= 9.0f) ? 0.0f : Math.exp(-0.5f * d2));
		// return (d2 >= 25.0f) ? 0.0f : exp(-0.5f*d2);
	}
	
	/**
	 * Helper method that returns whether a floating point number is
	 * equal to NaN floating point equivalent or not
	 * 
	 * @param value
	 * @return true or false
	 */
	private static boolean isNaN(final float value) {
		return (value != value);
	}
	
	/**
	 * Sets the minimum and maximum curvature boundaries for the mean and
	 * Gaussian curvature properties estimated on the model. In doing so
	 * it statistically removes the outliers of the mean and Gaussian data
	 * points such that the coloring is distinctive for the majority of the 
	 * points (since because the curvatures are estimated there are cases
	 * when the outliers are really tilting the data severely and this distorts
	 * the coloring).
	 * 
	 * @param curvatures
	 * 			of the vertices of the model under analysis
	 * @param model
	 * 			under analysis
	 */
	private static void setCurvatureBoundaryValues(HashMap<Vertex, Curvature> curvatures, Model model) {
		float avgMean = 0f, varMean = 0f;
		float avgGauss = 0f, varGauss = 0f;
		List<Float> curvMeanData = new ArrayList<Float>();
		List<Float> curvGaussData = new ArrayList<Float>();
		LOGGER.debug("Computing model estimated curvature statistics ...");
		for (int i = 0 ; i < model.getVertices().size() ; ++i) {
			Curvature c = curvatures.get(model.getVertices().get(i));
			c.setMeanCurvature();
			c.setGaussCurvature();
			c.setInitialHue(c.getHue());
			c.setInitialSaturation(c.getSaturation());
			float curvMean = c.getMeanCurvature();
			float curvGauss = c.getGaussCurvature();
			curvMeanData.add(curvMean);
			curvGaussData.add(curvGauss);
			avgMean += curvMean;
			varMean += curvMean * curvMean;
			avgGauss += curvGauss;
			varGauss += curvGauss * curvGauss;
			if (model.getLowMeanCurvature() > c.getMeanCurvature()) {
				model.setLowMeanCurvature(c.getMeanCurvature());
			}
			if (model.getHighMeanCurvature() < c.getMeanCurvature()) {
				model.setHighMeanCurvature(c.getMeanCurvature());
			}
			if (model.getLowGaussCurvature() > c.getGaussCurvature()) {
				model.setLowGaussCurvature(c.getGaussCurvature());
			}
			if (model.getHighGaussCurvature() < c.getGaussCurvature()) {
				model.setHighGaussCurvature(c.getGaussCurvature());
			}
		}
		avgMean /= model.getVertices().size();
		avgGauss /= model.getVertices().size();
		varMean /= model.getVertices().size() - avgMean * avgMean;
		varGauss /= model.getVertices().size() - avgGauss * avgGauss;
		LOGGER.debug("\n\tAvg. Mean Curv.: " + avgMean + " Var. Mean Curv.: " + varMean +
				"\nAvg. Gaussian Curv.: " + avgGauss + " Var. Gaussian Curv.: " + varGauss);
		LOGGER.debug("Analyzing curvatures and removing spurious points statistics for coloring scale ...");
		float[] statsMean = {avgMean, varMean};
		float[] statsGauss = {avgGauss, varGauss};
		removeSpuriousCurvatures(curvMeanData, statsMean);
		removeSpuriousCurvatures(curvGaussData, statsGauss);
		model.setAvgMeanCurvature(statsMean[0]);
		model.setAvgGaussCurvature(statsGauss[0]);
		model.setVarMeanCurvature(statsMean[1]);
		model.setVarGaussCurvature(statsGauss[1]);
		LOGGER.debug("Curvature statistical computation finished." +
				"\n\tAvg. Mean Curv.: " + model.getAvgMeanCurvature() + " Var. Mean Curv.: " + model.getVarMeanCurvature() +
				"\nAvg. Gaussian Curv.: " + model.getAvgGaussCurvature() + " Var. Gaussian Curv.: " + model.getVarGaussCurvature());
	}
	
	/**
	 * Removes the spurious data from the a curvature data set and recomputes its mean
	 * and variance without the outliers. The new values are stored in the {@code statsBeg}
	 * array which contains the starting mean and variance of the samples:
	 * {@code statsBeg[0] = } mean, {@code statsBeg[1] = } variance. Implementation makes use
	 * of the Chauvenet's criterion to do this. After the outliers are ignored at least 90%
	 * of the initial data is contained in the cleaned data set.
	 * 
	 * @param curvData
	 * 			curvature data set
	 * @param statsBeg
	 * 			array of the beginning mean and variance of the data set:
	 * 			statsBeg[0] = mean of samples;
	 * 			statsBeg[1] = variance of samples;
	 */
	private static void removeSpuriousCurvatures(List<Float> curvData, float[] statsBeg) {
		float avg = statsBeg[0];
		float var = statsBeg[1];
		float distChauv = 0f;
		int iter = 0;
		float N = (float)(curvData.size());
		int iterLim = (int)((float)(curvData.size()) * 0.1f);
		sortCurv(curvData);
		while (distChauv < 0.5f && iter < iterLim && N > 0) {
			float dataPoint;
			float std = (float)(Math.sqrt(var));
			float dFront = estimateProbability(((curvData.get(0) - avg) / std)) * N;
			float dBack = estimateProbability(((curvData.get(curvData.size()-1) - avg) / std)) * N;
			if (dFront < dBack) {
				distChauv = dFront;
				dataPoint = curvData.get(0);
			}
			else {
				distChauv = dBack;
				dataPoint = curvData.get((int)(N - 1f));
			}
			if (distChauv < 0.5) {
				curvData.remove(dataPoint);
				N--;
				float[] stats = getVarianceOfData(curvData);
				avg = stats[0];
				var = stats[1];
			}
			iter++;
		}
		statsBeg[0] = avg;
		statsBeg[1] = var;
		LOGGER.debug(iter);
	}
	
	/**
	 * Sorts ascendingly & in-line a data list of floating point numbers.
	 * @param curvData
	 * 			data list to be sorted
	 */
	private static void sortCurv(final List<Float> curvData) {
		if (curvData != null) {
			Collections.sort(curvData, new Comparator<Float>() {
				@Override
				public int compare(Float c1, Float c2) {
					return c1.compareTo(c2);
				}
			});
		}
	}
	
	/**
	 * Estimates roughly the cumulative probability of a point to find
	 * itself at a normally distributed (mean 0, stddev 1) value based on bins of width
	 * 0.5. The bins are considered left and right of the mean 0.
	 * 
	 * @param value
	 * 			of interest for the point
	 * @return cumulative probability to be in the 0.5 bin given the value.
	 */
	private static float estimateProbability(float value) {
		value = Math.abs(value);
		if (value <= 0.5f) {
			return 0.191f;
		}
		if (value <= 1f) {
			return 0.150f;
		}
		if (value <= 1.5f) {
			return 0.092f;
		}
		if (value <= 2f) {
			return 0.044f;
		}
		if (value <= 2.5f) {
			return 0.017f;
		}
		if (value <= 3f) {
			return 0.005f;
		}
		if (value <= 3.5f) {
			return 0.001f;
		}
		return 0f;
	}
	
	/**
	 * Computes the statistics of a list of floating point data samples.
	 * The values computed are the sample mean and the variance.
	 * 
	 * @param data
	 * 			list of floating point sample points
	 * @return array of the values containing which contains:
	 * 			statistics[0] = mean of sample points;
	 * 			statistics[1] = variance of sample points.
	 */
	private static float[] getVarianceOfData(List<Float> data) {
		float mean = 0f;
		float var = 0f;
		for (int i = 0 ; i < data.size() ; ++i) {
			mean += data.get(i);
			var += data.get(i) * data.get(i);
		}
		mean /= data.size();
		var /= data.size() - mean * mean;
		float[] statistics = {mean , var};
		return statistics;
	}
}
