/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Stefan Profanter - initial API and implementation, Year: 2012
 ******************************************************************************/
package edu.tum.cs.vis.model;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.Callable;

import javax.vecmath.Point3f;
import javax.vecmath.Vector3f;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import processing.core.PGraphics;
import edu.tum.cs.vis.model.util.Curvature;
import edu.tum.cs.vis.model.util.Group;
import edu.tum.cs.vis.model.util.Line;
import edu.tum.cs.vis.model.util.ThreadPool;
import edu.tum.cs.vis.model.util.Triangle;
import edu.tum.cs.vis.model.util.Vertex;

/**
 * @author Stefan Profanter
 * 
 */
public class Model {

	static void calculateVertexNormalsForTriangle(Triangle t) {
		Vertex p0 = t.getPosition()[0];
		Vertex p1 = t.getPosition()[1];
		Vertex p2 = t.getPosition()[2];

		Vector3f a = new Vector3f(p0);
		a.sub(p1);
		Vector3f b = new Vector3f(p1);
		b.sub(p2);
		Vector3f c = new Vector3f(p2);
		c.sub(p0);
		float l2a = a.lengthSquared(), l2b = b.lengthSquared(), l2c = c.lengthSquared();
		if (l2a == 0 || l2b == 0 || l2c == 0)
			return;

		Vector3f facenormal = t.getNormalVector();

		Vector3f normalP0 = (Vector3f) facenormal.clone();
		normalP0.scale(1.0f / (l2a * l2c));
		synchronized (p0.getNormalVector()) {
			p0.getNormalVector().add(normalP0);
		}

		Vector3f normalP1 = (Vector3f) facenormal.clone();
		normalP1.scale(1.0f / (l2b * l2a));
		synchronized (p1.getNormalVector()) {
			p1.getNormalVector().add(normalP1);
		}

		Vector3f normalP2 = (Vector3f) facenormal.clone();
		normalP2.scale(1.0f / (l2c * l2b));
		synchronized (p2.getNormalVector()) {
			p2.getNormalVector().add(normalP2);
		}
	}

	static void calculateVoronoiAreaForTriangle(Triangle t) {
		// Edges
		Vector3f e[] = new Vector3f[3];

		for (int i = 0; i < 3; i++) {
			e[i] = new Vector3f(t.getPosition()[(i + 2) % 3]);
			e[i].sub(t.getPosition()[(i + 1) % 3]);
		}

		// Compute corner weights
		Vector3f tmp = new Vector3f();
		tmp.cross(e[0], e[1]);
		float area = 0.5f * tmp.length();
		float l2[] = { e[0].lengthSquared(), e[1].lengthSquared(), e[2].lengthSquared() };
		float ew[] = { l2[0] * (l2[1] + l2[2] - l2[0]), l2[1] * (l2[2] + l2[0] - l2[1]),
				l2[2] * (l2[0] + l2[1] - l2[2]) };
		t.setCornerarea(new Vector3f());
		if (ew[0] <= 0.0f) {
			t.getCornerarea().y = -0.25f * l2[2] * area / (e[0].dot(e[2]));
			t.getCornerarea().z = -0.25f * l2[1] * area / (e[0].dot(e[1]));
			t.getCornerarea().x = area - t.getCornerarea().y - t.getCornerarea().z;
		} else if (ew[1] <= 0.0f) {
			t.getCornerarea().z = -0.25f * l2[0] * area / (e[1].dot(e[0]));
			t.getCornerarea().x = -0.25f * l2[2] * area / (e[1].dot(e[2]));
			t.getCornerarea().y = area - t.getCornerarea().z - t.getCornerarea().x;
		} else if (ew[2] <= 0.0f) {
			t.getCornerarea().x = -0.25f * l2[1] * area / (e[2].dot(e[1]));
			t.getCornerarea().y = -0.25f * l2[0] * area / (e[2].dot(e[0]));
			t.getCornerarea().z = area - t.getCornerarea().x - t.getCornerarea().y;
		} else {
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
	 * Rotate a coordinate system to be perpendicular to the given normal
	 */
	static void rot_coord_sys(final Vector3f old_u, final Vector3f old_v, final Vector3f new_norm,
			Vector3f new_u, Vector3f new_v) {
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

	private Group					group;

	private final ArrayList<Vertex>	vertices	= new ArrayList<Vertex>();

	final ArrayList<Triangle>		triangles	= new ArrayList<Triangle>();

	private final ArrayList<Line>	lines		= new ArrayList<Line>();

	/**
	 * @param line
	 */
	private void addLine(Line line) {
		for (int i = 0; i < line.getPosition().length; i++) {
			line.getPosition()[i] = checkOrAddVertex(line.getPosition()[i]);
		}
		lines.add(line);
	}

	/**
	 * @param tri
	 */
	private void addTriangle(Triangle tri) {
		for (int i = 0; i < tri.getPosition().length; i++) {
			tri.getPosition()[i] = checkOrAddVertex(tri.getPosition()[i]);
		}
		triangles.add(tri);

	}

	private void calculateCurvature() {

		for (int i = 0; i < vertices.size(); i++) {
			vertices.get(i).setCurvature(new Curvature());
		}

		// Set up an initial coordinate system per vertex
		for (int i = 0; i < triangles.size(); i++) {
			Triangle t = triangles.get(i);
			for (int j = 0; j < 3; j++) {
				Vertex v = t.getPosition()[j];
				v.getCurvature().setPrincipleDirectionMin(
						new Vector3f(t.getPosition()[(j + 1) % 3]));
				v.getCurvature().getPrincipleDirectionMin().sub(v);
			}
		}
		for (int i = 0; i < vertices.size(); i++) {

			Vector3f tmp = new Vector3f();
			tmp.cross(vertices.get(i).getCurvature().getPrincipleDirectionMin(), vertices.get(i)
					.getNormalVector());
			tmp.normalize();
			vertices.get(i).getCurvature().setPrincipleDirectionMin(tmp);

			tmp = new Vector3f();
			tmp.cross(vertices.get(i).getNormalVector(), vertices.get(i).getCurvature()
					.getPrincipleDirectionMin());
			vertices.get(i).getCurvature().setPrincipleDirectionMax(tmp);
		}

		// Compute curvature per-face

		for (int i = 0; i < triangles.size(); i++) {

			Triangle tri = triangles.get(i);
			// Edges
			Vector3f e[] = new Vector3f[3];

			for (int j = 0; j < 3; j++) {
				e[j] = new Vector3f(tri.getPosition()[(j + 2) % 3]);
				e[j].sub(tri.getPosition()[(j + 1) % 3]);
			}

			// N-T-B coordinate system per face
			Vector3f t = new Vector3f(e[0]);
			t.normalize();
			Vector3f n = new Vector3f();
			n.cross(e[0], e[1]);
			Vector3f b = new Vector3f();
			b.cross(n, t);
			b.normalize();

			// Estimate curvature based on variation of normals
			// along edges
			float m[] = { 0, 0, 0 };
			double w[][] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
			for (int j = 0; j < 3; j++) {

				float u = e[j].dot(t);
				float v = e[j].dot(b);
				w[0][0] += u * u;
				w[0][1] += u * v;
				// w[1][1] += v*v + u*u;
				// w[1][2] += u*v;
				w[2][2] += v * v;
				Vector3f dn = new Vector3f(tri.getPosition()[(j + 2) % 3].getNormalVector());
				dn.sub(tri.getPosition()[(j + 1) % 3].getNormalVector());
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
				System.out.println("Singular matrix");
				continue;
			}

			RealVector constants = new ArrayRealVector(new double[] { m[0], m[1], m[2] }, false);
			RealVector solution = solver.solve(constants);

			m[0] = (float) solution.getEntry(0);
			m[1] = (float) solution.getEntry(1);
			m[2] = (float) solution.getEntry(2);

			// Push it back out to the vertices
			for (int j = 0; j < 3; j++) {
				Vertex vj = tri.getPosition()[j];

				float c1, c12, c2;
				float ret[] = proj_curv(t, b, m[0], m[1], m[2], vj.getCurvature()
						.getPrincipleDirectionMin(), vj.getCurvature().getPrincipleDirectionMax());
				c1 = ret[0];
				c12 = ret[1];
				c2 = ret[2];

				float wt;
				if (j == 0)
					wt = tri.getCornerarea().x / vj.getPointarea();
				else if (j == 1)
					wt = tri.getCornerarea().y / vj.getPointarea();
				else
					wt = tri.getCornerarea().z / vj.getPointarea();

				vj.getCurvature().setCurvatureMin(vj.getCurvature().getCurvatureMin() + wt * c1);
				if (vj.getCurvature().getCurvatureMin() == Float.NaN) {
					System.out.println("Nan");
				}
				vj.getCurvature().setCurvatureMinMax(
						vj.getCurvature().getCurvatureMinMax() + wt * c12);
				vj.getCurvature().setCurvatureMax(vj.getCurvature().getCurvatureMax() + wt * c2);
			}
		}

		// Compute principal directions and curvatures at each vertex
		for (int i = 0; i < vertices.size(); i++) {
			Curvature c = vertices.get(i).getCurvature();
			Vector3f pdirRet[] = new Vector3f[2];
			float kRet[] = new float[2];
			diagonalize_curv(c.getPrincipleDirectionMin(), c.getPrincipleDirectionMax(),
					c.getCurvatureMin(), c.getCurvatureMinMax(), c.getCurvatureMax(),
					vertices.get(i).getNormalVector(), pdirRet, kRet);
			c.setPrincipleDirectionMin(pdirRet[0]);
			c.setPrincipleDirectionMax(pdirRet[1]);
			c.setCurvatureMin(kRet[0]);
			c.setCurvatureMax(kRet[1]);
		}
	}

	/**
	 * Ported from trimesh2 (Szymon Rusinkiewicz Princeton University)
	 * 
	 * uses average of per-face normals, weighted according to: Max, N.
	 * "Weights for Computing Vertex Normals from Facet Normals," Journal of Graphics Tools, Vol. 4,
	 * No. 2, 1999.
	 */
	private void calculateVertexNormals() {
		// Compute from faces

		List<Callable<Void>> threads = new LinkedList<Callable<Void>>();

		final int interval = 500;

		for (int start = 0; start < triangles.size(); start += interval) {
			final int st = start;
			threads.add(new Callable<Void>() {

				@Override
				public Void call() throws Exception {
					int end = Math.min(st + interval, triangles.size());
					for (int i = st; i < end; i++) {
						calculateVertexNormalsForTriangle(triangles.get(i));
					}
					return null;
				}

			});
		};

		ThreadPool.executeInPool(threads);

		for (Vertex v : vertices)
			v.getNormalVector().normalize();

	}

	/**
	 * Ported from trimesh2 (Szymon Rusinkiewicz Princeton University)
	 * 
	 * Compute the area "belonging" to each vertex or each corner of a triangle (defined as Voronoi
	 * area restricted to the 1-ring of a vertex, or to the triangle).
	 */
	private void calculateVoronoiArea() {
		// Compute from faces

		List<Callable<Void>> threads = new LinkedList<Callable<Void>>();

		final int interval = 500;

		for (Vertex v : vertices)
			v.setPointarea(0);

		for (int start = 0; start < triangles.size(); start++) {
			calculateVoronoiAreaForTriangle(triangles.get(start));
		};

		/*for (int start = 0; start < triangles.size(); start += interval) {
			final int st = start;
			threads.add(new Callable<Void>() {

				@Override
				public Void call() throws Exception {
					int end = Math.min(st + interval, triangles.size());
					for (int i = st; i < end; i++) {
						calculateVoronoiAreaForTriangle(triangles.get(i));
					}
					return null;
				}

			});
		};

		ThreadPool.executeInPool(threads);*/
	}

	private Vertex checkOrAddVertex(Vertex v) {

		for (Vertex p : vertices) {
			if (p.equals(v))
				return p;
		}
		vertices.add(v);
		return v;
	}

	void colorbycurv(float scale, float smooth) {
		/*	float smoothsigma = smooth;
			if (smoothsigma > 0.0f) {
				smoothsigma *= feature_size();
				diffuse_curv(mesh, smoothsigma);
			}*/
		float cscale = 10.0f * scale * typical_scale();
		System.out.println("Using scale = " + cscale);
		cscale = cscale * cscale;

		int nv = vertices.size();
		for (int i = 0; i < nv; i++) {
			Curvature c = vertices.get(i).getCurvature();
			float H = 0.5f * (c.getCurvatureMin() + c.getCurvatureMax());
			float K = c.getCurvatureMin() * c.getCurvatureMax();
			float h = (float) (4.0f / 3.0f * Math
					.abs(Math.atan2(H * H - K, H * H * Math.signum(H))));
			float s = (float) ((2 / Math.PI) * Math.atan((2.0f * H * H - K) * cscale));
			vertices.get(i).color = Color.getHSBColor(h, s, 1.0f);
		}
	}

	/**
	 * Given a curvature tensor, find principal directions and curvatures. Makes sure that pdir1 and
	 * pdir2 are perpendicular to normal
	 * 
	 * returns pdir1
	 */
	void diagonalize_curv(final Vector3f old_u, final Vector3f old_v, float ku, float kuv,
			float kv, final Vector3f new_norm, Vector3f pdir[], float k[]) {
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
			r_old_u.sub(r_old_v);
			pdir[0] = r_old_u;
		}
		pdir[1] = new Vector3f();
		pdir[1].cross(new_norm, pdir[0]);
	}

	/**
	 * @param g
	 * @param overrideColor
	 */
	public void draw(PGraphics g, Color overrideColor) {
		if (group != null)
			group.draw(g, overrideColor);
	}

	/**
	 * A characteristic "feature size" for the mesh. Computed as an approximation to the median edge
	 * length
	 * 
	 * @return
	 */
	float feature_size() {

		int nf = triangles.size();
		int nsamp = Math.min(nf / 2, 333);

		float samples[] = new float[nsamp * 3];
		int idx = 0;

		for (int i = 0; i < nsamp; i++) {
			// Quick 'n dirty portable random number generator
			int ind = (int) (Math.random() * nf);

			Point3f p0 = triangles.get(ind).getPosition()[0];
			Point3f p1 = triangles.get(ind).getPosition()[1];
			Point3f p2 = triangles.get(ind).getPosition()[2];

			Vector3f f1 = new Vector3f(p0);
			f1.sub(p1);
			samples[idx++] = f1.lengthSquared();

			Vector3f f2 = new Vector3f(p1);
			f2.sub(p2);
			samples[idx++] = f2.lengthSquared();

			Vector3f f3 = new Vector3f(p2);
			f3.sub(p0);
			samples[idx++] = f3.lengthSquared();
		}
		Arrays.sort(samples);
		return (float) Math.sqrt(samples[samples.length / 2]);
	}

	/**
	 * @return the group
	 */
	public Group getGroup() {
		return group;
	}

	/**
	 * @return the lines
	 */
	public ArrayList<Line> getLines() {
		return lines;
	}

	/**
	 * @return the triangles
	 */
	public ArrayList<Triangle> getTriangles() {
		return triangles;
	}

	/**
	 * @return the vertices
	 */
	public ArrayList<Vertex> getVertices() {
		return vertices;
	}

	/**
	 * 
	 */
	public void modelChanged() {
		vertices.clear();
		triangles.clear();
		lines.clear();
		processGroup(group);
		calculateVertexNormals();
		calculateVoronoiArea();
		calculateCurvature();
		colorbycurv(1, 0);
	}

	private void processGroup(Group g) {
		for (Line l : g.getMesh().getLines())
			addLine(l);
		for (Triangle t : g.getMesh().getTriangles())
			addTriangle(t);
		for (Group c : g.getChildren())
			processGroup(c);
	}

	/**
	 * Reproject a curvature tensor from the basis spanned by old_u and old_v (which are assumed to
	 * be unit-length and perpendicular) to the new_u, new_v basis. returns [new_ku, new_kuv,
	 * new_kv]
	 * 
	 * @param old_u
	 * @param old_v
	 * @param old_ku
	 * @param old_kuv
	 * @param old_kv
	 * @param new_u
	 * @param new_v
	 * @return
	 */
	private float[] proj_curv(final Vector3f old_u, final Vector3f old_v, float old_ku,
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
	 * @param group
	 *            the group to set
	 */
	public void setGroup(Group group) {
		this.group = group;
		modelChanged();
	}

	/**
	 * Compute a "typical scale" for the mesh: computed as 1% of the reciprocal of the 10-th
	 * percentile curvature
	 * 
	 * @return
	 */
	float typical_scale() {
		float frac = 0.1f;
		float mult = 0.01f;

		int nv = vertices.size();
		int nsamp = Math.min(nv, 500);
		if (nsamp == 0)
			return 0;

		float samples[] = new float[nsamp * 2];
		int idx = 0;
		for (int i = 0; i < nsamp; i++) {
			int ind = (int) (Math.random() * nv);

			samples[idx++] = Math.abs(vertices.get(ind).getCurvature().getCurvatureMin());
			samples[idx++] = Math.abs(vertices.get(ind).getCurvature().getCurvatureMax());
		}

		int which = (int) (frac * samples.length);
		Arrays.sort(samples);

		float f = 0;
		if (samples[which] == 0.0f) {
			// mesh->need_bsphere();
			// f = mult * mesh->bsphere.r;
			System.out.println("NOOOOOOOO Bsphere!!");
		} else {
			f = mult / samples[which];
		}
		System.out.println("Typical scale = " + f);
		return f;
	}
}
