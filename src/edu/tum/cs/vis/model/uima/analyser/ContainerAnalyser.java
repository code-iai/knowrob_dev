/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Stefan Profanter - initial API and implementation, Year: 2012
 ******************************************************************************/
package edu.tum.cs.vis.model.uima.analyser;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;

import javax.vecmath.Point3f;
import javax.vecmath.Vector3f;

import org.apache.log4j.Logger;

import edu.tum.cs.vis.model.uima.annotation.ContainerAnnotation;
import edu.tum.cs.vis.model.uima.annotation.MeshAnnotation;
import edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation;
import edu.tum.cs.vis.model.uima.annotation.primitive.ConeAnnotation;
import edu.tum.cs.vis.model.uima.annotation.primitive.PlaneAnnotation;
import edu.tum.cs.vis.model.uima.cas.MeshCas;
import edu.tum.cs.vis.model.util.ThreadPool;
import edu.tum.cs.vis.model.util.Triangle;

/**
 * @author Stefan Profanter
 * 
 */
public class ContainerAnalyser extends MeshAnalyser {

	/**
	 * Log4J Logger
	 */
	private static Logger	logger	= Logger.getLogger(ContainerAnalyser.class);

	static void checkConcavAnnotationCone(MeshCas cas, ConeAnnotation an,
			HashSet<PrimitiveAnnotation> annotations) {

		float areaOfCapSmall = (float) (Math.pow(an.getRadiusSmall(), 2) * Math.PI);
		float areaOfCapLarge = (float) (Math.pow(an.getRadiusLarge(), 2) * Math.PI);
		float heightHalf = an.getDirection().length();

		Point3f p1 = new Point3f(an.getCentroid());
		p1.add(an.getDirection());

		Point3f p2 = new Point3f(an.getCentroid());
		// p2.sub(d);

		boolean hasTopCap = false;
		boolean hasBottomCap = false;

		float prevPropAbs = Float.MAX_VALUE;

		PlaneAnnotation capAnnotation = null;

		for (MeshAnnotation<PlaneAnnotation> a : annotations) {
			if (a instanceof PlaneAnnotation) {
				PlaneAnnotation pa = (PlaneAnnotation) a;

				Point3f intersect = new Point3f();
				boolean isIntersect = false;
				for (Triangle t : pa.getMesh().getTriangles()) {
					if (t.intersectsRay(p1, p2, intersect)) {
						isIntersect = true;
						break;
					}
				}
				if (!isIntersect) {
					continue;
				}

				// intersection point found. Check distance:
				Vector3f midToIntersect = new Vector3f(intersect);
				midToIntersect.sub(an.getCentroid());
				float distToIntersect = midToIntersect.length();

				float proportion = distToIntersect / heightHalf;

				float propAbs = Math.abs(1f - proportion);

				// System.out.println("abs: " + propAbs + " prev: " + prevPropAbs);

				if (propAbs > 0.3 || propAbs > prevPropAbs)
					continue; // plane isn't at the top or bottom of the cylinder or is not as near
								// as previous found

				if (midToIntersect.dot(an.getDirection()) > 0) {
					// use top
					if (pa.getArea() / areaOfCapSmall > 0.7) {
						hasTopCap = true;
						capAnnotation = pa;
						prevPropAbs = propAbs;
					}
				} else {
					if (pa.getArea() / areaOfCapLarge > 0.7) {
						hasBottomCap = true;
						capAnnotation = pa;
						prevPropAbs = propAbs;
					}
				}

				if (hasBottomCap && hasTopCap)
					return; // isn't container. Is a closed box.
			}
		}
		if (hasBottomCap || hasTopCap) {
			ContainerAnnotation ca = new ContainerAnnotation(cas.getModel());
			ca.getMesh().getTriangles().addAll(an.getMesh().getTriangles());
			ca.getMesh().getTriangles().addAll(capAnnotation.getMesh().getTriangles());

			ca.setDirection((Vector3f) an.getDirection().clone());

			if (hasTopCap) {
				ca.getDirection().scale(-1);
			}

			ca.setVolume(an.getVolume());

			cas.addAnnotation(ca);
		}
	}

	/**
	 * Number of triangles already elaborated/processed. Used for indicating current process
	 */
	final AtomicInteger	itemsElaborated		= new AtomicInteger(0);

	private int			maxItemsToElaborate	= 0;

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
		return "Container";
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.analyser.MeshAnalyser#processStart(edu.tum.cs.vis.model.uima.cas.MeshCas)
	 */
	@Override
	public void processStart(final MeshCas cas) {
		List<Callable<Void>> threads = new LinkedList<Callable<Void>>();

		final HashSet<PrimitiveAnnotation> annotations = cas
				.findAnnotations(PrimitiveAnnotation.class);

		for (final MeshAnnotation a : annotations) {
			if (a instanceof ConeAnnotation && ((ConeAnnotation) a).isConcav()) {
				maxItemsToElaborate++;
				threads.add(new Callable<Void>() {

					@Override
					public Void call() throws Exception {
						checkConcavAnnotationCone(cas, (ConeAnnotation) a, annotations);

						itemsElaborated.incrementAndGet();
						return null;
					}
				});
			} // TODO Add container finding for half-spheres or rounded boxes or simple boxes
		}

		ThreadPool.executeInPool(threads);
		// increment, because if no containers, 1 indicates that
		// it has finished
		maxItemsToElaborate = 1;
		itemsElaborated.incrementAndGet();

	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.analyser.MeshAnalyser#updateProgress()
	 */
	@Override
	public void updateProgress() {
		if (maxItemsToElaborate > 0) {
			setProgress(itemsElaborated.get() / (float) (maxItemsToElaborate) * 100.0f);
		}

	}

}
