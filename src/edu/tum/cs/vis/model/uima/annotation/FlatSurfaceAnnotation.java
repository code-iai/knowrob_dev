/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v3.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors:
 *     Stefan Profanter - initial API and implementation, Year: 2012
 ******************************************************************************/
package edu.tum.cs.vis.model.uima.annotation;

import java.awt.Color;

import javax.vecmath.Vector3d;

import edu.tum.cs.vis.model.uima.analyzer.FlatSurfaceAnalyzer;
import edu.tum.cs.vis.model.uima.feature.Area;
import edu.tum.cs.vis.model.uima.feature.Dimension2D;
import edu.tum.cs.vis.model.uima.feature.NormalVector;
import edu.tum.cs.vis.model.util.Polygon;

/**
 * UIMA Annotation for a flat surface.
 * 
 * @author Stefan Profanter
 * @see FlatSurfaceAnalyzer
 * 
 */
public class FlatSurfaceAnnotation extends MeshAnnotation {

	/**
	 * auto generated
	 */
	private static final long	serialVersionUID	= 8767464136366533135L;

	/**
	 * Color for drawing this annotation
	 */
	// TODO change to static
	private Color				annotationColor		= new Color(255, 0, 0, 128);

	/**
	 * Total area of the flat surface
	 */
	private Area				area;
	/**
	 * Dimension2D of the surface in x and y coordinates. So normal vector will be (0,0,1)
	 */
	private Dimension2D			dimension;
	/**
	 * surface normal
	 */
	private NormalVector		normalVector;

	/**
	 * Default constructor
	 */
	// TODO remove
	public FlatSurfaceAnnotation() {
		annotationColor = new Color((float) Math.random(), (float) Math.random(),
				(float) Math.random());
	}

	@Override
	public Color getAnnotationColor() {
		return annotationColor;
	}

	/**
	 * Returns the area of the surface
	 * 
	 * @return the area
	 */
	public Area getArea() {
		return area;
	}

	/**
	 * Returns the dimension of the surface
	 * 
	 * @return the dimension in x and y coordinates. So normal vector will be (0,0,1)
	 */
	public Dimension2D getDimension() {
		return dimension;
	}

	/**
	 * Returns the surface normal
	 * 
	 * @return surface normal
	 */
	public Vector3d getNormalVector() {
		return normalVector;
	}

	/**
	 * Updates the members area, dimension, normalVector
	 */
	public void setFeatures() {
		if (mesh.getPolygons().size() == 0)
			return;
		normalVector = new NormalVector(mesh.getPolygons().get(0).getNormalVector());

		float ar = 0;

		for (Polygon p : mesh.getPolygons())
			ar += p.getArea();
		area = new Area();
		area.setSquareMM(ar);

		dimension = new Dimension2D();

		// TODO calculate dimension
	}
}
