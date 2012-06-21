package edu.tum.cs.vis.model.uima.annotation;

import java.awt.Color;

import javax.vecmath.Vector3d;

import edu.tum.cs.vis.model.uima.analyzer.FlatSurfaceAnalyzer;
import edu.tum.cs.vis.model.uima.feature.Area;
import edu.tum.cs.vis.model.uima.feature.Dimension2D;
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
	 * Color for drawing this annotation
	 */
	// TODO change to static
	private Color		annotationColor	= new Color(255, 0, 0, 128);

	/**
	 * Total area of the flat surface
	 */
	private Area		area;
	/**
	 * Dimension2D of the surface in x and y coordinates. So normal vector will be (0,0,1)
	 */
	private Dimension2D	dimension;
	/**
	 * surface normal
	 */
	private Vector3d	normalVector;

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

	public void setFeatures() {
		if (mesh.getPolygons().size() == 0)
			return;
		normalVector = mesh.getPolygons().get(0).getNormalVector();

		float ar = 0;

		for (Polygon p : mesh.getPolygons())
			ar += p.getArea();
		area = new Area();
		area.setSquareMM(ar);

		// TODO calculate dimension
	}
}
