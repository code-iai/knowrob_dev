package edu.tum.cs.vis.model.util;

import javax.vecmath.Point3f;

import processing.core.PApplet;

/**
 * A simple line of a model
 * 
 * @author Stefan Profanter
 *
 */
public class Line extends DrawObject {

	/**
	 * Color of the line
	 */
	public Appearance appearance;

	/**
	 * Default constructor
	 */
	public Line() {
		position = new Point3f[2];
	}

	/**
	 * Draw the line onto the applet.
	 * 
	 * @param applet
	 *            Applet to draw on
	 */
	public void draw(PApplet applet, int overrideColor) {
		setColor(applet, overrideColor);
		applet.line(position[0].x, position[0].y, position[0].z, position[1].x,
				position[1].y, position[1].z);
	}
}
