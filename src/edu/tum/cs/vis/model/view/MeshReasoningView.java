/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Stefan Profanter - initial API and implementation, Year: 2012
 ******************************************************************************/
package edu.tum.cs.vis.model.view;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.util.ArrayList;

import javax.swing.event.MouseInputListener;
import javax.vecmath.Point3f;
import javax.vecmath.Vector3f;

import peasy.PeasyCam;
import edu.tum.cs.uima.Annotation;
import edu.tum.cs.vis.model.uima.annotation.MeshAnnotation;
import edu.tum.cs.vis.model.uima.cas.MeshCas;
import edu.tum.cs.vis.model.util.Polygon;

/**
 * Viewing applet for showing the results of the reasoning process.
 * 
 * @author Stefan Profanter
 * 
 */
public final class MeshReasoningView extends PAppletSelection implements MouseInputListener {

	/**
	 * 
	 */
	private static final long				serialVersionUID	= 984696039698156574L;

	/**
	 * Cam for manipulating the view
	 */
	PeasyCam								cam;

	/**
	 * Background color
	 */
	private final Color						bgcolor				= new Color(10, 10, 10);

	/**
	 * Start point of mouse click ray
	 */
	private Point3f							rayStart			= new Point3f();
	/**
	 * End point of mouse click ray
	 */
	private Point3f							rayEnd				= new Point3f(1, 1, 1);
	/**
	 * List of all CASes which were manipulated with AnalysisEngines.
	 */
	private ArrayList<MeshCas>				casList				= new ArrayList<MeshCas>();

	/**
	 * List of selected polygons (polygons which intersect with mouse ray)
	 */
	private final ArrayList<Polygon>		selectedPolygons	= new ArrayList<Polygon>();
	/**
	 * List of selected annotations (annotations which contain one of selectedPolygons)
	 */
	private final ArrayList<MeshAnnotation>	selectedAnnotations	= new ArrayList<MeshAnnotation>();

	/**
	 * The controller for this view
	 */
	private MeshReasoningViewControl		control				= null;

	@Override
	public void draw() {

		scale(20);
		background(bgcolor.getRed(), bgcolor.getGreen(), bgcolor.getBlue());
		// draw axis
		noFill();
		strokeWeight(1);
		stroke(125, 0, 0);
		strokeWeight(1);
		line(0, 0, 0, width, 0, 0);
		stroke(0, 125, 0);
		line(0, 0, 0, 0, 0, width);
		stroke(0, 0, 125);
		line(0, 0, 0, 0, height, 0);

		/*strokeWeight(5);
		stroke(255, 255, 0);

		synchronized (rayStart) {
			line(rayStart.x, rayStart.y, rayStart.z, rayEnd.x, rayEnd.y, rayEnd.z);
		}
		pushMatrix();

		translate(intersect.x, intersect.y, intersect.z);
		noStroke();
		scale(0.1f);
		fill(0, 0, 255);
		sphere(1);
		popMatrix();

		fill(127);*/

		Vector3f camPos = new Vector3f(cam.getPosition());

		lights();
		pointLight(80f, 80f, 100f, camPos.x, camPos.y, camPos.z);

		// Must be called AFTER all scale, transform, rotate, ... calls
		captureViewMatrix();

		getSelectionGraphics().setDrawWithTransparency(selectedAnnotations.size() > 0);
		for (MeshCas c : casList) {
			c.draw(g);
			// c.getGroup().draw(g, null);
		}
		getSelectionGraphics().setDrawWithTransparency(false);

		synchronized (selectedAnnotations) {

			for (MeshAnnotation ma : selectedAnnotations) {
				ma.getMesh().drawPolygons(g, new Color(255, 125, 0));
			}
		}
	}

	/**
	 * Get list of all CAS objects (see UIMA Framework)
	 * 
	 * @return ArrayList of MeshCas
	 */
	public ArrayList<MeshCas> getCasList() {
		return casList;
	}

	@Override
	public void mouseClicked(MouseEvent e) {
		if (e.getButton() == 1) {
			calculatePickPoints(e.getX(), e.getY());
			if (!isMouseRayValid()) {
				System.out.println("Mouse ray not valid!!!");
				return;
			}
			rayStart = getMouseRayStart();
			rayEnd = getMouseRayEnd();

			boolean found = false;
			// Check if clicked on one of previous selected polygons
			for (Polygon p : selectedPolygons) {
				if (p.intersectsRay(rayEnd, rayStart, null)) {
					selectedPolygons.clear();
					selectedPolygons.add(p);
					found = true;
					break;
				}
			}

			if (!found) {
				// It is new selection
				selectedPolygons.clear();
				for (MeshCas c : casList) {
					c.getGroup().getIntersectedPolygons(rayEnd, rayStart, selectedPolygons);
				}
			}

			// Check if one of selected polygons is in already selected annotation
			ArrayList<Polygon> newSelected = new ArrayList<Polygon>();
			for (Polygon p : selectedPolygons) {
				synchronized (selectedAnnotations) {
					for (MeshAnnotation ma : selectedAnnotations)
						if (ma.meshContainsPolygon(p)) {
							newSelected.add(p);
						}
				}
			}
			if (newSelected.size() > 0) {
				// Currently selected was in one or more of the selected annotations, so select out
				// of current annotations
				selectedPolygons.clear();
				selectedPolygons.addAll(newSelected);
			}

			selectedPolygonsChanged();
		}

	}

	/**
	 * Called to update list of selectedAnnotations after list of selectedPolygons has changed
	 */
	private void selectedPolygonsChanged() {
		synchronized (selectedAnnotations) {
			selectedAnnotations.clear();
		}
		for (MeshCas c : casList) {
			for (Annotation a : c.getAnnotations()) {
				if (!(a instanceof MeshAnnotation))
					continue;
				MeshAnnotation ma = (MeshAnnotation) a;
				if (!ma.isDrawAnnotation())
					continue; // Skip not visible annotations
				if (selectedAnnotations.contains(ma))
					continue;
				for (Polygon p : selectedPolygons)
					if (ma.meshContainsPolygon(p)) {
						synchronized (selectedAnnotations) {
							selectedAnnotations.add(ma);
						}
						break;
					}
			}
		}
		control.showSelectedAnnotation(selectedAnnotations);
	}

	/**
	 * Set list of all CAS objects (see UIMA Framework)
	 * 
	 * @param casList
	 *            list to set
	 */
	public void setCasList(ArrayList<MeshCas> casList) {
		this.casList = casList;
	}

	/**
	 * @param control
	 *            the control to set
	 */
	public void setControl(MeshReasoningViewControl control) {
		this.control = control;
	}

	@Override
	public void setup() {
		size(1000, 1000, "edu.tum.cs.vis.model.view.PAppletSelectionGraphics");

		frameRate(30);
		cam = new PeasyCam(this, 0, 0, 0, 10);
		cam.setMinimumDistance(0.01);
		cam.setMaximumDistance(500);

		cam.setRightDragHandler(cam.getPanDragHandler());

		cam.setDistance(50);

		cam.rotateX((float) Math.PI / 2f);
		cam.rotateZ((float) Math.PI);

		captureViewMatrix();

		perspective();

		draw();
	}
}
