/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: 
 * 				Stefan Profanter - initial API and implementation, Year: 2012
 * 				Andrei Stoica - refactored implementation during Google Summer of Code 2014
 ******************************************************************************/
package org.knowrob.vis.model.view;

import java.awt.Color;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import javax.swing.event.MouseInputListener;
import javax.vecmath.Point3f;
import javax.vecmath.Vector3f;

import org.apache.log4j.Logger;

import peasy.PeasyCam;
import processing.core.PConstants;
import processing.core.PGraphics;
import org.knowrob.vis.model.uima.annotation.ComplexHandleAnnotation;
import org.knowrob.vis.model.uima.annotation.DrawableAnnotation;
import org.knowrob.vis.model.uima.annotation.PrimitiveAnnotation;
import org.knowrob.vis.model.uima.cas.MeshCas;
import org.knowrob.vis.model.util.Curvature;
import org.knowrob.vis.model.util.DrawSettings;
import org.knowrob.vis.model.util.DrawType;
import org.knowrob.vis.model.util.Edge;
import org.knowrob.vis.model.util.IntersectedTriangle;
import org.knowrob.vis.model.util.Region;
import org.knowrob.vis.model.util.Triangle;
import org.knowrob.vis.model.util.Vertex;
import org.knowrob.vis.tools.ImageGenerator.ImageGeneratorSettings;
import org.knowrob.vis.uima.Annotation;

/**
 * Class which implements the main GUI applet for showing the results of the reasoning process.
 * It implements the functionality that works in the backend to provide the model visualization
 * throughout the mesh reasoning process. It contains different view perspectives meant
 * for both debugging and showing the results of the processing. The possibility to save the
 * visualization at any time is also implemented.
 * 
 * @author Stefan Profanter
 * @author Andrei Stoica
 */
public final class MeshReasoningView extends PAppletSelection implements MouseInputListener,
		KeyListener {

	/**
	 * auto generated
	 */
	private static final long	serialVersionUID	= 984696039698156574L;

	/**
	 * Cam for manipulating the view
	 */
	private static PeasyCam							cam;
	
	/**
	 * Background color dark
	 */
	private final Color								bgcolorDark			= new Color(10, 10, 10);

	/**
	 * Background color white
	 */
	private final Color								bgcolorWhite		= new Color(255, 255, 255);

	/**
	 * White background flag
	 */
	private boolean									backgroundWhite		= false;
	
	/**
	 * Normals of vertices drawing flag
	 */
	private boolean									drawVertexNormals	= false;
	
	/**
	 * Normals of triangles drawing flag
	 */
	private boolean									drawTriangleNormals = false;
	
	/**
	 * Minimum normalized curvature direction of each vertex drawing flag
	 */
	private boolean									drawVertexCurvatureMin = false;
	
	/**
	 * Maximum normalized curvature direction of each vertex drawing flag
	 */
	private boolean 								drawVertexCurvatureMax = false;
	
	/**
	 * Minimum and maximum scaled curvature directions of each vertex drawing flag.
	 * The scaling is done with the associated curvature values.
	 */
	private boolean									drawVertexCurvature	= false;

	/**
	 * Voronoi area at each vertex drawing flag
	 */
	private boolean									drawVoronoiArea		= false;
	
	/**
	 * Sharp edges drawing flag
	 */
	private boolean									drawSharpEdges 		= false;
	
	/**
	 * Region edges drawing flag
	 */
	private boolean									drawRegionEdges		= false;
	
	/**
	 * Select only nearest triangle flag
	 */
	private boolean									selectNearestOnly	= true;
	
	/**
	 * Bounding box drawing flag
	 */
	private boolean									drawBoundingBox		= false;

	/**
	 * Start point of mouse click ray
	 */
	private Point3f									rayStart			= new Point3f();

	/**
	 * End point of mouse click ray
	 */
	private Point3f									rayEnd				= new Point3f(1, 1, 1);

	/**
	 * List of all CASes which were manipulated with AnalysisEngines.
	 */
	private ArrayList<MeshCas>						casList				= new ArrayList<MeshCas>();

	/**
	 * Current scale factor of the model in the MeshCas
	 */
	private float									modelScale			= 1f;

	/**
	 * User can manually scale model by pressing '+' or '-'
	 */
	private float									userScale			= 30f;

	/**
	 * Path where to save image if user clicks on "save image". Will be evaluated in draw method
	 */
	private String									imageSavePath		= null;

	/**
	 * List of selected triangles (triangles which intersect with mouse ray)
	 */
	private final ArrayList<IntersectedTriangle>	selectedTriangles	= new ArrayList<IntersectedTriangle>();

	/**
	 * List of selected annotations (annotations which contain one of selectedTriangles)
	 * */
	private final Map<DrawableAnnotation, Color>	selectedAnnotations	= new HashMap<DrawableAnnotation, Color>();

	/**
	 * If set to true, only the clicked triangle is selected. if set to false, the whole annotation
	 * is selected.
	 */
	private boolean									selectTrianglesOnly	= false;

	/**
	 * The controller for this view
	 */
	private MeshReasoningViewControl				control				= null;

	/**
	 * Draw settings for the view
	 */
	private final DrawSettings						drawSettings		= new DrawSettings();

	/**
	 * Image generator settings. If set to null, image generation is disabled.
	 */
	private ImageGeneratorSettings					imageGeneratorSettings;

	/**
	 * Adds annotation to selected annotations
	 * 
	 * @param a
	 *            annotation to add
	 */
	public void addSelectedAnnotation(DrawableAnnotation a) {
		addSelectedAnnotation(a, null);
	}

	/**
	 * Adds given annotation to selected annotations list which means the annotation is selected. You
	 * can additionally specify a color for the annotation.
	 * 
	 * @param a
	 *            annotation to select
	 * @param color
	 *            color of the selected annotation or null to use default
	 */
	public void addSelectedAnnotation(DrawableAnnotation a, Color color) {
		synchronized (selectedAnnotations) {
			selectedAnnotations.put(a, color);
		}
		control.showSelectedAnnotation(selectedAnnotations);
	}

	/**
	 * Removes / unselects all selected annotations
	 */
	public void clearSelectedAnnotations() {
		synchronized (selectedAnnotations) {
			selectedAnnotations.clear();
		}
		control.showSelectedAnnotation(selectedAnnotations);
	}

	@Override
	public void draw() {

		scale(modelScale * userScale);
		if (backgroundWhite
				|| (imageGeneratorSettings != null && imageGeneratorSettings.isWhiteBackground()))
			background(bgcolorWhite.getRed(), bgcolorWhite.getGreen(), bgcolorWhite.getBlue());
		else
			background(bgcolorDark.getRed(), bgcolorDark.getGreen(), bgcolorDark.getBlue());
		// draw axis

		if (imageGeneratorSettings == null || (imageGeneratorSettings.isDrawAxis())) {
			noFill();
			strokeWeight(1);
			stroke(125, 0, 0);
			strokeWeight(1);
			line(0, 0, 0, 1, 0, 0);
			stroke(0, 125, 0);
			line(0, 0, 0, 0, 1, 0);
			stroke(0, 0, 125);
			line(0, 0, 0, 0, 0, 1);
		}

		Vector3f camPos = new Vector3f(cam.getPosition());

		lights();
		pointLight(80f, 80f, 100f, camPos.x, camPos.y, camPos.z);

		// Must be called AFTER all scale, transform, rotate, ... calls
		captureViewMatrix();

		getSelectionGraphics().setDrawWithTransparency(
				(selectTrianglesOnly && selectedTriangles.size() > 0)
						|| (selectedAnnotations.size() > 0) || drawSharpEdges || drawRegionEdges);

		for (MeshCas c : casList) {
			if (c.getModel() == null)
				return;

			c.draw(g, drawSettings);
			if (drawBoundingBox) {
				g.noFill();
				g.stroke(255, 125, 0);
				c.getModel().getGroup().drawBoundingBox(g, true);
			}

		}

		getSelectionGraphics().setDrawWithTransparency(false);

		if (selectTrianglesOnly) {

			synchronized (selectedTriangles) {
				DrawSettings tmpSet = drawSettings
						.getTemporaryOverride(selectedTriangles.size() > 1 ? new Color(255, 125, 0,
								200) : new Color(255, 50, 0, 200));
				for (IntersectedTriangle t : selectedTriangles) {
					t.t.draw(g, tmpSet);
				}

				if (selectedTriangles.size() == 1) {
					IntersectedTriangle t = selectedTriangles.get(0);

					tmpSet = drawSettings.getTemporaryOverride(new Color(0, 50, 255, 200));
					for (Triangle n : t.t.getNeighbors()) {

						n.draw(g, tmpSet);
					}
				}

			}
		} else {
			synchronized (selectedAnnotations) {

				Color defaultColor = selectedAnnotations.size() > 1 ? new Color(255, 125, 0, 200)
						: new Color(255, 50, 0, 200);
				boolean drawPrimitives = (selectedAnnotations.size() == 1 || (imageGeneratorSettings != null && imageGeneratorSettings
						.isAlwaysDrawSelectedPrimitives()));
				for (DrawableAnnotation ma : selectedAnnotations.keySet()) {
					Color anColor = selectedAnnotations.get(ma);
					DrawSettings tmpSet = drawSettings
							.getTemporaryOverride(anColor == null ? defaultColor : anColor);
					ma.draw(g, tmpSet);

					if (drawPrimitives) {

						Color primColor = anColor != null ? new Color(tmpSet.getOverrideColor()
								.getRed(), tmpSet.getOverrideColor().getGreen(), tmpSet
								.getOverrideColor().getBlue(), 255) : null;
						if (ma instanceof PrimitiveAnnotation) {
							@SuppressWarnings("rawtypes")
							PrimitiveAnnotation an = (PrimitiveAnnotation) ma;

							an.drawPrimitiveAnnotation(g, primColor);

						} else if (ma instanceof ComplexHandleAnnotation) {
							ComplexHandleAnnotation an = (ComplexHandleAnnotation) ma;

							an.drawPrimitiveAnnotation(g, primColor);

						}
					}
				}
			}
		}

		if (drawVertexNormals || drawTriangleNormals || drawVertexCurvatureMin || drawVertexCurvatureMax 
				|| drawVertexCurvature || drawVoronoiArea || drawSharpEdges || drawRegionEdges) {
			g.strokeWeight(2f);
			for (MeshCas c : casList) {
				synchronized (c.getModel().getVertices()) {
					for (Vertex v : c.getModel().getVertices()) {
						if (drawVertexNormals || drawVoronoiArea) {
							g.stroke(41, 120, 40);
							Vector3f n = (Vector3f) v.getNormalVector().clone();
							n.scale(0.075f);
							g.line(v.x, v.y, v.z, v.x + n.x, v.y + n.y, v.z + n.z);
							g.fill(35, 150, 140);
							g.noStroke();
							g.sphereDetail(20);
							if (drawVoronoiArea) {
								g.pushMatrix();
								g.translate(v.x + n.x, v.y + n.y, v.z + n.z);
								g.sphere(v.getPointarea());
								g.popMatrix();
							}
						}
						if (drawVertexCurvature || drawVertexCurvatureMin || drawVertexCurvatureMax) {
							Curvature curv = c.getCurvature(v);
							if (curv == null) {
								continue;
							}
							g.strokeWeight(2f);
							if (drawVertexCurvatureMin) {
								g.stroke(50, 50, 255);
								Vector3f min = (Vector3f) curv.getPrincipleDirectionMin().clone();
								min.normalize();
								min.scale(1 / 35f);
								g.line(v.x - min.x/2, v.y - min.y/2, v.z - min.z/2, v.x + min.x/2, v.y + min.y/2, v.z + min.z/2);
							}
							if (drawVertexCurvatureMax) {
								g.stroke(255, 50, 50);
								Vector3f max = (Vector3f) curv.getPrincipleDirectionMax().clone();
								max.normalize();
								max.scale(1 / 35f);
								g.line(v.x - max.x/2, v.y - max.y/2, v.z - max.z/2, v.x + max.x/2, v.y + max.y/2, v.z + max.z/2);
							}
							g.strokeWeight(1f);
							if (drawVertexCurvature) {
								// min curvature direction - color code is BLUE
								g.stroke(50, 50, 255);
								Vector3f min = (Vector3f) curv.getPrincipleDirectionMin().clone();
								min.normalize();
								min.scale(curv.getCurvatureMin() / 20f);
	//							min.scale(c.getModel().getScale() / (8f * min.length()));
								g.line(v.x - min.x/2, v.y - min.y/2, v.z - min.z/2, v.x + min.x/2, v.y + min.y/2, v.z + min.z/2);
							
								// max curvature direction - color code is RED
								g.stroke(255, 50, 50);
								Vector3f max = (Vector3f) curv.getPrincipleDirectionMax().clone();
								max.normalize();
								max.scale(curv.getCurvatureMax() / 20f);
	//							max.scale(c.getModel().getScale() / (8f * max.length()));
								g.line(v.x - max.x/2, v.y - max.y/2, v.z - max.z/2, v.x + max.x/2, v.y + max.y/2, v.z + max.z/2);
							}
							g.strokeWeight(2f);
						}
					}
				}
				if (drawSharpEdges || drawRegionEdges) {
					g.strokeWeight(3f);
					if (drawSharpEdges && drawRegionEdges) {
						synchronized (c.getModel().getRegions()) {
							for (Region r : c.getModel().getRegions()) {
								for (int i = 0 ; i < r.getBoundaryEdges().size() ; ++i) {
									Edge edge = r.getBoundaryEdges().get(i);
									Vertex v = edge.getVerticesOfEdge()[0];
									Vector3f edgeVal = edge.getEdgeValue();
									if (edge.isSharpEdge()) {
										g.stroke(240, 140, 0);
									}
									else {
										g.stroke(0, 255, 0);
									}
									g.line(v.x, v.y, v.z, v.x - edgeVal.x, v.y - edgeVal.y, v.z - edgeVal.z);
								}
							}
						}
						g.stroke(240, 140, 0);
						synchronized (c.getModel().getTriangles()) {
							for (Triangle t : c.getModel().getTriangles()) {
								Edge[] edges = t.getEdges();
								for (int i = 0 ; i < edges.length ; ++i) {
									if (edges[i].isSharpEdge()) {
										Vertex v = edges[i].getVerticesOfEdge()[0];
										Vector3f edge = edges[i].getEdgeValue();
										g.line(v.x, v.y, v.z, v.x - edge.x, v.y - edge.y, v.z - edge.z);
									}
								}
							}
						}
					}
					else if (drawSharpEdges) {
						g.stroke(240, 140, 0);
						synchronized (c.getModel().getTriangles()) {
							for (Triangle t : c.getModel().getTriangles()) {
								Edge[] edges = t.getEdges();
								for (int i = 0 ; i < edges.length ; ++i) {
									if (edges[i].isSharpEdge()) {
										Vertex v = edges[i].getVerticesOfEdge()[0];
										Vector3f edge = edges[i].getEdgeValue();
										g.line(v.x, v.y, v.z, v.x - edge.x, v.y - edge.y, v.z - edge.z);
									}
								}
							}
						}
					}
					else if (drawRegionEdges) {
						g.stroke(0 , 255, 0);
						synchronized (c.getModel().getRegions()) {
							for (Region r : c.getModel().getRegions()) {
								for (int i = 0 ; i < r.getBoundaryEdges().size() ; ++i) {
									Edge edge = r.getBoundaryEdges().get(i);
									Vertex v = edge.getVerticesOfEdge()[0];
									Vector3f edgeVal = edge.getEdgeValue(); 
									g.line(v.x, v.y, v.z, v.x - edgeVal.x, v.y - edgeVal.y, v.z - edgeVal.z);
								}
							}
						}
					}
					g.strokeWeight(2f);
				}
				if (drawTriangleNormals) {
					g.strokeWeight(1f);
					for (Triangle t : c.getModel().getTriangles()) {
						Vector3f n = (Vector3f) t.getNormalVector().clone();
						n.scale(0.05f);
						g.stroke(35, 120, 200);
						g.line(t.getCentroid().x, t.getCentroid().y, t.getCentroid().z, t.getCentroid().x + n.x, t.getCentroid().y + n.y, t.getCentroid().z + n.z);
					}
					g.strokeWeight(2f);
				}
			}
		}

		// Check if user wants to save current view
		String imgGen = null;
		if (imageGeneratorSettings != null)
			imgGen = imageGeneratorSettings.getAndClearCurrentImageFile();
		if (imageSavePath != null || imgGen != null) {
			if (imgGen != null) {
				save(imgGen);
				imageGeneratorSettings.saveView(cam, false);
				Logger.getRootLogger().info("Image saved as: " + imgGen);
				imageGeneratorSettings.triggerSaved();
			} else {
				save(imageSavePath);

				synchronized (imageSavePath) {
					imageSavePath.notifyAll();
				}
				Logger.getRootLogger().info("Image saved as: " + imageSavePath);
				imageSavePath = null;
			}
		}

	}

	/**
	 * Draw a cone or cylinder with given properties on graphics context
	 * 
	 * @param g
	 *            graphics context
	 * @param sides
	 *            number of sides to use (indicates the level of detail to draw cone)
	 * @param r1
	 *            bottom radius
	 * @param r2
	 *            top radius
	 * @param h
	 *            height
	 * @param top
	 *            draw top cap
	 * @param bottom
	 *            draw bottom cap
	 */
	public static void drawCylinder(PGraphics g, int sides, float r1, float r2, float h,
			boolean top, boolean bottom) {
		float angle = (float) (2 * Math.PI / sides);
		float halfHeight = h / 2;
		if (top) {
			// top
			g.beginShape();
			for (int i = 0; i < sides; i++) {
				float x = (float) (Math.cos(i * angle) * r1);
				float y = (float) (Math.sin(i * angle) * r1);
				g.vertex(x, y, -halfHeight);
			}
			g.endShape(PConstants.CLOSE);
		}
		if (bottom) {
			// bottom
			g.beginShape();
			for (int i = 0; i < sides; i++) {
				float x = (float) (Math.cos(i * angle) * r2);
				float y = (float) (Math.sin(i * angle) * r2);
				g.vertex(x, y, halfHeight);
			}
			g.endShape(PConstants.CLOSE);
		}
		// draw body
		g.beginShape(PConstants.TRIANGLE_STRIP);
		for (int i = 0; i < sides + 1; i++) {
			float x1 = (float) (Math.cos(i * angle) * r1);
			float y1 = (float) (Math.sin(i * angle) * r1);
			float x2 = (float) (Math.cos(i * angle) * r2);
			float y2 = (float) (Math.sin(i * angle) * r2);
			g.vertex(x1, y1, -halfHeight);
			g.vertex(x2, y2, halfHeight);
		}
		g.endShape(PConstants.CLOSE);
	}

	
	/**
	 * Gets camera of mesh reasoning view.
	 * 
	 * @return the camera
	 */
	public PeasyCam getCam() {
		return cam;
	}

	/**
	 * Gets list of all CAS objects (see UIMA Framework)
	 * 
	 * @return ArrayList of MeshCas
	 */
	public ArrayList<MeshCas> getCasList() {
		return casList;
	}

	/**
	 * Gets the control handle of the mesh reasoning view
	 * 
	 * @return the control handle of the view
	 */
	public MeshReasoningViewControl getControl() {
		return control;
	}

	/**
	 * Gets the current model scale of the view
	 * 
	 * @return model scale
	 */
	public float getScale() {
		return modelScale;
	}

	/**
	 * Should background be drawn white?
	 * 
	 * @return the backgroundWhite
	 */
	public boolean isBackgroundWhite() {
		return backgroundWhite;
	}

	/**
	 * Should bounding box be drawn for each group?
	 * 
	 * @return the drawBoundingBox
	 */
	public boolean isDrawBoundingBox() {
		return drawBoundingBox;
	}
	
	/**
	 * Should the min curvature direction of each vertex be displayed?
	 * 
	 * @return the drawVertexCurvatureMin
	 */
	public boolean isDrawVertexCurvatureMin() {
		return drawVertexCurvatureMin;
	}
	
	/**
	 * Should the max curvature direction of each vertex be displayed?
	 * 
	 * @return the drawVertexCurvatureMax
	 */
	public boolean isDrawVertexCurvatureMax() {
		return drawVertexCurvatureMax;
	}

	/**
	 * Should curvature properties for each vertex be drawn?
	 * 
	 * @return the drawVertexCurvature
	 */
	public boolean isDrawVertexCurvature() {
		return drawVertexCurvature;
	}

	/**
	 * Should vertex normals be drawn?
	 * 
	 * @return the drawVertexNormals
	 */
	public boolean isDrawVertexNormals() {
		return drawVertexNormals;
	}
	
	/**
	 * Should triangle normals be drawn?
	 * 
	 * @return the drawTriangleNormals
	 */
	public boolean isDrawTriangleNormals() {
		return drawTriangleNormals;
	}

	/**
	 * Should voronoi area be drawn?
	 * 
	 * @return the drawVoronoiArea
	 */
	public boolean isDrawVoronoiArea() {
		return drawVoronoiArea;
	}
	
	/**
	 * Should the sharp edges be drawn?
	 * 
	 * @return the drawSharpEdges
	 */
	public boolean isDrawSharpEdges() {
		return drawSharpEdges;
	}
	
	/**
	 * Should the region edges be drawn?
	 * 
	 * @return the drawRegionEdges
	 */
	public boolean isDrawRegionEdges() {
		return drawRegionEdges;
	}

	/**
	 * Should only the nearest triangle or all triangles intersecting mouse ray be selected?
	 * 
	 * @return the selectNearestOnly
	 */
	public boolean isSelectNearestOnly() {
		return selectNearestOnly;
	}

	/**
	 * Gets the value of selectTrianglesOnly. True if only triangles will be selected, false if whole
	 * annotation is selected on mouse click.
	 * 
	 * @return selectTrianglesOnly value
	 */
	public boolean isSelectTrianglesOnly() {
		return selectTrianglesOnly;
	}

	@Override
	public void keyTyped(KeyEvent e) {
		char c = e.getKeyChar();
		if (c == '+') {
			userScale *= 1.5;
		} else if (c == '-') {
			userScale = userScale / 1.5f;
		} else if (c == ',') {
			cam.setDistance(cam.getDistance() * 0.8, 250);
		} else if (c == '.') {
			cam.setDistance(cam.getDistance() * 1.2, 250);
		} else if (c == '1') {
			setDrawType(DrawType.FILL);
		} else if (c == '2') {
			setDrawType(DrawType.LINES);
		} else if (c == '3') {
			setDrawType(DrawType.POINTS);
		} else if (c == 'w') {
			drawSettings.incLineWidth();
			draw();
		} else if (c == 'q') {
			drawSettings.decLineWidth();
			draw();
		} else if (imageGeneratorSettings != null && e.isControlDown() && c == 's') {
			// Ctrl+S
			imageGeneratorSettings.triggerViewInitialized();
		}
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
			// Check if clicked on one of previous selected triangles
			for (IntersectedTriangle p : selectedTriangles) {
				Point3f newIntersect = new Point3f();
				if (p.t.intersectsRay(rayEnd, rayStart, newIntersect)) {
					selectedTriangles.clear();
					p.intersection = newIntersect;
					selectedTriangles.add(p);
					found = true;
					break;
				}
			}

			if (!found) {

				// It is new selection
				selectedTriangles.clear();
				for (MeshCas c : casList) {
					if (c.getModel() == null)
						continue;
					ArrayList<IntersectedTriangle> toAdd = new ArrayList<IntersectedTriangle>();
					c.getModel().getGroup().getIntersectedTriangles(rayEnd, rayStart, toAdd);

					// Now check if all triangles in toAdd are currently visible (enabled by
					// checkbox). We only need to check, if main mesh is hidden.
					if (!c.isDrawMesh()) {
						ArrayList<IntersectedTriangle> toRemove = new ArrayList<IntersectedTriangle>();

						for (IntersectedTriangle p : toAdd) {

							boolean isVisible = false;
							for (Annotation a : c.getAnnotations()) {
								if (!(a instanceof DrawableAnnotation))
									continue;
								DrawableAnnotation ma = (DrawableAnnotation) a;
								if (!ma.isDrawAnnotation())
									continue; // Skip not visible annotations

								if (ma.containsTriangle(p.t)) {
									isVisible = true;
									break;
								}
							}
							if (!isVisible)
								toRemove.add(p);
						}
						toAdd.removeAll(toRemove);
					}

					selectedTriangles.addAll(toAdd);
				}
			}

			if (!selectTrianglesOnly) {
				// Check if one of selected triangles is in already selected annotation
				ArrayList<IntersectedTriangle> newSelected = new ArrayList<IntersectedTriangle>();
				for (IntersectedTriangle p : selectedTriangles) {
					synchronized (selectedAnnotations) {
						for (DrawableAnnotation ma : selectedAnnotations.keySet())
							if (ma.containsTriangle(p.t)) {
								newSelected.add(p);
							}
					}
				}
				if (newSelected.size() > 0) {
					// Currently selected was in one or more of the selected annotations, so select
					// out
					// of current annotations
					selectedTriangles.clear();
					selectedTriangles.addAll(newSelected);
				}
			}

			if (selectNearestOnly && selectedTriangles.size() > 1) {
				IntersectedTriangle nearest = null;
				float nearestDist = Float.MAX_VALUE;
				for (IntersectedTriangle p : selectedTriangles) {
					Vector3f camPos = new Vector3f(cam.getPosition());
					camPos.sub(p.intersection);
					float dist = camPos.lengthSquared();
					if (dist < nearestDist) {
						nearestDist = dist;
						nearest = p;
					}
				}
				selectedTriangles.clear();
				selectedTriangles.add(nearest);
				// uncomment below if visual debug of selected triangles is desired
//				System.out.println("========= New selection =========");
//				for (IntersectedTriangle t : selectedTriangles) {
//					System.out.println(t.t);
//					System.out.println("=========== Neighbors ===========");
//					for (Triangle n : t.t.getNeighbors()) {
//						System.out.println(n);
//						System.out.println("aon = " + (Math.toDegrees(t.t.getNormalVector().angle(n.getNormalVector())) + " aoninv = " + (Math.toDegrees(n.getNormalVector().angle(t.t.getNormalVector())))));
//					}
//				}
//				System.out.println("========= End selection =========\n");
			}

			selectedTrianglesChanged();
		}

	}

	/**
	 * Saves the current viewable context into a PNG image file with the given filename. If filename
	 * doesn't end with '.png', it will be added.
	 * 
	 * @param filename
	 *            Filename for the image to save
	 */
	public void saveImage(String filename) {
		saveImage(filename, false);

	}

	/**
	 * Saves the current view as PNG image.
	 * 
	 * @param filename
	 *            Filename of the image
	 * @param overwrite
	 *            overwrite existing image. If false and image with the same name already exists, a
	 *            new name is generated by adding a consecutive number.
	 */
	public void saveImage(String filename, boolean overwrite) {
		if (imageSavePath != null) {
			synchronized (imageSavePath) {
				try {
					imageSavePath.wait();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}

		String path = filename;
		if (!path.endsWith(".png"))
			path += ".png";

		if (!(new File(path).isAbsolute()))
			path = new File(sketchPath + "/images", path).getAbsolutePath();

		String basePath = path;
		int num = 0;

		while (!overwrite && new File(path).exists()) {
			num++;
			path = basePath.substring(0, basePath.lastIndexOf('.')) + "-" + num + ".png";
		}

		imageSavePath = path;
	}

	/**
	 * Updates list of selectedAnnotations after list of selectedTriangless has changed
	 */
	private void selectedTrianglesChanged() {
		if (!selectTrianglesOnly) {
			synchronized (selectedAnnotations) {
				selectedAnnotations.clear();
			}
			for (MeshCas c : casList) {
				for (Annotation a : c.getAnnotations()) {
					if (!(a instanceof DrawableAnnotation))
						continue;

					DrawableAnnotation ma = (DrawableAnnotation) a;
					if (!ma.isDrawAnnotation())
						continue; // Skip not visible annotations

					for (IntersectedTriangle p : selectedTriangles)
						if (ma.containsTriangle(p.t)) {
							synchronized (selectedAnnotations) {
								selectedAnnotations.put(ma, null);
							}
							break;
						}
				}
			}
			control.showSelectedAnnotation(selectedAnnotations);
		}

	}

	/**
	 * Sets background white flag to true if background should be white
	 * 
	 * @param backgroundWhite
	 *            the backgroundWhite to set
	 */
	public void setBackgroundWhite(boolean backgroundWhite) {
		this.backgroundWhite = backgroundWhite;
	}

	/**
	 * Sets list of all CAS objects (see UIMA Framework)
	 * 
	 * @param casList
	 *            list to set
	 */
	public void setCasList(ArrayList<MeshCas> casList) {
		this.casList = casList;
	}

	/**
	 * Sets control panel for this view
	 * 
	 * @param control
	 *            the control to set
	 */
	public void setControl(MeshReasoningViewControl control) {
		this.control = control;
	}

	/**
	 * Sets the bounding box drawing flag
	 * 
	 * @param drawBoundingBox
	 *            the drawBoundingBox to set
	 */
	public void setDrawBoundingBox(boolean drawBoundingBox) {
		this.drawBoundingBox = drawBoundingBox;
	}

	/**
	 * Sets the curvature coloring scheme. Three coloring schemes
	 * are available:
	 * 		Mean Curvature coloring (relative HSV scale converted to RGB),
	 * 		Gaussian Curvature coloring (relative HSV scale converted to RGB),
	 * 		Estimated Curvatures coloring (absolute HSV scale converted to RGB). 
	 * 
	 * @param whichColoring
	 *            setter for the coloring scheme:
	 *            		0 := Mean Curvature coloring
	 *            		1 := Gaussian Curvature coloring
	 *            		2 := Estimated Curvatures coloring
	 */
	public void setDrawCurvatureColor(int whichColoring) {

		for (MeshCas c : casList) {
			float offsetMean = 10;
			float offsetGauss = 20;
			for (Vertex v : c.getModel().getVertices()) {
				if (whichColoring == 0) {
					v.overrideColor = null;
				}
				else if (whichColoring == 1) {
					Curvature curv = c.getCurvature(v);
					v.overrideColor = (curv == null) ? null : curv.getMeanColor(c.getModel().getAvgMeanCurvature() - offsetMean,
							c.getModel().getAvgMeanCurvature() + offsetMean);
				} 
				else if (whichColoring == 2) {
					Curvature curv = c.getCurvature(v);
					v.overrideColor = (curv == null) ? null : curv.getGaussColor(c.getModel().getAvgGaussCurvature() - offsetGauss,
							c.getModel().getAvgGaussCurvature() + offsetGauss);
				}
				else if (whichColoring == 3) {
					Curvature curv = c.getCurvature(v);
					v.overrideColor = (curv == null) ? null : curv.getColor();
				}
			}
		}
	}

	/**
	 * Sets the draw type for the whole model and annotations. Available types are:
	 * 	 	FILLED = faces are filled with color, 
	 * 	 	LINES = faces are transparent, 
	 * 		POINTS = faces and lines are transparent.
	 * 
	 * @param t
	 *            new draw type
	 */
	public void setDrawType(DrawType t) {
		drawSettings.drawType = t;
		draw();
	}

	/**
	 * Sets the normalized minimum curvature direction drawing flag
	 * 
	 * @param drawVertexCurvatureMin
	 *            the drawVertexCurvatureMin to set
	 */
	public void setDrawVertexCurvatureMin(boolean drawVertexCurvatureMin) {
		this.drawVertexCurvatureMin = drawVertexCurvatureMin;
	}
	
	/**
	 * Sets the normalized maximum curvature direction drawing flag
	 * 
	 * @param drawVertexCurvatureMax
	 *            the drawVertexCurvatureMax to set
	 */
	public void setDrawVertexCurvatureMax(boolean drawVertexCurvatureMax) {
		this.drawVertexCurvatureMax = drawVertexCurvatureMax;
	}
	
	/**
	 * Sets the scaled minimum and maximum curvature directions drawing flag.
	 * The scaling performed is based on the associated minimum and maximum 
	 * curvature values.
	 * 
	 * @param drawVertexCurvature
	 *            the drawVertexCurvature to set
	 */
	public void setDrawVertexCurvature(boolean drawVertexCurvature) {
		this.drawVertexCurvature = drawVertexCurvature;
	}

	/**
	 * Sets the normal of each vertex drawing flag
	 * 
	 * @param drawVertexNormals
	 *            the drawVertexNormals to set
	 */
	public void setDrawVertexNormals(boolean drawVertexNormals) {
		this.drawVertexNormals = drawVertexNormals;
	}

	/**
	 * Sets the normal of each triangle drawing flag
	 * 
	 * @param drawTriangleNormals
	 * 			  the drwaTriangleNormals to set
	 */
	public void setDrawTriangleNormals(boolean drawTriangleNormals) {
		this.drawTriangleNormals = drawTriangleNormals;
	}
	
	/**
	 * Sets the Voronoi area associated to each vertex drawing flag
	 * 
	 * @param drawVoronoiArea
	 *            the drawVoronoiArea to set
	 */
	public void setDrawVoronoiArea(boolean drawVoronoiArea) {
		this.drawVoronoiArea = drawVoronoiArea;
	}
	
	/**
	 * Sets the sharp edges drawing flag
	 * 
	 * @param drawSharpEdges
	 * 			  the drawSharpEdges to set
	 */
	public void setDrawSharpEdges(boolean drawSharpEdges) {
		this.drawSharpEdges = drawSharpEdges;
	}
	
	/**
	 * Sets the region edges drawing flag
	 * 
	 * @param drawRegionEdges
	 * 			  the drawRegionEdges to set
	 */
	public void setDrawRegionEdges(boolean drawRegionEdges) {
		this.drawRegionEdges = drawRegionEdges;
	}

	/**
	 * Sets the image generator settings for the image generation. If null, image generation is disabled.
	 * 
	 * @param imageGeneratorSettings
	 *            new image generator settings
	 */
	public void setImageGeneratorSettings(ImageGeneratorSettings imageGeneratorSettings) {
		this.imageGeneratorSettings = imageGeneratorSettings;

	}

	/**
	 * Sets the manual rotation of the camera
	 * 
	 * @param pitch
	 *            rotation around x
	 * @param yaw
	 *            rotation around y
	 * @param roll
	 *            rotation around z
	 */
	public void setManualRotation(float pitch, float yaw, float roll) {
		cam.setRotations(pitch, yaw, roll);
	}

	/**
	 * Sets the model scale to scale model manually
	 * 
	 * @param modelScale
	 *            new scale factor
	 */
	public void setScale(float modelScale) {
		this.modelScale = modelScale;
	}

	/**
	 * Sets the select only the near triangle drawing flag
	 * 
	 * @param selectNearestOnly
	 *            the selectNearestOnly to set
	 */
	public void setSelectNearestOnly(boolean selectNearestOnly) {
		this.selectNearestOnly = selectNearestOnly;
	}

	/**
	 * Set to true if only clicked triangles should be selected, false if whole annotation should be
	 * selected.
	 * 
	 * @param select
	 *            new value for selectTrianglesOnly
	 */
	public void setSelectTrianglesOnly(boolean select) {
		selectTrianglesOnly = select;
	}

	@Override
	public void setup() {
		size(1000, 1000, "org.knowrob.vis.model.view.PAppletSelectionGraphics");

		frameRate(10);
		cam = new PeasyCam(this, 0, 0, 0, 10);
		cam.setMinimumDistance(0.01);
		cam.setMaximumDistance(500);

		cam.setRightDragHandler(cam.getPanDragHandler());

		cam.setDistance(50);
		cam.setWheelScale(0.25);

		cam.rotateX((float) Math.PI / 2f);
		cam.rotateZ((float) Math.PI);
		cam.feed();

		captureViewMatrix();

		perspective();

		draw();
		if (imageGeneratorSettings != null) {
			imageGeneratorSettings.triggerSetup();
		}
	}
}
