/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Stefan Profanter - initial API and implementation, Year: 2012
 ******************************************************************************/
package org.knowrob.vis.model.uima.annotation.primitive;

import java.awt.Color;
import java.util.HashMap;

import javax.vecmath.Matrix4d;
import javax.vecmath.Matrix4f;
import javax.vecmath.Tuple3f;
import javax.vecmath.Vector3d;
import javax.vecmath.Vector3f;

import org.semanticweb.owlapi.model.OWLClass;
import org.semanticweb.owlapi.model.OWLDataFactory;
import org.semanticweb.owlapi.model.OWLIndividual;
import org.semanticweb.owlapi.model.OWLNamedIndividual;
import org.semanticweb.owlapi.model.OWLObjectProperty;
import org.semanticweb.owlapi.model.OWLOntology;
import org.semanticweb.owlapi.model.OWLOntologyManager;
import org.semanticweb.owlapi.util.DefaultPrefixManager;

import processing.core.PGraphics;
import org.knowrob.owl.OWLThing;
import org.knowrob.owl.utils.OWLImportExport;
import org.knowrob.vis.model.Model;
import org.knowrob.vis.model.uima.annotation.PrimitiveAnnotation;
import org.knowrob.vis.model.util.Curvature;
import org.knowrob.vis.model.util.Vertex;

/**
 * @author Stefan Profanter
 * 
 */
public final class PlaneAnnotation extends PrimitiveAnnotation<PlaneAnnotation> {

	/**
	 * auto generated
	 */
	private static final long	serialVersionUID	= 7758656289829843165L;

	/**
	 * Plane representing the plane annotation
	 */
	private final Plane			plane;

	/**
	 * Creates a new plane annotation.
	 * 
	 * @param curvatures
	 *            Map of curvatures for vertices
	 * @param model
	 *            parent model
	 */
	public PlaneAnnotation(HashMap<Vertex, Curvature> curvatures, Model model) {
		super(PlaneAnnotation.class, curvatures, model, new Color(250, 200, 255));
		plane = new Plane();
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation#drawAnnotation(processing.core.PGraphics)
	 */
	@Override
	public void drawPrimitiveAnnotation(PGraphics g, Color color) {
		plane.draw(g, color == null ? getDrawColor() : color);

	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation#fitAnnotation()
	 */
	@Override
	public boolean fitAnnotation() {
		HashMap<Vertex, Float> vertices = new HashMap<Vertex, Float>();
		Vector3f centroid = getVerticesWithWeight(vertices);

		return plane.fit(centroid, vertices.keySet(), vertices, mesh.getTriangles());
	}

	/**
	 * Get centroid of plane at unscaled position
	 * 
	 * @return the centroid
	 */
	public Tuple3f getCentroid() {
		return model.getUnscaled(plane.getCentroid());
	}

	/**
	 * Get unsaled plane corners
	 * 
	 * @return the corner
	 */
	public Tuple3f[] getCorner() {

		return model.getUnscaled(plane.getCorner());
	}

	/**
	 * Get unscaled vector (length, direction) of long side.
	 * 
	 * @return the longSide
	 */
	public Vector3f getLongSide() {
		return model.getUnscaled(plane.getLongSide());
	}

	/**
	 * Get the unerlying plane.
	 * 
	 * WARNING: Sphere parameters are scaled. To get unscaled values put them through
	 * model.getUnscaled().
	 * 
	 * @return fitted plane for annotation
	 */
	public Plane getPlane() {
		return plane;
	}

	/**
	 * Get plane normal
	 * 
	 * @return the planeNormal
	 */
	public Vector3f getPlaneNormal() {
		return plane.getPlaneNormal();
	}

	/**
	 * Get pose matrix for plane.
	 * 
	 * @return 4x4 pose matrix of the plane relative to the object centroid
	 */
	@Override
	public Matrix4f getPoseMatrix() {
		Matrix4f mat = plane.getPoseMatrix();
		Vector3f pos = new Vector3f(); 
		mat.get(pos);
		mat.setTranslation(model.getUnscaled(pos));
		return mat;
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation#getPrimitiveArea()
	 */
	@Override
	public float getPrimitiveArea() {
		return model.getUnscaled(plane.getArea());
	}

	/**
	 * Get unscaled vector (length, direction) of long side.
	 * 
	 * @return the shortSide
	 */
	public Vector3f getShortSide() {
		return model.getUnscaled(plane.getShortSide());
	}


	/**
	 * Export this annotation to an OWL description
	 * 
	 * @param manager OWL ontology manager
	 * @param factory OWL data factory
	 * @param pm Prefix manager
	 * @param ontology Ontology to which the assertions shall be added
	 * @param scale 
	 * @return Reference to an OWLIndividual for this annotation
	 */
	public OWLIndividual writeToOWL(OWLIndividual obj_inst, OWLOntologyManager manager, OWLDataFactory factory, DefaultPrefixManager pm, OWLOntology ontology, float scale) {

		OWLClass part_class = factory.getOWLClass("knowrob:FlatPhysicalSurface", pm);
		OWLNamedIndividual part_inst = factory.getOWLNamedIndividual(OWLThing.getUniqueID("knowrob:FlatPhysicalSurface"), pm);
		manager.addAxiom(ontology, factory.getOWLClassAssertionAxiom(part_class, part_inst));

		// set as physicalPart of parent object
		OWLObjectProperty properPhysicalParts = factory.getOWLObjectProperty("knowrob:properPhysicalParts", pm);
		manager.addAxiom(ontology, factory.getOWLObjectPropertyAssertionAxiom(properPhysicalParts, obj_inst, part_inst));

		// set relative pose of this object part w.r.t. to the main object
		OWLIndividual mat_inst = OWLImportExport.createPoseInst(new Matrix4d(getPoseMatrix()), manager, factory, pm, ontology);
		manager.addAxiom(ontology, factory.getOWLObjectPropertyAssertionAxiom(
				factory.getOWLObjectProperty("knowrob:orientation",  pm), part_inst, mat_inst));
		manager.addAxiom(ontology, factory.getOWLObjectPropertyAssertionAxiom(
				factory.getOWLObjectProperty("knowrob:relativeTo",  pm), mat_inst, obj_inst));
		
		
		manager.addAxiom(ontology, factory.getOWLDataPropertyAssertionAxiom(
				factory.getOWLDataProperty("knowrob:areaCoverage",  pm), 
				part_inst, 
				getAreaCoverage()));

		manager.addAxiom(ontology, factory.getOWLDataPropertyAssertionAxiom(
				factory.getOWLDataProperty("knowrob:areaOfObject",  pm), 
				part_inst, 
				getArea()));


		// plane normal vector
		OWLIndividual vec_inst = OWLImportExport.createDirVector(new Vector3d(getPlaneNormal()), manager, factory, pm, ontology);
		
		manager.addAxiom(ontology, factory.getOWLObjectPropertyAssertionAxiom(
				factory.getOWLObjectProperty("knowrob:normalDirection",  pm), part_inst, vec_inst));

		// vector for long side of the plane
		vec_inst = OWLImportExport.createDirVector(new Vector3d(getLongSide()), manager, factory, pm, ontology);
		
		manager.addAxiom(ontology, factory.getOWLObjectPropertyAssertionAxiom(
				factory.getOWLObjectProperty("knowrob:objectLongSide",  pm), part_inst, vec_inst));

		// vector for short side of the plane
		vec_inst = OWLImportExport.createDirVector(new Vector3d(getShortSide()), manager, factory, pm, ontology);
		
		manager.addAxiom(ontology, factory.getOWLObjectPropertyAssertionAxiom(
				factory.getOWLObjectProperty("knowrob:objectShortSide",  pm), part_inst, vec_inst));
		
		return part_inst;
		
	}
}
