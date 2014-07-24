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
import java.util.HashSet;

import javax.vecmath.Matrix4d;
import javax.vecmath.Matrix4f;
import javax.vecmath.Tuple3f;
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
 * Primitive annotation for concave or convex sphere.
 * 
 * @author Stefan Profanter
 * 
 */
public final class SphereAnnotation extends PrimitiveAnnotation<SphereAnnotation> {

	/**
	 * auto generated
	 */
	private static final long	serialVersionUID	= 4870579150170536881L;

	/**
	 * Sphere representing annotation
	 */
	private final Sphere		sphere;

	/**
	 * Creates a new sphere annotation
	 * 
	 * @param curvatures
	 *            Map of curvatures for vertices
	 * @param model
	 *            parent model
	 * @param concave
	 *            is sphere concave or convex
	 */
	public SphereAnnotation(HashMap<Vertex, Curvature> curvatures, Model model, boolean concave) {
		super(SphereAnnotation.class, curvatures, model, concave ? new Color(0, 255, 0)
				: new Color(255, 0, 0));
		sphere = new Sphere(concave);
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation#drawAnnotation(processing.core.PGraphics)
	 */
	@Override
	public void drawPrimitiveAnnotation(PGraphics g, Color color) {

		sphere.draw(g, color == null ? getDrawColor() : color);
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation#fitAnnotation()
	 */
	@Override
	public boolean fitAnnotation() {

		HashMap<Vertex, Float> weights = new HashMap<Vertex, Float>();
		getVerticesWithWeight(weights);
		HashSet<Vertex> vert = new HashSet<Vertex>();
		Vector3f centroid = getVertices(vert);
		return sphere.fit(centroid, weights.keySet(), weights, mesh.getTriangles());
	}

	/**
	 * Get sphere center unscaled
	 * 
	 * @return the center
	 */
	public Tuple3f getCenter() {
		return model.getUnscaled(sphere.getCenter());
	}

	/**
	 * Since a sphere is rotationally symmetric, the orientation does not make much sense (is chosen
	 * to be unit orientation here), but for compatibility reasons, this method is nevertheless
	 * implemented.
	 * 
	 * @return 4x4 pose matrix of the plane relative to the object centroid
	 */
	@Override
	public Matrix4f getPoseMatrix() {
		Matrix4f mat = sphere.getPoseMatrix();
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
		return model.getUnscaled(sphere.getArea());
	}

	/**
	 * get unscaled radius of sphere
	 * 
	 * @return the radius
	 */
	public float getRadius() {
		return model.getUnscaled(sphere.getRadius());
	}

	/**
	 * Get fitted sphere for annotation.
	 * 
	 * WARNING: Sphere parameters are scaled. To get unscaled values put them through
	 * model.getUnscaled().
	 * 
	 * @return fitted sphere
	 */
	public Sphere getSphere() {
		return sphere;
	}

	/**
	 * Get unscaled sphere volume
	 * 
	 * @return the unscaled volume of the sphere
	 */
	public float getVolume() {

		return model.getUnscaled(sphere.getVolume());
	}

	/**
	 * Is sphere concave or convex?
	 * 
	 * @return true if concave
	 */
	public boolean isConcave() {
		return sphere.isConcave();
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

		OWLClass part_class = factory.getOWLClass("knowrob:Sphere", pm);
		OWLNamedIndividual part_inst = factory.getOWLNamedIndividual(OWLThing.getUniqueID("knowrob:Sphere"), pm);
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
		
		
		if(isConcave()) {
			OWLClass concave_class = factory.getOWLClass("knowrob:ConcaveTangibleObject", pm);
			manager.addAxiom(ontology, factory.getOWLClassAssertionAxiom(concave_class, part_inst));
		}

		manager.addAxiom(ontology, factory.getOWLDataPropertyAssertionAxiom(
				factory.getOWLDataProperty("knowrob:radius",  pm), 
				part_inst, 
				getRadius()));
		
		manager.addAxiom(ontology, factory.getOWLDataPropertyAssertionAxiom(
				factory.getOWLDataProperty("knowrob:volumeOfObject",  pm), 
				part_inst, 
				getVolume()));

		manager.addAxiom(ontology, factory.getOWLDataPropertyAssertionAxiom(
				factory.getOWLDataProperty("knowrob:areaOfObject",  pm), 
				part_inst, 
				getArea()));

		manager.addAxiom(ontology, factory.getOWLDataPropertyAssertionAxiom(
				factory.getOWLDataProperty("knowrob:areaCoverage",  pm), 
				part_inst, 
				getAreaCoverage()));
				
		return part_inst;
		
	}
}
