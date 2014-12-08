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
import java.util.LinkedHashMap;
import java.util.Set;

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
import org.knowrob.vis.model.uima.annotation.HandleAnnotation;
import org.knowrob.vis.model.uima.annotation.PrimitiveAnnotation;
import org.knowrob.vis.model.uima.cas.MeshCas;
import org.knowrob.vis.model.util.Curvature;
import org.knowrob.vis.model.util.Vertex;

/**
 * Annotation for primitive type: Cone (convex and concave)
 * 
 * @author Stefan Profanter
 * 
 */
public final class ConeAnnotation extends PrimitiveAnnotation<ConeAnnotation> implements
		HandleAnnotation {

	/**
	 * auto generated
	 */
	private static final long	serialVersionUID	= -7420446109108464883L;

	/**
	 * Cone representing the cone annotation.
	 */
	private final Cone			cone;

	/**
	 * Create new cone annotation.
	 * 
	 * @param curvatures
	 *            Map of curvatures for vertices
	 * @param model
	 *            parent model
	 * @param concave
	 *            is cone concave or convex
	 * 
	 */
	public ConeAnnotation(HashMap<Vertex, Curvature> curvatures, Model model, boolean concave) {
		super(ConeAnnotation.class, curvatures, model, concave ? new Color(0, 125, 125)
				: new Color(255, 255, 0));
		cone = new Cone(concave);
	}

	/**
	 * Create new cone annotation.
	 * 
	 * @param curvatures
	 *            Map of curvatures for vertices
	 * @param model
	 *            parent model
	 * @param concave
	 *            is cone concave or convex
	 * @param annotationColor
	 *            User defined color of the cone annotation.
	 */
	public ConeAnnotation(HashMap<Vertex, Curvature> curvatures, Model model, boolean concave,
			Color annotationColor) {
		super(ConeAnnotation.class, curvatures, model, annotationColor);
		cone = new Cone(concave);
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation#drawAnnotation(processing.core.PGraphics)
	 */
	@Override
	public void drawPrimitiveAnnotation(PGraphics g, Color color) {

		cone.draw(g, color == null ? getDrawColor() : color);

	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation#fitAnnotation()
	 */
	@Override
	public boolean fitAnnotation() {
		/*if (isConcave())
			return true;

		if (getArea() < 0.08)
			return true;*/

		LinkedHashMap<Vertex, Float> vertices = new LinkedHashMap<Vertex, Float>();
		Vector3f centroid = getVerticesWithWeight(vertices);
		return cone.fit(centroid, vertices.keySet(), vertices, mesh.getTriangles());
	}

	/**
	 * get centroid of cone at unscaled position
	 * 
	 * @return the centroid
	 */
	public Tuple3f getCentroid() {
		return model.getUnscaled(cone.getCentroid());
	}

	@Override
	public Cone getCone() {
		return cone;
	}

	/**
	 * 
	 * get direction (unscaled) of cone. Direction is aligned with generating line and shows from
	 * centroid to small radius cap. Length of direction is half height of the cone (center to one
	 * end).
	 * 
	 * @return the direction
	 */
	public Vector3f getDirection() {
		return new Vector3f(model.getUnscaled(cone.getDirection()));
	}

	/**
	 * Get total unscaled height of cone from bottom cap to top cap which is 2*directionUnscaled.
	 * 
	 * @return unscaled height
	 */
	public float getHeight() {
		return getDirection().length() * 2;
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.MeshAnnotation#getNeighborAnnotations(edu.tum.cs.vis.model.uima.cas.MeshCas)
	 */
	@Override
	public Set<ConeAnnotation> getNeighborAnnotations(MeshCas cas) {
		return getNeighborAnnotations(cas, ConeAnnotation.class);
	}

	/* (non-Javadoc)
	 * @see edu.tum.cs.vis.model.uima.annotation.PrologBaseAnnotation#getPoseMatrix()
	 */
	@Override
	public Matrix4f getPoseMatrix() {

		Matrix4f mat = cone.getPoseMatrix();
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
		return model.getUnscaled(cone.getArea());
	}

	/**
	 * Get average radius (unscaled) of cone which is the average between small and large radius
	 * 
	 * @return average radius unscaled
	 */
	public float getRadiusAvg() {
		return model.getUnscaled(cone.getRadiusAvg());
	}

	/**
	 * Get large radius (unscaled), which is at the bottom of cone.
	 * 
	 * @return the radiusLarge
	 */
	public float getRadiusLarge() {
		return model.getUnscaled(cone.getRadiusLarge());
	}

	/**
	 * Get small radius (unscaled), which is at the bottom of cone.
	 * 
	 * @return the radiusSmall
	 */
	public float getRadiusSmall() {
		return model.getUnscaled(cone.getRadiusSmall());
	}

	/**
	 * Get unscaled volume of cone.
	 * 
	 * @return the volume unscaled.
	 */
	public float getVolume() {

		return model.getUnscaled(cone.getVolume());
	}

	/**
	 * Is cone concave or convex?
	 * 
	 * @return the concave
	 */
	public boolean isConcave() {
		return cone.isConcave();
	}
	
	
	public OWLIndividual writeToOWL(OWLIndividual obj_inst, OWLOntologyManager manager, OWLDataFactory factory, DefaultPrefixManager pm, OWLOntology ontology, float scale) {

		OWLClass part_class = factory.getOWLClass("knowrob:Cone", pm);
		OWLNamedIndividual part_inst = factory.getOWLNamedIndividual(OWLThing.getUniqueID("knowrob:Cone"), pm);
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
				factory.getOWLDataProperty("knowrob:radius",  pm), 
				part_inst, 
				getRadiusAvg()));
		manager.addAxiom(ontology, factory.getOWLDataPropertyAssertionAxiom(
				factory.getOWLDataProperty("knowrob:minRadius",  pm), 
				part_inst, 
				getRadiusSmall()));
		manager.addAxiom(ontology, factory.getOWLDataPropertyAssertionAxiom(
				factory.getOWLDataProperty("knowrob:maxRadius",  pm), 
				part_inst, 
				getRadiusLarge()));

		manager.addAxiom(ontology, factory.getOWLDataPropertyAssertionAxiom(
				factory.getOWLDataProperty("knowrob:lengthOfObject",  pm), 
				part_inst, 
				getHeight()));


		manager.addAxiom(ontology, factory.getOWLDataPropertyAssertionAxiom(
				factory.getOWLDataProperty("knowrob:volumeOfObject",  pm), 
				part_inst, 
				getVolume()));

		manager.addAxiom(ontology, factory.getOWLDataPropertyAssertionAxiom(
				factory.getOWLDataProperty("knowrob:areaCoverage",  pm), 
				part_inst, 
				getAreaCoverage()));

		manager.addAxiom(ontology, factory.getOWLDataPropertyAssertionAxiom(
				factory.getOWLDataProperty("knowrob:areaOfObject",  pm), 
				part_inst, 
				getArea()));
		


		// create direction vector instance
		Vector3f d = getDirection();
		OWLIndividual vec_inst = OWLImportExport.createDirVector(new Vector3d(d.x, d.y, d.z), manager, factory, pm, ontology);
		
		manager.addAxiom(ontology, factory.getOWLObjectPropertyAssertionAxiom(
				factory.getOWLObjectProperty("knowrob:longitudinalDirection",  pm), part_inst, vec_inst));
		
		return part_inst;
		
	}
}
