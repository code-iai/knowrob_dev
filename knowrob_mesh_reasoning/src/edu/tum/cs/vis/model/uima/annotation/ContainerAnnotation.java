/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Stefan Profanter - initial API and implementation, Year: 2012
 ******************************************************************************/
package edu.tum.cs.vis.model.uima.annotation;

import java.awt.Color;

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

import edu.tum.cs.ias.knowrob.owl.OWLThing;
import edu.tum.cs.ias.knowrob.owl.utils.OWLImportExport;
import edu.tum.cs.vis.model.Model;

/**
 * Mesh annotation for container. A container is like a vessel where a bottom cap and walls are
 * found.
 * 
 * 
 * @author Stefan Profanter
 * 
 */
public class ContainerAnnotation extends MeshAnnotation<ContainerAnnotation> {

	/**
	 * auto generated
	 */
	private static final long	serialVersionUID	= 7483463728879308919L;

	/**
	 * Direction of container. Shows along generating line into direction, where the container is
	 * open. Length of this vector is exactly containers height.
	 */
	private Vector3f			direction;

	/**
	 * total volume of container
	 */
	private float				volume;

	/**
	 * Constructor for a container
	 * 
	 * @param model
	 *            parent model for annotation
	 */
	public ContainerAnnotation(Model model) {
		super(ContainerAnnotation.class, model, new Color(27, 93, 27));
	}

	/**
	 * Direction of container. Shows along generating line into direction, where the container is
	 * open. Length of this vector is exactly containers height.
	 * 
	 * @return the direction
	 */
	public Vector3f getDirection() {
		return direction;
	}

	/**
	 * Direction of container. Shows along generating line into direction, where the container is
	 * open. Length of this vector is exactly containers height as unscaled value.
	 * 
	 * @return the direction unscaled
	 */
	public Tuple3f getDirectionUnscaled() {
		return model.getUnscaled(direction);
	}

	/**
	 * Volume of container
	 * 
	 * @return the volume
	 */
	public float getVolume() {
		return volume;
	}

	/**
	 * Volume of container as unscaled value.
	 * 
	 * @return the volume
	 */
	public float getVolumeUnscaled() {
		return model.getUnscaled(volume);
	}

	/**
	 * Set direction of container. Length should be height of the container and vector should point
	 * into opening direction.
	 * 
	 * @param direction
	 *            the direction to set
	 */
	public void setDirection(Vector3f direction) {
		this.direction = direction;
	}

	/**
	 * Set volume of container
	 * 
	 * @param volume
	 *            the volume to set
	 */
	public void setVolume(float volume) {
		this.volume = volume;
	}


	/**
	 * Export this annotation to an OWL description
	 * 
	 * @param manager OWL ontology manager
	 * @param factory OWL data factory
	 * @param pm Prefix manager
	 * @param ontology Ontology to which the assertions shall be added
	 * @return Reference to an OWLIndividual for this annotation
	 */
	public OWLIndividual writeToOWL(OWLIndividual obj_inst, OWLOntologyManager manager, OWLDataFactory factory, DefaultPrefixManager pm, OWLOntology ontology) {

		OWLClass part_class = factory.getOWLClass("knowrob:Container", pm);
		OWLNamedIndividual part_inst = factory.getOWLNamedIndividual(OWLThing.getUniqueID("knowrob:Container"), pm);
		manager.addAxiom(ontology, factory.getOWLClassAssertionAxiom(part_class, part_inst));

		// set as physicalPart of parent object
		OWLObjectProperty properPhysicalParts = factory.getOWLObjectProperty("knowrob:properPhysicalParts", pm);
		manager.addAxiom(ontology, factory.getOWLObjectPropertyAssertionAxiom(properPhysicalParts, obj_inst, part_inst));

		// TODO: Containers do not have a pose -- once this is added, export it!

		manager.addAxiom(ontology, factory.getOWLDataPropertyAssertionAxiom(
				factory.getOWLDataProperty("knowrob:volumeOfObject",  pm), 
				part_inst, 
				getVolume()));
		
		
		// create direction vector instance
		OWLIndividual vec_inst = OWLImportExport.createDirVector(new Vector3d(getDirection()), manager, factory, pm, ontology);
		manager.addAxiom(ontology, factory.getOWLObjectPropertyAssertionAxiom(
				factory.getOWLObjectProperty("knowrob:longitudinalDirection",  pm), part_inst, vec_inst));
		
		
		return part_inst;
		
	}

}
