/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Stefan Profanter - initial API and implementation, Year: 2012
 ******************************************************************************/
package edu.tum.cs.vis.model.uima.cas;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.Serializable;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.log4j.Logger;
import org.semanticweb.owlapi.apibinding.OWLManager;
import org.semanticweb.owlapi.io.RDFXMLOntologyFormat;
import org.semanticweb.owlapi.model.AddImport;
import org.semanticweb.owlapi.model.IRI;
import org.semanticweb.owlapi.model.OWLClass;
import org.semanticweb.owlapi.model.OWLDataFactory;
import org.semanticweb.owlapi.model.OWLImportsDeclaration;
import org.semanticweb.owlapi.model.OWLNamedIndividual;
import org.semanticweb.owlapi.model.OWLOntology;
import org.semanticweb.owlapi.model.OWLOntologyManager;
import org.semanticweb.owlapi.util.DefaultPrefixManager;

import processing.core.PGraphics;
import edu.tum.cs.ias.knowrob.owl.OWLThing;
import edu.tum.cs.ias.knowrob.owl.utils.OWLFileUtils;
import edu.tum.cs.ias.knowrob.owl.utils.OWLImportExport;
import edu.tum.cs.uima.Annotation;
import edu.tum.cs.uima.JCas;
import edu.tum.cs.vis.model.Model;
import edu.tum.cs.vis.model.uima.annotation.ComplexHandleAnnotation;
import edu.tum.cs.vis.model.uima.annotation.ContainerAnnotation;
import edu.tum.cs.vis.model.uima.annotation.DrawableAnnotation;
import edu.tum.cs.vis.model.uima.annotation.MeshAnnotation;
import edu.tum.cs.vis.model.uima.annotation.primitive.ConeAnnotation;
import edu.tum.cs.vis.model.uima.annotation.primitive.PlaneAnnotation;
import edu.tum.cs.vis.model.uima.annotation.primitive.SphereAnnotation;
import edu.tum.cs.vis.model.util.Curvature;
import edu.tum.cs.vis.model.util.DrawSettings;
import edu.tum.cs.vis.model.util.Triangle;
import edu.tum.cs.vis.model.util.Vertex;

/**
 * UIMA CAS for 3D meshes.
 * 
 * @author Stefan Profanter
 * 
 */
public class MeshCas extends JCas implements Serializable {

	/**
	 * Log4J Logger
	 */
	@SuppressWarnings("unused")
	private static Logger						logger				= Logger.getLogger(MeshCas.class);

	/**
	 * Auto generated
	 */
	private static final long					serialVersionUID	= 4603166505444872760L;
	/**
	 * Group which represents the mesh and its child groups
	 */
	private Model								model;
	/**
	 * Filename of the model represented in this CAS
	 */
	private String								modelFile;

	/**
	 * Indicates if the original mesh should be drawn or not.
	 */
	private boolean								drawMesh			= true;

	/**
	 * Maps curvature property to each vertex of mesh
	 */
	private final HashMap<Vertex, Curvature>	curvatures			= new HashMap<Vertex, Curvature>();

	/**
	 * IRI for exported OWL files
	 */
	public final static String CAD = "http://ias.cs.tum.edu/kb/knowrob-cad.owl#";

	/**
	 * adds a new annotation to the annotations list
	 * 
	 * @param a
	 *            annotation to add
	 */
	public void addAnnotation(Annotation a) {
		synchronized (annotations) {
			annotations.add(a);
		}
	}

	/**
	 * Draw the original mesh and all the meshes of the annotations.
	 * 
	 * @param g
	 *            Applet to draw on
	 * 
	 * @param drawSettings
	 *            draw settings to draw annotation. Used e.g. to override annotation color.
	 */
	public void draw(PGraphics g, DrawSettings drawSettings) {
		if (model == null)
			return;
		if (drawMesh)
			model.draw(g, drawSettings);
		synchronized (annotations) {
			for (Annotation a : annotations) {
				if (!(a instanceof DrawableAnnotation))
					continue;
				DrawableAnnotation ma = (DrawableAnnotation) a;
				ma.draw(g, drawSettings);

			}
		}
	}

	/**
	 * Searches the annotation (type equals clazz) which contains the given triangle. If no
	 * annotation found, null be returned
	 * 
	 * @param clazz
	 *            Type of the annotation to find
	 * @param p
	 *            annotation must contain this triangle
	 * @return the found annotation or null
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public <T extends MeshAnnotation> T findAnnotation(Class<T> clazz, Triangle p) {
		for (Annotation a : getAnnotations()) {
			if (!(a instanceof MeshAnnotation))
				continue;
			MeshAnnotation ma = (MeshAnnotation) a;
			if ((ma.getClass() == clazz || clazz.isInstance(ma))
					&& ma.getMesh().getTriangles().contains(p))
				return (T) ma;
		}
		return null;
	}

	/**
	 * Searches all annotations which are an instance of the given class.
	 * 
	 * @param clazz
	 *            Type of the annotations to find
	 * @return the found annotations or empty
	 */
	@SuppressWarnings("unchecked")
	public <T extends DrawableAnnotation> HashSet<T> findAnnotations(Class<T> clazz) {
		HashSet<T> an = new HashSet<T>();
		for (Annotation a : getAnnotations()) {
			if (!(a instanceof DrawableAnnotation))
				continue;
			DrawableAnnotation ma = (DrawableAnnotation) a;
			if ((ma.getClass() == clazz || clazz.isInstance(ma))) {
				an.add((T) a);
			}
		}
		return an;
	}

	/**
	 * Get curvature property for vertex <tt>v</tt>
	 * 
	 * @param v
	 *            Vertex for which curvature is needed
	 * @return curvature properties or null if none defined
	 */
	public Curvature getCurvature(Vertex v) {
		return curvatures.get(v);
	}

	/**
	 * Get vertex curvature mapping.
	 * 
	 * @return map for vertex to curvature
	 */
	public HashMap<Vertex, Curvature> getCurvatures() {
		return curvatures;
	}

	/**
	 * Get parent model of mesh
	 * 
	 * @return the model
	 */
	public Model getModel() {
		return model;
	}

	/**
	 * Gets all triangles of mesh
	 * 
	 * @return array of triangles
	 */
	public Triangle[] getTriangles() {
		return model.getTriangles().toArray(new Triangle[0]);
	}

	/**
	 * Gets all vertices of mesh
	 * 
	 * @return array of vertices
	 */
	public Vertex[] getVertices() {
		return model.getVertices().toArray(new Vertex[0]);
	}

	/**
	 * Should mesh be drawn or not
	 * 
	 * @return the drawMesh
	 */
	public boolean isDrawMesh() {
		return drawMesh;
	}

	/**
	 * removes the annotation from annotations list
	 * 
	 * @param a
	 *            annotation to remove
	 */
	public void removeAnnotation(Annotation a) {
		synchronized (annotations) {
			annotations.remove(a);
		}
	}

	/**
	 * Set curvature property for given vertex. Overwrites existing one for vertex.
	 * 
	 * @param v
	 *            Vertex to which the curvature property should be assigned
	 * @param c
	 *            curvature property
	 */
	public void setCurvature(Vertex v, Curvature c) {
		curvatures.put(v, c);
	}

	/**
	 * Should mesh be drawn (be visible) or not
	 * 
	 * @param drawMesh
	 *            the drawMesh to set
	 */
	public void setDrawMesh(boolean drawMesh) {
		this.drawMesh = drawMesh;
	}

	/**
	 * Set parent model of mesh.
	 * 
	 * @param model
	 *            the model to set
	 */
	public void setModel(Model model) {
		this.model = model;
	}

	/**
	 * Get path of model file
	 * 
	 * @return Pathname of CAD model file
	 */
	public String getModelFile() {
		return modelFile;
	}

	/**
	 * Set model file path name
	 * 
	 * @param modelFile Pathname of model file
	 */
	public void setModelFile(String modelFile) {
		this.modelFile = modelFile;
	}


	/**
	 * Export the annotations from this CAS to an OWL file
	 * 
	 * @param filename OWL file to which the CAS content shall be exported
	 */
	public void writeToOWL(String filename) {

		OWLDataFactory factory;
		DefaultPrefixManager pm;

		OWLOntology ontology = null;

		try {

			// Create ontology manager and data factory
			OWLOntologyManager manager = OWLManager.createOWLOntologyManager();
			factory = manager.getOWLDataFactory();

			// Get prefix manager using the base IRI of the JoystickDrive ontology as default namespace
			pm = OWLImportExport.PREFIX_MANAGER;

			// Create empty OWL ontology
			ontology = manager.createOntology(IRI.create(CAD));
			OWLImportExport.PREFIX_MANAGER.setPrefix("cad:", CAD);
			manager.setOntologyFormat(ontology, new RDFXMLOntologyFormat());

			// Import KnowRob ontology
			OWLImportsDeclaration oid = factory.getOWLImportsDeclaration(IRI.create(OWLImportExport.KNOWROB));
			AddImport addImp = new AddImport(ontology,oid);
			manager.applyChange(addImp);


			// create object instance based on CAD file name
			OWLClass obj_class = factory.getOWLClass("knowrob:HumanScaleObject", pm);
			OWLNamedIndividual obj_inst = factory.getOWLNamedIndividual(OWLThing.getUniqueID("cad:" + new File(filename).getName().split("\\.")[0]), pm);
			manager.addAxiom(ontology, factory.getOWLClassAssertionAxiom(obj_class, obj_inst));


			// create instances for all object parts
			for(Annotation a : annotations) {

				if(a instanceof ComplexHandleAnnotation) {
					((ComplexHandleAnnotation) a).writeToOWL(obj_inst, manager, factory, pm, ontology, model.getScale());

				} else if(a instanceof ContainerAnnotation) {
					((ContainerAnnotation) a).writeToOWL(obj_inst, manager, factory, pm, ontology, model.getScale());

				} else if(a instanceof SphereAnnotation) {
					((SphereAnnotation) a).writeToOWL(obj_inst, manager, factory, pm, ontology, model.getScale());

				} else if(a instanceof PlaneAnnotation) {
					((PlaneAnnotation) a).writeToOWL(obj_inst, manager, factory, pm, ontology, model.getScale());

				} else if(a instanceof ConeAnnotation) {
					((ConeAnnotation) a).writeToOWL(obj_inst, manager, factory, pm, ontology, model.getScale());
				}
			}


			// export to OWL file
			String owldata = beautifyOWL(OWLFileUtils.saveOntologytoString(ontology, manager.getOntologyFormat(ontology)));

			
			BufferedWriter out = new BufferedWriter(new FileWriter(filename));
			out.write(owldata);
			out.close();
			

		} catch (Exception e) {
			ontology = null;
			e.printStackTrace();
		}
	}

	
	private String beautifyOWL(String owl_data) {
		
		String header = "\n\n" +
		"<!DOCTYPE rdf:RDF [\n" +
	    "    <!ENTITY local_path 'file://@OWL_PATH_PREFIX@/owl/'>\n" +
		"    <!ENTITY owl \"http://www.w3.org/2002/07/owl#\" >\n" +
		"    <!ENTITY xsd \"http://www.w3.org/2001/XMLSchema#\" >\n" +
		"    <!ENTITY knowrob \"http://ias.cs.tum.edu/kb/knowrob.owl#\" >\n" +
		"    <!ENTITY cad \"http://ias.cs.tum.edu/kb/knowrob-cad.owl#\" >\n" +
		"    <!ENTITY rdfs \"http://www.w3.org/2000/01/rdf-schema#\" >\n" +
		"    <!ENTITY rdf \"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" >\n" +
		"]>\n\n<rdf:RDF";
		
		owl_data = owl_data.replace("rdf:resource=\"http://ias.cs.tum.edu/kb/knowrob.owl#", 
				"rdf:resource=\"&knowrob;");
		owl_data = owl_data.replace("rdf:about=\"http://ias.cs.tum.edu/kb/knowrob.owl#", 
				"rdf:about=\"&knowrob;");

		owl_data = owl_data.replace("rdf:resource=\"http://ias.cs.tum.edu/kb/knowrob-cad.owl#", 
				"rdf:resource=\"&cad;");
		owl_data = owl_data.replace("rdf:about=\"http://ias.cs.tum.edu/kb/knowrob-cad.owl#", 
				"rdf:about=\"&cad;");

		owl_data = owl_data.replace("rdf:datatype=\"http://www.w3.org/2001/XMLSchema#", 
									"rdf:datatype=\"&xsd;");
		
		owl_data = owl_data.replace("<owl:imports rdf:resource=\"&knowrob;\"/>", 
				"<owl:imports rdf:resource=\"&local_path;knowrob.owl\"/>");
		
		
		
		owl_data = owl_data.replace("<rdf:RDF", header);
		return owl_data;
	}


}
