/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: 
 * 		Stefan Profanter - initial API and implementation, Year: 2012
 * 		Andrei Stoica - refactored implementation during Google Summer of Code 2014	
 ******************************************************************************/
package org.knowrob.vis.model;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.JFrame;

import org.apache.log4j.Logger;
import org.apache.log4j.xml.DOMConfigurator;

import org.knowrob.utils.ResourceRetriever;
import org.knowrob.vis.model.ItemModel;
import org.knowrob.vis.model.Model;
import org.knowrob.vis.model.uima.analyser.ComplexHandleAnalyser;
import org.knowrob.vis.model.uima.analyser.ContainerAnalyser;
import org.knowrob.vis.model.uima.analyser.EdgeAnalyser;
import org.knowrob.vis.model.uima.analyser.MeshAnalyser;
import org.knowrob.vis.model.uima.analyser.NeighborAnalyser;
import org.knowrob.vis.model.uima.analyser.PrimitiveAnalyser;
import org.knowrob.vis.model.uima.annotation.ComplexHandleAnnotation;
import org.knowrob.vis.model.uima.annotation.ContainerAnnotation;
import org.knowrob.vis.model.uima.annotation.DrawableAnnotation;
import org.knowrob.vis.model.uima.annotation.HandleAnnotation;
import org.knowrob.vis.model.uima.annotation.MeshAnnotation;
import org.knowrob.vis.model.uima.annotation.PrimitiveAnnotation;
import org.knowrob.vis.model.uima.annotation.primitive.ConeAnnotation;
import org.knowrob.vis.model.uima.annotation.primitive.PlaneAnnotation;
import org.knowrob.vis.model.uima.annotation.primitive.SphereAnnotation;
import org.knowrob.vis.model.uima.cas.MeshCas;
import org.knowrob.vis.model.util.ContainerAnnotationVolumeComparator;
import org.knowrob.vis.model.util.HandleComparator;
import org.knowrob.vis.model.util.PrimitiveAnnotationAreaComparator;
import org.knowrob.vis.model.util.Triangle;
import org.knowrob.vis.model.util.Vertex;
import org.knowrob.vis.model.util.algorithm.CasProcessing;
import org.knowrob.vis.model.util.algorithm.CurvatureCalculation;
import org.knowrob.vis.model.view.MeshReasoningView;
import org.knowrob.vis.model.view.MeshReasoningViewControl;
import org.knowrob.vis.model.view.PAppletSelectionGraphics;
import org.knowrob.vis.tools.ImageGenerator.ImageGeneratorAction;
import org.knowrob.vis.tools.ImageGenerator.ImageGeneratorSettings;
import org.knowrob.vis.uima.Annotation;
import org.knowrob.vis.util.PrintUtil;


/**
 * Main mesh reasoning class for parsing and analyzing CAD models. Provides methods for starting mesh
 * reasoning with and without GUI and implements the general workflow followed by the reasoning 
 * process on the CAD models. It implements both the backend functionality and frontend piece 
 * of the knowrob_mesh_reasoning add-on.
 * 
 * @author Stefan Profanter
 * @author Andrei Stoica (refactored the general processing method analyseByPath)
 */
public class MeshReasoning {

	/**
	 * Log4j logger
	 */
	private final static Logger	LOGGER	= Logger.getLogger(MeshReasoning.class);
	
	/**
	 * Settings for image generator. If null, image generation is disabled.
	 */
	private ImageGeneratorSettings		imageGeneratorSettings	= null;

	/**
	 * View for this mesh reasoning object
	 */
	private MeshReasoningView			mrv						= null;

	/**
	 * Mesh reasoning container
	 */
	private MeshCas						cas						= null;

	/**
	 * Main frame container for mesh reasoning view
	 */
	public JFrame						frame;

	/**
	 * View control for mesh reasoning view
	 */
	private MeshReasoningViewControl	control;

	
	/**
	 * Constructor for mesh reasoning object. Initializes object and creates mesh reasoning view if
	 * indicated.
	 * 
	 * @param withView
	 *            set to true if mesh reasoning view should be created and shown
	 * 
	 */
	public MeshReasoning(boolean withView) {
		this(withView, null);
	}

	/**
	 * Constructor for mesh reasoning object. Initializes object and creates mesh reasoning view if
	 * indicated.
	 * 
	 * @param withView
	 *            set to true if mesh reasoning view should be created and shown
	 * @param imageGenerator
	 *            ImageGeneratorSettings for image generation. Set to null if image generation not
	 *            used.
	 * 
	 */
	public MeshReasoning(boolean withView, ImageGeneratorSettings imageGenerator) {
		imageGeneratorSettings = imageGenerator;
		if (imageGeneratorSettings != null)
			initImageGenerator();
		// initImageGeneratorSpoon();
		cas = new MeshCas();

		if (withView) {
			frame = new JFrame();
			frame.setMinimumSize(new Dimension(800, 600));
			frame.setSize(1324, 768);
			frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			frame.setTitle("Mesh reasoning view");
			frame.setLocationRelativeTo(null);

			if (imageGeneratorSettings != null && imageGeneratorSettings.isRunBackground())
				frame.setState(java.awt.Frame.ICONIFIED);

			ArrayList<MeshAnalyser> analyser = new ArrayList<MeshAnalyser>();

			mrv = new MeshReasoningView();
			mrv.setImageGeneratorSettings(imageGeneratorSettings);
			control = new MeshReasoningViewControl(cas, analyser, mrv);
			mrv.setControl(control);
			mrv.init();

			mrv.getCasList().add(cas);

			frame.setLayout(new BorderLayout());
			frame.getContentPane().add(mrv, BorderLayout.CENTER);
			frame.getContentPane().add(control, BorderLayout.LINE_END);
			frame.setVisible(true);
		}

	}

	/**
	 * Main initialization method for creating mesh reasoning object. Constructs mesh reasoning
	 * object and initializes log4j logger.
	 * 
	 * @param withView
	 *            Also create GUI to visualize mesh reasoning
	 * @return new mesh reasoning object
	 */
	public static MeshReasoning initMeshReasoning(boolean withView) {
		return MeshReasoning.initMeshReasoning(withView, null);
	}

	/**
	 * Initialize mesh reasoning object for image generation.
	 * 
	 * @param withView
	 *            Also create GUI to visualize mesh reasoning. Should be true if image generation is
	 *            used.
	 * @param imageGenerator
	 *            ImageGeneratorSettings for image generation
	 * @return new Mesh Reasoning object.
	 */
	public static MeshReasoning initMeshReasoning(boolean withView,
			ImageGeneratorSettings imageGenerator) {
		DOMConfigurator.configureAndWatch("log4j.xml", 60 * 1000);
		return new MeshReasoning(withView, imageGenerator);
	}

	/**
	 * Start mesh reasoning on specified file path
	 * 
	 * @param path
	 *            path to CAD model. Can be physical file path or http://, ftp:// or even package://
	 *            which indicates a ros package
	 * 
	 */
	public void analyseByPath(String path) {

		if (imageGeneratorSettings != null) {
			// update file name
			imageGeneratorSettings.setCurrentModel(path.substring(imageGeneratorSettings
					.getInputBasePath().getAbsolutePath().length() + 1));
		}
		LOGGER.info("MeshReasoning started. Parsing model ...");
		LOGGER.debug("Path: " + path);
		long start = System.currentTimeMillis();

		// Load and parse model
		ItemModel itemModel = new ItemModel(path);

		if (itemModel.getParser() == null) {
			LOGGER.error("Couldn't parse model. Maybe path of model file is wrong.");
			return;
		}

		Model model = itemModel.getParser().getModel();
		LOGGER.debug("Model parsed. Took: "
				+ PrintUtil.prettyMillis(System.currentTimeMillis() - start) + " (Vertices: "
				+ model.getVertices().size() + ", Lines: " + model.getLines().size()
				+ ", Triangles: " + model.getTriangles().size() + ")");

		start = System.currentTimeMillis();

		// normalize model for further reasoning
		model.normalize();
		
		// list of current running analyzers used in mesh reasoning view
		ArrayList<MeshAnalyser> analyser;
		if (mrv != null) {
			analyser = mrv.getControl().getAnalyser();
		} else {
			cas = new MeshCas();
			analyser = new ArrayList<MeshAnalyser>(6);
		}
		
		// set model to MeshCas
		cas.setModel(model);
		
		// remember model path (e.g. for saving cache files)
		if (path.indexOf("://") <= 0) { // Is local file
			cas.setModelFile(path);	
		} else if (path.startsWith("package://")) {
			int st = path.indexOf('/') + 2;
			int end = path.indexOf('/', st);
			String serverName = path.substring(st, end);
			String filePath = path.substring(end + 1);
			String pkgPath = ResourceRetriever.findPackage(serverName);
			if (pkgPath != null)
				cas.setModelFile(pkgPath + "/" + filePath);
		}
		
		// analyze mesh by looking into its neighboring relation ships
		// and remove double sided or collinear triangles while updating
		// the vertex sharing of the model and the normal vectors of the
		// model
		NeighborAnalyser na = new NeighborAnalyser();
		analyser.add(na);
		Thread.yield();
		na.process(cas);
		
		// perform sharp edge detection on the model
		EdgeAnalyser ea = new EdgeAnalyser();
		analyser.add(ea);
		Thread.yield();
		ea.process(cas,imageGeneratorSettings);

		if (imageGeneratorSettings != null) {
			imageGeneratorSettings.waitSetup();

			if (imageGeneratorSettings.isInitViewFromFile()) {
				if (imageGeneratorSettings.initView(mrv.getCam())) {
					try {
						Thread.sleep(2500); // wait until model correctly shown
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}
			}

			if (imageGeneratorSettings.isSaveView() && !imageGeneratorSettings.isViewInitialized()) {
				// allow user to set a viewpoint
				LOGGER.info("Set the desired view and then press Ctrl+S to continue");
				imageGeneratorSettings.waitViewInitialized();
			}

			if (imageGeneratorSettings.isSavePlainModel()) {
				// wait until model is saved
				imageGeneratorSettings.waitSaved("plain");
			}
		}

		// estimate curvature values
		LOGGER.debug("Calculating curvature estimates ...");
		long curvatureStartTime = System.currentTimeMillis();
		CurvatureCalculation.calculateCurvatures(cas.getCurvatures(), cas.getModel());
		long curvatureDuration = System.currentTimeMillis() - curvatureStartTime;
		LOGGER.debug("Ended. Took: " + PrintUtil.prettyMillis(curvatureDuration));
	
		if (imageGeneratorSettings != null && imageGeneratorSettings.isSaveCurvatureColor()) {
			// draw estimated curvature coloring
			mrv.setDrawCurvatureColor(3);
			imageGeneratorSettings.waitSaved("curvature");
			mrv.setDrawCurvatureColor(0);
		}		
		
		// process the estimated curvature information
		CasProcessing casProcessor = new CasProcessing(cas);
		casProcessor.process();
		
		// create analysers and add them to the analyser list
		PrimitiveAnalyser pa = new PrimitiveAnalyser();
		analyser.add(pa);
		ContainerAnalyser ca = new ContainerAnalyser();
		analyser.add(ca);
		ComplexHandleAnalyser cha = new ComplexHandleAnalyser();
		analyser.add(cha);

		// optimize resources
		Thread.yield();
		// run the analysers
		pa.process(cas, imageGeneratorSettings);
		ca.process(cas, imageGeneratorSettings);
		cha.process(cas, imageGeneratorSettings);

		LOGGER.info("MeshReasoning completed. Took: " + PrintUtil.prettyMillis(System.currentTimeMillis() - start));
		if (imageGeneratorSettings != null && imageGeneratorSettings.isCloseAfterFinish()) {
			LOGGER.debug("Closing ...");
			System.exit(123);
		}
	}

	/**
	 * Clears all highlighted annotations in mesh reasoning view
	 */
	public void clearHightlight() {
		if (mrv != null)
			mrv.clearSelectedAnnotations();
	}

	/**
	 * Gets all annotations of given class.
	 * 
	 * @param clazz
	 *            Class of desired annotation type.
	 * @return HashSet of annotations found in model.
	 */
	@SuppressWarnings("rawtypes")
	public <T extends MeshAnnotation> HashSet<T> findAnnotations(Class<T> clazz) {
		return cas.findAnnotations(clazz);
	}

	/**
	 * Gets all the cone annotations of the CAD model
	 * 
	 * @return Set of cone annotations in model
	 */
	public HashSet<ConeAnnotation> findAnnotationsCone() {
		return cas.findAnnotations(ConeAnnotation.class);
	}

	/**
	 * Gets all container annotations of the CAD model
	 * 
	 * @return Set of container annotations in model
	 */
	public HashSet<ContainerAnnotation> findAnnotationsContainer() {
		return cas.findAnnotations(ContainerAnnotation.class);
	}

	/**
	 * Gets all plane annotations of the CAD model
	 * 
	 * @return Set of plane annotations in model
	 */
	public HashSet<PlaneAnnotation> findAnnotationsPlane() {
		return cas.findAnnotations(PlaneAnnotation.class);
	}

	/**
	 * Gets all sphere annotations of the CAD model
	 * 
	 * @return Set of sphere annotations in model
	 */
	public HashSet<SphereAnnotation> findAnnotationsSphere() {
		return cas.findAnnotations(SphereAnnotation.class);
	}

	/**
	 * Gets list of all found annotation types in mesh object
	 * 
	 * @return List of annotation names such as Sphere, Cone, Plane, ...
	 */
	public ArrayList<String> getAnnotationTypes() {
		HashSet<String> types = new HashSet<String>();
		for (Annotation a : cas.getAnnotations()) {
			if (!(a instanceof MeshAnnotation))
				continue;
			@SuppressWarnings("rawtypes")
			MeshAnnotation ma = (MeshAnnotation) a;
			String cl = ma.getClass().getSimpleName();
			if (cl.length() < 1)
				continue;
			int pos = cl.indexOf("Annotation");
			if (pos > 0)
				cl = cl.substring(0, pos);
			types.add(cl);
		}
		ArrayList<String> ret = new ArrayList<String>();
		ret.addAll(types);
		return ret;
	}

	/**
	 * Gets list of all handle annotations in model. List is ordered using HandleComparator.
	 * 
	 * @return Ordered list by handle probability of all handles found in model.
	 * 
	 * @see org.knowrob.vis.model.util.HandleComparator
	 */
	public HandleAnnotation[] getHandle() {
		return getHandle(-1f, -1f, -1f, -1f);
	}

	/**
	 * Gets list of all handle annotations in model. List is ordered using HandleComparator.
	 * 
	 * @param minRadius
	 *            handle minimum radius (in meters)
	 * @param maxRadius
	 *            handle maximum radius (in meters)
	 * 
	 * @return Ordered list by handle probability of all handles found in model.
	 * 
	 * @see org.knowrob.vis.model.util.HandleComparator
	 */
	public HandleAnnotation[] getHandle(double minRadius, double maxRadius) {
		return getHandle(minRadius, maxRadius, -1f, -1f);
	}

	/**
	 * Gets list of all handle annotations in model. List is ordered using HandleComparator.
	 * 
	 * @param minRadius
	 *            handle minimum radius (in meters)
	 * @param maxRadius
	 *            handle maximum radius (in meters)
	 * @param minLength
	 *            handle minimum length (in meters)
	 * @param maxLength
	 *            handle maximum length (in meters)
	 * 
	 * @return Ordered list by handle probability of all handles found in model.
	 * 
	 * @see org.knowrob.vis.model.util.HandleComparator
	 */
	public HandleAnnotation[] getHandle(double minRadius, double maxRadius, double minLength,
			double maxLength) {

		Set<ConeAnnotation> cones = findAnnotationsCone();
		Set<ComplexHandleAnnotation> complexHandle = cas.findAnnotations(ComplexHandleAnnotation.class);
				
		List<HandleAnnotation> allAnnotations = new ArrayList<HandleAnnotation>(cones.size() + complexHandle.size());
		
		for (ConeAnnotation c : cones) {
			if (c.isConcave())
				continue;
			allAnnotations.add(c);
		}
		allAnnotations.addAll(complexHandle);
		cones.clear();
		complexHandle.clear();
		
		HandleAnnotation[] handleList = allAnnotations.toArray(new HandleAnnotation[0]);
		Arrays.sort(handleList, new HandleComparator(cas.getModel(), minRadius, maxRadius, minLength, maxLength));
		
		ArrayList<HandleAnnotation> rets = new ArrayList<HandleAnnotation>();
		for (HandleAnnotation h : handleList) {
			
			if((HandleComparator.getHandleWeight(h, cas.getModel(),
					HandleComparator.DEFAULT_RADIUS_MIN, HandleComparator.DEFAULT_RADIUS_MAX,
					HandleComparator.DEFAULT_LENGTH_MIN, HandleComparator.DEFAULT_LENGTH_MAX) < HandleComparator.MIN_WEIGHT_FOR_HANDLE)) {
				continue;
			}

			if(cas.getModel().getUnscaled(h.getCone().getRadiusAvg()) > HandleComparator.DEFAULT_RADIUS_MAX) {
//				System.out.println("Over max radius: " + cas.getModel().getUnscaled(h.getCone().getRadiusAvg()) + " > " + HandleComparator.DEFAULT_RADIUS_MAX);
				continue;
			}
			
			if(cas.getModel().getUnscaled(h.getCone().getRadiusAvg()) < HandleComparator.DEFAULT_RADIUS_MIN) {
//				System.out.println("Under min radius: " + cas.getModel().getUnscaled(h.getCone().getRadiusAvg()) + " < " + HandleComparator.DEFAULT_RADIUS_MIN);
				continue;
			}

			rets.add(h);
		}
		return rets.toArray(new HandleAnnotation[0]);
	}

	/**
	 * Gets all triangles of the model.
	 * 
	 * @return Array of triangles
	 */
	public Triangle[] getTriangles() {
		return cas.getTriangles();
	}

	/**
	 * Gets all vertices of the model
	 * 
	 * @return Array of vertices
	 */
	public Vertex[] getVertices() {
		return cas.getVertices();
	}

	/**
	 * Highlight specified annotation in mesh reasoning view
	 * 
	 * @param a
	 *            Annotation to highlight
	 */
	public void highlightAnnotation(DrawableAnnotation a) {
		highlightAnnotation(a, null);
	}

	/**
	 * Highlight specified annotation in mesh reasoning view
	 * 
	 * @param a
	 *            Annotation to highlight
	 * @param color
	 *            highlight color. Set to null to use default color.
	 */
	public void highlightAnnotation(DrawableAnnotation a, Color color) {
		if (mrv == null)
			return;
		mrv.addSelectedAnnotation(a, color);
	}

	/**
	 * Initializes imageGeneratorSettings with predefined settings.
	 */
	private void initImageGenerator() {
		imageGeneratorSettings.setDrawAxis(false);
		imageGeneratorSettings.setWhiteBackground(true);
		// If there is a view settings file, load those settings
		imageGeneratorSettings.setInitViewFromFile(true);
		// If no view settings exist yet, allow to align model and then continue by pressing Ctrl+S
		imageGeneratorSettings.setSaveView(true);
		imageGeneratorSettings.setSavePlainModel(true);
		imageGeneratorSettings.setSaveCurvatureColor(true);
		imageGeneratorSettings.setCloseAfterFinish(true);
		imageGeneratorSettings.addAnalyserToSave(PrimitiveAnalyser.class, "segmented");
		// don't show window
		imageGeneratorSettings.setRunBackground(true);

		/**
		 * Analyser actions begin
		 */
		final MeshReasoning mr = this;

		imageGeneratorSettings.clearAnalyserActions();
		imageGeneratorSettings.addAnalyserAction(PrimitiveAnalyser.class,
				new ImageGeneratorAction() {

					@Override
					public void trigger(ImageGeneratorSettings localSettings) {
						mr.clearHightlight();
						@SuppressWarnings("rawtypes")
						PrimitiveAnnotation[] cones = mr.findAnnotationsCone().toArray(
								new PrimitiveAnnotation[0]);
						Arrays.sort(cones, new PrimitiveAnnotationAreaComparator());
						int max = 3;
						for (int i = 0; i < max && i < cones.length; i++) {
							mr.highlightAnnotation(cones[i]);
							try {
								Thread.sleep(200);
							} catch (InterruptedException e) {// do nothing
							}
							// wait til selected
							localSettings.waitSaved("cones" + (i + 1));
							mr.clearHightlight();
						}
					}
				});
		imageGeneratorSettings.addAnalyserAction(PrimitiveAnalyser.class,
				new ImageGeneratorAction() {

					@Override
					public void trigger(ImageGeneratorSettings localSettings) {
						mr.clearHightlight();
						@SuppressWarnings("rawtypes")
						PrimitiveAnnotation[] planes = mr.findAnnotationsPlane().toArray(
								new PrimitiveAnnotation[0]);
						Arrays.sort(planes, new PrimitiveAnnotationAreaComparator());
						int max = 3;
						for (int i = 0; i < max && i < planes.length; i++) {
							mr.highlightAnnotation(planes[i]);
							try {
								Thread.sleep(200);
							} catch (InterruptedException e) {// do nothing
							} // wait til selected
							localSettings.waitSaved("planes" + (i + 1));
							mr.clearHightlight();
						}
					}
				});
		imageGeneratorSettings.addAnalyserAction(PrimitiveAnalyser.class,
				new ImageGeneratorAction() {

					@Override
					public void trigger(ImageGeneratorSettings localSettings) {
						mr.clearHightlight();
						@SuppressWarnings("rawtypes")
						PrimitiveAnnotation[] spheres = mr.findAnnotationsSphere().toArray(
								new PrimitiveAnnotation[0]);
						Arrays.sort(spheres, new PrimitiveAnnotationAreaComparator());
						int max = 3;
						for (int i = 0; i < max && i < spheres.length; i++) {
							mr.highlightAnnotation(spheres[i]);
							try {
								Thread.sleep(200);
							} catch (InterruptedException e) {// do nothing
							} // wait til selected
							localSettings.waitSaved("spheres" + (i + 1));
							mr.clearHightlight();
						}
					}
				});

		imageGeneratorSettings.addAnalyserAction(ContainerAnalyser.class,
				new ImageGeneratorAction() {

					@Override
					public void trigger(ImageGeneratorSettings localSettings) {
						mr.clearHightlight();
						ContainerAnnotation[] container = mr.findAnnotations(
								ContainerAnnotation.class).toArray(new ContainerAnnotation[0]);
						Arrays.sort(container, new ContainerAnnotationVolumeComparator());
						int max = 3;
						for (int i = 0; i < max && i < container.length; i++) {
							mr.highlightAnnotation(container[i]);
							try {
								Thread.sleep(200);
							} catch (InterruptedException e) {// do nothing
							} // wait til selected
							localSettings.waitSaved("container" + (i + 1));
							mr.clearHightlight();
						}
					}
				});

		imageGeneratorSettings.addAnalyserAction(ComplexHandleAnalyser.class,
				new ImageGeneratorAction() {

					@Override
					public void trigger(ImageGeneratorSettings localSettings) {
						mr.clearHightlight();
						HandleAnnotation[] handles = mr.getHandle(
								HandleComparator.DEFAULT_RADIUS_MIN,
								HandleComparator.DEFAULT_RADIUS_MAX,
								HandleComparator.DEFAULT_LENGTH_MIN,
								HandleComparator.DEFAULT_LENGTH_MAX);
						int max = 3;
						for (int i = 0; i < max && i < handles.length; i++) {
							mr.highlightAnnotation((DrawableAnnotation) handles[i]);
							try {
								Thread.sleep(200);
							} catch (InterruptedException e) {// do nothing
							} // wait til selected
							localSettings.waitSaved("handle" + (i + 1));
							mr.clearHightlight();
						}
					}
				});

		/**
		 * Analyser actions end
		 */
	}

	/**
	 * Initializes image generator settings for spoon images where handle and spoon sphere are selected.
	 */
	@SuppressWarnings("unused")
	private void initImageGeneratorSpoon() {
		imageGeneratorSettings.setDrawAxis(false);
		imageGeneratorSettings.setWhiteBackground(true);
		// If there is a view settings file, load those settings
		imageGeneratorSettings.setInitViewFromFile(true);
		// If no view settings exist yet, allow to align model and then continue by pressing Ctrl+S
		imageGeneratorSettings.setSaveView(true);
		imageGeneratorSettings.setSavePlainModel(true);
		imageGeneratorSettings.setSaveCurvatureColor(true);
		imageGeneratorSettings.setCloseAfterFinish(true);
		// even draw fitted cone/sphere/plane if more than one annotation selected
		imageGeneratorSettings.setAlwaysDrawSelectedPrimitives(true);

		/**
		 * Analyser actions begin
		 */
		final MeshReasoning mr = this;

		imageGeneratorSettings.clearAnalyserActions();
		imageGeneratorSettings.addAnalyserAction(PrimitiveAnalyser.class,
				new ImageGeneratorAction() {

					@Override
					public void trigger(ImageGeneratorSettings localSettings) {
						mr.clearHightlight();
						@SuppressWarnings("rawtypes")
						PrimitiveAnnotation[] cones = mr.findAnnotationsCone().toArray(
								new PrimitiveAnnotation[0]);
						Arrays.sort(cones, new PrimitiveAnnotationAreaComparator());

						@SuppressWarnings("rawtypes")
						PrimitiveAnnotation[] spheres = mr.findAnnotationsSphere().toArray(
								new PrimitiveAnnotation[0]);
						Arrays.sort(spheres, new PrimitiveAnnotationAreaComparator());

						if (cones.length == 0) {
							LOGGER.warn("Skipping image because no cone found");
							return;
						}
						if (spheres.length == 0) {
							LOGGER.warn("Skipping image because no sphere found");
							return;
						}
						for (int i = 0; i < cones.length; i++) {
							if (!((ConeAnnotation) cones[i]).isConcave()) {
								mr.highlightAnnotation(cones[i]);
								break;
							}
						}
						for (int i = 0; i < spheres.length; i++) {
							if (((SphereAnnotation) spheres[i]).isConcave()) {
								mr.highlightAnnotation(spheres[i]);
								break;
							}
						}
						try {
							Thread.sleep(200);
						} catch (InterruptedException e) {// do nothing
						} // wait til selected
						localSettings.waitSaved("ConeSphere");
						mr.clearHightlight();

					}
				});

		/**
		 * Analyser actions end
		 */
	}

	/**
	 * Default image file name for saving current mesh reasoning view as a png image.
	 * 
	 * @param s
	 *            default file name
	 */
	public void setDefaultImageFilename(String s) {
		if (control != null)
			control.setDefaultImageFilename(s);
	}

	/**
	 * Sets the title of the main frame
	 * 
	 * @param title
	 *            new title
	 */
	public void setFrameTitle(String title) {
		frame.setTitle(title);
	}
}
