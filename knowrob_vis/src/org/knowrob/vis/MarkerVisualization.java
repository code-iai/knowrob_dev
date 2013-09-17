package org.knowrob.vis;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Vector;

import javax.vecmath.Matrix4d;
import javax.vecmath.Quat4d;

import ros.NodeHandle;
import ros.Publisher;
import ros.Ros;
import ros.RosException;
import ros.communication.Duration;
import ros.communication.Time;
import ros.pkg.std_msgs.msg.ColorRGBA;
import ros.pkg.visualization_msgs.msg.*;

import edu.tum.cs.ias.knowrob.owl.OWLThing;
import edu.tum.cs.ias.knowrob.prolog.PrologInterface;


/**
 * Visualization module for the KnowRob knowledge base 
 * 
 * The objects to be visualized are published as 'Markers'
 * on the visualization_marker topic and can be visualized
 * in the 'rviz' program.
 * 
 * @author tenorth@cs.uni-bremen.de
 *
 */
public class MarkerVisualization {

	static Ros ros;
	public static NodeHandle n;
	Thread markerPublisher;

	/**
	 * Store the markers to be published 
	 */
	protected Map<String, Marker> markers;


	/**
	 * Store the highlights and the original colors
	 */
	protected Map<String, ColorRGBA> highlighted;


	/**
	 * Counter for marker IDs 
	 */
	private static int id = 0;



	/**
	 * Constructor. Starts the marker publisher in a parallel thread. 
	 */
	public MarkerVisualization() {
		
		initRos();
//		PrologInterface.initJPLProlog("knowrob_vis");
		
		markers =  Collections.synchronizedMap(new HashMap<String, Marker>());
		highlighted = Collections.synchronizedMap(new HashMap<String, ColorRGBA>());

		// spawn new thread that publishes the markers in the HashMap
		markerPublisher = new Thread( new PublisherThread() );
		markerPublisher.start();
	}


	/**
	 * Add object 'identifier' to the visualization.
	 * 
	 * @param identifier OWL identifier of an object instance
	 */
	public void addObject(String identifier) {

		// read marker from Prolog
		Marker m = readMarkerFromProlog(identifier);

		// add marker to map
		synchronized (markers) {
			markers.put(identifier, m);
		}
	}



	/**
	 * Add object 'identifier' and all its parts to the visualization.
	 * 
	 * @param identifier OWL identifier of an object instance
	 */
	public void addObjectWithChildren(String identifier) {

		// add this object
		addObject(identifier);

		// read children and add them too
		for(String child : readChildren(identifier))
			addObject(child);
	}



	/**
	 * Remove object 'identifier' from the visualization.
	 * 
	 * @param identifier OWL identifier of an object instance
	 */
	public void removeObject(String identifier) {

		// remove the object from the list
		synchronized (markers) {
			markers.remove(identifier);
		}
	}

	/**
	 * Remove the object 'identifier' as well as its children
	 * 
	 * @param identifier OWL identifier of an object instance
	 */
	public void removeObjectWithChildren(String identifier) {

		// remove this object
		removeObject(identifier);

		// remove children and remove them too
		for(String child : readChildren(identifier))
			removeObject(child);
	}
	
	
	/**
	 * Remove all objects from the visualization
	 */
	public void clear() {
		synchronized (markers) {
			markers.clear();
		}
	}



	/**
	 * Highlight the object 'identifier' in red
	 * 
	 * @param identifier OWL identifier of an object instance
	 * @param highlight True to set, False to remove highlight
	 */
	public void highlight(String identifier, boolean highlight) {
		highlight(identifier, highlight, 200, 0, 0, 255);
	}

	
	
	/**
	 * Highlight the object 'identifier' in 'color'
	 * 
	 * @param identifier OWL identifier of an object instance
	 * @param highlight True to set, False to remove highlight
	 * @param color Integer of the form #AARRGGBB
	 */
	public void highlight(String identifier, boolean highlight, int color) {
		int a = color & 0xff000000 >> 24;
		int r = color & 0x00ff0000 >> 16;
		int g = color & 0x0000ff00 >> 8;
		int b = color & 0x000000ff;
		highlight(identifier, highlight, r, g, b, a);
	}
	
	
	/**
	 * Highlight the object 'identifier' with color RGBA
	 * 
	 * @param identifier OWL identifier of an object instance
	 * @param highlight True to set, False to remove highlight
	 * @param r Red value (0--255)
	 * @param g Green value (0--255)
	 * @param b Blue value (0--255)
	 * @param a Alpha value (0--255)
	 */
	public void highlight(String identifier, boolean highlight, int r, int g, int b, int a) {

		if(!highlight) {
			synchronized (markers) {				
				markers.get(identifier).color = highlighted.get(identifier);
			}
			
		} else {

			synchronized (highlighted) {
				highlighted.put(identifier, markers.get(identifier).color);
			}

			synchronized (markers) {
				markers.get(identifier).color = new ColorRGBA();
				markers.get(identifier).color.r = ((float) r)/255;
				markers.get(identifier).color.g = ((float) g)/255;
				markers.get(identifier).color.b = ((float) b)/255;
				markers.get(identifier).color.a = ((float) a)/255;
			}
		}
	}


	/**
	 * Highlight the object 'identifier' and its children in red
	 * 
	 * @param identifier OWL identifier of an object instance
	 * @param highlight True to set, False to remove highlight
	 */
	public void highlightWithChildren(String identifier, boolean highlight) {

		// remove this object
		highlight(identifier, highlight); 

		// remove children and highlight them too
		for(String child : readChildren(identifier))
			highlight(child, highlight); 
	}


	/**
	 * Clear all highlights
	 */
	public void clearHighlight() {

		// reset colors to cached original ones

		synchronized (markers) {
			for(String obj : highlighted.keySet()) {
				markers.get(obj).color = highlighted.get(obj);
			}
		}
	}



	// // // // // // // // // // // // // // // // // // // // // // // // // // // 
	//  
	// Helper methods: read data from Prolog, create data structures
	//  

	/**
	 * Thread-safe ROS initialization
	 */
	protected static void initRos() {

		ros = Ros.getInstance();

		if(!Ros.getInstance().isInitialized()) {
			ros.init("knowrob_vis");
		}
		n = ros.createNodeHandle();

	}


	/**
	 * Read children of an object instance that are either linked 
	 * by the 'properPhysicalParts' property or (inversely) by the
	 * 'describedInMap' property.
	 * 
	 * @param parent The parent whose children are to the returned, as OWL identifier
	 * @return Array of OWL identifiers of all children
	 */
	private String[] readChildren(String parent) {

		HashSet<String> children = new HashSet<String>();

		// if parent is map, read map objects; otherwise, the rdf_reachable returns the parent itself

		HashMap<String, Vector<String>> mapParts = PrologInterface.executeQuery(
				"rdf_reachable(PART, knowrob:describedInMap, '"+parent+"')");

		if(mapParts != null && mapParts.get("PART") != null) {
			for(int i=0;i<mapParts.get("PART").size();i++) {
				if(!mapParts.get("PART").get(i).equals(parent)) {

					// add object to children set
					children.add(OWLThing.removeSingleQuotes(mapParts.get("PART").get(i)));


					// read all physical parts of all child objects
					HashMap<String, Vector<String>> parts = PrologInterface.executeQuery(
							"rdf_reachable("+mapParts.get("PART").get(i)+", knowrob:properPhysicalParts, P)");

					if(parts != null && parts.get("P") != null) {

						for(int j=0;j<parts.get("P").size();j++) {
							if(!parts.get("P").get(j).toString().equals(mapParts.get("PART").get(i))) {
								children.add(OWLThing.removeSingleQuotes(parts.get("P").get(j)));
							}
						}
					}
				}
			}
		}
		return children.toArray(new String[]{});
	}


	/**
	 * Read object information from Prolog and create a marker from it
	 * 
	 * @param identifier
	 * @return
	 */
	Marker readMarkerFromProlog(String identifier) {

		Marker m = new Marker();

		m.header.frame_id = "/map";
		m.header.stamp = Time.now();
		m.ns = "knowrob_vis";
		m.id = id++;

		m.action = Marker.ADD;
		m.lifetime = new Duration();

		try {
			// read object pose
			HashMap<String, Vector<String>> res = PrologInterface.executeQuery(
					"current_object_pose('"+ identifier + "', [M00, M01, M02, M03, M10, M11, M12, M13, M20, M21, M22, M23, M30, M31, M32, M33])");

			if (res!=null && res.get("M00") != null && res.get("M00").size() > 0 && res.get("M00").get(0)!=null) {

				double[] p = new double[16];
				Matrix4d poseMat = new Matrix4d(p);

				for(int i=0;i<4;i++) {
					for(int j=0;j<4;j++) {
						poseMat.setElement(i, j, Double.valueOf(res.get("M"+i+j).get(0)));
					}
				}

				Quat4d q = new Quat4d();
				q.set(poseMat);

				m.pose.orientation.w = q.w;
				m.pose.orientation.x = q.x;
				m.pose.orientation.y = q.y;
				m.pose.orientation.z = q.z;

				m.pose.position.x = poseMat.m03;
				m.pose.position.y = poseMat.m13;
				m.pose.position.z = poseMat.m23;

			} else {
				m.type = Marker.CUBE;
			}

		} catch (Exception e) {
			e.printStackTrace();
		}


		try {

			// read object dimensions if available, default to a 10cm cube
			HashMap<String, Vector<String>> res = PrologInterface.executeQuery(
					"rdf_has('"+identifier+"', knowrob:depthOfObject,  literal(type(_, D))), " +
					"rdf_has('"+identifier+"', knowrob:widthOfObject,  literal(type(_, W))), " +
					"rdf_has('"+identifier+"', knowrob:heightOfObject, literal(type(_, H)))");

			if (res!=null && res.get("D") != null && res.get("D").size() > 0 && res.get("D").get(0)!=null) {
				m.scale.x = Double.valueOf(OWLThing.removeSingleQuotes(res.get("D").get(0)));
				m.scale.y = Double.valueOf(OWLThing.removeSingleQuotes(res.get("W").get(0)));
				m.scale.z = Double.valueOf(OWLThing.removeSingleQuotes(res.get("H").get(0)));

			} else {
				m.scale.x = 0.1;
				m.scale.y = 0.1;
				m.scale.z = 0.1;
			}


		} catch (Exception e) {
			e.printStackTrace();
		}


		try {

			// check if mesh is available for this object
			HashMap<String, Vector<String>> res = PrologInterface.executeQuery("get_model_path('"+ identifier + "', Path)");

			if (res!=null && res.get("Path") != null && res.get("Path").size() > 0 && res.get("Path").get(0)!=null) {

				m.type = Marker.MESH_RESOURCE;
				m.mesh_resource = OWLThing.removeSingleQuotes(res.get("Path").get(0));

			} else {
				m.type = Marker.CUBE;
			}

			// set light grey color by default
			m.color.r = 0.6f;
			m.color.g = 0.6f;
			m.color.b = 0.6f;
			m.color.a = 1.0f;

		} catch (Exception e) {
			e.printStackTrace();
		}

		return m;
	}


	/**
	 * Thread that publishes the current state of the marker set
	 * to the visualization_marker topic.
	 * 
	 * @author tenorth@cs.uni-bremen.de
	 *
	 */
	public class PublisherThread implements Runnable {

		@Override 
		public void run() {

			try {
				Publisher<Marker> pub = n.advertise("visualization_marker", new Marker(), 100);

				while(n.isValid()) {

					synchronized (markers) {
						for(String mrk : markers.keySet())
							pub.publish(markers.get(mrk));
					}
					n.spinOnce();
					Thread.sleep(500);
				}

				pub.shutdown();

			} catch (RosException e) {
				e.printStackTrace();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}

	
	public static void main(String args[]) {

		MarkerVisualization vis = new MarkerVisualization();
		vis.addObjectWithChildren("http://ias.cs.tum.edu/kb/ias_semantic_map.owl#SemanticEnvironmentMap0");
		vis.highlight("http://ias.cs.tum.edu/kb/knowrob.owl#Refrigerator67", true, 150, 0, 0, 180);

	}
}
