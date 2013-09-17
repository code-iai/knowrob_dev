package org.knowrob.vis;

import java.util.ArrayList;
import java.util.HashSet;

import javax.vecmath.Matrix4d;
import javax.vecmath.Quat4d;

import ros.NodeHandle;
import ros.Publisher;
import ros.Ros;
import ros.RosException;
import ros.communication.Duration;
import ros.communication.Time;
import ros.pkg.visualization_msgs.msg.*;

import edu.tum.cs.ias.knowrob.json_prolog.Prolog;
import edu.tum.cs.ias.knowrob.json_prolog.PrologBindings;
import edu.tum.cs.ias.knowrob.json_prolog.PrologQueryProxy;


public class MarkerVisualization {

	static Ros ros;
	public static NodeHandle n;


	public MarkerVisualization() {

		initRos();

	}



	public void publishMap(String map_id) {


		Prolog pl = new Prolog();
		PrologQueryProxy bdgs = pl.query("map_root_objects('"+map_id+"', Objs), member(Obj, Objs), map_object_info([Obj, Type, [M00, M01, M02, M03, M10, M11, M12, M13, M20, M21, M22, M23, M30, M31, M32, M33], [D, W, H]])");
		//PrologQueryProxy bdgs = pl.query("map_root_objects('"+map_id+"', Objs), member(Obj, Objs), map_object_info([Obj, Type, Pose, [D, W, H]])");

		ArrayList<Marker> markers = new ArrayList<Marker>();
		HashSet<String> objs = new HashSet<String>();

		int id = 0;
		for(PrologBindings bdg : bdgs) {

			String obj =  bdg.getBdgs_().get("Obj").getValue().toString();

			if(objs.contains(obj)) {
				continue;
			} else {
				objs.add(obj);
			}

			Marker m = new Marker();

			m.header.frame_id = "/map";
			m.header.stamp = Time.now();
			m.ns = "basic_shapes";
			m.id = id++;
			
			m.text = obj;
			
			m.action = Marker.ADD;
			m.lifetime = new Duration();

			m.type = Marker.CUBE;
			m.scale.x = Double.valueOf(bdg.getBdgs_().get("D").getValue().toString());
			m.scale.y = Double.valueOf(bdg.getBdgs_().get("W").getValue().toString());
			m.scale.z = Double.valueOf(bdg.getBdgs_().get("H").getValue().toString());


			double[] p = new double[16];
			Matrix4d poseMat = new Matrix4d(p);

			for(int i=0;i<4;i++) {
				for(int j=0;j<4;j++) {
					poseMat.setElement(i, j, Double.valueOf(bdg.getBdgs_().get("M"+i+j).getValue().toString()));
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
			
			m.color.r = 0.6f;
			m.color.g = 0.6f;
			m.color.b = 0.6f;
			m.color.a = 1.0f;

			
			markers.add(m);

		}


		Marker m = new Marker();

		m.header.frame_id = "/map";
		m.header.stamp = Time.now();
		m.ns = "basic_shapes";
		m.id = id++;
		
		m.action = Marker.ADD;
		m.lifetime = new Duration();

		m.type = Marker.MESH_RESOURCE;
		m.scale.x = 1;
		m.scale.y = 1;
		m.scale.z = 1;
		
		m.color.r = 0.6f;
		m.color.g = 0.6f;
		m.color.b = 0.6f;
		m.color.a = 1.0f;

//		m.mesh_resource = "package://knowrob_cad_models/models/electric-devices/iron_2.dae";
		m.mesh_resource = "package://pr2_description/meshes/base_v0/base.dae";
		
		markers.add(m);
		
		


		Publisher<Marker> pub;
		try {
			pub = n.advertise("visualization_marker", new Marker(), 100);

			while(n.isValid()) {

				for(Marker mrk : markers)
					pub.publish(mrk);
				
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



	/**
	 * Thread-safe ROS initialization
	 */
	protected static void initRos() {

		ros = Ros.getInstance();

		if(!Ros.getInstance().isInitialized()) {
			ros.init("knowrob_flowchart_vis");
		}
		n = ros.createNodeHandle();

	}

	public static void main(String args[]) {

		MarkerVisualization vis = new MarkerVisualization();
		vis.publishMap("http://ias.cs.tum.edu/kb/ias_semantic_map.owl#SemanticEnvironmentMap0");

	}
}
