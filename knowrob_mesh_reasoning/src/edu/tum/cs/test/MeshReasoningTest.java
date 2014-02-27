/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Stefan Profanter - initial API and implementation, Year: 2012
 ******************************************************************************/
package edu.tum.cs.test;

import edu.tum.cs.vis.model.MeshReasoning;

/**
 * Main test class for testing MeshReasoningTest
 * 
 * @author Stefan Profanter
 * 
 */
public class MeshReasoningTest {

	/**
	 * Main Method loading the mesh, drawing it and starting the analyzer
	 * 
	 * @param args
	 *            command line arguments not used
	 */
	public static void main(String[] args) {

//		String path = "/home/tenorth/work/ros/groovy/rosbuild_ws/stacks/knowrob_addons/knowrob_cad_models/models/drinking-vessels/cup2.dae";
//		String path = "/home/tenorth/work/ros/groovy/rosbuild_ws/stacks/knowrob_addons/knowrob_cad_models/models/kitchen-tools/spatula.kmz";
		String path = "/home/tenorth/work/ros/groovy/rosbuild_ws/stacks/knowrob_addons/knowrob_cad_models/models/electric-devices/pancake-maker.dae";
//		String path = "/home/tenorth/work/ros/groovy/rosbuild_ws/stacks/knowrob_addons/knowrob_cad_models/models/food-drinks/mondamin-pancake-mix.dae";
		

		MeshReasoning mr = MeshReasoning.initMeshReasoning(true);

		mr.analyseByPath(path);
	}
}
