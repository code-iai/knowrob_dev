/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: 
 * 				Stefan Profanter - initial API and implementation, Year: 2012
 * 				Andrei Stoica - improved mesh reasoning during Google Summer of Code 2014
 ******************************************************************************/
package org.knowrob.vis.test;

import org.knowrob.vis.model.MeshReasoning;

/**
 * Main test class for testing MeshReasoningTest with visualization
 * of the model under test and its reasoning process.
 * 
 * @author Stefan Profanter
 * @author Andrei Stoica - refactored main method called analyseByPath(String), {@link edu.tum.cs.vis.model.MeshReasoning#analyseByPath(String)}
 */
public class MeshReasoningTest {

	/**
	 * Main Method loading the mesh, drawing it and starting the analyzer
	 * for only one model. Modify path as needed.
	 * 
	 * @param args
	 *            command line arguments not used
	 */
	public static void main(String[] args) {
		// modify path as desired
		String path = "package://knowrob_cad_models/models/kitchen/drinking-vessels/cup2.dae";
		
		MeshReasoning mr = MeshReasoning.initMeshReasoning(true);

		mr.analyseByPath(path);
	}
}
