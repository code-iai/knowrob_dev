/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Stefan Profanter - initial API and implementation, Year: 2012
 ******************************************************************************/
package edu.tum.cs.vis.model;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.util.ArrayList;

import javax.swing.JFrame;

import org.apache.log4j.Logger;
import org.apache.log4j.xml.DOMConfigurator;

import edu.tum.cs.util.PrintUtil;
import edu.tum.cs.vis.model.uima.analyzer.ContainerAnalyzer;
import edu.tum.cs.vis.model.uima.analyzer.MeshAnalyzer;
import edu.tum.cs.vis.model.uima.analyzer.NeighborAnalyzer;
import edu.tum.cs.vis.model.uima.analyzer.PrimitiveAnalyzer;
import edu.tum.cs.vis.model.uima.cas.MeshCas;
import edu.tum.cs.vis.model.util.algorithm.CurvatureCalculation;
import edu.tum.cs.vis.model.view.MeshReasoningView;
import edu.tum.cs.vis.model.view.MeshReasoningViewControl;

/**
 * @author Stefan Profanter
 * 
 */
public class MeshReasoning {

	/**
	 * Log4j logger
	 */
	private static Logger	logger	= Logger.getRootLogger();

	public static MeshReasoning initMeshReasoning(boolean withView) {

		DOMConfigurator.configureAndWatch("log4j.xml", 60 * 1000);
		return new MeshReasoning(withView);
	}

	private MeshReasoningView	mrv	= null;

	public MeshReasoning(boolean withView) {
		MeshCas cas = new MeshCas();

		if (withView) {
			JFrame frame = new JFrame();
			frame.setMinimumSize(new Dimension(800, 600));
			frame.setSize(1324, 768);
			frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			frame.setTitle("Mesh reasoning view");
			frame.setLocationRelativeTo(null);

			ArrayList<MeshAnalyzer> analyzer = new ArrayList<MeshAnalyzer>();

			mrv = new MeshReasoningView();
			MeshReasoningViewControl control = new MeshReasoningViewControl(cas, analyzer, mrv);
			mrv.setControl(control);
			mrv.init();

			mrv.getCasList().add(cas);

			frame.setLayout(new BorderLayout());
			frame.getContentPane().add(mrv, BorderLayout.CENTER);
			frame.getContentPane().add(control, BorderLayout.LINE_END);
			frame.setVisible(true);
		}

	}

	public MeshCas analyzeByIdentifier(String prologIdentifier) {
		String path = Properties.getPropertyStringOfClassOrIndividual("knowrob:pathToCadModel",
				prologIdentifier);
		if (path != null)
			return analyzeByPath(path);
		return null;
	}

	public MeshCas analyzeByPath(String path) {

		logger.info("MeshReasoning started. Parsing model ...");
		long start = System.currentTimeMillis();

		ItemModel itemModel = new ItemModel(path);

		if (itemModel.getParser() == null) {
			logger.error("Couldn't parse model. Maybe path of model file is wrong.");
			return null;
		}

		Model model = itemModel.getParser().getModel();
		logger.debug("Model parsed. Took: "
				+ PrintUtil.prettyMillis(System.currentTimeMillis() - start) + " (Vertices: "
				+ model.getVertices().size() + ", Lines: " + model.getLines().size()
				+ ", Triangles: " + model.getTriangles().size() + ")");

		logger.debug("Calculating curvature ...");

		model.normalize();

		MeshCas cas;
		ArrayList<MeshAnalyzer> analyzer;
		if (mrv != null) {
			cas = mrv.getControl().getCas();
			analyzer = mrv.getControl().getAnalyzer();
		} else {
			cas = new MeshCas();
			analyzer = new ArrayList<MeshAnalyzer>();
		}
		cas.setModel(model);
		CurvatureCalculation.calculateCurvatures(cas.getCurvatures(), model);

		NeighborAnalyzer na = new NeighborAnalyzer();
		analyzer.add(na);
		PrimitiveAnalyzer pa = new PrimitiveAnalyzer();
		analyzer.add(pa);
		ContainerAnalyzer ca = new ContainerAnalyzer();
		analyzer.add(ca);

		Thread.yield();

		na.process(cas);
		pa.process(cas);
		ca.process(cas);

		return cas;
	}

}
