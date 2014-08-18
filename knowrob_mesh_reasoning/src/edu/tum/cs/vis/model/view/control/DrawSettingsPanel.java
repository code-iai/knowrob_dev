/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Stefan Profanter - initial API and implementation, Year: 2012
 ******************************************************************************/
package edu.tum.cs.vis.model.view.control;

import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import edu.tum.cs.vis.model.uima.cas.MeshCas;
import edu.tum.cs.vis.model.view.MeshReasoningView;

/**
 * Control panel for draw settings of the model.
 * Implements the buttons frontend which is featured in
 * the main applet of the mesh reasoning view.
 * 
 * @author Stefan Profanter
 * @author Andrei Stoica - implemented new buttons for an extended view
 */
public class DrawSettingsPanel extends JPanel implements ActionListener {

	/**
	 * auto generated
	 */
	private static final long		serialVersionUID	= -1936392151575823886L;

	/**
	 * cas for which this control is
	 */
	@SuppressWarnings("unused")
	private final MeshCas			cas;

	/**
	 * The view on which this control is used
	 */
	private final MeshReasoningView	view;

	/**
	 * Per vertex normals drawing check box button
	 */
	private final JCheckBox			cbxDrawVertexNormals;
	
	/**
	 * Per triangle normals drawing check box button
	 */
	private final JCheckBox			cbxDrawTriangleNormals;
	
	/**
	 * Per vertex normalized minimum curvature direction drawing check box button
	 */
	private final JCheckBox			cbxDrawVertexCurvatureMin;
	
	/**
	 * Per vertex normalized maximum curvature direction drawing check box button
	 */
	private final JCheckBox			cbxDrawVertexCurvatureMax;
	
	/**
	 * Per vertex scaled minimum and maximum curvature directions drawing check box button 
	 */
	private final JCheckBox			cbxDrawVertexCurvature;
	
	/**
	 * Model curvature coloring scheme drawing combo box button 
	 */
	@SuppressWarnings("rawtypes")
	private final JComboBox			cmbbxDrawCurvatureColor;
	
	/**
	 * Model curvature coloring scheme drawing combo box string options array
	 */
	private final String[]			cmbbxOptions;
	
	/**
	 * Region edges drawing check box button
	 */
	private final JCheckBox			cbxDrawRegionEdges;
	
	/**
	 * Sharp edges drawing check box button
	 */
	private final JCheckBox			cbxDrawSharpEdges;
	
	/**
	 * Select only triangles drawing type check box button
	 */
	private final JCheckBox			cbxSelectTriangles;
	
	/**
	 * Select only the nearest triangle drawing type check box button
	 */
	private final JCheckBox			cbxSelectNearestOnly;
	
	/**
	 * Per model group bounding box drawing check box button
	 */
	private final JCheckBox			cbxDrawBoundingBox;

	/**
	 * Per vertex Voronoi area sphere drawing check box button
	 */
	private final JCheckBox			cbxDrawVoronoiArea;
	
	/**
	 * White background check button
	 */
	private final JCheckBox			cbxWhiteBackground;
	
	/**
	 * Set view perspective (rotate camera to given parameters) check box button
	 */
	private final JButton			btnSetRotation;
	
	/**
	 * Creates a new draw settings panel
	 * 
	 * @param cas
	 *            main mesh cas
	 * @param view
	 *            parent mesh reasoning view
	 */
	@SuppressWarnings({ "rawtypes", "unchecked" })
	public DrawSettingsPanel(MeshCas cas, MeshReasoningView view) {
		this.view = view;
		this.cas = cas;
		// setPreferredSize(new Dimension(300, 300));

		GridLayout grid = new GridLayout(0, 2);
		setLayout(grid);
		cmbbxOptions = new String[4];
		cmbbxOptions[0] = new String("No curv. coloring");
		cmbbxOptions[1] = new String("Mean curv. coloring");
		cmbbxOptions[2] = new String("Gauss curv. coloring");
		cmbbxOptions[3] = new String("Est. curv. coloring");

		cbxDrawVertexNormals = new JCheckBox("Vertex normals");
		cbxDrawVertexNormals.addActionListener(this);
		cbxDrawVertexNormals.setSelected(false);
		this.add(cbxDrawVertexNormals);

		cbxDrawTriangleNormals = new JCheckBox("Triangle normals");
		cbxDrawTriangleNormals.addActionListener(this);
		cbxDrawTriangleNormals.setSelected(false);
		this.add(cbxDrawTriangleNormals);
		
		cbxDrawVertexCurvatureMin = new JCheckBox("Vertex min. curvature");
		cbxDrawVertexCurvatureMin.addActionListener(this);
		cbxDrawVertexCurvatureMin.setSelected(false);
		this.add(cbxDrawVertexCurvatureMin);
		
		cbxDrawVertexCurvatureMax = new JCheckBox("Vertex max. curvature");
		cbxDrawVertexCurvatureMax.addActionListener(this);
		cbxDrawVertexCurvatureMax.setSelected(false);
		this.add(cbxDrawVertexCurvatureMax);

		cbxDrawVertexCurvature = new JCheckBox("Vertex curvature");
		cbxDrawVertexCurvature.addActionListener(this);
		cbxDrawVertexCurvature.setSelected(false);
		this.add(cbxDrawVertexCurvature);
		
		// parameterization skipped because compatibility with Java JRE >1.6 was desired
		cmbbxDrawCurvatureColor = new JComboBox(cmbbxOptions);
		cmbbxDrawCurvatureColor.addActionListener(this);
		cmbbxDrawCurvatureColor.setSelectedIndex(0);
		this.add(cmbbxDrawCurvatureColor);

		cbxDrawRegionEdges = new JCheckBox("Region Edges");
		cbxDrawRegionEdges.addActionListener(this);
		cbxDrawRegionEdges.setSelected(false);
		this.add(cbxDrawRegionEdges);
		
		cbxDrawSharpEdges = new JCheckBox("Sharp Edges");
		cbxDrawSharpEdges.addActionListener(this);
		cbxDrawSharpEdges.setSelected(false);
		this.add(cbxDrawSharpEdges);
		
		cbxSelectTriangles = new JCheckBox("Select triangles");
		cbxSelectTriangles.addActionListener(this);
		cbxSelectTriangles.setSelected(view.isSelectTrianglesOnly());
		this.add(cbxSelectTriangles);
		
		cbxSelectNearestOnly = new JCheckBox("Select nearest");
		cbxSelectNearestOnly.addActionListener(this);
		cbxSelectNearestOnly.setSelected(view.isSelectNearestOnly());
		this.add(cbxSelectNearestOnly);

		cbxDrawVoronoiArea = new JCheckBox("Voronoi Area");
		cbxDrawVoronoiArea.addActionListener(this);
		cbxDrawVoronoiArea.setSelected(false);
		this.add(cbxDrawVoronoiArea);
		
		cbxDrawBoundingBox = new JCheckBox("Bounding box");
		cbxDrawBoundingBox.addActionListener(this);
		cbxDrawBoundingBox.setSelected(view.isDrawBoundingBox());
		this.add(cbxDrawBoundingBox);

		cbxWhiteBackground = new JCheckBox("White Background");
		cbxWhiteBackground.addActionListener(this);
		cbxWhiteBackground.setSelected(view.isBackgroundWhite());
		this.add(cbxWhiteBackground);

		btnSetRotation = new JButton("Set View");
		btnSetRotation.addActionListener(this);
		this.add(btnSetRotation);

	}

	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == cbxWhiteBackground) {
			view.setBackgroundWhite(cbxWhiteBackground.isSelected());
		}
		else if (e.getSource() == cbxDrawVertexNormals) {
			view.setDrawVertexNormals(cbxDrawVertexNormals.isSelected());
		}
		else if (e.getSource() == cbxDrawTriangleNormals) {
			view.setDrawTriangleNormals(cbxDrawTriangleNormals.isSelected());
		}
		else if (e.getSource() == cbxDrawVertexCurvatureMin) {
			view.setDrawVertexCurvatureMin(cbxDrawVertexCurvatureMin.isSelected()); 
		}
		else if (e.getSource() == cbxDrawVertexCurvatureMax) {
			view.setDrawVertexCurvatureMax(cbxDrawVertexCurvatureMax.isSelected());
		}
		else if (e.getSource() == cbxDrawVertexCurvature) {
			view.setDrawVertexCurvature(cbxDrawVertexCurvature.isSelected());
		}
		else if (e.getSource() == cmbbxDrawCurvatureColor) {
			if (cmbbxDrawCurvatureColor.getSelectedIndex() == 0) {
				view.setDrawCurvatureColor(0);
			}
			else if (cmbbxDrawCurvatureColor.getSelectedIndex() == 1) {
				view.setDrawCurvatureColor(1);
			}
			else if (cmbbxDrawCurvatureColor.getSelectedIndex() == 2) {
				view.setDrawCurvatureColor(2);
			}
			else if (cmbbxDrawCurvatureColor.getSelectedIndex() == 3) {
				view.setDrawCurvatureColor(3);
			}
		}
		else if (e.getSource() == cbxDrawVoronoiArea) {
			view.setDrawVoronoiArea(cbxDrawVoronoiArea.isSelected());
		}
		else if (e.getSource() == cbxDrawSharpEdges) {
			view.setDrawSharpEdges(cbxDrawSharpEdges.isSelected());
		}
		else if (e.getSource() == cbxDrawRegionEdges) {
			view.setDrawRegionEdges(cbxDrawRegionEdges.isSelected());
		}
		else if (e.getSource() == cbxDrawBoundingBox) {
			view.setDrawBoundingBox(cbxDrawBoundingBox.isSelected());
		}
		else if (e.getSource() == btnSetRotation) {
			String current = Math.round(view.getCam().getRotations()[0] * 180f / Math.PI) + ","
					+ Math.round(view.getCam().getRotations()[1] * 180f / Math.PI) + ","
					+ Math.round(view.getCam().getRotations()[2] * 180f / Math.PI);

			String s = (String) JOptionPane.showInputDialog(this,
					"Enter the view parameters (in degree) in the format 'pitch,yaw,roll'",
					"View parameters", JOptionPane.PLAIN_MESSAGE, null, null, current);

			// If a string was returned, say so.
			if ((s != null) && (s.length() > 0)) {
				String parts[] = s.split(",");
				if (parts.length != 3) {
					JOptionPane.showMessageDialog(this, "Invalid format.", "Input error",
							JOptionPane.ERROR_MESSAGE);
					return;
				}

				try {
					int pitch = Integer.parseInt(parts[0]);
					int yaw = Integer.parseInt(parts[1]);
					int roll = Integer.parseInt(parts[2]);
					view.setManualRotation((float) (pitch * Math.PI / 180f),
							(float) (yaw * Math.PI / 180f), (float) (roll * Math.PI / 180f));
				} catch (NumberFormatException ex) {
					JOptionPane.showMessageDialog(this, "Invalid format.", "Input error",
							JOptionPane.ERROR_MESSAGE);
					return;
				}
			}
		} else if (e.getSource() == cbxSelectNearestOnly)
			view.setSelectNearestOnly(cbxSelectNearestOnly.isSelected());
		else if (e.getSource() == cbxSelectTriangles)
			view.setSelectTrianglesOnly(cbxSelectTriangles.isSelected());
	}
}
