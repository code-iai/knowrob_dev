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
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import edu.tum.cs.vis.model.uima.cas.MeshCas;
import edu.tum.cs.vis.model.view.MeshReasoningView;

/**
 * Control panel for draw settings
 * 
 * @author Stefan Profanter
 * 
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
	 * The view for which this control is
	 */
	private final MeshReasoningView	view;

	/**
	 * draw vertex normals
	 */
	private final JCheckBox			cbxDrawVertexNormals;
	
	/**
	 * draw triangle normals
	 */
	private final JCheckBox			cbxDrawTriangleNormals;
	
	/**
	 * draw vertex normalized minimum curvature direction
	 */
	private final JCheckBox			cbxDrawVertexCurvatureMin;
	
	/**
	 * draw vertex normalized maximum curvature direction
	 */
	private final JCheckBox			cbxDrawVertexCurvatureMax;
	
	/**
	 * draw vertex together both min and max curvature directions
	 */
	private final JCheckBox			cbxDrawVertexCurvature;
	
	/**
	 * color model by curvature
	 */
	private final JCheckBox			cbxDrawCurvatureColor;
	
	/**
	 * Checkbox for displaying region edge boundaries
	 */
	private final JCheckBox			cbxDrawRegionEdges;
	
	/**
	 * draw sharp edges
	 */
	private final JCheckBox			cbxDrawSharpEdges;
	
	/**
	 * Checkbox for selectTrianglesOnly setting.
	 */
	private final JCheckBox			cbxSelectTriangles;
	
	/**
	 * select only nearest triangle or all triangles intersecting mouse ray
	 */
	private final JCheckBox			cbxSelectNearestOnly;
	
	/**
	 * draw bounding box for each group
	 */
	private final JCheckBox			cbxDrawBoundingBox;

	/**
	 * draw voronoi area
	 */
	private final JCheckBox			cbxDrawVoronoiArea;
	
	/**
	 * use white background
	 */
	private final JCheckBox			cbxWhiteBackground;
	
	/**
	 * button to set rotation of camera manually
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
	public DrawSettingsPanel(MeshCas cas, MeshReasoningView view) {
		this.view = view;
		this.cas = cas;
		// setPreferredSize(new Dimension(300, 300));

		GridLayout grid = new GridLayout(0, 2);
		setLayout(grid);

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
		
		cbxDrawCurvatureColor = new JCheckBox("Color by curvature");
		cbxDrawCurvatureColor.addActionListener(this);
		cbxDrawCurvatureColor.setSelected(false);
		this.add(cbxDrawCurvatureColor);

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
		if (e.getSource() == cbxWhiteBackground)
			view.setBackgroundWhite(cbxWhiteBackground.isSelected());
		else if (e.getSource() == cbxDrawVertexNormals)
			view.setDrawVertexNormals(cbxDrawVertexNormals.isSelected());
		else if (e.getSource() == cbxDrawTriangleNormals) 
			view.setDrawTriangleNormals(cbxDrawTriangleNormals.isSelected());
		else if (e.getSource() == cbxDrawVertexCurvatureMin)
			view.setDrawVertexCurvatureMin(cbxDrawVertexCurvatureMin.isSelected());
		else if (e.getSource() == cbxDrawVertexCurvatureMax)
			view.setDrawVertexCurvatureMax(cbxDrawVertexCurvatureMax.isSelected());
		else if (e.getSource() == cbxDrawVertexCurvature)
			view.setDrawVertexCurvature(cbxDrawVertexCurvature.isSelected());
		else if (e.getSource() == cbxDrawCurvatureColor)
			view.setDrawCurvatureColor(cbxDrawCurvatureColor.isSelected());
		else if (e.getSource() == cbxDrawVoronoiArea)
			view.setDrawVoronoiArea(cbxDrawVoronoiArea.isSelected());
		else if (e.getSource() == cbxDrawSharpEdges)
			view.setDrawSharpEdges(cbxDrawSharpEdges.isSelected());
		else if (e.getSource() == cbxDrawRegionEdges) 
			view.setDrawRegionEdges(cbxDrawRegionEdges.isSelected());
		else if (e.getSource() == cbxDrawBoundingBox)
			view.setDrawBoundingBox(cbxDrawBoundingBox.isSelected());
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
