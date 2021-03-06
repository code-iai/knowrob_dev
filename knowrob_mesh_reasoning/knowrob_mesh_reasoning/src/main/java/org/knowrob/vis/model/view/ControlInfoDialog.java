/*******************************************************************************
 * Copyright (c) 2013 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: 
 * 				Stefan Profanter - initial API and implementation, Year: 2013
 * 				Andrei Stoica - refactored implementation Google Summer of Code 2014
 ******************************************************************************/
package org.knowrob.vis.model.view;

import java.awt.BorderLayout;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

/**
 * Dialog window showing information on available shortcut keys.
 * 
 * @author Stefan Profanter
 * @author Andrei Stoica - corrected dialog info wrong options
 */
public class ControlInfoDialog extends JDialog implements ActionListener {

	/**
	 * auto generated
	 */
	private static final long		serialVersionUID	= -3897203470825310360L;

	/**
	 * Close button for dialog
	 */
	private final JButton			closeBtn;

	/**
	 * Main panel of dialog
	 */
	private JPanel					myPanel				= null;

	/**
	 * List of available shortcuts
	 */
	private final static String[][]	shortcuts			= new String[][] {
			{ "left click", "Camera rotation" },
			{ "right click", "Camera move" },
			{ "mouse wheel", "Camera zoom" },
			{ ",", "Zoom in" },
			{ ".", "Zoom out" },
			{ "+", "Increase scale" },
			{ "-", "Decrease scale" },
			{ "1", "Draw mode: Fill" },
			{ "2", "Draw mode: Lines" },
			{ "3", "Draw mode: Points" },
			{ "q", "Decrease line width" },
			{ "w", "Increase line width" },
			{ "<html>Shift+Click on<br />Annotation checkbox</html>",
			"<html>Toggle all<br />other annotations</html>" } };

	/**
	 * Main constructor for ControlInfoDialog
	 * 
	 * @param owner
	 *            Owner of dialog.
	 */
	public ControlInfoDialog(Frame owner) {
		super(owner, true);
		setTitle("Control info");
		myPanel = new JPanel(new BorderLayout());
		getContentPane().add(myPanel);

		JPanel pnlLeft = new JPanel(new GridBagLayout());
		JPanel pnlCenter = new JPanel(new GridBagLayout());
		myPanel.add(pnlLeft, BorderLayout.WEST);
		myPanel.add(pnlCenter, BorderLayout.CENTER);

		GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.HORIZONTAL;
		c.insets = new Insets(10, 10, 10, 10);

		for (int i = 0; i < shortcuts.length; i++) {
			c.gridx = 0;
			c.gridy = i;
			pnlLeft.add(new JLabel(shortcuts[i][0], SwingConstants.RIGHT), c);
			c.gridx = 0;
			c.gridy = i;
			pnlCenter.add(new JLabel(shortcuts[i][1]), c);
		}

		closeBtn = new JButton("close");
		closeBtn.addActionListener(this);
		myPanel.add(closeBtn, BorderLayout.SOUTH);
		pack();
		setLocationRelativeTo(owner);
		setVisible(true);

	}

	@Override
	public void actionPerformed(ActionEvent e) {
		if (closeBtn == e.getSource()) {
			setVisible(false);
		}
	}
}
