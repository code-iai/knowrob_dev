/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Stefan Profanter - initial API and implementation, Year: 2012
 ******************************************************************************/
package edu.tum.cs.uima;

/**
 * Dummy for UIMA Framework
 * 
 * @author Stefan Profanter
 * 
 */
public abstract class AnalysisEngine {

	/**
	 * Dummy for UIMA framework
	 * 
	 * @param aJCas
	 *            CAS Object to analyse
	 */
	public abstract void process(JCas aJCas);

}
