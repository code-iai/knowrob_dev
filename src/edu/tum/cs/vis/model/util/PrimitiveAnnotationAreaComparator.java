/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Stefan Profanter - initial API and implementation, Year: 2013
 ******************************************************************************/
package edu.tum.cs.vis.model.util;

import java.util.Comparator;

import edu.tum.cs.vis.model.uima.annotation.PrimitiveAnnotation;

/**
 * @author Stefan Profanter
 * 
 */
@SuppressWarnings("rawtypes")
public class PrimitiveAnnotationAreaComparator implements Comparator<PrimitiveAnnotation> {

	/* (non-Javadoc)
	 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
	 */
	@Override
	public int compare(PrimitiveAnnotation arg0, PrimitiveAnnotation arg1) {
		return Float.compare(arg0.getArea(), arg1.getArea());
	}

}
