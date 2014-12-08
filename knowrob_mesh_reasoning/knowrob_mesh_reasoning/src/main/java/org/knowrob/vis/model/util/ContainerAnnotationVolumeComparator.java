/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: 
 * 				Stefan Profanter - initial API and implementation, Year: 2013
 * 				Andrei Stoica - improved during Google Summer of Code 2014
 ******************************************************************************/
package edu.tum.cs.vis.model.util;

import java.util.Comparator;

import edu.tum.cs.vis.model.uima.annotation.ContainerAnnotation;

/**
 * Class that implements the Comparator abstract class of a ContainerAnnotation
 * It compares two container annotations based on their volume units. The comparison
 * returns the bigger ContainerAnnotation.
 * 
 * See:
 * {@link java.util.Comparator} and
 * {@link edu.tum.cs.vis.model.uima.annotation.ContainerAnnotation}
 * 
 * @author Stefan Profanter
 * @author Andrei Stoica - renamed (fixed naming misspell) and improved documentation
 */
public class ContainerAnnotationVolumeComparator implements Comparator<ContainerAnnotation> {

	/* (non-Javadoc)
	 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
	 */
	@Override
	public final int compare(ContainerAnnotation arg0, ContainerAnnotation arg1) {
		return Float.compare(arg0.getVolume(), arg1.getVolume()) * (-1);
	}

}
