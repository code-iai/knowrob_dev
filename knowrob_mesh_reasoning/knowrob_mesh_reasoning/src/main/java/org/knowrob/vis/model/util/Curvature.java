/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: 
 * 				Stefan Profanter - initial API and implementation, Year: 2012
 * 				Andrei Stoica - refactored during Google Summer of Code 2014
 ******************************************************************************/
package org.knowrob.vis.model.util;

import java.awt.Color;

import javax.vecmath.Vector3f;

import org.knowrob.vis.model.uima.annotation.primitive.PrimitiveType;

/**
 * Class which represents a curvature property for a vertex or a triangle.
 * It contains the curvature tensor information containing the principle
 * directions of the estimated curvature property along with their associated
 * scalar values, i.e. the minimum curvature value kMin and its principle direction
 * 3D vector pMin, and respectively, the maximum curvature value kMax and its principle
 * direction 3D vector pMax. Additional information such as the mean and Gaussian 
 * curvatures properties are also available. Besides two distance properties
 * used to classify and map the curvature to a coloring scheme are also available.
 * These two fields are the {@code hue} and the {@code saturation}. Given the 
 * averaging processing performed on the curvatures, initial values of these two 
 * properties are also available, since they are necessary for visualization 
 * purposes. 
 * 
 * @author Stefan Profanter
 * @author Andrei Stoica - extended and refactored the initial API and implementation
 */
public class Curvature {

	/**
	 * Curvature minimum principle direction
	 */
	private Vector3f		principleDirectionMin	= new Vector3f();

	/**
	 * Curvature maximum principle direction
	 */
	private Vector3f		principleDirectionMax	= new Vector3f();

	/**
	 * Minimum curvature floating point value
	 */
	private float			curvatureMin			= 0;
	
	
	/**
	 * Maximum curvature flaoting point value
	 */
	private float			curvatureMax			= 0;

	/**
	 * Value between minimum and maximum curvature as estimated from
	 * solving the second normal form through the LU decomposition
	 * in the curvature calculation process
	 * 
	 * {@link edu.tum.cs.vis.model.util.algorithm.CurvatureCalculation#calculateCurvatures(java.util.HashMap, edu.tum.cs.vis.model.Model)}
	 */
	private float			curvatureMinMax			= 0;
	
	/** 
	 * Mean curvature value := (kMin + kMax) / 2.
	 * Intended for storing the initial mean curvature value
	 */
	private float 			meanCurvature			= 0;

	/**
	 * Gaussian curvature value := kMin * kMax.
	 * Intended for storing the initial Gaussian curvature value 
	 */
	private float 			gaussCurvature			= 0;

	/**
	 * Primitive type given the estimated curvature properties
	 */
	private PrimitiveType	primitiveType;

	/**
	 * Currently computed hue from the curvature properties
	 */
	private float			hue;
	
	/**
	 * Currently computed saturation from curvature properties
	 */
	private float			saturation;

	/**
	 * Initially computed hue based on the initial curvature estimation
	 */
	private float 			initialHue;
	
	/**
	 * Initially computed saturation based on the initial curvature estimation
	 */
	private float			initialSaturation;
	
	/**
	 * Converts HSV to RGB color space according to:
	 * https://en.wikipedia.org/wiki/HSL_and_HSV#From_HSV
	 *
	 * The conversion is done from the scalar field of hue and
	 * saturation with their initial computed values to the
	 * RGB color scheme.
	 * 
	 * @param h
	 *            hue
	 * @param s
	 *            saturation
	 * @param v
	 *            value
	 * @return Color represented by HSV
	 */
	private static final Color hsv2rgb(float h, float s, float v) {
		float H = h, S = s, V = v;
		if (S <= 0.0f) {
			S = Math.abs(S);
			// return new Color(V, V, V);
		}
		H = (float) (H % (2 * Math.PI));
		if (H < 0.0f) {
			H += 2 * Math.PI;
		}
		// S and V is now between 0 and 1, H between 0 an 2*PI
		float hi = (float) (H * (Math.PI / 3)); // Divide by 60 degree
		int i = (int) Math.floor(hi);
		float f = hi - i;
		float p = V * (1.0f - S);
		float q = V * (1.0f - (S * f));
		float t = V * (1.0f - (S * (1.0f - f)));
		switch (i) {
			case 0:
				return new Color(V, t, p);
			case 1:
				return new Color(q, V, p);
			case 2:
				return new Color(p, V, t);
			case 3:
				return new Color(p, q, V);
			case 4:
				return new Color(t, p, V);
			default:
				return new Color(V, p, q);
		}
	}
	
	/**
	 * Converts HSV to RGB color space according to:
	 * https://en.wikipedia.org/wiki/HSL_and_HSV#From_HSV
	 * a scalar field of the curvature tensor
	 * 
	 * The conversion is done from the scalar field of the
	 * initially estimated mean or Gaussian curvature to the
	 * RGB color scheme. Minimum and maximum values of the
	 * scalars are needed for compressing the mapping all the
	 * components on the HSV scale which is then converted 
	 * to RGB.
	 * 
	 * @param H
	 * 				hue
	 * @param hMin
	 * 				minimum hue value for the scale
	 * @param hMax
	 * 				maximum hue value for the scale
	 * @param S
	 * 				saturation to be used (needs to be in between 0 and 1)
	 * @param V
	 * 				value to be used (needs to be in between 0 and 1)
	 * @return color
	 * 				RGB color
	 */
	private static final Color hsv2rgb(float H, float hMin, float hMax, float S, float V) {
		if (S < 0) {
			S = 0.2f;
		}
		if (S > 1) {
			S = 1;
		}
		if (H > hMax) {
			H = hMax;
		}
		if (H < hMin) {
			H = hMin;
		}
		H = 240f - 240f * (H - hMin) / (hMax - hMin);
		float finterval = (float) (H / 60f); // Divide by 60 degree
		int interval = (int) Math.floor(finterval);
		float f = finterval - interval;
		float p = V * (1.0f - S);
		float q = V * (1.0f - (S * f));
		float t = V * (1.0f - (S * (1.0f - f)));
		switch (interval) {
			case 0:
				return new Color(V, t, p);
			case 1:
				return new Color(q, V, p);
			case 2:
				return new Color(p, V, t);
			case 3:
				return new Color(p, q, V);
			case 4:
				return new Color(t, p, V);
			default:
				return new Color(V, p, q);
		}
	}
	
	/**
	 * Gets the estimated curvatures scalar field coloring using the
	 * initially computed hue and saturation scalar values.
	 * 
	 * @return RGB color
	 * 			according to the hue and saturation values
	 */
	public Color getColor() {
		return hsv2rgb(this.initialHue, this.initialSaturation, 1);
	}
	
	/**
	 * Gets the mean curvature coloring for the model.
	 * 
	 * @param lowMeanCurvature
	 * 			lowest mean curvature reference of the CAD model
	 * @param highMeanCurvature
	 * 			highest mean curvature reference of the CAD model
	 * @return color by mean curvature value
	 */
	public Color getMeanColor(final float lowMeanCurvature, final float highMeanCurvature) {
		return hsv2rgb(this.meanCurvature, lowMeanCurvature, highMeanCurvature, this.saturation, 1);
	}
	
	/**
	 * Gets the Gaussian curvature coloring for the model.
	 * 
	 * @param lowGaussCurvature
	 * 			lowest Gaussian curvature reference of the CAD model
	 * @param highGaussCurvature
	 * 			highest Gaussian curvature reference of the CAD model
	 * @return color by Gauss curvature value
	 */
	public Color getGaussColor(final float lowGaussCurvature, final float highGaussCurvature) {
		return hsv2rgb(this.gaussCurvature, lowGaussCurvature, highGaussCurvature, this.saturation, 1);
	}

	/**
	 * Gets the maximum curvature value
	 * 
	 * @return the maximum curvature - kMax
	 */
	public float getCurvatureMax() {
		return this.curvatureMax;
	}

	/**
	 * Gets minimum curvature value
	 * 
	 * @return the minimum curvature - kMin
	 */
	public float getCurvatureMin() {
		return this.curvatureMin;
	}

	/**
	 * Gets the min-max curvature value of the curvature calculation.
	 * 
	 * @return the min-max curvature property - kMinMax
	 */
	public float getCurvatureMinMax() {
		return this.curvatureMinMax;
	}
	
	/**
	 * Gets the mean curvature value
	 * 
	 * @return the mean curvature value
	 */
	public float getMeanCurvature() {
		return this.meanCurvature;
	}
	
	/**
	 * Gets the Gaussian curvature value
	 * 
	 * @return the Gaussian curvature value
	 */
	public float getGaussCurvature() {
		return this.gaussCurvature;
	}

	/**
	 * Gets the initially computed hue based on the estimated
	 * curvature values for the vertices of the CAD model
	 * 
	 * @return the initialHue
	 */
	public float getInitialHue() {
		return this.initialHue;
	}
	
	/**
	 * Gets the hue of the initial curvature tensor used for 
	 * coloring the model by curvature
	 * 
	 * @return the hue
	 */
	public float getHue() {
		return this.hue;
	}

	/**
	 * Gets the primitive type associated with the curvature properties
	 * 
	 * @return the primitiveType
	 */
	public PrimitiveType getPrimitiveType() {
		return this.primitiveType;
	}

	/**
	 * Gets the principle direction 3D vector corresponding to the 
	 * maximum curvature
	 * 
	 * @return the principleDirectionMax
	 */
	public Vector3f getPrincipleDirectionMax() {
		return this.principleDirectionMax;
	}

	/**
	 * Gets the principle direction 3D vector corresponding to the
	 * minimum curvature
	 * 
	 * @return the principleDirectionMin
	 */
	public Vector3f getPrincipleDirectionMin() {
		return this.principleDirectionMin;
	}

	/**
	 * Gets the initially computed saturation based on the estimated curvature
	 * values of the vertices of the CAD model
	 * 
	 * @return the initialSaturation
	 */
	public float getInitialSaturation() {
		return this.initialSaturation;
	}
	
	/**
	 * Get saturation for curvature values used for coloring model by curvature
	 * 
	 * @return the saturation
	 */
	public float getSaturation() {
		return this.saturation;
	}
	
	/**
	 * Gets the mean curvature value based on the min and max curvature values of the Curvature
	 */
	public void setMeanCurvature() {
		this.meanCurvature = (this.curvatureMin + this.curvatureMax) / 2f;
	}
	
	/**
	 * Gets the Gaussian curvature value based on the min and max curvature values of the Curvature
	 */
	public void setGaussCurvature() {
		this.gaussCurvature = this.curvatureMin * this.curvatureMax;
	}

	/**
	 * Set maximum curvature magnitude
	 * 
	 * @param curvatureMax
	 *            the curvatureMax to set
	 */
	public void setCurvatureMax(float curvatureMax) {
		this.curvatureMax = curvatureMax;
	}

	/**
	 * Set minimum curvature magnitude
	 * 
	 * @param curvatureMin
	 *            the curvatureMin to set
	 */
	public void setCurvatureMin(float curvatureMin) {
		this.curvatureMin = curvatureMin;
	}

	/**
	 * Set min max curvature magnitude
	 * 
	 * @param curvatureMinMax
	 *            the curvatureMinMax to set
	 */
	public void setCurvatureMinMax(float curvatureMinMax) {
		this.curvatureMinMax = curvatureMinMax;
	}

	/**
	 * Set initial hue value of the vertex curvature
	 * 
	 * @param hue
	 * 			the hue to be set
	 */
	public void setInitialHue(final float hue) {
		this.initialHue = hue;
	}
	
	/**
	 * set hue used for coloring model by curvature
	 * 
	 * @param hue
	 *            the hue to be set
	 */
	public void setHue(float hue) {
		this.hue = hue;
	}

	/**
	 * set primitive type which represents these curvature properties
	 * 
	 * @param primitiveType
	 *            the primitiveType to set
	 */
	public void setPrimitiveType(PrimitiveType primitiveType) {
		this.primitiveType = primitiveType;
	}

	/**
	 * Set max principle direction
	 * 
	 * @param principleDirectionMax
	 *            the principleDirectionMax to set
	 */
	public void setPrincipleDirectionMax(Vector3f principleDirectionMax) {
		this.principleDirectionMax = principleDirectionMax;
	}

	/**
	 * Set min principle direction
	 * 
	 * @param principleDirectionMin
	 *            the principleDirectionMin to set
	 */
	public void setPrincipleDirectionMin(Vector3f principleDirectionMin) {
		this.principleDirectionMin = principleDirectionMin;
	}

	/**
	 * Set initial saturation of the curvature vertex
	 * 
	 * @param saturation
	 * 			the saturation to be set
	 */
	public void setInitialSaturation(final float saturation) {
		this.initialSaturation = saturation;
	}
	
	/**
	 * Set hue used for coloring model by curvature
	 * 
	 * @param saturation
	 *            the saturation to set
	 */
	public void setSaturation(float saturation) {
		this.saturation = saturation;
	}
	
	@Override
	public String toString() {
		String ret = "";
		ret += "min: " + this.curvatureMin;
		ret += "\nmax: " + this.curvatureMax;
		ret += "\nminmax: " + this.curvatureMinMax;
		ret += "\ndirMin: " + this.principleDirectionMin.toString();
		ret += "\ndirMax: " + this.principleDirectionMax.toString();
		ret += "\nprimitiveType: " + this.primitiveType;
		return ret;
	}

}
