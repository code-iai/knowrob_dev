package edu.tum.cs.vis.model;

import java.util.HashMap;
import java.util.Vector;

import edu.tum.cs.vis.gui.applet.PrologUtil;

/**
 * Static class to obtain model properties such as path to model file,
 * width, height, depth.
 * 
 * 
 * @author Stefan Profanter
 * @see ItemModel
 */
public class Properties {
	
	/**
	 * Instantiate an ItemModel from prolog identifier of Instance or Class.
	 * The path to the model will be gathered with prolog an this model will then 
	 * be parsed.
	 * @param identifier Class or Instance identifier
	 * @return ItemModel instance if model parsed successfully. null otherwise.
	 */
	public static ItemModel getModelOfItem(String identifier) {
		try {
			if (!identifier.startsWith("'") || !identifier.endsWith("'"))
			{
				identifier = "'" + identifier + "'";
			} 
			HashMap<String, Vector<Object>> nfo = PrologUtil
					.executeQuery(
							"rdf_has(" + identifier	+ ",knowrob:pathToCadModel,literal(P)) ; "
							+ "owl_individual_of(" + identifier	+ ", Class), class_properties(Class, knowrob:pathToCadModel,literal(P))",
							null);
			
			if (nfo.get("P") != null && nfo.get("P").size() > 0) {
				String str = nfo.get("P").get(0).toString();
				
				if (str == null)
					return null;
				
				if (str.startsWith("'") && str.endsWith("'")) {
					str = str.substring(1, str.length()-1);
				}
				return new ItemModel(str);
				
			} else return null;
			
		} catch (Exception e) {
			System.err.println("Parsing of " + identifier  + " failed:");
			e.printStackTrace();
			return null;
		}
	}
	
	/**
	 * Get x-Dimension (width) of class or instance identified by ObjClassOrIdentifier from the CAD-Model.
	 * 
	 * @param ObjClassOrIdentifier ClassName or instance identifier
	 * @return x-Dimension or -1 if model not found
	 */
	public static float xDimOfObject(String ObjClassOrIdentifier) {
		ItemModel model = getModelOfItem(ObjClassOrIdentifier);
		if (model == null || model.getParser() == null)
			return -1;
		return model.getParser().getGroup().getTotalWidth();
	}
	

	
	/**
	 * Get y-Dimension (depth) of class or instance identified by ObjClassOrIdentifier from the CAD-Model.
	 * 
	 * @param ObjClassOrIdentifier ClassName or instance identifier
	 * @return y-Dimension or -1 if model not found
	 */
	public static float yDimOfObject(String ObjClassOrIdentifier) {
		ItemModel model = getModelOfItem(ObjClassOrIdentifier);
		if (model == null || model.getParser() == null)
			return -1;
		return model.getParser().getGroup().getTotalDepth();
	}
	

	/**
	 * Get z-Dimension (height) of class or instance identified by ObjClassOrIdentifier from the CAD-Model.
	 * 
	 * @param ObjClassOrIdentifier ClassName or instance identifier
	 * @return z-Dimension or -1 if model not found
	 */
	public static float zDimOfObject(String ObjClassOrIdentifier) {
		ItemModel model = getModelOfItem(ObjClassOrIdentifier);
		if (model == null || model.getParser() == null)
			return -1;
		return model.getParser().getGroup().getTotalHeight();
	}
	
	
}
