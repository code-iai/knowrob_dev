package edu.tum.cs.vis.model.util;

import java.util.ArrayList;

import javax.vecmath.Point3f;

import processing.core.PApplet;

/**
 * Represents a Mesh-Group with a Name.
 * A group can multiple child groups so it will present a mesh hierarchy.
 * 
 * @author Stefan Profanter
 *
 */
public class Group {
	
	/**
	 * Minimum position values of group
	 */
	protected Float minX = null;
	protected Float maxX = null;
	protected Float minY = null;
	protected Float maxY = null;
	protected Float minZ = null;
	protected Float maxZ = null;
	
	/**
	 * Holds the Mash for this group
	 */
	private Mesh mesh = new Mesh();
	
	/**
	 * Name of this group
	 */
	private String name;
	
	/**
	 * Child groups
	 */
	private ArrayList<Group> children = new ArrayList<Group>();

	/**
	 * Get the mesh containing triangles and lines
	 * @return the mesh
	 */
	public Mesh getMesh() {
		return mesh;
	}

	/**
	 * Set the mesh containing triangles and lines
	 * @param mesh the mesh
	 */
	public void setMesh(Mesh mesh) {
		this.mesh = mesh;
	}

	/**
	 * Get name of this group (normally parsed from model file)
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	/**
	 * Set name of this group (normally parsed from model file)
	 * @return the name
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * Get all direct child groups
	 * @return list with all child groups
	 */
	public ArrayList<Group> getChildren() {
		return children;
	}

	/**
	 * Set child groups of this group. Used to create a hierarchy
	 * @param children list of children
	 */
	public void setChildren(ArrayList<Group> children) {
		this.children = children;
	}
	
	/**
	 * Set minimum x coordinate of this group. No coordinate of all its children should be smaller than this value.
	 * Used for drawing a bounding box.
	 * @return minimum x
	 */
	public Float getMinX() {
		return minX;
	}

	/**
	 * Set maximum x coordinate of this group. No coordinate of all its children should be bigger than this value.
	 * Used for drawing a bounding box.
	 * @return maximum x
	 */
	public Float getMaxX() {
		return maxX;
	}
	
	/**
	 * Set minimum y coordinate of this group. No coordinate of all its children should be smaller than this value.
	 * Used for drawing a bounding box.
	 * @return minimum y
	 */
	public Float getMinY() {
		return minY;
	}

	/**
	 * Set maximum y coordinate of this group. No coordinate of all its children should be bigger than this value.
	 * Used for drawing a bounding box.
	 * @return maximum y
	 */
	public Float getMaxY() {
		return maxY;
	}

	/**
	 * Set minimum z coordinate of this group. No coordinate of all its children should be smaller than this value.
	 * Used for drawing a bounding box.
	 * @return minimum z
	 */
	public Float getMinZ() {
		return minZ;
	}

	/**
	 * Set maximum z coordinate of this group. No coordinate of all its children should be bigger than this value.
	 * Used for drawing a bounding box.
	 * @return maximum z
	 */
	public Float getMaxZ() {
		return maxZ;
	}
	
	/**
	 * Call this function only in the root group!
	 * It will center the whole Group with its children to the Point (0,0,0)
	 * which means the volumetric center will then be at point 0.
	 * 
	 */
	public void initialize(String textureBasePath) {
		initialize(textureBasePath, true);
	}
	
	/**
	 * Call this function only in the root group!
	 * It will center the whole Group with its children to the Point (0,0,0)
	 * which means the volumetric center will then be at point 0.
	 * 
	 */
	private void initialize(String textureBasePath,boolean isRoot) {
		//Initialize min and max values of mesh
		mesh.setTextureBasePath(textureBasePath);
		mesh.getWidth();
		mesh.getHeight();
		mesh.getDepth();
		for (Group g : children)
		{
			g.initialize(textureBasePath,false);
		}
		getTotalWidth();
		getTotalHeight();
		getTotalDepth();
		if (isRoot)
		{
			Point3f translation = new Point3f((minX + (getTotalWidth() / 2))*(-1), (minY
						+ (getTotalHeight() / 2))*(-1), (minZ + (getTotalDepth() / 2))*(-1));
			this.translate(translation);
		}
	}
	
	/**
	 * Translates (moves) the group and its children by the specified vector.
	 * @param translation Translation vector
	 */
	private void translate(Point3f translation) {
		minX += translation.x;
		maxX += translation.x;
		minY += translation.y;
		maxY += translation.y;
		minZ += translation.z;
		maxZ += translation.z;
		mesh.translate(translation);
		for (Group g : children)
		{
			g.translate(translation);
		}
	}
	
	/**
	 * Scales all coordinates by the given factor. Only used for initialization
	 * 
	 * @param factor
	 *            The scale factor
	 */
	public void scale(float factor) {
		mesh.scaleMesh(factor);

		for (Group g : children)
		{
			g.scale(factor);
		}
		minX *= factor;
		maxX *= factor;
		minY *= factor;
		maxY *= factor;
		minZ *= factor;
		maxZ *= factor;
	}
	
	
	/**
	 * Returns the total height of the group by searching the biggest distance on the
	 * y-axis between the vectors.
	 * 
	 * @return float as height of the group
	 */
	public float getTotalHeight() {
		if (minY == null || maxY == null)
		{
			minY = mesh.getMinY();
			maxY = mesh.getMaxY();
			
			for (Group g : children)
			{
				minY = Math.min(minY, g.getMinY());
				maxY = Math.max(maxY, g.getMaxY());
				
			}
		}
		return maxY - minY;
	}

	/**
	 * Returns the total width of the group by searching the biggest distance on the
	 * x-axis between the vectors.
	 * 
	 * @return float as width of the group
	 */
	public float getTotalWidth() {
		if (minX == null || maxX == null)
		{
			minX = mesh.getMinX();
			maxX = mesh.getMaxX();
			
			for (Group g : children)
			{
				minX = Math.min(minX, g.getMinX());
				maxX = Math.max(maxX, g.getMaxX());
				
			}
		}
		return maxX - minX;
	}

	/**
	 * Returns the total depth of the group by searching the biggest distance on the
	 * z-axis between the vectors.
	 * 
	 * @return float as depth of the group
	 */
	public float getTotalDepth() {
		if (minZ == null || maxZ == null)
		{
			minZ = mesh.getMinZ();
			maxZ = mesh.getMaxZ();
			
			for (Group g : children)
			{
				minZ = Math.min(minZ, g.getMinZ());
				maxZ = Math.max(maxZ, g.getMaxZ());
				
			}
		}
		return maxZ - minZ;
	}
	
	/**
	 * Draw method to draw the model on the applet.
	 * @param applet The applet to draw on.
	 * @param colorOverride override the draw color an texture. Draw whole object in the given color if != 0
	 */
	public void draw(PApplet applet, int overrideColor) {
		mesh.drawLines(applet, overrideColor);
		mesh.drawTriangles(applet, overrideColor);
		for (Group g : children)
		{
			g.draw(applet, overrideColor);			
		}
	}
	
	/**
	 * Draws the bounding box around the model with the current style
	 * 
	 * @param applet
	 *            Applet to draw on
	 *            
	 * @param recursive
	 * 			  Draw also Bounding-Box of children
	 */
	public void drawBoundingBox(PApplet applet, boolean recursive) {
		mesh.drawBoundingBox(applet);
		if (recursive)
		{
			for (Group g : children)
			{
				g.drawBoundingBox(applet, recursive);			
			}
		}
	}
	
	/**
	 * Add a child to this group
	 * @param g group to add as a child of this
	 */
	public void addChild(Group g)
	{
		children.add(g);
	}

	/**
	 * Mirrors whole group and children on the x coordinate.by setting each x to the inverse -x
	 */
	public void mirrorX() {
		mesh.mirrorX();
		minX *= (-1);
		float tmp = minX;
		minX = maxX * (-1);
		maxX = tmp;
		for (Group g : children)
		{
			g.mirrorX();			
		}	
	}
}
