/*
 * Created on 17:52:51 23-Jan-2008
 * Project: VIZ-GUI 
 * File: VizFrameData.java
 * 
 * @author stefan
 * 
 * TODO: 
 */
package uk.ac.ucl.chem.ccs.vizclient;

public class VizFrameData {

	int length = 0;
	
	int r[];
	int g[];
	int b[];
	
	int x[];
	int y[];
	
	
	/**
	 * @param length
	 */
	public VizFrameData(int length) {
		this.length = length;
		r = new int[length];
		 g = new int[length];
		 b = new int[length];
		 x = new int[length];
		 y = new int[length];
	}
	
	

	/**
	 * @return Returns the b.
	 */
	public int getB(int i) {
		return b[i];
	}
	/**
	 * @param b The b to set.
	 */
	public void setB(int i, int b) {
		this.b[i] = b;
	}
	/**
	 * @return Returns the g.
	 */
	public int getG(int i) {
		return g[i];
	}
	/**
	 * @param g The g to set.
	 */
	public void setG(int i, int g) {
		this.g[i] = g;
	}
	/**
	 * @return Returns the r.
	 */
	public int getR(int i) {
		return r[i];
	}
	/**
	 * @param r The r to set.
	 */
	public void setR(int i, int r) {
		this.r[i] = r;
	}
	/**
	 * @return Returns the x.
	 */
	public int getX(int i) {
		return x[i];
	}
	/**
	 * @param x The x to set.
	 */
	public void setX(int i, int x) {
		this.x[i] = x;
	}
	/**
	 * @return Returns the y.
	 */
	public int getY(int i) {
		return y[i];
	}
	/**
	 * @param y The y to set.
	 */
	public void setY(int i, int y) {
		this.y[i] = y;
	}
	/**
	 * @return Returns the length.
	 */
	public int getLength() {
		return length;
	}
}
