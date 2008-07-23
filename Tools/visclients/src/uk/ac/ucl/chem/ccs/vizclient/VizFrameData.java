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

	private long frameNo = 0;
	private long realFrameNo = 0;
	private int bufferSize = 0;
	private int length = 0;
	private double dataRate = 0.f;
	private double framePerSec = 0.f;
	private int r[][];
	private int g[][];
	private int b[][];
	
	private int x[];
	private int y[];
	

	/**
	 * @param length
	 */
	public VizFrameData(int length) {
		this.length = length;
		r = new int[length][4];
		 g = new int[length][4];
		 b = new int[length][4];
		 x = new int[length];
		 y = new int[length];
	}
	

	public double getDataRate() {
		return dataRate;
	}



	/**
	 * @return Returns the framePerSec.
	 */
	public double getFramePerSec() {
		return framePerSec;
	}
	/**
	 * @param framePerSec The framePerSec to set.
	 */
	public void setFramePerSec(double framePerSec) {
		this.framePerSec = framePerSec;
	}
	public void setDataRate(double dataRate) {
		this.dataRate = dataRate;
	}



	/**
	 * @return Returns the bufferSize.
	 */
	public int getBufferSize() {
		return bufferSize;
	}
	/**
	 * @param bufferSize The bufferSize to set.
	 */
	public void setBufferSize(int bufferSize) {
		this.bufferSize = bufferSize;
	}
	
	
	
	public long getRealFrameNo() {
		return realFrameNo;
	}


	public void setRealFrameNo(long realFrameNo) {
		this.realFrameNo = realFrameNo;
	}


	/**
	 * @return Returns the frameNo.
	 */
	public long getFrameNo() {
		return frameNo;
	}
	/**
	 * @param frameNo The frameNo to set.
	 */
	public void setFrameNo(long frame_no) {
		this.frameNo = frame_no;
	}
	/**
	 * @return Returns the b.
	 */
	public int getB(int i, int view) {
		return b[i][view];
	}
	/**
	 * @param b The b to set.
	 */
	public void setB(int i, int view, int b) {
		this.b[i][view] = b;
	}
	/**
	 * @return Returns the g.
	 */
	public int getG(int i, int view) {
		return g[i][view];
	}
	/**
	 * @param g The g to set.
	 */
	public void setG(int i, int view, int g) {
		this.g[i][view] = g;
	}
	/**
	 * @return Returns the r.
	 */
	public int getR(int i, int view) {
		return r[i][view];
	}
	/**
	 * @param r The r to set.
	 */
	public void setR(int i, int view, int r) {
		this.r[i][view] = r;
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
