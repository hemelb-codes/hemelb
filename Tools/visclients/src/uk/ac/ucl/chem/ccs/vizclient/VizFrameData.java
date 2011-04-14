/*
 * Created on 17:52:51 23-Jan-2008
 * Project: VIZ-GUI 
 * File: VizFrameData.java
 * 
 * @author Stefan Zasada
 * @author Robert Haines
 * 
 * TODO: 
 */
package uk.ac.ucl.chem.ccs.vizclient;

import java.nio.IntBuffer;

public class VizFrameData {

	//private long frameNo = 0;
	//private long realFrameNo = 0;
	private int bufferSize = 0;
	//private double dataRate = 0.f;
	//private double framePerSec = 0.f;
	
    private int vis_time_step;
    private double vis_time;
    private int vis_cycle;
    private int n_inlets;
    private double inlet_avg_vel[];
    private double vis_mouse_pressure;
    private double vis_stess_pressure;
	
    private int[][] pixels;
    private int width;
    private int height;
    
	/**
	 * @param length
	 */
    public VizFrameData(int width, int height) {
		this.width = width;
		this.height = height;
		 
		pixels = new int[4][width * height];
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < (width * height); j ++) {
				pixels[i][j] = (255 << 24) | (255 << 16) | (255 << 8);
			}
		}
	}
	
	public int getWidth() {
		return width;
	}

	public int getHeight() {
		return height;
	}

	public void setPixel(int index, int view, int red, int green, int blue) {
		int pixel = (red << 24) | (green << 16) | (blue << 8);
		pixels[view][index] = pixel;
	}
	
	public IntBuffer getBuffer(int view) {
		return IntBuffer.wrap(pixels[view]);
	}
	
//	public double getDataRate() {
//		return dataRate;
//	}



//	/**
//	 * @return Returns the framePerSec.
//	 */
//	public double getFramePerSec() {
//		return framePerSec;
//	}
	/**
//	 * @param framePerSec The framePerSec to set.
//	 */
//	public void setFramePerSec(double framePerSec) {
//		this.framePerSec = framePerSec;
//	}
//	public void setDataRate(double dataRate) {
//		this.dataRate = dataRate;
//	}



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
	

	public int getVis_time_step() {
		return vis_time_step;
	}


	public void setVis_time_step(int vis_time_step) {
		this.vis_time_step = vis_time_step;
	}


	public double getVis_time() {
		return vis_time;
	}


	public void setVis_time(double vis_time) {
		this.vis_time = vis_time;
	}


	public int getVis_cycle() {
		return vis_cycle;
	}


	public void setVis_cycle(int vis_cycle) {
		this.vis_cycle = vis_cycle;
	}


	public int getN_inlets() {
		return n_inlets;
	}


	public void setN_inlets(int n_inlets) {
		this.n_inlets = n_inlets;
		inlet_avg_vel = new double[n_inlets];
	}


	public double getInlet_avg_vel(int pos) {
		return inlet_avg_vel[pos];
	}


	public void setInlet_avg_vel(double inlet_vel, int pos) {
		this.inlet_avg_vel[pos] = inlet_vel;
	}


	public double getVis_mouse_pressure() {
		return vis_mouse_pressure;
	}


	public void setVis_mouse_pressure(double vis_mouse_pressure) {
		this.vis_mouse_pressure = vis_mouse_pressure;
	}


	public double getVis_stess_pressure() {
		return vis_stess_pressure;
	}


	public void setVis_stess_pressure(double vis_stess_pressure) {
		this.vis_stess_pressure = vis_stess_pressure;
	}


//	public long getRealFrameNo() {
//		return realFrameNo;
//	}


//	public void setRealFrameNo(long realFrameNo) {
//		this.realFrameNo = realFrameNo;
//	}


//	/**
//	 * @return Returns the frameNo.
//	 */
//	public long getFrameNo() {
//		return frameNo;
//	}
//	/**
//	 * @param frameNo The frameNo to set.
//	 */
//	public void setFrameNo(long frame_no) {
//		this.frameNo = frame_no;
//	}
}
