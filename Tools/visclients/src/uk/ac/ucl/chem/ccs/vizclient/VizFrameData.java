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

	//private long frameNo = 0;
	//private long realFrameNo = 0;
	private int bufferSize = 0;
	private int length = 0;
	//private double dataRate = 0.f;
	//private double framePerSec = 0.f;
	private int r[][];
	private int g[][];
	private int b[][];
	
	private int x[];
	private int y[];
	
	
    private double vis_pressure_min;
    private double vis_pressure_max;
    private double vis_velocity_min;
   private  double vis_velocity_max;
    private double vis_stress_min;
    private double vis_stress_max;
    private int vis_time_step;
    private double vis_time;
    private int vis_cycle;
    private int n_inlets;
    private double inlet_avg_vel[];
    private double vis_mouse_pressure;
    private double vis_stess_pressure;
	

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
	
	
	
	public double getVis_pressure_min() {
		return vis_pressure_min;
	}


	public void setVis_pressure_min(double vis_pressure_min) {
		this.vis_pressure_min = vis_pressure_min;
	}


	public double getVis_pressure_max() {
		return vis_pressure_max;
	}


	public void setVis_pressure_max(double vis_pressure_max) {
		this.vis_pressure_max = vis_pressure_max;
	}


	public double getVis_velocity_min() {
		return vis_velocity_min;
	}


	public void setVis_velocity_min(double vis_velocity_min) {
		this.vis_velocity_min = vis_velocity_min;
	}


	public double getVis_velocity_max() {
		return vis_velocity_max;
	}


	public void setVis_velocity_max(double vis_velocity_max) {
		this.vis_velocity_max = vis_velocity_max;
	}


	public double getVis_stress_min() {
		return vis_stress_min;
	}


	public void setVis_stress_min(double vis_stress_min) {
		this.vis_stress_min = vis_stress_min;
	}


	public double getVis_stress_max() {
		return vis_stress_max;
	}


	public void setVis_stress_max(double vis_stress_max) {
		this.vis_stress_max = vis_stress_max;
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
