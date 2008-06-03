package uk.ac.ucl.chem.ccs.vizclient;


public interface SteeringConnection {

	public VizFrameData getFrame ();
	
	public boolean isConnected();
	
	public boolean connect(int port, String hostname);
	
	public boolean connect();

	public boolean disconnect();

	public boolean magnify (int m);
	
	public boolean rotate (double dx, double dy);
	
}
