package uk.ac.ucl.chem.ccs.vizclient;

import java.io.*;
import java.net.*;
import java.util.Calendar;
public class DirectBiConnection implements SteeringConnection {

	private static final int BYTES_PER_PIXEL_DATA = 16;
	private DataInputStream d;
	private DataOutputStream dos;
	private Socket listenSocket;	
	
	// data per pixel are colour id and pixel id (2 * sizeof(int) * bytes)
	private int bits_per_char = 8;
	// r-, g-, b- pixel components need "bits_per_char" bits
	private int colour_mask = (1 << bits_per_char) - 1;
	// x-, y- pixel components need "bits_per_two_chars" bits
	//private int pixel_mask = (1 << bits_per_two_chars) - 1;		
	private int colour_data[] = {0, 0, 0};
	//private int pixel_data = 0;	
	private long frame_no = 0;
	private boolean connected = false;
	private String hostname;
	private int port, window;
	/**
	 * 
	 */
	public DirectBiConnection(int port, String hostname, int window) {
		// TODO Auto-generated method stub
		this.hostname = hostname;
		this.port = port;		
		this.window = window;
	}
	
	public DirectBiConnection(int port, String hostname) {
		this(port, hostname, 1024*1024);
	}
	
	public VizFrameData getFrame () {
		//check if connected?
		if (connected) {
		int frame_size=0;
		frame_no++;		
		Calendar cal = Calendar.getInstance();
		short x_data=0, y_data=0;
		
		long start_time = cal.getTimeInMillis();
		//System.err.println("Frame no " + frame_no);

			try {				
				frame_size = d.readInt();

				if (frame_size < 1) {
					return null;
				}
			} catch (IOException e) {
				connected = false;
				return null;
			}
	
		int col_pixels = frame_size/BYTES_PER_PIXEL_DATA;
		VizFrameData vizFrame = new VizFrameData(col_pixels);

		for (int i =0; i < col_pixels; i++) {

			try {
				
				x_data = d.readShort();
				y_data = d.readShort();
				
					colour_data[0] = d.readInt();
					colour_data[1] = d.readInt();
					colour_data[2] = d.readInt();

	
					
					//System.out.println("x= " + x_data + " y " + y_data);
				} catch (IOException e1) {
					connected = false;
					return null;
				}

			vizFrame.setR(i,0,(colour_data[0] >> (3*bits_per_char)) & colour_mask);
			vizFrame.setG(i,0,(colour_data[0] >> (2*bits_per_char)) & colour_mask);
			vizFrame.setB(i,0,(colour_data[0] >> (1*bits_per_char)) & colour_mask);
			vizFrame.setR(i,1,(colour_data[0] >> (0*bits_per_char)) & colour_mask);
			vizFrame.setG(i,1,(colour_data[1] >> (3*bits_per_char)) & colour_mask);
			vizFrame.setB(i,1,(colour_data[1] >> (2*bits_per_char)) & colour_mask);	    
			vizFrame.setR(i,2,(colour_data[1] >> (1*bits_per_char)) & colour_mask);
			vizFrame.setG(i,2,(colour_data[1] >> (0*bits_per_char)) & colour_mask);
			vizFrame.setB(i,2,(colour_data[2] >> (3*bits_per_char)) & colour_mask);
			vizFrame.setR(i,3,(colour_data[2] >> (2*bits_per_char)) & colour_mask);
			vizFrame.setG(i,3,(colour_data[2] >> (1*bits_per_char)) & colour_mask);
			vizFrame.setB(i,3,(colour_data[2] >> (0*bits_per_char)) & colour_mask);
			
			vizFrame.setX(i, x_data);
			vizFrame.setY(i, y_data);	    
		}

		Calendar cal2 = Calendar.getInstance();
		// PLZ CAN HAZ DATEZ PLS! k thnx bi.
		long end_time = cal2.getTimeInMillis();
		double total_time = (end_time - start_time)*1000;
		double data_rate = 1024.0*(frame_size / total_time);
				
		//System.err.println("bytes = " + frame_size + ", time = " + total_time + ", rate = " + data_rate + "KB/s");
		vizFrame.setDataRate(data_rate);
		vizFrame.setRealFrameNo(0);
		vizFrame.setFrameNo(frame_no);
		vizFrame.setBufferSize(frame_size);
		vizFrame.setFramePerSec(1/total_time);
		return vizFrame;
		} 
		return null;

	}
	
	public boolean isConnected() {
		return connected;
	}
	
	public boolean connect(int port, String hostname) {
		this.port = port;
		this.hostname = hostname;
		return connect();
	}
	
	public boolean connect() {
		if (!connected) {
		int s=0;
		try {
			listenSocket = new Socket();
			//set TCP buffer size
			listenSocket.setReceiveBufferSize(window);
			s=listenSocket.getReceiveBufferSize();
			listenSocket.connect(new InetSocketAddress(hostname, port));
			dos = new DataOutputStream(listenSocket.getOutputStream());
       		BufferedInputStream bufff = new BufferedInputStream(listenSocket.getInputStream());
			d = new DataInputStream(bufff);

			connected = true;
		} catch (UnknownHostException e) {
			System.err.println("can't connect to host: " + hostname);
		} catch (IOException e) {
			System.err.println("Couldn't get I/O for the connection to: " + hostname);
		}

		System.err.println("RecSize " + s);
		}
		return connected;
	}
	
	
	public boolean disconnect() {
		if (connected) {
		try {
			dos.close();
			d.close();
			listenSocket.close();
			connected = false;
		} catch (Exception e) {
			System.err.println("couldn't close connection");
		}
		}
		return !connected;
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#finalize()
	 */
	protected void finalize() throws Throwable {
		// TODO Auto-generated method stub
		super.finalize();
		disconnect();
	}



	public boolean send(SteeringData sd) {
		// TODO Auto-generated method stub
		try {
			dos.writeFloat(sd.getCtr_x());
			dos.writeFloat(sd.getCtr_y());
			dos.writeFloat(sd.getCtr_z());
			dos.writeFloat(sd.getLongitude());
			dos.writeFloat(sd.getLatitude());
			dos.writeFloat(sd.getZoom_factor());
			dos.writeFloat(sd.getVis_brightness());
			dos.writeFloat(sd.getVelocity_max());
			dos.writeFloat(sd.getStress_max());

			return true;
		} catch (Exception e) {
			
		}
		return false;
	}
	
}
