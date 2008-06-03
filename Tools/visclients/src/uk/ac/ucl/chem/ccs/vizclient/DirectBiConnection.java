package uk.ac.ucl.chem.ccs.vizclient;

import java.io.*;
import java.net.*;
import java.util.Calendar;
public class DirectBiConnection implements SteeringConnection {

	private static final int BYTES_PER_PIXEL_DATA = 8;
	private DataInputStream d;
	private DataOutputStream dos;
	private Socket listenSocket;	
	
	// data per pixel are colour id and pixel id (2 * sizeof(int) * bytes)
	private int bits_per_char = 8;
	private int bits_per_two_chars = 2 * bits_per_char;
	// r-, g-, b- pixel components need "bits_per_char" bits
	private int colour_mask = (1 << bits_per_char) - 1;
	// x-, y- pixel components need "bits_per_two_chars" bits
	//private int pixel_mask = (1 << bits_per_two_chars) - 1;		
	private int colour_data = 0;
	//private int pixel_data = 0;	
	private long frame_no = 0;
	private boolean connected = false;
	private String hostname;
	private int port;
	/**
	 * 
	 */
	public DirectBiConnection(int port, String hostname) {
		// TODO Auto-generated method stub
		this.hostname = hostname;
		this.port = port;		
	}

	
	public VizFrameData getFrame () {
		//check if connected?
		if (connected) {
		int frame_size=0;
		frame_no++;		
		Calendar cal = Calendar.getInstance();
		short x_data=0, y_data=0;
		
		long start_time = cal.getTimeInMillis();
		System.err.println("Frame no " + frame_no);

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
					colour_data = d.readInt();
					x_data = d.readShort();
					y_data = d.readShort();
				} catch (IOException e1) {
					connected = false;
					return null;
				}

			vizFrame.setR(i,(colour_data >> bits_per_two_chars) & colour_mask);
			vizFrame.setG(i,(colour_data >> bits_per_char     ) & colour_mask);
			vizFrame.setB(i,(colour_data                      ) & colour_mask);
	    
			
			vizFrame.setX(i, x_data);
			vizFrame.setY(i, y_data);	    
		}

		Calendar cal2 = Calendar.getInstance();
		// PLZ CAN HAZ DATEZ PLS! k thnx bi.
		long end_time = cal2.getTimeInMillis();
		double total_time = (end_time - start_time)*1000;
		double data_rate = 1024.0*(frame_size / total_time);
				
		System.err.println("bytes = " + frame_size + ", time = " + total_time +
				", rate = " + data_rate + "KB/s");
		vizFrame.setDataRate(data_rate);
		vizFrame.setRealFrameNo(0);
		vizFrame.setFrameNo(frame_no);
		vizFrame.setBufferSize(frame_size);
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
			listenSocket.setReceiveBufferSize(1024*1024);
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

	
	public boolean magnify(int m) {
		// TODO Auto-generated method stub
		try {
			dos.writeInt(2);
			dos.writeInt(m);
			return true;
		} catch (Exception e) {
			
		}
		return false;
	}

	public boolean rotate(double dx, double dy) {
		// TODO Auto-generated method stub
		try {
			dos.writeInt(1);
			dos.writeDouble(dx);
			dos.writeDouble(dy);
			return true;
		} catch (Exception e) {
			
		}
		return false;
	}
	
}
