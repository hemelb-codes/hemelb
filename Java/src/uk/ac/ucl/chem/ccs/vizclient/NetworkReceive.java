package uk.ac.ucl.chem.ccs.vizclient;

import java.io.*;
import java.net.*;
import java.awt.*;
import java.awt.image.*;

public class NetworkReceive {

	private static final int BYTES_PER_PIXEL_DATA = 8;
	private DataInputStream d;
	private Socket listenSocket;	
	
	private int size_x = 1024;
	private int size_y = 1024;
	
	/**
	 * 
	 */
	public NetworkReceive(int port, String hostname) {
		// TODO Auto-generated method stub

		
		//DataInputStream d = null;


		try {
			listenSocket = new Socket(hostname, port);
			d = new DataInputStream(listenSocket.getInputStream());
		} catch (UnknownHostException e) {
			System.err.println("can't connect to host: " + hostname);
			System.exit(1);
		} catch (IOException e) {
			System.err.println("Couldn't get I/O for the connection to: " + hostname);
			System.exit(1);
		}


		
	}

	
	public VizFrameData getFrame () {
		
		int frame_size=0;
		
		try {
			frame_size = d.readInt();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//check frame size
		if (frame_size == -1) {
			return null;
		}
		
		int col_pixels = frame_size/BYTES_PER_PIXEL_DATA;
	
		System.err.println("Buffer size is " + frame_size);
	
		VizFrameData vizFrame = new VizFrameData(col_pixels);
		
		// data per pixel are colour id and pixel id (2 * sizeof(int) * bytes)

		int bits_per_char = 8;
		int bits_per_two_chars = 2 * bits_per_char;

		// r-, g-, b- pixel components need "bits_per_char" bits
		int colour_mask = (1 << bits_per_char) - 1;

		// x-, y- pixel components need "bits_per_two_chars" bits
		int pixel_mask = (1 << bits_per_two_chars) - 1;
		
		
		int colour_data = 0;
		int pixel_data = 0;

		
		for (int i = 0; i < col_pixels; i++) {
			try {
				colour_data = d.readInt();
				pixel_data = d.readInt();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			

			
			vizFrame.setR(i,(colour_data >> bits_per_two_chars) & colour_mask);
			vizFrame.setG(i,(colour_data >> bits_per_char     ) & colour_mask);
			vizFrame.setB(i,(colour_data                      ) & colour_mask);
		  
			vizFrame.setX(i, (pixel_data >> bits_per_two_chars) & pixel_mask);
			vizFrame.setY(i, (pixel_data                      ) & pixel_mask);		    
		}
		return vizFrame;

	}
	
	
	/* (non-Javadoc)
	 * @see java.lang.Object#finalize()
	 */
	protected void finalize() throws Throwable {
		// TODO Auto-generated method stub
		super.finalize();
		try {
			d.close();
			listenSocket.close();
		} catch (Exception e) {
			System.err.println("couldn't close connection");
			e.printStackTrace();
		}
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		NetworkReceive nr = new NetworkReceive (Integer.parseInt(args[1]), args[0]);
		nr.getFrame();
	}
	
}
