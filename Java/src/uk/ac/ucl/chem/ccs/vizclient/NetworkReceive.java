package uk.ac.ucl.chem.ccs.vizclient;

import java.io.*;
import java.net.*;
import java.nio.ByteBuffer;

public class NetworkReceive {

	private static final int BYTES_PER_PIXEL_DATA = 8;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

		//set hostname and port
		int port = Integer.parseInt(args[1]);
		String hostname = args[0];

		Socket listenSocket = null;
		DataInputStream in = null;


		try {
			listenSocket = new Socket(hostname, port);
			in = new DataInputStream(listenSocket.getInputStream());
		} catch (UnknownHostException e) {
			System.err.println("can't connect to host: " + hostname);
			System.exit(1);
		} catch (IOException e) {
			System.err.println("Couldn't get I/O for the connection to: " + hostname);
			System.exit(1);
		}

		readData(in);


		try {
			in.close();
			listenSocket.close();
		} catch (Exception e) {
			System.err.println("couldn't close connection");
			e.printStackTrace();
		}
	}
	
	public static boolean readData (DataInputStream d) {
		
		int frame_size=0;
		
		try {
			frame_size = d.readInt();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//check frame size
		if (frame_size == -1) {
			return false;
		}
		
		int col_pixels = frame_size/BYTES_PER_PIXEL_DATA;
	
		System.err.println("Buffer size is " + frame_size);
	
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

		
			// data per pixel are colour id and pixel id (2 * sizeof(int) * bytes)

			int bits_per_char = 8;
			int bits_per_two_chars = 2 * bits_per_char;

			// r-, g-, b- pixel components need "bits_per_char" bits
			int colour_mask = (1 << bits_per_char) - 1;

			// x-, y- pixel components need "bits_per_two_chars" bits
			int pixel_mask = (1 << bits_per_two_chars) - 1;
			
			int r =((colour_data >> bits_per_two_chars) & colour_mask);
			int g = ((colour_data >> bits_per_char     ) & colour_mask);
			int b =  ((colour_data                      ) & colour_mask);
		  
		     int x = (pixel_data >> bits_per_two_chars) & pixel_mask;
		     int y = (pixel_data                      ) & pixel_mask;
			
			System.out.println("Colour r = " + r + " g = " + g + " b = " + b);
			System.out.println("Pixel x = " + x + " y = " + y);
			
		}
		
		return true;
		
	}
	
}
