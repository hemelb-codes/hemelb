package uk.ac.ucl.chem.ccs.vizclient;

import java.io.*;
import java.net.*;
import java.nio.ByteBuffer;
import java.util.Date;
import java.util.Calendar;
public class NetworkReceive  {

	private static final int BYTES_PER_PIXEL_DATA = 8;
	private DataInputStream d;
	private Socket listenSocket;	
	
	// data per pixel are colour id and pixel id (2 * sizeof(int) * bytes)
	private int bits_per_char = 8;
	private int bits_per_two_chars = 2 * bits_per_char;
	// r-, g-, b- pixel components need "bits_per_char" bits
	private int colour_mask = (1 << bits_per_char) - 1;
	// x-, y- pixel components need "bits_per_two_chars" bits
	private int pixel_mask = (1 << bits_per_two_chars) - 1;		
	private int colour_data = 0;
	private int pixel_data = 0;	
	private long frame_no = 0;
	
	/**
	 * 
	 */
	public NetworkReceive(int port, String hostname) {
		// TODO Auto-generated method stub
		int s=0;
		try {
			listenSocket = new Socket();
			//set TCP buffer size
			listenSocket.setReceiveBufferSize(1024*1024);
			s=listenSocket.getReceiveBufferSize();
			listenSocket.connect(new InetSocketAddress(hostname, port));
			
       		BufferedInputStream bufff = new BufferedInputStream(listenSocket.getInputStream());
			d = new DataInputStream(bufff);

		} catch (UnknownHostException e) {
			System.err.println("can't connect to host: " + hostname);
			System.exit(1);
		} catch (IOException e) {
			System.err.println("Couldn't get I/O for the connection to: " + hostname);
			System.exit(1);
		}

		System.err.println("RecSize " + s);
		
	}

	
	public VizFrameData getFrame () {
		
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
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
	
		int col_pixels = frame_size/BYTES_PER_PIXEL_DATA;
		VizFrameData vizFrame = new VizFrameData(col_pixels);

		for (int i =0; i < col_pixels; i++) {

			try {
					colour_data = d.readInt();
					x_data = d.readShort();
					y_data = d.readShort();
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}

			vizFrame.setR(i,(colour_data >> bits_per_two_chars) & colour_mask);
			vizFrame.setG(i,(colour_data >> bits_per_char     ) & colour_mask);
			vizFrame.setB(i,(colour_data                      ) & colour_mask);
		  
			//vizFrame.setX(i, (pixel_data >> bits_per_two_chars) & pixel_mask);
			//vizFrame.setY(i, (pixel_data                   ) & pixel_mask);	    
			
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
		
		vizFrame.setFrameNo(frame_no);
		vizFrame.setBufferSize(frame_size);
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
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		NetworkReceive nr = new NetworkReceive (Integer.parseInt(args[1]), args[0]);
		while (true) {
		nr.getFrame();
		}
	}
	
	public static final int byteArrayToInt(byte [] b) {
       // return (b[3] << 24)
        //        + ((b[2] & 0xFF) << 16)
          //      + ((b[1] & 0xFF) << 8)
            //    + (b[0] & 0xFF);
		
			int result = b[3] + 256*(b[2] + 256*(b[1] + 256*b[0]));
	    //System.err.println("b0 " + b[0] + " b1 " + b[1] + " b2 " +  b[2] + " b3 " + b[3]);
    return result;

}

	
}
