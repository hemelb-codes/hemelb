package uk.ac.ucl.chem.ccs.vizclient;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.net.InetSocketAddress;
import java.net.Socket;
import java.net.UnknownHostException;

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
	//private long frame_no = 0;
	private boolean connected = false;
	private String hostname;
	private int port, window;
	
	int frameno = 0;
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
			int width = 0;
			int height = 0;
			int frame_size=0;

			//frame_no++;		
			short x_data=0, y_data=0;

			//Calendar cal = Calendar.getInstance();
			//long start_time = cal.getTimeInMillis();
			//System.err.println("Frame no " + frame_no);

			try {			
				//System.err.println("Waiting for frame " + ++frameno);
				width = d.readInt();
				height = d.readInt();
				frame_size = d.readInt();
				//System.err.println("Got frame dimensions (" + width + ", " + height + ")");
				//System.err.println("Got frame size " + frame_size);
				if (frame_size < 1) {
					return null;
				}
			} catch (IOException e) {
				connected = false;
				return null;
			}

			int col_pixels = frame_size/BYTES_PER_PIXEL_DATA;
			VizFrameData vizFrame = new VizFrameData(width, height);

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

				int pixel = (y_data * width) + x_data;
				vizFrame.setPixel(pixel, 0,
						(colour_data[0] >> (3*bits_per_char)) & colour_mask,
						(colour_data[0] >> (2*bits_per_char)) & colour_mask,
						(colour_data[0] >> (1*bits_per_char)) & colour_mask);
				vizFrame.setPixel(pixel, 1,
						(colour_data[0] >> (0*bits_per_char)) & colour_mask,
						(colour_data[1] >> (3*bits_per_char)) & colour_mask,
						(colour_data[1] >> (2*bits_per_char)) & colour_mask);
				vizFrame.setPixel(pixel, 2,
						(colour_data[1] >> (1*bits_per_char)) & colour_mask,
						(colour_data[1] >> (0*bits_per_char)) & colour_mask,
						(colour_data[2] >> (3*bits_per_char)) & colour_mask);
				vizFrame.setPixel(pixel, 3,
						(colour_data[2] >> (2*bits_per_char)) & colour_mask,
						(colour_data[2] >> (1*bits_per_char)) & colour_mask,
						(colour_data[2] >> (0*bits_per_char)) & colour_mask);
			}

			// PLZ CAN HAZ DATEZ PLS! k thnx bi.
			//Calendar cal2 = Calendar.getInstance();
			//long end_time = cal2.getTimeInMillis();

			//double total_time = (end_time - start_time)/1000.f;
			//System.out.println("Total time "+ total_time);
			//double data_rate = (frame_size/1024.f) / total_time;

			//System.err.println("bytes = " + frame_size + ", time = " + total_time + ", rate = " + data_rate + "KB/s");

			try {
				
				vizFrame.setVis_time_step(d.readInt());
				vizFrame.setVis_time(d.readDouble());
				vizFrame.setVis_cycle(d.readInt());

				//read inlet vels
				int n_inlets = d.readInt();
			
				vizFrame.setN_inlets(n_inlets);
/*			for (int n=0; n< n_inlets; n++) {
				double inlet_avg_vel = d.readDouble();
				vizFrame.setInlet_avg_vel(inlet_avg_vel, n);
			}*/
			
				//mouse parameters
				vizFrame.setVis_mouse_pressure(d.readDouble());
				vizFrame.setVis_stess_pressure(d.readDouble());

				//System.err.println(vizFrame.getVis_pressure_min() + " : " + vizFrame.getVis_stess_pressure());

				//System.err.println("Finished getting frame");
			} catch (IOException e1) {
				connected = false;
				return null;
			}

			//vizFrame.setDataRate(0);
			//vizFrame.setRealFrameNo(0);
			//vizFrame.setFrameNo(frame_no);
			vizFrame.setBufferSize(frame_size);
			//vizFrame.setFramePerSec(1/total_time);
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
		// Write params in right order
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
			dos.writeFloat(sd.getPressure_min());
			dos.writeFloat(sd.getPressure_max());
			dos.writeFloat(sd.getVis_glyph_length());
			dos.writeFloat(sd.getPixels_x()*1.f);
			dos.writeFloat(sd.getPixels_y()*1.f);
			dos.writeFloat(sd.getVis_mouse_x()*1.f);
			dos.writeFloat(sd.getVis_mouse_y()*1.f);
			dos.writeFloat(sd.getKill()*1.f);
			dos.writeFloat(sd.getVis_mode()*1.f);
			dos.writeFloat(sd.getVis_streaklines_per_pulsatile_period());
			dos.writeFloat(sd.getVis_streakline_length());
			dos.writeFloat(10.0f); // Max framerate
			return true;
		} catch (Exception e) {

		}
		return false;
	}

}
