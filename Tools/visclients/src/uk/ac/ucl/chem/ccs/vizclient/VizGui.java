package uk.ac.ucl.chem.ccs.vizclient;

import info.clearthought.layout.TableLayout;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import com.sun.opengl.util.Animator;
import javax.swing.BorderFactory;
import java.awt.event.MouseWheelListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.event.MouseListener;
import java.awt.Point;
import java.awt.event.MouseEvent;
import javax.swing.WindowConstants;
import javax.swing.border.BevelBorder;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JButton;
import javax.swing.JTextArea;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import java.text.DecimalFormat;
import java.util.concurrent.ConcurrentLinkedQueue;
import javax.media.opengl.*;
import java.awt.event.KeyListener;
import java.awt.event.KeyEvent;


public class VizGui extends javax.swing.JPanel implements GLEventListener{
	private JTextArea notificationArea;
	static private JScrollPane jScrollPane1;
	private GLCanvas canvas1;
	private GLCapabilities cap;
	float scale_x;
	float scale_y;
	private String hostname;
	private int port, window;
    private NetThread thread=null; 
    private TimerThread timeFramesPerSec;
    private RotateThread rotateModel;
	private Animator animator = null;
    private boolean connected = false;
	private ConcurrentLinkedQueue queue;
	private SteeringConnection nr;
	private InfoPanel ifp = null;	
	private SteeringData sd;
	private int SCALE_X = 1024;
	private int SCALE_Y = 1024;
	
	private long framesThisSec = 0;
	private long framesLastSec = 0;
	private long framesPerSec = 0;
	
	private double bytesRec = 0.d;
	
	private double panel_width;
	private double panel_height;

	
	private int scale = 512;
	
	static private final int VIEW1 =0;
	static private final int VIEW2 =1;
	static private final int VIEW3 =2;
	static private final int VIEW4 =3;
	static private final int VIEWALL =4;
	private int view;
	
	public VizGui(int port, String hostname, int window) {
		super();
		queue = new ConcurrentLinkedQueue();
		this.port = port;
		this.hostname = hostname;
		this.window = window;
		initGUI();
		nr = new DirectBiConnection (port, hostname, window);
		sd = new SteeringData();
		view = VIEW1;
	}

	public VizGui () {
		this(1, "foohost", 1024*1024);
	}
	
	public void setHostPort (int port, String hostname) {
		this.port = port;
		this.hostname = hostname;
		notificationArea.append("Changed host " + hostname + ":" + port + "\n");
		  notificationArea.setCaretPosition(notificationArea.getDocument().getLength());
	}
	
	public boolean startReceive() {
		if (!connected) {
		nr = new DirectBiConnection (port, hostname, window);
		connected = nr.connect();
		if (!connected) {
			notificationArea.append("Connection error: Couldn't connect to host " + hostname + ":" + port +"\n");
		 notificationArea.setCaretPosition(notificationArea.getDocument().getLength());
			return connected;
		}
		notificationArea.append("Connection started to host " + hostname + ":" + port +"\n");
		  notificationArea.setCaretPosition(notificationArea.getDocument().getLength());
		thread = new NetThread();
	    thread.start();
	    
	    //start the timing thread
	    // timeFramesPerSec = new TimerThread();
	    // timeFramesPerSec.start();
	   
	    
		//start the animator thread
		animator = new Animator(canvas1);
		animator.start();
		ifp.updateSteeredParams();
		ifp.viewChanged();
		send();
		}
		return connected;
	}
	
	public boolean stopReceive() {
		if (connected){
			animator.stop();
			connected = !nr.disconnect();
			
			//thread.stopThread();
			notificationArea.append("Connection terminated to host " + hostname + ":" + port +"\n");
			  notificationArea.setCaretPosition(notificationArea.getDocument().getLength());
		}
		return !connected;
	}
	
	
	/** 
     * Executed exactly once to initialize the 
     * associated GLDrawable
     */ 

    public void init(GLAutoDrawable drawable) {

        // print every openGL call for debugging purposes
	//drawable.setGL(new TraceGL(drawable.getGL(), System.err));
	GL gl = drawable.getGL();
	
	gl.glDisable (GL.GL_DEPTH_TEST);
	gl.glDisable (GL.GL_BLEND);
	gl.glShadeModel (GL.GL_FLAT);
	gl.glDisable (GL.GL_DITHER);

	gl.glClearColor( 1.0f, 1.0f, 1.0f, 1.0f ); //white

	gl.glPointSize (1.F);
	
	int y = 512;
	int x = 512;
	scale_x = 1f / (float)x;
	scale_y = 1f / (float)y;
	
	
    }
    
    /** 
     * Executed if the associated GLDrawable is resized
     */ 
    public void reshape(GLAutoDrawable drawable, int x, int y, int width, int height) {
	GL gl = drawable.getGL();
	gl.glViewport( 0, 0, width, height );
	gl.glMatrixMode( GL.GL_PROJECTION ); 
	gl.glLoadIdentity();
   //  
	//gl.glOrtho(0, 400, 0, 300, -1, 1);
    	
    	
   	  float ortho_x = 0.5F * (float)width / (float)512;
      float ortho_y = 0.5F * (float)height / (float)512;
    	  
    //gl.glViewport(0, 0, width, height);
    	  
    	  //gl.glLoadIdentity ();
    	  gl.glOrtho(-ortho_x, ortho_x, -ortho_y, ortho_y, -1, 1);
    }

    /** This method handles the painting of the GLDrawable */
    
    public void display(GLAutoDrawable drawable) {
    VizFrameData vfd = (VizFrameData)queue.poll();
    
    if (vfd != null) {	
    GL gl = drawable.getGL();
	/** Clear the colour buffer */
	gl.glClear( GL.GL_COLOR_BUFFER_BIT );
	/** Draw some dots */ 

	if (view < VIEWALL) {
	
	gl.glBegin( GL.GL_POINTS );	
	  for (int i = 0; i < vfd.getLength(); i++)
	    {

	      gl.glColor3f (vfd.getR(i, view) * (1.0f / 255.0f),
	      				vfd.getG(i, view) * (1.0f / 255.0f),
						vfd.getB(i, view) * (1.0f / 255.0f));
	      
	      float a = -0.5f + (vfd.getX(i)/512f);
	      float b = -0.5f + (vfd.getY(i)/512f);
	      
	    //  System.err.println ("a = " + a + " b = " + b);
	      
	  gl.glVertex2f (a,b);
	    
	  
	  
	    }
	  gl.glEnd(); 
	  
	} else {
		gl.glBegin( GL.GL_POINTS );	

		
		  for (int i = 0; i < vfd.getLength(); i++)
		    {
//image 1
		      gl.glColor3f (vfd.getR(i, VIEW1) * (1.0f / 255.0f),
		      				vfd.getG(i, VIEW1) * (1.0f / 255.0f),
							vfd.getB(i, VIEW1) * (1.0f / 255.0f));
		      
		  gl.glVertex2f (-0.5f + scale_x * vfd.getX(i),-0.5f + scale_y * (vfd.getY(i)+ (SCALE_Y/2.0f)));
		  
		  
		//image 2
	      gl.glColor3f (vfd.getR(i, VIEW2) * (1.0f / 255.0f),
	      				vfd.getG(i, VIEW2) * (1.0f / 255.0f),
						vfd.getB(i, VIEW2) * (1.0f / 255.0f));
	      
	  gl.glVertex2f (-0.5f + scale_x * (vfd.getX(i)+ (SCALE_X/2.0f)),-0.5f + scale_y * (vfd.getY(i) + (SCALE_Y/2.0f)));
	  
	//image 3
      gl.glColor3f (vfd.getR(i, VIEW3) * (1.0f / 255.0f),
      				vfd.getG(i, VIEW3) * (1.0f / 255.0f),
					vfd.getB(i, VIEW3) * (1.0f / 255.0f));
      
      gl.glVertex2f (-0.5f + scale_x * vfd.getX(i) ,-0.5f + scale_y * vfd.getY(i));
  
  
  	//image 4
  	gl.glColor3f (vfd.getR(i, VIEW4) * (1.0f / 255.0f),
  				vfd.getG(i, VIEW4) * (1.0f / 255.0f),
				vfd.getB(i, VIEW4) * (1.0f / 255.0f));
  
  	gl.glVertex2f (-0.5f + scale_x * (vfd.getX(i)+ (SCALE_X/2.0f)),-0.5f + scale_y * vfd.getY(i));


		    
		    }
		  gl.glEnd();
	}
	  
	  if (ifp != null) {
		  ifp.updatePanel(vfd.getBufferSize(), vfd.getRealFrameNo(), vfd.getFrameNo());
		 // System.err.println(vfd.getFramePerSec());
	  }
	  
	 // notificationArea.append("Frame " + vfd.getFrameNo() + "  Buffer size " + vfd.getBufferSize() + "\n");
	 // notificationArea.setCaretPosition(notificationArea.getDocument().getLength());
	  
	
    }
    }
    /** This method handles things if display depth changes */
    public void displayChanged(GLAutoDrawable drawable,
			       boolean modeChanged,
			       boolean deviceChanged){
    }
	
	private void initGUI() {
		try {
			TableLayout thisLayout = new TableLayout(new double[][] {
					{ TableLayout.FILL },
					{ TableLayout.FILL, TableLayout.FILL, TableLayout.FILL,
							TableLayout.FILL, TableLayout.FILL,
							TableLayout.FILL, TableLayout.FILL,
							TableLayout.FILL, TableLayout.FILL,
							TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL,
							TableLayout.FILL, TableLayout.FILL,
							TableLayout.FILL, TableLayout.FILL,
							TableLayout.FILL, TableLayout.FILL,
							TableLayout.FILL } });
			this.setLayout(thisLayout);
			thisLayout.setHGap(5);
			thisLayout.setVGap(5);
			this.setPreferredSize(new java.awt.Dimension(570, 677));
			{
				cap = new GLCapabilities(); 
				canvas1 = new GLCanvas(cap);
				this.add(canvas1, "0, 3, 0, 19");
				canvas1.addGLEventListener(this);
			
				

				canvas1.addKeyListener(new KeyListener () {
					
					
				public void keyPressed(KeyEvent e) {
					
				}
				
				public void keyReleased(KeyEvent e)  {
					
				}
				
				public void keyTyped(KeyEvent e) {
					if (e.getKeyChar() == '1') {
						view = VIEW1;
						scale_x = 1.0f / (float)scale;
						scale_y = 1.0f / (float)scale;	
					} else if (e.getKeyChar() == '2') {
						view = VIEW2;
						scale_x = 1.0f / (float)scale;
						scale_y = 1.0f / (float)scale;	
					} else if (e.getKeyChar() == '3') {
						view = VIEW3;
						scale_x = 1.0f / (float)scale;
						scale_y = 1.0f / (float)scale;	
					} else if (e.getKeyChar() == '4') {
						view = VIEW4;
						scale_x = 1.0f / (float)scale;
						scale_y = 1.0f / (float)scale;	
					} else if (e.getKeyChar() == '5') {
						view = VIEWALL;

						scale_x = 1.0f / (float)SCALE_X;
						scale_y = 1.0f / (float)SCALE_Y;
					} else	if (e.getKeyChar() == 'b') {
						
						if (scale > 256) {
							scale = scale/2;
							scale_x = 1.0f / (float)scale;
							scale_y = 1.0f / (float)scale;	
						}
						
					} else if (e.getKeyChar() == 's') {
						if (scale < 1024) {
							scale = scale*2;
							scale_x = 1.0f / (float)scale;
							scale_y = 1.0f / (float)scale;	
						}					}
					ifp.viewChanged();
					//System.err.println("View " + view);
				}
				
				
				
				});
				
				canvas1.addMouseWheelListener(new MouseWheelListener () {
					
				    public void mouseWheelMoved(MouseWheelEvent e) {
				    	float mag = sd.getZoom_factor();	
				        int notches = (e.getWheelRotation());
				        mag = mag - (0.1f*notches);
				        sd.setZoom_factor(mag);
				        send();
				       // System.err.println("Mag foactor= " + sd.getZoom_factor());
				     }
					
				});
				canvas1.addMouseListener(new MouseListener () {
					Point start = null;
					
					public void mouseReleased(MouseEvent e) {
						//System.err.println("Button " + e.getButton());
						if (start != null && e.getButton() == 1) {
							Point end = e.getPoint();
							double dx = start.getX() - end.getX();
							double dy = start.getY() - end.getY();
							
							//dddd
							panel_width = canvas1.getWidth();
							panel_height = canvas1.getHeight();
							
							if (dx != 0) {
								//scale 
								dx = (90.d/panel_width)*dx;
								sd.updateLongitude((float)dx);
								//System.err.println("dx = " + dx + " panel " + panel_width);
							} if (dy != 0) {
								dy = -1*(90.d/panel_height)*dy;
								sd.updateLatitude((float)dy);
								//System.err.println("dy = " + dy);
							}
							send();
						} else if(start != null && e.getButton() == 3) {
							Point end = e.getPoint();
							double dx = start.getX() - end.getX();
							double dy = start.getY() - end.getY();
							panel_width = canvas1.getWidth();
							panel_height = canvas1.getHeight();
							
							sd.setCtr_x(100*(float)(dx/panel_width));
							sd.setCtr_y(100*(float)(dy/panel_height));
							
							System.err.println("Ctr x = " + sd.getCtr_x() + " crt y = " + sd.getCtr_y());
							send();
 						}
					}
					
					public void mouseClicked(MouseEvent e) {
					}
					
					public void mousePressed(MouseEvent e) {
						start = e.getPoint();
						//System.err.println("startx = " + start.getX() + " starty = " + start.getY());
					}
					
					public void mouseExited(MouseEvent e) {
					}
					
					public void mouseEntered(MouseEvent e) {
					}
				});
			}	
			{
				jScrollPane1 = new JScrollPane();
				this.add(jScrollPane1, "0, 0, 0, 2");
				jScrollPane1.setAutoscrolls(true);
				{
					notificationArea = new JTextArea();
					jScrollPane1.setViewportView(notificationArea);
					notificationArea.setText("-------Messages-------\n");
					notificationArea.setVisible(true);
					notificationArea.setBorder(BorderFactory
						.createBevelBorder(BevelBorder.LOWERED));
					notificationArea.setDoubleBuffered(true);
					notificationArea.setEditable(false);
				}
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void rotate (boolean rot) {
		if (rot) {
			//System.out.println("rotate");
			rotateModel = new RotateThread();
			rotateModel.start();
		} else {
			rotateModel.stopThread();
		}
	}
	
	public void resetSteering () {
		sd = new SteeringData();
		send();
		ifp.updateSteeredParams();
		view = VIEW1;
	}
	
	private void send() {
		if (nr != null && nr.isConnected()) {
			//System.err.println (" dX = "  +  dx + " dy = " +dy);
		nr.send(sd);
		}
	}
	
	public boolean isConnected() {
		return nr.isConnected();
	}
	
	public JPanel getInfoPanel () {
		ifp = new InfoPanel();
		return ifp;
	}
	
	//thread to rotate the model
	private class RotateThread extends Thread {
		boolean run;
		public RotateThread () {
			run=true;
		}
   
		
    	public void stopThread () {
    		//private boolean shouldirun = true;
    		run=false;
    	} 	
   	
    	public void run() {
    		while(nr.isConnected() && run==true) {
    			try {
    			Thread.sleep(250);    			
    			sd.updateLongitude(0.5f);
    			send();
    			
    			} catch (Exception e) {
    				
    			}
    		}
    			

    	}
	}
	
	//thread to time fps etc
	private class TimerThread extends Thread {
		private boolean run;
		
		public TimerThread () {
			run=true;
		}
   
		
    	public void stopThread () {
    		//private boolean shouldirun = true;
    		run=false;
    	} 	
   	
    	public void run() {
    		while(nr.isConnected()) {
    			try {
    			Thread.sleep(1000);    			
    			framesPerSec = (framesThisSec - framesLastSec);
    			framesLastSec = framesThisSec;
    			ifp.updateFPS(framesPerSec, bytesRec);
    			bytesRec = 0;
   			
    			} catch (Exception e) {
    				
    			}
    		}
    			

    	}
    
		
	}
	
	//thread to do the network receive
    private class NetThread extends Thread {
    	
    	//TimerThread timeFramesPerSec;

		public NetThread () {
			
    	}
    	
    	public void stopThread () {
    		nr.disconnect();
    		//timeFramesPerSec.stop();
    	}
    	

    	
    	public void run() {
    		//timeFramesPerSec.start();
    		
    		while(nr.isConnected()) {
    			VizFrameData vfd =nr.getFrame();
    			if (vfd != null) {
    			queue.offer(vfd);
    		//	framesThisSec=vfd.getFrameNo();
    		//	bytesRec = bytesRec + vfd.getBufferSize();
    		}
    		}
    			stopReceive();

    	}
    }
    
    //side info panel 
    private class InfoPanel extends JPanel {
    	
    	private JLabel jLabel4;
    	private JTextField dataRate;
    	private JLabel jLabel5;
    	private JLabel jLabel6;
    	private JTextField droppedFrames;
    	private JTextField viewing;
    	private JTextField frameSize;
    	private JTextField frameRecNo;
    	private JTextField frameNo;
    	private JLabel jLabel3;
    	private JLabel jLabel2;
    	private JLabel jLabel1;
    	private JTextField framePerSec;
    	private JLabel jLabel10;
    	
    	//steered params
    	private JLabel jLabel7;
    	private JLabel jLabel8;
    	private JLabel jLabel9;
    	private JTextField vizBrightness;
    	private JTextField velMax;
    	private JTextField stressMax;
    	private JButton updateParams;
    	

    	
    	
    	
    	public InfoPanel () {
    		super();
   
		TableLayout infoPanelLayout = new TableLayout(new double[][] {{TableLayout.FILL}, {TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}});
		infoPanelLayout.setHGap(5);
		infoPanelLayout.setVGap(5);
		this.setLayout(infoPanelLayout);
		this.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
		{
			jLabel1 = new JLabel();
			this.add(jLabel1, "0, 2");
			jLabel1.setText("Frame Rec No");
			jLabel1.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 0));
			jLabel1.setFont(new java.awt.Font("Dialog",1,12));
		}
		{
			jLabel2 = new JLabel();
			this.add(jLabel2, "0, 0");
			jLabel2.setText("Frame No");
			jLabel2.setFont(new java.awt.Font("Dialog",1,12));
		}
		{
			jLabel3 = new JLabel();
			this.add(jLabel3, "0, 4");
			jLabel3.setText("Frame Size (bits)");
			jLabel3.setFont(new java.awt.Font("Dialog",1,12));
		}
		{
			jLabel4 = new JLabel();
			this.add(jLabel4, "0, 8");
			jLabel4.setText("Data Rate (Kbits/sec)");
			jLabel4.setFont(new java.awt.Font("Dialog",1,12));
		}
		{
			frameNo = new JTextField();
			this.add(frameNo, "0, 1");
		}
		{
			frameRecNo = new JTextField();
			this.add(frameRecNo, "0, 3");
			frameRecNo.setText("");
		}
		{
			frameSize = new JTextField();
			this.add(frameSize, "0, 5");
			frameSize.setText("");
		}
		{
			droppedFrames = new JTextField();
			this.add(droppedFrames, "0, 7");
			droppedFrames.setText("");
		}
		{
			jLabel5 = new JLabel();
			this.add(jLabel5, "0, 6");
			jLabel5.setText("Frames Dropped");
			jLabel5.setFont(new java.awt.Font("Dialog",1,12));
		}
		{
			jLabel6 = new JLabel();
			this.add(jLabel6, "0, 10");
			jLabel6.setText("View");
			jLabel6.setFont(new java.awt.Font("Dialog",1,12));
		}
		{
			dataRate = new JTextField();
			this.add(dataRate, "0, 9");
		}
		{
			viewing = new JTextField();
			this.add(viewing, "0, 11");			
			viewing.setText("");
		}
		{
			framePerSec = new JTextField();
			this.add(framePerSec, "0, 13");
		}
		{
			jLabel10 = new JLabel();
			this.add(jLabel10, "0, 12");
			jLabel10.setText("Frames Per Second");
			jLabel10.setFont(new java.awt.Font("Dialog",1,12));
		}
		
		{
			jLabel7 = new JLabel();
			this.add(jLabel7, "0, 15");
			jLabel7.setText("Viz Brightness");
			jLabel7.setFont(new java.awt.Font("Dialog",1,12));
		}
		{
			vizBrightness = new JTextField();
			this.add(vizBrightness, "0, 16");
			//vizBrightness.setText(Float.toString(sd.getVis_brightness()));
		}
		{
			jLabel8 = new JLabel();
			this.add(jLabel8, "0, 17");
			jLabel8.setText("Velocity Max (m/s)");
			jLabel8.setFont(new java.awt.Font("Dialog",1,12));
		}
		{
			velMax = new JTextField();
			this.add(velMax, "0, 18");
			//velMax.setText(Float.toString(sd.getVelocity_max()));
		}
		{
			jLabel9 = new JLabel();
			this.add(jLabel9, "0, 19");
			jLabel9.setText("Stress Max (Pa)");
			jLabel9.setFont(new java.awt.Font("Dialog",1,12));
		}
		{
			stressMax = new JTextField();
			this.add(stressMax, "0, 20");
			//stressMax.setText(Float.toString(sd.getStress_max()));
		}
		{
			updateParams = new JButton("Update");
			this.add(updateParams, "0, 22");
			updateParams.addActionListener(new ActionListener () {
				public void actionPerformed(ActionEvent evt) {
				
					try {
						sd.setStress_max(Float.parseFloat(stressMax.getText()));
						sd.setVelocity_max(Float.parseFloat(velMax.getText()));
						sd.setVis_brightness(Float.parseFloat(vizBrightness.getText()));
					} catch (Exception e) {
						
					}
					
					send();
				}
			});
		}
		
    	}
    
    	public void updateSteeredParams() {
    		stressMax.setText(Float.toString(sd.getStress_max()));
    		velMax.setText(Float.toString(sd.getVelocity_max()));
    		vizBrightness.setText(Float.toString(sd.getVis_brightness()));
    	}
    	
    	public void updateFPS (long fps, double bps) {
    		framePerSec.setText(Long.toString(fps));
    		DecimalFormat df= new DecimalFormat("#####0.###"); 
    		dataRate.setText(df.format(bps));
    	}
    	
    	public void viewChanged () {
    		switch (view) {
    			case 0:
    				viewing.setText("Velocity Volume Rendering");
    				break;	
   				case 1:
   					viewing.setText("Stress Volume Rendering");
    				break;
    			case 2:
    				viewing.setText("Wall Pressure");
    				break;
    			case 3:
    				viewing.setText("Wall Stress");
    				break;
    			case 4:
    				viewing.setText("All Views");
    				break;
    		}
    		
    	}
    	
    	public void updatePanel (int fSize, long fNo, long fRecNo) {
    		frameNo.setText(Long.toString(fNo));
    		frameRecNo.setText(Long.toString(fRecNo));
    		frameSize.setText(Integer.toString(fSize));
    		droppedFrames.setText(Integer.toString(0));

    	

    	}
    	
    }

}
