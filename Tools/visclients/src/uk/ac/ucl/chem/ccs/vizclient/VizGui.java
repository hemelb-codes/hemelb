package uk.ac.ucl.chem.ccs.vizclient;

import info.clearthought.layout.TableLayout;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import com.sun.opengl.util.Animator;
import javax.swing.BorderFactory;
import java.awt.event.MouseWheelListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseListener;
import java.awt.Point;
import java.awt.event.MouseEvent;
import javax.swing.WindowConstants;
import javax.swing.border.BevelBorder;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.JScrollPane;
import javax.swing.JTextField;

import java.util.concurrent.ConcurrentLinkedQueue;
import javax.media.opengl.*;

public class VizGui extends javax.swing.JPanel implements GLEventListener{
	private JTextArea notificationArea;
	static private JScrollPane jScrollPane1;
	private GLCanvas canvas1;
	private GLCapabilities cap;
	float scale_x;
	float scale_y;
	private String hostname;
	private int port;
    private NetThread thread=null; 
	private Animator animator = null;
    private boolean connected = false;
	private ConcurrentLinkedQueue queue;
	private SteeringConnection nr;
	private InfoPanel ifp = null;

	
	public VizGui(int port, String hostname) {
		super();
		queue = new ConcurrentLinkedQueue();
		this.port = port;
		this.hostname = hostname;
		initGUI();
		nr = new DirectBiConnection (port, hostname);	
	}

	public VizGui () {
		this(1, "foohost");
	}
	
	public void setHostPort (int port, String hostname) {
		this.port = port;
		this.hostname = hostname;
		notificationArea.append("Changed host " + hostname + ":" + port + "\n");
		  notificationArea.setCaretPosition(notificationArea.getDocument().getLength());
	}
	
	public boolean startReceive() {
		if (!connected) {
		nr = new DirectBiConnection (port, hostname);
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

		//start the animator thread
		animator = new Animator(canvas1);
		animator.start();
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
	
	int y = 1024;
	int x = 1024;
	scale_x = 1.0f / (float)x;
	scale_y = 1.0f / (float)y;
	
	
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
    	
    	
   	 // float ortho_x = 0.5F * (float)w / (float)PIXELS_X;
      //float ortho_y = 0.5F * (float)h / (float)PIXELS_Y;
    	  
    //gl.glViewport(0, 0, width, height);
    	  
    	  //gl.glLoadIdentity ();
    	 // gl.glOrtho2D(-ortho_x, ortho_x, -ortho_y, ortho_y);*/
    }

    /** This method handles the painting of the GLDrawable */
    
    public void display(GLAutoDrawable drawable) {
    VizFrameData vfd = (VizFrameData)queue.poll();
    
    if (vfd != null) {	
    GL gl = drawable.getGL();
	/** Clear the colour buffer */
	gl.glClear( GL.GL_COLOR_BUFFER_BIT );
	/** Draw some dots */ 

	gl.glBegin( GL.GL_POINTS );	
	  for (int i = 0; i < vfd.getLength(); i++)
	    {

	      gl.glColor3f (vfd.getR(i) * (1.0f / 255.0f),
	      				vfd.getG(i) * (1.0f / 255.0f),
						vfd.getB(i) * (1.0f / 255.0f));
	      
	  gl.glVertex2f (-0.5f + scale_x * vfd.getX(i),-0.5f + scale_y * vfd.getY(i));
	    
	    }
	  gl.glEnd(); 
	  
	  if (ifp != null) {
		  ifp.updatePanel(vfd.getBufferSize(), vfd.getRealFrameNo(), vfd.getFrameNo(), vfd.getDataRate());
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
							TableLayout.FILL } });
			this.setLayout(thisLayout);
			thisLayout.setHGap(5);
			thisLayout.setVGap(5);
			this.setPreferredSize(new java.awt.Dimension(570, 677));
			{
				cap = new GLCapabilities(); 
				canvas1 = new GLCanvas(cap);
				this.add(canvas1, "0, 1, 0, 9");
				canvas1.addGLEventListener(this);
				canvas1.addMouseWheelListener(new MouseWheelListener () {
					
				    public void mouseWheelMoved(MouseWheelEvent e) {
				        String message;
				        String newline = "\n";
				        int notches = e.getWheelRotation();
				        if (notches < 0) {
				            message = "Mouse wheel moved UP "
				                         + -notches + " notch(es)" + newline;
				        } else {
				            message = "Mouse wheel moved DOWN "
				                         + notches + " notch(es)" + newline;
				        }
				        System.err.print(message);
				     }
					
				});
				canvas1.addMouseListener(new MouseListener () {
					Point start = null;
					
					public void mouseReleased(MouseEvent e) {
						if (start != null) {
							Point end = e.getPoint();
							double dx = start.getX() - end.getX();
							double dy = start.getY() - end.getY();
							
							if (dx != 0 && dy != 0) {
								//scale 
								dx = dx/50.d;
								dy = dy/50.d;
								
							if (nr != null && nr.isConnected()) {
								System.err.println (" dX = "  +  dx + " dy = " +dy);

							nr.rotate(dx, dy);
							}
							}
						}
					}
					
					public void mouseClicked(MouseEvent e) {
					}
					
					public void mousePressed(MouseEvent e) {
						start = e.getPoint();
					}
					
					public void mouseExited(MouseEvent e) {
					}
					
					public void mouseEntered(MouseEvent e) {
					}
				});
			}	
			{
				jScrollPane1 = new JScrollPane();
				this.add(jScrollPane1, "0, 0, 0, 0");
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
	
	public boolean isConnected() {
		return nr.isConnected();
	}
	
	public JPanel getInfoPanel () {
		ifp = new InfoPanel();
		return ifp;
	}
	
    private class NetThread extends Thread {
    	


		public NetThread () {
    
    	}
    	
    	public void stopThread () {
    		nr.disconnect();
    	}
    	

    	
    	public void run() {
    		while(nr.isConnected()) {
    			VizFrameData vfd =nr.getFrame();
    			if (vfd != null) {
    			queue.offer(vfd);
    		}
    		}
    			stopReceive();

    	}
    }
    
    private class InfoPanel extends JPanel {
    	
    	private JLabel jLabel4;
    	private JTextField dataRate;
    	private JLabel jLabel5;
    	private JTextField droppedFrames;
    	private JTextField frameSize;
    	private JTextField frameRecNo;
    	private JTextField frameNo;
    	private JLabel jLabel3;
    	private JLabel jLabel2;
    	private JLabel jLabel1;
    	
    	public InfoPanel () {
    		super();
   
		TableLayout infoPanelLayout = new TableLayout(new double[][] {{TableLayout.FILL}, {TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}});
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
			jLabel3.setText("Frame Size");
			jLabel3.setFont(new java.awt.Font("Dialog",1,12));
		}
		{
			jLabel4 = new JLabel();
			this.add(jLabel4, "0, 8");
			jLabel4.setText("Data Rate");
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
			dataRate = new JTextField();
			this.add(dataRate, "0, 9");
		}
    	}
    
    	public void updatePanel (int fSize, long fNo, long fRecNo, double dRate) {
    		frameNo.setText(Long.toString(fNo));
    		frameRecNo.setText(Long.toString(fRecNo));
    		frameSize.setText(Integer.toString(fSize));
    		droppedFrames.setText(Integer.toString(0));
    		dataRate.setText(Double.toString(dRate));
    	}
    	
    }

}
