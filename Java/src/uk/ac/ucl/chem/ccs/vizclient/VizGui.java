package uk.ac.ucl.chem.ccs.vizclient;

import info.clearthought.layout.TableLayout;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import com.sun.opengl.util.Animator;
import javax.swing.BorderFactory;

import javax.swing.WindowConstants;
import javax.swing.border.BevelBorder;
import javax.swing.JFrame;
import javax.swing.JTextPane;
import java.util.concurrent.ConcurrentLinkedQueue;
import javax.media.opengl.*;


public class VizGui extends javax.swing.JPanel implements GLEventListener{
	private JTextPane jTextPane1;
	private GLCanvas canvas1;
	private GLCapabilities cap;

	float scale_x = 1.0f / (float)1024;
	float scale_y = 1.0f / (float)1024;
	
	ConcurrentLinkedQueue queue;
	
	public VizGui(int port, String hostname) {
		super();
		queue = new ConcurrentLinkedQueue();
	    NetThread nt = new NetThread(port, hostname);
	    Thread thread = new Thread(nt);
	    thread.start();
		initGUI();
	}
	
	/**
	* Auto-generated main method to display this 
	* JPanel inside a new JFrame.
	*/
	public static void main(String[] args) {
		JFrame frame = new JFrame();
		VizGui vg = new VizGui(65250, "fermi.chem.ucl.ac.uk");
		frame.getContentPane().add(vg);
		frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
		frame.addWindowListener( new WindowAdapter() {
		    public void windowClosed(WindowEvent e){
				System.exit(0);
			    }
				public void windowClosing(WindowEvent e) {
				    windowClosed(e);
				}
			    }
		);
		Animator animator = new Animator(vg.canvas1);
		animator.start();

	}
	
	
	/** 
     * Executed exactly once to initialize the 
     * associated GLDrawable
     */ 

    public void init(GLAutoDrawable drawable) {

        // print every openGL call for debugging purposes
	drawable.setGL(new TraceGL(drawable.getGL(), System.err));
	GL gl = drawable.getGL();
	
	gl.glDisable (GL.GL_DEPTH_TEST);
	gl.glDisable (GL.GL_BLEND);
	gl.glShadeModel (GL.GL_FLAT);
	gl.glDisable (GL.GL_DITHER);
	
	/** 
	 * Set the background colour when the GLDrawable 
	 * is cleared 
	 */ 
	gl.glClearColor( 1.0f, 1.0f, 1.0f, 1.0f ); //white

	/** Set the drawing colour to black */
	//gl.glColor3f( 0.0f, 0.0f, 0.0f );
	//drawable.getGL().glPointSize(1.0f); //a 'dot' is 4 by 4 pixels
    }
    
    /** 
     * Executed if the associated GLDrawable is resized
     */ 
    public void reshape(GLAutoDrawable drawable, int x, int y, int width, int height) {
/*	GL gl = drawable.getGL();
	gl.glViewport( 0, 0, width, height );
	gl.glMatrixMode( GL.GL_PROJECTION ); 
	gl.glLoadIdentity();
        gl.glOrtho(0, 400, 0, 300, -1, 1);
    	
    	
   	  float ortho_x = 0.5F * (float)w / (float)PIXELS_X;
      float ortho_y = 0.5F * (float)h / (float)PIXELS_Y;
    	  
    gl.glViewport(0, 0, width, height);
    	  
    	  gl.glLoadIdentity ();
    	  gl.glOrtho2D(-ortho_x, ortho_x, -ortho_y, ortho_y);*/
    }

    /** This method handles the painting of the GLDrawable */
    
    public void display(GLAutoDrawable drawable) {
        VizFrameData vfd = (VizFrameData)queue.poll();
    
        if (queue != null) {	
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
				jTextPane1 = new JTextPane();
				this.add(jTextPane1, "0, 9, 0, 9");
				jTextPane1.setText("Messsages");
				jTextPane1.setEditable(false);
				jTextPane1.setVisible(true)
;				jTextPane1.setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED));
			}
			{
				cap = new GLCapabilities(); 
				canvas1 = new GLCanvas(cap);
				this.add(canvas1, "0, 0, 0, 8");
				canvas1.addGLEventListener(this);
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
    private class NetThread implements Runnable {
    	
    	NetworkReceive nr;
    	
    	public NetThread (int port, String hostname) {
    	    nr = new NetworkReceive (port, hostname);
    	}
    	
    	public void run() {
    		while(true) {
    			queue.offer(nr.getFrame());
    		}
    	}
    }

}
