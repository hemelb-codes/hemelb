package uk.ac.ucl.chem.ccs.vizclient;

import java.awt.*;
import java.awt.event.*;
import javax.media.opengl.*;
import javax.media.opengl.glu.*;

public class firstAttempt extends Frame implements GLEventListener{
	
	VizFrameData vfd;
	
    public firstAttempt() {
    	
    DirectBiConnection nr = new DirectBiConnection (65250, "fermi.chem.ucl.ac.uk");
    vfd = nr.getFrame();
    	
	GLCapabilities caps = new GLCapabilities();
	GLCanvas canvas = new GLCanvas(caps);

	canvas.addGLEventListener(this);
	
	add("Center", canvas);
	setSize(1024,1024);
	setVisible(true);
    }
    
    
    public static void main( String args[]) {
	firstAttempt frame = new firstAttempt();

	
	//exit if frame's close box is clicked
	frame.addWindowListener( new WindowAdapter() {
	    public void windowClosed(WindowEvent e){
		System.exit(0);
	    }
		public void windowClosing(WindowEvent e) {
		    windowClosed(e);
		}
	    }
				 );
	
    }

    /* The functions below are required because we are a GLEventListener.
       We could also have put them in another class and put that class in the
       addGLEventListener method above.
     */

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
	/*GL gl = drawable.getGL();
	gl.glViewport( 0, 0, width, height );
	gl.glMatrixMode( GL.GL_PROJECTION ); 
	gl.glLoadIdentity();
        gl.glOrtho(0, 400, 0, 300, -1, 1);*/
    }

    /** This method handles the painting of the GLDrawable */
    
    public void display(GLAutoDrawable drawable) {
	GL gl = drawable.getGL();
	/** Clear the colour buffer */
	gl.glClear( GL.GL_COLOR_BUFFER_BIT );
	/** Draw some dots */ 
	
	  float scale_x = 1.0f / (float)1024;
	  float scale_y = 1.0f / (float)1024;
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

    /** This method handles things if display depth changes */
    public void displayChanged(GLAutoDrawable drawable,
			       boolean modeChanged,
			       boolean deviceChanged){
    }
    

}

