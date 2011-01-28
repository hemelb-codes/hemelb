package uk.ac.ucl.chem.ccs.vizclient;

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;


import javax.swing.SwingUtilities;

/**
* This code was edited or generated using CloudGarden's Jigloo
* SWT/Swing GUI Builder, which is free for non-commercial
* use. If Jigloo is being used commercially (ie, by a corporation,
* company or business for any purpose whatever) then you
* should purchase a license for each developer using Jigloo.
* Please visit www.cloudgarden.com for details.
* Use of Jigloo implies acceptance of these licensing terms.
* A COMMERCIAL LICENSE HAS NOT BEEN PURCHASED FOR
* THIS MACHINE, SO JIGLOO OR THIS CODE CANNOT BE USED
* LEGALLY FOR ANY CORPORATE OR COMMERCIAL PURPOSE.
*/
public class VizStandalone extends javax.swing.JFrame {

	public VizStandalone (String resourceID, int w, String h, int p) {
		
		final String fresourceID = resourceID;
		final int fw = w;
		final String fh = h;
		final int fp = p;
	
		SwingUtilities.invokeLater(new Runnable() {
		    public void run() {
		VizSteererWindow vs;
		if (fresourceID != null) {
			vs = new VizSteererWindow(fresourceID ,fw, null);

		} else {
			vs = new VizSteererWindow(fh,fp,fw, null);
		}
		vs.setLocationRelativeTo(null);
		vs.setVisible(true);
		vs.showSideBars();
		vs.addWindowListener( new WindowAdapter() {
		    public void windowClosed(WindowEvent e){
					System.exit(0);
			    }
				public void windowClosing(WindowEvent e) {
				    windowClosed(e);
				}
			    }
		);
		    }
		});
	}
 	
	/**
	* Auto-generated main method to display this JFrame
	*/
	public static void main(String[] args) {
		String h = "localhost";
		int p = 65250;
		int w = 1024*1024;
		String resourceID = null;
		switch (args.length) {
		case 3:
			w = Integer.parseInt(args[2]);
		case 2:
			p = Integer.parseInt(args[1]);
		case 1:
			h = args[0];
			if (args[0].equals("-r")) {
				resourceID = args[1];
			}
		}
		
        // Use Apple Aqua L&F screen menu bar if available; set property before any frames created
	       try {
	           java.lang.System.setProperty("apple.laf.useScreenMenuBar", "true");
	       } catch (Exception e) {
	           // try the older menu bar property
	    	   try {
	           java.lang.System.setProperty("com.apple.macos.useScreenMenuBar", "true");
	    	   } catch (Exception e2) {
	    		   //cant set
	    	   }
	       }
		
		VizStandalone vs = new VizStandalone (resourceID, w, h, p);

	}
	
}