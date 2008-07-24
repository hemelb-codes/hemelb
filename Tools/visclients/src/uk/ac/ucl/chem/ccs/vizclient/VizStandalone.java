package uk.ac.ucl.chem.ccs.vizclient;

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;




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

	/**
	* Auto-generated main method to display this JFrame
	*/
	public static void main(String[] args) {
		String h = "localhost";
		int p = 65250;
		int w = 1024*1024;
		
		switch (args.length) {
		case 3:
			w = Integer.parseInt(args[2]);
		case 2:
			p = Integer.parseInt(args[1]);
		case 1:
			h = args[0];
		}
		
		VizSteererWindow vs = new VizSteererWindow(h,p,w, null);
		vs.setLocationRelativeTo(null);
		vs.setVisible(true);
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
	
}