package uk.ac.ucl.chem.ccs.vizclient;

import info.clearthought.layout.TableLayout;
import java.awt.Canvas;

import java.awt.Dimension;
import javax.swing.BorderFactory;

import javax.swing.WindowConstants;
import javax.swing.border.BevelBorder;
import javax.swing.JFrame;
import javax.swing.JInternalFrame;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.JSplitPane;
import javax.swing.JTextPane;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;

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
public class VizGui extends javax.swing.JPanel {
	private JTextPane jTextPane1;
	private GLCanvas canvas1;
	private GLCapabilities cap;
	
	/**
	* Auto-generated main method to display this 
	* JPanel inside a new JFrame.
	*/
	public static void main(String[] args) {
		JFrame frame = new JFrame();
		frame.getContentPane().add(new VizGui());
		frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
	}
	
	public VizGui() {
		super();
		initGUI();
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
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
