package uk.ac.ucl.chem.ccs.clinicalgui;

import java.awt.Dimension;
import java.awt.GridLayout;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.WindowConstants;

import uk.ac.ucl.chem.ccs.clinicalgui.res.ResPanel;


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
public class ReservationsPanel extends javax.swing.JPanel {
	private JPanel resPanel;

	/**
	* Auto-generated main method to display this 
	* JPanel inside a new JFrame.
	*/
	public static void main(String[] args) {
		JFrame frame = new JFrame();
		frame.getContentPane().add(new ReservationsPanel());
		frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
	}
	
	public ReservationsPanel() {
		super();
		initGUI();
	}
	
	private void initGUI() {
		try {
			GridLayout thisLayout = new GridLayout(1, 1);
			thisLayout.setColumns(1);
			thisLayout.setHgap(5);
			thisLayout.setVgap(5);
			this.setLayout(thisLayout);
			setPreferredSize(new Dimension(400, 300));
			{
				resPanel = new ResPanel(ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.resfile"));
				this.add(resPanel);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
