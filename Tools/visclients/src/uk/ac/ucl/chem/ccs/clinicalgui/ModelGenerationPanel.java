package uk.ac.ucl.chem.ccs.clinicalgui;
import ij.plugin.DICOM;
import info.clearthought.layout.TableLayout;

import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.WindowConstants;
import org.jdesktop.application.Application;
import java.io.*;

import org.dcm4che2.*;
import org.dcm4che.data.*;


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
public class ModelGenerationPanel extends javax.swing.JPanel {
	private JLabel ImageLabel;
	private JPanel jPanel1;
	private JButton LaunchGraphClientButton;
	public JButton PrevButton;
	public JButton NextButton;

	/**
	* Auto-generated main method to display this 
	* JPanel inside a new JFrame.
	*/
	public static void main(String[] args) {
		JFrame frame = new JFrame();
		frame.getContentPane().add(new ModelGenerationPanel());
		frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
	}
	
	public ModelGenerationPanel() {
		super();
		initGUI();
	}
	
	private void initGUI() {
		try {
			TableLayout thisLayout = new TableLayout(new double[][] {{TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}, {TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}});
			thisLayout.setHGap(5);
			thisLayout.setVGap(5);
			this.setLayout(thisLayout);
			setPreferredSize(new Dimension(400, 300));
			{
				//make image from DICOM file and add it to a label (to be displayed in GUI)
				DICOM d = new DICOM(new FileInputStream("C:/users/konstantin/desktop/dicomFiles/file6.dcm"));
				d.run("Name");
		        //d.show();
				ImageLabel = new JLabel(new ImageIcon(d.getImage()));
				this.add(ImageLabel, "2, 0, 3, 2");
				ImageLabel.setName("ImageLabel");
			}
			{
				jPanel1 = new JPanel();
				TableLayout jPanel1Layout = new TableLayout(new double[][] {{TableLayout.FILL}, {TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}});
				jPanel1Layout.setHGap(5);
				jPanel1Layout.setVGap(5);
				this.add(jPanel1, "0, 0, 1, 2");
				jPanel1.setLayout(jPanel1Layout);
			}
			{
				NextButton = new JButton();
				this.add(NextButton, "3, 3");
				NextButton.setName("NextButton");
			}
			{
				PrevButton = new JButton();
				this.add(PrevButton, "0, 3");
				PrevButton.setName("PrevButton");
			}
			{
				LaunchGraphClientButton = new JButton();
				this.add(LaunchGraphClientButton, "1, 3, 2, 3");
				LaunchGraphClientButton.setName("LaunchGraphClientButton");
			}
			Application.getInstance().getContext().getResourceMap(getClass()).injectComponents(this);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
