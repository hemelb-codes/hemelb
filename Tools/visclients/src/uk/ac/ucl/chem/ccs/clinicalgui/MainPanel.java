package uk.ac.ucl.chem.ccs.clinicalgui;

import java.awt.Dimension;
import java.awt.GridLayout;
import javax.swing.JFrame;
import javax.swing.JTabbedPane;
import javax.swing.WindowConstants;
import javax.swing.event.ChangeListener;
import javax.swing.event.ChangeEvent;
import uk.ac.ucl.chem.ccs.clinicalgui.res.ResPanel;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
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
public class MainPanel extends javax.swing.JPanel {
	private JTabbedPane jTabbedPane1;
	private ResPanel reservationPanel1;
	private ModelGenerationPanel modelGenerationPanel1;
	private SimulationMonitoringPanel simulationMonitoringPanel1;
	private SimulationLaunchPanel simulationLaunchPanel1;
	private int guiStage = 0;

	/**
	* Auto-generated main method to display this 
	* JPanel inside a new JFrame.
	*/
	public static void main(String[] args) {
		JFrame frame = new JFrame();
		frame.getContentPane().add(new MainPanel());
		frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
	}
	
	public MainPanel() {
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
				jTabbedPane1 = new JTabbedPane();
				this.add(jTabbedPane1);
				{
					reservationPanel1 = new ResPanel(ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.resfile"));

					jTabbedPane1.addTab("Reservation", null, reservationPanel1, null);
					reservationPanel1.next.addActionListener(new NextListener());
				}
				{
					modelGenerationPanel1 = new ModelGenerationPanel();
					jTabbedPane1.addTab("Model Generation", null, modelGenerationPanel1, null);
					modelGenerationPanel1.NextButton.addActionListener(new NextListener());
					modelGenerationPanel1.PrevButton.addActionListener(new PrevListener());
				}
				{
					simulationLaunchPanel1 = new SimulationLaunchPanel();
					jTabbedPane1.addTab("Simulation Launch", null, simulationLaunchPanel1, null);
				}
				{
					simulationMonitoringPanel1 = new SimulationMonitoringPanel();
					jTabbedPane1.addTab("Simulation Monitoring", null, simulationMonitoringPanel1, null);
				}
				jTabbedPane1.addChangeListener(new ChangeListener() {
					public void stateChanged(ChangeEvent e) {
						guiStage = jTabbedPane1.getSelectedIndex();
					}
				});
				
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private class NextListener implements ActionListener {
		public void actionPerformed (ActionEvent evt) {
			guiStage++;
			jTabbedPane1.setSelectedIndex(guiStage);
		}
		
	}
	
	private class PrevListener implements ActionListener {
		public void actionPerformed (ActionEvent evt) {
			guiStage--;
			jTabbedPane1.setSelectedIndex(guiStage);
		}
		
	}

}
