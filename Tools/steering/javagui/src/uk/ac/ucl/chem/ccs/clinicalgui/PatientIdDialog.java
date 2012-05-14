package uk.ac.ucl.chem.ccs.clinicalgui;
import info.clearthought.layout.TableLayout;
import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import org.jdesktop.application.Application;

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
public class PatientIdDialog extends javax.swing.JDialog {
	private JLabel PatientLabel;
	private JButton SubmitButton;
	private JButton CancelButton;
	private JTextField PatientTextField;
	private JLabel ChangePatientPrompt;

	/**
	* Auto-generated main method to display this JDialog
	*/
	public static void main(String[] args) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				JFrame frame = new JFrame();
				PatientIdDialog inst = new PatientIdDialog(frame);
				inst.setVisible(true);
			}
		});
	}
	
	public PatientIdDialog(JFrame frame) {
		super(frame);
		initGUI();
	}
	
	private void initGUI() {
		try {
			TableLayout thisLayout = new TableLayout(new double[][] {{TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}, {TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}});
			thisLayout.setHGap(5);
			thisLayout.setVGap(5);
			getContentPane().setLayout(thisLayout);
			    //add label displaying current patient ID
				PatientLabel = new JLabel();
				PatientLabel.setLayout(null);
				getContentPane().add(PatientLabel, "1, 2, 2, 2");
				PatientLabel.setText("selected patient: " + ClinicalGuiClient.PATIENT_ID);
				PatientLabel.setPreferredSize(new java.awt.Dimension(150, 150));
				//add label prompting user to enter new patient ID
				ChangePatientPrompt = new JLabel();
				getContentPane().add(ChangePatientPrompt, "0, 4");
				ChangePatientPrompt.setName("ChangePatientPrompt");
				//add submit and cancel buttons and text field for new patient ID
				getContentPane().add(getSubmitButton(),"2, 4");
				getContentPane().add(getCancelButton(), "3, 4");
				getContentPane().add(getPatientTextField(), "1, 4");
				SubmitButton.addActionListener(new ActionListener() {
					public void actionPerformed (ActionEvent evt) {
						ClinicalGuiClient.PATIENT_ID = PatientTextField.getText();
						PatientLabel.setText("selected patient: " + ClinicalGuiClient.PATIENT_ID);
						//PatientIdDialog.this.repaint();
					}
				});
				CancelButton.addActionListener(new ActionListener() {
					public void actionPerformed (ActionEvent evt) {
						PatientIdDialog.this.dispose();
					}
				});
			this.setSize(402, 216);
			Application.getInstance().getContext().getResourceMap(getClass()).injectComponents(getContentPane());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private JTextField getPatientTextField() {
		if(PatientTextField == null) {
			PatientTextField = new JTextField();
			PatientTextField.setName("PatientTextField");
		}
		return PatientTextField;
	}
	
	private JButton getSubmitButton() {
		if(SubmitButton == null) {
			SubmitButton = new JButton();
			SubmitButton.setName("SubmitButton");
		}
		return SubmitButton;
	}
	
	private JButton getCancelButton() {
		if(CancelButton == null) {
			CancelButton = new JButton();
			CancelButton.setName("CancelButton");
		}
		return CancelButton;
	}

}
