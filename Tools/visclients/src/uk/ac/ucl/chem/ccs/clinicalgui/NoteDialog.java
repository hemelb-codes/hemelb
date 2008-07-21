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
public class NoteDialog extends javax.swing.JDialog {
	private JLabel PromptLabel;
	public JButton SubmitButton;
	public JTextField NoteTextField;

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
	
	public NoteDialog(JFrame frame) {
		super(frame);
		initGUI();
	}
	
	private void initGUI() {
		try {
			TableLayout thisLayout = new TableLayout(new double[][] {{TableLayout.FILL}, {TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}});
			thisLayout.setHGap(5);
			thisLayout.setVGap(5);
			getContentPane().setLayout(thisLayout);
			    //add label displaying current patient ID
				PromptLabel = new JLabel();
				getContentPane().add(PromptLabel, "0,0");
				PromptLabel.setText("Enter a short note about this model below");
				getContentPane().add(getSubmitButton(),"0, 2");
				getContentPane().add(getNoteTextField(), "0, 1");
				SubmitButton.addActionListener(new ActionListener(){
					public void actionPerformed (ActionEvent evt) {
						ModelGenerationPanel.modelNote = NoteTextField.getText();
						NoteDialog.this.dispose();		
					}
			});
				
			this.setSize(402, 216);
			Application.getInstance().getContext().getResourceMap(getClass()).injectComponents(getContentPane());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private JTextField getNoteTextField() {
		if(NoteTextField == null) {
			NoteTextField = new JTextField();
			NoteTextField.setName("NoteTextField");
		}
		return NoteTextField;
	}
	
	private JButton getSubmitButton() {
		if(SubmitButton == null) {
			SubmitButton = new JButton();
			SubmitButton.setName("SubmitButton");
		}
		return SubmitButton;
	}
	

}
