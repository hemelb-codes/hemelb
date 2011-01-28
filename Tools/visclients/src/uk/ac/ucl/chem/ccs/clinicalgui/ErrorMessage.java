/*
 * AHE: Application Hosting Environment
 *
 * (C) Copyright 2006, University College London, United Kingdom
 * (C) Copyright 2006, University of Manchester, United Kingdom
 *
 * The Application Hosting Environment(AHE) comes with no warranty of
 * any kind. It is a copyrighted code distributed free of charge under
 * the terms of the GNU Public License (http://www.gnu.org/copyleft/gpl.html),
 * which is commonly known as "open source" distribution. This means that
 * anyone is free to use, modify, or extend AHE in any way they choose, but
 * if you distribute a modified version of AHE, it must remain open-source,
 * meaning you distribute it under the terms of the GPL. You should clearly
 * annotate such a code as a derivative version of AHE. If you release any code 
 * that includes AHE source code, then it must also be open-sourced, meaning 
 * you distribute it under the terms of the GPL.
 *
 */
/*
 * Project: AHE-GUI 
 *
 * @author stefan.zasada@ucl.ac.uk
 * 
 */

package uk.ac.ucl.chem.ccs.clinicalgui;
import java.awt.Container;
import java.awt.Dialog;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.SwingConstants;

public class ErrorMessage extends javax.swing.JDialog implements ActionListener {
	private JButton okButton;
	private JLabel message;


	
	public ErrorMessage(Container cont, String errormsg) {
		super((Frame)cont, true);
		initGUI(errormsg);
		this.setLocationRelativeTo((Frame)cont);
		this.setVisible(true);
	}
	
	public ErrorMessage(Frame frame, String errormsg) {
		super(frame, true);
		initGUI(errormsg);
		this.setLocationRelativeTo(frame);
		this.setVisible(true);
	}
	
	public ErrorMessage(Dialog frame, String errormsg) {
		super(frame, true);
		initGUI(errormsg);
		this.setLocationRelativeTo(frame);
		this.setVisible(true);
	}

	private void initGUI(String errormsg) {
		try {
			{
				GridBagLayout thisLayout = new GridBagLayout();
				thisLayout.rowWeights = new double[] {0.9, 0.1};
				thisLayout.rowHeights = new int[] {7, 3};
				thisLayout.columnWeights = new double[] {0.1};
				thisLayout.columnWidths = new int[] {7};
				getContentPane().setLayout(thisLayout);
				this.setTitle("AHE Client Error");
		        this.setResizable(false);
			}
			{
				okButton = new JButton();
				getContentPane().add(okButton, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.VERTICAL, new Insets(0, 0, 0, 0), 0, 0));
				okButton.setText("OK");
				okButton.setBounds(29, 140, 50, 22);
				okButton.addActionListener(this);
				okButton.setFont(new java.awt.Font("Sansserif",0,11));
				okButton.setPreferredSize(new java.awt.Dimension(76, 24));
				this.getRootPane().setDefaultButton(okButton);
			}
			{
				message = new JLabel();
				message.setText("<HTML>" + errormsg + "</HTML>");
				getContentPane().add(message, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
				message.setBounds(0, 7, 399, 105);
				message.setFont(new java.awt.Font("Sansserif",0,12));
				message.setHorizontalAlignment(SwingConstants.CENTER);
				message.setHorizontalTextPosition(SwingConstants.CENTER);
			}
			this.setSize(333, 187);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void actionPerformed(ActionEvent e) {

        this.setVisible(false);
        this.dispose();


	}

}
