package uk.ac.ucl.chem.ccs.vizclient;
import info.clearthought.layout.TableLayout;
import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

import javax.swing.WindowConstants;
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
	private JMenuBar jMenuBar1;
	private JMenu jMenu2;
	private JMenuItem jMenuItem1;

	private JMenuItem jMenuItem3;
	private JSeparator jSeparator1;
	private JMenuItem jMenuItem2;
	private VizGui vg;
	private String hostname;
	private int port;
	/**
	* Auto-generated main method to display this JFrame
	*/
	public static void main(String[] args) {
		
		VizStandalone vs = new VizStandalone(args[0],Integer.parseInt(args[1]));
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
	
	public VizStandalone(String hostname, int port) {
		super();
		this.hostname = hostname;
		this.port = port;
		initGUI();
	}
	
	private void initGUI() {
		try {
			setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
			{
				vg=new VizGui(port, hostname);
				getContentPane().add(vg);
			}
			{
				jMenuBar1 = new JMenuBar();
				setJMenuBar(jMenuBar1);
				{
					jMenu2 = new JMenu();
					jMenuBar1.add(jMenu2);
					jMenu2.setText("Connection");
					jMenu2.addActionListener(new ActionListener () {
						public void actionPerformed (ActionEvent e) {
							if (vg.isConnected() ){
								jMenuItem1.setEnabled(false);
								jMenuItem2.setEnabled(true);
							} else {
								jMenuItem1.setEnabled(true);
								jMenuItem2.setEnabled(false);
							}
						}
					});
					{
						jMenuItem3 = new JMenuItem();
						jMenu2.add(jMenuItem3);
						jMenuItem3.setText("Host");
						jMenuItem3.setEnabled(true);
						jMenuItem3.addActionListener(new ActionListener (){ 
							public void actionPerformed (ActionEvent e) {
								HostDialogue hostDia = new HostDialogue(VizStandalone.this);								
								if (hostDia.showDialogue()) {
									if (jMenuItem2.isEnabled()) {
										if (vg.stopReceive()) {
											jMenuItem1.setEnabled(true);
											jMenuItem2.setEnabled(false);
											}
									}
									VizStandalone.this.vg.setHostPort(port, hostname);
									//System.err.println("Using hostname " + hostname + ":" + port);
								}
							}
						});
					}
					{
						jSeparator1 = new JSeparator();
						jMenu2.add(jSeparator1);
					}
					{
						jMenuItem1 = new JMenuItem();
						jMenu2.add(jMenuItem1);
						jMenuItem1.setText("Connect");
						jMenuItem1.addActionListener(new ActionListener (){ 
							public void actionPerformed (ActionEvent e) {
								if (vg.startReceive()) {
								jMenuItem2.setEnabled(true);
								jMenuItem1.setEnabled(false);
								}
							}
						});
					}
					{
						jMenuItem2 = new JMenuItem();
						jMenu2.add(jMenuItem2);
						jMenuItem2.setText("Disconnect");
						jMenuItem2.setEnabled(false);
						jMenuItem2.addActionListener(new ActionListener (){ 
							public void actionPerformed (ActionEvent e) {
								if (vg.stopReceive()) {
								jMenuItem1.setEnabled(true);
								jMenuItem2.setEnabled(false);
								}
							}
						});
					}
				}
			}

			pack();
			//setSize(400, 300);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private class HostDialogue extends JDialog {
		private boolean changed = false;
		
		private JButton okButton;
		private JButton cancelButton;
		private JTextField portField;
		private JTextField hostField;
		private JLabel jLabel2;
		private JLabel jLabel3;
		
		public HostDialogue (JFrame f) {
			super(f);
			init();
			this.setLocationRelativeTo(f);
		}
		
		public boolean showDialogue() {
			this.setVisible (true);
			return changed;
		}
		
		private void init() {
		
			
			TableLayout jDialog1Layout = new TableLayout(new double[][] {{TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}, {TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}});
			jDialog1Layout.setHGap(5);
			jDialog1Layout.setVGap(5);
			this.getContentPane().setLayout(jDialog1Layout);
			this.setTitle("Change host details");
			this.setModal(true);
			{
				jLabel2 = new JLabel();
				this.getContentPane().add(jLabel2, "1, 1, 2, 2");
				jLabel2.setText("Host");
				jLabel2.setHorizontalAlignment(SwingConstants.RIGHT);
				jLabel2.setFont(new java.awt.Font("Dialog",1,12));
			}
			{
				jLabel3 = new JLabel();
				this.getContentPane().add(jLabel3, "1, 3, 2, 4");
				jLabel3.setText("Port");
				jLabel3.setHorizontalAlignment(SwingConstants.RIGHT);
				jLabel3.setFont(new java.awt.Font("Dialog",1,12));
			}
			{
				hostField = new JTextField();
				hostField.setText(hostname);
				this.getContentPane().add(hostField, "3, 1, 10, 2");
				
			}
			{
				portField = new JTextField();
				portField.setText(Integer.toString(port));
				this.getContentPane().add(portField, "3, 3, 4, 4");
			}
			{
				cancelButton = new JButton();
				this.getContentPane().add(cancelButton, "5, 3, 7, 4");
				cancelButton.setText("Cancel");
				cancelButton.addActionListener(new ActionListener () {
					public void actionPerformed (ActionEvent e) {
						HostDialogue.this.setVisible(false);
						HostDialogue.this.dispose();
					}
				}
				);
				
			}
			{
				okButton = new JButton();
				this.getContentPane().add(okButton, "8, 3, 10, 4");
				okButton.setText("OK");
				okButton.addActionListener(new ActionListener () {
					public void actionPerformed (ActionEvent e) {
						String newHostname = hostField.getText();
						int newPort =  Integer.parseInt(portField.getText());
						
						if (port != newPort || !hostname.equals(newHostname)) {
						hostname = newHostname;
						port = newPort;
						changed=true;
						}
						HostDialogue.this.setVisible(false);
						HostDialogue.this.dispose();
					}
				}
				);
			}
			this.setSize(386, 95);
		}
	

	}}
