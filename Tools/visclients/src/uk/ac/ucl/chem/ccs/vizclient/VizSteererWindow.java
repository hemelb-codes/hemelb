package uk.ac.ucl.chem.ccs.vizclient;
import info.clearthought.layout.TableLayout;
import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.event.MenuListener;
import javax.swing.event.MenuEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

import javax.swing.WindowConstants;
import javax.swing.border.BevelBorder;



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
public class VizSteererWindow extends javax.swing.JFrame {
	private JMenuBar mainMenuBar;
	private JMenu connectMenu;
	private JMenu steeringMenu;
	private JMenuItem resetMenuItem;
	private JMenuItem connectMenuItem;
	private JPanel jPanel1;
	private JPanel jPanel2;
	private JMenuItem quitMenuItem;
	private JSeparator jSeparator2;
	private JMenuItem hostMenuItem;
	private JSeparator jSeparator1;
	private JMenuItem disconnectMenuItem;
	private VizGui vg;
	private String hostname;
	private int port, window;
	/**
	* Auto-generated main method to display this JFrame
	*/
	public static void main(String[] args) {
		String h = "localhost";
		int p = 65250;
		int w = 1024*1024;
		
		VizSteererWindow vs = new VizSteererWindow(h,p,w, null);



	}
	
	public VizSteererWindow(String hostname, int port, int window, Component parent) {
		super();
		this.window = window;
		this.hostname = hostname;
		this.port = port;
		initGUI();
		this.setLocationRelativeTo(parent);
		this.setVisible(true);
	}
	
	private void initGUI() {
		try {
			setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
			{
				jPanel1 = new JPanel();
				GridBagLayout jPanel1Layout = new GridBagLayout();
				jPanel1Layout.rowWeights = new double[] {0.1, 0.1};
				jPanel1Layout.rowHeights = new int[] {7, 7};
				jPanel1Layout.columnWeights = new double[] {0.1, 0.1};
				jPanel1Layout.columnWidths = new int[] {7, 7};
				jPanel1.setLayout(jPanel1Layout);
				getContentPane().add(jPanel1, BorderLayout.CENTER);
				jPanel1.setSize(700, 700);
				{
					vg=new VizGui(port, hostname, window);
					jPanel1.add(vg, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
					vg.setBorder(BorderFactory.createEtchedBorder(BevelBorder.LOWERED));
					vg.getInfoPanel().setBorder(BorderFactory.createEtchedBorder(BevelBorder.LOWERED));
				}
				{
					jPanel2 = vg.getInfoPanel();
					jPanel1.add(jPanel2, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
				}

			}
			{
				mainMenuBar = new JMenuBar();
				setJMenuBar(mainMenuBar);
				{
					connectMenu = new JMenu();
					mainMenuBar.add(connectMenu);
					connectMenu.setText("Connection");
					connectMenu.addMenuListener(new MenuListener () {
						public void menuSelected (MenuEvent e) {
							if (vg.isConnected() ){
								connectMenuItem.setEnabled(false);
								disconnectMenuItem.setEnabled(true);
							} else {
								connectMenuItem.setEnabled(true);
								disconnectMenuItem.setEnabled(false);
							}
						}
						
						public void menuCanceled (MenuEvent e) {
							
						}
						
						public void menuDeselected (MenuEvent e) {
							
						}
						
					});
					{
						hostMenuItem = new JMenuItem();
						connectMenu.add(hostMenuItem);
						hostMenuItem.setText("Host");
						hostMenuItem.setEnabled(true);
						hostMenuItem.addActionListener(new ActionListener (){ 
							public void actionPerformed (ActionEvent e) {
								HostDialogue hostDia = new HostDialogue(VizSteererWindow.this);								
								if (hostDia.showDialogue()) {
									if (disconnectMenuItem.isEnabled()) {
										if (vg.stopReceive()) {
											connectMenuItem.setEnabled(true);
											disconnectMenuItem.setEnabled(false);
											}
									}
									VizSteererWindow.this.vg.setHostPort(port, hostname);
									//System.err.println("Using hostname " + hostname + ":" + port);
								}
							}
						});
					}
					{
						jSeparator1 = new JSeparator();
						connectMenu.add(jSeparator1);
					}
					{
						connectMenuItem = new JMenuItem();
						connectMenu.add(connectMenuItem);
						connectMenuItem.setText("Connect");
						connectMenuItem.addActionListener(new ActionListener (){ 
							public void actionPerformed (ActionEvent e) {
								if (vg.startReceive()) {
								disconnectMenuItem.setEnabled(true);
								connectMenuItem.setEnabled(false);
								}
							}
						});
					}
					{
						disconnectMenuItem = new JMenuItem();
						connectMenu.add(disconnectMenuItem);
						disconnectMenuItem.setText("Disconnect");
						disconnectMenuItem.setEnabled(false);
						disconnectMenuItem.addActionListener(new ActionListener (){ 
							public void actionPerformed (ActionEvent e) {
								if (vg.stopReceive()) {
								connectMenuItem.setEnabled(true);
								disconnectMenuItem.setEnabled(false);
								}
							}
						});
					}
					{
						jSeparator2 = new JSeparator();
						connectMenu.add(jSeparator2);
					}
					{
						quitMenuItem = new JMenuItem();
						connectMenu.add(quitMenuItem);
						quitMenuItem.setText("Quit");
						quitMenuItem.setEnabled(true);
						quitMenuItem.addActionListener(new ActionListener (){ 
							public void actionPerformed (ActionEvent e) {
								VizSteererWindow.this.dispose();
							}
						});
					}
					{
						steeringMenu = new JMenu();
						mainMenuBar.add(steeringMenu);
						steeringMenu.setText("Steering");
					}
					{	resetMenuItem = new JMenuItem();
						steeringMenu.add(resetMenuItem);
						resetMenuItem.setText("Reset");
						resetMenuItem.addActionListener(new ActionListener() {
							public void actionPerformed (ActionEvent e) {
								vg.resetSteering();
							}
						});
						
					}
					
				}
			}

			pack();
			this.setSize(1000, 800);
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
