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
import java.io.File;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.WindowConstants;
import javax.swing.border.BevelBorder;

import uk.ac.ucl.chem.ccs.vizclient.ExampleFileFilter;

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
	private JRadioButtonMenuItem view2MenuItem;
	private JRadioButtonMenuItem view1MenuItem;
	private JMenu viewMenu;
	private JPanel jPanel1;
	private JPanel jPanel2;
	private JMenuItem quitMenuItem;
	
	private JRadioButtonMenuItem viewpoint1;
	private JRadioButtonMenuItem viewpoint2;
	private JRadioButtonMenuItem viewpoint3;
	private JRadioButtonMenuItem viewpoint4;
	private JRadioButtonMenuItem viewpoint5;

	
	
	private JMenuItem killMenuItem;
	private JMenuItem loadMenuItem;
	private JMenuItem saveMenuItem;
	private JSeparator jSeparator4;
	private JCheckBoxMenuItem window3MenuItem;
	private JCheckBoxMenuItem window2MenuItem;
	private JRadioButtonMenuItem view3MenuItem;
	private JSeparator jSeparator2;
	private JSeparator jSeparator3;
	private JMenuItem hostMenuItem;
	private JSeparator jSeparator1;
	private JMenuItem disconnectMenuItem;
	private VizGui vg;
	private String hostname;
	private int port, window;
	private Component parent;
	private boolean rendezvousing = false;
	private JCheckBoxMenuItem rotateMenuItem;
	private JCheckBoxMenuItem window1MenuItem;
	private ButtonGroup buttonGroup1;
	private ButtonGroup buttonGroup2;

	/**
	* Auto-generated main method to display this JFrame
	*/
	public static void main(String[] args) {
		String h = "localhost";
		int p = 65250;
		int w = 1024*1024;
		
		VizSteererWindow vs = new VizSteererWindow(h,p,w, null);



	}

	public VizSteererWindow(String resourceID, int window, Component parent) {
		super("HemeLB Steering Client");
		this.window = window;
		initGUI();
		//this.parent = parent;
		this.setLocationRelativeTo(parent);
		//this.setVisible(true);
		rendezvousing = true;
		//poll for ID
		String service = "http://bunsen.chem.ucl.ac.uk:28080/ahe/test/rendezvous";
		HttpGetPoll hgp = new HttpGetPoll(service, resourceID);

		vg.appendNotification("Polling rendezvous service for connection data\n");
		
		if (hgp.pollService(5)) {
		port=hgp.getPort();
		hostname= hgp.getHost();
		vg.appendNotification("Found sim listening at " + hostname + ":" + port + "\n");
		vg.setHostPort(port, hostname);
		vg.startReceive();		
		} else {
			vg.appendNotification("Couldn't rendezvous with HemeLB simulation");
		}
		rendezvousing = false;
	}
	
	public VizSteererWindow(String hostname, int port, int window, Component parent) {
		super("HemeLB Steering Client");
		this.window = window;
		this.hostname = hostname;
		this.port = port;
		this.parent = parent;
		initGUI();
		this.setVisible(true);
		this.setLocationRelativeTo(parent);

	}
	
	public VizSteererWindow(String hostname, int port, int window, Component parent, boolean connectAtOnce) {
	
		this(hostname, port, window, parent);
		if (connectAtOnce) {
			vg.startReceive();
		}
	}
	
	public void showSideBars() {
		vg.showParamWindow1(true);
		vg.showParamWindow2(true);
		vg.showParamWindow3(true);
	}
	
	private void initGUI() {
		try {
			setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
			{
				jPanel1 = new JPanel();
				
				TableLayout tabPanelLayout = new TableLayout(new double[][] {{TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}, {TableLayout.FILL}});
				tabPanelLayout.setHGap(5);
				tabPanelLayout.setVGap(5);
				
				//GridBagLayout jPanel1Layout = new GridBagLayout();
				//jPanel1Layout.rowWeights = new double[] {0.1, 0.1};
				//jPanel1Layout.rowHeights = new int[] {7, 7};
				//jPanel1Layout.columnWeights = new double[] {0.1, 0.1};
				//jPanel1Layout.columnWidths = new int[] {7, 7};
				
				
				jPanel1.setLayout(tabPanelLayout);
				getContentPane().add(jPanel1, BorderLayout.CENTER);
				//jPanel1.setSize(700, 700);
				{
					vg=new VizGui(port, hostname, window, this);
//					jPanel1.add(vg, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));

					jPanel1.add(vg, "0, 0, 2, 0");

					vg.setBorder(BorderFactory.createEtchedBorder(BevelBorder.LOWERED));
				//	vg.getInfoPanel().setBorder(BorderFactory.createEtchedBorder(BevelBorder.LOWERED));
				}
				{
					jPanel2 = vg.getInfoPanel();
//					jPanel1.add(jPanel2, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
					jPanel1.add(jPanel2, "3, 0, 3, 0");

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
							if (vg.isConnected() && rendezvousing == false){
								connectMenuItem.setEnabled(false);
								disconnectMenuItem.setEnabled(true);
								hostMenuItem.setEnabled(true);
							} else if (!vg.isConnected() && rendezvousing == false) {
								connectMenuItem.setEnabled(true);
								disconnectMenuItem.setEnabled(false);
								hostMenuItem.setEnabled(true);
							} else {
								connectMenuItem.setEnabled(false);
								disconnectMenuItem.setEnabled(false);
								hostMenuItem.setEnabled(false);

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
						hostMenuItem.setText("Set HemeLB host");
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

						if (parent == null) {
							quitMenuItem.setText("Quit");
						} else {
							quitMenuItem.setText("Close");
						}
						quitMenuItem.setEnabled(true);
						quitMenuItem.addActionListener(new ActionListener (){ 
							public void actionPerformed (ActionEvent e) {
								vg.stopReceive();
								VizSteererWindow.this.dispose();
							}
						});
					}
					{
						steeringMenu = new JMenu();
						mainMenuBar.add(steeringMenu);
						steeringMenu.setText("Steering");
					}
					{
						viewMenu = new JMenu();
						mainMenuBar.add(viewMenu);
						viewMenu.setText("View");
						
					//TODO Add view menu listener, so right option is chosen.
						{
							view1MenuItem = new JRadioButtonMenuItem();
							viewMenu.add(view1MenuItem);
							view1MenuItem.setText("Vessel wall");
							view1MenuItem.setSelected(true);
							view1MenuItem.addActionListener(new ActionListener (){ 
								public void actionPerformed (ActionEvent e) {
									vg.changeVisMode(0);
								}
							});
						}
						{
							view2MenuItem = new JRadioButtonMenuItem();
							viewMenu.add(view2MenuItem);
							view2MenuItem.setText("Velocity glyphs");
							view2MenuItem.addActionListener(new ActionListener (){ 
								public void actionPerformed (ActionEvent e) {
									vg.changeVisMode(1);
								}
							});
						}
						{
							view3MenuItem = new JRadioButtonMenuItem();
							viewMenu.add(view3MenuItem);
							view3MenuItem.setText("Streaklines");
							view3MenuItem.addActionListener(new ActionListener (){ 
								public void actionPerformed (ActionEvent e) {
									vg.changeVisMode(2);
								}
							});
						}

						viewMenu.add(new JSeparator());
						
						viewpoint1 = new JRadioButtonMenuItem();
						viewMenu.add(viewpoint1);
						viewpoint1.setText("View 1");
						viewpoint1.addActionListener(new ActionListener (){ 
							public void actionPerformed (ActionEvent e) {
								vg.viewChanged(VizGui.VIEW1);
							}
						});
						
						viewpoint2 = new JRadioButtonMenuItem();
						viewMenu.add(viewpoint2);
						viewpoint2.setText("View 2");
						viewpoint2.addActionListener(new ActionListener (){ 
							public void actionPerformed (ActionEvent e) {
								vg.viewChanged(VizGui.VIEW2);
							}
						});
						
						viewpoint3 = new JRadioButtonMenuItem();
						viewMenu.add(viewpoint3);
						viewpoint3.setText("View 3");
						viewpoint3.addActionListener(new ActionListener (){ 
							public void actionPerformed (ActionEvent e) {
								vg.viewChanged(VizGui.VIEW3);
							}
						});
						
						viewpoint4 = new JRadioButtonMenuItem();
						viewMenu.add(viewpoint4);
						viewpoint4.setText("View 4");
						viewpoint4.addActionListener(new ActionListener (){ 
							public void actionPerformed (ActionEvent e) {
								vg.viewChanged(VizGui.VIEW4);
							}
						});

						viewpoint5 = new JRadioButtonMenuItem();
						viewMenu.add(viewpoint5);
						viewpoint5.setText("View 5");
						viewpoint5.addActionListener(new ActionListener (){ 
							public void actionPerformed (ActionEvent e) {
								vg.viewChanged(VizGui.VIEWALL);
							}
						});
						
						{
							buttonGroup1 = new ButtonGroup();
							buttonGroup1.add(view1MenuItem);
							buttonGroup1.add(view3MenuItem);
							buttonGroup1.add(view2MenuItem);
						}
						
						{
							buttonGroup2 = new ButtonGroup();
							buttonGroup2.add(viewpoint1);
							buttonGroup2.add(viewpoint2);
							buttonGroup2.add(viewpoint3);
							buttonGroup2.add(viewpoint4);
							buttonGroup2.add(viewpoint5);
						}
						
						{
							viewMenu.addMenuListener(new MenuListener () {
								public void menuSelected (MenuEvent e) {
										switch (vg.getVisMode()) {
										case 0:
											view1MenuItem.setSelected(true);
											break;
										case 1:
											view2MenuItem.setSelected(true);
											break;
										case 2:
											view3MenuItem.setSelected(true);
											break;
										}
										
										switch (vg.getView()) {
										case VizGui.VIEW1:
											viewpoint1.setSelected(true);
											break;
										case VizGui.VIEW2:
											viewpoint2.setSelected(true);
											break;
										case VizGui.VIEW3:
											viewpoint3.setSelected(true);
											break;
										case VizGui.VIEW4:
											viewpoint4.setSelected(true);
											break;
										case VizGui.VIEWALL:
											viewpoint5.setSelected(true);
											break;
										}
								}
			

							public void menuCanceled (MenuEvent e) {
								
							}
							
							public void menuDeselected (MenuEvent e) {
								
							}
							
						});
						}
						
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
					{	
					killMenuItem = new JMenuItem();
					steeringMenu.add(killMenuItem);
					killMenuItem.setText("Kill HemeLB");
					killMenuItem.addActionListener(new ActionListener() {
						public void actionPerformed (ActionEvent e) {
							vg.kill();
						}
					});
					
				} 
					{
						rotateMenuItem = new JCheckBoxMenuItem();
						rotateMenuItem.setText("Auto-rotate");
						steeringMenu.add(rotateMenuItem);
						rotateMenuItem.addActionListener(new ActionListener() {
							public void actionPerformed (ActionEvent e) {
	
									vg.rotate(rotateMenuItem.getState() );

							}
						});
					}
					{
						jSeparator3 = new JSeparator();
						steeringMenu.add(jSeparator3);
					}
					{
						window1MenuItem = new JCheckBoxMenuItem();
						window1MenuItem.setText("Scene parameters");
						steeringMenu.add(window1MenuItem);
						window1MenuItem.addActionListener(new ActionListener() {
							public void actionPerformed (ActionEvent e) {
									vg.showParamWindow1(window1MenuItem.getState());

							}
						});
					}
					{
						window2MenuItem = new JCheckBoxMenuItem();
						steeringMenu.add(window2MenuItem);
						window2MenuItem.setText("Resolution");
						window2MenuItem.addActionListener(new ActionListener() {
							public void actionPerformed (ActionEvent e) {
									vg.showParamWindow2(window2MenuItem.getState());

							}
						});
					}
					{
						window3MenuItem = new JCheckBoxMenuItem();
						steeringMenu.add(window3MenuItem);
						window3MenuItem.setText("Render parameters");
						window3MenuItem.addActionListener(new ActionListener() {
							public void actionPerformed (ActionEvent e) {
									vg.showParamWindow3(window3MenuItem.getState());

							}
						});
					}
					{
						jSeparator4 = new JSeparator();
						steeringMenu.add(jSeparator4);
					}
					{
						saveMenuItem = new JMenuItem();
						steeringMenu.add(saveMenuItem);
						saveMenuItem.setText("Save Parameters");
						saveMenuItem.addActionListener(new ActionListener() {
							public void actionPerformed (ActionEvent evt) {
								JFileChooser fc = new JFileChooser();	
								
							    ExampleFileFilter filter = new ExampleFileFilter();
							    filter.addExtension("hcd");
							    filter.setDescription("HemeLB Client Data");
							    fc.setFileFilter(filter);

								fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES );
								int returnVal = fc.showSaveDialog(VizSteererWindow.this);

							    if (returnVal == JFileChooser.APPROVE_OPTION) {
							    	String filename = fc.getSelectedFile().getAbsolutePath();
							    	if (!filename.endsWith(".hcd")) {
							    		filename = filename + ".hcd";
							    	}
							        vg.saveParameters(filename);
							        //System.out.println(fileLocation);
							    } 
							}
						});
					}
					{
						loadMenuItem = new JMenuItem();
						steeringMenu.add(loadMenuItem);
						loadMenuItem.setText("Load Paramters");
						loadMenuItem.addActionListener(new ActionListener() {
							public void actionPerformed (ActionEvent evt) {
								JFileChooser fc = new JFileChooser();	
								
							    ExampleFileFilter filter = new ExampleFileFilter();
							    filter.addExtension("hcd");
							    filter.setDescription("HemeLB Client Data");
							    fc.setFileFilter(filter);

								fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES );
								int returnVal = fc.showOpenDialog(VizSteererWindow.this);

							    if (returnVal == JFileChooser.APPROVE_OPTION) {
							        vg.loadParameters(fc.getSelectedFile().getAbsolutePath());
							        //System.out.println(fileLocation);
							    } 
							}
						});
					}
					{
						steeringMenu.addMenuListener(new MenuListener () {
							public void menuSelected (MenuEvent e) {
								window1MenuItem.setState(vg.paramWindow1Visible());
								window2MenuItem.setState(vg.paramWindow2Visible());
								window3MenuItem.setState(vg.paramWindow3Visible());

							}
						
						
						public void menuCanceled (MenuEvent e) {
							
						}
						
						public void menuDeselected (MenuEvent e) {
							
						}
						
					});
					}
					
				}
			}

			pack();
			this.setSize(1000, 800);
			
			//setSize(400, 300);
			//Application.getInstance().getContext().getResourceMap(getClass()).injectComponents(getContentPane());
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
