package uk.ac.ucl.chem.ccs.vizclient;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JSeparator;

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
					{
						jMenuItem3 = new JMenuItem();
						jMenu2.add(jMenuItem3);
						jMenuItem3.setText("Host");
						jMenuItem3.setEnabled(false);
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

}
