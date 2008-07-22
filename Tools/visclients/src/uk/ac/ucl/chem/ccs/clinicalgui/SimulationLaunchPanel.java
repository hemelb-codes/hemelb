package uk.ac.ucl.chem.ccs.clinicalgui;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.DICOM;
import info.clearthought.layout.TableLayout;

import java.awt.Dimension;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;

import javax.swing.ImageIcon;
import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.WindowConstants;
import javax.swing.JTree;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;

import org.globus.ftp.GridFTPSession;
import org.jdesktop.application.Application;
import java.io.*;
import java.util.Enumeration;
import java.util.Hashtable;

//import org.dcm4che2.imageio.*;
import org.dcm4che2.io.DicomInputStream;
//import org.dcm4che.data.*;
//import javax.imageio.*;

import uk.ac.ucl.chem.ccs.aheclient.res.AdvancedReservation;
import uk.ac.ucl.chem.ccs.clinicalgui.res.CreateReservation;
import uk.ac.ucl.chem.ccs.clinicalgui.res.ResPanel;
import uk.ac.ucl.chem.ccs.clinicalgui.res.SelectResource;

//import uk.ac.ucl.chem.ccs.clinicalgui.MainPanel.NextListener;


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
public class SimulationLaunchPanel extends javax.swing.JPanel {
	private JScrollPane treeView;
	private JLabel infoLabel;
	private JEditorPane noteLabel;
	private JButton LaunchHemelbButton;
	public JButton PrevButton;
	public JButton NextButton;
	private JPanel reservationInfoPanel;
	private JLabel resHeader;
	private JLabel resName;
	private JLabel resStart;
	private JLabel resEnd;
	private JLabel m1;
	private JLabel m2;
	private JLabel m3;
	private JLabel m4;


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
	
	public SimulationLaunchPanel() {
		super();
		initGUI(null);
	}
	
	private void updateInfo(String m){
		infoLabel.setText(m);
	}
	
	public void createDataTree(String[] highlightedPath){
		try{
			//DefaultMutableTreeNode top = GridServerInterface.populateServerData();
			JTree tree = GridServerInterface.populateServerData(highlightedPath);
			//tree.setExpandsSelectedPaths(true);
			tree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
			tree.addTreeSelectionListener(this.new myTreeListener());
			treeView = new JScrollPane(tree);
			this.add(treeView, "0, 0, 1, 3");
		}
		catch(Exception ex){
			updateInfo("could not load data from the grid");
		}
	}
	

	public class myTreeListener implements TreeSelectionListener{
		public void valueChanged(TreeSelectionEvent e) {
			String selectedPath = e.getPath().toString();
			String input = selectedPath.substring(1,selectedPath.length() - 1);
			String[] data = input.split(",");
			//System.out.println("length of path is " + data.length);
			if(data.length != 6)
				return;
			else if(!data[5].trim().equals("notes.txt"))
				return;
			else{
				try{
					noteLabel.setText("Note: " + GridServerInterface.getModelNote(data[1].trim(),data[2].trim(),data[3].trim(),data[4].trim()));
				}
				catch(Exception ex){
					updateInfo("retrieval of note unsuccessful");
				}
				updateInfo("successfully retrieved note");
			}
		}
	}
	
	public void updateReservationInfo(){
		if(ClinicalGuiClient.reservation != null){
			resName.setText("Name: " + ClinicalGuiClient.reservation.getResName());
			resStart.setText("Start: " + ClinicalGuiClient.reservation.getStart().toString());
			resEnd.setText("End: " + ClinicalGuiClient.reservation.getEnd().toString());
			m1.setText("Resources:");
			Hashtable res = ClinicalGuiClient.reservation.getRes();
			Enumeration myEnum = res.keys();
			if(myEnum.hasMoreElements())
				m2.setText((String) myEnum.nextElement());
			if(myEnum.hasMoreElements())
				m3.setText((String) myEnum.nextElement());
			if(res.size() > 3)
				m4.setText("...");
			else if(myEnum.hasMoreElements())
				m4.setText((String) myEnum.nextElement());
			//ClinicalGuiClient.reservation.getRes().
		}
	}
	
	public void initGUI(String[] highlightedPath) {
		try {
			
				TableLayout thisLayout = new TableLayout(new double[][] {{TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}, {TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}});
				thisLayout.setHGap(5);
				thisLayout.setVGap(5);
				this.setLayout(thisLayout);
				setPreferredSize(new Dimension(400, 300));
			    
			    createDataTree(highlightedPath);
				
				
				
				
			
			    infoLabel = new JLabel();
			    this.add(infoLabel, "0, 4, 3, 4");
			    infoLabel.setText("select a data set above");
			    
			    noteLabel = new JEditorPane();
			    this.add(noteLabel, "2, 0, 3, 1");
			    
			    reservationInfoPanel = new JPanel();
			    this.add(reservationInfoPanel, "2, 2, 3, 3");
			    reservationInfoPanel.setLayout(new TableLayout(new double[][] {{TableLayout.FILL}, {TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}}));
			    
			    resHeader = new JLabel();
			    reservationInfoPanel.add(resHeader,"0,0");
			    resHeader.setText("Reservation Details:");
			    
			    resName = new JLabel();
			    reservationInfoPanel.add(resName,"0,1");
			    
			    resStart = new JLabel();
			    reservationInfoPanel.add(resStart,"0,2");
			    
			    resEnd = new JLabel();
			    reservationInfoPanel.add(resEnd,"0,3");
			    
			    m1 = new JLabel();
			    reservationInfoPanel.add(m1,"0,4");
			    
			    m2= new JLabel();
			    reservationInfoPanel.add(m2,"0,5");
			    
			    m3 = new JLabel();
			    reservationInfoPanel.add(m3,"0,6");
			    
			    m4 = new JLabel();
			    reservationInfoPanel.add(m4,"0,7");
			    
			    
				
			    NextButton = new JButton();
				this.add(NextButton, "3, 5");
				NextButton.setName("NextButton");
			
			
				PrevButton = new JButton();
				this.add(PrevButton, "0, 5");
				PrevButton.setName("PrevButton");
				
			
				LaunchHemelbButton = new JButton();
				this.add(LaunchHemelbButton, "1, 5, 2, 5");
				LaunchHemelbButton.setName("LaunchHemelbButton");		      
				LaunchHemelbButton.addActionListener(new ActionListener() {
					public void actionPerformed (ActionEvent evt){
						if(ClinicalGuiClient.reservation == null){
							SelectResource sr = new SelectResource(SimulationLaunchPanel.this.getTopLevelAncestor());
							sr.showDialog();
						}
					}
				});
				
				
				updateReservationInfo();
				
				
			Application.getInstance().getContext().getResourceMap(getClass()).injectComponents(this);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
