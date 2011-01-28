package uk.ac.ucl.chem.ccs.clinicalgui;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.DICOM;
import info.clearthought.layout.TableLayout;

import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;
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

//import org.dcm4che2.imageio.*;
import org.dcm4che2.io.DicomInputStream;
//import org.dcm4che.data.*;
//import javax.imageio.*;



//import uk.ac.ucl.chem.ccs.clinicalgui.MainPanel.NextListener;


/**
 * @author Konstantin Voevodski
 *
 * implements the UI compoenents and functionality of the "Model Generation" panel
 */
public class ModelGenerationPanel extends javax.swing.JPanel {
	private JScrollPane treeView;
	//scrollabe window containing the DICOM data tree
	private JTree tree;
	//the tree object representing DICOM data available on server
	private JPanel attributePanel;
	//panel containing l1,l2,...,l18 labels
	private JLabel infoLabel;
	//displays the status of the GUI (error messages / status)
	private JButton LaunchGraphClientButton;
	//button to launch the segmentation tool
	private JButton GetDicomDataButton;
	//button to download DICOM data from server
	private JButton ViewDicomDataButton;
	//button to view the downloaded DICOM slices inside a scrollable window
	public JButton PrevButton;
	//button to go back to previous panel
	public JButton NextButton;
	//button to go forward to next panel
	private JLabel l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18;
	//labels listing the attributes of the DICOM data
	private SimulationLaunchPanel companionPanel;
	//the next panel whose display needs to be updated when user uploads hemlb input (seg tool output) data from this panel
	private MainPanel parentPanel;
	//the parent panel (needed to disable the tab panel of the parent)
	
	
	//****
	public static String modelNote;
	//holds the note entered for a particular model: set by NoteDialog
	private String selectedPath;
	//the currently selected path of the tree displaying the DICOM data
	private String patientId, studyId, seriesId;
	//the patient,study,series ids of currently selected data set
	public static ImageStack dicomImageStack;
	//the stack of dicom images used to create the scrollable window : set by DicomServerInterface
	private boolean dicomDataLoaded = false;
	//true if DICOM data has been downloaded
    //*****
	
	private String tmpdir = ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.dicomtmpdir");
	
	private String SEG_TOOL_PATH =  ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.segtool");
	//path to seg tool
	private String CONFIG_FILE_PATH = tmpdir +"/OutputConfig.dat";
	//path to output file of seg tool (used as input to hemelb)
	private String PARS_FILE_PATH = tmpdir +"/OutputPars.asc";
	//path of output file of seg tool (used as input to hemelb)
	private String CHECKPOINT_FILE_PATH = tmpdir +"/CheckPoint.dat";
	//path of output file of seg tool
	private String NOTE_FILE_PATH = tmpdir +"/notes.txt";
	//path to file containing note entered by user
	
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
	
	public void setParentPanel(MainPanel p){
		parentPanel = p;
	}
	public void setCompanionPanel(SimulationLaunchPanel p){
		companionPanel = p;
	}
	
	private void updateInfo(String m){
		infoLabel.setText(m);
	}
	
	private void disableUI(){
		PrevButton.setEnabled(false);
		GetDicomDataButton.setEnabled(false);
		ViewDicomDataButton.setEnabled(false);
		LaunchGraphClientButton.setEnabled(false);
		NextButton.setEnabled(false);
		parentPanel.disableTabPane();
		tree.setEnabled(false);
		System.out.println("disabled gui");
	}
	
	private void enableUI(){
		PrevButton.setEnabled(true);
		GetDicomDataButton.setEnabled(true);
		if(dicomDataLoaded){
			ViewDicomDataButton.setEnabled(true);
			LaunchGraphClientButton.setEnabled(true);
		}
		NextButton.setEnabled(true);
		parentPanel.enableTabPane();
		tree.setEnabled(true);
		System.out.println("enabled gui");
	}
	
	public class getDataListener implements ActionListener{
		public void actionPerformed (ActionEvent evt){
			if(selectedPath == null){
				updateInfo("you need to select a data set");
				return;
			}
			String input = selectedPath.substring(1,selectedPath.length() - 1);
			String[] data = input.split(",");
			if(data.length != 4){
				updateInfo("data set not fully specified: you need to select a series");
				return;
			}
			String patientId = data[1].trim();
			String studyId = data[2].trim();
			String seriesId = data[3].trim();
			if(dicomDataLoaded && ModelGenerationPanel.this.patientId.equals(patientId) && ModelGenerationPanel.this.studyId.equals(studyId)
		      && ModelGenerationPanel.this.seriesId.equals(seriesId)){
				updateInfo("Data already uploaded");
				return;
			}
			ModelGenerationPanel.this.patientId = patientId;
			ModelGenerationPanel.this.studyId = studyId;
			ModelGenerationPanel.this.seriesId = seriesId;
			updateInfo("Fetching DICOM data from server");
			ModelGenerationPanel.this.parentPanel.setCursor(new Cursor(Cursor.WAIT_CURSOR));
			disableUI();
			new Thread(){
				public void run() {
					int filesReceived = -1;
					try{
						filesReceived = DicomServerInterface.queryReceiveWrite(ModelGenerationPanel.this.patientId, ModelGenerationPanel.this.studyId, ModelGenerationPanel.this.seriesId);
					}
					catch(IOException ex){
						updateInfo("error writing received files");
						dicomDataLoaded = false;
					}
					if(filesReceived > 0){
						dicomDataLoaded = true;
						updateInfo("Successfully received " + filesReceived + " DICOM files");
					}
					else{
						updateInfo("The query did not retrieve any data");
						dicomDataLoaded = false;
					}
					enableUI();
					ModelGenerationPanel.this.parentPanel.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
			}}.start();		
			
		}
	}
	
	public class launchSegToolListener implements ActionListener{
		public void actionPerformed (ActionEvent evt){
			updateInfo("Launching segmentation tool");
			ModelGenerationPanel.this.parentPanel.setCursor(new Cursor(Cursor.WAIT_CURSOR));
			disableUI();
			new Thread(){
				public void run() {
					try{
						String cmd = SEG_TOOL_PATH + " " + DicomServerInterface.SLICES_FILE_PATH + " " + CONFIG_FILE_PATH + " " +  PARS_FILE_PATH + " " + CHECKPOINT_FILE_PATH + " 1.0 1.0 1";
					 	Process p = Runtime.getRuntime().exec(cmd);
					    p.waitFor();
					 }
					 catch(Exception ex){
					 	updateInfo("Error launching segmentation tool");
					 	return;
					 }
					 modelNote = "";
					 NoteDialog nt = new NoteDialog((JFrame) ModelGenerationPanel.this.getTopLevelAncestor());
				     nt.setModal(true);
					 nt.setVisible(true);
					 //when NoteDialog is done modelNote is note entered by user
					 ModelGenerationPanel.this.updateInfo("Uploading segmentation tool output to the grid");
					 try{
						 	BufferedWriter writer = new BufferedWriter(new FileWriter((NOTE_FILE_PATH)));
						 	writer.write(modelNote);
						 	writer.newLine();
						 	writer.close();
						 	String dateTime = new java.text.SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new java.util.Date());
						    GridServerInterface.uploadFiles(CONFIG_FILE_PATH, PARS_FILE_PATH, NOTE_FILE_PATH, patientId, studyId, seriesId, dateTime);
						    String[] selectedPath = new String[4];
						    selectedPath[0] = patientId;
						    selectedPath[1] = studyId;
						    selectedPath[2] = seriesId;
						    selectedPath[3] = dateTime;
							companionPanel.initGUI(selectedPath);
							parentPanel.addButtons(companionPanel);
							parentPanel.advance();
					 }
					 catch(Exception ex){
						 updateInfo("Error uploading segmentation tool output to the grid");
						 return;
					 }
					 updateInfo("Segmentation tool output has been successfully uploaded to the grid");
					 ModelGenerationPanel.this.parentPanel.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
					 enableUI();
			}}.start();		
			 
		}
	}
	
	private void initGUI() {
		try {
			
				TableLayout thisLayout = new TableLayout(new double[][] {{TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}, {TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}});
				thisLayout.setHGap(5);
				thisLayout.setVGap(5);
				this.setLayout(thisLayout);
				setPreferredSize(new Dimension(400, 300));
			
				attributePanel = new JPanel();
				TableLayout jPanelLayout = new TableLayout(new double[][] {{TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}, {TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}});
				jPanelLayout.setHGap(5);
				jPanelLayout.setVGap(5);
				this.add(attributePanel, "2, 0, 4, 3");
				attributePanel.setLayout(jPanelLayout);
				
				//these labels display the series attribute type 
				l1 = new JLabel();
				l1.setText("Patient Sex");
				l2 = new JLabel();
				l2.setText("Patient DOB");
				l3 = new JLabel();
				l3.setText("Study Description");
				l4 = new JLabel();
				l4.setText("Referring Physician");
				l5 = new JLabel();
				l5.setText("Series Date");
				l6 = new JLabel();
				l6.setText("Series Time");
				l7 = new JLabel();
				l7.setText("Series Description");
				l8 = new JLabel();
				l8.setText("Institution Name");
				l9 = new JLabel();
				l9.setText("Performing Physician");
				attributePanel.add(l1,"0,0");
				attributePanel.add(l2,"0,1");
				attributePanel.add(l3,"0,2");
				attributePanel.add(l4,"0,3");
				attributePanel.add(l5,"0,4");
				attributePanel.add(l6,"0,5");
				attributePanel.add(l7,"0,6");
				attributePanel.add(l8,"0,7");
				attributePanel.add(l9,"0,8");
				//these labels will display the attributes of the series selected by user
				l10 = new JLabel();
				l11 = new JLabel();
				l12 = new JLabel();
				l13 = new JLabel();
				l14 = new JLabel();
				l15 = new JLabel();
				l16 = new JLabel();
				l17 = new JLabel();
				l18 = new JLabel();
				attributePanel.add(l10,"1,0,2,0");
				attributePanel.add(l11,"1,1,2,1");
				attributePanel.add(l12,"1,2,2,2");
				attributePanel.add(l13,"1,3,2,3");
				attributePanel.add(l14,"1,4,2,4");
				attributePanel.add(l15,"1,5,2,5");
				attributePanel.add(l16,"1,6,2,6");
				attributePanel.add(l17,"1,7,2,7");
				attributePanel.add(l18,"1,8,2,8");
				
				tree = DicomServerInterface.populateServerData();
				treeView = new JScrollPane(tree);
				this.add(treeView, "0, 0, 1, 3");
				
				tree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
				tree.addTreeSelectionListener(new TreeSelectionListener(){
					public void valueChanged(TreeSelectionEvent e) {
						selectedPath = e.getPath().toString();
						String input = selectedPath.substring(1,selectedPath.length() - 1);
						String[] data = input.split(",");
						if(data.length == 2 ){
							String patientId = data[1].trim();
							String[] patientMetadata = DicomServerInterface.getPatientAttributes(patientId);
							l10.setText(patientMetadata[0]);
							l11.setText(patientMetadata[1]);
							l12.setText("");
							l13.setText("");
							l14.setText("");
							l15.setText("");
							l16.setText("");
							l17.setText("");
							l18.setText("");
						}
						else if(data.length == 3){
							String patientId = data[1].trim();
							String studyId = data[2].trim();
							String[] patientMetadata = DicomServerInterface.getPatientAttributes(patientId);
							l10.setText(patientMetadata[0]);
							l11.setText(patientMetadata[1]);
							String[] studyMetadata = DicomServerInterface.getStudyAttributes(patientId,studyId);
							l12.setText(studyMetadata[0]);
							l13.setText(studyMetadata[1]);
							l14.setText("");
							l15.setText("");
							l16.setText("");
							l17.setText("");
							l18.setText("");
						}
						else if(data.length == 4){
							String patientId = data[1].trim();
							String studyId = data[2].trim();
							String seriesId = data[3].trim();
							String[] patientMetadata = DicomServerInterface.getPatientAttributes(patientId);
							l10.setText(patientMetadata[0]);
							l11.setText(patientMetadata[1]);
							String[] studyMetadata = DicomServerInterface.getStudyAttributes(patientId,studyId);
							l12.setText(studyMetadata[0]);
							l13.setText(studyMetadata[1]);
							String[] seriesMetadata = DicomServerInterface.getSeriesAttributes(patientId,studyId,seriesId);
							l14.setText(seriesMetadata[0]);
							l15.setText(seriesMetadata[1]);
							l16.setText(seriesMetadata[2]);
							l17.setText(seriesMetadata[3]);
							l18.setText(seriesMetadata[4]);
						}
					}
				});
			
			    infoLabel = new JLabel();
			    this.add(infoLabel, "0, 4, 4, 4");
			    infoLabel.setText("select a data set above");
				
			    NextButton = new JButton();
				this.add(NextButton, "4, 5");
				NextButton.setName("NextButton");
			
			
				PrevButton = new JButton();
				this.add(PrevButton, "0, 5");
				PrevButton.setName("PrevButton");
				
				GetDicomDataButton = new JButton();
				this.add(GetDicomDataButton, "1, 5");
				GetDicomDataButton.setName("LaunchDataGenerationButton");		      
				GetDicomDataButton.addActionListener(this.new getDataListener());
				
				ViewDicomDataButton = new JButton();
				this.add(ViewDicomDataButton, "2, 5");
				ViewDicomDataButton.setName("ViewDicomDataButton");	
				ViewDicomDataButton.setEnabled(false);
				ViewDicomDataButton.addActionListener(new ActionListener() {
					public void actionPerformed (ActionEvent evt){
							new ImagePlus("Image Slices",dicomImageStack).show();
							 
					}
				});
		
			
				LaunchGraphClientButton = new JButton();
				this.add(LaunchGraphClientButton, "3, 5");
				LaunchGraphClientButton.setName("LaunchGraphClientButton");		      
				LaunchGraphClientButton.setEnabled(false);
				LaunchGraphClientButton.addActionListener(this.new launchSegToolListener());
				
				
				
				
			Application.getInstance().getContext().getResourceMap(getClass()).injectComponents(this);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
