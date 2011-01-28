package uk.ac.ucl.chem.ccs.clinicalgui;

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

import info.clearthought.layout.TableLayout;

import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.Rectangle;
import java.util.Hashtable;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.io.File;
import javax.swing.table.AbstractTableModel;
import javax.swing.JFileChooser;
import javax.swing.Timer;
import javax.swing.border.LineBorder;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.ProgressMonitor;
import javax.swing.border.TitledBorder;
import javax.swing.table.TableModel;
import javax.swing.BorderFactory;
//import org.realitygrid.steerer.Steerer;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import uk.ac.ucl.chem.ccs.aheclient.util.GridSAMStateInfo;
import uk.ac.ucl.chem.ccs.aheclient.util.JobFileElement;
import uk.ac.ucl.chem.ccs.aheclient.util.JobRegistryElement;
import uk.ac.ucl.chem.ccs.aheclient.util.StageFilesIn;
import uk.ac.ucl.chem.ccs.aheclient.util.Tools;
import java.util.Vector;
import java.util.Iterator;
import java.lang.Integer;
import uk.ac.ucl.chem.ccs.aheclient.util.AHEJobObject;
import uk.ac.ucl.chem.ccs.aheclient.wsrf.MonitorSimCall;
import uk.ac.ucl.chem.ccs.aheclient.wsrf.StartCall;
import uk.ac.ucl.chem.ccs.aheclient.wsrf.TerminateSimCall;
import uk.ac.ucl.chem.ccs.vizclient.HttpGetPoll;
import uk.ac.ucl.chem.ccs.vizclient.VizSteererWindow;

import java.awt.GraphicsEnvironment;



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
public class DisplayJobPanel extends javax.swing.JPanel {
        private JPanel jPanel1;
        private JPanel gridsamStatus;
        private JPanel controls;
        private JLabel jLabel4;
        private JLabel jLabel8;
        private JButton changeLocationButton;
        private JButton downloadButton;
        private JTable outputFilesTable;
        private JScrollPane filesScrollPane;
        private JTextField jobEPR;
        private JTextField jobStdin;
        private JTextField jobStderr;
        private JTextField jobSdtout;
        private JTextField jobArgs;
        private JTextField jobConf;
        private JTextField jobCpus;
        private JTextField rm;
        private JLabel jobStatus;
        private JTextField jobType;
        private JTextField resourceID;
        private JTextField jobName;
        private JLabel jLabel13;
        private JLabel jLabel12;
        private JLabel jLabel11;
        private JLabel jLabel10;
        private JLabel jLabel9;
        private JLabel jLabel7;
        private JLabel jLabel6;
        private JLabel jLabel5;
        private JLabel jLabel3;
        private JLabel jLabel2;
        private JPanel jPanel3;
        private JScrollPane jobDetailsSP;
        private JCheckBox deleteFiles;
        private JButton vizButton;
        private JButton teminateJob;
        private JButton updateStatus;
        private JButton pollingButton;
        private JLabel jLabel1;
        private JSlider jSlider1;
        private JPanel polling;
        private JLabel time1;
        private JLabel jLabel14;
        private JLabel time;
        private JList detailsList;
        private JPanel dtails;
        private JTextArea gridsamStatusResults;
        private JScrollPane jScrollPane1;
        private JPanel stagedFiles;
        private JTabbedPane jTabbedPane1;
        private JPanel jPanel2;
        private AHEJobObject ajo;
        private Timer pollTimer;
        private Vector downloadFiles;
        private Rectangle closeBut;
        //private JobRegistryElement jre = null;
        private JPanel regSteering;
        private String fileLocation = null;
    	private boolean steeredApp = false;
    	private VizSteererWindow vs;
    	private JButton steer;
    	private String rID;
    	private JTextField steerERP;
    	private String h = "localhost";
    	private int p = 65250;
    	private int w = 1024*1024;
        
        private static Log cat = LogFactory.getLog(DisplayJobPanel.class);

       
        public DisplayJobPanel(AHEJobObject ajo) {
                super();
                this.ajo = ajo;
                //this.jre = jre;
                downloadFiles = new Vector();
                
                //hemelb should use a string
        		rID=ajo.getResourceID();
        		rID = rID.substring(0, 4);
        		if (rID.startsWith("0")) {
        			rID = 1+rID;
        		}
                
                initGUI();
                
                //look for registered app
        		if (steeredApp) {
        			PollThread pt = new PollThread();
        			pt.start();
        		}
        }
        
        public void setJobObject(AHEJobObject job){
        	ajo = job;
        	initGUI();
        }
       
        private void initGUI() {
        	    if(ajo == null){
        	    	try {
        				setPreferredSize(new Dimension(400, 300));
        			} catch (Exception e) {
        				e.printStackTrace();
        			}
        			//JLabel l = new JLabel("No simulation running");
        			//this.add(l);
        			//this.setEnabled(false);
        			return;
        	    }
        	    System.out.println("I am not null: drawing the display job panel");
                try {
                        TableLayout thisLayout = new TableLayout(new double[][] {
                                        { TableLayout.FILL },
                                        { TableLayout.PREFERRED, TableLayout.PREFERRED,
                                                        TableLayout.PREFERRED, TableLayout.FILL,
                                                        TableLayout.FILL } });
                        thisLayout.setHGap(5);
                        thisLayout.setVGap(5);
                        this.setLayout(thisLayout);

                       
                        {
                                jPanel1 = new JPanel();
                                TableLayout jPanel1Layout = new TableLayout(new double[][] {
                                                { TableLayout.FILL, TableLayout.FILL, TableLayout.FILL,
                                                                TableLayout.PREFERRED, TableLayout.PREFERRED },
                                                { TableLayout.FILL, TableLayout.FILL, TableLayout.FILL,
                                                                TableLayout.FILL, TableLayout.FILL } });

                                jPanel1Layout.setHGap(5);
                                jPanel1Layout.setVGap(5);
                                jPanel1.setLayout(jPanel1Layout);
                                jPanel1.setLayout(jPanel1Layout);

                                jPanel1.setBorder(BorderFactory.createEtchedBorder());
                                this.add(jPanel1, "0, 0, 0, 2");
                                jPanel1.setPreferredSize(new java.awt.Dimension(630, 305));
                                jPanel1.setSize(630, 305);
                                {
                                        dtails = new JPanel();
                                        GridLayout dtailsLayout = new GridLayout(1, 1);
                                        dtailsLayout.setColumns(1);
                                        dtailsLayout.setHgap(5);
                                        dtailsLayout.setVgap(5);
                                        dtails.setBorder(BorderFactory.createTitledBorder("Job Details"));
                                        jPanel1.add(dtails, "0,  0,  2,  4");
                                        dtails.setLayout(dtailsLayout);
                                        {
                                                jobDetailsSP = new JScrollPane();
                                                dtails.add(jobDetailsSP);
                                                {
                                                        jPanel3 = new JPanel();
                                                        TableLayout jPanel3Layout = new TableLayout(
                                                                new double[][] {
                                                                                { TableLayout.PREFERRED, TableLayout.FILL },
                                                                                { TableLayout.FILL, TableLayout.FILL,
                                                                                                TableLayout.FILL,
                                                                                                TableLayout.FILL,
                                                                                                TableLayout.FILL,
                                                                                                TableLayout.FILL,
                                                                                                TableLayout.FILL,
                                                                                                TableLayout.FILL,
                                                                                                TableLayout.FILL,
                                                                                                TableLayout.FILL,
                                                                                                TableLayout.FILL,
                                                                                                TableLayout.FILL } });
                                                        jPanel3Layout.setHGap(5);
                                                        jPanel3Layout.setVGap(5);
                                                        jPanel3.setLayout(jPanel3Layout);
                                                        jobDetailsSP.setViewportView(jPanel3);
                                                        jPanel3.setBackground(new java.awt.Color(156,199,219));
                                                        {
                                                                jLabel2 = new JLabel();
                                                                jPanel3.add(jLabel2, "0, 0");
                                                                jLabel2.setText("Job Start Time");
                                                        }
                                                        {
                                                                jLabel3 = new JLabel();
                                                                jPanel3.add(jLabel3, "0, 1");
                                                                jLabel3.setText("Resource ID");
                                                        }
                                                        {
                                                                jLabel4 = new JLabel();
                                                                jPanel3.add(jLabel4, "0, 2");
                                                                jLabel4.setText("Job Type");
                                                        }
                                                        {
                                                                jLabel5 = new JLabel();
                                                                jPanel3.add(jLabel5, "0, 3");
                                                                jLabel5.setText("Status");
                                                        }
                                                        {
                                                                jLabel6 = new JLabel();
                                                                jPanel3.add(jLabel6, "0, 4");
                                                                jLabel6.setText("Machine");
                                                        }
                                                        {
                                                                jLabel7 = new JLabel();
                                                                jPanel3.add(jLabel7, "0, 5");
                                                                jLabel7.setText("CPUs Requested");
                                                        }
                                                        {
                                                                jLabel8 = new JLabel();
                                                                jPanel3.add(jLabel8, "0, 6");
                                                                jLabel8.setText("Configuration File");
                                                        }
                                                        {
                                                                jLabel9 = new JLabel();
                                                                jPanel3.add(jLabel9, "0, 7");
                                                                jLabel9.setText("Job Arguments");
                                                        }
                                                        {
                                                                jLabel10 = new JLabel();
                                                                jPanel3.add(jLabel10, "0, 8");
                                                                jLabel10.setText("Job Stdout");
                                                        }
                                                        {
                                                                jLabel11 = new JLabel();
                                                                jPanel3.add(jLabel11, "0, 9");
                                                                jLabel11.setText("Job Stderr");
                                                        }
                                                        {
                                                                jLabel12 = new JLabel();
                                                                jPanel3.add(jLabel12, "0, 10");
                                                                jLabel12.setText("Job Stdin");
                                                        }
                                                        {
                                                                jLabel13 = new JLabel();
                                                                jPanel3.add(jLabel13, "0, 11");
                                                                jLabel13.setText("Resource Endpoint");
                                                        }
                                                        {
                                                                jobName = new JTextField();
                                                                jPanel3.add(jobName, "1, 0");
                                                                jobName.setBorder(new LineBorder(new java.awt.Color(0,0,0),1,false));
                                                                jobName.setOpaque(true);
                                                                jobName.setBackground(new java.awt.Color(255,255,255));
                                                                jobName.setEditable(false);
                                               
                                                        }
                                                        {
                                                                resourceID = new JTextField();
                                                                jPanel3.add(resourceID, "1, 1");
                                                                resourceID.setBorder(new LineBorder(new java.awt.Color(0,0,0),1,false));
                                                                resourceID.setOpaque(true);
                                                                resourceID.setBackground(new java.awt.Color(255,255,255));
                                                                resourceID.setEditable(false);
                                                               
                                                        }
                                                        {
                                                                jobType = new JTextField();
                                                                jPanel3.add(jobType, "1, 2");
                                                                jobType.setBorder(new LineBorder(new java.awt.Color(0,0,0),1,false));
                                                                jobType.setBackground(new java.awt.Color(255,255,255));
                                                                jobType.setOpaque(true);
                                                                jobType.setEditable(false);
                                                               
                                                        }
                                                        {
                                                                jobStatus = new JLabel();
                                                                jPanel3.add(jobStatus, "1, 3");
                                                                jobStatus.setOpaque(true);
                                                                jobStatus.setBorder(new LineBorder(new java.awt.Color(0,0,0), 1, false));
                                                                jobStatus.setBackground(new java.awt.Color(255,255,255));
                                                        }
                                                        {
                                                                rm = new JTextField();
                                                                jPanel3.add(rm, "1, 4");
                                                                rm.setBackground(new java.awt.Color(255,255,255));
                                                                rm.setOpaque(true);
                                                                rm.setBorder(new LineBorder(new java.awt.Color(0,0,0),1,false));
                                                                rm.setEditable(false);
                                                               
                                                        }
                                                        {
                                                                jobCpus = new JTextField();
                                                                jPanel3.add(jobCpus, "1, 5");
                                                                jobCpus.setBorder(new LineBorder(new java.awt.Color(0,0,0),1,false));
                                                                jobCpus.setOpaque(true);
                                                                jobCpus.setBackground(new java.awt.Color(255,255,255));
                                                                jobCpus.setEditable(false);
                                                       
                                                        }
                                                        {
                                                                jobConf = new JTextField();
                                                                jPanel3.add(jobConf, "1, 6");
                                                                jobConf.setBackground(new java.awt.Color(255,255,255));
                                                                jobConf.setOpaque(true);
                                                                jobConf.setBorder(new LineBorder(new java.awt.Color(0,0,0),1,false));
                                                                jobConf.setEditable(false);
                                                               
                                                        }
                                                        {
                                                                jobArgs = new JTextField();
                                                                jPanel3.add(jobArgs, "1, 7");
                                                                jobArgs.setBackground(new java.awt.Color(255,255,255));
                                                                jobArgs.setOpaque(true);
                                                                jobArgs.setBorder(new LineBorder(new java.awt.Color(0,0,0),1,false));
                                                                jobArgs.setEditable(false);
                                                               
                                                        }
                                                        {
                                                                jobSdtout = new JTextField();
                                                                jPanel3.add(jobSdtout, "1, 8");
                                                                jobSdtout.setBackground(new java.awt.Color(255,255,255));
                                                                jobSdtout.setOpaque(true);
                                                                jobSdtout.setBorder(new LineBorder(new java.awt.Color(0,0,0),1,false));
                                                                jobSdtout.setEditable(false);
                                                               
                                                        }
                                                        {
                                                                jobStderr = new JTextField();
                                                                jPanel3.add(jobStderr, "1, 9");
                                                                jobStderr.setBackground(new java.awt.Color(255,255,255));
                                                                jobStderr.setOpaque(true);
                                                                jobStderr.setBorder(new LineBorder(new java.awt.Color(0,0,0),1,false));
                                                                jobStderr.setEditable(false);
                                                               
                                                        }
                                                        {
                                                                jobStdin = new JTextField();
                                                                jPanel3.add(jobStdin, "1, 10");
                                                                jobStdin.setBorder(new LineBorder(new java.awt.Color(0,0,0),1,false));
                                                                jobStdin.setBackground(new java.awt.Color(255,255,255));
                                                                jobStdin.setOpaque(true);
                                                                jobStdin.setEditable(false);
                                                               
                                                        }
                                                        {
                                                                jobEPR = new JTextField();
                                                                jPanel3.add(jobEPR, "1, 11");
                                                                jobEPR.setOpaque(true);
                                                                jobEPR.setBorder(new LineBorder(new java.awt.Color(0,0,0),1,false));
                                                                jobEPR.setBackground(new java.awt.Color(255,255,255));
                                                                jobEPR.setEditable(false);
                                                       
                                                        }
                                                }

                                        }
                                }
                                {
                                        controls = new JPanel();
                                        TableLayout controlsLayout = new TableLayout(
                                                new double[][] {
                                                                { TableLayout.FILL },
                                                                { TableLayout.FILL, TableLayout.FILL,
                                                                                TableLayout.FILL, TableLayout.FILL } });
                                        controlsLayout.setHGap(5);
                                        controlsLayout.setVGap(5);
                                        controls.setLayout(controlsLayout);
                                        controls.setBorder(BorderFactory.createTitledBorder("Operations"));
                                        jPanel1.add(controls, "3,  0,  4,  2");
                                        {
                                                updateStatus = new JButton();
                                                controls.add(updateStatus, "0, 0");
                                                updateStatus.setText("Update Job Status");
                                                updateStatus.addActionListener(new ActionListener() {
                                                        public void actionPerformed (ActionEvent evt) {
                                                                if (updateStatus.getText().equals("Start Job")) {
                                                                StartCall sc = new StartCall(ajo,
                                                                                ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.myproxy-lifetime"),
                                                                                ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.myproxy-port"),
                                                                                ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.myproxy-dn"),
                                                                                ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.myproxy-server"),
                                                                                ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.myproxy-pw"),
                                                                                ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.myproxy-un"));
                                                                ajo = sc.makeCall();
                                                                updatePanel();
                                                        } else {
                                                                pollJobState();
                                                        }
                                                        }
                                                });
                                        }
                                        {
                                                teminateJob = new JButton();
                                                controls.add(teminateJob, "0, 1");
                                                teminateJob.setText("Terminate Job");
                                                teminateJob.addActionListener(new ActionListener() {
                                                        public void actionPerformed (ActionEvent evt) {
                                                                TerminateSimCall tsc = new TerminateSimCall (ajo.getEndPoint());
                                                                DisplayJobPanel.this.setCursor(new Cursor(Cursor.WAIT_CURSOR));
                                                                boolean tcsstatus = tsc.makeCall();
                                                                DisplayJobPanel.this.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
                                                                if (tcsstatus) {
                                                                        ajo.setState(AHEJobObject.GRIDSAM_TERMINATING);
                                                                        updateState ();
                                                                } else {
                                                                        ErrorMessage em = new ErrorMessage(DisplayJobPanel.this, "Error terminating job. Check log for details");;
                                                                }
                                                        }
                                                });
                                        }
                                        {
                                        		vizButton = new JButton();
                                        			 controls.add(vizButton, "0, 2");
                                                     vizButton.setText("Visualize");
                                                     
                                                     vizButton.addActionListener(new ActionListener() {
                                                         public void actionPerformed (ActionEvent evt) {
                                                        		String h = "localhost";
                                                        		int p = 65250;
                                                        		int w = 1024*1024;
                                                        		
                                                        		VizSteererWindow vs = new VizSteererWindow(h,p,w, (JFrame) DisplayJobPanel.this.getTopLevelAncestor());
                                                         }
                                                     });
                                        }
                                
                                        
                                        {
                                                deleteFiles = new JCheckBox();
                                                controls.add(deleteFiles, "0, 3");
                                                deleteFiles.setText("Delete staged files when destroying job");
                                                deleteFiles.setFont(new java.awt.Font("Sansserif",0,11));
                                                deleteFiles.setSelected(true);
                                        }
                                }
                                {
                                        polling = new JPanel();
                                        GridBagLayout pollingLayout = new GridBagLayout();
                                        pollingLayout.rowWeights = new double[] {0.1, 0.1, 0.1, 0.1};
                                        pollingLayout.rowHeights = new int[] {7, 7, 7, 7};
                                        pollingLayout.columnWeights = new double[] {0.0, 0.1};
                                        pollingLayout.columnWidths = new int[] {109, 7};
                                        polling.setBorder(BorderFactory.createTitledBorder("Status Polling"));
                                        jPanel1.add(polling, "3,  3,  4,  4");
                                        polling.setLayout(pollingLayout);
                                        {
                                                jLabel1 = new JLabel();
                                                polling.add(jLabel1, new GridBagConstraints(0, 0, 2, 1, 0.0, 0.0, GridBagConstraints.NORTHWEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
                                                jLabel1.setText("Set the polling interval ");
                                                jLabel1.setFont(new java.awt.Font("Sansserif",0,11));
                                        }
                                        {
                                                jSlider1 = new JSlider();
                                                polling.add(jSlider1, new GridBagConstraints(0, 1, 2, 2, 0.0, 0.0, GridBagConstraints.NORTHWEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
                                                jSlider1.setMaximum(60);
                                                jSlider1.setValue(0);
                                                //jSlider1.setMinorTickSpacing(1);
                                                //jSlider1.createStandardLabels(5);
                                                Hashtable lab = new Hashtable ();
                                                lab.put(new Integer(0), new JLabel("0"));
                                                lab.put(new Integer(20), new JLabel("10"));
                                                lab.put(new Integer(40), new JLabel("20"));
                                                lab.put(new Integer(60), new JLabel("30"));
                                                jSlider1.setLabelTable(lab);
                                               
                                                jSlider1.setPaintTicks(true);
                                                jSlider1.setPaintLabels(true);
                                                jSlider1.setSnapToTicks(false);
                                                jSlider1.setMajorTickSpacing(2);
                                                jSlider1.setFont(new java.awt.Font("Sansserif",0,11));
                                                jSlider1.addChangeListener(new ChangeListener () {
                                                        public void stateChanged(ChangeEvent e) {
                                                                if (jSlider1.getValue() != 0) {
                                                                        Integer i = new Integer(jSlider1.getValue());
                                                                time1.setText(Float.toString(i.floatValue()/2));
                                                                if(pollingButton.getText().equals("Stop Polling")) {
                                                                        pollTimer.stop();                                                                       
                                                                       
                                                                        pollTimer.setInitialDelay(jSlider1.getValue() * 30000);
                                                                        pollTimer.setDelay(jSlider1.getValue() * 30000);
                                                                        pollTimer.start();
                                                                }
                                                        } else {
                                                                if(pollingButton.getText().equals("Stop Polling")) {
                                                                        pollTimer.stop();
                                                                        pollingButton.setText("Start Polling");
                                                                }
                                                                time1.setText("0.0");
                                                        }
                                                        }
                                                });
                                        }
                                        {
                                                pollingButton = new JButton();
                                                polling.add(pollingButton, new GridBagConstraints(1, 3, 1, 1, 0.0, 0.0, GridBagConstraints.SOUTHEAST, GridBagConstraints.VERTICAL, new Insets(0, 0, 0, 0), 0, 0));
                                                pollingButton.setText("Start Polling");
                                                pollingButton.addActionListener(new ActionListener() {
                                                        public void actionPerformed (ActionEvent evt) {
                                                                if (jSlider1.getValue() > 0) {
                                                                if(pollingButton.getText().equals("Start Polling")) {
                                                                        pollTimer = new Timer(jSlider1.getValue() * 30000, new ActionListener() {
                                                                                public void actionPerformed (ActionEvent evt) {
                                                                                        pollJobState();
                                                                                }
                                                                                });
                                                                        pollTimer.setInitialDelay(1);
                                                                        pollTimer.start();
                                                                        pollingButton.setText("Stop Polling");
                                                                } else {
                                                                        pollTimer.stop();
                                                                        pollingButton.setText("Start Polling");
                                                                }
                                                                }
                                                                }
                                                });
                                               
                                        }
                                        {
                                                time = new JLabel();
                                                polling.add(time, new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.VERTICAL, new Insets(0, 0, 0, 0), 0, 0));
                                                time.setBackground(new java.awt.Color(255,255,255));
                                                time.setText("Every");
                                        }
                                        {
                                                jLabel14 = new JLabel();
                                                polling.add(jLabel14, new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0, GridBagConstraints.EAST, GridBagConstraints.VERTICAL, new Insets(0, 0, 0, 0), 0, 0));
                                                jLabel14.setText("mins");
                                        }
                                        {
                                                time1 = new JLabel();
                                                polling.add(time1, new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.VERTICAL, new Insets(0, 0, 0, 0), 0, 0));
                                                time1.setText("0.0");;
                                        }
                                }
                        }
                        {
                                jPanel2 = new JPanel();
                                GridLayout jPanel2Layout = new GridLayout(1, 1);
                                jPanel2Layout.setColumns(1);
                                jPanel2Layout.setHgap(5);
                                jPanel2Layout.setVgap(5);
                                jPanel2.setLayout(jPanel2Layout);
                                TitledBorder title2;
                                title2 = BorderFactory.createTitledBorder("Job Output");
                                jPanel2.setBorder(title2);
                                this.add(jPanel2, "0, 3, 0, 4");
                                jPanel2.setPreferredSize(new java.awt.Dimension(630, 254));
                                {
                                        jTabbedPane1 = new JTabbedPane();
                                        jPanel2.add(jTabbedPane1);

                                        {
                                                gridsamStatus = new JPanel();
                                                GridLayout gridsamStatusLayout = new GridLayout(1, 1);
                                                gridsamStatusLayout.setColumns(1);
                                                gridsamStatusLayout.setHgap(5);
                                                gridsamStatusLayout.setVgap(5);
                                                gridsamStatus.setLayout(gridsamStatusLayout);
                                                jTabbedPane1.addTab("AHE Job Status", null, gridsamStatus, null);
                                                {
                                                        jScrollPane1 = new JScrollPane();
                                                        gridsamStatus.add(jScrollPane1);
                                                        {
                                                                gridsamStatusResults = new JTextArea();
                                                                jScrollPane1
                                                                        .setViewportView(gridsamStatusResults);
                                                                gridsamStatusResults.setFont(new java.awt.Font("Monospaced",0,12));
                                                        }
                                                }
                                        }
                                        {
                                                stagedFiles = new JPanel();
                                                TableLayout stagedFilesLayout = new TableLayout(
                                                        new double[][] {
                                                                        { TableLayout.FILL, TableLayout.FILL,
                                                                                        TableLayout.FILL, TableLayout.FILL,
                                                                                        TableLayout.PREFERRED,
                                                                                        TableLayout.PREFERRED },
                                                                        { TableLayout.FILL, TableLayout.FILL,
                                                                                        TableLayout.FILL, TableLayout.FILL,
                                                                                        TableLayout.FILL,
                                                                                        TableLayout.PREFERRED } });
                                                stagedFilesLayout.setHGap(5);
                                                stagedFilesLayout.setVGap(5);
                                                stagedFiles.setLayout(stagedFilesLayout);
                                                jTabbedPane1.addTab(
                                                        "Staged Files",
                                                        null,
                                                        stagedFiles,
                                                        null);
                                                {
                                                        filesScrollPane = new JScrollPane();
                                                        stagedFiles.add(filesScrollPane, "0, 0, 5, 4");
                                                        {
                                                               
                                                                outputFilesTable = new JTable();
                                                               
                                                                int col1 = 0, col2 = 0;
                                                                int fsize = outputFilesTable.getFont().getSize() - 5;
                                                                Object data[][] = new Object [ajo.getOutfiles().size()+ ajo.getInfiles().size()][3];
                                                               
       
                                                               
                                                                int i = 0;
                                                                if (ajo.getOutfiles() != null) {
                                                                        Iterator it = ajo.getOutfiles().iterator();
                                                                        while (it.hasNext ()) {
                                                                                JobFileElement je = (JobFileElement)it.next();
                                                                                data [i][0] = new Boolean (true);
                                                                                data [i][1] = je.getName();
                                                                                if (je.getName().length() > col1){
                                                                                        col1 = je.getName().length();
                                                                                }
                                                                                String url = Tools.getUrlNoUP(je.getRemotepath());
                                                                                data [i][2] = url;
                                                                                if (url.length() > col2){
                                                                                        col2 = url.length();
                                                                                }
                                                                                i++;
                                                                        }
                                                                }
                                                               
                                                               
                                                                if (ajo.getInfiles() != null) {
                                                                        Iterator it = ajo.getInfiles().iterator();
                                                                        while (it.hasNext ()) {
                                                                                JobFileElement je = (JobFileElement)it.next();
                                                                                data [i][0] = new Boolean (false);
                                                                                data [i][1] = je.getName();
                                                                                if (je.getName().length() > col1){
                                                                                        col1 = je.getName().length();
                                                                                }
                                                                                String url = Tools.getUrlNoUP(je.getRemotepath());
                                                                                data [i][2] = url;
                                                                                if (url.length() > col2){
                                                                                        col2 = url.length();
                                                                                }
                                                                                i++;
                                                                        }
                                                                }
                                                               
                                                                String colNames [] = { "Download", "File Name",  "File Location"};
                                                               
                                                                TableModel outputFilesTableModel = new MyTableModel(data, colNames);
                                                                outputFilesTable.setIntercellSpacing(new Dimension (3,3));
                                                                outputFilesTable.setModel(outputFilesTableModel);
                                                                outputFilesTable.getColumnModel().getColumn(0).setPreferredWidth(70);
                                                                outputFilesTable.getColumnModel().getColumn(1).setPreferredWidth(col1*fsize);
                                                                outputFilesTable.getColumnModel().getColumn(2).setPreferredWidth(col2*fsize);
                                                                outputFilesTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
                                                                filesScrollPane.setViewportView(outputFilesTable);
                                                                this.addComponentListener(new ComponentAdapter () {
                                                                        public void componentResized(ComponentEvent e) {
                                                                                if (outputFilesTable.getWidth() < filesScrollPane.getWidth()) {
                                                                                        outputFilesTable.setAutoResizeMode(JTable.AUTO_RESIZE_LAST_COLUMN );
                                                                                } else {
                                                                                        outputFilesTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
                                                                                }
                                                                        }
                                                                });
                                                               
                                                        }
                                               
                                                }
                                                {
                                                        downloadButton = new JButton();
                                                        stagedFiles.add(downloadButton, "5, 5");
                                                        downloadButton.setText("Download");
                                                        downloadButton.addActionListener(new ActionListener() {
                                                                public void actionPerformed (ActionEvent evt) {
                                                                        int outFileSize = ajo.getOutfiles().size();
                                                                        for (int row = 0; row < outputFilesTable.getRowCount(); row++){
                                                                                if (((Boolean)outputFilesTable.getValueAt(row, 0)).booleanValue() == true) {
                                                                                       
                                                                                        if (row < outFileSize) {
                                                                                                JobFileElement je = (JobFileElement)ajo.getOutfiles().elementAt(row);
                                                                                               
                                                                                                if (fileLocation != null) {
                                                                                                        je.setLocalpath(Tools.checkURL(fileLocation) + je.getName());
                                                                                                }
                                                                                               
                                                                                                downloadFiles.add(je);
                                                                                        } else {
                                                                                                JobFileElement je = (JobFileElement)ajo.getInfiles().elementAt(row - outFileSize);
                                                                                               
                                                                                                if (fileLocation != null) {
                                                                                                        je.setLocalpath(Tools.checkURL(fileLocation) + je.getName());
                                                                                                }
                                                                                               
                                                                                                downloadFiles.add(je);
                                                                                        }
                                                                                }
                                                                               
                                                                        }
                                                                       
                                                                        StageFilesIn task = new StageFilesIn(ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.ahedavserver"),
                                                                                        ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.ahedavuser"),
                                                                                        ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.ahedavpasswd"));
                                                                        task.init(downloadFiles);
                                                                       
                                                                        ProgressMonitor progressMonitor = new ProgressMonitor(DisplayJobPanel.this,
                                                                        "Downloading Files",
                                                                        null, 0, task.getLength());
                                                                        //progressMonitor.setMillisToDecideToPopup(1);
                                                                        progressMonitor.setMillisToPopup(100);
                                                                        //jProgressBar1.setMaximum(task.getLength());
                                                                        //jProgressBar1.setValue(0);

                                                                        while(task.filesToStage()) {
                                                                                if (task.stageNext()) {
                                                                                        progressMonitor.setProgress(task.getCurrent());
                                                                                } else {
                                                                                        cat.error(task.getError());
                                                                                       
                                                                                }
                                                                               
                                                                        }
                                                                       
                                                                }
                                                        });
                                                                               
                                                }
                                                {
                                                        changeLocationButton = new JButton();
                                                        stagedFiles.add(changeLocationButton, "4, 5");
                                                        changeLocationButton.setText("Local Dir");
                                                        changeLocationButton.addActionListener(new ActionListener() {
                                                                public void actionPerformed (ActionEvent evt) {
                                                                        JFileChooser fc = new JFileChooser();       
                                                                        fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
                                                                        int returnVal = fc.showOpenDialog(DisplayJobPanel.this);

                                                                    if (returnVal == JFileChooser.APPROVE_OPTION) {
                                                                        File file = fc.getSelectedFile();
                                                                        fileLocation = file.getAbsolutePath();
                                                                        //System.out.println(fileLocation);
                                                                    }
                                                                }
                                                        });
                                                                               
                                                }
                                        }
                                        {
                                                if (ajo.getReGSWSEPR() != null) {
           
                                					regSteering = new JPanel();
                            						
                            						TableLayout steerLayout = new TableLayout(new double[][] {
                            								{ TableLayout.FILL },
                            								{ TableLayout.FILL, TableLayout.FILL,
                            										TableLayout.FILL, TableLayout.FILL,
                            										TableLayout.FILL } });
                            						regSteering.setLayout(steerLayout);
                            						steeredApp = true;
                            						jTabbedPane1.addTab("ReG Steering", null, regSteering, null);
                            						{
                            							JLabel look = new JLabel("Steering address");
                            							steerERP = new JTextField();
                            							steer = new JButton("Start Steerer");
                            							steer.setEnabled(false);	
                            							steer.addActionListener(new ActionListener() {
                            								public void actionPerformed (ActionEvent evt) {
                            									vs = new VizSteererWindow(h,p,w, DisplayJobPanel.this.getTopLevelAncestor());
                            								}});
                            								
                            							regSteering.add(look, "0,1");
                            							regSteering.add(steerERP, "0,2");
                            									regSteering.add(steer, "0,3");
                            							
                            							
                            						}
                                                	
                                                }
                                        }
                                       
                                }
                        }
                        updatePanel ();
                        this.setPreferredSize(new java.awt.Dimension(630, 605));
                        this.setSize(630, 605);
                        this.setOpaque(false);

                } catch (Exception e) {
                        e.printStackTrace();
                }
        }
       
        public void updatePanel () {
                if (ajo != null) {
                jobName.setText(ajo.getGridsam_Time());
                resourceID.setText(ajo.getResourceID());
                jobType.setText(ajo.getJobType());
                rm.setText(ajo.getRm());
                jobCpus.setText(Integer.toString(ajo.getCpuCount()));
                jobConf.setText(ajo.getConfFile());
                jobArgs.setText(ajo.getArgument());
                jobEPR.setText(ajo.getEndPoint());
                jobStdin.setText(ajo.getStdin());
                jobStderr.setText(ajo.getStderr());
                jobSdtout.setText(ajo.getStdout());
               
                updateState();
                updateGridSAMPanel ();
                } else {
                        ErrorMessage em = new ErrorMessage (this.getTopLevelAncestor(), "Error updating job. Check logs for details");
                }
                }
       
        private void updateState () {
               
                switch(ajo.getState()) {
                        case AHEJobObject.AHE_PREPARING:
                        jobStatus.setForeground(Color.BLUE);
                        jobStatus.setText("AHE PREPARING");
                        break;
                       
                        case AHEJobObject.AHE_FILES_STAGED:
                        jobStatus.setForeground(Color.ORANGE);
                        jobStatus.setText("AHE FILES STAGED");
                        break;
                       
                        case AHEJobObject.AHE_JOB_BUILT:
                        jobStatus.setForeground(Color.YELLOW);
                        jobStatus.setText("AHE JOB BUILT");
                        break;
                       
                        case AHEJobObject.GRIDSAM_PENDING:
                        jobStatus.setForeground(Color.MAGENTA);
                        jobStatus.setText("AHE PENDING");
                        break;
                       
                        case AHEJobObject.GRIDSAM_STAGING_IN:
                        jobStatus.setForeground(new Color(92, 201, 151));
                        jobStatus.setText("AHE STAGING IN");
                        break;
                       
                        case AHEJobObject.GRIDSAM_STAGED_IN:
                        jobStatus.setForeground(new Color(1, 240, 242));
                        jobStatus.setText("AHE STAGED IN");
                        break;
                       
                        case AHEJobObject.GRIDSAM_STAGING_OUT:
                        jobStatus.setForeground(new Color(188, 140, 217));
                        jobStatus.setText("AHE STAGING OUT");
                        break;
                       
                        case AHEJobObject.GRIDSAM_STAGED_OUT:
                        jobStatus.setForeground(new Color(135, 90, 133));
                        jobStatus.setText("AHE PREPARING");
                        break;
                       
                        case AHEJobObject.GRIDSAM_ACTIVE:
                        jobStatus.setForeground(new Color(207, 165, 92));
                        jobStatus.setText("AHE ACTIVE");
                        break;
                       
                        case AHEJobObject.GRIDSAM_EXECUTED:
                        jobStatus.setForeground(new Color(127, 124, 133));
                        jobStatus.setText("AHE EXECUTED");
                        break;
                       
                        case AHEJobObject.GRIDSAM_FAILED:
                        jobStatus.setForeground(Color.RED);
                        jobStatus.setText("AHE FAILED");
                        if(pollingButton.getText().equals("Stop Polling")) {
                                pollTimer.stop();
                                pollingButton.setText("Start Polling");
                        }
                        break;
                       
                        case AHEJobObject.GRIDSAM_DONE:
                        jobStatus.setForeground(Color.GREEN);
                        jobStatus.setText("AHE DONE");
                        if(pollingButton.getText().equals("Stop Polling")) {
                                pollTimer.stop();
                                pollingButton.setText("Start Polling");
                        }
                       
                        break;
                       
                        case AHEJobObject.GRIDSAM_TERMINATING:
                        jobStatus.setForeground(new Color(242, 52, 154));
                        jobStatus.setText("AHE TERMINATING");
                        break;
                                       
                        case AHEJobObject.GRIDSAM_TERMINATED:
                        jobStatus.setForeground(new Color(76, 67, 20));
                        jobStatus.setText("AHE TERMINATED");
                        if(pollingButton.getText().equals("Stop Polling")) {
                                pollTimer.stop();
                                pollingButton.setText("Start Polling");
                        }
                        break;
                       
                        case AHEJobObject.GRIDSAM_UNDEFINED:
                        jobStatus.setForeground(new Color(89, 9, 54));
                        jobStatus.setText("AHE UNDEFINED");
                        if(pollingButton.getText().equals("Stop Polling")) {
                                pollTimer.stop();
                                pollingButton.setText("Start Polling");
                        }
                        break;
                       
                       
                }
               
                if (AHEJobObject.GRIDSAM_PENDING <= ajo.getState() && ajo.getState() <= AHEJobObject.GRIDSAM_EXECUTED) {
                        teminateJob.setEnabled(true);
                } else {
                        teminateJob.setEnabled(false);
                }
               
                if (AHEJobObject.GRIDSAM_FAILED <= ajo.getState() && ajo.getState() <= AHEJobObject.GRIDSAM_UNDEFINED) {
                        updateStatus.setEnabled(false);
                        pollingButton.setText("Start Polling");
                        pollingButton.setEnabled(false);
                } else {
                        updateStatus.setEnabled(true);
                }
               
                if (ajo.getState() == AHEJobObject.AHE_JOB_BUILT) {
                        updateStatus.setText("Start Job");
                } else {
                        updateStatus.setText("Update Job");
                }
               
                if (ajo.getState() == AHEJobObject.GRIDSAM_DONE || ajo.getState() == AHEJobObject.GRIDSAM_STAGED_OUT) {
                        downloadButton.setEnabled(true);
                } else {
                        downloadButton.setEnabled(false);
                }
               
               
        }
       
        private void updateGridSAMPanel () {
                Vector vec = ajo.getGridSAMOutput();
               
                gridsamStatusResults.setText("-------------AHE Job Status-------------\n");
                                                       
                if (vec != null) {
                       
                        Iterator it = vec.iterator ();
                        while (it.hasNext ()) {
                                GridSAMStateInfo gssi = (GridSAMStateInfo)it.next();
                                gridsamStatusResults.append("\n" + gssi.toString());
                                gridsamStatusResults.append("\n----------------------------------------\n");
                        }
                }
        }

        private void pollJobState () {
                MonitorSimCall gpc = new MonitorSimCall (ajo.getEndPoint());
                ajo.setGridSAMOutput(gpc.makeCall());
                updateState();
                updateGridSAMPanel ();       
        }
       
       
        private class MyTableModel extends AbstractTableModel {
               
               
           private String[] columnNames;
           private Object[][] data;
               
           public MyTableModel (Object[][] data, String[] columnNames) {
                   this.columnNames = columnNames;
                   this.data = data;
           }
               
            public int getColumnCount() {
                return columnNames.length;
            }

            public int getRowCount() {
                return data.length;
            }

            public String getColumnName(int col) {
                return columnNames[col];
            }

            public Object getValueAt(int row, int col) {
                return data[row][col];
            }

            public Class getColumnClass(int c) {
                return getValueAt(0, c).getClass();
            }

            /*
             * Don't need to implement this method unless your table's
             * editable.
             */
            public boolean isCellEditable(int row, int col) {
                //Note that the data/cell address is constant,
                //no matter where the cell appears onscreen.
                if (col > 0) {
                    return false;
                } else {
                    return true;
                }
            }

            /*
             * Don't need to implement this method unless your table's
             * data can change.
             */
            public void setValueAt(Object value, int row, int col) {
                data[row][col] = value;
                fireTableCellUpdated(row, col);
            }
          
        }

        public boolean getDeleteFiles () {
                return deleteFiles.isSelected();
        }
       
        public AHEJobObject getJobObject() {
                return ajo;
        }

        public Rectangle getCloseBut() {
                return closeBut;
        }

        public void setCloseBut(Rectangle closeBut) {
                this.closeBut = closeBut;
        }
 /*
        public JobRegistryElement getJre() {
                return jre;
        }

        public void setJre(JobRegistryElement jre) {
                this.jre = jre;
        }
       */
       
        
        private class PollThread extends Thread {
    		private boolean run;

    		public PollThread () {
    			run=true;

    		}


    		public void stopThread () {
    			//private boolean shouldirun = true;
    			run=false;

    		} 	

    		public void run() {
    			while(run) {
    				String service = ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.rendezvous");
    				System.err.println("Polling: " + service +" for " + rID);
    				HttpGetPoll hgp = new HttpGetPoll(service, rID);

    				
    				if (hgp.pollService(5)) {
    				p=hgp.getPort();
    				h= hgp.getHost();
    				
    				steer.setEnabled(true);
    				steerERP.setText(h+":"+p);
    				System.err.println("Found sim listening at " + h + ":" + p + "\n");

    				} else {
    					System.err.println("Couldn't rendezvous with HemeLB simulation");
    				}
    				run = false;
    			}


    		}


    	}
       
}
