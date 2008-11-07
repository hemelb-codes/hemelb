/*
 * Created on 30-Jun-2008
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package uk.ac.ucl.chem.ccs.clinicalgui;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.DICOM;

import java.awt.Image;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JTree;
import javax.swing.tree.DefaultMutableTreeNode;

import org.dcm4che2.data.DicomObject;
import org.dcm4che2.data.Tag;


/**
 * @author Konstantin Voevodski
 *
 * implements methods for querying and obtaining data from the DICOM server
 */
public class DicomServerInterface {
	/*
	private static String SERVER_LOC = "GENIUS@bunsen.chem.ucl.ac.uk:11112";
	private static String CLIENT_LOC = "DCMRCV@128.40.177.220:11113";
	private static String CLIENT_NAME = "DCMRCV";
	private static String TEMP_DIRECTORY = "/tmp/dcmQueryOutput";
	public static String SLICES_FILE_PATH = "/tmp/SlicesFile.txt";
	*/
	
	private static String SERVER_LOC = ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.dicomservername") + "@" +
		ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.dicomserverhost") + ":" + 
		ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.dicomserverport");
	
	private static String CLIENT_LOC = ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.dicomclientname") + "@" +
		ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.dicomclienthost") + ":" + 
		ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.dicomclientport");
	
	private static String CLIENT_NAME = ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.dicomclientname");
	private static String TEMP_DIRECTORY = ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.dicomtmpdir") + "/dcmQueryOutput";
	public static String SLICES_FILE_PATH = ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.dicomtmpdir") + "/SlicesFile.txt";
	
	private static File tempDirectory = new File(TEMP_DIRECTORY);
	private static Pattern p = Pattern.compile("\\[.*\\]");
	private static DefaultMutableTreeNode top = new DefaultMutableTreeNode("Available Data");
	
	public static int queryReceiveWrite(String patientId, String studyId, String seriesId) throws IOException{
	//queries, receives and writes DICOM files corresponding to a particular dataset
    //returns number of DICOM files received
		tempDirectory.mkdir();
		String[] a = new String[3];
		a[0] = CLIENT_LOC;
		a[1] = "-dest";
		a[2] = TEMP_DIRECTORY;
		DcmRcv receiver = new DcmRcv();
		receiver.run(a);
		a = new String[7];
		a[0] = SERVER_LOC;
		//a[1] = "-qPatientID=" + patientId;
		//a[2] = "-qStudyID=" + studyId;
		//a[3] = "-q0020000E=" + seriesId;
		//a[4] = "-rSliceThickness";
		//a[5] = "-rPixelSpacing";
		//a[4] = "-S";
		a[1] = "-cmove";
		a[2] = CLIENT_NAME;
		a[3] = "-qPatientID=" + patientId;
		a[4] = "-qStudyID=" + studyId;
		a[5] = "-q0020000E=" + seriesId;
		a[6] = "-I";
		//double[] imageData = {-1,-1};
		List<DicomObject> results = DcmQR.run(a);
		/*
		if(result.size() == 0)
			return imageData;
		DicomObject dcmObj = result.get(0);
		//we assume that slice thickness, pixel spacing is the same for all the fetched files
		String s1 = dcmObj.get(Tag.SliceThickness).toString();
		String s2 = dcmObj.get(Tag.PixelSpacing).toString();
		System.out.println("slice thickness attribute: " + s1);
		System.out.println("pixel spacing attribute " + s2);
		imageData[0] = Double.parseDouble(parseBrackets(s1));
		imageData[1] = Double.parseDouble(parseBrackets(s2));
		*/
		receiver.stop();
		if(results.size() == 0)
			return 0;
		writeSlicesFile();
		String[] filenames = tempDirectory.list();
		for(String filename: filenames){
			File tempFile = new File(TEMP_DIRECTORY + "/" + filename);
			tempFile.delete();
		}
		tempDirectory.delete();
		return results.size();
		//return imageData;
	}
	
	public static String[] getPatientAttributes(String patientId){
	//returns Patient Sex, Patient DOB of the data set given by (patientId)
		
		String[] a = new String[4];
		a[0] = SERVER_LOC;
		a[1] = "-qPatientID=" + patientId;
		a[2] = "-r00100040";
		a[3] = "-r00100030";
		List<DicomObject> results = DcmQR.run(a);
		DicomObject dcmObj = results.get(0);
		String[] output = new String[2];
		output[0] = parseBrackets(dcmObj.get(Tag.PatientSex).toString());
		output[1] = parseBrackets(dcmObj.get(Tag.PatientBirthDate).toString());
		return output;
	}
	
	public static String[] getStudyAttributes(String patientId, String studyId){
	//returns Study Description, Referring Physician of the data set given by (patientId, studyId)
		
		String[] a = new String[5];
		a[0] = SERVER_LOC;
		a[1] = "-qPatientID=" + patientId;
		a[2] = "-qStudyID=" + studyId;
		a[3] = "-r00081030";
		a[4] = "-r00080090";
		List<DicomObject> results = DcmQR.run(a);
		DicomObject dcmObj = results.get(0);
		String[] output = new String[2];
		//System.out.println("the size of results is " + results.size());
		output[0] = parseBrackets(dcmObj.get(Tag.StudyDescription).toString());
		output[1] = parseBrackets(dcmObj.get(Tag.ReferringPhysicianName).toString());
		return output;
	}
	
	public static String[] getSeriesAttributes(String patientId, String studyId, String seriesId){
	//returns Date, Time, Series Description, Institution Name, and Performing Physician name of the data set given by (patientId, studyId, seriesId)
		
		String[] a = new String[10];
		a[0] = SERVER_LOC;
		a[1] = "-qPatientID=" + patientId;
		a[2] = "-qStudyID=" + studyId;
		a[3] = "-q0020000E=" + seriesId;
		a[4] = "-I";
		a[5] = "-r00080021";
		a[6] = "-r00080031";
		a[7] = "-r0008103E";
		a[8] = "-r00080080";
		a[9] = "-r00081050";
		List<DicomObject> results = DcmQR.run(a);
		DicomObject dcmObj = results.get(0);
		String[] output = new String[5];
		output[0] = parseBrackets(dcmObj.get(Tag.SeriesDate).toString());
		output[1] = parseBrackets(dcmObj.get(Tag.SeriesTime).toString());
		output[2] = parseBrackets(dcmObj.get(Tag.SeriesDescription).toString());
		output[3] = parseBrackets(dcmObj.get(Tag.InstitutionName).toString());
		output[4] = parseBrackets(dcmObj.get(Tag.PerformingPhysicianName).toString());
		return output;
	}
	
	public static JTree populateServerData() {
	//returns a tree object representing all data available on DICOM server
	//tree organized by (patient,study,series)
		String[] a = new String[5];
		a[0] = SERVER_LOC;
		a[1] = "-rPatientID";
		a[2] = "-rStudyID";
		a[3] = "-r0020000E";
		a[4] = "-S";
		List<DicomObject> result = DcmQR.run(a);
		Hashtable<String,HashSet<String>> Patients = new Hashtable<String,HashSet<String>>();
		//map each patient id to a set of study ids
		Hashtable<String,HashSet<String>> Studies = new Hashtable<String,HashSet<String>>();
		//map each study id to a set of series ids
		System.out.println("there are " + result.size() + " objects returned");
		String patient, study, series;
        for(DicomObject dcmObj : result){
        	patient = parseBrackets(dcmObj.get(Tag.PatientID).toString());
        	study = parseBrackets(dcmObj.get(Tag.StudyID).toString());
            series = parseBrackets(dcmObj.get(Tag.SeriesInstanceUID).toString());
            System.out.println(patient);
            System.out.println(study);
            System.out.println(series);
            augmentTable(Patients,patient,study);
            augmentTable(Studies,study,series);
        }
        Enumeration myEnum = Patients.keys();
        while(myEnum.hasMoreElements()){
        	patient = (String) myEnum.nextElement();
        	DefaultMutableTreeNode patientNode = new DefaultMutableTreeNode(patient);
        	top.add(patientNode);
        	HashSet<String> studiesSet = (HashSet<String>) Patients.get(patient);
        	for(String s: studiesSet){
        		DefaultMutableTreeNode studyNode = new DefaultMutableTreeNode(s);
        		patientNode.add(studyNode);
        		HashSet<String> seriesSet = (HashSet<String>) Studies.get(s);
        		for(String s2: seriesSet){
        			DefaultMutableTreeNode seriesNode = new DefaultMutableTreeNode(s2);
        			studyNode.add(seriesNode);
        		}
        		
        	}
        }
        return new JTree(top);
	}
	
	private static String parseBrackets(String s){
	//parses out outermost brackets in s, returning what is inside: given [text_here] returns text_here
		 Matcher m = p.matcher(s);
		 m.find();
		 String temp = m.group();
		 return(temp.substring(1,temp.length()-1));
	}
	
	private static void augmentTable(Hashtable table, String a, String b){
	//adds b to a's hashset: used for keeping track of the DICOM data hierarchy, for example,
	//keeping track of all studies associated with a particular patient
		HashSet s;
		if(table.containsKey(a))
			s = (HashSet) table.get(a);
		else
			s = new HashSet();
		s.add(b);
		table.put(a,s);
			
	}

	
	private static void writeSlicesFile() throws IOException{
	//subroutine of queryReceiveWrite: reads slice files received from server and creates one file
	//containing the pixel data of all slices, which is used as input to seg tool
		String[] filenames = tempDirectory.list();
		Arrays.sort(filenames);
		BufferedWriter writer = new BufferedWriter(new FileWriter(SLICES_FILE_PATH));
		if(filenames.length == 0)
			return;
		DICOM d = new DICOM(new FileInputStream(TEMP_DIRECTORY + "/" + filenames[0]));
		d.run("Name");
		int width = d.getWidth();
		int height = d.getHeight();
		writer.write(width + "\t" + height);
		writer.newLine();
		ImageStack s = new ImageStack(width,height);
		int counter = 1;
		for(String filename : filenames){
			System.out.println(filename);
			d = new DICOM(new FileInputStream(TEMP_DIRECTORY + "/" + filename));
			d.run("Name");
			Image i = d.getImage();
			for(int x = 1; x <= height; x++){
				for(int y = 1; y <= width; y++){
					int[] pixels = d.getPixel(x,y);
					writer.write(Integer.toString(pixels[0]));
					if(x != height || y != width)
						writer.write("\t");
				}
			}
			s.addSlice("slice" + Integer.toString(counter++), d.getProcessor());
			writer.newLine();
		}
		writer.close();
		ModelGenerationPanel.dicomImageStack = s;
	}

}
