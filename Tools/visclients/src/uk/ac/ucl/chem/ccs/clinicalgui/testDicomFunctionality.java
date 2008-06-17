package uk.ac.ucl.chem.ccs.clinicalgui;

import ij.ImagePlus;
import ij.plugin.DICOM;

import java.awt.Image;
import java.io.*;

public class testDicomFunctionality {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		DICOM d = new DICOM(new FileInputStream("C:/users/konstantin/desktop/dicomFiles/file4.dcm"));
		d.run("Name");
        d.show();

	}

}
