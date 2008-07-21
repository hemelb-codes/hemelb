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
		DICOM d = new DICOM(new FileInputStream("/tmp/dcmQueryOutput/1.3.12.2.1107.5.99.2.11737.30000008032011100623400000250"));
		d.run("Name");
		Image i = d.getImage();
        //d.show();
        int width = d.getWidth();
		int height = d.getHeight();
		BufferedWriter writer = new BufferedWriter(new FileWriter("/tmp/dcmQueryOutput/pixels.txt"));
		writer.write(width + "\t" + height + "\t");
		//int[] pixels;
		//pixels = d.getPixel(1,1);
		//System.out.println(pixels[0] + " " + pixels[1] + " " + pixels[2] + " " + pixels[3]);
		//System.out.println(pixels.length);
		for(int x = 1; x <= height; x++){
			for(int y = 1; y <= width; y++){
				int[] pixels = d.getPixel(x,y);
				//System.out.println(pixels.length);
				if(pixels[0] != 0)
					System.out.println("non-zero pixel");
				writer.write(Integer.toString(pixels[0]));
				if(x != height || y != width)
					writer.write("\t");
			}
		}
		writer.newLine();
		writer.close();

	}

}
