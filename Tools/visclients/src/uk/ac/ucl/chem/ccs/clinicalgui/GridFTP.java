package uk.ac.ucl.chem.ccs.clinicalgui;
import java.io.File;
import java.util.Vector;

import org.globus.ftp.FileInfo;
import org.globus.ftp.GridFTPClient;
import org.globus.ftp.GridFTPSession;

/**
* testing GridFTP class
* 
* 
*/

public class GridFTP{
	
public static void main(String[] args) throws Exception{
	 	String host = "bunsen.chem.ucl.ac.uk";
        String remoteFile = "/home/konstantin/dir/pars.asc";
        String localFile = "/tmp/downloadTest.asc";
        MyGridFtp o = new MyGridFtp(host, 2811);
        /*
        o.changeDir("/home/konstantin/dir");
        Vector<FileInfo> v = o.list();
        for(FileInfo f : v){
        	System.out.println(f.getName());
        	System.out.println(f.isFile());
        	System.out.println(f.isDirectory());
        }
        */
        o.download(localFile,remoteFile);
        o.close();
}
}
