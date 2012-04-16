package uk.ac.ucl.chem.ccs.clinicalgui;

import java.io.*;
import java.util.Vector;

import org.globus.ftp.*;
import org.globus.gsi.*;
import org.globus.gsi.gssapi.*;
import org.ietf.jgss.GSSCredential;

/**
 * @author Konstantin Voevodski
 *
 * implements methods for getting/putting files and managing directories on the grid ftp server
 */
public class MyGridFtp {
        private GSSCredential cred = null;
        private GridFTPClient client = null;
        private GlobusCredential globusCred = null;
  
        
        public MyGridFtp(String host, int port) {
        	try {
               globusCred = new GlobusCredential(ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.proxyfile"));
               cred = new GlobusGSSCredentialImpl(globusCred,GSSCredential.DEFAULT_LIFETIME);
               client = new GridFTPClient(host, port);
               client.setClientWaitParams(100000, 1000);
               client.authenticate(cred);
               client.setPassiveMode(true);

        	} catch (Exception e) {
        	//	System.err.println("got here");
        		e.printStackTrace();
        	}
        }
        
        public void upload(String localFile, String remoteFile) throws Exception{
        		client.put(new File(localFile), remoteFile, false); 
        }
        
        public void download(String localFile, String remoteFile) throws Exception{
    		client.get(remoteFile, new File(localFile));
        }
        
        public void makeDir(String path) throws Exception{
        	client.makeDir(path);
        }
        
        public Vector list() throws Exception{
        Vector v = new Vector();
        	try {
         v =client.list();
        	} catch (Exception e) {
        		System.err.println("exception on list");
        		e.printStackTrace();
        	}
        	
        	return v;
        	}
        
        public void changeDir(String path) throws Exception{
        	//System.out.println("changing dir to " + path);
        	try{
        	client.changeDir(path);
        	} catch (Exception e) {
        		System.err.println("exception on cd");

        	}
        	}
        
        public void close() throws Exception{
                client.close();
                cred.dispose();
        }

}
