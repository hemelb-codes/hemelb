package uk.ac.ucl.chem.ccs.clinicalgui;

import java.io.*;
import java.util.Vector;

import org.globus.ftp.*;
import org.globus.gsi.*;
import org.globus.gsi.gssapi.*;
import org.ietf.jgss.GSSCredential;

/**
 * Sample Java class to perform a transfer from a remote GRIDFTP host
 */
public class MyGridFtp {
        private GSSCredential cred = null;
        private GridFTPClient client = null;
        private GlobusCredential globusCred = null;
  
        
        public MyGridFtp(String host, int port) throws Exception {
               globusCred = new GlobusCredential("/home/konstantin/x509up_u8012");
               cred = new GlobusGSSCredentialImpl(globusCred,GSSCredential.DEFAULT_LIFETIME);
               client = new GridFTPClient(host , port);
               client.authenticate(cred);
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
        	return client.list();
        }
        
        public void changeDir(String path) throws Exception{
        	client.changeDir(path);
        }
        
        public void close() throws Exception{
                client.close();
                cred.dispose();
        }

}
