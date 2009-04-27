package uk.ac.ucl.chem.ccs.vizclient;

import java.net.URL;
import java.io.IOException;
import java.lang.Integer;
import org.apache.commons.httpclient.HttpClient;
import org.apache.commons.httpclient.HttpException;
import org.apache.commons.httpclient.HttpMethod;
import org.apache.commons.httpclient.methods.GetMethod;

public class HttpGetPoll {
 
	private URL target;

	private int port = 0;
	private String host = null;
	
	public HttpGetPoll (String serviceAddress, String jobID) {
	
		try {
			if (!serviceAddress.endsWith("/")) {
				serviceAddress = serviceAddress + "/";
			}
			target = new URL (serviceAddress + jobID);
		} catch (Exception e) {
			target = null;
		}
	}
	
	public HttpGetPoll (String steeringServiceEPR) {
		
		try {
			target = new URL (steeringServiceEPR);
		} catch (Exception e) {
			target = null;
		}
	}

	public boolean pollService (int freq) {
		return pollService (freq, 0);
	}

	public boolean pollService (int freq, int timeout) {
		
		if (target == null) {
			return false;
		}
		
		int exectTime = 0;
		boolean httpOk = false;
        String responseBody = null;

    	
        //create a singular HttpClient object
        HttpClient client = new HttpClient();

        //establish a connection within 5 seconds
        client.getHttpConnectionManager().getParams().setConnectionTimeout(5000);
        
		while (!httpOk) {

	        HttpMethod method = new GetMethod(target.toExternalForm());
	        method.setFollowRedirects(true);

	        //execute the method
	        try{
	            if (client.executeMethod(method) == 200) {
	//            if (method.getStatusCode() == 200) {
	            	httpOk = true;
		            responseBody = method.getResponseBodyAsString();
	            }
	            	            
	        } catch (HttpException he) {
	            System.err.println("Http error connecting to '" + target + "'");
	            System.err.println(he.getMessage());
	            return false;
	        } catch (IOException ioe){
	            System.err.println("Unable to connect to '" + target + "'");
	            return false;
	        }
	        //clean up the connection resources
	        method.releaseConnection();

	        //pause between retry
	        try {
	            Thread.currentThread().sleep(1000*freq);
	        }
	          catch (InterruptedException e) {
	        }
	       
	          if (timeout != 0) {
	        	  exectTime = exectTime + freq;
	          }
	        
	        if (exectTime > timeout) {
	        	return false;
	        }
	        
		}
		
		if (responseBody != null) {
			String ans[] = responseBody.split(":");
			host = ans[0];
			port = Integer.parseInt(ans[1]);
		}		
		
		return true;
	}

	public int getPort() {
		return port;
	}

	public String getHost() {
		return host;
	}
	
	public static void main(String[] args){
		HttpGetPoll hgp = new HttpGetPoll("http://bunsen.chem.ucl.ac.uk:28080/ahe/test/rendezvous", "123456578");
		hgp.pollService(5);
		
		System.out.println("Host = " + hgp.getHost() + "\nport = "+ hgp.getPort());
	}
	
	
}
