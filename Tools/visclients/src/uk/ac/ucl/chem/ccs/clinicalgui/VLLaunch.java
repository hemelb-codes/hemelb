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
 * Created on 2:30:45 PM Mar 20, 2006
 * Project: AHE-GUI
 * File: AHEPrepare.java
 *
 * @author stefan.zasada@ucl.ac.uk
 *
 * TODO:
 */
package uk.ac.ucl.chem.ccs.clinicalgui;

import java.io.File;
import java.util.Iterator;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.log4j.PropertyConfigurator;
import java.util.Vector;
import uk.ac.ucl.chem.ccs.aheclient.res.AdvancedReservation;
import uk.ac.ucl.chem.ccs.aheclient.wsrf.PrepareCall;
import uk.ac.ucl.chem.ccs.aheclient.wsrf.StartCall; //import uk.ac.ucl.chem.ccs.aheclient.confparser.AHEConfParser;
//import uk.ac.ucl.chem.ccs.aheclient.confparser.DefaultConfParser;
import uk.ac.ucl.chem.ccs.aheclient.util.AHEJobObject;
import uk.ac.ucl.chem.ccs.aheclient.util.GridSAMStateInfo;
import uk.ac.ucl.chem.ccs.aheclient.util.JobFileElement;
import uk.ac.ucl.chem.ccs.aheclient.util.JobRegistryElement;
import uk.ac.ucl.chem.ccs.aheclient.util.PrepareResponse;
import uk.ac.ucl.chem.ccs.aheclient.util.StageFilesIn;
import uk.ac.ucl.chem.ccs.aheclient.util.StageFilesOut;
import uk.ac.ucl.chem.ccs.aheclient.util.Tools;

/**
 * Prepare, start, monitor job and download output.
 */
public class VLLaunch
{
  private static Log         cat;
  private PrepareCall        pc = null;
  private PrepareResponse    pr = null;
  private JobRegistryElement jre;
  private String             rtParsPath;
  private String             parsPath;
  private String             configPath;

  // private AHEJobObject ajo = null;

  /**
   * Initialize
   */
  public VLLaunch(String s1, String s2, String s3)
  {
    System.out.println(s1);
    System.out.println(s2);
    System.out.println(s3);
    rtParsPath = s1;
    parsPath = s2;
    configPath = s3;
  }

  /**
   * Method to make the PrepareCall, once the constructor has initialised the properties and
   * keystore etc.
   * 
   * @param rmCPUCount
   * @param simulationName
   * @param app
   * @return PrepareResponse
   */
  public PrepareResponse prepare(String simulationName)
  {
    String app = ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.appname");
    int rmCPUCount = 1;
    // Note, most of the parameters to this call are used for resource matching and aren't needed,
    // hence the empty strings
    // Only call once
    if (pr == null)
    {
      pc = new PrepareCall(ClinicalGuiClient.prop
          .getProperty("uk.ac.ucl.chem.ccs.aheclient.jobregepr"), rmCPUCount, simulationName,
          "single", app, "", 12 * 60, "", "", "", "", "", "", "");

      pr = pc.makeCall();
      String appWSRes = pr.getAppWSResource();
      String appServiceGroupEntry = pr.getAppServiceGroupEntry();

      jre = new JobRegistryElement(appWSRes, appServiceGroupEntry, app, "", "", "", "");
    }
    return pr;
  }

  /**
   * Method to make the StartCall, once the constructor has initialised the properties and keystore
   * etc. Prepare must be performed first
   * 
   * @param app
   * @param config
   * @param rm
   * @param rmCPUCount
   * @return
   */
  public AHEJobObject start(String rm, int rmCPUCount, AdvancedReservation ad)
  {
    String app = "hemelb";

    AHEJobObject ajo = null;

    if (pr != null)
    {
      // The Start call needs lots of the properties from the config file
      String stagePath = ClinicalGuiClient.prop
          .getProperty("uk.ac.ucl.chem.ccs.aheclient.ahedavserver");
      String davun = ClinicalGuiClient.prop.getProperty("uk.ac.ucl.chem.ccs.aheclient.ahedavuser");
      String davpw = ClinicalGuiClient.prop
          .getProperty("uk.ac.ucl.chem.ccs.aheclient.ahedavpasswd");
      String myproxyserver = ClinicalGuiClient.prop
          .getProperty("uk.ac.ucl.chem.ccs.aheclient.myproxy-server");
      String myproxydn = ClinicalGuiClient.prop
          .getProperty("uk.ac.ucl.chem.ccs.aheclient.myproxy-dn");
      String myproxyport = ClinicalGuiClient.prop
          .getProperty("uk.ac.ucl.chem.ccs.aheclient.myproxy-port");
      String myproxylifetime = ClinicalGuiClient.prop
          .getProperty("uk.ac.ucl.chem.ccs.aheclient.myproxy-lifetime");
      String myproxyun = ClinicalGuiClient.prop
          .getProperty("uk.ac.ucl.chem.ccs.aheclient.myproxy-un");
      String myproxypw = ClinicalGuiClient.prop
          .getProperty("uk.ac.ucl.chem.ccs.aheclient.myproxy-pw");

      // vectors of input and output files
      Vector inFiles = new Vector();
      Vector outFiles = new Vector();

      // get resource id
      String resourceID = jre.getResourceID();

      // set input files
      inFiles.add(new JobFileElement("par1", "rtPars", "foo", rtParsPath));
      inFiles.add(new JobFileElement("par2", "pars", "foo", parsPath));
      inFiles.add(new JobFileElement("par3", "config", "foo", configPath));

      // set up the rendezvous id
      String args = resourceID.substring(0, 4);
      System.err.print("resource id is: " + args);
      if (args.startsWith("0"))
      {
        args = "1" + args;
      }

      inFiles.add(new JobFileElement("argument", args, null, null));
      outFiles.add(new JobFileElement("file1", "file1", "file1", "file1"));

      // set up AJO with job details
      ajo = new AHEJobObject(app, resourceID, jre.getMemberServiceEPR(), inFiles, outFiles, rm,
          "pars.asc", rmCPUCount);
      ajo.setAppServiceGroup(jre.getServiceGroupEntryEPR());
      ajo.setSimName(jre.getComponentTaskDescription());
      ajo.setAdvancedReservation(ad);
      
      // make call
      StartCall sc = new StartCall(ajo, myproxylifetime, myproxyport, myproxydn, myproxyserver,
          myproxypw, myproxyun);

      ajo = sc.makeCall();

    }
    return ajo;
  }

}
