package uk.ac.ucl.chem.ccs.clinicalgui;

/* ***** BEGIN LICENSE BLOCK *****
 * Version: MPL 1.1/GPL 2.0/LGPL 2.1
 *
 * The contents of this file are subject to the Mozilla Public License Version
 * 1.1 (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
 * for the specific language governing rights and limitations under the
 * License.
 *
 * The Original Code is part of dcm4che, an implementation of DICOM(TM) in
 * Java(TM), hosted at http://sourceforge.net/projects/dcm4che.
 *
 * The Initial Developer of the Original Code is
 * Gunter Zeilinger, Huetteldorferstr. 24/10, 1150 Vienna/Austria/Europe.
 * Portions created by the Initial Developer are Copyright (C) 2002-2005
 * the Initial Developer. All Rights Reserved.
 *
 * Contributor(s):
 * Gunter Zeilinger <gunterze@gmail.com>
 *
 * Alternatively, the contents of this file may be used under the terms of
 * either the GNU General Public License Version 2 or later (the "GPL"), or
 * the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
 * in which case the provisions of the GPL or the LGPL are applicable instead
 * of those above. If you wish to allow use of your version of this file only
 * under the terms of either the GPL or the LGPL, and not to allow others to
 * use your version of this file under the terms of the MPL, indicate your
 * decision by deleting the provisions above and replace them with the notice
 * and other provisions required by the GPL or the LGPL. If you do not delete
 * the provisions above, a recipient may use your version of this file under
 * the terms of any one of the MPL, the GPL or the LGPL.
 *
 * ***** END LICENSE BLOCK ***** */

//package org.dcm4che2.tool.dcmrcv;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.security.GeneralSecurityException;
import java.security.KeyStore;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.dcm4che2.data.BasicDicomObject;
import org.dcm4che2.data.DicomObject;
import org.dcm4che2.data.Tag;
import org.dcm4che2.data.UID;
import org.dcm4che2.io.DicomOutputStream;
import org.dcm4che2.net.Association;
import org.dcm4che2.net.CommandUtils;
import org.dcm4che2.net.Device;
import org.dcm4che2.net.DicomServiceException;
import org.dcm4che2.net.Executor;
import org.dcm4che2.net.NetworkApplicationEntity;
import org.dcm4che2.net.NetworkConnection;
import org.dcm4che2.net.NewThreadExecutor;
import org.dcm4che2.net.PDVInputStream;
import org.dcm4che2.net.Status;
import org.dcm4che2.net.TransferCapability;
import org.dcm4che2.net.service.StorageService;
import org.dcm4che2.net.service.VerificationService;

/**
 * @author gunter zeilinger(gunterze@gmail.com)
 * @version $Revision$ $Date$
 * @since Oct 13, 2005
 */
public class DcmRcv extends StorageService {

    private static final int KB = 1024;

    private static final String USAGE = "dcmrcv [Options] [<aet>[@<ip>]:]<port>";

    private static final String DESCRIPTION = "DICOM Server listening on specified <port> for incoming association "
            + "requests. If no local IP address of the network interface is specified "
            + "connections on any/all local addresses are accepted. If <aet> is "
            + "specified, only requests with matching called AE title will be "
            + "accepted.\n" + "Options:";

    private static final String EXAMPLE = "\nExample: dcmrcv DCMRCV:11112 -dest /tmp \n"
            + "=> Starts server listening on port 11112, accepting association "
            + "requests with DCMRCV as called AE title. Received objects "
            + "are stored to /tmp.";

    private static char[] SECRET = { 's', 'e', 'c', 'r', 'e', 't' };
    
    private static final String[] ONLY_DEF_TS = { UID.ImplicitVRLittleEndian };

    private static final String[] NATIVE_TS = { UID.ExplicitVRLittleEndian,
            UID.ExplicitVRBigEndian, UID.ImplicitVRLittleEndian };

    private static final String[] NATIVE_LE_TS = { UID.ExplicitVRLittleEndian,
            UID.ImplicitVRLittleEndian };

    private static final String[] NON_RETIRED_TS = { UID.JPEGLSLossless,
            UID.JPEGLossless, UID.JPEGLosslessNonHierarchical14,
            UID.JPEG2000LosslessOnly, UID.DeflatedExplicitVRLittleEndian,
            UID.RLELossless, UID.ExplicitVRLittleEndian,
            UID.ExplicitVRBigEndian, UID.ImplicitVRLittleEndian,
            UID.JPEGBaseline1, UID.JPEGExtended24, UID.JPEGLSLossyNearLossless,
            UID.JPEG2000, UID.MPEG2, };

    private static final String[] NON_RETIRED_LE_TS = { UID.JPEGLSLossless,
            UID.JPEGLossless, UID.JPEGLosslessNonHierarchical14,
            UID.JPEG2000LosslessOnly, UID.DeflatedExplicitVRLittleEndian,
            UID.RLELossless, UID.ExplicitVRLittleEndian,
            UID.ImplicitVRLittleEndian, UID.JPEGBaseline1, UID.JPEGExtended24,
            UID.JPEGLSLossyNearLossless, UID.JPEG2000, UID.MPEG2, };

    private static final String[] CUIDS = {
            UID.BasicStudyContentNotificationSOPClassRetired,
            UID.StoredPrintStorageSOPClassRetired,
            UID.HardcopyGrayscaleImageStorageSOPClassRetired,
            UID.HardcopyColorImageStorageSOPClassRetired,
            UID.ComputedRadiographyImageStorage,
            UID.DigitalXRayImageStorageForPresentation,
            UID.DigitalXRayImageStorageForProcessing,
            UID.DigitalMammographyXRayImageStorageForPresentation,
            UID.DigitalMammographyXRayImageStorageForProcessing,
            UID.DigitalIntraoralXRayImageStorageForPresentation,
            UID.DigitalIntraoralXRayImageStorageForProcessing,
            UID.StandaloneModalityLUTStorageRetired,
            UID.EncapsulatedPDFStorage, UID.StandaloneVOILUTStorageRetired,
            UID.GrayscaleSoftcopyPresentationStateStorageSOPClass,
            UID.ColorSoftcopyPresentationStateStorageSOPClass,
            UID.PseudoColorSoftcopyPresentationStateStorageSOPClass,
            UID.BlendingSoftcopyPresentationStateStorageSOPClass,
            UID.XRayAngiographicImageStorage, UID.EnhancedXAImageStorage,
            UID.XRayRadiofluoroscopicImageStorage, UID.EnhancedXRFImageStorage,
            UID.XRayAngiographicBiPlaneImageStorageRetired,
            UID.PositronEmissionTomographyImageStorage,
            UID.StandalonePETCurveStorageRetired, UID.CTImageStorage,
            UID.EnhancedCTImageStorage, UID.NuclearMedicineImageStorage,
            UID.UltrasoundMultiframeImageStorageRetired,
            UID.UltrasoundMultiframeImageStorage, UID.MRImageStorage,
            UID.EnhancedMRImageStorage, UID.MRSpectroscopyStorage,
            UID.RTImageStorage, UID.RTDoseStorage, UID.RTStructureSetStorage,
            UID.RTBeamsTreatmentRecordStorage, UID.RTPlanStorage,
            UID.RTBrachyTreatmentRecordStorage,
            UID.RTTreatmentSummaryRecordStorage,
            UID.NuclearMedicineImageStorageRetired,
            UID.UltrasoundImageStorageRetired, UID.UltrasoundImageStorage,
            UID.RawDataStorage, UID.SpatialRegistrationStorage,
            UID.SpatialFiducialsStorage, UID.RealWorldValueMappingStorage,
            UID.SecondaryCaptureImageStorage,
            UID.MultiframeSingleBitSecondaryCaptureImageStorage,
            UID.MultiframeGrayscaleByteSecondaryCaptureImageStorage,
            UID.MultiframeGrayscaleWordSecondaryCaptureImageStorage,
            UID.MultiframeTrueColorSecondaryCaptureImageStorage,
            UID.VLImageStorageRetired, UID.VLEndoscopicImageStorage,
            UID.VideoEndoscopicImageStorage, UID.VLMicroscopicImageStorage,
            UID.VideoMicroscopicImageStorage,
            UID.VLSlideCoordinatesMicroscopicImageStorage,
            UID.VLPhotographicImageStorage, UID.VideoPhotographicImageStorage,
            UID.OphthalmicPhotography8BitImageStorage,
            UID.OphthalmicPhotography16BitImageStorage,
            UID.StereometricRelationshipStorage,
            UID.VLMultiframeImageStorageRetired,
            UID.StandaloneOverlayStorageRetired, UID.BasicTextSR,
            UID.EnhancedSR, UID.ComprehensiveSR, UID.ProcedureLogStorage,
            UID.MammographyCADSR, UID.KeyObjectSelectionDocument,
            UID.ChestCADSR, UID.StandaloneCurveStorageRetired,
            UID._12leadECGWaveformStorage, UID.GeneralECGWaveformStorage,
            UID.AmbulatoryECGWaveformStorage, UID.HemodynamicWaveformStorage,
            UID.CardiacElectrophysiologyWaveformStorage,
            UID.BasicVoiceAudioWaveformStorage, UID.HangingProtocolStorage,
            UID.SiemensCSANonImageStorage };

    private Executor executor = new NewThreadExecutor("DCMRCV");

    private Device device = new Device("DCMRCV");

    private NetworkApplicationEntity ae = new NetworkApplicationEntity();

    private NetworkConnection nc = new NetworkConnection();

    private String[] tsuids = NON_RETIRED_LE_TS;

    private File destination;

    private boolean devnull;

    private int fileBufferSize = 256;

    private int rspdelay = 0;

    private String keyStoreURL = "resource:tls/test_sys_2.p12";
    
    private char[] keyStorePassword = SECRET; 

    private char[] keyPassword; 
    
    private String trustStoreURL = "resource:tls/mesa_certs.jks";
    
    private char[] trustStorePassword = SECRET; 
    
    public DcmRcv() {
        super(CUIDS);
        device.setNetworkApplicationEntity(ae);
        device.setNetworkConnection(nc);
        ae.setNetworkConnection(nc);
        ae.setAssociationAcceptor(true);
        ae.register(new VerificationService());
        ae.register(this);
    }

    public final void setAEtitle(String aet) {
        ae.setAETitle(aet);
    }

    public final void setHostname(String hostname) {
        nc.setHostname(hostname);
    }

    public final void setPort(int port) {
        nc.setPort(port);
    }

    public final void setTlsWithoutEncyrption() {
        nc.setTlsWithoutEncyrption();
    }

    public final void setTls3DES_EDE_CBC() {
        nc.setTls3DES_EDE_CBC();
    }

    public final void setTlsAES_128_CBC() {
        nc.setTlsAES_128_CBC();
    }
    
    public final void disableSSLv2Hello() {
        nc.disableSSLv2Hello();
    }
    
    public final void setTlsNeedClientAuth(boolean needClientAuth) {
        nc.setTlsNeedClientAuth(needClientAuth);
    }
    
    public final void setKeyStoreURL(String url) {
        keyStoreURL = url;
    }
    
    public final void setKeyStorePassword(String pw) {
        keyStorePassword = pw.toCharArray();
    }
    
    public final void setKeyPassword(String pw) {
        keyPassword = pw.toCharArray();
    }
    
    public final void setTrustStorePassword(String pw) {
        trustStorePassword = pw.toCharArray();
    }
    
    public final void setTrustStoreURL(String url) {
        trustStoreURL = url;
    }
        
    public final void setPackPDV(boolean packPDV) {
        ae.setPackPDV(packPDV);
    }

    public final void setAssociationReaperPeriod(int period) {
        device.setAssociationReaperPeriod(period);
    }

    public final void setTcpNoDelay(boolean tcpNoDelay) {
        nc.setTcpNoDelay(tcpNoDelay);
    }

    public final void setRequestTimeout(int timeout) {
        nc.setRequestTimeout(timeout);
    }

    public final void setReleaseTimeout(int timeout) {
        nc.setReleaseTimeout(timeout);
    }

    public final void setSocketCloseDelay(int delay) {
        nc.setSocketCloseDelay(delay);
    }

    public final void setIdleTimeout(int timeout) {
        ae.setIdleTimeout(timeout);
    }

    public final void setDimseRspTimeout(int timeout) {
        ae.setDimseRspTimeout(timeout);
    }

    public final void setMaxPDULengthSend(int maxLength) {
        ae.setMaxPDULengthSend(maxLength);
    }

    public void setMaxPDULengthReceive(int maxLength) {
        ae.setMaxPDULengthReceive(maxLength);
    }

    public final void setReceiveBufferSize(int bufferSize) {
        nc.setReceiveBufferSize(bufferSize);
    }

    public final void setSendBufferSize(int bufferSize) {
        nc.setSendBufferSize(bufferSize);
    }

    private void setDimseRspDelay(int delay) {
        rspdelay = delay;
    }

    private static CommandLine parse(String[] args) {
        Options opts = new Options();
        
        OptionBuilder.withArgName("NULL|3DES|AES");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "enable TLS connection without, 3DES or AES encryption");
        opts.addOption(OptionBuilder.create("tls"));
        
        opts.addOption("nossl2", false, "disable acceptance of SSLv2Hello TLS handshake");
        opts.addOption("noclientauth", false, "disable client authentification for TLS");        

        OptionBuilder.withArgName("file|url");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "file path or URL of P12 or JKS keystore, resource:tls/test_sys_2.p12 by default");
        opts.addOption(OptionBuilder.create("keystore"));

        OptionBuilder.withArgName("password");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "password for keystore file, 'secret' by default");
        opts.addOption(OptionBuilder.create("keystorepw"));

        OptionBuilder.withArgName("password");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "password for accessing the key in the keystore, keystore password by default");
        opts.addOption(OptionBuilder.create("keypw"));

        OptionBuilder.withArgName("file|url");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "file path or URL of JKS truststore, resource:tls/mesa_certs.jks by default");
        opts.addOption(OptionBuilder.create("truststore"));

        OptionBuilder.withArgName("password");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "password for truststore file, 'secret' by default");
        opts.addOption(OptionBuilder.create("truststorepw"));
        
        OptionBuilder.withArgName("dir");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "store received objects into files in specified directory <dir>."
                        + " Do not store received objects by default.");
        opts.addOption(OptionBuilder.create("dest"));

        opts.addOption("defts", false, "accept only default transfer syntax.");
        opts.addOption("bigendian", false,
                "accept also Explict VR Big Endian transfer syntax.");
        opts.addOption("native", false,
                "accept only transfer syntax with uncompressed pixel data.");

        OptionBuilder.withArgName("maxops");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "maximum number of outstanding operations performed "
                        + "asynchronously, unlimited by default.");
        opts.addOption(OptionBuilder.create("async"));

        opts.addOption("pdv1", false, "send only one PDV in one P-Data-TF PDU, "
                + "pack command and data PDV in one P-DATA-TF PDU by default.");
        opts.addOption("tcpdelay", false,
                "set TCP_NODELAY socket option to false, true by default");

        OptionBuilder.withArgName("ms");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "delay in ms for Socket close after sending A-ABORT, 50ms by default");
        opts.addOption(OptionBuilder.create("soclosedelay"));

        OptionBuilder.withArgName("ms");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "delay in ms for DIMSE-RSP; useful for testing asynchronous mode");
        opts.addOption(OptionBuilder.create("rspdelay"));

        OptionBuilder.withArgName("ms");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "timeout in ms for receiving -ASSOCIATE-RQ, 5s by default");
        opts.addOption(OptionBuilder.create("requestTO"));

        OptionBuilder.withArgName("ms");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "timeout in ms for receiving A-RELEASE-RP, 5s by default");
        opts.addOption(OptionBuilder.create("releaseTO"));

        OptionBuilder.withArgName("ms");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "period in ms to check for outstanding DIMSE-RSP, 10s by default");
        opts.addOption(OptionBuilder.create("reaper"));

        OptionBuilder.withArgName("ms");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "timeout in ms for receiving DIMSE-RQ, 60s by default");
        opts.addOption(OptionBuilder.create("idleTO"));

        OptionBuilder.withArgName("KB");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "maximal length in KB of received P-DATA-TF PDUs, 16KB by default");
        opts.addOption(OptionBuilder.create("rcvpdulen"));

        OptionBuilder.withArgName("KB");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "maximal length in KB of sent P-DATA-TF PDUs, 16KB by default");
        opts.addOption(OptionBuilder.create("sndpdulen"));

        OptionBuilder.withArgName("KB");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "set SO_RCVBUF socket option to specified value in KB");
        opts.addOption(OptionBuilder.create("sorcvbuf"));

        OptionBuilder.withArgName("KB");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "set SO_SNDBUF socket option to specified value in KB");
        opts.addOption(OptionBuilder.create("sosndbuf"));

        OptionBuilder.withArgName("KB");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "minimal buffer size to write received object to file, 1KB by default");
        opts.addOption(OptionBuilder.create("bufsize"));

        opts.addOption("h", "help", false, "print this message");
        opts.addOption("V", "version", false,
                "print the version information and exit");
        CommandLine cl = null;
        try {
            cl = new GnuParser().parse(opts, args);
        } catch (ParseException e) {
            exit("dcmrcv: " + e.getMessage());
            throw new RuntimeException("unreachable");
        }
        if (cl.hasOption("V")) {
            Package p = DcmRcv.class.getPackage();
            System.out.println("dcmrcv v" + p.getImplementationVersion());
            System.exit(0);
        }
        if (cl.hasOption("h") || cl.getArgList().size() == 0) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(USAGE, DESCRIPTION, opts, EXAMPLE);
            System.exit(0);
        }
        return cl;
    }

    public void run(String[] args) {
    	System.out.println(args[0] + " " + args[1] + " " + args[2]);
        CommandLine cl = parse(args);
        DcmRcv dcmrcv = new DcmRcv();
        final List argList = cl.getArgList();
        String port = (String) argList.get(0);
        String[] aetPort = split(port, ':', 1);
        dcmrcv.setPort(parseInt(aetPort[1], "illegal port number", 1, 0xffff));
        if (aetPort[0] != null) {
            String[] aetHost = split(aetPort[0], '@', 0);
            dcmrcv.setAEtitle(aetHost[0]);
            dcmrcv.setHostname(aetHost[1]);
        }

        if (cl.hasOption("dest"))
            dcmrcv.setDestination(cl.getOptionValue("dest"));
        if (cl.hasOption("defts"))
            dcmrcv.setTransferSyntax(ONLY_DEF_TS);
        else if (cl.hasOption("native"))
            dcmrcv.setTransferSyntax(cl.hasOption("bigendian") ? NATIVE_TS
                    : NATIVE_LE_TS);
        else if (cl.hasOption("bigendian"))
            dcmrcv.setTransferSyntax(NON_RETIRED_TS);
        if (cl.hasOption("reaper"))
            dcmrcv.setAssociationReaperPeriod(parseInt(cl
                            .getOptionValue("reaper"),
                            "illegal argument of option -reaper", 1,
                            Integer.MAX_VALUE));
        if (cl.hasOption("idleTO"))
            dcmrcv.setIdleTimeout(parseInt(cl.getOptionValue("idleTO"),
                            "illegal argument of option -idleTO", 1,
                            Integer.MAX_VALUE));
        if (cl.hasOption("requestTO"))
            dcmrcv.setRequestTimeout(parseInt(cl.getOptionValue("requestTO"),
                    "illegal argument of option -requestTO", 1,
                    Integer.MAX_VALUE));
        if (cl.hasOption("releaseTO"))
            dcmrcv.setReleaseTimeout(parseInt(cl.getOptionValue("releaseTO"),
                    "illegal argument of option -releaseTO", 1,
                    Integer.MAX_VALUE));
        if (cl.hasOption("soclosedelay"))
            dcmrcv.setSocketCloseDelay(parseInt(cl
                    .getOptionValue("soclosedelay"),
                    "illegal argument of option -soclosedelay", 1, 10000));
        if (cl.hasOption("rspdelay"))
            dcmrcv.setDimseRspDelay(parseInt(cl.getOptionValue("rspdelay"),
                    "illegal argument of option -rspdelay", 0, 10000));
        if (cl.hasOption("rcvpdulen"))
            dcmrcv.setMaxPDULengthReceive(parseInt(cl
                    .getOptionValue("rcvpdulen"),
                    "illegal argument of option -rcvpdulen", 1, 10000)
                    * KB);
        if (cl.hasOption("sndpdulen"))
            dcmrcv.setMaxPDULengthSend(parseInt(cl.getOptionValue("sndpdulen"),
                    "illegal argument of option -sndpdulen", 1, 10000)
                    * KB);
        if (cl.hasOption("sosndbuf"))
            dcmrcv.setSendBufferSize(parseInt(cl.getOptionValue("sosndbuf"),
                    "illegal argument of option -sosndbuf", 1, 10000)
                    * KB);
        if (cl.hasOption("sorcvbuf"))
            dcmrcv.setReceiveBufferSize(parseInt(cl.getOptionValue("sorcvbuf"),
                    "illegal argument of option -sorcvbuf", 1, 10000)
                    * KB);
        if (cl.hasOption("bufsize"))
            dcmrcv.setFileBufferSize(parseInt(cl.getOptionValue("bufsize"),
                    "illegal argument of option -bufsize", 1, 10000)
                    * KB);

        dcmrcv.setPackPDV(!cl.hasOption("pdv1"));
        dcmrcv.setTcpNoDelay(!cl.hasOption("tcpdelay"));
        if (cl.hasOption("async"))
            dcmrcv.setMaxOpsPerformed(parseInt(cl.getOptionValue("async"),
                    "illegal argument of option -async", 0, 0xffff));
        dcmrcv.initTransferCapability();
        if (cl.hasOption("tls")) {
            String cipher = cl.getOptionValue("tls");
            if ("NULL".equalsIgnoreCase(cipher)) {
                dcmrcv.setTlsWithoutEncyrption();
            } else if ("3DES".equalsIgnoreCase(cipher)) {
                dcmrcv.setTls3DES_EDE_CBC();
            } else if ("AES".equalsIgnoreCase(cipher)) {
                dcmrcv.setTlsAES_128_CBC();
            } else {
                exit("Invalid parameter for option -tls: " + cipher);
            }
            if (cl.hasOption("nossl2")) {
                dcmrcv.disableSSLv2Hello();
            }
            dcmrcv.setTlsNeedClientAuth(!cl.hasOption("noclientauth"));

            if (cl.hasOption("keystore")) {
                dcmrcv.setKeyStoreURL(cl.getOptionValue("keystore"));
            }
            if (cl.hasOption("keystorepw")) {
                dcmrcv.setKeyStorePassword(
                        cl.getOptionValue("keystorepw"));
            }
            if (cl.hasOption("keypw")) {
                dcmrcv.setKeyPassword(cl.getOptionValue("keypw"));
            }
            if (cl.hasOption("truststore")) {
                dcmrcv.setTrustStoreURL(
                        cl.getOptionValue("truststore"));
            }
            if (cl.hasOption("truststorepw")) {
                dcmrcv.setTrustStorePassword(
                        cl.getOptionValue("truststorepw"));
            }
            try {
                dcmrcv.initTLS();
            } catch (Exception e) {
                System.err.println("ERROR: Failed to initialize TLS context:"
                        + e.getMessage());
                System.exit(2);
            }
        }        
        try {
            dcmrcv.start();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void setTransferSyntax(String[] tsuids) {
        this.tsuids = tsuids;
    }

    private void initTransferCapability() {
        TransferCapability[] tc = new TransferCapability[CUIDS.length + 1];
        tc[0] = new TransferCapability(UID.VerificationSOPClass, ONLY_DEF_TS,
                TransferCapability.SCP);
        for (int i = 0; i < CUIDS.length; i++)
            tc[i + 1] = new TransferCapability(CUIDS[i], tsuids,
                    TransferCapability.SCP);
        ae.setTransferCapability(tc);
    }

    private void setFileBufferSize(int size) {
        fileBufferSize = size;
    }

    private void setMaxOpsPerformed(int maxOps) {
        ae.setMaxOpsPerformed(maxOps);
    }

    private void setDestination(String filePath) {
        this.destination = new File(filePath);
        this.devnull = "/dev/null".equals(filePath);
        if (!devnull)
            destination.mkdir();
    }

    public void initTLS() throws GeneralSecurityException, IOException {
        KeyStore keyStore = loadKeyStore(keyStoreURL, keyStorePassword);
        KeyStore trustStore = loadKeyStore(trustStoreURL, trustStorePassword);
        device.initTLS(keyStore,
                keyPassword != null ? keyPassword : keyStorePassword,
                trustStore);
    }
    
    private static KeyStore loadKeyStore(String url, char[] password)
            throws GeneralSecurityException, IOException {
        KeyStore key = KeyStore.getInstance(toKeyStoreType(url));
        InputStream in = openFileOrURL(url);
        try {
            key.load(in, password);
        } finally {
            in.close();
        }
        return key;
    }

    private static InputStream openFileOrURL(String url) throws IOException {
        if (url.startsWith("resource:")) {
            return DcmRcv.class.getClassLoader().getResourceAsStream(
                    url.substring(9));
        }
        try {
            return new URL(url).openStream();
        } catch (MalformedURLException e) {
            return new FileInputStream(url);
        }
    }

    private static String toKeyStoreType(String fname) {
        return fname.endsWith(".p12") || fname.endsWith(".P12")
                 ? "PKCS12" : "JKS";
    }
    
    public void start() throws IOException {
        device.startListening(executor);
        System.out.println("Start Server listening on port " + nc.getPort());
    }
    
    public void stop(){
    	//device.stopListening();
    	/*
    	NetworkConnection[] connections = device.getNetworkConnection();
    	for(NetworkConnection connection : connections)
    		connection.unbind();
    		*/
    }

    private static String[] split(String s, char delim, int defPos) {
        String[] s2 = new String[2];
        s2[defPos] = s;
        int pos = s.indexOf(delim);
        if (pos != -1) {
            s2[0] = s.substring(0, pos);
            s2[1] = s.substring(pos + 1);
        }
        return s2;
    }

    private static void exit(String msg) {
        System.err.println(msg);
        System.err.println("Try 'dcmrcv -h' for more information.");
        System.exit(1);
    }

    private static int parseInt(String s, String errPrompt, int min, int max) {
        try {
            int i = Integer.parseInt(s);
            if (i >= min && i <= max)
                return i;
        } catch (NumberFormatException e) {
            // parameter is not a valid integer; fall through to exit
        }
        exit(errPrompt);
        throw new RuntimeException();
    }

    /** Overwrite {@link StorageService#cstore} to send delayed C-STORE RSP 
     * by separate Thread, so reading of following received C-STORE RQs from
     * the open association is not blocked.
     */
    @Override
    public void cstore(final Association as, final int pcid, DicomObject rq, 
            PDVInputStream dataStream, String tsuid) 
            throws DicomServiceException, IOException {
        final DicomObject rsp = CommandUtils.mkRSP(rq, CommandUtils.SUCCESS);
        onCStoreRQ(as, pcid, rq, dataStream, tsuid, rsp);
        if (rspdelay > 0) {
            executor.execute(new Runnable() {
                public void run() {
                    try {
                        Thread.sleep(rspdelay);
                        as.writeDimseRSP(pcid, rsp);
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }
            });
        } else {
            as.writeDimseRSP(pcid, rsp);
        }
        onCStoreRSP(as, pcid, rq, dataStream, tsuid, rsp);
    }

    
    @Override
    protected void onCStoreRQ(Association as, int pcid, DicomObject rq,
            PDVInputStream dataStream, String tsuid, DicomObject rsp)
            throws IOException, DicomServiceException {
        if (destination == null) {
            super.onCStoreRQ(as, pcid, rq, dataStream, tsuid, rsp);
        } else {
            try {
                String cuid = rq.getString(Tag.AffectedSOPClassUID);
                String iuid = rq.getString(Tag.AffectedSOPInstanceUID);
                BasicDicomObject fmi = new BasicDicomObject();
                fmi.initFileMetaInformation(cuid, iuid, tsuid);
                File file = devnull ? destination : new File(destination, iuid);
                FileOutputStream fos = new FileOutputStream(file);
                BufferedOutputStream bos = new BufferedOutputStream(fos,
                        fileBufferSize);
                DicomOutputStream dos = new DicomOutputStream(bos);
                dos.writeFileMetaInformation(fmi);
                dataStream.copyTo(dos);
                dos.close();
            } catch (IOException e) {
                throw new DicomServiceException(rq, Status.ProcessingFailure, e
                        .getMessage());
            }
        }
    }

}

