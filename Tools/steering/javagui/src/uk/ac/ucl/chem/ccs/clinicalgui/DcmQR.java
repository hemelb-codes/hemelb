package uk.ac.ucl.chem.ccs.clinicalgui;

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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.dcm4che2.data.BasicDicomObject;
import org.dcm4che2.data.DicomElement;
import org.dcm4che2.data.DicomObject;
import org.dcm4che2.data.Tag;
import org.dcm4che2.data.UID;
import org.dcm4che2.data.UIDDictionary;
import org.dcm4che2.data.VR;
import org.dcm4che2.io.DicomOutputStream;
import org.dcm4che2.net.Association;
import org.dcm4che2.net.CommandUtils;
import org.dcm4che2.net.ConfigurationException;
import org.dcm4che2.net.Device;
import org.dcm4che2.net.DicomServiceException;
import org.dcm4che2.net.DimseRSP;
import org.dcm4che2.net.DimseRSPHandler;
import org.dcm4che2.net.Executor;
import org.dcm4che2.net.ExtQueryTransferCapability;
import org.dcm4che2.net.ExtRetrieveTransferCapability;
import org.dcm4che2.net.NetworkApplicationEntity;
import org.dcm4che2.net.NetworkConnection;
import org.dcm4che2.net.NewThreadExecutor;
import org.dcm4che2.net.NoPresentationContextException;
import org.dcm4che2.net.PDVInputStream;
import org.dcm4che2.net.Status;
import org.dcm4che2.net.TransferCapability;
import org.dcm4che2.net.UserIdentity;
import org.dcm4che2.net.service.DicomService;
import org.dcm4che2.net.service.StorageService;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author gunter zeilinger(gunterze@gmail.com)
 * @version $Revision$ $Date$
 * @since Jan, 2006
 */
public class DcmQR {
    private static Logger LOG = LoggerFactory.getLogger(DcmQR.class);

    private static final int KB = 1024;

    private static final String USAGE = "dcmqr <aet>[@<host>[:<port>]] [Options]";

    private static final String DESCRIPTION = 
        "Query specified remote Application Entity (=Query/Retrieve SCP) "
        + "and optional (s. option -cget/-cmove) retrieve instances of "
        + "matching entities. If <port> is not specified, DICOM default port "
        + "104 is assumed. If also no <host> is specified localhost is assumed. "
        + "Also Storage Services can be provided (s. option -cstore) to receive "
        + "retrieved instances. For receiving objects retrieved by C-MOVE in a "
        + "separate association, a local listening port must be specified "
        + "(s.option -L).\n"
        + "Options:";
    
    private static final String EXAMPLE = 
        "\nExample: dcmqr -L QRSCU:11113 QRSCP@localhost:11112 -cmove QRSCU " +
        "-qStudyDate=20060204 -qModalitiesInStudy=CT -cstore CT -cstore PR:LE " +
        "-cstoredest /tmp\n"
        + "=> Query Application Entity QRSCP listening on local port 11112 for "
        + "CT studies from Feb 4, 2006 and retrieve matching studies by C-MOVE "
        + "to own Application Entity QRSCU listing on local port 11113, "
        + "storing received CT images and Grayscale Softcopy Presentation "
        + "states to /tmp.";

    private static char[] SECRET = { 's', 'e', 'c', 'r', 'e', 't' };

    
    
    private static enum QueryRetrieveLevel {
        PATIENT("PATIENT", PATIENT_RETURN_KEYS, PATIENT_LEVEL_FIND_CUID,
                PATIENT_LEVEL_GET_CUID, PATIENT_LEVEL_MOVE_CUID),
        STUDY("STUDY", STUDY_RETURN_KEYS, STUDY_LEVEL_FIND_CUID,
                STUDY_LEVEL_GET_CUID, STUDY_LEVEL_MOVE_CUID),
        SERIES("SERIES", SERIES_RETURN_KEYS, SERIES_LEVEL_FIND_CUID,
                SERIES_LEVEL_GET_CUID, SERIES_LEVEL_MOVE_CUID),
        IMAGE("IMAGE", INSTANCE_RETURN_KEYS, SERIES_LEVEL_FIND_CUID,
                SERIES_LEVEL_GET_CUID, SERIES_LEVEL_MOVE_CUID);
        
        private final String code;
        private final int[] returnKeys;
        private final String[] findClassUids;
        private final String[] getClassUids;
        private final String[] moveClassUids;

        private QueryRetrieveLevel(String code, int[] returnKeys,
                String[] findClassUids, String[] getClassUids,
                String[] moveClassUids) {
            this.code = code;
            this.returnKeys = returnKeys;
            this.findClassUids = findClassUids;
            this.getClassUids = getClassUids;
            this.moveClassUids = moveClassUids;
        }
        
        public String getCode() {
            return code;
        }

        public int[] getReturnKeys() {
            return returnKeys;
        }
        
        public String[] getFindClassUids() {
            return findClassUids;
        }
        
        public String[] getGetClassUids() {
            return getClassUids;
        }
        
        public String[] getMoveClassUids() {
            return moveClassUids;
        }
    }

    private static final String[] PATIENT_LEVEL_FIND_CUID = {
        UID.PatientRootQueryRetrieveInformationModelFIND,
        UID.PatientStudyOnlyQueryRetrieveInformationModelFINDRetired };

    private static final String[] STUDY_LEVEL_FIND_CUID = {
        UID.StudyRootQueryRetrieveInformationModelFIND,
        UID.PatientRootQueryRetrieveInformationModelFIND,
        UID.PatientStudyOnlyQueryRetrieveInformationModelFINDRetired };

    private static final String[] SERIES_LEVEL_FIND_CUID = {
        UID.StudyRootQueryRetrieveInformationModelFIND,
        UID.PatientRootQueryRetrieveInformationModelFIND, };

    private static final String[] PATIENT_LEVEL_GET_CUID = {
        UID.PatientRootQueryRetrieveInformationModelGET,
        UID.PatientStudyOnlyQueryRetrieveInformationModelGETRetired };
    
    private static final String[] STUDY_LEVEL_GET_CUID = {
        UID.StudyRootQueryRetrieveInformationModelGET,
        UID.PatientRootQueryRetrieveInformationModelGET,
        UID.PatientStudyOnlyQueryRetrieveInformationModelGETRetired };

    private static final String[] SERIES_LEVEL_GET_CUID = {
        UID.StudyRootQueryRetrieveInformationModelGET,
        UID.PatientRootQueryRetrieveInformationModelGET };

    private static final String[] PATIENT_LEVEL_MOVE_CUID = {
        UID.PatientRootQueryRetrieveInformationModelMOVE,
        UID.PatientStudyOnlyQueryRetrieveInformationModelMOVERetired };

    private static final String[] STUDY_LEVEL_MOVE_CUID = {
        UID.StudyRootQueryRetrieveInformationModelMOVE,
        UID.PatientRootQueryRetrieveInformationModelMOVE,
        UID.PatientStudyOnlyQueryRetrieveInformationModelMOVERetired };

    private static final String[] SERIES_LEVEL_MOVE_CUID = {
        UID.StudyRootQueryRetrieveInformationModelMOVE,
        UID.PatientRootQueryRetrieveInformationModelMOVE };

    private static final int[] PATIENT_RETURN_KEYS = {
        Tag.PatientName,
        Tag.PatientID,
        Tag.PatientBirthDate,
        Tag.PatientSex,
        Tag.NumberOfPatientRelatedStudies,
        Tag.NumberOfPatientRelatedSeries,
        Tag.NumberOfPatientRelatedInstances };

    private static final int[] PATIENT_MATCHING_KEYS = { 
        Tag.PatientName,
        Tag.PatientID,
        Tag.IssuerOfPatientID,
        Tag.PatientBirthDate,
        Tag.PatientSex };

    private static final int[] STUDY_RETURN_KEYS = {
        Tag.StudyDate,
        Tag.StudyTime,
        Tag.AccessionNumber,
        Tag.StudyID,
        Tag.StudyInstanceUID,
        Tag.NumberOfStudyRelatedSeries,
        Tag.NumberOfStudyRelatedInstances };

    private static final int[] STUDY_MATCHING_KEYS = {
        Tag.StudyDate,
        Tag.StudyTime,
        Tag.AccessionNumber,
        Tag.ModalitiesInStudy,
        Tag.ReferringPhysicianName,
        Tag.StudyID,
        Tag.StudyInstanceUID };

    private static final int[] PATIENT_STUDY_MATCHING_KEYS = {
        Tag.StudyDate,
        Tag.StudyTime,
        Tag.AccessionNumber,
        Tag.ModalitiesInStudy,
        Tag.ReferringPhysicianName,
        Tag.PatientName,
        Tag.PatientID,
        Tag.IssuerOfPatientID,
        Tag.PatientBirthDate,
        Tag.PatientSex,
        Tag.StudyID,
        Tag.StudyInstanceUID };

    private static final int[] SERIES_RETURN_KEYS = {
        Tag.Modality,
        Tag.SeriesNumber,
        Tag.SeriesInstanceUID,
        Tag.NumberOfSeriesRelatedInstances };

    private static final int[] SERIES_MATCHING_KEYS = {
        Tag.Modality,
        Tag.SeriesNumber,
        Tag.SeriesInstanceUID,
        Tag.RequestAttributesSequence
    };

    private static final int[] INSTANCE_RETURN_KEYS = {
        Tag.InstanceNumber,
        Tag.SOPClassUID,
        Tag.SOPInstanceUID, };

    private static final int[] MOVE_KEYS = {
        Tag.QueryRetrieveLevel,
        Tag.PatientID,
        Tag.StudyInstanceUID,
        Tag.SeriesInstanceUID,
        Tag.SOPInstanceUID, };

    private static final String[] IVRLE_TS = {
        UID.ImplicitVRLittleEndian };

    private static final String[] NATIVE_LE_TS = {
        UID.ExplicitVRLittleEndian,
        UID.ImplicitVRLittleEndian};

    private static final String[] NATIVE_BE_TS = {
        UID.ExplicitVRBigEndian,
        UID.ImplicitVRLittleEndian};

    private static final String[] DEFLATED_TS = {
        UID.DeflatedExplicitVRLittleEndian,
        UID.ExplicitVRLittleEndian,
        UID.ImplicitVRLittleEndian};

    private static final String[] NOPX_TS = {
        UID.NoPixelData,
        UID.ExplicitVRLittleEndian,
        UID.ImplicitVRLittleEndian};

    private static final String[] NOPXDEFL_TS = {
        UID.NoPixelDataDeflate,
        UID.NoPixelData,
        UID.ExplicitVRLittleEndian,
        UID.ImplicitVRLittleEndian};

    private static final String[] JPLL_TS = {
        UID.JPEGLossless,
        UID.JPEGLosslessNonHierarchical14,
        UID.JPEGLSLossless,
        UID.JPEG2000LosslessOnly,
        UID.ExplicitVRLittleEndian,
        UID.ImplicitVRLittleEndian};

    private static final String[] JPLY_TS = {
        UID.JPEGBaseline1,
        UID.JPEGExtended24,
        UID.JPEGLSLossyNearLossless,
        UID.JPEG2000,
        UID.ExplicitVRLittleEndian,
        UID.ImplicitVRLittleEndian};

    private static final String[] MPEG2_TS = { UID.MPEG2 };

    private static final String[] DEF_TS = {
        UID.JPEGLossless,
        UID.JPEGLosslessNonHierarchical14,
        UID.JPEGLSLossless,
        UID.JPEGLSLossyNearLossless,
        UID.JPEG2000LosslessOnly,
        UID.JPEG2000,
        UID.JPEGBaseline1,
        UID.JPEGExtended24,
        UID.MPEG2,
        UID.DeflatedExplicitVRLittleEndian,
        UID.ExplicitVRBigEndian,
        UID.ExplicitVRLittleEndian,
        UID.ImplicitVRLittleEndian};

    private static enum TS {
        IVLE(IVRLE_TS),
        LE(NATIVE_LE_TS),
        BE(NATIVE_BE_TS),
        DEFL(DEFLATED_TS),
        JPLL(JPLL_TS),
        JPLY(JPLY_TS),
        MPEG2(MPEG2_TS),
        NOPX(NOPX_TS),
        NOPXD(NOPXDEFL_TS);
        
        final String[] uids;
        TS(String[] uids) { this.uids = uids; }
    }

    private static enum CUID {
        CR(UID.ComputedRadiographyImageStorage),
        CT(UID.CTImageStorage),
        MR(UID.MRImageStorage),
        US(UID.UltrasoundImageStorage),
        NM(UID.NuclearMedicineImageStorage),
        SC(UID.SecondaryCaptureImageStorage),
        XA(UID.XRayAngiographicImageStorage),
        XRF(UID.XRayRadiofluoroscopicImageStorage),
        DX(UID.DigitalXRayImageStorageForPresentation),
        MG(UID.DigitalMammographyXRayImageStorageForPresentation),
        PR(UID.GrayscaleSoftcopyPresentationStateStorageSOPClass),
        KO(UID.KeyObjectSelectionDocument),
        SR(UID.BasicTextSR);

        final String uid;
        CUID(String uid) { this.uid = uid; }
        
    }

    private static final String[] EMPTY_STRING = {};

    private Executor executor = new NewThreadExecutor("DCMQR");

    private NetworkApplicationEntity remoteAE = new NetworkApplicationEntity();

    private NetworkConnection remoteConn = new NetworkConnection();

    private Device device = new Device("DCMQR");

    private NetworkApplicationEntity ae = new NetworkApplicationEntity();

    private NetworkConnection conn = new NetworkConnection();

    private Association assoc;

    private int priority = 0;

    private boolean cget; 

    private String moveDest;

    private File storeDest;

    private boolean devnull;

    private int fileBufferSize = 256;

    private boolean evalRetrieveAET = false;

    private QueryRetrieveLevel qrlevel = QueryRetrieveLevel.STUDY;

    private List<String> privateFind = new ArrayList<String>();

    private final List<TransferCapability> storeTransferCapability =
            new ArrayList<TransferCapability>(8);

    private DicomObject keys = new BasicDicomObject();

    private int cancelAfter = Integer.MAX_VALUE;

    private int completed;

    private int warning;

    private int failed;

    private boolean relationQR;

    private boolean dateTimeMatching;

    private boolean fuzzySemanticPersonNameMatching;

    private boolean noExtNegotiation;

    private String keyStoreURL = "resource:tls/test_sys_1.p12";
    
    private char[] keyStorePassword = SECRET; 

    private char[] keyPassword; 
    
    private String trustStoreURL = "resource:tls/mesa_certs.jks";
    
    private char[] trustStorePassword = SECRET;
    
    public DcmQR() {
        remoteAE.setInstalled(true);
        remoteAE.setAssociationAcceptor(true);
        remoteAE.setNetworkConnection(new NetworkConnection[] { remoteConn });

        device.setNetworkApplicationEntity(ae);
        device.setNetworkConnection(conn);
        ae.setNetworkConnection(conn);
        ae.setAssociationInitiator(true);
        ae.setAssociationAcceptor(true);
        ae.setAETitle("DCMQR");
    }

    public final void setLocalHost(String hostname) {
        conn.setHostname(hostname);
    }

    public final void setLocalPort(int port) {
        conn.setPort(port);
    }

    public final void setRemoteHost(String hostname) {
        remoteConn.setHostname(hostname);
    }

    public final void setRemotePort(int port) {
        remoteConn.setPort(port);
    }

    public final void setTlsWithoutEncyrption() {
        conn.setTlsWithoutEncyrption();
        remoteConn.setTlsWithoutEncyrption();
    }

    public final void setTls3DES_EDE_CBC() {
        conn.setTls3DES_EDE_CBC();
        remoteConn.setTls3DES_EDE_CBC();
    }

    public final void setTlsAES_128_CBC() {
        conn.setTlsAES_128_CBC();
        remoteConn.setTlsAES_128_CBC();
    }

    public final void disableSSLv2Hello() {
        conn.disableSSLv2Hello();
    }

    public final void setTlsNeedClientAuth(boolean needClientAuth) {
        conn.setTlsNeedClientAuth(needClientAuth);
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

    public final void setCalledAET(String called, boolean reuse) {
        remoteAE.setAETitle(called);
        if (reuse)
            ae.setReuseAssocationToAETitle(new String[] { called });
    }

    public final void setCalling(String calling) {
        ae.setAETitle(calling);
    }

    public final void setUserIdentity(UserIdentity userIdentity) {
        ae.setUserIdentity(userIdentity);
    }

    public final void setPriority(int priority) {
        this.priority = priority;
    }

    public final void setConnectTimeout(int connectTimeout) {
        conn.setConnectTimeout(connectTimeout);
    }

    public final void setMaxPDULengthReceive(int maxPDULength) {
        ae.setMaxPDULengthReceive(maxPDULength);
    }

    public final void setMaxOpsInvoked(int maxOpsInvoked) {
        ae.setMaxOpsInvoked(maxOpsInvoked);
    }

    public final void setMaxOpsPerformed(int maxOps) {
        ae.setMaxOpsPerformed(maxOps);
    }

    public final void setPackPDV(boolean packPDV) {
        ae.setPackPDV(packPDV);
    }

    public final void setAssociationReaperPeriod(int period) {
        device.setAssociationReaperPeriod(period);
    }

    public final void setDimseRspTimeout(int timeout) {
        ae.setDimseRspTimeout(timeout);
    }

    public final void setTcpNoDelay(boolean tcpNoDelay) {
        conn.setTcpNoDelay(tcpNoDelay);
    }

    public final void setAcceptTimeout(int timeout) {
        conn.setAcceptTimeout(timeout);
    }

    public final void setReleaseTimeout(int timeout) {
        conn.setReleaseTimeout(timeout);
    }

    public final void setSocketCloseDelay(int timeout) {
        conn.setSocketCloseDelay(timeout);
    }

    public final void setMaxPDULengthSend(int maxPDULength) {
        ae.setMaxPDULengthSend(maxPDULength);
    }

    public final void setReceiveBufferSize(int bufferSize) {
        conn.setReceiveBufferSize(bufferSize);
    }

    public final void setSendBufferSize(int bufferSize) {
        conn.setSendBufferSize(bufferSize);
    }

    public final void setFileBufferSize(int size) {
        fileBufferSize = size;
    }

    private static CommandLine parse(String[] args) {
        Options opts = new Options();
        OptionBuilder.withArgName("aet[@host][:port]");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "set AET, local address and listening port of local Application Entity");
        opts.addOption(OptionBuilder.create("L"));

        OptionBuilder.withArgName("username");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "enable User Identity Negotiation with specified username and "
                + " optional passcode");
        opts.addOption(OptionBuilder.create("username"));

        OptionBuilder.withArgName("passcode");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "optional passcode for User Identity Negotiation, "
                + "only effective with option -username");
        opts.addOption(OptionBuilder.create("passcode"));

        opts.addOption("uidnegrsp", false,
                "request positive User Identity Negotation response, "
                + "only effective with option -username");
        
        OptionBuilder.withArgName("NULL|3DES|AES");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "enable TLS connection without, 3DES or AES encryption");
        opts.addOption(OptionBuilder.create("tls"));

        opts.addOption("nossl2", false, 
                "disable SSLv2Hello TLS handshake");
        opts.addOption("noclientauth", false, 
                "disable client authentification for TLS");

        OptionBuilder.withArgName("file|url");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "file path or URL of P12 or JKS keystore, resource:tls/test_sys_1.p12 by default");
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
        
        OptionBuilder.withArgName("aet");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "retrieve instances of matching entities by C-MOVE to specified destination.");
        opts.addOption(OptionBuilder.create("cmove"));

        opts.addOption("cget", false, "retrieve instances of matching entities by C-GET.");

        OptionBuilder.withArgName("cuid[:ts]");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription(
                "negotiate support of specified Storage SOP Class and Transfer "
                + "Syntaxes. The Storage SOP\nClass may be specified by its UID "
                + "or by one\nof following key words:\n"
                + "CR  - Computed Radiography Image Storage\n"
                + "CT  - CT Image Storage\n"
                + "MR  - MRImageStorage\n"
                + "US  - Ultrasound Image Storage\n"
                + "NM  - Nuclear Medicine Image Storage\n"
                + "SC  - Secondary Capture Image Storage\n"
                + "XA  - XRay Angiographic Image Storage\n"
                + "XRF - XRay Radiofluoroscopic Image Storage\n"
                + "DX  - Digital X-Ray Image Storage for Presentation\n"
                + "                            MG  - Digital Mammography X-Ray Image Storage\n"
                + "for Presentation\n"
                + "PR  - Grayscale Softcopy Presentation State Storage\n"
                + "                            KO  - Key Object Selection Document Storage\n"
                + "SR  - Basic Text Structured Report Document Storage\n"
                + "                            The Transfer Syntaxes may be specified by a comma\n"
                + "                            separated list of UIDs or by one of following key\n"
                + "                            words:\n"
                + "                            IVRLE - offer only Implicit VR Little Endian\n"
                + "                            Transfer Syntax\n"
                + "                            LE - offer Explicit and Implicit VR Little Endian\n"
                + "                            Transfer Syntax\n"
                + "                            BE - offer Explicit VR Big Endian Transfer Syntax\n"
                + "                            DEFL - offer Deflated Explicit VR Little\n"
                + "                            Endian Transfer Syntax\n"
                + "                            JPLL - offer JEPG Loss Less Transfer Syntaxes\n"
                + "                            JPLY - offer JEPG Lossy Transfer Syntaxes\n"
                + "                            MPEG2 - offer MPEG2 Transfer Syntax\n"
                + "                            NOPX - offer No Pixel Data Transfer Syntax\n"
                + "                            NOPXD - offer No Pixel Data Delflate Transfer Syntax\n"
                + "                            If only the Storage SOP Class is specified, all\n"
                + "                            Transfer Syntaxes listed above except No Pixel Data\n"
                + "                            and No Pixel Data Delflate Transfer Syntax are\n"
                + "                            offered.");
        opts.addOption(OptionBuilder.create("cstore"));

        OptionBuilder.withArgName("dir");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "store received objects into files in specified directory <dir>."
                        + " Do not store received objects\nby default.");
        opts.addOption(OptionBuilder.create("cstoredest"));

        opts.addOption("ivrle", false, "offer only Implicit VR Little Endian Transfer Syntax.");

        OptionBuilder.withArgName("maxops");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("maximum number of outstanding C-MOVE-RQ " +
                "it may invoke asynchronously, 1 by default.");
        opts.addOption(OptionBuilder.create("async"));

        OptionBuilder.withArgName("maxops");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "maximum number of outstanding storage operations performed "
                        + "asynchronously, unlimited by\n                            default.");
        opts.addOption(OptionBuilder.create("storeasync"));

        opts.addOption("noextneg", false, "disable extended negotiation.");
        opts.addOption("rel", false,
                "negotiate support of relational queries and retrieval.");
        opts.addOption("datetime", false,
                "negotiate support of combined date and time attribute range matching.");
        opts.addOption("fuzzy", false, 
                "negotiate support of fuzzy semantic person name attribute matching.");

        opts.addOption("retall", false, "negotiate private FIND SOP Classes " +
                "to fetch all available attributes of matching entities.");
        opts.addOption("blocked", false, "negotiate private FIND SOP Classes " +
                "to return attributes of several matching entities per FIND\n" +
                "                            response.");
        opts.addOption("vmf", false, "negotiate private FIND SOP Classes to " +
                "return attributes of legacy CT/MR images of one series as\n" +
                "                           virtual multiframe object.");
        opts.addOption("pdv1", false,
                "send only one PDV in one P-Data-TF PDU, pack command and data "
                + "PDV in one P-DATA-TF PDU\n"
                + "                           by default.");
        opts.addOption("tcpdelay", false,
                "set TCP_NODELAY socket option to false, true by default");

        OptionBuilder.withArgName("ms");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "timeout in ms for TCP connect, no timeout by default");
        opts.addOption(OptionBuilder.create("connectTO"));

        OptionBuilder.withArgName("ms");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "delay in ms for Socket close after sending A-ABORT, 50ms by default");
        opts.addOption(OptionBuilder.create("soclosedelay"));

        OptionBuilder.withArgName("ms");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "period in ms to check for outstanding DIMSE-RSP, 10s by default");
        opts.addOption(OptionBuilder.create("reaper"));

        OptionBuilder.withArgName("ms");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "timeout in ms for receiving DIMSE-RSP, 60s by default");
        opts.addOption(OptionBuilder.create("rspTO"));

        OptionBuilder.withArgName("ms");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "timeout in ms for receiving A-ASSOCIATE-AC, 5s by default");
        opts.addOption(OptionBuilder.create("acceptTO"));

        OptionBuilder.withArgName("ms");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "timeout in ms for receiving A-RELEASE-RP, 5s by default");
        opts.addOption(OptionBuilder.create("releaseTO"));

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
        opts.addOption(OptionBuilder.create("filebuf"));

        OptionGroup qrlevel = new OptionGroup();

        OptionBuilder.withDescription("perform patient level query, multiple "
                + "exclusive with -S and -I, perform study level query\n"
                + "                            by default.");
        OptionBuilder.withLongOpt("patient");
        opts.addOption(OptionBuilder.create("P"));

        OptionBuilder.withDescription("perform series level query, multiple "
                + "exclusive with -P and -I, perform study level query\n"
                + "                            by default.");
        OptionBuilder.withLongOpt("series");
        opts.addOption(OptionBuilder.create("S"));

        OptionBuilder.withDescription("perform instance level query, multiple "
                + "exclusive with -P and -S, perform study level query\n"
                + "                            by default.");
        OptionBuilder.withLongOpt("image");
        opts.addOption(OptionBuilder.create("I"));

        opts.addOptionGroup(qrlevel);

        OptionBuilder.withArgName("[seq/]attr=value");
        OptionBuilder.hasArgs();
        OptionBuilder.withValueSeparator('=');
        OptionBuilder.withDescription("specify matching key. attr can be " +
                "specified by name or tag value (in hex), e.g. PatientName\n" +
                "or 00100010. Attributes in nested Datasets can\n" +
                "be specified by including the name/tag value of\n" +
                "                            the sequence attribute, e.g. 00400275/00400009\n" +
                "for Scheduled Procedure Step ID in the Request\n" +
                "Attributes Sequence");
        opts.addOption(OptionBuilder.create("q"));

        OptionBuilder.withArgName("attr");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("specify additional return key. attr can " +
                "be specified by name or tag value (in hex).");
        opts.addOption(OptionBuilder.create("r"));

        OptionBuilder.withArgName("num");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("cancel query after receive of specified " +
                "number of responses, no cancel by default");
        opts.addOption(OptionBuilder.create("C"));

        OptionBuilder.withArgName("aet");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("retrieve matching objects to specified " +
                "move destination.");
        opts.addOption(OptionBuilder.create("cmove"));
        
        opts.addOption("evalRetrieveAET", false,
                "Only Move studies not allready stored on destination AET");
        opts.addOption("lowprior", false,
                "LOW priority of the C-FIND/C-MOVE operation, MEDIUM by default");
        opts.addOption("highprior", false,
                "HIGH priority of the C-FIND/C-MOVE operation, MEDIUM by default");

        OptionBuilder.withArgName("num");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("repeat query (and retrieve) several times");
        opts.addOption(OptionBuilder.create("repeat"));

        OptionBuilder.withArgName("ms");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription(
                "delay in ms between repeated query (and retrieve), no delay by default");
        opts.addOption(OptionBuilder.create("repeatdelay"));

        opts.addOption("reuseassoc", false,
                "Reuse association for repeated query (and retrieve)");
        opts.addOption("closeassoc", false,
                "Close association between repeated query (and retrieve)");

        opts.addOption("h", "help", false, "print this message");
        opts.addOption("V", "version", false,
                "print the version information and exit");
        CommandLine cl = null;
        try {
            cl = new GnuParser().parse(opts, args);
        } catch (ParseException e) {
            exit("dcmqr: " + e.getMessage());
            throw new RuntimeException("unreachable");
        }
        if (cl.hasOption('V')) {
            Package p = DcmQR.class.getPackage();
            System.out.println("dcmqr v" + p.getImplementationVersion());
            System.exit(0);
        }
        if (cl.hasOption('h') || cl.getArgList().size() != 1) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(USAGE, DESCRIPTION, opts, EXAMPLE);
            System.exit(0);
        }

        return cl;
    }

    @SuppressWarnings("unchecked")
    public static List<DicomObject> run(String[] args) {
    	//System.out.println(args[0] + " " + args[1] + " " + args[2] + " " + args[3]);
        CommandLine cl = parse(args);
        DcmQR dcmqr = new DcmQR();
        final List<String> argList = cl.getArgList();
        String remoteAE = argList.get(0);
        String[] calledAETAddress = split(remoteAE, '@');
        dcmqr.setCalledAET(calledAETAddress[0], cl.hasOption("reuseassoc"));
        if (calledAETAddress[1] == null) {
            dcmqr.setRemoteHost("127.0.0.1");
            dcmqr.setRemotePort(104);
        } else {
            String[] hostPort = split(calledAETAddress[1], ':');
            dcmqr.setRemoteHost(hostPort[0]);
            dcmqr.setRemotePort(toPort(hostPort[1]));
        }
        if (cl.hasOption("L")) {
            String localAE = cl.getOptionValue("L");
            String[] localPort = split(localAE, ':');
            if (localPort[1] != null) {
                dcmqr.setLocalPort(toPort(localPort[1]));                
            }
            String[] callingAETHost = split(localPort[0], '@');
            dcmqr.setCalling(callingAETHost[0]);
            if (callingAETHost[1] != null) {
                dcmqr.setLocalHost(callingAETHost[0]);
            }
        }
        if (cl.hasOption("username")) {
            String username = cl.getOptionValue("username");
            UserIdentity userId;
            if (cl.hasOption("passcode")) {
                String passcode = cl.getOptionValue("passcode");
                userId = new UserIdentity.UsernamePasscode(username,
                        passcode.toCharArray());
            } else {
                userId = new UserIdentity.Username(username);
            }
            userId.setPositiveResponseRequested(cl.hasOption("uidnegrsp"));
            dcmqr.setUserIdentity(userId);
        }
        if (cl.hasOption("connectTO"))
            dcmqr.setConnectTimeout(parseInt(cl.getOptionValue("connectTO"),
                    "illegal argument of option -connectTO", 1,
                    Integer.MAX_VALUE));
        if (cl.hasOption("reaper"))
            dcmqr.setAssociationReaperPeriod(parseInt(cl.getOptionValue("reaper"),
                            "illegal argument of option -reaper", 1,
                            Integer.MAX_VALUE));
        if (cl.hasOption("rspTO"))
            dcmqr.setDimseRspTimeout(parseInt(cl.getOptionValue("rspTO"),
                    "illegal argument of option -rspTO", 1, Integer.MAX_VALUE));
        if (cl.hasOption("acceptTO"))
            dcmqr.setAcceptTimeout(parseInt(cl.getOptionValue("acceptTO"),
                    "illegal argument of option -acceptTO", 1,
                    Integer.MAX_VALUE));
        if (cl.hasOption("releaseTO"))
            dcmqr.setReleaseTimeout(parseInt(cl.getOptionValue("releaseTO"),
                    "illegal argument of option -releaseTO", 1,
                    Integer.MAX_VALUE));
        if (cl.hasOption("soclosedelay"))
            dcmqr.setSocketCloseDelay(parseInt(cl
                    .getOptionValue("soclosedelay"),
                    "illegal argument of option -soclosedelay", 1, 10000));
        if (cl.hasOption("rcvpdulen"))
            dcmqr.setMaxPDULengthReceive(parseInt(cl
                    .getOptionValue("rcvpdulen"),
                    "illegal argument of option -rcvpdulen", 1, 10000)
                    * KB);
        if (cl.hasOption("sndpdulen"))
            dcmqr.setMaxPDULengthSend(parseInt(cl.getOptionValue("sndpdulen"),
                    "illegal argument of option -sndpdulen", 1, 10000)
                    * KB);
        if (cl.hasOption("sosndbuf"))
            dcmqr.setSendBufferSize(parseInt(cl.getOptionValue("sosndbuf"),
                    "illegal argument of option -sosndbuf", 1, 10000)
                    * KB);
        if (cl.hasOption("sorcvbuf"))
            dcmqr.setReceiveBufferSize(parseInt(cl.getOptionValue("sorcvbuf"),
                    "illegal argument of option -sorcvbuf", 1, 10000)
                    * KB);
        if (cl.hasOption("filebuf"))
            dcmqr.setFileBufferSize(parseInt(cl.getOptionValue("filebuf"),
                    "illegal argument of option -filebuf", 1, 10000)
                    * KB);
        dcmqr.setPackPDV(!cl.hasOption("pdv1"));
        dcmqr.setTcpNoDelay(!cl.hasOption("tcpdelay"));
        dcmqr.setMaxOpsInvoked(cl.hasOption("async") ? parseInt(cl
                .getOptionValue("async"), "illegal argument of option -async",
                0, 0xffff) : 1);
        dcmqr.setMaxOpsPerformed(cl.hasOption("cstoreasync") ? parseInt(cl
                .getOptionValue("cstoreasync"), "illegal argument of option -cstoreasync",
                0, 0xffff) : 0);
        if (cl.hasOption("C"))
            dcmqr.setCancelAfter(parseInt(cl.getOptionValue("C"),
                    "illegal argument of option -C", 1, Integer.MAX_VALUE));
        if (cl.hasOption("lowprior"))
            dcmqr.setPriority(CommandUtils.LOW);
        if (cl.hasOption("highprior"))
            dcmqr.setPriority(CommandUtils.HIGH);
        if (cl.hasOption("cstore")) {
            String[] storeTCs = cl.getOptionValues("cstore");
            for (String storeTC : storeTCs) {
                String cuid;
                String[] tsuids;
                int colon = storeTC.indexOf(':');
                if (colon == -1) {
                    cuid = storeTC;
                    tsuids = DEF_TS;
                } else {
                    cuid = storeTC.substring(0, colon);
                    String ts = storeTC.substring(colon+1);
                    try {
                        tsuids = TS.valueOf(ts).uids;
                    } catch (IllegalArgumentException e) {
                        tsuids = ts.split(",");
                    }
                }
                try {
                    cuid = CUID.valueOf(cuid).uid;
                } catch (IllegalArgumentException e) {
                    // assume cuid already contains UID
                }
                dcmqr.addStoreTransferCapability(cuid, tsuids);
            }
            if (cl.hasOption("cstoredest"))
                dcmqr.setStoreDestination(cl.getOptionValue("cstoredest"));
        }
        dcmqr.setCGet(cl.hasOption("cget"));
        if (cl.hasOption("cmove"))
            dcmqr.setMoveDest(cl.getOptionValue("cmove"));
        if (cl.hasOption("evalRetrieveAET"))
            dcmqr.setEvalRetrieveAET(true);
        if (cl.hasOption("P"))
            dcmqr.setQueryLevel(QueryRetrieveLevel.PATIENT);
        else if (cl.hasOption("S"))
            dcmqr.setQueryLevel(QueryRetrieveLevel.SERIES);
        else if (cl.hasOption("I"))
            dcmqr.setQueryLevel(QueryRetrieveLevel.IMAGE);
        else
            dcmqr.setQueryLevel(QueryRetrieveLevel.STUDY);
        if (cl.hasOption("noextneg"))
            dcmqr.setNoExtNegotiation(true);
        if (cl.hasOption("rel"))
            dcmqr.setRelationQR(true);
        if (cl.hasOption("datetime"))
            dcmqr.setDateTimeMatching(true);
        if (cl.hasOption("fuzzy"))
            dcmqr.setFuzzySemanticPersonNameMatching(true);
        if (!cl.hasOption("P")) {
            if (cl.hasOption("retall"))
                dcmqr.addPrivate(
                        UID.PrivateStudyRootQueryRetrieveInformationModelFIND);
            if (cl.hasOption("blocked"))
                dcmqr.addPrivate(
                        UID.PrivateBlockedStudyRootQueryRetrieveInformationModelFIND);
            if (cl.hasOption("vmf"))
                dcmqr.addPrivate(
                        UID.PrivateVirtualMultiframeStudyRootQueryRetrieveInformationModelFIND);
        }
        if (cl.hasOption("q")) {
            String[] matchingKeys = cl.getOptionValues("q");
            for (int i = 1; i < matchingKeys.length; i++, i++)
                dcmqr.addMatchingKey(Tag.toTagPath(matchingKeys[i - 1]), matchingKeys[i]);
        }
        if (cl.hasOption("r")) {
            String[] returnKeys = cl.getOptionValues("r");
            for (int i = 0; i < returnKeys.length; i++)
                dcmqr.addReturnKey(Tag.toTagPath(returnKeys[i]));
        }

        dcmqr.configureTransferCapability(cl.hasOption("ivrle"));
        

        int repeat = cl.hasOption("repeat") ? parseInt(cl
                .getOptionValue("repeat"),
                "illegal argument of option -repeat", 1, Integer.MAX_VALUE) : 0;
        int interval = cl.hasOption("repeatdelay") ? parseInt(cl
                .getOptionValue("repeatdelay"),
                "illegal argument of option -repeatdelay", 1, Integer.MAX_VALUE)
                : 0;
        boolean closeAssoc = cl.hasOption("closeassoc");

        if (cl.hasOption("tls")) {
            String cipher = cl.getOptionValue("tls");
            if ("NULL".equalsIgnoreCase(cipher)) {
                dcmqr.setTlsWithoutEncyrption();
            } else if ("3DES".equalsIgnoreCase(cipher)) {
                dcmqr.setTls3DES_EDE_CBC();
            } else if ("AES".equalsIgnoreCase(cipher)) {
                dcmqr.setTlsAES_128_CBC();
            } else {
                exit("Invalid parameter for option -tls: " + cipher);
            }
            if (cl.hasOption("nossl2")) {
                dcmqr.disableSSLv2Hello();
            }
            dcmqr.setTlsNeedClientAuth(!cl.hasOption("noclientauth"));
            if (cl.hasOption("keystore")) {
                dcmqr.setKeyStoreURL(cl.getOptionValue("keystore"));
            }
            if (cl.hasOption("keystorepw")) {
                dcmqr.setKeyStorePassword(
                        cl.getOptionValue("keystorepw"));
            }
            if (cl.hasOption("keypw")) {
                dcmqr.setKeyPassword(cl.getOptionValue("keypw"));
            }
            if (cl.hasOption("truststore")) {
                dcmqr.setTrustStoreURL(
                        cl.getOptionValue("truststore"));
            }
            if (cl.hasOption("truststorepw")) {
                dcmqr.setTrustStorePassword(
                        cl.getOptionValue("truststorepw"));
            }
            long t1 = System.currentTimeMillis();
            try {
                dcmqr.initTLS();
            } catch (Exception e) {
                System.err.println("ERROR: Failed to initialize TLS context:"
                        + e.getMessage());
                System.exit(2);
            }
            long t2 = System.currentTimeMillis();
            LOG.info("Initialize TLS context in {} s", Float.valueOf((t2 - t1) / 1000f));
        }
        try {
            dcmqr.start();
        } catch (Exception e) {
            System.err.println("ERROR: Failed to start server for receiving " +
                    "Storage Commitment results:" + e.getMessage());
            System.exit(2);
        }
        List<DicomObject> result;
        try {
            long t1 = System.currentTimeMillis();
            try {
                dcmqr.open();
            } catch (Exception e) {
                LOG.error("Failed to establish association:", e);
                System.exit(2);
            }
            long t2 = System.currentTimeMillis();
            LOG.info("Connected to {} in {} s", remoteAE, Float.valueOf((t2 - t1) / 1000f));
    
            for (;;) {
                try {
                	//************************************************************************************
                    result = dcmqr.query();
                    //********************************************************************************
                    long t3 = System.currentTimeMillis();
                    LOG.info("Received {} matching entries in {} s", Integer.valueOf(result.size()),
                            Float.valueOf((t3 - t2) / 1000f));
                    if (dcmqr.isCMove() || dcmqr.isCGet()) {
                        if (dcmqr.isCMove())
                            dcmqr.move(result);
                        else
                            dcmqr.get(result);
                        long t4 = System.currentTimeMillis();
                        LOG.info("Retrieved {} objects (warning: {}, failed: {}) in {}s",
                                        new Object[] {
                                                Integer.valueOf(dcmqr
                                                        .getTotalRetrieved()),
                                                Integer.valueOf(dcmqr.getWarning()),
                                                Integer.valueOf(dcmqr.getFailed()),
                                                Float.valueOf((t4 - t3) / 1000f) });
                    }
                    if (repeat == 0 || closeAssoc) {
                        try {
                            dcmqr.close();
                        } catch (InterruptedException e) {
                            LOG.error(e.getMessage(), e);
                        }
                        LOG.info("Released connection to {}",remoteAE);
                    }
                    if (repeat-- == 0)
                        break;
                    Thread.sleep(interval);
                    long t4 = System.currentTimeMillis();
                    dcmqr.open();
                    t2 = System.currentTimeMillis();
                    LOG.info("Reconnect or reuse connection to {} in {} s",
                            remoteAE, Float.valueOf((t2 - t4) / 1000f));
                } catch (IOException e) {
                    LOG.error(e.getMessage(), e);
                } catch (InterruptedException e) {
                    LOG.error(e.getMessage(), e);
                } catch (ConfigurationException e) {
                    LOG.error(e.getMessage(), e);
                }
            }
        } finally {
            dcmqr.stop();
        }
        return result;
    }

    private void addStoreTransferCapability(String cuid, String[] tsuids) {
        storeTransferCapability.add(new TransferCapability(
                cuid, tsuids, TransferCapability.SCP));
    }

    private void setEvalRetrieveAET(boolean evalRetrieveAET) {
       this.evalRetrieveAET = evalRetrieveAET;
    }

    private boolean isEvalRetrieveAET() {
        return evalRetrieveAET;
    }

    private void setNoExtNegotiation(boolean b) {
        this.noExtNegotiation = b;
    }

    private void setFuzzySemanticPersonNameMatching(boolean b) {
        this.fuzzySemanticPersonNameMatching = b;
    }

    private void setDateTimeMatching(boolean b) {
        this.dateTimeMatching = b;
    }

    private void setRelationQR(boolean b) {
        this.relationQR = b;
    }

    public final int getFailed() {
        return failed;
    }

    public final int getWarning() {
        return warning;
    }

    private final int getTotalRetrieved() {
        return completed + warning;
    }

    private void setCancelAfter(int limit) {
        this.cancelAfter = limit;
    }

    private void addMatchingKey(int[] tagPath, String value) {
        keys.putString(tagPath, null, value);
    }

    private void addReturnKey(int[] tagPath) {
        keys.putNull(tagPath, null);
    }

    private void configureTransferCapability(boolean ivrle) {
        String[] findcuids = qrlevel.getFindClassUids();
        String[] movecuids = moveDest != null ? qrlevel.getMoveClassUids()
                : EMPTY_STRING;
        String[] getcuids = cget ? qrlevel.getGetClassUids()
                : EMPTY_STRING;
        TransferCapability[] tcs = new TransferCapability[findcuids.length
                + privateFind.size() + movecuids.length + getcuids.length
                + storeTransferCapability.size()];
        int i = 0;
        for (String cuid : findcuids)
            tcs[i++] = mkFindTC(cuid, ivrle ? IVRLE_TS : NATIVE_LE_TS);
        for (String cuid : privateFind)
            tcs[i++] = mkFindTC(cuid, ivrle ? IVRLE_TS : DEFLATED_TS);
        for (String cuid : movecuids)
            tcs[i++] = mkRetrieveTC(cuid, ivrle ? IVRLE_TS : NATIVE_LE_TS);
        for (String cuid : getcuids)
            tcs[i++] = mkRetrieveTC(cuid, ivrle ? IVRLE_TS : NATIVE_LE_TS);
        for (TransferCapability tc : storeTransferCapability) {
            tcs[i++] = tc;
        }
        ae.setTransferCapability(tcs);
        if (!storeTransferCapability.isEmpty()) {
            ae.register(createStorageService());
        }
    }

    private DicomService createStorageService() {
        String[] cuids = new String[storeTransferCapability.size()];
        int i = 0;
        for (TransferCapability tc : storeTransferCapability) {
            cuids[i++] = tc.getSopClass();
        }
        return new StorageService(cuids) {
            @Override
            protected void onCStoreRQ(Association as, int pcid, DicomObject rq,
                    PDVInputStream dataStream, String tsuid, DicomObject rsp)
                    throws IOException, DicomServiceException {
                if (storeDest == null) {
                    super.onCStoreRQ(as, pcid, rq, dataStream, tsuid, rsp);
                } else {
                    try {
                        String cuid = rq.getString(Tag.AffectedSOPClassUID);
                        String iuid = rq.getString(Tag.AffectedSOPInstanceUID);
                        BasicDicomObject fmi = new BasicDicomObject();
                        fmi.initFileMetaInformation(cuid, iuid, tsuid);
                        File file = devnull ? storeDest : new File(storeDest, iuid);
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
            
        };
    }

    private TransferCapability mkRetrieveTC(String cuid, String[] ts) {
        ExtRetrieveTransferCapability tc = new ExtRetrieveTransferCapability(
                cuid, ts, TransferCapability.SCU);
        tc.setExtInfoBoolean(
                ExtRetrieveTransferCapability.RELATIONAL_RETRIEVAL, relationQR);
        if (noExtNegotiation)
            tc.setExtInfo(null);
        return tc;
    }

    private TransferCapability mkFindTC(String cuid, String[] ts) {
        ExtQueryTransferCapability tc = new ExtQueryTransferCapability(cuid,
                ts, TransferCapability.SCU);
        tc.setExtInfoBoolean(ExtQueryTransferCapability.RELATIONAL_QUERIES,
                relationQR);
        tc.setExtInfoBoolean(ExtQueryTransferCapability.DATE_TIME_MATCHING,
                dateTimeMatching);
        tc.setExtInfoBoolean(ExtQueryTransferCapability.FUZZY_SEMANTIC_PN_MATCHING,
                fuzzySemanticPersonNameMatching);
        if (noExtNegotiation)
            tc.setExtInfo(null);
        return tc;
    }

    private void setQueryLevel(QueryRetrieveLevel qrlevel) {
        this.qrlevel = qrlevel;
        keys.putString(Tag.QueryRetrieveLevel, VR.CS, qrlevel.getCode());
        for (int tag : qrlevel.getReturnKeys()) {
            keys.putNull(tag, null);
        }
    }

    public final void addPrivate(String cuid) {
        privateFind.add(cuid);
    }

    private void setCGet(boolean cget) {
        this.cget = cget;
    }

    private boolean isCGet() {
        return cget;
    }

    private void setMoveDest(String aet) {
        moveDest = aet;
    }

    private boolean isCMove() {
        return moveDest != null;
    }

    private static int toPort(String port) {
        return port != null ? parseInt(port, "illegal port number", 1, 0xffff)
                : 104;
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

    private static String[] split(String s, char delim) {
        String[] s2 = { s, null };
        int pos = s.indexOf(delim);
        if (pos != -1) {
            s2[0] = s.substring(0, pos);
            s2[1] = s.substring(pos + 1);
        }
        return s2;
    }

    private static void exit(String msg) {
        System.err.println(msg);
        System.err.println("Try 'dcmqr -h' for more information.");
        System.exit(1);
    }

    public void start() throws IOException { 
        if (conn.isListening()) {
            conn.bind(executor );
            System.out.println("Start Server listening on port " + conn.getPort());
        }
    }

    public void stop() {
        if (conn.isListening()) {
            conn.unbind();
        }
    }

    public void open() throws IOException, ConfigurationException,
            InterruptedException {
        assoc = ae.connect(remoteAE, executor);
    }

    public List<DicomObject> query() throws IOException, InterruptedException {
        List<DicomObject> result = new ArrayList<DicomObject>();
        TransferCapability tc = selectFindTransferCapability();
        String cuid = tc.getSopClass();
        String tsuid = selectTransferSyntax(tc);
        if (tc.getExtInfoBoolean(ExtQueryTransferCapability.RELATIONAL_QUERIES)
                || containsUpperLevelUIDs(cuid)) {
            LOG.info("Send Query Request using {}:\n{}",
                    UIDDictionary.getDictionary().prompt(cuid), keys);
            DimseRSP rsp = assoc.cfind(cuid, priority, keys, tsuid, cancelAfter);
            while (rsp.next()) {
                DicomObject cmd = rsp.getCommand();
                if (CommandUtils.isPending(cmd)) {
                    DicomObject data = rsp.getDataset();
                    result.add(data);
                    LOG.info("Query Response #{}:\n{}", Integer.valueOf(result.size()), data);
                }
            }
        } else {
            List<DicomObject> upperLevelUIDs = queryUpperLevelUIDs(cuid, tsuid);
            List<DimseRSP> rspList = new ArrayList<DimseRSP>(upperLevelUIDs.size());
            for (int i = 0, n = upperLevelUIDs.size(); i < n; i++) {
                upperLevelUIDs.get(i).copyTo(keys);
                LOG.info("Send Query Request #{}/{} using {}:\n{}",
                        new Object[] {
                            Integer.valueOf(i+1),
                            Integer.valueOf(n),
                            UIDDictionary.getDictionary().prompt(cuid),
                            keys
                        });
                rspList.add(assoc.cfind(cuid, priority, keys, tsuid, cancelAfter));
            }
            for (int i = 0, n = rspList.size(); i < n; i++) {
                DimseRSP rsp = rspList.get(i);
                for (int j = 0; rsp.next(); ++j) {
                    DicomObject cmd = rsp.getCommand();
                    if (CommandUtils.isPending(cmd)) {
                        DicomObject data = rsp.getDataset();
                        result.add(data);
                        LOG.info("Query Response #{} for Query Request #{}/{}:\n{}",
                                new Object[]{ Integer.valueOf(j+1), Integer.valueOf(i+1), Integer.valueOf(n), data });
                    }
                }
            }
        }
        return result;
    }

    @SuppressWarnings("fallthrough")
    private boolean containsUpperLevelUIDs(String cuid) {
        switch (qrlevel) {
        case IMAGE:
            if (!keys.containsValue(Tag.SeriesInstanceUID)) {
                return false;
            }
            // fall through
        case SERIES:
            if (!keys.containsValue(Tag.StudyInstanceUID)) {
                return false;
            }
            // fall through
        case STUDY:
            if (Arrays.asList(PATIENT_LEVEL_FIND_CUID).contains(cuid)
                    && !keys.containsValue(Tag.PatientID)) {
                return false;
            }
            // fall through
        case PATIENT:
            // fall through
        }
        return true;
    }

    private List<DicomObject> queryUpperLevelUIDs(String cuid, String tsuid)
            throws IOException, InterruptedException {
        List<DicomObject> keylist = new ArrayList<DicomObject>();
        if (Arrays.asList(PATIENT_LEVEL_FIND_CUID).contains(cuid)) {
            queryPatientIDs(cuid, tsuid, keylist);
            if (qrlevel == QueryRetrieveLevel.STUDY) {
                return keylist;
            }
            keylist = queryStudyOrSeriesIUIDs(cuid, tsuid, keylist,
                    Tag.StudyInstanceUID, STUDY_MATCHING_KEYS, QueryRetrieveLevel.STUDY);
        } else {
            keylist.add(new BasicDicomObject());
            keylist = queryStudyOrSeriesIUIDs(cuid, tsuid, keylist,
                    Tag.StudyInstanceUID, PATIENT_STUDY_MATCHING_KEYS, QueryRetrieveLevel.STUDY);
        }
        if (qrlevel == QueryRetrieveLevel.IMAGE) {
            keylist = queryStudyOrSeriesIUIDs(cuid, tsuid, keylist,
                    Tag.SeriesInstanceUID, SERIES_MATCHING_KEYS, QueryRetrieveLevel.SERIES);
        }
        return keylist;
    }

    private void queryPatientIDs(String cuid, String tsuid,
            List<DicomObject> keylist) throws IOException, InterruptedException {
        String patID = keys.getString(Tag.PatientID);
        String issuer = keys.getString(Tag.IssuerOfPatientID);
        if (patID != null) {
            DicomObject patIdKeys = new BasicDicomObject();
            patIdKeys.putString(Tag.PatientID, VR.LO, patID);
            if (issuer != null) {
                patIdKeys.putString(Tag.IssuerOfPatientID, VR.LO, issuer);
            }
            keylist.add(patIdKeys);
        } else {
            DicomObject patLevelQuery = new BasicDicomObject();
            keys.subSet(PATIENT_MATCHING_KEYS).copyTo(patLevelQuery);
            patLevelQuery.putNull(Tag.PatientID, VR.LO);
            patLevelQuery.putNull(Tag.IssuerOfPatientID, VR.LO);
            patLevelQuery.putString(Tag.QueryRetrieveLevel, VR.CS, "PATIENT");
            LOG.info("Send Query Request using {}:\n{}",
                    UIDDictionary.getDictionary().prompt(cuid), patLevelQuery);
            DimseRSP rsp = assoc.cfind(cuid, priority, patLevelQuery, tsuid,
                    Integer.MAX_VALUE);
            for (int i = 0; rsp.next(); ++i) {
                DicomObject cmd = rsp.getCommand();
                if (CommandUtils.isPending(cmd)) {
                    DicomObject data = rsp.getDataset();
                    LOG.info("Query Response #{}:\n{}", Integer.valueOf(i+1), data);
                    DicomObject patIdKeys = new BasicDicomObject();
                    patIdKeys.putString(Tag.PatientID, VR.LO,
                            data.getString(Tag.PatientID));
                    issuer = keys.getString(Tag.IssuerOfPatientID);
                    if (issuer != null) {
                        patIdKeys.putString(Tag.IssuerOfPatientID, VR.LO,
                                issuer);
                    }
                    keylist.add(patIdKeys);
                }
            }
        }
    }

    private List<DicomObject> queryStudyOrSeriesIUIDs(String cuid, String tsuid,
            List<DicomObject> upperLevelIDs, int uidTag, int[] matchingKeys,
            QueryRetrieveLevel qrLevel) throws IOException,
            InterruptedException {
        List<DicomObject> keylist = new ArrayList<DicomObject>();
        String uid = keys.getString(uidTag);
        for (DicomObject upperLevelID : upperLevelIDs) {
            if (uid != null) {
                DicomObject suidKey = new BasicDicomObject();
                upperLevelID.copyTo(suidKey);
                suidKey.putString(uidTag, VR.UI, uid);
                keylist.add(suidKey);
            } else {
                DicomObject keys2 = new BasicDicomObject();
                keys.subSet(matchingKeys).copyTo(keys2);
                upperLevelID.copyTo(keys2);
                keys2.putNull(uidTag, VR.UI);
                keys2.putString(Tag.QueryRetrieveLevel, VR.CS, qrLevel.getCode());
                LOG.info("Send Query Request using {}:\n{}",
                        UIDDictionary.getDictionary().prompt(cuid), keys2);
                DimseRSP rsp = assoc.cfind(cuid, priority, keys2,
                        tsuid, Integer.MAX_VALUE);
                for (int i = 0; rsp.next(); ++i) {
                    DicomObject cmd = rsp.getCommand();
                    if (CommandUtils.isPending(cmd)) {
                        DicomObject data = rsp.getDataset();
                        LOG.info("Query Response #{}:\n{}", Integer.valueOf(i+1), data);
                        DicomObject suidKey = new BasicDicomObject();
                        upperLevelID.copyTo(suidKey);
                        suidKey.putString(uidTag, VR.UI, data.getString(uidTag));
                        keylist.add(suidKey);
                    }
                }
            }
        }
        return keylist;
    }

    private TransferCapability selectFindTransferCapability()
            throws NoPresentationContextException {
        TransferCapability tc;
        if ((tc = selectTransferCapability(privateFind)) != null)
            return tc;
        if ((tc = selectTransferCapability(qrlevel.getFindClassUids())) != null)
            return tc;
        throw new NoPresentationContextException(UIDDictionary.getDictionary()
                .prompt(qrlevel.getFindClassUids()[0])
                + " not supported by" + remoteAE.getAETitle());
    }

    private String selectTransferSyntax(TransferCapability tc) {
        String[] tcuids = tc.getTransferSyntax();
        if (Arrays.asList(tcuids).indexOf(UID.DeflatedExplicitVRLittleEndian) != -1)
            return UID.DeflatedExplicitVRLittleEndian;
        return tcuids[0];
    }

    public void move(List<DicomObject> findResults)
            throws IOException, InterruptedException {
        if (moveDest == null)
            throw new IllegalStateException("moveDest == null");
        TransferCapability tc = selectTransferCapability(qrlevel.getMoveClassUids());
        if (tc == null)
            throw new NoPresentationContextException(UIDDictionary
                    .getDictionary().prompt(qrlevel.getMoveClassUids()[0])
                    + " not supported by" + remoteAE.getAETitle());
        String cuid = tc.getSopClass();
        String tsuid = selectTransferSyntax(tc);
        for (int i = 0, n = Math.min(findResults.size(), cancelAfter); i < n; ++i) {
            DicomObject keys = findResults.get(i).subSet(MOVE_KEYS);
            if (isEvalRetrieveAET()
                    && moveDest.equals(findResults.get(i)
                            .getString(Tag.RetrieveAETitle))) {
                LOG.info("Skipping {}:\n{}",
                        UIDDictionary.getDictionary().prompt(cuid), keys);
            } else {
                LOG.info("Send Retrieve Request using {}:\n{}",
                        UIDDictionary.getDictionary().prompt(cuid), keys);
                assoc.cmove(cuid, priority, keys, tsuid, moveDest, rspHandler);
            }
        }
        assoc.waitForDimseRSP();
    }


    public void get(List<DicomObject> findResults)
            throws IOException, InterruptedException {
        TransferCapability tc = selectTransferCapability(qrlevel.getGetClassUids());
        if (tc == null)
            throw new NoPresentationContextException(UIDDictionary
                    .getDictionary().prompt(qrlevel.getGetClassUids()[0])
                    + " not supported by" + remoteAE.getAETitle());
        String cuid = tc.getSopClass();
        String tsuid = selectTransferSyntax(tc);
        for (int i = 0, n = Math.min(findResults.size(), cancelAfter); i < n; ++i) {
            DicomObject keys = findResults.get(i).subSet(MOVE_KEYS);
            LOG.info("Send Retrieve Request using {}:\n{}",
                    UIDDictionary.getDictionary().prompt(cuid), keys);
            assoc.cget(cuid, priority, keys, tsuid, rspHandler);
        }
        assoc.waitForDimseRSP();
    }
    
    private final DimseRSPHandler rspHandler = new DimseRSPHandler() {
        @Override
        public void onDimseRSP(Association as, DicomObject cmd,
                DicomObject data) {
            DcmQR.this.onMoveRSP(as, cmd, data);
        }
    };

    protected void onMoveRSP(Association as, DicomObject cmd, DicomObject data) {
        if (!CommandUtils.isPending(cmd)) {
            completed += cmd.getInt(Tag.NumberOfCompletedSuboperations);
            warning += cmd.getInt(Tag.NumberOfWarningSuboperations);
            failed += cmd.getInt(Tag.NumberOfFailedSuboperations);
        }

    }

    private TransferCapability selectTransferCapability(String[] cuid) {
        TransferCapability tc;
        for (int i = 0; i < cuid.length; i++) {
            tc = assoc.getTransferCapabilityAsSCU(cuid[i]);
            if (tc != null)
                return tc;
        }
        return null;
    }

    private TransferCapability selectTransferCapability(List<String> cuid) {
        TransferCapability tc;
        for (int i = 0, n = cuid.size(); i < n; i++) {
            tc = assoc.getTransferCapabilityAsSCU(cuid.get(i));
            if (tc != null)
                return tc;
        }
        return null;
    }

    public void close() throws InterruptedException {
        assoc.release(true);
    }

    private void setStoreDestination(String filePath) {
        this.storeDest = new File(filePath);
        this.devnull = "/dev/null".equals(filePath);
        if (!devnull)
            storeDest.mkdir();
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
            return DcmQR.class.getClassLoader().getResourceAsStream(
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
}
