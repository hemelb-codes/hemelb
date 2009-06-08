package uk.ac.ucl.chem.ccs.vizclient;
import info.clearthought.layout.TableLayout;
import javax.swing.WindowConstants;
import java.awt.Container;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

/**
* This code was edited or generated using CloudGarden's Jigloo
* SWT/Swing GUI Builder, which is free for non-commercial
* use. If Jigloo is being used commercially (ie, by a corporation,
* company or business for any purpose whatever) then you
* should purchase a license for each developer using Jigloo.
* Please visit www.cloudgarden.com for details.
* Use of Jigloo implies acceptance of these licensing terms.
* A COMMERCIAL LICENSE HAS NOT BEEN PURCHASED FOR
* THIS MACHINE, SO JIGLOO OR THIS CODE CANNOT BE USED
* LEGALLY FOR ANY CORPORATE OR COMMERCIAL PURPOSE.
*/
public class ParamWindow1 extends javax.swing.JDialog {
	private JLabel jLabel1;
	private JPanel jPanel1;
	private JLabel jLabel4;
	public JTextField zoomField;
	public JTextField latitudeField;
	public JTextField longtitudeField;
	public JTextField vis_ctr_zField;
	public JTextField vis_ctr_yField;
	public JTextField vis_ctr_xField;
	public JButton updateButton;
	private JLabel jLabel6;
	private JLabel jLabel5;
	private JLabel jLabel3;
	private JLabel jLabel2;
	private SteeringData sd;
	
	/**
	* Auto-generated main method to display this JFrame
	*/
	public static void main(String[] args) {
	
		ParamWindow1 inst = new ParamWindow1();
		//inst.setVisible(true);
	}

	
	public ParamWindow1() {
		super();
		//this.sd = sd;
		initGUI();
		}
	
	public ParamWindow1(Container parent) {
		super();
		//this.sd = sd;
		initGUI();
//		int x = parent.getX() + parent.getWidth();
//		//int y = parent.getLocationOnScreen().y;
//		this.setLocation(x,y);
//	 if (parent == null){
//		 System.err.println("null");
//	 }
//	 System.err.println(x);

		this.setLocationRelativeTo(parent);
		}	
		
	public boolean showDialog() {
		this.setModal(false);
		this.setVisible(true);
		return true;
	}
	
	private void initGUI() {
		try {
			setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
			TableLayout thisLayout = new TableLayout(new double[][] {{TableLayout.FILL}, {TableLayout.FILL}});
			thisLayout.setHGap(5);
			thisLayout.setVGap(5);
			getContentPane().setLayout(thisLayout);
			this.setTitle("Scene parameters");
			this.setResizable(false);
			this.setUndecorated(false);
			{
				jPanel1 = new JPanel();
				TableLayout jPanel1Layout = new TableLayout(new double[][] {{TableLayout.FILL}, {TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}});
				jPanel1Layout.setHGap(5);
				jPanel1Layout.setVGap(5);
				jPanel1.setLayout(jPanel1Layout);
				getContentPane().add(jPanel1, "0, 0");
				jPanel1.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
				{
					jLabel1 = new JLabel();
					jPanel1.add(jLabel1, "0, 0");
					jLabel1.setText("Scene centre x");
				}
				{
					jLabel2 = new JLabel();
					jPanel1.add(jLabel2, "0, 2");
					jLabel2.setText("Scene centre y");
				}
				{
					jLabel3 = new JLabel();
					jPanel1.add(jLabel3, "0, 4");
					jLabel3.setText("Scene centre z");
				}
				{
					jLabel4 = new JLabel();
					jPanel1.add(jLabel4, "0, 6");
					jLabel4.setText("Longitude (deg)");
				}
				{
					jLabel5 = new JLabel();
					jPanel1.add(jLabel5, "0, 8");
					jLabel5.setText("Latitude (deg)");
				}
				{
					jLabel6 = new JLabel();
					jPanel1.add(jLabel6, "0, 10");
					jLabel6.setText("Zoom");
				}
				{
					updateButton = new JButton();
					jPanel1.add(updateButton, "0, 12");
					updateButton.setText("Update");
				}
				{
					vis_ctr_xField = new JTextField();
					jPanel1.add(vis_ctr_xField, "0, 1");
				}
				{
					vis_ctr_yField = new JTextField();
					jPanel1.add(vis_ctr_yField, "0, 3");
				}
				{
					vis_ctr_zField = new JTextField();
					jPanel1.add(vis_ctr_zField, "0, 5");
				}
				{
					longtitudeField = new JTextField();
					jPanel1.add(longtitudeField, "0, 7");
				}
				{
					latitudeField = new JTextField();
					jPanel1.add(latitudeField, "0, 9");
				}
				{
					zoomField = new JTextField();
					jPanel1.add(zoomField, "0, 11");
				}
			}

			pack();
			this.setSize(150, 447);
		//	Application.getInstance().getContext().getResourceMap(getClass()).injectComponents(getContentPane());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void update (SteeringData sd) {
		this.sd = sd;
		updateText();
	}
	
	private void updateText () {
		zoomField.setText(Float.toString(sd.getZoom_factor()));
		latitudeField.setText(Float.toString(sd.getLatitude()));
		longtitudeField.setText(Float.toString(sd.getLongitude()));
		vis_ctr_zField.setText(Float.toString(sd.getCtr_z()));
		vis_ctr_yField.setText(Float.toString(sd.getCtr_y()));
		vis_ctr_xField.setText(Float.toString(sd.getCtr_x()));
	}
	
	public void updateLon (SteeringData sd) {
		this.sd = sd;
		longtitudeField.setText(Float.toString(sd.getLongitude()));
	}

	private void cleanup() {
		this.dispose();
	}
	
}
