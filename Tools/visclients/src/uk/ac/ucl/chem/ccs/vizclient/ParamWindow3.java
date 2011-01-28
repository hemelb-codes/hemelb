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
public class ParamWindow3 extends javax.swing.JDialog {
	private JLabel jLabel1;
	private JPanel jPanel1;
	private JLabel jLabel4;
	private JLabel jLabel6;
	private JLabel jLabel5;
	private JLabel jLabel3;
	private JLabel jLabel2;
	private JLabel jLabel7;
	private JLabel jLabel8;
	public JTextField viz_brightnessField;
	public JTextField velocity_maxField;
	public JTextField pressure_maxField;
	public JTextField pressure_minField;
	public JTextField stress_maxField;
	public JTextField vis_glyph_lengthField;
	public JTextField vis_streaklines_per_pulsatile_periodField;
	public JTextField vis_streakline_lengthField;
	public JButton updateButton;

	private SteeringData sd;
	
	/**
	* Auto-generated main method to display this JFrame
	*/
	public static void main(String[] args) {
	
		ParamWindow3 inst = new ParamWindow3();
		inst.setVisible(true);
	}

	
	public ParamWindow3() {
		super();
		//this.sd = sd;
		initGUI();
		}
	
	public ParamWindow3(Container parent) {
		super();
		//this.sd = sd;
		initGUI();
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
			this.setTitle("Render parameters");
			this.setResizable(false);
			this.setUndecorated(false);
			{
				jPanel1 = new JPanel();
				TableLayout jPanel1Layout = new TableLayout(new double[][] {{TableLayout.FILL}, {TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}});
				jPanel1Layout.setHGap(5);
				jPanel1Layout.setVGap(5);
				jPanel1.setLayout(jPanel1Layout);
				getContentPane().add(jPanel1, "0, 0");
				jPanel1.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
				{
					jLabel1 = new JLabel();
					jPanel1.add(jLabel1, "0, 0");
					jLabel1.setText("Brightness");
				}
				{
					jLabel2 = new JLabel();
					jPanel1.add(jLabel2, "0, 2");
					jLabel2.setText("Max velocity (m/s)");
				}
				{
					jLabel3 = new JLabel();
					jPanel1.add(jLabel3, "0, 6");
					jLabel3.setText("Max pressure (mm.Hg)");
				}
				{
					jLabel4 = new JLabel();
					jPanel1.add(jLabel4, "0, 4");
					jLabel4.setText("Min pressure (mm.Hg)");
				}
				{
					jLabel5 = new JLabel();
					jPanel1.add(jLabel5, "0, 8");
					jLabel5.setText("Max stress (Pa)");
				}
				{
					jLabel6 = new JLabel();
					jPanel1.add(jLabel6, "0, 10");
					jLabel6.setText("Glyph length");
				}
				{
					jLabel7 = new JLabel();
					jPanel1.add(jLabel7, "0, 12");
					jLabel7.setText("Streaklines/period");
				}
				{
					jLabel8 = new JLabel();
					jPanel1.add(jLabel8, "0, 14");
					jLabel8.setText("Streakline length");
				}
				{
					updateButton = new JButton();
					jPanel1.add(updateButton, "0, 16");
					updateButton.setText("Update");
				}
				{
					vis_glyph_lengthField = new JTextField();
					jPanel1.add(vis_glyph_lengthField, "0, 11");
				}
				{
					pressure_maxField = new JTextField();
					jPanel1.add(pressure_maxField, "0, 7");
				}
				{
					pressure_minField = new JTextField();
					jPanel1.add(pressure_minField, "0, 5");
				}
				{
					stress_maxField = new JTextField();
					jPanel1.add(stress_maxField, "0, 9");
				}
				{
					velocity_maxField = new JTextField();
					jPanel1.add(velocity_maxField, "0, 3");
				}
				{
					viz_brightnessField = new JTextField();
					jPanel1.add(viz_brightnessField, "0, 1");
				}
				{
					vis_streaklines_per_pulsatile_periodField = new JTextField();
					jPanel1.add(vis_streaklines_per_pulsatile_periodField, "0, 13");
				}
				{
					vis_streakline_lengthField = new JTextField();
					jPanel1.add(vis_streakline_lengthField, "0, 15");
				}
			}

			pack();
			this.setSize(150, 585);
			//Application.getInstance().getContext().getResourceMap(getClass()).injectComponents(getContentPane());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void update (SteeringData sd) {
		this.sd = sd;
		updateText();
	}
	
	private void updateText () {
		viz_brightnessField.setText(Float.toString(sd.getVis_brightness()));
		velocity_maxField.setText(Float.toString(sd.getVelocity_max()));
		pressure_maxField.setText(Float.toString(sd.getPressure_max()));
		pressure_minField.setText(Float.toString(sd.getPressure_min()));
		stress_maxField.setText(Float.toString(sd.getStress_max()));
		vis_glyph_lengthField.setText(Float.toString(sd.getVis_glyph_length()));
		vis_streaklines_per_pulsatile_periodField.setText(Float.toString(sd.getVis_streaklines_per_pulsatile_period()));
		vis_streakline_lengthField.setText(Float.toString(sd.getVis_streakline_length()));
	}
	
	

	private void cleanup() {
		this.dispose();
	}
	
}
