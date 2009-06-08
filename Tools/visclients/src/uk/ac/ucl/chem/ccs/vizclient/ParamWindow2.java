package uk.ac.ucl.chem.ccs.vizclient;
import info.clearthought.layout.TableLayout;
import javax.swing.WindowConstants;
//import org.jdesktop.application.Application;
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
public class ParamWindow2 extends javax.swing.JDialog {
	private JLabel jLabel1;
	private JPanel jPanel1;
	public JTextField pixels_xField;
	public JTextField pixels_yField;
	public JButton updateButton;
	private JLabel jLabel2;
	private SteeringData sd;
	
	/**
	* Auto-generated main method to display this JFrame
	*/
	public static void main(String[] args) {
	
		ParamWindow2 inst = new ParamWindow2();
		inst.setVisible(true);
	}

	
	public ParamWindow2() {
		super();
		//this.sd = sd;
		initGUI();
		}
	
	public ParamWindow2(Container parent) {
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
			this.setTitle("Resolution");
			this.setResizable(false);
			this.setUndecorated(false);
			{
				jPanel1 = new JPanel();
				TableLayout jPanel1Layout = new TableLayout(new double[][] {{TableLayout.FILL}, {TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}});
				jPanel1Layout.setHGap(5);
				jPanel1Layout.setVGap(5);
				jPanel1.setLayout(jPanel1Layout);
				getContentPane().add(jPanel1, "0, 0");
				jPanel1.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
				{
					jLabel1 = new JLabel();
					jPanel1.add(jLabel1, "0, 0");
					jLabel1.setText("Render width (px)");
				}
				{
					jLabel2 = new JLabel();
					jPanel1.add(jLabel2, "0, 2");
					jLabel2.setText("Render height (px)");
				}
				{
					updateButton = new JButton();
					jPanel1.add(updateButton, "0, 4");
					updateButton.setText("Update");
				}
				{
					pixels_xField = new JTextField();
					jPanel1.add(pixels_xField, "0, 1");
				}
				{
					pixels_yField = new JTextField();
					jPanel1.add(pixels_yField, "0, 3");
				}
			}

			pack();
			this.setSize(150, 172);
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
		pixels_yField.setText(Integer.toString(sd.getPixels_y()));
		pixels_xField.setText(Integer.toString(sd.getPixels_x()));
	}
	
	

	private void cleanup() {
		this.dispose();
	}
	
}
