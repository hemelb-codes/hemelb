package uk.ac.ucl.chem.ccs.clinicalgui.res;
import uk.ac.ucl.chem.ccs.aheclient.res.*;
import uk.ac.ucl.chem.ccs.aheclient.util.AHEJobObject;
import uk.ac.ucl.chem.ccs.aheclient.util.PrepareResponse;
import uk.ac.ucl.chem.ccs.aheclient.util.ResourceElement;
import uk.ac.ucl.chem.ccs.clinicalgui.*;
import info.clearthought.layout.TableLayout;
import java.awt.BorderLayout;
import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.ComboBoxModel;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import javax.swing.WindowConstants;
import java.awt.Container;
import javax.swing.table.TableModel;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import edu.lsu.cct.cosched.clientAPI.*;
import java.util.Vector;
import java.util.Hashtable;
import java.util.Iterator;
import javax.swing.table.AbstractTableModel;

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
public class SelectResource extends javax.swing.JDialog {
	private JPanel holder;
	private JButton cancelButton;
	private JButton finishButton;
	private JTable resourceTable;
	private JScrollPane jScrollPane1;
	private VLLaunch v;
	String patientId;
	String dateTime;
	DisplayJobPanel djp;

	/**
	* Auto-generated main method to display this JFrame
	*/
	public static void main(String[] args) {
	
		SelectResource inst = new SelectResource();
		inst.setVisible(true);
	}

	
	public SelectResource() {
		super();
		initGUI();
		}
	
	public SelectResource(Container parent,String patientId, String studyId, String seriesId, String dateTime, DisplayJobPanel djp) {
		super();
		this.djp = djp;
		this.patientId = patientId;
		this.dateTime = dateTime;
		v = new VLLaunch(GridServerInterface.getParamsPath(), GridServerInterface.getRootPath() + "/" + patientId + "/" + studyId + "/" + seriesId + "/" + dateTime + "/pars.asc",
				GridServerInterface.getRootPath() + "/" + patientId + "/" + studyId + "/" + seriesId + "/" + dateTime + "/config.dat");
		initGUI();
		this.setLocationRelativeTo(parent);
		}	
	
	public void showDialog() {
		this.setModal(true);
		this.setVisible(true);
		//return new AdvancedReservation();
	}
	
	private void initGUI() {
		try {
			setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
			this.setTitle("Select Resource");
				holder = new JPanel();
				getContentPane().add(holder, BorderLayout.CENTER);
				TableLayout holderLayout = new TableLayout(new double[][] {
						{ TableLayout.FILL, TableLayout.FILL, TableLayout.FILL,
								TableLayout.FILL },
						{ TableLayout.FILL, TableLayout.FILL, TableLayout.FILL,
								TableLayout.FILL, TableLayout.FILL,
								TableLayout.FILL, TableLayout.FILL,
								TableLayout.FILL, TableLayout.FILL,
								TableLayout.FILL, TableLayout.FILL,
								TableLayout.FILL } });
				holderLayout.setHGap(5);
				holderLayout.setVGap(5);
				holder.setLayout(holderLayout);
				holder.setPreferredSize(new java.awt.Dimension(410, 525));
			
							finishButton = new JButton();
							holder.add(finishButton, "2, 10, 3, 11");
							finishButton.setText("Launch Hemelb");
							finishButton.addActionListener(new ActionListener() {
								public void actionPerformed(ActionEvent evt) {
									finish();
								}
							});
					
							cancelButton = new JButton();
							holder.add(cancelButton, "0, 10, 1, 11");
							cancelButton.setText("Cancel");
							cancelButton
								.addActionListener(new ActionListener() {
									public void actionPerformed(ActionEvent evt) {
										cleanup();
									}
								});
		
		           
								//Build a table model with the RMs in
								CoallocatorFactory.loadProperties(true);
								try {
									Coallocator co = CoallocatorFactory
										.getCoallocator();
								} catch (Exception e) {
									//foo
								}
								//Vector allComp = Resource
									//.getKnownResourcesOfType(SimpleComputeResource.class);
								PrepareResponse pr = v.prepare(patientId + "|" + dateTime);
								Vector allComp = pr.getResources();
						
								Object[][] rms = new Object[allComp.size()][4];
								int x = 0;

								//Add compute resource
								Iterator it = allComp.iterator();
								while (it.hasNext()) {
									rms[x][0] = new Boolean(false);
									rms[x][1] = ((ResourceElement) it.next()).getCommonName();
									rms[x][2] = "Compute";
									rms[x][3] = new Integer(0);
									x++;
								}
        
								TableModel resourceTableModel = new MyTableModel(
									rms,
									new String[] { "Select", "Resource",
											"Type", "Procs" });
								
								
								
								resourceTable = new JTable(resourceTableModel);
								jScrollPane1 = new JScrollPane();
								jScrollPane1.setViewportView(resourceTable);
								holder.add(jScrollPane1, "0, 0, 3, 9");
								

					resourceTable.getColumnModel().getColumn(0).setPreferredWidth(5);
					resourceTable.getColumnModel().getColumn(1).setPreferredWidth(150);
					resourceTable.getColumnModel().getColumn(2).setPreferredWidth(5);
					resourceTable.getColumnModel().getColumn(3).setPreferredWidth(5);
				
					resourceTable.setRowHeight(20);
			//}
			pack();
			this.setSize(416, 531);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private void finish() {
		
		String resource = null;
		int proc = -1;
		
		for (int k=0; k < resourceTable.getRowCount(); k++) {
			Boolean bool = (Boolean)resourceTable.getValueAt(k, 0);
			if (bool.booleanValue() == true) {
				resource = (String)resourceTable.getValueAt(k, 1);
				proc = ((Integer)resourceTable.getValueAt(k, 3)).intValue();
				if (proc < 1) {
					ErrorMessage ems = new ErrorMessage(this, "Processor count must be > 0");
					return;
				}
				break;			
			}
		}
		if(resource == null){
			ErrorMessage ems = new ErrorMessage(this, "You need to select a resource");
			return;
		}
		djp.setJobObject(v.start(resource,proc, null));
		System.out.println("reset the simulation monitoring panel");
		cleanup();
	}
	
	private void cleanup() {
		this.dispose();
	}
	
	class MyTableModel extends AbstractTableModel {
	    private String[] columnNames;
	    private Object[][] data;
	    
	    public MyTableModel (Object[][] data, String[] columnNames) {
	    	super();
	    	this.data = data;
	    	this.columnNames = columnNames;
	    }

	    public int getColumnCount() {
	        return columnNames.length;
	    }

	    public int getRowCount() {
	        return data.length;
	    }

	    public String getColumnName(int col) {
	        return columnNames[col];
	    }

	    public Object getValueAt(int row, int col) {
	        return data[row][col];
	    }

	    public Class getColumnClass(int c) {
	        return getValueAt(0, c).getClass();
	    }

	    /*
	     * Don't need to implement this method unless your table's
	     * editable.
	     */
	    public boolean isCellEditable(int row, int col) {
	        //Note that the data/cell address is constant,
	        //no matter where the cell appears onscreen.
	        if (col == 0 || col == 3) {
	            return true;
	        } else {
	            return false;
	        }
	    }

	    /*
	     * Don't need to implement this method unless your table's
	     * data can change.
	     */
	    public void setValueAt(Object value, int row, int col) {
	        data[row][col] = value;
	        fireTableCellUpdated(row, col);
	    }
	    
	}
}
