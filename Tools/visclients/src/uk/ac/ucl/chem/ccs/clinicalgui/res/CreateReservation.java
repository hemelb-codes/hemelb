package uk.ac.ucl.chem.ccs.clinicalgui.res;
import uk.ac.ucl.chem.ccs.aheclient.res.*;
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
public class CreateReservation extends javax.swing.JDialog {
	private JPanel holder;
	private JPanel harcContainer;
	private JPanel Controls;
	private JTextField nameField;
	private JLabel jLabel2;
	private JPanel middlePanel;
	private JButton cancelButton;
	private JButton finishButton;
	private JPanel buttonPanel;
	private JTable resourceTable;
	private JScrollPane jScrollPane1;
	private JTextField durationField;
	private JTextField timeField;
	private JComboBox yearCombo;
	private JComboBox monthComb;
	private JComboBox dayCombo;
	private JLabel duration;
	private JLabel time;
	private JLabel date;
	private ButtonGroup buttonGroup1;
	private JRadioButton gur;
	private JRadioButton harc;
	private JLabel jLabel1;
	private Calendar today;  
	private AdvancedReservation adv; 

	/**
	* Auto-generated main method to display this JFrame
	*/
	public static void main(String[] args) {
	
		CreateReservation inst = new CreateReservation();
		inst.setVisible(true);
	}

	
	public CreateReservation() {
		super();
		adv = null;
		today = Calendar.getInstance();
		initGUI();
		}
	
	public CreateReservation(Container parent) {
		super();
		today = Calendar.getInstance();
		initGUI();
		this.setLocationRelativeTo(parent);
		}	
	
	public AdvancedReservation getReservation() {
		return adv;
	}
	
	public AdvancedReservation showDialog() {
		this.setModal(true);
		this.setVisible(true);
		return new AdvancedReservation();
	}
	
	private void initGUI() {
		try {
			setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
			this.setTitle("Create New Reservation");
			{
				{
					buttonGroup1 = new ButtonGroup();
				}
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
				{
					jLabel1 = new JLabel();
					holder.add(jLabel1, "0, 0");
					jLabel1.setText("Type");
					jLabel1.setBorder(BorderFactory.createEmptyBorder(
						0,
						10,
						0,
						0));
				}
				{
					harc = new JRadioButton();
					holder.add(harc, "1, 0");
					harc.setText("HARC");
					buttonGroup1.add(harc);
					harc.setSelected(true);
				}
				{
					gur = new JRadioButton();
					holder.add(gur, "2, 0");
					gur.setText("GUR");
					buttonGroup1.add(gur);
					gur.setEnabled(false);
				}
				{
					harcContainer = new JPanel();
					TableLayout harcContainerLayout = new TableLayout(
						new double[][] {
								{ TableLayout.FILL, TableLayout.FILL,
										TableLayout.FILL, TableLayout.FILL },
								{ TableLayout.FILL, TableLayout.FILL,
										TableLayout.FILL, TableLayout.FILL,
										TableLayout.FILL, TableLayout.FILL,
										TableLayout.FILL, TableLayout.FILL,
										TableLayout.FILL, TableLayout.FILL,
										TableLayout.FILL, TableLayout.FILL } });
					harcContainerLayout.setHGap(5);
					harcContainerLayout.setVGap(5);
					harcContainer.setLayout(harcContainerLayout);
					holder.add(harcContainer, "0, 1, 3, 11");
					{
						Controls = new JPanel();
						TableLayout ControlsLayout = new TableLayout(
							new double[][] {
									{ TableLayout.FILL, TableLayout.FILL,
											TableLayout.FILL, TableLayout.FILL,
											TableLayout.FILL, TableLayout.FILL,
											TableLayout.FILL, TableLayout.FILL },
									{ TableLayout.FILL, TableLayout.FILL,
											TableLayout.FILL } });
						ControlsLayout.setHGap(5);
						ControlsLayout.setVGap(5);
						Controls.setLayout(ControlsLayout);
						harcContainer.add(Controls, "0, 0, 3, 2");
						Controls.setBorder(BorderFactory.createEmptyBorder(
							10,
							10,
							10,
							10));
						{
							time = new JLabel();
							Controls.add(time, "0, 1");
							time.setText("Time");
						}
						{
							date = new JLabel();
							Controls.add(date, "0, 0");
							date.setText("Date");
						}
						{
							duration = new JLabel();
							Controls.add(duration, "4, 1, 5, 1");
							duration.setText("Duration");
						}
						{
							String[] days = new String[] { "1", "2", "3", "4",
									"5", "6", "7", "8", "9", "10", "11", "12",
									"13", "14", "15", "16", "17", "18", "19",
									"20", "21", "22", "23", "24", "25", "26",
									"27", "28", "29", "30", "31" };
							ComboBoxModel dayComboModel = new DefaultComboBoxModel(
								days);
							dayCombo = new JComboBox();
							Controls.add(dayCombo, "2, 0");
							dayComboModel.setSelectedItem(days[today
								.get(Calendar.DAY_OF_MONTH) - 1]);
							dayCombo.setModel(dayComboModel);
							dayCombo.setEditable(true);

						}
						{
							String[] months = new String[] { "January",
									"February", "March", "April", "May",
									"June", "July", "August", "September",
									"October", "November", "December" };
							ComboBoxModel monthCombModel = new DefaultComboBoxModel(
								months);
							monthComb = new JComboBox();
							Controls.add(monthComb, "3, 0, 5, 0");
							monthCombModel.setSelectedItem(months[today
								.get(Calendar.MONTH)]);
							monthComb.setModel(monthCombModel);
							monthComb.setBackground(new java.awt.Color(255,255,255));
						}
						{
							int year = today.get(Calendar.YEAR);
							String[] years = new String[10];
							for (int i = 0; i < 10; i++) {
								years[i] = Integer.toString(year + i);
							}
							ComboBoxModel yearComboModel = new DefaultComboBoxModel(
								years);
							yearCombo = new JComboBox();
							Controls.add(yearCombo, "6, 0, 7, 0");
							yearCombo.setModel(yearComboModel);
							yearCombo.setBackground(new java.awt.Color(255,255,255));
						}
						{
							timeField = new JTextField();
							Controls.add(timeField, "2, 1, 3, 1");
							String minute = new String();
							int min = today.get(Calendar.MINUTE);
							if (min < 10) {
								minute = "0" + Integer.toString(min);
							} else {
								minute = Integer.toString(min);
							}
							timeField.setText(Integer.toString(today
								.get(Calendar.HOUR_OF_DAY))
								+ ":"
								+ minute);
						}
						{
							durationField = new JTextField();
							Controls.add(durationField, "6, 1, 7, 1");
							durationField.setText("00:00");
						}
						{
							jLabel2 = new JLabel();
							Controls.add(jLabel2, "0, 2");
							jLabel2.setText("Name");
						}
						{
							nameField = new JTextField();
							Controls.add(nameField, "2, 2, 7, 2");
						}
					}
					{
						buttonPanel = new JPanel();
						TableLayout buttonPanelLayout = new TableLayout(
							new double[][] {
									{ TableLayout.FILL, TableLayout.FILL,
											TableLayout.FILL },
									{ TableLayout.FILL, TableLayout.FILL } });
						buttonPanelLayout.setHGap(5);
						buttonPanelLayout.setVGap(5);
						buttonPanel.setLayout(buttonPanelLayout);
						harcContainer.add(buttonPanel, "0, 10, 3, 11");
						buttonPanel.setBorder(BorderFactory.createEmptyBorder(
							20,
							10,
							20,
							10));
						{
							finishButton = new JButton();
							buttonPanel.add(finishButton, "2, 0, 2, 1");
							finishButton.setText("Finish");
							finishButton.addActionListener(new ActionListener() {
								public void actionPerformed(ActionEvent evt) {
									finish();
								}
							});
						}
						{
							cancelButton = new JButton();
							buttonPanel.add(cancelButton, "1, 0, 1, 1");
							cancelButton.setText("Cancel");
							cancelButton
								.addActionListener(new ActionListener() {
									public void actionPerformed(ActionEvent evt) {
										cleanup();
									}
								});
						}
					}
					{
						middlePanel = new JPanel();
						TableLayout middlePanelLayout = new TableLayout(
							new double[][] { { TableLayout.FILL },
									{ TableLayout.FILL } });
						middlePanelLayout.setHGap(5);
						middlePanelLayout.setVGap(5);
						middlePanel.setLayout(middlePanelLayout);
						harcContainer.add(middlePanel, "0, 3, 3, 9");
						middlePanel.setBorder(BorderFactory.createEmptyBorder(
							0,
							10,
							0,
							10));
						{
							jScrollPane1 = new JScrollPane();
							middlePanel.add(jScrollPane1, "0,  0");
							{
								//Build a table model with the RMs in
								CoallocatorFactory.loadProperties(true);
								try {
									Coallocator co = CoallocatorFactory
										.getCoallocator();
								} catch (Exception e) {
									//foo
								}
								Vector allComp = Resource
									.getKnownResourcesOfType(SimpleComputeResource.class);
					

								Object[][] rms = new Object[allComp.size()][4];
								int x = 0;

								//Add compute resource
								Iterator it = allComp.iterator();
								while (it.hasNext()) {
									rms[x][0] = new Boolean(false);
									rms[x][1] = ((SimpleComputeResource) it
										.next()).name();
									rms[x][2] = "Compute";
									rms[x][3] = new Integer(0);
									x++;
								}

								TableModel resourceTableModel = new MyTableModel(
									rms,
									new String[] { "Select", "Resource",
											"Type", "Procs" });
								
								
								
								
								resourceTable = new JTable(resourceTableModel);
								jScrollPane1.setViewportView(resourceTable);
								

					resourceTable.getColumnModel().getColumn(0).setPreferredWidth(5);
					resourceTable.getColumnModel().getColumn(1).setPreferredWidth(150);
					resourceTable.getColumnModel().getColumn(2).setPreferredWidth(5);
					resourceTable.getColumnModel().getColumn(3).setPreferredWidth(5);
				
					resourceTable.setRowHeight(20);

							}
						}
					}
				}
			}
			pack();
			this.setSize(416, 531);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private void finish() {
		boolean worked = true;
		ReservationManager rm = new HARCManager();
		String day=Integer.toString(dayCombo.getSelectedIndex() + 1);
		String month=Integer.toString(monthComb.getSelectedIndex() + 1);
		String year = (String)yearCombo.getSelectedItem();
		String startTime = timeField.getText();
		String duration = durationField.getText();
		
		
		Date start = null;
		DateFormat df = new SimpleDateFormat("dd/MM/yyyy HH:mm");
		try {
			System.err.println("....");
			start = df.parse(day + "/" + month + "/" + year + " " + startTime);

		} catch (Exception e) {
			ErrorMessage ems = new ErrorMessage(this, "Error parsing reservation start time");
			worked = false;
		}

		
		String[] durationComponents = duration.split(":");
		
		if (durationComponents.length < 2) {
			ErrorMessage ems = new ErrorMessage(this, "Error parsing reservation duration");
			worked = false;
		}
		
		 //work out end time
		 int durHour = 0, durMin =0;
		try {
		durHour = Integer.parseInt(durationComponents[0]);
		 durMin = Integer.parseInt(durationComponents[1]);
		} catch (Exception e) {
			ErrorMessage ems = new ErrorMessage(this, "Error parsing reservation duration");
			worked = false;
		}
		
		Calendar dur = Calendar.getInstance();
		dur.setTime(start);
		dur.add(Calendar.HOUR, durHour);
		dur.add(Calendar.MINUTE, durMin);
		
		Date end = dur.getTime();

		
		//check dates
		if (end.before(start)) {
			ErrorMessage ems = new ErrorMessage(this, "End date must be after start date");
			worked = false;
		}
		
		 if (worked) {
		//parse table selections
		Hashtable rms = new Hashtable();

		
		for (int k=0; k < resourceTable.getRowCount(); k++) {
			Boolean bool = (Boolean)resourceTable.getValueAt(k, 0);
			if (bool.booleanValue() == true) {
				String resource = (String)resourceTable.getValueAt(k, 1);
				int proc = ((Integer)resourceTable.getValueAt(k, 3)).intValue();
				if (proc < 1) {
					ErrorMessage ems = new ErrorMessage(this, "Processor count must be > 0");
					return;
				}
				
				ReservationElement re = new ReservationElement(resource,
						proc);
				rms.put(resource, re);				
			}
		}
		
		adv = rm.create(start, end, rms);
		if (adv == null) {
			ErrorMessage ems = new ErrorMessage(this, rm.getError());
		} else {
			adv.setResName(nameField.getText());
		}

		 }

	if (worked) {
		cleanup();
	}
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
