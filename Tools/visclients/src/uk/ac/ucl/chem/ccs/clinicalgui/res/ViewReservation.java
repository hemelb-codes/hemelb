package uk.ac.ucl.chem.ccs.clinicalgui.res;
import info.clearthought.layout.TableLayout;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.DateFormat;
import java.util.Calendar;
import java.util.Enumeration;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.JTextPane;
import javax.swing.WindowConstants;
import javax.swing.border.BevelBorder;
import javax.swing.border.TitledBorder;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;

import uk.ac.ucl.chem.ccs.aheclient.gui.ErrorMessage;
import uk.ac.ucl.chem.ccs.aheclient.res.AdvancedReservation;
import uk.ac.ucl.chem.ccs.aheclient.res.HARCManager;
import uk.ac.ucl.chem.ccs.aheclient.res.ReservationElement;

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
public class ViewReservation extends javax.swing.JDialog {
	private JPanel holder;
	private JPanel harcContainer;
	private JPanel Controls;
	private JTextField dateField;
	private JTextField nameField;
	private JLabel jLabel2;
	private JPanel middlePanel;
	private JButton statusButton;
	private JButton closeButton;
	private JPanel buttonPanel;
	private JTable resourceTable;
	private JScrollPane jScrollPane1;
	private JTextField timeField;
	private JLabel time;
	private JTextPane statusPane;
	private JScrollPane jScrollPane2;
	private JPanel statusPanel;
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
	
		ViewReservation inst = new ViewReservation(new AdvancedReservation());
		inst.setVisible(true);
	}

	
	public ViewReservation(AdvancedReservation adv) {
		super();
		this.adv = adv;
		initGUI();
		}
	
	public ViewReservation(Container parent, AdvancedReservation adv) {
		super();
		this.adv = adv;
		initGUI();
		this.setLocationRelativeTo(parent);
		}	
		
	public boolean showDialog() {
		this.setModal(true);
		this.setVisible(true);
		return true;
	}
	
	private void initGUI() {
		try {
			setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
			this.setTitle(adv.getResName());
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
						statusPanel = new JPanel();
						harcContainer.add(statusPanel, "0, 7, 3, 9");
						TableLayout statusPanelLayout = new TableLayout(
							new double[][] { { TableLayout.FILL },
									{ TableLayout.FILL } });
						statusPanelLayout.setHGap(5);
						statusPanelLayout.setVGap(5);
						statusPanel.setLayout(statusPanelLayout);
						statusPanel.setBorder(BorderFactory.createEmptyBorder(0, 10, 0, 10));
						{
							jScrollPane2 = new JScrollPane();
							statusPanel.add(jScrollPane2, "0, 0");
							{
								statusPane = new JTextPane();
								jScrollPane2.setViewportView(statusPane);
								statusPane.setEditable(false);
								statusPane.setText("Status...");
							}
						}
					}
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
							time.setText("End");
						}
						{
							date = new JLabel();
							Controls.add(date, "0, 0");
							date.setText("Start");
						}
						{
							timeField = new JTextField();
							Controls.add(timeField, "2, 1, 7, 1");
							timeField.setEditable(false);
							timeField.setBackground(java.awt.Color.WHITE);
							timeField.setBorder(BorderFactory.createEtchedBorder(BevelBorder.LOWERED));
							timeField.setText(adv.getEnd().toString());
						}
						{
							jLabel2 = new JLabel();
							Controls.add(jLabel2, "0, 2");
							jLabel2.setText("Name");
						}
						{
							nameField = new JTextField();
							nameField.setBackground(java.awt.Color.WHITE);
							Controls.add(nameField, "2, 2, 7, 2");
							nameField.setEditable(false);
							nameField.setText(adv.getResName());
							nameField.setBorder(BorderFactory.createEtchedBorder(BevelBorder.LOWERED));
						}
						{
							dateField = new JTextField();
							Controls.add(dateField, "2, 0, 7, 0");
							dateField.setEditable(false);
							dateField.setBackground(java.awt.Color.WHITE);
							dateField.setBorder(BorderFactory.createEtchedBorder(BevelBorder.LOWERED));
							dateField.setText(adv.getStart().toString());
							
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
							closeButton = new JButton();
							buttonPanel.add(closeButton, "2, 0, 2, 1");
							closeButton.setText("Close");
							closeButton.addActionListener(new ActionListener() {
								public void actionPerformed(ActionEvent evt) {
									cleanup();
								}
							});
						}
						{
							statusButton = new JButton();
							buttonPanel.add(statusButton, "1, 0, 1, 1");
							statusButton.setText("Get Status");
							statusButton
								.addActionListener(new ActionListener() {
									public void actionPerformed(ActionEvent evt) {
										updateStatus();
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
						harcContainer.add(middlePanel, "0, 3, 3, 6");
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
								int size = adv.getRes().size();
								int y = 0;
								Enumeration en = adv.getRes().elements();
								Object[][] rms = new Object [size][4];
	
								while (en.hasMoreElements()) {
									ReservationElement re = (ReservationElement)en.nextElement();
									rms[y][0] = re.getRmName();
									rms[y][1] = re.getResID();
									if (re.getType() == ReservationElement.COMP) {
									rms[y][2] = "Computer";
									} else {
										rms[y][2] = "Network";
									}
									rms[y][3] = new Integer(re.getProcCount());
									y++;
								}

								TableModel resourceTableModel = new MyTableModel(
									rms,
									new String[] {"Resource","Res ID", 
											"Type", "Procs" });
								
								resourceTable = new JTable(resourceTableModel);
								jScrollPane1.setViewportView(resourceTable);
								

					resourceTable.getColumnModel().getColumn(0).setPreferredWidth(100);
					resourceTable.getColumnModel().getColumn(1).setPreferredWidth(50);
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

	private void updateStatus() {
		HARCManager hm = new HARCManager();
		String status = hm.status(adv);
		if (status != null){
			statusPane.setText(status);
		} else {
			ErrorMessage ems = new ErrorMessage(this, hm.getError());
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
	            return false;
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
