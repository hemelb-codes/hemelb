package uk.ac.ucl.chem.ccs.clinicalgui.res;
import info.clearthought.layout.TableLayout;
import uk.ac.ucl.chem.ccs.clinicalgui.*;
import uk.ac.ucl.chem.ccs.aheclient.res.*;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.SimpleDateFormat;
import java.util.Hashtable;
import java.util.Vector;
import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import java.util.Iterator;
import java.util.Enumeration;
import javax.swing.WindowConstants;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableModel;
import org.jdesktop.application.Application;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTable;


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
public class ResPanel extends javax.swing.JPanel {
	private JLabel jLabel1;
	private JPanel holder;
	private JButton view;
	public JButton next;
	private JTable resTable;
	private JPanel jPanel1;
	private JButton newRes;
	private JButton delete;
	private JButton edit;
	private JButton select;
	private JRadioButton gur;
	private JRadioButton hard;
	private JRadioButton add;
	private ButtonGroup buttonGroup1;
	private JPanel controlPanel;
	private JScrollPane resList;
	private Vector tableData;
	private String filename;
	private SimulationLaunchPanel simPanel;
	
	
	/**
	* Auto-generated main method to display this 
	* JPanel inside a new JFrame.
	*/

	
	public ResPanel(String filename) {
		super();
		tableData = new Vector();
		this.filename = filename;
		loadData();
		initGUI();
	}
	
	public void setSimPanel(SimulationLaunchPanel p){
		simPanel = p;
	}
	
	private void initGUI() {
		try {
			this.setLayout(null);
			{
				{
				}
				{
					buttonGroup1 = new ButtonGroup();
				}
				holder = new JPanel();
				TableLayout holderLayout = new TableLayout(new double[][] {
						{ TableLayout.FILL, TableLayout.FILL, TableLayout.FILL,
								TableLayout.FILL },
						{ TableLayout.FILL, TableLayout.FILL, TableLayout.FILL,
								TableLayout.FILL, TableLayout.FILL,
								TableLayout.FILL, TableLayout.FILL,
								TableLayout.FILL, TableLayout.FILL,
								TableLayout.FILL, TableLayout.FILL,
								TableLayout.FILL, TableLayout.FILL,
								TableLayout.FILL } });
				holderLayout.setHGap(5);
				holderLayout.setVGap(5);
				holder.setLayout(holderLayout);
				this.add(holder);
				holder.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
				holder.setBounds(0, 0, 553, 322);
				holder.setSize(630, 613);
				{
					jLabel1 = new JLabel();
					holder.add(jLabel1, "0,   0,   1,   0");
					TableLayout jLabel1Layout = new TableLayout(new double[][] {
							{ TableLayout.FILL, TableLayout.FILL,
									TableLayout.FILL, TableLayout.FILL },
							{ TableLayout.FILL, TableLayout.FILL,
									TableLayout.FILL, TableLayout.FILL } });
					jLabel1Layout.setHGap(5);
					jLabel1Layout.setVGap(5);
					jLabel1.setLayout(null);
					jLabel1.setText("Advanced Reservations");
					jLabel1.setFont(new java.awt.Font("SansSerif",0,18));
					jLabel1.setOpaque(true);
				}
				{
					controlPanel = new JPanel();
					TableLayout controlPanelLayout = new TableLayout(
						new double[][] {
								{ TableLayout.FILL, TableLayout.FILL,
										TableLayout.FILL, TableLayout.FILL },
								{ TableLayout.FILL, TableLayout.FILL,
										TableLayout.FILL, TableLayout.FILL } });
					controlPanelLayout.setHGap(5);
					controlPanelLayout.setVGap(5);
					controlPanel.setLayout(controlPanelLayout);
					holder.add(controlPanel, "0, 11, 3, 13");
					{
						add = new JRadioButton();
						controlPanel.add(add, "0, 0");
						add.setText("All");
						add.setBorder(BorderFactory.createEmptyBorder(0, 10, 0, 0));
						buttonGroup1.add(add);
						add.setSelected(true);
					}
					{
						hard = new JRadioButton();
						controlPanel.add(hard, "0, 1");
						hard.setText("HARC");
						hard.setBorder(BorderFactory.createEmptyBorder(0, 10, 0, 0));
						buttonGroup1.add(hard);
					}
					{
						gur = new JRadioButton();
						controlPanel.add(gur, "0, 2");
						gur.setText("GUR");
						gur.setEnabled(false);
						gur.setBorder(BorderFactory.createEmptyBorder(0, 10, 0, 0));
						buttonGroup1.add(gur);
					}
					{
						view = new JButton();
						controlPanel.add(view, "2, 0");
						view.setText("View");
						view.addActionListener(new ActionListener() {
							public void actionPerformed (ActionEvent evt) {
								viewDetails();
							}
						});
						
					}
					    select = new JButton();
					    controlPanel.add(select, "2,1");
					    select.setText("Select");
					    select.addActionListener(new ActionListener() {
							public void actionPerformed (ActionEvent evt) {
								int selectedRow = resTable.getSelectedRow();
								System.out.println("Selected row is " + selectedRow);
								if (selectedRow > -1){
									System.out.println("So this is called");
									AdvancedReservation selected = (AdvancedReservation)tableData.get(selectedRow);
									ClinicalGuiClient.reservation = selected;
									simPanel.updateReservationInfo();
								}
							}
						});
					    
					{
						edit = new JButton();
						controlPanel.add(edit, "3, 0");
						edit.setText("Edit");
						edit.setEnabled(false);
					}
					{
						delete = new JButton();
						controlPanel.add(delete, "3, 1");
						delete.setText("Delete");
						delete.addActionListener(new ActionListener() {
							public void actionPerformed (ActionEvent evt) {
								deleteAction();
							}
						});
					}
					{
						newRes = new JButton();
						controlPanel.add(newRes, "1, 0");
						newRes.setText("Create");
						newRes.addActionListener(new ActionListener() {
							public void actionPerformed (ActionEvent evt) {
								CreateReservation newRes = new CreateReservation(ResPanel.this.getTopLevelAncestor());
								newRes.showDialog();
								AdvancedReservation res = newRes.getReservation();
								if (res != null){
									tableData.add(res);
									updateTable();
									saveData();
								}
							}
						});
					}
					{
						next = new JButton();
						controlPanel.add(next, "3, 3");
						next.setName("next");
					}
				}
				{
					jPanel1 = new JPanel();
					TableLayout jPanel1Layout = new TableLayout(new double[][] {
							{ TableLayout.FILL }, { TableLayout.FILL } });
					jPanel1Layout.setHGap(5);
					jPanel1Layout.setVGap(5);
					jPanel1.setLayout(jPanel1Layout);
					holder.add(jPanel1, "0, 1, 3, 10");
					jPanel1.setBackground(new java.awt.Color(255,255,255));
					jPanel1.setForeground(new java.awt.Color(255,255,255));
					{
						resList = new JScrollPane();
						jPanel1.add(resList, "0,  0");
						resList
							.setBackground(new java.awt.Color(255, 255, 255));
						resList
							.setForeground(new java.awt.Color(255, 255, 255));
						resList.setOpaque(false);
						{
							
							resTable = new JTable();
							resList.setViewportView(resTable);
							resTable.setFocusable(false);
							resTable.setGridColor(new java.awt.Color(133,114,114));
							resTable.setIntercellSpacing(new java.awt.Dimension(0, 1));
							updateTable();
							

							
						}
					}
				}
			}
			this.setPreferredSize(new java.awt.Dimension(630, 623));
			this.setSize(630, 623);
			Application.getInstance().getContext().getResourceMap(getClass()).injectComponents(this);
		} catch (Exception e) {
			e.printStackTrace();
		}
		

	}

	private void updateTable() {
		Object [][] data = new Object[tableData.size()][5];
		
		Iterator it = tableData.iterator();
		int x = 0;
		while (it.hasNext()) {
			SimpleDateFormat df = new SimpleDateFormat ("yyyy/MM/dd HH:mm z");
			AdvancedReservation ar = (AdvancedReservation)it.next();
			data[x][0] = ar.getResName();
			data[x][1] = df.format(ar.getStart());
			data[x][2] = df.format(ar.getEnd());
			if (ar.getType() == AdvancedReservation.HARC){
			data[x][3] = "HARC";
			} else {
				data[x][3] = "GUR";
			}
			data[x][4] = new Integer(ar.getRes().size());
			x++;
		}
		
		TableModel resTableModel = new MyTableModel(
				data,
				new String[] { "Name", "Start", "End", "Type", "Resources" });
		resTable.setModel(resTableModel);
		resTable.getColumnModel().getColumn(0).setPreferredWidth(10);
		resTable.getColumnModel().getColumn(1).setPreferredWidth(150);
		resTable.getColumnModel().getColumn(2).setPreferredWidth(150);
		resTable.getColumnModel().getColumn(3).setPreferredWidth(5);
		resTable.getColumnModel().getColumn(4).setPreferredWidth(5);
		resTable.setRowHeight(20);
	}
	
	private void loadData() {
		tableData = ReservationPersistance.read(filename);
	}

	private void saveData () {
		ReservationPersistance.write(filename, tableData);
	}
	
	private void clickedNext () {
		int selectedRow = resTable.getSelectedRow();
		System.out.println("Selected row is " + selectedRow);
		if (selectedRow > -1){
			System.out.println("So this is called");
				ClinicalGuiClient.reservation = (AdvancedReservation)tableData.get(selectedRow);
		}
	}
	
	private void viewDetails () {
		int selectedRow = resTable.getSelectedRow();
		System.out.println("Selected row is " + selectedRow);
		if (selectedRow > -1){
			System.out.println("So this is called");
			AdvancedReservation selected = (AdvancedReservation)tableData.get(selectedRow);
			ViewReservation vr = new ViewReservation(ResPanel.this.getTopLevelAncestor(), selected);
			vr.showDialog();
		}
	}
	
	private void deleteAction() {
		HARCManager hm = new HARCManager();
		int selectedRows [] = resTable.getSelectedRows();
		if (selectedRows.length > 0){
			boolean ok = false;
			int g = selectedRows.length-1;
			System.out.println(g);
			while (g >= 0) {
				AdvancedReservation ar = (AdvancedReservation)tableData.get(selectedRows[g]);
				ok = hm.delete(ar);
				if (ok) {
					tableData.remove(selectedRows[g]);
				} else {
					ErrorMessage ems = new ErrorMessage(ResPanel.this.getTopLevelAncestor(), hm.getError());
				}
				g--;
			}
			saveData();
			updateTable();
		}
		
		
	}

	class MyTableModel extends AbstractTableModel {
		 protected boolean isSortAsc = true;
		 protected int sortCol = 0;
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

