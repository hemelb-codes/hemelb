package uk.ac.ucl.chem.ccs.clinicalgui;

import info.clearthought.layout.TableLayout;
import info.clearthought.layout.TableLayoutConstants;

import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Enumeration;
import java.util.Hashtable;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.WindowConstants;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.TreeSelectionModel;

import org.jdesktop.application.Application;

import uk.ac.ucl.chem.ccs.clinicalgui.res.SelectResource;

/**
 * @author Konstantin Voevodski
 * 
 *         implements the UI compoenents and functionality of the "Simulation Launch" panel
 */
public class SimulationLaunchPanel extends javax.swing.JPanel
{
  /**
   * 
   */
  private static final long serialVersionUID = 1L;
  private JScrollPane       treeView;
  // contains the tree object representing the model data available on the grid ftp server
  private JEditorPane       infoLabel;
  // displays the status of the GUI (error messages/status)
  private JEditorPane       noteLabel;
  // displays the
  private JButton           LaunchHemelbButton;
  // button to launch hemelb (on model selected in the tree)
  public JButton            PrevButton;
  // go back to prev panel
  public JButton            NextButton;
  // go forward to next panel
  private JPanel            reservationInfoPanel;
  // panel which contains reservation info (if a reservation is selected in ResPanel
  // *****
  // components of reservationPanel:
  private JLabel            resHeader;
  private JLabel            resName;
  private JLabel            resStart;
  private JLabel            resEnd;
  private JLabel            m1;
  private JLabel            m2;
  private JLabel            m3;
  private JLabel            m4;
  // ********
  private String            selectedPath;
  // the selected path in the model tree
  private JPanel            reservationTypePanel;
  // lets the user select between reservation types (have reservation, don't have reservation,
  // SPRUCE)
  // ****
  // components of reservationTypePanel:
  private JLabel            resTypeHeader;
  private JRadioButton      b1;
  private JRadioButton      b2;
  private JRadioButton      b3;
  private ButtonGroup       bGroup;
  // *************
  private DisplayJobPanel   companionPanel;
  // the next panel (simulation monitoring tab: DisplayJobPanel)
  private MainPanel         parentPanel;

  // the parent panel (needed to advance to next panel when HemeLB is launched)

  /**
   * Auto-generated main method to display this JPanel inside a new JFrame.
   */
  public static void main(String[] args)
  {
    JFrame frame = new JFrame();
    frame.getContentPane().add(new ModelGenerationPanel());
    frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
    frame.pack();
    frame.setVisible(true);
  }

  public SimulationLaunchPanel()
  {
    super();
    initGUI(null);
  }

  public void setParentPanel(MainPanel p)
  {
    parentPanel = p;
  }

  public void setCompanionPanel(DisplayJobPanel djp)
  {
    companionPanel = djp;
  }

  private void updateInfo(String m)
  {
    infoLabel.setText(m);
  }

  public void createDataTree(String[] highlightedPath)
  {
    try
    {
      // DefaultMutableTreeNode top = GridServerInterface.populateServerData();
      JTree tree = GridServerInterface.populateServerData(highlightedPath);
      // tree.setExpandsSelectedPaths(true);
      tree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
      tree.addTreeSelectionListener(this.new myTreeListener());
      treeView = new JScrollPane(tree);
      this.add(treeView, "0, 0, 1, 3");
    } catch (Exception ex)
    {
      ex.printStackTrace();
      updateInfo("could not load data from the grid");
    }
  }

  public class myTreeListener implements TreeSelectionListener
  {
    public void valueChanged(TreeSelectionEvent e)
    {
      selectedPath = e.getPath().toString();
      String input = selectedPath.substring(1, selectedPath.length() - 1);
      String[] data = input.split(",");
      // System.out.println("length of path is " + data.length);
      if (data.length != 6)
      {
        return;
      }
      else if (!data[5].trim().equals("notes.txt"))
      {
        return;
      }
      else
      {
        try
        {
          noteLabel.setText("Note: "
              + GridServerInterface.getModelNote(data[1].trim(), data[2].trim(), data[3].trim(),
                  data[4].trim()));
        } catch (Exception ex)
        {
          updateInfo("retrieval of note unsuccessful");
        }
        updateInfo("successfully retrieved note");
      }
    }
  }

  public void updateReservationInfo()
  {
    if (ClinicalGuiClient.reservation != null)
    {
      b1.setSelected(true);
      resName.setText("Name: " + ClinicalGuiClient.reservation.getResName());
      resStart.setText("Start: " + ClinicalGuiClient.reservation.getStart().toString());
      resEnd.setText("End: " + ClinicalGuiClient.reservation.getEnd().toString());
      m1.setText("Resources:");
      Hashtable res = ClinicalGuiClient.reservation.getRes();
      Enumeration myEnum = res.keys();
      if (myEnum.hasMoreElements())
      {
        m2.setText((String) myEnum.nextElement());
      }
      if (myEnum.hasMoreElements())
      {
        m3.setText((String) myEnum.nextElement());
      }
      if (res.size() > 3)
      {
        m4.setText("...");
      }
      else if (myEnum.hasMoreElements())
      {
        m4.setText((String) myEnum.nextElement());
        // ClinicalGuiClient.reservation.getRes().
      }
    }
    else
    {
      b2.setSelected(true);
    }
  }

  public void initGUI(String[] highlightedPath)
  {
    try
    {

      TableLayout thisLayout = new TableLayout(new double[][] {
          { TableLayoutConstants.FILL, TableLayoutConstants.FILL, TableLayoutConstants.FILL,
              TableLayoutConstants.FILL },
          { TableLayoutConstants.FILL, TableLayoutConstants.FILL, TableLayoutConstants.FILL,
              TableLayoutConstants.FILL, TableLayoutConstants.FILL, TableLayoutConstants.FILL } });
      thisLayout.setHGap(5);
      thisLayout.setVGap(5);
      this.setLayout(thisLayout);
      setPreferredSize(new Dimension(400, 300));

      infoLabel = new JEditorPane();
      this.add(infoLabel, "0, 4, 1, 4");
      infoLabel.setText("select a data set above");
      infoLabel.setBorder(BorderFactory.createEtchedBorder());

      createDataTree(highlightedPath);

      noteLabel = new JEditorPane();
      this.add(noteLabel, "2, 0, 3, 1");
      noteLabel.setBorder(BorderFactory.createEtchedBorder());

      reservationInfoPanel = new JPanel();
      this.add(reservationInfoPanel, "2, 2, 3, 3");
      reservationInfoPanel.setLayout(new TableLayout(new double[][] {
          { TableLayoutConstants.FILL },
          { TableLayoutConstants.FILL, TableLayoutConstants.FILL, TableLayoutConstants.FILL,
              TableLayoutConstants.FILL, TableLayoutConstants.FILL, TableLayoutConstants.FILL,
              TableLayoutConstants.FILL, TableLayoutConstants.FILL } }));
      reservationInfoPanel.setBorder(BorderFactory.createEtchedBorder());

      resHeader = new JLabel();
      reservationInfoPanel.add(resHeader, "0,0");
      resHeader.setText("Reservation Details:");

      resName = new JLabel();
      reservationInfoPanel.add(resName, "0,1");

      resStart = new JLabel();
      reservationInfoPanel.add(resStart, "0,2");

      resEnd = new JLabel();
      reservationInfoPanel.add(resEnd, "0,3");

      m1 = new JLabel();
      reservationInfoPanel.add(m1, "0,4");

      m2 = new JLabel();
      reservationInfoPanel.add(m2, "0,5");

      m3 = new JLabel();
      reservationInfoPanel.add(m3, "0,6");

      m4 = new JLabel();
      reservationInfoPanel.add(m4, "0,7");

      reservationTypePanel = new JPanel();
      this.add(reservationTypePanel, "2,4,3,4");
      reservationTypePanel.setLayout(new TableLayout(new double[][] {
          { TableLayoutConstants.FILL, TableLayoutConstants.FILL, TableLayoutConstants.FILL },
          { TableLayoutConstants.FILL, TableLayoutConstants.FILL } }));
      b1 = new JRadioButton("Use res");
      b2 = new JRadioButton("No res");
      b3 = new JRadioButton("SPRUCE");
      bGroup = new ButtonGroup();
      bGroup.add(b1);
      bGroup.add(b2);
      bGroup.add(b3);

      resTypeHeader = new JLabel("Reservation type:");
      reservationTypePanel.add(resTypeHeader, "0,0,1,0");
      reservationTypePanel.add(b1, "0,1");
      reservationTypePanel.add(b2, "1,1");
      reservationTypePanel.add(b3, "2,1");
      reservationTypePanel.setBorder(BorderFactory.createEtchedBorder());

      NextButton = new JButton();
      this.add(NextButton, "3, 5");
      NextButton.setName("NextButton");

      PrevButton = new JButton();
      this.add(PrevButton, "0, 5");
      PrevButton.setName("PrevButton");

      LaunchHemelbButton = new JButton();
      this.add(LaunchHemelbButton, "1, 5, 2, 5");
      LaunchHemelbButton.setName("LaunchHemelbButton");
      LaunchHemelbButton.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent evt)
        {
          if (selectedPath == null)
          {
            updateInfo("you need to select a model first");
            return;
          }
          String input = selectedPath.substring(1, selectedPath.length() - 1);
          String[] data = input.split(",");
          if (data.length < 5)
          {
            updateInfo("model not fully specified");
            return;
          }
          if (b1.isSelected())
          {
            if (ClinicalGuiClient.reservation == null)
            {
              updateInfo("you don't have a reservation");
              return;
            }
            else
            {
              VLLaunch l = new VLLaunch(GridServerInterface.getParamsPath(), GridServerInterface
                  .getRootPath()
                  + "/"
                  + data[1].trim()
                  + "/"
                  + data[2].trim()
                  + "/"
                  + data[3].trim()
                  + "/"
                  + data[4].trim() + "/pars.asc", GridServerInterface.getRootPath() + "/"
                  + data[1].trim() + "/" + data[2].trim() + "/" + data[3].trim() + "/"
                  + data[4].trim() + "/config.dat");
              l.prepare(data[1].trim() + "|" + data[4].trim());
              companionPanel.setJobObject(l.start("", 0, ClinicalGuiClient.reservation));
              parentPanel.advance();
            }
          }
          if (b2.isSelected())
          {
            SelectResource sr = new SelectResource(
                SimulationLaunchPanel.this.getTopLevelAncestor(), data[1].trim(), data[2].trim(),
                data[3].trim(), data[4].trim(), companionPanel);
            sr.showDialog();
            parentPanel.advance();
          }
          else if (b3.isSelected())
          {
            updateInfo("Spruce refers to trees of the genus Picea, a genus of about 35 "
                + "species of coniferous evergreen trees in the Family Pinaceae.");
            return;
          }
          else
          {
            updateInfo("you need to select a reservation type first");
            return;
          }
        }
      });

      updateReservationInfo();
      Application.getInstance().getContext().getResourceMap(getClass()).injectComponents(this);
    } catch (Exception e)
    {
      e.printStackTrace();
    }
  }

}
