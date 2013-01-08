

import javax.swing.JFrame;
import javax.swing.JPanel;
import java.awt.BorderLayout;
import java.awt.Dimension;
import javax.swing.JTabbedPane;
import javax.swing.JDesktopPane;
import javax.swing.JButton;
import java.awt.GridBagLayout;
import java.awt.GridBagConstraints;
import javax.swing.JInternalFrame;
import javax.swing.JSplitPane;
import javax.swing.JTextField;
import javax.swing.JLabel;
import javax.swing.JCheckBox;
public class classSetup {

	private JFrame jFrame = null;  //  @jve:decl-index=0:visual-constraint="86,19"
	private JPanel jContentPane = null;
	private JPanel jPanel = null;
	private JTabbedPane jTabbedPane = null;
	private JPanel Results = null;
	private JPanel Movie = null;
	/**
	 * This method initializes jFrame	
	 * 	
	 * @return javax.swing.JFrame	
	 */
	private JFrame getJFrame() {
		if (jFrame == null) {
			jFrame = new JFrame();
			jFrame.setSize(new Dimension(723, 261));
			jFrame.setTitle("testing");
			jFrame.setContentPane(getJContentPane());
		}
		return jFrame;
	}

	/**
	 * This method initializes jContentPane	
	 * 	
	 * @return javax.swing.JPanel	
	 */
	private JPanel getJContentPane() {
		if (jContentPane == null) {
			jContentPane = new JPanel();
			jContentPane.setLayout(new BorderLayout());
			jContentPane.add(getJPanel(), BorderLayout.CENTER);
			
		}
		return jContentPane;
	}

	/**
	 * This method initializes jPanel	
	 * 	
	 * @return javax.swing.JPanel	
	 */
	private JPanel getJPanel() {
		if (jPanel == null) {
			GridBagConstraints gridBagConstraints9 = new GridBagConstraints();
			gridBagConstraints9.fill = GridBagConstraints.BOTH;
			gridBagConstraints9.gridy = 0;
			gridBagConstraints9.weightx = 1.0;
			gridBagConstraints9.weighty = 1.0;
			gridBagConstraints9.gridx = 0;
			jPanel = new JPanel();
			jPanel.setLayout(new GridBagLayout());
			jPanel.add(getJTabbedPane(), gridBagConstraints9);
		}
		return jPanel;
	}

	/**
	 * This method initializes jTabbedPane	
	 * 	
	 * @return javax.swing.JTabbedPane	
	 */
	private JTabbedPane getJTabbedPane() {
		if (jTabbedPane == null) {
			jTabbedPane = new JTabbedPane();
			jTabbedPane.addTab("Results", null, getResults(), null);
			jTabbedPane.addTab("Movie", null, getMovie(), null);
		}
		return jTabbedPane;
	}

	/**
	 * This method initializes Results	
	 * 	
	 * @return javax.swing.JPanel	
	 */
	private JPanel getResults() {
		if (Results == null) {
			Results = new JPanel();
			Results.setLayout(new GridBagLayout());
		}
		return Results;
	}

	/**
	 * This method initializes Movie	
	 * 	
	 * @return javax.swing.JPanel	
	 */
	private JPanel getMovie() {
		if (Movie == null) {
			Movie = new JPanel();
			Movie.setLayout(new GridBagLayout());
		}
		return Movie;
	}

}
