

import javax.swing.JTextField;
import java.awt.Dimension;
import javax.swing.JFrame;
import javax.swing.JPanel;
import java.awt.BorderLayout;
import java.awt.GridBagLayout;
import java.awt.GridBagConstraints;
import javax.swing.JLabel;
import javax.swing.JCheckBox;
import javax.swing.JButton;
public class mainframe {

	private JFrame jFrame = null;  //  @jve:decl-index=0:visual-constraint="62,0"
	private JPanel mainframe = null;
	private JPanel Center = null;
	private JTextField Sequence11 = null;
	private JLabel Seq1 = null;
	private JLabel testing = null;
	private JPanel jPanel = null;
	private JLabel Structure = null;
	private JTextField Structure2 = null;
	private JCheckBox structureCheckbox = null;
	private JLabel Structure1 = null;
	private JCheckBox AlignmentCheckbox1 = null;
	private JCheckBox contraFoldCheckbox11 = null;
	private JLabel contrafold = null;
	private JButton getInteraction = null;
	private JPanel Right = null;
	private JTextField Sequence111 = null;
	private JLabel Seq11 = null;
	private JLabel testing1 = null;
	private JLabel Structure3 = null;
	private JTextField Structure21 = null;
	private JCheckBox structureCheckbox1 = null;
	private JLabel Structure11 = null;
	private JCheckBox AlignmentCheckbox11 = null;
	private JCheckBox contraFoldCheckbox111 = null;
	private JLabel contrafold1 = null;
	/**
	 * This method initializes jFrame	
	 * 	
	 * @return javax.swing.JFrame	
	 */
	private JFrame getJFrame() {
		if (jFrame == null) {
			jFrame = new JFrame();
			jFrame.setSize(new Dimension(795, 307));
			jFrame.setTitle("Antisense RNA");
			jFrame.setContentPane(getMainframe());
		}
		return jFrame;
	}

	/**
	 * This method initializes mainframe	
	 * 	
	 * @return javax.swing.JPanel	
	 */
	private JPanel getMainframe() {
		if (mainframe == null) {
			mainframe = new JPanel();
			mainframe.setLayout(new BorderLayout());
			mainframe.add(getsRNA(), BorderLayout.WEST);
			mainframe.add(getJPanel(), BorderLayout.NORTH);
			mainframe.add(getGetInteraction(), BorderLayout.SOUTH);
			mainframe.add(getRight(), BorderLayout.EAST);
		}
		return mainframe;
	}

	/**
	 * This method initializes Center	
	 * 	
	 * @return javax.swing.JPanel	
	 */
	private JPanel getsRNA() {
		if (Center == null) {
			GridBagConstraints gridBagConstraints8 = new GridBagConstraints();
			gridBagConstraints8.gridx = 2;
			gridBagConstraints8.gridy = 10;
			contrafold = new JLabel();
			contrafold.setText("contraFOLD");
			GridBagConstraints gridBagConstraints3 = new GridBagConstraints();
			gridBagConstraints3.gridx = 0;
			gridBagConstraints3.gridy = 7;
			GridBagConstraints gridBagConstraints = new GridBagConstraints();
			gridBagConstraints.gridx = 0;
			gridBagConstraints.gridy = 6;
			GridBagConstraints gridBagConstraints7 = new GridBagConstraints();
			gridBagConstraints7.gridx = 2;
			gridBagConstraints7.gridy = 6;
			Structure1 = new JLabel();
			Structure1.setText("Alignment");
			GridBagConstraints gridBagConstraints6 = new GridBagConstraints();
			gridBagConstraints6.gridx = 0;
			gridBagConstraints6.gridy = 4;
			GridBagConstraints gridBagConstraints5 = new GridBagConstraints();
			gridBagConstraints5.fill = GridBagConstraints.VERTICAL;
			gridBagConstraints5.gridy = 5;
			gridBagConstraints5.weightx = 3.0;
			gridBagConstraints5.gridx = 2;
			GridBagConstraints gridBagConstraints4 = new GridBagConstraints();
			gridBagConstraints4.gridx = 2;
			gridBagConstraints4.gridy = 4;
			Structure = new JLabel();
			Structure.setText("Structure");
			GridBagConstraints gridBagConstraints2 = new GridBagConstraints();
			gridBagConstraints2.gridx = 2;
			gridBagConstraints2.gridwidth = 3;
			gridBagConstraints2.gridy = 0;
			testing = new JLabel();
			testing.setText("sRNA");
			GridBagConstraints gridBagConstraints11 = new GridBagConstraints();
			gridBagConstraints11.gridx = 2;
			gridBagConstraints11.gridheight = 2;
			gridBagConstraints11.gridy = 1;
			Seq1 = new JLabel();
			Seq1.setText("Sequence");
			GridBagConstraints gridBagConstraints1 = new GridBagConstraints();
			gridBagConstraints1.fill = GridBagConstraints.VERTICAL;
			gridBagConstraints1.gridx = 2;
			gridBagConstraints1.gridy = 3;
			gridBagConstraints1.gridwidth = 2;
			gridBagConstraints1.weightx = 1.0;
			Center = new JPanel();
			Center.setLayout(new GridBagLayout());
			Center.add(getSequence11(), gridBagConstraints1);
			Center.add(Seq1, gridBagConstraints11);
			Center.add(testing, gridBagConstraints2);
			Center.add(Structure, gridBagConstraints4);
			Center.add(getStructure2(), gridBagConstraints5);
			Center.add(getStructureCheckbox(), gridBagConstraints6);
			Center.add(Structure1, gridBagConstraints7);
			Center.add(getAlignmentCheckbox1(), gridBagConstraints);
			Center.add(getContraFoldCheckbox11(), gridBagConstraints3);
			Center.add(contrafold, gridBagConstraints8);
		}
		return Center;
	}

	/**
	 * This method initializes Sequence11	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getSequence11() {
		if (Sequence11 == null) {
			Sequence11 = new JTextField(4);
			Sequence11.setText("testingdf");
			//Sequence11.setSize(300,10);
		}
		return Sequence11;
	}

	/**
	 * This method initializes jPanel	
	 * 	
	 * @return javax.swing.JPanel	
	 */
	private JPanel getJPanel() {
		if (jPanel == null) {
			jPanel = new JPanel();
			jPanel.setLayout(new GridBagLayout());
		}
		return jPanel;
	}

	/**
	 * This method initializes Structure2	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getStructure2() {
		if (Structure2 == null) {
			Structure2 = new JTextField();
		}
		return Structure2;
	}

	/**
	 * This method initializes structureCheckbox	
	 * 	
	 * @return javax.swing.JCheckBox	
	 */
	private JCheckBox getStructureCheckbox() {
		if (structureCheckbox == null) {
			structureCheckbox = new JCheckBox();
		}
		return structureCheckbox;
	}

	/**
	 * This method initializes AlignmentCheckbox1	
	 * 	
	 * @return javax.swing.JCheckBox	
	 */
	private JCheckBox getAlignmentCheckbox1() {
		if (AlignmentCheckbox1 == null) {
			AlignmentCheckbox1 = new JCheckBox();
		}
		return AlignmentCheckbox1;
	}

	/**
	 * This method initializes contraFoldCheckbox11	
	 * 	
	 * @return javax.swing.JCheckBox	
	 */
	private JCheckBox getContraFoldCheckbox11() {
		if (contraFoldCheckbox11 == null) {
			contraFoldCheckbox11 = new JCheckBox();
		}
		return contraFoldCheckbox11;
	}

	/**
	 * This method initializes getInteraction	
	 * 	
	 * @return javax.swing.JButton	
	 */
	private JButton getGetInteraction() {
		if (getInteraction == null) {
			getInteraction = new JButton();
			getInteraction.setText("getInteraction");
		}
		return getInteraction;
	}

	/**
	 * This method initializes Right	
	 * 	
	 * @return javax.swing.JPanel	
	 */
	private JPanel getRight() {
		if (Right == null) {
			GridBagConstraints gridBagConstraints81 = new GridBagConstraints();
			gridBagConstraints81.gridx = 2;
			gridBagConstraints81.gridy = 4;
			contrafold1 = new JLabel();
			contrafold1.setText("contraFOLD");
			GridBagConstraints gridBagConstraints31 = new GridBagConstraints();
			gridBagConstraints31.gridx = 0;
			gridBagConstraints31.gridy = 4;
			GridBagConstraints gridBagConstraints9 = new GridBagConstraints();
			gridBagConstraints9.gridx = 0;
			gridBagConstraints9.gridy = 3;
			GridBagConstraints gridBagConstraints71 = new GridBagConstraints();
			gridBagConstraints71.gridx = 2;
			gridBagConstraints71.gridy = 3;
			Structure11 = new JLabel();
			Structure11.setText("Alignment");
			GridBagConstraints gridBagConstraints61 = new GridBagConstraints();
			gridBagConstraints61.gridx = 0;
			gridBagConstraints61.gridy = 2;
			GridBagConstraints gridBagConstraints51 = new GridBagConstraints();
			gridBagConstraints51.fill = GridBagConstraints.VERTICAL;
			gridBagConstraints51.gridy = 2;
			gridBagConstraints51.weightx = 1.0;
			gridBagConstraints51.gridx = 3;
			GridBagConstraints gridBagConstraints41 = new GridBagConstraints();
			gridBagConstraints41.gridx = 2;
			gridBagConstraints41.gridy = 2;
			Structure3 = new JLabel();
			Structure3.setText("Structure");
			GridBagConstraints gridBagConstraints21 = new GridBagConstraints();
			gridBagConstraints21.gridx = 3;
			gridBagConstraints21.gridy = 0;
			testing1 = new JLabel();
			testing1.setText("JLabel");
			GridBagConstraints gridBagConstraints111 = new GridBagConstraints();
			gridBagConstraints111.gridheight = 2;
			gridBagConstraints111.gridy = 0;
			gridBagConstraints111.gridx = 2;
			Seq11 = new JLabel();
			Seq11.setText("Sequence");
			GridBagConstraints gridBagConstraints12 = new GridBagConstraints();
			gridBagConstraints12.fill = GridBagConstraints.VERTICAL;
			gridBagConstraints12.gridy = 1;
			gridBagConstraints12.weightx = 1.0;
			gridBagConstraints12.gridx = 3;
			Right = new JPanel();
			Right.setLayout(new GridBagLayout());
			Right.add(getSequence111(), gridBagConstraints12);
			Right.add(Seq11, gridBagConstraints111);
			Right.add(testing1, gridBagConstraints21);
			Right.add(Structure3, gridBagConstraints41);
			Right.add(getStructure21(), gridBagConstraints51);
			Right.add(getStructureCheckbox1(), gridBagConstraints61);
			Right.add(Structure11, gridBagConstraints71);
			Right.add(getAlignmentCheckbox11(), gridBagConstraints9);
			Right.add(getContraFoldCheckbox111(), gridBagConstraints31);
			Right.add(contrafold1, gridBagConstraints81);
		}
		return Right;
	}

	/**
	 * This method initializes Sequence111	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getSequence111() {
		if (Sequence111 == null) {
			Sequence111 = new JTextField();
		}
		return Sequence111;
	}

	/**
	 * This method initializes Structure21	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getStructure21() {
		if (Structure21 == null) {
			Structure21 = new JTextField();
		}
		return Structure21;
	}

	/**
	 * This method initializes structureCheckbox1	
	 * 	
	 * @return javax.swing.JCheckBox	
	 */
	private JCheckBox getStructureCheckbox1() {
		if (structureCheckbox1 == null) {
			structureCheckbox1 = new JCheckBox();
		}
		return structureCheckbox1;
	}

	/**
	 * This method initializes AlignmentCheckbox11	
	 * 	
	 * @return javax.swing.JCheckBox	
	 */
	private JCheckBox getAlignmentCheckbox11() {
		if (AlignmentCheckbox11 == null) {
			AlignmentCheckbox11 = new JCheckBox();
		}
		return AlignmentCheckbox11;
	}

	/**
	 * This method initializes contraFoldCheckbox111	
	 * 	
	 * @return javax.swing.JCheckBox	
	 */
	private JCheckBox getContraFoldCheckbox111() {
		if (contraFoldCheckbox111 == null) {
			contraFoldCheckbox111 = new JCheckBox();
		}
		return contraFoldCheckbox111;
	}

}
