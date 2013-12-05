// SimpleNode.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package pal.tree;

import pal.misc.*;
import java.io.*;
import pal.io.*;


/**
 * data structure for a node (includes branch) in a binary/non-binary
 * rooted/unrooted tree
 *
 * @version $Id: SimpleNode.java,v 1.9 2001/07/13 14:39:13 korbinian Exp $
 *
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 */
public class SimpleNode implements Node
{	
	/** parent node */
	private Node parent;
	
	/** number of child nodes */
	private int numChilds;

	/** number of node as displayed */
	private int number;
	
	/** sequences associated with node */
	private byte[] sequence;
	
	/** partial likelihood associated with node */
	private double[][][] partial;
	
	/** length of branch to parent node */
	private double length;
	
	/** standard error of length of branch to parent node */ 
	private double lengthSE;
	
	/** height of this node */
	private double height;

	/** height SE of this node */
	private double heightSE;

	
	/** identifier of node/associated branch */
	private Identifier identifier;



	// the following constructors should eventually become
	// "friendly" to prevent anyone calling them directly.
	// Instead, the NodeFactory should be used!
	
	/** constructor default node */
	public SimpleNode()
	{
		parent = null;
		child = null;
		numChilds = 0;
		length = 0.0;
		lengthSE = 0.0;
		height = 0.0;
		identifier = Identifier.ANONYMOUS;
		
		number = 0;
		sequence = null;
		partial = null;
	}	

	/** constructor used to clone a node and all children */
	public SimpleNode(Node n)
	{
		this(n, true);
	}	



	public void reset()
	{
		parent = null;
		child = null;
		numChilds = 0;
		length = 0.0;
		lengthSE = 0.0;
		height = 0.0;
		identifier = Identifier.ANONYMOUS;
		
		number = 0;
		sequence = null;
		partial = null;
	}

	public SimpleNode(Node n, boolean keepIds) {
		init(n, keepIds);
		for (int i = 0; i < n.getChildCount(); i++) {
			addChild(new SimpleNode(n.getChild(i)));
		}
	}


	protected void init(Node n) {
		init(n, true);
	}

	/**
	 * Initialized node instance variables based on given Node.
	 * children are ignored.
	 */
	protected void init(Node n, boolean keepId) {
		parent = null;
		length = n.getBranchLength();
		lengthSE = n.getBranchLengthSE();
		height = n.getNodeHeight();
		if (keepId) {
			identifier = n.getIdentifier();
		} else { identifier = Identifier.ANONYMOUS; }
		
		number = n.getNumber();
		sequence = n.getSequence();
		
		child = null;
		numChilds = 0;
	}

	/**
	 * Returns the parent node of this node.
	 */
	public final Node getParent() {
		return parent;
	}

	/** Set the parent node of this node. */
	public void setParent(Node node)
	{
		parent = node;
	}

	/**
	 * removes parent.
	 */
	public final void removeParent() {
		parent = null;
	}

	/**
	 * Returns the sequence at this node, in the form of a String.
	 */
	public String getSequenceString() {
		return new String(sequence);
	}

	/**
	 * Returns the sequence at this node, in the form of an array of bytes.
	 */
	public byte[] getSequence() {
		return sequence;
	}

	/**
	 * Sets the sequence at this node, in the form of an array of bytes.
	 */
	public void setSequence(byte[] s) {
		sequence = s;
	}

	/**
	 * Get the length of the branch attaching this node to its parent.
	 */
	public final double getBranchLength() {
		return length;
	}

	/**
	 * Set the length of the branch attaching this node to its parent.
	 */
	public final void setBranchLength(double value) {
		length = value;
	}

	/**
	 * Get the length SE of the branch attaching this node to its parent.
	 */
	public final double getBranchLengthSE() {
		return lengthSE;
	}

	/**
	 * Set the length SE of the branch attaching this node to its parent.
	 */
	public final void setBranchLengthSE(double value) {
		lengthSE = value;
	}


	/**
	 * Get the height of this node relative to the most recent node.
	 */
	public final double getNodeHeight() {
		return height;
	}

	/**
	 * Set the height of this node relative to the most recent node.
	 */
	public final void setNodeHeight(double value) {
		height = Math.abs(value);
	}

	/**
	 * Set the height SE of this node relative to the most recent node.
	 */
	public final void setNodeHeightSE(double value) {
		heightSE = value;
	}

	/**
	 * Set the height SE of this node relative to the most recent node.
	 */
	public final double getNodeHeightSE() {
		return heightSE;
	}

	/**
	 * Returns the identifier for this node.
	 */
	public final Identifier getIdentifier() {
		return identifier;
	}

	/**
	 * Set identifier for this node.
	 */
	public final Identifier setIdentifier(Identifier id) {
		identifier = id;
		return identifier;
	}

	public void setNumber(int n) {
		number = n;
	}

	public int getNumber() {
		return number;
	}

	

	/**
	 * get child node
	 *
	 * @param n number of child
	 *
	 * @return child node
	 */ 
	public Node getChild(int n)
	{
		return child[n];
	}
	
	/**
	 * set child node
	 *
	 * @param n number
	 * @node node new child node
	 */
	public void setChild(int n, Node node)
	{
		child[n] = node;
		child[n].setParent(this);
	}
	
	/**
	 * check whether this node is an internal node
	 *
	 * @return result (true or false)
	 */
	public boolean hasChildren()
	{
		return !isLeaf();
	}

	/**
	 * check whether this node is an external node
	 *
	 * @return result (true or false)
	 */
	public boolean isLeaf()
	{
		if (numChilds == 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	/**
	 * check whether this node is a root node 
	 *
	 * @return result (true or false)
	 */
	public boolean isRoot()
	{
		if (parent == null)
		{
			return true;
		}
		else
		{
			return false;
		}
	}


	/**
	 * add new child node
	 *
	 * @param n new child node
	 */
	public void addChild(Node n)
	{
		insertChild(n, numChilds);
	}

	/**
	 * add new child node (insertion at a specific position)
	 *
	 * @param n new child node
	 + @param pos position
	 */
	public void insertChild(Node n, int pos)
	{
		SimpleNode c = (SimpleNode) n;
	
		Node[] newChild = new Node[numChilds + 1];
		
		for (int i = 0; i < pos; i++)
		{
			newChild[i] = child[i];
		}
		newChild[pos] = c;
		for (int i = pos; i < numChilds; i++)
		{
			newChild[i+1] = child[i];
		}
		
		child = newChild;

		numChilds++;
		
		c.setParent(this);	
	}

	
	/**
	 * remove child
	 *
	 * @param n number of child to be removed
	 */
	public void removeChild(int n)
	{
		if (n >= numChilds)
		{
			throw new IllegalArgumentException("Nonexistent child");
		}
		Node[] newChild = new Node[numChilds-1];
		
		for (int i = 0; i < n; i++)
		{
			newChild[i] = child[i];
		}
		
		for (int i = n; i < numChilds-1; i++)
		{
			newChild[i] = child[i+1];
		}
	
		//remove parent link from removed child!
		((Node)child[n]).setParent(null);
	
		child = newChild;
		numChilds--;
	}

	/**
	 * determines the height of this node and its descendants
	 * from branch lengths, assuming contemporaneous tips.
	 */
	public void lengths2HeightsContemp()
	{
		double largestHeight = 0.0;
		
		if (!isLeaf())
		{
			for (int i = 0; i < numChilds; i++)
			{
				NodeUtils.lengths2Heights(getChild(i));
				
				double newHeight = 
					getChild(i).getNodeHeight() + getChild(i).getBranchLength();
				
				if (newHeight > largestHeight)
				{
					largestHeight = newHeight;
				}
			}
		}
		
		setNodeHeight(largestHeight);
	}
	
	/**
	 * Returns the number of children this node has.
	 */
	public int getChildCount() {
		return numChilds;
	}
	  



	public String toString() {
	
		StringWriter sw = new StringWriter();
		NodeUtils.printNH(new PrintWriter(sw), this, true, false);
		return sw.toString();
	}
	
	
	//
	// Private stuff
	//
	
	private Node[] child;
}

