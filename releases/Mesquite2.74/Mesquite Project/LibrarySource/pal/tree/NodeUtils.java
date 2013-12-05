// NodeUtils.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package pal.tree;

import pal.misc.*;
import java.io.*;
import pal.io.*;
import pal.util.*;

/**
 * Helper routines for dealing with nodes. 
 *
 * @version $Id: NodeUtils.java,v 1.14 2001/07/13 14:39:13 korbinian Exp $
 *
 * @author Alexei Drummond
 * @author Korbinian Strimmer
 */
public class NodeUtils {

	/**
	 * Converts lengths to heights, *without* assuming contemporaneous
	 * tips.
	 */
	public static void lengths2Heights(Node root) {
	
		lengths2Heights(root, getGreatestDistance(root));
	}

	/**
	 * Converts lengths to heights, but maintains tip heights.
	 */
	public static void lengths2HeightsKeepTips(Node node, boolean useMax) {
		
		if (!node.isLeaf()) {
			for (int i = 0; i < node.getChildCount(); i++) {
				lengths2HeightsKeepTips(node.getChild(i), useMax);
			}

			double totalHL = 0.0;
			double maxHL = 0.0;
			double hl = 0.0;
			double maxH = 0.0;
			double h = 0.0;
			for (int i = 0; i < node.getChildCount(); i++) {
				h = node.getChild(i).getNodeHeight();
				hl = node.getChild(i).getBranchLength() + h;
				if (hl > maxHL) maxHL = hl;
				if (h > maxH) maxH = h;
				totalHL += hl;
			}
			if (useMax) {
				hl = maxHL; // set parent height to maximum parent height implied by children
			} else {
				hl = totalHL /  node.getChildCount(); // get mean parent height
				if (hl < maxH) hl = maxHL; // if mean parent height is not greater than all children height, fall back on max parent height.
			}
			node.setNodeHeight(hl); // set new parent height
			
			// change lengths in children to reflect changes.
			for (int i = 0; i < node.getChildCount(); i++) {
				h = node.getChild(i).getNodeHeight();
				node.getChild(i).setBranchLength(hl - h);
			}
		}
	}
	

	/**
	 * sets this nodes height value to newHeight and all children's
	 * height values based on length of branches.
	 */
	private static void lengths2Heights(Node node, double newHeight) {
	
		if (!node.isRoot()) {
			newHeight -= node.getBranchLength();
			node.setNodeHeight(newHeight);
		} else {
			node.setNodeHeight(newHeight);
		}
		
		for (int i = 0; i < node.getChildCount(); i++) {
			lengths2Heights(node.getChild(i), newHeight);
		} 
	}
  
  	/**
	 * Exchange field info between two nodes. Specifically
	 * identifiers, branch lengths, node heights and branch length
	 * SEs.
	 */
	public static void exchangeInfo(Node node1, Node node2) {
		
		Identifier swaps;
		double swapd;
		
		swaps = node1.getIdentifier();
		node1.setIdentifier(node2.getIdentifier());
		node2.setIdentifier(swaps);
		
		swapd = node1.getBranchLength();
		node1.setBranchLength(node2.getBranchLength());
		node2.setBranchLength(swapd);

		swapd = node1.getNodeHeight();
		node1.setNodeHeight(node2.getNodeHeight());
		node2.setNodeHeight(swapd);

		swapd = node1.getBranchLengthSE();
		node1.setBranchLengthSE(node2.getBranchLengthSE());
		node2.setBranchLengthSE(swapd);
	}

	/**
	 * determines branch lengths of this and all descendent nodes
	 * from heights
	 */
	public static void heights2Lengths(Node node) {
		heights2Lengths(node, true); //respect minimum
	}

	/**
	 * determines branch lengths of this and all descendent nodes
	 * from heights
	 */
	public static void heights2Lengths(Node node, boolean respectMinimum) {
	
		for (int i = 0; i < node.getChildCount(); i++) {
			heights2Lengths(node.getChild(i));
		}
		
		if (node.isRoot()) {
			node.setBranchLength(0.0);
		}
		else {
			node.setBranchLength(node.getParent().getNodeHeight() - node.getNodeHeight());
			if (respectMinimum && (node.getBranchLength() < BranchLimits.MINARC))
			{
				node.setBranchLength(BranchLimits.MINARC);
			}
		}
	}

	/**
	 * determines branch lengths of this node and its immediate descendent nodes
	 * from heights.
	 */
	public static void localHeights2Lengths(Node node, boolean respectMinimum) {
	
		for (int i = 0; i < node.getChildCount(); i++) {
			Node child = node.getChild(i);
			
			child.setBranchLength(node.getNodeHeight() - child.getNodeHeight());
		}
		
		if (node.isRoot()) {
			node.setBranchLength(0.0);
		}
		else {
			node.setBranchLength(node.getParent().getNodeHeight() - node.getNodeHeight());
			if (respectMinimum && (node.getBranchLength() < BranchLimits.MINARC))
			{
				node.setBranchLength(BranchLimits.MINARC);
			}
		}
	}

   
	/**
	 * Get the distance to furthest leaf from this nodes parent.
	 */
	private static double getGreatestDistance(Node node) {
		
		double distance = 0.0;
		if (!node.isLeaf()) {
			if (!node.isRoot()) {
				distance = node.getBranchLength();
			}
			double max = getGreatestDistance(node.getChild(0));
			double posmax = 0.0;
			for (int i = 1; i < node.getChildCount(); i++) {
				posmax = getGreatestDistance(node.getChild(i));
				if (posmax > max) max = posmax;
			}
			distance += max;
	    
	    		return distance;
		} else {
	    		return node.getBranchLength();
		}
	}

	/**
	 * Finds the largest child (in terms of node height).
	 */
	public static double findLargestChild(Node node) {
		// find child with largest height
		double max = node.getChild(0).getNodeHeight();
		for (int j = 1; j < node.getChildCount(); j++)
		{
			if (node.getChild(j).getNodeHeight() > max)
			{
				max = node.getChild(j).getNodeHeight();
			}
		}
		return max;
	}

	/**
	 * remove child
	 *
	 * @param node child node to be removed
	 */
	public static void removeChild(Node parent, Node child)
	{
		int rm = -1;
		for (int i = 0; i < parent.getChildCount(); i++)
		{
			if (child == parent.getChild(i))
			{
				rm = i;
				break;
			}
		}
		
		parent.removeChild(rm);
	}
	
	/**
	 * remove internal branch (collapse node with its parent)
	 *
	 * @param node node associated with internal branch
	 */
	public static void removeBranch(Node node)
	{
		if (node.isRoot() || node.isLeaf())
		{
			throw new IllegalArgumentException("INTERNAL NODE REQUIRED (NOT ROOT)");
		}
		
		Node parent = node.getParent();
		
		// add childs of node to parent
		// (node still contains the link to childs
		// to allow later restoration)
		int numChilds = node.getChildCount();
		for (int i = 0; i < numChilds; i++)
		{
			parent.addChild(node.getChild(i));
		}
		
		// remove node from parent
		// (link to parent is restored and the 
		// position is stored)
		int rm = -1;
		for (int i = 0; i < parent.getChildCount(); i++)
		{
			if (node == parent.getChild(i))
			{
				rm = i;
				break;
			}
		}
		parent.removeChild(rm);
		node.setParent(parent);
		node.setNumber(rm);
	}

	/**
	 * restore internal branch 
	 *
	 * @param node node associated with internal branch
	 */
	public static void restoreBranch(Node node)
	{
		if (node.isRoot() || node.isLeaf())
		{
			throw new IllegalArgumentException("INTERNAL NODE REQUIRED (NOT ROOT)");
		}
		
		Node parent = node.getParent();
		
		// remove childs of node from parent and make node their parent
		int numChilds = node.getChildCount();
		for (int i = 0; i < numChilds; i++)
		{
			Node c = node.getChild(i);
			removeChild(parent, c);
			c.setParent(node);
		}
		
		// insert node into parent
		parent.insertChild(node, node.getNumber());
	}

	
	
	/**
	 * join two childs, introducing a new node/branch in the tree
	 * that replaces the first child
	 * 
	 * @param n1 number of first child
	 * @param n2 number of second child 
	 */
	public static void joinChilds(Node node, int n1, int n2) {
	
		if (n1 == n2) {
			throw new IllegalArgumentException("CHILDREN MUST BE DIFFERENT");
		}
				
		int c1, c2;
		if (n2 < n1)
		{
			c1 = n2;
			c2 = n1;
		}
		else
		{
			c1 = n1;
			c2 = n2;
		}

		Node newNode = NodeFactory.createNode();

		Node child1 = node.getChild(c1);
		Node child2 = node.getChild(c2);
		
		node.setChild(c1, newNode);
		newNode.setParent(node);
		node.removeChild(c2); // now parent of child2 = null
				
		newNode.addChild(child1);
		newNode.addChild(child2);
	}

	/**
	 * determine preorder successor of this node
	 *
	 * @return next node
	 */ 
	public static Node preorderSuccessor(Node node) {
		
		Node next = null;
		
		if (node.isLeaf()) {
			Node cn = node, ln = null; // Current and last node
			
			// Go up
			do
			{
				if (cn.isRoot())
				{
					next = cn;
					break;
				}
				ln = cn;
				cn = cn.getParent();
			}
			while (cn.getChild(cn.getChildCount()-1) == ln);
				
			// Determine next node
			if (next == null)
			{
				// Go down one node
				for (int i = 0; i < cn.getChildCount()-1; i++)
				{
					if (cn.getChild(i) == ln)
					{
						next = cn.getChild(i+1);
						break;
					}
				}
			}
		}
		else
		{
			next = node.getChild(0);
		}
		
		return next;
	}

	/**
	 * determine postorder successor of a node
	 *
	 * @return next node
	 */ 
	public static Node postorderSuccessor(Node node) {
		
		Node cn = null;
		Node parent = node.getParent();
		
		if (node.isRoot())
		{
			cn = node;
		}
		else
		{
			
			// Go up one node
			if (parent.getChild(parent.getChildCount()-1) == node) {
				return parent;
			}
			// Go down one node
			for (int i = 0; i < parent.getChildCount()-1; i++)
			{
				if (parent.getChild(i) == node)
				{
					cn = parent.getChild(i+1);
					break;
				}
			}
		}
		
		// Go down until leaf
		while (cn.getChildCount() > 0)
		{
				
			cn = cn.getChild(0);
		}
		
		return cn;
	}

	/**
	 * prints node in New Hamshire format.
	 */
	static void printNH(PrintWriter out, Node node, 
		boolean printLengths, boolean printInternalLabels) {
		
		printNH(out, node, printLengths, printInternalLabels, 0);
	}

	
	private static int printNH(PrintWriter out, Node node,
		boolean printLengths, boolean printInternalLabels, int column) {
		
		column = breakLine(out, column);

		if (!node.isLeaf())
		{
			out.print("(");
			column++;
			
			for (int i = 0; i < node.getChildCount(); i++)
			{
				if (i != 0)
				{
					out.print(",");
					column++;
				}
				
				column = printNH(out, node.getChild(i), printLengths, printInternalLabels, column);
			}
			
			out.print(")");
			column++;
		}
		
		if (!node.isRoot())
		{
			if (node.isLeaf() || printInternalLabels)
			{
				column = breakLine(out, column);
		
				String id = node.getIdentifier().toString();
				out.print(id);
				column += id.length();
			}

			if (printLengths)
			{
				out.print(":");
				column++;
		
				column = breakLine(out, column);
		
				column += FormattedOutput.getInstance().displayDecimal(out, node.getBranchLength(), 7);
			}
		}
		
		return column;
	}

	private static int breakLine(PrintWriter out, int column)
	{
		if (column > 70)
		{
			out.println();
			column = 0;
		}
		
		return column;
	}

	/**
	 * Returns the first node in this tree that has the
	 * required identifier.
	 */
	public static Node findByIdentifier(Node node, Identifier identifier) {

		Log.getDefaultLogger().debug("node identifier = " + node.getIdentifier());
		Log.getDefaultLogger().debug("target identifier = " + identifier);

		if (node.getIdentifier().getName().equals(identifier.getName())) {
			return node;
		} else {
			Node pos = null;
			for (int i = 0; i < node.getChildCount(); i++) {
				pos = findByIdentifier(node.getChild(i), identifier);
				if (pos != null) return pos;
			}
			//if (pos == null && !node.isRoot()) {
			//	pos = findByIdentifier(node.getParent(), identifier);
			//}
			if (pos != null) return pos;
			return null;
		}
	}

	/**
	 * Root tree at this node.
	 */
	public static Node root(Node node) {
	
		if (!node.isRoot()) {
	  
			Node myParent = node.getParent();
			removeChild(myParent, node);
		
			root(myParent);
		
			while (myParent.getChildCount() == 1) {
				myParent = myParent.getChild(0);
			}
		
			node.addChild(myParent);
			lengths2Heights(node);
		}
		return node;
	}

	/**
	 * Root the tree above the node with this identifier.
	 */
	public static Node rootAbove(Identifier id, Node root) {
		return rootAbove(findByIdentifier(root, id));
	}

	/**
	 * Root tree above this node;
	 */
	public static Node rootAbove(Node node) {
		
		if (!node.isRoot()) {

			Node root = NodeFactory.createNode();
		
			Node myParent = node.getParent();
			removeChild(myParent, node);
		
			Log.getDefaultLogger().debug("Before root() call");
			root(myParent);
			Log.getDefaultLogger().debug("After root() call");
			
			while (myParent.getChildCount() == 1) {
				myParent = myParent.getChild(0);
			}
		
			root.addChild(myParent);
			root.addChild(node);

			lengths2Heights(root);

			return root;
		
		} else return node;
	}

	/**
	 * determine distance to root
	 *
	 * @return distance to root
	 */
	public static double getDistanceToRoot(Node node)
	{
		if (node.isRoot())
		{
			return 0.0;
		}
		else
		{
			return node.getBranchLength() + getDistanceToRoot(node.getParent());
		}
	}

	/**
	 * Return the number of terminal leaves below this node or 1 if this is
	 * a terminal leaf.
	 */
	public static int getLeafCount(Node node) {
		
		int count = 0;
		if (!node.isLeaf()) {
			for (int i = 0; i < node.getChildCount(); i++) {
				count += getLeafCount(node.getChild(i));
			}
		} else {
			count = 1;
		}
		return count;
	}

		/** returns number of branches centered around an internal node in an unrooted tree */
	public static final int getUnrootedBranchCount(Node center) {
		if (center.isRoot()) 	{
			return center.getChildCount();
		}
		else {
			return center.getChildCount()+1;
		}
	}
}

