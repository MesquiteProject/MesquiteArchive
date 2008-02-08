/* Mesquite source code.  Copyright 1997-2007 W. Maddison and D. Maddison.
Version 2.01, December 2007.
Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. 
The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.
Perhaps with your help we can be more than a few, and make Mesquite better.

Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
Mesquite's web site is http://mesquiteproject.org

This source code and its compiled class files are free and modifiable under the terms of 
GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)
 */

package mesquite.consensus.lib;

import java.util.Vector;

import mesquite.consensus.lib.Bipartition;
import mesquite.lib.*;

public class BipartitionVector extends Vector {
	boolean rooted = true;
	int numTaxa =  0;
	int numTreesTotal = 0;
	Bits[] nodes;
	int princess = 0;
	Taxa taxa;
	Bits allTaxa;
	Bits bits1;
	Bits bits2;
	public static final int MAJRULEMODE = 0;
	public static final int STRICTMODE = 1;
	int mode = STRICTMODE;

	public void setTaxa(Taxa taxa){
		numTreesTotal = 0;
		this.taxa = taxa;
		numTaxa = taxa.getNumTaxa();
		nodes = new Bits[MesquiteTree.standardNumNodeSpaces(taxa)];
		for (int i= 0; i<nodes.length; i++)
			nodes[i] = new Bits(numTaxa);
		allTaxa = new Bits(numTaxa);
		bits1 = new Bits(numTaxa);
		bits2 = new Bits(numTaxa);

	}
	public void zeroFrequencies(){
		numTreesTotal = 0;
		allTaxa.clearAllBits();
		for (int i=0; i<size(); i++){
			Bipartition b = getBipart(i);
			b.freq = 0;
		}
	}
	public Bipartition getBipart(int i){
		if (i<size() && i>=0)
			return (Bipartition)elementAt(i);
		return null;
	}

	public void dump(){
		sort(true);
		MesquiteMessage.println("\nBipartition frequencies");
		for (int i=0; i<size(); i++){
			Bipartition stored = getBipart(i);
			double freq = stored.freq*1.0/numTreesTotal;
			if (freq>0.05) {
				String s = MesquiteDouble.toPrecisionSpecifiedString(freq,3);
				MesquiteMessage.println(stored.bits.toAsteriskString() + "   " +s);
			}
		}

	}
	public void tradeValues(int i, int j){
		long freq1 = ((Bipartition)elementAt(i)).getFreq();
		Bits bits1 = ((Bipartition)elementAt(i)).getBits();
		long freq2 = ((Bipartition)elementAt(j)).getFreq();
		Bits bits2 = ((Bipartition)elementAt(j)).getBits();
		((Bipartition)elementAt(i)).setFreq(freq2);
		((Bipartition)elementAt(j)).setFreq(freq1);
		((Bipartition)elementAt(i)).setBits(bits2);
		((Bipartition)elementAt(j)).setBits(bits1);
	}

	public void sort(boolean descending){
		for (int i=0; i<size(); i++) {
			for (int j= i+1; j<size(); j++) {
				long freq1 = ((Bipartition)elementAt(i)).getFreq();
				long freq2 = ((Bipartition)elementAt(j)).getFreq();
				if ((freq1<freq2 && descending) || (freq1>freq2 && !descending))
					tradeValues(i,j);
			}
		}
	}

	/** removes all bipartitions with frequencies lower that that specified */
	public void removeBipartitions(int lowerFrequencyLimit) {
		for (int i=size()-1; i>=0; i--) {
			if (((Bipartition)elementAt(i)).getFreq() < lowerFrequencyLimit)
				remove(i);
		}

	}

	/** Returns a new BipartitionVector that contains only those Bipartitions that are in a frequency greater than or equal to the lowerFrequencyLimit */
	public BipartitionVector getFilteredBipartitions (int lowerFrequencyLimit) {
		BipartitionVector bpv = new BipartitionVector();
		bpv.setTaxa(taxa);
		for (int i=0; i<size(); i++){
			Bipartition bp = ((Bipartition)elementAt(i));
			if (bp.getFreq() >= lowerFrequencyLimit) {
				Bipartition newbp = new Bipartition(numTaxa);
				newbp.bits.setBits(bp.getBits());
				newbp.setFreq(bp.getFreq());
				bpv.addElement(newbp);
			}
		}
		return bpv;
	}

	public boolean compatible(Bipartition bp, Bits bits){
		bits1.setBits(bp.getBits());
		bits2.setBits(bits);
		bits1.andBits(bits2);
		if (bits1.anyBitsOn()){   // intersection has something in it
			bits1.setBits(bp.getBits());
			bits1.invertAllBits();
			bits1.andBits(bits2);  // intersection of bits1 complement and bits2
			if (bits1.anyBitsOn()){   // intersection has something in it
				bits1.setBits(bp.getBits());
				bits2.invertAllBits();
				bits1.andBits(bits2);  // intersection of bits2 complement and bits1
				if (bits1.anyBitsOn()){   // intersection has something in it
					bits1.setBits(bp.getBits());
					bits1.invertAllBits();
					bits2.invertAllBits();
					bits1.andBits(bits2);  // intersection of bits2 complement and bits1 complement
					if (bits1.anyBitsOn())   // intersection has something in it
							return false;
				}

			}
		}
		return true;
	}


	private Bipartition addBipart(Bits bits){
		switch (mode) {
		case STRICTMODE: 
			if (numTreesTotal==0) {  // first tree; just add it
				Bipartition bipart = new Bipartition(numTaxa);
				bipart.bits.setBits(bits);
				addElement(bipart);
				return bipart;
			} else {
				boolean foundConflict = false;
				int identical = -1;
				for (int i=size()-1; i>=0; i--){
					Bipartition stored = getBipart(i);
					if (compatible(stored,bits)){ //then we are ok
						if (stored.equals(bits, rooted))  // record this in case we need to delete it later
							identical=i;
					} else {
						foundConflict=true;
						remove(i);
					}
				}
				if (foundConflict && identical>0)
					remove(identical);
			}
			return null;
		case MAJRULEMODE: 
			for (int i=0; i<size(); i++){
				Bipartition stored = getBipart(i);

				if (stored.equals(bits, rooted)){
					stored.increment();
					return stored;
				}
			}
			Bipartition bipart = new Bipartition(numTaxa);
			bipart.bits.setBits(bits);
			bipart.increment();
			addElement(bipart);
			return bipart;
		default:
			return null;
		}
	}




	private void getPartitions(Tree tree, int node){
		if (tree.nodeIsTerminal(node)){
			nodes[node].setBit(tree.taxonNumberOfNode(node));
			return;
		}
		nodes[node].clearAllBits();
		for (int daughter=tree.firstDaughterOfNode(node); tree.nodeExists(daughter); daughter = tree.nextSisterOfNode(daughter) ) {
			getPartitions(tree, daughter);
			nodes[node].orBits(nodes[daughter]);
		}
		if (node != princess || rooted)
			addBipart(nodes[node]);
	}

	/** adds tree to existing */
	public void addTree(Tree tree){
		princess = tree.firstDaughterOfNode(tree.getRoot());
		for (int i= 0; i<nodes.length; i++)
			nodes[i].clearAllBits();
		getPartitions(tree, tree.getRoot());
		numTreesTotal++;
		allTaxa.orBits(nodes[tree.getRoot()]);
	}

	NameReference freqRef = NameReference.getNameReference("consensusFrequency");
	private void resolveByBypartition(MesquiteTree tree, Bipartition stored){
		int newNode = tree.makeClade(stored.bits);
		if (mode==MAJRULEMODE) {
			double prop = stored.freq*1.0/numTreesTotal;
			tree.setAssociatedDouble(freqRef, newNode, prop);
			tree.setNodeLabel(MesquiteDouble.toStringDigitsSpecified(prop, 3), newNode);
		}
	}
	public Tree makeTree(double minFreq){
		MesquiteTree tree = new MesquiteTree(taxa);
		tree.setToDefaultBush(allTaxa.numBitsOn(), false);
		for (int i=0; i<size(); i++){
			Bipartition stored = getBipart(i);
			if (!rooted)
				stored.bits.standardizeComplement(0);
			if (mode==STRICTMODE) {
				resolveByBypartition(tree, stored);
			} else {
				if (minFreq==0.5) {
					if (stored.freq*1.0/numTreesTotal>minFreq)
						resolveByBypartition(tree, stored);
				}
				else	if (stored.freq*1.0/numTreesTotal>=minFreq){
					resolveByBypartition(tree, stored);
				}
			}
		}

		return tree;
	}
	public Tree makeTree(){
		return makeTree(0.5);
	}

	public int getMode() {
		return mode;
	}
	public void setMode(int mode) {
		this.mode = mode;
	}
}
