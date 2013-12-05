// Alignment.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package pal.alignment;

import pal.io.*;
import pal.datatype.*;
import pal.misc.*;

import java.io.*;


/**
 * interface for any alignment data.
 *
 * @version $Id: Alignment.java,v 1.11 2001/07/13 14:39:12 korbinian Exp $
 *
 * @author Alexei Drummond
 * @author Korbinian Strimmer
 */
public interface Alignment extends Serializable, IdGroup
{
	//
	// Public stuff
	//

	/** character used to designate gaps */
	char GAP = '-';

	// Abstract method

	/** sequence alignment at (sequence, site) */
	char getData(int seq, int site);	
	
   	/**
	 * Return number of sites for each sequence in this alignment
	 * @note for people who like accessor methods over public instance variables...
	 */
	int getSiteCount();	

	/**
	 * Return number of sequences in this alignment
	 */
	int getSequenceCount();

	/**
	 * Return DataType of this alignment.
	 */
	DataType getDataType();

	/**
	 * Sets the dataType of this alignment.
	 */
	void setDataType(DataType dataType);

	/**
	 * Returns string representation of single sequence in 
	 * alignment with gap characters included.
	 */
	String getAlignedSequenceString(int sequence);

	/**
	 * Returns frequency of character states.
	 */
	double[] getFrequency();

	/**
	 * Sets frequency of character states.
	 */
	void setFrequency(double[] frequencies);
}	

