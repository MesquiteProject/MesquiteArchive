// F81.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package pal.substmodel;

import pal.misc.*;

import java.io.*;


/**
 * Felsenstein 1981 model of nucleotide evolution
 *
 * @version $Id: F81.java,v 1.5 2001/07/13 14:39:13 korbinian Exp $
 *
 * @author Korbinian Strimmer
 */
public class F81 extends NucleotideModel implements Serializable
{
	/**
	 * constructor
	 *
	 * @param freq nucleotide frequencies
	 */
	public F81(double[] freq)
	{
		super(freq);
		
		makeF81();
		fromQToR();
	}
 
	// Get numerical code describing the model type
	public int getModelID()
	{
		return 4;
	}
 
 	// interface Report
 
	public void report(PrintWriter out)
	{
		out.println("Model of substitution: F81 (Felsenstein 1981)");
		printFrequencies(out);
		printRatios(out);
	}	
	
	// interface Parameterized
	
	public int getNumParameters()
	{
		return 0;
	}
	
	public void setParameter(double param, int n)
	{
		return;
	}

	public double getParameter(int n)
	{
		return 0.0;
	}

	public void setParameterSE(double paramSE, int n)
	{
		return;
	}

	public double getLowerLimit(int n)
	{
		return 0.0;
	}
	
	public double getUpperLimit(int n)
	{
		return 0.0;
	}
	
	public double getDefaultValue(int n)
	{
		return 0.0;
	}


	//
	// Private stuff
	// 

	// Make F81 model
	private void makeF81()
	{
		// Q matrix
		rate[0][1] = 1; rate[0][2] = 1.0; rate[0][3] = 1;
		rate[1][2] = 1; rate[1][3] = 1.0;
		rate[2][3] = 1;
	}
}

