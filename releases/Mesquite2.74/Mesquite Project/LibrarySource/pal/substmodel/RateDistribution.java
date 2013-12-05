// RateDistribution.java
//
// (c) 1999-2000 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package pal.substmodel;

import pal.misc.*;
import pal.io.*;

import java.io.*;


/**
 * abstract base class for models of rate variation over sites
 * employing a discrete rate distribution
 *
 * @version $Id: RateDistribution.java,v 1.7 2001/07/13 14:39:13 korbinian Exp $
 *
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 */

public abstract class RateDistribution implements Parameterized, Report, Cloneable, Serializable
{
	//
	// Public stuff
	//
		
	/** number of rate categories*/
	public int numRates;
	
	/** rates of each rate category */
	public double[] rate;
	
	/** probability of each rate */
	public double[] probability;

	/**
	 * construct discrete distribution
	 *
	 *  @param n number of rate categories
	 */
	public RateDistribution(int n)
	{
		format = FormattedOutput.getInstance();
		
		numRates = n;
		rate = new double[n];
		probability = new double[n];
	}
	
	// interface Report (remains abstract)
	
	// interface Parameterized (remains abstract)

	//
	// Protected stuff
	//
	
	protected FormattedOutput format;
	
	protected void printRates(PrintWriter out)
	{
		out.println("Relative rates and their probabilities:\n");
		format.displayIntegerWhite(out, numRates);
		out.println("   Rate      Probability");
		for (int i = 0; i < numRates; i++)
		{
			format.displayInteger(out, i+1, numRates);
			out.print("   ");
			format.displayDecimal(out, rate[i], 5);
			out.print("   ");
			format.displayDecimal(out, probability[i], 5);
			out.println();
		}
	}

	public Object clone() {
		try {
			RateDistribution rd = (RateDistribution)super.clone();
			return rd;
		} catch (CloneNotSupportedException e) {
			// this shouldn't happen
			throw new InternalError();
		}
	}
}

