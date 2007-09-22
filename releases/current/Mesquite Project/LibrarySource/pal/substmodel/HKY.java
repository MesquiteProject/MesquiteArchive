// HKY.java//// (c) 1999-2001 PAL Development Core Team//// This package may be distributed under the// terms of the Lesser GNU General Public License (LGPL)package pal.substmodel;import pal.misc.*;import java.io.*;/** * Hasegawa-Kishino-Yano model of nucleotide evolution * * @version $Id: HKY.java,v 1.5 2001/07/13 14:39:13 korbinian Exp $ * * @author Korbinian Strimmer */public class HKY extends NucleotideModel implements Serializable{	/**	 * Constructor 1	 *	 * @param kappa transition/transversion rate ratio	 * @param freq nucleotide frequencies	 */	public HKY(double kappa, double[] freq)	{		super(freq);				this.kappa = kappa;				makeHKY();		fromQToR();				showSE = false;	} 	/**	 * Constructor 2	 *	 * @param params parameter list	 * @param freq nucleotide frequencies	 */	public HKY(double[] params, double[] freq)	{		this(params[0], freq);	}	// Get numerical code describing the model type	public int getModelID()	{		return 2;	}   	// interface Report 	public void report(PrintWriter out)	{		out.println("Model of substitution: HKY (Hasegawa et al. 1985)");		out.print("Transition/transversion rate ratio kappa: ");		format.displayDecimal(out, kappa, 2);		if (showSE)		{			out.print("  (S.E. ");			format.displayDecimal(out, kappaSE, 2);			out.print(")");		}		out.println();		out.println();		printFrequencies(out);		printRatios(out);	}		// interface Parameterized	public int getNumParameters()	{		return 1;	}		public void setParameter(double param, int n)	{		kappa = param;				makeHKY();		fromQToR();	}	public double getParameter(int n)	{		return kappa;	}	public void setParameterSE(double paramSE, int n)	{		kappaSE = paramSE;			showSE = true;	}	public double getLowerLimit(int n)	{		return 0.0001;	}		public double getUpperLimit(int n)	{		return 100.0;	}		public double getDefaultValue(int n)	{		return 4.0;	}	//	// Private stuff	// 	private boolean showSE;	private double kappa, kappaSE;	// Make HKY model	private void makeHKY()	{		// Q matrix		rate[0][1] = 1; rate[0][2] = kappa; rate[0][3] = 1;		rate[1][2] = 1; rate[1][3] = kappa;		rate[2][3] = 1;	}}