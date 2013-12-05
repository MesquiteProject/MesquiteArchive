// WindowedMutationRate.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
 
package pal.mep;

import pal.math.*;
import pal.misc.*;
import pal.io.*;

import java.io.*;

/**
 * This class models a windowed mutation rate
 * (parameter: mu = mutation rate). <BR>
 *
 * @version $Id: WindowedMutationRate.java,v 1.2 2001/07/13 14:39:13 korbinian Exp $
 *
 * @author Alexei Drummond
 */
public class WindowedMutationRate extends MutationRateModel implements Report, Summarizable, Parameterized, Serializable
{	
	//
	// Public stuff
	//

	/** mutation rates */
	public double muBackground;
	public double muWindow;

	/** mutation rate SEs */
	public double muBackgroundSE;
	public double muWindowSE;
	
	/** mutation rate change times */
	public double windowCenter;
	public double windowWidth;

	public boolean backgroundFixed = false;

	/**
	 * Construct demographic model with default settings
	 */
	public WindowedMutationRate(double windowCenter, double windowWidth, int units) {
	
		super();
	
		setUnits(units);

		this.windowCenter = windowCenter;
		this.windowWidth = windowWidth;

		muBackground = getDefaultValue(0);
		muWindow = getDefaultValue(0);	
	}


	/**
	 * Construct mutation rate model of a give rate in given units.
	 */
	public WindowedMutationRate(double muBackground, 
		double windowCenter, double windowWidth, int units) {
	
		super();

		this.muBackground = muBackground;
		backgroundFixed = true;

		this.windowCenter = windowCenter;
		this.windowWidth = windowWidth;

		muWindow = getDefaultValue(0);	
		setUnits(units);
	}

	/**
	 * Construct mutation rate model of a give rate in given units.
	 */
	public WindowedMutationRate(double muWindow, double muBackground, 
		double windowCenter, double windowWidth, int units, boolean fixedb) {
	
		super();

		this.muWindow = muWindow;
		this.muBackground = muBackground;
		backgroundFixed = fixedb;

		this.windowCenter = windowCenter;
		this.windowWidth = windowWidth;
	
		setUnits(units);
	}

	/**
	 * Construct mutation rate model of a give rate in given units.
	 */
	public WindowedMutationRate(double muWindow, double muBackground, 
		double windowCenter, double windowWidth, int units) {
	
		this(muWindow, muBackground, windowCenter, windowWidth, units, false);	
	}


	public Object clone()
	{
		return new WindowedMutationRate(muWindow, muBackground, 
			windowCenter, windowWidth, getUnits(), backgroundFixed); 
	}

	public String[] getSummaryTypes() {
		if (summaryTypes == null) {
			summaryTypes = new String[4];
			summaryTypes[0] = "window mu";
			summaryTypes[1] = "background mu";
			summaryTypes[2] = "window center";
			summaryTypes[3] = "window width";
		}
		return summaryTypes;
	}

	public double getSummaryValue(int summaryType) {
		
		switch (summaryType) {
			case 0: return muWindow;
			case 1: return muBackground;
			case 2: return windowCenter;
			case 3: return windowWidth;
		}
		throw new RuntimeException("Assertion error: unknown summary type :"+summaryType);
	}

	/**
	 * returns current day mutation rate.
	 */
	public double getMu()
	{
		return getMutationRate(0.0);
	}

		
	// Implementation of abstract methods
	
	public final double getMutationRate(double t)
	{
		if ((t > windowCenter - (windowWidth / 2.0)) && 
			(t <= windowCenter + (windowWidth / 2.0))) {
		
			return muWindow;
		} 

		return muBackground;
	}

	/**
	 * Window must not span zero!
	 */
	public final double getExpectedSubstitutions(double time)
	{
		double height = 0.0;
	
		// bit before window
		double totalTime = windowCenter - (windowWidth / 2.0);
		if (totalTime > time) return muBackground * time;
		if (totalTime >= 0.0) 
			height += muBackground * totalTime;	
		else System.err.println("Mutation window spans time zero!");
	
		// window
		if ((totalTime + windowWidth) > time) { 
			return height + (muWindow * (time - totalTime));
		}
		height += muWindow * windowWidth;
		
		totalTime += windowWidth;
		
		// bit after window
		return height + (muBackground * (time - totalTime));
	}

	/**
	 * Window must not span zero!
	 */
	public final double getTime(double expected)
	{
		//NOTE: VERY APPROXIMATE!!! SHOULD BE IMPLEMENTED PROPERLY
	
		return expected / muBackground;	
	}

	/**
	 * Linearly scales this mutation rate model.
	 * @param scale getExpectedSubstitutions should return scale instead of 1.0 at time t.
	 */
	public final void scale(double scale) {
		muBackground *= scale;
		muWindow *= scale;
	}
	
	// Parameterized interface

	public int getNumParameters()
	{
		if (backgroundFixed) return 1;
		return 2;
	}
	
	public double getParameter(int k)
	{
		switch (k) {
			case 0: return muWindow;
			case 1: return muBackground;
		}
		return muWindow;
	}

	public double getUpperLimit(int k)
	{
		return 1e12;
	}

	public double getLowerLimit(int k)
	{
		return 1e-12;
	}

	public double getDefaultValue(int k)
	{
		//arbitrary default values
		if (getUnits() == GENERATIONS) {
			return 1e-6;
		} else {
			return 1e-6;
		}
	}

	public void setParameter(double value, int k)
	{
		switch (k) {
			case 0: muWindow = value; break;
			case 1: muBackground = value; break;
		}
	}

	public void setParameterSE(double value, int k) {
		switch (k) {
			case 0: muWindowSE = value; break;
			case 1: muBackgroundSE = value; break;
		}
	}

	public String toString()
	{
		OutputTarget out = OutputTarget.openString();
		report(out);
		out.close();
		
		return out.getString();
	}
	
	public void report(PrintWriter out)
	{
		out.println("Mutation rate model: windowed mutation rate ");
			
		out.print("Unit of time: ");
		if (getUnits() == GENERATIONS)
		{
			out.print("generations");
		}
		else
		{
			out.print("expected substitutions");
		}
		out.println();
		out.println();
		out.println("Parameters of demographic function:");
		out.print("window = ");
		fo.displayDecimal(out, windowCenter - (windowWidth / 2.0), 6);
		out.print(" - ");
		fo.displayDecimal(out, windowCenter + (windowWidth / 2.0), 6);
		out.println();
		out.print("window mutation rate = ");
		fo.displayDecimal(out, muWindow, 9);
		out.println();
		out.print("background mutation rate = ");
		fo.displayDecimal(out, muBackground, 9);
		out.println();
		if (backgroundFixed) {
			out.println("background mutation rate fixed.");
		} else {
			out.println("background mutation rate free to vary.");
		}
	}

	public String toSingleLine() {
		String line = "";
		line += "win mu\t" + muWindow + "\t";
		line += "bg mu\t" + muBackground + "\t";
		line += "win cen\t" + windowCenter + "\t";
		line += "win wid\t" + windowWidth + "\t";
		return line;
	}



	String[] summaryTypes = null;
}

