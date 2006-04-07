/* Mesquite source code.  Copyright 1997-2005 W. Maddison and D. Maddison. Version 1.06, August 2005.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.lib;import java.awt.*;import mesquite.lib.duties.*;/* ======================================================================== *//** This establishes and curates color tables for character states for use in shading etc.Eventually facilitities should be provided for changing the colors, and distinct tables for categorical and continuousstates, and perhaps even distinct tables could be allowed at each character*/public  class MesquiteColorTable  {	static final int maxNumStates = 64;	public static final int DEFAULT = -1;	public static final int GRAYSCALE = 0;	public static final int COLORS = 1;	public static final int COLORS_NO_BW = 2;	public static final int BIGCOLORS = 3;	public static final int GREENSCALE = 4;	public static final int BLUESCALE = 5;	public static final int REDSCALE = 6;	private static Color[][] defaultColorTable = null;	private static Color[] defaultBIGColorTable = null;	private static Color[][] defaultGrayTable = null;	private static Color[][] defaultGreenTable = null;	private static Color[][] defaultBlueTable = null;	private static Color[][] defaultRedTable = null;	private Color[][] colorTable;	private int mode = COLORS;	static {		defaultColorTable = new Color[maxNumStates][maxNumStates];		for (int i=0; i<maxNumStates; i++) {			for (int j=0; j<maxNumStates; j++) {				defaultColorTable[i][j] = Color.white;			}		}		for (int maxState=1; maxState<maxNumStates; maxState++) {			defaultColorTable[maxState][0] = Color.white;			defaultColorTable[maxState][maxState] = Color.black;			for (int theState=1; theState<maxState; theState++) {				defaultColorTable[maxState][theState] = new Color(Color.HSBtoRGB((float)((maxState-theState) * 0.8 /maxState),(float)1.0,(float)1.0));				//defaultColorTable[maxState][theState] = new Color(Color.HSBtoRGB((float)(theState * 0.8 /maxState),(float)1.0,(float)1.0));			}		}		defaultBIGColorTable = new Color[1000];		//defaultBIGColorTable[0] = Color.white;		//defaultBIGColorTable[999] = Color.black;		for (int theState=0; theState<1000; theState++) {			defaultBIGColorTable[theState] = new Color(Color.HSBtoRGB((float)((1000-theState) * 0.8 /1000),(float)1.0,(float)1.0));		}		defaultGrayTable = new Color[maxNumStates][maxNumStates];		for (int i=0; i<maxNumStates; i++) {			for (int j=0; j<maxNumStates; j++) {				defaultGrayTable[i][j] = Color.white;			}		}		for (int maxState=1; maxState<maxNumStates; maxState++) {			defaultGrayTable[maxState][0] = Color.white;			defaultGrayTable[maxState][maxState] = Color.black;			for (int theState=1; theState<maxState; theState++) {				float c = (float)((maxState-theState) * 1.0 /maxState);				defaultGrayTable[maxState][theState] = new Color(c,c,c);			}		}				defaultGreenTable = new Color[maxNumStates][maxNumStates];		for (int i=0; i<maxNumStates; i++) {			for (int j=0; j<maxNumStates; j++) {				defaultGreenTable[i][j] = Color.white;			}		}		for (int maxState=1; maxState<maxNumStates; maxState++) {			defaultGreenTable[maxState][0] = Color.white;			defaultGreenTable[maxState][maxState] = Color.black;			for (int theState=1; theState<maxState; theState++) {				float c = (float)((maxState-theState) * 1.0 /maxState);				defaultGreenTable[maxState][theState] = new Color(c, 1, c);			}		}		defaultBlueTable = new Color[maxNumStates][maxNumStates];		for (int i=0; i<maxNumStates; i++) {			for (int j=0; j<maxNumStates; j++) {				defaultBlueTable[i][j] = Color.white;			}		}		for (int maxState=1; maxState<maxNumStates; maxState++) {			defaultBlueTable[maxState][0] = Color.white;			defaultBlueTable[maxState][maxState] = Color.black;			for (int theState=1; theState<maxState; theState++) {				float c = (float)((maxState-theState) * 1.0 /maxState);				defaultBlueTable[maxState][theState] = new Color(c, c, 1);			}		}		defaultRedTable = new Color[maxNumStates][maxNumStates];		for (int i=0; i<maxNumStates; i++) {			for (int j=0; j<maxNumStates; j++) {				defaultRedTable[i][j] = Color.white;			}		}		for (int maxState=1; maxState<maxNumStates; maxState++) {			defaultRedTable[maxState][0] = Color.white;			defaultRedTable[maxState][maxState] = Color.black;			for (int theState=1; theState<maxState; theState++) {				float c = (float)((maxState-theState) * 1.0 /maxState);				defaultRedTable[maxState][theState] = new Color(1, c, c);			}		}	}	public MesquiteColorTable() {		colorTable = new Color[maxNumStates][maxNumStates];		for (int i=0; i<maxNumStates; i++) {			for (int j=0; j<maxNumStates; j++) {				colorTable[i][j] = defaultColorTable[i][j];			}		}	}	/**sets whether to use grayscale or color*/	public void setMode(int mode){		this.mode = mode;	}	/**gets whether to use grayscale or color*/	public int getMode(){		return mode;	}	/** gets color for state i with given maximum state possible, from default color table*/	public static Color getDefaultColor(int maxState, int i, int mode) {		if (maxState<0 || i<0)			return Color.white;		else if (mode == GRAYSCALE){			if ((maxState)<defaultGrayTable.length && (i)<defaultGrayTable[0].length)				return defaultGrayTable[maxState][i];			else				return defaultGrayTable[maxNumStates-1][maxNumStates-1];					}		else if (mode == GREENSCALE){			if ((maxState)<defaultGreenTable.length && (i)<defaultGreenTable[0].length)				return defaultGreenTable[maxState][i];			else				return defaultGreenTable[maxNumStates-1][maxNumStates-1];					}		else if (mode == BLUESCALE){			if ((maxState)<defaultBlueTable.length && (i)<defaultBlueTable[0].length)				return defaultBlueTable[maxState][i];			else				return defaultBlueTable[maxNumStates-1][maxNumStates-1];					}		else if (mode == REDSCALE){			if ((maxState)<defaultRedTable.length && (i)<defaultRedTable[0].length)				return defaultRedTable[maxState][i];			else				return defaultRedTable[maxNumStates-1][maxNumStates-1];					}		else if (mode == COLORS_NO_BW && (maxState+2)<defaultColorTable.length && (i+1)<defaultColorTable[0].length)			return defaultColorTable[maxState+2][i+1]; //TODO: guard that in bounds!		else if ((maxState)<defaultColorTable.length && (i)<defaultColorTable[0].length)			return defaultColorTable[maxState][i];		else 			return defaultColorTable[maxNumStates-1][maxNumStates-1];	}	/** gets color for state i with given maximum state possible, from color table*/	public Color getColor(int maxState, int i) {		if (maxState<0 || i<0)			return Color.white;		else if (colorTable!=null) //TODO: guard that in bounds!			return colorTable[maxState][i];		else			return defaultColorTable[maxState][i];	}	/** gets color for double d given range*/	public static Color getColor(double state, double min, double max) {		if (state>= min && state <= max && min !=max) {			int i =(int)(1000*(state-min)/(max-min));			if (i==1000)				i=999;			return defaultBIGColorTable[i];		}		else			return Color.white;	}	/** gets green value for double "state" given range*/	public static Color getGreenScale(double state, double min, double max, boolean log, int power) {		if (state>= min && state <= max && min !=max) {			double fraction = (state-min)/(max-min);   // the fraction, with 1 being black, 0 white			if (power>1) {				double origFraction=fraction;				for (int i = 1; i<=3; i++) {					fraction= fraction*origFraction;				}			}			if (log) {				fraction = (Math.exp(fraction)-1)/(Math.exp(1)-1);   			}			fraction = 1.0-fraction;   // the fraction, with 0 being black, 1 white			return new Color((float)fraction, (float)1, (float)fraction);		}		else			return Color.white;	}	public static Color getGreenScale(double state, double min, double max, boolean log) {		return getGreenScale(state,min,max,log,1);	}	/** gets green value for double "state" given range*/	public static Color getYellowScale(double state, double min, double max, boolean log, int power) {		if (state>= min && state <= max && min !=max) {			double fraction = (state-min)/(max-min);   // the fraction, with 1 being black, 0 white			if (power>1) {				double origFraction=fraction;				for (int i = 1; i<=3; i++) {					fraction= fraction*origFraction;				}			}			if (log) {				fraction = (Math.exp(fraction)-1)/(Math.exp(1)-1);   			}			fraction = 1.0-fraction;   // the fraction, with 0 being black, 1 white			return new Color((float)1, (float)1, (float)fraction);		}		else			return Color.white;	}	public static Color getYellowScale(double state, double min, double max, boolean log) {		return getYellowScale(state,min,max,log,1);	}	/** gets red value for double "state" given range*/	public static Color getRedScale(double state, double min, double max, boolean log, int power) {		if (state>= min && state <= max && min !=max) {			double fraction = (state-min)/(max-min);   // the fraction, with 1 being black, 0 white			if (power>1) {				double origFraction=fraction;				for (int i = 1; i<=3; i++) {					fraction= fraction*origFraction;				}			}			if (log) {				fraction = (Math.exp(fraction)-1)/(Math.exp(1)-1);   			}			fraction = 1.0-fraction;   // the fraction, with 0 being black, 1 white			return new Color((float)1, (float)fraction, (float)fraction);		}		else			return Color.white;	}	public static Color getRedScale(double state, double min, double max, boolean log) {		return getRedScale(state,min,max,log,1);	}	/** gets blue value for double "state" given range*/	public static Color getBlueScale(double state, double min, double max, boolean log, int power) {		if (state>= min && state <= max && min !=max) {			double fraction = (state-min)/(max-min);   // the fraction, with 1 being black, 0 white			if (power>1) {				double origFraction=fraction;				for (int i = 1; i<=3; i++) {					fraction= fraction*origFraction;				}			}			if (log) {				fraction = (Math.exp(fraction)-1)/(Math.exp(1)-1);   			}			fraction = 1.0-fraction;   // the fraction, with 0 being black, 1 white			return new Color((float)fraction, (float)fraction, (float)1);		}		else			return Color.white;	}	public static Color getBlueScale(double state, double min, double max, boolean log) {		return getBlueScale(state,min,max,log,1);	}	/** gets gray value for double "state" given range*/	public static Color getGrayScale(double state, double min, double max, boolean log, int power) {		if (state>= min && state <= max && min !=max) {			double fraction = (state-min)/(max-min);   // the fraction, with 1 being black, 0 white			if (power>1) {				double origFraction=fraction;				for (int i = 1; i<=3; i++) {					fraction= fraction*origFraction;				}			}			if (log) {				fraction = (Math.exp(fraction)-1)/(Math.exp(1)-1);   			}			fraction = 1.0-fraction;   // the fraction, with 0 being black, 1 white			return new Color((float)fraction, (float)fraction, (float)fraction);		}		else			return Color.white;	}	public static Color getGrayScale(double state, double min, double max, boolean log) {		return getGrayScale(state,min,max,log,1);	}	public static Color getGrayScale(double state, double min, double max) {		return getGrayScale(state,min,max,false,1);	}	/** sets color in color table */	public void setColor(int maxState, int i, Color color) {		if (colorTable!=null)			colorTable[maxState][i] = color;	}}