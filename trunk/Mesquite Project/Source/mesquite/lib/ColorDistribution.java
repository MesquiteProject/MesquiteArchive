/* Mesquite source code.  Copyright 1997-2007 W. Maddison and D. Maddison.Version 2.01, December 2007.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.lib;import java.awt.*;import mesquite.lib.duties.*;/* ======================================================================== *//**A class to indicate how something is to be colored, either an array of colors to be shown in equal pieces, or same with distribution of weights for pie slices.  A maximum number of MAXCOLORS colors may be shown.*/public class ColorDistribution {	Color[] colors, colorsDimmed;	int numColors=0;	double[] weights;	boolean sequential = false;	static final int MAXCOLORS = 64;		public static int numberOfRed = 5;	public static int numberOfGreen = 11;	public static int numberOfBlue = 14;	public static Color lightGreen, veryLightGreen, darkGreen, lightBlue, veryLightBlue, violetBlue, veryLightGray, veryVeryLightGray, veryVeryVeryLightGray, lightRed, darkRed, veryVeryLightGreen;	public static Color darkBrown, brown, straw, lightYellow, veryLightYellow, tabLineBrown, mesquiteBrown, darkMesquiteBrown, veryDarkMesquiteBrown, lightMesquiteBrown, brightMesquiteBrown;	public static Color uneditable;	public static Color unassigned;	public static Color inapplicable;		public static Color [] codPosMedium, codPosDark;	public static Color spinDark, spinLight;	public static Color[] projectLight, projectDark; //pale, light, medium, dark, project, 	public final static int numColorSchemes = 4;	public static Color burlyWood, navajoWhite, bisque, sienna, paleGoldenRod, veryPaleGoldenRod;	public static NameReference colorNameReference;	public static StringArray standardColorNames;	static ObjectArray standardColors, standardColorsDimmed;	public static double dimmingConstant = 0.3;	public static int NO_COLOR = 18;	static {		spinLight = new Color((float)0.3, (float)0.6, (float)0.99);		spinDark = new Color((float)0.1, (float)0.1, (float)0.70);		unassigned = new Color(230, 230, 230);		inapplicable = veryVeryLightGray;;				veryLightGray = brighter(Color.lightGray, 0.5);		veryVeryLightGray = brighter(veryLightGray, 0.5);		veryVeryVeryLightGray = new Color((float)0.98, (float)0.98, (float)0.98);		darkRed = new Color((float)0.5, (float)0.2, (float)0.1);		lightRed = new Color((float)0.9, (float)0.48, (float)0.35);		darkGreen = new Color((float)0.1, (float)0.5, (float)0.2);		lightGreen = new Color((float)0.35, (float)0.9, (float)0.48);		veryLightGreen = new Color((float)0.55, (float)0.99, (float)0.68);		veryVeryLightGreen = new Color((float)0.70, (float)0.99, (float)0.83);		lightBlue = new Color((float)0.35, (float)0.48, (float)0.9);		veryLightBlue = new Color((float)0.55, (float)0.68, (float)0.99);				violetBlue = new Color((float)0.55, (float)0.40, (float)0.89);		tabLineBrown = new Color((float)0.47, (float)0.41, (float)0.26);		darkBrown = new Color((float)0.45, (float)0.40, (float)0.15);		brown = new Color((float)0.65, (float)0.58, (float)0.25);		straw = new Color((float)0.85, (float)0.80, (float)0.38);		lightYellow = new Color((float)0.95, (float)0.95, (float)0.64);		veryLightYellow = new Color((float)0.99, (float)0.99, (float)0.78);				burlyWood = new Color((float)0.87, (float)0.7216, (float)0.5294); //222, 184, 135  DEB887; medium  SHOULD BE 0.87, (float)0.7216, (float)0.5294		navajoWhite =  new Color((float)1.0, (float)0.87, (float)0.6784); //FFDEAD; light		bisque =  new Color((float)1.0, (float)0.894, (float)0.7686); // FFE4C4; pale		sienna =  new Color((float)0.6275, (float)0.3216, (float)0.1765); // A0522D; dark		paleGoldenRod = new Color((float)0.9333, (float)0.9398, (float)0.6667); //EEE8AA  green should be 0.9333, 0.9098, 0.66667		veryPaleGoldenRod = brighter(paleGoldenRod,0.5);		lightMesquiteBrown = new Color(188, 168, 122);		brightMesquiteBrown = new Color(228, 200, 132);		mesquiteBrown = new Color(108, 98, 82);		darkMesquiteBrown = new Color(92, 82, 70);		veryDarkMesquiteBrown = new Color(78, 68, 55);		//spinLight = new Color((float)0.6, (float)0.9, (float)0.6);		//spinDark = new Color((float)0.05, (float)0.5, (float)0.05);//		darkMesquiteBrown = new Color(82, 72, 60);		//mesquiteBrown = new Color(88, 78, 62);		//mesquiteBrown = new Color(77, 65, 47);		codPosMedium = new Color[4];		codPosDark = new Color[4];		codPosDark[0] = new Color((float)0.1, (float)0.2, (float)0.5);   // first positions, blue		//codPosMedium[0] = new Color((float)0.35, (float)0.48, (float)0.9);		codPosMedium[0] = new Color((float)0.45, (float)0.55, (float)0.94);		codPosDark[1] = new Color((float)0.1, (float)0.5, (float)0.2);   // second positions, green		codPosMedium[1] = new Color((float)0.35, (float)0.9, (float)0.48);		codPosDark[2] = new Color((float)0.5, (float)0.2, (float)0.1);   // third positions, red		codPosMedium[2] = new Color((float)0.9, (float)0.48, (float)0.35);		codPosDark[3] = Color.gray;       // noncoding, unspecified		codPosMedium[3] = Color.lightGray;		/* traditional theme		pale = veryLightGray;		light = Color.lightGray;		project = Color.lightGray;		medium = Color.gray;		dark = Color.darkGray;		/**/		uneditable = lightYellow;		/* autumn theme*/		/*pale = new Color[numColorSchemes];		light = new Color[numColorSchemes];		medium = new Color[numColorSchemes];		dark = new Color[numColorSchemes];		project = new Color[numColorSchemes];*/		projectLight = new Color[numColorSchemes];		projectDark = new Color[numColorSchemes];		/*pale[0] = bisque; 		light[0] = paleGoldenRod;  		medium[0] = burlyWood; 		project[0] = burlyWood;		dark[0] = sienna;  /**/		projectLight[0] = paleGoldenRod;		projectDark[0] = sienna;		colorNameReference = NameReference.getNameReference("Color");		standardColors = new ObjectArray(18);		standardColors.setValue(0, Color.black);		standardColors.setValue(1, Color.darkGray);		standardColors.setValue(2, Color.gray);		standardColors.setValue(3, Color.lightGray);		standardColors.setValue(4, Color.white);		standardColors.setValue(5, Color.red); //	public static int numberOfRed = 5;		standardColors.setValue(6, Color.orange);		standardColors.setValue(7, Color.yellow);		standardColors.setValue(8, paleGoldenRod);		standardColors.setValue(9, burlyWood);		standardColors.setValue(10, sienna);		standardColors.setValue(11, Color.green);//	public static int numberOfGreen = 11;		standardColors.setValue(12, lightGreen);		standardColors.setValue(13, Color.cyan);		standardColors.setValue(14, Color.blue); // public static int numberOfBlue = 14;		standardColors.setValue(15, lightBlue);		standardColors.setValue(16, Color.magenta);		standardColors.setValue(17, Color.pink);		//DO NOT ASSIGN A COLOR TO 18		standardColorNames = new StringArray(18);		standardColorNames.setValue(0, "Black");		standardColorNames.setValue(1, "Dark Gray");		standardColorNames.setValue(2, "Gray");		standardColorNames.setValue(3, "Light Gray");		standardColorNames.setValue(4, "White");		standardColorNames.setValue(5, "Red");		standardColorNames.setValue(6, "Orange");		standardColorNames.setValue(7, "Yellow");		standardColorNames.setValue(8, "Goldenrod");		standardColorNames.setValue(9, "Wood");		standardColorNames.setValue(10, "Sienna");		standardColorNames.setValue(11, "Green");		standardColorNames.setValue(12, "Light Green");		standardColorNames.setValue(13, "Cyan");		standardColorNames.setValue(14, "Blue");		standardColorNames.setValue(15, "Light Blue");		standardColorNames.setValue(16, "Magenta");		standardColorNames.setValue(17, "Pink");		//DO NOT ASSIGN A COLOR TO 18		standardColorsDimmed = new ObjectArray(18);		standardColorsDimmed.setValue(0, Color.gray);		standardColorsDimmed.setValue(1, brighter(Color.darkGray, dimmingConstant));		standardColorsDimmed.setValue(2, brighter(Color.gray, dimmingConstant));		standardColorsDimmed.setValue(3, brighter(Color.lightGray, dimmingConstant));		standardColorsDimmed.setValue(4, Color.white);		standardColorsDimmed.setValue(5, brighter(Color.red, dimmingConstant));		standardColorsDimmed.setValue(6, brighter(Color.orange, dimmingConstant));		standardColorsDimmed.setValue(7, brighter(Color.yellow, dimmingConstant));		standardColorsDimmed.setValue(8, brighter(paleGoldenRod, dimmingConstant));		standardColorsDimmed.setValue(9, brighter(burlyWood, dimmingConstant));		standardColorsDimmed.setValue(10, brighter(sienna, dimmingConstant));		standardColorsDimmed.setValue(11, brighter(Color.green, dimmingConstant));		standardColorsDimmed.setValue(12, brighter(lightGreen, dimmingConstant));		standardColorsDimmed.setValue(13, brighter(Color.cyan, dimmingConstant));		standardColorsDimmed.setValue(14, brighter(Color.blue, dimmingConstant));		standardColorsDimmed.setValue(15, brighter(lightBlue, dimmingConstant));		standardColorsDimmed.setValue(16, brighter(Color.magenta, dimmingConstant));		standardColorsDimmed.setValue(17, brighter(Color.pink, dimmingConstant));		//DO NOT ASSIGN A COLOR TO 18	}	public ColorDistribution() {		colors = new Color[MAXCOLORS];		colorsDimmed = new Color[MAXCOLORS];		weights = new double[MAXCOLORS];			}		public static Color getContentBackground(){		return new Color(235,235,235);	}	public static Color getContentBackgroundPale(){		return new Color(240,240,240);	}	public static Color getContentElement(){		return new Color(225,225,225);	}	public static Color getContentDarkElement(){		return new Color(150,150,150);	}	public static Color getContentEdgePale(){		return new Color(240,240,240);	}	public static Color getContentEdgeDark(){		return new Color(150,150,150);	}	//THESE SHOULD ALL BE REPLACED BY references to static colors, rather than having to instantiate a colour each time		//EXTERNAL INTERFACE AREA (main tab panel, project resources panel)	public static Color getExtInterfaceBackground(){  //general background to main tabs and project panel		return new Color(74, 68, 59);//chocolate theme		//return new Color(168,168,179); //slate theme	}	public static Color getExtInterfaceElement(){ // slightly contrasting color for unselected main tabs and some parts of project panel		return new Color(88, 82, 74);//chocolate theme	//return new Color(188,188,200);//slate theme	}	public static Color getExtInterfaceElementContrast(){   //element that is contrasted, e.g. the highlighted one, in the main tabs/project area		return new Color(102, 98, 86); //chocolate theme		//return new Color(221,221,232);//slate theme	}		public static Color getExtInterfaceEdgeContrast(){			return new Color(190,185,175);//chocolate theme			//return new Color(80,80,90);//slate theme	}	public static Color getExtInterfaceTextContrast(){  //high contrast text in external interface area, e.g. highlighted tab		return new Color(238,224,210);//chocolate theme		//return new Color(28,28,33); //slate theme	}	public static Color getExtInterfaceTextMedium(){ //medium contrast text in external interface area, e.g. project resources text		return new Color(210,190,154);//chocolate theme		//new Color(48,48,56); //slate theme	}	public static Color getExtInterfaceTextMuted(){  //muted text in external interface area,e.g. unselected tabs		return new Color(190,178,130);//chocolate theme	}		//INTERNAL INTERFACE AREA (tool palettes, graphics/text/etc tabs)	public static Color getInterfaceBackground(){  //background of tool palettes and small tabs bar		return new Color(216,204,172);  // chocolate theme  		//return new Color(196,192,192);  // slate theme  	}	public static Color getInterfaceElement(){ //unselected tool buttons and small tab		return new Color(234, 228,212); // chocolate theme 		//return new Color(216,212,212);  // slate theme 	}	public static Color getInterfaceElementContrast(){ //selected small tab		return new Color(236,230,222); // chocolate theme 		//return new Color(228,228,228);  // slate theme 	}	public static Color getInterfaceEdgeNegative(){  //edge to unselected small tab and tool button		return new Color(226,222,222); // chocolate theme }	public static Color getInterfaceEdgePositive(){  //edge to selected small tab		return new Color(114,100,90); // chocolate theme 		//return new Color(100,100,104);  // slate theme 	}	public static Color getInterfaceTextContrast(){  //contrasting text; must be opposite to element contrast; e.g. for selected tab		return new Color(44,40,40);  //chocolate theme		//return new Color(40,40,44);  // slate theme 	}	public static Color getInterfaceTextMuted(){  //muted text; e.g. for unselected tabs		return new Color(64,60,60); // chocolate theme 		//return new Color(60,60,64);  // slate theme 	}		public static Color getProjectLight(int i){  //currently used for selected tool button		return new Color(242, 224, 185);	}	public static Color getProjectDark(int i){ //currently used for selected tool button edge		return new Color(160, 82, 45);   //just brown	}			private static float brighten(int v, double percent){		float b = (float)((255-(255-v)*percent)/255);		if (b<0)			b=0;		else if (b>1)			b=1;		return b;	}		public static void setTransparentGraphics(Graphics g, float f) {		if (g!=null && (g instanceof Graphics2D)) {			if (f>0.0f && f<1.0f)				((Graphics2D)g).setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, f)); //			else if (f==0.0f)//				((Graphics2D)g).setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.00001f)); 		}	}	public static void setTransparentGraphics(Graphics g) {		setTransparentGraphics(g,0.5f); 	}	public static void setOpaqueGraphics(Graphics g) {		if (g!=null  && (g instanceof Graphics2D))			((Graphics2D)g).setComposite(AlphaComposite.getInstance(AlphaComposite.SRC, 1));   	}	public static Color brighter(Color c, double percent){		if (c==null)			return null;		int green = c.getGreen();		int red = c.getRed();		int blue = c.getBlue();		return new Color(brighten(red, percent), brighten(green, percent), brighten(blue, percent));	}	private static float darken(int v, double percent){		float b = (float)((v*percent)/255);		if (b<0)			b=0;		else if (b>1)			b=1;		return b;	}		public static Color darker(Color c, double percent){		if (c==null)			return null;		int green = c.getGreen();		int red = c.getRed();		int blue = c.getBlue();		return new Color(darken(red, percent), darken(green, percent), darken(blue, percent));	}	public static int getStandardColorNumber(String name){		int ci = standardColorNames.indexOf(name);		return ci;	}	public static int getStandardColorNumber(Color color){		int ci = standardColors.indexOf(color);		return ci;	}	public static Color getStandardColor(String name){		int ci = standardColorNames.indexOf(name);		if (ci<0)			return null;		return (Color)standardColors.getValue(ci);	}	public static String getStandardColorName(Color color){		if (color==null)			return null;		int ic = standardColors.indexOf(color);		if (ic>=0)			return standardColorNames.getValue(ic);		else			return null;	}	public static String getStandardColorName(int ci){		if (ci<0)			return null;		return (String)standardColorNames.getValue(ci);	}	public static Color getStandardColor(int ci){		if (ci<0)			return null;		return (Color)standardColors.getValue(ci);	}	public static Color getStandardColorDimmed(int ci){		if (ci<0)			return null;		return (Color)standardColorsDimmed.getValue(ci);	}	/** Initialize colors by setting the number of colors to 0, the weights to 0, and the colors to null*/	public void initialize() {		numColors=0;		for (int i=0; i<MAXCOLORS; i++) {			colors[i]=null;			colorsDimmed[i]=null;			weights[i]=0;		}	}	/** set color for state (or other unit) i to the given color.*/	public void setColor(int i, Color color) {		if (i>=0 && i<MAXCOLORS) {			if (colors[i]==null)				numColors++;			else if (color==null)				numColors--;			colors[i]=color;			Color dimmed = color;			if (color !=null){				if ( color.getGreen()==0 && color.getRed() ==0 && color.getBlue()==0)					dimmed = Color.gray;				else					dimmed = ColorDistribution.brighter(color, dimmingConstant);			}			colorsDimmed[i] = dimmed;		}	}	/** finds index of color.*/	public int indexOf(Color color) {		if (color == null)			return -1;		for (int i=0; i<colors.length; i++)			if (colors[i] != null && color.equals(colors[i])) 				return i;		return -1;	}		public static Color getColorFromArguments(String arguments, MesquiteInteger pos) {		int red =  MesquiteInteger.fromString(arguments, pos);		int green =  MesquiteInteger.fromString(arguments, pos);		int blue =  MesquiteInteger.fromString(arguments, pos);		return new Color(red,green,blue);	}		public static String getColorStringForSnapshot(Color color) {		return color.getRed() + " " + color.getGreen() + " " + color.getBlue();	}	/** gets the color for unit i*/	public Color getColor(int i) {		return colors[i];	}	/** gets the color for unit i*/	public Color getColor(int i, boolean regularStrength) {		if (regularStrength)			return colors[i];		else			return colorsDimmed[i];	}	/** sets the weight (for pie charts) for unit i*/	public void setWeight(int i, double weight) {		if (i>=0 && i<MAXCOLORS) {			weights[i]=weight;		}	}	/** get the weight for state (or other unit) i*/	public double getWeight(int i) {		return weights[i];	}	/** sets whether or not the colors are sequential, e.g. as for stochastic character mapping. 	 In this case the weights are the point at which the change occurs, e.g. proportional position on branch*/	public void setSequential(boolean s){		sequential = s;	}	/** returns whether the colors are to be considered sequential */	public boolean getSequential(){		return sequential;	}	/** get the total number of colors assigned*/	public int getNumColors() {		return numColors;	}	public static int getColorScheme(MesquiteModule mod){		if (mod == null || mod.getProject()==null)			return 0;		else			return mod.getProject().getProjectColor();	}	public static int getColorScheme(MesquiteProject mp){		if (mp == null)			return 0;		else			return mp.getProjectColor();	}	public String toString(){		String s = "ColorDistribution: ";		if (colors == null)			return s + " NO COLORS";		for (int i=0; i <colors.length && i<numColors; i++)			s += " " + i + ": " + colors[i] + " (weight " + weights[i] + ");";		return s;	}}