/* Mesquite source code.  Copyright 1997-2006 W. Maddison and D. Maddison.Version 1.1, May 2006.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.trees.BallsNSticks;import java.util.*;import java.awt.*;import mesquite.lib.*;import mesquite.lib.duties.*;import java.awt.geom.*;/* ======================================================================== */public class BallsNSticks extends DrawTree {	NodeLocsVH nodeLocsTask;	MesquiteCommand edgeWidthCommand;	MesquiteString orientationName, lineStyleName;	Vector drawings;	int oldEdgeWidth = 2;	int oldSpotSize = 22;	int ornt, style;	static final int DIAGONAL = 0;	static final int SQUARE = 1;	static final int CURVED = 2;	public MesquiteBoolean cosmic = new MesquiteBoolean(false);//	public MesquiteBoolean useArc = new MesquiteBoolean(false);	MesquiteSubmenuSpec orientationSubmenu;	MesquiteSubmenuSpec lineStyleSubmenu;	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) {		nodeLocsTask= (NodeLocsVH)hireEmployee(commandRec, NodeLocsVH.class, "Calculator of node locations");		if (nodeLocsTask == null)			return sorry(commandRec, getName() + " couldn't start because no node location module was obtained.");		drawings = new Vector();		ornt = TreeDisplay.UP; 		//addMenuItem("Display", "Edge width...", makeCommand("edgeWidth")); 		defineMenus(false);		orientationName = new MesquiteString("Up");		orientationSubmenu.setSelected(orientationName);		lineStyleName = new MesquiteString("Diagonal");		style = DIAGONAL;		lineStyleSubmenu.setSelected(lineStyleName);		return true; 	 }  	   	 public void employeeQuit(MesquiteModule m){  	 	iQuit();  	 }	public void defineMenus(boolean accumulating){  		orientationSubmenu = addSubmenu(null, "Orientation");		addItemToSubmenu(null, orientationSubmenu, "Up", makeCommand("orientUp",  this)); 		addItemToSubmenu(null, orientationSubmenu, "Right", makeCommand("orientRight",  this)); 		addItemToSubmenu(null, orientationSubmenu, "Down", makeCommand("orientDown",  this)); 		addItemToSubmenu(null, orientationSubmenu, "Left", makeCommand("orientLeft",  this)); 		lineStyleSubmenu = addSubmenu(null, "Line Style");		addItemToSubmenu(null, lineStyleSubmenu, "Diagonal", makeCommand("useDiagonal",  this)); 		addItemToSubmenu(null, lineStyleSubmenu, "Square", makeCommand("useSquare",  this)); 		addItemToSubmenu(null, lineStyleSubmenu, "Curved", makeCommand("useCurved",  this));//		addCheckMenuItem( null, "Curved Lines", makeCommand("toggleArc",  this), useArc); 		addMenuItem( "Line Width...", makeCommand("setEdgeWidth",  this)); 		addMenuItem( "Preferred Spot Size...", makeCommand("setSpotDiameter",  this));		addCheckMenuItem( null, "Cosmic", makeCommand("toggleCosmic",  this), cosmic);	}	public   TreeDrawing createTreeDrawing(TreeDisplay treeDisplay, int numTaxa) {		BallsNSticksDrawing treeDrawing =  new BallsNSticksDrawing (treeDisplay, numTaxa, this);		if (legalOrientation(treeDisplay.getOrientation())){			orientationName.setValue(orient(treeDisplay.getOrientation()));			ornt = treeDisplay.getOrientation();		}		else			treeDisplay.setOrientation(ornt);		drawings.addElement(treeDrawing);		return treeDrawing;	}	public boolean legalOrientation (int orientation){		return (orientation == TreeDisplay.UP || orientation == TreeDisplay.DOWN || orientation == TreeDisplay.RIGHT || orientation == TreeDisplay.LEFT);	}	/*................................................................................................................. 	public void endJob() { 		if (MesquiteTrunk.trackActivity) logln ("MesquiteModule " + getName() + "  closing down "); 		//dtd.disposePolys(coordTask.treeDisplay.getTree(), coordTask.treeDisplay.getTree().getRoot());		closeDownAllEmployees (this); 		employees.removeElement(nodeLocsTask); 		nodeLocsTask= null;   	 }   	/*.................................................................................................................*/	public String orient (int orientation){		if (orientation == TreeDisplay.UP)			return "Up";		else if (orientation == TreeDisplay.DOWN)			return "Down";		else if (orientation == TreeDisplay.RIGHT)			return "Right";		else if (orientation == TreeDisplay.LEFT)			return "Left";		else return "other";	}	/*.................................................................................................................*/  	 public Snapshot getSnapshot(MesquiteFile file) {    	 	Snapshot temp = new Snapshot();  	 	temp.addLine("setSpotDiameter " + oldSpotSize);   	 	temp.addLine("setEdgeWidth " + oldEdgeWidth);   	 	if (ornt== TreeDisplay.UP)  	 		temp.addLine("orientUp");   	 	else if (ornt== TreeDisplay.DOWN)  	 		temp.addLine("orientDown");   	 	else if (ornt== TreeDisplay.LEFT)  	 		temp.addLine("orientLeft");   	 	else if (ornt== TreeDisplay.RIGHT)  	 		temp.addLine("orientRight");   	 	if (style== DIAGONAL)  	 		temp.addLine("useDiagonal");   	 	else if (style== SQUARE)  	 		temp.addLine("useSquare");   	 	else if (style== CURVED)  	 		temp.addLine("useCurved"); 		temp.addLine("toggleCosmic " + cosmic.toOffOnString()); 	 	return temp;  	 }	MesquiteInteger pos = new MesquiteInteger();	/*.................................................................................................................*/	 public Object doCommand(String commandName, String arguments, CommandRecord commandRec, CommandChecker checker) {    	 	if (checker.compare(this.getClass(), "Sets the diameter of the spots", "[diameter in pixels]", commandName, "setSpotDiameter")) {			int newDiameter= MesquiteInteger.fromFirstToken(arguments, pos);			if (!MesquiteInteger.isCombinable(newDiameter))				newDiameter = MesquiteInteger.queryInteger(containerOfModule(), "Set spot diameter", "Enter preferred diameter of spots at nodes.  This sets the preferred spot size; if there is not room in the drawing for spots so large, then the actual spot size may be smaller.", oldSpotSize, 6, 100);    	 		if (newDiameter>=6 && newDiameter<100 && newDiameter!=oldSpotSize) {				Enumeration e = drawings.elements();				oldSpotSize = newDiameter;				while (e.hasMoreElements()) {					Object obj = e.nextElement();					BallsNSticksDrawing treeDrawing = (BallsNSticksDrawing)obj;    	 				treeDrawing.spotSize=newDiameter;					treeDrawing.preferredSpotSize = newDiameter;    	 				treeDrawing.treeDisplay.setMinimumTaxonNameDistance(0, treeDrawing.spotSize/2  + 4);    	 			}	 			parametersChanged(null, commandRec);    	 		}    	 		    	 	}    	 	else if (checker.compare(this.getClass(), "Sets how wide the branches of the tree are drawn", "[width in pixels]", commandName, "setEdgeWidth")) {			int newWidth= MesquiteInteger.fromFirstToken(arguments, pos);			if (!MesquiteInteger.isCombinable(newWidth))				newWidth = MesquiteInteger.queryInteger(containerOfModule(), "Set edge width", "Edge Width:", oldEdgeWidth, 1, 24);    	 		if (newWidth>0 && newWidth<24 && newWidth!=oldEdgeWidth) {				Enumeration e = drawings.elements();				while (e.hasMoreElements()) {					Object obj = e.nextElement();					BallsNSticksDrawing treeDrawing = (BallsNSticksDrawing)obj;    	 				treeDrawing.setEdgeWidth(newWidth);    	 			}    	 			oldEdgeWidth = newWidth;				if (commandRec != null && !commandRec.scripting()) parametersChanged(null, commandRec);    	 		}    	 	}    	 	else if (checker.compare(this.getClass(), "Sets whether or not \"cosmic\" mode is on", "[on or off]", commandName, "toggleCosmic")) {    	 		boolean current = cosmic.getValue();    	 		cosmic.toggleValue(parser.getFirstToken(arguments));    	 		if (current!=cosmic.getValue())    	 			parametersChanged(null, commandRec);    	 	}     	 	else if (checker.compare(this.getClass(), "Sets whether or not arcs are to be used", "[on or off]", commandName, "toggleArc")) {    	 		MesquiteBoolean useArc = new MesquiteBoolean(style == CURVED);  //here for compatibility with 1. 06 scripts    	 		boolean current = useArc.getValue();    	 		useArc.toggleValue(parser.getFirstToken(arguments));    	 		if (useArc.getValue()){    	 			style = CURVED;    	 		lineStyleName.setValue("Curved");     	 	}   	 		if (current!=useArc.getValue())    	 			parametersChanged(null, commandRec);    	 	}    	 	else if (checker.compare(this.getClass(), "Sets line style", null, commandName, "useDiagonal")) {    	 		int current = style;    	 		style = DIAGONAL;       			lineStyleName.setValue("Diagonal");   	 		if (current!=style)    	 			parametersChanged(null, commandRec);    	 	}       	 	else if (checker.compare(this.getClass(), "Sets line style", null, commandName, "useSquare")) {    	 		int current = style;    	 		style = SQUARE;       			lineStyleName.setValue("Square");   	 		if (current!=style)    	 			parametersChanged(null, commandRec);    	 	}       	 	else if (checker.compare(this.getClass(), "Sets line style", null, commandName, "useCurved")) {    	 		int current = style;    	 		style = CURVED;    			lineStyleName.setValue("Curved");    			if (current!=style)    	 			parametersChanged(null, commandRec);    	 	}   	 	else if (checker.compare(this.getClass(), "Returns the module employed to set the node locations", null, commandName, "getNodeLocsEmployee")) {    	 		return nodeLocsTask;    	 	}    	 	else if (checker.compare(this.getClass(), "Orients the tree so that the terminals are pointing up", null, commandName, "orientUp")) {			Enumeration e = drawings.elements();			ornt = 0;			while (e.hasMoreElements()) {				Object obj = e.nextElement();				BallsNSticksDrawing treeDrawing = (BallsNSticksDrawing)obj;		    	 	treeDrawing.reorient(TreeDisplay.UP);			    	 ornt = treeDrawing.treeDisplay.getOrientation();		    	 }			orientationName.setValue(orient(ornt));			parametersChanged(null, commandRec);    	 	}    	 	else if (checker.compare(this.getClass(), "Orients the tree so that the terminals are pointing down", null, commandName, "orientDown")) {			Enumeration e = drawings.elements();			ornt = 0;			while (e.hasMoreElements()) {				Object obj = e.nextElement();				BallsNSticksDrawing treeDrawing = (BallsNSticksDrawing)obj;		    	 	treeDrawing.reorient(TreeDisplay.DOWN);			    	 ornt = treeDrawing.treeDisplay.getOrientation();		    	 }			orientationName.setValue(orient(ornt)); 			parametersChanged(null, commandRec);   	 	}    	 	else if (checker.compare(this.getClass(), "Orients the tree drawing so that the terminal taxa are at right", null, commandName, "orientRight")) {			Enumeration e = drawings.elements();			ornt = 0;			while (e.hasMoreElements()) {				Object obj = e.nextElement();				BallsNSticksDrawing treeDrawing = (BallsNSticksDrawing)obj;		    	 	treeDrawing.reorient(TreeDisplay.RIGHT);			    	 ornt = treeDrawing.treeDisplay.getOrientation();		    	 }			orientationName.setValue(orient(ornt));			parametersChanged(null, commandRec);    	 	}    	 	else if (checker.compare(this.getClass(), "Orients the tree drawing so that the terminal taxa are at left", null, commandName, "orientLeft")) {			Enumeration e = drawings.elements();			ornt = 0;			while (e.hasMoreElements()) {				Object obj = e.nextElement();				BallsNSticksDrawing treeDrawing = (BallsNSticksDrawing)obj;		    	 	treeDrawing.reorient(TreeDisplay.LEFT);			    	 ornt = treeDrawing.treeDisplay.getOrientation();		    	 }			orientationName.setValue(orient(ornt));			parametersChanged(null, commandRec);    	 	} 		else { 			return super.doCommand(commandName, arguments, commandRec, checker);		}		return null;   	 }	/*.................................................................................................................*/    	 public String getName() {		return "Balls & Sticks";   	 }   	 	/*.................................................................................................................*/	/** returns whether this module is requesting to appear as a primary choice */   	public boolean requestPrimaryChoice(){   		return true;     	}	/*.................................................................................................................*/   	  	/** returns an explanation of what the module does.*/ 	public String getExplanation() { 		return "Draws trees with spots and the nodes and thin lines for branches." ;   	 }}/* ======================================================================== *//* ======================================================================== */class BallsNSticksDrawing extends TreeDrawing  {	public RotatedRectangle[] branchPoly;	public BallsNSticks ownerModule;	int spotSize;	int edgeWidth;	public int preferredSpotSize = 22;	int oldNumTaxa = 0; 	public static final int inset=1;	private boolean ready=false;	private int foundBranch;	NameReference triangleNameRef;	BasicStroke defaultStroke;		public BallsNSticksDrawing (TreeDisplay treeDisplay, int numTaxa, BallsNSticks ownerModule) {		super(treeDisplay, MesquiteTree.standardNumNodeSpaces(numTaxa));		triangleNameRef = NameReference.getNameReference("triangled");		this.ownerModule = ownerModule;		spotSize = ownerModule.oldSpotSize;		edgeWidth = ownerModule.oldEdgeWidth;		if (ownerModule.cosmic.getValue())			edgeWidth = 8;		try{			defaultStroke = new BasicStroke();		}		catch (Throwable t){		}	    	treeDisplay.setMinimumTaxonNameDistance(0, spotSize/2+ 4);		this.treeDisplay = treeDisplay;		oldNumTaxa = numTaxa;		ready = true;	}	public void resetNumNodes(int numNodes){		super.resetNumNodes(numNodes);		branchPoly= new RotatedRectangle[numNodes];		for (int i=0; i<numNodes; i++) {			branchPoly[i] = new RotatedRectangle();		}	}	/*_________________________________________________*/	private void UPdefinePoly(RotatedRectangle poly, boolean internalNode, int Nx, int Ny, int mNx, int mNy) {		poly.setShape(Nx, Ny, mNx, mNy, edgeWidth, Nx<mNx, RotatedRectangle.RECTANGLE); //FLATHORIZONTAL	}	/*_________________________________________________*/	private void UPCalcBranchPolys(Tree tree, int node)	{		if (!tree.getAssociatedBit(triangleNameRef,node))		for (int d = tree.firstDaughterOfNode(node); tree.nodeExists(d); d = tree.nextSisterOfNode(d))			UPCalcBranchPolys(tree, d);		UPdefinePoly(branchPoly[node], tree.nodeIsInternal(node), x[node],y[node], x[tree.motherOfNode(node)], y[tree.motherOfNode(node)]);	}	/*_________________________________________________*/	private void DOWNdefinePoly(RotatedRectangle poly, boolean internalNode, int Nx, int Ny, int mNx, int mNy) {		poly.setShape(Nx, Ny, mNx, mNy, edgeWidth, Nx>mNx, RotatedRectangle.RECTANGLE);	}	/*_________________________________________________*/	private void DOWNCalcBranchPolys(Tree tree, int node)	{		if (!tree.getAssociatedBit(triangleNameRef,node))			for (int d = tree.firstDaughterOfNode(node); tree.nodeExists(d); d = tree.nextSisterOfNode(d))					DOWNCalcBranchPolys(tree, d);		DOWNdefinePoly(branchPoly[node], tree.nodeIsInternal(node),x[node],y[node], x[tree.motherOfNode(node)], y[tree.motherOfNode(node)]);	}	/*_________________________________________________*/	/*_________________________________________________*/	private void RIGHTdefinePoly(RotatedRectangle poly, boolean internalNode, int Nx, int Ny, int mNx, int mNy) {		poly.setShape(Nx, Ny, mNx, mNy, edgeWidth, false, RotatedRectangle.RECTANGLE);	}	/*_________________________________________________*/	private void RIGHTCalcBranchPolys(Tree tree, int node)	{		if (!tree.getAssociatedBit(triangleNameRef,node))		for (int d = tree.firstDaughterOfNode(node); tree.nodeExists(d); d = tree.nextSisterOfNode(d))			RIGHTCalcBranchPolys(tree, d);		RIGHTdefinePoly(branchPoly[node], tree.nodeIsInternal(node),x[node],y[node], x[tree.motherOfNode(node)], y[tree.motherOfNode(node)]);	}	/*_________________________________________________*/	private void LEFTdefinePoly(RotatedRectangle poly, boolean internalNode, int Nx, int Ny, int mNx, int mNy) {		poly.setShape( Nx, Ny, mNx, mNy, edgeWidth, false, RotatedRectangle.RECTANGLE);	}	/*_________________________________________________*/	private void LEFTCalcBranchPolys(Tree tree, int node)	{		if (!tree.getAssociatedBit(triangleNameRef,node))		for (int d = tree.firstDaughterOfNode(node); tree.nodeExists(d); d = tree.nextSisterOfNode(d))			LEFTCalcBranchPolys(tree, d);		LEFTdefinePoly(branchPoly[node], tree.nodeIsInternal(node),x[node],y[node], x[tree.motherOfNode(node)], y[tree.motherOfNode(node)]);	}	/*_________________________________________________*/	private void calculateLines(Tree tree, int node) {		for (int d = tree.firstDaughterOfNode(node); tree.nodeExists(d); d = tree.nextSisterOfNode(d))			calculateLines( tree, d);		lineTipY[node]=y[node];		lineTipX[node]=x[node];		lineBaseY[node]=y[tree.motherOfNode(node)];		lineBaseX[node]=x[tree.motherOfNode(node)];	}	/*_________________________________________________*/	private void calcBranchPolys(Tree tree, int drawnRoot, CommandRecord commandRec) {		if (ownerModule==null) {MesquiteTrunk.mesquiteTrunk.logln("ownerModule null"); return;}		if (ownerModule.nodeLocsTask==null) {ownerModule.logln("nodelocs task null"); return;}		if (treeDisplay==null) {ownerModule.logln("treeDisplay null"); return;}		if (tree==null) { ownerModule.logln("tree null"); return;}				ownerModule.nodeLocsTask.calculateNodeLocs(treeDisplay,  tree, drawnRoot,  treeDisplay.getField(), commandRec); //Graphics g removed as parameter May 02		calculateLines(tree, drawnRoot);		spotSize = preferredSpotSize;		if (treeDisplay.getTaxonSpacing()<spotSize+2) {			spotSize= treeDisplay.getTaxonSpacing()-2;			if (spotSize<2)				spotSize=2;		}	    	treeDisplay.setMinimumTaxonNameDistance(0, spotSize/2+ 4);		if (treeDisplay.getTaxonSpacing()<edgeWidth+2) {			edgeWidth = treeDisplay.getTaxonSpacing()-2;			if (edgeWidth<2)				edgeWidth =2;		}		if (treeDisplay.getOrientation()==TreeDisplay.UP) {			UPCalcBranchPolys(tree, drawnRoot);		}		else if (treeDisplay.getOrientation()==TreeDisplay.DOWN){			DOWNCalcBranchPolys(tree, drawnRoot);		}		else  if (treeDisplay.getOrientation()==TreeDisplay.RIGHT) {			RIGHTCalcBranchPolys(tree, drawnRoot);		}		else  if (treeDisplay.getOrientation()==TreeDisplay.LEFT){			LEFTCalcBranchPolys(tree, drawnRoot);		}	}		/*_________________________________________________*/	/** Draw highlight for branch node with current color of graphics context */	public void drawHighlight(Tree tree, int node, Graphics g, boolean flip){		Color tC = g.getColor();		if (flip)			g.setColor(Color.red);		else			g.setColor(Color.blue);		for (int i=1; i<4; i++)			g.drawOval( x[node]- spotSize/2 - 2 - i, y[node]- spotSize/2 - 2 - i, spotSize + 3 + i + i, spotSize + 3 + i + i);		g.setColor(tC);	}	/*_________________________________________________*/	private void drawArc(Tree tree, Graphics g, int node, int width) {		if (tree.nodeExists(node)) {			int nM = tree.motherOfNode(node);			int xN=x[node];			int xnM = x[nM];			int yN =y[node];			int ynM = y[nM];			boolean done = false;			try{				if (MesquiteWindow.Java2Davailable && g instanceof Graphics2D) {					if (treeDisplay.getOrientation()==TreeDisplay.UP) {						if (xnM>xN){ //leans left							xN += width/2;							xnM += width/2;							ynM += edgeWidth - width/2;							yN += width/2;						}						else {							xN += width/2;							xnM += width/2;							ynM += width/2 ;							yN += width/2;						}											}					else if (treeDisplay.getOrientation()==TreeDisplay.DOWN){ //����						if (xnM>xN){ //leans left							xN += width/2;							xnM += width/2;							ynM -= edgeWidth - width/2;							yN -= width/2;						}						else {							xN += width/2;							xnM += width/2;							ynM -= width/2 ;							yN -= width/2;						}					}					else  if (treeDisplay.getOrientation()==TreeDisplay.RIGHT) {						if (ynM>yN){ //leans left							yN += width/2;							ynM += width/2;							xnM -= edgeWidth - width/2;							xN -= width/2;						}						else {							yN += width/2;							ynM += width/2;							xnM -= width/2 ;							xN -= width/2;						}											}					else  if (treeDisplay.getOrientation()==TreeDisplay.LEFT){  //����						if (ynM>yN){ //leans right							yN += width/2;							ynM += width/2;							xnM += edgeWidth - width/2;							xN += width/2;						}						else {							yN += width/2;							ynM += width/2;							xnM += width/2 ;							xN += width/2;						}					}					Arc2D.Double arc = null;					if (treeDisplay.getOrientation()==TreeDisplay.UP) {						if (xnM>xN) {  //leans left							//g.setColor(Color.blue);							arc = new Arc2D.Double(xN, yN-(ynM-yN), (xnM-xN)*2,  (ynM - yN)*2, 180, 90, Arc2D.OPEN); // left							//g.drawRect(xN, yN-(ynM-yN), (xnM-xN)*2,  (ynM - yN)*2);						}						else {							//g.setColor(Color.green);							arc = new Arc2D.Double(xnM-(xN-xnM), yN - (ynM - yN), (xN-xnM)*2,  (ynM - yN)*2, 0, -90, Arc2D.OPEN); //right							//g.drawRect(xnM-(xN-xnM), yN - (ynM - yN), (xN-xnM)*2,  (ynM - yN)*2);						}					}					else if (treeDisplay.getOrientation()==TreeDisplay.DOWN){//����						if (xnM>xN) {  //leans right							//g.setColor(Color.blue);							arc = new Arc2D.Double(xN, ynM, (xnM-xN)*2,  -(ynM - yN)*2, 90, 90, Arc2D.OPEN); // left							//g.drawRect(xN, yN-(ynM-yN), (xnM-xN)*2,  (ynM - yN)*2);						}						else {							//g.setColor(Color.green);							arc = new Arc2D.Double(xnM-(xN-xnM), ynM, (xN-xnM)*2,  -(ynM - yN)*2, 0, 90, Arc2D.OPEN); //right							//g.drawRect(xnM-(xN-xnM), yN - (ynM - yN), (xN-xnM)*2,  (ynM - yN)*2);						}					}					else  if (treeDisplay.getOrientation()==TreeDisplay.RIGHT) {						if (ynM>yN) { //leans left							//g.setColor(Color.blue);							arc = new Arc2D.Double(xnM, yN, (xN-xnM)*2,  (ynM - yN)*2, 90, 90, Arc2D.OPEN); // left							//g.drawRect(xN, yN-(ynM-yN), (xnM-xN)*2,  (ynM - yN)*2);						}						else {							//g.setColor(Color.green);							arc = new Arc2D.Double(xnM, ynM + (ynM - yN), (xN-xnM)*2,  -(ynM - yN)*2, 180,90, Arc2D.OPEN); //right							//g.drawRect(xnM-(xN-xnM), yN - (ynM - yN), (xN-xnM)*2,  (ynM - yN)*2);						}					}					else  if (treeDisplay.getOrientation()==TreeDisplay.LEFT){ //����						if (ynM>yN) { //leans right							//g.setColor(Color.blue);							arc = new Arc2D.Double(xN - (xnM-xN), yN, -(xN-xnM)*2,  (ynM - yN)*2, 0, 90, Arc2D.OPEN); 							//g.drawRect(xN, yN-(ynM-yN), (xnM-xN)*2,  (ynM - yN)*2);						}						else {							//g.setColor(Color.green);							arc = new Arc2D.Double(xN - (xnM-xN), ynM + (ynM - yN), -(xN-xnM)*2, - (ynM - yN)*2, 0,-90, Arc2D.OPEN); 							//g.drawRect(xnM-(xN-xnM), yN - (ynM - yN), (xN-xnM)*2,  (ynM - yN)*2);						}					}					if (arc!=null) {						BasicStroke wideStroke = new BasicStroke(width);						Graphics2D g2 = (Graphics2D)g;						g2.setStroke(wideStroke);						g2.draw(arc);						done  = true;						g2.setStroke(defaultStroke);					}				}								}			catch (Throwable t){			}			if (!done){				if (treeDisplay.getOrientation()==TreeDisplay.UP) {					if (xnM > xN)  ynM += edgeWidth-1;				}				else if (treeDisplay.getOrientation()==TreeDisplay.DOWN){ //����					if (xnM > xN)  ynM -= edgeWidth-1;				}				else  if (treeDisplay.getOrientation()==TreeDisplay.RIGHT) {					if (ynM > yN)  xnM -= edgeWidth-1;				}				else  if (treeDisplay.getOrientation()==TreeDisplay.LEFT){  //����					if (ynM > yN) xnM += edgeWidth-1;				}				else					System.out.println("Error: wrong tree orientation in BallsNSticks");				for (int i=0; i<width; i++) {					if (treeDisplay.getOrientation()==TreeDisplay.UP) {						if (xnM>xN) {							g.drawArc(xN, yN - (ynM - yN), (xnM-xN)*2,  (ynM - yN)*2, 180, 90); // left							ynM--;						}						else {							g.drawArc(xnM-(xN-xnM), yN - (ynM - yN), (xN-xnM)*2,  (ynM - yN)*2, 0, -90); //right							ynM++; 						}						xN++;					}					else if (treeDisplay.getOrientation()==TreeDisplay.DOWN){//����						if (xnM>xN) {							g.drawArc(xN,ynM, (xnM-xN)*2,  (yN -ynM)*2, 90, 90); //right							ynM++;						}						else {							g.drawArc(xnM-(xN-xnM),ynM, (xN-xnM)*2,   (yN -ynM)*2, 0, 90); //left 							ynM--;  						}						xN++;					}					else  if (treeDisplay.getOrientation()==TreeDisplay.RIGHT) {						if (ynM>yN) {							g.drawArc(xnM, yN, (xN-xnM)*2,  (ynM - yN)*2, 90, 90);  //left							xnM++;						}						else {							g.drawArc(xnM,ynM - (yN -ynM), (xN-xnM)*2,  (yN -ynM)*2, 180, 90);  //right 							xnM--;  						}						yN++;					}					else  if (treeDisplay.getOrientation()==TreeDisplay.LEFT){ //����						if (ynM>yN) {							g.drawArc(xN - (xnM-xN), yN, (xnM-xN)*2,  (ynM - yN)*2, 0, 90);  //right							xnM--;						}						else {							g.drawArc(xN - (xnM-xN),ynM - (yN -ynM), (xnM-xN)*2,  (yN -ynM)*2, 0, -90);  //left 							xnM++;  						}						yN++;					}									}			}		}	}	/*_________________________________________________*/	private void drawSquare(Tree tree, Graphics g, int node, int width) {		if (tree.nodeExists(node)) {			int nM = tree.motherOfNode(node);			int xN=x[node];			int xnM = x[nM];			int yN =y[node];			int ynM = y[nM];			if (treeDisplay.getOrientation()==TreeDisplay.UP) {				if (xnM>xN){ //leans left					xN += width/2;					xnM += width/2;					ynM += edgeWidth - width/2;					yN += width/2;					g.fillRect(xN, ynM, xnM-xN, edgeWidth);				}				else {					xN += width/2;					xnM += width/2;					ynM += width/2 ;					yN += width/2;					g.fillRect(xnM, ynM, xN - xnM+edgeWidth, edgeWidth);				}				g.fillRect(xN, yN, edgeWidth, ynM-yN);			}			else if (treeDisplay.getOrientation()==TreeDisplay.DOWN){ //����				if (xnM>xN){ //leans left					xN += width/2;					xnM += width/2;					ynM -= edgeWidth - width/2;					yN -= width/2;					g.fillRect(xN, ynM, xnM-xN, edgeWidth);				}				else {					xN += width/2;					xnM += width/2;					ynM -= width/2 ;					yN -= width/2;					g.fillRect(xnM, ynM, xN - xnM+edgeWidth, edgeWidth);				}				g.fillRect(xN, ynM, edgeWidth, yN-ynM);			}			else  if (treeDisplay.getOrientation()==TreeDisplay.RIGHT) {				if (ynM>yN){ //leans left					yN += width/2;					ynM += width/2;					xnM -= edgeWidth - width/2;					xN -= width/2;					g.fillRect(xnM, yN, edgeWidth, ynM-yN);				}				else {					yN += width/2;					ynM += width/2;					xnM -= width/2 ;					xN -= width/2;					g.fillRect(xnM, ynM, edgeWidth, yN - ynM);				}				g.fillRect(xnM, yN, xN-xnM, edgeWidth);			}			else  if (treeDisplay.getOrientation()==TreeDisplay.LEFT){  //����				if (ynM>yN){ //leans right					yN += width/2;					ynM += width/2;					xnM += edgeWidth - width/2;					xN += width/2;					g.fillRect(xnM, yN, edgeWidth, ynM-yN);				}				else {					yN += width/2;					ynM += width/2;					xnM += width/2 ;					xN += width/2;					g.fillRect(xnM, ynM, edgeWidth, yN - ynM);				}				g.fillRect(xN, yN, xnM - xN+ edgeWidth, edgeWidth);			}		}		}	/*_________________________________________________*/	private   void drawOneBranch(Tree tree, Graphics g, int node, int drawnRoot) {		if (tree.nodeExists(node)) {			if (!tree.getAssociatedBit(triangleNameRef,node))				for (int d = tree.firstDaughterOfNode(node); tree.nodeExists(d); d = tree.nextSisterOfNode(d))						drawOneBranch( tree, g, d, drawnRoot);			if (tree.getRooted() || tree.getRoot()!=node) {				g.setColor(treeDisplay.getBranchColor(node));				if (ownerModule.style == BallsNSticks.SQUARE)					drawSquare(tree, g, node, edgeWidth);								else if (ownerModule.style == BallsNSticks.CURVED)					drawArc(tree, g, node, edgeWidth);								else {					branchPoly[node].fill(g, ownerModule.cosmic.getValue()); 					// if (drawnRoot==node && tree.getRoot()!=node)					if (tree.numberOfParentsOfNode(node)>1) { //for reticulate trees						for (int i=1; i<=tree.numberOfParentsOfNode(node); i++) {							int anc =tree.parentOfNode(node, i);							if (anc!= tree.motherOfNode(node)) {								g.drawLine(x[node],y[node], x[tree.parentOfNode(node, i)],y[tree.parentOfNode(node, i)]);								g.drawLine(x[node]+1,y[node], x[tree.parentOfNode(node, i)]+1,y[tree.parentOfNode(node, i)]);								g.drawLine(x[node],y[node]+1, x[tree.parentOfNode(node, i)],y[tree.parentOfNode(node, i)]+1);								g.drawLine(x[node]+1,y[node]+1, x[tree.parentOfNode(node, i)]+1,y[tree.parentOfNode(node, i)]+1);							}						}					}				}				if (tree.getAssociatedBit(triangleNameRef,node)) {					for (int j=0; j<2; j++)					for (int i=0; i<2; i++) {						g.drawLine(x[node]+i,y[node]+j, x[tree.leftmostTerminalOfNode(node)]+i,y[tree.leftmostTerminalOfNode(node)]+j);						g.drawLine(x[tree.leftmostTerminalOfNode(node)]+i,y[tree.leftmostTerminalOfNode(node)]+j, x[tree.rightmostTerminalOfNode(node)]+i,y[tree.rightmostTerminalOfNode(node)]+j);						g.drawLine(x[node]+i,y[node]+j, x[tree.rightmostTerminalOfNode(node)]+i,y[tree.rightmostTerminalOfNode(node)]+j);					}				}				drawSpot( g, node);			}		}	}	/*_________________________________________________*/	public   void drawTree(Tree tree, int drawnRoot, Graphics g) {	        if (MesquiteTree.OK(tree)) {	        	if (tree.getNumNodeSpaces()!=numNodes)	        		resetNumNodes(tree.getNumNodeSpaces());	        	g.setColor(treeDisplay.branchColor);	       	 	drawOneBranch(tree, g, drawnRoot, drawnRoot);  	       	 }	   }	/*_________________________________________________*/	public   void recalculatePositions(Tree tree, CommandRecord commandRec) {	        if (MesquiteTree.OK(tree)) {	        	if (tree.getNumNodeSpaces()!=numNodes)	        		resetNumNodes(tree.getNumNodeSpaces());	        	if (!tree.nodeExists(getDrawnRoot()))	        		setDrawnRoot(tree.getRoot());	        	calcBranchPolys(tree, getDrawnRoot(), commandRec);		}	}		/*_________________________________________________*/	private boolean inSpot(int node, int h, int v){		if ((h-x[node])*(h-x[node]) + (v-y[node])*(v-y[node]) < spotSize*spotSize/4) //use radius			return true;		else			return false;		// ask if x, y is in node's spot    g.fillOval( x[node]- spotSize/2, y[node]- spotSize/2, spotSize, spotSize);	}	/*_________________________________________________*/	private void drawSpot(Graphics g, int node){		GraphicsUtil.fillOval(g, x[node]- spotSize/2, y[node]- spotSize/2, spotSize, spotSize, ownerModule.cosmic.getValue());	}	/*_________________________________________________*/	private void fillSpot(Graphics g, int node){		GraphicsUtil.fillOval(g, x[node]- spotSize/2 + 2, y[node]- spotSize/2 + 2, spotSize - 4, spotSize - 4, ownerModule.cosmic.getValue());	}	/*_________________________________________________*/	public  void fillTerminalBox(Tree tree, int node, Graphics g) {	}	/*_________________________________________________*/	public  void fillTerminalBoxWithColors(Tree tree, int node, ColorDistribution colors, Graphics g){	}	/*_________________________________________________*/	public  int findTerminalBox(Tree tree, int drawnRoot, int x, int y){		return -1;	}	private boolean ancestorIsTriangled(Tree tree, int node) {		if (tree.getAssociatedBit(triangleNameRef, tree.motherOfNode(node)))			return true;		if (tree.getRoot() == node || tree.getSubRoot() == node)			return false;		return ancestorIsTriangled(tree, tree.motherOfNode(node));	}		/*_________________________________________________*/	public void fillBranchWithColors(Tree tree, int node, ColorDistribution colors, Graphics g) {		if (node>0 && (tree.getRooted() || tree.getRoot()!=node) && !ancestorIsTriangled(tree, node)) {			Color c = g.getColor();			int numColors = colors.getNumColors();			if (numColors==1){				g.setColor(colors.getColor(0, !tree.anySelected()|| tree.getSelected(node)));				fillSpot(g,node);			}			else if (numColors>0) {				int startAngle=90;//was 270				double totalFreq=0;				for (int i=0; i<numColors; i++) totalFreq += colors.getWeight(i);				int arcAngle = 0;				for (int i=0; i<numColors; i++) {					Color color;					if ((color = colors.getColor(i, !tree.anySelected()|| tree.getSelected(node)))!=null)						g.setColor(color);										arcAngle = (int)((colors.getWeight(i)/totalFreq)*360);					GraphicsUtil.fillArc(g, x[node]- spotSize/2 + 2, y[node]- spotSize/2 + 2, spotSize - 4, spotSize - 4, startAngle, arcAngle, ownerModule.cosmic.getValue());					startAngle+=arcAngle; 				}			}			if (c!=null) g.setColor(c);		}	}	/*_________________________________________________*/	public   void fillBranch(Tree tree, int node, Graphics g) {		if (node>0 && (tree.getRooted() || tree.getRoot()!=node) && !ancestorIsTriangled(tree, node)) {			fillSpot(g,node);		}	}	   	/*_________________________________________________*/	private void ScanBranches(Tree tree, int node, int x, int y)	{		if (foundBranch==0) {			if (branchPoly[node].contains(x, y) || inSpot(node, x, y))				foundBranch = node;			if (!tree.getAssociatedBit(triangleNameRef, node)) {				for (int d = tree.firstDaughterOfNode(node); tree.nodeExists(d); d = tree.nextSisterOfNode(d))					ScanBranches(tree, d, x, y);			}		}	}	/*_________________________________________________*/	public   int findBranch(Tree tree, int drawnRoot, int x, int y) {	        if (MesquiteTree.OK(tree) && ready) {	        	foundBranch=0;	       		 ScanBranches(tree, drawnRoot, x, y);	       		 if (foundBranch == tree.getRoot() && !tree.getRooted())	       		 	return 0;	       		 else	       		 return foundBranch;	       	}	       	return 0;	}		/*_________________________________________________*/	public void reorient(int orientation) {		treeDisplay.setOrientation(orientation);		treeDisplay.pleaseUpdate(true, null);	}	/*_________________________________________________*/	public void setEdgeWidth(int edw) {		edgeWidth = edw;	}	/*_________________________________________________*/	public   void dispose() { 		for (int i=0; i<numNodes; i++) {			branchPoly[i] = null;		}		super.dispose();	}}	