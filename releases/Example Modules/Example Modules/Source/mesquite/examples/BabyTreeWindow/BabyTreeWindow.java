package mesquite.examples.BabyTreeWindow;import java.awt.*;import mesquite.lib.*;import mesquite.lib.duties.*;/*========================================*//*A simple example module that displays a tree window that depends on a tree in a standard tree window*/public class BabyTreeWindow extends TreeWindowAssistantN {	DrawTreeCoordinator treeDrawCoordTask;	NumberForTree numberTask;	BTreeWindow bTreeWindow;	  	 /*--------------------------------------*/	/*The basic substitute for a constructor for modules  <b>(overrides method of MesquiteModule)</b>*/	public boolean startJob(String arguments, Object condition, boolean hiredByName) {		treeDrawCoordTask= (DrawTreeCoordinator)hireEmployee(DrawTreeCoordinator.class, null);		if (treeDrawCoordTask == null)			return sorry(getName() + " couldn't start because no tree draw coordinating module was obtained.");		numberTask= (NumberForTree)hireEmployee(NumberForTree.class, "Number to calculate"); 		makeMenu("Baby");		addSubmenu(null, "Number to calculate", makeCommand("setNumberTask",  this), NumberForTree.class);		addMenuItem( "-", null); 		bTreeWindow= new BTreeWindow( this, treeDrawCoordTask, numberTask); 		setModuleWindow(bTreeWindow); 		if (!MesquiteThread.isScripting()) 			bTreeWindow.setVisible(true); 		resetContainingMenuBar();		resetAllWindowsMenus(); 		return true;  	 }	/*.................................................................................................................*/	public boolean isSubstantive(){		return false;	}  	 /*--------------------------------------*/  	 public void employeeQuit(MesquiteModule m){  	 	if (m==numberTask) {  	 		numberTask = null;			bTreeWindow.setNumberTask(null);  	 		resetContainingMenuBar();  	 	}  	 	else  	 		iQuit();  	 }  	 /*--------------------------------------*/	/*A method necessary with modules of subclass TreeWindowAssistantN; allows standard tree window	on which this module depends to indicate tree has changed  <b>(overrides method of TreeWindowAssistantN)</b>*/	public   void setTree(Tree tree) {		bTreeWindow.setTree(tree);	}  	 /*--------------------------------------*/	/*Makes the module shut down when the go-away box of the window is touched  <b>(overrides method of MesquiteModule)</b>*/ 	public void windowGoAway(MesquiteWindow whichWindow) {		whichWindow.hide();		whichWindow.dispose();		iQuit();	}  	 /*--------------------------------------*/	/*Returns the snapshot necessary to get this module back to the current state.  Note that it incorporates	a snapshot from its window  <b>(overrides method of MesquiteModule)</b>*/  	 public Snapshot getSnapshot(MesquiteFile file) {  	 	if (bTreeWindow ==null)  	 		return null;  	 	Snapshot fromWindow = bTreeWindow.getSnapshot(file);  	 	if (fromWindow == null || fromWindow.getNumLines() ==0)  	 		return null;   	 	Snapshot sn = new Snapshot();		sn.addLine("getWindow");		sn.addLine("tell It");		sn.incorporate(fromWindow, true);		sn.addLine("endTell");		sn.addLine("getTreeDrawCoordinator", treeDrawCoordTask);		sn.addLine("setNumberTask", numberTask);    	 	sn.addLine("showWindow");  	 	return sn;  	 }  	 /*--------------------------------------*/	/*The standard method for Commandable interface; receives commands either for snapshotting purposes	or from menu actions  <b>(overrides method of MesquiteModule)</b>*/    	 public Object doCommand(String commandName, String arguments, CommandChecker checker) {    	 	if (checker.compare(this.getClass(), "Returns the module serving as the window's draw tree coordinator", null, commandName, "getTreeDrawCoordinator"))    	 		return treeDrawCoordTask;    	 	else if (checker.compare(this.getClass(), "Sets which module class should calculate a number for the tree", "[name of module]", commandName, "setNumberTask")) {    	 		NumberForTree temp= (NumberForTree)replaceEmployee(NumberForTree.class, arguments, null, numberTask);			if (temp!=null) {				numberTask = temp;				bTreeWindow.setNumberTask(numberTask);				resetContainingMenuBar();				return numberTask;			}    	 	}    	 	else    	 		return super.doCommand(commandName, arguments, checker);		return null;   	 }  	 /*--------------------------------------*/	/*Receives message from employees that their parameters have changed an recalculation may be needed  <b>(overrides method of MesquiteModule)</b>*/ 	public void employeeParametersChanged(MesquiteModule employee, MesquiteModule source, Notification notification) { 		if (employee== numberTask) 			bTreeWindow.recalculate(); 	}  	 /*--------------------------------------*/	/*Indicates to the name of this module for purposes of menu listings and documentation.  <b>(overrides method of MesquiteModule)</b>*/    	 public String getName() {		return "Baby Tree Window (example module)";   	 }  	 /*--------------------------------------*/ 	/*Returns an explanation of what the module does.  <b>(overrides method of MesquiteModule)</b>*/ 	public String getExplanation() { 		return "Displays a single tree (the same as in a tree window).  This is an example module to show how to program Mesquite." ;   	 }  	 /*--------------------------------------*/	/*Returns the authors of the module.  <b>(overrides method of MesquiteModule)</b>*/  	 public String getAuthors() {		return "Wayne Maddison & David Maddison";   	 }}	/*========================================*//*The window (Frame) itself shown on the screen, containing a tree that is the same as the one in the standard tree window*/class BTreeWindow extends MesquiteWindow  {	TreeDisplay treeDisplay;	DrawTreeCoordinator treeDrawCoordTask;	TextField p;	MesquiteNumber num = new MesquiteNumber();	MesquiteString resultString = new MesquiteString();	NumberForTree numberTask;	MesquiteTree tree;		public BTreeWindow (BabyTreeWindow ownerModule, DrawTreeCoordinator treeDrawCoordTask, NumberForTree numberTask){		super(ownerModule, true);   		this.treeDrawCoordTask = treeDrawCoordTask;  		this.numberTask = numberTask;      		setWindowSize(500,400);		setBackground(Color.white);		p = new TextField();		p.setBackground(Color.yellow);		addToWindow(p);		resetTitle();	}  	 /*--------------------------------------*/	/* Used to get the title for window <b>(overrides abstract method of MesquiteWindow)</b>*/	public void resetTitle(){		setTitle("Baby Tree Window"); 	}  	 /*--------------------------------------*/	/* Resize the tree display and other components.*/	public void sizeDisplays(){		if (treeDisplay==null) return;		int totalWidth = getWidth();		int totalHeight = getHeight() - 30;		treeDisplay.setSize(totalWidth,totalHeight);		treeDisplay.setFieldSize(totalWidth,totalHeight);		if (p==null)			return;		p.setBounds(0,totalHeight, totalWidth, 30);	}  	 /*--------------------------------------*/	/*Sets the tree to be shown in the window.*/	public void setTree(Tree newTree){		if (treeDisplay == null) {			Taxa taxa = newTree.getTaxa();			treeDisplay =treeDrawCoordTask.createOneTreeDisplay(taxa, this); 			addToWindow(treeDisplay);			treeDisplay.setLocation(0,0);			sizeDisplays();		}				if (treeDisplay.getTree()!=null)			treeDisplay.getTree().dispose();		if (newTree!=null) {			tree = newTree.cloneTree();			treeDisplay.setTree(tree);			recalculate();			treeDisplay.suppressDrawing(false);			treeDisplay.setVisible(true);			treeDisplay.repaint();		}	}  	 /*--------------------------------------*/	/*Sets what module is to be used for calculating the number for the tree*/	public void setNumberTask(NumberForTree numTask){		numberTask = numTask;		recalculate();	}  	 /*--------------------------------------*/	/*Recalculates the number for the tree*/	public void recalculate(){		if (numberTask!=null && tree !=null){			resultString.setValue("");			num.setToUnassigned();			numberTask.calculateNumber(tree, num, resultString);			p.setText(resultString.toString());		}		else {			p.setText("");		}	}  	 /*--------------------------------------*/	/*Sets the size of the window (setSize and setBounds should not be used!!!>  <b>(overrides method of MesquiteWindow)</b>*/	public void setWindowSize(int w, int h){		super.setWindowSize(w,h);		sizeDisplays();	}  	 /*--------------------------------------*/	/* Called when the window has been resized, e.g. by user. <b>(overrides method of MesquiteWindow)</b>*/	public void windowResized(){		sizeDisplays();	}  	 /*--------------------------------------*/	/*Disposes of the window*/	public void dispose(){		if (treeDisplay!=null){			if (treeDisplay.getTree()!=null)				treeDisplay.getTree().dispose();			treeDisplay.dispose();		}		super.dispose();	}}