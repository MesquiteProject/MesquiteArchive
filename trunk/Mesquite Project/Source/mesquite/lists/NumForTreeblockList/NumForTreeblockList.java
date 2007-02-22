/* Mesquite source code.  Copyright 1997-2006 W. Maddison and D. Maddison. Version 1.11, June 2006.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html) */package mesquite.lists.NumForTreeblockList;/*~~  */import mesquite.lists.lib.*;import java.util.*;import java.awt.*;import mesquite.lib.*;import mesquite.lib.duties.*;import mesquite.lib.table.*;/* ======================================================================== */public class NumForTreeblockList extends TreeblocksListAssistant implements MesquiteListener {	/*.................................................................................................................*/	public String getName() {		return "Number for Tree Block (in List of Tree Blocks window)";	}	public String getNameForMenuItem() {		return "Number for Tree Block";	}	public String getExplanation() {		return "Supplies numbers for tree blocks for a tree blocks list window." ;	}	public void getEmployeeNeeds(){  //This gets called on startup to harvest information; override this and inside, call registerEmployeeNeed		EmployeeNeed e = registerEmployeeNeed(NumberForTreeBlock.class, "The List of Tree Blocks window can display a column showing a value for each of the tree blocks.",		"You can select a value to show in the Number For Tree Blocks submenu of the Columns menu of the List of Tree Blocks Window. ");	}	/*.................................................................................................................*/	ListableVector treeBlocks=null;	NumberForTreeBlock numberTask;	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) {		if (arguments !=null) {			numberTask = (NumberForTreeBlock)hireNamedEmployee(commandRec, NumberForTreeBlock.class, arguments);			if (numberTask==null) {				return sorry(commandRec, "Number for tree block (for list) can't start because the requested module was not successfully hired");			}			return true;		}		numberTask = (NumberForTreeBlock)hireEmployee(commandRec, NumberForTreeBlock.class, "Value to calculate for tree block (for tree block list)");		if (numberTask==null) {			return sorry(commandRec, "Number for tree block (for list) can't start because the no calculating module was successfully hired");		}		return true;	}	/** Returns whether or not it's appropriate for an employer to hire more than one instance of this module.   	If false then is hired only once; second attempt fails.*/	public boolean canHireMoreThanOnce(){		return true;	}	/*.................................................................................................................*/	public void employeeQuit(MesquiteModule m){		iQuit();	}	public void setTableAndObject(MesquiteTable table, Object obj, CommandRecord commandRec){		if (treeBlocks !=null)			treeBlocks.removeListener(this);		if (obj instanceof ListableVector)			this.treeBlocks = (ListableVector)obj;		treeBlocks.addListener(this);		doCalcs(commandRec);	}	/*.................................................................................................................*/	/** passes which object is being disposed (from MesquiteListener interface)*/	public void disposing(Object obj){		//TODO: respond	}	/*.................................................................................................................*/	/** passes which object is being disposed (from MesquiteListener interface)*/	public boolean okToDispose(Object obj, int queryUser){		return true;  //TODO: respond	}	public void changed(Object caller, Object obj, Notification notification, CommandRecord commandRec){		if (Notification.appearsCosmetic(notification))			return;		doCalcs(commandRec);		parametersChanged(notification, commandRec);	}	/*.................................................................................................................*/	public Class getHireSubchoice(){		return NumberForTreeBlock.class;	}	/*.................................................................................................................*/	public Snapshot getSnapshot(MesquiteFile file) { 		Snapshot temp = new Snapshot();		temp.addLine("setValueTask ", numberTask); 		return temp;	}	/*.................................................................................................................*/	public Object doCommand(String commandName, String arguments, CommandRecord commandRec, CommandChecker checker) {		if (checker.compare(this.getClass(), "Sets module that calculates a number for a block of trees", "[name of module]", commandName, "setValueTask")) {			NumberForTreeBlock temp= (NumberForTreeBlock)hireNamedEmployee(commandRec, NumberForTreeBlock.class, arguments);			if (temp!=null) {				numberTask = temp;				return temp;			}		}		else			return super.doCommand(commandName, arguments, commandRec, checker);		return null;	}	public String getTitle() {		if (numberTask==null)			return "";		return numberTask.getVeryShortName();	}	/*.................................................................................................................*/	public void employeeParametersChanged(MesquiteModule employee, MesquiteModule source, Notification notification, CommandRecord commandRec) {		doCalcs(commandRec);		parametersChanged(notification, commandRec);	}	/*.................................................................................................................*/	NumberArray na = new NumberArray(0);	StringArray explArray = new StringArray(0);	/*.................................................................................................................*/	public void doCalcs(CommandRecord commandRec){		if (numberTask==null || treeBlocks == null)			return;		outputInvalid(commandRec);		int numBlocks = treeBlocks.size();		explArray.resetSize(numBlocks);		MesquiteString expl = new MesquiteString();		na.resetSize(numBlocks);		MesquiteNumber mn = new MesquiteNumber();		for (int ic=0; ic<numBlocks; ic++) {			TreeVector trees = (TreeVector)treeBlocks.elementAt(ic);			mn.setToUnassigned();			numberTask.calculateNumber(trees, mn, expl, commandRec);			na.setValue(ic, mn);			explArray.setValue(ic, expl.getValue());		}	}	public String getExplanationForRow(int ic){		if (explArray == null || explArray.getSize() <= ic)			return null;		return explArray.getValue(ic);	}	public String getStringForRow(int ic){		if (na==null)			return "";		return na.toString(ic);	}	public String getWidestString(){		if (numberTask==null)			return "888888";		return numberTask.getVeryShortName()+"   ";	}	/*.................................................................................................................*/	/** returns whether this module is requesting to appear as a primary choice */	public boolean requestPrimaryChoice(){		return true;  	}	/*.................................................................................................................*/	public boolean isPrerelease(){		return false;  	}	public void endJob() {		if (treeBlocks !=null)			treeBlocks.removeListener(this);		super.endJob();	}}