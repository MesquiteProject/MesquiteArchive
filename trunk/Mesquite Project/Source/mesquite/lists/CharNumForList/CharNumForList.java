/* Mesquite source code.  Copyright 1997-2006 W. Maddison and D. Maddison. Version 1.11, June 2006.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html) */package mesquite.lists.CharNumForList;/*~~  */import mesquite.lists.lib.*;import java.util.*;import java.awt.*;import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.lib.duties.*;import mesquite.lib.table.*;/* ======================================================================== */public class CharNumForList extends CharListAssistant implements MesquiteListener {	/*.................................................................................................................*/	public String getName() {		return "Number for Character (in List of Characters window)";	}	public String getNameForMenuItem() {		return "Number for Character";	}	public String getExplanation() {		return "Supplies numbers for characters for a character list window." ;	}	public void getEmployeeNeeds(){  //This gets called on startup to harvest information; override this and inside, call registerEmployeeNeed		EmployeeNeed e = registerEmployeeNeed(NumberForCharacter.class, getName() + " needs a method to calculate a value for each of the characters.",		"You can select a value to show in the Number For Characters submenu of the Columns menu of the List of Characters Window. ");	}	/*.................................................................................................................*/	CharacterData data=null;	NumberForCharacter numberTask;	Object dataCondition;	MesquiteBoolean shadeCells;	boolean suppressed = false;	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) {		if (arguments!=null) {			numberTask = (NumberForCharacter)hireNamedEmployee(commandRec, NumberForCharacter.class, arguments);			if (numberTask==null) {				return sorry(commandRec, "Number for character (for character list) could not start because the requested calculator module was not obtained");			}		}		if (!commandRec.scripting() && numberTask == null) {			numberTask = (NumberForCharacter)hireEmployee(commandRec, NumberForCharacter.class, "Value to calculate for characters (for Character List)");			if (numberTask==null) {				return sorry(commandRec, "Number for character (for character list) could not start because no calculator module was obtained");			}		}		shadeCells = new MesquiteBoolean(false);		addCheckMenuItem(null, "Color Cells", makeCommand("toggleShadeCells",  this), shadeCells);		return true;	}	/*.................................................................................................................*/	/** Returns whether or not it's appropriate for an employer to hire more than one instance of this module.   	If false then is hired only once; second attempt fails.*/	public boolean canHireMoreThanOnce(){		return true;	}	/*.................................................................................................................*/	public void employeeQuit(MesquiteModule m){		iQuit();	}	/*.................................................................................................................*/	public Class getHireSubchoice(){		return NumberForCharacter.class;	}	/*.................................................................................................................*/	public Snapshot getSnapshot(MesquiteFile file) { 		Snapshot temp = new Snapshot();		temp.addLine("suppress"); 		temp.addLine("setValueTask ", numberTask); 		temp.addLine("toggleShadeCells " + shadeCells.toOffOnString());		temp.addLine("desuppress"); 		return temp;	}	/*.................................................................................................................*/	public Object doCommand(String commandName, String arguments, CommandRecord commandRec, CommandChecker checker) {		if (checker.compare(this.getClass(), "Sets the module to calculate a number for each character", "[name of module]", commandName, "setValueTask")) {			NumberForCharacter temp= (NumberForCharacter)replaceEmployee(commandRec, NumberForCharacter.class, arguments, "Number for a character", numberTask);			if (temp!=null) {				numberTask = temp;				if (!suppressed){					doCalcs(commandRec);					parametersChanged(null, commandRec);				}				return temp;			}		}		else if (checker.compare(this.getClass(), "Sets whether or not to color cells", "[on or off]", commandName, "toggleShadeCells")) {			boolean current = shadeCells.getValue();			shadeCells.toggleValue(parser.getFirstToken(arguments));			if (current!=shadeCells.getValue()) {				outputInvalid(commandRec);				parametersChanged(null, commandRec);			}		}		else if (checker.compare(this.getClass(), "Suppresses calculation", null, commandName, "suppress")) {			suppressed = true;		}		else if (checker.compare(this.getClass(), "Releases suppression of calculation", null, commandName, "desuppress")) {			if (suppressed){				suppressed = false;				outputInvalid(commandRec);				doCalcs(commandRec);				parametersChanged(null, commandRec);			}		}		else			return super.doCommand(commandName, arguments, commandRec, checker);		return null;	}	public String getTitle() {		if (numberTask==null)			return "";		return numberTask.getVeryShortName();	}	/*.................................................................................................................*/	/** passes which object is being disposed (from MesquiteListener interface)*/	public void disposing(Object obj){		if (obj == data)			data=null;	}	/*.................................................................................................................*/	/** passes which object is being disposed (from MesquiteListener interface)*/	public boolean okToDispose(Object obj, int queryUser){		return true;  //TODO: respond	}	/*.................................................................................................................*/	public void changed(Object caller, Object obj, Notification notification, CommandRecord commandRec){		int code = Notification.getCode(notification);		if (code == AssociableWithSpecs.SPECSSET_CHANGED || code == MesquiteListener.VALUE_CHANGED ||code == MesquiteListener.DATA_CHANGED || code == MesquiteListener.PARTS_CHANGED || code == MesquiteListener.PARTS_ADDED || code == MesquiteListener.PARTS_DELETED || code == MesquiteListener.PARTS_MOVED){			if (!suppressed){				outputInvalid(commandRec);				doCalcs(commandRec);				parametersChanged(notification, commandRec);			}		}	}	/*.................................................................................................................*/	public void employeeParametersChanged(MesquiteModule employee, MesquiteModule source, Notification notification, CommandRecord commandRec) {		if (!suppressed){			outputInvalid(commandRec);			doCalcs(commandRec);			parametersChanged(notification, commandRec);		}	}	/*.................................................................................................................*/	public void setTableAndData(MesquiteTable table, CharacterData data, CommandRecord commandRec){		if (data==null)			return;		if (this.data !=null)			this.data.removeListener(this);		this.data = data;		data.addListener(this);		if (!suppressed)			doCalcs(commandRec);		//table would be used if selection needed	}	/*.................................................................................................................*/	NumberArray na = new NumberArray(0);	StringArray explArray = new StringArray(0);	MesquiteNumber min = new MesquiteNumber();	MesquiteNumber max = new MesquiteNumber();	/*.................................................................................................................*/	public void doCalcs(CommandRecord commandRec){		if (numberTask==null || data == null || suppressed)			return;		int numChars = data.getNumChars();		na.resetSize(numChars);		explArray.resetSize(numChars);		MesquiteString expl = new MesquiteString();		na.deassignArrayToInteger();		MesquiteNumber mn = new MesquiteNumber();		for (int ic=0; ic<numChars; ic++) {			commandRec.tick("Number for character in character list; examining character " + ic);			CharacterDistribution states = data.getCharacterDistribution(ic);			mn.setToUnassigned();			numberTask.calculateNumber(states, mn, expl, commandRec);			na.setValue(ic, mn);			explArray.setValue(ic, expl.getValue());		}		na.placeMinimumValue(min);		na.placeMaximumValue(max);	}	/** Gets background color for cell for row ic.  Override it if you want to change the color from the default. */	public Color getBackgroundColorOfCell(int ic, boolean selected){		if (!shadeCells.getValue())			return null;		if (min.isCombinable() && max.isCombinable() && na != null && na.isCombinable(ic)){			return MesquiteColorTable.getGreenScale(na.getDouble(ic), min.getDoubleValue(), max.getDoubleValue(), false);		}		return null;	}	public String getExplanationForRow(int ic){		if (explArray == null || explArray.getSize() <= ic)			return null;		return explArray.getValue(ic);	}	public boolean useString(int ic){		return true;	}	public String getStringForCharacter(int ic){		if (na==null)			return "";		return na.toString(ic); 	}	public String getWidestString(){		if (numberTask==null)			return "888888";		return numberTask.getVeryShortName()+"   ";	}	/*.................................................................................................................*/	/** returns whether this module is requesting to appear as a primary choice */	public boolean requestPrimaryChoice(){		return true;  	}	/*.................................................................................................................*/	public boolean isPrerelease(){		return false;  	}	public void endJob() {		if (data !=null)			data.removeListener(this);		super.endJob();	}}