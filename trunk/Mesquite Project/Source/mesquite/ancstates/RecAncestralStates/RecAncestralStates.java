/* Mesquite source code.  Copyright 1997-2006 W. Maddison and D. Maddison.Version 1.11, June 2006.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html) */package mesquite.ancstates.RecAncestralStates;import java.util.*;import java.awt.*;import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.lib.duties.*;/* ======================================================================== */public class RecAncestralStates extends CharHistorySource {	public void getEmployeeNeeds(){  //This gets called on startup to harvest information; override this and inside, call registerEmployeeNeed		EmployeeNeed e = registerEmployeeNeed(CharSourceCoordObed.class, getName() + " needs a source of characters whose ancestral states are to be reconstructed.",		"The source of characters for reconstruction can be chosen initially or later in the Source of Characters or Character Source submenus");		e.setPriority(1);		EmployeeNeed e2 = registerEmployeeNeed(CharMapper.class, getName() + " needs a method by which to reconstruct ancestral states.",		"The reconstruction method can be chosen initially or in the Reconstrction Method submenu");	}	CharMapper assignTask;	public CharSourceCoordObed characterSourceTask;	private CharacterDistribution currentObservedStates=null;	private Taxa currentTaxa=null;	int currentChar=0;	int lastCharRetrieved = -1;	long oldTreeVersion = 0;	long oldTreeID = 0;	private String currentCharacterName;	MesquiteString assignTaskName;	MesquiteCommand atC;	private CharacterHistory recon;	Object hiringCondition;	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, boolean hiredByName) {		hiringCondition = condition;		if (condition!=null) 			characterSourceTask = (CharSourceCoordObed)hireCompatibleEmployee(CharSourceCoordObed.class, condition, "Source of characters to reconstruct");		else			characterSourceTask = (CharSourceCoordObed)hireEmployee(CharSourceCoordObed.class, "Source of characters to reconstruct");		if (characterSourceTask == null) {			return sorry(getName() + " couldn't start because no source of characters obtained.");		}		assignTask = (CharMapper)hireEmployee(CharMapper.class, "Ancestral state reconstruction method");		if (assignTask == null) {			return sorry(getName() + " couldn't start because no reconstruction module obtained.");		}		assignTask.setOneCharacterAtATime(true);		atC = makeCommand("setMethod",  this);		assignTask.setHiringCommand(atC);		assignTaskName = new MesquiteString(assignTask.getName());		if (numModulesAvailable(CharMapper.class)>1){			MesquiteSubmenuSpec mss = addSubmenu(null, "Reconstruction Method", atC, CharMapper.class);			mss.setSelected(assignTaskName);		}		return true;	}	/*.................................................................................................................*/	/** Generated by an employee who quit.  The MesquiteModule should act accordingly. */	public void employeeQuit(MesquiteModule employee) {		if (employee == characterSourceTask || employee == assignTask)  // character source quit and none rehired automatically			iQuit();	}	public boolean allowsStateWeightChoice(){		if (assignTask == null)			return false;		else			return assignTask.allowsStateWeightChoice();	}	/*.................................................................................................................*/	public Snapshot getSnapshot(MesquiteFile file) {		Snapshot temp = new Snapshot();		temp.addLine("getCharacterSource ",characterSourceTask);		temp.addLine("setMethod ",assignTask);		return temp;	}	MesquiteInteger pos = new MesquiteInteger();	/*.................................................................................................................*/	public Object doCommand(String commandName, String arguments, CommandChecker checker) {		if (checker.compare(this.getClass(), "Sets module used to reconstruct ancestral states", "[name of module]", commandName, "setMethod")) {			CharMapper temp=  (CharMapper)replaceEmployee(CharMapper.class, arguments, "Reconstruction method", assignTask);			if (temp!=null) {				assignTask= temp;				assignTask.setHiringCommand(atC);				assignTaskName.setValue(assignTask.getName());				assignTask.setOneCharacterAtATime(true);				parametersChanged();			}			return assignTask;		}		else if (checker.compare(this.getClass(), "Sets module supplying characters", "[name of module]", commandName, "setCharacterSource")) { //temporary while data files still exist that don't use coordinators			return characterSourceTask.doCommand(commandName, arguments, checker);		}		else if (checker.compare(this.getClass(), "Returns module supplying characters", null, commandName, "getCharacterSource")) {			return characterSourceTask;		}		else			return  super.doCommand(commandName, arguments, checker);	}	Tree tree = null;	public  void prepareHistory(Tree tree, int ic){		this.tree = tree;		if (tree == null) {			return;		}		if (tree.getTaxa() != currentTaxa|| (characterSourceTask.usesTree() && (tree.getID() != oldTreeID || tree.getVersionNumber() != oldTreeVersion)) || ic != lastCharRetrieved || currentObservedStates == null ) {			int maxnum = characterSourceTask.getNumberOfCharacters(tree);			if (ic>= maxnum)				ic = maxnum-1;			currentObservedStates = characterSourceTask.getCharacter(tree, ic);			currentTaxa = tree.getTaxa();			oldTreeVersion = tree.getVersionNumber();			oldTreeID = tree.getID();			lastCharRetrieved = ic;		}		if (currentObservedStates == null)			return;		assignTask.setObservedStates(tree, currentObservedStates);	}	/*.................................................................................................................*/	public CharacterHistory getMapping(int im, CharacterHistory history, MesquiteString resultString) {		if (tree == null){			if (!MesquiteThread.isScripting())				MesquiteMessage.warnUser("Error: tree not obtained for reconstruct ancestral states");			if (resultString!=null)				resultString.setValue("No character history could be reconstructed because no tree was supplied for the calculations");			return null;		}		if (currentObservedStates == null) {			if (!MesquiteThread.isScripting())				MesquiteMessage.warnUser("Error: character distribution not obtained for reconstruct ancestral states");			if (resultString!=null)				resultString.setValue("No character history could be reconstructed because no character distribution was supplied for the calculations");			return null;		}		//preparing to receive mapping		history =currentObservedStates.adjustHistorySize(tree, history);		assignTask.getMapping(im, history, resultString);		if (resultString!=null)			resultString.prepend("Character: " + currentObservedStates.getName() + "\n");		if (history == null) { 			if (resultString!=null)				resultString.setValue("Ancestral states not successfully reconstructed at nodes.");			MesquiteMessage.warnProgrammer("Error: states not successfully reconstructed at nodes (Reconstruct Ancestral States)");			return null;		}		else  {			history.setObservedStates(currentObservedStates);			return history;		}	}	/** returns the name of history ic*/	public String getHistoryName(Taxa taxa, int ic){		return characterSourceTask.getCharacterName(taxa, ic);	}	/** returns the name of history ic*/	public String getHistoryName(Tree tree, int ic){		return characterSourceTask.getCharacterName(tree, ic);	}	public int getNumberOfHistories(Tree tree){		if (characterSourceTask == null || tree == null)			return 0;		return characterSourceTask.getNumberOfCharacters(tree.getTaxa());	}	public int getNumberOfHistories(Taxa taxa){		if (characterSourceTask == null)			return 0;		return characterSourceTask.getNumberOfCharacters(taxa);	}	public int getNumberOfMappings(Tree tree,  int ic){		if (assignTask == null)			return 0;		return assignTask.getNumberOfMappings();	}	public int getNumberOfMappings(Taxa taxa,  int ic){		if (assignTask == null)			return 0;		return assignTask.getNumberOfMappings();	}	/** returns the name of history ic and mapping im*/	public String getMappingName(Taxa taxa, int ic, int im){		return getHistoryName(taxa, ic);	}	/** returns the name of history ic and mapping im*/	public String getMappingName(Tree tree, int ic, int im){		return getHistoryName(tree, ic);	}	/** returns the name of histories for menu items, e.g. if each history represents a character, return "Character"*/	public  String getHistoryTypeName(){		return "Character";	}	/*.................................................................................................................*/	public void employeeParametersChanged(MesquiteModule employee, MesquiteModule source, Notification notification) {		if (employee==characterSourceTask || employee == assignTask) {			currentObservedStates = null;  // to force retrieval of new observedStates			parametersChanged(notification);		}	}	/*.................................................................................................................*/	public String getName() {		return "Reconstruct Ancestral States";	}	/*.................................................................................................................*/	public String getNameAndParameters() {		if (assignTask == null)			return getName();		return "Reconstructed Ancestral States (" + assignTask.getName() + ")";	}	public boolean showCitation(){		return true;	}	/*.................................................................................................................*/	public boolean isPrerelease(){		return false;	}	/*.................................................................................................................*/	public boolean isSubstantive(){		return true;	}	/*.................................................................................................................*/	/** returns whether this module is requesting to appear as a primary choice */	public boolean requestPrimaryChoice(){		return true;  	}	/*.................................................................................................................*/	public String getExplanation() {		return "Reconstructs ancestral states on the nodes of a tree.";	}}