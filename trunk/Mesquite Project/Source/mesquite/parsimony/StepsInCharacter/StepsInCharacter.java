/* Mesquite source code.  Copyright 1997-2006 W. Maddison and D. Maddison. Version 1.11, June 2006.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.parsimony.StepsInCharacter;/*~~  */import java.util.*;import java.awt.*;import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.lib.duties.*;import mesquite.parsimony.lib.*;/* ======================================================================== */public class StepsInCharacter extends NumberForTree implements Incrementable {	MesquiteNumber steps;	CharSourceCoordObed characterSourceTask;	CharacterSteps charStepsTask;	Taxa oldTaxa = null;	CharacterDistribution charStates;	int currentChar = 0;	long oldTreeID;	long oldTreeVersion;	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) {    	 	steps=new MesquiteNumber(); 		characterSourceTask = (CharSourceCoordObed)hireEmployee(commandRec, CharSourceCoordObed.class, "Source of characters for parsimony steps"); 		if (characterSourceTask == null) 			return sorry(commandRec, getName() + " couldn't start because no source of characters was obtained."); 		charStepsTask = (CharacterSteps)hireEmployee(commandRec, CharacterSteps.class, null); 		if (charStepsTask == null) 			return sorry(commandRec, getName() + " couldn't start because no steps counting module was obtained.");		addMenuItem( "Next Character", makeCommand("nextCharacter",  this));		addMenuItem( "Previous Character", makeCommand("previousCharacter",  this));		addMenuItem( "Choose Character...", makeCommand("chooseCharacter",  this));  		return true;  	 }  	 	/*.................................................................................................................*/	/** returns whether this module is requesting to appear as a primary choice */   	public boolean requestPrimaryChoice(){   		return true;     	}  	 public void employeeQuit(MesquiteModule m){  	 	if (m==charStepsTask || m == characterSourceTask)  	 		iQuit();  	 }   	/** Called to provoke any necessary initialization.  This helps prevent the module's intialization queries to the user from   	happening at inopportune times (e.g., while a long chart calculation is in mid-progress)*/   	public void initialize(Tree tree, CommandRecord commandRec){   		if (tree == null)   			return;   		oldTaxa = tree.getTaxa();   		characterSourceTask.initialize(tree.getTaxa(), commandRec);   	}	/*.................................................................................................................*/	public void calculateNumber(Tree tree, MesquiteNumber result, MesquiteString resultString, CommandRecord commandRec) {    	 	if (result==null)    	 		return;		result.setToUnassigned();		if (resultString!=null)			resultString.setValue("");    	 	steps.setValue((int)0);    	 	Taxa taxa = tree.getTaxa(); 		if (taxa != oldTaxa && oldTaxa!=null || (characterSourceTask.usesTree() && (tree.getID() != oldTreeID || tree.getVersionNumber() != oldTreeVersion))) { 			currentChar = 0; 			charStates = characterSourceTask.getCharacter(tree, currentChar, commandRec);			oldTreeID = tree.getID();			oldTreeVersion = tree.getVersionNumber(); 		} 		else if (charStates == null) {			if (currentChar<0 || currentChar>=characterSourceTask.getNumberOfCharacters(tree, commandRec))				currentChar = 0; 			charStates = characterSourceTask.getCharacter(tree, currentChar, commandRec);		}		result.setToUnassigned();		charStepsTask.calculateNumber(tree, charStates, result, resultString, commandRec);		steps.setValue(result);		oldTaxa = taxa;		result.setValue(steps);		if (resultString!=null)			resultString.append(" (char. " + CharacterStates.toExternal(currentChar) + ")");	}	/*.................................................................................................................*/	public boolean biggerIsBetter() {		return false;	}	/*.................................................................................................................*/   	 public void employeeParametersChanged(MesquiteModule employee, MesquiteModule source, Notification notification, CommandRecord commandRec) {   	 	if (employee==characterSourceTask) {   			currentChar = 0;   			charStates = null;			parametersChanged(notification, commandRec);   	 	}   	 	else   	 		super.employeeParametersChanged(employee, source, notification, commandRec);   	 }	/*.................................................................................................................*/   	 public void setCurrent(long i, CommandRecord commandRec){ 		if (characterSourceTask==null || oldTaxa==null) 			return; 		if ((i>=0) && (i<=characterSourceTask.getNumberOfCharacters(oldTaxa, commandRec)-1)) { 			currentChar = (int)i;			charStates=null;			//parametersChanged(null, commandRec);		}   	 } 	public String getItemTypeName(){ 		return "Character"; 	}	/*.................................................................................................................*/ 	public long toInternal(long i){ 		return(CharacterStates.toInternal((int)i)); 	}	/*.................................................................................................................*/ 	public long toExternal(long i){ 		return(CharacterStates.toExternal((int)i)); 	}	/*.................................................................................................................*/   	 public long getCurrent(CommandRecord commandRec){   	 	return currentChar;   	 }	/*.................................................................................................................*/ 	public long getMin(CommandRecord commandRec){ 		return 0; 	}	/*.................................................................................................................*/ 	public long getMax(CommandRecord commandRec){ 		if (characterSourceTask==null || oldTaxa==null) 			return 0; 		return characterSourceTask.getNumberOfCharacters(oldTaxa, commandRec)-1; 	}	/*.................................................................................................................*/  	 public Snapshot getSnapshot(MesquiteFile file) {   	 	Snapshot temp = new Snapshot();  	 	temp.addLine("getCharacterSource ",characterSourceTask);  	 	temp.addLine("setCharacter " + CharacterStates.toExternal(currentChar));  	 	return temp;  	 }	/*.................................................................................................................*/    	 public Object doCommand(String commandName, String arguments, CommandRecord commandRec, CommandChecker checker) {    	 	if (checker.compare(this.getClass(), "Sets module supplying characters", "[name of module]", commandName, "setCharacterSource")) {//temporary, for data files using old system without coordinators 			return characterSourceTask.doCommand(commandName, arguments, commandRec, checker);    	 	}    	 	else if (checker.compare(this.getClass(), "Returns module supplying characters", null, commandName, "getCharacterSource")) { 			return characterSourceTask;    	 	}    	 	else if (checker.compare(this.getClass(), "Goes to next character", null, commandName, "nextCharacter")) {    	 		if (currentChar>=characterSourceTask.getNumberOfCharacters(oldTaxa, commandRec)-1)    	 			currentChar=0;    	 		else    	 			currentChar++;			charStates = null; 			parametersChanged(null, commandRec);     	 	}    	 	else if (checker.compare(this.getClass(), "Goes to previous character", null, commandName, "previousCharacter")) {    	 		if (currentChar<=0)    	 			currentChar=characterSourceTask.getNumberOfCharacters(oldTaxa, commandRec)-1;    	 		else    	 			currentChar--;   			charStates = null; 			parametersChanged(null, commandRec);     	 	}    	 	else if (checker.compare(this.getClass(), "Queries the user about what character to use", null, commandName, "chooseCharacter")) {    	 		int ic=characterSourceTask.queryUserChoose(oldTaxa, " to count steps ", commandRec);    	 		if (MesquiteInteger.isCombinable(ic)) {	   			currentChar = ic;	   			charStates = null;	 			parametersChanged(null, commandRec);  			}    	 	}    	 	else if (checker.compare(this.getClass(), "Sets the character to use", "[character number]", commandName, "setCharacter")) {    	 		int icNum = MesquiteInteger.fromFirstToken(arguments, stringPos);    	 		if (!MesquiteInteger.isCombinable(icNum))    	 			return null;    	 		int ic = CharacterStates.toInternal(icNum);    	 		if ((ic>=0) && characterSourceTask.getNumberOfCharacters(oldTaxa, commandRec)==0) {    	 			currentChar = ic;	   			charStates = null;    	 		}    	 		else if ((ic>=0) && (ic<=characterSourceTask.getNumberOfCharacters(oldTaxa, commandRec)-1)) {    	 			currentChar = ic;	   			charStates = null;	 			parametersChanged(null, commandRec);  			}    	 	}    	 	else    	 		return super.doCommand(commandName, arguments, commandRec, checker);		return null;   	 }   	/*.................................................................................................................*/    	 public String getParameters() {    	 	if (characterSourceTask == null)    	 		return null;		return "Current character: " + CharacterStates.toExternal(currentChar) + " from: " + characterSourceTask.getNameAndParameters() ;   	 }	/*.................................................................................................................*/    	 public String getName() {		return "Steps in Character";   	 }	/*.................................................................................................................*/    	 public String getVeryShortName() {		return "Steps";   	 }    	 public boolean showCitation(){   	 	return true;   	 }  	 public boolean isPrerelease(){   	 	return false;   	 }	/*.................................................................................................................*/  	 public String getExplanation() {		return "Calculates the parsimony steps in a character for a given tree.";   	 }}