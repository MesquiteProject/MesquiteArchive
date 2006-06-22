/* Mesquite source code.  Copyright 1997-2006 W. Maddison and D. Maddison.Version 1.11, June 2006.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.cont.MapContinuous;/*~~  */import java.util.*;import java.awt.*;import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.lib.duties.*;import mesquite.cont.lib.*;/* ======================================================================== */public class MapContinuous extends NumbersForNodesIncr implements Incrementable {	CharStatesForNodes assignTask;	CharacterDistribution observedStates;	ContinuousHistory reconstructedStates, recon;	CharSourceCoordObed characterSourceTask;	NumberArray numbers;	Taxa currentTaxa;	Tree tree;	int currentChar=0;	long oldTreeVersion = 0;	long oldTreeID = 0;	boolean characterSet = false;	int lastCharRetrieved = -1;		//choice of what item to show	int currentItem=0;	MesquiteMenuItemSpec itemItem;	MesquiteCommand itemChoiceCommand;	String itemName=null;	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) { 		characterSourceTask = (CharSourceCoordObed)hireCompatibleEmployee(commandRec, CharSourceCoordObed.class, ContinuousState.class, "Source of characters (for Map Continuous)"); 		if (characterSourceTask == null) { 			return sorry(commandRec, getName() + " couldn't start because no source of characters obtained."); 		}				//todo: allow subsequentChoice 		assignTask = (CharStatesForNodes)hireNamedEmployee(commandRec, CharStatesForNodes.class, StringUtil.tokenize("Parsimony Ancestral States")); 		if (assignTask == null)			assignTask = (CharStatesForNodes)hireEmployee(commandRec, CharStatesForNodes.class, "Reconstruction method"); 		if (assignTask == null) 			return sorry(commandRec, getName() + " couldn't start because no reconstructing module obtained.");		/**/		if (!(NumbersForNodesIncr.class.isAssignableFrom(getHiredAs()))){			addMenuItem( "Next Character", makeCommand("nextCharacter",  this));			addMenuItem( "Previous Character", makeCommand("previousCharacter",  this));			addMenuItem( "Choose Character", makeCommand("chooseCharacter",  this));		}		/**/		itemChoiceCommand = MesquiteModule.makeCommand("setItem",  this);		itemItem = addMenuItem("Item to Map...", itemChoiceCommand);		return true;   	} 		/*.................................................................................................................*/	/** Generated by an employee who quit.  The MesquiteModule should act accordingly. */ 	public void employeeQuit(MesquiteModule employee) { 		if (employee == characterSourceTask || employee == assignTask)  // character source quit and none rehired automatically 			iQuit();	}	/*.................................................................................................................*/	private String itemString(){ 		if (itemName!=null)			return " (item " + itemName + ")"; 		else			return "";  	} 	/** returns current parameters, for logging etc..*/ 	public String getParameters() { 			return assignTask.getName() + " of " + characterSourceTask.getNameAndParameters() + itemString();   	 }	/*.................................................................................................................*/   	 public void setCurrent(long i, CommandRecord commandRec){ 		if (characterSourceTask==null || currentTaxa==null){ 			currentChar = (int)i;			characterSet=true; 		} 		else if ((i>=0) && (i<=characterSourceTask.getNumberOfCharacters(currentTaxa, commandRec)-1)) { 			currentChar = (int)i;			characterSet=true;		}   	 } 	public String getItemTypeName(){ 		return "Character"; 	}	/*.................................................................................................................*/ 	public long toInternal(long i){ 		return(CharacterStates.toInternal((int)i)); 	}	/*.................................................................................................................*/ 	public long toExternal(long i){ 		return(CharacterStates.toExternal((int)i)); 	}	/*.................................................................................................................*/   	 public long getCurrent(CommandRecord commandRec){   	 	return currentChar;   	 }	/*.................................................................................................................*/ 	public long getMin(CommandRecord commandRec){ 		return 0; 	}	/*.................................................................................................................*/ 	public long getMax(CommandRecord commandRec){ 		if (characterSourceTask==null || currentTaxa==null) 			return 0; 		return characterSourceTask.getNumberOfCharacters(currentTaxa, commandRec)-1; 	}	/*.................................................................................................................*/  	 public Snapshot getSnapshot(MesquiteFile file) {  	 //TODO: allow change in assignTask, etc.   	 	Snapshot temp = new Snapshot();  	 	temp.addLine( "getCharacterSource " , characterSourceTask);  	 	temp.addLine("setCharacter " + CharacterStates.toExternal(currentChar));  	 	temp.addLine("setItem " + (currentItem));  	 	return temp;  	 }	/*.................................................................................................................*/    	 public Object doCommand(String commandName, String arguments, CommandRecord commandRec, CommandChecker checker) {    	 	if (checker.compare(this.getClass(), "Sets module supplying characters", "[name of module]", commandName, "setCharacterSource")) {//temporary, for data files using old system without coordinators 			return characterSourceTask.doCommand(commandName, arguments, commandRec, checker);    	 	}    	 	else if (checker.compare(this.getClass(), "Returns module supplying characters", null, commandName, "getCharacterSource")) { 			return characterSourceTask;    	 	}    	 	else if (checker.compare(this.getClass(), "Goes to next character", null, commandName, "nextCharacter")) {    	 		if (currentChar>=characterSourceTask.getNumberOfCharacters(currentTaxa, commandRec)-1)    	 			currentChar=0;    	 		else    	 			currentChar++;				parametersChanged(null, commandRec);    	 	}    	 	else if (checker.compare(this.getClass(), "Goes to next character", null, commandName, "previousCharacter")) {    	 		if (currentChar<=0)    	 			currentChar=characterSourceTask.getNumberOfCharacters(currentTaxa, commandRec)-1;    	 		else    	 			currentChar--;				parametersChanged(null, commandRec);    	 	}    	 	else if (checker.compare(this.getClass(), "Queries the user which character to use", null, commandName, "chooseCharacter")) {    	 		int ic=characterSourceTask.queryUserChoose(currentTaxa, " to map " + whatIsMyPurpose(), commandRec);    	 		if (MesquiteInteger.isCombinable(ic)) {	   			currentChar = ic;	 			characterSet=true;	 			parametersChanged(null, commandRec); //? 			}    	 	}    	 	else if (checker.compare(this.getClass(), "Sets the character to use", "[number of character]", commandName, "setCharacter")) {    	 		int ic = CharacterStates.toInternal(MesquiteInteger.fromString(arguments));   			if (currentTaxa==null) {    	 			currentChar = ic;	 			characterSet=true;   			}    	 		if ((ic>=0) && (ic<=characterSourceTask.getNumberOfCharacters(currentTaxa, commandRec)-1)) {    	 			currentChar = ic;	 			characterSet=true;				parametersChanged(null, commandRec); 			}    	 	}    	 	else if (checker.compare(this.getClass(), "Sets the item to use (in a multi-item continuous data matrix)", "[item number]", commandName, "setItem")) {    	 		int ic = MesquiteInteger.fromString(arguments);    	 		if (!MesquiteInteger.isCombinable(ic) && reconstructedStates!=null){				ic = reconstructedStates.userQueryItem("Select item to map", this);    	 		}   			if (!MesquiteInteger.isCombinable(ic))   				return null;   			if (currentTaxa==null) {    	 			currentItem = ic;   			} 			else if (reconstructedStates !=null && reconstructedStates instanceof ContinuousHistory) {	   	 		if ((ic>=0) && (ic<=reconstructedStates.getNumItems()-1)) {	    	 			currentItem = ic;					parametersChanged(null, commandRec);	 			} 			}    	 	}    	 	else    	 		return super.doCommand(commandName, arguments, commandRec, checker);		return null;   	 }	/*.................................................................................................................*/ 	public void employeeParametersChanged(MesquiteModule employee, MesquiteModule source, Notification notification, CommandRecord commandRec) { 		observedStates = null;//to force recalculation 		super.employeeParametersChanged(module, source, notification, commandRec);   	 }   	/** Called to provoke any necessary initialization.  This helps prevent the module's intialization queries to the user from   	happening at inopportune times (e.g., while a long chart calculation is in mid-progress)*/   	public void initialize(Tree tree, CommandRecord commandRec){   		currentTaxa = tree.getTaxa();   		characterSourceTask.initialize(currentTaxa, commandRec);   	}   int times = 0;	/*.................................................................................................................*/	public  void calculateNumbers(Tree tree, NumberArray result, MesquiteString resultString, CommandRecord commandRec) {		if (tree==null)			return;		this.tree = tree;		itemName=null;		Taxa taxa = tree.getTaxa();		result.deassignArray();		if (taxa != currentTaxa || (characterSourceTask.usesTree() && (tree.getID() != oldTreeID || tree.getVersionNumber() != oldTreeVersion)) || currentChar != lastCharRetrieved || observedStates == null ) { 			int maxnum = characterSourceTask.getNumberOfCharacters(tree, commandRec);			if (currentChar>= maxnum)				currentChar = maxnum-1;			observedStates = characterSourceTask.getCharacter(tree, currentChar, commandRec);			currentTaxa = taxa;			oldTreeVersion = tree.getVersionNumber();			oldTreeID = tree.getID();			lastCharRetrieved = currentChar;		}		if (observedStates ==null)			return;		recon = (ContinuousHistory)observedStates.adjustHistorySize(tree, recon);		assignTask.calculateStates(tree, observedStates, recon, null, commandRec);		if (recon!=null)			reconstructedStates=(ContinuousHistory)recon.clone(reconstructedStates);				if (reconstructedStates != null) {			int numItems = reconstructedStates.getNumItems();			if (numItems>1){				if (!itemItem.isEnabled()){					itemItem.setEnabled(true);					MesquiteTrunk.mesquiteTrunk.resetMenuItemEnabling();				}				if (currentItem>= numItems)					currentItem=0;				itemName = reconstructedStates.getItemName(currentItem);							}			else {				currentItem = 0;				if (itemItem.isEnabled()){					itemItem.setEnabled(false);					MesquiteTrunk.mesquiteTrunk.resetMenuItemEnabling();				}			} 					}		if (resultString !=null)			resultString.setValue(itemString());					if (reconstructedStates !=null && reconstructedStates instanceof ContinuousHistory) {			DoubleArray dub = reconstructedStates.getItem(currentItem);			result.setValues(dub.getMatrix());		}	}	/*.................................................................................................................*/    	 public String getName() {		return "Map Continuous";   	 }	/*.................................................................................................................*/	public boolean showCitation(){		return true;	}	/*.................................................................................................................*/	public boolean isPrerelease(){		return false;	}	/*.................................................................................................................*/	/** returns whether this module is requesting to appear as a primary choice */   	public boolean requestPrimaryChoice(){   		return true;     	}	/*.................................................................................................................*/ 	/** returns an explanation of what the module does.*/ 	public String getExplanation() { 		return "Supplies a reconstruction of continuous-valued ancestral states on a tree."  		+ " Differs from the basic reconstruction methods in that it supplies simply numbers for nodes, instead of a CharacterHistory";   	 }   	 }