/* Mesquite module ~~*/package mesquite.speciation.StoredCharHistories;/*~~  */import java.applet.*;import java.util.*;import java.awt.*;import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.lib.duties.*;import mesquite.categ.lib.*;/** ======================================================================== */public class StoredCharHistories extends CharHistorySource {	private Taxa currentTaxa=null;	int currentChar=0;	int lastCharRetrieved = -1;	long oldTreeVersion = 0;	long oldTreeID = 0;		/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) {		if (getProject().getNumberOfFileElements(StoredCharacterHistory.class)==0){			return sorry(commandRec, "There are no stored character histories.");		}		return true; 	}	public boolean allowsStateWeightChoice(){		return false;	}	/*.................................................................................................................*/  	 public Snapshot getSnapshot(MesquiteFile file) {   	 	Snapshot temp = new Snapshot();  	 	//temp.addLine("setMethod ",assignTask);  	 	//temp.addLine("setCharacterSource ",characterSourceTask);  	 	return temp;  	 }	MesquiteInteger pos = new MesquiteInteger();	/*.................................................................................................................*/    	 public Object doCommand(String commandName, String arguments, CommandRecord commandRec, CommandChecker checker) {    	 	if (checker.compare(this.getClass(), "Sets module used to reconstruct ancestral states", "[name of module]", commandName, "setMethod")) {    	 	}    	 	else if (checker.compare(this.getClass(), "Sets module supplying characters", "[name of module]", commandName, "setCharacterSource")) {    	 	}    	 	else    	 		return super.doCommand(commandName, arguments, commandRec, checker);		return null;   	 }   	CharacterHistory getBlankHistory(Tree tree){   		return new CategoricalHistory(tree.getTaxa(), tree.getNumNodeSpaces());   	}   	Tree tree = null;	public  void prepareHistory(Tree tree, int ic, CommandRecord commandRec){		this.tree = null;		if (tree == null) {			return;		}  		if (ic<0)   			return;   		this.currentChar =ic;	}	/*.................................................................................................................*/	public CharacterHistory getMapping(int im, CharacterHistory history, MesquiteString resultString, CommandRecord commandRec) {   		if (currentChar<0)   			return null; 		if (tree == null) {			return null;		} 		resultString.setValue("");//check to see if tree is same tree!!!!		if (currentChar >= getProject().getNumberOfFileElements(StoredCharacterHistory.class)) {			if (resultString!=null)				resultString.setValue("Sorry, stored history was not found");			return getBlankHistory(tree);		}		StoredCharacterHistory sch = ((StoredCharacterHistory)getProject().getFileElement(StoredCharacterHistory.class, currentChar));		if (sch==null) {			if (resultString!=null)				resultString.setValue("Sorry, stored history was not found");			return getBlankHistory(tree);		}		if (tree == sch.getTree()) {			if (resultString!=null)				resultString.setValue("Stored History: \n" + sch.getName());			return sch.getHistory();		}		if (tree instanceof MesquiteTree && sch.getTree() instanceof MesquiteTree){			if (!((MesquiteTree)tree).equalsCoreArrays((MesquiteTree)sch.getTree())) {				if (resultString!=null)					resultString.setValue("Sorry, stored history doesn't apply to same tree (history: " + sch.getName() + ")");				return getBlankHistory(tree);			}		}		else			if (tree!=sch.getTree()) {				MesquiteMessage.warnUser("Sorry, trees don't match in character history chosen");				if (resultString!=null)					resultString.setValue("Sorry, stored history doesn't apply to same tree (history: " + sch.getName() + ")");				return getBlankHistory(tree);			}				if (resultString!=null)			resultString.setValue("Stored History: \n" + sch.getName());		return sch.getHistory();	}   	/** returns the name of history ic*/   	public String getHistoryName(Taxa taxa, int ic, CommandRecord commandRec){   		return "Stored character history " + ic;   	}   	/** returns the name of history ic*/   	public String getHistoryName(Tree tree, int ic, CommandRecord commandRec){   		return "Stored character history " + ic;   	}	public int getNumberOfHistories(Tree tree,  CommandRecord commandRec){		if (tree == null)			return 0;		return getProject().getNumberOfFileElements(StoredCharacterHistory.class);	}	public int getNumberOfHistories(Taxa taxa,  CommandRecord commandRec){		return getProject().getNumberOfFileElements(StoredCharacterHistory.class);	}	public int getNumberOfMappings(Tree tree,  int ic, CommandRecord commandRec){		return 1;	}	public int getNumberOfMappings(Taxa taxa,  int ic, CommandRecord commandRec){		return 1;	}  	/** returns the name of history ic and mapping im*/   	public String getMappingName(Taxa taxa, int ic, int im, CommandRecord commandRec){   		return getHistoryName(taxa, ic, commandRec);   	}   	/** returns the name of history ic and mapping im*/   	public String getMappingName(Tree tree, int ic, int im, CommandRecord commandRec){   		return getHistoryName(tree, ic, commandRec);   	}  	/** returns the name of histories for menu items, e.g. if each history represents a character, return "Character"*/   	public  String getHistoryTypeName(){   		return "Character";   	}	/*.................................................................................................................*/    	 public String getName() {		return "Stored character histories";   	 }	/*.................................................................................................................*/  	 public String getExplanation() { 		return "Supplies character histories from data files (as opposed to simulated character histories, for example)." ;   	 }}