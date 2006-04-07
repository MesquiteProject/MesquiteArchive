/* Mesquite source code.  Copyright 1997-2005 W. Maddison and D. Maddison. Version 1.06, August 2005.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.charMatrices.AlterData; import java.util.*;import java.awt.*;import java.awt.image.*;import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.lib.duties.*;import mesquite.lib.table.*;/* ======================================================================== */public class AlterData extends DataWindowAssistantI {	MesquiteTable table;	CharacterData data;	MesquiteSubmenuSpec mss= null;	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) {		mss = addSubmenu(null, "Alter/Transform", makeCommand("doAlter",  this));		mss.setList(DataAlterer.class);		return true;	}	/*.................................................................................................................*/	/** returns whether this module is requesting to appear as a primary choice */   	public boolean requestPrimaryChoice(){   		return true;     	}	/*.................................................................................................................*/   	 public boolean isPrerelease(){   	 	return false;   	 }	/*.................................................................................................................*/	public void setTableAndData(MesquiteTable table, CharacterData data, CommandRecord commandRec){		this.table = table;		this.data = data;		mss.setCompatibilityCheck(data.getStateClass());		resetContainingMenuBar();			}		/* possible alterers:		recode (as in MacClade)		reverse (sequences) �		shuffle �		random fill �		fill (via dialog) �		search and replace			DNA only:		complement �			continuous only:		log transform �		scale �		standardize (mean 0, variance 1)	*/	/*.................................................................................................................*/   	 public boolean isSubstantive(){   	 	return false;   	 }	/*.................................................................................................................*/    	 public Object doCommand(String commandName, String arguments, CommandRecord commandRec, CommandChecker checker) {    	 	if (checker.compare(this.getClass(), "Hires module to alter the data matrix", "[name of module]", commandName, "doAlter")) {   	 		if (table!=null && data !=null){	    	 	if (data.getEditorInhibition()){	    	 		discreetAlert(commandRec, "This matrix is marked as locked against editing.");	    	 		return null;	    	 	}	    	 	DataAlterer tda= (DataAlterer)hireNamedEmployee(commandRec, DataAlterer.class, arguments);				if (tda!=null) {					boolean a = tda.alterData(data, table, commandRec);	 	   			if (a) {	 	   				table.repaintAll();						data.notifyListeners(this, new Notification(MesquiteListener.DATA_CHANGED), commandRec);					}					if (!tda.pleaseLeaveMeOn()) {						fireEmployee(tda);					}				}			}    	 	}    	 	else    	 		return super.doCommand(commandName, arguments, commandRec, checker);	return null;   	 }	/*.................................................................................................................*/    	 public String getName() {		return "Alter Data";   	 }	/*.................................................................................................................*/ 	/** returns an explanation of what the module does.*/ 	public String getExplanation() { 		return "Manages data-transforming modules." ;   	 }   	 }