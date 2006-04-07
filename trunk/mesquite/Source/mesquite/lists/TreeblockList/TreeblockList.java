/* Mesquite (package mesquite.lists).  Copyright 1997-2005 W. Maddison and D. Maddison. Version 1.06, August 2005.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.lists.TreeblockList;/*~~  */import mesquite.lists.lib.*;import java.util.*;import java.awt.*;import mesquite.lib.*;import mesquite.lib.duties.*;import mesquite.lib.table.*;/* ======================================================================== */public class TreeblockList extends ListLVModule {	TreesManager manager;	ListableVector treeBlocks;	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) {		manager = (TreesManager)findElementManager(TreeVector.class);		if (manager==null)			return sorry(commandRec, getName() + " couldn't start because no tree manager module was found");		return true;  	 }  	 public boolean showing(Object obj){    	 		return (getModuleWindow()!=null && obj == treeBlocks);  	 }  	 	public boolean rowsDeletable(){		return true;	}	public boolean resetMenusOnNameChange(){		return true;	}	public boolean deleteRow(int row){		if (row<0 || row>= getNumberOfRows())			return false;		TreeVector trees = (TreeVector)treeBlocks.elementAt(row);		getProject().removeFileElement(trees);//must remove first, before disposing		trees.dispose();		return true;	}	/*.................................................................................................................*/    	 public Object doCommand(String commandName, String arguments, CommandRecord commandRec, CommandChecker checker) {    	 	if (checker.compare(this.getClass(), "Hires utility module to operate on the trees blocks", "[name of module]", commandName, "doUtility")) {    	 		TreeBlockListUtility tda= (TreeBlockListUtility)hireNamedEmployee(commandRec, TreeBlockListUtility.class, arguments);			if (tda!=null) {				int count = 0;				for (int i=0; i<treeBlocks.size() && i<getNumberOfRows(); i++)					if (isRowSelected(i))						count++;				boolean doAll = false;				if (count == 0) {					count = treeBlocks.size();					doAll = true;				}				TreeVector[] selected = new TreeVector[count];				count = 0;				for (int i=0; i<treeBlocks.size() && i<getNumberOfRows(); i++)					if (doAll || isRowSelected(i))						selected[count++] = (TreeVector)treeBlocks.elementAt(i);				boolean a = tda.operateOnTreeBlocks(selected, commandRec);				if (!tda.pleaseLeaveMeOn())					fireEmployee(tda);			}    	 	}    	 	else    	 		return super.doCommand(commandName, arguments, commandRec, checker);		return null;   	 }  	 public void showListWindow(Object obj, CommandRecord commandRec){		setModuleWindow( new ListableVectorWindow(this));		treeBlocks = manager.getTreeBlockVector();     	 	((ListableVectorWindow)getModuleWindow()).setObject(treeBlocks, commandRec);		//makeMenu("Tree_Blocks");		makeMenu("List");		MesquiteSubmenuSpec mss2 = addSubmenu(null, "Utilities", MesquiteModule.makeCommand("doUtility",  this));		mss2.setList(TreeBlockListUtility.class);				if (!commandRec.scripting()){			TreeblocksListAssistant assistant = (TreeblocksListAssistant)hireNamedEmployee(commandRec, TreeblocksListAssistant.class, "#TreeblocksListNumber");	 		if (assistant!= null){	 			((ListableVectorWindow)getModuleWindow()).addListAssistant(assistant, commandRec);				assistant.setUseMenubar(false);			}			assistant = (TreeblocksListAssistant)hireNamedEmployee(commandRec, TreeblocksListAssistant.class, "#TreeblocksListTaxa");	 		if (assistant!= null){	 			((ListableVectorWindow)getModuleWindow()).addListAssistant(assistant, commandRec);				assistant.setUseMenubar(false);			}		} 		resetContainingMenuBar();		resetAllWindowsMenus();  	 }	/*.................................................................................................................*/	/* following required by ListModule*/  	 public Object getMainObject(){  	 	return treeBlocks;  	 }  	 public int getNumberOfRows(){  	 	if (treeBlocks==null)  	 		return 0;  	 	else  	 		return treeBlocks.size();  	 }	public Class getAssistantClass(){		return TreeblocksListAssistant.class;	}	public String getItemTypeName(){		return "Tree block";	}	public String getItemTypeNamePlural(){		return "Tree blocks";	}	/*.................................................................................................................*/  	 public Snapshot getSnapshot(MesquiteFile file) {    	 	if (getModuleWindow()==null || !getModuleWindow().isVisible())   	 		return null;   	 	Snapshot temp = new Snapshot();      	 	if (getModuleWindow()!=null)			getModuleWindow().incorporateSnapshot(temp, file);  	 	temp.addLine("showWindow");  	 	return temp;  	 }	/*.................................................................................................................*/    	 public String getName() {		return "Tree Blocks List";   	 }   	 	/*.................................................................................................................*/   	  	/** returns an explanation of what the module does.*/ 	public String getExplanation() { 		return "Makes windows listing tree blocks and information about them." ;   	 }	/*.................................................................................................................*/	/** Requests a window to close.  In the process, subclasses of MesquiteWindow might close down their owning MesquiteModules etc.*/ 	public void windowGoAway(MesquiteWindow whichWindow) {			whichWindow.hide();	}   	    	 }