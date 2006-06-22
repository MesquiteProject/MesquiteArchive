/* Mesquite source code.  Copyright 1997-2006 W. Maddison and D. Maddison.Version 1.11, June 2006.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.trees.StoredTreeBlocks;/*~~  */import java.util.*;import java.awt.*;import mesquite.lib.*;import mesquite.lib.duties.*;/** Supplies tree blocks stored in the projects.*/public class StoredTreeBlocks extends TreeBlockSource implements MesquiteListener {	int currentTreeBlockIndex=MesquiteInteger.unassigned;	TreeVector currentTreeBlock = null;	TreeVector lastUsedTreeBlock = null;	TreesManager manager;	Taxa preferredTaxa =null;	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) {		manager = (TreesManager)findElementManager(TreeVector.class);		if (manager==null)			return sorry(commandRec, getName() + " couldn't start because no tree manager module was found.");		if (manager.getNumberTreeBlocks()==0 && !commandRec.scripting())			return sorry(commandRec, "No stored blocks of trees are available.");		return true;  	 }	/*.................................................................................................................*/   	 public boolean isSubstantive(){   	 	return true;   	 }	/*.................................................................................................................*/   	 public boolean isPrerelease(){   	 	return false;   	 }	/*.................................................................................................................*/	/** returns whether this module is requesting to appear as a primary choice */   	public boolean requestPrimaryChoice(){   		return true;     	}	/*.................................................................................................................*/	/** passes which object changed*/	public void changed(Object caller, Object obj, Notification notification, CommandRecord commandRec){		int code = Notification.getCode(notification);		if (!doomed && code != MesquiteListener.SELECTION_CHANGED && code != MesquiteListener.ANNOTATION_CHANGED && code != MesquiteListener.ANNOTATION_DELETED && code != MesquiteListener.ANNOTATION_ADDED)			parametersChanged(notification, commandRec);	}		/*.................................................................................................................*/	/** passes which object disposed*/	public void disposing(Object obj){		if (obj == preferredTaxa) {			setHiringCommand(null); //since there is no rehiring			iQuit();		}	}	/** Asks whether it's ok to delete the object as far as the listener is concerned (e.g., is it in use?)*/	public boolean okToDispose(Object obj, int queryUser){		return true;	}	/*.................................................................................................................*/  	public void setPreferredTaxa(Taxa taxa) {  		if (preferredTaxa != taxa){	  		if (preferredTaxa != null)	  			preferredTaxa.removeListener(this);	  		taxa.addListener(this);  		}  		preferredTaxa = taxa;  	}   	/** Called to provoke any necessary initialization.  This helps prevent the module's intialization queries to the user from   	happening at inopportune times (e.g., while a long chart calculation is in mid-progress)*/   	public void initialize(Taxa taxa, CommandRecord commandRec){   	}	/*.................................................................................................................*/   	public TreeVector getFirstBlock(Taxa taxa, CommandRecord commandRec) {   		currentTreeBlockIndex=0;   		return getCurrentBlock(taxa, commandRec);   	}	/*.................................................................................................................*/   	public TreeVector getBlock(Taxa taxa, int ic, CommandRecord commandRec) {   		currentTreeBlockIndex=ic;   		return getCurrentBlock(taxa, commandRec);   	}   	private void checkBlock(Taxa taxa, CommandRecord commandRec){		int nt = manager.getNumberTreeBlocks(taxa);		setPreferredTaxa(taxa);		if ((!MesquiteInteger.isCombinable(currentTreeBlockIndex) || currentTreeBlockIndex>=nt || currentTreeBlockIndex<0)) {			if (commandRec.scripting())				currentTreeBlockIndex = 0;			else if (nt<=1)				currentTreeBlockIndex = 0;			else {				String[] list = new String[nt];				for (int i=0; i< nt; i++)					list[i]=manager.getTreeBlock(taxa, i).getName();				currentTreeBlockIndex = ListDialog.queryList(containerOfModule(), "Use which tree block?", "Use which tree block? \n(for " + employer.getName() + ")", MesquiteString.helpString,list, 0);				if (!MesquiteInteger.isCombinable(currentTreeBlockIndex))					currentTreeBlockIndex = 0;			}		} 		currentTreeBlock = manager.getTreeBlock(taxa, currentTreeBlockIndex);   		if (currentTreeBlock!=lastUsedTreeBlock) {   			currentTreeBlock.addListener(this);   			if (lastUsedTreeBlock!=null)   				lastUsedTreeBlock.removeListener(this);   			lastUsedTreeBlock = currentTreeBlock;   		}   	}	/*.................................................................................................................*/   	public TreeVector getCurrentBlock(Taxa taxa, CommandRecord commandRec) {   		if (currentTreeBlockIndex>getNumberOfTreeBlocks(taxa, commandRec) || currentTreeBlockIndex<0)   			return  null; 		checkBlock(taxa, commandRec);   		return currentTreeBlock;   	}	/*.................................................................................................................*/   	public TreeVector getNextBlock(Taxa taxa, CommandRecord commandRec) {   		currentTreeBlockIndex++;   		return getCurrentBlock(taxa, commandRec);   	}	/*.................................................................................................................*/   	public int getNumberOfTreeBlocks(Taxa taxa, CommandRecord commandRec) {		return manager.getNumberTreeBlocks(taxa);   	}   	/*.................................................................................................................*/   	public String getTreeBlockNameString(Taxa taxa, int index, CommandRecord commandRec) {		return manager.getTreeBlock(taxa, index).getName();   	}	/*.................................................................................................................*/   	public String getCurrentTreeBlockNameString(Taxa taxa, CommandRecord commandRec) { 		checkBlock(taxa, commandRec); 		return currentTreeBlock.getName();  	}	/*.................................................................................................................*/    	 public String getName() {		return "Stored Tree Blocks";   	 }   	 	/*.................................................................................................................*/  	 public String getExplanation() {		return "Supplies lists of trees stored, for instance in a file.";   	 }}