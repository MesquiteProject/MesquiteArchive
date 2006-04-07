/* Mesquite (package mesquite.parsimony).  Copyright 1997-2005 W. Maddison and D. Maddison. Version 1.06, August 2005.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.parsimony.TreelengthForMatrix;/*~~  */import java.util.*;import java.awt.*;import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.lib.duties.*;import mesquite.parsimony.lib.*;/* ======================================================================== */public class TreelengthForMatrix extends NumberForMatrixAndTree {	MesquiteNumber treelength;	CharacterSteps charStepsTask;	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) {    	 	treelength=new MesquiteNumber(); 		charStepsTask = (CharacterSteps)hireEmployee(commandRec, CharacterSteps.class, null); 		if (charStepsTask ==null) 			return sorry(commandRec, getName() + " couldn't start because no step counting module was obtained."); 		return true;  	 }  	 public void employeeQuit(MesquiteModule m){  	 		iQuit();  	 }  	 	/*.................................................................................................................*/	/** returns whether this module is requesting to appear as a primary choice */   	public boolean requestPrimaryChoice(){   		return true;     	}   	/** Called to provoke any necessary initialization.  This helps prevent the module's intialization queries to the user from   	happening at inopportune times (e.g., while a long chart calculation is in mid-progress)*/   	public void initialize(Tree tree, MCharactersDistribution matrix, CommandRecord commandRec){   	}	/*.................................................................................................................*/	public void calculateNumber(Tree tree, MCharactersDistribution matrix, MesquiteNumber result, MesquiteString resultString, CommandRecord commandRec) {    	 	if (result==null || tree == null || matrix == null)    	 		return;    	 	if (resultString !=null)    	 		resultString.setValue("");		result.setToUnassigned();		    	 	treelength.setValue((int)0);    	 	int count=0;     	 	Taxa taxa = tree.getTaxa(); 		MesquiteNumber cNum = new MesquiteNumber(); 		MesquiteNumber wt = new MesquiteNumber(); 		boolean someExcluded = false;		CharWeightSet weightSet = null;		CharacterData data = matrix.getParentData();		if (data !=null)			weightSet = (CharWeightSet)data.getCurrentSpecsSet(CharWeightSet.class);    	 	for (int ic=0;  ic<matrix.getNumChars(); ic++) {			if (matrix.isCurrentlyIncluded(ic)){				CharacterDistribution charStates = matrix.getCharacterDistribution(ic);				if (charStates!=null){					cNum.setToUnassigned();					charStepsTask.calculateNumber(tree, charStates, cNum, null, commandRec);					if (cNum.isCombinable()){						if (weightSet!=null) {							weightSet.placeValue(ic, wt);							if (wt.isCombinable())								cNum.multiplyBy(wt);						}						treelength.add(cNum);						count++;					}					else						someExcluded = true;				}			}			else				someExcluded = true;		}	    	String exclString = null;		if (someExcluded)			exclString = " (" + count + " characters included)";		else			exclString = "";		result.setValue(treelength);		if (resultString!=null)			resultString.setValue("Treelength: "+ result.toString() + exclString);	}	/*.................................................................................................................*/    	 public String getName() {		return "Treelength ";   	 }	/*.................................................................................................................*/  	 public boolean isPrerelease() {		return false;   	 }   	 public boolean showCitation(){   	 	return true;   	 }	/*.................................................................................................................*/  	 public String getExplanation() {		return "Calculates the parsimony treelength of a given tree and matrix.";   	 }}