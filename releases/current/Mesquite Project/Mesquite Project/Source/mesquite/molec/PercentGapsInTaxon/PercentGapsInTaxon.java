/* Mesquite source code.  Copyright 1997-2007 W. Maddison and D. Maddison.Version 2.0, September 2007.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html) */package mesquite.molec.PercentGapsInTaxon;/*~~  */import java.util.*;import java.awt.*;import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.lib.duties.*;/* ======================================================================== */public class PercentGapsInTaxon extends NumberForTaxon {	public void getEmployeeNeeds(){  //This gets called on startup to harvest information; override this and inside, call registerEmployeeNeed		EmployeeNeed e = registerEmployeeNeed(MatrixSourceCoord.class, getName() + "  needs a source of characters.",		"The source of characters is arranged initially");	}	MatrixSourceCoord matrixSourceTask;	Taxa currentTaxa = null;	MCharactersDistribution observedStates =null;	MesquiteBoolean countEdges;	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, boolean hiredByName) {		matrixSourceTask = (MatrixSourceCoord)hireEmployee(MatrixSourceCoord.class, "Source of character matrix (for percent gaps)"); 		if (matrixSourceTask==null)			return sorry(getName() + " couldn't start because no source of character matrices was obtained.");		countEdges = new MesquiteBoolean(true);		addCheckMenuItem(null, "Count Gaps at Edges", makeCommand("toggleEdges", this), countEdges);		return true;	}	/*.................................................................................................................*/	/** Generated by an employee who quit.  The MesquiteModule should act accordingly. */	public void employeeQuit(MesquiteModule employee) {		if (employee == matrixSourceTask)  // character source quit and none rehired automatically			iQuit();	}	/*.................................................................................................................*/	/** returns whether this module is requesting to appear as a primary choice */	public boolean requestPrimaryChoice(){		return true;  	}	/*.................................................................................................................*/	public Snapshot getSnapshot(MesquiteFile file) { 		Snapshot temp = new Snapshot();		temp.addLine("getMatrixSource", matrixSourceTask);		temp.addLine("toggleEdges " + countEdges.toOffOnString());		return temp;	}	MesquiteInteger pos = new MesquiteInteger();	/*.................................................................................................................*/	public Object doCommand(String commandName, String arguments, CommandChecker checker) {		if (checker.compare(this.getClass(), "Sets whether or not to count gaps at the edges of the matrix", "[on = cut; off]", commandName, "toggleEdges")) {			countEdges.toggleValue(parser.getFirstToken(arguments));			parametersChanged();		}		else if (checker.compare(this.getClass(), "Returns the matrix source", null, commandName, "getMatrixSource")) {			return matrixSourceTask;		}		else return  super.doCommand(commandName, arguments, checker);		return null;	}	/** Called to provoke any necessary initialization.  This helps prevent the module's intialization queries to the user from   	happening at inopportune times (e.g., while a long chart calculation is in mid-progress)*/	public void initialize(Taxa taxa){		currentTaxa = taxa;		matrixSourceTask.initialize(currentTaxa);	}	public void calculateNumber(Taxon taxon, MesquiteNumber result, MesquiteString resultString){		if (result==null)			return;	   	clearResultAndLastResult(result);		Taxa taxa = taxon.getTaxa();		int it = taxa.whichTaxonNumber(taxon);		if (taxa != currentTaxa || observedStates == null ) {			observedStates = matrixSourceTask.getCurrentMatrix(taxa);			currentTaxa = taxa;		}		if (observedStates==null)			return;		CharacterData data = observedStates.getParentData();		CharInclusionSet incl = null;		if (data !=null)			incl = (CharInclusionSet)data.getCurrentSpecsSet(CharInclusionSet.class);		int numChars = observedStates.getNumChars();		int charExc = 0;		boolean edges = countEdges.getValue();		if (numChars != 0) {			CharacterState cs = null;			boolean hitNonGap = edges;			int lastBase = numChars-1;			for (int ic=numChars-1; ic>=0 && !hitNonGap; ic--) {				if (incl == null || incl.isSelected(ic)){					cs = observedStates.getCharacterState(cs, ic, it);					if (!cs.isInapplicable()){						lastBase = ic;						hitNonGap = true;					}				}			}			numChars = lastBase +1;			int numGaps = 0;			int useNumChars = 0;			hitNonGap = edges;			for (int ic=0; ic<numChars; ic++) {				if (incl == null || incl.isSelected(ic)){					cs = observedStates.getCharacterState(cs, ic, it);					if (cs.isInapplicable()){						if (hitNonGap || edges){							numGaps++;							useNumChars++;						}					}					else {						useNumChars++;						hitNonGap = true;					}				}			}			if (useNumChars == 0)				result.setValue(1.0);			else				result.setValue(((double)numGaps)/useNumChars);		}			String exs = "";		if (charExc > 0)			exs = " (" + Integer.toString(charExc) + " characters excluded)";		if (resultString!=null)			resultString.setValue("Proportion of inapplicable codings (gaps) in matrix "+ observedStates.getName() + exs + ": " + result.toString());		saveLastResult(result);		saveLastResultString(resultString);	}	/*.................................................................................................................*/	public void employeeParametersChanged(MesquiteModule employee, MesquiteModule source, Notification notification) {		observedStates = null;		super.employeeParametersChanged(employee, source, notification);	}	/*.................................................................................................................*/	public String getName() {		return "Proportion Gaps in Taxon";  	}	/*.................................................................................................................*/	public boolean isPrerelease() {		return false;	}	public String getParameters() {		return "Proportion gaps in taxon in matrix from: " + matrixSourceTask.getParameters();	}	/*.................................................................................................................*/	/** returns an explanation of what the module does.*/	public String getExplanation() {		return "Reports the proportion of gaps (inapplicable codings) in a taxon for a data matrix." ;	}}