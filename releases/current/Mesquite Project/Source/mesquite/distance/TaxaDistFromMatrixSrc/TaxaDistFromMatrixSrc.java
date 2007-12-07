/* Mesquite source code.  Copyright 1997-2007 W. Maddison and D. Maddison.Version 2.01, December 2007.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html) */package mesquite.distance.TaxaDistFromMatrixSrc;/*~~  */import java.util.*;import java.awt.*;import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.lib.duties.*;import mesquite.cont.lib.*;import mesquite.categ.lib.*;import mesquite.distance.lib.*;/* ======================================================================== *//* incrementable, with each being based on a different matrix */public class TaxaDistFromMatrixSrc extends IncTaxaDistanceSource {	public void getEmployeeNeeds(){  //This gets called on startup to harvest information; override this and inside, call registerEmployeeNeed		EmployeeNeed e = registerEmployeeNeed(MatrixSourceCoord.class, getName() + "  needs a source of matrices.",		"The source of matrices is arranged initially");		EmployeeNeed e2 = registerEmployeeNeed(MatrixSourceCoordObed.class, getName() + "  needs a source of matrices.",		"The source of matrices is arranged initially");		EmployeeNeed e3 = registerEmployeeNeed(TaxaDistFromMatrix.class, getName() + "  needs a method to calculate distances.",		"The method to calculate distances can be selected initially or in the Distance Calculator submenu");	}	MesquiteNumber num;	MatrixSourceCoord matrixSourceTask;	MatrixSourceCoordObed matrixSourceObedTask;	TaxaDistFromMatrix distanceTask;	Taxa currentTaxa;	boolean matrixSet = false;	int currentMatrix = 0;	MCharactersDistribution observedStates =null;	MesquiteString distanceTaskName;	MesquiteCommand mesqC;	MesquiteSubmenuSpec mss;	Class currentStateClass = null;	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, boolean hiredByName) {		distanceTask = (TaxaDistFromMatrix)hireEmployee(TaxaDistFromMatrix.class, "Distance calculator");		if (distanceTask == null) {			return sorry(getName() + " couldn't start because no distance calculator was obtained.");		}		mesqC = makeCommand("setDistanceTask",  this);		distanceTask.setHiringCommand(mesqC);		distanceTaskName = new MesquiteString(distanceTask.getName());		if (numModulesAvailable(TaxaDistFromMatrix.class)>1){			mss = addSubmenu(null, "Distance Calculator", mesqC, TaxaDistFromMatrix.class);			mss.setSelected(distanceTaskName);		}		Class stateClass = distanceTask.getRequiredStateClass();		if (stateClass == null){			if (getHiredAs() == IncTaxaDistanceSource.class)				matrixSourceObedTask = (MatrixSourceCoordObed)hireEmployee(MatrixSourceCoordObed.class, "Source of matrices (for " + getName() + ")");			else				matrixSourceTask = (MatrixSourceCoord)hireEmployee(MatrixSourceCoord.class, "Source of matrices (for " + getName() + ")");		}		else {			if (getHiredAs() == IncTaxaDistanceSource.class)				matrixSourceObedTask = (MatrixSourceCoordObed)hireCompatibleEmployee(MatrixSourceCoordObed.class, stateClass, "Source of matrices (for " + getName() + ")");			else				matrixSourceTask = (MatrixSourceCoord)hireCompatibleEmployee(MatrixSourceCoord.class, stateClass, "Source of matrices (for " + getName() + ")");		}		if (matrixSourceTask == null && matrixSourceObedTask == null) {			return sorry(getName() + " couldn't start because no source of a character matrix was obtained.");		}		return true;	}	/*.................................................................................................................*/	/** Generated by an employee who quit.  The MesquiteModule should act accordingly. */	public void employeeQuit(MesquiteModule employee) {		iQuit();	}	/*.................................................................................................................*/	public Snapshot getSnapshot(MesquiteFile file) {		Snapshot temp = new Snapshot();		if (matrixSourceTask != null)			temp.addLine( "getMatrixSource " , matrixSourceTask);		else if (matrixSourceObedTask != null)			temp.addLine( "getMatrixSource " , matrixSourceObedTask);		temp.addLine( "setDistanceTask " , distanceTask);		return temp;	}	/*.................................................................................................................*/	public Object doCommand(String commandName, String arguments, CommandChecker checker) {		if (checker.compare(this.getClass(), "Returns the module supplying matrices", null, commandName, "getMatrixSource")) {			if (matrixSourceTask !=null)				return matrixSourceTask;			return matrixSourceObedTask;		}		else if (checker.compare(this.getClass(), "Sets the module doing the distance calculation from a matrix", null, commandName, "setDistanceTask")) {			Class currentStateClass = null;			if (observedStates!=null)				currentStateClass = observedStates.getStateClass();			TaxaDistFromMatrix temp=  (TaxaDistFromMatrix)replaceCompatibleEmployee(TaxaDistFromMatrix.class, arguments, distanceTask,currentStateClass);			if (temp!=null) {				distanceTask= temp;				distanceTask.setHiringCommand(mesqC);				distanceTaskName.setValue(distanceTask.getName());				if (!MesquiteThread.isScripting())					parametersChanged(); 			}			return distanceTask;		}		else			return  super.doCommand(commandName, arguments, checker);	}	/** Called to provoke any necessary initialization.  This helps prevent the module's initialization queries to the user from   	happening at inopportune times (e.g., while a long chart calculation is in mid-progress)*/	public void initialize(Taxa taxa){		currentTaxa = taxa;		if (matrixSourceTask!=null)			matrixSourceTask.initialize(currentTaxa);		else if (matrixSourceObedTask!=null)			matrixSourceObedTask.initialize(currentTaxa);	}	public TaxaDistance getTaxaDistance(Taxa taxa){		if (matrixSourceTask!=null)			observedStates =matrixSourceTask.getCurrentMatrix(taxa);		else if (matrixSourceObedTask!=null)			observedStates =matrixSourceObedTask.getMatrix(taxa, currentMatrix);		currentTaxa = taxa;		if (observedStates==null) {			MesquiteMessage.warnProgrammer("Observed states null in " + getName() + " for taxa " + taxa); 			return null;		}		if (observedStates.getStateClass() !=currentStateClass){			currentStateClass = observedStates.getStateClass();			mss.setListableFilter(currentStateClass);			//mss.setListableFilter(null);			resetContainingMenuBar();		}		return distanceTask.getTaxaDistance(taxa, observedStates);	}	/*.................................................................................................................*/	public void setCurrent(long i){  //SHOULD NOT notify (e.g., parametersChanged)		currentMatrix = (int)i;	}	public long getCurrent(){		return currentMatrix;	}	public String getItemTypeName(){		return "Matrix";	}	public long getMin(){		return 0;	}	public long getMax(){		if (matrixSourceObedTask == null)			return 0;		return matrixSourceObedTask.getNumberOfMatrices(currentTaxa)-1;	}	public long toInternal(long i){		return i-1;	}	public long toExternal(long i){ //return whether 0 based or 1 based counting		return i+1;	}	/*.................................................................................................................*/	public void employeeParametersChanged(MesquiteModule employee, MesquiteModule source, Notification notification) {		observedStates=null;		super.employeeParametersChanged(this, source, notification);	}	/*.................................................................................................................*/	public String getName() {		return "Distances from Character Matrix";  	}	/*.................................................................................................................*/	public String getNameAndParameters() {		String s ="Distances: " + distanceTask.getName()+ " -- Matrix: " ;  		if (matrixSourceTask !=null)			s += matrixSourceTask.getParameters();		else if (matrixSourceObedTask !=null)			s += matrixSourceObedTask.getParameters();		return s;	}	public String getParameters() {		String s = "";		if (distanceTask !=null)			s += "Distances calculated by " + distanceTask.getNameAndParameters();		if (matrixSourceTask !=null)			s += " Matrix source: " + matrixSourceTask.getParameters();		else if (matrixSourceObedTask !=null)			s += " Matrix source: " + matrixSourceObedTask.getParameters();		return s;	}	/*.................................................................................................................*/	/** returns an explanation of what the module does.*/	public String getExplanation() {		return "Distances calculated from a character matrix." ;	}	public boolean isPrerelease(){		return false;	}	/*.................................................................................................................*/	public boolean showCitation(){		return false;	}}