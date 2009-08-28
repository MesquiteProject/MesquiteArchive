/* Mesquite source code.  Copyright 1997-2009 W. Maddison and D. Maddison.Version 2.7, August 2009.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html) */package mesquite.stochchar.AsymmMkExplorer;import java.util.*;import java.awt.*;import java.awt.event.*;import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.lib.duties.*;import mesquite.stochchar.lib.AsymmModel;import mesquite.stochchar.lib.MargLikelihoodForModel;import mesquite.categ.lib.*;/* ======================================================================== */public class AsymmMkExplorer extends TreeWindowAssistantOA implements ParametersExplorable {	public void getEmployeeNeeds(){  //This gets called on startup to harvest information; override this and inside, call registerEmployeeNeed		EmployeeNeed e = registerEmployeeNeed(mesquite.stochchar.zMargLikeCateg.zMargLikeCateg.class, getName() + "  needs a module to calculate likelihoods.",		"The module to calculate likelihoods is arranged automatically");		EmployeeNeed e2 = registerEmployeeNeed(CharSourceCoordObed.class, getName() + "  needs a source of characters.",		"The source of characters is arranged initially");	}	//Vector traces;	CharSourceCoordObed characterSourceTask;	ParametersExplorer explorer;	int currentChar = 0;	Taxa taxa;	MargLikelihoodForModel likelihoodTask = null;	MesquiteBoolean showAsRateBias = new MesquiteBoolean(true);	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, boolean hiredByName) {		//traces = new Vector();		explorer = (ParametersExplorer)hireEmployee(ParametersExplorer.class, "Parameters explorer");		if (explorer == null)			return sorry(getName() + " couldn't start because no parameters explorer module was obtained.");		likelihoodTask = (MargLikelihoodForModel)hireNamedEmployee(MargLikelihoodForModel.class, "#zMargLikeCateg");		if (likelihoodTask == null)			return sorry(getName() + " couldn't start because no likelihood calculator obtained");		characterSourceTask = (CharSourceCoordObed)hireCompatibleEmployee(CharSourceCoordObed.class, new RequiresExactlyCategoricalData(), "Source of characters (for " + getName() + ")");		if (characterSourceTask == null)			return sorry(getName() + " couldn't start because no source of characters was obtained.");		explorer.addCheckMenuItem( null, "Use Rate-Bias Notation", makeCommand("toggleRateBias",  this), showAsRateBias);		explorer.addMenuItem( "Next Character", makeCommand("nextCharacter",  this));		explorer.addMenuItem( "Previous Character", makeCommand("previousCharacter",  this));		explorer.addMenuItem( "Choose Character...", makeCommand("chooseCharacter",  this));		resetContainingMenuBar();		if (!MesquiteThread.isScripting()){			currentChar = characterSourceTask.queryUserChoose(taxa, "Choose character for likelihood surface exploration");			if (currentChar <0 || !MesquiteInteger.isCombinable(currentChar))				currentChar = 0;		}		asymmModel = new AsymmModel("Estimating Asymm", CategoricalState.class);		asymmModel.setUseRateBiasNotation(true);		param0 = new MesquiteParameter("Rate", "Rate of state change", MesquiteDouble.unassigned, 0, MesquiteDouble.infinite, 0.000, 0.1);		param1 = new MesquiteParameter("Bias", "Bias of gains (0 to 1) versus loses (1 to 0)", MesquiteDouble.unassigned, 0, MesquiteDouble.infinite, 0.1, 10.0);		parameters = new MesquiteParameter[]{param0, param1};		explorer.setExplorable(this);		return true;	}	/*.................................................................................................................*/	/** Generated by an employee who quit.  The MesquiteModule should act accordingly. */	public void employeeQuit(MesquiteModule employee) {		iQuit();	}		/*.................................................................................................................*/	public Snapshot getSnapshot(MesquiteFile file) {		Snapshot temp = new Snapshot();		temp.addLine("getCharacterSource ",characterSourceTask);		temp.addLine("setCharacter " + CharacterStates.toExternal(currentChar));		temp.addLine("toggleRateBias " + showAsRateBias.toOffOnString());		return temp;	}	boolean suspendCommandReceived = false;	MesquiteInteger pos = new MesquiteInteger();	/*.................................................................................................................*/	public Object doCommand(String commandName, String arguments, CommandChecker checker) {		if (checker.compare(this.getClass(), "Returns module supplying characters", null, commandName, "getCharacterSource")) {			return characterSourceTask;		}		else if (checker.compare(this.getClass(), "Goes to next character", null, commandName, "toggleRateBias")) {			showAsRateBias.toggleValue(parser.getFirstToken(arguments));			setNotation(showAsRateBias.getValue());				if (!MesquiteThread.isScripting())					parametersChanged();		}		else if (checker.compare(this.getClass(), "Goes to next character", null, commandName, "nextCharacter")) {			if (currentChar>=characterSourceTask.getNumberOfCharacters(taxa)-1)				currentChar=0;			else				currentChar++;			//charStates = null;			parametersChanged(); //?		}		else if (checker.compare(this.getClass(), "Goes to previous character", null, commandName, "previousCharacter")) {			if (currentChar<=0)				currentChar=characterSourceTask.getNumberOfCharacters(taxa)-1;			else				currentChar--;			//charStates = null;			parametersChanged(); //?		}		else if (checker.compare(this.getClass(), "Queries the user about what character to use", null, commandName, "chooseCharacter")) {			int ic=characterSourceTask.queryUserChoose(taxa, " to calculate value for tree ");			if (MesquiteInteger.isCombinable(ic)) {				currentChar = ic;				//charStates = null;				parametersChanged(); //?			}		}		else if (checker.compare(this.getClass(), "Sets the character to use", "[character number]", commandName, "setCharacter")) {			int icNum = MesquiteInteger.fromFirstToken(arguments, stringPos);			if (!MesquiteInteger.isCombinable(icNum))				return null;			int ic = CharacterStates.toInternal(icNum);			if ((ic>=0) && characterSourceTask.getNumberOfCharacters(taxa)==0) {				currentChar = ic;				//charStates = null;			}			else if ((ic>=0) && (ic<=characterSourceTask.getNumberOfCharacters(taxa)-1)) {				currentChar = ic;				//charStates = null;				parametersChanged(); //?			}		}		else			return  super.doCommand(commandName, arguments, checker);		return null;	}	/*.................................................................................................................*/	public void employeeParametersChanged(MesquiteModule employee, MesquiteModule source, Notification notification) {		doCalculations(true);	}	AsymmModel asymmModel;	MesquiteParameter param0, param1;	MesquiteParameter[] parameters;	Tree tree;	boolean useRateBias = true;;	CharacterDistribution observedStates;	MesquiteNumber likelihood = new MesquiteNumber();	void setNotation(boolean useRateBias){		if (this.useRateBias == useRateBias)			return;		if (useRateBias){			param0.setName("Rate");			param0.setExplanation("Rate of state change");			param1.setName("Bias");			param1.setExplanation("Bias of gains (0 to 1) versus loses (1 to 0)");		}		else {			param0.setName("Forward Rate");			param0.setExplanation("Rate of 0 to 1 change");			param1.setName("Backward Rate");			param1.setExplanation("Rate of 1 to 0 change");		}		asymmModel.setUseRateBiasNotation(useRateBias);		this.useRateBias = useRateBias;	}	/*.................................................................................................................*/	public  void setTree(Tree tree){		taxa = tree.getTaxa();		this.tree = tree;		doCalculations(true);	}	public MesquiteParameter[] getExplorableParameters(){		return parameters;	}	public double calculate(MesquiteString resultString){		asymmModel.deassignParameters();		asymmModel.setParam0(param0.getValue());		asymmModel.setParam1(param1.getValue());		likelihoodTask.calculateLogProbability( tree,  observedStates, asymmModel, resultString, likelihood);		return likelihood.getDoubleValue();	}	public void restoreAfterExploration(){	}	MesquiteString resultString = new MesquiteString();	/*.................................................................................................................*/	public void doCalculations(boolean doPreps) {		observedStates = characterSourceTask.getCharacter(tree, currentChar);		explorer.explorableChanged(this);	}	/*.................................................................................................................*/	public void endJob() {		setModuleWindow(null);		super.endJob();	}	/*.................................................................................................................*/	public String getName() {		return "Likelihood Surface AsymmMk Model";	}	/*.................................................................................................................*/	/** returns the version number at which this module was first released.  If 0, then no version number is claimed.  If a POSITIVE integer	 * then the number refers to the Mesquite version.  This should be used only by modules part of the core release of Mesquite.	 * If a NEGATIVE integer, then the number refers to the local version of the package, e.g. a third party package*/	public int getVersionOfFirstRelease(){		return 200;  	}	/*.................................................................................................................*/	public boolean isPrerelease(){		return false;	}	/*.................................................................................................................*/	/** returns whether this module is requesting to appear as a primary choice */	public boolean requestPrimaryChoice(){		return true;  	}	/*.................................................................................................................*/	public String getExplanation() {		return "Shows the likelihood surface for the AsymmMk Model.";	}}