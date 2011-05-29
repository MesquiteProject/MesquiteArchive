/* Mesquite source code.  Copyright 1997-2010 W. Maddison and D. Maddison. 
Version 2.73, July 2010.
Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. 
The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.
Perhaps with your help we can be more than a few, and make Mesquite better.

Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
Mesquite's web site is http://mesquiteproject.org

This source code and its compiled class files are free and modifiable under the terms of 
GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)
 */package mesquite.align.lib;

import mesquite.align.lib.*;
import mesquite.categ.lib.*;
import mesquite.lib.*;
import mesquite.lib.characters.*;
import mesquite.lib.duties.*;

/* TODO: 
 *   allow user to alter pairwise alignment parameters
 */

public abstract class AlignScoreForTaxonGen extends NumberForTaxon {
	protected MatrixSourceCoord matrixSourceTask;
	protected Taxa currentTaxa = null;
	protected MCharactersDistribution observedStates =null;
	protected PairwiseAligner aligner;
	protected int alphabetLength;
	protected MesquiteInteger comparisonTaxon = new MesquiteInteger(0);
	//MesquiteTimer timer;
	/*.................................................................................................................*/
	public boolean startJob(String arguments, Object condition, boolean hiredByName) {
		matrixSourceTask = (MatrixSourceCoord)hireEmployee(MatrixSourceCoord.class, "Source of character matrix (for number of stops)"); 
		if (matrixSourceTask==null)
			return sorry(getName() + " couldn't start because no source of character matrices was obtained.");
		//timer =  new MesquiteTimer();
		addMenuItem("Reference Taxon...", MesquiteModule.makeCommand("setReferenceTaxon", this));
		/*
		 pairwiseTask = (TwoSequenceAligner)hireEmployee(TwoSequenceAligner.class, "Pairwise Aligner");
		 if (pairwiseTask == null)
		 return sorry(getName() + " couldn't start because no pairwise aligner obtained.");
		 */		return true;
	}
	
	/*.................................................................................................................*/
	public Snapshot getSnapshot(MesquiteFile file) {
		Snapshot temp = new Snapshot();
		if (comparisonTaxon.getValue()!=0)
			temp.addLine("setReferenceTaxon " + comparisonTaxon.getValue());
		return temp;
	}
	/*.................................................................................................................*/
	/** Generated by an employee who quit.  The MesquiteModule should act accordingly. */
	public void employeeQuit(MesquiteModule employee) {
		if (employee == matrixSourceTask)  // character source quit and none rehired automatically
			iQuit();
	}
	/*.................................................................................................................*/
	/** Override if one wishes to modify the alignment costs away from the default. */
	public int[][] modifyAlignmentCosts(int[][] defaultSubs) {
		return defaultSubs;
	}
	/*.................................................................................................................*/
	public void initAligner() {
		MesquiteInteger gapOpen = new MesquiteInteger();
		MesquiteInteger gapExtend = new MesquiteInteger();
		MesquiteInteger gapOpenTerminal = new MesquiteInteger();
		MesquiteInteger gapExtendTerminal = new MesquiteInteger();
		AlignUtil.getDefaultGapCosts(gapOpen, gapExtend, gapOpenTerminal, gapExtendTerminal);  
		int subs[][] = AlignUtil.getDefaultSubstitutionCosts(alphabetLength);  
		subs = modifyAlignmentCosts(subs);
		aligner = new PairwiseAligner(false,subs,gapOpen.getValue(), gapExtend.getValue(), gapOpenTerminal.getValue(), gapExtendTerminal.getValue(), alphabetLength);
		aligner.setUseLowMem(true);
	}
	/*.................................................................................................................*/
	/** returns whether this module is requesting to appear as a primary choice */
	public boolean requestPrimaryChoice(){
		return false;  
	}
	
	/** Called to provoke any necessary initialization.  This helps prevent the module's intialization queries to the user from
	 happening at inopportune times (e.g., while a long chart calculation is in mid-progress)*/
	public void initialize(Taxa taxa){
		currentTaxa = taxa;
		matrixSourceTask.initialize(currentTaxa);
		
	}
	
	public void calculateNumber(Taxon taxon, MesquiteNumber result, MesquiteString resultString){
		if (result==null)
			return;
		clearResultAndLastResult(result);
		if (taxon == null || comparisonTaxon == null)
			return;
		Taxa taxa = taxon.getTaxa();
		if (comparisonTaxon.getValue()<0 || comparisonTaxon.getValue()>=taxa.getNumTaxa())
			comparisonTaxon.setValue(0);
		int it = taxa.whichTaxonNumber(taxon);
		if (taxa != currentTaxa || observedStates == null ) {
			observedStates = matrixSourceTask.getCurrentMatrix(taxa);
			currentTaxa = taxa;
		}
		if (observedStates==null || !(observedStates.getParentData() instanceof DNAData))
			return;
		DNAData data = (DNAData)observedStates.getParentData();
		if (data == null)
			return;
		MesquiteNumber score = new MesquiteNumber();
		
		
		if (aligner==null) {
			alphabetLength = ((CategoricalState)data.makeCharacterState()).getMaxPossibleState()+1;	  
			initAligner();
		}
		
		getAlignmentScore(data, (MCategoricalDistribution)observedStates,comparisonTaxon.getValue(),it,score);
		
		
		if (result !=null)
			result.setValue(score);
		if (resultString!=null)
			resultString.setValue(getScoreName() + " of sequence in matrix "+ observedStates.getName() + ": " + score.getIntValue());
		saveLastResult(result);
		saveLastResultString(resultString);
	}
	/*.................................................................................................................*/
	protected abstract void getAlignmentScore(DNAData data, MCategoricalDistribution observedStates, int it1, int it2, MesquiteNumber score) ;
	
	/*.................................................................................................................*/
	public boolean queryReferenceTaxon() {
		MesquiteInteger buttonPressed = new MesquiteInteger(1);
		ExtensibleDialog dialog = new ExtensibleDialog(containerOfModule(), "Reference Taxon",buttonPressed);  //MesquiteTrunk.mesquiteTrunk.containerOfModule()
		
		int maxNum = MesquiteInteger.unassigned;
		if (currentTaxa!=null)
			maxNum = currentTaxa.getNumTaxa();
		IntegerField refTaxonField = dialog.addIntegerField("Reference Taxon", comparisonTaxon.getValue()+1, 8, 1, maxNum);
		
		dialog.completeAndShowDialog(true);
		if (buttonPressed.getValue()==0)  {
			comparisonTaxon.setValue(refTaxonField.getValue()-1);
			storePreferences();
		}
		dialog.dispose();
		return (buttonPressed.getValue()==0);
	}
	/*.................................................................................................................*/
	public Object doCommand(String commandName, String arguments, CommandChecker checker) {
		if (checker.compare(this.getClass(), "Allows one to specify which taxon is the reference taxon.", "[taxon number]", commandName, "setReferenceTaxon")) {
			MesquiteInteger io = new MesquiteInteger(0);
			int newRefTaxon = MesquiteInteger.fromString(arguments, io);
			if (newRefTaxon<0 ||  !MesquiteInteger.isCombinable(newRefTaxon)){
				queryReferenceTaxon();
			}
			else{
				comparisonTaxon.setValue(newRefTaxon);
			}
			parametersChanged();
			
		}
		else
			return  super.doCommand(commandName, arguments, checker);
		return null;
	}
	/*.................................................................................................................*/
	/** Returns CompatibilityTest so other modules know if this is compatible with some object. */
	public CompatibilityTest getCompatibilityTest(){
		return new RequiresAnyMolecularData();
	}
	/*.................................................................................................................*/
	public void employeeParametersChanged(MesquiteModule employee, MesquiteModule source, Notification notification) {
		observedStates = null;
		super.employeeParametersChanged(employee, source, notification);
	}
	/*.................................................................................................................*/
	public abstract String getScoreName() ;
	/*.................................................................................................................*/
	public boolean isSubstantive(){
		return true;
	}
	
	public String getParameters() {
		return getScoreName() + " of sequence in matrix from: " + matrixSourceTask.getParameters();
	}
	/*.................................................................................................................*/
	
	/** returns an explanation of what the module does.*/
	public String getExplanation() {
		return "Reports the " + getScoreName() + " for a taxon in a data matrix." ;
	}
	
}


