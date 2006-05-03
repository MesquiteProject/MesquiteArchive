/* Mesquite source code.  Copyright 1997-2005 W. Maddison and D. Maddison. Version 1.06, August 2005.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.distance.lib;/*~~  */import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.categ.lib.*;import mesquite.cont.lib.GeographicStateTest;import mesquite.distance.lib.*;/* ======================================================================== *//* incrementable, with each being based on a different matrix */public abstract class DNATaxaDistFromMatrix extends TaxaDistFromMatrix {	MesquiteBoolean estimateAmbiguityDifferences = new MesquiteBoolean(true);	/*.................................................................................................................*/	public boolean superStartJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) {		super.superStartJob(arguments, condition, commandRec, hiredByName);		addCheckMenuItemToSubmenu(null, distParamSubmenu, "Estimate Ambiguity Differences", MesquiteModule.makeCommand("toggleEstimateAmbiguityDifferences", this), estimateAmbiguityDifferences);		return true;  	 }	 	/*.................................................................................................................*/	 public Snapshot getSnapshot(MesquiteFile file) {	 	Snapshot snapshot = new Snapshot();	 	snapshot.addLine("toggleEstimateAmbiguityDifferences  " + estimateAmbiguityDifferences.toOffOnString());	 	return snapshot;	 }	 /*.................................................................................................................*/	 public Object doCommand(String commandName, String arguments, CommandRecord commandRec, CommandChecker checker) {	 	if (checker.compare(this.getClass(), "Sets whether sites with missing data, gaps, or ambiguities have differences estimated by other sites or are ignored.", "[on; off]", commandName, "toggleEstimateAmbiguityDifferences")) {	 		estimateAmbiguityDifferences.toggleValue(new Parser().getFirstToken(arguments));	 		parametersChanged(null, commandRec);	 }		else	 		return super.doCommand(commandName, arguments, commandRec, checker);	 	return null;	 }	/*.................................................................................................................*/	 public boolean getEstimateAmbiguityDifferences() {		 return estimateAmbiguityDifferences.getValue();  	 } 	public Class getRequiredStateClass(){		return DNAState.class;	}	public CompatibilityTest getCompatibilityTest(){		return new DNAStateOnlyTest();	} 	public boolean isSubstantive(){   		return true;   	} }