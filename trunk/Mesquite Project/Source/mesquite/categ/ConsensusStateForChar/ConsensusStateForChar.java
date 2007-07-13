package mesquite.categ.ConsensusStateForChar;


import mesquite.categ.lib.*;
import mesquite.lib.*;
import mesquite.lib.table.*;

public class ConsensusStateForChar extends CategStateForCharacter {
	MesquiteBoolean showOnlyModal = new MesquiteBoolean(false);
	double nonModalThreshold = 0.2;

	MesquiteMenuItemSpec menuItem1, menuItem2;

	public boolean startJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) {
		menuItem1= addCheckMenuItem(null,"Show Modal Value Only", makeCommand("toggleShowOnlyModal", this), showOnlyModal);
		menuItem2= addMenuItem(null,"Threshold for Non-Modal States...", makeCommand("setNonModalThreshold", this));
		return true;
	}
	/*.................................................................................................................*/
	public Snapshot getSnapshot(MesquiteFile file) {
		Snapshot temp = new Snapshot();

		temp.addLine("toggleShowOnlyModal " + showOnlyModal.toOffOnString());
		temp.addLine("setNonModalThreshold " + nonModalThreshold);

		return temp;
	}

	MesquiteInteger pos = new MesquiteInteger();
	/*.................................................................................................................*/
	public Object doCommand(String commandName, String arguments, CommandRecord commandRec, CommandChecker checker) {
		if (checker.compare(this.getClass(), "Sets the minimum frequency among applicable states required to include state in consensus", "[threshold value]", commandName, "setNonModalThreshold")) {
			double currentValue = nonModalThreshold;
			pos.setValue(0);
			double t = MesquiteDouble.fromString(arguments, pos);

			if (!MesquiteDouble.isCombinable(t) && !MesquiteThread.isScripting()) 
				t = MesquiteDouble.queryDouble(containerOfModule(), "Threshold for Non-Modal States", "Minimum frequency among applicable states required to include state in consensus:", nonModalThreshold);
			if (MesquiteDouble.isCombinable(t))
				nonModalThreshold = t;
			if (currentValue!=nonModalThreshold)
				parametersChanged(null, commandRec);
		}
		
		else if (checker.compare(this.getClass(), "Sets whether or not only to show only the modal value.", "[on or off]", commandName, "toggleShowOnlyModal")) {
			boolean current = showOnlyModal.getValue();
			showOnlyModal.toggleValue(parser.getFirstToken(arguments));
			if (current!=showOnlyModal.getValue()) {
				parametersChanged(null, commandRec);
			}
		}
		else
			return super.doCommand(commandName, arguments, commandRec, checker);
		return null;
	}
	/*.................................................................................................................*/
	public long getConsensusState(CategoricalData data, int ic, MesquiteTable table){
		int maxState = data.getMaxState(ic);
		double[] freq = new double[CategoricalState.getMaxPossibleStateStatic()];
		int numTaxa = data.getNumTaxa();
		int numTaxaWithData = 0;
		for (int it=0; it<numTaxa; it++) {
			if (!getSelectedOnly() || table!=null && (table.isRowSelected(it) || !table.anyRowSelected())) {
				long s= data.getState(ic,it);
				if (CategoricalState.isInapplicable(s))
					continue;
				int card = CategoricalState.cardinality(s);
				if (card==1)
					freq[CategoricalState.getOnlyElement(s)]+=1.0;
				else
					for (int iState = 0; iState<=maxState; iState++){
						if (CategoricalState.isElement(s,iState))
							freq[iState]+=1.0/card;
					}
				numTaxaWithData++;
			}
		}
//		now let's figure out the frequency of the most common state
		double max = 0.0;
		for (int iState = 0; iState<=maxState; iState++){
			if (freq[iState]>max) {
				max = freq[iState];
			}
		}

//		now let's figure out which states are close enough to the most common state.
		long consensus = CategoricalState.emptySet();
		if (showOnlyModal.getValue()) {
			for (int iState = 0; iState<=maxState; iState++){
				if (freq[iState]>=max) {
					consensus = CategoricalState.addToSet(consensus,iState);
				}
			}
		}
		else 
			for (int iState = 0; iState<=maxState; iState++){
				if (freq[iState]>nonModalThreshold*numTaxaWithData) {
					consensus = CategoricalState.addToSet(consensus,iState);
				}
			}
		if (CategoricalState.cardinality(consensus) > 1)
			consensus = CategoricalState.setUncertainty(consensus,true);
		else  if (CategoricalState.cardinality(consensus) ==0)
			consensus = CategoricalState.inapplicable;
		return consensus;
	}


	public void calculateState(CategoricalData data, int ic,MesquiteTable table, CategoricalState resultState, MesquiteString resultString) {
		if (data==null || resultState==null)
			return;
		long s = getConsensusState(data, ic, table);
		resultState.setValue(s);
		resultString.setValue(resultState.toString());
	}

	public String getName() {
		return "Consensus State for Character";
	}

}
