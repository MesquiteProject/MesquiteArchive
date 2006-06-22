/* Mesquite source code.  Copyright 1997-2006 W. Maddison and D. Maddison.Version 1.11, June 2006.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.charts.TreeBlockValues;/*~~  */import java.awt.*;import java.util.*;import mesquite.lib.*;import mesquite.lib.duties.*;/* ======================================================================== *//**=== Class TreeBlocksChart.  ===*/public class TreeBlockValues extends FileAssistantCH {	public NumberForTreeBlock numberTask;	private TreeBlockSource treeBlockSourceTask;	private Taxa taxa;	MesquiteString treeBlockSourceName;	MesquiteString numberTaskName;	MesquiteSubmenuSpec msNT;	ChartWindow cWindow;	ItemsCharter chartWindowTask;	MesquiteBoolean live = new MesquiteBoolean(true);	MesquiteCommand tbstC, ntC;	static int numMade = 0;	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) {		chartWindowTask = (ItemsCharter)hireEmployee(commandRec, ItemsCharter.class, null);		if (chartWindowTask == null)			return sorry(commandRec, getName() + " couldn't start because no charting module was obtained.");		makeMenu("Chart");						//Tree block source %%%%%%%%		treeBlockSourceName = new MesquiteString();		tbstC =makeCommand("setTreeBlockSource",  this);		if (numModulesAvailable(TreeBlockSource.class)>1) {			MesquiteSubmenuSpec mss = addSubmenu(null, "Tree Block Source", tbstC, TreeBlockSource.class);			mss.setSelected(treeBlockSourceName);		}		taxa = getProject().chooseTaxa(containerOfModule(), "For which block of taxa do you want to show a chart of values for tree blocks?", commandRec); 		if (taxa ==null) 			return sorry(commandRec, getName() + " couldn't start because no taxa block was obtained."); 		taxa.addListener(this);				treeBlockSourceTask= (TreeBlockSource)hireEmployee(commandRec, TreeBlockSource.class, "Source of tree blocks");		if (treeBlockSourceTask == null)			return sorry(commandRec, getName() + " couldn't start because no source of tree blocks was obtained.");		treeBlockSourceTask.setHiringCommand(tbstC);		treeBlockSourceName.setValue(treeBlockSourceTask.getName());		treeBlockSourceTask.setPreferredTaxa(taxa);						//values etc.  %%%%%%%%		ntC =makeCommand("setCalculator",  this);		msNT = addSubmenu(null, "Values", ntC, NumberForTreeBlock.class);		numberTask=(NumberForTreeBlock)hireEmployee(commandRec, NumberForTreeBlock.class, "Value to calculate for tree blocks");		if (numberTask == null)			return sorry(commandRec, getName() + " couldn't start because no calculator module was obtained.");		numberTask.setHiringCommand(ntC);		numberTaskName = new MesquiteString(numberTask.getName());		msNT.setSelected(numberTaskName);			    	 			cWindow = chartWindowTask.makeChartWindow(this, commandRec);		setModuleWindow( cWindow);		chartWindowTask.setTaxa(taxa);		chartWindowTask.setNumberTask(numberTask);		chartWindowTask.setItemsSource(treeBlockSourceTask);		cWindow.setChartTitle("Tree Blocks Chart " + (++numMade));		cWindow.resetTitle();		if (!commandRec.scripting()){			cWindow.setChartVisible();			cWindow.setVisible(true);			chartWindowTask.doCounts(commandRec);		} 		resetContainingMenuBar();		resetAllWindowsMenus();		//window.setTitle("Trees: " + treeSourceTask.getSourceName());		return true; 	}	/*.................................................................................................................*/	/** returns whether this module is requesting to appear as a primary choice */   	public boolean requestPrimaryChoice(){   		return true;     	}	/*.................................................................................................................*/	public boolean isSubstantive(){		return false;	} 	public void employeeQuit(MesquiteModule m){ 		if (m==chartWindowTask) 			iQuit(); 	}	/*.................................................................................................................*/	public void endJob(){			if (taxa!=null)				taxa.removeListener(this);			super.endJob();	}	/*.................................................................................................................*/	/** passes which object is being disposed (from MesquiteListener interface)*/	public void disposing(Object obj){		if (obj instanceof Taxa && (Taxa)obj == taxa) {			iQuit();		}	}	/*.................................................................................................................*/ 	public void windowGoAway(MesquiteWindow whichWindow) {			whichWindow.hide();			whichWindow.dispose();			iQuit();	}	/*.................................................................................................................*/  	 public Snapshot getSnapshot(MesquiteFile file) {    	 	Snapshot temp = new Snapshot();  	 	temp.addLine("suspendCalculations");   	 	temp.addLine("setTaxa " + getProject().getTaxaReferenceExternal(taxa));  	 	temp.addLine("setTreeBlockSource ", treeBlockSourceTask);    	 	temp.addLine("setCalculator ", numberTask);   	 	temp.addLine("getCharter", chartWindowTask);   	 	temp.addLine("setChartVisible");   	 	temp.addLine("resumeCalculations");   	 	temp.addLine("doCounts");   	 	temp.addLine("showWindow");  	 	return temp;  	 }	MesquiteInteger pos = new MesquiteInteger();	/*.................................................................................................................*/    	 public Object doCommand(String commandName, String arguments, CommandRecord commandRec, CommandChecker checker) {     	 	 if (checker.compare(getClass(), "Sets block of taxa to use", "[block reference, number, or name]", commandName, "setTaxa")){   	 		Taxa t = getProject().getTaxa(checker.getFile(), parser.getFirstToken(arguments));   	 		if (t!=null){	   	 		if (taxa!=null)	   	 			taxa.removeListener(this);	   	 		taxa = t;	   	 		if (taxa!=null)	   	 			taxa.addListener(this);	   	 		if (taxa!=null && treeBlockSourceTask!=null)	   	 			treeBlockSourceTask.setPreferredTaxa(taxa);	   	 		return taxa;   	 		}      	 	 }      	 	else if (checker.compare(this.getClass(), "Requests that chart values be recalculated", null, commandName, "doCounts")){			if (cWindow!=null)				chartWindowTask.doCounts(commandRec);     	 	}     	 	else if (checker.compare(this.getClass(), "Suspends calculations", null, commandName, "suspendCalculations")){			chartWindowTask.incrementSuspension(commandRec);     	 	}     	 	else if (checker.compare(this.getClass(), "Resumes calculations", null, commandName, "resumeCalculations")){			chartWindowTask.decrementSuspension(commandRec);     	 	}    	 	else if (checker.compare(this.getClass(), "Sets whether or not the chart is \"live\" (that is, whether its calculatons are updated automatically when parameters or data change)", "[on or off]", commandName, "setLive")) {    	 		live.toggleValue(parser.getFirstToken(arguments));    	 	}    	 	else if (checker.compare(this.getClass(), "Returns the chart drawing module", null, commandName, "getCharter")) {       			return chartWindowTask;    	 	}    	 	else if (checker.compare(this.getClass(), "Sets the chart to visible", null, commandName, "setChartVisible")) {			if (cWindow!=null)				cWindow.setChartVisible();    	 		    	 	}    	 	else if (checker.compare(getClass(), "Sets the module supplying tree blocks", "[name of module]", commandName, "setTreeBlockSource")) {    	 		TreeBlockSource temp =   (TreeBlockSource)replaceEmployee(commandRec, TreeBlockSource.class, arguments, "Source of trees for chart", treeBlockSourceTask); 			if (temp!=null) { 				treeBlockSourceTask = temp;				treeBlockSourceTask.setHiringCommand(tbstC);				treeBlockSourceName.setValue(treeBlockSourceTask.getName());	   	 		if (taxa!=null)	   	 			treeBlockSourceTask.setPreferredTaxa(taxa);				if (cWindow!=null){					chartWindowTask.setTaxa(taxa);					chartWindowTask.setNumberTask(numberTask);					chartWindowTask.setItemsSource(treeBlockSourceTask);					if (!commandRec.scripting()){						chartWindowTask.doCounts(commandRec);	 				} 				}	 			return treeBlockSourceTask; 			}    	 		//treeValues -- resize matrix according to number of trees, or if infinite, to chosen number    	 	}   	 	else if (checker.compare(getClass(), "Sets the module calculating numbers for the tree blocks", "[name of module]", commandName, "setCalculator")) {    	 		NumberForTreeBlock temp =  (NumberForTreeBlock)replaceEmployee(commandRec, NumberForTreeBlock.class, arguments, "Value to calculate for trees", numberTask);    	 		//((TreeBlockValuesWindow)window).setNumberTask(numberTask); 			if (temp!=null) { 				numberTask = temp;				numberTask.setHiringCommand(ntC);				numberTaskName.setValue(numberTask.getName());				if (cWindow!=null){					chartWindowTask.setTaxa(taxa);					chartWindowTask.setNumberTask(numberTask);					chartWindowTask.setItemsSource(treeBlockSourceTask);					if (!commandRec.scripting()){						chartWindowTask.doCounts(commandRec);	 				} 				}	 			return numberTask;	 		}    	 	}    	 	else    	 		return super.doCommand(commandName, arguments, commandRec, checker);    	 	return null;    	 }	/*.................................................................................................................*/	public void employeeParametersChanged(MesquiteModule employee, MesquiteModule source, Notification notification, CommandRecord commandRec) {			if (cWindow!=null && !commandRec.scripting() && live.getValue())				if (Notification.getCode(notification) == MesquiteListener.PARTS_ADDED && source == treeBlockSourceTask) {					int[] notifParam = notification.getParameters();					if (notifParam!=null && notifParam.length==2)						chartWindowTask.doCounts(commandRec, notifParam[0]+1, notifParam[0]+notifParam[1], false);					else						chartWindowTask.doCounts(commandRec);				}				else chartWindowTask.doCounts(commandRec); //perhaps put intemschart in charge liveness and pass here whether through ePC	}	/*.................................................................................................................*/    	 public String getName() {		return "Tree Block Values";   	 }	/*.................................................................................................................*/    	 public String getNameForMenuItem() {		return "Tree Blocks";   	 }	/*.................................................................................................................*/  	 public String getExplanation() {		return "Makes a chart showing some value for each of many tree blocks.";   	 }   	 }