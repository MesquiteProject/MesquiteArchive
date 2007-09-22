/* Mesquite source code.  Copyright 1997-2007 W. Maddison and D. Maddison.Version 2.0, September 2007.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html) */package mesquite.align.ClustalAlign;/*~~  */import java.util.*;import java.lang.*;import java.io.*;import java.awt.*;import java.awt.event.*;import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.lib.duties.*;import mesquite.categ.lib.*;import mesquite.lib.table.*;import mesquite.align.lib.*;/* ======================================================================== */public class ClustalAlign extends MultipleSequenceAligner implements ActionListener{	String clustalPath;	SingleLineTextField clustalPathField =  null;	boolean preferencesSet = false;	String clustalOptions = " -gapopen=15.0  -gapext=6.66 " ;	Random rng;	public static int runs = 0;	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, boolean hiredByName) {		loadPreferences();		rng = new Random(System.currentTimeMillis());		return true;	}	/*.................................................................................................................*/	/** returns whether this module is requesting to appear as a primary choice */	public boolean requestPrimaryChoice(){		return true;  	}	/*.................................................................................................................*/	public void processSingleXMLPreference (String tag, String content) {		if ("clustalPath".equalsIgnoreCase(tag)) {			clustalPath = StringUtil.cleanXMLEscapeCharacters(content);			if (!StringUtil.blank(clustalPath)) {				if (!clustalPath.endsWith(MesquiteFile.fileSeparator))					clustalPath+=MesquiteFile.fileSeparator;			}		}		else if ("clustalOptions".equalsIgnoreCase(tag))			clustalOptions = StringUtil.cleanXMLEscapeCharacters(content);		preferencesSet = true;	}	/*.................................................................................................................*/	public String preparePreferencesForXML () {		StringBuffer buffer = new StringBuffer(200);		StringUtil.appendXMLTag(buffer, 2, "clustalPath", clustalPath);  		StringUtil.appendXMLTag(buffer, 2, "clustalOptions", clustalOptions);  		preferencesSet = true;		return buffer.toString();	}	/*.................................................................................................................*/	public boolean queryOptions() {		MesquiteInteger buttonPressed = new MesquiteInteger(1);		ExtensibleDialog queryFilesDialog = new ExtensibleDialog(containerOfModule(), "ClustalW Locations & Options",buttonPressed);  //MesquiteTrunk.mesquiteTrunk.containerOfModule()		queryFilesDialog.addLabel("ClustalW - File Locations & Options");		clustalPathField = queryFilesDialog.addTextField("Directory containing ClustalW:", clustalPath, 40);		Button clustalBrowseButton = queryFilesDialog.addAListenedButton("Browse...",null, this);		clustalBrowseButton.setActionCommand("clustalBrowse");		SingleLineTextField clustalOptionsField = queryFilesDialog.addTextField("Clustal options:", clustalOptions, 26, true);		queryFilesDialog.completeAndShowDialog(true);		if (buttonPressed.getValue()==0)  {			clustalPath = clustalPathField.getText();			clustalOptions = clustalOptionsField.getText();			storePreferences();		}		queryFilesDialog.dispose();		return (buttonPressed.getValue()==0);	}	/*.................................................................................................................*/	private void saveNBRFFile(MolecularData data, String directoryPath, String fileName, boolean[] taxaToAlign, int firstSite, int lastSite) {		if (data==null)			return;				runs++;		String path = createSupportDirectory() + MesquiteFile.fileSeparator + fileName;  // place files in support directory for module		incrementMenuResetSuppression();		Taxa taxa = data.getTaxa();		FileCoordinator coord = getFileCoordinator();		MesquiteFile tempDataFile = (MesquiteFile)coord.doCommand("newLinkedFile", StringUtil.tokenize(path), CommandChecker.defaultChecker); //TODO: never scripting???		TaxaManager taxaManager = (TaxaManager)findElementManager(Taxa.class);		Taxa newTaxa =taxa.cloneTaxa(taxaToAlign);		newTaxa.addToFile(tempDataFile, null, taxaManager);		//rename taxa so clustal doesn't screw around with names		for (int it=0; it<newTaxa.getNumTaxa(); it++)			newTaxa.setTaxonName(it, "t" + it);		CharMatrixManager matrixManager = data.getMatrixManager();		int numNewChars=0;		int firstChar = -1;		for (int ic=0; ic<data.getNumChars(); ic++){			if (data.getSelected(ic) || (firstSite>=0 && MesquiteInteger.isCombinable(firstSite) && ic>= firstSite && lastSite<data.getNumChars() && MesquiteInteger.isCombinable(lastSite) && ic<= lastSite)){				numNewChars++;				if (firstChar<0) //first one found					firstChar=ic;			}		}		MolecularData newData = (MolecularData)matrixManager.getNewData(newTaxa, numNewChars);		for (int ic=0; ic<newData.getNumChars(); ic++){			if (taxaToAlign!=null) {				int count=0;				for (int i  =  0; i<taxaToAlign.length; i++) 					if (taxaToAlign[i]) {						newData.setState(ic, count, data.getState(ic+firstChar, i));						count++;					}			}			else for (int it=0; it<newTaxa.getNumTaxa(); it++)				newData.setState(ic, it, data.getState(ic+firstChar, it));		}		//newData = data.cloneData(); ???				newData.setName(data.getName());		newData.addToFile(tempDataFile, getProject(), null);		FileInterpreterI exporter=null;		if (data instanceof DNAData)			exporter = (FileInterpreterI)coord.findEmployeeWithName("#InterpretNBRFDNA");		else if (data instanceof ProteinData)			exporter = (FileInterpreterI)coord.findEmployeeWithName("#InterpretNBRFProtein");		if (exporter!=null) {			String ext = exporter.preferredDataFileExtension();			if (StringUtil.blank(ext))				ext = "";			else				ext = "." + ext;			String s = "file = " + StringUtil.tokenize(fileName + ext) + " directory = " + StringUtil.tokenize(directoryPath) + " noTrees suppressAllGapTaxa";			//if (data.anySelected()) 			//	s += " writeOnlySelectedData";			s+= " usePrevious";			coord.export(exporter, tempDataFile, s);		}		tempDataFile.close();		decrementMenuResetSuppression();	}	/*.................................................................................................................*/	public long[][] alignSequences(MCategoricalDistribution matrix, boolean[] taxaToAlign, int firstSite, int lastSite, int firstTaxon, int lastTaxon) {		if (!queryOptions())			return null;		if (!(matrix.getParentData() != null && matrix.getParentData() instanceof MolecularData)){			discreetAlert( "Sorry, Clustal Aligner works only if given a full MolecularData object");			return null;		}		MolecularData data = (MolecularData)matrix.getParentData();		boolean isProtein = data instanceof ProteinData;		boolean pleaseStorePref = false;		if (!preferencesSet || StringUtil.blank(clustalPath)) {			clustalPath = MesquiteFile.chooseDirectory("Choose directory containing clustalW: ");			if (StringUtil.blank(clustalPath))				return null;			if (!clustalPath.endsWith(MesquiteFile.fileSeparator))				clustalPath+=MesquiteFile.fileSeparator;			pleaseStorePref = true;		}		mesquiteTrunk.incrementProjectBrowserRefreshSuppression();		if (pleaseStorePref)			storePreferences();		data.setEditorInhibition(true);		String unique = MesquiteTrunk.getUniqueIDBase() + Math.abs(rng.nextInt());		String rootDir = createSupportDirectory() + MesquiteFile.fileSeparator;  //replace this with current directory of file//		StringBuffer fileBuffer = getFileInBuffer(data);		String fileName = "tempAlign" + MesquiteFile.massageStringToFilePathSafe(unique) + ".nbrf";   //replace this with actual file name?		String filePath = rootDir +  fileName;		if (taxaToAlign!=null)			saveNBRFFile(data, rootDir, fileName, taxaToAlign, firstSite, lastSite);		else if (!(firstTaxon==0 && lastTaxon==matrix.getNumTaxa())) {  // we are doing something other than all taxa.			boolean[] taxaToAlignLocal = new boolean[matrix.getNumTaxa()];			for (int it = 0; it<matrix.getNumTaxa(); it++)				taxaToAlignLocal[it] =  (it>=firstTaxon && it<= lastTaxon);			saveNBRFFile(data, rootDir, fileName, taxaToAlignLocal, firstSite, lastSite);		}		else			saveNBRFFile(data, rootDir, fileName, null, -1, -1);		String runningFilePath = rootDir + "running" + MesquiteFile.massageStringToFilePathSafe(unique);		String outFilePath = rootDir + "alignedFile" + MesquiteFile.massageStringToFilePathSafe(unique) + ".nbrf";//		MesquiteFile.putFileContents(filePath, fileBuffer.toString(), true);		StringBuffer shellScript = new StringBuffer(1000);		if (!MesquiteTrunk.isWindows())			shellScript.append(getClustalCommand() + "  -infile=" + StringUtil.protectForUnix(filePath) + " -outfile=" + StringUtil.protectForUnix(outFilePath) + " -align -output=pir ");		else			shellScript.append(getClustalCommand() + " \\ -infile=" + StringUtil.protectForWindows(filePath) + " -outfile=" + StringUtil.protectForWindows(outFilePath) + " -align -output=pir ");		if (isProtein)			shellScript.append("-type=protein ");		else			shellScript.append("-type=dna ");		shellScript.append(clustalOptions + StringUtil.lineEnding());//		shellScript.append(ShellScriptUtil.getRemoveCommand(runningFilePath));		String scriptPath = rootDir + "clustscript" + MesquiteFile.massageStringToFilePathSafe(unique) + ".bat";		MesquiteFile.putFileContents(scriptPath, shellScript.toString(), true);		boolean success = ShellScriptUtil.executeAndWaitForShell(scriptPath, runningFilePath, null, true, getName());		if (success){			FileCoordinator coord = getFileCoordinator();			MesquiteFile tempDataFile = null;			CommandRecord oldCR = MesquiteThread.getCurrentCommandRecord();			CommandRecord scr = new CommandRecord(true);			MesquiteThread.setCurrentCommandRecord(scr);			if (data instanceof DNAData)				tempDataFile = (MesquiteFile)coord.doCommand("linkFile", StringUtil.tokenize(outFilePath) + " " + StringUtil.tokenize("#InterpretNBRFDNA") + " suppressImportFileSave ", CommandChecker.defaultChecker); //TODO: never scripting???			else				tempDataFile = (MesquiteFile)coord.doCommand("linkFile", StringUtil.tokenize(outFilePath) + " " + StringUtil.tokenize("#InterpretNBRFProtein") + " suppressImportFileSave ", CommandChecker.defaultChecker); //TODO: never scripting???			MesquiteThread.setCurrentCommandRecord(oldCR);			CharacterData alignedData = getProject().getCharacterMatrix(tempDataFile,  0);			long[][] aligned = null;			Taxa alignedTaxa =  alignedData.getTaxa();			Taxa originalTaxa =  data.getTaxa();			if (alignedData!=null) {				int numChars = alignedData.getNumChars();				//sorting to get taxon names in correct order				int[] keys = new int[alignedData.getNumTaxa()];				for (int it = 0; it<alignedData.getNumTaxa(); it++){					String name = alignedTaxa.getTaxonName(it);					keys[it] = MesquiteInteger.fromString(name.substring(1, name.length()));  //this is original taxon number					if (!MesquiteInteger.isCombinable(keys[it])) {						success=false;						break;					}				}				if (success) {					for (int i=1; i<alignedTaxa.getNumTaxa(); i++) {						for (int j= i-1; j>=0 && keys[j]>keys[j+1]; j--) {							alignedTaxa.swapParts(j, j+1);							int kj = keys[j];							keys[j] = keys[j+1];							keys[j+1] = kj;							//alignedData.swapTaxa(j, j+1);						}					}					alignedData.changed(this, alignedTaxa, new Notification(MesquiteListener.PARTS_MOVED));					if (alignedData instanceof MolecularData){						aligned = new long[alignedData.getNumChars()][originalTaxa.getNumTaxa()];						for (int ic = 0; ic<alignedData.getNumChars(); ic++)							for (int it = 0; it<alignedData.getNumTaxa(); it++){								//	String name = alignedTaxa.getTaxonName(it);								//	int iCurrent = MesquiteInteger.fromString(name.substring(1, name.length()));  //this is original taxon number								//	int iTaxon = IntegerArray.indexOf(keys, iCurrent);								//if (iTaxon>=0 && MesquiteInteger.isCombinable(iTaxon))								aligned[ic][keys[it]] = ((MolecularData)alignedData).getState(ic, it);							}					}				}			}			if (tempDataFile!=null)				tempDataFile.close();			mesquiteTrunk.decrementProjectBrowserRefreshSuppression();			if (runs == 1)				deleteSupportDirectory();			runs--;			data.setEditorInhibition(false);			if (success) 				return aligned;			return null;		}		if (runs == 1)			deleteSupportDirectory();		runs--;		mesquiteTrunk.decrementProjectBrowserRefreshSuppression();		data.setEditorInhibition(false);		return null;	}		/*.................................................................................................................*	boolean executeAndWaitForClustal(String scriptPath, String runningFilePath, String outFilePath){		try{			ShellScriptUtil.setScriptFileToBeExecutable(scriptPath);			MesquiteFile.putFileContents(runningFilePath, "Clustal is running...", true);			Process proc = ShellScriptUtil.executeScript(scriptPath);			if (proc==null) {				return false;			}			else				while (MesquiteFile.fileExists(runningFilePath)){					try {						Thread.sleep(200);					}					catch (InterruptedException e){						MesquiteMessage.notifyProgrammer("InterruptedException in ClustalAlign");						return false;					}				}		}		catch (IOException e){			MesquiteMessage.notifyProgrammer("IOException in ClustalAlign");//			warning			return false;		}		return true;	}	/*.................................................................................................................*/	public boolean recoverClustalResults(MolecularData data, String outFilePath){		mesquiteTrunk.incrementProjectBrowserRefreshSuppression();		//reading aligned file		FileCoordinator coord = getFileCoordinator();		MesquiteFile tempDataFile = null;		if (data instanceof DNAData)			tempDataFile = (MesquiteFile)coord.doCommand("linkFile", StringUtil.tokenize(outFilePath) + " " + StringUtil.tokenize("#InterpretNBRFDNA") + " suppressImportFileSave ", CommandChecker.defaultChecker); //TODO: never scripting???		else			tempDataFile = (MesquiteFile)coord.doCommand("linkFile", StringUtil.tokenize(outFilePath) + " " + StringUtil.tokenize("#InterpretNBRFProtein") + " suppressImportFileSave ", CommandChecker.defaultChecker); //TODO: never scripting???		CharacterData alignedData = getProject().getCharacterMatrix(tempDataFile,  0);		return true;	}		/*------------------------------------------------------*    	/*.................................................................................................................*/	String getClustalCommand(){		if (MesquiteTrunk.isWindows())			return StringUtil.protectForWindows(clustalPath + "clustalw.exe");		else			return StringUtil.protectForUnix(clustalPath + "clustalw");	}	/*.................................................................................................................*/	public boolean showCitation() {		return true;	}	/*.................................................................................................................*/	public boolean isSubstantive(){		return true;	}	/*.................................................................................................................*/	public boolean isPrerelease(){		return false;	}	/*.................................................................................................................*/	public String getName() {		return "Clustal Align";	}	/*.................................................................................................................*/	public String getNameForMenuItem() {		return "Clustal Align...";	}	/*.................................................................................................................*/	/** returns an explanation of what the module does.*/	public String getExplanation() {		return "Sends the selected sequence to Clustal to align." ;	}	/*.................................................................................................................*/	/** returns the version number at which this module was first released.  If 0, then no version number is claimed.  If a POSITIVE integer	 * then the number refers to the Mesquite version.  This should be used only by modules part of the core release of Mesquite.	 * If a NEGATIVE integer, then the number refers to the local version of the package, e.g. a third party package*/	public int getVersionOfFirstRelease(){		return -100;  	}	/*.................................................................................................................*/	public  void actionPerformed(ActionEvent e) {		if (e.getActionCommand().equalsIgnoreCase("clustalBrowse")) {			clustalPath = MesquiteFile.chooseDirectory("Choose directory containing clustal: ");			if (!StringUtil.blank(clustalPath)) {				if (!clustalPath.endsWith(MesquiteFile.fileSeparator))					clustalPath+=MesquiteFile.fileSeparator;				clustalPathField.setText(clustalPath);			}		}	}}