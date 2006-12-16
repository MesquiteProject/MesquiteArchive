/* Mesquite source code.  Copyright 1997-2004 W. Maddison and D. Maddison. Version 1.05, September 2004.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.paup.PAUPNJ;/*~~  */import java.io.IOException;import java.util.*;import java.awt.*;import java.awt.event.*;import mesquite.categ.lib.DNAData;import mesquite.categ.lib.MCategoricalDistribution;import mesquite.categ.lib.MolecularData;import mesquite.categ.lib.ProteinData;import mesquite.lib.*;import mesquite.lib.characters.CharacterData;import mesquite.lib.characters.MCharactersDistribution;import mesquite.lib.duties.*;/** Supplies trees from tree blocks in a file.  Reads trees only when needed; hence suitable for files with too many trees to be held in memory at once, but slower than StoredTrees.*/public class PAUPNJ extends TreeSource implements MesquiteListener, ActionListener {	private MatrixSourceCoord matrixSourceTask;	int bestTree = -1;	String scoresPath = null;	boolean loaded = false;	String treename = null;	Random rng;	boolean preferencesSet = false;	String paupPath;	String paupOptions = " nj " ;	SingleLineTextField paupPathField =  null;	Taxa taxa;	protected MCharactersDistribution observedStates;	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) {		matrixSourceTask = (MatrixSourceCoord)hireCompatibleEmployee(commandRec, MatrixSourceCoord.class, getCharacterClass(), "Source of matrix (for " + getName() + ")");		if (matrixSourceTask == null)			return sorry(commandRec, getName() + " couldn't start because no source of matrix (for " + getName() + ") was obtained");		loadPreferences();		rng = new Random(System.currentTimeMillis());				return true;  	 }  	/** Called to provoke any necessary initialization.  This helps prevent the module's intialization queries to the user from	 happening at inopportune times (e.g., while a long chart calculation is in mid-progress)*/	public void initialize(Taxa taxa, CommandRecord commandRec){		this.taxa = taxa;		if (matrixSourceTask!=null) {			matrixSourceTask.initialize(taxa, commandRec);			if (observedStates ==null)				observedStates = matrixSourceTask.getCurrentMatrix(taxa, commandRec);		}	}	/*.................................................................................................................*/	public Class getCharacterClass() {		return null;	}	/*.................................................................................................................*/	public void processSingleXMLPreference (String tag, String content) {		if ("paupPath".equalsIgnoreCase(tag)) {			paupPath = StringUtil.cleanXMLEscapeCharacters(content);			if (!StringUtil.blank(paupPath)) {				if (!paupPath.endsWith(MesquiteFile.fileSeparator))					paupPath+=MesquiteFile.fileSeparator;			}		}		else if ("paupOptions".equalsIgnoreCase(tag))			paupOptions = StringUtil.cleanXMLEscapeCharacters(content);		preferencesSet = true;	}	/*.................................................................................................................*/	public String preparePreferencesForXML () {		StringBuffer buffer = new StringBuffer(200);		StringUtil.appendXMLTag(buffer, 2, "paupPath", paupPath);  		StringUtil.appendXMLTag(buffer, 2, "paupOptions", paupOptions);  		preferencesSet = true;		return buffer.toString();	}	/*.................................................................................................................*/	public boolean queryOptions() {		MesquiteInteger buttonPressed = new MesquiteInteger(1);		ExtensibleDialog queryFilesDialog = new ExtensibleDialog(containerOfModule(), "PAUP* Locations & Options",buttonPressed);  //MesquiteTrunk.mesquiteTrunk.containerOfModule()		queryFilesDialog.addLabel("PAUP* - File Locations & Options");		paupPathField = queryFilesDialog.addTextField("Directory containing PAUP*:", paupPath, 40);		Button paupBrowseButton = queryFilesDialog.addAListenedButton("Browse...",null, this);		paupBrowseButton.setActionCommand("paupBrowse");		SingleLineTextField paupOptionsField = queryFilesDialog.addTextField("PAUP* options:", paupOptions, 26, true);		queryFilesDialog.completeAndShowDialog(true);		if (buttonPressed.getValue()==0)  {			paupPath = paupPathField.getText();			paupOptions = paupOptionsField.getText();			storePreferences();		}		queryFilesDialog.dispose();		return (buttonPressed.getValue()==0);	}	/*.................................................................................................................*/   	 public boolean isSubstantive(){   	 	return true;   	 }	/*.................................................................................................................*/   	 public boolean isPrerelease(){   	 	return true;   	 }   	 	/*.................................................................................................................*/  	public void setPreferredTaxa(Taxa taxa) {  	}  	  		/*.................................................................................................................*/   	public Tree getCurrentTree(Taxa taxa, CommandRecord commandRec) {   		if (observedStates ==null)   			initialize(taxa, commandRec);   		Tree t = doPAUP(observedStates,commandRec);   		return t;    	}	/*.................................................................................................................*/   	public Tree getTree(Taxa taxa, int itree, CommandRecord commandRec) {  		return getCurrentTree(taxa, commandRec);   	}	/*.................................................................................................................*/   	public int getNumberOfTrees(Taxa taxa, CommandRecord commandRec) {		return 1;    	}   	/*.................................................................................................................*/   	public String getTreeNameString(Taxa taxa, int itree, CommandRecord commandRec) {   		return treename;   	}	/*.................................................................................................................*/   	public String getCurrentTreeNameString(Taxa taxa, CommandRecord commandRec) {   		return treename;   	}   	   	   	/*.................................................................................................................*/   	public void writeNEXUSFile(Taxa taxa, String dir, String fileName, String path, MolecularData data, CommandRecord commandRec) {   		if (path != null) {   			MesquiteFile f = MesquiteFile.newFile(dir, fileName);   			f.openWriting(true);   			f.writeLine("#NEXUS" + StringUtil.lineEnding());   			f.writeLine(((TaxaManager)findElementManager(Taxa.class)).getTaxaBlock(data.getTaxa(), null));   			data.getMatrixManager().writeCharactersBlock(data, null, f, null);   			f.closeWriting();   		}   	}		/*.................................................................................................................*/	public Tree doPAUP(MCharactersDistribution matrix, CommandRecord commandRec) {		if (!queryOptions())			return null;		if (matrix==null)			return null;		if (!(matrix.getParentData() != null && matrix.getParentData() instanceof MolecularData)){			discreetAlert(commandRec, "Sorry, paup Aligner works only if given a full MolecularData object");			return null;		}		MolecularData data = (MolecularData)matrix.getParentData();		boolean isProtein = data instanceof ProteinData;			boolean pleaseStorePref = false;		if (!preferencesSet || StringUtil.blank(paupPath)) {			paupPath = MesquiteFile.chooseDirectory("Choose directory containing PAUP*: ");			if (StringUtil.blank(paupPath))				return null;			if (!paupPath.endsWith(MesquiteFile.fileSeparator))				paupPath+=MesquiteFile.fileSeparator;			pleaseStorePref = true;		}		mesquiteTrunk.incrementProjectBrowserRefreshSuppression();		if (pleaseStorePref)			storePreferences();		data.setEditorInhibition(true);		String unique = MesquiteTrunk.getUniqueIDBase() + Math.abs(rng.nextInt());		String rootDir = createSupportDirectory() + MesquiteFile.fileSeparator;  //replace this with current directory of file		String fileName = "tempData" + MesquiteFile.massageStringToFilePathSafe(unique) + ".nex";   //replace this with actual file name?		String filePath = rootDir +  fileName;	  	writeNEXUSFile(taxa,  rootDir,  fileName,  filePath,  data,  commandRec);	  		 		String runningFilePath = rootDir + "running" + MesquiteFile.massageStringToFilePathSafe(unique);		String outFilePath = rootDir + "tempTree" + MesquiteFile.massageStringToFilePathSafe(unique) + ".tre";		StringBuffer shellScript = new StringBuffer(1000);				String commandFileName =  "paupCmds.txt";		String treeFileName = "tree.txt";		String commandFilePath = rootDir + commandFileName;		String treeFilePath = rootDir + treeFileName;		StringBuffer sb = new StringBuffer();		sb.append("#NEXUS" + StringUtil.lineEnding() + StringUtil.lineEnding() + "BEGIN PAUP;" + StringUtil.lineEnding());		sb.append("\t" + "set torder=right tcompress outroot = monophyl nowarntree nowarntsave nowarnroot;" + StringUtil.lineEnding());		sb.append("\t" + "exec " + filePath + ";" + StringUtil.lineEnding());		if (!isProtein)			sb.append("\t" + "dset distance=hky85;" + StringUtil.lineEnding());		sb.append("\t" + "nj;" + StringUtil.lineEnding());		sb.append("\t" + "savetree brlens=yes replace=yes file='" + rootDir + treeFileName + "' ;" + StringUtil.lineEnding());		sb.append("\t" + "quit;" + StringUtil.lineEnding());		sb.append("ENDBLOCK;" + StringUtil.lineEnding());		MesquiteFile.putFileContents(commandFilePath, sb.toString(), true);				shellScript.append(getPAUPCommand()+ " " + commandFilePath + StringUtil.lineEnding());			shellScript.append(ShellScriptUtil.getRemoveCommand(runningFilePath));		String scriptPath = rootDir + "paupScript" + MesquiteFile.massageStringToFilePathSafe(unique) + ".bat";		MesquiteFile.putFileContents(scriptPath, shellScript.toString(), true);		boolean success = executeAndWaitForPAUP(scriptPath, runningFilePath, outFilePath);		if (success){			FileCoordinator coord = getFileCoordinator();			MesquiteFile tempDataFile = null;			tempDataFile = (MesquiteFile)coord.doCommand("includeTreeFile", StringUtil.tokenize(treeFilePath) + " " + StringUtil.tokenize("#InterpretNEXUS") + " suppressImportFileSave ", CommandRecord.scriptingRecord, CommandChecker.defaultChecker); //TODO: never scripting???			TreesManager manager = (TreesManager)findElementManager(TreeVector.class);			Tree t =null;			int numTB = manager.getNumberTreeBlocks(taxa);			TreeVector tv = manager.getTreeBlock(taxa,numTB-1);			if (tv!=null) {				t = tv.getTree(0);				if (t!=null)					success=true;			}						mesquiteTrunk.decrementProjectBrowserRefreshSuppression();			if (tempDataFile!=null)				tempDataFile.close();			deleteSupportDirectory();			data.setEditorInhibition(false);			if (success) 				return t;			return null;		}		deleteSupportDirectory();		mesquiteTrunk.decrementProjectBrowserRefreshSuppression();		data.setEditorInhibition(false);		return null;	}					/*.................................................................................................................*/	boolean executeAndWaitForPAUP(String scriptPath, String runningFilePath, String outFilePath){		try{			ShellScriptUtil.setScriptFileToBeExecutable(scriptPath);			MesquiteFile.putFileContents(runningFilePath, "PAUP is running...", true);			Process proc = ShellScriptUtil.executeScript(scriptPath);			if (proc==null) {				return false;			}			else				while (MesquiteFile.fileExists(runningFilePath)){					try {						Thread.sleep(200);					}					catch (InterruptedException e){						return false;					}				}		}		catch (IOException e){//			warning			return false;		}		return true;	}	/*.................................................................................................................*	public boolean recoverPAUPResults(MolecularData data, String outFilePath, CommandRecord commandRec){		mesquiteTrunk.incrementProjectBrowserRefreshSuppression();			return true;	}		/*------------------------------------------------------*    	/*.................................................................................................................*/	String getPAUPCommand(){		if (MesquiteTrunk.isWindows())			return StringUtil.protectForWindows(paupPath + "paup.exe");		else			return StringUtil.protectForUnix(paupPath + "paup");	}				/*.................................................................................................................*/ 	/** returns the version number at which this module was first released.  If 0, then no version number is claimed.  If a POSITIVE integer 	 * then the number refers to the Mesquite version.  This should be used only by modules part of the core release of Mesquite. 	 * If a NEGATIVE integer, then the number refers to the local version of the package, e.g. a third party package*/    	public int getVersionOfFirstRelease(){    		return 110;      	}	/*.................................................................................................................*/    	 public String getName() {		return "PAUP NJ";   	 }	/*.................................................................................................................*/	/** returns whether this module is requesting to appear as a primary choice */   	public boolean requestPrimaryChoice(){   		return true;     	}	/*.................................................................................................................*/  	 public String getExplanation() {		return "Supplies tree from PAUP's NJ";   	 }	/*.................................................................................................................*/   	public String getParameters() {		return "Trees obtained from file " + scoresPath;   	}	/*.................................................................................................................*/	public  void actionPerformed(ActionEvent e) {		if (e.getActionCommand().equalsIgnoreCase("paupBrowse")) {			paupPath = MesquiteFile.chooseDirectory("Choose directory containing paup: ");			if (!StringUtil.blank(paupPath)) {				if (!paupPath.endsWith(MesquiteFile.fileSeparator))					paupPath+=MesquiteFile.fileSeparator;				paupPathField.setText(paupPath);			}		}	}}