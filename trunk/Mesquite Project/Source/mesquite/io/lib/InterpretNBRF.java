/* Mesquite (package mesquite.io).  Copyright 2000-2006 D. Maddison and W. Maddison. Version 1.11, June 2006.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.io.lib;/*~~  */import java.util.*;import java.awt.*;import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.lib.duties.*;import mesquite.categ.lib.*;/* - make more things protected- document classes, add doc target to project (use mesquite.lib for example).- scriptinggeneral exporting stuff:	- tokenize properlygeneral importing stuff	- set punctuation on read in, including single quotes	- set reading to not consider [ ], or to consider different comment tokens*//* Post version 1.0:	- verbose file reading	- on reading, store comment line as footnote*//* ============  a file interpreter for NBRF files ============*//** This is the class for interpreting NBRF/PIR files.  It is subclassed to make interpreters specifically forDNA and Protein files. */public abstract class InterpretNBRF extends FileInterpreterI {	Class[] acceptedClasses;/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) {		acceptedClasses = new Class[] {ProteinState.class, DNAState.class}; 		return true;  //make this depend on taxa reader being found?)  	 }  	 /*.................................................................................................................*/	public boolean canExportEver() {  		 return true;  //	}/*.................................................................................................................*/	public boolean canExportProject(MesquiteProject project) {  		 return project.getNumberCharMatrices(acceptedClasses) > 0;  //	}	/*.................................................................................................................*/	public boolean canExportData(Class dataClass) {  		 for (int i = 0; i<acceptedClasses.length; i++)		 	if (dataClass==acceptedClasses[i])		 		return true;		 return false; 	}/*.................................................................................................................*/	public boolean canImport() {  		 return true;	}	/** Returns whether the module can read (import) files considering the passed argument string (e.g., fuse) */	public boolean canImport(String arguments){		return true;	}/*.................................................................................................................*/	public abstract void setNBRFState(CharacterData data, int ic, int it, char c);/*.................................................................................................................*/	public abstract CharacterData createData(CharactersManager charTask, Taxa taxa);/*.................................................................................................................*/	public void readFile(MesquiteProject mf, MesquiteFile file, String arguments, CommandRecord commandRec) {		incrementMenuResetSuppression();		ProgressIndicator progIndicator = new ProgressIndicator(mf,"Importing File "+ file.getName(), file.existingLength());		progIndicator.start();		file.linkProgressIndicator(progIndicator);		boolean fuse = false;		String fRA = parser.getFirstToken(arguments);		while (!StringUtil.blank(fRA)) {			if (fRA.equalsIgnoreCase(StringUtil.argumentMarker + "fuseTaxaCharBlocks"))				fuse = true;			fRA = parser.getNextToken();		}		if (file.openReading()) {			TaxaManager taxaTask = (TaxaManager)findElementManager(Taxa.class);			 CharactersManager charTask = (CharactersManager)findElementManager(CharacterData.class);			Taxa taxa = null;			if (fuse)				taxa = getProject().chooseTaxa(containerOfModule(), "To which block of taxa do you want to fuse the taxa from the file \"" + file.getName() + " being read in?  If you choose cancel, a new taxa block will be created instead.", true,commandRec);			if (taxa == null){				taxa = taxaTask.makeNewTaxa(getProject().getTaxas().getUniqueName("Untitled Block of Taxa"), 0, false);				taxa.addToFile(file, getProject(), taxaTask);			}			CategoricalData data = null;			if (fuse)				data = (CategoricalData)getProject().chooseData(containerOfModule(), null, taxa, CategoricalState.class, "Select matrix with which to fuse the matrix from the file \"" + file.getName() + "  being read.   If you choose cancel, a new matrix will be created instead.",  true, commandRec);			if (data == null){				data =(CategoricalData)createData(charTask,taxa);				data.addToFile(file, getProject(), null);			}			boolean wassave = data.saveChangeHistory;			data.saveChangeHistory = false;						//file.readLine() reads to next CR, LF, or CRLF.  To read to a target string, use file.readLine(targetString);						/*Each module is automatically given its own Parser object named "parser".  				To change its punctuation, use, e.g.:					parser.setPunctuationString(".,:;*&^%$#@!");				To change its whitespace, use, e.g.					parser.setWhitespaceString(" \t\r\n");				To return to default NEXUS punctuation & whitespace, use					parser.setPunctuationString(null);					parser.setWhitespaceString(null);			*/					//	int numChars=0;			int numTaxa = 0;			if (fuse)				numTaxa = taxa.getNumTaxa();			StringBuffer sb = new StringBuffer(1000);			file.readLine(sb);			String line = sb.toString();						boolean abort = false;			while (!StringUtil.blank(line) && !abort) {				parser.setString(line); //sets the string to be used by the parser to "line" and sets the pos to 0								parser.setPunctuationString(";");				String token = parser.getFirstToken(line); //should be >DL				char c = parser.nextDarkChar();      // should be a semicolon				parser.setPunctuationString(null);								token = parser.getNextToken();  //taxon Name				taxa.addTaxa(numTaxa-1, 1, true);				Taxon t = taxa.getTaxon(numTaxa);								if (t!=null) {					t.setName(token);					progIndicator.setText("Reading taxon: "+token);					file.readLine(sb);  // skip over comment line					line = sb.toString();					if (line==null) break;					line = file.readLine("*");  // pull in sequence					if (line==null) break;					parser.setString(line); 					int ic = 0;					while (parser.getPosition()<line.length()) {						c=parser.nextDarkChar();						if (c!='*' && c!= '\0') {							if (data.getNumChars() <= ic) {								data.addCharacters(data.getNumChars()-1, 1, false);   // add a character if needed							}							//data.setState(ic, numTaxa, c);    // setting state to that specified by character c							setNBRFState(data,ic, numTaxa, c);    // setting state to that specified by character c						}						ic += 1;					}				}				numTaxa++;	//			file.readLine(sb);				line = file.readNextDarkLine();		// added 1.01				if (file.getFileAborted()) {					abort = true;				}			}			data.saveChangeHistory = wassave;			data.resetChangedSinceSave();			finishImport(progIndicator, file, abort, commandRec);		}		decrementMenuResetSuppression();	}	/* ============================  exporting ============================*/	boolean includeGaps = false;	boolean excludeEmpty = false;	/*.................................................................................................................*/		public boolean getExportOptions(boolean dataSelected, boolean taxaSelected){		MesquiteInteger buttonPressed = new MesquiteInteger(1);		ExporterDialog exportDialog = new ExporterDialog(this,containerOfModule(), "Export NBRF Options", buttonPressed);				Checkbox includeGapsCheckBox = exportDialog.addCheckBox("include gaps", includeGaps);		Checkbox excludeEmptyCheckBox = exportDialog.addCheckBox("exclude empty sequences", suppressAllGapTaxa);				exportDialog.completeAndShowDialog(dataSelected, taxaSelected);					boolean ok = (exportDialog.query(dataSelected, taxaSelected)==0);		suppressAllGapTaxa = excludeEmptyCheckBox.getState();				includeGaps = includeGapsCheckBox.getState();		exportDialog.dispose();		return ok;	}		/*.................................................................................................................*/		public boolean getExportOptionsSimple(boolean dataSelected, boolean taxaSelected){   // an example of a simple query, that only proved line delimiter choice; not used here		return (ExporterDialog.query(this,containerOfModule(), "Export NBRF Options")==0);	}			/*.................................................................................................................*/	public String getLineStart(){  		return ">DL; ";	}		/*.................................................................................................................*/	public abstract CharacterData findDataToExport(MesquiteFile file, String arguments, CommandRecord commandRec);	/*.................................................................................................................*/	public StringBuffer getFile(CharacterData data) {		if (data==null)			return null;		Taxa taxa = data.getTaxa();		int numTaxa = taxa.getNumTaxa();		int numChars = data.getNumChars();		StringBuffer outputBuffer = new StringBuffer(numTaxa*(20 + numChars));				int counter = 1;		for (int it = 0; it<numTaxa; it++){			if ((!writeOnlySelectedTaxa || (taxa.getSelected(it))) && (!suppressAllGapTaxa || !taxonEntirelyInapplicable(it, data))) {			// TO DO: also have the option of only writing taxa with data in them				counter = 1;				outputBuffer.append(getLineStart());				outputBuffer.append(ParseUtil.tokenize(taxa.getTaxonName(it)) + getLineEnding());				outputBuffer.append(ParseUtil.tokenize(taxa.getTaxonName(it)) + getLineEnding());				for (int ic = 0; ic<numChars; ic++) {					if (!writeOnlySelectedData || (data.getSelected(ic))){						int currentSize = outputBuffer.length();						if (includeGaps || (!data.isInapplicable(ic,it))) {								data.statesIntoStringBuffer(ic, it, outputBuffer, false);								counter ++;						}						if (outputBuffer.length()-currentSize>1) {							alert("Sorry, this data matrix can't be exported to this format (some character states aren't represented by a single symbol [char. " + CharacterStates.toExternal(ic) + ", taxon " + Taxon.toExternal(it) + "])");							return null;						}						if ((counter % 50 == 1) && (counter > 1)) {    // modulo							outputBuffer.append(getLineEnding());						}					}				}				outputBuffer.append("*" + getLineEnding());			}		}		return outputBuffer;	}   	/** returns whether the character ic is entirely inapplicable codings*/   	public boolean taxonEntirelyInapplicable(int it, CharacterData data){   		for (int ic = 0; ic< data.getNumChars(); ic++)   			if (!data.isInapplicable(ic, it))   				return false;   		return true;   	}	boolean suppressAllGapTaxa = false;	/*.................................................................................................................*/	public void exportFile(MesquiteFile file, String arguments, CommandRecord commandRec) { //if file is null, consider whole project open to export		Arguments args = new Arguments(new Parser(arguments), true);		boolean usePrevious = args.parameterExists("usePrevious");		suppressAllGapTaxa = args.parameterExists("suppressAllGapTaxa");		CharacterData data = findDataToExport(file, arguments, commandRec);		if (data ==null) {			showLogWindow(true);			logln("WARNING: No suitable data available for export to a file of format \"" + getName() + "\".  The file will not be written.\n");			return;		}		Taxa taxa = data.getTaxa();		writeOnlySelectedData = args.parameterExists("writeOnlySelectedData");		if (!commandRec.scripting() && !usePrevious)			if (!getExportOptions(data.anySelected(), taxa.anySelected()))				return;				StringBuffer outputBuffer = getFile(data);				saveExportedFileWithExtension(outputBuffer, arguments, "nbrf", commandRec);		suppressAllGapTaxa = false;	}	/*.................................................................................................................*/    	 public String getName() {		return "NBRF/PIR file";   	 }	/*.................................................................................................................*/   	  	/** returns an explanation of what the module does.*/ 	public String getExplanation() { 		return "Imports and exports NBRF files that consist of molecular sequence data." ;   	 }	/*.................................................................................................................*/   	    	 }	