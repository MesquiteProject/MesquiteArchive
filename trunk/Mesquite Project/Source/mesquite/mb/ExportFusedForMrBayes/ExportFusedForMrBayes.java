/* Mesquite (package mesquite.io).  Copyright 2000-2005 D. Maddison and W. Maddison. Version 1.06, August 2005.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html) */package mesquite.mb.ExportFusedForMrBayes;/*~~  */import java.util.*;import java.awt.*;import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.lib.characters.CharacterData;import mesquite.lib.duties.*;import mesquite.assoc.lib.*;import mesquite.categ.lib.*;import mesquite.cont.lib.*;public class ExportFusedForMrBayes extends FileInterpreterI {	AssociationSource associationTask;	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) {		return true;  //make this depend on taxa reader being found?)	}	public boolean isPrerelease(){		return true;	}	public boolean isSubstantive(){		return true;	}	/*.................................................................................................................*/	public String preferredDataFileExtension() {		return "nex";	}	/*.................................................................................................................*/	public boolean canExportEver() {  		return true;  //	}	/*.................................................................................................................*/	public boolean canExportProject(MesquiteProject project) {  		return (project.getNumberCharMatrices(CategoricalState.class) > 0) ;	}	/*.................................................................................................................*/	public boolean canExportData(Class dataClass) {  		return CategoricalData.class.isAssignableFrom(dataClass);	}	/*.................................................................................................................*/	public boolean canImport() {  		return false;	}	/*.................................................................................................................*/	public void readFile(MesquiteProject mf, MesquiteFile file, String arguments, CommandRecord commandRec) {	}	/* ============================  exporting ============================*/	/*.................................................................................................................*/	boolean convertAmbiguities = false;	boolean useData = true;	String addendum = "";	String fileName = "untitled.nex";		public boolean getExportOptions(CharacterData data, boolean dataSelected, boolean taxaSelected){		MesquiteInteger buttonPressed = new MesquiteInteger(1);		ExporterDialog exportDialog = new ExporterDialog(this,containerOfModule(), "Export Fused NEXUS For MrBayes", buttonPressed);		exportDialog.setSuppressLineEndQuery(true);		exportDialog.setDefaultButton(null);//		Checkbox convertToMissing = exportDialog.addCheckBox("convert partial ambiguities to missing", convertAmbiguities);		exportDialog.addLabel("MrBayes block: ");				addendum = getMrBayesBlock(data);				TextArea fsText =exportDialog.addTextAreaSmallFont(addendum,16);						exportDialog.completeAndShowDialog(dataSelected, taxaSelected);					boolean ok = (exportDialog.query(dataSelected, taxaSelected)==0);		//		convertAmbiguities = convertToMissing.getState();		addendum = fsText.getText();		exportDialog.dispose();		return ok;	}			private String basicBlock(){		String sT = "begin mrbayes;\n\tset autoclose=yes nowarn=yes;";		sT +="\n\tlset nst=6 rates=invgamma;\n\tunlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all); \n\tprset applyto=(all) ratepr=variable;\n\tmcmcp ngen= 10000000 printfreq=1000  samplefreq=1000 nchains=4 savebrlens=yes;\n\tmcmc;\nend;";		return sT;	}	private String getMrBayesBlock(CharacterData data){		CharacterPartition characterPartition = (CharacterPartition)data.getCurrentSpecsSet(CharacterPartition.class);		if (characterPartition == null)			return basicBlock();		CharactersGroup[] parts = characterPartition.getGroups();		if (parts==null) {			return basicBlock();		}		else {			String sT = "begin mrbayes;\n\tset autoclose=yes nowarn=yes;  ";			int numCharSets = 0;			//charsets			for (int i=0; i<parts.length; i++) {				String q = ListableVector.getListOfMatches((Listable[])characterPartition.getProperties(), parts[i], CharacterStates.toExternal(0));				if (q != null) {					sT +=  "\n\tcharset " + StringUtil.tokenize(parts[i].getName()) + " = " + q + ";";					numCharSets++;				}			}			if (numCharSets <=1)				return "";			sT += "\n\tpartition currentPartition = " + numCharSets + ": ";			boolean firstTime = true;			String nums = "";			int num = 0;			for (int i=0; i<parts.length; i++) {				String q = ListableVector.getListOfMatches((Listable[])characterPartition.getProperties(), parts[i], CharacterStates.toExternal(0));				if (q != null) {					if (!firstTime){						sT += ", ";						nums += ", ";					}					firstTime = false;					sT += StringUtil.tokenize(parts[i].getName());					nums += Integer.toString(num+1);					num++;				}			}			sT +=";\n\tset partition = currentPartition;\n\tlset applyto=(" + nums + ") nst=6 rates=invgamma;\n\tunlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all); \n\tprset applyto=(all) ratepr=variable;\n\tmcmcp ngen= 10000000 printfreq=1000  samplefreq=1000 nchains=4 savebrlens=yes;\n\tmcmc;\nend;";			return sT;		}	}	/*.................................................................................................................*/	/*.................................................................................................................*/	boolean heterogeneous = false;	Class previousClass = null;//datatype=mixed(dna:1-100,protein:101-200)	int addDataTypes(Taxa taxa, StringBuffer buffer, StringBuffer partitionBuffer,  int totNumChars, Class dataSuperclass){		int numMatrices = getProject().getNumberCharMatrices( taxa, dataSuperclass);		for (int iM = 0; iM < numMatrices; iM++){			CharacterData data = getProject().getCharacterMatrix( taxa, iM, dataSuperclass);			if (data != null) { 				if (needsComma)					buffer.append(',');				needsComma = true;				if (data instanceof RNAData){					if (previousClass != null && previousClass != RNAData.class)						heterogeneous = true;					previousClass = RNAData.class;					buffer.append(" rna");				}				else if (data instanceof DNAData){					if (previousClass != null && previousClass != DNAData.class)						heterogeneous = true;					previousClass = DNAData.class;					buffer.append(" dna");				}				else if (data instanceof ProteinData){					if (previousClass != null && previousClass != ProteinData.class)						heterogeneous = true;					previousClass = ProteinData.class;					buffer.append(" protein");				}				else if (data instanceof CategoricalData){					if (previousClass != null && previousClass != CategoricalData.class)						heterogeneous = true;					previousClass = CategoricalData.class;					buffer.append(" standard");				}				else if (data instanceof ContinuousData){					if (previousClass != null && previousClass != ContinuousData.class)						heterogeneous = true;					previousClass = ContinuousData.class;					buffer.append(" continuous");				}				buffer.append(": " + (totNumChars+1) + "-" + (totNumChars + data.getNumChars()));				totNumChars += data.getNumChars();			}		}		return totNumChars;	}	/*.................................................................................................................*/	void composeForTaxon(Taxon[] components, int iTaxaBlock, StringBuffer buffer, boolean master, Class dataSuperclass){		if (!master && (components == null || components.length == 0)){ //no representatives in this taxa block; write gaps for those matrices			Taxa taxa = getProject().getTaxa(iTaxaBlock);			int numMatrices = getProject().getNumberCharMatrices( taxa, dataSuperclass);			for (int iM = 0; iM < numMatrices; iM++){				CharacterData data = getProject().getCharacterMatrix( taxa, iM, dataSuperclass);				if (data != null){					for (int ic=0; ic<data.getNumChars(); ic++)						buffer.append('-');				}			}			return;		}		Taxon first = components[0];  //this will tell taxa block;		Taxa taxa = first.getTaxa();		int numMatrices = getProject().getNumberCharMatrices( taxa, dataSuperclass);		for (int iM = 0; iM < numMatrices; iM++){			CharacterData data = getProject().getCharacterMatrix( taxa, iM, dataSuperclass);			if (data != null){				for (int ic=0; ic<data.getNumChars(); ic++){					int it = components[0].getNumber();  //Todo: this should fuse all components; at the moment it only takes it from the first taxon.					data.statesIntoNEXUSStringBuffer(ic, it, buffer);				}			} 		}	}	/*.................................................................................................................*/	void composeForMasterTaxon(Taxa masterTaxa, int it, Taxon[][] components, StringBuffer buffer, Class dataSuperclass, CommandRecord commandRec){		buffer.append(StringUtil.tokenize(masterTaxa.getTaxonName(it)) + "   ");		composeForTaxon(new Taxon[]{masterTaxa.getTaxon(it)}, -1, buffer, true, dataSuperclass);		if (components != null)			for (int iTaxaBlock = 0; iTaxaBlock<components.length; iTaxaBlock++){				if (iTaxaBlock != getProject().getTaxaNumber(masterTaxa))					composeForTaxon(components[iTaxaBlock], iTaxaBlock, buffer, false, dataSuperclass);			}	}	boolean needsComma= false;	/*	begin data;	   dimensions ntax=4 nchar=10;	   format datatype=dna gap=-;	   matrix	   taxon_1 AACGATTCGT	   taxon_2 AAGGAT--CA	   taxon_3 AACGACTCCT	   taxon_4 AAGGATTCCT	   ;	end;	/*.................................................................................................................*/	void fillMatrix(Taxa masterTaxa, StringBuffer buffer, Class dataSuperclass, CommandRecord commandRec){		if (associationTask == null){			associationTask = (AssociationSource)hireEmployee(commandRec, AssociationSource.class, "Source of taxon associations");		}		buffer.append("#NEXUS\n\nbegin data;\n");		StringBuffer[] taxaStrings = new StringBuffer[masterTaxa.getNumTaxa()];		StringBuffer dataTypesBuffer = new StringBuffer();		StringBuffer partitionBuffer = new StringBuffer();		for (int it=0; it< masterTaxa.getNumTaxa(); it++)			taxaStrings[it] = new StringBuffer(100);		Taxon[][][] components = new Taxon[masterTaxa.getNumTaxa()][getProject().getNumberTaxas()][];		needsComma = false;		heterogeneous = false;		previousClass = null;		int totNumChars =  addDataTypes(masterTaxa, dataTypesBuffer, partitionBuffer, 0, dataSuperclass);		if (associationTask != null){			for (int iTaxaBlock = 0; iTaxaBlock < getProject().getNumberTaxas(); iTaxaBlock++){				Taxa aTaxa = getProject().getTaxa(iTaxaBlock);				if (aTaxa != masterTaxa){					TaxaAssociation assoc = associationTask.getAssociation(masterTaxa, aTaxa, 0, commandRec);  //Todo: permit choice					if (assoc != null){						totNumChars = addDataTypes(assoc.getOtherTaxa(masterTaxa), dataTypesBuffer,partitionBuffer,  totNumChars, dataSuperclass);						for (int it = 0; it<masterTaxa.getNumTaxa(); it++){							components[it][iTaxaBlock] = assoc.getAssociates(masterTaxa.getTaxon(it));						}					}				}			}		}		buffer.append("dimensions ntax = " + masterTaxa.getNumTaxa() + " nchar = " + totNumChars + ";\n");		buffer.append("format datatype = ");		if (heterogeneous)			buffer.append("mixed(" + dataTypesBuffer + ")");		else if (previousClass == RNAData.class)			buffer.append("rna");		else if (previousClass == DNAData.class)			buffer.append("dna");		else if (previousClass == ProteinData.class)			buffer.append("protein");		else if (previousClass == CategoricalData.class)			buffer.append("standard");		else if (previousClass == ContinuousData.class)			buffer.append("continuous");					buffer.append(" gap = - missing =?;\nmatrix\n");		for (int it=0; it< masterTaxa.getNumTaxa(); it++)			composeForMasterTaxon(masterTaxa, it, components[it], taxaStrings[it], dataSuperclass, commandRec);		buffer.append('\n');		for (int i=0; i< masterTaxa.getNumTaxa(); i++){			buffer.append(taxaStrings[i]);			buffer.append('\n');			taxaStrings[i].setLength(0);		}		buffer.append("\n;\nend;\n");	}	/*.................................................................................................................*/	public void exportFile(MesquiteFile file, String arguments, CommandRecord commandRec) { //if file is null, consider whole project open to export		Arguments args = new Arguments(new Parser(arguments), true);		boolean usePrevious = args.parameterExists("usePrevious");		Taxa masterTaxa = (Taxa)getProject().chooseTaxa(containerOfModule(), "Select master block of taxa", false, commandRec);		StringBuffer buffer = new StringBuffer(500);		fillMatrix(masterTaxa, buffer, CategoricalState.class, commandRec);/*		if (!commandRec.scripting() && !usePrevious)			if (!getExportOptions())				return;*/		saveExportedFileWithExtension(buffer, arguments, "nex", commandRec);	}	/*.................................................................................................................*/	/** returns the version number at which this module was first released.  If 0, then no version number is claimed.  If a POSITIVE integer	 * then the number refers to the Mesquite version.  This should be used only by modules part of the core release of Mesquite.	 * If a NEGATIVE integer, then the number refers to the local version of the package, e.g. a third party package*/	public int getVersionOfFirstRelease(){		return 107;  	}	/*.................................................................................................................*/	public String getName() {		return "Fused Export NEXUS for MrBayes";	}	/*.................................................................................................................*/	/** returns an explanation of what the module does.*/	public String getExplanation() {		return "Exports NEXUS files for use by MrBayes with matrices fused." ;	}	/*.................................................................................................................*/}